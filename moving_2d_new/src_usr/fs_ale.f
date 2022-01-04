!----------------------------------------------------------------------
!     Author: Prabal Negi
!     Description: Routines for free surface mesh movement.
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      subroutine gen_mapping_mvb

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'MASS'          ! bm1, temporary
      include 'PARALLEL'

      include 'FS_ALE'

      real xsort(lx1*ly1*lz1*lelv)
      real xlocal_u(lx1*lelv)
      integer local_num(lx1*ly1*lz1*lelv) ! local numbering
      integer sortind(lx1*ly1*lz1*lelv)
      integer i,j,k,n,n1

      integer lelxx
      parameter (lelxx=500)
      real gxsort1(lx1*lelxx),gxsort2(lx1*lelxx)
      integer gind(lx1*lelxx),gind2(lx1*lelxx),gind3(lx1*lelxx)
      integer gl_num_u(lx1*lelxx)         ! unique global numbers

      integer nunq                        ! Unique local points
      integer n2,n3
      integer glnos
     
      real xlast,tol

      integer igl_running_sum       ! function
      integer iglsum

      real xmean,ymean,vol
      character cb*3
      real temp1(lx1,ly1,lz1,lelv)
      integer kx1,kx2,ky1,ky2,kz1,kz2,nx,ny,nz,iface,ie,ieg
      integer ix,iy,iz,nfaces
      real vlsum
      integer ifield, nxyz
      real temp0(lx1,ly1,lz1)

      call mntr_log(fs_id,fs_log,'Generating Global Mapping')

      n = lx1*ly1*lz1*nelv

!     Keep original mesh location      
      call copy3(xm0_fs,ym0_fs,zm0_fs,xm1,ym1,zm1,n)

!     Local sorting      
      call copy(xsort,ym1,n)
      call sort(xsort,sortind,n)

!     Local unique elements      
      xlast = -9999.0
      nunq  = 0
      tol   = 1.0e-12         ! This is arbitrary. 
                              ! Should use a better measure
      do i=1,n
        if (abs(xsort(i)-xlast).gt.tol) then
          nunq = nunq+1
          xlocal_u(nunq) = xsort(i)
          xlast = xsort(i)
        endif
        j = sortind(i)
        local_num(j) = nunq
      enddo  

      n2 = igl_running_sum(nunq)
      n3 = iglsum(nunq,1)

!     Global sorting                                          
      call rzero(gxsort1,lx1*lelxx)
      call copy(gxsort1(n2-nunq+1),xlocal_u,nunq)

      call gop(gxsort1,gxsort2,'+  ',n3)

      call copy(gxsort2,gxsort1,n3)       ! in principle we shouldn't
                                          ! need this

      call sort(gxsort1,gind,n3)

      xlast = -9999.0
      glnos = 0
      do i=1,n3
        if (abs(gxsort1(i)-xlast).gt.tol) then
          glnos = glnos+1
          xlast = gxsort1(i)
        endif
        gind2(i) = glnos
      enddo  

      call nekgsync()

      call i8zero(fs_gl_num,lx1*lelxx)
      do i=1,n3
        j = gind(i)
        gind3(j) = gind2(i)
      enddo

      do i=1,nunq
        j = n2-nunq+i
        gl_num_u(i) = gind3(j)
      enddo

!     Set global numbering in the whole field      
      xlast = -9999.0
      n1    = 0
      do i=1,n
        j = local_num(i)
        fs_gl_num(i) = gl_num_u(j)
      enddo  

!     Gather-Scatter Setup      
      call fs_setupds(fs_gs_handle,n,fs_gl_num)

      call rzero(fs_mask,n)

!     Multiplicity (Inverse) (vmult)
      call rzero(fs_vmult,n)
      nfaces = 2*ndim
      ifield = 1
      nxyz   = lx1*ly1*lz1
      nx     = lx1
      ny     = ly1
      nz     = lz1
      fs_nel = 0
      do ie=1,nelv
        call copy(temp0,xm1(1,1,1,ie),nxyz)
        call col2(temp0,bm1(1,1,1,ie),nxyz)
        xmean = vlsum(temp0,nxyz)
        vol   = vlsum(bm1(1,1,1,ie),nxyz)
        xmean = xmean/vol
        call copy(temp0,ym1(1,1,1,ie),nxyz)
        call col2(temp0,bm1(1,1,1,ie),nxyz)
        ymean = vlsum(temp0,nxyz)
        ymean = ymean/vol
        do iface=1,nfaces
          cb  = cbc(iface,ie,ifield)
          ieg = lglel(ie)
          call facind (kx1,kx2,ky1,ky2,kz1,kz2,nx,ny,nz,iface)

          fs_cbc(iface,ie) = '   '
!         Change these conditions for actual case            
          if (cb.eq.'O  '.and.(xmean.lt.0.0)) then
             fs_nel = fs_nel + 1
             fs_cbc(iface,ie) = 'INT'
             fs_elno(fs_nel)  = ie
             fs_iface(fs_nel) = iface
             do iz=kz1,kz2
             do iy=ky1,ky2
             do ix=kx1,kx2
               fs_vmult(ix,iy,iz,ie) = 1.0
               temp1(ix,iy,iz,ie) = ((ym1(ix,iy,iz,ie)-3.0)/2.0)**2 !ymean
               fs_mask(ix,iy,iz,ie)  = 1.0
             enddo
             enddo
             enddo
          endif                     
        enddo
      enddo      

      call fgslib_gs_op(fs_gs_handle,fs_vmult,1,1,0)  ! 1 ==> +
      call invcol1(fs_vmult,n)

      call fgslib_gs_op(fs_gs_handle,temp1,1,1,0)  ! 1 ==> +
      call col2(temp1,fs_vmult,n)

      call fs_symo_save       ! save SYM/O corner points

      call mntr_log(fs_id,fs_log,'Global Mapping Done.')

      return
      end subroutine gen_mapping_mvb
!----------------------------------------------------------------------
      subroutine fs_gen_damping

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'GEOM'

      include 'FS_ALE'

      real dummy
      common /scrcg/ dummy(lx1,ly1,lz1,lelt)

      real x,x0,mu
      integer i,n

      call mntr_log(fs_id,fs_log,'Generating Damping Function')


      n = lx1*ly1*lz1*nelv

!     Find interface position      
      call copy(dummy,xm1,n)
      call col2(dummy,fs_mask,n)
      call fgslib_gs_op(fs_gs_handle,dummy,1,1,0)     ! 1 ==> +
      call col2(dummy,fs_vmult,n)


!     Create damping function      
      mu = 0.80
      do i=1,n
        x0               = dummy(i,1,1,1)
        x                = xm1(i,1,1,1)
        if (abs(x-x0).lt.fs_ofst) then
          fs_damp(i,1,1,1) = 1.0
        else  
          fs_damp(i,1,1,1) = exp(-((x-x0-fs_ofst)/mu)**2)
        endif  
      enddo

      call col2(fs_damp,v1mask,n)

      ifto = .true.
      call outpost(fs_mask,fs_vmult,fs_damp,pr,fs_damp,'gln')


      return
      end subroutine fs_gen_damping        

!---------------------------------------------------------------------- 
      subroutine fs_setupds(gs_handle,ntot,glo_num)

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
!      include 'NONCON'

      integer gs_handle
      integer*8 glo_num(1)

      integer mid,mp,nekcomm,nekgroup,nekreal      
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      integer ntot
      real t0,t1
      real dnekclock          ! function

      t0 = dnekclock()

c     Initialize gather-scatter code
!      ntot      = lx1*ly1*lz1*nelv
      call fgslib_gs_setup(gs_handle,glo_num,ntot,nekcomm,mp)

c     call gs_chkr(glo_num)

      t1 = dnekclock() - t0
      if (nio.eq.0) then
         write(6,1) t1,gs_handle
    1    format('FS setupds time',1pe11.4,' seconds ',i3)
      endif
c
      return
      end
!-----------------------------------------------------------------------

      subroutine fs_mvmesh()

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'GEOM'
      include 'MVGEOM'
      include 'INPUT'
      include 'FS_ALE'

      integer i,ntot1

      if (.not.ifmvbd) return
      if (.not.ifusermv) return

      param(59) = 1.       ! ifdef = .true.

      ntot1 = lx1*ly1*lz1*nelv

!     Zero out everything except the free surface      
      call opcopy(wx,wy,wz,vx,vy,vz)

      if (fs_ifgrid) then
!       Directly get smooth normals/tangentials by global
!       approximation of the interface.      
        call fs_smoothmv(wx,wy,wz)
      else
!       Move mesh in the normal direction.
!       With possible tangential correction.
!       And projection to smooth field
        call fs_smooth_meshmv(wx,wy,wz)
      endif  

      call fs_int_project(wx,wy,wz)

!     Spatially damp out the mesh velocity      
      call col2(wx,fs_damp,ntot1)
      call col2(wy,fs_damp,ntot1)
      if (ndim.eq.3) then
        call col2(wz,fs_damp,ntot1)
      else  
!       We don't move mesh along Z
        call rzero(wz,ntot1)
      endif

!     Use Gordon Hall for mesh motion
      if (fs_ifgh) call fs_gh_correction(wx,wy,wz) 

      return
      end subroutine fs_mvmesh        

!-----------------------------------------------------------------------

      subroutine fs_mvmeshn(ux,uy,uz)
!     Only in 2D for now
!     Project velocities on to Normal directions

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'MASS'

      include 'FS_ALE'

      real ux(lx1,ly1,lz1,lelv)
      real uy(lx1,ly1,lz1,lelv)
      real uz(lx1,ly1,lz1,lelv)

      integer i,n,nface
      integer js1,js2,jf1,jf2,jskip1,jskip2,ifc,e
      integer j1,j2,j3,nxyz
      integer ifld

      real rnor,rtn1

      character cb*3

      real dummy1,dummy2,dummy3
      common /scrsf/ dummy1(lx1,ly1,lz1,lelt),
     $               dummy2(lx1,ly1,lz1,lelt),
     $               dummy3(lx1,ly1,lz1,lelt)

      integer nsave
      integer ies(2),ixs(2),iys(2)
      save ies,ixs,iys,nsave
      real wallvx(lx1*lelt),wallvy(lx1*lelt)
      real tol
      integer icalld
      save icalld
      data icalld /0/

      ifld  = 1
      nxyz  = lx1*ly1*lz1
      nface = 2*ndim

      do i = 1,fs_nsymo
        e  = fs_ie(i)
        j1 = fs_ix(i)
        j2 = fs_iy(i)
        j3 = fs_iz(i)
        wallvx(i) = ux(j1,j2,j3,e)
        wallvy(i) = 0.0
      enddo  
        

      do 200 e=1,nelv
      do 200 ifc=1,nface
        cb  = fs_cbc(ifc,e)
        if (cb.eq.'INT') then
          call facind2 (js1,jf1,jskip1,js2,jf2,jskip2,ifc)
          i = 0
          do 220 j2=js2,jf2,jskip2
          do 220 j1=js1,jf1,jskip1
             i = i + 1
!            normal component         
             rnor = ( ux(j1,j2,1,e)*unx(i,1,ifc,e) +
     $                uy(j1,j2,1,e)*uny(i,1,ifc,e) )
!            remove tangential component
             ux(j1,j2,1,e) = rnor*unx(i,1,ifc,e)
             uy(j1,j2,1,e) = rnor*uny(i,1,ifc,e)


  220      continue
        endif                 
  200 continue

      call dsavg(ux)
      call dsavg(uy)

      do i=1,fs_nsymo
        e  = fs_ie(i)
        j1 = fs_ix(i)
        j2 = fs_iy(i)
        j3 = fs_iz(i)

        ux(j1,j2,j3,e) = wallvx(i)
        uy(j1,j2,j3,e) = wallvy(i)
      enddo   


      return
      end subroutine fs_mvmeshn        
!---------------------------------------------------------------------- 

      subroutine fs_projt(ux,uy,uz,tdir)
!     Project velocities on to tangential directions

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'MASS'

      include 'FS_ALE'

      real ux(lx1,ly1,lz1,lelv)
      real uy(lx1,ly1,lz1,lelv)
      real uz(lx1,ly1,lz1,lelv)

      integer i,n,nface
      integer js1,js2,jf1,jf2,jskip1,jskip2,ifc,e
      integer j1,j2,j3,nxyz
      integer ifld

      real rnor,rtn1,rtn2
      integer tdir

      character cb*3

      real dummy1,dummy2,dummy3
      common /scrsf/ dummy1(lx1,ly1,lz1,lelt),
     $               dummy2(lx1,ly1,lz1,lelt),
     $               dummy3(lx1,ly1,lz1,lelt)

      integer nsave
      integer ies(lx1),ixs(lx1),iys(lx1)
      save ies,ixs,iys,nsave
      real wallvx(lx1),wallvy(lx1)
      real tol
      integer icalld
      save icalld
      data icalld /0/

      nxyz  = lx1*ly1*lz1
      nface = 2*ndim

      do i = 1,nsave
        wallvx(i) = 0.0
        wallvy(i) = 0.0
      enddo  
        

      do 200 e=1,nelv
        do 200 ifc=1,nface
          cb  = fs_cbc(ifc,e)
          if (cb.eq.'INT') then
            call facind2 (js1,jf1,jskip1,js2,jf2,jskip2,ifc)
            i = 0
            do 220 j2=js2,jf2,jskip2
            do 220 j1=js1,jf1,jskip1
               i = i + 1
!              normal component         
               rnor = ( ux(j1,j2,1,e)*unx(i,1,ifc,e) +
     $                  uy(j1,j2,1,e)*uny(i,1,ifc,e) )
!              tangential compnent            
               rtn1 = ( ux(j1,j2,1,e)*t1x(i,1,ifc,e) +
     $                  uy(j1,j2,1,e)*t1y(i,1,ifc,e) )
               rtn2 = ( ux(j1,j2,1,e)*t2x(i,1,ifc,e) +
     $                  uy(j1,j2,1,e)*t2y(i,1,ifc,e) )
!              project to tangential direction
               if (tdir.eq.1) then
                 ux(j1,j2,1,e) = rtn1*t1x(i,1,ifc,e)
                 uy(j1,j2,1,e) = rtn1*t1y(i,1,ifc,e)
                 if (ndim.eq.3) uz(j1,j2,1,e) = rtn1*t1z(i,1,ifc,e)
               elseif (tdir.eq.2) then
                 ux(j1,j2,1,e) = rtn2*t2x(i,1,ifc,e)
                 uy(j1,j2,1,e) = rtn2*t2y(i,1,ifc,e)
                 if (ndim.eq.3) uz(j1,j2,1,e) = rtn2*t2z(i,1,ifc,e)
               else
                 if (nio.eq.0) write(6,*) 'unknown tdir:', tdir
                 call exitt  
               endif
  220        continue
          endif                 
  200   continue

      call dsavg(ux)
      call dsavg(uy)

      do i=1,nsave
        e  = fs_ie(i)
        j1 = fs_ix(i)
        j2 = fs_iy(i)
        j3 = fs_iz(i)

        ux(j1,j2,j3,e) = wallvx(i)
        uy(j1,j2,j3,e) = wallvy(i)
      enddo   

      return
      end subroutine fs_projt 
!---------------------------------------------------------------------- 

      subroutine fs_int_project(wx,wy,wz) 

!     project surface values to the interior of the domain        

      implicit none

      include 'SIZE'
      include 'FS_ALE'

      real wx(lx1,ly1,lz1,lelv)
      real wy(lx1,ly1,lz1,lelv)
      real wz(lx1,ly1,lz1,lelv)

      integer ntot1

      ntot1 = lx1*ly1*lz1*nelv

      call col2(wx,fs_mask,ntot1)
      call col2(wy,fs_mask,ntot1)
      if (ndim.eq.3) call col2(wz,fs_mask,ntot1)

!     Extend the velocity to the interior of the domain
!     using the custom gather-scatter operator
      call fgslib_gs_op(fs_gs_handle,wx,1,1,0)  ! 1 ==> +
      call fgslib_gs_op(fs_gs_handle,wy,1,1,0)  ! 1 ==> +
      if (ndim.eq.3) call fgslib_gs_op(fs_gs_handle,wz,1,1,0)  ! 1 ==> +

!     Take care of multiplicity 
      call col2(wx,fs_vmult,ntot1)
      call col2(wy,fs_vmult,ntot1)
      if (ndim.eq.3) call col2(wz,fs_vmult,ntot1)

      return
      end subroutine fs_int_project
!----------------------------------------------------------------------

      subroutine fs_symo_save

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'FS_ALE'

      integer nface,nsave
      integer js1,js2,jf1,jf2,jskip1,jskip2,ifc,e
      integer j1,j2,j3
      integer ifld

      character cb*3
      real tol
      real dummy1,dummy2,dummy3
      common /scrsf/ dummy1(lx1,ly1,lz1,lelt),
     $               dummy2(lx1,ly1,lz1,lelt),
     $               dummy3(lx1,ly1,lz1,lelt)


!     Store corner point locations for SYM/O conditions
      ifld  = 1 
      tol   = 1.0e-14
      nsave = 0
      nface = 2*ndim
      fs_nsymo = 0
      call rzero(dummy1,lx1*ly1*lz1*nelv)
      do e=1,nelv
      do ifc=1,nface
        cb  = cbc(ifc,e,ifld)
!       Change these conditions for actual case            
        if (cb.eq.'SYM'.or.cb.eq.'O  ') then
          call facind2 (js1,jf1,jskip1,js2,jf2,jskip2,ifc)
          do j2=js2,jf2,jskip2
          do j1=js1,jf1,jskip1
            dummy1(j1,j2,1,e) = dummy1(j1,j2,1,e)+1.0
            if (abs(dummy1(j1,j2,1,e)-2.0).lt.tol) then
              nsave         = nsave + 1
              fs_ie(nsave)  = e
              fs_ix(nsave)  = j1
              fs_iy(nsave)  = j2
              fs_iz(nsave)  = 1
              fs_nsymo      = nsave
            endif
          enddo ! j2
          enddo ! j1
        endif                 
      enddo
      enddo

      return
      end subroutine fs_symo_save
!----------------------------------------------------------------------
      subroutine fs_gh_correction(wx,wy,wz)

      implicit none

      include 'SIZE'
      include 'INPUT'   ! cbc
      include 'TSTEP'

      include 'FS_ALE'

!     Mesh velocities      
      real wx(lx1,ly1,lz1,lelv)
      real wy(lx1,ly1,lz1,lelv)
      real wz(lx1,ly1,lz1,lelv)

      integer e,f,nfaces,ntot1,ifld
      real dti
      integer i,j,k

      character cb*3

!     mask for interior points
      real maski
      common /scrcg/ maski(lx1,ly1,lz1,lelt)


      call mntr_log(fs_id,fs_log,'Gordon Hall Correction')

      ntot1  = lx1*ly1*lz1*nelv
      nfaces = 2*ndim
      ifld   = 1

!     create a mask for interior points.
      call rzero(maski,ntot1)

      do e=1,nelv     
      do f=1,nfaces   
         cb = cbc(f,e,ifld)
!        fill up edges with ones            
         call facev(maski,e,f,1.0,nx1,ny1,nz1)
!        Put ones on the interface 
!         if (cb.eq.'O  ') call facev(maski,e,f,1.0,nx1,ny1,nz1)
!        Don't move Periodic points         
         if (cb.eq.'P  ') call facev(maski,e,f,0.0,nx1,ny1,nz1)
!        Keeping wall fixed            
         if (cb.eq.'W  ') call facev(maski,e,f,0.0,nx1,ny1,nz1)
      enddo
      enddo

!     Only keep the edge values            
      call col2(wx,maski,ntot1)
      call col2(wy,maski,ntot1)
      if (ndim.eq.3) call col2(wz,maski,ntot1)

!     At first order, the displacement is just wx*dt, etc
!     We can do a better approximation if this works.
!     (Although better approximation is probably not needed)      
      call opcmult(wx,wy,wz,dt)           ! This is now displacement
      
!     Project displacement into the interior of the domain      
      call fs_gordonhall(wx,wy,wz)

      dti = 1.0/dt
      call opcmult(wx,wy,wz,dti)           ! This is now velocities

      return
      end subroutine fs_gh_correction        
!----------------------------------------------------------------------       
      subroutine fs_gordonhall(mvx,mvy,mvz) 
!     fix up geometry irregularities
!     routine has been modified from fix_geom
!     to move internal points of the domain.

      implicit none

      include 'SIZE'
      include 'INPUT'   ! cbc
      include 'GEOM'
      include 'TSTEP'   ! ifield
      include 'WZ'

      integer lt
      parameter (lt = lx1*ly1*lz1)
      real xb,yb,zb,tmsk,tmlt,w1,w2
      common /scrns/ xb(lt,lelt),yb(lt,lelt),zb(lt,lelt),
     $               tmsk(lt,lelt),tmlt(lt,lelt),w1(lt),w2(lt)

      real xtmp,ytmp,ztmp
      common /scrsf/ xtmp(lx1,ly1,lz1,lelt),
     $               ytmp(lx1,ly1,lz1,lelt),
     $               ztmp(lx1,ly1,lz1,lelt)


      real mvx(lt,lelt),mvy(lt,lelt),mvz(lt,lelt) ! displaced boundary points
      
      integer i,e,kpass,n,f,nxyz,nfaces
      character*3 cb

      real s,xm,ym,zm,xx,yy,zz
      real glamax
      real glmax

      n      = nx1*ny1*nz1*nelt
      nxyz   = nx1*ny1*nz1
      nfaces = 2*ndim
      ifield = 1                   ! velocity field
!      if (ifheat) ifield = 2       ! temperature field

!     Use the temporary arrays for x,y,z      
      call opcopy(xtmp,ytmp,ztmp,xm1,ym1,zm1) 

      call rone  (tmlt,n)
      call dssum (tmlt,nx1,ny1,nz1)  ! denominator

      call rone  (tmsk,n)
      do e=1,nelfld(ifield)      ! fill mask where bc is periodic
      do f=1,nfaces              ! so we don't translate periodic bcs (z only)
         cb =cbc(f,e,ifield)
         if (cb.eq.'P  ') call facev (tmsk,e,f,0.0,nx1,ny1,nz1)
!        Keeping wall fixed            
         if (cb.eq.'W  ') call facev (tmsk,e,f,0.0,nx1,ny1,nz1)
      enddo
      enddo

      do kpass = 1,ndim+1     ! extra pass is just to test convergence

         call copy(xb,xtmp,n)
         call copy(yb,ytmp,n)
         if (ndim.eq.3) call copy(zb,ztmp,n)
           
         call dssum(xb,nx1,ny1,nz1)
         call dssum(yb,nx1,ny1,nz1)
         if (ndim.eq.3) call dssum(zb,nx1,ny1,nz1)

         xm = 0.
         ym = 0.
         zm = 0.

         do e=1,nelt
            do i=1,nxyz                       ! compute averages of geometry
               s     = 1./tmlt(i,e)
               xb(i,e) = s*xb(i,e)
               yb(i,e) = s*yb(i,e)
               if (ndim.eq.3) zb(i,e) = s*zb(i,e)

!              if there are tears between elemnts   
               xb(i,e) = xb(i,e) - xtmp(i,1,1,e)   ! local displacements
               yb(i,e) = yb(i,e) - ytmp(i,1,1,e)
               if (ndim.eq.3) zb(i,e) = zb(i,e) - ztmp(i,1,1,e)

               xb(i,e) = xb(i,e)*tmsk(i,e)
               yb(i,e) = yb(i,e)*tmsk(i,e)
               if (ndim.eq.3) zb(i,e) = zb(i,e)*tmsk(i,e)

!              add in the boundary displacements in the first pass   
               if (kpass.le.ndim) then
                 xb(i,e)=xb(i,e)+mvx(i,e)
                 yb(i,e)=yb(i,e)+mvy(i,e)
                 if (ndim.eq.3) zb(i,e)=zb(i,e)+mvz(i,e)
               endif   

               xm = max(xm,abs(xb(i,e)))
               ym = max(ym,abs(yb(i,e)))
               if (ndim.eq.3) zm = max(zm,abs(zb(i,e)))
            enddo

            if (kpass.le.ndim) then
               call gh_face_extend(xb(1,e),zgm1,nx1,kpass,w1,w2)
               call gh_face_extend(yb(1,e),zgm1,nx1,kpass,w1,w2)
               if (ndim.eq.3)
     $           call gh_face_extend(zb(1,e),zgm1,nx1,kpass,w1,w2)
            endif ! kpass.le.ndim

         enddo

         if (kpass.le.ndim) then
           call add2(xtmp,xb,n)
           call add2(ytmp,yb,n)
           if (ndim.eq.3) call add2(ztmp,zb,n)

           call sub2(mvx,xb,n)
           call sub2(mvy,yb,n)
           if (ndim.eq.3) call sub2(mvz,zb,n)
         endif
        
         xx = glamax(xb,n)
         yy = glamax(yb,n)
         if (ndim.eq.3) zz = glamax(zb,n)

         xm = glmax(xm,1)
         ym = glmax(ym,1)
         if (ndim.eq.3) zm = glmax(zm,1)

!         if (nio.eq.0) write(6,1) xm,ym,zm,xx,yy,zz,kpass
!    1    format(1p6e12.4,' xyz repair',i2)

      enddo

      call opcopy(mvx,mvy,mvz,xtmp,ytmp,ztmp)
      call opsub2(mvx,mvy,mvz,xm1,ym1,zm1)
!      call opcopy(xm1,ym1,zm1,xtmp,ytmp,ztmp)
      param(59) = 1.       ! ifdef = .true.
      
      return
      end subroutine fs_gordonhall
c-----------------------------------------------------------------------

