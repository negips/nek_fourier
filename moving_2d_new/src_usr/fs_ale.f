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
      real xlocal_u(lx1*nelv)
      integer local_num(lx1*ly1*lz1*lelv) ! local numbering
      integer sortind(lx1*ly1*lz1*lelv)
      integer i,j,k,n,n1

      integer nelxx
      parameter (nelxx=50)
      real gxsort1(lx1*nelxx),gxsort2(lx1*nelxx)
      integer gind(lx1*nelxx),gind2(lx1*nelxx),gind3(lx1*nelxx)
      integer gl_num_u(lx1*nelxx)         ! unique global numbers

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


      n = lx1*ly1*lz1*nelv

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

!      write(6,*) 'nid,nunq', nid,nunq
!      if (nid.eq.0) write(6,*) 'xlocal_u', (xlocal_u(i), i=1,nunq)

      n2 = igl_running_sum(nunq)
      n3 = iglsum(nunq,1)

!     Global sorting                                          
      call rzero(gxsort1,lx1*nelxx)
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

      call i8zero(fs_gl_num,lx1*nelxx)
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

!     Multiplicity (vmult)
      call rzero(fs_vmult,n)
      nfaces = 2*ndim
      ifield = 1
      nxyz   = lx1*ly1*lz1
      nx     = lx1
      ny     = ly1
      nz     = lz1
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
             fs_cbc(iface,ie) = 'INT'
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

      real x,x0,mu
      integer i,n

      n = lx1*ly1*lz1*nelv

!     Create damping function      
      x0 = -1.0
      mu = 0.80
      do i=1,n
        x                = xm1(i,1,1,1)
        fs_damp(i,1,1,1) = exp(-((x-x0)/mu)**2)
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

      ntot1 = lx1*ly1*lz1*nelv

!     Zero out everything except the free surface      
      call opcopy(wx,wy,wz,vx,vy,vz)
!      call fs_mvmeshn(wx,wy,wz)
!      call col2(wy,v2mask,ntot1)          ! normal vel at 'SYM' = 0.0

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

!     Spatially damp out the mesh velocity      
      call col2(wx,fs_damp,ntot1)
      call col2(wy,fs_damp,ntot1)
      if (ndim.eq.3) then
        call col2(wz,fs_damp,ntot1)
      else  
!       We don't move mesh along Z
        call rzero(wz,ntot1)
      endif  
   
!     prabal      
!      call rzero3(wx,wy,wz,ntot1) 

      return
      end subroutine fs_mvmesh        

!-----------------------------------------------------------------------

      subroutine fs_mvmeshn(ux,uy,uz)
!     Only in 2D for now

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
      integer j1,j2,nxyz
      integer ifld

      real rnor,rtn1

      character cb*3
c      COMMON /CTMP1/ DUMMY1(LCTMP1)
c      COMMON /CTMP0/ DUMMY0(LCTMP0)
c
c      COMMON /SCRNS/ DUMMY2(LX1,LY1,LZ1,LELT,7)
c      COMMON /SCRUZ/ DUMMY3(LX1,LY1,LZ1,LELT,4)
c      COMMON /SCREV/ DUMMY4(LX1,LY1,LZ1,LELT,2)
c      COMMON /SCRVH/ DUMMY5(LX1,LY1,LZ1,LELT,2)
c      COMMON /SCRMG/ DUMMY6(LX1,LY1,LZ1,LELT,4)
c      COMMON /SCRCH/ DUMMY7(LX1,LY1,LZ1,LELT,2)
c      COMMON /SCRSF/ DUMMY8(LX1,LY1,LZ1,LELT,3)
c      COMMON /SCRCG/ DUMM10(LX1,LY1,LZ1,LELT,1)

      real dummy1,dummy2,dummy3
      common /scrsf/ dummy1(lx1,ly1,lz1,lelt),
     $               dummy2(lx1,ly1,lz1,lelt),
     $               dummy3(lx1,ly1,lz1,lelt)

      integer nsave
      integer ies(2),ixs(2),iys(2)
      save ies,ixs,iys,nsave
      real wallvx(2)
      real tol
      integer icalld
      save icalld
      data icalld /0/

      ifld  = 1
      nxyz  = lx1*ly1*lz1
      nface = 2*ndim

      tol   = 1.0e-14
      if (icalld.eq.0) then
        nsave = 0
        call rzero(dummy1,nxyz*nelv)
        do 100 e=1,nelv
          do 100 ifc=1,nface
            cb  = cbc(ifc,e,ifld)
!           Change these conditions for actual case            
            if (cb.eq.'SYM'.or.cb.eq.'O  ') then
              call facind2 (js1,jf1,jskip1,js2,jf2,jskip2,ifc)
              do 120 j2=js2,jf2,jskip2
              do 120 j1=js1,jf1,jskip1
                dummy1(j1,j2,1,e) = dummy1(j1,j2,1,e)+1.0
                if (abs(dummy1(j1,j2,1,e)-2.0).lt.tol) then
                  nsave = nsave + 1
                  ies(nsave)    = e
                  ixs(nsave)    = j1
                  iys(nsave)    = j2
                endif
  120         continue
            endif                 
  100     continue
        icalld = icalld + 1    
      endif       ! icalld

      if (nsave.eq.2) then
        wallvx(1) = ux(ixs(1),iys(1),1,ies(1))
        wallvx(2) = ux(ixs(2),iys(2),1,ies(2))
      else
        if (nio.eq.0) write(6,*) 'Did not find 2 corner points'
        call exitt
      endif  

      do 200 e=1,nelv
        do 200 ifc=1,nface
          cb  = fs_cbc(ifc,e)
          if (cb.eq.'INT') then
            call facind2 (js1,jf1,jskip1,js2,jf2,jskip2,ifc)
            do 220 j2=js2,jf2,jskip2
            do 220 j1=js1,jf1,jskip1
!              normal component         
               rnor = ( ux(j1,j2,1,e)*vnx(j1,j2,1,e) +
     $                  uy(j1,j2,1,e)*vny(j1,j2,1,e) )
!              tangential compnent            
               rtn1 = ( ux(j1,j2,1,e)*v1x(j1,j2,1,e) +
     $                  uy(j1,j2,1,e)*v1y(j1,j2,1,e) )
!              remove tangential component
               ux(j1,j2,1,e) = ux(j1,j2,1,e) - rtn1*v1x(j1,j2,1,e)
               uy(j1,j2,1,e) = uy(j1,j2,1,e) - rtn1*v1y(j1,j2,1,e)
!               ux(j1,j2,1,e) = rnor*vnx(j1,j2,1,e)
!               uy(j1,j2,1,e) = rnor*vny(j1,j2,1,e)
              
  220        continue
          endif                 
  200   continue

      call dsavg(ux)
      call dsavg(uy)

      ux(ixs(1),iys(1),1,ies(1)) = wallvx(1) 
      ux(ixs(2),iys(2),1,ies(2)) = wallvx(2)
      uy(ixs(1),iys(1),1,ies(1)) = 0.0 
      uy(ixs(2),iys(2),1,ies(2)) = 0.0

      return
      end subroutine fs_mvmeshn        
!---------------------------------------------------------------------- 




