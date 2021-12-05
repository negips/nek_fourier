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

      include 'TEST'

      real xsort(lx1*ly1*lz1*lelv)
      real xlocal_u(lx1*nelv)
      integer local_num(lx1*ly1*lz1*lelv) ! local numbering
      integer sortind(lx1*ly1*lz1*lelv)
      integer i,j,k,n,n1

      integer nelxx
      parameter (nelxx=10)
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

      real x,x0,mu

      n = lx1*ly1*lz1*nelv

!     Local sorting      
      call copy(xsort,ym1,n)
      call sort(xsort,sortind,n)

!     Local unique elements      
      xlast = -9999.0
      nunq  = 0
      tol   = 1.0e-08         ! This is arbitrary. 
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

!      write(6,*) 'nid,glnos', nid,glnos
!      do i=1,n3
!        if (nid.eq.0) write(6,*) gxsort1(i),gxsort2(i),gind(i),gind2(i)
!      enddo  

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
!        vx(i,1,1,1) = fs_gl_num(i)+0.0
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

!         Change these conditions for actual case            
          if (cb.eq.'O  '.and.(xmean.lt.0.0)) then
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

      ifto = .true.
      do i=1,n
        tmp1(i,1,1,1) = fs_gl_num(i)+0.0
        tmp2(i,1,1,1) = fs_vmult(i,1,1,1)+0.0
      enddo  


!     Create damping function      
      x0 = -1.0
      mu = 0.80
      do i=1,n
        x              = xm1(i,1,1,1)
        fs_damp(i,1,1,1) = exp(-((x-x0)/mu)**2)
      enddo

      call col2(fs_damp,v1mask,n)

      call outpost(fs_mask,fs_vmult,vz,pr,fs_damp,'gln')

      return
      end subroutine gen_mapping_mvb

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
      include 'FS_ALE'

      integer i,ntot1

      ntot1 = lx1*ly1*lz1*nelv

!     Zero out everything except the free surface      
      call opcopy(wx,wy,wz,vx,vy,vz)
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
      if (ndim.eq.3) call col2(wz,fs_damp,ntot1)
     

      return
      end subroutine fs_mvmesh        

!-----------------------------------------------------------------------








