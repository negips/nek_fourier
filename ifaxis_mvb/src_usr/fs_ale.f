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

      
      integer*8 fs_gl_num(lx1*ly1*lz1*lelv)  ! final global numbers
      integer fs_gs_handle    ! free surface gather-scatter handle
      real fs_vmult(lx1,ly1,lz1,lelv)
      common /fs_gsh/ fs_gs_handle, fs_gl_num, fs_vmult

      real xmean,ymean,vol
      character cb*3
      real tmp1(lx1,ly1,lz1,lelv)
      integer kx1,kx2,ky1,ky2,kz1,kz2,nx,ny,nz,iface,ie,ieg
      integer ix,iy,iz,nfaces
      real vlsum
      integer ifield, nxyz
      real tmp0(lx1,ly1,lz1)


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

!      if (nid.eq.0) write(6,*) 'Global Numbering'
!      do i=1,nunq
!        if (nid.eq.0) write(6,*) nid,xlocal_u(i),gl_num_u(i)
!      enddo  
!      call nekgsync()
!      do i=1,nunq
!        if (nid.eq.1) write(6,*) nid,xlocal_u(i),gl_num_u(i)
!      enddo  
!      call nekgsync()

!     Set global numbering in the whole field      
      xlast = -9999.0
      n1    = 0
      do i=1,n
        j = local_num(i)
        fs_gl_num(i) = gl_num_u(j)
        vx(i,1,1,1) = fs_gl_num(i)+0.0
      enddo  

!     Gather-Scatter Setup      
      call fs_setupds(fs_gs_handle,n,fs_gl_num)

!     Multiplicity (vmult)
      call rzero(fs_vmult,n)
      nfaces = 2*ndim
      ifield = 1
      nxyz   = lx1*ly1*lz1
      nx     = lx1
      ny     = ly1
      nz     = lz1
      do ie=1,nelv
        call copy(tmp0,xm1(1,1,1,ie),nxyz)
        call col2(tmp0,bm1(1,1,1,ie),nxyz)
        xmean = vlsum(tmp0,nxyz)
        vol   = vlsum(bm1(1,1,1,ie),nxyz)
        xmean = xmean/vol
        call copy(tmp0,ym1(1,1,1,ie),nxyz)
        call col2(tmp0,bm1(1,1,1,ie),nxyz)
        ymean = vlsum(tmp0,nxyz)
        ymean = ymean/vol
        do iface=1,nfaces
          cb  = cbc(iface,ie,ifield)
          ieg = lglel(ie)
          call facind (kx1,kx2,ky1,ky2,kz1,kz2,nx,ny,nz,iface)

!         Change these conditions for actual case            
          if (cb.eq.'P  '.and.(xmean.lt.0.0)) then
             do iz=kz1,kz2
             do iy=ky1,ky2
             do ix=kx1,kx2
               fs_vmult(ix,iy,iz,ie) = 1.0
               tmp1(ix,iy,iz,ie) = ((ym1(ix,iy,iz,ie)-3.0)/2.0)**2 !ymean
             enddo
             enddo
             enddo
          endif                     
        enddo
      enddo      

      call fgslib_gs_op(fs_gs_handle,fs_vmult,1,1,0)  ! 1 ==> +
      call invcol1(fs_vmult,n)
      call copy(vy,fs_vmult,n)

      call fgslib_gs_op(fs_gs_handle,tmp1,1,1,0)  ! 1 ==> +
      call col2(tmp1,fs_vmult,n)
      call copy(t,tmp1,n)

      ifto = .true.
      call outpost(vx,vy,vz,pr,t,'gln')

      call exitt

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
c-----------------------------------------------------------------------












