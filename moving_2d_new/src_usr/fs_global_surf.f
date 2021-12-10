!----------------------------------------------------------------------
!     Author: Prabal Negi
!     Description: Routines for Generating global basis function.
!                  for surface representation
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      subroutine gen_global_basis

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'MASS'          ! bm1, temporary
      include 'PARALLEL'

      include 'FS_ALE'

      include 'GFLDR'

      integer lxfs,lyfs
      parameter (lxfs = 50)
      parameter (lyfs = 3)

      integer i,j,e,n

      real xg_fs(lxfs),wx_fs(lxfs)
      real yg_fs(lyfs),wy_fs(lyfs)
      real w2_fs(lxfs,lyfs)

      real xm1_fs(lxfs,lyfs,1)
      real ym1_fs(lxfs,lyfs,1)
      real fld_fs(lxfs,lyfs,1)

!     testing interpolation
      integer nxf,nyf,nzf
      integer nhash,nmax
      integer ldim2
      integer lzfs
      real zg_fs(1)
      real xin,yin,zin,fldout

      integer inth_test       ! handle 

      integer nintp           ! no of interpolation points
      integer*8 nfail
      integer*8 nfail_sum
      integer*8 i8glsum

      real toldist
      
      integer nidd,npp,nekcomm,nekgroup,nekreal
      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      call rzero(xg_fs,lxfs)
      call rzero(yg_fs,lxfs)
      call rzero(wx_fs,lxfs)
      call rzero(wy_fs,lxfs)

      ndim = 3
      call zwgll (xg_fs,wx_fs,lxfs)
      if (ndim.eq.3) call zwgll (yg_fs,wy_fs,lyfs)
      
      if (ndim.eq.3) then
        do i=1,lxfs
        do j=1,lyfs
          w2_fs(i,j) = wx_fs(i)*wy_fs(j)
        enddo
        enddo
      else
        call copy(w2_fs,wx_fs,lxfs)
      endif 

      if (ndim.eq.3) then
        do e=1,1
        do j=1,lyfs
        do i=1,lxfs
          xm1_fs(i,j,1)=xg_fs(i)
          ym1_fs(i,j,1)=yg_fs(j)
          fld_fs(i,j,1)=xg_fs(i)*xg_fs(i)
        enddo
        enddo
        enddo
      endif   
      
      ndim = 2

!     testing
!     initialize interpolation tool using source mesh
      nxf   = 2*lxfs
      nyf   = 2*lyfs
      nzf   = 2*1
      nhash = 1*lxfs*lyfs*1
      nmax  = 128
      nels  = 1
      ldim2 = 2

      nxf   = 1
      nyf   = 1
      nzf   = 1

      lzfs  = 1

      call fgslib_findpts_setup(inth_test,nekcomm,np,ldim2,
     &                          xm1_fs,ym1_fs,zg_fs,lxfs,lyfs,lzfs,
     &                          nels,nxf,nyf,nzf,bb_t,
     &                          nhash,nhash,nmax,tol)
    
!      if (iffpts) then ! locate points (iel,iproc,r,s,t)
         nfail = 0
         toldist = 5e-6
         toldist = 5e-14

!                    
         nintp  = 1
         xin    = 0.15
         yin    = 0.5
         zin    = 0.
         call fgslib_findpts(inth_test,
     &                       grcode,1,
     &                       gproc,1,
     &                       gelid,1,
     &                       grst,ldim2,
     &                       gdist,1,
     &                       xin,1,
     &                       yin,1,
     &                       zin,1,nintp)

         do i=1,nintp
            if(grcode(i).eq.1 .and. sqrt(gdist(i)).gt.toldist)
     &        nfail = nfail + 1
            if(grcode(i).eq.2) nfail = nfail + 1
         enddo

         nfail_sum = i8glsum(nfail,1)
         if(nfail_sum.gt.0) then
           if(nio.eq.0) write(6,*)
     &       ' WARNING: Unable to find all mesh points in source fld ',
     &       nfail_sum
         endif
!      endif

!      write(6,*) 'ids', grcode(1),gproc(1),gelid(1),(grst(i),i=1,3)

      ! evaluate inut field at given points
      call fgslib_findpts_eval(inth_test,
     &                         fldout,1,
     &                         grcode,1,
     &                         gproc,1,
     &                         gelid,1,
     &                         grst,ldim2,nintp,
     &                         fld_fs)

      write(6,*) fldout


      call exitt

      return
      end subroutine gen_global_basis
!----------------------------------------------------------------------
      subroutine fs_basis_setupds(gs_handle,ntot,glo_num)

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

      t1 = dnekclock() - t0
      if (nio.eq.0) then
         write(6,1) t1,gs_handle
    1    format('FS Glob. Surf. setupds time',1pe11.4,' seconds ',i3)
      endif
c
      return
      end subroutine fs_basis_setupds
!-----------------------------------------------------------------------





