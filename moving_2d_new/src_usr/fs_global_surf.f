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

      integer e,ifc,n,ne,nfaces

      integer kx1,kx2,ky1,ky2,kz1,kz2
      integer ix,iy,iz,ia

      character cb*3

      call rzero(zx_fs,lxfs)
      call rzero(zy_fs,lxfs)
      call rzero(wx_fs,lxfs)
      call rzero(wy_fs,lxfs)

      call zwgll (zx_fs,wx_fs,lxfs)
      call zwgll (zy_fs,wy_fs,lyfs)
      
      if (ndim.eq.3) then
        do ix=1,lxfs
        do iy=1,lyfs
          w2_fs(ix,iy) = wx_fs(ix)*wy_fs(iy)
        enddo
        enddo
      else
        call copy(w2_fs,wx_fs,lxfs)
      endif 

      if (ndim.eq.3) then
        do e=1,1
        do iy=1,lyfs
        do ix=1,lxfs
          xg_fs(ix,iy)=zx_fs(ix)
          yg_fs(ix,iy)=zy_fs(iy)
          fld_fs(ix,iy)=zx_fs(ix)*zx_fs(ix)
        enddo
        enddo
        enddo
      endif   

!     Get the surface x,y,z
      nfaces = 2*ndim
      ne     = 0              ! total number of interface elements
      do e=1,nelv
      do ifc=1,nfaces
        cb  = fs_cbc(ifc,e)
        if (cb.eq.'INT') then
          ia = 0
          ne = ne+1
          call facind (kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,ifc)
          do iz=kz1,kz2
          do iy=ky1,ky2
          do ix=kx1,kx2
            ia = ia + 1
            xm1_fs(ia,1,ne) = xm1(ix,iy,iz,e)
            ym1_fs(ia,1,ne) = ym1(ix,iy,iz,e)
            if (ndim.eq.3) zm1_fs(ia,1,ne) = zm1(ix,iy,iz,e)
          enddo
          enddo
          enddo
        endif
      enddo
      enddo  



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
      subroutine temporary_1

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'GFLDR'
      include 'PARALLEL'

      include 'FS_ALE'

      integer i

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

!     testing
!     initialize interpolation tool using global mesh
      nxf   = 2*lxfs
      nyf   = 2*lyfs
      nzf   = 2*1
      nhash = 1*lxfs*lyfs*1
      nmax  = 128
!     We do all Global calculations on nid 0
      if (nid.eq.0) then
        nels  = 1
      else
        nels  = 0
      endif  

      ldim2 = 2

      nxf   = 1
      nyf   = 1
      nzf   = 1

      lzfs  = 1

!     Interpolation handle for Global mesh.      
      call fgslib_findpts_setup(inth_test,nekcomm,np,ldim2,
     &                          xg_fs,yg_fs,yg_fs,lxfs,lyfs,lzfs,
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
         call fgslib_findpts(intgh_fs,
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
      call fgslib_findpts_eval(intgh_fs,
     &                         fldout,1,
     &                         grcode,1,
     &                         gproc,1,
     &                         gelid,1,
     &                         grst,ldim2,nintp,
     &                         fld_fs)

      write(6,*) fldout


      return
      end subroutine temporary_1        
