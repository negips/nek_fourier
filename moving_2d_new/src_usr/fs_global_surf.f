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
      include 'WZ'

      include 'FS_ALE'

      include 'GFLDR'

      integer e,ifc,n,ne,nfaces

      integer kx1,kx2,ky1,ky2,kz1,kz2
      integer ix,iy,iz,ia,ii

      character cb*3
      real xm1_min,xm1_max,ym1_min,ym1_max
      real glmin,glmax

      call rzero(zx_fs,lxfs)
      call rzero(zy_fs,lxfs)
      call rzero(wx_fs,lxfs)
      call rzero(wy_fs,lxfs)

      call zwgll (zx_fs,wx_fs,lxfs)
      call zwgll (zy_fs,wy_fs,lyfs)

      xm1_min = glmin(xm1,lx1*ly1*lz1*nelv)
      xm1_max = glmax(xm1,lx1*ly1*lz1*nelv)
      ym1_min = glmin(ym1,lx1*ly1*lz1*nelv)
      ym1_max = glmax(ym1,lx1*ly1*lz1*nelv)
    
      if (ndim.eq.3) then
        do ix=1,lxfs
        do iy=1,lyfs
          w2_fs(ix,iy) = wx_fs(ix)*wy_fs(iy)
        enddo
        enddo
      else
        call copy(w2_fs,wx_fs,lxfs)
      endif 

      call rzero(xg_fs,lxfs*lyfs*2)
      call rzero(yg_fs,lxfs*lyfs*2)
      call rzero(zg_fs,lxfs*lyfs*2)

      do iy=1,lyfs
        do ix=1,lxfs
          xg_fs(ix,iy,1)=(zx_fs(ix)+1.0)*(xm1_max-xm1_min)/2.0 + xm1_min
          yg_fs(ix,iy,1)=(zy_fs(iy)+1.0)*(ym1_max-ym1_min)/2.0 + ym1_min
          if (ndim.eq.2) xg_fs(ix,iy,1) = zx_fs(ix)
          fld_fs(ix,iy,1)=yg_fs(ix,iy,1)**2
        enddo
      enddo


!     Get the surface x,y,z
      nfaces = 2*ndim
      ne     = 0              ! total number of interface elements
      ii     = 0
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
            ii = ii + 1
            xm1_fs(ii,1,1) = xm1(ix,iy,iz,e)
            ym1_fs(ii,1,1) = ym1(ix,iy,iz,e)
            if (ndim.eq.3) zm1_fs(ia,1,ne) = zm1(ix,iy,iz,e)
            if (ndim.eq.2) xm1_fs(ii,1,ne) = zgm1(iy,1)
          enddo
          enddo
          enddo
        endif
      enddo
      enddo  

!!     In 2d our interface is along the y-dir      
!      if (ndim.eq.2) then
!        do e=1,ne
!        do ix=1,lx1
!        do iy=1,ly1
!          ym1_fs(ix,iy,e)=ym1_fs(ix,1,e)
!          xm1_fs(ix,iy,e)=zgm1(iy,1)
!        enddo
!        enddo
!        enddo
!      endif


      return
      end subroutine gen_global_basis
!----------------------------------------------------------------------
      subroutine fs_intp_setup

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

      integer nintp           ! no of interpolation points
      
      integer nidd,npp,nekcomm,nekgroup,nekreal
      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

!     initialize interpolation tool using global mesh
      nxf   = 2*lx1
      nyf   = 2*ly1
      nzf   = 2*1
      nhash = fs_nel*lx1*ly1
      if (nhash.eq.0) then
        nhash = lx1*ly1
      endif  
      nmax  = 128
!     We do all Global calculations on nid 0
      if (nid.eq.0) then
        nels  = 1
      else
        nels  = 0
      endif  

      ldim2 = 2

!     Interpolation handle for Global mesh.      
      call fgslib_findpts_setup(intgh_fs,nekcomm,np,ldim2,
     &                          xg_fs,yg_fs,zg_fs,lxfs,lyfs,lzfs,
     &                          nels,nxf,nyf,nzf,bb_t,
     &                          nhash,nhash,nmax,tol)


!!     initialize interpolation tool using local sem mesh
!      nxf   = 2*lx1
!      nyf   = 2*ly1
!      nzf   = 2*1
!      nhash = fs_nel*lx1*ly1
!      nmax  = 128
!!     We do all Global calculations on nid 0
!      nels  = fs_nel
!
!      ldim2 = 2
!
!      
!!     Interpolation handle for SEM surface mesh. 
!      call fgslib_findpts_setup(intlh_fs,nekcomm,np,ldim2,
!     &                          xm1_fs,ym1_fs,zm1_fs,lx1,ly1,1,
!     &                          nels,nxf,nyf,nzf,bb_t,
!     &                          nhash,nhash,nmax,tol)

      if (nio.eq.0) write(6,*) 'FS: Interpolation Setup: Done'

      return
      end subroutine fs_intp_setup
!---------------------------------------------------------------------- 
      subroutine fs_getpts

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
      real xin,yin,zin,fldout

      integer nintp           ! no of interpolation points
      integer*8 nfail
      integer*8 nfail_sum
      integer*8 i8glsum

      real toldist

      integer ix,iy

      nfail = 0
      toldist = 5e-14
!
      if (ndim.eq.2) then      
        nintp  = fs_nel*lx1
      else
        nintp  = fs_nel*lx1*ly1
      endif
      ldim2  = 2

      call fgslib_findpts(intgh_fs,
     &                    grcode,1,
     &                    gproc,1,
     &                    gelid,1,
     &                    grst,ldim2,
     &                    gdist,1,
     &                    xm1_fs,1,
     &                    ym1_fs,1,
     &                    zm1_fs,1,nintp)


      do i=1,nintp
         if(grcode(i).eq.1 .and. sqrt(gdist(i)).gt.toldist)
     &     nfail = nfail + 1
         if(grcode(i).eq.2) nfail = nfail + 1
      enddo

      nfail_sum = i8glsum(nfail,1)
      if(nfail_sum.gt.0) then
        if(nio.eq.0) write(6,*)
     &    ' WARNING: Unable to find all mesh points in source fld ',
     &    nfail_sum
      endif

      call rzero(fs_gfldout,lx1*ly1*lelv)
!     evaluate inut field at given points
      call fgslib_findpts_eval(intgh_fs,
     &                         fs_gfldout,1,
     &                         grcode,1,
     &                         gproc,1,
     &                         gelid,1,
     &                         grst,ldim2,nintp,
     &                         fld_fs)

      if (nid.eq.1) then
        write(6,*) 'xm1_fs', nintp
        write(6,*) (xm1_fs(i,1,1), i= 10,13)
        write(6,*) 'ym1_fs'
        write(6,*) (ym1_fs(i,1,1), i= 10,13)
        write(6,*) 'fld_fs'
        write(6,*) (fs_gfldout(i), i= 10,13)
      endif  


!      if (nio.eq.0) then
!        write(6,*) 'xg_fs'
!        do iy=1,lyfs
!          write(6,*) (xg_fs(ix,iy,1), ix=1,lxfs)
!        enddo  
!        write(6,*) 'yg_fs'
!        do iy=1,lyfs
!          write(6,*) (yg_fs(ix,iy,1), ix=1,lxfs)
!        enddo
!        write(6,*) 'fld_fs'
!        do iy=1,lyfs
!          write(6,*) (fld_fs(ix,iy,1), ix=1,lxfs)
!        enddo
!      endif  

      if (nio.eq.0) write(6,*) 'FS: Getpts: Done'

      return
      end subroutine fs_getpts
!----------------------------------------------------------------------

      subroutine fs_restore_int

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'MASS'          ! bm1, temporary
      include 'PARALLEL'
      include 'WZ'

      include 'FS_ALE'

      integer e,ifc,n,ne,nfaces

      integer kx1,kx2,ky1,ky2,kz1,kz2
      integer ix,iy,iz,ia,ii

      character cb*3
      real xm1_min,xm1_max,ym1_min,ym1_max
      real glmin,glmax

      call opzero(vx,vy,vz)

!     Get the surface x,y,z
      nfaces = 2*ndim

      ii  = 0
      do e=1,nelv
      do ifc=1,nfaces
        cb  = fs_cbc(ifc,e)
        if (cb.eq.'INT') then
          call facind (kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,ifc)
          do iz=kz1,kz2
          do iy=ky1,ky2
          do ix=kx1,kx2
            ii = ii + 1
            vx(ix,iy,iz,e) = fs_gfldout(ii)
          enddo
          enddo
          enddo
        endif
      enddo
      enddo  

      call outpost(vx,vy,vz,pr,t,'int')

      return
      end subroutine fs_restore_int
!----------------------------------------------------------------------



