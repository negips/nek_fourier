!----------------------------------------------------------------------
!     Author: Prabal Negi
!     Description: Routines for Generating global basis function.
!                  for surface representation
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      subroutine fs_smooth_meshmv(wx,wy,wz)

      implicit none

      include 'SIZE'
      include 'FS_ALE'

      real wx(lx1,ly1,lz1,lelv),wy(lx1,ly1,lz1,lelv)
      real wz(lx1,ly1,lz1,lelv)

      integer icalld
      save icalld
      data icalld /0/

      if (.not.fs_ifgsm) return


      if (icalld.eq.0) then
        call fs_global_basis
        icalld = icalld+1
      endif

      call fs_gllo_xyz
      call fs_gllo_flds(wx,wy,wz)
      call fs_intp_setup
      call fs_get_localpts
      call fs_get_globalpts
      call fs_restore_int(wx,wy,wz,'Norm')

!     Correction for tangential movement to minimize
!     mesh deformation
      if (fs_iftc) then
        call fs_tang_corr(wx,wy,wz)
      endif  

!     Free the handles      
      call fgslib_findpts_free(intgh_fs)
      call fgslib_findpts_free(intlh_fs)

      call mntr_log(fs_id,fs_log,'Interface Smoothening Done')

      return
      end subroutine fs_smooth_meshmv
!----------------------------------------------------------------------

      subroutine fs_tang_corr(wx,wy,wz)

      implicit none

      include 'SIZE'
      include 'FS_ALE'
      include 'GEOM'

      include 'TEST'
      include 'SOLN'

      real wx(lx1,ly1,lz1,lelv),wy(lx1,ly1,lz1,lelv)
      real wz(lx1,ly1,lz1,lelv)
     
      integer i,n
      integer e,j1,j2,j3
      real p,s          ! sign
      real mu           ! exponential decay
      real dy, tol

      integer ifc,nface
      integer js1,js2,jf1,jf2,jskip1,jskip2

      character cb*3

!     Temporary arrays
      real wk1,wk2,wk3,wk4,wk5,wk6,wk7
      common /scrns3/ wk1(lx1,ly1,lz1,lelt),    ! Tangengial directions     
     $               wk2(lx1,ly1,lz1,lelt),     ! Tangengial directions
     $               wk3(lx1,ly1,lz1,lelt),     ! Tangengial directions
     $               wk4(lx1,ly1,lz1,lelt),     ! Displacement x 
     $               wk5(lx1,ly1,lz1,lelt),     ! Displacement y
     $               wk6(lx1,ly1,lz1,lelt),     ! Displacement z
     $               wk7(lx1,ly1,lz1,lelt)


!      call fs_gllo_xyz
!      call fs_intp_setup

      n = lx1*ly1*lz1*nelv

      call opzero(wk1,wk2,wk3)
      nface = 2*ndim
      do e=1,nelv
      do ifc=1,nface
        cb  = fs_cbc(ifc,e)
!       Change these conditions for actual case            
        if (cb.eq.'INT') then
          call facind2 (js1,jf1,jskip1,js2,jf2,jskip2,ifc)
          i = 0
          do j2=js2,jf2,jskip2
          do j1=js1,jf1,jskip1
            i = i + 1 
            wk1(j1,j2,1,e) = t1x(i,1,ifc,e)
            wk2(j1,j2,1,e) = t1y(i,1,ifc,e)
            if (ndim.eq.3) wk3(j1,j2,1,e) = t1z(i,1,ifc,e)
          enddo ! j2
          enddo ! j1
        endif                 
      enddo
      enddo
      call dsavg(wk1)
      call dsavg(wk2)
      if (ndim.eq.3) call dsavg(wk3)
      call fs_int_project(wk1,wk2,wk3)

!      call outpost(wk1,wk2,wk3,pr,wk3,'tst')

      call fs_gllo_flds(wk1,wk2,wk3)
      call fs_get_localpts
      call fs_get_globalpts
      call fs_restore_int(wk1,wk2,wk3,'Tang')    ! Tangential velocities

!      call outpost(wk1,wk2,wk3,pr,wk3,'tst')

      call opcopy(wk4,wk5,wk6,xm1,ym1,zm1)
      call opsub2(wk4,wk5,wk6,xm0_fs,ym0_fs,zm0_fs) ! dx,dy,dz

     
!     No correction for SYM/O points
      do i=1,fs_nsymo
        e  = fs_ie(i)
        j1 = fs_ix(i)
        j2 = fs_iy(i)
        j3 = fs_iz(i)
        wk1(j1,j2,j3,e) = 0.
        wk2(j1,j2,j3,e) = 0.
        wk3(j1,j2,j3,e) = 0.
      enddo  

      mu = 0.1
      tol = 1.0e-14
      do i=1,n
        dy          = wk5(i,1,1,1)
        p           = wk2(i,1,1,1)
        if (abs(p).gt.tol) then
          s           = - p/abs(p)
        else
          s           = 0.0
        endif  
        wx(i,1,1,1) = wx(i,1,1,1) + s*mu*dy*wk1(i,1,1,1)
        wy(i,1,1,1) = wy(i,1,1,1) + s*mu*dy*wk2(i,1,1,1)
        if (ndim.eq.3) wz(i,1,1,1) = wz(i,1,1,1) + s*mu*dy*wk3(i,1,1,1)
      enddo  

!      call outpost(wx,wy,wz,pr,wz,'tst')
!     Add v2x,v2y...             

!!     Free the handles      
!      call fgslib_findpts_free(intgh_fs)
!      call fgslib_findpts_free(intlh_fs)


      return
      end subroutine fs_tang_corr
!----------------------------------------------------------------------
      subroutine fs_global_basis

      implicit none

      include 'SIZE'
      include 'FS_ALE'

      integer ix,iy,iz

      call rzero(zx_fs,lxfs)
      call rzero(zy_fs,lxfs)
      call rzero(wx_fs,lxfs)
      call rzero(wy_fs,lxfs)

      call zwgll(zx_fs,wx_fs,lxfs)
      call zwgll(zy_fs,wy_fs,lyfs)
   
      if (ndim.eq.3) then
        do ix=1,lxfs
        do iy=1,lyfs
          w2_fs(ix,iy) = wx_fs(ix)*wy_fs(iy)
        enddo
        enddo
      else
        call copy(w2_fs,wx_fs,lxfs)
      endif

      return
      end subroutine fs_global_basis
!----------------------------------------------------------------------
      subroutine fs_gllo_xyz

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

!      real wx(lx1,ly1,lz1,lelv),wy(lx1,ly1,lz1,lelv)
!      real wz(lx1,ly1,lz1,lelv)

      integer e,ifc,n,ne,nfaces

      integer kx1,kx2,ky1,ky2,kz1,kz2
      integer ix,iy,iz,ia,ii

      character cb*3
      real xm1_min,xm1_max,ym1_min,ym1_max
      real glmin,glmax

      real s(3)   ! surface normals
      real vn     ! normal velocity      

      xm1_min = glmin(xm1,lx1*ly1*lz1*nelv)
      xm1_max = glmax(xm1,lx1*ly1*lz1*nelv)
      ym1_min = glmin(ym1,lx1*ly1*lz1*nelv)
      ym1_max = glmax(ym1,lx1*ly1*lz1*nelv)
 
      call rzero(xg_fs,lxfs*lyfs*2)
      call rzero(yg_fs,lxfs*lyfs*2)
      call rzero(zg_fs,lxfs*lyfs*2)

      do iy=1,lyfs
        do ix=1,lxfs
          xg_fs(ix,iy,1)=(zx_fs(ix)+1.0)*(xm1_max-xm1_min)/2.0 + xm1_min
          yg_fs(ix,iy,1)=(zy_fs(iy)+1.0)*(ym1_max-ym1_min)/2.0 + ym1_min
          if (ndim.eq.2) xg_fs(ix,iy,1) = zx_fs(ix)
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
          ne = ne+1
          call facind (kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,ifc)
          do iz=kz1,kz2
          do iy=ky1,ky2
          do ix=kx1,kx2
!            call getSnormal(s,ix,iy,iz,ifc,e)   ! surface normals
            ii = ii + 1
            xi_fs(ii) = xm1(ix,iy,iz,e)
            yi_fs(ii) = ym1(ix,iy,iz,e)
            if (ndim.eq.3) zi_fs(ii) = zm1(ix,iy,iz,e)
            if (ndim.eq.2) xi_fs(ii) = zgm1(ix,1)

            xm1_fs(ix,iy,ne) = xm1(ix,iy,iz,e)
            ym1_fs(ix,iy,ne) = ym1(ix,iy,iz,e)
            if (ndim.eq.3) zm1_fs(ia,1,ne) = zm1(ix,iy,iz,e)
!           ndim.eq.2 needs very special treatment
            if (ndim.eq.2) then
              if (kx1.eq.kx2) then
                do ia = 1,lx1
                  xm1_fs(ia,iy,ne)  = zgm1(ia,1)
                  ym1_fs(ia,iy,ne)  = ym1(ix,iy,iz,e)
                enddo
              elseif (ky1.eq.ky2) then
                do ia = 1,ly1
                  xm1_fs(ix,ia,ne)  = zgm1(ia,1)
                  ym1_fs(ix,ia,ne)  = ym1(ix,iy,iz,e)
                enddo
              endif     ! kx1.eq.kx2
            endif       ! ndim.eq.2
          enddo         ! ix
          enddo         ! iy
          enddo         ! iz
        endif           ! cb.eq.INT
      enddo             ! ifc
      enddo             ! e

      return
      end subroutine fs_gllo_xyz
!----------------------------------------------------------------------
       subroutine fs_gllo_flds(wx,wy,wz)

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

      real wx(lx1,ly1,lz1,lelv),wy(lx1,ly1,lz1,lelv)
      real wz(lx1,ly1,lz1,lelv)

      integer e,ifc,n,ne,nfaces

      integer kx1,kx2,ky1,ky2,kz1,kz2
      integer ix,iy,iz,ia,ii

      character cb*3

      real s(3)   ! surface normals
      real vn     ! normal velocity      

!     Get the surface x,y,z
      nfaces = 2*ndim
      ne     = 0              ! total number of interface elements
      ii     = 0
      do e=1,nelv
      do ifc=1,nfaces
        cb  = fs_cbc(ifc,e)
        if (cb.eq.'INT') then
          ne = ne+1
          call facind (kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,ifc)
          do iz=kz1,kz2
          do iy=ky1,ky2
          do ix=kx1,kx2
!            call getSnormal(s,ix,iy,iz,ifc,e)   ! surface normals
            lfld_fs(ix,iy,ne,1)= wx(ix,iy,iz,e) ! ym1_fs(ix,iy,ne)**2
            lfld_fs(ix,iy,ne,2)= wy(ix,iy,iz,e) ! ym1_fs(ix,iy,ne)**2
            if (ndim.eq.3) lfld_fs(ix,iy,ne,1)= wz(ix,iy,iz,e)
!           ndim.eq.2 needs very special treatment
            if (ndim.eq.2) then
              if (kx1.eq.kx2) then
                do ia = 1,lx1
                  lfld_fs(ia,iy,ne,1) = wx(ix,iy,iz,e)
                  lfld_fs(ia,iy,ne,2) = wy(ix,iy,iz,e)
                enddo
              elseif (ky1.eq.ky2) then
                do ia = 1,ly1
                  lfld_fs(ix,ia,ne,1) = wx(ix,iy,iz,e)
                  lfld_fs(ix,ia,ne,2) = wy(ix,iy,iz,e)
                enddo
              endif     ! kx1.eq.kx2
            endif       ! ndim.eq.2
          enddo         ! ix
          enddo         ! iy
          enddo         ! iz
        endif           ! cb.eq.INT
      enddo             ! ifc
      enddo             ! e

      return
      end subroutine fs_gllo_flds
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

      call mntr_log(fs_id,fs_log,'Global Interpolation Setup Done')

!     initialize interpolation tool using local sem mesh
      nxf   = 2*lx1
      nyf   = 2*ly1
      nzf   = 2*1
      nhash = fs_nel*lx1*ly1
      nmax  = 128
!     We do all Global calculations on nid 0
      nels  = fs_nel

      
!     Interpolation handle for SEM surface mesh. 
      call fgslib_findpts_setup(intlh_fs,nekcomm,np,ldim2,
     &                          xm1_fs,ym1_fs,zm1_fs,lx1,ly1,1,
     &                          nels,nxf,nyf,nzf,bb_t,
     &                          nhash,nhash,nmax,tol)

      call mntr_log(fs_id,fs_log,'Local Interpolation Setup Done')

      return
      end subroutine fs_intp_setup
!---------------------------------------------------------------------- 
      subroutine fs_get_globalpts

!     Find the local sem surface points within the global mesh.

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
        nintp  = fs_nel*lx1*ly1
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

      call rzero(fs_gfldout,lx1*ly1*lelv*3)
!     evaluate input field at sem points
      do i=1,ndim
        call fgslib_findpts_eval(intgh_fs,
     &                           fs_gfldout(1,1,1,i),1,
     &                           grcode,1,
     &                           gproc,1,
     &                           gelid,1,
     &                           grst,ldim2,nintp,
     &                           fld_fs(1,1,i))
      enddo  

      call mntr_log(fs_id,fs_log,'Get Global pts: Done')

      return
      end subroutine fs_get_globalpts
!----------------------------------------------------------------------
      subroutine fs_get_localpts

!     Find the Global mesh points within the local sem surface mesh.

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'GFLDR'
      include 'PARALLEL'

      include 'FS_ALE'

      integer i,j

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

      nfail = 0
      toldist = 5e-14

      if (nid.eq.0) then
        nintp  = lxfs*lyfs
      else
        nintp  = 0  
      endif
      ldim2  = 2

      call fgslib_findpts(intlh_fs,
     &                    grcode,1,
     &                    gproc,1,
     &                    gelid,1,
     &                    grst,ldim2,
     &                    gdist,1,
     &                    xg_fs,1,
     &                    yg_fs,1,
     &                    zg_fs,1,nintp)


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

!     Evaluate fields at global mesh points      
      do i=1,ndim
!       evaluate inut field at given points
        call fgslib_findpts_eval(intlh_fs,
     &                           fs_lfldout(1,1,i),1,
     &                           grcode,1,
     &                           gproc,1,
     &                           gelid,1,
     &                           grst,ldim2,nintp,
     &                           lfld_fs(1,1,1,i))

!        write(6,*) 'Fldout', (fs_lfldout(j,1,i),j=1,nintp)

!       This is now the globally smooth field from which we interpolate
!       the local sem points 
        if (nintp.gt.0) then
          call copy(fld_fs(1,1,i),fs_lfldout(1,1,i),nintp)
        else
          call rzero(fld_fs(1,1,i),lxfs*lyfs)
        endif  
      enddo  

      call mntr_log(fs_id,fs_log,'Get Local pts: Done')

      return
      end subroutine fs_get_localpts
!----------------------------------------------------------------------

      subroutine fs_restore_int(wx,wy,wz,dirc)

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'MASS'          ! bm1, temporary
      include 'PARALLEL'
      include 'TSTEP'

      include 'FS_ALE'

      real wx(lx1,ly1,lz1,lelv),wy(lx1,ly1,lz1,lelv)
      real wz(lx1,ly1,lz1,lelv)

      integer e,ifc,n,ne,nfaces

      integer kx1,kx2,ky1,ky2,kz1,kz2
      integer ix,iy,iz,ia,ii

      character cb*3
      real xm1_min,xm1_max,ym1_min,ym1_max
      real glmin,glmax,glsum

      real erx,ery,erz,err         ! errors due to smoothening
      real ar,arsum                ! area

      character str*120
      integer loglev

      character dirc*4

!     Get the surface x,y,z
      nfaces = 2*ndim

      ar  = 0.0
      ne  = 0
      erx = 0.0
      ery = 0.0
      erz = 0.0
      do e=1,nelv
      do ifc=1,nfaces
        cb  = fs_cbc(ifc,e)
        if (cb.eq.'INT') then
          ne = ne + 1
          call facind (kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,ifc)
          ii = 0
          do iz=kz1,kz2
          do iy=ky1,ky2
          do ix=kx1,kx2
            ii = ii+1
            erx = erx + 
     $      area(ii,1,ifc,e)*(fs_gfldout(ix,iy,ne,1)-wx(ix,iy,iz,e))**2
            ery = ery + 
     $      area(ii,1,ifc,e)*(fs_gfldout(ix,iy,ne,2)-wy(ix,iy,iz,e))**2
            if (ndim.eq.3)
     $        erz = erz + 
     $      area(ii,1,ifc,e)*(fs_gfldout(ix,iy,ne,3)-wz(ix,iy,iz,e))**2
          
            ar = ar +  area(ii,1,ifc,e)

            wx(ix,iy,iz,e) = fs_gfldout(ix,iy,ne,1)
            wy(ix,iy,iz,e) = fs_gfldout(ix,iy,ne,2)
            if (ndim.eq.3) wz(ix,iy,iz,e) = fs_gfldout(ix,iy,ne,3)
          enddo
          enddo
          enddo
        endif
      enddo
      enddo  

      err = (erx + ery + erz)
      erx = glsum(err,1)      ! now contains total error
      arsum = glsum(ar,1)
      err = sqrt(erx/arsum)

      write(str,'(A4,1x,A24,1x,E11.4E2)') dirc,
     $      'Smooth projection error:',err

!     We (almost) always want to see this error      
      loglev = 7
      call mntr_log(fs_id,loglev,str)

      return
      end subroutine fs_restore_int
!----------------------------------------------------------------------



