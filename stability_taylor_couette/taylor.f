c- constants -----------------------------------------------------------

! #define tSTATSTART uparam(1) /* start time for averaging */
! #define tSTATFREQ  uparam(2) /* output frequency for statistics */

c data extraction along wall normal direction
! #define INTP_NMAX 200 /* number of sample points */
! #define XCINT 1.0     /* x coordinate of 1D line*/
! #define ZCINT 1.0     /* z coordinate of 1D line */

c mesh dimensions
! #define BETAM 2.4     /* wall normal stretching parameter */
! #define PI (4.*atan(1.))
! #define XLEN (2.*PI)
! #define ZLEN PI
! #define NUMBER_ELEMENTS_X 16
! #define NUMBER_ELEMENTS_Y 12
! #define NUMBER_ELEMENTS_Z 8

c-----------------------------------------------------------------------
      subroutine uservp (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      integer e

      utrans = 1.
      udiff  = param(2)

      if (ifield .eq. 2) then
         e = gllel(ieg)
         udiff = param(8)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine userf  (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      ffx = 0.0 
      ffy = 0.0
      ffz = 0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine userq  (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      qvol =  0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine userchk

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'SOLN'
      include 'INPUT'
      include 'TSTEP'
!      include 'TOTAL'
      include 'MASS'

      include 'F3D'

      include 'TEST'

      integer ntot1,ntot2
      integer i,j

      integer igeom


      if (istep.eq.0) then
        call frame_start
      endif  

      call frame_monitor

!      call slp_mark_faces()

!      ifaxis = .false.
      ifto = .true.

      ntot1 = lx1*ly1*lz1*nelv
      ntot2 = lx2*ly2*lz2*lelv
      call copy(t,vz,ntot1)

!      call outpost(tmp1,tmp2,tmp3,pr,tmp3,'tmp')

      if (istep.eq.0) then
        call outpost(vx,vy,vz,pr,t,'   ')
        call initp_f3d

        if (.not.iff3d) then
          call init_pertfld_f3d
        endif  

        call rzero(zm1,ntot1)
        call outpost(vx,vy,vz,pr,vz,'ini')
      endif

!!     for k_f3d=0, we seem to get squire modes
!      if (iff3d.and.(abs(k_f3d).lt.1.0e-12)) then
!        call rzero(vzp(1,1),ntot1)
!        call rzero(vzp(1,2),ntot1)
!      endif 

!      if (mod(istep,iostep).eq.0) then
!!      if (istep.le.10) then
!        i = 1
!        call outpost(vxp(1,i),vyp(1,i),vzp(1,i),
!     $               prp(1,i),vzp(1,i),'ptr')
!        i = 2
!        call outpost(vxp(1,i),vyp(1,i),vzp(1,i),
!     $               prp(1,i),vzp(1,i),'pti')
!      endif

!     Call time stepper      
      call tst_solve()

      if (istep.eq.nsteps.or.lastep.eq.1) then
        call frame_end
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine userbc (ix,iy,iz,iside,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      ux = 0.0
      uy = 0.0
      uz = 0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine useric (ix,iy,iz,ieg)

      implicit none

      include 'SIZE'
      include 'INPUT'         ! if3d
      include 'PARALLEL'
      include 'NEKUSE'
      include 'GEOM'

      include 'F3D'

      integer ix,iy,iz,ieg
      real pi

      integer jp
      common /ppointr/ jp

      real fcoeff(3)
      real xl(3)
      real mth_ran_dst

      logical ifcouette
      logical ifpoiseuille
      logical iftaylor

      real r1,r2,omega1,omega2
      real a1,a2
      save r1,r2,a1,a2,omega1,omega2

      integer isave
      save isave
      data isave /0/

      real glmin,glmax
      integer n


      ifcouette         = .false.
      ifpoiseuille      = .false.
      iftaylor          = .true.

      pi = 4.0*atan(1.0)

      if (isave.eq.0) then
        n = lx1*ly1*lz1*nelv
        r1 = glmin(ym1,n)
        r2 = glmax(ym1,n)
        omega1 = 1.0/r1
        omega2 = 0.0
        a1 = (omega2*(r2**2) - omega1*r1*r1)/(r2**2 - r1**2)
        a2 = (omega1 - omega2)*(r1**2)*(r2**2)/(r2**2 - r1**2)
        isave = 1
      endif


      if (jp.eq.0) then

        if (ifpoiseuille) then
          ux = 1.0 - y**2
          uy = 0.
          uz = 0.0 + 0.0
        elseif (ifcouette) then
          ux = 0.0 + 1.0*y
          uy = 0.
          uz = 0.0 + 0.0
        elseif (iftaylor) then
          ux = 0.0 + 0.0
          uy = 0.
          uz = a1*y + a2/y
        endif  
      else

!       perturbation; white noise
        xl(1) = X
        xl(2) = Y
        if (IF3D.or.iff3d) xl(3) = Y+X
        
        if (jp.eq.1) then
          fcoeff(1)=  3.0e4
          fcoeff(2)= -1.5e3
          fcoeff(3)=  0.5e5
        else
          fcoeff(1)=  9.0e4
          fcoeff(2)=  1.5e3
          fcoeff(3)= -2.5e5
        endif          
        ux=UPARAM(1)*mth_ran_dst(ix,iy,iz,ieg,xl,fcoeff)
        if (jp.eq.1) then
          fcoeff(1)=  2.3e4
          fcoeff(2)=  2.3e3
          fcoeff(3)= -2.0e5
        else
          fcoeff(1)=  1.3e4
          fcoeff(2)= -5.8e3
          fcoeff(3)= -1.9e5
        endif
        uy=UPARAM(1)*mth_ran_dst(ix,iy,iz,ieg,xl,fcoeff)
        if (IF3D.or.iff3d) then
           if (jp.eq.1) then           
             fcoeff(1)= 2.e4
             fcoeff(2)= 1.e3
             fcoeff(3)= 1.e5
           else
             fcoeff(1)= -1.9e4
             fcoeff(2)= -8.0e3
             fcoeff(3)=  3.2e5
           endif 
           uz=UPARAM(1)*mth_ran_dst(ix,iy,iz,ieg,xl,fcoeff)
        else
           uz = 0.0
        endif

      endif


      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat   ! This routine to modify element vertices

      implicit none  
  
      include 'SIZE'      ! _before_ mesh is generated, which 
      include 'INPUT'
      include 'GEOM'
      include 'TSTEP'
!      include 'TOTAL'     ! guarantees GLL mapping of mesh.

      integer n,i,j
      real r0

!      ifaxis = .true.   ! just for initialization
      param(42)=1       ! 0: GMRES (nonsymmetric), 1: PCG w/o weights
      param(43)=1       ! 0: Additive multilevel (param 42=0), 1: Original 2 level
      param(44)=0       ! 0: E based Schwartz, 1: A based Schwartz

      n = nelv * 2**ldim
!      xmin = glmin(xc,n)
!      xmax = glmax(xc,n)
!      ymin = glmin(yc,n)
!      ymax = glmax(yc,n)
!      zmin = glmin(zc,n)
!      zmax = glmax(zc,n)
!
!      xscale = XLEN/(xmax-xmin)
!      yscale = 1./(ymax-ymin)
!      zscale = ZLEN/(zmax-zmin)

      pi = 4.*atan(1.0)

      if (abs(uparam(3)).gt.1.0e-6) then
        r0   = abs(uparam(3))+1
      endif

      if (nio.eq.0) write(6,*) 'R0:', r0

      do j=1,nelv
      do i=1,2**ldim
!         xc(i,j) = 2.0*pi*(xc(i,j))*alphai
         yc(i,j) = yc(i,j) + r0
!         yc(i,1) = tanh(BETAM*(2*yc(i,1)-1))/tanh(BETAM)
!         zc(i,1) = zscale*zc(i,1)
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat2   ! This routine to modify mesh coordinates


      include 'SIZE'
      include 'TOTAL'


!      call outpost(vx,vy,vz,pr,t,'   ')

!      do iel=1,nelt
!      do ifc=1,2*ndim
!         if (cbc(ifc,iel,1) .eq. 'W  ') boundaryID(ifc,iel) = 1 
!         cbc(ifc,iel,2) = cbc(ifc,iel,1) 
!         if (cbc(ifc,iel,1) .eq. 'W  ') cbc(ifc,iel,2) = 't  '
!      enddo
!      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat3
      include 'SIZE'
      include 'TOTAL'

!      param(54) = -1  ! use >0 for const flowrate or <0 bulk vel
                      ! flow direction is given by (1=x, 2=y, 3=z) 
!      param(55) = 1.0 ! flowrate/bulk-velocity 

      return
      end
c-----------------------------------------------------------------------
!=======================================================================
!> @brief Register user specified modules
      subroutine frame_usr_register
      implicit none

      include 'SIZE'
      include 'FRAMELP'

!     register modules
      call io_register
      call chkpt_register
      call frame_register_f3d
      call tst_register

      return
      end subroutine
!======================================================================
!> @brief Initialise user specified modules
      subroutine frame_usr_init
      implicit none

      include 'SIZE'
      include 'FRAMELP'

!     initialise modules
      call chkpt_init
      call frame_get_param_f3d
      call tst_init

      return
      end subroutine
!======================================================================
!> @brief Finalise user specified modules
      subroutine frame_usr_end
      implicit none

      include 'SIZE'
      include 'FRAMELP'

      
      return
      end subroutine

!-----------------------------------------------------------------------

      subroutine check_vbasea

      implicit none

      include 'SIZE'
      include 'ARN_ARPD'
      include 'TSTEPPERD'
      include 'F3D'
      include 'TSTEP'
      include 'SOLN'
      include 'INPUT'

      integer i,j

      do j=1,2
        i = 1
        call copy(vxp(1,2),vbasea(i,j),tst_nv)
        i = i + tst_nv
        call copy(vyp(1,2),vbasea(i,j),tst_nv)
        i = i + tst_nv
        if (if3d.or.iff3d) then
          call copy(vzp(1,2),vbasea(i,j),tst_nv)
          i = i + tst_nv
        endif
        if (arna_ifpr) then
          call copy(prp(1,2),vbasea(i,j),tst_np)
          i = i + tst_np
        endif  

        ifto = .true.
        call outpost(vxp(1,2),vyp(1,2),vzp(1,2),prp(1,2),
     $               vzp(1,2),'vba')
      enddo 

      call exitt

      return
      endsubroutine check_vbasea
!---------------------------------------------------------------------- 



c automatically added by makenek
      subroutine usrdat0() 

      return
      end

c automatically added by makenek
      subroutine usrsetvert(glo_num,nel,nx,ny,nz) ! to modify glo_num
      integer*8 glo_num(1)

      return
      end

c automatically added by makenek
      subroutine userqtl

      call userqtl_scig

      return
      end
