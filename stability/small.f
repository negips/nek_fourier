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

!       prabal. Temporarily initializing pressure
        do i=1,ntot1
          tmp1(i,1,1,1) = sin(2.0*pi*ym1(i,1,1,1))
          tmp2(i,1,1,1) = sin(2.0*pi*xm1(i,1,1,1))
        enddo

        call rzero(zm1,ntot1)
        call outpost(vx,vy,vz,pr,vz,'ini')
      endif

      if (mod(istep,iostep).eq.0) then
!      if (istep.le.10) then
        i = 1
        call outpost(vxp(1,i),vyp(1,i),vzp(1,i),
     $               prp(1,i),vzp(1,i),'ptr')
        i = 2
        call outpost(vxp(1,i),vyp(1,i),vzp(1,i),
     $               prp(1,i),vzp(1,i),'pti')
      endif

      if (istep.eq.0) then
        call outpost(tmp1,tmp2,tmp3,
     $               tmp4,tmp3,'tmp')
        call outpost(tmp5,tmp6,tmp7,
     $               tmp8,tmp7,'tmp')
        call outpost(tmp9,tmp10,tmp11,
     $               tmp12,tmp11,'tmp')
      endif  

!      call exitt

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
      include 'PARALLEL'
!      include 'TOTAL'
      include 'NEKUSE'

      include 'F3D'

      integer ix,iy,iz,ieg
      real pi

      integer jp
      common /ppointr/ jp

      pi = 4.0*atan(1.0)

      if (jp.eq.0) then
        ux = 0.0 + 1.0*y
        uy = 0.
        uz = 0.0 + 1.0*y
      else
        ux = 0.0 + (1.0e-0)*sin(x)*sin(2*pi*(y-1.0)/3)
        uy = 0.0 + (2.0e-1)*sin(2*pi*(y-1.0)/3)*sin(jp + x + y)
!        uz = -1.0 + (2.0e-0)*rand()
        uz = ux*sin(jp + x*y+0.)

!        ux = cos(x)
!        uy = sin(y)
!      
!        if (jp.eq.1) then
!          uz = 1.0/k_3dsp*(sin(x) - cos(y))
!        elseif (jp.eq.2) then
!          uz = 1.0/k_3dsp*(-sin(x) + cos(y))
!        endif

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

!      ifaxis = .true.   ! just for initialization
      param(42)=1
      param(43)=1
      param(44)=1

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

      write(6,*) n,nelv,lelv

      do j=1,nelv
      do i=1,2**ldim
!         xc(i,j) = pi*(xc(i,j)+2.0)
         yc(i,j) = yc(i,j) + 1.0
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
      call frame_get_param_f3d

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
