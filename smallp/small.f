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

      include '3DS'

      include 'TEST'

      integer ntot1,ntot2
      integer i,j

      integer igeom

!      call slp_mark_faces()

!      ifaxis = .false.
      ifto = .true.

      ntot1 = lx1*ly1*lz1*nelv
      ntot2 = lx2*ly2*lz2*lelv
      call copy(t,vz,ntot1)

!      call outpost(tmp1,tmp2,tmp3,pr,tmp3,'tmp')

      if (istep.eq.0) then
        call outpost(vx,vy,vz,pr,t,'   ')
        call initp_3ds

!       prabal. Temporarily initializing pressure
        do i=1,ntot1
          tmp1(i,1,1,1) = sin(2.0*pi*ym1(i,1,1,1))
          tmp2(i,1,1,1) = sin(2.0*pi*xm1(i,1,1,1))
        enddo
!        call map_pm1_to_pr(tmp1,1)
!        call copy(prp(1,1),pr,ntot2)
!        call map_pm1_to_pr(tmp2,1)
!        call copy(prp(1,2),pr,ntot2)
!        call rzero(pr,ntot2)

!        param(44) = 0

        write(6,*) 'param40', param(40)
        write(6,*) 'param41', param(41)
        write(6,*) 'param42', param(42)
        write(6,*) 'param43', param(43)
        write(6,*) 'param44', param(44)

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

!      call chkdiv_3ds

!      call rzero3(tmp1,tmp2,tmp3,ntot1)
!      call rzero3(tmp5,tmp6,tmp7,ntot1)
!      call rzero3(tmp9,tmp10,tmp11,ntot1)
!
!      call copy(tmp1,vxp(1,1),ntot1)
!      call opdiv_3ds(tmp4,tmp1,tmp2,tmp3) ! du/dx
!      call col2(tmp4,bm2inv,ntot2)
!
!      call copy(tmp5,vyp(1,1),ntot1)
!      call opdiv_3ds(tmp8,tmp6,tmp5,tmp7) ! dv/dy
!      call col2(tmp8,bm2inv,ntot2)
!
!      call copy(tmp9,vzp(1,1),ntot1)
!      call opdiv_3ds(tmp8,tmp11,tmp10,tmp9) ! dw/dz
!!      call col2(tmp12,bm2inv,ntot2)
!
!      if3d_3ds = .false.
!      call opdiv_3ds(tmp4,vxp(1,1),vyp(1,1),tmp3)
!
!      call add3(tmp12,tmp8,tmp4,ntot2)

!      call col2(tmp9,bm1,ntot1)

!      if3d_3ds = .true.
!      call opdiv_3ds(tmp12,vxp(1,1),vyp(1,1),vzp(1,2))

!      call map12_all_3ds(prp,vxp(1,1))
!      call map21_all_3ds(vz,prp)
!      call sub3(tmp2,vxp,vz,ntot1)
!      call map21_weak(vz,prp)
!      call col2(vz,binvm1,ntot1)
!      call sub3(tmp3,vxp,vz,ntot1)


      if (istep.eq.0) then
        call outpost(tmp1,tmp2,tmp3,
     $               tmp4,tmp3,'tmp')
        call outpost(tmp5,tmp6,tmp7,
     $               tmp8,tmp7,'tmp')
        call outpost(tmp9,tmp10,tmp11,
     $               tmp12,tmp11,'tmp')
      endif  

!      call exitt

!      if (istep.eq.1) then
!        
!        call init_pertfld_3ds()
!
!        istep = 1
!        igeom = 2
!        ifield = 1
!        jp = 1
!
!        call setsolv
!        call comment
!        
!        call settime
!        call incomprp_cyl(igeom)
!
!        jp = 1
!        call velpr_update_3ds(igeom)
!
!        call exitt
!      endif  

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

      include '3DS'

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
