!======================================================================
!     Routines for introducing third component in a 2d simulation
!     Author: Prabal S. Negi
!     Right now we assume the third component is homogeneous.
!     Later, a fourier dependence can be added for the linearized solve.
!
!====================================================================== 
!-----------------------------------------------------------------------

      subroutine init_f3d()

      implicit none

      include 'SIZE'
      include 'SOLN'

      include 'F3D'

      integer icalld
      save icalld
      data icalld /0/

      integer nxyz,ntot1
      

      nxyz  = lx1*ly1*lz1
      ntot1 = nxyz*nelv

      iff3d = .true.

!     Need to initialize some variables
!     V3MASK

      call copy(v3mask,v1mask,ntot1)

!     Velocities can be initialized from useric. 


      return
      end subroutine init_f3d
!----------------------------------------------------------------------

      subroutine makef_f3d

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MASS'
      include 'TSTEP'   ! ifield

      include 'F3D'


      integer nxyz,ntot1

      real ta1,ta2,ta3
      common /scruz/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)


      nxyz = lx1*ly1*lz1
      ntot1 = nxyz*nelv

      if (.not.iff3d) then
        call rzero(bfz_f3d,ntot1)    
        return
      endif


!     Build user defined forcing for uz
      call makeuf_f3d

      if3d = .true.
!      if (filterType.eq.2) call make_hpf
!     hpf field stored in ta3
!      call xaddcol3(bfz_f3d,ta3,bm1,ntot1)
      if3d = .false. 

!     Put vx,vy,vz on Dealiased grid (rst form)
!     Standard nek routine. Don't need to change anything (yet) 
!      call set_convect_new(vxd,vyd,vzd,vx,vy,vz)
      call advab_f3d

!     Just leaving it here but we don't need this.
!      call admeshv      ! subroutine not defined yet

      if (iftran) call makeabf_f3d

      if ((iftran.and..not.ifchar).or.
     $    (iftran.and..not.ifnav.and.ifchar)) call makebdf_f3d


      return
      end subroutine makef_f3d

!----------------------------------------------------------------------

      subroutine makeuf_f3d

!     Compute and add: (1) user specified forcing function (FX,FY,FZ)

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
      include 'PARALLEL'
      include 'NEKUSE'

      include 'F3D'

      real ta1,ta2,ta3
      common /scruz/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)

      integer ntot1,iel,i,j,k,ielg


      ntot1 = lx1*ly1*lz1*nelv

      time = time-dt
      call rzero(bfz_f3d,ntot1)

      do 100 iel=1,nelv
         ielg = lglel(iel)
         do 100 k=1,lz1
         do 100 j=1,ly1
         do 100 i=1,lx1
            call nekasgn (i,j,k,iel)
            call userf   (i,j,k,ielg)
            bfz_f3d(i,j,k,iel) = ffz
 100  continue

      call col2  (bfz_f3d,bm1,ntot1)
      time = time+dt

      return
      end subroutine makeuf_f3d

!----------------------------------------------------------------------
      subroutine advab_f3d

!     Eulerian scheme, add convection term to forcing function 
!     at current time step.

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'

      include 'F3D'

      include 'TEST'

      real ta1,ta2,ta3
      common /scruz/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)

      integer ntot1

      ntot1 = lx1*ly1*lz1*nelv

!      call setup_convect(2)

      call convop  (ta3,ta1)
      call subcol3 (bfz_f3d,ta3,bm1,ntot1)


      return
      end subroutine advab_f3d
c-----------------------------------------------------------------------

      subroutine makeabf_f3d
!
!     Eulerian scheme, add convection term to forcing function 
!     at current time step.

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
      include 'INPUT'

      include 'F3D'


      real ta1,ta2,ta3
      common /scruz/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)

      integer ntot1
      real ab0,ab1,ab2

      ntot1 = lx1*ly1*lz1*nelv


      ab0 = ab(1)
      ab1 = ab(2)
      ab2 = ab(3)
      call add3s2 (ta3,abz1,abz2,ab1,ab2,ntot1)
      call copy   (abz2,abz1,ntot1)
      call copy   (abz1,bfz,ntot1)
      call add2s1 (bfz_f3d,ta3,ab0,ntot1)
      if (.not.iflomach) call col2 (bfz_f3d,vtrans,ntot1)

      return
      end subroutine makeabf_f3d

!-----------------------------------------------------------------------
      subroutine makebdf_f3d

!     Add contributions to F from lagged BD terms.

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'GEOM'
      include 'INPUT'
      include 'TSTEP'

      include 'F3D'

      real ta1,ta2,ta3,tb1,tb2,tb3,h2
      common /scrns/ ta1(lx1,ly1,lz1,lelv)
     $ ,             ta2(lx1,ly1,lz1,lelv)
     $ ,             ta3(lx1,ly1,lz1,lelv)
     $ ,             tb1(lx1,ly1,lz1,lelv)
     $ ,             tb2(lx1,ly1,lz1,lelv)
     $ ,             tb3(lx1,ly1,lz1,lelv)
     $ ,             h2 (lx1,ly1,lz1,lelv)


      integer ilag,ntot1
      real const

      ntot1 = lx1*ly1*lz1*nelv
      const = 1./dt

      if(iflomach) then
        call cfill(h2,const,ntot1)
      else
        call cmult2(h2,vtrans(1,1,1,1,ifield),const,ntot1)
      endif

      call col3(tb3,vz,bm1,ntot1)
      call cmult(tb3,bd(2),ntot1)

      do 100 ilag=2,nbd
         if (ifgeom) then

            call col3(ta3,vzlag(1,1,1,1,ilag-1),bm1lag(1,1,1,1,ilag-1),
     $                                          ntot1)
            call cmult(ta3,bd(ilag+1),ntot1)
         else
            call col3(ta3,vzlag(1,1,1,1,ilag-1),bm1,ntot1)
            call cmult(ta3,bd(ilag+1),ntot1)
         endif
         call add2 (tb3,ta3,ntot1)
 100  continue
      call xaddcol3(bfz_f3d,tb3,h2,ntot1)     

      return
      end subroutine makebdf_f3d
!-----------------------------------------------------------------------

      subroutine lagvel_f3d

!     Keep old velocity field(s)

      implicit none 

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'TSTEP'

      integer ilag,ntot1

      ntot1 = lx1*ly1*lz1*nelv

      do 100 ilag=3-1,2,-1
        call copy (vzlag (1,1,1,1,ilag),vzlag (1,1,1,1,ilag-1),ntot1)
 100  continue

      call copy (vzlag,vz,ntot1)

      return
      end subroutine lagvel_f3d
!----------------------------------------------------------------------
      subroutine ophx_f3d (out1,out2,out3,inp1,inp2,inp3,h1,h2)

!     OUT = (H1*A+H2*B) * INP  

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'

      include 'F3D'


      real out1 (lx1,ly1,lz1,1)
      real out2 (lx1,ly1,lz1,1)
      real out3 (lx1,ly1,lz1,1)
      real inp1 (lx1,ly1,lz1,1)
      real inp2 (lx1,ly1,lz1,1)
      real inp3 (lx1,ly1,lz1,1)
      real h1   (lx1,ly1,lz1,1)
      real h2   (lx1,ly1,lz1,1)

      integer imesh,matmod
      

      imesh = 1

      if (ifstrs) then
         matmod = 0
         call axhmsf (out1,out2,out3,inp1,inp2,inp3,h1,h2,matmod)
      else

!        the numbers are only needed for axis-symmetric formulation
!        need to come back to this later. 
         call axhelm (out1,inp1,h1,h2,imesh,1)
         call axhelm (out2,inp2,h1,h2,imesh,2)
         call axhelm (out3,inp3,h1,h2,imesh,3)
      endif

      return
      end subroutine ophx_f3d
!-----------------------------------------------------------------------
      subroutine cresvif_f3d (resv1,resv2,resv3,h1,h2)

!     Compute startresidual/right-hand-side in the velocity solver

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      

!      include 'TOTAL'

      include 'F3D'

      real           resv1 (lx1,ly1,lz1,1)
      real           resv2 (lx1,ly1,lz1,1)
      real           resv3 (lx1,ly1,lz1,1)
      real           h1    (lx1,ly1,lz1,1)
      real           h2    (lx1,ly1,lz1,1)
      common /scruz/ w1    (lx1,ly1,lz1,lelv)
     $ ,             w2    (lx1,ly1,lz1,lelv)
     $ ,             w3    (lx1,ly1,lz1,lelv)

      integer igeom
      common /cgeom/ igeom

      integer ntot1,ntot2
      

      ntot1 = lx1*ly1*lz1*nelv
      ntot2 = lx2*ly2*lz2*nelv

      if (igeom.eq.2) call lagvel

!     prabal
!     for 3d solve
      if (igeom.eq.2) call lagvel_f3d

      if3d = .true.
      call bcdirvc (vx,vy,vz,v1mask,v2mask,v3mask)
      if3d = .false.

!     prabal. We don't care about traction conditions for now.
!     Maybe need to look at it if added stiffness terms are needed
!     Or if surface tension is needed
      call bcneutr

      call extrapp (pr,prlag)
      call opgradt (resv1,resv2,resv3,pr)
!      call opzero(resv1,resv2,resv3)
      call rzero(resv3,ntot1)             ! homogeneous in z

      call copy(bfz,bfz_f3d,ntot1)
      call opadd2(resv1,resv2,resv3,bfx,bfy,bfz)
      call add2(resv3,bfz,ntot1)

!     prabal
      call ophx_f3d(w1,w2,w3,vx,vy,vz,h1,h2)
      call opsub2(resv1,resv2,resv3,w1,w2,w3)
      call sub2(resv3,w3,ntot1)

      return
      end subroutine cresvif_f3d
!-----------------------------------------------------------------------
      subroutine plan3_f3d (igeom)

!     Compute pressure and velocity using consistent approximation spaces.     
!     Operator splitting technique.

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'EIGEN'
      include 'SOLN'
      include 'TSTEP'

      include 'F3D'

      real resv1,resv2,resv3
      real dv1,dv2,dv3
      real h1,h2
      common /scrns/  resv1 (lx1,ly1,lz1,lelv)
     $ ,              resv2 (lx1,ly1,lz1,lelv)
     $ ,              resv3 (lx1,ly1,lz1,lelv)
     $ ,              dv1   (lx1,ly1,lz1,lelv)
     $ ,              dv2   (lx1,ly1,lz1,lelv)
     $ ,              dv3   (lx1,ly1,lz1,lelv)
      common /scrvh/  h1    (lx1,ly1,lz1,lelv)
     $ ,              h2    (lx1,ly1,lz1,lelv)

      real ut1(lx1,ly1,lz1,lelt)
      real ut2(lx1,ly1,lz1,lelt)
      real ut3(lx1,ly1,lz1,lelt)
      integer flowtype(lelt)
      common /testvel1/ ut1,ut2,ut3,flowtype

      real ut4(lx1,ly1,lz1,lelt)
      real ut5(lx1,ly1,lz1,lelt)
      real ut6(lx1,ly1,lz1,lelt)

      common /testvel2/ ut4,ut5,ut6

      integer intype
      integer igeom
      integer ntot1



      ntot1 = lx1*ly1*lz1*nelv   

      if (igeom.eq.1) then

!        old geometry

         call makef

         call makef_f3d

      else

!        new geometry, new b.c.

         intype = -1
         call sethlm  (h1,h2,intype)

!-------------------------------------------------- 
!!        prabal
!!        for introducing slip across element faces   
!!        ut1, ut2, ut3 ---> contains slip velocity arrays   
!         call opzero  (ut4,ut5,ut6)
!         call ophx    (ut4,ut5,ut6,ut1,ut2,ut3,h1,h2)       ! Ax
!         call opsub2  (bfx,bfy,bfz,ut4,ut5,ut6)             ! bfx - Ax
!-------------------------------------------------- 

         call cresvif_f3d (resv1,resv2,resv3,h1,h2)
     
         if3d = .true. 
         call ophinv(dv1,dv2,dv3,resv1,resv2,resv3,h1,h2,tolhv,nmxv)
         if3d = .false.   


         call opadd2(vx,vy,vz,dv1,dv2,dv3)
         call add2(vz,dv3,ntot1) 

         call incomprn(vx,vy,vz,pr)

      endif

      return
      end subroutine plan3_f3d

!----------------------------------------------------------------------



!----------------------------------------------------------------------





