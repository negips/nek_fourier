!======================================================================
!     Routines for Weak Laplacian/Weak Stress divergence 
!     in Cylindrical Coordinates.
!     Assuming Fourier in the 3rd (\theta) direction      
!     Author: Prabal S. Negi
!
!====================================================================== 
!-----------------------------------------------------------------------

      subroutine stnrate_cyl (u1r,u2r,u3r,u1i,u2i,u3i,nel,matmod)

      implicit none

!     Compute strainrates

!     CAUTION : Stresses and strainrates share the same scratch commons

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'          ! ym1 == R
      include 'TSTEP'

      include '3DS'

!     Real Variables      
      real erxt,errt,erxx,erxr,errr,ertt
      common /ctmp0/ erxt(lx1*ly1*lz1*lelt)      ! Er_x\theta
     $             , errt(lx1*ly1*lz1*lelt)      ! Er_Rt
      common /ctmp1/ erxx(lx1*ly1*lz1*lelt)      ! Er_xx
     $             , erxr(lx1*ly1*lz1*lelt)      ! Er_xR
     $             , errr(lx1*ly1*lz1*lelt)      ! Er_RR
     $             , ertt(lx1*ly1*lz1*lelt)      ! Er_\theta\theta

!     Made a new scratch array here.      
!     Imaginary Variables
      real eixt,eirt,eixx,eixr,eirr,eitt
      common /scrns/   eixt(lx1*ly1*lz1*lelt)      ! Ei_x\theta
     $               , eirt(lx1*ly1*lz1*lelt)      ! Ei_Rt
     $               , eixx(lx1*ly1*lz1*lelt)      ! Ei_xx
     $               , eixr(lx1*ly1*lz1*lelt)      ! Ei_xR
     $               , eirr(lx1*ly1*lz1*lelt)      ! Ei_RR
     $               , eitt(lx1*ly1*lz1*lelt)      ! Ei_\theta\theta


!     Gradients of real part of velocities      
      real ur1x,ur1r,ur1t,ur2x,ur2r,ur2t,ur3x,ur3r,ur3t
      common /scrpsn1/ ur1x(lx1*ly1*lz1*lelt)      ! du1dx_r
     $               , ur1r(lx1*ly1*lz1*lelt)      ! du1dr_r
     $               , ur2x(lx1*ly1*lz1*lelt)      ! du2dx_r
     $               , ur2r(lx1*ly1*lz1*lelt)      ! du2dr_r
     $               , ur3x(lx1*ly1*lz1*lelt)      ! du3dx_r
     $               , ur3r(lx1*ly1*lz1*lelt)      ! du3dr_r


!     Gradient of imaginary part of velocities      
      real ui1x,ui1r,ui1t,ui2x,ui2r,ui2t,ui3x,ui3r,ui3t
      common /scrpsn2/ ui1x(lx1*ly1*lz1*lelt)      ! du1dx_i
     $               , ui1r(lx1*ly1*lz1*lelt)      ! du1dr_i
     $               , ui2x(lx1*ly1*lz1*lelt)      ! du2dx_i
     $               , ui2r(lx1*ly1*lz1*lelt)      ! du2dr_i
     $               , ui3x(lx1*ly1*lz1*lelt)      ! du3dx_i
     $               , ui3r(lx1*ly1*lz1*lelt)      ! du3dr_i


      real u1r,u2r,u3r
      dimension u1r(lx1,ly1,lz1,lelv)
     $        , u2r(lx1,ly1,lz1,lelv)
     $        , u3r(lx1,ly1,lz1,lelv)

      real u1i,u2i,u3i
      dimension u1i(lx1,ly1,lz1,lelv)
     $        , u2i(lx1,ly1,lz1,lelv)
     $        , u3i(lx1,ly1,lz1,lelv)

      real rinv(lx1*ly1*lz1*lelv),wk1(lx1*ly1*lz1*lelv)
      real wk2(lx1*ly1*lz1*lelv)
      common /scrsf/ rinv,wk1,wk2


      ntot1 = lx1*ly1*lz1*nel

!     rzero3 zeros all 3 components
!     without checking dimensionality

!     Zero real parts      
      call rzero (ur1x,ntot1)
      call rzero (ur1r,ntot1)

      call rzero (ur2x,ntot1)
      call rzero (ur2r,ntot1)

      call rzero (ur3x,ntot1)
      call rzero (ur3r,ntot1)


!     Zero Imaginary parts      
      call rzero (ui1x,ntot1)
      call rzero (ui1r,ntot1)

      call rzero (ui2x,ntot1)
      call rzero (ui2r,ntot1)

      call rzero (ui3x,ntot1)
      call rzero (ui3r,ntot1)

!     Zero Real parts  
      call rzero3 (erxx,errr,ertt,ntot1)
      call rzero3 (erxr,erxt,errt,ntot1)

!     Zero imaginary parts  
      call rzero3 (eixx,eirr,eitt,ntot1)
      call rzero3 (eixr,eixt,eirt,ntot1)

!     uxyz does not zero out variables.
!     values are just added on
!     Derivatives of Real parts      
      call uxyz  (u1r,ur1x,ur1r,wk1,nel)
      call uxyz  (u2r,ur2x,ur2r,wk1,nel)
      call uxyz  (u3r,ur3x,ur3r,wk1,nel)

!     Since the division by the Jacobian is missing
      call invcol2 (ur1x,jacm1,ntot1)
      call invcol2 (ur1r,jacm1,ntot1)
      call invcol2 (ur2x,jacm1,ntot1)
      call invcol2 (ur2r,jacm1,ntot1)
      call invcol2 (ur3x,jacm1,ntot1)
      call invcol2 (ur3r,jacm1,ntot1)

!     Derivatives of Imaginary parts
      call uxyz  (u1i,ui1x,ui1r,wk1,nel)
      call uxyz  (u2i,ui2x,ui2r,wk1,nel)
      call uxyz  (u3i,ui3x,ui3r,wk1,nel)

!     Since the division by the Jacobian is missing
      call invcol2 (ui1x,jacm1,ntot1)
      call invcol2 (ui1r,jacm1,ntot1)
      call invcol2 (ui2x,jacm1,ntot1)
      call invcol2 (ui2r,jacm1,ntot1)
      call invcol2 (ui3x,jacm1,ntot1)
      call invcol2 (ui3r,jacm1,ntot1)

!     rinv = 1/R
      call invers2(rinv,ym1,ntot1)

!!    Er_xx 
!     (Real)
      call copy(erxx,ur1x,ntot1)    ! du1/dx

!!    Ei_xx 
!     (Imaginary)
      call copy(eixx,ui1x,ntot1)    ! du1(im)/dx

!!    Er_Rx
!     (Real)      
      call copy(erxr,ur1r,ntot1)    ! du1/dr
      call add2(erxr,ur2x,ntot1)    ! + du2/dx
      call cmult(erxr,0.5,ntot1)    ! 1/2*[du1/dr + du2/dx]

!!    Ei_Rx
!     (Imaginary)      
      call copy(eixr,ui1r,ntot1)    ! du1(im)/dr
      call add2(eixr,ui2x,ntot1)    ! + du2(im)/dx
      call cmult(eixr,0.5,ntot1)    ! 1/2*[du1(im)/dr + du2(im)/dx]


!!    Er_\thetax
!     (Real)      
      call col3(erxt,u1i,rinv,ntot1)       ! u1(im)/R
      call cmult(erxt,-k_3dsp,ntot1)       ! -k/R*u1(im)
      call add2(erxt,ur3x,ntot1)           ! + du3/dx
      call cmult(erxt,0.5,ntot1)           ! 1/2*[du3/dx - k/R*u1(im)]

!!    Ei_\thetax
!     (Imaginary) 
      call col3(eixt,u1r,rinv,ntot1)       ! u1/R
      call cmult(eixt,k_3dsp,ntot1)        ! k/R*u1
      call add2(eixt,ui3x,ntot1)           ! + du3(im)/dx
      call cmult(eixt,0.5,ntot1)           ! 1/2*[du3(im)/dx + k/R*u1]


!!    Er_RR
!     (Real)      
      call copy(errr,ur2r,ntot1)           ! du2/dr

!!    Ei_RR
!     (Imaginary)      
      call copy(eirr,ui2r,ntot1)           ! du2(im)/dr

      
!!    Er_\thetaR
!     (Real)      
      call copy(errt,ur3r,ntot1)           ! du3/dr
      call col3(wk1,u3r,rinv,ntot1)        ! u3/R
      call sub2(errt,wk1,ntot1)            ! du3/dr - u3/R
      call col3(wk2,u2i,rinv,ntot1)        ! u2(im)/R
      call add2s2(errt,wk2,-k_3dsp,ntot1)  ! du3/dr - u3/R - k/R*u2(im)
      call cmult(errt,0.5,ntot1)           ! 0.5*[du3/dr - u3/R - k/R*u2(im)]

!!    Ei_\thetaR
!     (Imaginary)      
      call copy(eirt,ui3r,ntot1)           ! du3(im)/dr
      call col3(wk1,u3i,rinv,ntot1)        ! u3(im)/R
      call sub2(eirt,wk1,ntot1)            ! du3(im)/dr - u3(im)/R
      call col3(wk2,u2r,rinv,ntot1)        ! u2/R
      call add2s2(eirt,wk2,k_3dsp,ntot1)   ! du3(im)/dr - u3(im)/R + k/R*u2
      call cmult(eirt,0.5,ntot1)           ! 0.5*[du3(im)/dr - u3(im)/R + k/R*u2]


!!    Er_\theta\theta
!     (Real)      
      call col3(ertt,u3i,rinv,ntot1)       ! u3(im)/R      
      call cmult(ertt,-k_3dsp,ntot1)       ! -k/R*u3(im)
      call xaddcol3(ertt,u2r,rinv,ntot1)   ! -k/R*u3(im) + u2/R

!!    Ei_\theta\theta
!     (Imaginary)      
      call col3(eitt,u3r,rinv,ntot1)       ! u3/R      
      call cmult(eitt,k_3dsp,ntot1)        ! k/R*u3
      call xaddcol3(eitt,u2i,rinv,ntot1)   ! k/R*u3 + u2(im)/R


      return
      end
c-----------------------------------------------------------------------
      subroutine stress (h1,h2,nel,matmod,ifaxis)
C
C     MATMOD.GE.0        Fluid material models
C     MATMOD.LT.0        Solid material models
C
C     CAUTION : Stresses and strainrates share the same scratch commons
C
      include 'SIZE'
      common /ctmp1/ txx(lx1,ly1,lz1,lelt)
     $             , txy(lx1,ly1,lz1,lelt)
     $             , tyy(lx1,ly1,lz1,lelt)
     $             , tzz(lx1,ly1,lz1,lelt)
      common /ctmp0/ txz(lx1,ly1,lz1,lelt)
     $             , tyz(lx1,ly1,lz1,lelt)
      common /scrsf/ t11(lx1,ly1,lz1,lelt)
     $             , t22(lx1,ly1,lz1,lelt)
     $             , t33(lx1,ly1,lz1,lelt)
     $             , hii(lx1,ly1,lz1,lelt)
C
      DIMENSION H1(LX1,LY1,LZ1,1),H2(LX1,LY1,LZ1,1)
      LOGICAL IFAXIS

      NTOT1 = lx1*ly1*lz1*NEL

      IF (MATMOD.EQ.0) THEN

C        Newtonian fluids

         CONST = 2.0
         CALL CMULT2 (HII,H1,CONST,NTOT1)
         CALL COL2   (TXX,HII,NTOT1)
         CALL COL2   (TXY,H1 ,NTOT1)
         CALL COL2   (TYY,HII,NTOT1)
         IF (IFAXIS .OR. ldim.EQ.3) CALL COL2 (TZZ,HII,NTOT1)
         IF (ldim.EQ.3) THEN
            CALL COL2 (TXZ,H1 ,NTOT1)
            CALL COL2 (TYZ,H1 ,NTOT1)
         endif
C
      ELSEIF (MATMOD.EQ.-1) THEN
C
C        Elastic solids
C
         CONST = 2.0
         CALL ADD3S   (HII,H1,H2,CONST,NTOT1)
         CALL COPY    (T11,TXX,NTOT1)
         CALL COPY    (T22,TYY,NTOT1)
         CALL COL3    (TXX,HII,T11,NTOT1)
         CALL ADDCOL3 (TXX,H1 ,T22,NTOT1)
         CALL COL3    (TYY,H1 ,T11,NTOT1)
         CALL ADDCOL3 (TYY,HII,T22,NTOT1)
         CALL COL2    (TXY,H2     ,NTOT1)
         IF (IFAXIS .OR. ldim.EQ.3) THEN
            CALL COPY (T33,TZZ,NTOT1)
            CALL COL3    (TZZ,H1 ,T11,NTOT1)
            CALL ADDCOL3 (TZZ,H1 ,T22,NTOT1)
            CALL ADDCOL3 (TZZ,HII,T33,NTOT1)
            CALL ADDCOL3 (TXX,H1 ,T33,NTOT1)
            CALL ADDCOL3 (TYY,H1 ,T33,NTOT1)
         endif
         IF (ldim.EQ.3) THEN
            CALL COL2 (TXZ,H2     ,NTOT1)
            CALL COL2 (TYZ,H2     ,NTOT1)
         endif
C
      endif
C
      return
      end
c-----------------------------------------------------------------------
      subroutine aijuj (au1,au2,au3,nel,ifaxis)
C
      include 'SIZE'
      common /ctmp1/ txx(lx1,ly1,lz1,lelt)
     $             , txy(lx1,ly1,lz1,lelt)
     $             , tyy(lx1,ly1,lz1,lelt)
     $             , tzz(lx1,ly1,lz1,lelt)
      common /ctmp0/ txz(lx1,ly1,lz1,lelt)
     $             , tyz(lx1,ly1,lz1,lelt)
C
      DIMENSION AU1(LX1,LY1,LZ1,1)
     $        , AU2(LX1,LY1,LZ1,1)
     $        , AU3(LX1,LY1,LZ1,1)
      LOGICAL IFAXIS
C
      CALL TTXYZ (AU1,TXX,TXY,TXZ,NEL)
      CALL TTXYZ (AU2,TXY,TYY,TYZ,NEL)
      IF (IFAXIS)    CALL AXITZZ (AU2,TZZ,NEL)
      IF (ldim.EQ.3) CALL TTXYZ  (AU3,TXZ,TYZ,TZZ,NEL)
C
      return
      end
c-----------------------------------------------------------------------

