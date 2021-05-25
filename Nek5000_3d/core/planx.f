      SUBROUTINE PLAN3 (IGEOM)
C-----------------------------------------------------------------------
C
C     Compute pressure and velocity using consistent approximation spaces.     
C     Operator splitting technique.
C
C-----------------------------------------------------------------------
      include 'SIZE'
      include 'INPUT'
      include 'EIGEN'
      include 'SOLN'
      include 'TSTEP'

      include 'F3D'
C
      COMMON /SCRNS/  RESV1 (LX1,LY1,LZ1,LELV)
     $ ,              RESV2 (LX1,LY1,LZ1,LELV)
     $ ,              RESV3 (LX1,LY1,LZ1,LELV)
     $ ,              DV1   (LX1,LY1,LZ1,LELV)
     $ ,              DV2   (LX1,LY1,LZ1,LELV)
     $ ,              DV3   (LX1,LY1,LZ1,LELV)
      COMMON /SCRVH/  H1    (LX1,LY1,LZ1,LELV)
     $ ,              H2    (LX1,LY1,LZ1,LELV)

!     prabal
      real ut1(lx1,ly1,lz1,lelt)
      real ut2(lx1,ly1,lz1,lelt)
      real ut3(lx1,ly1,lz1,lelt)
      integer flowtype(lelt)
      common /testvel1/ ut1,ut2,ut3,flowtype

      real ut4(lx1,ly1,lz1,lelt)
      real ut5(lx1,ly1,lz1,lelt)
      real ut6(lx1,ly1,lz1,lelt)

      common /testvel2/ ut4,ut5,ut6


!     prabal
      if (iff3d) then
        call plan3_f3d(igeom)
        return
      endif    


C
      IF (IGEOM.EQ.1) THEN
C
C        Old geometry
C
         CALL MAKEF
C
      ELSE
C
C        New geometry, new b.c.
C
         intype = -1
         call sethlm  (h1,h2,intype)


!-------------------------------------------------- 
!!        prabal
!!        For introducing slip across element faces   
!!        ut1, ut2, ut3 ---> Contains slip velocity arrays   
!         call opzero  (ut4,ut5,ut6)
!         call ophx    (ut4,ut5,ut6,ut1,ut2,ut3,h1,h2)       ! Ax
!         call opsub2  (bfx,bfy,bfz,ut4,ut5,ut6)             ! bfx - Ax
!-------------------------------------------------- 

         call cresvif (resv1,resv2,resv3,h1,h2)
         call ophinv  (dv1,dv2,dv3,resv1,resv2,resv3,h1,h2,tolhv,nmxv)
         call opadd2  (vx,vy,vz,dv1,dv2,dv3)

!!        prabal
!         call opadd2(vx,vy,vz,ut1,ut2,ut3)            ! Add slip velocity back
c
         call incomprn(vx,vy,vz,pr)
C
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE LAGPRES 
C--------------------------------------------------------------------
C
C     Keep old pressure values
C
C--------------------------------------------------------------------
      include 'SIZE'
      include 'SOLN'
      include 'TSTEP'

      common /cgeom/ igeom

      IF (NBDINP.EQ.3.and.igeom.le.2) THEN
         NTOT2 = lx2*ly2*lz2*NELV
         CALL COPY (PRLAG,PR,NTOT2)
      ENDIF
      RETURN
      END
C
      subroutine cresvif (resv1,resv2,resv3,h1,h2)
C---------------------------------------------------------------------
C
C     Compute startresidual/right-hand-side in the velocity solver
C
C---------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      REAL           RESV1 (LX1,LY1,LZ1,1)
      REAL           RESV2 (LX1,LY1,LZ1,1)
      REAL           RESV3 (LX1,LY1,LZ1,1)
      REAL           H1    (LX1,LY1,LZ1,1)
      REAL           H2    (LX1,LY1,LZ1,1)
      COMMON /SCRUZ/ W1    (LX1,LY1,LZ1,LELV)
     $ ,             W2    (LX1,LY1,LZ1,LELV)
     $ ,             W3    (LX1,LY1,LZ1,LELV)

      common /cgeom/ igeom

!     prabal
      real ut1(lx1,ly1,lz1,lelt)
      real ut2(lx1,ly1,lz1,lelt)
      real ut3(lx1,ly1,lz1,lelt)
      integer flowtype(lelt)
      common /testvel1/ ut1,ut2,ut3,flowtype

      real ut4(lx1,ly1,lz1,lelt)
      real ut5(lx1,ly1,lz1,lelt)
      real ut6(lx1,ly1,lz1,lelt)

      common /testvel2/ ut4,ut5,ut6


      NTOT1 = lx1*ly1*lz1*NELV
      NTOT2 = lx2*ly2*lz2*NELV
      if (igeom.eq.2) CALL LAGVEL

!     prabal
!      call opzero(vx,vy,vz)   ! zero out velocity field
!                              ! (before BCs are applied).

      CALL BCDIRVC (VX,VY,VZ,v1mask,v2mask,v3mask)
      CALL BCNEUTR

      call extrapp (pr,prlag)
      call opgradt (resv1,resv2,resv3,pr)
      CALL OPADD2  (RESV1,RESV2,RESV3,BFX,BFY,BFZ)

      CALL OPHX    (W1,W2,W3,VX,VY,VZ,H1,H2)
      CALL OPSUB2  (RESV1,RESV2,RESV3,W1,W2,W3)


      RETURN
      END
c-----------------------------------------------------------------------


