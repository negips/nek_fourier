!======================================================================
!     Solver Routines to take care of complex variables 
!     in Cylindrical Coordinates.
!     Author: Prabal S. Negi
!
!====================================================================== 
!-----------------------------------------------------------------------

      subroutine ophinv_cyl(or1,or2,or3,oi1,oi2,oi3,
     $                      ir1,ir2,ir3,ii1,ii2,ii3,
     $                      h1,h2,tolh,nmxhi)

C     Ok = (H1*A+H2*B)-1 * Ik  (implicit)

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'SOLN'
      include 'VPROJ'
      include 'TSTEP'

!     Real outputs
      real or1(lx1,ly1,lz1,1),or2(lx1,ly1,lz1,1),or3(lx1,ly1,lz1,1)
!     Imaginary outputs      
      real oi1(lx1,ly1,lz1,1),oi2(lx1,ly1,lz1,1),oi3(lx1,ly1,lz1,1)
     
!     Real Inputs      
      real ir1(lx1,ly1,lz1,1),ir2(lx1,ly1,lz1,1),ir3(lx1,ly1,lz1,1)
!     Imaginary Outputs      
      real ii1(lx1,ly1,lz1,1),ii2(lx1,ly1,lz1,1),ii3(lx1,ly1,lz1,1)
     
      real h1 (lx1,ly1,lz1,1) , h2 (lx1,ly1,lz1,1)

      integer mtmp,i,matmod,nmxhi
      real tolh
 
c      ifproj = .false.
c      if (param(94).gt.0)    ifproj = .true.
c      if (ifprojfld(ifield)) ifproj = .true.
c 
c      if (.not.ifproj) then
c         if (ifield.eq.1) call ophinv
c     $      (o1,o2,o3,i1,i2,i3,h1,h2,tolh,nmxhi)
c         if (ifield.eq.ifldmhd) call ophinv
c     $      (o1,o2,o3,i1,i2,i3,h1,h2,tolh,nmxhi)
c         return
c      endif
 
      mtmp = param(93)
      do i=1,2*ldim
         ivproj(1,i) = min(mxprev,mtmp) - 1
      enddo
 
      imesh = 1
 
      if (ifstrs) then
         matmod = 0
         call hmhzsf_cyl('NOMG',or1,or2,or3,oi1,oi2,oi3,
     $                  ir1,ir2,ir3,ii1,ii2,ii3,h1,h2,
     $                  v1mask,v2mask,v3mask,vmult,
     $                  tolh,nmxhi,matmod)
      else
         if (ifield.eq.1) then
           if (nid.eq.0) then 
             write(6,*) 'IFSTRS = ', ifstrs
             write(6,*) 'Cylindrical Solver needs coupled solve'
             write(6,*) 'Set IFSTRS = TRUE'
           endif  

           call exitt

         elseif (ifield.eq.ifldmhd) then  ! B-field

           if (nid.eq.0) 
     $       write(6,*) 'Cylindrical Solver not implemented for MHD'

           call exitt
         endif
      endif
C
      return
      end subroutine ophinv_cyl 
c--------------------------------------------------------------------

      subroutine hmhzsf_cyl(name,ur1,ur2,ur3,ui1,ui2,ui3,
     $                   rr1,rr2,rr3,ri1,ri2,ri3,h1,h2,
     $                   rmask1,rmask2,rmask3,rmult,
     $                   tol,maxit,matmod)

!     Solve coupled Helmholtz equations (stress formulation)

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'SOLN'   ! For outpost diagnostic call
      include 'TSTEP'
      include 'ORTHOSTRS'
      include 'CTIMER'

      real ur1(1),ur2(1),ur3(1),ui1(1),ui2(1),ui3(1)
      real rr1(1),rr2(1),rr3(1),ri1(1),ri2(1),ri3(1)
      real rmask1(1),rmask2(1),rmask3(1),rmult(1)
      real h1(1),h2(1)
      character name*4

      integer iproj,nel,vol,n,maxit,matmod
      real tol

#ifdef TIMER
      nhmhz = nhmhz + 1
      etime1 = dnekclock()
#endif

      nel = nelfld(ifield)
      vol = volfld(ifield)
      n   = lx1*ly1*lz1*nel

      napproxstrs(1) = 0
      iproj = 0
      if (ifprojfld(ifield)) iproj = param(94)
      if (iproj.gt.0.and.istep.ge.iproj) napproxstrs(1)=param(93)
      napproxstrs(1)=min(napproxstrs(1),mxprev)

      call rmask   (rr1,rr2,rr3,nel)
      call rmask   (ri1,ri2,ri3,nel)
     
      call opdssum (rr1,rr2,rr3)
      call opdssum (ri1,ri2,ri3)

      call rzero3  (ur1,ur2,ur3,n)
      call rzero3  (ui1,ui2,ui3,n)

      if (imesh.eq.1) then
!         call chktcgs (r1,r2,r3,rmask1,rmask2,rmask3,rmult,binvm1
!     $                ,vol,tol,nel)

!         call strs_project_a(r1,r2,r3,h1,h2,rmult,ifield,ierr,matmod)

!         call cggosf  (u1,u2,u3,r1,r2,r3,h1,h2,rmult,binvm1
!     $                ,vol,tol,maxit,matmod)

!         call strs_project_b(u1,u2,u3,h1,h2,rmult,ifield,ierr)

      else

!         call chktcgs (r1,r2,r3,rmask1,rmask2,rmask3,rmult,bintm1
!     $                ,vol,tol,nel)
!         call cggosf  (u1,u2,u3,r1,r2,r3,h1,h2,rmult,bintm1
!     $                ,vol,tol,maxit,matmod)

      endif

#ifdef TIMER
      thmhz=thmhz+(dnekclock()-etime1)
#endif

      return
      end subroutine hmhzsf_cyl
!----------------------------------------------------------------------       

      subroutine cggosf_cyl(u1r,u2r,u3r,u1i,u2i,u3i,
     $                      r1r,r2r,r3r,r1i,r2i,r3i,h1,h2,
     $                      rmult,binv,vol,tin,maxit,matmod)

!     Conjugate gradient iteration for solution of coupled 
!     Helmholtz equations, with complex entires
!     Actual entiries treated with real variables 

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
!      include 'TOTAL'
      include 'MASS'
      include 'DOMAIN'
      include 'FDMH1'

      real dpc,p1r,p2r,p3r
      common /screv/  dpc(lx1*ly1*lz1*lelt)
     $     ,          p1r (lx1*ly1*lz1*lelt)
      common /scrch/  p2r (lx1*ly1*lz1*lelt)
     $     ,          p3r (lx1*ly1*lz1*lelt)

      real qq1,qq2,qq3,pp1r,pp2r,pp3r,wa
      common /scrsl/  qq1(lx1*ly1*lz1*lelt)
     $     ,          qq2(lx1*ly1*lz1*lelt)
     $     ,          qq3(lx1*ly1*lz1*lelt)
      common /scrmg/  pp1r(lx1*ly1*lz1*lelt)
     $     ,          pp2r(lx1*ly1*lz1*lelt)
     $     ,          pp3r(lx1*ly1*lz1*lelt)
     $     ,          wa (lx1*ly1*lz1*lelt)
      real ap1r(1),ap2r(1),ap3r(1)
      real ap1i(1),ap2i(1),ap3i(1)

!     Why do we do this equivalence?
      equivalence (ap1r,pp1r),(ap2r,pp2r),(ap3r,pp3r)
      equivalence (ap1i,pp1i),(ap2i,pp2i),(ap3i,pp3i)

      common /fastmd/ ifdfrm(lelt), iffast(lelt), ifh2, ifsolv
      common /cprint/ ifprint
      logical ifdfrm, iffast, ifh2, ifsolv, ifprint

      real u1r(1),u2r(1),u3r(1),u1i(1),u2i(1),u3i(1),
     $     r1r(1),r2r(1),r3r(1),r1i(1),r2i(1),r3i(1),
     $     h1(1),h2(1),rmult(1),binv(1)

      real p1i,p2i,p3i,pp1i
      common /scruz/ p1i(lx1*ly1*lz1*lelt),
     $               p2i(lx1*ly1*lz1*lelt),
     $               p3i(lx1*ly1*lz1*lelt),
     $               pp1i(lx1*ly1*lz1*lelt)

      real pp2i,pp3i
      common /scrvh/ pp2i(lx1*ly1*lz1*lelt),
     $               pp3i(lx1*ly1*lz1*lelt)


      integer maxit,matmod,nel,nxyz,n,iter,ifin
      logical iffdm,ifcrsl

      real tin,tol,vol
      real alphar,alphai,betar,betai,papr,papi
      real rpp1r,rpp1i,rpp2r,rpp2i

      real r0,rbnorm

      real glsc2,glsc3,op_glsc2_wt
      real opnorm2_wt_comp

!     This common block is used in weak laplacian calculations
!     In this routine this is only used as a work array
!     Since the entries are not required after the weak laplacian
!     has been evaluated      
      real wk1r,wk2r,wk3r,wk1i,wk2i,wk3i,wk4
      common /scrns/   wk1r(lx1*ly1*lz1*lelt)
     $               , wk2r(lx1*ly1*lz1*lelt)
     $               , wk3r(lx1*ly1*lz1*lelt)
     $               , wk1i(lx1*ly1*lz1*lelt)
     $               , wk2i(lx1*ly1*lz1*lelt)
     $               , wk3i(lx1*ly1*lz1*lelt)
     $               , wk4(lx1*ly1*lz1*lelt) 




!     No Fast Diagonalization Method            
      iffdm  = .false.
!     No Coarse grid
      ifcrsl = .false.

      nel   = nelfld(ifield)
      nxyz  = lx1*ly1*lz1
      n     = nxyz*nel

      if (istep.le.1.and.iffdm) call set_fdm_prec_h1A

      tol  = tin

c     overrule input tolerance
      if (restol(ifield).ne.0) tol=restol(ifield)
      if (ifcrsl) call set_up_h1_crs_strs(h1,h2,ifield,matmod)

      if ( .not.ifsolv ) then           !     Set logical flags
         call setfast (h1,h2,imesh)
         ifsolv = .true.
      endif

!      call opdot (wa,r1r,r2r,r3r,r1r,r2r,r3r,n)
!      rbnorm = glsc3(wa,binv,rmult,n)
!      rbnorm = sqrt ( rbnorm / vol )
      call col3(wa,binv,rmult,n)
      rbnorm  = opnorm2_wt_comp(r1r,r2r,r3r,r1i,r2i,r3i,wa,n)
      rbnorm = sqrt(rbnorm/vol) 

      if (rbnorm .lt. tol**2) then
         iter = 0
         r0 = rbnorm
c        if ( .not.ifprint )  goto 9999
         if (matmod.ge.0.and.nio.eq.0) write (6,3000) 
     $                                 istep,iter,rbnorm,r0,tol
         if (matmod.lt.0.and.nio.eq.0) write (6,3010) 
     $                                 istep,iter,rbnorm,r0,tol
         goto 9999
      endif

!!     prabal. Commenting for now      
!!     Evaluate diagional pre-conidtioner for fluid solve
!      call setprec (dpc,h1,h2,imesh,1)
!      call setprec (wa ,h1,h2,imesh,2)
!      call add2    (dpc,wa,n)
!
!      if (ldim.eq.3) then
!         call setprec (wa,h1,h2,imesh,3)
!         call add2    (dpc,wa,n)
!      endif
      call rone(dpc,n)        ! prabal. remove this

      if (iffdm) then
!         call set_fdm_prec_h1b(dpc,h1,h2,nel)
!         call fdm_h1a (pp1,r1,dpc,nel,ktype(1,1,1),wa)
!         call fdm_h1a (pp2,r2,dpc,nel,ktype(1,1,2),wa)
!         call fdm_h1a (pp3,r3,dpc,nel,ktype(1,1,3),wa)
!         call rmask   (pp1,pp2,pp3,nel)
!         call opdssum (pp1,pp2,pp3)
      else
!!        Real           
!         call col3 (pp1,dpc,r1r,n)
!         call col3 (pp2,dpc,r2r,n)
!         call col3 (pp3,dpc,r3r,n)
!
!!        Imaginary         
!         call col3 (pp1,dpc,r1i,n)
!         call col3 (pp2,dpc,r2i,n)
!         call col3 (pp3,dpc,r3i,n)

        call copy3(pp1r,pp2r,pp3r,r1r,r2r,r3r,n)
        call copy3(pp1i,pp2i,pp3i,r1i,r2i,r3i,n)
      endif

      if (ifcrsl) then
!         call crs_strs(p1,p2,p3,r1,r2,r3)
!         call rmask   (p1,p2,p3,nel)
      else
         call rzero3(p1r,p2r,p3r,n)
         call rzero3(p1i,p2i,p3i,n)
      endif

      call opadd2(p1r,p2r,p3r,pp1r,pp2r,pp3r)
      call opadd2(p1i,p2i,p3i,pp1i,pp2i,pp3i)

!      rpp1 = op_glsc2_wt(p1,p2,p3,r1r,r2r,r3r,rmult)
      call opglsc2_wt_comp(rpp1r,rpp1i,p1r,p2r,p3r,p1i,p2i,p3i,
     $                        r1r,r2r,r3r,r1i,r2i,r3i,rmult,n)

      maxit=200
      do 1000 iter=1,maxit
!         call axhmsf_cyl(ap1,ap2,ap3,p1,p2,p3,h1,h2,matmod)
         call axhmsf_cyl(Ap1r,Ap2r,Ap3r,Ap1i,Ap2i,Ap3i,
     $                      p1r,p2r,p3r,p1i,p2i,p3i,h1,h2)

!        prabal. Need to change rmask         
         call rmask   (ap1r,ap2r,ap3r,nel)
         call rmask   (ap1i,ap2i,ap3i,nel)

         call dssum3(Ap1r,Ap2r,Ap3r)
         call dssum3(Ap1i,Ap2i,Ap3i)

!         pap   = op_glsc2_wt(p1,p2,p3,ap1,ap2,ap3,rmult)

!        pAp = pn-1*Apn-1            
         call opglsc2_wt_comp(pApr,pApi,p1r,p2r,p3r,p1i,p2i,p3i,
     $                        Ap1r,Ap2r,Ap3r,Ap1i,Ap2i,Ap3i,rmult,n)

!        \alpha = (rn-1*rn-1)/(pn-1*Apn-1) 
         alphar = (rpp1r*pApr + rpp1i*pApi)/(pApr**2 + pApi**2)
         alphai = (rpp1i*pApr - rpp1r*pApi)/(pApr**2 + pApi**2)

!        ur = ur + \alphar*pr - \alphai*pi
         call opadds3(u1r,u2r,u3r,p1r,p2r,p3r,alphar,n,2)
         call opadds3(u1r,u2r,u3r,p1r,p2r,p3r,-alphai,n,2)

!        ui = ui + \alphar*pi + \alphai*pr
         call opadds3(u1i,u2i,u3i,p1i,p2i,p3i,alphar,n,2)
         call opadds3(u1i,u2i,u3i,p1r,p2r,p3r,alphai,n,2)

!        rr = rr - \alphar*Apr + \alphai*Api
         call opadds3(r1r,r2r,r3r,Ap1r,Ap2r,Ap3r,-alphar,n,2)
         call opadds3(r1r,r2r,r3r,Ap1i,Ap2i,Ap3i,alphai,n,2)

!        ri = ri - \alphar*Api - \alphai*Apr
         call opadds3(r1i,r2i,r3i,Ap1i,Ap2i,Ap3i,-alphar,n,2)
         call opadds3(r1i,r2i,r3i,Ap1r,Ap2r,Ap3r,-alphai,n,2)


!         call opdot  (wa,r1,r2,r3,r1,r2,r3,n)
         call col3(wa,binv,rmult,n)
         rbnorm = opnorm2_wt_comp(r1r,r2r,r3r,r1i,r2i,r3i,wa,n)
         rbnorm = sqrt (rbnorm/vol)

         if (iter.eq.1) r0 = rbnorm

         if (rbnorm.lt.tol) then
            ifin = iter
            if (nio.eq.0) then
               if (matmod.ge.0) write(6,3000) istep,ifin,rbnorm,r0,tol
               if (matmod.lt.0) write(6,3010) istep,ifin,rbnorm,r0,tol
            endif
            goto 9999
         endif

         if (iffdm) then
!            call fdm_h1a (pp1,r1,dpc,nel,ktype(1,1,1),wa)
!            call fdm_h1a (pp2,r2,dpc,nel,ktype(1,1,2),wa)
!            call fdm_h1a (pp3,r3,dpc,nel,ktype(1,1,3),wa)
!            call rmask   (pp1,pp2,pp3,nel)
!            call opdssum (pp1,pp2,pp3)
         else
!            call col3 (pp1,dpc,r1r,n)
!            call col3 (pp2,dpc,r2r,n)
!            if (if3d) call col3 (pp3,dpc,r3r,n)

!           prabal. preconditioner application 
            call copy3(pp1r,pp2r,pp3r,r1r,r2r,r3r,n)
            call copy3(pp1i,pp2i,pp3i,r1i,r2i,r3i,n)
           
         endif

         if (ifcrsl) then
!           call crs_strs(qq1,qq2,qq3,r1,r2,r3)
!           call rmask   (qq1,qq2,qq3,nel)
!           call opadd2  (pp1,pp2,pp3,qq1,qq2,qq3)
         endif


         call opglsc2_wt_comp(rpp1r,rpp1i,p1r,p2r,p3r,p1i,p2i,p3i,
     $                        r1r,r2r,r3r,r1i,r2i,r3i,rmult,n)

         rpp2r = rpp1r
         rpp2i = rpp1i

!        \beta = (rn*rn)/(rn-1*rn-1)
         betar = rpp1r*rpp2r + rpp1i*rpp2i
         betai = rpp1i*rpp2r - rpp1r*rpp2i


         call copy3(wk1r,wk2r,wk3r,p1r,p2r,p3r,n)
         call copy3(wk1i,wk2i,wk3i,p1i,p2i,p3i,n)
         
!        pr = ppr + \betar*pr - \betai*pi
         call copy3(p1r,p2r,p3r,pp1r,pp2r,pp3r,n)
         call opadds3(p1r,p2r,p3r,wk1r,wk2r,wk3r,betar,n,2)   
         call opadds3(p1r,p2r,p3r,wk1i,wk2i,wk3i,-betai,n,2)   

!        pi = ppi + \betar*pi + \betai*pr
         call copy3(p1i,p2i,p3i,pp1i,pp2i,pp3i,n)
         call opadds3(p1i,p2i,p3i,wk1i,wk2i,wk3i,betar,n,2)   
         call opadds3(p1i,p2i,p3i,wk1r,wk2r,wk3r,betai,n,2)   

 1000 continue
      if (matmod.ge.0.and.nio.eq.0) write (6,3001) 
     $                              istep,iter,rbnorm,r0,tol
      if (matmod.lt.0.and.nio.eq.0) write (6,3011) 
     $                              istep,iter,rbnorm,r0,tol

 9999 continue
      ifsolv = .false.


 3000 format(i11,'  Helmh3 fluid  ',I6,1p3E13.4)
 3010 format(i11,'  Helmh3 mesh   ',I6,1p3E13.4)
 3001 format(i11,'  Helmh3 fluid unconverged! ',I6,1p3E13.4)
 3011 format(i11,'  Helmh3 mesh unconverged! ',I6,1p3E13.4)

      return
      end subroutine cggosf_cyl
!-----------------------------------------------------------------------

      real function opnorm2_wt_comp(u1r,u2r,u3r,u1i,u2i,u3i,wt,n)

!                _       
!     Calculate: u *Wt* u
!     Assuming the Weight itself is real        

      implicit none

!      real opnorm2_wt_comp

      integer i,n

      real u1r(n),u2r(n),u3r(n),u1i(n),u2i(n),u3i(n)
      real wt(n)

      real a1,a2,wk

      a1 = 0.
      a2 = 0.
      do i=1,n
        a1 = a1 + wt(i)*(u1r(i)*u1r(i) + u2r(i)*u2r(i) + u3r(i)*u3r(i))
        a2 = a2 + wt(i)*(u1i(i)*u1i(i) + u2i(i)*u2i(i) + u3i(i)*u3i(i))
      enddo
     
      wk = a1 + a2

      call glsum(wk,1)

      opnorm2_wt_comp = wk

      return
      end function opnorm2_wt_comp

!-----------------------------------------------------------------------

      subroutine opglsc2_wt_comp(scr,sci,u1r,u2r,u3r,u1i,u2i,u3i,
     $                              v1r,v2r,v3r,v1i,v2i,v3i,wt,n)

!                _       
!     Calculate: u *Wt* v
!     Assuming the Weight itself is real        

      implicit none

      integer i,n

      real u1r(n),u2r(n),u3r(n),u1i(n),u2i(n),u3i(n)
      real v1r(n),v2r(n),v3r(n),v1i(n),v2i(n),v3i(n)
      real wt(n)

      real scr,sci,a1,a2

      a1 = 0.
      a2 = 0.
      do i=1,n
        a1 = a1 + wt(i)*(u1r(i)*v1r(i) + u2r(i)*v2r(i) + u3r(i)*v3r(i))
        a1 = a1 + wt(i)*(u1i(i)*v1i(i) + u2i(i)*v2i(i) + u3i(i)*v3i(i))

        a2 = a2 + wt(i)*(u1r(i)*v1i(i) + u2r(i)*v2i(i) + u3r(i)*v3i(i))
        a2 = a2 - wt(i)*(u1i(i)*v1r(i) + u2i(i)*v2r(i) + u3i(i)*v3r(i))
      enddo

!     a1 = ur*Wt*vr + ui*Wt*vi
!     a2 = ur*Wt*vi - ui*Wt*vr
!     sc = a1 + i*a2      

      scr = a1
      sci = a2

!     Sum over all processors 
      call glsum(scr,1)
      call glsum(sci,1)


      return
      end subroutine opglsc2_wt_comp

!-----------------------------------------------------------------------

      subroutine copy3(a1,a2,a3,b1,b2,b3,n)

      implicit none  

      integer i,n

      real a1(1),a2(1),a3(1)
      real b1(1),b2(1),b3(1)

      do i=1,n
         a1(i)=b1(i)
         a2(i)=b2(i)
         a3(i)=b3(i)
      enddo


      return
      end subroutine copy3
!-----------------------------------------------------------------------

      subroutine dssum3(a1,a2,a3)

      implicit none

      include 'SIZE'

      real a1(1),a2(1),a3(1)

      call dssum(a1,nx1,ny1,nz1)
      call dssum(a2,nx1,ny1,nz1)      
      call dssum(a3,nx1,ny1,nz1)

      return  
      end subroutine dssum3        
!-----------------------------------------------------------------------

      subroutine opadds3 (a1,a2,a3,b1,b2,b3,const,n,isc)

      implicit none

      integer n,isc
      real const

      real a1(1),a2(1),a3(1),b1(1),b2(1),b3(1)

      if (isc.eq.1) then
         call add2s1 (a1,b1,const,n)
         call add2s1 (a2,b2,const,n)
         call add2s1 (a3,b3,const,n)
      elseif (isc.eq.2) then
         call add2s2 (a1,b1,const,n)
         call add2s2 (a2,b2,const,n)
         call add2s2 (a3,b3,const,n)
      endif

      return
      end subroutine opadds3
!-----------------------------------------------------------------------








