!======================================================================
!     Routines for introducing Fourier basis in the 3rd direction
!     Author: Prabal S. Negi
!
!====================================================================== 
!-----------------------------------------------------------------------
      subroutine init_fou()

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'

      include '3DS'

      integer icalld
      save icalld
      data icalld /0/

      integer nxyz,ntot1
      

      iffou = .true.

      if (.not.iffou) return

      if (nio.eq.0) then
        write(6,*) 'Initializing zFourier', iffou
      endif

      if (iffou.and.ndim.eq.3) then
        write(6,'(A6,1x,I2)') 'Ndim =', ndim
        write(6,*) 'Ensure ndim = 2'
        call exitt
      endif

      nxyz  = lx1*ly1*lz1
      ntot1 = nxyz*nelv

!     Initialize forward/backward ffts
      call init_plans_fou()

!     Initial field     
      call init_fld_fou()

!     Need to initialize some variables
!     V3MASK
      call copy(v3mask,v1mask,ntot1)

      if (nio.eq.0) then
        write(6,*) 'Initializing zfourier', iffou
      endif

      return
      end subroutine init_fou
!----------------------------------------------------------------------

      subroutine init_plans_fou

      implicit none

      include 'fftw3.f'

      include 'SIZE'
      include 'FOURIER'
      
      integer datalen

!     FFTW_MEASURE  - Instructs FFTW to run and measure the execution time of several FFTs
!                     in order to find the best way to compute the transform of size n. 
!                     This process takes some time (usually a few seconds), depending on the 
!                     machine and on the size of the transform.
!     FFTW_ESTIMATE - Does not run any computation and just builds a reasonable plan for the transform,
!                     which is probably sub optimal.
!     Other options with longer initialization times are
!     FFTW_PATIENT
!     FFTW_EXHAUSTIVE 



      datalen = nfmodes
     
!     Forward plan (r --> c) 
      call dfftw_plan_dft_r2c_1d(fou_planf,datalen,fou_datain,
     $                           fou_fftamp,FFTW_ESHAUSTIVE)


      datalen = nfmodes
!     Backward plan (c --> r)
      call dfftw_plan_dft_c2r_1d(fou_planb,datalen,fou_fftamp,
     $                           fou_datain,FFTW_EXHAUSTIVE)


!     Dealiased Plans
      datalen = nfmodes_d

!     Forward plan (r --> c) 
      call dfftw_plan_dft_r2c_1d(fou_planf_d,datalen,fou_vxd,
     $                           fou_amp_vxd,FFTW_ESHAUSTIVE)


      datalen = nfmodes_d
!     Backward plan (c --> r)
      call dfftw_plan_dft_c2r_1d(fou_planb_d,datalen,fou_amp_vxd,
     $                           fou_vxd,FFTW_EXHAUSTIVE)

      return
      end subroutine init_plans_fou
!-----------------------------------------------------------------------

      subroutine fluid_fou (igeom)

      implicit none

      include 'SIZE'
      include 'DEALIAS'
      include 'INPUT'
      include 'SOLN'
      include 'TSTEP'

      include 'FOURIER'

      real*8 ts, dnekclock

      if (.not.iffou) return 
 
      ifield = 1
      imesh  = 1
      call unorm
      call settolv

      ts = dnekclock() 

      if(nio.eq.0 .and. igeom.eq.2) 
     &   write(*,'(13x,a)') 'Solving for fluid'

      if (ifsplit) then
        write(6,*) 'Fourier basis not implemented for Pn-Pn'
        call exitt    
        
      elseif (iftran) then

        call plan3_fou(igeom)    !  Same as PLAN 1 w/o nested iteration
                                 !  Std. NEKTON time stepper  !

!         if (igeom.ge.2) call chkptol         ! check pressure tolerance
         if (igeom.eq.ngeom) then 
           if (ifneknekc) then
!              call vol_flow_ms    ! check for fixed flow rate
           else
              call vol_flow       ! check for fixed flow rate
           endif
         endif

      else   !  steady Stokes, non-split

!        call plan1 (igeom) ! The NEKTON "Classic".

        write(6,*) 'Fourier basis not implemented for Steady-State yet'
        call exitt    

      endif

      if(nio.eq.0 .and. igeom.ge.2) 
     &   write(*,'(4x,i7,a,1p2e12.4)') 
     &   istep,'  Fluid_fou done',time,dnekclock()-ts

      return
      end subroutine fluid_fou
c-----------------------------------------------------------------------

      subroutine plan3_3ds (igeom)

!     Compute pressure and velocity using consistent approximation spaces.     
!     Operator splitting technique.

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'EIGEN'
      include 'SOLN'
      include 'TSTEP'

      include 'FOURIER'

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

         call makef_fou

      else

!        new geometry, new b.c.

         intype = -1
         call sethlm  (h1,h2,intype)

         call cresvif_3ds (resv1,resv2,resv3,h1,h2)

!         debugging  
!         call opcopy(vx,vy,vz,resv1,resv2,resv3)
!         call copy(vz,resv3,ntot1)      
!         call outpost(vx,vy,vz,pr,vz,'   ')
!         call exitt

!         debugging  
!         call opcopy(vx,vy,vz,bfx,bfy,bfz)
!         call copy(vz,bfz,ntot1)      
!         call outpost(vx,vy,vz,pr,vz,'   ')
!         call exitt

!         call copy(resv1,resv3,ntot1)    
 
     
         if3d = .true. 
         call ophinv(dv1,dv2,dv3,resv1,resv2,resv3,h1,h2,tolhv,nmxv)
         if3d = .false.   


         call opadd2(vx,vy,vz,dv1,dv2,dv3)
         call add2(vz,dv3,ntot1) 

!!        prabal
!         call opadd2(vx,vy,vz,ut1,ut2,ut3)            ! add slip velocity back

!         debugging  
!         call opcopy(vx,vy,vz,dv1,dv2,dv3)
!         call copy(vz,dv3,ntot1)      
!         call outpost(vx,vy,vz,pr,vz,'   ')   
!         call exitt   

         call incomprn(vx,vy,vz,pr)

      endif

      return
      end subroutine plan3_fou

!----------------------------------------------------------------------

      subroutine makef_fou

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'INPUT'
      include 'TSTEP'
      include 'CTIMER'
      include 'MVGEOM'

      integer i
      real dnekclock
      real etime1

      etime1 = dnekclock()

!     Phsical space Operations
      call makeuf_fou
      call advab_fou 

!     Perform FFT on the dealiased grid


!     Fourier space operations

      do k=1,nfmodes      
        call makeabf_fou(k)
        call makebdf_fou(k)
      enddo

c     Adding this call allows prescribed pressure bc for PnPn-2
c     if (.not.ifsplit.and..not.ifstrs)         call bcneutr
      
      tmakf=tmakf+(dnekclock()-etime1)

      return
      end subroutine makef_fou

!---------------------------------------------------------------------- 

      subroutine fluidp_3ds (igeom)

!     Driver for perturbation velocity

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'SOLN'

      integer igeom
      integer jp0,jp2
      integer i,j

      jp2 = npert/2

      jp = 0
      do i = 1,jp2
        jp0 = jp    
        do j = 1,2
          jp = jp0 + j    
          if (nio.eq.0.and.igeom.eq.2) write(6,1) istep,time,jp
   1      format(i9,1pe14.7,' Perturbation Solve (Momentum):',i5)

          call perturbv_mom_3ds (igeom)

!!         prabal  
!          if (igeom.eq.1)  
!     $      call outpost(bfxp(1,jp),bfyp(1,jp),bfzp(1,jp),
!     $                  prp(1,jp),bfzp(1,jp),'bfp')

!!         prabal  
!          if (igeom.eq.1)  
!     $      call outpost(vx,vy,vz,pr,vz,'ch1')


        enddo

        do j=1,2
          jp = jp0 + j 
          if (nio.eq.0.and.igeom.eq.2) write(6,2) istep,time,jp
   2      format(i9,1pe14.7,' Perturbation Solve (Pressure):',i5)

!          call perturbv_prp_3ds ()
           call incomprp_3ds(igeom) 

!!         prabal  
!          if (igeom.eq.2)  
!     $      call outpost(vx,vy,vz,pr,vz,'ch2')

        enddo

        do j = 1,2
          jp = jp0 + j    
          call velp_update_3ds (igeom)

!!         prabal  
!          if (igeom.eq.2)  
!     $      call outpost(vx,vy,vz,pr,vz,'ch3')

        enddo
      enddo ! i=1,jp2


      jp=0   ! set jp to zero, for baseline flow

      return
      end subroutine fluidp_3ds
c-----------------------------------------------------------------------

      subroutine perturbv_mom_3ds (igeom)

      implicit none

!     Solve the convection-diffusion equation for the perturbation field, 
!     with projection onto a div-free space.


      include 'SIZE'
      include 'INPUT'
      include 'EIGEN'
      include 'SOLN'
      include 'TSTEP'
      include 'MASS'

      include 'TEST'

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

      integer intype,igeom
      integer ntot1


      ifield = 1

      ntot1 = lx1*ly1*lz1*nelv

      if (igeom.eq.1) then

!        Old geometry, old velocity

         call makefp_3ds
         call lagfieldp_3ds

!        Add third component of convective term   
!         call advab_w_3ds   

      else
c
c        New geometry, new velocity
c
         intype = -1
         call sethlm_3dsp(h1,h2,intype)
         call cresvipp_3ds(resv1,resv2,resv3,h1,h2)

         if3d = .true.   
         call ophinv   (dv1,dv2,dv3,resv1,resv2,resv3,h1,h2,tolhv,nmxv)
         if3d = .false.

         call add2(vxp(1,jp),dv1,ntot1)
         call add2(vyp(1,jp),dv2,ntot1)
         call add2(vzp(1,jp),dv3,ntot1)

      endif

      return
      end subroutine perturbv_mom_3ds
!-----------------------------------------------------------------------

!      subroutine perturbv_prp_3ds (igeom)
!
!      implicit none
!
!!     Solve the convection-diffusion equation for the perturbation field, 
!!     with projection onto a div-free space.
!
!
!      include 'SIZE'
!      include 'INPUT'
!      include 'EIGEN'
!      include 'SOLN'
!      include 'TSTEP'
!      include 'MASS'
!
!      include 'TEST'
!
!
!      ifield = 1
!
!      call incomprp_3ds ()
!
!
!      return
!      end subroutine perturbv_prp_3ds


!-----------------------------------------------------------------------
      subroutine init_fld_fou

!     Initialize starting field.
!     Assume its initialized in Physical space

      implicit none

      include 'SIZE'
      include 'SOLN'    ! jp
      include 'TSTEP'   ! ifield

      include 'FOURIER'

      include 'INPUT'
      include 'TSTEP'
      include 'PARALLEL'
      include 'NEKUSE'

      integer i,j,k
      integer e,eg,nz
      integer ntot1

      nel   = nelfld(ifield)

      ntot1 = nx1*ny1*nz1*nelv

      nz    = nfmodes

      do k=1,nz
        do e=1,nel
          eg = lglel(e)
          do j=1,ly1
            do i=1,lx1
              call nekasgn (i,j,1,e)
              z = fou_zm1(k) 
              call useric  (i,j,k,eg)
              if (jp.eq.0) then
                if (ifield.eq.1) then
                  fou_icvx(i,j,e,k) = ux
                  fou_icvy(i,j,e,k) = uy
                  fou_icvz(i,j,e,k) = uz
                elseif (ifield.eq.ifldmhd .and. ifmhd) then
!!                   prabal. MHD not implemented yet
                else
!                   prabal. Passive scalar not implemented yet   
!                    t(i,j,k,e,ifield-1) = temp
                endif
              else
!                 prabal. We'll come to perturbation mode some other time
              endif

            enddo     ! i
          enddo       ! j
        enddo         ! e
        call dsavg(fou_icvx)
        call dsavg(fou_icvy)
        call dsavg(fou_icvz)
      enddo           ! k

!     prabal. Perform a fourier transform
      if (jp.eq.0) then
        if (ifield.eq.1) then

          call phy_to_fou(fou_rvx,fou_ivx,fou_icvx,fou_datain,
     $                    fou_fftamp,fou_planf,nz)  

          call phy_to_fou(fou_rvy,fou_ivy,fou_icvy,fou_datain,
     $                    fou_fftamp,fou_planf,nz)  

          call phy_to_fou(fou_rvz,fou_ivz,fou_icvz,fou_datain,
     $                    fou_fftamp,fou_planf,nz)  

        elseif (ifield.eq.ifldmhd .and. ifmhd) then
!!           prabal. MHD not implemented yet
        else
!           prabal. Passive scalar not implemented yet   
!            t(i,j,k,e,ifield-1) = temp
        endif
      else
!         prabal. We'll come to perturbation mode some other time
      endif


      return
      end subroutine init_fld_fou
!----------------------------------------------------------------------
      subroutine makefp_3ds

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MASS'
      include 'TSTEP'   ! ifield

      include '3DS'


      integer nxyz,ntot1

      real ta1,ta2,ta3
      common /scruz/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)


      nxyz = lx1*ly1*lz1
      ntot1 = nxyz*nelv

!     Build user defined forcing
      call makeufp_3ds


!      if3d = .true.
!      if (filterType.eq.2) call make_hpf
!     hpf field stored in ta3
!      call xaddcol3(bfz_3ds,ta3,bm1,ntot1)
!      if3d = .false. 
!          call outpost(bfxp(1,jp),bfyp(1,jp),bfzp(1,jp),
!     $                  prp(1,jp),bfzp(1,jp),'bf1')

      call advabp_3ds
!      if (ifnav.and.(.not.ifchar).and.(ifadj)) call advabp_adjoint_3dsp


      if (iftran) call makextp_3ds
      call makebdfp_3ds


      return
      end subroutine makefp_3ds

!----------------------------------------------------------------------

      subroutine makeufp_3ds

!     Compute and add: (1) user specified forcing function (FX,FY,FZ)

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
      include 'PARALLEL'
      include 'NEKUSE'

      include '3DS'

      real ta1,ta2,ta3
      common /scruz/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)

      integer ntot1,iel,i,j,k,ielg
      integer ijke


      ntot1 = lx1*ly1*lz1*nelv

      time = time-dt
      call rzero(bfxp(1,jp),ntot1)
      call rzero(bfyp(1,jp),ntot1)
      call rzero(bfzp(1,jp),ntot1)

      do 100 iel=1,nelv
         ielg = lglel(iel)
         do 100 k=1,lz1
         do 100 j=1,ly1
         do 100 i=1,lx1
            call nekasgn (i,j,k,iel)
            call userf   (i,j,k,ielg)
            ijke = i+lx1*((j-1)+ly1*((k-1) + lz1*(iel-1)))
            bfxp(ijke,jp) = ffx
            bfyp(ijke,jp) = ffy
            bfzp(ijke,jp) = ffz
 100  continue

!     Not sure why we multiply by density
!      call col2(bfzp(1,jp),vtrans(1,1,1,1,ifield),nx1*ny1*nz1*nelv)

      call col2  (bfxp(1,jp),bm1,ntot1)
      call col2  (bfyp(1,jp),bm1,ntot1)
      call col2  (bfzp(1,jp),bm1,ntot1)
      time = time+dt

      return
      end subroutine makeufp_3ds

!-----------------------------------------------------------------------
      subroutine advabp_3ds

!     Eulerian scheme, add convection term to forcing function
!     at current time step.

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'

      include '3DS'

      real ta1,ta2,ta3
      real tb1,tb2,tb3
      common /scrns/ ta1 (lx1*ly1*lz1*lelv)
     $ ,             ta2 (lx1*ly1*lz1*lelv)
     $ ,             ta3 (lx1*ly1*lz1*lelv)
     $ ,             tb1 (lx1*ly1*lz1*lelv)
     $ ,             tb2 (lx1*ly1*lz1*lelv)
     $ ,             tb3 (lx1*ly1*lz1*lelv)

      integer i,ntot1
      real tmp


      ntot1 = lx1*ly1*lz1*nelv

      if (if3d.or.if3d_3ds) then
         call copy  (tb1,vx,ntot1)                   ! Save velocity
         call copy  (tb2,vy,ntot1)                   ! Save velocity
         call copy  (tb3,vz,ntot1)                   ! Save velocity

!        U <-- dU
         call copy  (vx,vxp(1,jp),ntot1)                   ! Save velocity
         call copy  (vy,vxp(1,jp),ntot1)                   ! Save velocity
         call copy  (vz,vxp(1,jp),ntot1)                   ! Save velocity

         call convop  (ta1,tb1)                                ! du.grad U
         call convop  (ta2,tb2)
         call convop  (ta3,tb3)

!        Restore velocity
         call copy  (vx,tb1,ntot1)
         call copy  (vy,tb2,ntot1)
         call copy  (vz,tb3,ntot1)

         do i=1,ntot1
            tmp = bm1(i,1,1,1)*vtrans(i,1,1,1,ifield)
            bfxp(i,jp) = bfxp(i,jp)-tmp*ta1(i)
            bfyp(i,jp) = bfyp(i,jp)-tmp*ta2(i)
            bfzp(i,jp) = bfzp(i,jp)-tmp*ta3(i)
         enddo

         call convop  (ta1,vxp(1,jp))       !  U.grad dU
         call convop  (ta2,vyp(1,jp))
         call convop  (ta3,vzp(1,jp))

         do i=1,ntot1
            tmp = bm1(i,1,1,1)*vtrans(i,1,1,1,ifield)
            bfxp(i,jp) = bfxp(i,jp)-tmp*ta1(i)
            bfyp(i,jp) = bfyp(i,jp)-tmp*ta2(i)
            bfzp(i,jp) = bfzp(i,jp)-tmp*ta3(i)
         enddo

!        Add z convection for all components
!        Assuming dU/dz, dV/dz, dW/dz = 0 of course
         if (mod(jp,2).eq.1) then
            call convect_w_3ds(ta1,vxp(1,jp+1),vz)           
            call convect_w_3ds(ta2,vyp(1,jp+1),vz)           
            call convect_w_3ds(ta3,vzp(1,jp+1),vz)           

            do i=1,ntot1
               tmp = -k_3dsp*bm1(i,1,1,1)*vtrans(i,1,1,1,ifield)
               bfxp(i,jp) = bfxp(i,jp)-tmp*ta1(i)
               bfyp(i,jp) = bfyp(i,jp)-tmp*ta2(i)
               bfzp(i,jp) = bfzp(i,jp)-tmp*ta3(i)
            enddo

         else
            call convect_w_3ds(ta1,vxp(1,jp-1),vz)           
            call convect_w_3ds(ta2,vyp(1,jp-1),vz)           
            call convect_w_3ds(ta3,vzp(1,jp-1),vz)           

            do i=1,ntot1
               tmp = k_3dsp*bm1(i,1,1,1)*vtrans(i,1,1,1,ifield)
               bfxp(i,jp) = bfxp(i,jp)-tmp*ta1(i)
               bfyp(i,jp) = bfyp(i,jp)-tmp*ta2(i)
               bfzp(i,jp) = bfzp(i,jp)-tmp*ta3(i)
            enddo

         endif

      else  ! 2D without fourier third component

         call opcopy  (tb1,tb2,tb3,vx,vy,vz)                   ! Save velocity
         call opcopy  (vx,vy,vz,vxp(1,jp),vyp(1,jp),vzp(1,jp)) ! U <-- dU
         call convop  (ta1,tb1)                                ! du.grad U
         call convop  (ta2,tb2)
         call opcopy  (vx,vy,vz,tb1,tb2,tb3)  ! Restore velocity

         do i=1,ntot1
            tmp = bm1(i,1,1,1)*vtrans(i,1,1,1,ifield)
            bfxp(i,jp) = bfxp(i,jp)-tmp*ta1(i)
            bfyp(i,jp) = bfyp(i,jp)-tmp*ta2(i)
         enddo

         call convop  (ta1,vxp(1,jp))       !  U.grad dU
         call convop  (ta2,vyp(1,jp))

         do i=1,ntot1
            tmp = bm1(i,1,1,1)*vtrans(i,1,1,1,ifield)
            bfxp(i,jp) = bfxp(i,jp)-tmp*ta1(i)
            bfyp(i,jp) = bfyp(i,jp)-tmp*ta2(i)
         enddo

      endif


      return
      end subroutine advabp_3ds
c--------------------------------------------------------------------


      subroutine makextp_3ds

!     Add extrapolation terms to perturbation source terms

!     (nek5 equivalent for velocity is "makeabf")

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'

      include '3DS'

      common /scrns/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)

      ntot1 = lx1*ly1*lz1*nelv

      ab0 = ab(1)
      ab1 = ab(2)
      ab2 = ab(3)
      call add3s2 (ta1,exx1p(1,jp),exx2p(1,jp),ab1,ab2,ntot1)
      call add3s2 (ta2,exy1p(1,jp),exy2p(1,jp),ab1,ab2,ntot1)
      call copy   (exx2p(1,jp),exx1p(1,jp),ntot1)
      call copy   (exy2p(1,jp),exy1p(1,jp),ntot1)
      call copy   (exx1p(1,jp),bfxp (1,jp),ntot1)
      call copy   (exy1p(1,jp),bfyp (1,jp),ntot1)
      call add2s1 (bfxp(1,jp),ta1,ab0,ntot1)
      call add2s1 (bfyp(1,jp),ta2,ab0,ntot1)
      if (if3d.or.if3d_3ds) then
         call add3s2 (ta3,exz1p(1,jp),exz2p(1,jp),ab1,ab2,ntot1)
         call copy   (exz2p(1,jp),exz1p(1,jp),ntot1)
         call copy   (exz1p(1,jp),bfzp (1,jp),ntot1)
         call add2s1 (bfzp(1,jp),ta3,ab0,ntot1)
      endif
c
      return
      end subroutine makextp_3ds
!-----------------------------------------------------------------------

      subroutine makebdfp_3ds

!     Add contributions to perturbation source from lagged BD terms.

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'GEOM'
      include 'INPUT'
      include 'TSTEP'

      include '3DS'

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
      call cmult2(h2,vtrans(1,1,1,1,ifield),const,ntot1)
      call opcolv3c (tb1,tb2,tb3
     $              ,vxp(1,jp),vyp(1,jp),vzp(1,jp),bm1,bd(2))

      if (if3d_3ds) then
        call col3(tb3,vzp(1,jp),bm1,ntot1)
        call cmult(tb3,bd(2),ntot1)
      endif

!     Add contribution from lag terms
      do ilag=2,nbd
         if (ifgeom) then
            call opcolv3c(ta1,ta2,ta3,vxlagp(1,ilag-1,jp),
     $                                vylagp(1,ilag-1,jp),
     $                                vzlagp(1,ilag-1,jp),
     $                                bm1lag(1,1,1,1,ilag-1),bd(ilag+1))
         else
            call opcolv3c(ta1,ta2,ta3,vxlagp(1,ilag-1,jp),
     $                                vylagp(1,ilag-1,jp),
     $                                vzlagp(1,ilag-1,jp),
     $                                bm1                   ,bd(ilag+1))
         endif
         call opadd2  (tb1,tb2,tb3,ta1,ta2,ta3)
      enddo
      call opadd2col (bfxp(1,jp),bfyp(1,jp),bfzp(1,jp),tb1,tb2,tb3,h2)

!     Add contribution of lag terms to vzp
      if (if3d_3ds) then
        do ilag=2,nbd
           if (ifgeom) then
              call col3(ta3,vzlagp(1,ilag-1,jp),bm1lag(1,1,1,1,ilag-1),
     $                                ntot1) 
              call cmult(ta3,bd(ilag+1),ntot1)           
           else
              call col3(ta3,vzlagp(1,ilag-1,jp),bm1,ntot1) 
              call cmult(ta3,bd(ilag+1),ntot1)           
           endif
           call add2(tb3,ta3,ntot1)
        enddo
        call add2col2(bfzp(1,jp),tb3,h2,ntot1)
      endif       ! if3d_3ds


      return
      end subroutine makebdfp_3ds
!-----------------------------------------------------------------------

      subroutine lagfieldp_3ds

!     Keep old velocity field(s)

      implicit none 

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'TSTEP'

      include '3DS'

      integer ilag,ntot1

      ntot1 = lx1*ly1*lz1*nelv

      do ilag=nbdinp-1,2,-1
         call opcopy
     $     (vxlagp(1,ilag  ,jp),vylagp(1,ilag  ,jp),vzlagp(1,ilag  ,jp)
     $     ,vxlagp(1,ilag-1,jp),vylagp(1,ilag-1,jp),vzlagp(1,ilag-1,jp))
      enddo
      call opcopy(vxlagp(1,1,jp),vylagp(1,1,jp),vzlagp(1,1,jp)
     $           ,vxp   (1,jp)  ,vyp   (1,jp)  ,vzp   (1,jp) )

!     if3d_3ds
      if (if3d_3ds) then
        do ilag=nbdinp-1,2,-1
          call copy(vzlagp(1,ilag,jp),vzlagp(1,ilag-1,jp),ntot1)
        enddo
        call copy(vzlagp(1,1,jp),vzp(1,jp),ntot1)
      endif


      return
      end subroutine lagfieldp_3ds
!----------------------------------------------------------------------
      subroutine ophx_3ds (out1,out2,out3,inp1,inp2,inp3,h1,h2)

!     OUT = (H1*A+H2*B) * INP  

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'

      include '3DS'


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
      end subroutine ophx_3ds
!-----------------------------------------------------------------------
      subroutine cresvipp_3ds (resv1,resv2,resv3,h1,h2)

!     Compute startresidual/right-hand-side in the velocity solver

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'    ! v?mask
      include 'MASS'    ! bm1

      include '3DS'

      real           resv1 (lx1,ly1,lz1,lelv)
      real           resv2 (lx1,ly1,lz1,lelv)
      real           resv3 (lx1,ly1,lz1,lelv)
      real           h1    (lx1,ly1,lz1,lelv)
      real           h2    (lx1,ly1,lz1,lelv)

      real w1,w2,w3
      common /scruz/ w1    (lx1,ly1,lz1,lelv)
     $ ,             w2    (lx1,ly1,lz1,lelv)
     $ ,             w3    (lx1,ly1,lz1,lelv)

      integer igeom
      common /cgeom/ igeom

      integer ntot1,ntot2
      real const
      

      ntot1 = lx1*ly1*lz1*nelv
      ntot2 = lx2*ly2*lz2*nelv

      if3d = .true.
      call bcdirvc (vxp(1,jp),vyp(1,jp),vzp(1,jp),
     $              v1mask,v2mask,v3mask)
      if3d = .false.

!     prabal. We don't care about traction conditions for now.
!     Maybe need to look at it if added stiffness terms are needed
!     Or if surface tension is needed
!      call bcneutr

      if (mod(jp,2).eq.1) then

!       Need to Extrapolate both pressures at the same time    
        call extrapprp (prextr_3ds(1,1))
        call opgradt (resv1,resv2,resv3,prextr_3ds(1,1))


        jp = jp + 1
        call extrapprp (prextr_3ds(1,2))
        jp = jp - 1

        call map21_all_3ds(resv3,prextr_3ds(1,2))
        const = k_3dsp
        call cmult(resv3,const,ntot1)
        call col2(resv3,bm1,ntot1)
      else
        call extrapprp (prextr_3ds(1,2))
        call opgradt (resv1,resv2,resv3,prextr_3ds(1,2))

        jp = jp - 1
        call extrapprp (prextr_3ds(1,1))
        jp = jp + 1

        call map21_all_3ds(resv3,prextr_3ds(1,1))
        const = -k_3dsp
        call cmult(resv3,const,ntot1)
        call col2(resv3,bm1,ntot1)
      endif    

      call add2(resv1,bfxp(1,jp),ntot1)
      call add2(resv2,bfyp(1,jp),ntot1)
      call add2(resv3,bfzp(1,jp),ntot1)

!     prabal
      call ophx_3ds(w1,w2,w3,vxp(1,jp),vyp(1,jp),
     $              vzp(1,jp),h1,h2)
      call sub2(resv1,w1,ntot1)
      call sub2(resv2,w2,ntot1)
      call sub2(resv3,w3,ntot1)

      return
      end subroutine cresvipp_3ds
!-----------------------------------------------------------------------

      subroutine sethlm_3dsp(h1,h2,intype)

      implicit none

c     Set the variable property arrays H1 and H2
c     in the Helmholtz equation.
c     (associated with variable IFIELD)
c     INTYPE =      integration type

      include 'SIZE'
      include 'INPUT'
      include 'SOLN' 
      include 'TSTEP'   ! ifield

      include '3DS'

      real h1(1),h2(1)

      integer intype
      integer nel,ntot1

      real dtbd

      real k2

      nel   = nelfld(ifield)
      ntot1 = lx1*ly1*lz1*nel

      k2    = k_3dsp**2

      if (iftran) then
         dtbd = bd(1)/dt
         call copy  (h1,vdiff (1,1,1,1,ifield),ntot1)
         if (intype.eq.0) then
            call rzero (h2,ntot1)
         else
            call cmult2(h2,vtrans(1,1,1,1,ifield),dtbd,ntot1)

!           Add second derivative of the 3rd direction to the operator
            call add2s2(h2,vdiff(1,1,1,1,ifield),k2,ntot1)
         endif
      else
         call copy  (h1,vdiff (1,1,1,1,ifield),ntot1)
         call rzero (h2,ntot1)

!        Add second derivative of the 3rd direction to the operator
         call add2s2(h2,vdiff(1,1,1,1,ifield),k2,ntot1)
      endif


      return
      end subroutine sethlm_3dsp

!----------------------------------------------------------------------

      subroutine incomprp_3ds (igeom)
c
c     Project U onto the closest incompressible field
c
c     Output: updated values of U, iproj, proj; and
c             up    := pressure correction req'd to impose div U = 0
c
c
c     Dependencies: ifield ==> which "density" (vtrans) is used.
c

      implicit none

      include 'SIZE'
      include 'SOLN'          ! vxp,vyp,vzp,prp,jp
      include 'INPUT'         ! npert
      include 'TSTEP'         ! dt,ifield
      include 'CTIMER'

      include '3DS'

      real w1,w2,w3
      real dv1,dv2,dv3,dp
      common /scrns/ w1    (lx1,ly1,lz1,lelv)
     $ ,             w2    (lx1,ly1,lz1,lelv)
     $ ,             w3    (lx1,ly1,lz1,lelv)
     $ ,             dv1   (lx1,ly1,lz1,lelv)
     $ ,             dv2   (lx1,ly1,lz1,lelv)
     $ ,             dv3   (lx1,ly1,lz1,lelv)
     $ ,             dp    (lx2,ly2,lz2,lelv)

      real h1,h2,h2inv
      common /scrvh/ h1    (lx1,ly1,lz1,lelv)
     $ ,             h2    (lx1,ly1,lz1,lelv)
      common /scrhi/ h2inv (lx1,ly1,lz1,lelv)

      real dp2
      common /scrch/ dp2(lx2,ly2,lz2,lelv)
      logical ifprjp

      integer ntot1,ntot2,intype,istart

      real dtbd,const

      integer jpi

      integer igeom

      if (igeom.eq.1) return

      if (icalld.eq.0) tpres=0.0
      icalld=icalld+1
      npres=icalld

      jpi = 1
      if (mod(jp,2).eq.0) jpi = 2

      ntot1  = lx1*ly1*lz1*nelv
      ntot2  = lx2*ly2*lz2*nelv
      dtbd   = bd(1)/dt

      call rzero   (h1,ntot1)
      call copy    (h2,vtrans(1,1,1,1,ifield),ntot1)
      call cmult2  (h2,vtrans(1,1,1,1,ifield),dtbd,ntot1)
      call invers2 (h2inv,h2,ntot1)

      call opdiv   (prcorr_3ds(1,jpi),vxp(1,jp),vyp(1,jp),vzp(1,jp))

      if (jpi.eq.1) then
        call map12_all_3ds(dp2,vzp(1,jpi+1))
        const = -k_3dsp
      else
        call map12_all_3ds(dp2,vzp(1,jpi-1))
        const =  k_3dsp
      endif

      call cmult(dp2,const,ntot2)

!      call add2s2(prcorr_3ds(1,jpi),dp2,const,ntot2)

      call chsign  (prcorr_3ds(1,jpi),ntot2)
      call ortho   (prcorr_3ds(1,jpi))


      ifprjp=.false.    ! project out previous pressure solutions?
      istart=param(95)  
      if (istep.ge.istart.and.istart.ne.0) ifprjp=.true.

      ! Most likely, the following can be commented out. (pff, 1/6/2010)
      if (npert.gt.1.or.ifbase)            ifprjp=.false.

      intype =  2             ! Changing integration type here.
                              ! Need to modify cdabdtp accordingly
                              ! Also need to modify uzprec

      call esolver (prcorr_3ds(1,jpi),h1,h2,h2inv,intype)


      return
      end subroutine incomprp_3ds
!------------------------------------------------------------------------
      subroutine velp_update_3ds (igeom)

!     Update Pressure and velocities based on pressure correction


      implicit none

      include 'SIZE'
      include 'SOLN'          ! vxp,vyp,vzp,prp,jp
      include 'INPUT'         ! npert
      include 'TSTEP'         ! dt,ifield
      include 'MASS'
      include 'CTIMER'

      include '3DS'

      include 'TEST'

      real w1,w2,w3
      real dv1,dv2,dv3,dp
      common /scrns/ w1    (lx1,ly1,lz1,lelv)
     $ ,             w2    (lx1,ly1,lz1,lelv)
     $ ,             w3    (lx1,ly1,lz1,lelv)
     $ ,             dv1   (lx1,ly1,lz1,lelv)
     $ ,             dv2   (lx1,ly1,lz1,lelv)
     $ ,             dv3   (lx1,ly1,lz1,lelv)
     $ ,             dp    (lx2,ly2,lz2,lelv)

      real h1,h2,h2inv
      common /scrvh/ h1    (lx1,ly1,lz1,lelv)
     $ ,             h2    (lx1,ly1,lz1,lelv)
      common /scrhi/ h2inv (lx1,ly1,lz1,lelv)

      real dp2
      common /scrch/ dp2(lx2,ly2,lz2,lelv)
      logical ifprjp

      integer ntot1,ntot2

      real dtbd,const

      integer jpi

      integer igeom

      if (igeom.eq.1) return

      jpi = 1
      if (mod(jp,2).eq.0) jpi = 2

      ntot1  = lx1*ly1*lz1*nelv
      ntot2  = lx2*ly2*lz2*nelv
      dtbd   = bd(1)/dt

!     In principle these should still be preserved from the last routine.
!     But I calculate it again for clarity
      call rzero   (h1,ntot1)
      call copy    (h2,vtrans(1,1,1,1,ifield),ntot1)
      call cmult2  (h2,vtrans(1,1,1,1,ifield),dtbd,ntot1)
      call invers2 (h2inv,h2,ntot1)

      call opgradt (w1 ,w2 ,w3 ,prcorr_3ds(1,jpi))

      ! Map pressure to velcity mesh
      if (jpi.eq.1) then

!        call ortho(prcorr_3ds(1,jpi+1))
        call map21_all_3ds(w3,prcorr_3ds(1,jpi+1)) 
        const = -k_3dsp

      else
!        call ortho(prcorr_3ds(1,jpi-1))
        call map21_all_3ds(w3,prcorr_3ds(1,jpi-1)) 
        const = k_3dsp
      endif

      call cmult(w3,const,ntot1)
      call col2(w3,bm1,ntot1)

      if3d = .true.
      call opbinv_3ds(dv1,dv2,dv3,w1 ,w2 ,w3 ,h2inv)
      if3d = .false.

      call add2(vxp(1,jp),dv1,ntot1)
      call add2(vyp(1,jp),dv2,ntot1)
      if (if3d.or.if3d_3ds) call add2(vzp(1,jp),dv3,ntot1)

      call extrapprp(prextr_3ds(1,jpi))
      call lagpresp
      call add3(prp(1,jp),prextr_3ds(1,jpi),prcorr_3ds(1,jpi),ntot2)

      return
      end subroutine velp_update_3ds
!------------------------------------------------------------------------
      subroutine map12_all_3ds(pm2,pm1)

      implicit none

      include 'SIZE'

      real pm1(lx1,ly1,lz1,lelv)
      real pm2(lx2,ly2,lz2,lelv)
      integer e

      do e=1,nelv
         call map12 (pm2(1,1,1,e),pm1(1,1,1,e),e)
      enddo
   
      return
      end subroutine map12_all_3ds
!-----------------------------------------------------------------------

      subroutine map21_all_3ds (y,x)

!     Map X from mesh M2 to mesh M1 (Y)
      
      implicit none

      INCLUDE 'SIZE'

      real x(lx2,ly2,lz2,lelv)
      real y(lx1,ly1,lz1,lelv)

      integer e

      do e=1,nelv
        call map21t(y(1,1,1,e),x(1,1,1,e),e)
      enddo

      return
      end subroutine map21_all_3ds

!-----------------------------------------------------------------------

      subroutine opbinv_3ds (out1,out2,out3,inp1,inp2,inp3,h2inv)

!     Compute OUT = (H2*B)-1 * INP   (explicit)


      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'SOLN'
C
      real out1  (1)
      real out2  (1)
      real out3  (1)
      real inp1  (1)
      real inp2  (1)
      real inp3  (1)
      real h2inv (1)

      real tmp
      integer i,isbcnt,ntot

      include 'OPCTR'
C
#ifdef TIMER
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'opbinv'
      endif
#endif


      ntot=lx1*ly1*lz1*nelv

!      call opmask  (inp1,inp2,inp3)
!      call opdssum (inp1,inp2,inp3)
      call col2 (inp1,v1mask,ntot)
      call col2 (inp2,v2mask,ntot)
      call col2 (inp3,v3mask,ntot)

      call dssum(inp1,lx1,ly1,lz1)
      call dssum(inp2,lx1,ly1,lz1)
      call dssum(inp3,lx1,ly1,lz1)



#ifdef TIMER
      isbcnt = ntot*(1+ldim)
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
#endif

      call invcol3 (out1,bm1,h2inv,ntot)  ! this is expensive and should
      call dssum   (out1,lx1,ly1,lz1)     ! be changed (pff, 3/18/09)
      if (if3d) then
         do i=1,ntot
            tmp = 1./out1(i)
            out1(i)=inp1(i)*tmp
            out2(i)=inp2(i)*tmp
            out3(i)=inp3(i)*tmp
         enddo
      else
         do i=1,ntot
            tmp = 1./out1(i)
            out1(i)=inp1(i)*tmp
            out2(i)=inp2(i)*tmp
         enddo
      endif

      return
      end subroutine opbinv_3ds
c-----------------------------------------------------------------------

      subroutine convect_w_3ds(cku,u,Cz)

!     Compute dealiased form:  J^T Bf *JCz .Ju w/ correct Jacobians

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
!      include 'TOTAL'

      real cku(1),u(1),cx(1),cy(1),cz(1)
      logical ifuf,ifcf            ! u and/or c already on fine mesh?

      integer lxy,ltd
      parameter (lxy=lx1*ly1*lz1,ltd=lxd*lyd*lzd)

      real fx,fy,fz
      real ur,us,ut
      real tr,uf,wd2,jacm1d
      common /scrcv/ fx(ltd),fy(ltd),fz(ltd)
     $             , ur(ltd),us(ltd),ut(ltd)
     $             , jacm1d(ltd),tr(ltd)
     $             , wd2(ltd),uf(ltd)

      integer e,k
      integer iu,ic,ijc,ick

      integer nxyz1,nxyzc,nxyzd,nxyzu,nxyzj

      real zd,wd
      common /dealias1/ zd(lxd),wd(lxd)

      integer i,j,l

      nxyz1 = lx1*ly1*lz1
      nxyzd = lxd*lyd*lzd

      nxyzu = nxyz1
!      if (ifuf) nxyzu = nxyzd

      nxyzc = nxyz1
!      if (ifcf) nxyzc = nxyzd

      nxyzj = nxyz1


      iu  = 1    ! pointer to scalar field u
      ic  = 1    ! pointer to vector field C
      ijc = 1    ! pointer to scalar JACM1
      ick = 1    ! pointer to scalar cku 

      call zwgl (zd,wd,lxd)  ! zwgl -- NOT zwgll!

      if (if3d) then
        do k=1,lzd
        do j=1,lyd
        do i=1,lxd
           l = (k-1)*lyd*lxd + (j-1)*lxd + i 
           wd2(l) = wd(i)*wd(j)*wd(k)
        enddo
        enddo
        enddo
      else
        do j=1,lyd
        do i=1,lxd
           l = (j-1)*lxd + i 
           wd2(l) = wd(i)*wd(j)
        enddo
        enddo

      endif


      do e=1,nelv

!       Interpolate Convecting Field   
        call intp_rstd(fz,cz(ic),lx1,lxd,if3d,0) ! 0 --> forward

!       Interpolate Convected Field   
        call intp_rstd(uf,u(iu),lx1,lxd,if3d,0) ! 0 --> forward

!       Interpolate Jacobian (Probably only needs to be done once) 
        call intp_rstd(jacm1d,jacm1(ijc,1,1,1),lx1,lxd,if3d,0) ! 0 --> forward

        do i=1,nxyzd ! mass matrix included, per DFM (4.8.5)
           tr(i) = wd2(i)*jacm1d(i)*uf(i)*fz(i)
        enddo
        call intp_rstd(cku(ick),tr,lx1,lxd,if3d,1) ! Project back to coarse

        ic  = ic  + nxyzc
        iu  = iu  + nxyzu
        ijc = ijc + nxyzj
        ick = ick + nxyz1

      enddo

      return
      end subroutine convect_w_3ds
!-----------------------------------------------------------------------
      subroutine cdabdtp_3ds(ap,wp,h1,h2,h2inv,intype)

!     INTYPE= 0  Compute the matrix-vector product    DA(-1)DT*p
!     INTYPE= 1  Compute the matrix-vector product    D(B/DT)(-1)DT*p
!     INTYPE=-1  Compute the matrix-vector product    D(A+B/DT)(-1)DT*p
!     INTYPE= 2  Compute the matrix-vector product    D(B/DT)(-1)DT*p
!                  with fourier in 3rd component 

      implicit none

      include 'SIZE'
      include 'INPUT'         ! if3d
      include 'MASS'
      include '3DS'

      include 'TEST'

!      include 'TOTAL'
      real           ap    (lx2,ly2,lz2,lelv)
      real           wp    (lx2,ly2,lz2,lelv)
      real           h1    (lx1,ly1,lz1,lelv)
      real           h2    (lx1,ly1,lz1,lelv)
      real           h2inv (lx1,ly1,lz1,lelv)

      real ta1,ta2,ta3,tb1,tb2,tb3
      common /scrns/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)
     $ ,             tb1 (lx1,ly1,lz1,lelv)
     $ ,             tb2 (lx1,ly1,lz1,lelv)
     $ ,             tb3 (lx1,ly1,lz1,lelv)

      real ttmp2(lx2,ly2,lz2,lelv)         ! lazy work. Should use a scratch array

      integer ntot1,ntot2,intype

      real const

      ntot1 = nx1*ny1*nz1*nelv
      ntot2 = nx2*ny2*nz2*nelv

      const = k_3dsp**2

      call opgradt (ta1,ta2,ta3,wp)
      call map21_all_3ds(ta3,wp)
      call col2(ta3,bm1,ntot1)            ! opgradt includes a mass matrix
      call cmult  (ta3,const,ntot1)

      if3d = .true.
      call opbinv_3ds (tb1,tb2,tb3,ta1,ta2,ta3,h2inv)
      if3d = .false.

!      call copy(tmp1,ta3,ntot1)
!      call copy(tmp2,tb3,ntot1)

      call map12_all_3ds(ttmp2,tb3)
!      call cmult  (tmp2,const,ntot2)

      call opdiv  (ap,tb1,tb2,tb3)
      call add2   (ap,tmp2,ntot2)


!      call copy(tmp7,ttmp2,ntot2)
!      write(6,*) 'nx1', nx1,ny1,nz1,nelv,ntot1,const
!      write(6,*) 'nx2', nx2,ny2,nz2,nelv,ntot2,const

      return
      end subroutine cdabdtp_3ds
!-----------------------------------------------------------------------

      subroutine copy_fou(dataout,datain,nz)

      implicit none

      include 'SIZE'

      real datain(1),dataout(1)
      integer i,skip,nz
      
      skip = lx1*ly1*1*lelv

      do i=1,nz
        j = (i-1)*skip + 1    
        dataout(i) = datain(j)
      enddo

      return
      end subroutine copy_fou
!-----------------------------------------------------------------------

      subroutine phy_to_fou(fldout_r,flout_i,fldin,wkr,wkc,plan,nz)

      implicit none
      
      include 'SIZE'

      integer e,j,i,k,nz

      real fldin(lx1,ly1,lelv,nz)
      real fldout_r(lx1,ly1,lelv,nz)
      real fldout_i(lx1,ly1,lelv,nz)
      real wkr(nz)
      complex wkc(nz)
      integer plan

!     prabal. Perform a fourier transform
      do e=1,nel
        do j=1,ly1
          do i=1,lx1
            call copy_fou(wkr,fldin(i,j,e,1),nz)    ! copy data to a continuous array
            call dfftw_execute_dft_r2c(plan,wkr,wkc)

            do k=1,nz
              fldout_r(i,j,e,k)=real(wkc(k))   
              fldout_i(i,j,e,k)=aimag(wkc(k))
            enddo  

          enddo
        enddo
      enddo

      return
      end subroutine phy_to_fou_all

!-----------------------------------------------------------------------



