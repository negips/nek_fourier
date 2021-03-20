!======================================================================
! Name        : fft
! Author      : Prabal S. Negi
! Version     : last modification 2019.03.28
! Copyright   : GPL
! Description : Calculate the FFT transform
!               for input field snapshots
!======================================================================      
      subroutine fft_make(maxfld,prefix,bmask,iskip,istart)
!
!     Main driver routine for FFT calculations using windowing
!
!     maxfld   : max number of snapshots to be used
!     prefix   : prefix for the snapshots
!     bmask    : user-defined binary mask used to select regions of interest
!     iskip    : read every iskip:th field
!     istart   : file index to start reading sequence with
!    
!======================================================================
!----------------------------------------------------------------------
     
      implicit none

      include 'fftw3.f'

      include 'SIZE_DEF'
      include 'SIZE'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'TSTEP_DEF'
      include 'TSTEP'


!      include 'DMD_POD'

      integer maxfld,iskip,istart,iend
      character*3 prefix
      real bmask(lx1,ly1,lz1,lelv)

      integer fftl
      parameter (fftl=LKRYL)         ! make sure it is divisible by 2

      integer fftl_in
      parameter (fftl_in=fftl)

      integer fftl_out
      parameter (fftl_out=fftl)     ! ideally should be (fftl_out=fftl/2+1)

      real data_in(fftl_in)
      complex data_out(fftl_out)
      complex data_out2(fftl_out)

      real Amp_avg(LKRYL)

      integer*8 plan,plan2
      integer data_len

      integer i,j,nt,nmodes,istore
!      real pi

      real t0,tn,Tlen,deltaT,Freq,Afreq
      real AmpCorr            ! Amplitude correction when using windowing
      real wwts(LKRYL)        ! windowing weights
      integer wtype           ! windowing type
      character windowing*32  ! name
      

      pi = 4.*atan(1.)

      nmodes=10000
      wtype = 2
      call blank(windowing,32)
      if (wtype.eq.1) windowing='Rectangular window'
      if (wtype.eq.2) windowing='Hanning window'
      if (wtype.eq.3) windowing='Flat-top window'

!     read fields into VKRYL
      istore=0
      do i=istart,istart+(maxfld-1)*iskip,iskip
        istore = istore+1
        call load_kryl_fld(prefix,i,bmask,istore)
        if (i.eq.istart) t0=TIME    ! start time
      enddo
      tn=TIME           ! end time
      Tlen=tn-t0
      data_len = istore
      deltaT = Tlen/data_len
      Freq = 1./deltaT
      Afreq=2.*pi*Freq

      if (nmodes.gt.data_len) nmodes=(data_len/2+1)

!     Write to logfile
      if (nid.eq.0) then
        write(6,*) 'Evaluating FFT for ',istore, 'fields'
        write(6,*) 'Start time=',t0
        write(6,*) 'End time=',tn
        write(6,*) 'Time Period=',Tlen
        write(6,*) 'Sampling Freq=',Freq
        write(6,*) 'Windowing: ',windowing
        write(6,*) 'No. of FFT modes: ',nmodes
      endif        


!     Forward plan      
      call dfftw_plan_dft_r2c_1d(plan,data_len,data_in,data_out,
     $                           FFTW_ESTIMATE)

!     Backward plan      
      call dfftw_plan_dft_c2r_1d(plan2,data_len,data_out,data_in,
     $                           FFTW_ESTIMATE)
 
      call set_window_wts(wwts,AmpCorr,data_len,wtype,nid)

      call rzero(Amp_avg,LKRYL)
      do i=1,LG
        do j=1,data_len
          data_in(j) = wwts(j)*VKRYL(i,j)
        enddo  

!       FFTW computes the un-normalized transform [See doc. pg 4]
!       Therefore need to divide amplitudes by n      
        call dfftw_execute_dft_r2c(plan,data_in,data_out)

        call copy(data_out2,data_out,2*data_len)   ! twice since its complex

        do j=1,nmodes
          call rzero(data_out,2*data_len)
          call rzero(data_in,data_len)

          data_out(j) = data_out2(j)/AmpCorr
          Amp_avg(j)  = Amp_avg(j) + abs(data_out(j))

          call dfftw_execute_dft_c2r(plan2,data_out,data_in)

          D_MODES(i,j) = data_in(1)       ! take the first time instant
                                          ! of all modes
        enddo 

      enddo  ! i=1,LG
      call dfftw_destroy_plan(plan)
      call dfftw_destroy_plan(plan2)

      call cmult(Amp_avg,1./(LG+0.),nmodes)

      if (nid.eq.0) then
        open(unit=10101,file='fft_amp.out',status='unknown',
     $       form='formatted')
        write(10101,'(A4,1x,2(A24,1x))'),'i','Freq','Amp'
      endif  

      nt = nx1*ny1*nz1*nelv
      do i=1,nmodes
        call opcopy(vx,vy,vz,D_MODES(1,i),
     $              D_MODES(1+nt,i),D_MODES(1+2*nt,i))
        istep=i
        time=Freq*(i-1.)/(2.*data_len)
        call outpost(vx,vy,vz,pr,t,'fft')
        if (nid.eq.0) then
          write(10101,'(I4,1x,2(E24.16E3,1x))'),i,time,Amp_avg(i)
          flush(10101)
        endif  
       
      enddo  
      close(10101)

      return
      end subroutine fft_make
!---------------------------------------------------------------------- 

      subroutine set_window_wts(wwts,AmpCorr,n,wtype,nid)

      implicit none

      integer i,j,n1
      integer nid             ! only for write statement

      integer n               ! data size
      real wwts(n)            ! calculated weights
      real AmpCorr            ! Amplitude correction
      integer wtype           ! window type:
                              ! 1: Rectangular
                              ! 2: Hann
                              ! 3: Flat-top

                              
      real a0,a1,a2,a3,a4     ! Flat-top coefficients
      real pi
      real vlsum              ! function

      pi = 4.0*atan(1.0)

      n1 = n-1
      AmpCorr = 0.
      call rzero(wwts,n)      ! initialize weights

      if (wtype.eq.1) then
!       Rectangular            
        do i=1,n
          wwts(i)=1.
        enddo
      elseif (wtype.eq.2) then
!       Hann/Hanning window            
        do i=1,n
          j=i-1
          wwts(i)=0.5*(1. - cos(2.*pi*j/n1))
        enddo
      elseif (wtype.eq.3) then
!       Flat-top window
!       Coefficients used from Matlab website            
        a0=0.21557895
        a1=0.41663158
        a2=0.277263158
        a3=0.083578947
        a4=0.006947368 
        do i=1,n
          j=i-1
          wwts(i)=a0 - a1*cos(2.*pi*j/n1) + a2*cos(4.*pi*j/n1) 
     $               - a3*cos(6.*pi*j/n1) + a4*cos(8.*pi*j/n1) 
        enddo
      else
        if (nid.eq.0) write(6,*) 'Unknown windowing mode'
        call exitt
      endif

      AmpCorr = vlsum(wwts,n) 
      if (nid.eq.0) write(6,*) 'Amp Corr:', AmpCorr

      return
      end subroutine set_window_wts
!---------------------------------------------------------------------- 









