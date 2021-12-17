!======================================================================
!     Description: Routines for free surface simulation
!     Author: Prabal S. Negi
!     Subroutines to make the module compatible with Framework
!      
!====================================================================== 
!-----------------------------------------------------------------------
!> @brief Register FS_ALE module
!! @note This routine should be called in frame_usr_register

      subroutine frame_register_fs()     
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'FRAMELP'
      include 'FS_ALE'

      integer lpmid
      real ltim

      real dnekclock

!     timing
      ltim = dnekclock()

!     Check if the current module was already registered
      call mntr_mod_is_name_reg(lpmid,fs_name)
      if (lpmid.gt.0) then
         call mntr_warn(lpmid,
     $        'module ['//trim(fs_name)//'] already registered')
         return
      endif

!     Find parent module
      call mntr_mod_is_name_reg(lpmid,'FRAME')
      if (lpmid.le.0) then
         lpmid = 1
         call mntr_abort(lpmid,
     $        'Parent module ['//'FRAME'//'] not registered')
      endif

!     Register module
      call mntr_mod_reg(fs_id,lpmid,fs_name,
     $          'Free Surface ALE.')

!     Register timers
      call mntr_tmr_is_name_reg(lpmid,'FRM_TOT')
!     Total time
      call mntr_tmr_reg(fs_tmr_tot_id,lpmid,fs_id,
     $     'FS_TOT','FS module total time',.false.)
!     Initialisation time
      call mntr_tmr_reg(fs_tmr_ini_id,lpmid,fs_id,
     $     'FS_INI','FS initialisation time',.true.)
       
!     Register and set active section
      call rprm_sec_reg(fs_sec_id,fs_id,
     $     '_'//adjustl(fs_name),
     $     'Runtime parameter section for FS module')
      call rprm_sec_set_act(.true.,fs_sec_id)

!     Register parameters
!     IFFS
      call rprm_rp_reg(fs_iffs_id,fs_sec_id,'FS_IFFS',
     $     'Enable FS? ',
     $     rpar_log,0,0.0,.false.,'  ')
     
!     FS_OFST
      call rprm_rp_reg(fs_ofst_id,fs_sec_id,'FS_OFST',
     $     'Damping Offset ',
     $     rpar_real,0,0.0,.false.,'  ')

!!     FS_SLIPL
!      call rprm_rp_reg(fs_slipl_id,fs_sec_id,'FS_SLIPL',
!     $     'Slip Length ',
!     $     rpar_real,0,0.0,.false.,'  ')
!
!!     FS_BLENDL
!      call rprm_rp_reg(fs_blendl_id,fs_sec_id,'FS_BLENDL',
!     $     'Blending Length ',
!     $     rpar_real,0,1.0,.false.,'  ')
     
      ! set initialisation flag
!      otd_ifinit=.false.

      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(fs_tmr_tot_id,1,ltim)

      return
      end subroutine frame_register_fs

!---------------------------------------------------------------------- 

      subroutine frame_get_param_fs()

      implicit none

      include 'SIZE'
      include 'INPUT'               ! param(59), initc
      include 'TSTEP'               ! time
      include 'FRAMELP'
      include 'FS_ALE'

      ! local variables
      integer       itmp
      real          rtmp, ltim
      logical       ltmp
      character*20  ctmp
      character*2   str1, str2
      character*200 lstring

      ! functions
      real dnekclock
!-----------------------------------------------------------------------
      ! timing
      ltim = dnekclock()

      ! get runtime parameters
!     iffs
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,fs_iffs_id,rpar_log)
      fs_iffs = ltmp
!     fs_ofst
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,fs_ofst_id,rpar_real)
      fs_ofst = rtmp
!!     fs_slipl
!      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,fs_slipl_id,rpar_real)
!      fs_slipl = rtmp
!!     fs_blendl
!      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,fs_blendl_id,rpar_real)
!      fs_blendl = rtmp

      return
      end subroutine frame_get_param_fs 
!---------------------------------------------------------------------- 

