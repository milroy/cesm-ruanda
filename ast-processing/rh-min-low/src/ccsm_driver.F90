program ccsm_driver

!-------------------------------------------------------------------------------
!
! Purpose: Main program for NCAR CCSM4/cpl7. Can have different
!          land, sea-ice, and ocean models plugged in at compile-time.
!          These models can be either: stub, dead, data, or active
!          components or some combination of the above.
!
!               stub -------- Do nothing.
!               dead -------- Send analytic data back.
!               data -------- Send data back interpolated from input files.
!               active ------ Prognostically simulate the given component.
!
! Method: Call appropriate initialization, run (time-stepping), and 
!         finalization routines.
! 
!-------------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! share code & libs
   !----------------------------------------------------------------------------
   use shr_sys_mod,       only: shr_sys_abort
   use perf_mod
   use ESMF
   use ccsm_comp_mod

   implicit none

   !--------------------------------------------------------------------------
   ! Local Variables
   !--------------------------------------------------------------------------
   integer                    :: localrc






   !--------------------------------------------------------------------------
   ! Setup and initialize the communications and logging.  
   !--------------------------------------------------------------------------
   call ccsm_pre_init1()

   !--------------------------------------------------------------------------
   ! Initialize ESMF.  This is done outside of the ESMF_INTERFACE ifdef
   ! because it is needed for the time manager, even if the ESMF_INTERFACE
   ! is not used.
   !--------------------------------------------------------------------------
   call ESMF_Initialize()

   !--------------------------------------------------------------------------
   ! Read in the configuration information and initialize the time manager.
   !--------------------------------------------------------------------------
   call ccsm_pre_init2()


   !--------------------------------------------------------------------------
   ! If ESMF is not defined, then just call the initialize, run and finalize
   ! routines directly.
   !--------------------------------------------------------------------------
   call ccsm_init()
   call ccsm_run()
   call ccsm_final()


   !--------------------------------------------------------------------------
   ! Clean-up
   !--------------------------------------------------------------------------
   call ESMF_Finalize( )


end program ccsm_driver
