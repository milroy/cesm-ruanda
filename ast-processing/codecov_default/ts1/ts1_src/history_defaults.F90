module history_defaults
!----------------------------------------------------------------------- 
! 
! Purpose: contains calls to setup default history stuff that has not found
!          a proper home yet. Shouldn't really exist.
!
! Public functions/subroutines:
!   bldfld
! 
! Author: B.A. Boville from code in cam_history.F90
!-----------------------------------------------------------------------
  use shr_kind_mod, only: r8 => shr_kind_r8, r4 => shr_kind_r4
  use constituents, only: pcnst, cnst_name
  use ppgrid,       only: pver, pverp
  use pmgrid,       only: plev, plevp
  use dycore,       only: dycore_is

  use cam_history,  only: phys_decomp, dyn_decomp, addfld, add_default
  implicit none

  PRIVATE

  public :: bldfld


CONTAINS


!#######################################################################
  subroutine bldfld ()
!
!----------------------------------------------------------------------- 
! 
! Purpose: 
!
! Build Master Field List of all possible fields in a history file.  Each field has 
! associated with it a "long_name" netcdf attribute that describes what the field is, 
! and a "units" attribute.
! 
! Method: Call a subroutine to add each field
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
!
! Local workspace
!
    integer m                     ! Index

!
! Call addfld to add each field to the Master Field List.
!
    call addfld ('SGH     ','m       ',1,    'I','Standard deviation of orography',phys_decomp)
    call addfld ('SGH30   ','m       ',1,    'I','Standard deviation of 30s orography',phys_decomp)


!jt
!jt Maybe add this to scam specific initialization
!jt


    call addfld ('DQP     ','kg/kg/s ',pver, 'A','Specific humidity tendency due to precipitation',phys_decomp)

  end subroutine bldfld

!#######################################################################

end module history_defaults
