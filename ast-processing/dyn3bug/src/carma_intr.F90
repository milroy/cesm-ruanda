!! This module is stub for a coupler between the 1 model and the Community Aerosol
!! and Radiation Model for Atmospheres (CARMA) microphysics model. It is used when
!! CARMA is not being used, so that the 1 code that calls CARMA does not need to
!! be changed. The real version of this routine exists in the directory
!! physics/carma/cam. A CARMA model can be activated by using configure with the
!! option:
!!
!!  -carma <carma_pkg>
!!
!! where carma_pkg is the name for a particular microphysical model.
!!
!! @author  Chuck Bardeen
!! @version May 2009
module carma_intr

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use pmgrid,         only: plat, plev, plevp, plon
  use ppgrid,         only: pcols, pver, pverp
  use constituents,   only: pcnst
  use physics_types,  only: physics_state, physics_ptend, physics_ptend_init
  use physics_buffer, only: physics_buffer_desc


  implicit none
  
  private
  save

  ! Public interfaces
  
  ! 1 Physics Interface
  public carma_register                 ! register consituents
  public carma_is_active                ! retrns true if this package is active (microphysics = .true.)
  public carma_implements_cnst          ! returns true if consituent is implemented by this package
  public carma_init_cnst                ! initialize constituent mixing ratios, if not read from initial file
  public carma_init                     ! initialize timestep independent variables
  public carma_final                    ! finalize the CARMA module
  public carma_timestep_init            ! initialize timestep dependent variables
  public carma_timestep_tend            ! interface to tendency computation
  public carma_accumulate_stats         ! collect stats from all MPI tasks
  
  ! Other Microphysics
  public carma_emission_tend            ! calculate tendency from emission source function
  public carma_wetdep_tend              ! calculate tendency from wet deposition
  
contains


  subroutine carma_register
    implicit none

    return
  end subroutine carma_register


! Subprogram not used   function carma_is_active()
! Subprogram not used     implicit none
! Subprogram not used   
! Subprogram not used     logical :: carma_is_active
! Subprogram not used   
! Subprogram not used     carma_is_active = .false.
! Subprogram not used     
! Subprogram not used     return
! Subprogram not used   end function carma_is_active


! Subprogram not used   function carma_implements_cnst(name)
! Subprogram not used     implicit none
! Subprogram not used     
! Subprogram not used     character(len=*), intent(in) :: name   !! constituent name
! Subprogram not used     logical :: carma_implements_cnst       ! return value
! Subprogram not used 
! Subprogram not used     carma_implements_cnst = .false.
! Subprogram not used     
! Subprogram not used     return
! Subprogram not used   end function carma_implements_cnst
  

  subroutine carma_init
    implicit none
    
    return
  end subroutine carma_init


  subroutine carma_final
    implicit none
        
    return
  end subroutine carma_final
  

  subroutine carma_timestep_init
    implicit none

    return
  end subroutine carma_timestep_init


  subroutine carma_timestep_tend(state, cam_in, cam_out, ptend, dt, pbuf, dlf, rliq, prec_str, snow_str, &
    prec_sed, snow_sed, ustar, obklen)
    use hycoef,           only: hyai, hybi, hyam, hybm
    use time_manager,     only: get_nstep, get_step_size, is_first_step
    use camsrfexch,       only: cam_in_t, cam_out_t
    use scamMod,          only: single_column
 
    implicit none

    type(physics_state), intent(inout) :: state                 !! physics state variables
    type(cam_in_t), intent(in)         :: cam_in                !! surface inputs
    type(cam_out_t), intent(inout)     :: cam_out               !! cam output to surface models
    type(physics_ptend), intent(out)   :: ptend                 !! constituent tendencies
    real(r8), intent(in)               :: dt                    !! time step (s)
    type(physics_buffer_desc), pointer :: pbuf(:)               !! physics buffer
    real(r8), intent(in), optional     :: dlf(pcols,pver)       !! Detraining cld H20 from convection (kg/kg/s)
    real(r8), intent(inout), optional  :: rliq(pcols)           !! vertical integral of liquid not yet in q(ixcldliq)
    real(r8), intent(out), optional    :: prec_str(pcols)       !! [Total] sfc flux of precip from stratiform (m/s) 
    real(r8), intent(out), optional    :: snow_str(pcols)       !! [Total] sfc flux of snow from stratiform   (m/s)
    real(r8), intent(out), optional    :: prec_sed(pcols)       !! total precip from cloud sedimentation (m/s)
    real(r8), intent(out), optional    :: snow_sed(pcols)       !! snow from cloud ice sedimentation (m/s)
    real(r8), intent(in), optional     :: ustar(pcols)          !! friction velocity (m/s)
    real(r8), intent(in), optional     :: obklen(pcols)         !! Obukhov length [ m ]
    
    call physics_ptend_init(ptend,state%psetcols,'none') !Initialize an empty ptend for use with physics_update

    if (present(prec_str))  prec_str(:)    = 0._r8
    if (present(snow_str))  snow_str(:)    = 0._r8
    if (present(prec_sed))  prec_sed(:)    = 0._r8
    if (present(snow_sed))  snow_sed(:)    = 0._r8

    return
  end subroutine carma_timestep_tend


! Subprogram not used   subroutine carma_init_cnst(name, q, gcid)
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     character(len=*), intent(in) :: name               !! constituent name
! Subprogram not used     real(r8), intent(out)        :: q(plon,plev,plat)  !! mass mixing ratio
! Subprogram not used     integer, intent(in)          :: gcid(:)            !! global column id
! Subprogram not used     
! Subprogram not used     if (name == "carma") then
! Subprogram not used       q = 0._r8
! Subprogram not used     end if 
! Subprogram not used     
! Subprogram not used     return
! Subprogram not used   end subroutine carma_init_cnst


! Subprogram not used   subroutine carma_emission_tend(state, ptend, cam_in, dt)
! Subprogram not used     use camsrfexch,       only: cam_in_t
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used     
! Subprogram not used     type(physics_state), intent(in )    :: state                !! physics state
! Subprogram not used     type(physics_ptend), intent(inout)  :: ptend                !! physics state tendencies
! Subprogram not used     type(cam_in_t),      intent(inout)  :: cam_in               !! surface inputs
! Subprogram not used     real(r8),            intent(in)     :: dt                   !! time step (s)
! Subprogram not used 
! Subprogram not used     return
! Subprogram not used   end subroutine carma_emission_tend 


! Subprogram not used   subroutine carma_wetdep_tend(state, ptend, dt,  pbuf, dlf, cam_out)
! Subprogram not used     use camsrfexch,       only: cam_out_t
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     real(r8),             intent(in)    :: dt             !! time step (s)
! Subprogram not used     type(physics_state),  intent(in )   :: state          !! physics state
! Subprogram not used     type(physics_ptend),  intent(inout) :: ptend          !! physics state tendencies
! Subprogram not used     type(physics_buffer_desc), pointer  :: pbuf(:)        !! physics buffer
! Subprogram not used     real(r8), intent(in)                :: dlf(pcols,pver)       !! Detraining cld H20 from convection (kg/kg/s)
! Subprogram not used     type(cam_out_t),      intent(inout) :: cam_out        !! cam output to surface models
! Subprogram not used 
! Subprogram not used     return
! Subprogram not used   end subroutine carma_wetdep_tend


  subroutine carma_accumulate_stats()
    implicit none

  end subroutine carma_accumulate_stats
end module carma_intr
