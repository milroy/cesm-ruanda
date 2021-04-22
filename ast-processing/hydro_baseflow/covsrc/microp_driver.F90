module microp_driver

!-------------------------------------------------------------------------------------------------------
!
! Driver for 1 microphysics parameterizations
!
!-------------------------------------------------------------------------------------------------------

use shr_kind_mod,  only: r8 => shr_kind_r8
use ppgrid,        only: pver
use physics_types, only: physics_state, physics_ptend, physics_tend,  &
                         physics_ptend_copy, physics_ptend_sum
use physics_buffer,only: pbuf_get_index, pbuf_get_field, physics_buffer_desc
use phys_control,  only: phys_getopts

use cldwat2m_macro,only: ini_macro
use micro_mg_cam,  only: micro_mg_cam_readnl, micro_mg_cam_register, &
                         micro_mg_cam_implements_cnst, micro_mg_cam_init_cnst, &
                         micro_mg_cam_init, micro_mg_cam_tend
use cam_logfile,   only: iulog
use abortutils,    only: endrun
use perf_mod,      only: t_startf, t_stopf

implicit none
private
save

public :: &
   microp_driver_readnl,          &
   microp_driver_register,        &
   microp_driver_init_cnst,       &
   microp_driver_implements_cnst, &
   microp_driver_init,            &
   microp_driver_tend

character(len=16)  :: microp_scheme   ! Microphysics scheme

!===============================================================================
contains
!===============================================================================

subroutine microp_driver_readnl(nlfile)

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Read in namelist for microphysics scheme
   !-----------------------------------------------------------------------

   call phys_getopts(microp_scheme_out=microp_scheme)

   select case (microp_scheme)
   case ('MG')
      call micro_mg_cam_readnl(nlfile)
   case ('RK')
      ! microp_driver doesn't handle this one
      continue
   case default
      call endrun('microp_driver_readnl:: unrecognized microp_scheme')
   end select

end subroutine microp_driver_readnl

subroutine microp_driver_register

   ! Register microphysics constituents and fields in the physics buffer.
   !-----------------------------------------------------------------------


   select case (microp_scheme)
   case ('MG')
      call micro_mg_cam_register()
   case ('RK')
      ! microp_driver doesn't handle this one
      continue
   case default
      call endrun('microp_driver_register:: unrecognized microp_scheme')
   end select

end subroutine microp_driver_register

!===============================================================================

! Subprogram not used function microp_driver_implements_cnst(name)
! Subprogram not used 
! Subprogram not used    ! Return true if specified constituent is implemented by the
! Subprogram not used    ! microphysics package
! Subprogram not used 
! Subprogram not used    character(len=*), intent(in) :: name        ! constituent name
! Subprogram not used    logical :: microp_driver_implements_cnst    ! return value
! Subprogram not used 
! Subprogram not used    ! Local workspace
! Subprogram not used    integer :: m
! Subprogram not used    !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    microp_driver_implements_cnst = .false.
! Subprogram not used 
! Subprogram not used    select case (microp_scheme)
! Subprogram not used    case ('MG')
! Subprogram not used       microp_driver_implements_cnst = micro_mg_cam_implements_cnst(name)
! Subprogram not used    case ('RK')
! Subprogram not used       ! microp_driver doesn't handle this one
! Subprogram not used       continue
! Subprogram not used    case default
! Subprogram not used       call endrun('microp_driver_implements_cnst:: unrecognized microp_scheme')
! Subprogram not used    end select
! Subprogram not used 
! Subprogram not used end function microp_driver_implements_cnst

!===============================================================================

! Subprogram not used subroutine microp_driver_init_cnst(name, q, gcid)
! Subprogram not used 
! Subprogram not used    ! Initialize the microphysics constituents, if they are
! Subprogram not used    ! not read from the initial file.
! Subprogram not used 
! Subprogram not used    character(len=*), intent(in)  :: name     ! constituent name
! Subprogram not used    real(r8),         intent(out) :: q(:,:)   ! mass mixing ratio (gcol, plev)
! Subprogram not used    integer,          intent(in)  :: gcid(:)  ! global column id
! Subprogram not used    !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    select case (microp_scheme)
! Subprogram not used    case ('MG')
! Subprogram not used       call micro_mg_cam_init_cnst(name, q, gcid)
! Subprogram not used    case ('RK')
! Subprogram not used       ! microp_driver doesn't handle this one
! Subprogram not used       continue
! Subprogram not used    case default
! Subprogram not used       call endrun('microp_driver_init_cnst:: unrecognized microp_scheme')
! Subprogram not used    end select
! Subprogram not used 
! Subprogram not used end subroutine microp_driver_init_cnst

!===============================================================================

subroutine microp_driver_init(pbuf2d)

   type(physics_buffer_desc), pointer :: pbuf2d(:,:)

   ! Initialize the microphysics parameterizations
   !-----------------------------------------------------------------------

   call ini_macro()

   select case (microp_scheme)
   case ('MG')
      call micro_mg_cam_init(pbuf2d)
   case ('RK')
      ! microp_driver doesn't handle this one
      continue
   case default
      call endrun('microp_driver_init:: unrecognized microp_scheme')
   end select


end subroutine microp_driver_init

!===============================================================================

subroutine microp_driver_tend(state, ptend, dtime, pbuf)

   ! Call the microphysics parameterization run methods.

   ! Input arguments

   type(physics_state), intent(in)    :: state       ! State variables
   type(physics_ptend), intent(out)   :: ptend       ! Package tendencies
   type(physics_buffer_desc), pointer :: pbuf(:)

   real(r8), intent(in)  :: dtime                    ! Timestep

   ! Local variables

   integer :: lchnk
   integer :: ncol

   !======================================================================

   lchnk = state%lchnk
   ncol  = state%ncol

   ! Call MG Microphysics

   select case (microp_scheme)
   case ('MG')
      call t_startf('microp_mg_tend')
      call micro_mg_cam_tend(state, ptend, dtime, pbuf)
      call t_stopf('microp_mg_tend')
   case ('RK')
      ! microp_driver doesn't handle this one
      continue
   case default
      call endrun('microp_driver_tend:: unrecognized microp_scheme')
   end select

end subroutine microp_driver_tend

end module microp_driver
