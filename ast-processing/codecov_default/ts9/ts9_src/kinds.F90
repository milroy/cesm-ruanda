module kinds
      use shr_kind_mod, only : SHR_KIND_I4, SHR_KIND_R8, SHR_KIND_I8, SHR_KIND_CL
      use cam_logfile, only : iulog ! _EXTERNAL
implicit none
private
!
!  most floating point variables should be of type real_kind = real*8
!  For higher precision, we also have quad_kind = real*16, but this
!  is only supported on IBM systems
! 
  integer (kind=4), public, parameter::  &
  int_kind     = SHR_KIND_I4,            &
  log_kind     = kind(.true.),           &
  long_kind    = SHR_KIND_I8,            &
  real_kind    = SHR_KIND_R8,            &
  longdouble_kind    = 8
  public :: shr_kind_cl, iulog


end module kinds

