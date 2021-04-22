module seq_cdata_mod

  use shr_kind_mod     , only: r8=> shr_kind_r8
  use shr_sys_mod      , only: shr_sys_flush
  use shr_sys_mod      , only: shr_sys_abort
  use seq_infodata_mod , only: seq_infodata_type
  use mct_mod
  use seq_comm_mct

  implicit none
  save
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------
  public :: seq_cdata_setptrs
  public :: seq_cdata_init    ! only used by xxx_comp_esm.F90 for data models

  !--------------------------------------------------------------------------
  ! Public data
  !--------------------------------------------------------------------------
  ! in general, this type just groups together related data via pointers

  type seq_cdata
     character(len=16)                :: name               ! user defined name
     integer                          :: ID                 ! component id
     integer                          :: mpicom             ! mpi communicator
     type(mct_gGrid)         ,pointer :: dom => null()      ! domain info
     type(mct_gsMap)         ,pointer :: gsMap => null()    ! decomp info
     type(seq_infodata_type) ,pointer :: infodata => null() ! Input init object
  end type seq_cdata
  
  public seq_cdata

!==============================================================================
contains
!==============================================================================

  subroutine seq_cdata_setptrs(cdata, ID, mpicom, dom, gsMap, infodata, name)

    !-----------------------------------------------------------------------
    !
    ! Arguments
    type(seq_cdata)         ,intent(in)       :: cdata      ! input
    integer                 ,optional         :: ID         ! component id
    integer                 ,optional         :: mpicom     ! mpi comm
    type(mct_gGrid)         ,optional,pointer :: dom        ! domain
    type(mct_gsMap)         ,optional,pointer :: gsMap      ! decomp
    type(seq_infodata_type) ,optional,pointer :: infodata   ! INIT object
    character(len=*)        ,optional         :: name       ! name
    !
    ! Local variables
    character(*),parameter :: subName = '(seq_cdata_setptrs) '
    !-----------------------------------------------------------------------

    if (present(name     )) name     =  cdata%name
    if (present(ID       )) ID       =  cdata%ID
    if (present(mpicom   )) mpicom   =  cdata%mpicom
    if (present(dom      )) dom      => cdata%dom
    if (present(gsMap    )) gsMap    => cdata%gsMap
    if (present(infodata )) infodata => cdata%infodata

  end subroutine seq_cdata_setptrs

  !===============================================================================

! Subprogram not used   subroutine seq_cdata_init(cdata,ID,dom,gsMap,infodata,name)
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ! Description
! Subprogram not used     ! This is here only for backwards compatibility with current data model
! Subprogram not used     ! xxx_comp_esmf.F90 interfaces
! Subprogram not used     !
! Subprogram not used     ! Arguments
! Subprogram not used     implicit none
! Subprogram not used     type(seq_cdata)         ,intent(inout)       :: cdata      ! initialized
! Subprogram not used     integer                 ,intent(in)          :: ID         ! component id
! Subprogram not used     type(mct_gGrid)         ,intent(in),target   :: dom        ! domain
! Subprogram not used     type(mct_gsMap)         ,intent(in),target   :: gsMap      ! decomp
! Subprogram not used     type(seq_infodata_type) ,intent(in),target   :: infodata   ! INIT object
! Subprogram not used     character(len=*)        ,intent(in),optional :: name       ! user defined name
! Subprogram not used     !
! Subprogram not used     ! Local variables
! Subprogram not used     !
! Subprogram not used     integer :: mpicom     ! mpi communicator
! Subprogram not used     character(*),parameter :: subName = '(seq_cdata_init) '
! Subprogram not used     logical :: iamroot    ! iamroot
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     call seq_comm_setptrs(ID, mpicom=mpicom, iamroot=iamroot)
! Subprogram not used 
! Subprogram not used     if (present(name)) then
! Subprogram not used       cdata%name   =  name
! Subprogram not used     else
! Subprogram not used       cdata%name   =  'undefined'
! Subprogram not used     endif
! Subprogram not used     cdata%ID       =  ID
! Subprogram not used     cdata%mpicom   =  mpicom
! Subprogram not used     cdata%dom      => dom
! Subprogram not used     cdata%gsMap    => gsMap
! Subprogram not used     cdata%infodata => infodata
! Subprogram not used 
! Subprogram not used   end subroutine seq_cdata_init

end module seq_cdata_mod
