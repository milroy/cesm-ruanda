module prep_wav_mod

  use shr_kind_mod    , only: r8 => SHR_KIND_R8 
  use shr_kind_mod    , only: cs => SHR_KIND_CS
  use shr_kind_mod    , only: cl => SHR_KIND_CL
  use shr_sys_mod     , only: shr_sys_abort, shr_sys_flush
  use seq_comm_mct    , only: num_inst_atm, num_inst_ice, num_inst_ocn 
  use seq_comm_mct    , only: num_inst_wav, num_inst_frc
  use seq_comm_mct    , only: CPLID, WAVID, logunit
  use seq_comm_mct    , only: seq_comm_getdata=>seq_comm_setptrs 
  use seq_infodata_mod, only: seq_infodata_getdata, seq_infodata_type   
  use seq_map_type_mod 
  use seq_map_mod
  use seq_flds_mod
  use t_drv_timers_mod
  use mct_mod
  use perf_mod
  use component_type_mod, only: component_get_x2c_cx, component_get_c2x_cx
  use component_type_mod, only: wav, ocn, ice, atm

  implicit none
  save
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: prep_wav_init
  public :: prep_wav_mrg

  public :: prep_wav_calc_a2x_wx
  public :: prep_wav_calc_o2x_wx
  public :: prep_wav_calc_i2x_wx

  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

  private :: prep_wav_merge

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  ! mappers
  type(seq_map), pointer :: mapper_sa2w
  type(seq_map), pointer :: mapper_so2w
  type(seq_map), pointer :: mapper_si2w

  ! attribute vectors 
  type(mct_aVect), pointer :: o2x_wx(:) ! Ocn export, wav grid, cpl pes 
  type(mct_aVect), pointer :: i2x_wx(:) ! Ice export, wav grid, cpl pes 
  type(mct_aVect), pointer :: a2x_wx(:) ! Atm export, wav grid, cpl pes 

  ! accumulation variables
  ! none at this time

  ! seq_comm_getData variables
  integer :: mpicom_CPLID                     ! MPI cpl communicator
  !================================================================================================

contains

  !================================================================================================

  subroutine prep_wav_init(infodata, atm_c2_wav, ocn_c2_wav, ice_c2_wav)
    
    !---------------------------------------------------------------
    ! Description
    ! Initialize module attribute vectors and all other non-mapping
    ! module variables
    !
    ! Arguments
    type(seq_infodata_type) , intent(in)    :: infodata
    logical                 , intent(in)    :: atm_c2_wav ! .true.  => atm to wav coupling on
    logical                 , intent(in)    :: ocn_c2_wav ! .true.  => ocn to wav coupling on
    logical                 , intent(in)    :: ice_c2_wav ! .true.  => ocn to wav coupling on
    !
    ! Local Variables
    integer                     :: eai , eoi, eii, ewi
    integer                     :: lsize_w
    logical                     :: samegrid_ow   ! samegrid ocean and wave
    logical                     :: samegrid_aw   ! samegrid atm and wave
    logical                     :: iamroot_CPLID ! .true. => CPLID masterproc
    logical                     :: esmf_map_flag ! .true. => use esmf for mapping
    logical                     :: wav_present   ! .true. => wav is present
    character(CL)               :: atm_gnam      ! atm grid
    character(CL)               :: ocn_gnam      ! ocn grid
    character(CL)               :: wav_gnam      ! wav grid
    type(mct_avect) , pointer   :: w2x_wx
    character(*)    , parameter :: subname = '(prep_wav_init)'
    character(*)    , parameter :: F00 = "('"//subname//" : ', 4A )"
    !---------------------------------------------------------------

    call seq_infodata_getData(infodata, &
         wav_present=wav_present      , &
         ocn_gnam=ocn_gnam            , &
         wav_gnam=wav_gnam            , &
         esmf_map_flag=esmf_map_flag  )

    allocate(mapper_sa2w)
    allocate(mapper_so2w)
    allocate(mapper_si2w)

    if (wav_present) then

       call seq_comm_getData(CPLID, mpicom=mpicom_CPLID, iamroot=iamroot_CPLID)

       w2x_wx => component_get_c2x_cx(wav(1)) 
       lsize_w = mct_aVect_lsize(w2x_wx)

       allocate(a2x_wx(num_inst_atm))
       do eai = 1,num_inst_atm
          call mct_aVect_init(a2x_wx(eai), rList=seq_flds_a2x_fields, lsize=lsize_w)
          call mct_aVect_zero(a2x_wx(eai))
       enddo
       allocate(o2x_wx(num_inst_ocn))
       do eoi = 1,num_inst_ocn
          call mct_aVect_init(o2x_wx(eoi), rList=seq_flds_o2x_fields, lsize=lsize_w)
          call mct_aVect_zero(o2x_wx(eoi))
       enddo
       allocate(i2x_wx(num_inst_ice))
       do eii = 1,num_inst_ice
          call mct_aVect_init(i2x_wx(eii), rList=seq_flds_i2x_fields, lsize=lsize_w)
          call mct_aVect_zero(i2x_wx(eii))
       enddo
     
       samegrid_ow = .true.
       samegrid_aw = .true.
       if (trim(ocn_gnam) /= trim(wav_gnam)) samegrid_ow = .false.
       if (trim(atm_gnam) /= trim(wav_gnam)) samegrid_aw = .false.

       if (atm_c2_wav) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Sa2w'
          end if
          call seq_map_init_rcfile(mapper_Sa2w, atm(1), wav(1), &
               'seq_maps.rc','atm2wav_smapname:','atm2wav_smaptype:',samegrid_aw, &
               'mapper_Sa2w initialization')
       endif
       call shr_sys_flush(logunit)
       if (ocn_c2_wav) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_So2w'
          end if
          call seq_map_init_rcfile(mapper_So2w, ocn(1), wav(1), &
               'seq_maps.rc','ocn2wav_smapname:','ocn2wav_smaptype:',samegrid_ow, &
               'mapper_So2w initialization')
       endif
       call shr_sys_flush(logunit)  !TODO ??? is this in Tony's code
       if (ice_c2_wav) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Si2w'
          end if
          call seq_map_init_rcfile(mapper_Si2w, ice(1), wav(1), &
               'seq_maps.rc','ice2wav_smapname:','ice2wav_smaptype:',samegrid_ow, &
               'mapper_Si2w initialization')
       endif
       call shr_sys_flush(logunit)

    end if

  end subroutine prep_wav_init

  !================================================================================================

! Subprogram not used   subroutine prep_wav_mrg(infodata, fractions_wx, timer_mrg)
! Subprogram not used 
! Subprogram not used     !---------------------------------------------------------------
! Subprogram not used     ! Description
! Subprogram not used     ! Merge all wav inputs
! Subprogram not used     !
! Subprogram not used     ! Arguments
! Subprogram not used     type(seq_infodata_type) , intent(in)    :: infodata
! Subprogram not used     type(mct_aVect)         , intent(in)    :: fractions_wx(:)
! Subprogram not used     character(len=*)        , intent(in)    :: timer_mrg
! Subprogram not used     !
! Subprogram not used     ! Local Variables
! Subprogram not used     integer                  :: eai, eoi, eii, ewi, efi 
! Subprogram not used     type(mct_avect), pointer :: x2w_wx
! Subprogram not used     character(*), parameter  :: subname = '(prep_wav_mrg)'
! Subprogram not used     !---------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     call t_drvstartf (trim(timer_mrg),barrier=mpicom_CPLID)
! Subprogram not used     do ewi = 1,num_inst_wav
! Subprogram not used        ! Use fortran mod to address ensembles in merge
! Subprogram not used        eai = mod((ewi-1),num_inst_atm) + 1
! Subprogram not used        eoi = mod((ewi-1),num_inst_ocn) + 1
! Subprogram not used        eii = mod((ewi-1),num_inst_ice) + 1
! Subprogram not used        efi = mod((ewi-1),num_inst_frc) + 1
! Subprogram not used 
! Subprogram not used        x2w_wx => component_get_x2c_cx(wav(ewi)) 
! Subprogram not used        call prep_wav_merge(a2x_wx(eai), o2x_wx(eoi), i2x_wx(eii), fractions_wx(efi), x2w_wx)
! Subprogram not used     enddo
! Subprogram not used     call t_drvstopf  (trim(timer_mrg))
! Subprogram not used 
! Subprogram not used   end subroutine prep_wav_mrg

  !================================================================================================

! Subprogram not used   subroutine prep_wav_merge(a2x_w, o2x_w, i2x_w, frac_w, x2w_w)
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     ! Arguments
! Subprogram not used     type(mct_aVect), intent(in)    :: a2x_w  ! input
! Subprogram not used     type(mct_aVect), intent(in)    :: o2x_w  ! input
! Subprogram not used     type(mct_aVect), intent(in)    :: i2x_w  ! input
! Subprogram not used     type(mct_aVect), intent(in)    :: frac_w ! input
! Subprogram not used     type(mct_aVect), intent(inout) :: x2w_w  ! output
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     integer       :: nflds,i,i1,o1
! Subprogram not used     logical       :: iamroot
! Subprogram not used     logical, save :: first_time = .true.
! Subprogram not used     character(CL),allocatable :: mrgstr(:)   ! temporary string
! Subprogram not used     character(CL) :: field   ! string converted to char
! Subprogram not used     type(mct_aVect_sharedindices),save :: a2x_sharedindices
! Subprogram not used     type(mct_aVect_sharedindices),save :: o2x_sharedindices
! Subprogram not used     type(mct_aVect_sharedindices),save :: i2x_sharedindices
! Subprogram not used     character(*), parameter   :: subname = '(prep_wav_merge) '
! Subprogram not used 
! Subprogram not used     !----------------------------------------------------------------------- 
! Subprogram not used 
! Subprogram not used     call seq_comm_getdata(CPLID, iamroot=iamroot)
! Subprogram not used 
! Subprogram not used     if (first_time) then
! Subprogram not used        nflds = mct_aVect_nRattr(x2w_w)
! Subprogram not used 
! Subprogram not used        allocate(mrgstr(nflds))
! Subprogram not used        do i = 1,nflds
! Subprogram not used           field = mct_aVect_getRList2c(i, x2w_w)
! Subprogram not used           mrgstr(i) = subname//'x2w%'//trim(field)//' ='
! Subprogram not used        enddo
! Subprogram not used 
! Subprogram not used        call mct_aVect_setSharedIndices(a2x_w, x2w_w, a2x_SharedIndices)
! Subprogram not used        call mct_aVect_setSharedIndices(o2x_w, x2w_w, a2x_SharedIndices)
! Subprogram not used        call mct_aVect_setSharedIndices(i2x_w, x2w_w, a2x_SharedIndices)
! Subprogram not used 
! Subprogram not used        !--- document copy operations ---
! Subprogram not used        do i=1,a2x_SharedIndices%shared_real%num_indices
! Subprogram not used           i1=a2x_SharedIndices%shared_real%aVindices1(i)
! Subprogram not used           o1=a2x_SharedIndices%shared_real%aVindices2(i)
! Subprogram not used           field = mct_aVect_getRList2c(i1, a2x_w)
! Subprogram not used           mrgstr(o1) = trim(mrgstr(o1))//' = a2x%'//trim(field)
! Subprogram not used        enddo
! Subprogram not used        do i=1,o2x_SharedIndices%shared_real%num_indices
! Subprogram not used           i1=o2x_SharedIndices%shared_real%aVindices1(i)
! Subprogram not used           o1=o2x_SharedIndices%shared_real%aVindices2(i)
! Subprogram not used           field = mct_aVect_getRList2c(i1, o2x_w)
! Subprogram not used           mrgstr(o1) = trim(mrgstr(o1))//' = o2x%'//trim(field)
! Subprogram not used        enddo
! Subprogram not used        do i=1,i2x_SharedIndices%shared_real%num_indices
! Subprogram not used           i1=i2x_SharedIndices%shared_real%aVindices1(i)
! Subprogram not used           o1=i2x_SharedIndices%shared_real%aVindices2(i)
! Subprogram not used           field = mct_aVect_getRList2c(i1, i2x_w)
! Subprogram not used           mrgstr(o1) = trim(mrgstr(o1))//' = i2x%'//trim(field)
! Subprogram not used        enddo
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     ! Create input wave state directly from atm, ocn, ice output state
! Subprogram not used 
! Subprogram not used     call mct_aVect_copy(aVin=a2x_w, aVout=x2w_w, vector=mct_usevector, sharedIndices=a2x_SharedIndices)
! Subprogram not used     call mct_aVect_copy(aVin=o2x_w, aVout=x2w_w, vector=mct_usevector, sharedIndices=o2x_SharedIndices)
! Subprogram not used     call mct_aVect_copy(aVin=i2x_w, aVout=x2w_w, vector=mct_usevector, sharedIndices=i2x_SharedIndices)
! Subprogram not used 
! Subprogram not used     if (first_time) then
! Subprogram not used        if (iamroot) then
! Subprogram not used           write(logunit,'(A)') subname//' Summary:'
! Subprogram not used           do i = 1,nflds
! Subprogram not used              write(logunit,'(A)') trim(mrgstr(i))
! Subprogram not used           enddo
! Subprogram not used        endif
! Subprogram not used        deallocate(mrgstr)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     first_time = .false.
! Subprogram not used 
! Subprogram not used   end subroutine prep_wav_merge

  !================================================================================================

! Subprogram not used   subroutine prep_wav_calc_a2x_wx(timer)
! Subprogram not used     !---------------------------------------------------------------
! Subprogram not used     ! Description
! Subprogram not used     ! Create a2x_wx (note that a2x_wx is a local module variable)
! Subprogram not used     !
! Subprogram not used     ! Arguments
! Subprogram not used     character(len=*), intent(in) :: timer
! Subprogram not used     !
! Subprogram not used     ! Local Variables
! Subprogram not used     integer :: eai
! Subprogram not used     type(mct_aVect), pointer :: a2x_ax
! Subprogram not used     character(*), parameter  :: subname = '(prep_wav_calc_a2x_wx)'
! Subprogram not used     !---------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
! Subprogram not used     do eai = 1,num_inst_atm
! Subprogram not used        a2x_ax => component_get_c2x_cx(atm(eai))
! Subprogram not used        call seq_map_map(mapper_Sa2w, a2x_ax, a2x_wx(eai), norm=.true.)
! Subprogram not used     enddo
! Subprogram not used     call t_drvstopf  (trim(timer))
! Subprogram not used   end subroutine prep_wav_calc_a2x_wx

  !================================================================================================

! Subprogram not used   subroutine prep_wav_calc_o2x_wx(timer)
! Subprogram not used     !---------------------------------------------------------------
! Subprogram not used     ! Description
! Subprogram not used     ! Create o2x_wx (note that o2x_wx is a local module variable)
! Subprogram not used     !
! Subprogram not used     ! Arguments
! Subprogram not used     character(len=*), intent(in) :: timer
! Subprogram not used     !
! Subprogram not used     ! Local Variables
! Subprogram not used     integer :: eoi
! Subprogram not used     type(mct_aVect), pointer :: o2x_ox
! Subprogram not used     character(*), parameter  :: subname = '(prep_wav_calc_o2x_wx)'
! Subprogram not used     !---------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
! Subprogram not used     do eoi = 1,num_inst_ocn
! Subprogram not used        o2x_ox => component_get_c2x_cx(ocn(eoi))
! Subprogram not used        call seq_map_map(mapper_So2w, o2x_ox, o2x_wx(eoi), norm=.true.)
! Subprogram not used     enddo
! Subprogram not used     call t_drvstopf  (trim(timer))
! Subprogram not used   end subroutine prep_wav_calc_o2x_wx

  !================================================================================================

! Subprogram not used   subroutine prep_wav_calc_i2x_wx(timer)
! Subprogram not used     !---------------------------------------------------------------
! Subprogram not used     ! Description
! Subprogram not used     ! Create i2x_wx (note that i2x_wx is a local module variable)
! Subprogram not used     !
! Subprogram not used     ! Arguments
! Subprogram not used     character(len=*), intent(in) :: timer
! Subprogram not used     !
! Subprogram not used     ! Local Variables
! Subprogram not used     integer :: eii
! Subprogram not used     type(mct_aVect), pointer :: i2x_ix
! Subprogram not used     character(*), parameter   :: subname = '(prep_wav_calc_i2x_wx)'
! Subprogram not used     !---------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
! Subprogram not used     do eii = 1,num_inst_ice
! Subprogram not used        i2x_ix => component_get_c2x_cx(ice(eii))
! Subprogram not used        call seq_map_map(mapper_Si2w, i2x_ix, i2x_wx(eii), norm=.true.)
! Subprogram not used     enddo
! Subprogram not used     call t_drvstopf  (trim(timer))
! Subprogram not used   end subroutine prep_wav_calc_i2x_wx

end module prep_wav_mod
