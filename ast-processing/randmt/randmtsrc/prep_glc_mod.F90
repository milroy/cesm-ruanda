module prep_glc_mod

  use shr_kind_mod    , only: r8 => SHR_KIND_R8 
  use shr_kind_mod    , only: cs => SHR_KIND_CS
  use shr_kind_mod    , only: cl => SHR_KIND_CL
  use shr_sys_mod     , only: shr_sys_abort, shr_sys_flush
  use seq_comm_mct    , only: num_inst_glc, num_inst_lnd, num_inst_frc
  use seq_comm_mct    , only: CPLID, GLCID, logunit
  use seq_comm_mct    , only: seq_comm_getData=>seq_comm_setptrs 
  use seq_infodata_mod, only: seq_infodata_type, seq_infodata_getdata  
  use seq_map_type_mod 
  use seq_map_mod
  use seq_flds_mod
  use t_drv_timers_mod
  use mct_mod
  use perf_mod
  use component_type_mod, only: component_get_x2c_cx, component_get_c2x_cx
  use component_type_mod, only: glc, lnd

  implicit none
  save
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: prep_glc_init
  public :: prep_glc_mrg

  public :: prep_glc_accum
  public :: prep_glc_accum_avg

  public :: prep_glc_calc_l2x_gx

  public :: prep_glc_get_l2x_gx
  public :: prep_glc_get_l2gacc_lx
  public :: prep_glc_get_l2gacc_lx_cnt
  public :: prep_glc_get_mapper_SFl2g

  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

  private :: prep_glc_merge

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  ! mappers
  type(seq_map), pointer :: mapper_SFl2g

  ! attribute vectors 
  type(mct_aVect), pointer :: l2x_gx(:) ! Lnd export, glc grid, cpl pes - allocated in driver

  ! accumulation variables
  type(mct_aVect), pointer :: l2gacc_lx(:) ! Lnd export, lnd grid, cpl pes - allocated in driver
  integer        , target :: l2gacc_lx_cnt ! l2gacc_lx: number of time samples accumulated

  ! other module variables
  integer :: mpicom_CPLID  ! MPI cpl communicator
  !================================================================================================

contains

  !================================================================================================

  subroutine prep_glc_init(infodata, lnd_c2_glc)

    !---------------------------------------------------------------
    ! Description
    ! Initialize module attribute vectors and mapping variables
    !
    ! Arguments
    type (seq_infodata_type) , intent(inout) :: infodata
    logical                  , intent(in)    :: lnd_c2_glc ! .true.  => lnd to glc coupling on
    !
    ! Local Variables
    integer                          :: eli, egi
    integer                          :: lsize_l
    integer                          :: lsize_g
    logical                          :: esmf_map_flag ! .true. => use esmf for mapping
    logical                          :: iamroot_CPLID ! .true. => CPLID masterproc
    logical                          :: glc_present   ! .true. => glc is present
    character(CL)                    :: lnd_gnam      ! lnd grid
    character(CL)                    :: glc_gnam      ! glc grid
    type(mct_avect), pointer         :: l2x_lx
    type(mct_avect), pointer         :: x2g_gx
    character(*), parameter          :: subname = '(prep_glc_init)'
    character(*), parameter          :: F00 = "('"//subname//" : ', 4A )"
    !---------------------------------------------------------------

    call seq_infodata_getData(infodata , &
         esmf_map_flag=esmf_map_flag   , &
         glc_present=glc_present       , &
         lnd_gnam=lnd_gnam             , &
         glc_gnam=glc_gnam)

    allocate(mapper_SFl2g)

    if (glc_present) then

       call seq_comm_getData(CPLID, &
            mpicom=mpicom_CPLID, iamroot=iamroot_CPLID)

       l2x_lx => component_get_c2x_cx(lnd(1))
       lsize_l = mct_aVect_lsize(l2x_lx)
       
       x2g_gx => component_get_x2c_cx(glc(1))
       lsize_g = mct_aVect_lsize(x2g_gx)
       
       allocate(l2x_gx(num_inst_lnd))
       allocate(l2gacc_lx(num_inst_lnd))
       do eli = 1,num_inst_lnd
          call mct_aVect_initSharedFields(l2x_lx, x2g_gx, l2x_gx(eli) ,lsize=lsize_g)
          call mct_aVect_zero(l2x_gx(eli))
          
          call mct_aVect_initSharedFields(l2x_lx, x2g_gx, l2gacc_lx(eli), lsize=lsize_l)
          call mct_aVect_zero(l2gacc_lx(eli))
       enddo
       l2gacc_lx_cnt = 0
       
       if (iamroot_CPLID) then
          write(logunit,*) ' '
          write(logunit,F00) 'Initializing mapper_SFl2g'
       end if
       call seq_map_init_rearrolap(mapper_SFl2g, lnd(1), glc(1), 'mapper_SFl2g')
       call shr_sys_flush(logunit)

    end if

  end subroutine prep_glc_init

  !================================================================================================

! Subprogram not used   subroutine prep_glc_accum(timer)
! Subprogram not used 
! Subprogram not used     !---------------------------------------------------------------
! Subprogram not used     ! Description
! Subprogram not used     ! Accumulate glc inputs
! Subprogram not used     !
! Subprogram not used     ! Arguments
! Subprogram not used     character(len=*), intent(in) :: timer
! Subprogram not used     !
! Subprogram not used     ! Local Variables
! Subprogram not used     integer :: eli
! Subprogram not used     type(mct_avect), pointer :: l2x_lx
! Subprogram not used     character(*), parameter :: subname = '(prep_glc_accum)'
! Subprogram not used     !---------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
! Subprogram not used     do eli = 1,num_inst_lnd
! Subprogram not used        l2x_lx => component_get_c2x_cx(lnd(eli))
! Subprogram not used        if (l2gacc_lx_cnt == 0) then
! Subprogram not used           call mct_avect_copy(l2x_lx, l2gacc_lx(eli))
! Subprogram not used        else
! Subprogram not used           call mct_avect_accum(l2x_lx, l2gacc_lx(eli))
! Subprogram not used        endif
! Subprogram not used     end do
! Subprogram not used     l2gacc_lx_cnt = l2gacc_lx_cnt + 1
! Subprogram not used     call t_drvstopf  (trim(timer))
! Subprogram not used 
! Subprogram not used   end subroutine prep_glc_accum

  !================================================================================================

! Subprogram not used   subroutine prep_glc_accum_avg(timer)
! Subprogram not used 
! Subprogram not used     !---------------------------------------------------------------
! Subprogram not used     ! Description
! Subprogram not used     ! Finalize accumulation of glc inputs
! Subprogram not used     !
! Subprogram not used     ! Arguments
! Subprogram not used     character(len=*), intent(in) :: timer
! Subprogram not used     !
! Subprogram not used     ! Local Variables
! Subprogram not used     integer :: eli
! Subprogram not used     character(*), parameter :: subname = '(prep_glc_accum_avg)'
! Subprogram not used     !---------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
! Subprogram not used     if (l2gacc_lx_cnt > 1) then
! Subprogram not used        do eli = 1,num_inst_lnd
! Subprogram not used           call mct_avect_avg(l2gacc_lx(eli), l2gacc_lx_cnt)
! Subprogram not used        end do
! Subprogram not used     end if
! Subprogram not used     l2gacc_lx_cnt = 0
! Subprogram not used     call t_drvstopf  (trim(timer))
! Subprogram not used     
! Subprogram not used   end subroutine prep_glc_accum_avg

  !================================================================================================
  
! Subprogram not used   subroutine prep_glc_mrg(infodata, timer_mrg) 
! Subprogram not used 
! Subprogram not used     !---------------------------------------------------------------
! Subprogram not used     ! Description
! Subprogram not used     ! Merge glc inputs
! Subprogram not used     !
! Subprogram not used     ! Arguments
! Subprogram not used     type(seq_infodata_type) , intent(in)    :: infodata
! Subprogram not used     character(len=*)        , intent(in)    :: timer_mrg
! Subprogram not used     !
! Subprogram not used     ! Local Variables
! Subprogram not used     integer :: egi, eli
! Subprogram not used     type(mct_avect), pointer :: x2g_gx
! Subprogram not used     character(*), parameter  :: subname = '(prep_glc_mrg)'
! Subprogram not used     !---------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     call t_drvstartf (trim(timer_mrg),barrier=mpicom_CPLID)
! Subprogram not used     do egi = 1,num_inst_glc
! Subprogram not used        ! Use fortran mod to address ensembles in merge
! Subprogram not used        eli = mod((egi-1),num_inst_lnd) + 1
! Subprogram not used 
! Subprogram not used        x2g_gx => component_get_x2c_cx(glc(egi)) 
! Subprogram not used        call prep_glc_merge(l2x_gx(eli), x2g_gx)
! Subprogram not used     enddo
! Subprogram not used     call t_drvstopf  (trim(timer_mrg))
! Subprogram not used 
! Subprogram not used   end subroutine prep_glc_mrg

  !================================================================================================

! Subprogram not used   subroutine prep_glc_merge( s2x_g, x2g_g )
! Subprogram not used 
! Subprogram not used     !----------------------------------------------------------------------- 
! Subprogram not used     ! Arguments
! Subprogram not used     type(mct_aVect), intent(inout)  :: s2x_g  ! input
! Subprogram not used     type(mct_aVect), intent(inout)  :: x2g_g  ! output
! Subprogram not used     !----------------------------------------------------------------------- 
! Subprogram not used 
! Subprogram not used     integer       :: nflds,i,i1,o1
! Subprogram not used     logical       :: iamroot
! Subprogram not used     logical, save :: first_time = .true.
! Subprogram not used     character(CL),allocatable :: mrgstr(:)   ! temporary string
! Subprogram not used     character(CL) :: field   ! string converted to char
! Subprogram not used     type(mct_aVect_sharedindices),save :: s2x_sharedindices
! Subprogram not used     character(*), parameter   :: subname = '(prep_glc_merge) '
! Subprogram not used 
! Subprogram not used     !----------------------------------------------------------------------- 
! Subprogram not used 
! Subprogram not used     call seq_comm_getdata(CPLID, iamroot=iamroot)
! Subprogram not used 
! Subprogram not used     if (first_time) then
! Subprogram not used        nflds = mct_aVect_nRattr(x2g_g)
! Subprogram not used 
! Subprogram not used        allocate(mrgstr(nflds))
! Subprogram not used        do i = 1,nflds
! Subprogram not used           field = mct_aVect_getRList2c(i, x2g_g)
! Subprogram not used           mrgstr(i) = subname//'x2g%'//trim(field)//' ='
! Subprogram not used        enddo
! Subprogram not used 
! Subprogram not used        call mct_aVect_setSharedIndices(s2x_g, x2g_g, s2x_SharedIndices)
! Subprogram not used 
! Subprogram not used        !--- document copy operations ---
! Subprogram not used        do i=1,s2x_SharedIndices%shared_real%num_indices
! Subprogram not used           i1=s2x_SharedIndices%shared_real%aVindices1(i)
! Subprogram not used           o1=s2x_SharedIndices%shared_real%aVindices2(i)
! Subprogram not used           field = mct_aVect_getRList2c(i1, s2x_g)
! Subprogram not used           mrgstr(o1) = trim(mrgstr(o1))//' = s2x%'//trim(field)
! Subprogram not used        enddo
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     ! Create input glc state directly from land snow output state
! Subprogram not used     call mct_aVect_copy(aVin=s2x_g, aVout=x2g_g, vector=mct_usevector, sharedIndices=s2x_SharedIndices)
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
! Subprogram not used   end subroutine prep_glc_merge

  !================================================================================================

! Subprogram not used   subroutine prep_glc_calc_l2x_gx(timer)
! Subprogram not used     !---------------------------------------------------------------
! Subprogram not used     ! Description
! Subprogram not used     ! Create l2x_gx (note that l2x_gx is a local module variable)
! Subprogram not used     ! Also l2x_gx is really the accumulated l2xacc_lx mapped to l2x_gx
! Subprogram not used     !
! Subprogram not used     ! Arguments
! Subprogram not used     character(len=*), intent(in) :: timer
! Subprogram not used     !
! Subprogram not used     ! Local Variables
! Subprogram not used     integer :: eli
! Subprogram not used     character(*), parameter :: subname = '(prep_glc_calc_l2x_gx)'
! Subprogram not used     !---------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
! Subprogram not used     do eli = 1,num_inst_lnd
! Subprogram not used        call seq_map_map(mapper_SFl2g, l2gacc_lx(eli), l2x_gx(eli), norm=.true.)
! Subprogram not used     enddo
! Subprogram not used     call t_drvstopf  (trim(timer))
! Subprogram not used   end subroutine prep_glc_calc_l2x_gx

  !================================================================================================

! Subprogram not used   function prep_glc_get_l2x_gx()
! Subprogram not used     type(mct_aVect), pointer :: prep_glc_get_l2x_gx(:)
! Subprogram not used     prep_glc_get_l2x_gx => l2x_gx(:)   
! Subprogram not used   end function prep_glc_get_l2x_gx

! Subprogram not used   function prep_glc_get_l2gacc_lx()
! Subprogram not used     type(mct_aVect), pointer :: prep_glc_get_l2gacc_lx(:)
! Subprogram not used     prep_glc_get_l2gacc_lx => l2gacc_lx(:)   
! Subprogram not used   end function prep_glc_get_l2gacc_lx

! Subprogram not used   function prep_glc_get_l2gacc_lx_cnt()
! Subprogram not used     integer, pointer :: prep_glc_get_l2gacc_lx_cnt
! Subprogram not used     prep_glc_get_l2gacc_lx_cnt => l2gacc_lx_cnt
! Subprogram not used   end function prep_glc_get_l2gacc_lx_cnt

  function prep_glc_get_mapper_SFl2g()
    type(seq_map), pointer :: prep_glc_get_mapper_SFl2g
    prep_glc_get_mapper_SFl2g => mapper_SFl2g  
  end function prep_glc_get_mapper_SFl2g

end module prep_glc_mod
