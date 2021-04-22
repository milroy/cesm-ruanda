module pftdynMod

!---------------------------------------------------------------------------
!BOP
!
! !MODULE: pftdynMod
!
! !USES:
  use spmdMod
  use clmtype
  use decompMod   , only : get_proc_bounds
  use clm_varsur  , only : pctspec
  use clm_varpar  , only : max_pft_per_col
  use clm_varctl  , only : iulog, use_c13, use_cn, use_cndv
  use shr_sys_mod , only : shr_sys_flush
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils  , only : endrun
  use ncdio_pio   , only : file_desc_t, ncd_pio_openfile, ncd_inqdid, ncd_inqdlen, ncd_io, check_dim
!
! !DESCRIPTION:
! Determine pft weights at current time using dynamic landuse datasets.
! ASSUMES that only have one dynamic landuse dataset.
!
! !PUBLIC TYPES:
  implicit none
  private
  save
  public :: pftdyn_init
  public :: pftdyn_interp
  public :: pftdyn_wbal_init
  public :: pftdyn_wbal
  public :: pftdyn_cnbal
  public :: pftwt_init
  public :: pftwt_interp
  public :: CNHarvest
  public :: CNHarvestPftToColumn
!
! !REVISION HISTORY:
! Created by Peter Thornton
! slevis modified to handle CNDV and crop model
! 19 May 2009: PET - modified to handle harvest fluxes
!
!EOP
!
! ! PRIVATE TYPES
  integer , pointer   :: yearspft(:)
  real(r8), pointer   :: wtpft1(:,:)   
  real(r8), pointer   :: wtpft2(:,:)
  real(r8), pointer   :: harvest(:)   
  real(r8), pointer   :: wtcol_old(:)
  integer :: nt1
  integer :: nt2
  integer :: ntimes
  logical :: do_harvest
  type(file_desc_t)  :: ncid   ! netcdf id
!---------------------------------------------------------------------------

contains
  
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: pftdyn_init
!
! !INTERFACE:
! Subprogram not used   subroutine pftdyn_init()
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used ! Initialize dynamic landuse dataset (position it to the right time samples
! Subprogram not used ! that bound the initial model date)
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used     use clm_time_manager, only : get_curr_date
! Subprogram not used     use clm_varctl  , only : flanduse_timeseries
! Subprogram not used     use clm_varpar  , only : numpft, maxpatch_pft
! Subprogram not used     use fileutils   , only : getfil
! Subprogram not used !
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! !LOCAL VARIABLES:
! Subprogram not used !EOP
! Subprogram not used     integer  :: i,j,m,n,g                       ! indices
! Subprogram not used     real(r8) :: sumpct                          ! sum for error check
! Subprogram not used     integer  :: varid                           ! netcdf ids
! Subprogram not used     integer  :: year                            ! year (0, ...) for nstep+1
! Subprogram not used     integer  :: mon                             ! month (1, ..., 12) for nstep+1
! Subprogram not used     integer  :: day                             ! day of month (1, ..., 31) for nstep+1
! Subprogram not used     integer  :: sec                             ! seconds into current date for nstep+1
! Subprogram not used     integer  :: ier, ret                        ! error status
! Subprogram not used     logical  :: found                           ! true => input dataset bounding dates found
! Subprogram not used     logical  :: readvar	                        ! true => variable is on input dataset
! Subprogram not used     integer  :: begg,endg                       ! beg/end indices for land gridcells
! Subprogram not used     integer  :: begl,endl                       ! beg/end indices for land landunits
! Subprogram not used     integer  :: begc,endc                       ! beg/end indices for land columns
! Subprogram not used     integer  :: begp,endp                       ! beg/end indices for land pfts
! Subprogram not used     real(r8), pointer :: pctgla(:)          ! percent of gcell is glacier
! Subprogram not used     real(r8), pointer :: pctlak(:)          ! percent of gcell is lake
! Subprogram not used     real(r8), pointer :: pctwet(:)          ! percent of gcell is wetland
! Subprogram not used     real(r8), pointer :: pcturb(:)          ! percent of gcell is urbanized
! Subprogram not used     type(gridcell_type), pointer :: gptr        ! pointer to gridcell derived subtype
! Subprogram not used     character(len=256) :: locfn                 ! local file name
! Subprogram not used     character(len= 32) :: subname='pftdyn_init' ! subroutine name
! Subprogram not used  !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)
! Subprogram not used 
! Subprogram not used     ! Error check
! Subprogram not used 
! Subprogram not used     if ( maxpatch_pft /= numpft+1 )then
! Subprogram not used        call endrun( subname//' maxpatch_pft does NOT equal numpft+1 -- this is invalid for dynamic PFT case' )
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     allocate(pctgla(begg:endg),pctlak(begg:endg))
! Subprogram not used     allocate(pctwet(begg:endg),pcturb(begg:endg))
! Subprogram not used 
! Subprogram not used     ! Set pointers into derived type
! Subprogram not used 
! Subprogram not used     gptr => grc
! Subprogram not used 
! Subprogram not used     ! pctspec must be saved between time samples
! Subprogram not used     ! position to first time sample - assume that first time sample must match starting date
! Subprogram not used     ! check consistency -  special landunits, grid, frac and mask
! Subprogram not used     ! only do this once
! Subprogram not used 
! Subprogram not used     ! read data PCT_PFT corresponding to correct year
! Subprogram not used 
! Subprogram not used     allocate(wtpft1(begg:endg,0:numpft), wtpft2(begg:endg,0:numpft), stat=ier)
! Subprogram not used     if (ier /= 0) then
! Subprogram not used        call endrun( subname//' allocation error for wtpft1, wtpft2' )
! Subprogram not used     end if
! Subprogram not used     
! Subprogram not used     allocate(harvest(begg:endg),stat=ier)
! Subprogram not used     if (ier /= 0) then
! Subprogram not used        call endrun( subname//' allocation error for harvest')
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     allocate(wtcol_old(begp:endp),stat=ier)
! Subprogram not used     if (ier /= 0) then
! Subprogram not used        call endrun( subname//' allocation error for wtcol_old' )
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     if (masterproc) then
! Subprogram not used        write(iulog,*) 'Attempting to read pft dynamic landuse data .....'
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! Obtain file
! Subprogram not used     call getfil (flanduse_timeseries, locfn, 0)
! Subprogram not used     call ncd_pio_openfile (ncid, locfn, 0)
! Subprogram not used 
! Subprogram not used     ! Obtain pft years from dynamic landuse file
! Subprogram not used     
! Subprogram not used     call ncd_inqdid(ncid, 'time', varid)
! Subprogram not used     call ncd_inqdlen(ncid, varid, ntimes)
! Subprogram not used 
! Subprogram not used     ! Consistency check
! Subprogram not used     
! Subprogram not used     call check_dim(ncid, 'lsmpft', numpft+1)
! Subprogram not used 
! Subprogram not used     allocate (yearspft(ntimes), stat=ier)
! Subprogram not used     if (ier /= 0) then
! Subprogram not used        write(iulog,*) subname//' allocation error for yearspft'; call endrun()
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     call ncd_io(ncid=ncid, varname='YEAR', flag='read', data=yearspft)
! Subprogram not used 
! Subprogram not used     call ncd_io(ncid=ncid, varname='PCT_WETLAND', flag='read', data=pctwet, &
! Subprogram not used          dim1name=grlnd, readvar=readvar)
! Subprogram not used     if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_WETLAND NOT on landuse_timeseries file' )
! Subprogram not used 
! Subprogram not used     call ncd_io(ncid=ncid, varname= 'PCT_LAKE', flag='read', data=pctlak, &
! Subprogram not used          dim1name=grlnd, readvar=readvar)
! Subprogram not used     if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_LAKE NOT on landuse_timeseries file' )
! Subprogram not used 
! Subprogram not used     call ncd_io(ncid=ncid, varname= 'PCT_GLACIER', flag='read', data=pctgla, &
! Subprogram not used          dim1name=grlnd, readvar=readvar)
! Subprogram not used     if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_GLACIER NOT on landuse_timeseries file' )
! Subprogram not used 
! Subprogram not used     call ncd_io(ncid=ncid, varname= 'PCT_URBAN'  , flag='read', data=pcturb, &
! Subprogram not used          dim1name=grlnd, readvar=readvar)
! Subprogram not used     if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_URBAN NOT on landuse_timeseries file' )
! Subprogram not used 
! Subprogram not used     ! Consistency check
! Subprogram not used     do g = begg,endg
! Subprogram not used     !   this was causing a fail, even though values are the same to within 1e-15
! Subprogram not used     !   if (pctlak(g)+pctwet(g)+pcturb(g)+pctgla(g) /= pctspec(g)) then 
! Subprogram not used        if (abs((pctlak(g)+pctwet(g)+pcturb(g)+pctgla(g))-pctspec(g)) > 1e-13_r8) then 
! Subprogram not used           write(iulog,*) subname//'mismatch between input pctspec = ',&
! Subprogram not used                      pctlak(g)+pctwet(g)+pcturb(g)+pctgla(g),&
! Subprogram not used                     ' and that obtained from surface dataset ', pctspec(g),' at g= ',g
! Subprogram not used            call endrun()
! Subprogram not used        end if
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     ! Determine if current date spans the years
! Subprogram not used     ! If current year is less than first dynamic PFT timeseries year,
! Subprogram not used     ! then use the first year from dynamic pft file for both nt1 and nt2,
! Subprogram not used     ! forcing constant weights until the model year enters the dynamic
! Subprogram not used     ! pft dataset timeseries range.
! Subprogram not used     ! If current year is equal to or greater than the last dynamic pft
! Subprogram not used     ! timeseries year, then use the last year for both nt1 and nt2, 
! Subprogram not used     ! forcing constant weights for the remainder of the simulation.
! Subprogram not used     ! This mechanism permits the introduction of a dynamic pft period in the middle
! Subprogram not used     ! of a simulation, with constant weights before and after the dynamic period.
! Subprogram not used     ! PET: harvest - since harvest is specified as a rate for each year, this
! Subprogram not used     ! approach will not work. Instead, need to seta flag that indicates harvest is
! Subprogram not used     ! zero for the period before the beginning and after the end of the dynpft timeseries.
! Subprogram not used 
! Subprogram not used     call get_curr_date(year, mon, day, sec)
! Subprogram not used 
! Subprogram not used     if (year < yearspft(1)) then
! Subprogram not used        nt1 = 1
! Subprogram not used        nt2 = 1
! Subprogram not used        do_harvest = .false.
! Subprogram not used     else if (year >= yearspft(ntimes)) then
! Subprogram not used        nt1 = ntimes
! Subprogram not used        nt2 = ntimes
! Subprogram not used        do_harvest = .false.
! Subprogram not used     else
! Subprogram not used        found = .false.
! Subprogram not used        do n = 1,ntimes-1 
! Subprogram not used           if (year == yearspft(n)) then
! Subprogram not used              nt1 = n
! Subprogram not used              nt2 = nt1 + 1
! Subprogram not used              found = .true.
! Subprogram not used              do_harvest = .true.
! Subprogram not used           end if   
! Subprogram not used        end do
! Subprogram not used        if (.not. found) then
! Subprogram not used           write(iulog,*) subname//' error: model year not found in landuse_timeseries file'
! Subprogram not used           write(iulog,*)'model year = ',year
! Subprogram not used           call endrun()
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! Get pctpft time samples bracketing the current time
! Subprogram not used 
! Subprogram not used     if (masterproc) then
! Subprogram not used        write(iulog,*) 'Get PFTDYN data for year: ', yearspft(nt1)
! Subprogram not used     end if
! Subprogram not used     call pftdyn_getdata(nt1, wtpft1, begg,endg,0,numpft)
! Subprogram not used     if (masterproc) then
! Subprogram not used        write(iulog,*) 'Get PFTDYN data for year: ', yearspft(nt2)
! Subprogram not used     end if
! Subprogram not used     call pftdyn_getdata(nt2, wtpft2, begg,endg,0,numpft)
! Subprogram not used     
! Subprogram not used     if (use_cn) then
! Subprogram not used        ! Get harvest rate at the nt1 time
! Subprogram not used        call pftdyn_getharvest(nt1,begg,endg)
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! convert weights from percent to proportion
! Subprogram not used     do m = 0,numpft
! Subprogram not used        do g = begg,endg
! Subprogram not used           wtpft1(g,m) = wtpft1(g,m)/100._r8
! Subprogram not used           wtpft2(g,m) = wtpft2(g,m)/100._r8
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used        
! Subprogram not used     deallocate(pctgla,pctlak,pctwet,pcturb)
! Subprogram not used 
! Subprogram not used   end subroutine pftdyn_init

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: pftdyn_interp
!
! !INTERFACE:
! Subprogram not used   subroutine pftdyn_interp()
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used ! Time interpolate dynamic landuse data to get pft weights for model time
! Subprogram not used ! Note that harvest data are stored as rates (not weights) and so time interpolation is 
! Subprogram not used ! not necessary - the harvest rate is held constant through the year.  This is consistent with
! Subprogram not used ! the treatment of changing PFT weights, where interpolation of the annual endpoint weights leads to 
! Subprogram not used ! a constant rate of change in PFT weight through the year, with abrupt changes in the rate at
! Subprogram not used ! annual boundaries. This routine is still used to get the next harvest time slice, when needed.
! Subprogram not used ! This routine is also used to turn off the harvest switch when the model year runs past the end of
! Subprogram not used ! the dynpft time series.
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used     use clm_time_manager, only : get_curr_date, get_curr_calday, &
! Subprogram not used                                  get_days_per_year
! Subprogram not used     use clm_varcon      , only : istsoil
! Subprogram not used     use clm_varcon      , only : istcrop
! Subprogram not used     use clm_varpar      , only : numpft
! Subprogram not used     implicit none
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! !LOCAL VARIABLES:
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used     integer  :: begg,endg        ! beg/end indices for land gridcells
! Subprogram not used     integer  :: begl,endl        ! beg/end indices for land landunits
! Subprogram not used     integer  :: begc,endc        ! beg/end indices for land columns
! Subprogram not used     integer  :: begp,endp        ! beg/end indices for land pfts
! Subprogram not used     integer  :: i,j,m,p,l,g,c    ! indices
! Subprogram not used     integer  :: year             ! year (0, ...) for nstep+1
! Subprogram not used     integer  :: mon              ! month (1, ..., 12) for nstep+1
! Subprogram not used     integer  :: day              ! day of month (1, ..., 31) for nstep+1
! Subprogram not used     integer  :: sec              ! seconds into current date for nstep+1
! Subprogram not used     real(r8) :: cday             ! current calendar day (1.0 = 0Z on Jan 1)
! Subprogram not used     real(r8) :: days_per_year    ! days per year
! Subprogram not used     integer  :: ier              ! error status
! Subprogram not used     integer  :: lbc,ubc
! Subprogram not used     real(r8) :: wt1              ! time interpolation weights
! Subprogram not used     real(r8), pointer :: wtpfttot1(:)            ! summation of pft weights for renormalization
! Subprogram not used     real(r8), pointer :: wtpfttot2(:)            ! summation of pft weights for renormalization
! Subprogram not used     real(r8), parameter :: wtpfttol = 1.e-10     ! tolerance for pft weight renormalization
! Subprogram not used     type(gridcell_type), pointer :: gptr         ! pointer to gridcell derived subtype
! Subprogram not used     type(landunit_type), pointer :: lptr         ! pointer to landunit derived subtype
! Subprogram not used     type(pft_type)     , pointer :: pptr         ! pointer to pft derived subtype
! Subprogram not used     character(len=32) :: subname='pftdyn_interp' ! subroutine name
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)
! Subprogram not used 
! Subprogram not used     ! Set pointers into derived type
! Subprogram not used 
! Subprogram not used     gptr => grc
! Subprogram not used     lptr => lun
! Subprogram not used     pptr => pft
! Subprogram not used 
! Subprogram not used     allocate(wtpfttot1(begc:endc),wtpfttot2(begc:endc))
! Subprogram not used     wtpfttot1(:) = 0._r8
! Subprogram not used     wtpfttot2(:) = 0._r8
! Subprogram not used 
! Subprogram not used     ! Interpolat pctpft to current time step - output in pctpft_intp
! Subprogram not used     ! Map interpolated pctpft to subgrid weights
! Subprogram not used     ! assumes that maxpatch_pft = numpft + 1, that each landunit has only 1 column, 
! Subprogram not used     ! SCAM and CNDV have not been defined, and create_croplandunit = .false.
! Subprogram not used 
! Subprogram not used     ! If necessary, obtain new time sample
! Subprogram not used 
! Subprogram not used     ! Get current date
! Subprogram not used 
! Subprogram not used     call get_curr_date(year, mon, day, sec)
! Subprogram not used 
! Subprogram not used     ! Obtain new time sample if necessary.
! Subprogram not used     ! The first condition is the regular crossing of a year boundary
! Subprogram not used     ! when within the dynpft timeseries range. The second condition is
! Subprogram not used     ! the case of the first entry into the dynpft timeseries range from
! Subprogram not used     ! an earlier period of constant weights.
! Subprogram not used 
! Subprogram not used     if (year > yearspft(nt1) .or. (nt1 == 1 .and. nt2 == 1 .and. year == yearspft(1))) then
! Subprogram not used 
! Subprogram not used        if (year >= yearspft(ntimes)) then
! Subprogram not used           nt1 = ntimes
! Subprogram not used           nt2 = ntimes
! Subprogram not used        else
! Subprogram not used           nt1        = nt2
! Subprogram not used           nt2        = nt2 + 1
! Subprogram not used           do_harvest = .true.
! Subprogram not used        end if
! Subprogram not used        
! Subprogram not used        if (year > yearspft(ntimes)) then
! Subprogram not used           do_harvest = .false.
! Subprogram not used        endif
! Subprogram not used        
! Subprogram not used        if (nt2 > ntimes .and. masterproc) then
! Subprogram not used           write(iulog,*)subname,' error - current year is past input data boundary'
! Subprogram not used        end if
! Subprogram not used        
! Subprogram not used        do m = 0,numpft
! Subprogram not used           do g = begg,endg
! Subprogram not used              wtpft1(g,m) = wtpft2(g,m)
! Subprogram not used           end do
! Subprogram not used        end do
! Subprogram not used 
! Subprogram not used        if (masterproc) then
! Subprogram not used           write(iulog,*) 'Get PFTDYN data for year: ', yearspft(nt2)
! Subprogram not used        end if
! Subprogram not used        call pftdyn_getdata(nt2, wtpft2, begg,endg,0,numpft)
! Subprogram not used        if (use_cn) then
! Subprogram not used           call pftdyn_getharvest(nt1,begg,endg)
! Subprogram not used        end if
! Subprogram not used 
! Subprogram not used        do m = 0,numpft
! Subprogram not used           do g = begg,endg
! Subprogram not used              wtpft2(g,m) = wtpft2(g,m)/100._r8
! Subprogram not used           end do
! Subprogram not used        end do
! Subprogram not used     
! Subprogram not used     end if  ! end of need new data if-block 
! Subprogram not used 
! Subprogram not used     ! Interpolate pft weight to current time
! Subprogram not used 
! Subprogram not used     cday          = get_curr_calday() 
! Subprogram not used     days_per_year = get_days_per_year()
! Subprogram not used 
! Subprogram not used     wt1 = ((days_per_year + 1._r8) - cday)/days_per_year
! Subprogram not used 
! Subprogram not used     do p = begp,endp
! Subprogram not used        c = pptr%column(p)
! Subprogram not used        g = pptr%gridcell(p)
! Subprogram not used        l = pptr%landunit(p)
! Subprogram not used        if (lptr%itype(l) == istsoil .or. lptr%itype(l) == istcrop) then
! Subprogram not used           m = pptr%itype(p)
! Subprogram not used           wtcol_old(p)      = pptr%wtcol(p)
! Subprogram not used !         --- recoded for roundoff performance, tcraig 3/07 from k.lindsay
! Subprogram not used !         pptr%wtgcell(p)   = wtpft1(g,m)*wt1 + wtpft2(g,m)*wt2
! Subprogram not used           wtpfttot1(c) = wtpfttot1(c)+pptr%wtgcell(p)    
! Subprogram not used           pptr%wtgcell(p)   = wtpft2(g,m) + wt1*(wtpft1(g,m)-wtpft2(g,m))
! Subprogram not used           pptr%wtlunit(p)   = pptr%wtgcell(p) / lptr%wtgcell(l)
! Subprogram not used           pptr%wtcol(p)     = pptr%wtgcell(p) / lptr%wtgcell(l)
! Subprogram not used           wtpfttot2(c) = wtpfttot2(c)+pptr%wtgcell(p)
! Subprogram not used        end if
! Subprogram not used 
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used !   Renormalize pft weights so that sum of pft weights relative to grid cell 
! Subprogram not used !   remain constant even as land cover changes.  Doing this eliminates 
! Subprogram not used !   soil balance error warnings.  (DML, 4/8/2009)
! Subprogram not used     do p = begp,endp
! Subprogram not used        c = pptr%column(p)
! Subprogram not used        g = pptr%gridcell(p)
! Subprogram not used        l = pptr%landunit(p)
! Subprogram not used        if (lptr%itype(l) == istsoil .or. lptr%itype(l) == istcrop) then
! Subprogram not used           if (wtpfttot2(c) /= 0 .and. &
! Subprogram not used               abs(wtpfttot1(c)-wtpfttot2(c)) > wtpfttol) then
! Subprogram not used              pptr%wtgcell(p)   = (wtpfttot1(c)/wtpfttot2(c))*pptr%wtgcell(p)
! Subprogram not used              pptr%wtlunit(p)   = pptr%wtgcell(p) / lptr%wtgcell(l)
! Subprogram not used              pptr%wtcol(p)     = pptr%wtgcell(p) / lptr%wtgcell(l)
! Subprogram not used           end if
! Subprogram not used        end if
! Subprogram not used 
! Subprogram not used     end do
! Subprogram not used    
! Subprogram not used     deallocate(wtpfttot1,wtpfttot2) 
! Subprogram not used     
! Subprogram not used   end subroutine pftdyn_interp

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: pftdyn_getdata
!
! !INTERFACE:
! Subprogram not used   subroutine pftdyn_getdata(ntime, pctpft, begg, endg, pft0, maxpft)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used ! Obtain dynamic landuse data (pctpft) and make sure that
! Subprogram not used ! percentage of PFTs sum to 100% cover for vegetated landunit
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used     use clm_varpar  , only : numpft
! Subprogram not used !
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     integer , intent(in)  :: ntime
! Subprogram not used     integer , intent(in)  :: begg,endg,pft0,maxpft
! Subprogram not used     real(r8), intent(out) :: pctpft(begg:endg,pft0:maxpft)
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! !LOCAL VARIABLES:
! Subprogram not used !EOP
! Subprogram not used     integer  :: i,j,m,n
! Subprogram not used     integer  :: err, ierr, ret
! Subprogram not used     real(r8) :: sumpct,sumerr                     ! temporary
! Subprogram not used     real(r8), pointer :: arrayl(:,:)              ! temporary array
! Subprogram not used     logical  :: readvar
! Subprogram not used    
! Subprogram not used     character(len=32) :: subname='pftdyn_getdata' ! subroutine name
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used     
! Subprogram not used     allocate(arrayl(begg:endg,pft0:maxpft))	
! Subprogram not used     call ncd_io(ncid=ncid, varname= 'PCT_PFT', flag='read', data=arrayl, &
! Subprogram not used          dim1name=grlnd, nt=ntime, readvar=readvar)
! Subprogram not used     pctpft(begg:endg,pft0:maxpft) = arrayl(begg:endg,pft0:maxpft)
! Subprogram not used     deallocate(arrayl)		
! Subprogram not used     if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_PFT NOT on pftdyn file' )
! Subprogram not used 
! Subprogram not used     err = 0
! Subprogram not used     do n = begg,endg
! Subprogram not used        if (pctspec(n) < 100._r8) then
! Subprogram not used           sumpct = 0._r8
! Subprogram not used           do m = 0, numpft
! Subprogram not used              sumpct = sumpct + pctpft(n,m) * 100._r8/(100._r8-pctspec(n))
! Subprogram not used           end do
! Subprogram not used           if (abs(sumpct - 100._r8) > 0.1_r8) then
! Subprogram not used              err = 1; ierr = n; sumerr = sumpct
! Subprogram not used           end if
! Subprogram not used           if (sumpct < -0.000001_r8) then
! Subprogram not used              err = 2; ierr = n; sumerr = sumpct
! Subprogram not used           end if
! Subprogram not used        end if
! Subprogram not used     end do
! Subprogram not used     if (err == 1) then
! Subprogram not used        write(iulog,*) subname,' error: sum(pct) over numpft+1 is not = 100.',sumerr,ierr,pctspec(ierr),pctpft(ierr,:)
! Subprogram not used        call endrun()
! Subprogram not used     else if (err == 2) then
! Subprogram not used        write(iulog,*)subname,' error: sum(pct) over numpft+1 is < 0.',sumerr,ierr,pctspec(ierr),pctpft(ierr,:)
! Subprogram not used        call endrun()
! Subprogram not used     end if
! Subprogram not used     
! Subprogram not used   end subroutine pftdyn_getdata

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: pftdyn_getharvest
!
! !INTERFACE:
! Subprogram not used   subroutine pftdyn_getharvest(ntime, begg, endg)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used ! Obtain harvest data 
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     integer , intent(in)  :: ntime
! Subprogram not used     integer , intent(IN)  :: begg     ! beg indices for land gridcells
! Subprogram not used     integer , intent(IN)  :: endg     ! end indices for land gridcells
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! !LOCAL VARIABLES:
! Subprogram not used !EOP
! Subprogram not used     real(r8), pointer :: arrayl(:)                   ! temporary array
! Subprogram not used     logical :: readvar 
! Subprogram not used     character(len=32) :: subname='pftdyn_getharvest' ! subroutine name
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used     
! Subprogram not used     allocate(arrayl(begg:endg))
! Subprogram not used     
! Subprogram not used     call ncd_io(ncid=ncid, varname= 'HARVEST_VH1', flag='read', data=arrayl, dim1name=grlnd, &
! Subprogram not used          nt=ntime, readvar=readvar)
! Subprogram not used     if (.not. readvar) call endrun( trim(subname)//' ERROR: HARVEST_VH1 not on landuse_timeseries file' )
! Subprogram not used     harvest(begg:endg) = arrayl(begg:endg)
! Subprogram not used     
! Subprogram not used     call ncd_io(ncid=ncid, varname= 'HARVEST_VH2', flag='read', data=arrayl, dim1name=grlnd, &
! Subprogram not used          nt=ntime, readvar=readvar)
! Subprogram not used     if (.not. readvar) call endrun( trim(subname)//' ERROR: HARVEST_VH2 not on landuse_timeseries file' )
! Subprogram not used     harvest(begg:endg) = harvest(begg:endg) + arrayl(begg:endg)
! Subprogram not used     
! Subprogram not used     call ncd_io(ncid=ncid, varname= 'HARVEST_SH1', flag='read', data=arrayl, dim1name=grlnd, &
! Subprogram not used          nt=ntime, readvar=readvar)
! Subprogram not used     if (.not. readvar) call endrun( trim(subname)//' ERROR: HARVEST_SH1 not on landuse_timeseries file' )
! Subprogram not used     harvest(begg:endg) = harvest(begg:endg) + arrayl(begg:endg)
! Subprogram not used     
! Subprogram not used     call ncd_io(ncid=ncid, varname= 'HARVEST_SH2', flag='read', data=arrayl, dim1name=grlnd, &
! Subprogram not used          nt=ntime, readvar=readvar)
! Subprogram not used     if (.not. readvar) call endrun( trim(subname)//' ERROR: HARVEST_SH2 not on landuse_timeseries file' )
! Subprogram not used     harvest(begg:endg) = harvest(begg:endg) + arrayl(begg:endg)
! Subprogram not used     
! Subprogram not used     call ncd_io(ncid=ncid, varname= 'HARVEST_SH3', flag='read', data=arrayl, dim1name=grlnd, &
! Subprogram not used          nt=ntime, readvar=readvar)
! Subprogram not used     if (.not. readvar) call endrun( trim(subname)//' ERROR: HARVEST_SH3 not on landuse_timeseries file' )
! Subprogram not used     harvest(begg:endg) = harvest(begg:endg) + arrayl(begg:endg)
! Subprogram not used 
! Subprogram not used     deallocate(arrayl)
! Subprogram not used 
! Subprogram not used   end subroutine pftdyn_getharvest

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: pftdyn_wbal_init
!
! !INTERFACE:
  subroutine pftdyn_wbal_init( begc, endc )
!
! !DESCRIPTION:
! initialize the column-level mass-balance correction term.
! Called in every timestep.
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    integer, intent(IN)  :: begc, endc    ! proc beginning and ending column indices
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: c             ! indices
    type(column_type),   pointer :: cptr         ! pointer to column derived subtype
!-----------------------------------------------------------------------

    ! Set pointers into derived type

    cptr => col

    ! set column-level canopy water mass balance correction flux
    ! term to 0 at the beginning of every timestep
    
    do c = begc,endc
       cwf%h2ocan_loss(c) = 0._r8
    end do
    
  end subroutine pftdyn_wbal_init

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: pftdyn_wbal
!
! !INTERFACE:
! Subprogram not used   subroutine pftdyn_wbal( begg, endg, begc, endc, begp, endp )
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used ! modify pft-level state and flux variables to maintain water balance with
! Subprogram not used ! dynamic pft-weights.
! Subprogram not used ! Canopy water balance does not need to consider harvest fluxes, since pft weights are
! Subprogram not used ! not affected by harvest.
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used     use clm_varcon  , only : istsoil
! Subprogram not used     use clm_varcon  , only : istcrop
! Subprogram not used     use clm_time_manager, only : get_step_size
! Subprogram not used !
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     integer, intent(IN)  :: begg     ! beg indices for land gridcells
! Subprogram not used     integer, intent(IN)  :: endg     ! end indices for land gridcells
! Subprogram not used     integer, intent(IN)  :: begc     ! beg indices for land columns
! Subprogram not used     integer, intent(IN)  :: endc     ! end indices for land columns
! Subprogram not used     integer, intent(IN)  :: begp     ! beg indices for land plant function types
! Subprogram not used     integer, intent(IN)  :: endp     ! end indices for land plant function types
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! !LOCAL VARIABLES:
! Subprogram not used !EOP
! Subprogram not used     integer  :: pi,p,c,l,g    ! indices
! Subprogram not used     integer  :: ier           ! error code
! Subprogram not used     real(r8) :: dtime         ! land model time step (sec)
! Subprogram not used     real(r8) :: dwt           ! change in pft weight (relative to column)
! Subprogram not used     real(r8) :: init_h2ocan   ! initial canopy water mass
! Subprogram not used     real(r8) :: new_h2ocan    ! canopy water mass after weight shift
! Subprogram not used     real(r8), allocatable :: loss_h2ocan(:) ! canopy water mass loss due to weight shift
! Subprogram not used     type(landunit_type), pointer :: lptr         ! pointer to landunit derived subtype
! Subprogram not used     type(column_type),   pointer :: cptr         ! pointer to column derived subtype
! Subprogram not used     type(pft_type)   ,   pointer :: pptr         ! pointer to pft derived subtype
! Subprogram not used     character(len=32) :: subname='pftdyn_wbal' ! subroutine name
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     ! Set pointers into derived type
! Subprogram not used 
! Subprogram not used     lptr => lun
! Subprogram not used     cptr => col
! Subprogram not used     pptr => pft
! Subprogram not used 
! Subprogram not used     ! Allocate loss_h2ocan
! Subprogram not used     allocate(loss_h2ocan(begp:endp), stat=ier)
! Subprogram not used     if (ier /= 0) then
! Subprogram not used           write(iulog,*)subname,' allocation error for loss_h2ocan'; call endrun()
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! Get time step
! Subprogram not used 
! Subprogram not used     dtime = get_step_size()
! Subprogram not used 
! Subprogram not used     ! set column-level canopy water mass balance correction flux
! Subprogram not used     ! term to 0 at the beginning of every weight-shifting timestep
! Subprogram not used 
! Subprogram not used     do c = begc,endc
! Subprogram not used        cwf%h2ocan_loss(c) = 0._r8 ! is this OR pftdyn_wbal_init redundant?
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     do p = begp,endp
! Subprogram not used        l = pptr%landunit(p)
! Subprogram not used        loss_h2ocan(p) = 0._r8
! Subprogram not used 
! Subprogram not used        if (lptr%itype(l) == istsoil .or. lptr%itype(l) == istcrop) then
! Subprogram not used 
! Subprogram not used           ! calculate the change in weight for the timestep
! Subprogram not used           dwt = pptr%wtcol(p)-wtcol_old(p)
! Subprogram not used   
! Subprogram not used           if (dwt > 0._r8) then
! Subprogram not used           
! Subprogram not used              ! if the pft gained weight, then the 
! Subprogram not used              ! initial canopy water state is redistributed over the
! Subprogram not used              ! new (larger) area, conserving mass.
! Subprogram not used 
! Subprogram not used              pws%h2ocan(p) = pws%h2ocan(p) * (wtcol_old(p)/pptr%wtcol(p))
! Subprogram not used           
! Subprogram not used           else if (dwt < 0._r8) then
! Subprogram not used           
! Subprogram not used              ! if the pft lost weight on the timestep, then the canopy water
! Subprogram not used              ! mass associated with the lost weight is directed to a 
! Subprogram not used              ! column-level flux term that gets added to the precip flux
! Subprogram not used              ! for every pft calculation in Hydrology1()
! Subprogram not used              
! Subprogram not used              init_h2ocan = pws%h2ocan(p) * wtcol_old(p)
! Subprogram not used              loss_h2ocan(p) = pws%h2ocan(p) * (-dwt)
! Subprogram not used              new_h2ocan = init_h2ocan - loss_h2ocan(p)
! Subprogram not used              if (abs(new_h2ocan) < 1e-8_r8) then
! Subprogram not used                 new_h2ocan = 0._r8
! Subprogram not used                 loss_h2ocan(p) = init_h2ocan
! Subprogram not used              end if
! Subprogram not used              if (pptr%wtcol(p) /= 0._r8) then  
! Subprogram not used                 pws%h2ocan(p) = new_h2ocan/pptr%wtcol(p)
! Subprogram not used              else
! Subprogram not used                 pws%h2ocan(p) = 0._r8
! Subprogram not used                 loss_h2ocan(p) = init_h2ocan
! Subprogram not used              end if 
! Subprogram not used        
! Subprogram not used 
! Subprogram not used           end if
! Subprogram not used 
! Subprogram not used        end if
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     do pi = 1,max_pft_per_col
! Subprogram not used        do c = begc,endc
! Subprogram not used           if (pi <= cptr%npfts(c)) then
! Subprogram not used              p = cptr%pfti(c) + pi - 1
! Subprogram not used              cwf%h2ocan_loss(c) = cwf%h2ocan_loss(c) + loss_h2ocan(p)/dtime
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     ! Deallocate loss_h2ocan
! Subprogram not used     deallocate(loss_h2ocan)
! Subprogram not used     
! Subprogram not used   end subroutine pftdyn_wbal
  
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: pftdyn_cnbal
!
! !INTERFACE:
! Subprogram not used   subroutine pftdyn_cnbal( begc, endc, begp, endp )
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used ! modify pft-level state and flux variables to maintain carbon and nitrogen balance with
! Subprogram not used ! dynamic pft-weights.
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used     use shr_kind_mod, only : r8 => shr_kind_r8
! Subprogram not used     use shr_const_mod,only : SHR_CONST_PDB
! Subprogram not used     use decompMod   , only : get_proc_bounds
! Subprogram not used     use clm_varcon  , only : istsoil
! Subprogram not used     use clm_varcon  , only : istcrop
! Subprogram not used     use clm_varpar  , only : numveg, numpft
! Subprogram not used     use pftvarcon   , only : pconv, pprod10, pprod100
! Subprogram not used     use clm_varcon  , only : c13ratio
! Subprogram not used     use clm_time_manager, only : get_step_size
! Subprogram not used !
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     integer, intent(IN)  :: begp, endp    ! proc beginning and ending pft indices
! Subprogram not used     integer, intent(IN)  :: begc, endc    ! proc beginning and ending column indices
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! !LOCAL VARIABLES:
! Subprogram not used !EOP
! Subprogram not used     integer  :: pi,p,c,l,g    ! indices
! Subprogram not used     integer  :: ier           ! error code
! Subprogram not used     real(r8) :: dwt           ! change in pft weight (relative to column)
! Subprogram not used     real(r8) :: dt            ! land model time step (sec)
! Subprogram not used     real(r8) :: init_h2ocan   ! initial canopy water mass
! Subprogram not used     real(r8) :: new_h2ocan    ! canopy water mass after weight shift
! Subprogram not used     real(r8), allocatable :: dwt_leafc_seed(:)       ! pft-level mass gain due to seeding of new area
! Subprogram not used     real(r8), allocatable :: dwt_leafn_seed(:)       ! pft-level mass gain due to seeding of new area
! Subprogram not used     real(r8), allocatable :: dwt_leafc13_seed(:)     ! pft-level mass gain due to seeding of new area
! Subprogram not used     real(r8), allocatable :: dwt_deadstemc_seed(:)       ! pft-level mass gain due to seeding of new area
! Subprogram not used     real(r8), allocatable :: dwt_deadstemn_seed(:)       ! pft-level mass gain due to seeding of new area
! Subprogram not used     real(r8), allocatable :: dwt_deadstemc13_seed(:)     ! pft-level mass gain due to seeding of new area
! Subprogram not used     real(r8), allocatable :: dwt_frootc_to_litter(:)       ! pft-level mass loss due to weight shift
! Subprogram not used     real(r8), allocatable :: dwt_livecrootc_to_litter(:)   ! pft-level mass loss due to weight shift
! Subprogram not used     real(r8), allocatable :: dwt_deadcrootc_to_litter(:)   ! pft-level mass loss due to weight shift
! Subprogram not used     real(r8), allocatable, target :: dwt_frootc13_to_litter(:)     ! pft-level mass loss due to weight shift
! Subprogram not used     real(r8), allocatable, target :: dwt_livecrootc13_to_litter(:) ! pft-level mass loss due to weight shift
! Subprogram not used     real(r8), allocatable, target :: dwt_deadcrootc13_to_litter(:) ! pft-level mass loss due to weight shift
! Subprogram not used     real(r8), allocatable, target :: dwt_frootn_to_litter(:)       ! pft-level mass loss due to weight shift
! Subprogram not used     real(r8), allocatable, target :: dwt_livecrootn_to_litter(:)   ! pft-level mass loss due to weight shift
! Subprogram not used     real(r8), allocatable, target :: dwt_deadcrootn_to_litter(:)   ! pft-level mass loss due to weight shift
! Subprogram not used     real(r8), allocatable :: conv_cflux(:)         ! pft-level mass loss due to weight shift
! Subprogram not used     real(r8), allocatable :: prod10_cflux(:)       ! pft-level mass loss due to weight shift
! Subprogram not used     real(r8), allocatable :: prod100_cflux(:)      ! pft-level mass loss due to weight shift
! Subprogram not used     real(r8), allocatable, target :: conv_c13flux(:)       ! pft-level mass loss due to weight shift
! Subprogram not used     real(r8), allocatable, target :: prod10_c13flux(:)     ! pft-level mass loss due to weight shift
! Subprogram not used     real(r8), allocatable, target :: prod100_c13flux(:)    ! pft-level mass loss due to weight shift
! Subprogram not used     real(r8), allocatable, target :: conv_nflux(:)         ! pft-level mass loss due to weight shift
! Subprogram not used     real(r8), allocatable, target :: prod10_nflux(:)       ! pft-level mass loss due to weight shift
! Subprogram not used     real(r8), allocatable, target :: prod100_nflux(:)      ! pft-level mass loss due to weight shift
! Subprogram not used     real(r8) :: c3_del13c     ! typical del13C for C3 photosynthesis (permil, relative to PDB)
! Subprogram not used     real(r8) :: c4_del13c     ! typical del13C for C4 photosynthesis (permil, relative to PDB)
! Subprogram not used     real(r8) :: c3_r1         ! isotope ratio (13c/12c) for C3 photosynthesis
! Subprogram not used     real(r8) :: c4_r1         ! isotope ratio (13c/12c) for C4 photosynthesis
! Subprogram not used     real(r8) :: c3_r2         ! isotope ratio (13c/[12c+13c]) for C3 photosynthesis
! Subprogram not used     real(r8) :: c4_r2         ! isotope ratio (13c/[12c+13c]) for C4 photosynthesis
! Subprogram not used     real(r8) :: t1,t2,wt_new,wt_old
! Subprogram not used     real(r8) :: init_state, change_state, new_state
! Subprogram not used     real(r8) :: tot_leaf, pleaf, pstor, pxfer
! Subprogram not used     real(r8) :: leafc_seed, leafn_seed
! Subprogram not used     real(r8) :: deadstemc_seed, deadstemn_seed
! Subprogram not used     real(r8) :: leafc13_seed, deadstemc13_seed
! Subprogram not used     real(r8), pointer :: dwt_ptr0, dwt_ptr1, dwt_ptr2, dwt_ptr3, ptr
! Subprogram not used     type(landunit_type), pointer :: lptr         ! pointer to landunit derived subtype
! Subprogram not used     type(column_type),   pointer :: cptr         ! pointer to column derived subtype
! Subprogram not used     type(pft_type)   ,   pointer :: pptr         ! pointer to pft derived subtype
! Subprogram not used     integer          , pointer   :: pcolumn(:)   ! column of corresponding pft
! Subprogram not used     character(len=32) :: subname='pftdyn_cbal' ! subroutine name
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used     
! Subprogram not used     ! Set pointers into derived type
! Subprogram not used 
! Subprogram not used     lptr    => lun
! Subprogram not used     cptr    => col
! Subprogram not used     pptr    => pft
! Subprogram not used     pcolumn => pptr%column
! Subprogram not used 
! Subprogram not used     ! Allocate pft-level mass loss arrays
! Subprogram not used     allocate(dwt_leafc_seed(begp:endp), stat=ier)
! Subprogram not used     if (ier /= 0) then
! Subprogram not used           write(iulog,*)subname,' allocation error for dwt_leafc_seed'; call endrun()
! Subprogram not used     end if
! Subprogram not used     allocate(dwt_leafn_seed(begp:endp), stat=ier)
! Subprogram not used     if (ier /= 0) then
! Subprogram not used           write(iulog,*)subname,' allocation error for dwt_leafn_seed'; call endrun()
! Subprogram not used     end if
! Subprogram not used     if (use_c13) then
! Subprogram not used        allocate(dwt_leafc13_seed(begp:endp), stat=ier)
! Subprogram not used        if (ier /= 0) then
! Subprogram not used           write(iulog,*)subname,' allocation error for dwt_leafc13_seed'; call endrun()
! Subprogram not used        end if
! Subprogram not used     endif
! Subprogram not used     allocate(dwt_deadstemc_seed(begp:endp), stat=ier)
! Subprogram not used     if (ier /= 0) then
! Subprogram not used           write(iulog,*)subname,' allocation error for dwt_deadstemc_seed'; call endrun()
! Subprogram not used     end if
! Subprogram not used     allocate(dwt_deadstemn_seed(begp:endp), stat=ier)
! Subprogram not used     if (ier /= 0) then
! Subprogram not used           write(iulog,*)subname,' allocation error for dwt_deadstemn_seed'; call endrun()
! Subprogram not used     end if
! Subprogram not used     if (use_c13) then
! Subprogram not used        allocate(dwt_deadstemc13_seed(begp:endp), stat=ier)
! Subprogram not used        if (ier /= 0) then
! Subprogram not used           write(iulog,*)subname,' allocation error for dwt_deadstemc13_seed'; call endrun()
! Subprogram not used        end if
! Subprogram not used     endif
! Subprogram not used     allocate(dwt_frootc_to_litter(begp:endp), stat=ier)
! Subprogram not used     if (ier /= 0) then
! Subprogram not used           write(iulog,*)subname,' allocation error for dwt_frootc_to_litter'; call endrun()
! Subprogram not used     end if
! Subprogram not used     allocate(dwt_livecrootc_to_litter(begp:endp), stat=ier)
! Subprogram not used     if (ier /= 0) then
! Subprogram not used           write(iulog,*)subname,' allocation error for dwt_livecrootc_to_litter'; call endrun()
! Subprogram not used     end if
! Subprogram not used     allocate(dwt_deadcrootc_to_litter(begp:endp), stat=ier)
! Subprogram not used     if (ier /= 0) then
! Subprogram not used        write(iulog,*)subname,' allocation error for dwt_deadcrootc_to_litter'; call endrun()
! Subprogram not used     end if
! Subprogram not used     if (use_c13) then
! Subprogram not used        allocate(dwt_frootc13_to_litter(begp:endp), stat=ier)
! Subprogram not used        if (ier /= 0) then
! Subprogram not used           write(iulog,*)subname,' allocation error for dwt_frootc13_to_litter'; call endrun()
! Subprogram not used        end if
! Subprogram not used        allocate(dwt_livecrootc13_to_litter(begp:endp), stat=ier)
! Subprogram not used        if (ier /= 0) then
! Subprogram not used           write(iulog,*)subname,' allocation error for dwt_livecrootc13_to_litter'; call endrun()
! Subprogram not used        end if
! Subprogram not used        allocate(dwt_deadcrootc13_to_litter(begp:endp), stat=ier)
! Subprogram not used        if (ier /= 0) then
! Subprogram not used           write(iulog,*)subname,' allocation error for dwt_deadcrootc13_to_litter'; call endrun()
! Subprogram not used        end if
! Subprogram not used     endif
! Subprogram not used     allocate(dwt_frootn_to_litter(begp:endp), stat=ier)
! Subprogram not used     if (ier /= 0) then
! Subprogram not used           write(iulog,*)subname,' allocation error for dwt_frootn_to_litter'; call endrun()
! Subprogram not used     end if
! Subprogram not used     allocate(dwt_livecrootn_to_litter(begp:endp), stat=ier)
! Subprogram not used     if (ier /= 0) then
! Subprogram not used           write(iulog,*)subname,' allocation error for dwt_livecrootn_to_litter'; call endrun()
! Subprogram not used     end if
! Subprogram not used     allocate(dwt_deadcrootn_to_litter(begp:endp), stat=ier)
! Subprogram not used     if (ier /= 0) then
! Subprogram not used           write(iulog,*)subname,' allocation error for dwt_deadcrootn_to_litter'; call endrun()
! Subprogram not used     end if
! Subprogram not used     allocate(conv_cflux(begp:endp), stat=ier)
! Subprogram not used     if (ier /= 0) then
! Subprogram not used           write(iulog,*)subname,' allocation error for conv_cflux'; call endrun()
! Subprogram not used     end if
! Subprogram not used     allocate(prod10_cflux(begp:endp), stat=ier)
! Subprogram not used     if (ier /= 0) then
! Subprogram not used           write(iulog,*)subname,' allocation error for prod10_cflux'; call endrun()
! Subprogram not used     end if
! Subprogram not used     allocate(prod100_cflux(begp:endp), stat=ier)
! Subprogram not used     if (ier /= 0) then
! Subprogram not used           write(iulog,*)subname,' allocation error for prod100_cflux'; call endrun()
! Subprogram not used     end if
! Subprogram not used     if (use_c13) then
! Subprogram not used        allocate(conv_c13flux(begp:endp), stat=ier)
! Subprogram not used        if (ier /= 0) then
! Subprogram not used           write(iulog,*)subname,' allocation error for conv_c13flux'; call endrun()
! Subprogram not used        end if
! Subprogram not used        allocate(prod10_c13flux(begp:endp), stat=ier)
! Subprogram not used        if (ier /= 0) then
! Subprogram not used           write(iulog,*)subname,' allocation error for prod10_c13flux'; call endrun()
! Subprogram not used        end if
! Subprogram not used        allocate(prod100_c13flux(begp:endp), stat=ier)
! Subprogram not used        if (ier /= 0) then
! Subprogram not used           write(iulog,*)subname,' allocation error for prod100_c13flux'; call endrun()
! Subprogram not used        end if
! Subprogram not used     endif
! Subprogram not used     allocate(conv_nflux(begp:endp), stat=ier)
! Subprogram not used     if (ier /= 0) then
! Subprogram not used           write(iulog,*)subname,' allocation error for conv_nflux'; call endrun()
! Subprogram not used     end if
! Subprogram not used     allocate(prod10_nflux(begp:endp), stat=ier)
! Subprogram not used     if (ier /= 0) then
! Subprogram not used           write(iulog,*)subname,' allocation error for prod10_nflux'; call endrun()
! Subprogram not used     end if
! Subprogram not used     allocate(prod100_nflux(begp:endp), stat=ier)
! Subprogram not used     if (ier /= 0) then
! Subprogram not used           write(iulog,*)subname,' allocation error for prod100_nflux'; call endrun()
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! Get time step
! Subprogram not used     dt = real( get_step_size(), r8 )
! Subprogram not used 
! Subprogram not used 	do p = begp,endp
! Subprogram not used                 c = pcolumn(p)
! Subprogram not used 		! initialize all the pft-level local flux arrays
! Subprogram not used 		dwt_leafc_seed(p) = 0._r8
! Subprogram not used 		dwt_leafn_seed(p) = 0._r8
! Subprogram not used                 if (use_c13) then
! Subprogram not used                    dwt_leafc13_seed(p) = 0._r8
! Subprogram not used                 endif
! Subprogram not used 		dwt_deadstemc_seed(p) = 0._r8
! Subprogram not used 		dwt_deadstemn_seed(p) = 0._r8
! Subprogram not used                 if (use_c13) then
! Subprogram not used                    dwt_deadstemc13_seed(p) = 0._r8
! Subprogram not used                 endif
! Subprogram not used 		dwt_frootc_to_litter(p) = 0._r8
! Subprogram not used 		dwt_livecrootc_to_litter(p) = 0._r8
! Subprogram not used 		dwt_deadcrootc_to_litter(p) = 0._r8
! Subprogram not used                 if (use_c13) then
! Subprogram not used                    dwt_frootc13_to_litter(p) = 0._r8
! Subprogram not used                    dwt_livecrootc13_to_litter(p) = 0._r8
! Subprogram not used                    dwt_deadcrootc13_to_litter(p) = 0._r8
! Subprogram not used                 endif
! Subprogram not used 		dwt_frootn_to_litter(p) = 0._r8
! Subprogram not used 		dwt_livecrootn_to_litter(p) = 0._r8
! Subprogram not used 		dwt_deadcrootn_to_litter(p) = 0._r8
! Subprogram not used 		conv_cflux(p) = 0._r8
! Subprogram not used 		prod10_cflux(p) = 0._r8
! Subprogram not used 		prod100_cflux(p) = 0._r8
! Subprogram not used                 if (use_c13) then
! Subprogram not used                    conv_c13flux(p) = 0._r8
! Subprogram not used                    prod10_c13flux(p) = 0._r8
! Subprogram not used                    prod100_c13flux(p) = 0._r8
! Subprogram not used                 endif
! Subprogram not used 		conv_nflux(p) = 0._r8
! Subprogram not used 		prod10_nflux(p) = 0._r8
! Subprogram not used 		prod100_nflux(p) = 0._r8
! Subprogram not used        
! Subprogram not used 		l = pptr%landunit(p)
! Subprogram not used 		if (lptr%itype(l) == istsoil .or. lptr%itype(l) == istcrop) then
! Subprogram not used 
! Subprogram not used 			! calculate the change in weight for the timestep
! Subprogram not used 			dwt = pptr%wtcol(p)-wtcol_old(p)
! Subprogram not used 
! Subprogram not used 			! PFTs for which weight increases on this timestep
! Subprogram not used 			if (dwt > 0._r8) then
! Subprogram not used 
! Subprogram not used 				! first identify PFTs that are initiating on this timestep
! Subprogram not used 				! and set all the necessary state and flux variables
! Subprogram not used 				if (wtcol_old(p) == 0._r8) then
! Subprogram not used 
! Subprogram not used 					! set initial conditions for PFT that is being initiated
! Subprogram not used 					! in this time step.  Based on the settings in cnIniTimeVar.
! Subprogram not used 
! Subprogram not used 					! pft-level carbon state variables
! Subprogram not used 					pcs%leafc(p)              = 0._r8
! Subprogram not used 					pcs%leafc_storage(p)      = 0._r8
! Subprogram not used 					pcs%leafc_xfer(p)         = 0._r8
! Subprogram not used 					pcs%frootc(p)             = 0._r8
! Subprogram not used 					pcs%frootc_storage(p)     = 0._r8
! Subprogram not used 					pcs%frootc_xfer(p)        = 0._r8
! Subprogram not used 					pcs%livestemc(p)          = 0._r8
! Subprogram not used 					pcs%livestemc_storage(p)  = 0._r8
! Subprogram not used 					pcs%livestemc_xfer(p)     = 0._r8
! Subprogram not used 					pcs%deadstemc(p)          = 0._r8
! Subprogram not used 					pcs%deadstemc_storage(p)  = 0._r8
! Subprogram not used 					pcs%deadstemc_xfer(p)     = 0._r8
! Subprogram not used 					pcs%livecrootc(p)         = 0._r8
! Subprogram not used 					pcs%livecrootc_storage(p) = 0._r8
! Subprogram not used 					pcs%livecrootc_xfer(p)    = 0._r8
! Subprogram not used 					pcs%deadcrootc(p)         = 0._r8
! Subprogram not used 					pcs%deadcrootc_storage(p) = 0._r8
! Subprogram not used 					pcs%deadcrootc_xfer(p)    = 0._r8
! Subprogram not used 					pcs%gresp_storage(p)      = 0._r8
! Subprogram not used 					pcs%gresp_xfer(p)         = 0._r8
! Subprogram not used 					pcs%cpool(p)              = 0._r8
! Subprogram not used 					pcs%xsmrpool(p)           = 0._r8
! Subprogram not used 					pcs%pft_ctrunc(p)         = 0._r8
! Subprogram not used 					pcs%dispvegc(p)           = 0._r8
! Subprogram not used 					pcs%storvegc(p)           = 0._r8
! Subprogram not used 					pcs%totvegc(p)            = 0._r8
! Subprogram not used 					pcs%totpftc(p)            = 0._r8
! Subprogram not used 
! Subprogram not used 					! pft-level carbon-13 state variables
! Subprogram not used                                         if (use_c13) then
! Subprogram not used                                            pc13s%leafc(p)              = 0._r8
! Subprogram not used                                            pc13s%leafc_storage(p)      = 0._r8
! Subprogram not used                                            pc13s%leafc_xfer(p)         = 0._r8
! Subprogram not used                                            pc13s%frootc(p)             = 0._r8
! Subprogram not used                                            pc13s%frootc_storage(p)     = 0._r8
! Subprogram not used                                            pc13s%frootc_xfer(p)        = 0._r8
! Subprogram not used                                            pc13s%livestemc(p)          = 0._r8
! Subprogram not used                                            pc13s%livestemc_storage(p)  = 0._r8
! Subprogram not used                                            pc13s%livestemc_xfer(p)     = 0._r8
! Subprogram not used                                            pc13s%deadstemc(p)          = 0._r8
! Subprogram not used                                            pc13s%deadstemc_storage(p)  = 0._r8
! Subprogram not used                                            pc13s%deadstemc_xfer(p)     = 0._r8
! Subprogram not used                                            pc13s%livecrootc(p)         = 0._r8
! Subprogram not used                                            pc13s%livecrootc_storage(p) = 0._r8
! Subprogram not used                                            pc13s%livecrootc_xfer(p)    = 0._r8
! Subprogram not used                                            pc13s%deadcrootc(p)         = 0._r8
! Subprogram not used                                            pc13s%deadcrootc_storage(p) = 0._r8
! Subprogram not used                                            pc13s%deadcrootc_xfer(p)    = 0._r8
! Subprogram not used                                            pc13s%gresp_storage(p)      = 0._r8
! Subprogram not used                                            pc13s%gresp_xfer(p)         = 0._r8
! Subprogram not used                                            pc13s%cpool(p)              = 0._r8
! Subprogram not used                                            pc13s%xsmrpool(p)           = 0._r8
! Subprogram not used                                            pc13s%pft_ctrunc(p)         = 0._r8
! Subprogram not used                                            pc13s%dispvegc(p)           = 0._r8
! Subprogram not used                                            pc13s%storvegc(p)           = 0._r8
! Subprogram not used                                            pc13s%totvegc(p)            = 0._r8
! Subprogram not used                                            pc13s%totpftc(p)            = 0._r8
! Subprogram not used                                         endif
! Subprogram not used 
! Subprogram not used 					! pft-level nitrogen state variables
! Subprogram not used 					pns%leafn(p)	           = 0._r8
! Subprogram not used 					pns%leafn_storage(p)      = 0._r8
! Subprogram not used 					pns%leafn_xfer(p)         = 0._r8
! Subprogram not used 					pns%frootn(p)	           = 0._r8
! Subprogram not used 					pns%frootn_storage(p)     = 0._r8
! Subprogram not used 					pns%frootn_xfer(p)        = 0._r8
! Subprogram not used 					pns%livestemn(p)	       = 0._r8
! Subprogram not used 					pns%livestemn_storage(p)  = 0._r8
! Subprogram not used 					pns%livestemn_xfer(p)     = 0._r8
! Subprogram not used 					pns%deadstemn(p)	       = 0._r8
! Subprogram not used 					pns%deadstemn_storage(p)  = 0._r8
! Subprogram not used 					pns%deadstemn_xfer(p)     = 0._r8
! Subprogram not used 					pns%livecrootn(p)         = 0._r8
! Subprogram not used 					pns%livecrootn_storage(p) = 0._r8
! Subprogram not used 					pns%livecrootn_xfer(p)    = 0._r8
! Subprogram not used 					pns%deadcrootn(p)         = 0._r8
! Subprogram not used 					pns%deadcrootn_storage(p) = 0._r8
! Subprogram not used 					pns%deadcrootn_xfer(p)    = 0._r8
! Subprogram not used 					pns%retransn(p)	       = 0._r8
! Subprogram not used 					pns%npool(p)	           = 0._r8
! Subprogram not used 					pns%pft_ntrunc(p)         = 0._r8
! Subprogram not used 					pns%dispvegn(p)           = 0._r8
! Subprogram not used 					pns%storvegn(p)           = 0._r8
! Subprogram not used 					pns%totvegn(p)            = 0._r8
! Subprogram not used 					pns%totpftn (p)           = 0._r8
! Subprogram not used 
! Subprogram not used 					! initialize same flux and epv variables that are set
! Subprogram not used 					! in CNiniTimeVar
! Subprogram not used 					pcf%psnsun(p) = 0._r8
! Subprogram not used 					pcf%psnsha(p) = 0._r8
! Subprogram not used                                         if (use_c13) then
! Subprogram not used                                            pc13f%psnsun(p) = 0._r8
! Subprogram not used                                            pc13f%psnsha(p) = 0._r8
! Subprogram not used                                         endif
! Subprogram not used 					pps%laisun(p) = 0._r8
! Subprogram not used 					pps%laisha(p) = 0._r8
! Subprogram not used 					pps%lncsun(p) = 0._r8
! Subprogram not used 					pps%lncsha(p) = 0._r8
! Subprogram not used 					pps%vcmxsun(p) = 0._r8
! Subprogram not used 					pps%vcmxsha(p) = 0._r8
! Subprogram not used                                         if (use_c13) then
! Subprogram not used                                            pps%alphapsnsun(p) = 0._r8
! Subprogram not used                                            pps%alphapsnsha(p) = 0._r8
! Subprogram not used                                         endif
! Subprogram not used 
! Subprogram not used 					pepv%dormant_flag(p) = 1._r8
! Subprogram not used 					pepv%days_active(p) = 0._r8
! Subprogram not used 					pepv%onset_flag(p) = 0._r8
! Subprogram not used 					pepv%onset_counter(p) = 0._r8
! Subprogram not used 					pepv%onset_gddflag(p) = 0._r8
! Subprogram not used 					pepv%onset_fdd(p) = 0._r8
! Subprogram not used 					pepv%onset_gdd(p) = 0._r8
! Subprogram not used 					pepv%onset_swi(p) = 0.0_r8
! Subprogram not used 					pepv%offset_flag(p) = 0._r8
! Subprogram not used 					pepv%offset_counter(p) = 0._r8
! Subprogram not used 					pepv%offset_fdd(p) = 0._r8
! Subprogram not used 					pepv%offset_swi(p) = 0._r8
! Subprogram not used 					pepv%lgsf(p) = 0._r8
! Subprogram not used 					pepv%bglfr(p) = 0._r8
! Subprogram not used 					pepv%bgtr(p) = 0._r8
! Subprogram not used 					! difference from CNiniTimeVar: using column-level
! Subprogram not used 					! information to initialize annavg_t2m.
! Subprogram not used 					pepv%annavg_t2m(p) = cps%cannavg_t2m(c)
! Subprogram not used 					pepv%tempavg_t2m(p) = 0._r8
! Subprogram not used 					pepv%gpp(p) = 0._r8
! Subprogram not used 					pepv%availc(p) = 0._r8
! Subprogram not used 					pepv%xsmrpool_recover(p) = 0._r8
! Subprogram not used                                         if (use_c13) then
! Subprogram not used                                            pepv%xsmrpool_c13ratio(p) = c13ratio
! Subprogram not used                                         endif
! Subprogram not used 					pepv%alloc_pnow(p) = 1._r8
! Subprogram not used 					pepv%c_allometry(p) = 0._r8
! Subprogram not used 					pepv%n_allometry(p) = 0._r8
! Subprogram not used 					pepv%plant_ndemand(p) = 0._r8
! Subprogram not used 					pepv%tempsum_potential_gpp(p) = 0._r8
! Subprogram not used 					pepv%annsum_potential_gpp(p) = 0._r8
! Subprogram not used 					pepv%tempmax_retransn(p) = 0._r8
! Subprogram not used 					pepv%annmax_retransn(p) = 0._r8
! Subprogram not used 					pepv%avail_retransn(p) = 0._r8
! Subprogram not used 					pepv%plant_nalloc(p) = 0._r8
! Subprogram not used 					pepv%plant_calloc(p) = 0._r8
! Subprogram not used 					pepv%excess_cflux(p) = 0._r8
! Subprogram not used 					pepv%downreg(p) = 0._r8
! Subprogram not used 					pepv%prev_leafc_to_litter(p) = 0._r8
! Subprogram not used 					pepv%prev_frootc_to_litter(p) = 0._r8
! Subprogram not used 					pepv%tempsum_npp(p) = 0._r8
! Subprogram not used 					pepv%annsum_npp(p) = 0._r8
! Subprogram not used                                         if (use_c13) then
! Subprogram not used                                            pepv%rc13_canair(p) = 0._r8
! Subprogram not used                                            pepv%rc13_psnsun(p) = 0._r8
! Subprogram not used                                            pepv%rc13_psnsha(p) = 0._r8
! Subprogram not used                                         endif
! Subprogram not used 
! Subprogram not used 				end if  ! end initialization of new pft
! Subprogram not used 
! Subprogram not used 				! (still in dwt > 0 block)
! Subprogram not used 
! Subprogram not used 				! set the seed sources for leaf and deadstem
! Subprogram not used 				! leaf source is split later between leaf, leaf_storage, leaf_xfer
! Subprogram not used 				leafc_seed   = 0._r8
! Subprogram not used 				leafn_seed   = 0._r8
! Subprogram not used                                 if (use_c13) then
! Subprogram not used                                    leafc13_seed = 0._r8
! Subprogram not used                                 endif
! Subprogram not used 				deadstemc_seed   = 0._r8
! Subprogram not used 				deadstemn_seed   = 0._r8
! Subprogram not used                                 if (use_c13) then
! Subprogram not used                                    deadstemc13_seed = 0._r8
! Subprogram not used                                 endif
! Subprogram not used 				if (pptr%itype(p) /= 0) then
! Subprogram not used 					leafc_seed = 1._r8
! Subprogram not used 					leafn_seed  = leafc_seed / pftcon%leafcn(pptr%itype(p))
! Subprogram not used 					if (pftcon%woody(pptr%itype(p)) == 1._r8) then
! Subprogram not used 						deadstemc_seed = 0.1_r8
! Subprogram not used 						deadstemn_seed = deadstemc_seed / pftcon%deadwdcn(pptr%itype(p))
! Subprogram not used 					end if
! Subprogram not used 
! Subprogram not used                                         if (use_c13) then
! Subprogram not used                                            ! 13c state is initialized assuming del13c = -28 permil for C3, and -13 permil for C4.
! Subprogram not used                                            ! That translates to ratios of (13c/(12c+13c)) of 0.01080455 for C3, and 0.01096945 for C4
! Subprogram not used                                            ! based on the following formulae: 
! Subprogram not used                                            ! r1 (13/12) = PDB + (del13c * PDB)/1000.0
! Subprogram not used                                            ! r2 (13/(13+12)) = r1/(1+r1)
! Subprogram not used                                            ! PDB = 0.0112372_R8  (ratio of 13C/12C in Pee Dee Belemnite, C isotope standard)
! Subprogram not used                                            c3_del13c = -28._r8
! Subprogram not used                                            c4_del13c = -13._r8
! Subprogram not used                                            c3_r1 = SHR_CONST_PDB + ((c3_del13c*SHR_CONST_PDB)/1000._r8)
! Subprogram not used                                            c3_r2 = c3_r1/(1._r8 + c3_r1)
! Subprogram not used                                            c4_r1 = SHR_CONST_PDB + ((c4_del13c*SHR_CONST_PDB)/1000._r8)
! Subprogram not used                                            c4_r2 = c4_r1/(1._r8 + c4_r1)
! Subprogram not used                                            
! Subprogram not used                                            if (pftcon%c3psn(pptr%itype(p)) == 1._r8) then
! Subprogram not used                                               leafc13_seed     = leafc_seed     * c3_r2
! Subprogram not used                                               deadstemc13_seed = deadstemc_seed * c3_r2
! Subprogram not used                                            else
! Subprogram not used                                               leafc13_seed     = leafc_seed     * c4_r2
! Subprogram not used                                               deadstemc13_seed = deadstemc_seed * c4_r2
! Subprogram not used                                            end if
! Subprogram not used                                         endif
! Subprogram not used                                      end if
! Subprogram not used 
! Subprogram not used 				! When PFT area expands (dwt > 0), the pft-level mass density 
! Subprogram not used 				! is modified to conserve the original pft mass distributed
! Subprogram not used 				! over the new (larger) area, plus a term to account for the 
! Subprogram not used 				! introduction of new seed source for leaf and deadstem
! Subprogram not used 				t1 = wtcol_old(p)/pptr%wtcol(p)
! Subprogram not used 				t2 = dwt/pptr%wtcol(p)
! Subprogram not used 
! Subprogram not used 				tot_leaf = pcs%leafc(p) + pcs%leafc_storage(p) + pcs%leafc_xfer(p)
! Subprogram not used 				pleaf = 0._r8
! Subprogram not used 				pstor = 0._r8
! Subprogram not used 				pxfer = 0._r8
! Subprogram not used 				if (tot_leaf /= 0._r8) then
! Subprogram not used 					! when adding seed source to non-zero leaf state, use current proportions
! Subprogram not used 					pleaf = pcs%leafc(p)/tot_leaf
! Subprogram not used 					pstor = pcs%leafc_storage(p)/tot_leaf
! Subprogram not used 					pxfer = pcs%leafc_xfer(p)/tot_leaf
! Subprogram not used 				else
! Subprogram not used 					! when initiating from zero leaf state, use evergreen flag to set proportions
! Subprogram not used 					if (pftcon%evergreen(pptr%itype(p)) == 1._r8) then
! Subprogram not used 						pleaf = 1._r8
! Subprogram not used 					else
! Subprogram not used 						pstor = 1._r8
! Subprogram not used 					end if
! Subprogram not used 				end if 
! Subprogram not used 				pcs%leafc(p)         = pcs%leafc(p)*t1         + leafc_seed*pleaf*t2
! Subprogram not used 				pcs%leafc_storage(p) = pcs%leafc_storage(p)*t1 + leafc_seed*pstor*t2
! Subprogram not used 				pcs%leafc_xfer(p)    = pcs%leafc_xfer(p)*t1    + leafc_seed*pxfer*t2
! Subprogram not used 				pcs%frootc(p)  		   = pcs%frootc(p) 			* t1
! Subprogram not used 				pcs%frootc_storage(p)     = pcs%frootc_storage(p) 	* t1
! Subprogram not used 				pcs%frootc_xfer(p) 	   = pcs%frootc_xfer(p)		* t1
! Subprogram not used 				pcs%livestemc(p)		   = pcs%livestemc(p)  		* t1
! Subprogram not used 				pcs%livestemc_storage(p)  = pcs%livestemc_storage(p)  * t1
! Subprogram not used 				pcs%livestemc_xfer(p)     = pcs%livestemc_xfer(p) 	* t1
! Subprogram not used 				pcs%deadstemc(p)     = pcs%deadstemc(p)*t1     + deadstemc_seed*t2
! Subprogram not used 				pcs%deadstemc_storage(p)  = pcs%deadstemc_storage(p)  * t1
! Subprogram not used 				pcs%deadstemc_xfer(p)     = pcs%deadstemc_xfer(p) 	* t1
! Subprogram not used 				pcs%livecrootc(p)  	   = pcs%livecrootc(p) 		* t1
! Subprogram not used 				pcs%livecrootc_storage(p) = pcs%livecrootc_storage(p) * t1
! Subprogram not used 				pcs%livecrootc_xfer(p)    = pcs%livecrootc_xfer(p)	* t1
! Subprogram not used 				pcs%deadcrootc(p)  	   = pcs%deadcrootc(p) 		* t1
! Subprogram not used 				pcs%deadcrootc_storage(p) = pcs%deadcrootc_storage(p) * t1
! Subprogram not used 				pcs%deadcrootc_xfer(p)    = pcs%deadcrootc_xfer(p)	* t1
! Subprogram not used 				pcs%gresp_storage(p)	   = pcs%gresp_storage(p)  	* t1
! Subprogram not used 				pcs%gresp_xfer(p)  	   = pcs%gresp_xfer(p) 		* t1
! Subprogram not used 				pcs%cpool(p)			   = pcs%cpool(p)  			* t1
! Subprogram not used 				pcs%xsmrpool(p)		   = pcs%xsmrpool(p)			* t1
! Subprogram not used 				pcs%pft_ctrunc(p)  	   = pcs%pft_ctrunc(p) 		* t1
! Subprogram not used 				pcs%dispvegc(p)		   = pcs%dispvegc(p)			* t1
! Subprogram not used 				pcs%storvegc(p)		   = pcs%storvegc(p)			* t1
! Subprogram not used 				pcs%totvegc(p) 		   = pcs%totvegc(p)			* t1
! Subprogram not used 				pcs%totpftc(p) 		   = pcs%totpftc(p)			* t1
! Subprogram not used 
! Subprogram not used 				! pft-level carbon-13 state variables 
! Subprogram not used                                 if (use_c13) then
! Subprogram not used                                    tot_leaf = pc13s%leafc(p) + pc13s%leafc_storage(p) + pc13s%leafc_xfer(p)
! Subprogram not used                                    pleaf = 0._r8
! Subprogram not used                                    pstor = 0._r8
! Subprogram not used                                    pxfer = 0._r8
! Subprogram not used                                    if (tot_leaf /= 0._r8) then
! Subprogram not used                                       pleaf = pc13s%leafc(p)/tot_leaf
! Subprogram not used                                       pstor = pc13s%leafc_storage(p)/tot_leaf
! Subprogram not used                                       pxfer = pc13s%leafc_xfer(p)/tot_leaf
! Subprogram not used                                    else
! Subprogram not used                                       ! when initiating from zero leaf state, use evergreen flag to set proportions
! Subprogram not used                                       if (pftcon%evergreen(pptr%itype(p)) == 1._r8) then
! Subprogram not used                                          pleaf = 1._r8
! Subprogram not used                                       else
! Subprogram not used                                          pstor = 1._r8
! Subprogram not used                                       end if
! Subprogram not used                                    end if
! Subprogram not used                                    pc13s%leafc(p)         = pc13s%leafc(p)*t1         + leafc13_seed*pleaf*t2
! Subprogram not used                                    pc13s%leafc_storage(p) = pc13s%leafc_storage(p)*t1 + leafc13_seed*pstor*t2
! Subprogram not used                                    pc13s%leafc_xfer(p)    = pc13s%leafc_xfer(p)*t1    + leafc13_seed*pxfer*t2
! Subprogram not used                                    pc13s%frootc(p)			 = pc13s%frootc(p) 		* t1
! Subprogram not used                                    pc13s%frootc_storage(p)	         = pc13s%frootc_storage(p) 	* t1
! Subprogram not used                                    pc13s%frootc_xfer(p)		 = pc13s%frootc_xfer(p)		* t1
! Subprogram not used                                    pc13s%livestemc(p) 		 = pc13s%livestemc(p)  		* t1
! Subprogram not used                                    pc13s%livestemc_storage(p)       = pc13s%livestemc_storage(p)      * t1
! Subprogram not used                                    pc13s%livestemc_xfer(p)	         = pc13s%livestemc_xfer(p) 	* t1
! Subprogram not used                                    pc13s%deadstemc(p)               = pc13s%deadstemc(p)*t1     + deadstemc13_seed*t2
! Subprogram not used                                    pc13s%deadstemc_storage(p)       = pc13s%deadstemc_storage(p)      * t1
! Subprogram not used                                    pc13s%deadstemc_xfer(p)	         = pc13s%deadstemc_xfer(p) 	* t1
! Subprogram not used                                    pc13s%livecrootc(p)		 = pc13s%livecrootc(p) 		* t1
! Subprogram not used                                    pc13s%livecrootc_storage(p)      = pc13s%livecrootc_storage(p)     * t1
! Subprogram not used                                    pc13s%livecrootc_xfer(p)         = pc13s%livecrootc_xfer(p)	* t1
! Subprogram not used                                    pc13s%deadcrootc(p)		 = pc13s%deadcrootc(p) 		* t1
! Subprogram not used                                    pc13s%deadcrootc_storage(p)      = pc13s%deadcrootc_storage(p)     * t1
! Subprogram not used                                    pc13s%deadcrootc_xfer(p)         = pc13s%deadcrootc_xfer(p)	* t1
! Subprogram not used                                    pc13s%gresp_storage(p) 	         = pc13s%gresp_storage(p)  	* t1
! Subprogram not used                                    pc13s%gresp_xfer(p)		 = pc13s%gresp_xfer(p) 		* t1
! Subprogram not used                                    pc13s%cpool(p) 			 = pc13s%cpool(p)  		* t1
! Subprogram not used                                    pc13s%xsmrpool(p)  		 = pc13s%xsmrpool(p)		* t1
! Subprogram not used                                    pc13s%pft_ctrunc(p)		 = pc13s%pft_ctrunc(p) 		* t1
! Subprogram not used                                    pc13s%dispvegc(p)  		 = pc13s%dispvegc(p)		* t1
! Subprogram not used                                    pc13s%storvegc(p)  		 = pc13s%storvegc(p)		* t1
! Subprogram not used                                    pc13s%totvegc(p)		 = pc13s%totvegc(p)		* t1
! Subprogram not used                                    pc13s%totpftc(p)		 = pc13s%totpftc(p)		* t1
! Subprogram not used                                 endif
! Subprogram not used 
! Subprogram not used 				tot_leaf = pns%leafn(p) + pns%leafn_storage(p) + pns%leafn_xfer(p)
! Subprogram not used 				pleaf = 0._r8
! Subprogram not used 				pstor = 0._r8
! Subprogram not used 				pxfer = 0._r8
! Subprogram not used 				if (tot_leaf /= 0._r8) then
! Subprogram not used 					pleaf = pns%leafn(p)/tot_leaf
! Subprogram not used 					pstor = pns%leafn_storage(p)/tot_leaf
! Subprogram not used 					pxfer = pns%leafn_xfer(p)/tot_leaf
! Subprogram not used 				else
! Subprogram not used 					! when initiating from zero leaf state, use evergreen flag to set proportions
! Subprogram not used 					if (pftcon%evergreen(pptr%itype(p)) == 1._r8) then
! Subprogram not used 						pleaf = 1._r8
! Subprogram not used 					else
! Subprogram not used 						pstor = 1._r8
! Subprogram not used 					end if
! Subprogram not used 				end if 
! Subprogram not used 				! pft-level nitrogen state variables
! Subprogram not used 				pns%leafn(p)         = pns%leafn(p)*t1         + leafn_seed*pleaf*t2
! Subprogram not used 				pns%leafn_storage(p) = pns%leafn_storage(p)*t1 + leafn_seed*pstor*t2
! Subprogram not used 				pns%leafn_xfer(p)    = pns%leafn_xfer(p)*t1    + leafn_seed*pxfer*t2
! Subprogram not used 				pns%frootn(p)  		   = pns%frootn(p) 		* t1
! Subprogram not used 				pns%frootn_storage(p)         = pns%frootn_storage(p) 	* t1
! Subprogram not used 				pns%frootn_xfer(p) 	   = pns%frootn_xfer(p)		* t1
! Subprogram not used 				pns%livestemn(p)		   = pns%livestemn(p)  		* t1
! Subprogram not used 				pns%livestemn_storage(p)      = pns%livestemn_storage(p)      * t1
! Subprogram not used 				pns%livestemn_xfer(p)         = pns%livestemn_xfer(p) 	* t1
! Subprogram not used 				pns%deadstemn(p)              = pns%deadstemn(p)*t1     + deadstemn_seed*t2
! Subprogram not used 				pns%deadstemn_storage(p)      = pns%deadstemn_storage(p)      * t1
! Subprogram not used 				pns%deadstemn_xfer(p)         = pns%deadstemn_xfer(p) 	* t1
! Subprogram not used 				pns%livecrootn(p)  	   = pns%livecrootn(p) 		* t1
! Subprogram not used 				pns%livecrootn_storage(p)     = pns%livecrootn_storage(p)     * t1
! Subprogram not used 				pns%livecrootn_xfer(p)        = pns%livecrootn_xfer(p)	* t1
! Subprogram not used 				pns%deadcrootn(p)  	   = pns%deadcrootn(p) 		* t1
! Subprogram not used 				pns%deadcrootn_storage(p)     = pns%deadcrootn_storage(p)     * t1
! Subprogram not used 				pns%deadcrootn_xfer(p)        = pns%deadcrootn_xfer(p)        * t1
! Subprogram not used 				pns%retransn(p)		   = pns%retransn(p)		* t1
! Subprogram not used 				pns%npool(p)		   = pns%npool(p)  		* t1
! Subprogram not used 				pns%pft_ntrunc(p)  	   = pns%pft_ntrunc(p)        	* t1
! Subprogram not used 				pns%dispvegn(p)		   = pns%dispvegn(p)		* t1
! Subprogram not used 				pns%storvegn(p)		   = pns%storvegn(p)		* t1
! Subprogram not used 				pns%totvegn(p) 		   = pns%totvegn(p)		* t1
! Subprogram not used 				pns%totpftn(p) 		   = pns%totpftn(p)		* t1
! Subprogram not used 
! Subprogram not used 				! update temporary seed source arrays
! Subprogram not used 				! These are calculated in terms of the required contributions from
! Subprogram not used 				! column-level seed source
! Subprogram not used 				dwt_leafc_seed(p)   = leafc_seed   * dwt
! Subprogram not used                                 if (use_c13) then
! Subprogram not used                                    dwt_leafc13_seed(p) = leafc13_seed * dwt
! Subprogram not used                                 endif
! Subprogram not used 				dwt_leafn_seed(p)   = leafn_seed   * dwt
! Subprogram not used 				dwt_deadstemc_seed(p)   = deadstemc_seed   * dwt
! Subprogram not used                                 if (use_c13) then
! Subprogram not used                                    dwt_deadstemc13_seed(p) = deadstemc13_seed * dwt
! Subprogram not used                                 endif
! Subprogram not used 				dwt_deadstemn_seed(p)   = deadstemn_seed   * dwt
! Subprogram not used 
! Subprogram not used 			else if (dwt < 0._r8) then
! Subprogram not used 
! Subprogram not used 				! if the pft lost weight on the timestep, then the carbon and nitrogen state
! Subprogram not used 				! variables are directed to litter, CWD, and wood product pools.
! Subprogram not used 
! Subprogram not used 				! N.B. : the conv_cflux, prod10_cflux, and prod100_cflux fluxes are accumulated
! Subprogram not used 				! as negative values, but the fluxes for pft-to-litter are accumulated as 
! Subprogram not used 				! positive values
! Subprogram not used 
! Subprogram not used 				! set local weight variables for this pft
! Subprogram not used 				wt_new = pptr%wtcol(p)
! Subprogram not used 				wt_old = wtcol_old(p)
! Subprogram not used 
! Subprogram not used 				!---------------
! Subprogram not used 				! C state update
! Subprogram not used 				!---------------
! Subprogram not used 
! Subprogram not used 				! leafc 
! Subprogram not used 				ptr => pcs%leafc(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) - init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! leafc_storage 
! Subprogram not used 				ptr => pcs%leafc_storage(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) - init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! leafc_xfer 
! Subprogram not used 				ptr => pcs%leafc_xfer(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) - init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! frootc 
! Subprogram not used 				ptr => pcs%frootc(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					dwt_frootc_to_litter(p) = dwt_frootc_to_litter(p) - change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					dwt_frootc_to_litter(p) = dwt_frootc_to_litter(p) + init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! frootc_storage 
! Subprogram not used 				ptr => pcs%frootc_storage(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) - init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! frootc_xfer 
! Subprogram not used 				ptr => pcs%frootc_xfer(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) - init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! livestemc 
! Subprogram not used 				ptr => pcs%livestemc(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) - init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! livestemc_storage 
! Subprogram not used 				ptr => pcs%livestemc_storage(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) - init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! livestemc_xfer 
! Subprogram not used 				ptr => pcs%livestemc_xfer(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) - init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! deadstemc 
! Subprogram not used 				ptr => pcs%deadstemc(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) + change_state*pconv(pptr%itype(p))
! Subprogram not used 					prod10_cflux(p) = prod10_cflux(p) + change_state*pprod10(pptr%itype(p))
! Subprogram not used 					prod100_cflux(p) = prod100_cflux(p) + change_state*pprod100(pptr%itype(p))
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) - init_state*pconv(pptr%itype(p))
! Subprogram not used 					prod10_cflux(p) = prod10_cflux(p) - init_state*pprod10(pptr%itype(p))
! Subprogram not used 					prod100_cflux(p) = prod100_cflux(p) - init_state*pprod100(pptr%itype(p))
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! deadstemc_storage 
! Subprogram not used 				ptr => pcs%deadstemc_storage(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) - init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! deadstemc_xfer 
! Subprogram not used 				ptr => pcs%deadstemc_xfer(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) - init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! livecrootc 
! Subprogram not used 				ptr => pcs%livecrootc(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					dwt_livecrootc_to_litter(p) = dwt_livecrootc_to_litter(p) - change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					dwt_livecrootc_to_litter(p) = dwt_livecrootc_to_litter(p) + init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! livecrootc_storage 
! Subprogram not used 				ptr => pcs%livecrootc_storage(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) - init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! livecrootc_xfer 
! Subprogram not used 				ptr => pcs%livecrootc_xfer(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) - init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! deadcrootc 
! Subprogram not used 				ptr => pcs%deadcrootc(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					dwt_deadcrootc_to_litter(p) = dwt_deadcrootc_to_litter(p) - change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					dwt_deadcrootc_to_litter(p) = dwt_deadcrootc_to_litter(p) + init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! deadcrootc_storage 
! Subprogram not used 				ptr => pcs%deadcrootc_storage(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) - init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! deadcrootc_xfer 
! Subprogram not used 				ptr => pcs%deadcrootc_xfer(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) - init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! gresp_storage 
! Subprogram not used 				ptr => pcs%gresp_storage(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) - init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! gresp_xfer 
! Subprogram not used 				ptr => pcs%gresp_xfer(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) - init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! cpool 
! Subprogram not used 				ptr => pcs%cpool(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) - init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! xsmrpool 
! Subprogram not used 				ptr => pcs%xsmrpool(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) - init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! pft_ctrunc 
! Subprogram not used 				ptr => pcs%pft_ctrunc(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					conv_cflux(p) = conv_cflux(p) - init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used                                 if (use_c13) then
! Subprogram not used                                    !-----------------
! Subprogram not used                                    ! C13 state update
! Subprogram not used                                    !-----------------
! Subprogram not used                                    
! Subprogram not used                                    ! set pointers to the conversion and product pool fluxes for this pft
! Subprogram not used                                    ! dwt_ptr0 is reserved for local assignment to dwt_xxx_to_litter fluxes
! Subprogram not used                                    dwt_ptr1 => conv_c13flux(p)
! Subprogram not used                                    dwt_ptr2 => prod10_c13flux(p)
! Subprogram not used                                    dwt_ptr3 => prod100_c13flux(p)
! Subprogram not used                                    
! Subprogram not used                                    ! leafc 
! Subprogram not used                                    ptr => pc13s%leafc(p)
! Subprogram not used                                    init_state = ptr*wt_old
! Subprogram not used                                    change_state = ptr*dwt
! Subprogram not used                                    new_state = init_state+change_state
! Subprogram not used                                    if (wt_new /= 0._r8) then
! Subprogram not used                                       ptr = new_state/wt_new
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 + change_state
! Subprogram not used                                    else
! Subprogram not used                                       ptr = 0._r8
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 - init_state
! Subprogram not used                                    end if
! Subprogram not used                                    
! Subprogram not used                                    ! leafc_storage 
! Subprogram not used                                    ptr => pc13s%leafc_storage(p)
! Subprogram not used                                    init_state = ptr*wt_old
! Subprogram not used                                    change_state = ptr*dwt
! Subprogram not used                                    new_state = init_state+change_state
! Subprogram not used                                    if (wt_new /= 0._r8) then
! Subprogram not used                                       ptr = new_state/wt_new
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 + change_state
! Subprogram not used                                    else
! Subprogram not used                                       ptr = 0._r8
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 - init_state
! Subprogram not used                                    end if
! Subprogram not used                                    
! Subprogram not used                                    ! leafc_xfer 
! Subprogram not used                                    ptr => pc13s%leafc_xfer(p)
! Subprogram not used                                    init_state = ptr*wt_old
! Subprogram not used                                    change_state = ptr*dwt
! Subprogram not used                                    new_state = init_state+change_state
! Subprogram not used                                    if (wt_new /= 0._r8) then
! Subprogram not used                                       ptr = new_state/wt_new
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 + change_state
! Subprogram not used                                    else
! Subprogram not used                                       ptr = 0._r8
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 - init_state
! Subprogram not used                                    end if
! Subprogram not used                                    
! Subprogram not used                                    ! frootc 
! Subprogram not used                                    ptr => pc13s%frootc(p)
! Subprogram not used                                    dwt_ptr0 => dwt_frootc13_to_litter(p)
! Subprogram not used                                    init_state = ptr*wt_old
! Subprogram not used                                    change_state = ptr*dwt
! Subprogram not used                                    new_state = init_state+change_state
! Subprogram not used                                    if (wt_new /= 0._r8) then
! Subprogram not used                                       ptr = new_state/wt_new
! Subprogram not used                                       dwt_ptr0 = dwt_ptr0 - change_state
! Subprogram not used                                    else
! Subprogram not used                                       ptr = 0._r8
! Subprogram not used                                       dwt_ptr0 = dwt_ptr0 + init_state
! Subprogram not used                                    end if
! Subprogram not used                                    
! Subprogram not used                                    ! frootc_storage 
! Subprogram not used                                    ptr => pc13s%frootc_storage(p)
! Subprogram not used                                    init_state = ptr*wt_old
! Subprogram not used                                    change_state = ptr*dwt
! Subprogram not used                                    new_state = init_state+change_state
! Subprogram not used                                    if (wt_new /= 0._r8) then
! Subprogram not used                                       ptr = new_state/wt_new
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 + change_state
! Subprogram not used                                    else
! Subprogram not used                                       ptr = 0._r8
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 - init_state
! Subprogram not used                                    end if
! Subprogram not used 
! Subprogram not used                                    ! frootc_xfer 
! Subprogram not used                                    ptr => pc13s%frootc_xfer(p)
! Subprogram not used                                    init_state = ptr*wt_old
! Subprogram not used                                    change_state = ptr*dwt
! Subprogram not used                                    new_state = init_state+change_state
! Subprogram not used                                    if (wt_new /= 0._r8) then
! Subprogram not used                                       ptr = new_state/wt_new
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 + change_state
! Subprogram not used                                    else
! Subprogram not used                                       ptr = 0._r8
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 - init_state
! Subprogram not used                                    end if
! Subprogram not used                                    
! Subprogram not used                                    ! livestemc 
! Subprogram not used                                    ptr => pc13s%livestemc(p)
! Subprogram not used                                    init_state = ptr*wt_old
! Subprogram not used                                    change_state = ptr*dwt
! Subprogram not used                                    new_state = init_state+change_state
! Subprogram not used                                    if (wt_new /= 0._r8) then
! Subprogram not used                                       ptr = new_state/wt_new
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 + change_state
! Subprogram not used                                    else
! Subprogram not used                                       ptr = 0._r8
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 - init_state
! Subprogram not used                                    end if
! Subprogram not used                                    
! Subprogram not used                                    ! livestemc_storage 
! Subprogram not used                                    ptr => pc13s%livestemc_storage(p)
! Subprogram not used                                    init_state = ptr*wt_old
! Subprogram not used                                    change_state = ptr*dwt
! Subprogram not used                                    new_state = init_state+change_state
! Subprogram not used                                    if (wt_new /= 0._r8) then
! Subprogram not used                                       ptr = new_state/wt_new
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 + change_state
! Subprogram not used                                    else
! Subprogram not used                                       ptr = 0._r8
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 - init_state
! Subprogram not used                                    end if
! Subprogram not used                                    
! Subprogram not used                                    ! livestemc_xfer 
! Subprogram not used                                    ptr => pc13s%livestemc_xfer(p)
! Subprogram not used                                    init_state = ptr*wt_old
! Subprogram not used                                    change_state = ptr*dwt
! Subprogram not used                                    new_state = init_state+change_state
! Subprogram not used                                    if (wt_new /= 0._r8) then
! Subprogram not used                                       ptr = new_state/wt_new
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 + change_state
! Subprogram not used                                    else
! Subprogram not used                                       ptr = 0._r8
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 - init_state
! Subprogram not used                                    end if
! Subprogram not used 
! Subprogram not used                                    ! deadstemc 
! Subprogram not used                                    ptr => pc13s%deadstemc(p)
! Subprogram not used                                    init_state = ptr*wt_old
! Subprogram not used                                    change_state = ptr*dwt
! Subprogram not used                                    new_state = init_state+change_state
! Subprogram not used                                    if (wt_new /= 0._r8) then
! Subprogram not used                                       ptr = new_state/wt_new
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 + change_state*pconv(pptr%itype(p))
! Subprogram not used                                       dwt_ptr2 = dwt_ptr2 + change_state*pprod10(pptr%itype(p))
! Subprogram not used                                       dwt_ptr3 = dwt_ptr3 + change_state*pprod100(pptr%itype(p))
! Subprogram not used                                    else
! Subprogram not used                                       ptr = 0._r8
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 - init_state*pconv(pptr%itype(p))
! Subprogram not used                                       dwt_ptr2 = dwt_ptr2 - init_state*pprod10(pptr%itype(p))
! Subprogram not used                                       dwt_ptr3 = dwt_ptr3 - init_state*pprod100(pptr%itype(p))
! Subprogram not used                                    end if
! Subprogram not used                                    
! Subprogram not used                                    ! deadstemc_storage 
! Subprogram not used                                    ptr => pc13s%deadstemc_storage(p)
! Subprogram not used                                    init_state = ptr*wt_old
! Subprogram not used                                    change_state = ptr*dwt
! Subprogram not used                                    new_state = init_state+change_state
! Subprogram not used                                    if (wt_new /= 0._r8) then
! Subprogram not used                                       ptr = new_state/wt_new
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 + change_state
! Subprogram not used                                    else
! Subprogram not used                                       ptr = 0._r8
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 - init_state
! Subprogram not used                                    end if
! Subprogram not used                                    
! Subprogram not used                                    ! deadstemc_xfer 
! Subprogram not used                                    ptr => pc13s%deadstemc_xfer(p)
! Subprogram not used                                    init_state = ptr*wt_old
! Subprogram not used                                    change_state = ptr*dwt
! Subprogram not used                                    new_state = init_state+change_state
! Subprogram not used                                    if (wt_new /= 0._r8) then
! Subprogram not used                                       ptr = new_state/wt_new
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 + change_state
! Subprogram not used                                    else
! Subprogram not used                                       ptr = 0._r8
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 - init_state
! Subprogram not used                                    end if
! Subprogram not used                                    
! Subprogram not used                                    ! livecrootc 
! Subprogram not used                                    ptr => pc13s%livecrootc(p)
! Subprogram not used                                    dwt_ptr0 => dwt_livecrootc13_to_litter(p)
! Subprogram not used                                    init_state = ptr*wt_old
! Subprogram not used                                    change_state = ptr*dwt
! Subprogram not used                                    new_state = init_state+change_state
! Subprogram not used                                    if (wt_new /= 0._r8) then
! Subprogram not used                                       ptr = new_state/wt_new
! Subprogram not used                                       dwt_ptr0 = dwt_ptr0 - change_state
! Subprogram not used                                    else
! Subprogram not used                                       ptr = 0._r8
! Subprogram not used                                       dwt_ptr0 = dwt_ptr0 + init_state
! Subprogram not used                                    end if
! Subprogram not used                                    
! Subprogram not used                                    ! livecrootc_storage 
! Subprogram not used                                    ptr => pc13s%livecrootc_storage(p)
! Subprogram not used                                    init_state = ptr*wt_old
! Subprogram not used                                    change_state = ptr*dwt
! Subprogram not used                                    new_state = init_state+change_state
! Subprogram not used                                    if (wt_new /= 0._r8) then
! Subprogram not used                                       ptr = new_state/wt_new
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 + change_state
! Subprogram not used                                    else
! Subprogram not used                                       ptr = 0._r8
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 - init_state
! Subprogram not used                                    end if
! Subprogram not used                                    
! Subprogram not used                                    ! livecrootc_xfer 
! Subprogram not used                                    ptr => pc13s%livecrootc_xfer(p)
! Subprogram not used                                    init_state = ptr*wt_old
! Subprogram not used                                    change_state = ptr*dwt
! Subprogram not used                                    new_state = init_state+change_state
! Subprogram not used                                    if (wt_new /= 0._r8) then
! Subprogram not used                                       ptr = new_state/wt_new
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 + change_state
! Subprogram not used                                    else
! Subprogram not used                                       ptr = 0._r8
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 - init_state
! Subprogram not used                                    end if
! Subprogram not used                                    
! Subprogram not used                                    ! deadcrootc 
! Subprogram not used                                    ptr => pc13s%deadcrootc(p)
! Subprogram not used                                    dwt_ptr0 => dwt_deadcrootc13_to_litter(p)
! Subprogram not used                                    init_state = ptr*wt_old
! Subprogram not used                                    change_state = ptr*dwt
! Subprogram not used                                    new_state = init_state+change_state
! Subprogram not used                                    if (wt_new /= 0._r8) then
! Subprogram not used                                       ptr = new_state/wt_new
! Subprogram not used                                       dwt_ptr0 = dwt_ptr0 - change_state
! Subprogram not used                                    else
! Subprogram not used                                       ptr = 0._r8
! Subprogram not used                                       dwt_ptr0 = dwt_ptr0 + init_state
! Subprogram not used                                    end if
! Subprogram not used                                    
! Subprogram not used                                    ! deadcrootc_storage 
! Subprogram not used                                    ptr => pc13s%deadcrootc_storage(p)
! Subprogram not used                                    init_state = ptr*wt_old
! Subprogram not used                                    change_state = ptr*dwt
! Subprogram not used                                    new_state = init_state+change_state
! Subprogram not used                                    if (wt_new /= 0._r8) then
! Subprogram not used                                       ptr = new_state/wt_new
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 + change_state
! Subprogram not used                                    else
! Subprogram not used                                       ptr = 0._r8
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 - init_state
! Subprogram not used                                    end if
! Subprogram not used                                    
! Subprogram not used                                    ! deadcrootc_xfer 
! Subprogram not used                                    ptr => pc13s%deadcrootc_xfer(p)
! Subprogram not used                                    init_state = ptr*wt_old
! Subprogram not used                                    change_state = ptr*dwt
! Subprogram not used                                    new_state = init_state+change_state
! Subprogram not used                                    if (wt_new /= 0._r8) then
! Subprogram not used                                       ptr = new_state/wt_new
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 + change_state
! Subprogram not used                                    else
! Subprogram not used                                       ptr = 0._r8
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 - init_state
! Subprogram not used                                    end if
! Subprogram not used 
! Subprogram not used                                    ! gresp_storage 
! Subprogram not used                                    ptr => pc13s%gresp_storage(p)
! Subprogram not used                                    init_state = ptr*wt_old
! Subprogram not used                                    change_state = ptr*dwt
! Subprogram not used                                    new_state = init_state+change_state
! Subprogram not used                                    if (wt_new /= 0._r8) then
! Subprogram not used                                       ptr = new_state/wt_new
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 + change_state
! Subprogram not used                                    else
! Subprogram not used                                       ptr = 0._r8
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 - init_state
! Subprogram not used                                    end if
! Subprogram not used 
! Subprogram not used                                    ! gresp_xfer 
! Subprogram not used                                    ptr => pc13s%gresp_xfer(p)
! Subprogram not used                                    init_state = ptr*wt_old
! Subprogram not used                                    change_state = ptr*dwt
! Subprogram not used                                    new_state = init_state+change_state
! Subprogram not used                                    if (wt_new /= 0._r8) then
! Subprogram not used                                       ptr = new_state/wt_new
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 + change_state
! Subprogram not used                                    else
! Subprogram not used                                       ptr = 0._r8
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 - init_state
! Subprogram not used                                    end if
! Subprogram not used 
! Subprogram not used                                    ! cpool 
! Subprogram not used                                    ptr => pc13s%cpool(p)
! Subprogram not used                                    init_state = ptr*wt_old
! Subprogram not used                                    change_state = ptr*dwt
! Subprogram not used                                    new_state = init_state+change_state
! Subprogram not used                                    if (wt_new /= 0._r8) then
! Subprogram not used                                       ptr = new_state/wt_new
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 + change_state
! Subprogram not used                                    else
! Subprogram not used                                       ptr = 0._r8
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 - init_state
! Subprogram not used                                    end if
! Subprogram not used 
! Subprogram not used                                    ! pft_ctrunc 
! Subprogram not used                                    ptr => pc13s%pft_ctrunc(p)
! Subprogram not used                                    init_state = ptr*wt_old
! Subprogram not used                                    change_state = ptr*dwt
! Subprogram not used                                    new_state = init_state+change_state
! Subprogram not used                                    if (wt_new /= 0._r8) then
! Subprogram not used                                       ptr = new_state/wt_new
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 + change_state
! Subprogram not used                                    else
! Subprogram not used                                       ptr = 0._r8
! Subprogram not used                                       dwt_ptr1 = dwt_ptr1 - init_state
! Subprogram not used                                    end if
! Subprogram not used                                 endif
! Subprogram not used 
! Subprogram not used 				!---------------
! Subprogram not used 				! N state update
! Subprogram not used 				!---------------
! Subprogram not used 
! Subprogram not used 				! set pointers to the conversion and product pool fluxes for this pft
! Subprogram not used 				! dwt_ptr0 is reserved for local assignment to dwt_xxx_to_litter fluxes
! Subprogram not used 				dwt_ptr1 => conv_nflux(p)
! Subprogram not used 				dwt_ptr2 => prod10_nflux(p)
! Subprogram not used 				dwt_ptr3 => prod100_nflux(p)
! Subprogram not used 
! Subprogram not used 				! leafn 
! Subprogram not used 				ptr => pns%leafn(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 - init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! leafn_storage  
! Subprogram not used 				ptr => pns%leafn_storage(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 - init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! leafn_xfer  
! Subprogram not used 				ptr => pns%leafn_xfer(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 - init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! frootn 
! Subprogram not used 				ptr => pns%frootn(p)
! Subprogram not used 				dwt_ptr0 => dwt_frootn_to_litter(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					dwt_ptr0 = dwt_ptr0 - change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					dwt_ptr0 = dwt_ptr0 + init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! frootn_storage 
! Subprogram not used 				ptr => pns%frootn_storage(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 - init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! frootn_xfer  
! Subprogram not used 				ptr => pns%frootn_xfer(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 - init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! livestemn  
! Subprogram not used 				ptr => pns%livestemn(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 - init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! livestemn_storage 
! Subprogram not used 				ptr => pns%livestemn_storage(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 - init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! livestemn_xfer 
! Subprogram not used 				ptr => pns%livestemn_xfer(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 - init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! deadstemn 
! Subprogram not used 				ptr => pns%deadstemn(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 + change_state*pconv(pptr%itype(p))
! Subprogram not used 					dwt_ptr2 = dwt_ptr2 + change_state*pprod10(pptr%itype(p))
! Subprogram not used 					dwt_ptr3 = dwt_ptr3 + change_state*pprod100(pptr%itype(p))
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 - init_state*pconv(pptr%itype(p))
! Subprogram not used 					dwt_ptr2 = dwt_ptr2 - init_state*pprod10(pptr%itype(p))
! Subprogram not used 					dwt_ptr3 = dwt_ptr3 - init_state*pprod100(pptr%itype(p))
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! deadstemn_storage 
! Subprogram not used 				ptr => pns%deadstemn_storage(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 - init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! deadstemn_xfer 
! Subprogram not used 				ptr => pns%deadstemn_xfer(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 - init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! livecrootn 
! Subprogram not used 				ptr => pns%livecrootn(p)
! Subprogram not used 				dwt_ptr0 => dwt_livecrootn_to_litter(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					dwt_ptr0 = dwt_ptr0 - change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					dwt_ptr0 = dwt_ptr0 + init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! livecrootn_storage  
! Subprogram not used 				ptr => pns%livecrootn_storage(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 - init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! livecrootn_xfer  
! Subprogram not used 				ptr => pns%livecrootn_xfer(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 - init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! deadcrootn 
! Subprogram not used 				ptr => pns%deadcrootn(p)
! Subprogram not used 				dwt_ptr0 => dwt_deadcrootn_to_litter(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					dwt_ptr0 = dwt_ptr0 - change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					dwt_ptr0 = dwt_ptr0 + init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! deadcrootn_storage  
! Subprogram not used 				ptr => pns%deadcrootn_storage(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 - init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! deadcrootn_xfer  
! Subprogram not used 				ptr => pns%deadcrootn_xfer(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 - init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! retransn  
! Subprogram not used 				ptr => pns%retransn(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 - init_state
! Subprogram not used 				end if
! Subprogram not used 
! Subprogram not used 				! npool  
! Subprogram not used 				ptr => pns%npool(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 - init_state
! Subprogram not used 				end if
! Subprogram not used 				
! Subprogram not used 				! pft_ntrunc  
! Subprogram not used 				ptr => pns%pft_ntrunc(p)
! Subprogram not used 				init_state = ptr*wt_old
! Subprogram not used 				change_state = ptr*dwt
! Subprogram not used 				new_state = init_state+change_state
! Subprogram not used 				if (wt_new /= 0._r8) then
! Subprogram not used 					ptr = new_state/wt_new
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 + change_state
! Subprogram not used 				else
! Subprogram not used 					ptr = 0._r8
! Subprogram not used 					dwt_ptr1 = dwt_ptr1 - init_state
! Subprogram not used 				end if
! Subprogram not used 				
! Subprogram not used 			end if       ! weight decreasing
! Subprogram not used 		end if           ! is soil
! Subprogram not used 	end do               ! pft loop
! Subprogram not used     
! Subprogram not used 	! calculate column-level seeding fluxes
! Subprogram not used 	do pi = 1,max_pft_per_col
! Subprogram not used 		do c = begc, endc
! Subprogram not used 			if ( pi <=  cptr%npfts(c) ) then
! Subprogram not used 				p = cptr%pfti(c) + pi - 1
! Subprogram not used 				
! Subprogram not used 				! C fluxes
! Subprogram not used 				ccf%dwt_seedc_to_leaf(c) = ccf%dwt_seedc_to_leaf(c) + dwt_leafc_seed(p)/dt
! Subprogram not used 				ccf%dwt_seedc_to_deadstem(c) = ccf%dwt_seedc_to_deadstem(c) &
! Subprogram not used                                                                     + dwt_deadstemc_seed(p)/dt
! Subprogram not used 				
! Subprogram not used                                 ! C13 fluxes
! Subprogram not used                                 if (use_c13) then
! Subprogram not used                                    cc13f%dwt_seedc_to_leaf(c) = cc13f%dwt_seedc_to_leaf(c) + dwt_leafc13_seed(p)/dt
! Subprogram not used                                    cc13f%dwt_seedc_to_deadstem(c) = cc13f%dwt_seedc_to_deadstem(c) &
! Subprogram not used                                                                          + dwt_deadstemc13_seed(p)/dt
! Subprogram not used                                 endif
! Subprogram not used 				
! Subprogram not used 				! N fluxes
! Subprogram not used 				cnf%dwt_seedn_to_leaf(c) = cnf%dwt_seedn_to_leaf(c) + dwt_leafn_seed(p)/dt
! Subprogram not used 				cnf%dwt_seedn_to_deadstem(c) = cnf%dwt_seedn_to_deadstem(c) &
! Subprogram not used                                                                     + dwt_deadstemn_seed(p)/dt
! Subprogram not used 			end if
! Subprogram not used 		end do
! Subprogram not used 	end do
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 	! calculate pft-to-column for fluxes into litter and CWD pools
! Subprogram not used 	do pi = 1,max_pft_per_col
! Subprogram not used 		do c = begc, endc
! Subprogram not used 			if ( pi <=  cptr%npfts(c) ) then
! Subprogram not used 				p = cptr%pfti(c) + pi - 1
! Subprogram not used 
! Subprogram not used 				! fine root litter carbon fluxes
! Subprogram not used 				ccf%dwt_frootc_to_litr1c(c) = ccf%dwt_frootc_to_litr1c(c) + &
! Subprogram not used                                                             (dwt_frootc_to_litter(p)*pftcon%fr_flab(pptr%itype(p)))/dt
! Subprogram not used 				ccf%dwt_frootc_to_litr2c(c) = ccf%dwt_frootc_to_litr2c(c) + &
! Subprogram not used                                                             (dwt_frootc_to_litter(p)*pftcon%fr_fcel(pptr%itype(p)))/dt
! Subprogram not used 				ccf%dwt_frootc_to_litr3c(c) = ccf%dwt_frootc_to_litr3c(c) + &
! Subprogram not used                                                             (dwt_frootc_to_litter(p)*pftcon%fr_flig(pptr%itype(p)))/dt
! Subprogram not used 
! Subprogram not used 				! fine root litter C13 fluxes
! Subprogram not used                                 if (use_c13) then
! Subprogram not used                                    cc13f%dwt_frootc_to_litr1c(c) = cc13f%dwt_frootc_to_litr1c(c) + &
! Subprogram not used                                                                (dwt_frootc13_to_litter(p)*pftcon%fr_flab(pptr%itype(p)))/dt
! Subprogram not used                                    cc13f%dwt_frootc_to_litr2c(c) = cc13f%dwt_frootc_to_litr2c(c) + &
! Subprogram not used                                                                (dwt_frootc13_to_litter(p)*pftcon%fr_fcel(pptr%itype(p)))/dt
! Subprogram not used                                    cc13f%dwt_frootc_to_litr3c(c) = cc13f%dwt_frootc_to_litr3c(c) + &
! Subprogram not used                                                                (dwt_frootc13_to_litter(p)*pftcon%fr_flig(pptr%itype(p)))/dt
! Subprogram not used                                 endif
! Subprogram not used 
! Subprogram not used 				! fine root litter nitrogen fluxes
! Subprogram not used 				cnf%dwt_frootn_to_litr1n(c) = cnf%dwt_frootn_to_litr1n(c) + &
! Subprogram not used                                                             (dwt_frootn_to_litter(p)*pftcon%fr_flab(pptr%itype(p)))/dt
! Subprogram not used 				cnf%dwt_frootn_to_litr2n(c) = cnf%dwt_frootn_to_litr2n(c) + &
! Subprogram not used                                                             (dwt_frootn_to_litter(p)*pftcon%fr_fcel(pptr%itype(p)))/dt
! Subprogram not used 				cnf%dwt_frootn_to_litr3n(c) = cnf%dwt_frootn_to_litr3n(c) + &
! Subprogram not used                                                             (dwt_frootn_to_litter(p)*pftcon%fr_flig(pptr%itype(p)))/dt
! Subprogram not used 
! Subprogram not used 				! livecroot fluxes to cwd
! Subprogram not used 				ccf%dwt_livecrootc_to_cwdc(c) = ccf%dwt_livecrootc_to_cwdc(c) + &
! Subprogram not used                                                             (dwt_livecrootc_to_litter(p))/dt
! Subprogram not used                                 if (use_c13) then
! Subprogram not used                                    cc13f%dwt_livecrootc_to_cwdc(c) = cc13f%dwt_livecrootc_to_cwdc(c) + &
! Subprogram not used                                                                (dwt_livecrootc13_to_litter(p))/dt
! Subprogram not used                                 endif
! Subprogram not used 				cnf%dwt_livecrootn_to_cwdn(c) = cnf%dwt_livecrootn_to_cwdn(c) + &
! Subprogram not used                                                             (dwt_livecrootn_to_litter(p))/dt
! Subprogram not used 
! Subprogram not used 				! deadcroot fluxes to cwd
! Subprogram not used 				ccf%dwt_deadcrootc_to_cwdc(c) = ccf%dwt_deadcrootc_to_cwdc(c) + &
! Subprogram not used                                                             (dwt_deadcrootc_to_litter(p))/dt
! Subprogram not used                                 if (use_c13) then
! Subprogram not used                                    cc13f%dwt_deadcrootc_to_cwdc(c) = cc13f%dwt_deadcrootc_to_cwdc(c) + &
! Subprogram not used                                                                (dwt_deadcrootc13_to_litter(p))/dt
! Subprogram not used                                 endif
! Subprogram not used 				cnf%dwt_deadcrootn_to_cwdn(c) = cnf%dwt_deadcrootn_to_cwdn(c) + &
! Subprogram not used                                                             (dwt_deadcrootn_to_litter(p))/dt
! Subprogram not used 			end if
! Subprogram not used 		end do
! Subprogram not used 	end do
! Subprogram not used 
! Subprogram not used 	! calculate pft-to-column for fluxes into product pools and conversion flux
! Subprogram not used 	do pi = 1,max_pft_per_col
! Subprogram not used 		do c = begc,endc
! Subprogram not used 			if (pi <= cptr%npfts(c)) then
! Subprogram not used 				p = cptr%pfti(c) + pi - 1
! Subprogram not used 
! Subprogram not used 				! column-level fluxes are accumulated as positive fluxes.
! Subprogram not used 				! column-level C flux updates
! Subprogram not used 				ccf%dwt_conv_cflux(c) = ccf%dwt_conv_cflux(c) - conv_cflux(p)/dt
! Subprogram not used 				ccf%dwt_prod10c_gain(c) = ccf%dwt_prod10c_gain(c) - prod10_cflux(p)/dt
! Subprogram not used 				ccf%dwt_prod100c_gain(c) = ccf%dwt_prod100c_gain(c) - prod100_cflux(p)/dt
! Subprogram not used 
! Subprogram not used 				! column-level C13 flux updates
! Subprogram not used                                 if (use_c13) then
! Subprogram not used                                    cc13f%dwt_conv_cflux(c) = cc13f%dwt_conv_cflux(c) - conv_c13flux(p)/dt
! Subprogram not used                                    cc13f%dwt_prod10c_gain(c) = cc13f%dwt_prod10c_gain(c) - prod10_c13flux(p)/dt
! Subprogram not used                                    cc13f%dwt_prod100c_gain(c) = cc13f%dwt_prod100c_gain(c) - prod100_c13flux(p)/dt
! Subprogram not used                                 endif
! Subprogram not used 
! Subprogram not used 				! column-level N flux updates
! Subprogram not used 				cnf%dwt_conv_nflux(c) = cnf%dwt_conv_nflux(c) - conv_nflux(p)/dt
! Subprogram not used 				cnf%dwt_prod10n_gain(c) = cnf%dwt_prod10n_gain(c) - prod10_nflux(p)/dt
! Subprogram not used 				cnf%dwt_prod100n_gain(c) = cnf%dwt_prod100n_gain(c) - prod100_nflux(p)/dt
! Subprogram not used 
! Subprogram not used 			end if
! Subprogram not used 		end do
! Subprogram not used 	end do
! Subprogram not used 
! Subprogram not used 	! Deallocate pft-level flux arrays
! Subprogram not used         deallocate(dwt_leafc_seed)
! Subprogram not used         deallocate(dwt_leafn_seed)
! Subprogram not used         if (use_c13) then
! Subprogram not used            deallocate(dwt_leafc13_seed)
! Subprogram not used         endif
! Subprogram not used         deallocate(dwt_deadstemc_seed)
! Subprogram not used         deallocate(dwt_deadstemn_seed)
! Subprogram not used         if (use_c13) then
! Subprogram not used            deallocate(dwt_deadstemc13_seed)
! Subprogram not used         endif
! Subprogram not used 	deallocate(dwt_frootc_to_litter)
! Subprogram not used 	deallocate(dwt_livecrootc_to_litter)
! Subprogram not used 	deallocate(dwt_deadcrootc_to_litter)
! Subprogram not used         if (use_c13) then
! Subprogram not used            deallocate(dwt_frootc13_to_litter)
! Subprogram not used            deallocate(dwt_livecrootc13_to_litter)
! Subprogram not used            deallocate(dwt_deadcrootc13_to_litter)
! Subprogram not used         endif
! Subprogram not used 	deallocate(dwt_frootn_to_litter)
! Subprogram not used 	deallocate(dwt_livecrootn_to_litter)
! Subprogram not used 	deallocate(dwt_deadcrootn_to_litter)
! Subprogram not used 	deallocate(conv_cflux)
! Subprogram not used 	deallocate(prod10_cflux)
! Subprogram not used 	deallocate(prod100_cflux)
! Subprogram not used         if (use_c13) then
! Subprogram not used            deallocate(conv_c13flux)
! Subprogram not used            deallocate(prod10_c13flux)
! Subprogram not used            deallocate(prod100_c13flux)
! Subprogram not used         endif
! Subprogram not used 	deallocate(conv_nflux)
! Subprogram not used 	deallocate(prod10_nflux)
! Subprogram not used 	deallocate(prod100_nflux)
! Subprogram not used     
! Subprogram not used end subroutine pftdyn_cnbal

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: pftwt_init
!
! !INTERFACE:
! Subprogram not used   subroutine pftwt_init()
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used ! Initialize time interpolation of cndv pft weights from annual to time step
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used   use clm_varctl, only : nsrest, nsrStartup
! Subprogram not used !
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used ! !LOCAL VARIABLES:
! Subprogram not used     integer  :: ier, p                        ! error status, do-loop index
! Subprogram not used     integer  :: begp,endp                     ! beg/end indices for land pfts
! Subprogram not used     character(len=32) :: subname='pftwt_init' ! subroutine name
! Subprogram not used     type(pft_type), pointer :: pptr           ! ponter to pft derived subtype
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     pptr => pft
! Subprogram not used 
! Subprogram not used     call get_proc_bounds(begp=begp,endp=endp)
! Subprogram not used 
! Subprogram not used     allocate(wtcol_old(begp:endp),stat=ier)
! Subprogram not used     if (ier /= 0) then
! Subprogram not used        call endrun( subname//'::ERROR: pftwt_init allocation error for wtcol_old')
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     if (nsrest == nsrStartup) then
! Subprogram not used        do p = begp,endp
! Subprogram not used           pdgvs%fpcgrid(p) = pptr%wtcol(p)
! Subprogram not used           pdgvs%fpcgridold(p) = pptr%wtcol(p)
! Subprogram not used           wtcol_old(p) = pptr%wtcol(p)
! Subprogram not used        end do
! Subprogram not used      else
! Subprogram not used        do p = begp,endp
! Subprogram not used           wtcol_old(p) = pptr%wtcol(p)
! Subprogram not used        end do
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used   end subroutine pftwt_init

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: pftwt_interp
!
! !INTERFACE:
! Subprogram not used   subroutine pftwt_interp( begp, endp )
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used ! Time interpolate cndv pft weights from annual to time step
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used     use clm_time_manager, only : get_curr_calday, get_curr_date, &
! Subprogram not used                                  get_days_per_year
! Subprogram not used     use clm_time_manager, only : get_step_size, get_nstep
! Subprogram not used     use clm_varcon      , only : istsoil ! CNDV incompatible with dynLU
! Subprogram not used     use clm_varctl      , only : finidat
! Subprogram not used !
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     integer, intent(IN)  :: begp,endp                ! beg/end indices for land pfts
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used ! !LOCAL VARIABLES:
! Subprogram not used     integer  :: c,g,l,p            ! indices
! Subprogram not used     real(r8) :: cday               ! current calendar day (1.0 = 0Z on Jan 1)
! Subprogram not used     real(r8) :: wt1                ! time interpolation weights
! Subprogram not used     real(r8) :: dtime              ! model time step
! Subprogram not used     real(r8) :: days_per_year      ! days per year
! Subprogram not used     integer  :: nstep              ! time step number
! Subprogram not used     integer  :: year               ! year (0, ...) at nstep + 1
! Subprogram not used     integer  :: mon                ! month (1, ..., 12) at nstep + 1
! Subprogram not used     integer  :: day                ! day of month (1, ..., 31) at nstep + 1
! Subprogram not used     integer  :: sec                ! seconds into current date at nstep + 1
! Subprogram not used     type(landunit_type), pointer :: lptr ! pointer to landunit derived subtype
! Subprogram not used     type(pft_type)     , pointer :: pptr ! ...     to pft derived subtype
! Subprogram not used     character(len=32) :: subname='pftwt_interp' ! subroutine name
! Subprogram not used 
! Subprogram not used ! !CALLED FROM:
! Subprogram not used !  subr. driver
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     ! Set pointers into derived type
! Subprogram not used 
! Subprogram not used     lptr => lun
! Subprogram not used     pptr => pft
! Subprogram not used 
! Subprogram not used     ! Interpolate pft weight to current time step
! Subprogram not used     ! Map interpolated pctpft to subgrid weights
! Subprogram not used     ! assumes maxpatch_pft = numpft + 1, each landunit has 1 column, 
! Subprogram not used     ! SCAM not defined and create_croplandunit = .false.
! Subprogram not used 
! Subprogram not used     nstep         = get_nstep()
! Subprogram not used     dtime         = get_step_size()
! Subprogram not used     cday          = get_curr_calday(offset=-int(dtime))
! Subprogram not used     days_per_year = get_days_per_year()
! Subprogram not used 
! Subprogram not used     wt1 = ((days_per_year + 1._r8) - cday)/days_per_year
! Subprogram not used 
! Subprogram not used     call get_curr_date(year, mon, day, sec, offset=int(dtime))
! Subprogram not used 
! Subprogram not used     do p = begp,endp
! Subprogram not used        g = pptr%gridcell(p)
! Subprogram not used        l = pptr%landunit(p)
! Subprogram not used 
! Subprogram not used        if (lptr%itype(l) == istsoil .and. lptr%wtgcell(l) > 0._r8) then ! CNDV incompatible with dynLU
! Subprogram not used           wtcol_old(p)    = pptr%wtcol(p)
! Subprogram not used           pptr%wtcol(p)   = pdgvs%fpcgrid(p) + &
! Subprogram not used                      wt1 * (pdgvs%fpcgridold(p) - pdgvs%fpcgrid(p))
! Subprogram not used           pptr%wtlunit(p) = pptr%wtcol(p)
! Subprogram not used           pptr%wtgcell(p) = pptr%wtcol(p) * lptr%wtgcell(l)
! Subprogram not used 
! Subprogram not used           if (mon==1 .and. day==1 .and. sec==dtime .and. nstep>0) then
! Subprogram not used              pdgvs%fpcgridold(p) = pdgvs%fpcgrid(p)
! Subprogram not used           end if
! Subprogram not used        end if
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used   end subroutine pftwt_interp

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNHarvest
!
! !INTERFACE:
! Subprogram not used subroutine CNHarvest (num_soilc, filter_soilc, num_soilp, filter_soilp)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used ! Harvest mortality routine for coupled carbon-nitrogen code (CN)
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used    use clmtype
! Subprogram not used    use pftvarcon       , only : noveg, nbrdlf_evr_shrub, pprodharv10
! Subprogram not used    use clm_varcon      , only : secspday
! Subprogram not used    use clm_time_manager, only : get_days_per_year
! Subprogram not used !
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used    implicit none
! Subprogram not used    integer, intent(in) :: num_soilc       ! number of soil columns in filter
! Subprogram not used    integer, intent(in) :: filter_soilc(:) ! column filter for soil points
! Subprogram not used    integer, intent(in) :: num_soilp       ! number of soil pfts in filter
! Subprogram not used    integer, intent(in) :: filter_soilp(:) ! pft filter for soil points
! Subprogram not used !
! Subprogram not used ! !CALLED FROM:
! Subprogram not used ! subroutine CNEcosystemDyn
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 3/29/04: Created by Peter Thornton
! Subprogram not used !
! Subprogram not used ! !LOCAL VARIABLES:
! Subprogram not used !
! Subprogram not used ! local pointers to implicit in arrays
! Subprogram not used    integer , pointer :: pgridcell(:)   ! pft-level index into gridcell-level quantities
! Subprogram not used    integer , pointer :: ivt(:)         ! pft vegetation type
! Subprogram not used 
! Subprogram not used    real(r8), pointer :: leafc(:)              ! (gC/m2) leaf C
! Subprogram not used    real(r8), pointer :: frootc(:)             ! (gC/m2) fine root C
! Subprogram not used    real(r8), pointer :: livestemc(:)          ! (gC/m2) live stem C
! Subprogram not used    real(r8), pointer :: deadstemc(:)          ! (gC/m2) dead stem C
! Subprogram not used    real(r8), pointer :: livecrootc(:)         ! (gC/m2) live coarse root C
! Subprogram not used    real(r8), pointer :: deadcrootc(:)         ! (gC/m2) dead coarse root C
! Subprogram not used    real(r8), pointer :: xsmrpool(:)           ! (gC/m2) abstract C pool to meet excess MR demand
! Subprogram not used    real(r8), pointer :: leafc_storage(:)      ! (gC/m2) leaf C storage
! Subprogram not used    real(r8), pointer :: frootc_storage(:)     ! (gC/m2) fine root C storage
! Subprogram not used    real(r8), pointer :: livestemc_storage(:)  ! (gC/m2) live stem C storage
! Subprogram not used    real(r8), pointer :: deadstemc_storage(:)  ! (gC/m2) dead stem C storage
! Subprogram not used    real(r8), pointer :: livecrootc_storage(:) ! (gC/m2) live coarse root C storage
! Subprogram not used    real(r8), pointer :: deadcrootc_storage(:) ! (gC/m2) dead coarse root C storage
! Subprogram not used    real(r8), pointer :: gresp_storage(:)      ! (gC/m2) growth respiration storage
! Subprogram not used    real(r8), pointer :: leafc_xfer(:)         ! (gC/m2) leaf C transfer
! Subprogram not used    real(r8), pointer :: frootc_xfer(:)        ! (gC/m2) fine root C transfer
! Subprogram not used    real(r8), pointer :: livestemc_xfer(:)     ! (gC/m2) live stem C transfer
! Subprogram not used    real(r8), pointer :: deadstemc_xfer(:)     ! (gC/m2) dead stem C transfer
! Subprogram not used    real(r8), pointer :: livecrootc_xfer(:)    ! (gC/m2) live coarse root C transfer
! Subprogram not used    real(r8), pointer :: deadcrootc_xfer(:)    ! (gC/m2) dead coarse root C transfer
! Subprogram not used    real(r8), pointer :: gresp_xfer(:)         ! (gC/m2) growth respiration transfer
! Subprogram not used    real(r8), pointer :: leafn(:)              ! (gN/m2) leaf N
! Subprogram not used    real(r8), pointer :: frootn(:)             ! (gN/m2) fine root N
! Subprogram not used    real(r8), pointer :: livestemn(:)          ! (gN/m2) live stem N
! Subprogram not used    real(r8), pointer :: deadstemn(:)          ! (gN/m2) dead stem N
! Subprogram not used    real(r8), pointer :: livecrootn(:)         ! (gN/m2) live coarse root N
! Subprogram not used    real(r8), pointer :: deadcrootn(:)         ! (gN/m2) dead coarse root N
! Subprogram not used    real(r8), pointer :: retransn(:)           ! (gN/m2) plant pool of retranslocated N
! Subprogram not used    real(r8), pointer :: leafn_storage(:)      ! (gN/m2) leaf N storage
! Subprogram not used    real(r8), pointer :: frootn_storage(:)     ! (gN/m2) fine root N storage
! Subprogram not used    real(r8), pointer :: livestemn_storage(:)  ! (gN/m2) live stem N storage
! Subprogram not used    real(r8), pointer :: deadstemn_storage(:)  ! (gN/m2) dead stem N storage
! Subprogram not used    real(r8), pointer :: livecrootn_storage(:) ! (gN/m2) live coarse root N storage
! Subprogram not used    real(r8), pointer :: deadcrootn_storage(:) ! (gN/m2) dead coarse root N storage
! Subprogram not used    real(r8), pointer :: leafn_xfer(:)         ! (gN/m2) leaf N transfer
! Subprogram not used    real(r8), pointer :: frootn_xfer(:)        ! (gN/m2) fine root N transfer
! Subprogram not used    real(r8), pointer :: livestemn_xfer(:)     ! (gN/m2) live stem N transfer
! Subprogram not used    real(r8), pointer :: deadstemn_xfer(:)     ! (gN/m2) dead stem N transfer
! Subprogram not used    real(r8), pointer :: livecrootn_xfer(:)    ! (gN/m2) live coarse root N transfer
! Subprogram not used    real(r8), pointer :: deadcrootn_xfer(:)    ! (gN/m2) dead coarse root N transfer
! Subprogram not used !
! Subprogram not used ! local pointers to implicit in/out arrays
! Subprogram not used !
! Subprogram not used ! local pointers to implicit out arrays
! Subprogram not used    real(r8), pointer :: hrv_leafc_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_frootc_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_livestemc_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_deadstemc_to_prod10c(:)
! Subprogram not used    real(r8), pointer :: hrv_deadstemc_to_prod100c(:)
! Subprogram not used    real(r8), pointer :: hrv_livecrootc_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_deadcrootc_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_xsmrpool_to_atm(:)
! Subprogram not used    real(r8), pointer :: hrv_leafc_storage_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_frootc_storage_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_livestemc_storage_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_deadstemc_storage_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_livecrootc_storage_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_deadcrootc_storage_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_gresp_storage_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_leafc_xfer_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_frootc_xfer_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_livestemc_xfer_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_deadstemc_xfer_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_livecrootc_xfer_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_deadcrootc_xfer_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_gresp_xfer_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_leafn_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_frootn_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_livestemn_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_deadstemn_to_prod10n(:)
! Subprogram not used    real(r8), pointer :: hrv_deadstemn_to_prod100n(:)
! Subprogram not used    real(r8), pointer :: hrv_livecrootn_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_deadcrootn_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_retransn_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_leafn_storage_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_frootn_storage_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_livestemn_storage_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_deadstemn_storage_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_livecrootn_storage_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_deadcrootn_storage_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_leafn_xfer_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_frootn_xfer_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_livestemn_xfer_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_deadstemn_xfer_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_livecrootn_xfer_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_deadcrootn_xfer_to_litter(:)
! Subprogram not used !
! Subprogram not used ! !OTHER LOCAL VARIABLES:
! Subprogram not used    integer :: p                         ! pft index
! Subprogram not used    integer :: g                         ! gridcell index
! Subprogram not used    integer :: fp                        ! pft filter index
! Subprogram not used    real(r8):: am                        ! rate for fractional harvest mortality (1/yr)
! Subprogram not used    real(r8):: m                         ! rate for fractional harvest mortality (1/s)
! Subprogram not used    real(r8):: days_per_year             ! days per year
! Subprogram not used !EOP
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    ! assign local pointers to pft-level arrays
! Subprogram not used    pgridcell                      => pft%gridcell
! Subprogram not used    
! Subprogram not used    ivt                            => pft%itype
! Subprogram not used    leafc                          => pcs%leafc
! Subprogram not used    frootc                         => pcs%frootc
! Subprogram not used    livestemc                      => pcs%livestemc
! Subprogram not used    deadstemc                      => pcs%deadstemc
! Subprogram not used    livecrootc                     => pcs%livecrootc
! Subprogram not used    deadcrootc                     => pcs%deadcrootc
! Subprogram not used    xsmrpool                       => pcs%xsmrpool
! Subprogram not used    leafc_storage                  => pcs%leafc_storage
! Subprogram not used    frootc_storage                 => pcs%frootc_storage
! Subprogram not used    livestemc_storage              => pcs%livestemc_storage
! Subprogram not used    deadstemc_storage              => pcs%deadstemc_storage
! Subprogram not used    livecrootc_storage             => pcs%livecrootc_storage
! Subprogram not used    deadcrootc_storage             => pcs%deadcrootc_storage
! Subprogram not used    gresp_storage                  => pcs%gresp_storage
! Subprogram not used    leafc_xfer                     => pcs%leafc_xfer
! Subprogram not used    frootc_xfer                    => pcs%frootc_xfer
! Subprogram not used    livestemc_xfer                 => pcs%livestemc_xfer
! Subprogram not used    deadstemc_xfer                 => pcs%deadstemc_xfer
! Subprogram not used    livecrootc_xfer                => pcs%livecrootc_xfer
! Subprogram not used    deadcrootc_xfer                => pcs%deadcrootc_xfer
! Subprogram not used    gresp_xfer                     => pcs%gresp_xfer
! Subprogram not used    leafn                          => pns%leafn
! Subprogram not used    frootn                         => pns%frootn
! Subprogram not used    livestemn                      => pns%livestemn
! Subprogram not used    deadstemn                      => pns%deadstemn
! Subprogram not used    livecrootn                     => pns%livecrootn
! Subprogram not used    deadcrootn                     => pns%deadcrootn
! Subprogram not used    retransn                       => pns%retransn
! Subprogram not used    leafn_storage                  => pns%leafn_storage
! Subprogram not used    frootn_storage                 => pns%frootn_storage
! Subprogram not used    livestemn_storage              => pns%livestemn_storage
! Subprogram not used    deadstemn_storage              => pns%deadstemn_storage
! Subprogram not used    livecrootn_storage             => pns%livecrootn_storage
! Subprogram not used    deadcrootn_storage             => pns%deadcrootn_storage
! Subprogram not used    leafn_xfer                     => pns%leafn_xfer
! Subprogram not used    frootn_xfer                    => pns%frootn_xfer
! Subprogram not used    livestemn_xfer                 => pns%livestemn_xfer
! Subprogram not used    deadstemn_xfer                 => pns%deadstemn_xfer
! Subprogram not used    livecrootn_xfer                => pns%livecrootn_xfer
! Subprogram not used    deadcrootn_xfer                => pns%deadcrootn_xfer
! Subprogram not used    hrv_leafc_to_litter              => pcf%hrv_leafc_to_litter
! Subprogram not used    hrv_frootc_to_litter             => pcf%hrv_frootc_to_litter
! Subprogram not used    hrv_livestemc_to_litter          => pcf%hrv_livestemc_to_litter
! Subprogram not used    hrv_deadstemc_to_prod10c         => pcf%hrv_deadstemc_to_prod10c
! Subprogram not used    hrv_deadstemc_to_prod100c        => pcf%hrv_deadstemc_to_prod100c
! Subprogram not used    hrv_livecrootc_to_litter         => pcf%hrv_livecrootc_to_litter
! Subprogram not used    hrv_deadcrootc_to_litter         => pcf%hrv_deadcrootc_to_litter
! Subprogram not used    hrv_xsmrpool_to_atm              => pcf%hrv_xsmrpool_to_atm
! Subprogram not used    hrv_leafc_storage_to_litter      => pcf%hrv_leafc_storage_to_litter
! Subprogram not used    hrv_frootc_storage_to_litter     => pcf%hrv_frootc_storage_to_litter
! Subprogram not used    hrv_livestemc_storage_to_litter  => pcf%hrv_livestemc_storage_to_litter
! Subprogram not used    hrv_deadstemc_storage_to_litter  => pcf%hrv_deadstemc_storage_to_litter
! Subprogram not used    hrv_livecrootc_storage_to_litter => pcf%hrv_livecrootc_storage_to_litter
! Subprogram not used    hrv_deadcrootc_storage_to_litter => pcf%hrv_deadcrootc_storage_to_litter
! Subprogram not used    hrv_gresp_storage_to_litter      => pcf%hrv_gresp_storage_to_litter
! Subprogram not used    hrv_leafc_xfer_to_litter         => pcf%hrv_leafc_xfer_to_litter
! Subprogram not used    hrv_frootc_xfer_to_litter        => pcf%hrv_frootc_xfer_to_litter
! Subprogram not used    hrv_livestemc_xfer_to_litter     => pcf%hrv_livestemc_xfer_to_litter
! Subprogram not used    hrv_deadstemc_xfer_to_litter     => pcf%hrv_deadstemc_xfer_to_litter
! Subprogram not used    hrv_livecrootc_xfer_to_litter    => pcf%hrv_livecrootc_xfer_to_litter
! Subprogram not used    hrv_deadcrootc_xfer_to_litter    => pcf%hrv_deadcrootc_xfer_to_litter
! Subprogram not used    hrv_gresp_xfer_to_litter         => pcf%hrv_gresp_xfer_to_litter
! Subprogram not used    hrv_leafn_to_litter              => pnf%hrv_leafn_to_litter
! Subprogram not used    hrv_frootn_to_litter             => pnf%hrv_frootn_to_litter
! Subprogram not used    hrv_livestemn_to_litter          => pnf%hrv_livestemn_to_litter
! Subprogram not used    hrv_deadstemn_to_prod10n         => pnf%hrv_deadstemn_to_prod10n
! Subprogram not used    hrv_deadstemn_to_prod100n        => pnf%hrv_deadstemn_to_prod100n
! Subprogram not used    hrv_livecrootn_to_litter         => pnf%hrv_livecrootn_to_litter
! Subprogram not used    hrv_deadcrootn_to_litter         => pnf%hrv_deadcrootn_to_litter
! Subprogram not used    hrv_retransn_to_litter           => pnf%hrv_retransn_to_litter
! Subprogram not used    hrv_leafn_storage_to_litter      => pnf%hrv_leafn_storage_to_litter
! Subprogram not used    hrv_frootn_storage_to_litter     => pnf%hrv_frootn_storage_to_litter
! Subprogram not used    hrv_livestemn_storage_to_litter  => pnf%hrv_livestemn_storage_to_litter
! Subprogram not used    hrv_deadstemn_storage_to_litter  => pnf%hrv_deadstemn_storage_to_litter
! Subprogram not used    hrv_livecrootn_storage_to_litter => pnf%hrv_livecrootn_storage_to_litter
! Subprogram not used    hrv_deadcrootn_storage_to_litter => pnf%hrv_deadcrootn_storage_to_litter
! Subprogram not used    hrv_leafn_xfer_to_litter         => pnf%hrv_leafn_xfer_to_litter
! Subprogram not used    hrv_frootn_xfer_to_litter        => pnf%hrv_frootn_xfer_to_litter
! Subprogram not used    hrv_livestemn_xfer_to_litter     => pnf%hrv_livestemn_xfer_to_litter
! Subprogram not used    hrv_deadstemn_xfer_to_litter     => pnf%hrv_deadstemn_xfer_to_litter
! Subprogram not used    hrv_livecrootn_xfer_to_litter    => pnf%hrv_livecrootn_xfer_to_litter
! Subprogram not used    hrv_deadcrootn_xfer_to_litter    => pnf%hrv_deadcrootn_xfer_to_litter
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    days_per_year = get_days_per_year()
! Subprogram not used 
! Subprogram not used    ! pft loop
! Subprogram not used    do fp = 1,num_soilp
! Subprogram not used       p = filter_soilp(fp)
! Subprogram not used       g = pgridcell(p)
! Subprogram not used       
! Subprogram not used       ! If this is a tree pft, then
! Subprogram not used       ! get the annual harvest "mortality" rate (am) from harvest array
! Subprogram not used       ! and convert to rate per second
! Subprogram not used       if (ivt(p) > noveg .and. ivt(p) < nbrdlf_evr_shrub) then
! Subprogram not used 
! Subprogram not used          if (do_harvest) then
! Subprogram not used             am = harvest(g)
! Subprogram not used             m  = am/(days_per_year * secspday)
! Subprogram not used          else
! Subprogram not used             m = 0._r8
! Subprogram not used          end if   
! Subprogram not used 
! Subprogram not used          ! pft-level harvest carbon fluxes
! Subprogram not used          ! displayed pools
! Subprogram not used          hrv_leafc_to_litter(p)               = leafc(p)               * m
! Subprogram not used          hrv_frootc_to_litter(p)              = frootc(p)              * m
! Subprogram not used          hrv_livestemc_to_litter(p)           = livestemc(p)           * m
! Subprogram not used          hrv_deadstemc_to_prod10c(p)          = deadstemc(p)           * m * &
! Subprogram not used                                                 pprodharv10(ivt(p))
! Subprogram not used          hrv_deadstemc_to_prod100c(p)         = deadstemc(p)           * m * &
! Subprogram not used                                                 (1.0_r8 - pprodharv10(ivt(p)))
! Subprogram not used          hrv_livecrootc_to_litter(p)          = livecrootc(p)          * m
! Subprogram not used          hrv_deadcrootc_to_litter(p)          = deadcrootc(p)          * m
! Subprogram not used          hrv_xsmrpool_to_atm(p)               = xsmrpool(p)            * m
! Subprogram not used 
! Subprogram not used          ! storage pools
! Subprogram not used          hrv_leafc_storage_to_litter(p)       = leafc_storage(p)       * m
! Subprogram not used          hrv_frootc_storage_to_litter(p)      = frootc_storage(p)      * m
! Subprogram not used          hrv_livestemc_storage_to_litter(p)   = livestemc_storage(p)   * m
! Subprogram not used          hrv_deadstemc_storage_to_litter(p)   = deadstemc_storage(p)   * m
! Subprogram not used          hrv_livecrootc_storage_to_litter(p)  = livecrootc_storage(p)  * m
! Subprogram not used          hrv_deadcrootc_storage_to_litter(p)  = deadcrootc_storage(p)  * m
! Subprogram not used          hrv_gresp_storage_to_litter(p)       = gresp_storage(p)       * m
! Subprogram not used 
! Subprogram not used          ! transfer pools
! Subprogram not used          hrv_leafc_xfer_to_litter(p)          = leafc_xfer(p)          * m
! Subprogram not used          hrv_frootc_xfer_to_litter(p)         = frootc_xfer(p)         * m
! Subprogram not used          hrv_livestemc_xfer_to_litter(p)      = livestemc_xfer(p)      * m
! Subprogram not used          hrv_deadstemc_xfer_to_litter(p)      = deadstemc_xfer(p)      * m
! Subprogram not used          hrv_livecrootc_xfer_to_litter(p)     = livecrootc_xfer(p)     * m
! Subprogram not used          hrv_deadcrootc_xfer_to_litter(p)     = deadcrootc_xfer(p)     * m
! Subprogram not used          hrv_gresp_xfer_to_litter(p)          = gresp_xfer(p)          * m
! Subprogram not used 
! Subprogram not used          ! pft-level harvest mortality nitrogen fluxes
! Subprogram not used          ! displayed pools
! Subprogram not used          hrv_leafn_to_litter(p)               = leafn(p)               * m
! Subprogram not used          hrv_frootn_to_litter(p)              = frootn(p)              * m
! Subprogram not used          hrv_livestemn_to_litter(p)           = livestemn(p)           * m
! Subprogram not used          hrv_deadstemn_to_prod10n(p)          = deadstemn(p)           * m * &
! Subprogram not used                                                 pprodharv10(ivt(p))
! Subprogram not used          hrv_deadstemn_to_prod100n(p)         = deadstemn(p)           * m * &
! Subprogram not used                                                 (1.0_r8 - pprodharv10(ivt(p)))
! Subprogram not used          hrv_livecrootn_to_litter(p)          = livecrootn(p)          * m
! Subprogram not used          hrv_deadcrootn_to_litter(p)          = deadcrootn(p)          * m
! Subprogram not used          hrv_retransn_to_litter(p)            = retransn(p)            * m
! Subprogram not used 
! Subprogram not used          ! storage pools
! Subprogram not used          hrv_leafn_storage_to_litter(p)       = leafn_storage(p)       * m
! Subprogram not used          hrv_frootn_storage_to_litter(p)      = frootn_storage(p)      * m
! Subprogram not used          hrv_livestemn_storage_to_litter(p)   = livestemn_storage(p)   * m
! Subprogram not used          hrv_deadstemn_storage_to_litter(p)   = deadstemn_storage(p)   * m
! Subprogram not used          hrv_livecrootn_storage_to_litter(p)  = livecrootn_storage(p)  * m
! Subprogram not used          hrv_deadcrootn_storage_to_litter(p)  = deadcrootn_storage(p)  * m
! Subprogram not used 
! Subprogram not used          ! transfer pools
! Subprogram not used          hrv_leafn_xfer_to_litter(p)          = leafn_xfer(p)          * m
! Subprogram not used          hrv_frootn_xfer_to_litter(p)         = frootn_xfer(p)         * m
! Subprogram not used          hrv_livestemn_xfer_to_litter(p)      = livestemn_xfer(p)      * m
! Subprogram not used          hrv_deadstemn_xfer_to_litter(p)      = deadstemn_xfer(p)      * m
! Subprogram not used          hrv_livecrootn_xfer_to_litter(p)     = livecrootn_xfer(p)     * m
! Subprogram not used          hrv_deadcrootn_xfer_to_litter(p)     = deadcrootn_xfer(p)     * m
! Subprogram not used          
! Subprogram not used       end if  ! end tree block
! Subprogram not used 
! Subprogram not used    end do ! end of pft loop
! Subprogram not used 
! Subprogram not used    ! gather all pft-level litterfall fluxes from harvest to the column
! Subprogram not used    ! for litter C and N inputs
! Subprogram not used 
! Subprogram not used    call CNHarvestPftToColumn(num_soilc, filter_soilc)
! Subprogram not used 
! Subprogram not used end subroutine CNHarvest
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNHarvestPftToColumn
!
! !INTERFACE:
! Subprogram not used subroutine CNHarvestPftToColumn (num_soilc, filter_soilc)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used ! called at the end of CNHarvest to gather all pft-level harvest litterfall fluxes
! Subprogram not used ! to the column level and assign them to the three litter pools
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used   use clmtype
! Subprogram not used   use clm_varpar, only : max_pft_per_col, maxpatch_pft
! Subprogram not used !
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used   implicit none
! Subprogram not used   integer, intent(in) :: num_soilc       ! number of soil columns in filter
! Subprogram not used   integer, intent(in) :: filter_soilc(:) ! soil column filter
! Subprogram not used !
! Subprogram not used ! !CALLED FROM:
! Subprogram not used ! subroutine CNphenology
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 9/8/03: Created by Peter Thornton
! Subprogram not used !
! Subprogram not used ! !LOCAL VARIABLES:
! Subprogram not used !
! Subprogram not used ! local pointers to implicit in scalars
! Subprogram not used    integer , pointer :: ivt(:)      ! pft vegetation type
! Subprogram not used    real(r8), pointer :: wtcol(:)    ! pft weight relative to column (0-1)
! Subprogram not used    real(r8), pointer :: pwtgcell(:) ! weight of pft relative to corresponding gridcell
! Subprogram not used    real(r8), pointer :: lf_flab(:)  ! leaf litter labile fraction
! Subprogram not used    real(r8), pointer :: lf_fcel(:)  ! leaf litter cellulose fraction
! Subprogram not used    real(r8), pointer :: lf_flig(:)  ! leaf litter lignin fraction
! Subprogram not used    real(r8), pointer :: fr_flab(:)  ! fine root litter labile fraction
! Subprogram not used    real(r8), pointer :: fr_fcel(:)  ! fine root litter cellulose fraction
! Subprogram not used    real(r8), pointer :: fr_flig(:)  ! fine root litter lignin fraction
! Subprogram not used    integer , pointer :: npfts(:)    ! number of pfts for each column
! Subprogram not used    integer , pointer :: pfti(:)     ! beginning pft index for each column
! Subprogram not used    real(r8), pointer :: hrv_leafc_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_frootc_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_livestemc_to_litter(:)
! Subprogram not used    real(r8), pointer :: phrv_deadstemc_to_prod10c(:)
! Subprogram not used    real(r8), pointer :: phrv_deadstemc_to_prod100c(:)
! Subprogram not used    real(r8), pointer :: hrv_livecrootc_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_deadcrootc_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_leafc_storage_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_frootc_storage_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_livestemc_storage_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_deadstemc_storage_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_livecrootc_storage_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_deadcrootc_storage_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_gresp_storage_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_leafc_xfer_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_frootc_xfer_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_livestemc_xfer_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_deadstemc_xfer_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_livecrootc_xfer_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_deadcrootc_xfer_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_gresp_xfer_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_leafn_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_frootn_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_livestemn_to_litter(:)
! Subprogram not used    real(r8), pointer :: phrv_deadstemn_to_prod10n(:)
! Subprogram not used    real(r8), pointer :: phrv_deadstemn_to_prod100n(:)
! Subprogram not used    real(r8), pointer :: hrv_livecrootn_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_deadcrootn_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_retransn_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_leafn_storage_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_frootn_storage_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_livestemn_storage_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_deadstemn_storage_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_livecrootn_storage_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_deadcrootn_storage_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_leafn_xfer_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_frootn_xfer_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_livestemn_xfer_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_deadstemn_xfer_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_livecrootn_xfer_to_litter(:)
! Subprogram not used    real(r8), pointer :: hrv_deadcrootn_xfer_to_litter(:)
! Subprogram not used !
! Subprogram not used ! local pointers to implicit in/out arrays
! Subprogram not used    real(r8), pointer :: hrv_leafc_to_litr1c(:)
! Subprogram not used    real(r8), pointer :: hrv_leafc_to_litr2c(:)
! Subprogram not used    real(r8), pointer :: hrv_leafc_to_litr3c(:)
! Subprogram not used    real(r8), pointer :: hrv_frootc_to_litr1c(:)
! Subprogram not used    real(r8), pointer :: hrv_frootc_to_litr2c(:)
! Subprogram not used    real(r8), pointer :: hrv_frootc_to_litr3c(:)
! Subprogram not used    real(r8), pointer :: hrv_livestemc_to_cwdc(:)
! Subprogram not used    real(r8), pointer :: chrv_deadstemc_to_prod10c(:)
! Subprogram not used    real(r8), pointer :: chrv_deadstemc_to_prod100c(:)
! Subprogram not used    real(r8), pointer :: hrv_livecrootc_to_cwdc(:)
! Subprogram not used    real(r8), pointer :: hrv_deadcrootc_to_cwdc(:)
! Subprogram not used    real(r8), pointer :: hrv_leafc_storage_to_litr1c(:)
! Subprogram not used    real(r8), pointer :: hrv_frootc_storage_to_litr1c(:)
! Subprogram not used    real(r8), pointer :: hrv_livestemc_storage_to_litr1c(:)
! Subprogram not used    real(r8), pointer :: hrv_deadstemc_storage_to_litr1c(:)
! Subprogram not used    real(r8), pointer :: hrv_livecrootc_storage_to_litr1c(:)
! Subprogram not used    real(r8), pointer :: hrv_deadcrootc_storage_to_litr1c(:)
! Subprogram not used    real(r8), pointer :: hrv_gresp_storage_to_litr1c(:)
! Subprogram not used    real(r8), pointer :: hrv_leafc_xfer_to_litr1c(:)
! Subprogram not used    real(r8), pointer :: hrv_frootc_xfer_to_litr1c(:)
! Subprogram not used    real(r8), pointer :: hrv_livestemc_xfer_to_litr1c(:)
! Subprogram not used    real(r8), pointer :: hrv_deadstemc_xfer_to_litr1c(:)
! Subprogram not used    real(r8), pointer :: hrv_livecrootc_xfer_to_litr1c(:)
! Subprogram not used    real(r8), pointer :: hrv_deadcrootc_xfer_to_litr1c(:)
! Subprogram not used    real(r8), pointer :: hrv_gresp_xfer_to_litr1c(:)
! Subprogram not used    real(r8), pointer :: hrv_leafn_to_litr1n(:)
! Subprogram not used    real(r8), pointer :: hrv_leafn_to_litr2n(:)
! Subprogram not used    real(r8), pointer :: hrv_leafn_to_litr3n(:)
! Subprogram not used    real(r8), pointer :: hrv_frootn_to_litr1n(:)
! Subprogram not used    real(r8), pointer :: hrv_frootn_to_litr2n(:)
! Subprogram not used    real(r8), pointer :: hrv_frootn_to_litr3n(:)
! Subprogram not used    real(r8), pointer :: hrv_livestemn_to_cwdn(:)
! Subprogram not used    real(r8), pointer :: chrv_deadstemn_to_prod10n(:)
! Subprogram not used    real(r8), pointer :: chrv_deadstemn_to_prod100n(:)
! Subprogram not used    real(r8), pointer :: hrv_livecrootn_to_cwdn(:)
! Subprogram not used    real(r8), pointer :: hrv_deadcrootn_to_cwdn(:)
! Subprogram not used    real(r8), pointer :: hrv_retransn_to_litr1n(:)
! Subprogram not used    real(r8), pointer :: hrv_leafn_storage_to_litr1n(:)
! Subprogram not used    real(r8), pointer :: hrv_frootn_storage_to_litr1n(:)
! Subprogram not used    real(r8), pointer :: hrv_livestemn_storage_to_litr1n(:)
! Subprogram not used    real(r8), pointer :: hrv_deadstemn_storage_to_litr1n(:)
! Subprogram not used    real(r8), pointer :: hrv_livecrootn_storage_to_litr1n(:)
! Subprogram not used    real(r8), pointer :: hrv_deadcrootn_storage_to_litr1n(:)
! Subprogram not used    real(r8), pointer :: hrv_leafn_xfer_to_litr1n(:)
! Subprogram not used    real(r8), pointer :: hrv_frootn_xfer_to_litr1n(:)
! Subprogram not used    real(r8), pointer :: hrv_livestemn_xfer_to_litr1n(:)
! Subprogram not used    real(r8), pointer :: hrv_deadstemn_xfer_to_litr1n(:)
! Subprogram not used    real(r8), pointer :: hrv_livecrootn_xfer_to_litr1n(:)
! Subprogram not used    real(r8), pointer :: hrv_deadcrootn_xfer_to_litr1n(:)
! Subprogram not used !
! Subprogram not used ! local pointers to implicit out arrays
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! !OTHER LOCAL VARIABLES:
! Subprogram not used    integer :: fc,c,pi,p               ! indices
! Subprogram not used !EOP
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    ! assign local pointers
! Subprogram not used    lf_flab                        => pftcon%lf_flab
! Subprogram not used    lf_fcel                        => pftcon%lf_fcel
! Subprogram not used    lf_flig                        => pftcon%lf_flig
! Subprogram not used    fr_flab                        => pftcon%fr_flab
! Subprogram not used    fr_fcel                        => pftcon%fr_fcel
! Subprogram not used    fr_flig                        => pftcon%fr_flig
! Subprogram not used 
! Subprogram not used    ! assign local pointers to column-level arrays
! Subprogram not used    npfts                          => col%npfts
! Subprogram not used    pfti                           => col%pfti
! Subprogram not used    hrv_leafc_to_litr1c              => ccf%hrv_leafc_to_litr1c
! Subprogram not used    hrv_leafc_to_litr2c              => ccf%hrv_leafc_to_litr2c
! Subprogram not used    hrv_leafc_to_litr3c              => ccf%hrv_leafc_to_litr3c
! Subprogram not used    hrv_frootc_to_litr1c             => ccf%hrv_frootc_to_litr1c
! Subprogram not used    hrv_frootc_to_litr2c             => ccf%hrv_frootc_to_litr2c
! Subprogram not used    hrv_frootc_to_litr3c             => ccf%hrv_frootc_to_litr3c
! Subprogram not used    hrv_livestemc_to_cwdc            => ccf%hrv_livestemc_to_cwdc
! Subprogram not used    chrv_deadstemc_to_prod10c        => ccf%hrv_deadstemc_to_prod10c
! Subprogram not used    chrv_deadstemc_to_prod100c       => ccf%hrv_deadstemc_to_prod100c
! Subprogram not used    hrv_livecrootc_to_cwdc           => ccf%hrv_livecrootc_to_cwdc
! Subprogram not used    hrv_deadcrootc_to_cwdc           => ccf%hrv_deadcrootc_to_cwdc
! Subprogram not used    hrv_leafc_storage_to_litr1c      => ccf%hrv_leafc_storage_to_litr1c
! Subprogram not used    hrv_frootc_storage_to_litr1c     => ccf%hrv_frootc_storage_to_litr1c
! Subprogram not used    hrv_livestemc_storage_to_litr1c  => ccf%hrv_livestemc_storage_to_litr1c
! Subprogram not used    hrv_deadstemc_storage_to_litr1c  => ccf%hrv_deadstemc_storage_to_litr1c
! Subprogram not used    hrv_livecrootc_storage_to_litr1c => ccf%hrv_livecrootc_storage_to_litr1c
! Subprogram not used    hrv_deadcrootc_storage_to_litr1c => ccf%hrv_deadcrootc_storage_to_litr1c
! Subprogram not used    hrv_gresp_storage_to_litr1c      => ccf%hrv_gresp_storage_to_litr1c
! Subprogram not used    hrv_leafc_xfer_to_litr1c         => ccf%hrv_leafc_xfer_to_litr1c
! Subprogram not used    hrv_frootc_xfer_to_litr1c        => ccf%hrv_frootc_xfer_to_litr1c
! Subprogram not used    hrv_livestemc_xfer_to_litr1c     => ccf%hrv_livestemc_xfer_to_litr1c
! Subprogram not used    hrv_deadstemc_xfer_to_litr1c     => ccf%hrv_deadstemc_xfer_to_litr1c
! Subprogram not used    hrv_livecrootc_xfer_to_litr1c    => ccf%hrv_livecrootc_xfer_to_litr1c
! Subprogram not used    hrv_deadcrootc_xfer_to_litr1c    => ccf%hrv_deadcrootc_xfer_to_litr1c
! Subprogram not used    hrv_gresp_xfer_to_litr1c         => ccf%hrv_gresp_xfer_to_litr1c
! Subprogram not used    hrv_leafn_to_litr1n              => cnf%hrv_leafn_to_litr1n
! Subprogram not used    hrv_leafn_to_litr2n              => cnf%hrv_leafn_to_litr2n
! Subprogram not used    hrv_leafn_to_litr3n              => cnf%hrv_leafn_to_litr3n
! Subprogram not used    hrv_frootn_to_litr1n             => cnf%hrv_frootn_to_litr1n
! Subprogram not used    hrv_frootn_to_litr2n             => cnf%hrv_frootn_to_litr2n
! Subprogram not used    hrv_frootn_to_litr3n             => cnf%hrv_frootn_to_litr3n
! Subprogram not used    hrv_livestemn_to_cwdn            => cnf%hrv_livestemn_to_cwdn
! Subprogram not used    chrv_deadstemn_to_prod10n        => cnf%hrv_deadstemn_to_prod10n
! Subprogram not used    chrv_deadstemn_to_prod100n       => cnf%hrv_deadstemn_to_prod100n
! Subprogram not used    hrv_livecrootn_to_cwdn           => cnf%hrv_livecrootn_to_cwdn
! Subprogram not used    hrv_deadcrootn_to_cwdn           => cnf%hrv_deadcrootn_to_cwdn
! Subprogram not used    hrv_retransn_to_litr1n           => cnf%hrv_retransn_to_litr1n
! Subprogram not used    hrv_leafn_storage_to_litr1n      => cnf%hrv_leafn_storage_to_litr1n
! Subprogram not used    hrv_frootn_storage_to_litr1n     => cnf%hrv_frootn_storage_to_litr1n
! Subprogram not used    hrv_livestemn_storage_to_litr1n  => cnf%hrv_livestemn_storage_to_litr1n
! Subprogram not used    hrv_deadstemn_storage_to_litr1n  => cnf%hrv_deadstemn_storage_to_litr1n
! Subprogram not used    hrv_livecrootn_storage_to_litr1n => cnf%hrv_livecrootn_storage_to_litr1n
! Subprogram not used    hrv_deadcrootn_storage_to_litr1n => cnf%hrv_deadcrootn_storage_to_litr1n
! Subprogram not used    hrv_leafn_xfer_to_litr1n         => cnf%hrv_leafn_xfer_to_litr1n
! Subprogram not used    hrv_frootn_xfer_to_litr1n        => cnf%hrv_frootn_xfer_to_litr1n
! Subprogram not used    hrv_livestemn_xfer_to_litr1n     => cnf%hrv_livestemn_xfer_to_litr1n
! Subprogram not used    hrv_deadstemn_xfer_to_litr1n     => cnf%hrv_deadstemn_xfer_to_litr1n
! Subprogram not used    hrv_livecrootn_xfer_to_litr1n    => cnf%hrv_livecrootn_xfer_to_litr1n
! Subprogram not used    hrv_deadcrootn_xfer_to_litr1n    => cnf%hrv_deadcrootn_xfer_to_litr1n
! Subprogram not used 
! Subprogram not used    ! assign local pointers to pft-level arrays
! Subprogram not used    ivt                            => pft%itype
! Subprogram not used    wtcol                          => pft%wtcol
! Subprogram not used    pwtgcell                       => pft%wtgcell  
! Subprogram not used    hrv_leafc_to_litter              => pcf%hrv_leafc_to_litter
! Subprogram not used    hrv_frootc_to_litter             => pcf%hrv_frootc_to_litter
! Subprogram not used    hrv_livestemc_to_litter          => pcf%hrv_livestemc_to_litter
! Subprogram not used    phrv_deadstemc_to_prod10c        => pcf%hrv_deadstemc_to_prod10c
! Subprogram not used    phrv_deadstemc_to_prod100c       => pcf%hrv_deadstemc_to_prod100c
! Subprogram not used    hrv_livecrootc_to_litter         => pcf%hrv_livecrootc_to_litter
! Subprogram not used    hrv_deadcrootc_to_litter         => pcf%hrv_deadcrootc_to_litter
! Subprogram not used    hrv_leafc_storage_to_litter      => pcf%hrv_leafc_storage_to_litter
! Subprogram not used    hrv_frootc_storage_to_litter     => pcf%hrv_frootc_storage_to_litter
! Subprogram not used    hrv_livestemc_storage_to_litter  => pcf%hrv_livestemc_storage_to_litter
! Subprogram not used    hrv_deadstemc_storage_to_litter  => pcf%hrv_deadstemc_storage_to_litter
! Subprogram not used    hrv_livecrootc_storage_to_litter => pcf%hrv_livecrootc_storage_to_litter
! Subprogram not used    hrv_deadcrootc_storage_to_litter => pcf%hrv_deadcrootc_storage_to_litter
! Subprogram not used    hrv_gresp_storage_to_litter      => pcf%hrv_gresp_storage_to_litter
! Subprogram not used    hrv_leafc_xfer_to_litter         => pcf%hrv_leafc_xfer_to_litter
! Subprogram not used    hrv_frootc_xfer_to_litter        => pcf%hrv_frootc_xfer_to_litter
! Subprogram not used    hrv_livestemc_xfer_to_litter     => pcf%hrv_livestemc_xfer_to_litter
! Subprogram not used    hrv_deadstemc_xfer_to_litter     => pcf%hrv_deadstemc_xfer_to_litter
! Subprogram not used    hrv_livecrootc_xfer_to_litter    => pcf%hrv_livecrootc_xfer_to_litter
! Subprogram not used    hrv_deadcrootc_xfer_to_litter    => pcf%hrv_deadcrootc_xfer_to_litter
! Subprogram not used    hrv_gresp_xfer_to_litter         => pcf%hrv_gresp_xfer_to_litter
! Subprogram not used    hrv_leafn_to_litter              => pnf%hrv_leafn_to_litter
! Subprogram not used    hrv_frootn_to_litter             => pnf%hrv_frootn_to_litter
! Subprogram not used    hrv_livestemn_to_litter          => pnf%hrv_livestemn_to_litter
! Subprogram not used    phrv_deadstemn_to_prod10n        => pnf%hrv_deadstemn_to_prod10n
! Subprogram not used    phrv_deadstemn_to_prod100n       => pnf%hrv_deadstemn_to_prod100n
! Subprogram not used    hrv_livecrootn_to_litter         => pnf%hrv_livecrootn_to_litter
! Subprogram not used    hrv_deadcrootn_to_litter         => pnf%hrv_deadcrootn_to_litter
! Subprogram not used    hrv_retransn_to_litter           => pnf%hrv_retransn_to_litter
! Subprogram not used    hrv_leafn_storage_to_litter      => pnf%hrv_leafn_storage_to_litter
! Subprogram not used    hrv_frootn_storage_to_litter     => pnf%hrv_frootn_storage_to_litter
! Subprogram not used    hrv_livestemn_storage_to_litter  => pnf%hrv_livestemn_storage_to_litter
! Subprogram not used    hrv_deadstemn_storage_to_litter  => pnf%hrv_deadstemn_storage_to_litter
! Subprogram not used    hrv_livecrootn_storage_to_litter => pnf%hrv_livecrootn_storage_to_litter
! Subprogram not used    hrv_deadcrootn_storage_to_litter => pnf%hrv_deadcrootn_storage_to_litter
! Subprogram not used    hrv_leafn_xfer_to_litter         => pnf%hrv_leafn_xfer_to_litter
! Subprogram not used    hrv_frootn_xfer_to_litter        => pnf%hrv_frootn_xfer_to_litter
! Subprogram not used    hrv_livestemn_xfer_to_litter     => pnf%hrv_livestemn_xfer_to_litter
! Subprogram not used    hrv_deadstemn_xfer_to_litter     => pnf%hrv_deadstemn_xfer_to_litter
! Subprogram not used    hrv_livecrootn_xfer_to_litter    => pnf%hrv_livecrootn_xfer_to_litter
! Subprogram not used    hrv_deadcrootn_xfer_to_litter    => pnf%hrv_deadcrootn_xfer_to_litter
! Subprogram not used 
! Subprogram not used    do pi = 1,maxpatch_pft
! Subprogram not used       do fc = 1,num_soilc
! Subprogram not used          c = filter_soilc(fc)
! Subprogram not used 
! Subprogram not used          if (pi <=  npfts(c)) then
! Subprogram not used             p = pfti(c) + pi - 1
! Subprogram not used 
! Subprogram not used             if (pwtgcell(p)>0._r8) then
! Subprogram not used 
! Subprogram not used                ! leaf harvest mortality carbon fluxes
! Subprogram not used                hrv_leafc_to_litr1c(c) = hrv_leafc_to_litr1c(c) + &
! Subprogram not used                   hrv_leafc_to_litter(p) * lf_flab(ivt(p)) * wtcol(p)
! Subprogram not used                hrv_leafc_to_litr2c(c) = hrv_leafc_to_litr2c(c) + &
! Subprogram not used                   hrv_leafc_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p)
! Subprogram not used                hrv_leafc_to_litr3c(c) = hrv_leafc_to_litr3c(c) + &
! Subprogram not used                   hrv_leafc_to_litter(p) * lf_flig(ivt(p)) * wtcol(p)
! Subprogram not used 
! Subprogram not used                ! fine root harvest mortality carbon fluxes
! Subprogram not used                hrv_frootc_to_litr1c(c) = hrv_frootc_to_litr1c(c) + &
! Subprogram not used                   hrv_frootc_to_litter(p) * fr_flab(ivt(p)) * wtcol(p)
! Subprogram not used                hrv_frootc_to_litr2c(c) = hrv_frootc_to_litr2c(c) + &
! Subprogram not used                   hrv_frootc_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p)
! Subprogram not used                hrv_frootc_to_litr3c(c) = hrv_frootc_to_litr3c(c) + &
! Subprogram not used                   hrv_frootc_to_litter(p) * fr_flig(ivt(p)) * wtcol(p)
! Subprogram not used 
! Subprogram not used                ! wood harvest mortality carbon fluxes
! Subprogram not used                hrv_livestemc_to_cwdc(c)  = hrv_livestemc_to_cwdc(c)  + &
! Subprogram not used                   hrv_livestemc_to_litter(p)  * wtcol(p)
! Subprogram not used                chrv_deadstemc_to_prod10c(c)  = chrv_deadstemc_to_prod10c(c)  + &
! Subprogram not used                   phrv_deadstemc_to_prod10c(p)  * wtcol(p)
! Subprogram not used                chrv_deadstemc_to_prod100c(c)  = chrv_deadstemc_to_prod100c(c)  + &
! Subprogram not used                   phrv_deadstemc_to_prod100c(p)  * wtcol(p)
! Subprogram not used                hrv_livecrootc_to_cwdc(c) = hrv_livecrootc_to_cwdc(c) + &
! Subprogram not used                   hrv_livecrootc_to_litter(p) * wtcol(p)
! Subprogram not used                hrv_deadcrootc_to_cwdc(c) = hrv_deadcrootc_to_cwdc(c) + &
! Subprogram not used                   hrv_deadcrootc_to_litter(p) * wtcol(p)
! Subprogram not used 
! Subprogram not used                ! storage harvest mortality carbon fluxes
! Subprogram not used                hrv_leafc_storage_to_litr1c(c)      = hrv_leafc_storage_to_litr1c(c)      + &
! Subprogram not used                   hrv_leafc_storage_to_litter(p)      * wtcol(p)
! Subprogram not used                hrv_frootc_storage_to_litr1c(c)     = hrv_frootc_storage_to_litr1c(c)     + &
! Subprogram not used                   hrv_frootc_storage_to_litter(p)     * wtcol(p)
! Subprogram not used                hrv_livestemc_storage_to_litr1c(c)  = hrv_livestemc_storage_to_litr1c(c)  + &
! Subprogram not used                   hrv_livestemc_storage_to_litter(p)  * wtcol(p)
! Subprogram not used                hrv_deadstemc_storage_to_litr1c(c)  = hrv_deadstemc_storage_to_litr1c(c)  + &
! Subprogram not used                   hrv_deadstemc_storage_to_litter(p)  * wtcol(p)
! Subprogram not used                hrv_livecrootc_storage_to_litr1c(c) = hrv_livecrootc_storage_to_litr1c(c) + &
! Subprogram not used                   hrv_livecrootc_storage_to_litter(p) * wtcol(p)
! Subprogram not used                hrv_deadcrootc_storage_to_litr1c(c) = hrv_deadcrootc_storage_to_litr1c(c) + &
! Subprogram not used                   hrv_deadcrootc_storage_to_litter(p) * wtcol(p)
! Subprogram not used                hrv_gresp_storage_to_litr1c(c)      = hrv_gresp_storage_to_litr1c(c)      + &
! Subprogram not used                   hrv_gresp_storage_to_litter(p)      * wtcol(p)
! Subprogram not used 
! Subprogram not used                ! transfer harvest mortality carbon fluxes
! Subprogram not used                hrv_leafc_xfer_to_litr1c(c)      = hrv_leafc_xfer_to_litr1c(c)      + &
! Subprogram not used                   hrv_leafc_xfer_to_litter(p)      * wtcol(p)
! Subprogram not used                hrv_frootc_xfer_to_litr1c(c)     = hrv_frootc_xfer_to_litr1c(c)     + &
! Subprogram not used                   hrv_frootc_xfer_to_litter(p)     * wtcol(p)
! Subprogram not used                hrv_livestemc_xfer_to_litr1c(c)  = hrv_livestemc_xfer_to_litr1c(c)  + &
! Subprogram not used                   hrv_livestemc_xfer_to_litter(p)  * wtcol(p)
! Subprogram not used                hrv_deadstemc_xfer_to_litr1c(c)  = hrv_deadstemc_xfer_to_litr1c(c)  + &
! Subprogram not used                   hrv_deadstemc_xfer_to_litter(p)  * wtcol(p)
! Subprogram not used                hrv_livecrootc_xfer_to_litr1c(c) = hrv_livecrootc_xfer_to_litr1c(c) + &
! Subprogram not used                   hrv_livecrootc_xfer_to_litter(p) * wtcol(p)
! Subprogram not used                hrv_deadcrootc_xfer_to_litr1c(c) = hrv_deadcrootc_xfer_to_litr1c(c) + &
! Subprogram not used                   hrv_deadcrootc_xfer_to_litter(p) * wtcol(p)
! Subprogram not used                hrv_gresp_xfer_to_litr1c(c)      = hrv_gresp_xfer_to_litr1c(c)      + &
! Subprogram not used                   hrv_gresp_xfer_to_litter(p)      * wtcol(p)
! Subprogram not used 
! Subprogram not used                ! leaf harvest mortality nitrogen fluxes
! Subprogram not used                hrv_leafn_to_litr1n(c) = hrv_leafn_to_litr1n(c) + &
! Subprogram not used                   hrv_leafn_to_litter(p) * lf_flab(ivt(p)) * wtcol(p)
! Subprogram not used                hrv_leafn_to_litr2n(c) = hrv_leafn_to_litr2n(c) + &
! Subprogram not used                   hrv_leafn_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p)
! Subprogram not used                hrv_leafn_to_litr3n(c) = hrv_leafn_to_litr3n(c) + &
! Subprogram not used                   hrv_leafn_to_litter(p) * lf_flig(ivt(p)) * wtcol(p)
! Subprogram not used 
! Subprogram not used                ! fine root litter nitrogen fluxes
! Subprogram not used                hrv_frootn_to_litr1n(c) = hrv_frootn_to_litr1n(c) + &
! Subprogram not used                   hrv_frootn_to_litter(p) * fr_flab(ivt(p)) * wtcol(p)
! Subprogram not used                hrv_frootn_to_litr2n(c) = hrv_frootn_to_litr2n(c) + &
! Subprogram not used                   hrv_frootn_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p)
! Subprogram not used                hrv_frootn_to_litr3n(c) = hrv_frootn_to_litr3n(c) + &
! Subprogram not used                   hrv_frootn_to_litter(p) * fr_flig(ivt(p)) * wtcol(p)
! Subprogram not used 
! Subprogram not used                ! wood harvest mortality nitrogen fluxes
! Subprogram not used                hrv_livestemn_to_cwdn(c)  = hrv_livestemn_to_cwdn(c)  + &
! Subprogram not used                   hrv_livestemn_to_litter(p)  * wtcol(p)
! Subprogram not used                chrv_deadstemn_to_prod10n(c)  = chrv_deadstemn_to_prod10n(c)  + &
! Subprogram not used                   phrv_deadstemn_to_prod10n(p)  * wtcol(p)
! Subprogram not used                chrv_deadstemn_to_prod100n(c)  = chrv_deadstemn_to_prod100n(c)  + &
! Subprogram not used                   phrv_deadstemn_to_prod100n(p)  * wtcol(p)
! Subprogram not used                hrv_livecrootn_to_cwdn(c) = hrv_livecrootn_to_cwdn(c) + &
! Subprogram not used                   hrv_livecrootn_to_litter(p) * wtcol(p)
! Subprogram not used                hrv_deadcrootn_to_cwdn(c) = hrv_deadcrootn_to_cwdn(c) + &
! Subprogram not used                   hrv_deadcrootn_to_litter(p) * wtcol(p)
! Subprogram not used 
! Subprogram not used                ! retranslocated N pool harvest mortality fluxes
! Subprogram not used                hrv_retransn_to_litr1n(c) = hrv_retransn_to_litr1n(c) + &
! Subprogram not used                   hrv_retransn_to_litter(p) * wtcol(p)
! Subprogram not used 
! Subprogram not used                ! storage harvest mortality nitrogen fluxes
! Subprogram not used                hrv_leafn_storage_to_litr1n(c)      = hrv_leafn_storage_to_litr1n(c)      + &
! Subprogram not used                   hrv_leafn_storage_to_litter(p)      * wtcol(p)
! Subprogram not used                hrv_frootn_storage_to_litr1n(c)     = hrv_frootn_storage_to_litr1n(c)     + &
! Subprogram not used                   hrv_frootn_storage_to_litter(p)     * wtcol(p)
! Subprogram not used                hrv_livestemn_storage_to_litr1n(c)  = hrv_livestemn_storage_to_litr1n(c)  + &
! Subprogram not used                   hrv_livestemn_storage_to_litter(p)  * wtcol(p)
! Subprogram not used                hrv_deadstemn_storage_to_litr1n(c)  = hrv_deadstemn_storage_to_litr1n(c)  + &
! Subprogram not used                   hrv_deadstemn_storage_to_litter(p)  * wtcol(p)
! Subprogram not used                hrv_livecrootn_storage_to_litr1n(c) = hrv_livecrootn_storage_to_litr1n(c) + &
! Subprogram not used                   hrv_livecrootn_storage_to_litter(p) * wtcol(p)
! Subprogram not used                hrv_deadcrootn_storage_to_litr1n(c) = hrv_deadcrootn_storage_to_litr1n(c) + &
! Subprogram not used                   hrv_deadcrootn_storage_to_litter(p) * wtcol(p)
! Subprogram not used 
! Subprogram not used                ! transfer harvest mortality nitrogen fluxes
! Subprogram not used                hrv_leafn_xfer_to_litr1n(c)      = hrv_leafn_xfer_to_litr1n(c)      + &
! Subprogram not used                   hrv_leafn_xfer_to_litter(p)      * wtcol(p)
! Subprogram not used                hrv_frootn_xfer_to_litr1n(c)     = hrv_frootn_xfer_to_litr1n(c)     + &
! Subprogram not used                   hrv_frootn_xfer_to_litter(p)     * wtcol(p)
! Subprogram not used                hrv_livestemn_xfer_to_litr1n(c)  = hrv_livestemn_xfer_to_litr1n(c)  + &
! Subprogram not used                   hrv_livestemn_xfer_to_litter(p)  * wtcol(p)
! Subprogram not used                hrv_deadstemn_xfer_to_litr1n(c)  = hrv_deadstemn_xfer_to_litr1n(c)  + &
! Subprogram not used                   hrv_deadstemn_xfer_to_litter(p)  * wtcol(p)
! Subprogram not used                hrv_livecrootn_xfer_to_litr1n(c) = hrv_livecrootn_xfer_to_litr1n(c) + &
! Subprogram not used                   hrv_livecrootn_xfer_to_litter(p) * wtcol(p)
! Subprogram not used                hrv_deadcrootn_xfer_to_litr1n(c) = hrv_deadcrootn_xfer_to_litr1n(c) + &
! Subprogram not used                   hrv_deadcrootn_xfer_to_litter(p) * wtcol(p)
! Subprogram not used 
! Subprogram not used             end if
! Subprogram not used          end if
! Subprogram not used 
! Subprogram not used       end do
! Subprogram not used 
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used end subroutine CNHarvestPftToColumn
!-----------------------------------------------------------------------

end module pftdynMod
