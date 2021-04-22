module CropRestMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: CropRestMod
! 
! !DESCRIPTION: 
! Read/Write to/from Crop info to CLM restart file. 
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use spmdMod     , only : masterproc
  use abortutils  , only : endrun
!
! !PUBLIC TYPES:
  implicit none
  private
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: CropRest        ! Restart prognostic crop model
  public :: CropRestYear    ! Get the number of years crop has spunup
  public :: CropRestIncYear ! Increment the crop spinup years
!
! !REVISION HISTORY:
! Module created by slevis following CNRestMod by Peter Thornton
!

! !PRIVATE DATA MEMBERS:
  integer :: restyear = 0         ! Restart year from the initial conditions file, incremented as time elapses

!EOP
!----------------------------------------------------------------------- 

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CropRest
!
! !INTERFACE:
! Subprogram not used   subroutine CropRest ( ncid, flag )
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION: 
! Subprogram not used ! Read/write Crop restart data
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used     use clmtype
! Subprogram not used     use clm_atmlnd      , only : clm_a2l
! Subprogram not used     use clm_varpar      , only : numrad
! Subprogram not used     use decompMod       , only : get_proc_bounds
! Subprogram not used     use clm_time_manager, only : is_restart
! Subprogram not used     use ncdio_pio
! Subprogram not used !
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     type(file_desc_t)  :: ncid             ! netcdf id
! Subprogram not used     character(len=*), intent(in) :: flag   !'read' or 'write'
! Subprogram not used !
! Subprogram not used ! !CALLED FROM:
! Subprogram not used ! subroutine restart in module restFileMod
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! Author: slevis
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used ! !LOCAL VARIABLES:
! Subprogram not used     integer :: c,p,j                      ! indices 
! Subprogram not used     integer :: begp, endp                 ! per-proc beginning and ending pft indices
! Subprogram not used     integer :: begc, endc                 ! per-proc beginning and ending column indices 
! Subprogram not used     integer :: begl, endl                 ! per-proc beginning and ending landunit indices
! Subprogram not used     integer :: begg, endg                 ! per-proc gridcell ending gridcell indices
! Subprogram not used     real(r8):: m                          ! multiplier for the exit_spinup code
! Subprogram not used     logical :: readvar                    ! determine if variable is on initial file
! Subprogram not used     character(len=128) :: varname         ! temporary
! Subprogram not used     type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived subtype
! Subprogram not used     type(landunit_type), pointer :: lptr  ! pointer to landunit derived subtype
! Subprogram not used     type(column_type)  , pointer :: cptr  ! pointer to column derived subtype
! Subprogram not used     type(pft_type)     , pointer :: pptr  ! pointer to pft derived subtype
! Subprogram not used     integer , pointer :: iptemp(:)        ! pointer to memory to be allocated
! Subprogram not used     integer :: ier                        ! error status
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     ! Prognostic crop restart year
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='restyear', xtype=ncd_int,  &
! Subprogram not used             long_name='Number of years prognostic crop ran', units="years")
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='restyear', data=restyear, &
! Subprogram not used             ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' )then
! Subprogram not used            if ( readvar ) then
! Subprogram not used               call checkDates( )
! Subprogram not used            else
! Subprogram not used               if ( is_restart()) call endrun
! Subprogram not used            end if
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! Set pointers into derived type
! Subprogram not used 
! Subprogram not used     gptr => grc
! Subprogram not used     lptr => lun
! Subprogram not used     cptr => col
! Subprogram not used     pptr => pft
! Subprogram not used 
! Subprogram not used     ! Determine necessary subgrid bounds
! Subprogram not used 
! Subprogram not used     call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
! Subprogram not used 
! Subprogram not used     !--------------------------------
! Subprogram not used     ! pft physical state variables 
! Subprogram not used     !--------------------------------
! Subprogram not used 
! Subprogram not used     ! peaklai
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='peaklai', xtype=ncd_int,  &
! Subprogram not used             dim1name='pft',long_name='Flag if at max allowed LAI or not', &
! Subprogram not used             flag_values=(/0,1/), nvalid_range=(/0,1/),                    &
! Subprogram not used             flag_meanings=(/'NOT-at-peak', 'AT_peak-LAI' /) )
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='peaklai', data=pps%peaklai, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! idop
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='idop', xtype=ncd_int,  &
! Subprogram not used             dim1name='pft',long_name='Date of planting',units='jday', &
! Subprogram not used             nvalid_range=(/1,366/) )
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='idop', data=pps%idop, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! aleaf
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='aleaf', xtype=ncd_double,  &
! Subprogram not used             dim1name='pft',long_name='leaf allocation coefficient',units='')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='aleaf', data=pps%aleaf, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! aleafi
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='aleafi', xtype=ncd_double,  &
! Subprogram not used             dim1name='pft',long_name='Saved leaf allocation coefficient from phase 2', &
! Subprogram not used             units='')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='aleafi', data=pps%aleafi, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! astem
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='astem', xtype=ncd_double,  &
! Subprogram not used             dim1name='pft',long_name='stem allocation coefficient',units='')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='astem', data=pps%astem, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! astemi
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='astemi', xtype=ncd_double,  &
! Subprogram not used             dim1name='pft',long_name='Saved stem allocation coefficient from phase 2',&
! Subprogram not used             units='')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='astemi', data=pps%astemi, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! htmx 
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='htmx', xtype=ncd_double,  &
! Subprogram not used             dim1name='pft',long_name='max height attained by a crop during year',&
! Subprogram not used             units='m')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='htmx', data=pps%htmx, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! hdidx
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='hdidx', xtype=ncd_double,  &
! Subprogram not used             dim1name='pft',long_name='cold hardening index',units='')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='hdidx', data=pps%hdidx, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! vf
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='vf', xtype=ncd_double,  &
! Subprogram not used             dim1name='pft',long_name='vernalization factor',units='')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='vf', data=pps%vf, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! cumvd
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='cumvd', xtype=ncd_double,  &
! Subprogram not used             dim1name='pft',long_name='cumulative vernalization d',units='')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='cumvd', data=pps%cumvd, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! croplive
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='croplive', xtype=ncd_log,  &
! Subprogram not used             dim1name='pft',long_name='Flag that crop is alive, but not harvested')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='croplive', data=pps%croplive, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! cropplant
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='cropplant', xtype=ncd_log,  &
! Subprogram not used             dim1name='pft',long_name='Flag that crop is planted, but not harvested' )
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='cropplant', data=pps%cropplant, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! harvdate
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='harvdate', xtype=ncd_int,  &
! Subprogram not used             dim1name='pft',long_name='harvest date',units='jday', &
! Subprogram not used             nvalid_range=(/1,366/) )
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='harvdate', data=pps%harvdate, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! gdd1020
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='gdd1020', xtype=ncd_double,  &
! Subprogram not used             dim1name='pft', &
! Subprogram not used             long_name='20 year average of growing degree-days base 10C from planting', &
! Subprogram not used             units='ddays')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='gdd1020', data=pps%gdd1020, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! gdd820
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='gdd820', xtype=ncd_double,  &
! Subprogram not used             dim1name='pft', &
! Subprogram not used             long_name='20 year average of growing degree-days base 8C from planting', &
! Subprogram not used             units='ddays')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='gdd820', data=pps%gdd820, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! gdd020
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='gdd020', xtype=ncd_double,  &
! Subprogram not used             dim1name='pft', &
! Subprogram not used             long_name='20 year average of growing degree-days base 0C from planting', &
! Subprogram not used             units='ddays')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='gdd020', data=pps%gdd020, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! gddmaturity
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='gddmaturity', xtype=ncd_double,  &
! Subprogram not used             dim1name='pft',long_name='Growing degree days needed to harvest',units='ddays')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='gddmaturity', data=pps%gddmaturity, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! huileaf
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='huileaf', xtype=ncd_double,  &
! Subprogram not used             dim1name='pft', &
! Subprogram not used             long_name='heat unit index needed from planting to leaf emergence',units='')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='huileaf', data=pps%huileaf, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! huigrain
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='huigrain', xtype=ncd_double,  &
! Subprogram not used             dim1name='pft',long_name='heat unit index needed to reach vegetative maturity', &
! Subprogram not used             units='')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='huigrain', data=pps%huigrain, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! grainc
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='grainc', xtype=ncd_double,  &
! Subprogram not used             dim1name='pft',long_name='grain C',units='gC/m2')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='grainc', data=pcs%grainc, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! grainc_storage
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='grainc_storage', xtype=ncd_double,  &
! Subprogram not used             dim1name='pft',long_name='grain C storage',units='gC/m2')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='grainc_storage', data=pcs%grainc_storage, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! grainc_xfer
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='grainc_xfer', xtype=ncd_double,  &
! Subprogram not used             dim1name='pft',long_name='grain C transfer',units='gC/m2')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='grainc_xfer', data=pcs%grainc_xfer, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! grainn
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='grainn', xtype=ncd_double,  &
! Subprogram not used             dim1name='pft',long_name='grain N',units='gN/m2')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='grainn', data=pns%grainn, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! grainn_storage
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='grainn_storage', xtype=ncd_double,  &
! Subprogram not used             dim1name='pft',long_name='grain N storage',units='gN/m2')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='grainn_storage', data=pns%grainn_storage, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! grainn_xfer
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='grainn_xfer', xtype=ncd_double,  &
! Subprogram not used             dim1name='pft',long_name='grain N transfer',units='gN/m2')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='grainn_xfer', data=pns%grainn_xfer, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! grainc_xfer_to_grainc
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='grainc_xfer_to_grainc', xtype=ncd_double,  &
! Subprogram not used             dim1name='pft',long_name='grain C growth from storage',units='gC/m2/s')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='grainc_xfer_to_grainc', data=pcf%grainc_xfer_to_grainc, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! livestemc_to_litter
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='livestemc_to_litter', xtype=ncd_double,  &
! Subprogram not used             dim1name='pft',long_name='live stem C litterfall',units='gC/m2/s')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='livestemc_to_litter', data=pcf%livestemc_to_litter, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! grainc_to_food
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='grainc_to_food', xtype=ncd_double,  &
! Subprogram not used             dim1name='pft',long_name='grain C to food',units='gC/m2/s')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='grainc_to_food', data=pcf%grainc_to_food, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! grainn_xfer_to_grainn
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='grainn_xfer_to_grainn', xtype=ncd_double,  &
! Subprogram not used             dim1name='pft',long_name='grain N growth from storage',units='gN/m2/s')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='grainn_xfer_to_grainn', data=pnf%grainn_xfer_to_grainn, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! livestemn_to_litter
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='livestemn_to_litter', xtype=ncd_double,  &
! Subprogram not used             dim1name='pft',long_name='livestem N to litter',units='gN/m2/s')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='livestemn_to_litter', data=pnf%livestemn_to_litter, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! grainn_to_food
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='grainn_to_food', xtype=ncd_double,  &
! Subprogram not used             dim1name='pft',long_name='grain N to food',units='gN/m2/s')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='grainn_to_food', data=pnf%grainn_to_food, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! cpool_to_grainc
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='cpool_to_grainc', xtype=ncd_double,  &
! Subprogram not used             dim1name='pft',long_name='allocation to grain C',units='gC/m2/s')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='cpool_to_grainc', data=pcf%cpool_to_grainc, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! cpool_to_grainc_storage
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='cpool_to_grainc_storage', xtype=ncd_double,  &
! Subprogram not used             dim1name='pft',long_name='allocation to grain C storage',units='gC/m2/s')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='cpool_to_grainc_storage', data=pcf%cpool_to_grainc_storage, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! npool_to_grainn
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='npool_to_grainn', xtype=ncd_double,  &
! Subprogram not used             dim1name='pft',long_name='allocation to grain N',units='gN/m2/s')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='npool_to_grainn', data=pnf%npool_to_grainn, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! npool_to_grainn_storage
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='npool_to_grainn_storage', xtype=ncd_double,  &
! Subprogram not used             dim1name='pft',long_name='allocation to grain N storage',units='gN/m2/s')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='npool_to_grainn_storage', data=pnf%npool_to_grainn_storage, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! cpool_grain_gr
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='cpool_grain_gr', xtype=ncd_double,  &
! Subprogram not used             dim1name='pft',long_name='grain growth respiration',units='gC/m2/s')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='cpool_grain_gr', data=pcf%cpool_grain_gr, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! cpool_grain_storage_gr
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='cpool_grain_storage_gr', xtype=ncd_double,  &
! Subprogram not used             dim1name='pft',long_name='grain growth respiration to storage',units='gC/m2/s')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='cpool_grain_storage_gr', data=pcf%cpool_grain_storage_gr, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! transfer_grain_gr
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='transfer_grain_gr', xtype=ncd_double,  &
! Subprogram not used             dim1name='pft',long_name='grain growth respiration from storage',units='gC/m2/s')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='transfer_grain_gr', data=pcf%transfer_grain_gr, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! grainc_storage_to_xfer
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='grainc_storage_to_xfer', xtype=ncd_double,  &
! Subprogram not used             dim1name='pft',long_name='grain C shift storage to transfer',units='gC/m2/s')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='grainc_storage_to_xfer', data=pcf%grainc_storage_to_xfer, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! grainn_storage_to_xfer
! Subprogram not used     if (flag == 'define') then
! Subprogram not used        call ncd_defvar(ncid=ncid, varname='grainn_storage_to_xfer', xtype=ncd_double,  &
! Subprogram not used             dim1name='pft',long_name='grain N shift storage to transfer',units='gN/m2/s')
! Subprogram not used     else if (flag == 'read' .or. flag == 'write') then
! Subprogram not used        call ncd_io(varname='grainn_storage_to_xfer', data=pnf%grainn_storage_to_xfer, &
! Subprogram not used             dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
! Subprogram not used        if (flag=='read' .and. .not. readvar) then
! Subprogram not used            if (is_restart()) call endrun
! Subprogram not used        end if       
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used   end subroutine CropRest

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CropRestYear
!
! !INTERFACE:
! Subprogram not used   integer function CropRestYear ( )
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION: 
! Subprogram not used ! Return the restart year for prognostic crop
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! Author: Erik Kluzek
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used ! !LOCAL VARIABLES:
! Subprogram not used      CropRestYear = restyear
! Subprogram not used   end function CropRestYear

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CropRestIncYear
!
! !INTERFACE:
  subroutine CropRestIncYear ()
!
! !DESCRIPTION: 
! Increment the crop restart year, if appropriate
!
! This routine should be called every time step, but only once per clump (to avoid
! inadvertently updating nyrs multiple times)
!
! !USES:
    use surfrdMod        , only : crop_prog
    use clm_time_manager , only : get_curr_date, is_first_step
    implicit none
!
! !LOCAL VARIABLES:
    integer kyr                     ! current year
    integer kmo                     !         month of year  (1, ..., 12)
    integer kda                     !         day of month   (1, ..., 31)
    integer mcsec                   !         seconds of day (0, ..., seconds/day)
!-----------------------------------------------------------------------

    ! Update restyear only when running with prognostic crop
    if ( crop_prog )then
       ! Update restyear when it's the start of a new year - but don't do that at the
       ! very start of the run
       call get_curr_date (   kyr, kmo, kda, mcsec)
       if ((kmo == 1 .and. kda == 1 .and. mcsec == 0) .and. .not. is_first_step()) then
          restyear = restyear + 1
       end if
    end if

  end subroutine CropRestIncYear

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: checkDates
!
! !INTERFACE:
! Subprogram not used   subroutine checkDates( )
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION: 
! Subprogram not used ! Make sure the dates are compatible. The date given to startup the model
! Subprogram not used ! and the date on the restart file must be the same although years can be
! Subprogram not used ! different. The dates need to be checked when the restart file is being
! Subprogram not used ! read in for a startup or branch case (they are NOT allowed to be different
! Subprogram not used ! for a restart case).
! Subprogram not used !
! Subprogram not used ! For the prognostic crop model the date of planting is tracked and growing
! Subprogram not used ! degree days is tracked (with a 20 year mean) -- so shifting the start dates
! Subprogram not used ! messes up these bits of saved information.
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used     use clm_time_manager, only : get_driver_start_ymd, get_start_date
! Subprogram not used     use clm_varctl      , only : iulog
! Subprogram not used     use clm_varctl      , only : nsrest, nsrBranch, nsrStartup
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! Author: Erik Kluzek
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used ! !LOCAL VARIABLES:
! Subprogram not used     integer :: stymd       ! Start date YYYYMMDD from driver
! Subprogram not used     integer :: styr        ! Start year from driver
! Subprogram not used     integer :: stmon_day   ! Start date MMDD from driver
! Subprogram not used     integer :: rsmon_day   ! Restart date MMDD from restart file
! Subprogram not used     integer :: rsyr        ! Restart year from restart file
! Subprogram not used     integer :: rsmon       ! Restart month from restart file
! Subprogram not used     integer :: rsday       ! Restart day from restart file
! Subprogram not used     integer :: tod         ! Restart time of day from restart file
! Subprogram not used     character(len=*), parameter :: formDate = '(A,i4.4,"/",i2.2,"/",i2.2)' ! log output format
! Subprogram not used     character(len=32) :: subname = 'CropRest::checkDates'
! Subprogram not used     !
! Subprogram not used     ! If branch or startup make sure the startdate is compatible with the date
! Subprogram not used     ! on the restart file.
! Subprogram not used     !
! Subprogram not used     if ( nsrest == nsrBranch .or. nsrest == nsrStartup )then
! Subprogram not used        stymd       = get_driver_start_ymd()
! Subprogram not used        styr        = stymd / 10000
! Subprogram not used        stmon_day   = stymd - styr*10000
! Subprogram not used        call get_start_date( rsyr, rsmon, rsday, tod )
! Subprogram not used        rsmon_day = rsmon*100 + rsday
! Subprogram not used        if ( masterproc ) &
! Subprogram not used        write(iulog,formDate) 'Date on the restart file is: ', rsyr, rsmon, rsday
! Subprogram not used        if ( stmon_day /= rsmon_day )then
! Subprogram not used           write(iulog,formDate) 'Start date is: ', styr, stmon_day/100, &
! Subprogram not used                                  (stmon_day - stmon_day/100)
! Subprogram not used           call endrun( trim(subname)// &
! Subprogram not used           ' ERROR: For prognostic crop to work correctly, the start date (month and day)'// &
! Subprogram not used           ' and the date on the restart file needs to match (years can be different)' )
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used   end subroutine checkDates

end module CropRestMod

