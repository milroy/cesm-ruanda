!-------------------------------------------------------------------
! manages reading and interpolation of linoz data
! Created by: Francis Vitt
!-------------------------------------------------------------------
module linoz_data

  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils,   only : endrun
  use spmd_utils,   only : masterproc
  use tracer_data,  only : trfld,trfile
  use cam_logfile,  only : iulog

  implicit none

  private  ! all unless made public
  save 

  public :: fields
  public :: linoz_data_init
  public :: linoz_data_adv
  public :: init_linoz_data_restart
  public :: write_linoz_data_restart
  public :: read_linoz_data_restart
  public :: has_linoz_data
  public :: linoz_data_defaultopts
  public :: linoz_data_setopts

  type(trfld), pointer :: fields(:) => null()
  type(trfile) :: file

  logical :: has_linoz_data = .false.
  integer, parameter, public :: N_FLDS = 8
  integer :: number_flds

  character(len=256) :: filename = ''
  character(len=256) :: filelist = ''
  character(len=256) :: datapath = ''
  character(len=32)  :: datatype = 'CYCLICAL'
  logical            :: rmv_file = .false.
  integer            :: cycle_yr  = 0
  integer            :: fixed_ymd = 0
  integer            :: fixed_tod = 0

  character(len=16), dimension(N_FLDS), parameter :: fld_names = & ! data field names
       (/'o3_clim         ','t_clim          ','o3col_clim      ','PmL_clim        ', &
         'dPmL_dO3        ','dPmL_dT         ','dPmL_dO3col     ','cariolle_pscs   '/)

  character(len=16), dimension(N_FLDS), parameter :: fld_units = & ! data field names
       (/'vmr             ','K               ','Dobson Units    ','mr/s            ', &
         '/s              ','mr/K            ','mr/DU           ','/s              '/)

  integer :: index_map(N_FLDS)

  integer, public, parameter :: o3_clim_ndx = 1
  integer, public, parameter :: t_clim_ndx = 2
  integer, public, parameter :: o3col_clim_ndx = 3
  integer, public, parameter :: PmL_clim_ndx = 4

  integer, public, parameter :: dPmL_dO3_ndx = 5
  integer, public, parameter :: dPmL_dT_ndx = 6
  integer, public, parameter :: dPmL_dO3col_ndx = 7
  integer, public, parameter :: cariolle_pscs_ndx = 8

contains

!-------------------------------------------------------------------
!-------------------------------------------------------------------
! Subprogram not used   subroutine linoz_data_init()
! Subprogram not used 
! Subprogram not used     use tracer_data, only : trcdata_init
! Subprogram not used     use cam_history, only : addfld, phys_decomp
! Subprogram not used     use ppgrid,      only : pver
! Subprogram not used     use error_messages, only: handle_err
! Subprogram not used     use ppgrid,         only: pcols, pver, begchunk, endchunk
! Subprogram not used     use physics_buffer, only : physics_buffer_desc
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     integer :: ndx, istat, i
! Subprogram not used     
! Subprogram not used     if ( has_linoz_data ) then
! Subprogram not used        if ( masterproc ) then
! Subprogram not used           write(iulog,*) 'linoz_data_ini: linoz data :'//trim(filename)
! Subprogram not used        endif
! Subprogram not used     else
! Subprogram not used        return
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     allocate(file%in_pbuf(size(fld_names)))
! Subprogram not used     file%in_pbuf(:) = .false.
! Subprogram not used     call trcdata_init( fld_names, filename, filelist, datapath, fields, file, &
! Subprogram not used                        rmv_file, cycle_yr, fixed_ymd, fixed_tod, datatype)
! Subprogram not used         
! Subprogram not used     number_flds = 0
! Subprogram not used     if (associated(fields)) number_flds = size( fields )
! Subprogram not used 
! Subprogram not used     if( number_flds < 1 ) then
! Subprogram not used        if ( masterproc ) then
! Subprogram not used           write(iulog,*) 'linoz_data_init: There are no linoz data'
! Subprogram not used           write(iulog,*) ' '
! Subprogram not used        endif
! Subprogram not used        return
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     do i = 1,number_flds
! Subprogram not used        ndx = get_ndx( fields(i)%fldnam )
! Subprogram not used        index_map(i) = ndx
! Subprogram not used 
! Subprogram not used        if (ndx < 1) then
! Subprogram not used           call endrun('linoz_data_init: '//trim(fields(i)%fldnam)//' is not one of the named linoz data fields ')
! Subprogram not used        endif
! Subprogram not used        call addfld(fld_names(i), fld_units(i), pver, 'I', 'linoz data', phys_decomp )
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   end subroutine linoz_data_init

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine linoz_data_setopts(&
       linoz_data_file_in,      &
       linoz_data_filelist_in,  &
       linoz_data_path_in,      &
       linoz_data_type_in,      &
       linoz_data_rmfile_in,    &
       linoz_data_cycle_yr_in,  &
       linoz_data_fixed_ymd_in, &
       linoz_data_fixed_tod_in  &
       )

    implicit none

    character(len=*), intent(in), optional :: linoz_data_file_in
    character(len=*), intent(in), optional :: linoz_data_filelist_in
    character(len=*), intent(in), optional :: linoz_data_path_in
    character(len=*), intent(in), optional :: linoz_data_type_in
    logical,          intent(in), optional :: linoz_data_rmfile_in
    integer,          intent(in), optional :: linoz_data_cycle_yr_in
    integer,          intent(in), optional :: linoz_data_fixed_ymd_in
    integer,          intent(in), optional :: linoz_data_fixed_tod_in

    if ( present(linoz_data_file_in) ) then
       filename = linoz_data_file_in
    endif
    if ( present(linoz_data_filelist_in) ) then
       filelist = linoz_data_filelist_in
    endif
    if ( present(linoz_data_path_in) ) then
       datapath = linoz_data_path_in
    endif
    if ( present(linoz_data_type_in) ) then
       datatype = linoz_data_type_in
    endif
    if ( present(linoz_data_rmfile_in) ) then
       rmv_file = linoz_data_rmfile_in
    endif
    if ( present(linoz_data_cycle_yr_in) ) then
       cycle_yr = linoz_data_cycle_yr_in
    endif
    if ( present(linoz_data_fixed_ymd_in) ) then
       fixed_ymd = linoz_data_fixed_ymd_in
    endif
    if ( present(linoz_data_fixed_tod_in) ) then
       fixed_tod = linoz_data_fixed_tod_in
    endif

    if (len_trim(filename) > 0 ) has_linoz_data = .true.

  endsubroutine linoz_data_setopts

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine linoz_data_defaultopts(   &
       linoz_data_file_out,     &
       linoz_data_filelist_out, &
       linoz_data_path_out,     &
       linoz_data_type_out,     &
       linoz_data_rmfile_out,   &
       linoz_data_cycle_yr_out, &
       linoz_data_fixed_ymd_out,&
       linoz_data_fixed_tod_out &
       ) 

    implicit none

    character(len=*), intent(out), optional :: linoz_data_file_out
    character(len=*), intent(out), optional :: linoz_data_filelist_out
    character(len=*), intent(out), optional :: linoz_data_path_out
    character(len=*), intent(out), optional :: linoz_data_type_out
    logical,          intent(out), optional :: linoz_data_rmfile_out
    integer,          intent(out), optional :: linoz_data_cycle_yr_out
    integer,          intent(out), optional :: linoz_data_fixed_ymd_out
    integer,          intent(out), optional :: linoz_data_fixed_tod_out

    if ( present(linoz_data_file_out) ) then
       linoz_data_file_out = filename
    endif
    if ( present(linoz_data_filelist_out) ) then
       linoz_data_filelist_out = filelist
    endif
    if ( present(linoz_data_path_out) ) then
       linoz_data_path_out = datapath
    endif
    if ( present(linoz_data_type_out) ) then
       linoz_data_type_out = datatype
    endif
    if ( present(linoz_data_rmfile_out) ) then
       linoz_data_rmfile_out = rmv_file
    endif
    if ( present(linoz_data_cycle_yr_out) ) then
       linoz_data_cycle_yr_out = cycle_yr
    endif
    if ( present(linoz_data_fixed_ymd_out) ) then
       linoz_data_fixed_ymd_out = fixed_ymd
    endif
    if ( present(linoz_data_fixed_tod_out) ) then
       linoz_data_fixed_tod_out = fixed_tod
    endif

  endsubroutine linoz_data_defaultopts

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine linoz_data_adv( pbuf2d, state )

    use tracer_data,  only : advance_trcdata
    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk
    use ppgrid,       only : pcols, pver
    use string_utils, only : to_lower, GLC
    use cam_history,  only : outfld
    use physconst,    only : boltz                ! J/K/molecule
    use physics_buffer, only : physics_buffer_desc

    implicit none

  ! args
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    type(physics_state), intent(in):: state(begchunk:endchunk)                 

  ! local vars
    integer :: ind,c,ncol,i
    real(r8) :: to_mmr(pcols,pver)

    if( .not. has_linoz_data ) return

    call advance_trcdata( fields, file, state, pbuf2d  )
    
    ! set the tracer fields with the correct units
    do i = 1,number_flds
       ind = index_map(i)
       do c = begchunk,endchunk
          ncol = state(c)%ncol
          call outfld( fields(i)%fldnam, fields(i)%data(:ncol,:,c), ncol, state(c)%lchnk )
       enddo
    enddo

  end subroutine linoz_data_adv

!-------------------------------------------------------------------
!-------------------------------------------------------------------
! Subprogram not used   subroutine init_linoz_data_restart( piofile )
! Subprogram not used     use pio, only : file_desc_t
! Subprogram not used     use tracer_data, only : init_trc_restart
! Subprogram not used     implicit none
! Subprogram not used     type(file_desc_t),intent(inout) :: piofile     ! pio File pointer
! Subprogram not used 
! Subprogram not used     call init_trc_restart( 'linoz_data', piofile, file )
! Subprogram not used 
! Subprogram not used   end subroutine init_linoz_data_restart
!-------------------------------------------------------------------
! Subprogram not used   subroutine write_linoz_data_restart( PioFile )
! Subprogram not used     use tracer_data, only : write_trc_restart
! Subprogram not used     use pio, only : file_desc_t
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     type(file_desc_T) :: piofile
! Subprogram not used 
! Subprogram not used     call write_trc_restart( piofile, file )
! Subprogram not used 
! Subprogram not used   end subroutine write_linoz_data_restart

!-------------------------------------------------------------------
!-------------------------------------------------------------------
! Subprogram not used   subroutine read_linoz_data_restart( PioFile )
! Subprogram not used     use tracer_data, only : read_trc_restart
! Subprogram not used     use pio, only : file_desc_t
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     type(file_desc_T) :: piofile
! Subprogram not used 
! Subprogram not used     call read_trc_restart( 'linoz_data', piofile, file )
! Subprogram not used 
! Subprogram not used   end subroutine read_linoz_data_restart

!-------------------------------------------------------------------
!-------------------------------------------------------------------
! Subprogram not used   integer function get_ndx( name )
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used     character(len=*), intent(in) :: name
! Subprogram not used 
! Subprogram not used     integer :: i
! Subprogram not used 
! Subprogram not used     get_ndx = 0
! Subprogram not used     do i = 1,N_FLDS
! Subprogram not used       if ( trim(name) == trim(fld_names(i)) ) then
! Subprogram not used         get_ndx = i
! Subprogram not used         return
! Subprogram not used       endif
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used   end function get_ndx

end module linoz_data
