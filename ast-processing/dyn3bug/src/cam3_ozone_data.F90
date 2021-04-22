module cam3_ozone_data

!----------------------------------------------------------------------- 
! Purpose:
!
! Interpolates zonal ozone datasets used by CAM3 and puts the field 'O3' into
! the physics buffer.
! 
! Revision history:
! 2004-07-31  B. Eaton       Assemble module from comozp.F90, oznini.F90, oznint.F90, radozn.F90
! 2004-08-19  B. Eaton       Modify ozone_data_vert_interp to return mass mixing ratio.
! 2004-08-30  B. Eaton       Add ozone_data_get_cnst method.
! 2008 June   B. Eaton       Change name to cam3_ozone_data to support backwards compatibility
!                            for reading the CAM3 ozone data.  Add *_readnl method so module
!                            reads its own namelist.  Add cam3_ozone_data_on variable to
!                            turn the module on from the namelist.  By default it's off.
!-----------------------------------------------------------------------

use shr_kind_mod,   only: r8 => shr_kind_r8
use spmd_utils,     only: masterproc
use ppgrid,         only: begchunk, endchunk, pcols, pver
use abortutils,     only: endrun
use cam_logfile,    only: iulog
use physics_types,  only: physics_state
use boundarydata,   only: boundarydata_type, boundarydata_init, boundarydata_update, &
                          boundarydata_vert_interp
use mpishorthand

implicit none
private
save

! Public methods
public ::&
   cam3_ozone_data_readnl,        &! get namelist input
   cam3_ozone_data_register,      &! register ozone with physics buffer
   cam3_ozone_data_init,          &! open dataset and spatially interpolate data bounding initial time
   cam3_ozone_data_timestep_init   ! interpolate to current time

! Namelist variables
logical, public    :: cam3_ozone_data_on = .false. ! switch to turn module on/off
logical            :: ozncyc = .true.  ! .true. => assume annual cycle ozone data
character(len=256) :: bndtvo = ' '     ! full pathname for time-variant ozone dataset

! Local
integer            :: oz_idx           ! index into phys_buffer for ozone

type(boundarydata_type) :: ozonedata
character(len=6), parameter, dimension(1) :: nc_name = (/'OZONE '/) ! constituent names

!================================================================================================
contains
!================================================================================================

subroutine cam3_ozone_data_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'cam3_ozone_data_readnl'

   namelist /cam3_ozone_data_nl/ cam3_ozone_data_on, bndtvo, ozncyc
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'cam3_ozone_data_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, cam3_ozone_data_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

   ! Broadcast namelist variables
   call mpibcast(cam3_ozone_data_on, 1, mpilog, 0, mpicom)
   call mpibcast(bndtvo, len(bndtvo), mpichar, 0, mpicom)
   call mpibcast(ozncyc, 1, mpilog, 0, mpicom)

end subroutine cam3_ozone_data_readnl

!================================================================================================

! Subprogram not used subroutine cam3_ozone_data_register()
! Subprogram not used    use physics_buffer, only : pbuf_add_field, dtype_r8
! Subprogram not used 
! Subprogram not used    call pbuf_add_field('O3','physpkg',dtype_r8,(/pcols,pver/),oz_idx)
! Subprogram not used 
! Subprogram not used end subroutine cam3_ozone_data_register

!================================================================================================

! Subprogram not used subroutine cam3_ozone_data_init(phys_state)
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: Do initial read of time-variant ozone boundary dataset, containing
! Subprogram not used !          ozone mixing ratios as a function of latitude and pressure.  Read two
! Subprogram not used !          consecutive months between which the current date lies.  Routine
! Subprogram not used !          RADOZ2 then evaluates the two path length integrals (with and without
! Subprogram not used !          pressure weighting) from zero to the interfaces between the input
! Subprogram not used !          levels.  It also stores the contribution to the integral from each
! Subprogram not used !          layer.
! Subprogram not used ! 
! Subprogram not used ! Method: Call appropriate netcdf wrapper routines and interpolate to model grid
! Subprogram not used ! 
! Subprogram not used ! Author: CCM Core Group
! Subprogram not used ! Modified: P. Worley, August 2003, for chunking and performance optimization
! Subprogram not used !           J. Edwards, Dec 2005, functionality now performed by zonalbndrydata
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    use cam_history,      only: addfld, phys_decomp
! Subprogram not used 
! Subprogram not used    type(physics_state), intent(in) :: phys_state(begchunk:endchunk) 
! Subprogram not used    !-----------------------------------------------------------------------
! Subprogram not used     
! Subprogram not used    call addfld ('O3VMR', 'm3/m3', pver, 'A', 'Ozone volume mixing ratio', phys_decomp, sampling_seq='rad_lwsw')
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    ! Initialize for one field (arg_4=1) and do not vertically interpolate (arg_6=3)
! Subprogram not used    call boundarydata_init(bndtvo, phys_state, nc_name, 1, ozonedata, 3)
! Subprogram not used 
! Subprogram not used    if (masterproc) then
! Subprogram not used       write(iulog,*)'cam3_ozone_data_init: Initializing CAM3 prescribed ozone'
! Subprogram not used       write(iulog,*)'Time-variant boundary dataset (ozone) is: ', trim(bndtvo)
! Subprogram not used       if (ozncyc) then
! Subprogram not used          write(iulog,*)'OZONE dataset will be reused for each model year'
! Subprogram not used       else
! Subprogram not used          write(iulog,*)'OZONE dataset will not be cycled'
! Subprogram not used       end if
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used end subroutine cam3_ozone_data_init

!================================================================================================

! Subprogram not used subroutine cam3_ozone_data_timestep_init(pbuf2d,  phys_state)
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: Interpolate ozone mixing ratios to current time, reading in new monthly
! Subprogram not used !          data if necessary, and spatially interpolating it.
! Subprogram not used ! 
! Subprogram not used ! Method: Find next month of ozone data to interpolate.  Linearly interpolate 
! Subprogram not used !         vertically and horizontally
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    
! Subprogram not used    use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_get_chunk
! Subprogram not used 
! Subprogram not used    
! Subprogram not used    type(physics_state), intent(in) :: phys_state(begchunk:endchunk) 
! Subprogram not used    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
! Subprogram not used    real(r8),pointer :: tmpptr(:,:)
! Subprogram not used 
! Subprogram not used    integer lchnk
! Subprogram not used     
! Subprogram not used    call boundarydata_update(phys_state, ozonedata)
! Subprogram not used 
! Subprogram not used    do lchnk = begchunk, endchunk
! Subprogram not used       call pbuf_get_field(pbuf_get_chunk(pbuf2d, lchnk), oz_idx, tmpptr)
! Subprogram not used       call ozone_data_get_cnst(phys_state(lchnk), tmpptr)
! Subprogram not used    enddo
! Subprogram not used 
! Subprogram not used end subroutine cam3_ozone_data_timestep_init

!================================================================================================

! Subprogram not used subroutine ozone_data_get_cnst(state, q)
! Subprogram not used 
! Subprogram not used    use cam_history, only: outfld
! Subprogram not used    use physconst,   only: mwo3
! Subprogram not used 
! Subprogram not used    type(physics_state),  intent(in) :: state
! Subprogram not used    real(r8)                         :: q(:,:)     ! constituent mass mixing ratio
! Subprogram not used 
! Subprogram not used    ! local variables
! Subprogram not used    integer :: lchnk            ! chunk identifier
! Subprogram not used    integer :: i, k
! Subprogram not used    real(r8) :: ozmixin(pcols,ozonedata%levsiz)
! Subprogram not used    ! *** N.B. this hardwired mw of dry air needs to be changed to the share value
! Subprogram not used    real(r8), parameter :: mwdry = 28.9644_r8  ! Effective molecular weight of dry air (g/mol)
! Subprogram not used    real(r8), parameter :: mwr =  mwo3/mwdry   ! convert from the dataset values of vmr to mmr
! Subprogram not used    !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    lchnk = state%lchnk
! Subprogram not used 
! Subprogram not used    ozmixin=0._r8
! Subprogram not used    do k=1,ozonedata%levsiz
! Subprogram not used       do i=1,state%ncol
! Subprogram not used          ozmixin(i,k) = ozonedata%datainst(state%latmapback(i),k,lchnk,1)
! Subprogram not used       end do
! Subprogram not used    end do
! Subprogram not used    call boundarydata_vert_interp(lchnk, state%ncol, ozonedata%levsiz, &
! Subprogram not used                                  1, ozonedata%pin, state%pmid, ozmixin , q)
! Subprogram not used 
! Subprogram not used    call outfld('O3VMR', q, pcols, lchnk)
! Subprogram not used 
! Subprogram not used    do k=1,pver
! Subprogram not used       do i=1,state%ncol
! Subprogram not used          q(i,k) = mwr*q(i,k)
! Subprogram not used       end do
! Subprogram not used    end do
! Subprogram not used     
! Subprogram not used end subroutine ozone_data_get_cnst

!================================================================================================

end module cam3_ozone_data

