module cam_history_buffers
  use shr_kind_mod, only: r8 => shr_kind_r8
  use cam_history_support, only: max_chars, field_info, hentry
  use abortutils, only : endrun
  use pio, only : var_desc_t
  
  implicit none
!
! dim_index_2d, dim_index_3d: 2-D & 3-D dimension index lower & upper bounds
!
  type dim_index_2d                   ! 2-D dimension index
     sequence
     integer :: beg1, end1            ! lower & upper bounds of 1st dimension
     integer :: beg2, end2            ! lower & upper bounds of 2nd dimension
  end type dim_index_2d
  
  type dim_index_3d                   ! 3-D dimension index
     integer :: beg1, end1            ! lower & upper bounds of 1st dimension
     integer :: beg2, end2            ! lower & upper bounds of 2nd dimension
     integer :: beg3, end3            ! lower & upper bounds of 3rd dimension
  end type dim_index_3d


contains
  subroutine hbuf_accum_inst (buf8, field, nacs, dimind, idim, flag_xyfill, fillvalue)
    !
    !-----------------------------------------------------------------------
    !
    ! Purpose: Accumulate instantaneous values of field in 2-D hbuf.
    !          Set accumulation counter to 1.
    !
    !-----------------------------------------------------------------------
    !
    real(r8), pointer :: buf8(:,:)    ! 2-D history buffer
    integer, pointer                 :: nacs(:) ! accumulation counter
    integer, intent(in)              :: idim    ! Longitude dimension of field array
    logical, intent(in)              :: flag_xyfill ! non-applicable xy points flagged with fillvalue
    real(r8),          intent(in ) :: field(idim,*)   ! real*8 array
    type (dim_index_2d), intent(in ) :: dimind  ! 2-D dimension index
    real(r8), intent(in) :: fillvalue
    !
    ! Local indices
    !
    integer :: ib, ie    ! beginning and ending indices of first dimension
    integer :: jb, je    ! beginning and ending indices of second dimension
    integer :: ieu, jeu  ! number of elements in each dimension
    integer :: i, k      ! loop indices

    logical :: bad       ! flag indicates input field fillvalues not applied consistently
    ! with vertical level

    ib = dimind%beg1
    ie = dimind%end1
    jb = dimind%beg2
    je = dimind%end2

    ieu = ie-ib+1
    jeu = je-jb+1

    do k=1,jeu
       do i=1,ieu
          buf8(i,k) = field(i,k)
       end do
    end do

    if (flag_xyfill) then
       do i=1,ieu
          if (field(i,1) == fillvalue) then
             nacs(i) = 0
          else
             nacs(i) = 1
          end if
       end do
    else
       nacs(1) = 1
    end if

    return
  end subroutine hbuf_accum_inst
  !#######################################################################

! Subprogram not used   subroutine hbuf_accum_add (buf8, field, nacs, dimind, idim, flag_xyfill, fillvalue)
! Subprogram not used     !
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     ! Purpose: Add the values of field to 2-D hbuf.
! Subprogram not used     !          Increment accumulation counter by 1.
! Subprogram not used     !
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     real(r8), pointer :: buf8(:,:)    ! 2-D history buffer
! Subprogram not used     integer, pointer                 :: nacs(:) ! accumulation counter
! Subprogram not used     integer, intent(in) :: idim           ! Longitude dimension of field array
! Subprogram not used     logical, intent(in)              :: flag_xyfill ! non-applicable xy points flagged with fillvalue
! Subprogram not used     real(r8),          intent(in ) :: field(idim,*)   ! real*8 array
! Subprogram not used     type (dim_index_2d), intent(in ) :: dimind  ! 2-D dimension index
! Subprogram not used     real(r8), intent(in) :: fillvalue
! Subprogram not used     !
! Subprogram not used     ! Local indices
! Subprogram not used     !
! Subprogram not used     integer :: ib, ie    ! beginning and ending indices of first dimension
! Subprogram not used     integer :: jb, je    ! beginning and ending indices of second dimension
! Subprogram not used     integer :: ieu, jeu  ! number of elements in each dimension
! Subprogram not used     integer :: i,k       ! indices
! Subprogram not used 
! Subprogram not used     ib = dimind%beg1
! Subprogram not used     ie = dimind%end1
! Subprogram not used     jb = dimind%beg2
! Subprogram not used     je = dimind%end2
! Subprogram not used 
! Subprogram not used     ieu = ie-ib+1
! Subprogram not used     jeu = je-jb+1
! Subprogram not used 
! Subprogram not used     if (flag_xyfill) then
! Subprogram not used        do k=1,jeu
! Subprogram not used           do i=1,ieu
! Subprogram not used              if (field(i,k) /= fillvalue) then
! Subprogram not used                 buf8(i,k) = buf8(i,k) + field(i,k)
! Subprogram not used              end if
! Subprogram not used           end do
! Subprogram not used        end do
! Subprogram not used        !
! Subprogram not used        ! Ensure input field has fillvalue defined invariant in the z-direction, then increment nacs
! Subprogram not used        !
! Subprogram not used        call check_accum (field, idim, ieu, jeu, fillvalue)
! Subprogram not used        do i=1,ieu
! Subprogram not used           if (field(i,1) /= fillvalue) then
! Subprogram not used              nacs(i) = nacs(i) + 1
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used     else
! Subprogram not used        do k=1,jeu
! Subprogram not used           do i=1,ieu
! Subprogram not used              buf8(i,k) = buf8(i,k) + field(i,k)
! Subprogram not used           end do
! Subprogram not used        end do
! Subprogram not used        nacs(1) = nacs(1) + 1
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     return
! Subprogram not used   end subroutine hbuf_accum_add

  !#######################################################################

! Subprogram not used   subroutine hbuf_accum_add00z (buf8, field, nacs, dimind, idim, flag_xyfill, fillvalue)
! Subprogram not used     !
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     ! Purpose: Add the values of field to 2-D hbuf, only of time is 00z.
! Subprogram not used     !          Increment accumulation counter by 1.
! Subprogram not used     !
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     use time_manager, only:  get_curr_date
! Subprogram not used 
! Subprogram not used     real(r8), pointer :: buf8(:,:)    ! 2-D history buffer
! Subprogram not used     integer, pointer                 :: nacs(:) ! accumulation counter
! Subprogram not used     integer, intent(in) :: idim           ! Longitude dimension of field array
! Subprogram not used     logical, intent(in)              :: flag_xyfill ! non-applicable xy points flagged with fillvalue
! Subprogram not used     real(r8),          intent(in ) :: field(idim,*)   ! real*8 array
! Subprogram not used     type (dim_index_2d), intent(in ) :: dimind  ! 2-D dimension index
! Subprogram not used     real(r8), intent(in) :: fillvalue
! Subprogram not used     !
! Subprogram not used     ! Local indices
! Subprogram not used     !
! Subprogram not used     integer :: ib, ie    ! beginning and ending indices of first dimension
! Subprogram not used     integer :: jb, je    ! beginning and ending indices of second dimension
! Subprogram not used     integer :: ieu, jeu  ! number of elements in each dimension
! Subprogram not used     integer :: i,k       ! indices
! Subprogram not used     integer :: yr, mon, day, tod
! Subprogram not used 
! Subprogram not used     ! get the time of day, return if not 00z
! Subprogram not used     call get_curr_date (yr,mon,day,tod)
! Subprogram not used     if (tod /= 0) return
! Subprogram not used 
! Subprogram not used     ib = dimind%beg1
! Subprogram not used     ie = dimind%end1
! Subprogram not used     jb = dimind%beg2
! Subprogram not used     je = dimind%end2
! Subprogram not used 
! Subprogram not used     ieu = ie-ib+1
! Subprogram not used     jeu = je-jb+1
! Subprogram not used 
! Subprogram not used     if (flag_xyfill) then
! Subprogram not used        do k=1,jeu
! Subprogram not used           do i=1,ieu
! Subprogram not used              if (field(i,k) /= fillvalue) then
! Subprogram not used                 buf8(i,k) = buf8(i,k) + field(i,k)
! Subprogram not used              end if
! Subprogram not used           end do
! Subprogram not used        end do
! Subprogram not used        !
! Subprogram not used        ! Ensure input field has fillvalue defined invariant in the z-direction, then increment nacs
! Subprogram not used        !
! Subprogram not used        call check_accum (field, idim, ieu, jeu, fillvalue)
! Subprogram not used        do i=1,ieu
! Subprogram not used           if (field(i,1) /= fillvalue) then
! Subprogram not used              nacs(i) = nacs(i) + 1
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used     else
! Subprogram not used        do k=1,jeu
! Subprogram not used           do i=1,ieu
! Subprogram not used              buf8(i,k) = buf8(i,k) + field(i,k)
! Subprogram not used           end do
! Subprogram not used        end do
! Subprogram not used        nacs(1) = nacs(1) + 1
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     return
! Subprogram not used   end subroutine hbuf_accum_add00z

  !#######################################################################

! Subprogram not used   subroutine hbuf_accum_max (buf8, field, nacs, dimind, idim, flag_xyfill, fillvalue)
! Subprogram not used     !
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     ! Purpose: Accumulate the maximum values of field in 2-D hbuf
! Subprogram not used     !          Set accumulation counter to 1.
! Subprogram not used     !
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     real(r8), pointer :: buf8(:,:)    ! 2-D history buffer
! Subprogram not used     integer, pointer                 :: nacs(:) ! accumulation counter
! Subprogram not used     integer, intent(in) :: idim           ! Longitude dimension of field array
! Subprogram not used     logical, intent(in)              :: flag_xyfill ! non-applicable xy points flagged with fillvalue
! Subprogram not used     real(r8),          intent(in ) :: field(idim,*)   ! real*8 array
! Subprogram not used     type (dim_index_2d), intent(in ) :: dimind  ! 2-D dimension index
! Subprogram not used     real(r8), intent(in) :: fillvalue
! Subprogram not used     !
! Subprogram not used     ! Local indices
! Subprogram not used     !
! Subprogram not used     integer :: ib, ie    ! beginning and ending indices of first dimension
! Subprogram not used     integer :: jb, je    ! beginning and ending indices of second dimension
! Subprogram not used     integer :: ieu, jeu  ! number of elements in each dimension
! Subprogram not used     integer :: i, k
! Subprogram not used 
! Subprogram not used     ib = dimind%beg1
! Subprogram not used     ie = dimind%end1
! Subprogram not used     jb = dimind%beg2
! Subprogram not used     je = dimind%end2
! Subprogram not used 
! Subprogram not used     ieu = ie-ib+1
! Subprogram not used     jeu = je-jb+1
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     if (flag_xyfill) then
! Subprogram not used        do k=1,jeu
! Subprogram not used           do i=1,ieu
! Subprogram not used              if (nacs(i) == 0) then
! Subprogram not used                 buf8(i,k) = -huge (buf8)
! Subprogram not used              end if
! Subprogram not used              if (field(i,k) > buf8(i,k) .and. field(i,k) /= fillvalue) then
! Subprogram not used                 buf8(i,k) = field(i,k)
! Subprogram not used              end if
! Subprogram not used           end do
! Subprogram not used        end do
! Subprogram not used     else
! Subprogram not used        do k=1,jeu
! Subprogram not used           do i=1,ieu
! Subprogram not used              if (nacs(1) == 0) then
! Subprogram not used                 buf8(i,k) = -huge (buf8)
! Subprogram not used              end if
! Subprogram not used              if (field(i,k) > buf8(i,k)) then
! Subprogram not used                 buf8(i,k) = field(i,k)
! Subprogram not used              end if
! Subprogram not used           end do
! Subprogram not used        end do
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     if (flag_xyfill) then
! Subprogram not used        call check_accum (field, idim, ieu, jeu,fillvalue)
! Subprogram not used        do i=1,ieu
! Subprogram not used           if (field(i,1) /= fillvalue) then
! Subprogram not used              nacs(i) = 1
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used     else
! Subprogram not used        nacs(1) = 1
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     return
! Subprogram not used   end subroutine hbuf_accum_max

  !#######################################################################

! Subprogram not used   subroutine hbuf_accum_min (buf8, field, nacs, dimind, idim, flag_xyfill, fillvalue)
! Subprogram not used     !
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     ! Purpose: Accumulate the minimum values of field in 2-D hbuf
! Subprogram not used     !          Set accumulation counter to 1.
! Subprogram not used     !
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     real(r8), pointer :: buf8(:,:)    ! 2-D history buffer
! Subprogram not used     integer, pointer                 :: nacs(:) ! accumulation counter
! Subprogram not used     integer, intent(in) :: idim           ! Longitude dimension of field array
! Subprogram not used     logical, intent(in)              :: flag_xyfill ! non-applicable xy points flagged with fillvalue
! Subprogram not used     real(r8),          intent(in ) :: field(idim,*)   ! real*8 array
! Subprogram not used     type (dim_index_2d), intent(in ) :: dimind  ! 2-D dimension index
! Subprogram not used     real(r8), intent(in) :: fillvalue
! Subprogram not used     !
! Subprogram not used     ! Local indices
! Subprogram not used     !
! Subprogram not used     integer :: ib, ie    ! beginning and ending indices of first dimension
! Subprogram not used     integer :: jb, je    ! beginning and ending indices of second dimension
! Subprogram not used     integer :: ieu, jeu  ! number of elements in each dimension
! Subprogram not used     integer :: i, k
! Subprogram not used 
! Subprogram not used     ib = dimind%beg1
! Subprogram not used     ie = dimind%end1
! Subprogram not used     jb = dimind%beg2
! Subprogram not used     je = dimind%end2
! Subprogram not used 
! Subprogram not used     ieu = ie-ib+1
! Subprogram not used     jeu = je-jb+1
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     if (flag_xyfill) then
! Subprogram not used        do k=1,jeu
! Subprogram not used           do i=1,ieu
! Subprogram not used              if (nacs(i) == 0) then
! Subprogram not used                 buf8(i,k) = +huge (buf8)
! Subprogram not used              end if
! Subprogram not used              if (field(i,k) < buf8(i,k) .and. field(i,k) /= fillvalue) then
! Subprogram not used                 buf8(i,k) = field(i,k)
! Subprogram not used              end if
! Subprogram not used           end do
! Subprogram not used        end do
! Subprogram not used     else
! Subprogram not used        do k=1,jeu
! Subprogram not used           do i=1,ieu
! Subprogram not used              if (nacs(1) == 0) then
! Subprogram not used                 buf8(i,k) = +huge (buf8)
! Subprogram not used              end if
! Subprogram not used              if (field(i,k) < buf8(i,k)) then
! Subprogram not used                 buf8(i,k) = field(i,k)
! Subprogram not used              end if
! Subprogram not used           end do
! Subprogram not used        end do
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     if (flag_xyfill) then
! Subprogram not used        call check_accum (field, idim, ieu, jeu, fillvalue)
! Subprogram not used        do i=1,ieu
! Subprogram not used           if (field(i,1) /= fillvalue) then
! Subprogram not used              nacs(i) = 1
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used     else
! Subprogram not used        nacs(1) = 1
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     return
! Subprogram not used   end subroutine hbuf_accum_min

! Subprogram not used   subroutine hbuf_accum_addlcltime (buf8, field, nacs, dimind, idim, flag_xyfill, fillvalue, c , decomp_type,&
! Subprogram not used        lcltod_start, lcltod_stop)
! Subprogram not used     !
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     ! Purpose: Add the values of field to 2-D hbuf, only if the local time
! Subprogram not used     !          is in the range specified.
! Subprogram not used     !          Increment accumulation counter by 1.
! Subprogram not used     !
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     use time_manager, only:  get_curr_date
! Subprogram not used     use phys_grid,    only:  get_rlon_all_p
! Subprogram not used     use physconst,    only:  pi
! Subprogram not used     use phys_grid,     only: get_ncols_p, pcols
! Subprogram not used     use cam_pio_utils, only: phys_decomp
! Subprogram not used     use dyn_grid,      only: get_horiz_grid_dim_d, get_horiz_grid_d, dyn_grid_get_elem_coords
! Subprogram not used 
! Subprogram not used     type (dim_index_2d), intent(in ) :: dimind  ! 2-D dimension index
! Subprogram not used     real(r8), pointer    :: buf8(:,:)    ! 2-D history buffer
! Subprogram not used     integer,  pointer    :: nacs(:) ! accumulation counter
! Subprogram not used     integer,  intent(in) :: idim           ! Longitude dimension of field array
! Subprogram not used     logical,  intent(in) :: flag_xyfill ! non-applicable xy points flagged with fillvalue
! Subprogram not used     real(r8), intent(in) :: field(idim,*)   ! real*8 array
! Subprogram not used     integer,  intent(in) :: c              ! chunk (physics) or latitude (dynamics) index
! Subprogram not used 
! Subprogram not used     integer, intent(in)  :: decomp_type, lcltod_start, lcltod_stop
! Subprogram not used     real(r8), intent(in) :: fillvalue
! Subprogram not used     
! Subprogram not used     !
! Subprogram not used     ! Local indices
! Subprogram not used     !
! Subprogram not used     integer  :: ib, ie    ! beginning and ending indices of first dimension
! Subprogram not used     integer  :: jb, je    ! beginning and ending indices of second dimension
! Subprogram not used     integer  :: ieu, jeu  ! number of elements in each dimension
! Subprogram not used     integer  :: i,k       ! indices
! Subprogram not used     integer  :: yr, mon, day, tod
! Subprogram not used     integer  :: ncols, hdim1_d, hdim2_d
! Subprogram not used 
! Subprogram not used     integer,  allocatable :: lcltod(:)     ! local time of day (secs)
! Subprogram not used     logical,  allocatable :: inavg(:)      ! is the column in the desired local time range?
! Subprogram not used     real(r8), allocatable :: rlon(:)       ! column longitude (radians)
! Subprogram not used     integer,  allocatable :: cdex(:)       ! global column index 
! Subprogram not used 
! Subprogram not used     ib = dimind%beg1
! Subprogram not used     ie = dimind%end1
! Subprogram not used     jb = dimind%beg2
! Subprogram not used     je = dimind%end2
! Subprogram not used 
! Subprogram not used     ieu = ie-ib+1
! Subprogram not used     jeu = je-jb+1
! Subprogram not used 
! Subprogram not used     allocate( inavg(1:ieu) , lcltod(1:ieu) )
! Subprogram not used     lcltod(:) = 0
! Subprogram not used 
! Subprogram not used     !
! Subprogram not used     ! Get the time of day and longitude and compute the local time.
! Subprogram not used     !
! Subprogram not used     call get_curr_date (yr,mon,day,tod)      
! Subprogram not used 
! Subprogram not used     if ( decomp_type == phys_decomp ) then 
! Subprogram not used 
! Subprogram not used        ncols = get_ncols_p(c)
! Subprogram not used        ie = ncols
! Subprogram not used        ieu = ncols
! Subprogram not used        allocate( rlon(ncols) )
! Subprogram not used        call get_rlon_all_p(c, ncols, rlon)
! Subprogram not used        lcltod(1:ncols) = mod((tod) + (nint(86400_r8 * rlon(1:ncols) / 2._r8 / pi)), 86400)
! Subprogram not used 
! Subprogram not used     else 
! Subprogram not used 
! Subprogram not used        ncols = ieu
! Subprogram not used        allocate(rlon(ncols),cdex(ncols))
! Subprogram not used        call dyn_grid_get_elem_coords( c, rlon=rlon, cdex=cdex )
! Subprogram not used        lcltod(:) = -999999
! Subprogram not used        where( cdex(:)>0 ) lcltod(1:ieu) = mod((tod) + (nint(86400_r8 * rlon(1:ncols) / 2._r8 / pi)), 86400)
! Subprogram not used 
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     !
! Subprogram not used     ! Set a flag to indicate that the column is in the requested local time range.      
! Subprogram not used     ! If lcltod_stop is less than lcltod_stop, then the time is wrapping around 24 hours.
! Subprogram not used     !
! Subprogram not used     inavg(:)  = .false.
! Subprogram not used 
! Subprogram not used     if (lcltod_stop < lcltod_start) then
! Subprogram not used        ! the ".and.(lcltod(:)>0" condition was added to exclude the undifined (-999999) columns
! Subprogram not used        where((lcltod(:) >= lcltod_start) .or. ((lcltod(:) < lcltod_stop).and.(lcltod(:)>0))) inavg(:) = .true.
! Subprogram not used     else if (lcltod_stop == lcltod_start) then
! Subprogram not used        where(lcltod(ib:ie) == lcltod_start) inavg(1:ieu) = .true.
! Subprogram not used     else
! Subprogram not used        where((lcltod(:) >= lcltod_start) .and. (lcltod(:) < lcltod_stop)) inavg(:) = .true.
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     if (flag_xyfill) then
! Subprogram not used        do k=1,jeu
! Subprogram not used           do i=1,ieu
! Subprogram not used              if (inavg(i) .and. (field(i,k) /= fillvalue)) then
! Subprogram not used                 buf8(i,k) = buf8(i,k) + field(i,k)
! Subprogram not used              end if
! Subprogram not used           end do
! Subprogram not used        end do
! Subprogram not used        !
! Subprogram not used        ! Ensure input field has fillvalue defined invariant in the z-direction, then increment nacs
! Subprogram not used        !
! Subprogram not used        call check_accum (field, idim, ieu, jeu, fillvalue)
! Subprogram not used        do i=1,ieu
! Subprogram not used           if (inavg(i) .and. (field(i,1) /= fillvalue)) then
! Subprogram not used              nacs(i) = nacs(i) + 1
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used     else
! Subprogram not used        do k=1,jeu
! Subprogram not used           do i=1,ieu
! Subprogram not used              if (inavg(i)) then
! Subprogram not used                 buf8(i,k) = buf8(i,k) + field(i,k)
! Subprogram not used              end if
! Subprogram not used           end do
! Subprogram not used        end do
! Subprogram not used        !
! Subprogram not used        ! NOTE: Because of the local time, some columns in the chunk may not be used in the
! Subprogram not used        ! average, so nacs need to be dimension with the number of columns unlike the other 
! Subprogram not used        ! averaging techniques.
! Subprogram not used        !
! Subprogram not used        do i=1,ieu
! Subprogram not used           if (inavg(i)) nacs(i) = nacs(i) + 1
! Subprogram not used        end do
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     deallocate( inavg , rlon, lcltod )
! Subprogram not used     if (allocated(cdex)) deallocate(cdex)
! Subprogram not used 
! Subprogram not used     return
! Subprogram not used   end subroutine hbuf_accum_addlcltime

!#######################################################################


  !#######################################################################

! Subprogram not used   subroutine check_accum (field, idim, ieu, jeu, fillvalue)
! Subprogram not used     !
! Subprogram not used     integer, intent(in)  :: idim
! Subprogram not used     real(r8), intent(in) :: field(idim,*)   ! real*8 array
! Subprogram not used     integer, intent(in)  :: ieu,jeu         ! loop ranges
! Subprogram not used 
! Subprogram not used     logical :: bad
! Subprogram not used     integer :: i,k
! Subprogram not used     real(r8), intent(in) :: fillvalue
! Subprogram not used     !
! Subprogram not used     ! For multilevel fields ensure that all levels have fillvalue applied consistently
! Subprogram not used     !
! Subprogram not used     bad = .false.
! Subprogram not used     do k=2,jeu
! Subprogram not used        do i=1,ieu
! Subprogram not used           if (field(i,1) == fillvalue .and. field(i,k) /= fillvalue .or. &
! Subprogram not used                field(i,1) /= fillvalue .and. field(i,k) == fillvalue) then
! Subprogram not used              bad = .true.
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     if (bad) then
! Subprogram not used        call endrun ('CHECK_ACCUM: inconsistent level values')
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     return
! Subprogram not used   end subroutine check_accum

end module cam_history_buffers
