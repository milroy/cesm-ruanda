module interpolate_data
! Description:
!   Routines for interpolation of data in latitude, longitude and time.
!
! Author: Gathered from a number of places and put into the current format by Jim Edwards
!
! Modules Used:
!
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use abortutils,     only: endrun
  use cam_logfile,    only: iulog
  implicit none
  private
!
! Public Methods:
!

  public :: interp_type, lininterp, vertinterp, bilin, get_timeinterp_factors
  public :: lininterp_init, lininterp_finish
  type interp_type
     real(r8), pointer :: wgts(:)
     real(r8), pointer :: wgtn(:)
     integer, pointer  :: jjm(:)
     integer, pointer  :: jjp(:)
  end type interp_type
  interface lininterp
     module procedure lininterp_original
     module procedure lininterp_full1d
     module procedure lininterp1d
     module procedure lininterp2d2d
     module procedure lininterp2d1d
     module procedure lininterp3d2d
  end interface

integer, parameter, public :: extrap_method_zero = 0
integer, parameter, public :: extrap_method_bndry = 1
integer, parameter, public :: extrap_method_cycle = 2

contains
! Subprogram not used   subroutine lininterp_full1d (arrin, yin, nin, arrout, yout, nout)
! Subprogram not used     integer, intent(in) :: nin, nout
! Subprogram not used     real(r8), intent(in) :: arrin(nin), yin(nin), yout(nout)
! Subprogram not used     real(r8), intent(out) :: arrout(nout)
! Subprogram not used     type (interp_type) :: interp_wgts
! Subprogram not used 
! Subprogram not used     call lininterp_init(yin, nin, yout, nout, extrap_method_bndry, interp_wgts)
! Subprogram not used     call lininterp1d(arrin, nin, arrout, nout, interp_wgts)
! Subprogram not used     call lininterp_finish(interp_wgts)
! Subprogram not used 
! Subprogram not used   end subroutine lininterp_full1d

  subroutine lininterp_init(yin, nin, yout, nout, extrap_method, interp_wgts, &
       cyclicmin, cyclicmax)
!
! Description:
!   Initialize a variable of type(interp_type) with weights for linear interpolation.
!       this variable can then be used in calls to lininterp1d and lininterp2d.
!   yin is a 1d array of length nin of locations to interpolate from - this array must
!       be monotonic but can be increasing or decreasing
!   yout is a 1d array of length nout of locations to interpolate to, this array need
!       not be ordered
!   extrap_method determines how to handle yout points beyond the bounds of yin
!       if 0 set values outside output grid to 0
!       if 1 set to boundary value
!       if 2 set to cyclic boundaries
!         optional values cyclicmin and cyclicmax can be used to set the bounds of the
!         cyclic mapping - these default to 0 and 360.
!

    integer, intent(in) :: nin
    integer, intent(in) :: nout
    real(r8), intent(in) :: yin(:)           ! input mesh
    real(r8), intent(in) :: yout(:)         ! output mesh
    integer, intent(in) :: extrap_method       ! if 0 set values outside output grid to 0
                                               ! if 1 set to boundary value
                                               ! if 2 set to cyclic boundaries
    real(r8), intent(in), optional :: cyclicmin, cyclicmax

    type (interp_type), intent(out) :: interp_wgts
    real(r8) :: cmin, cmax
    real(r8) :: extrap
    real(r8) :: dyinwrap
    real(r8) :: ratio
    real(r8) :: avgdyin
    integer :: i, j, icount
    integer :: jj
    real(r8), pointer :: wgts(:)
    real(r8), pointer :: wgtn(:)
    integer, pointer :: jjm(:)
    integer, pointer :: jjp(:)
    logical :: increasing
    !
    ! Check validity of input coordinate arrays: must be monotonically increasing,
    ! and have a total of at least 2 elements
    !
    if (nin.lt.2) then
       call endrun('LININTERP: Must have at least 2 input points for interpolation')
    end if
    if(present(cyclicmin)) then
       cmin=cyclicmin
    else
       cmin=0_r8
    end if
    if(present(cyclicmax)) then
       cmax=cyclicmax
    else
       cmax=360_r8
    end if
    if(cmax<=cmin) then
       call endrun('LININTERP: cyclic min value must be < max value')
    end if
    increasing=.true.
    icount = 0
    do j=1,nin-1
       if (yin(j).gt.yin(j+1)) icount = icount + 1
    end do
    if(icount.eq.nin-1) then
       increasing = .false.
       icount=0
    endif
    if (icount.gt.0) then
       call endrun('LININTERP: Non-monotonic input coordinate array found')
    end if
    allocate(interp_wgts%jjm(nout), &
         interp_wgts%jjp(nout), &
         interp_wgts%wgts(nout), &
         interp_wgts%wgtn(nout))

    jjm => interp_wgts%jjm
    jjp => interp_wgts%jjp
    wgts =>  interp_wgts%wgts
    wgtn =>  interp_wgts%wgtn

    !
    ! Initialize index arrays for later checking
    !
    jjm = 0
    jjp = 0

    extrap = 0._r8
    select case (extrap_method)
    case (extrap_method_zero)
       !
       ! For values which extend beyond boundaries, set weights
       ! such that values will be 0.
       !
       do j=1,nout
          if(increasing) then
             if (yout(j).lt.yin(1)) then
                jjm(j) = 1
                jjp(j) = 1
                wgts(j) = 0._r8
                wgtn(j) = 0._r8
                extrap = extrap + 1._r8
             else if (yout(j).gt.yin(nin)) then
                jjm(j) = nin
                jjp(j) = nin
                wgts(j) = 0._r8
                wgtn(j) = 0._r8
                extrap = extrap + 1._r8
             end if
          else
             if (yout(j).gt.yin(1)) then
                jjm(j) = 1
                jjp(j) = 1
                wgts(j) = 0._r8
                wgtn(j) = 0._r8
                extrap = extrap + 1._r8
             else if (yout(j).lt.yin(nin)) then
                jjm(j) = nin
                jjp(j) = nin
                wgts(j) = 0._r8
                wgtn(j) = 0._r8
                extrap = extrap + 1._r8
             end if
          end if
       end do
    case (extrap_method_bndry)
       !
       ! For values which extend beyond boundaries, set weights
       ! such that values will just be copied.
       !
       do j=1,nout
          if(increasing) then
             if (yout(j).le.yin(1)) then
                jjm(j) = 1
                jjp(j) = 1
                wgts(j) = 1._r8
                wgtn(j) = 0._r8
                extrap = extrap + 1._r8
             else if (yout(j).gt.yin(nin)) then
                jjm(j) = nin
                jjp(j) = nin
                wgts(j) = 1._r8
                wgtn(j) = 0._r8
                extrap = extrap + 1._r8
             end if
          else
             if (yout(j).gt.yin(1)) then
                jjm(j) = 1
                jjp(j) = 1
                wgts(j) = 1._r8
                wgtn(j) = 0._r8
                extrap = extrap + 1._r8
             else if (yout(j).le.yin(nin)) then
                jjm(j) = nin
                jjp(j) = nin
                wgts(j) = 1._r8
                wgtn(j) = 0._r8
                extrap = extrap + 1._r8
             end if
          end if
       end do
    case (extrap_method_cycle)
       !
       ! For values which extend beyond boundaries, set weights
       ! for circular boundaries
       !
       dyinwrap = yin(1) + (cmax-cmin) - yin(nin)
       avgdyin = abs(yin(nin)-yin(1))/(nin-1._r8)
       ratio = dyinwrap/avgdyin
       if (ratio < 0.9_r8 .or. ratio > 1.1_r8) then
          write(iulog,*) 'Lininterp: Bad dyinwrap value =',dyinwrap,&
               ' avg=', avgdyin, yin(1),yin(nin)
          call endrun('interpolate_data')
       end if

       do j=1,nout
          if(increasing) then
             if (yout(j) <= yin(1)) then
                jjm(j) = nin
                jjp(j) = 1
                wgts(j) = (yin(1)-yout(j))/dyinwrap
                wgtn(j) = (yout(j)+(cmax-cmin) - yin(nin))/dyinwrap
             else if (yout(j) > yin(nin)) then
                jjm(j) = nin
                jjp(j) = 1
                wgts(j) = (yin(1)+(cmax-cmin)-yout(j))/dyinwrap
                wgtn(j) = (yout(j)-yin(nin))/dyinwrap
             end if
          else
             if (yout(j) > yin(1)) then
                jjm(j) = nin
                jjp(j) = 1
                wgts(j) = (yin(1)-yout(j))/dyinwrap
                wgtn(j) = (yout(j)+(cmax-cmin) - yin(nin))/dyinwrap
             else if (yout(j) <= yin(nin)) then
                jjm(j) = nin
                jjp(j) = 1
                wgts(j) = (yin(1)+(cmax-cmin)-yout(j))/dyinwrap
                wgtn(j) = (yout(j)+(cmax-cmin)-yin(nin))/dyinwrap
             end if

          endif
       end do
    end select

    !
    ! Loop though output indices finding input indices and weights
    !
    if(increasing) then
       do j=1,nout
          do jj=1,nin-1
             if (yout(j).gt.yin(jj) .and. yout(j).le.yin(jj+1)) then
                jjm(j) = jj
                jjp(j) = jj + 1
                wgts(j) = (yin(jj+1)-yout(j))/(yin(jj+1)-yin(jj))
                wgtn(j) = (yout(j)-yin(jj))/(yin(jj+1)-yin(jj))
                exit
             end if
          end do
       end do
    else
       do j=1,nout
          do jj=1,nin-1
             if (yout(j).le.yin(jj) .and. yout(j).gt.yin(jj+1)) then
                jjm(j) = jj
                jjp(j) = jj + 1
                wgts(j) = (yin(jj+1)-yout(j))/(yin(jj+1)-yin(jj))
                wgtn(j) = (yout(j)-yin(jj))/(yin(jj+1)-yin(jj))
                exit
             end if
          end do
       end do
    end if

    !
    ! Check that interp/extrap points have been found for all outputs
    !
    icount = 0
    do j=1,nout
       if (jjm(j).eq.0 .or. jjp(j).eq.0) icount = icount + 1
       ratio=wgts(j)+wgtn(j)
       if((ratio<0.9_r8.or.ratio>1.1_r8).and.extrap_method.ne.0) then
          write(iulog,*) j, wgts(j),wgtn(j),jjm(j),jjp(j), increasing,extrap_method
          call endrun('Bad weight computed in LININTERP_init')
       end if
    end do
    if (icount.gt.0) then
       call endrun('LININTERP: Point found without interp indices')
    end if

  end subroutine lininterp_init

  subroutine lininterp1d (arrin, n1, arrout, m1, interp_wgts)
    !-----------------------------------------------------------------------
    !
    ! Purpose: Do a linear interpolation from input mesh to output
    !          mesh with weights as set in lininterp_init.
    !
    !
    ! Author: Jim Edwards
    !
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    implicit none
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    integer, intent(in) :: n1                 ! number of input latitudes
    integer, intent(in) :: m1                ! number of output latitudes

    real(r8), intent(in) :: arrin(n1)    ! input array of values to interpolate
    type(interp_type), intent(in) :: interp_wgts
    real(r8), intent(out) :: arrout(m1) ! interpolated array

    !
    ! Local workspace
    !
    integer j                ! latitude indices
    integer, pointer :: jjm(:)
    integer, pointer :: jjp(:)

    real(r8), pointer :: wgts(:)
    real(r8), pointer :: wgtn(:)


    jjm => interp_wgts%jjm
    jjp => interp_wgts%jjp
    wgts =>  interp_wgts%wgts
    wgtn =>  interp_wgts%wgtn

    !
    ! Do the interpolation
    !
    do j=1,m1
      arrout(j) = arrin(jjm(j))*wgts(j) + arrin(jjp(j))*wgtn(j)
    end do

    return
  end subroutine lininterp1d

! Subprogram not used   subroutine lininterp2d2d(arrin, n1, n2, arrout, m1, m2, wgt1, wgt2)
! Subprogram not used     implicit none
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     ! Arguments
! Subprogram not used     !
! Subprogram not used     integer, intent(in) :: n1, n2, m1, m2
! Subprogram not used     real(r8), intent(in) :: arrin(n1,n2)    ! input array of values to interpolate
! Subprogram not used     type(interp_type), intent(in) :: wgt1, wgt2
! Subprogram not used     real(r8), intent(out) :: arrout(m1,m2) ! interpolated array
! Subprogram not used     !
! Subprogram not used     ! locals
! Subprogram not used     !
! Subprogram not used     integer i,j                ! indices
! Subprogram not used     integer, pointer :: iim(:), jjm(:)
! Subprogram not used     integer, pointer :: iip(:), jjp(:)
! Subprogram not used 
! Subprogram not used     real(r8), pointer :: wgts1(:), wgts2(:)
! Subprogram not used     real(r8), pointer :: wgtn1(:), wgtn2(:)
! Subprogram not used 
! Subprogram not used     real(r8) :: arrtmp(n1,m2)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     jjm => wgt2%jjm
! Subprogram not used     jjp => wgt2%jjp
! Subprogram not used     wgts2 => wgt2%wgts
! Subprogram not used     wgtn2 => wgt2%wgtn
! Subprogram not used 
! Subprogram not used     iim => wgt1%jjm
! Subprogram not used     iip => wgt1%jjp
! Subprogram not used     wgts1 => wgt1%wgts
! Subprogram not used     wgtn1 => wgt1%wgtn
! Subprogram not used 
! Subprogram not used     do j=1,m2
! Subprogram not used       do i=1,n1
! Subprogram not used         arrtmp(i,j) = arrin(i,jjm(j))*wgts2(j) + arrin(i,jjp(j))*wgtn2(j)
! Subprogram not used       end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     do j=1,m2
! Subprogram not used       do i=1,m1
! Subprogram not used         arrout(i,j) = arrtmp(iim(i),j)*wgts1(i) + arrtmp(iip(i),j)*wgtn1(i)
! Subprogram not used       end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used   end subroutine lininterp2d2d
  subroutine lininterp2d1d(arrin, n1, n2, arrout, m1, wgt1, wgt2, fldname)
    implicit none
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    integer, intent(in) :: n1, n2, m1
    real(r8), intent(in) :: arrin(n1,n2)    ! input array of values to interpolate
    type(interp_type), intent(in) :: wgt1, wgt2
    real(r8), intent(out) :: arrout(m1) ! interpolated array
    character(len=*), intent(in), optional :: fldname(:)
    !
    ! locals
    !
    integer i                ! indices
    integer, pointer :: iim(:), jjm(:)
    integer, pointer :: iip(:), jjp(:)

    real(r8), pointer :: wgts(:), wgte(:)
    real(r8), pointer :: wgtn(:), wgtw(:)

    jjm => wgt2%jjm
    jjp => wgt2%jjp
    wgts => wgt2%wgts
    wgtn => wgt2%wgtn

    iim => wgt1%jjm
    iip => wgt1%jjp
    wgtw => wgt1%wgts
    wgte => wgt1%wgtn

    do i=1,m1
       arrout(i) = arrin(iim(i),jjm(i))*wgtw(i)*wgts(i)+arrin(iip(i),jjm(i))*wgte(i)*wgts(i) + &
                   arrin(iim(i),jjp(i))*wgtw(i)*wgtn(i)+arrin(iip(i),jjp(i))*wgte(i)*wgtn(i)
    end do


  end subroutine lininterp2d1d
  subroutine lininterp3d2d(arrin, n1, n2, n3, arrout, m1, len1, wgt1, wgt2)
    implicit none
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    integer, intent(in) :: n1, n2, n3, m1, len1   ! m1 is to len1 as ncols is to pcols
    real(r8), intent(in) :: arrin(n1,n2,n3)    ! input array of values to interpolate
    type(interp_type), intent(in) :: wgt1, wgt2
    real(r8), intent(out) :: arrout(len1, n3) ! interpolated array

    !
    ! locals
    !
    integer i, k               ! indices
    integer, pointer :: iim(:), jjm(:)
    integer, pointer :: iip(:), jjp(:)

    real(r8), pointer :: wgts(:), wgte(:)
    real(r8), pointer :: wgtn(:), wgtw(:)

    jjm => wgt2%jjm
    jjp => wgt2%jjp
    wgts => wgt2%wgts
    wgtn => wgt2%wgtn

    iim => wgt1%jjm
    iip => wgt1%jjp
    wgtw => wgt1%wgts
    wgte => wgt1%wgtn

    do k=1,n3
       do i=1,m1
          arrout(i,k) = arrin(iim(i),jjm(i),k)*wgtw(i)*wgts(i)+arrin(iip(i),jjm(i),k)*wgte(i)*wgts(i) + &
               arrin(iim(i),jjp(i),k)*wgtw(i)*wgtn(i)+arrin(iip(i),jjp(i),k)*wgte(i)*wgtn(i)
       end do
    end do

  end subroutine lininterp3d2d




  subroutine lininterp_finish(interp_wgts)
    type(interp_type) :: interp_wgts

    deallocate(interp_wgts%jjm, &
         interp_wgts%jjp, &
         interp_wgts%wgts, &
         interp_wgts%wgtn)

    nullify(interp_wgts%jjm, &
         interp_wgts%jjp, &
         interp_wgts%wgts, &
         interp_wgts%wgtn)
  end subroutine lininterp_finish

! Subprogram not used   subroutine lininterp_original (arrin, yin, nlev, nlatin, arrout, &
! Subprogram not used        yout, nlatout)
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     ! Purpose: Do a linear interpolation from input mesh defined by yin to output
! Subprogram not used     !          mesh defined by yout.  Where extrapolation is necessary, values will
! Subprogram not used     !          be copied from the extreme edge of the input grid.  Vectorization is over
! Subprogram not used     !          the vertical (nlev) dimension.
! Subprogram not used     !
! Subprogram not used     ! Method: Check validity of input, then determine weights, then do the N-S interpolation.
! Subprogram not used     !
! Subprogram not used     ! Author: Jim Rosinski
! Subprogram not used     ! Modified: Jim Edwards so that there is no requirement of monotonacity on the yout array
! Subprogram not used     !
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     implicit none
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     ! Arguments
! Subprogram not used     !
! Subprogram not used     integer, intent(in) :: nlev                   ! number of vertical levels
! Subprogram not used     integer, intent(in) :: nlatin                 ! number of input latitudes
! Subprogram not used     integer, intent(in) :: nlatout                ! number of output latitudes
! Subprogram not used 
! Subprogram not used     real(r8), intent(in) :: arrin(nlev,nlatin)    ! input array of values to interpolate
! Subprogram not used     real(r8), intent(in) :: yin(nlatin)           ! input mesh
! Subprogram not used     real(r8), intent(in) :: yout(nlatout)         ! output mesh
! Subprogram not used 
! Subprogram not used     real(r8), intent(out) :: arrout(nlev,nlatout) ! interpolated array
! Subprogram not used     !
! Subprogram not used     ! Local workspace
! Subprogram not used     !
! Subprogram not used     integer j, jj              ! latitude indices
! Subprogram not used     integer jjprev             ! latitude indices
! Subprogram not used     integer k                  ! level index
! Subprogram not used     integer icount             ! number of values
! Subprogram not used 
! Subprogram not used     real(r8) extrap            ! percent grid non-overlap
! Subprogram not used     !
! Subprogram not used     ! Dynamic
! Subprogram not used     !
! Subprogram not used     integer :: jjm(nlatout)
! Subprogram not used     integer :: jjp(nlatout)
! Subprogram not used 
! Subprogram not used     real(r8) :: wgts(nlatout)
! Subprogram not used     real(r8) :: wgtn(nlatout)
! Subprogram not used     !
! Subprogram not used     ! Check validity of input coordinate arrays: must be monotonically increasing,
! Subprogram not used     ! and have a total of at least 2 elements
! Subprogram not used     !
! Subprogram not used     if (nlatin.lt.2) then
! Subprogram not used        call endrun('LININTERP: Must have at least 2 input points for interpolation')
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     icount = 0
! Subprogram not used     do j=1,nlatin-1
! Subprogram not used        if (yin(j).gt.yin(j+1)) icount = icount + 1
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     if (icount.gt.0) then
! Subprogram not used        call endrun('LININTERP: Non-monotonic coordinate array(s) found')
! Subprogram not used     end if
! Subprogram not used     !
! Subprogram not used     ! Initialize index arrays for later checking
! Subprogram not used     !
! Subprogram not used     do j=1,nlatout
! Subprogram not used        jjm(j) = 0
! Subprogram not used        jjp(j) = 0
! Subprogram not used     end do
! Subprogram not used     !
! Subprogram not used     ! For values which extend beyond N and S boundaries, set weights
! Subprogram not used     ! such that values will just be copied.
! Subprogram not used     !
! Subprogram not used     extrap = 0._r8
! Subprogram not used 
! Subprogram not used     do j=1,nlatout
! Subprogram not used        if (yout(j).le.yin(1)) then
! Subprogram not used           jjm(j) = 1
! Subprogram not used           jjp(j) = 1
! Subprogram not used           wgts(j) = 1._r8
! Subprogram not used           wgtn(j) = 0._r8
! Subprogram not used           extrap=extrap+1._r8
! Subprogram not used        else if (yout(j).gt.yin(nlatin)) then
! Subprogram not used           jjm(j) = nlatin
! Subprogram not used           jjp(j) = nlatin
! Subprogram not used           wgts(j) = 1._r8
! Subprogram not used           wgtn(j) = 0._r8
! Subprogram not used           extrap=extrap+1._r8
! Subprogram not used        endif
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     !
! Subprogram not used     ! Loop though output indices finding input indices and weights
! Subprogram not used     !
! Subprogram not used     do j=1,nlatout
! Subprogram not used        do jj=1,nlatin-1
! Subprogram not used           if (yout(j).gt.yin(jj) .and. yout(j).le.yin(jj+1)) then
! Subprogram not used              jjm(j) = jj
! Subprogram not used              jjp(j) = jj + 1
! Subprogram not used              wgts(j) = (yin(jj+1)-yout(j))/(yin(jj+1)-yin(jj))
! Subprogram not used              wgtn(j) = (yout(j)-yin(jj))/(yin(jj+1)-yin(jj))
! Subprogram not used              exit
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used     !
! Subprogram not used     ! Check that interp/extrap points have been found for all outputs
! Subprogram not used     !
! Subprogram not used     icount = 0
! Subprogram not used     do j=1,nlatout
! Subprogram not used        if (jjm(j).eq.0 .or. jjp(j).eq.0) then
! Subprogram not used           icount = icount + 1
! Subprogram not used        end if
! Subprogram not used     end do
! Subprogram not used     if (icount.gt.0) then
! Subprogram not used        call endrun('LININTERP: Point found without interp indices')
! Subprogram not used     end if
! Subprogram not used     !
! Subprogram not used     ! Do the interpolation
! Subprogram not used     !
! Subprogram not used     do j=1,nlatout
! Subprogram not used        do k=1,nlev
! Subprogram not used           arrout(k,j) = arrin(k,jjm(j))*wgts(j) + arrin(k,jjp(j))*wgtn(j)
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     return
! Subprogram not used   end subroutine lininterp_original


! Subprogram not used   subroutine bilin (arrin, xin, yin, nlondin, nlonin, &
! Subprogram not used        nlevdin, nlev, nlatin, arrout, xout, &
! Subprogram not used        yout, nlondout, nlonout, nlevdout, nlatout)
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     ! Purpose:
! Subprogram not used     !
! Subprogram not used     ! Do a bilinear interpolation from input mesh defined by xin, yin to output
! Subprogram not used     ! mesh defined by xout, yout.  Circularity is assumed in the x-direction so
! Subprogram not used     ! input x-grid must be in degrees east and must start from Greenwich.  When
! Subprogram not used     ! extrapolation is necessary in the N-S direction, values will be copied
! Subprogram not used     ! from the extreme edge of the input grid.  Vectorization is over the
! Subprogram not used     ! longitude dimension.  Input grid is assumed rectangular. Output grid
! Subprogram not used     ! is assumed ragged, i.e. xout is a 2-d mesh.
! Subprogram not used     !
! Subprogram not used     ! Method: Interpolate first in longitude, then in latitude.
! Subprogram not used     !
! Subprogram not used     ! Author: Jim Rosinski
! Subprogram not used     !
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     use shr_kind_mod, only: r8 => shr_kind_r8
! Subprogram not used     use abortutils,   only: endrun
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     implicit none
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     ! Input arguments
! Subprogram not used     !
! Subprogram not used     integer, intent(in) :: nlondin                        ! longitude dimension of input grid
! Subprogram not used     integer, intent(in) :: nlonin                         ! number of real longitudes (input)
! Subprogram not used     integer, intent(in) :: nlevdin                        ! vertical dimension of input grid
! Subprogram not used     integer, intent(in) :: nlev                           ! number of vertical levels
! Subprogram not used     integer, intent(in) :: nlatin                         ! number of input latitudes
! Subprogram not used     integer, intent(in) :: nlatout                        ! number of output latitudes
! Subprogram not used     integer, intent(in) :: nlondout                       ! longitude dimension of output grid
! Subprogram not used     integer, intent(in) :: nlonout(nlatout)               ! number of output longitudes per lat
! Subprogram not used     integer, intent(in) :: nlevdout                       ! vertical dimension of output grid
! Subprogram not used 
! Subprogram not used     real(r8), intent(in) :: arrin(nlondin,nlevdin,nlatin) ! input array of values to interpolate
! Subprogram not used     real(r8), intent(in) :: xin(nlondin)                  ! input x mesh
! Subprogram not used     real(r8), intent(in) :: yin(nlatin)                   ! input y mesh
! Subprogram not used     real(r8), intent(in) :: xout(nlondout,nlatout)        ! output x mesh
! Subprogram not used     real(r8), intent(in) :: yout(nlatout)                 ! output y mesh
! Subprogram not used     !
! Subprogram not used     ! Output arguments
! Subprogram not used     !
! Subprogram not used     real(r8), intent(out) :: arrout(nlondout,nlevdout,nlatout) ! interpolated array
! Subprogram not used     !
! Subprogram not used     ! Local workspace
! Subprogram not used     !
! Subprogram not used     integer :: i, ii, iw, ie, iiprev ! longitude indices
! Subprogram not used     integer :: j, jj, js, jn, jjprev ! latitude indices
! Subprogram not used     integer :: k                     ! level index
! Subprogram not used     integer :: icount                ! number of bad values
! Subprogram not used 
! Subprogram not used     real(r8) :: extrap               ! percent grid non-overlap
! Subprogram not used     real(r8) :: dxinwrap             ! delta-x on input grid for 2-pi
! Subprogram not used     real(r8) :: avgdxin              ! avg input delta-x
! Subprogram not used     real(r8) :: ratio                ! compare dxinwrap to avgdxin
! Subprogram not used     real(r8) :: sum                  ! sum of weights (used for testing)
! Subprogram not used     !
! Subprogram not used     ! Dynamic
! Subprogram not used     !
! Subprogram not used     integer :: iim(nlondout)         ! interpolation index to the left
! Subprogram not used     integer :: iip(nlondout)         ! interpolation index to the right
! Subprogram not used     integer :: jjm(nlatout)          ! interpolation index to the south
! Subprogram not used     integer :: jjp(nlatout)          ! interpolation index to the north
! Subprogram not used 
! Subprogram not used     real(r8) :: wgts(nlatout)        ! interpolation weight to the north
! Subprogram not used     real(r8) :: wgtn(nlatout)        ! interpolation weight to the north
! Subprogram not used     real(r8) :: wgte(nlondout)       ! interpolation weight to the north
! Subprogram not used     real(r8) :: wgtw(nlondout)       ! interpolation weight to the north
! Subprogram not used     real(r8) :: igrid(nlonin)        ! interpolation weight to the north
! Subprogram not used     !
! Subprogram not used     ! Check validity of input coordinate arrays: must be monotonically increasing,
! Subprogram not used     ! and have a total of at least 2 elements
! Subprogram not used     !
! Subprogram not used     if (nlonin < 2 .or. nlatin < 2) then
! Subprogram not used        call endrun ('BILIN: Must have at least 2 input points for interpolation')
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     if (xin(1) < 0._r8 .or. xin(nlonin) > 360._r8) then
! Subprogram not used        call endrun ('BILIN: Input x-grid must be between 0 and 360')
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     icount = 0
! Subprogram not used     do i=1,nlonin-1
! Subprogram not used        if (xin(i) >= xin(i+1)) icount = icount + 1
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     do j=1,nlatin-1
! Subprogram not used        if (yin(j) >= yin(j+1)) icount = icount + 1
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     do j=1,nlatout-1
! Subprogram not used        if (yout(j) >= yout(j+1)) icount = icount + 1
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     do j=1,nlatout
! Subprogram not used        do i=1,nlonout(j)-1
! Subprogram not used           if (xout(i,j) >= xout(i+1,j)) icount = icount + 1
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     if (icount > 0) then
! Subprogram not used        call endrun ('BILIN: Non-monotonic coordinate array(s) found')
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     if (yout(nlatout) <= yin(1) .or. yout(1) >= yin(nlatin)) then
! Subprogram not used        call endrun ('BILIN: No overlap in y-coordinates')
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     do j=1,nlatout
! Subprogram not used        if (xout(1,j) < 0._r8 .or. xout(nlonout(j),j) > 360._r8) then
! Subprogram not used           call endrun ('BILIN: Output x-grid must be between 0 and 360')
! Subprogram not used        end if
! Subprogram not used 
! Subprogram not used        if (xout(nlonout(j),j) <= xin(1) .or.  &
! Subprogram not used             xout(1,j)          >= xin(nlonin)) then
! Subprogram not used           call endrun ('BILIN: No overlap in x-coordinates')
! Subprogram not used        end if
! Subprogram not used     end do
! Subprogram not used     !
! Subprogram not used     ! Initialize index arrays for later checking
! Subprogram not used     !
! Subprogram not used     do j=1,nlatout
! Subprogram not used        jjm(j) = 0
! Subprogram not used        jjp(j) = 0
! Subprogram not used     end do
! Subprogram not used     !
! Subprogram not used     ! For values which extend beyond N and S boundaries, set weights
! Subprogram not used     ! such that values will just be copied.
! Subprogram not used     !
! Subprogram not used     do js=1,nlatout
! Subprogram not used        if (yout(js) > yin(1)) exit
! Subprogram not used        jjm(js) = 1
! Subprogram not used        jjp(js) = 1
! Subprogram not used        wgts(js) = 1._r8
! Subprogram not used        wgtn(js) = 0._r8
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     do jn=nlatout,1,-1
! Subprogram not used        if (yout(jn) <= yin(nlatin)) exit
! Subprogram not used        jjm(jn) = nlatin
! Subprogram not used        jjp(jn) = nlatin
! Subprogram not used        wgts(jn) = 1._r8
! Subprogram not used        wgtn(jn) = 0._r8
! Subprogram not used     end do
! Subprogram not used     !
! Subprogram not used     ! Loop though output indices finding input indices and weights
! Subprogram not used     !
! Subprogram not used     jjprev = 1
! Subprogram not used     do j=js,jn
! Subprogram not used        do jj=jjprev,nlatin-1
! Subprogram not used           if (yout(j) > yin(jj) .and. yout(j) <= yin(jj+1)) then
! Subprogram not used              jjm(j) = jj
! Subprogram not used              jjp(j) = jj + 1
! Subprogram not used              wgts(j) = (yin(jj+1) - yout(j)) / (yin(jj+1) - yin(jj))
! Subprogram not used              wgtn(j) = (yout(j)   - yin(jj)) / (yin(jj+1) - yin(jj))
! Subprogram not used              goto 30
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used        call endrun ('BILIN: Failed to find interp values')
! Subprogram not used 30     jjprev = jj
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     dxinwrap = xin(1) + 360._r8 - xin(nlonin)
! Subprogram not used     !
! Subprogram not used     ! Check for sane dxinwrap values.  Allow to differ no more than 10% from avg
! Subprogram not used     !
! Subprogram not used     avgdxin = (xin(nlonin)-xin(1))/(nlonin-1._r8)
! Subprogram not used     ratio = dxinwrap/avgdxin
! Subprogram not used     if (ratio < 0.9_r8 .or. ratio > 1.1_r8) then
! Subprogram not used        write(iulog,*)'BILIN: Insane dxinwrap value =',dxinwrap,' avg=', avgdxin
! Subprogram not used        call endrun
! Subprogram not used     end if
! Subprogram not used     !
! Subprogram not used     ! Check that interp/extrap points have been found for all outputs, and that
! Subprogram not used     ! they are reasonable (i.e. within 32-bit roundoff)
! Subprogram not used     !
! Subprogram not used     icount = 0
! Subprogram not used     do j=1,nlatout
! Subprogram not used        if (jjm(j) == 0 .or. jjp(j) == 0) icount = icount + 1
! Subprogram not used        sum = wgts(j) + wgtn(j)
! Subprogram not used        if (sum < 0.99999_r8 .or. sum > 1.00001_r8) icount = icount + 1
! Subprogram not used        if (wgts(j) < 0._r8 .or. wgts(j) > 1._r8) icount = icount + 1
! Subprogram not used        if (wgtn(j) < 0._r8 .or. wgtn(j) > 1._r8) icount = icount + 1
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     if (icount > 0) then
! Subprogram not used        call endrun ('BILIN: something bad in latitude indices or weights')
! Subprogram not used     end if
! Subprogram not used     !
! Subprogram not used     ! Do the bilinear interpolation
! Subprogram not used     !
! Subprogram not used     do j=1,nlatout
! Subprogram not used        !
! Subprogram not used        ! Initialize index arrays for later checking
! Subprogram not used        !
! Subprogram not used        do i=1,nlondout
! Subprogram not used           iim(i) = 0
! Subprogram not used           iip(i) = 0
! Subprogram not used        end do
! Subprogram not used        !
! Subprogram not used        ! For values which extend beyond E and W boundaries, set weights such that
! Subprogram not used        ! values will be interpolated between E and W edges of input grid.
! Subprogram not used        !
! Subprogram not used        do iw=1,nlonout(j)
! Subprogram not used           if (xout(iw,j) > xin(1)) exit
! Subprogram not used           iim(iw) = nlonin
! Subprogram not used           iip(iw) = 1
! Subprogram not used           wgtw(iw) = (xin(1)        - xout(iw,j))   /dxinwrap
! Subprogram not used           wgte(iw) = (xout(iw,j)+360._r8 - xin(nlonin))/dxinwrap
! Subprogram not used        end do
! Subprogram not used 
! Subprogram not used        do ie=nlonout(j),1,-1
! Subprogram not used           if (xout(ie,j) <= xin(nlonin)) exit
! Subprogram not used           iim(ie) = nlonin
! Subprogram not used           iip(ie) = 1
! Subprogram not used           wgtw(ie) = (xin(1)+360._r8 - xout(ie,j))   /dxinwrap
! Subprogram not used           wgte(ie) = (xout(ie,j)    - xin(nlonin))/dxinwrap
! Subprogram not used        end do
! Subprogram not used        !
! Subprogram not used        ! Loop though output indices finding input indices and weights
! Subprogram not used        !
! Subprogram not used        iiprev = 1
! Subprogram not used        do i=iw,ie
! Subprogram not used           do ii=iiprev,nlonin-1
! Subprogram not used              if (xout(i,j) > xin(ii) .and. xout(i,j) <= xin(ii+1)) then
! Subprogram not used                 iim(i) = ii
! Subprogram not used                 iip(i) = ii + 1
! Subprogram not used                 wgtw(i) = (xin(ii+1) - xout(i,j)) / (xin(ii+1) - xin(ii))
! Subprogram not used                 wgte(i) = (xout(i,j) - xin(ii))   / (xin(ii+1) - xin(ii))
! Subprogram not used                 goto 60
! Subprogram not used              end if
! Subprogram not used           end do
! Subprogram not used           call endrun ('BILIN: Failed to find interp values')
! Subprogram not used 60        iiprev = ii
! Subprogram not used        end do
! Subprogram not used 
! Subprogram not used        icount = 0
! Subprogram not used        do i=1,nlonout(j)
! Subprogram not used           if (iim(i) == 0 .or. iip(i) == 0) icount = icount + 1
! Subprogram not used           sum = wgtw(i) + wgte(i)
! Subprogram not used           if (sum < 0.99999_r8 .or. sum > 1.00001_r8) icount = icount + 1
! Subprogram not used           if (wgtw(i) < 0._r8 .or. wgtw(i) > 1._r8) icount = icount + 1
! Subprogram not used           if (wgte(i) < 0._r8 .or. wgte(i) > 1._r8) icount = icount + 1
! Subprogram not used        end do
! Subprogram not used 
! Subprogram not used        if (icount > 0) then
! Subprogram not used           write(iulog,*)'BILIN: j=',j,' Something bad in longitude indices or weights'
! Subprogram not used           call endrun
! Subprogram not used        end if
! Subprogram not used        !
! Subprogram not used        ! Do the interpolation, 1st in longitude then latitude
! Subprogram not used        !
! Subprogram not used        do k=1,nlev
! Subprogram not used           do i=1,nlonin
! Subprogram not used              igrid(i) = arrin(i,k,jjm(j))*wgts(j) + arrin(i,k,jjp(j))*wgtn(j)
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used           do i=1,nlonout(j)
! Subprogram not used              arrout(i,k,j) = igrid(iim(i))*wgtw(i) + igrid(iip(i))*wgte(i)
! Subprogram not used           end do
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     return
! Subprogram not used   end subroutine bilin

  subroutine vertinterp(ncol, ncold, nlev, pmid, pout, arrin, arrout)

    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Vertically interpolate input array to output pressure level
    ! Copy values at boundaries.
    !
    ! Method:
    !
    ! Author:
    !
    !-----------------------------------------------------------------------

    implicit none

    !------------------------------Arguments--------------------------------
    integer , intent(in)  :: ncol              ! column dimension
    integer , intent(in)  :: ncold             ! declared column dimension
    integer , intent(in)  :: nlev              ! vertical dimension
    real(r8), intent(in)  :: pmid(ncold,nlev)  ! input level pressure levels
    real(r8), intent(in)  :: pout              ! output pressure level
    real(r8), intent(in)  :: arrin(ncold,nlev) ! input  array
    real(r8), intent(out) :: arrout(ncold)     ! output array (interpolated)
    !--------------------------------------------------------------------------

    !---------------------------Local variables-----------------------------
    integer i,k               ! indices
    integer kupper(ncold)     ! Level indices for interpolation
    real(r8) dpu              ! upper level pressure difference
    real(r8) dpl              ! lower level pressure difference
    logical found(ncold)      ! true if input levels found
    logical error             ! error flag
    !-----------------------------------------------------------------
    !
    ! Initialize index array and logical flags
    !
    do i=1,ncol
       found(i)  = .false.
       kupper(i) = 1
    end do
    error = .false.
    !
    ! Store level indices for interpolation.
    ! If all indices for this level have been found,
    ! do the interpolation
    !
    do k=1,nlev-1
       do i=1,ncol
          if ((.not. found(i)) .and. pmid(i,k)<pout .and. pout<=pmid(i,k+1)) then
             found(i) = .true.
             kupper(i) = k
          end if
       end do
    end do
    !
    ! If we've fallen through the k=1,nlev-1 loop, we cannot interpolate and
    ! must extrapolate from the bottom or top data level for at least some
    ! of the longitude points.
    !
    do i=1,ncol
       if (pout <= pmid(i,1)) then
          arrout(i) = arrin(i,1)
       else if (pout >= pmid(i,nlev)) then
          arrout(i) = arrin(i,nlev)
       else if (found(i)) then
          dpu = pout - pmid(i,kupper(i))
          dpl = pmid(i,kupper(i)+1) - pout
          arrout(i) = (arrin(i,kupper(i)  )*dpl + arrin(i,kupper(i)+1)*dpu)/(dpl + dpu)
       else
          error = .true.
       end if
    end do
    !
    ! Error check
    !
    if (error) then
       call endrun ('VERTINTERP: ERROR FLAG')
    end if

    return
  end subroutine vertinterp

! Subprogram not used   subroutine get_timeinterp_factors (cycflag, np1, cdayminus, cdayplus, cday, &
! Subprogram not used        fact1, fact2, str)
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     ! Purpose: Determine time interpolation factors (normally for a boundary dataset)
! Subprogram not used     !          for linear interpolation.
! Subprogram not used     !
! Subprogram not used     ! Method:  Assume 365 days per year.  Output variable fact1 will be the weight to
! Subprogram not used     !          apply to data at calendar time "cdayminus", and fact2 the weight to apply
! Subprogram not used     !          to data at time "cdayplus".  Combining these values will produce a result
! Subprogram not used     !          valid at time "cday".  Output arguments fact1 and fact2 will be between
! Subprogram not used     !          0 and 1, and fact1 + fact2 = 1 to roundoff.
! Subprogram not used     !
! Subprogram not used     ! Author:  Jim Rosinski
! Subprogram not used     !
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     implicit none
! Subprogram not used     !
! Subprogram not used     ! Arguments
! Subprogram not used     !
! Subprogram not used     logical, intent(in) :: cycflag             ! flag indicates whether dataset is being cycled yearly
! Subprogram not used 
! Subprogram not used     integer, intent(in) :: np1                 ! index points to forward time slice matching cdayplus
! Subprogram not used 
! Subprogram not used     real(r8), intent(in) :: cdayminus          ! calendar day of rearward time slice
! Subprogram not used     real(r8), intent(in) :: cdayplus           ! calendar day of forward time slice
! Subprogram not used     real(r8), intent(in) :: cday               ! calenar day to be interpolated to
! Subprogram not used     real(r8), intent(out) :: fact1             ! time interpolation factor to apply to rearward time slice
! Subprogram not used     real(r8), intent(out) :: fact2             ! time interpolation factor to apply to forward time slice
! Subprogram not used 
! Subprogram not used     character(len=*), intent(in) :: str        ! string to be added to print in case of error (normally the callers name)
! Subprogram not used     !
! Subprogram not used     ! Local workspace
! Subprogram not used     !
! Subprogram not used     real(r8) :: deltat                         ! time difference (days) between cdayminus and cdayplus
! Subprogram not used     real(r8), parameter :: daysperyear = 365._r8  ! number of days in a year
! Subprogram not used     !
! Subprogram not used     ! Initial sanity checks
! Subprogram not used     !
! Subprogram not used     if (np1 == 1 .and. .not. cycflag) then
! Subprogram not used        call endrun ('GETFACTORS:'//str//' cycflag false and forward month index = Jan. not allowed')
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     if (np1 < 1) then
! Subprogram not used        call endrun ('GETFACTORS:'//str//' input arg np1 must be > 0')
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     if (cycflag) then
! Subprogram not used        if ((cday < 1._r8) .or. (cday > (daysperyear+1._r8))) then
! Subprogram not used           write(iulog,*) 'GETFACTORS:', str, ' bad cday=',cday
! Subprogram not used           call endrun ()
! Subprogram not used        end if
! Subprogram not used     else
! Subprogram not used        if (cday < 1._r8) then
! Subprogram not used           write(iulog,*) 'GETFACTORS:', str, ' bad cday=',cday
! Subprogram not used           call endrun ()
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used     !
! Subprogram not used     ! Determine time interpolation factors.  Account for December-January
! Subprogram not used     ! interpolation if dataset is being cycled yearly.
! Subprogram not used     !
! Subprogram not used     if (cycflag .and. np1 == 1) then                     ! Dec-Jan interpolation
! Subprogram not used        deltat = cdayplus + daysperyear - cdayminus
! Subprogram not used        if (cday > cdayplus) then                         ! We are in December
! Subprogram not used           fact1 = (cdayplus + daysperyear - cday)/deltat
! Subprogram not used           fact2 = (cday - cdayminus)/deltat
! Subprogram not used        else                                              ! We are in January
! Subprogram not used           fact1 = (cdayplus - cday)/deltat
! Subprogram not used           fact2 = (cday + daysperyear - cdayminus)/deltat
! Subprogram not used        end if
! Subprogram not used     else
! Subprogram not used        deltat = cdayplus - cdayminus
! Subprogram not used        fact1 = (cdayplus - cday)/deltat
! Subprogram not used        fact2 = (cday - cdayminus)/deltat
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     if (.not. valid_timeinterp_factors (fact1, fact2)) then
! Subprogram not used        write(iulog,*) 'GETFACTORS: ', str, ' bad fact1 and/or fact2=', fact1, fact2
! Subprogram not used        call endrun ()
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     return
! Subprogram not used   end subroutine get_timeinterp_factors

! Subprogram not used   logical function valid_timeinterp_factors (fact1, fact2)
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     ! Purpose: check sanity of time interpolation factors to within 32-bit roundoff
! Subprogram not used     !
! Subprogram not used     !---------------------------------------------------------------------------
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     real(r8), intent(in) :: fact1, fact2           ! time interpolation factors
! Subprogram not used 
! Subprogram not used     valid_timeinterp_factors = .true.
! Subprogram not used 
! Subprogram not used     ! The fact1 .ne. fact1 and fact2 .ne. fact2 comparisons are to detect NaNs.
! Subprogram not used     if (abs(fact1+fact2-1._r8) > 1.e-6_r8 .or. &
! Subprogram not used          fact1 > 1.000001_r8 .or. fact1 < -1.e-6_r8 .or. &
! Subprogram not used          fact2 > 1.000001_r8 .or. fact2 < -1.e-6_r8 .or. &
! Subprogram not used          fact1 .ne. fact1 .or. fact2 .ne. fact2) then
! Subprogram not used 
! Subprogram not used        valid_timeinterp_factors = .false.
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     return
! Subprogram not used   end function valid_timeinterp_factors

end module interpolate_data
