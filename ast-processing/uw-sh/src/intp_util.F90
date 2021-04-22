
! Subprogram not used subroutine chktime( time, nrec )
! Subprogram not used 
! Subprogram not used    !----------------------------------------------------------------------- 
! Subprogram not used    ! Purpose: 
! Subprogram not used    ! Make sure the time coordinate looks like calander day, and is increasing.
! Subprogram not used    ! Calendar day can either start with 1 Jan 0Z = day 1.0  or
! Subprogram not used    !  1 Jan 0Z = day 0.0
! Subprogram not used    !
! Subprogram not used    ! Author: B. Eaton
! Subprogram not used    !----------------------------------------------------------------------- 
! Subprogram not used 
! Subprogram not used    use shr_kind_mod, only: r8 => shr_kind_r8
! Subprogram not used    use abortutils,   only: endrun
! Subprogram not used    use cam_logfile,  only: iulog
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used    integer, intent(in) :: nrec                 ! size of time array
! Subprogram not used    real(r8), intent(in) :: time(nrec)           ! time coordinate expressed as calendar day.
! Subprogram not used    
! Subprogram not used    ! Local varibles:
! Subprogram not used    integer :: i
! Subprogram not used    !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    if ( time(1) .lt. 0._r8 .or. time(1) .ge. 366._r8 ) then
! Subprogram not used       write(iulog,*)'chktime: illegal time coordinate ',time(1)
! Subprogram not used       call endrun
! Subprogram not used    end if
! Subprogram not used    do i = 2, nrec
! Subprogram not used       if ( time(i) .lt. 0._r8 .or. time(i) .ge. 366._r8 ) then
! Subprogram not used          write(iulog,*)'chktime: illegal time coordinate ', time(i)
! Subprogram not used          call endrun
! Subprogram not used       end if
! Subprogram not used       if ( time(i) .le. time(i-1) ) then
! Subprogram not used          write(iulog,*)'chktime: time not increasing ', time(i-1), time(i)
! Subprogram not used          call endrun
! Subprogram not used       end if
! Subprogram not used    end do
! Subprogram not used 
! Subprogram not used    return
! Subprogram not used 
! Subprogram not used end subroutine chktime

!#######################################################################

subroutine findplb( x, nx, xval, index )

   !----------------------------------------------------------------------- 
   ! Purpose: 
   ! "find periodic lower bound"
   ! Search the input array for the lower bound of the interval that
   ! contains the input value.  The returned index satifies:
   ! x(index) .le. xval .lt. x(index+1)
   ! Assume the array represents values in one cycle of a periodic coordinate.
   ! So, if xval .lt. x(1), or xval .ge. x(nx), then the index returned is nx.
   !
   ! Author: B. Eaton
   !----------------------------------------------------------------------- 

   use shr_kind_mod, only: r8 => shr_kind_r8

   implicit none

   integer, intent(in) ::   nx         ! size of x
   real(r8), intent(in) ::  x(nx)      ! strictly increasing array
   real(r8), intent(in) ::  xval       ! value to be searched for in x
   
   integer, intent(out) ::  index

   ! Local variables:
   integer i
   !-----------------------------------------------------------------------

   if ( xval .lt. x(1) .or. xval .ge. x(nx) ) then
      index = nx
      return
   end if

   do i = 2, nx
      if ( xval .lt. x(i) ) then
         index = i-1
         return
      end if
   end do

end subroutine findplb

!#######################################################################

! Subprogram not used subroutine linintp( npts, t1, t2, tint, f1, f2, fint )
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    !-----------------------------------------------------------------------
! Subprogram not used    ! Purpose:
! Subprogram not used    ! Linearly interpolate between f1(t1) and f2(t2) to fint(tint),
! Subprogram not used    ! where f1, f2, and f3 are chunked data structures
! Subprogram not used    !
! Subprogram not used    ! Author: B. Eaton
! Subprogram not used    ! Chunked by P. Worley
! Subprogram not used    ! Un-Chunked by P. Rasch (it wasnt right)
! Subprogram not used    !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used ! Linearly interpolate between f1(t1) and f2(t2) to fint(tint).
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    use shr_kind_mod, only: r8 => shr_kind_r8
! Subprogram not used    use ppgrid
! Subprogram not used    use phys_grid, only: get_ncols_p
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! Input arguments:
! Subprogram not used       integer,   intent(in)::   npts
! Subprogram not used       real(r8),  intent(in)::   t1               ! time level of f1
! Subprogram not used       real(r8),  intent(in)::   t2               ! time level of f2
! Subprogram not used       real(r8),  intent(in)::   tint             ! interpolant time
! Subprogram not used       real(r8),  intent(in)::   f1(npts)         ! field at time t1
! Subprogram not used       real(r8),  intent(in)::   f2(npts)         ! field at time t2
! Subprogram not used 
! Subprogram not used ! Output arguments:
! Subprogram not used       real(r8), intent(out)::  fint(npts)       ! field at time tint
! Subprogram not used 
! Subprogram not used ! Local variables:
! Subprogram not used       integer i
! Subprogram not used       real(r8) factor
! Subprogram not used !------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used !      call t_startf ('linintp')
! Subprogram not used 
! Subprogram not used       factor = ( tint - t1 )/( t2 - t1)
! Subprogram not used 
! Subprogram not used       do i = 1, npts
! Subprogram not used          fint(i) = f1(i) + ( f2(i) - f1(i) )*factor
! Subprogram not used       end do
! Subprogram not used 
! Subprogram not used !      call t_stopf ('linintp')
! Subprogram not used 
! Subprogram not used       return
! Subprogram not used end subroutine linintp

!#######################################################################
