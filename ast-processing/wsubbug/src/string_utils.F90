module string_utils


   implicit none
   private

! Public interface methods

   public ::&
      to_upper, &   ! Convert character string to upper case
      to_lower, &   ! Convert character string to lower case
      INCSTR, &     ! increments a string
      GLC           ! Position of last significant character in string

contains

function to_upper(str)

!----------------------------------------------------------------------- 
! Purpose: 
! Convert character string to upper case.
! 
! Method: 
! Use achar and iachar intrinsics to ensure use of ascii collating sequence.
!
! Author:  B. Eaton, July 2001
!     
! $Id$
!----------------------------------------------------------------------- 
   implicit none

   character(len=*), intent(in) :: str      ! String to convert to upper case
   character(len=len(str))      :: to_upper

! Local variables

   integer :: i                ! Index
   integer :: aseq             ! ascii collating sequence
   integer :: lower_to_upper   ! integer to convert case
   character(len=1) :: ctmp    ! Character temporary
!-----------------------------------------------------------------------
   lower_to_upper = iachar("A") - iachar("a")

   do i = 1, len(str)
      ctmp = str(i:i)
      aseq = iachar(ctmp)
      if ( aseq >= iachar("a") .and. aseq <= iachar("z") ) &
           ctmp = achar(aseq + lower_to_upper)
      to_upper(i:i) = ctmp
   end do

end function to_upper

function to_lower(str)

!----------------------------------------------------------------------- 
! Purpose: 
! Convert character string to lower case.
! 
! Method: 
! Use achar and iachar intrinsics to ensure use of ascii collating sequence.
!
! Author:  B. Eaton, July 2001
!     
! $Id$
!----------------------------------------------------------------------- 
   implicit none

   character(len=*), intent(in) :: str      ! String to convert to lower case
   character(len=len(str))      :: to_lower

! Local variables

   integer :: i                ! Index
   integer :: aseq             ! ascii collating sequence
   integer :: upper_to_lower   ! integer to convert case
   character(len=1) :: ctmp    ! Character temporary
!-----------------------------------------------------------------------
   upper_to_lower = iachar("a") - iachar("A")

   do i = 1, len(str)
      ctmp = str(i:i)
      aseq = iachar(ctmp)
      if ( aseq >= iachar("A") .and. aseq <= iachar("Z") ) &
           ctmp = achar(aseq + upper_to_lower)
      to_lower(i:i) = ctmp
   end do

end function to_lower

! Subprogram not used integer function INCSTR( s, inc )
! Subprogram not used   !-----------------------------------------------------------------------
! Subprogram not used   ! 	... Increment a string whose ending characters are digits.
! Subprogram not used   !           The incremented integer must be in the range [0 - (10**n)-1]
! Subprogram not used   !           where n is the number of trailing digits.
! Subprogram not used   !           Return values:
! Subprogram not used   !
! Subprogram not used   !            0 success
! Subprogram not used   !           -1 error: no trailing digits in string
! Subprogram not used   !           -2 error: incremented integer is out of range
! Subprogram not used   !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used   implicit none
! Subprogram not used 
! Subprogram not used   !-----------------------------------------------------------------------
! Subprogram not used   ! 	... Dummy variables
! Subprogram not used   !-----------------------------------------------------------------------
! Subprogram not used   integer, intent(in) :: &
! Subprogram not used        inc                                       ! value to increment string (may be negative)
! Subprogram not used   character(len=*), intent(inout) :: &
! Subprogram not used        s                                         ! string with trailing digits
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   !-----------------------------------------------------------------------
! Subprogram not used   ! 	... Local variables
! Subprogram not used   !-----------------------------------------------------------------------
! Subprogram not used   integer :: &
! Subprogram not used        i, &                          ! index
! Subprogram not used        lstr, &                       ! number of significant characters in string
! Subprogram not used        lnd, &                        ! position of last non-digit
! Subprogram not used        ndigit, &                     ! number of trailing digits
! Subprogram not used        ival, &                       ! integer value of trailing digits
! Subprogram not used        pow, &                        ! power of ten
! Subprogram not used        digit                         ! integer value of a single digit
! Subprogram not used 
! Subprogram not used   lstr   = GLC( s )
! Subprogram not used   lnd    = LASTND( s )
! Subprogram not used   ndigit = lstr - lnd
! Subprogram not used 
! Subprogram not used   if( ndigit == 0 ) then
! Subprogram not used      INCSTR = -1
! Subprogram not used      return
! Subprogram not used   end if
! Subprogram not used 
! Subprogram not used   !-----------------------------------------------------------------------
! Subprogram not used   !     	... Calculate integer corresponding to trailing digits.
! Subprogram not used   !-----------------------------------------------------------------------
! Subprogram not used   ival = 0
! Subprogram not used   pow  = 0
! Subprogram not used   do i = lstr,lnd+1,-1
! Subprogram not used      digit = ICHAR(s(i:i)) - ICHAR('0')
! Subprogram not used      ival  = ival + digit * 10**pow
! Subprogram not used      pow   = pow + 1
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used   !-----------------------------------------------------------------------
! Subprogram not used   !     	... Increment the integer
! Subprogram not used   !-----------------------------------------------------------------------
! Subprogram not used   ival = ival + inc
! Subprogram not used   if( ival < 0 .or. ival > 10**ndigit-1 ) then
! Subprogram not used      INCSTR = -2
! Subprogram not used      return
! Subprogram not used   end if
! Subprogram not used 
! Subprogram not used   !-----------------------------------------------------------------------
! Subprogram not used   !     	... Overwrite trailing digits
! Subprogram not used   !-----------------------------------------------------------------------
! Subprogram not used   pow = ndigit
! Subprogram not used   do i = lnd+1,lstr
! Subprogram not used      digit  = MOD( ival,10**pow ) / 10**(pow-1)
! Subprogram not used      s(i:i) = CHAR( ICHAR('0') + digit )
! Subprogram not used      pow    = pow - 1
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used   INCSTR = 0
! Subprogram not used 
! Subprogram not used end function INCSTR

! Subprogram not used integer function LASTND( cs )
! Subprogram not used   !-----------------------------------------------------------------------
! Subprogram not used   ! 	... Position of last non-digit in the first input token.
! Subprogram not used   ! 	    Return values:
! Subprogram not used   !     	    > 0  => position of last non-digit
! Subprogram not used   !     	    = 0  => token is all digits (or empty)
! Subprogram not used   !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used   implicit none
! Subprogram not used 
! Subprogram not used   !-----------------------------------------------------------------------
! Subprogram not used   ! 	... Dummy arguments
! Subprogram not used   !-----------------------------------------------------------------------
! Subprogram not used   character(len=*), intent(in) :: cs       !  Input character string
! Subprogram not used 
! Subprogram not used   !-----------------------------------------------------------------------
! Subprogram not used   ! 	... Local variables
! Subprogram not used   !-----------------------------------------------------------------------
! Subprogram not used   integer :: n, nn, digit
! Subprogram not used 
! Subprogram not used   n = GLC( cs )
! Subprogram not used   if( n == 0 ) then     ! empty string
! Subprogram not used      LASTND = 0
! Subprogram not used      return
! Subprogram not used   end if
! Subprogram not used 
! Subprogram not used   do nn = n,1,-1
! Subprogram not used      digit = ICHAR( cs(nn:nn) ) - ICHAR('0')
! Subprogram not used      if( digit < 0 .or. digit > 9 ) then
! Subprogram not used         LASTND = nn
! Subprogram not used         return
! Subprogram not used      end if
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used   LASTND = 0    ! all characters are digits
! Subprogram not used 
! Subprogram not used end function LASTND

integer function GLC( cs )
  !-----------------------------------------------------------------------
  ! 	... Position of last significant character in string. 
  !           Here significant means non-blank or non-null.
  !           Return values:
  !               > 0  => position of last significant character
  !               = 0  => no significant characters in string
  !-----------------------------------------------------------------------

  implicit none

  !-----------------------------------------------------------------------
  ! 	... Dummy arguments
  !-----------------------------------------------------------------------
  character(len=*), intent(in) :: cs       !  Input character string

  !-----------------------------------------------------------------------
  ! 	... Local variables
  !-----------------------------------------------------------------------
  integer :: l, n

  l = LEN( cs )
  if( l == 0 ) then
     GLC = 0
     return
  end if

  do n = l,1,-1
     if( cs(n:n) /= ' ' .and. cs(n:n) /= CHAR(0) ) then
        exit
     end if
  end do
  GLC = n

end function GLC

end module string_utils
