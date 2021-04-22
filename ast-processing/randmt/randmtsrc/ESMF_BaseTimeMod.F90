! $Id$
!
! Earth System Modeling Framework
! Copyright 2002-2003, University Corporation for Atmospheric Research,
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics
! Laboratory, University of Michigan, National Centers for Environmental
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
! NASA Goddard Space Flight Center.
! Licensed under the GPL.
!
!==============================================================================
!
!     ESMF BaseTime Module
      module ESMF_BaseTimeMod
!
!==============================================================================
!
! This file contains the BaseTime class definition and all BaseTime class
! methods.
!
!------------------------------------------------------------------------------
! INCLUDES











! Note that MAX_ALARMS must match MAX_WRF_ALARMS defined in 
! ../../frame/module_domain.F !!!  Eliminate this dependence with 
! grow-as-you-go AlarmList in ESMF_Clock...  

!
!===============================================================================
!BOPI
! !MODULE: ESMF_BaseTimeMod - Base ESMF time definition 
!
! !DESCRIPTION:
! Part of Time Manager F90 API wrapper of C++ implemenation
!
! This module serves only as the common Time definition inherited
! by {\tt ESMF\_TimeInterval} and {\tt ESMF\_Time}
!
! See {\tt ../include/ESMC\_BaseTime.h} for complete description
!
!------------------------------------------------------------------------------
! !USES:
      use ESMF_BaseMod    ! ESMF Base class
      implicit none
!
!------------------------------------------------------------------------------
! !PRIVATE TYPES:
      private
!------------------------------------------------------------------------------
!     ! ESMF_BaseTime
!
!     ! Base class type to match C++ BaseTime class in size only;
!     !  all dereferencing within class is performed by C++ implementation

      type ESMF_BaseTime
        integer(ESMF_KIND_I8) :: S   ! whole seconds
        integer(ESMF_KIND_I8) :: Sn  ! fractional seconds, numerator
        integer(ESMF_KIND_I8) :: Sd  ! fractional seconds, denominator
      end type

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
      public ESMF_BaseTime
!------------------------------------------------------------------------------
!
! !PUBLIC MEMBER FUNCTIONS:
!
! overloaded operators
      public seccmp
      public normalize_basetime
      public operator(+)
      private ESMF_BaseTimeSum
      public operator(-)
      private ESMF_BaseTimeDifference
      public operator(/)
      private ESMF_BaseTimeQuotI
      private ESMF_BaseTimeQuotI8
      public operator(.EQ.)
      private ESMF_BaseTimeEQ
      public operator(.NE.)
      private ESMF_BaseTimeNE
      public operator(.LT.)
      private ESMF_BaseTimeLT
      public operator(.GT.)
      private ESMF_BaseTimeGT
      public operator(.LE.)
      private ESMF_BaseTimeLE
      public operator(.GE.)
      private ESMF_BaseTimeGE

!==============================================================================
!
! INTERFACE BLOCKS
!
!==============================================================================
      interface operator(+)
        module procedure ESMF_BaseTimeSum
      end interface
      interface operator(-)
        module procedure ESMF_BaseTimeDifference
      end interface
      interface operator(/)
        module procedure ESMF_BaseTimeQuotI,ESMF_BaseTimeQuotI8
      end interface
      interface operator(.EQ.)
        module procedure ESMF_BaseTimeEQ
      end interface
      interface operator(.NE.)
        module procedure ESMF_BaseTimeNE
      end interface
      interface operator(.LT.)
        module procedure ESMF_BaseTimeLT
      end interface
      interface operator(.GT.)
        module procedure ESMF_BaseTimeGT
      end interface
      interface operator(.LE.)
        module procedure ESMF_BaseTimeLE
      end interface
      interface operator(.GE.)
        module procedure ESMF_BaseTimeGE
      end interface


!==============================================================================

      contains

!==============================================================================

SUBROUTINE normalize_basetime( basetime )
  ! Factor so abs(Sn) < Sd and ensure that signs of S and Sn match.
  ! Also, enforce consistency.
  ! YR and MM fields are ignored.
  IMPLICIT NONE
  TYPE(ESMF_BaseTime), INTENT(INOUT) :: basetime

  !PRINT *,'DEBUG:  BEGIN normalize_basetime()'
  ! Consistency check...
  IF ( basetime%Sd < 0 ) THEN
    CALL wrf_error_fatal( &
      'normalize_basetime:  denominator of seconds cannot be negative' )
  ENDIF
  IF ( ( basetime%Sd == 0 ) .AND. ( basetime%Sn .NE. 0 ) ) THEN
    CALL wrf_error_fatal( &
      'normalize_basetime:  denominator of seconds cannot be zero when numerator is non-zero' )
  ENDIF
  ! factor so abs(Sn) < Sd
  IF ( basetime%Sd > 0 ) THEN
    IF ( ABS( basetime%Sn ) .GE. basetime%Sd ) THEN
      !PRINT *,'DEBUG:  normalize_basetime() A1:  S,Sn,Sd = ',basetime%S,basetime%Sn,basetime%Sd
      basetime%S = basetime%S + ( basetime%Sn / basetime%Sd )
      basetime%Sn = mod( basetime%Sn, basetime%Sd )
      !PRINT *,'DEBUG:  normalize_basetime() A2:  S,Sn,Sd = ',basetime%S,basetime%Sn,basetime%Sd
    ENDIF
    ! change sign of Sn if it does not match S
    IF ( ( basetime%S > 0 ) .AND. ( basetime%Sn < 0 ) ) THEN
      !PRINT *,'DEBUG:  normalize_basetime() B1:  S,Sn,Sd = ',basetime%S,basetime%Sn,basetime%Sd
      basetime%S = basetime%S - 1_ESMF_KIND_I8
      basetime%Sn = basetime%Sn + basetime%Sd
      !PRINT *,'DEBUG:  normalize_basetime() B2:  S,Sn,Sd = ',basetime%S,basetime%Sn,basetime%Sd
    ENDIF
    IF ( ( basetime%S < 0 ) .AND. ( basetime%Sn > 0 ) ) THEN
      !PRINT *,'DEBUG:  normalize_basetime() C1:  S,Sn,Sd = ',basetime%S,basetime%Sn,basetime%Sd
      basetime%S = basetime%S + 1_ESMF_KIND_I8
      basetime%Sn = basetime%Sn - basetime%Sd
      !PRINT *,'DEBUG:  normalize_basetime() C2:  S,Sn,Sd = ',basetime%S,basetime%Sn,basetime%Sd
    ENDIF
  ENDIF
  !PRINT *,'DEBUG:  END normalize_basetime()'
END SUBROUTINE normalize_basetime

!==============================================================================

! Add two basetimes
      FUNCTION ESMF_BaseTimeSum( basetime1, basetime2 )
        TYPE(ESMF_BaseTime) :: ESMF_BaseTimeSum
        TYPE(ESMF_BaseTime), INTENT(IN) :: basetime1
        TYPE(ESMF_BaseTime), INTENT(IN) :: basetime2
        ! locals
        INTEGER (ESMF_KIND_I8) :: Sn1, Sd1, Sn2, Sd2, lcd
!  PRINT *,'DEBUG:  BEGIN ESMF_BaseTimeSum()'
!  PRINT *,'DEBUG:  ESMF_BaseTimeSum():  basetime1%S = ',basetime1%S
!  PRINT *,'DEBUG:  ESMF_BaseTimeSum():  basetime1%Sn = ',basetime1%Sn
!  PRINT *,'DEBUG:  ESMF_BaseTimeSum():  basetime1%Sd = ',basetime1%Sd
!  PRINT *,'DEBUG:  ESMF_BaseTimeSum():  basetime2%S = ',basetime2%S
!  PRINT *,'DEBUG:  ESMF_BaseTimeSum():  basetime2%Sn = ',basetime2%Sn
!  PRINT *,'DEBUG:  ESMF_BaseTimeSum():  basetime2%Sd = ',basetime2%Sd
        ESMF_BaseTimeSum   = basetime1
        ESMF_BaseTimeSum%S = ESMF_BaseTimeSum%S + basetime2%S
        Sn1 = basetime1%Sn
        Sd1 = basetime1%Sd
        Sn2 = basetime2%Sn
        Sd2 = basetime2%Sd
!  PRINT *,'DEBUG:  ESMF_BaseTimeSum():  Sn1 = ',Sn1
!  PRINT *,'DEBUG:  ESMF_BaseTimeSum():  Sd1 = ',Sd1
!  PRINT *,'DEBUG:  ESMF_BaseTimeSum():  Sn2 = ',Sn2
!  PRINT *,'DEBUG:  ESMF_BaseTimeSum():  Sd2 = ',Sd2
        IF      ( ( Sd1 .EQ. 0 ) .AND. ( Sd2 .EQ. 0 ) ) THEN
!  PRINT *,'DEBUG:  ESMF_BaseTimeSum():  no fractions'
          ESMF_BaseTimeSum%Sn = 0
          ESMF_BaseTimeSum%Sd = 0
        ELSE IF ( ( Sd1 .NE. 0 ) .AND. ( Sd2 .EQ. 0 ) ) THEN
          ESMF_BaseTimeSum%Sn = Sn1
          ESMF_BaseTimeSum%Sd = Sd1
        ELSE IF ( ( Sd1 .EQ. 0 ) .AND. ( Sd2 .NE. 0 ) ) THEN
          ESMF_BaseTimeSum%Sn = Sn2
          ESMF_BaseTimeSum%Sd = Sd2
        ELSE IF ( ( Sd1 .NE. 0 ) .AND. ( Sd2 .NE. 0 ) ) THEN
          CALL compute_lcd( Sd1 , Sd2 , lcd )
          ESMF_BaseTimeSum%Sd = lcd
          ESMF_BaseTimeSum%Sn = (Sn1 * lcd / Sd1) + (Sn2 * lcd / Sd2)
        ENDIF
!  PRINT *,'DEBUG:  ESMF_BaseTimeSum():  ESMF_BaseTimeSum%S = ',ESMF_BaseTimeSum%S
!  PRINT *,'DEBUG:  ESMF_BaseTimeSum():  ESMF_BaseTimeSum%Sn = ',ESMF_BaseTimeSum%Sn
!  PRINT *,'DEBUG:  ESMF_BaseTimeSum():  ESMF_BaseTimeSum%Sd = ',ESMF_BaseTimeSum%Sd
        CALL normalize_basetime( ESMF_BaseTimeSum )
!  PRINT *,'DEBUG:  END ESMF_BaseTimeSum()'
      END FUNCTION ESMF_BaseTimeSum


! Subtract two basetimes
      FUNCTION ESMF_BaseTimeDifference( basetime1, basetime2 )
        TYPE(ESMF_BaseTime) :: ESMF_BaseTimeDifference
        TYPE(ESMF_BaseTime), INTENT(IN) :: basetime1
        TYPE(ESMF_BaseTime), INTENT(IN) :: basetime2
        ! locals
        TYPE(ESMF_BaseTime) :: neg2

        neg2%S  = -basetime2%S
        neg2%Sn = -basetime2%Sn
        neg2%Sd =  basetime2%Sd

        ESMF_BaseTimeDifference = basetime1 + neg2

      END FUNCTION ESMF_BaseTimeDifference


! Divide basetime by 8-byte integer
! Subprogram not used       FUNCTION ESMF_BaseTimeQuotI8( basetime, divisor )
! Subprogram not used         TYPE(ESMF_BaseTime) :: ESMF_BaseTimeQuotI8
! Subprogram not used         TYPE(ESMF_BaseTime), INTENT(IN) :: basetime
! Subprogram not used         INTEGER(ESMF_KIND_I8), INTENT(IN) :: divisor
! Subprogram not used         ! locals
! Subprogram not used         INTEGER(ESMF_KIND_I8) :: d, n, dinit
! Subprogram not used 
! Subprogram not used !PRINT *,'DEBUG ESMF_BaseTimeQuotI8() A:  S,Sn,Sd = ', &
! Subprogram not used !  basetime%S,basetime%Sn,basetime%Sd
! Subprogram not used !PRINT *,'DEBUG ESMF_BaseTimeQuotI8() A:  divisor = ', divisor
! Subprogram not used         IF ( divisor == 0_ESMF_KIND_I8 ) THEN
! Subprogram not used           CALL wrf_error_fatal( 'ESMF_BaseTimeQuotI8:  divide by zero' )
! Subprogram not used         ENDIF
! Subprogram not used 
! Subprogram not used !$$$ move to default constructor
! Subprogram not used         ESMF_BaseTimeQuotI8%S  = 0
! Subprogram not used         ESMF_BaseTimeQuotI8%Sn = 0
! Subprogram not used         ESMF_BaseTimeQuotI8%Sd = 0
! Subprogram not used 
! Subprogram not used         ! convert to a fraction and divide by multipling the denonminator by 
! Subprogram not used         ! the divisor
! Subprogram not used         IF ( basetime%Sd == 0 ) THEN
! Subprogram not used           dinit = 1_ESMF_KIND_I8
! Subprogram not used         ELSE
! Subprogram not used           dinit = basetime%Sd
! Subprogram not used         ENDIF
! Subprogram not used         n = basetime%S * dinit + basetime%Sn
! Subprogram not used         d = dinit * divisor
! Subprogram not used !PRINT *,'DEBUG ESMF_BaseTimeQuotI8() B:  n,d = ',n,d
! Subprogram not used         CALL simplify( n, d, ESMF_BaseTimeQuotI8%Sn, ESMF_BaseTimeQuotI8%Sd )
! Subprogram not used !PRINT *,'DEBUG ESMF_BaseTimeQuotI8() C:  S,Sn,Sd = ', &
! Subprogram not used !  ESMF_BaseTimeQuotI8%S,ESMF_BaseTimeQuotI8%Sn,ESMF_BaseTimeQuotI8%Sd
! Subprogram not used         CALL normalize_basetime( ESMF_BaseTimeQuotI8 )
! Subprogram not used !PRINT *,'DEBUG ESMF_BaseTimeQuotI8() D:  S,Sn,Sd = ', &
! Subprogram not used !  ESMF_BaseTimeQuotI8%S,ESMF_BaseTimeQuotI8%Sn,ESMF_BaseTimeQuotI8%Sd
! Subprogram not used       END FUNCTION ESMF_BaseTimeQuotI8

! Divide basetime by integer
! Subprogram not used       FUNCTION ESMF_BaseTimeQuotI( basetime, divisor )
! Subprogram not used         TYPE(ESMF_BaseTime) :: ESMF_BaseTimeQuotI
! Subprogram not used         TYPE(ESMF_BaseTime), INTENT(IN) :: basetime
! Subprogram not used         INTEGER, INTENT(IN) :: divisor
! Subprogram not used         IF ( divisor == 0 ) THEN
! Subprogram not used           CALL wrf_error_fatal( 'ESMF_BaseTimeQuotI:  divide by zero' )
! Subprogram not used         ENDIF
! Subprogram not used         ESMF_BaseTimeQuotI = basetime / INT( divisor, ESMF_KIND_I8 )
! Subprogram not used       END FUNCTION ESMF_BaseTimeQuotI


! .EQ. for two basetimes
! Subprogram not used       FUNCTION ESMF_BaseTimeEQ( basetime1, basetime2 )
! Subprogram not used         LOGICAL :: ESMF_BaseTimeEQ
! Subprogram not used         TYPE(ESMF_BaseTime), INTENT(IN) :: basetime1
! Subprogram not used         TYPE(ESMF_BaseTime), INTENT(IN) :: basetime2
! Subprogram not used         INTEGER :: retval
! Subprogram not used         CALL seccmp( basetime1%S, basetime1%Sn, basetime1%Sd, &
! Subprogram not used                      basetime2%S, basetime2%Sn, basetime2%Sd, &
! Subprogram not used                      retval )
! Subprogram not used         ESMF_BaseTimeEQ = ( retval .EQ. 0 )
! Subprogram not used       END FUNCTION ESMF_BaseTimeEQ


! .NE. for two basetimes
! Subprogram not used       FUNCTION ESMF_BaseTimeNE( basetime1, basetime2 )
! Subprogram not used         LOGICAL :: ESMF_BaseTimeNE
! Subprogram not used         TYPE(ESMF_BaseTime), INTENT(IN) :: basetime1
! Subprogram not used         TYPE(ESMF_BaseTime), INTENT(IN) :: basetime2
! Subprogram not used         INTEGER :: retval
! Subprogram not used         CALL seccmp( basetime1%S, basetime1%Sn, basetime1%Sd, &
! Subprogram not used                      basetime2%S, basetime2%Sn, basetime2%Sd, &
! Subprogram not used                      retval )
! Subprogram not used         ESMF_BaseTimeNE = ( retval .NE. 0 )
! Subprogram not used       END FUNCTION ESMF_BaseTimeNE


! .LT. for two basetimes
      FUNCTION ESMF_BaseTimeLT( basetime1, basetime2 )
        LOGICAL :: ESMF_BaseTimeLT
        TYPE(ESMF_BaseTime), INTENT(IN) :: basetime1
        TYPE(ESMF_BaseTime), INTENT(IN) :: basetime2
        INTEGER :: retval
        CALL seccmp( basetime1%S, basetime1%Sn, basetime1%Sd, &
                     basetime2%S, basetime2%Sn, basetime2%Sd, &
                     retval )
        ESMF_BaseTimeLT = ( retval .LT. 0 )
      END FUNCTION ESMF_BaseTimeLT


! .GT. for two basetimes
! Subprogram not used       FUNCTION ESMF_BaseTimeGT( basetime1, basetime2 )
! Subprogram not used         LOGICAL :: ESMF_BaseTimeGT
! Subprogram not used         TYPE(ESMF_BaseTime), INTENT(IN) :: basetime1
! Subprogram not used         TYPE(ESMF_BaseTime), INTENT(IN) :: basetime2
! Subprogram not used         INTEGER :: retval
! Subprogram not used         CALL seccmp( basetime1%S, basetime1%Sn, basetime1%Sd, &
! Subprogram not used                      basetime2%S, basetime2%Sn, basetime2%Sd, &
! Subprogram not used                      retval )
! Subprogram not used         ESMF_BaseTimeGT = ( retval .GT. 0 )
! Subprogram not used       END FUNCTION ESMF_BaseTimeGT


! .LE. for two basetimes
! Subprogram not used       FUNCTION ESMF_BaseTimeLE( basetime1, basetime2 )
! Subprogram not used         LOGICAL :: ESMF_BaseTimeLE
! Subprogram not used         TYPE(ESMF_BaseTime), INTENT(IN) :: basetime1
! Subprogram not used         TYPE(ESMF_BaseTime), INTENT(IN) :: basetime2
! Subprogram not used         INTEGER :: retval
! Subprogram not used         CALL seccmp( basetime1%S, basetime1%Sn, basetime1%Sd, &
! Subprogram not used                      basetime2%S, basetime2%Sn, basetime2%Sd, &
! Subprogram not used                      retval )
! Subprogram not used         ESMF_BaseTimeLE = ( retval .LE. 0 )
! Subprogram not used       END FUNCTION ESMF_BaseTimeLE


! .GE. for two basetimes
      FUNCTION ESMF_BaseTimeGE( basetime1, basetime2 )
        LOGICAL :: ESMF_BaseTimeGE
        TYPE(ESMF_BaseTime), INTENT(IN) :: basetime1
        TYPE(ESMF_BaseTime), INTENT(IN) :: basetime2
        INTEGER :: retval
        CALL seccmp( basetime1%S, basetime1%Sn, basetime1%Sd, &
                     basetime2%S, basetime2%Sn, basetime2%Sd, &
                     retval )
        ESMF_BaseTimeGE = ( retval .GE. 0 )
      END FUNCTION ESMF_BaseTimeGE

!==============================================================================

! Subprogram not used SUBROUTINE compute_lcd( e1, e2, lcd )
! Subprogram not used       IMPLICIT NONE
! Subprogram not used       INTEGER(ESMF_KIND_I8), INTENT(IN) :: e1, e2
! Subprogram not used       INTEGER(ESMF_KIND_I8), INTENT(OUT) :: lcd
! Subprogram not used       INTEGER, PARAMETER ::  nprimes = 9
! Subprogram not used       INTEGER(ESMF_KIND_I8), DIMENSION(nprimes), PARAMETER :: primes = (/2,3,5,7,11,13,17,19,23/)
! Subprogram not used       INTEGER i
! Subprogram not used       INTEGER(ESMF_KIND_I8) d1, d2, p
! Subprogram not used 
! Subprogram not used       d1 = e1 ; d2 = e2
! Subprogram not used       IF ( d1 .EQ. 0 .AND. d2 .EQ. 0 ) THEN ; lcd = 1 ; RETURN ; ENDIF
! Subprogram not used       IF ( d1 .EQ. 0 ) d1 = d2 
! Subprogram not used       IF ( d2 .EQ. 0 ) d2 = d1 
! Subprogram not used       IF ( d1 .EQ. d2 ) THEN ; lcd = d1 ; RETURN ; ENDIF
! Subprogram not used       lcd = d1 * d2
! Subprogram not used       DO i = 1, nprimes
! Subprogram not used         p = primes(i)
! Subprogram not used         DO WHILE (lcd/p .NE. 0 .AND. &
! Subprogram not used           mod(lcd/p,d1) .EQ. 0 .AND. mod(lcd/p,d2) .EQ. 0) 
! Subprogram not used           lcd = lcd / p 
! Subprogram not used         END DO
! Subprogram not used       ENDDO
! Subprogram not used END SUBROUTINE compute_lcd

!==============================================================================

! Subprogram not used SUBROUTINE simplify( ni, di, no, do ) 
! Subprogram not used     IMPLICIT NONE
! Subprogram not used     INTEGER(ESMF_KIND_I8), INTENT(IN)  :: ni, di
! Subprogram not used     INTEGER(ESMF_KIND_I8), INTENT(OUT) :: no, do
! Subprogram not used     INTEGER, PARAMETER ::  nprimes = 9
! Subprogram not used     INTEGER(ESMF_KIND_I8), DIMENSION(nprimes), PARAMETER :: primes = (/2,3,5,7,11,13,17,19,23/)
! Subprogram not used     INTEGER(ESMF_KIND_I8) :: pr, d, n
! Subprogram not used     INTEGER :: np
! Subprogram not used     LOGICAL keepgoing
! Subprogram not used     IF ( ni .EQ. 0 ) THEN
! Subprogram not used       do = 1
! Subprogram not used       no = 0
! Subprogram not used       RETURN
! Subprogram not used     ENDIF
! Subprogram not used     IF ( mod( di , ni ) .EQ. 0 ) THEN
! Subprogram not used       do = di / ni
! Subprogram not used       no = 1
! Subprogram not used       RETURN
! Subprogram not used     ENDIF
! Subprogram not used     d = di
! Subprogram not used     n = ni
! Subprogram not used     DO np = 1, nprimes
! Subprogram not used       pr = primes(np)
! Subprogram not used       keepgoing = .TRUE.
! Subprogram not used       DO WHILE ( keepgoing )
! Subprogram not used         keepgoing = .FALSE.
! Subprogram not used         IF ( d/pr .NE. 0 .AND. n/pr .NE. 0 .AND. MOD(d,pr) .EQ. 0 .AND. MOD(n,pr) .EQ. 0 ) THEN
! Subprogram not used           d = d / pr
! Subprogram not used           n = n / pr
! Subprogram not used           keepgoing = .TRUE.
! Subprogram not used         ENDIF
! Subprogram not used       ENDDO
! Subprogram not used     ENDDO
! Subprogram not used     do = d
! Subprogram not used     no = n
! Subprogram not used     RETURN
! Subprogram not used END SUBROUTINE simplify

!==============================================================================

! spaceship operator for seconds + Sn/Sd
SUBROUTINE seccmp(S1, Sn1, Sd1, S2, Sn2, Sd2, retval )
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: retval
!
! !ARGUMENTS:
  INTEGER(ESMF_KIND_I8), INTENT(IN) :: S1, Sn1, Sd1
  INTEGER(ESMF_KIND_I8), INTENT(IN) :: S2, Sn2, Sd2
! local
  INTEGER(ESMF_KIND_I8) :: lcd, n1, n2

  n1 = Sn1
  n2 = Sn2
  if ( ( n1 .ne. 0 ) .or. ( n2 .ne. 0 ) ) then
    CALL compute_lcd( Sd1, Sd2, lcd )
    if ( Sd1 .ne. 0 ) n1 = n1 * ( lcd / Sd1 )
    if ( Sd2 .ne. 0 ) n2 = n2 * ( lcd / Sd2 )
  endif

  if ( S1 .GT. S2 ) retval = 1
  if ( S1 .LT. S2 ) retval = -1
  IF ( S1 .EQ. S2 ) THEN
    IF (n1 .GT. n2) retval = 1
    IF (n1 .LT. n2) retval = -1
    IF (n1 .EQ. n2) retval = 0
  ENDIF
END SUBROUTINE seccmp

!==============================================================================

      end module ESMF_BaseTimeMod
