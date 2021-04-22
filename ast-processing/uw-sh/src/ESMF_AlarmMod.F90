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
!     ESMF Alarm Module
      module ESMF_AlarmMod
!
!==============================================================================
!
! This file contains the Alarm class definition and all Alarm class 
! methods.
!
!------------------------------------------------------------------------------
! INCLUDES


!===============================================================================
!BOPI
!
! !MODULE: ESMF_AlarmMod
!
! !DESCRIPTION:
! Part of Time Manager F90 API wrapper of C++ implemenation
!
! Defines F90 wrapper entry points for corresponding
! C++ class {\tt ESMC\_Alarm}
!
! See {\tt ../include/ESMC\_Alarm.h} for complete description
!
!------------------------------------------------------------------------------
! !USES:
      ! inherit from ESMF base class
      use ESMF_BaseMod

      ! associated derived types
      use ESMF_TimeIntervalMod
      use ESMF_TimeMod

      implicit none

!------------------------------------------------------------------------------
! !PRIVATE TYPES:
     private
!------------------------------------------------------------------------------
!     ! ESMF_Alarm
!
!     ! F90 class type to match C++ Alarm class in size only;
!     !  all dereferencing within class is performed by C++ implementation

! internals for ESMF_Alarm
      type ESMF_AlarmInt
        character(len=256) :: name = " "
        type(ESMF_TimeInterval) :: RingInterval
        type(ESMF_Time)  :: RingTime
        type(ESMF_Time)  :: PrevRingTime
        type(ESMF_Time)  :: StopTime
        integer :: ID
        integer :: AlarmMutex
        logical :: Ringing
        logical :: Enabled
        logical :: RingTimeSet
        logical :: RingIntervalSet
        logical :: StopTimeSet
      end type

! Actual public type:  this bit allows easy mimic of "deep" ESMF_AlarmCreate
! in ESMF 2.1.0+.  Note that ESMF_AlarmCreate is in a separate module to avoid 
! cyclic dependence.  
! NOTE:  DO NOT ADD NON-POINTER STATE TO THIS DATA TYPE.  It emulates ESMF 
!        shallow-copy-masquerading-as-reference-copy insanity.  
      type ESMF_Alarm
        type(ESMF_AlarmInt), pointer :: alarmint => null()
      end type

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
      public ESMF_Alarm
      public ESMF_AlarmInt   ! needed on AIX but not PGI
!------------------------------------------------------------------------------

! !PUBLIC MEMBER FUNCTIONS:
      public ESMF_AlarmDestroy
      public ESMF_AlarmSet
      public ESMF_AlarmGet
!      public ESMF_AlarmGetRingInterval
!      public ESMF_AlarmSetRingInterval
!      public ESMF_AlarmGetRingTime
!      public ESMF_AlarmSetRingTime
!      public ESMF_AlarmGetPrevRingTime
!      public ESMF_AlarmSetPrevRingTime
!      public ESMF_AlarmGetStopTime
!      public ESMF_AlarmSetStopTime
      public ESMF_AlarmEnable
      public ESMF_AlarmDisable
      public ESMF_AlarmRingerOn
      public ESMF_AlarmRingerOff
      public ESMF_AlarmIsRinging
!      public ESMF_AlarmCheckRingTime
      public operator(==)
 
! Required inherited and overridden ESMF_Base class methods

!      public ESMF_AlarmRead
!      public ESMF_AlarmWrite
      public ESMF_AlarmValidate
      public ESMF_AlarmPrint

! !PRIVATE MEMBER FUNCTIONS:
      private ESMF_AlarmEQ
!EOPI

!------------------------------------------------------------------------------
! The following line turns the CVS identifier string into a printable variable.
      character(*), parameter, private :: version = &
      '$Id$'

!==============================================================================
!
! INTERFACE BLOCKS
!
!==============================================================================
!BOP
! !INTERFACE:
      interface operator(==)

! !PRIVATE MEMBER FUNCTIONS:
      module procedure ESMF_AlarmEQ

! !DESCRIPTION:
!     This interface overloads the == operator for the {\tt ESMF\_Alarm} class
!
!EOP
      end interface
!
!------------------------------------------------------------------------------

!==============================================================================

      contains

!==============================================================================

!------------------------------------------------------------------------------
!
! This section includes the Set methods.
!
!------------------------------------------------------------------------------
!BOP
! !IROUTINE: ESMF_AlarmSet - Initializes an alarm

! !INTERFACE:
      subroutine ESMF_AlarmSet(alarm, name, RingTime, RingInterval, &
                               StopTime, Enabled, rc)

! !ARGUMENTS:
      type(ESMF_Alarm), intent(inout) :: alarm
      character(len=*), intent(in), optional :: name
      type(ESMF_Time), intent(in), optional :: RingTime
      type(ESMF_TimeInterval), intent(in), optional :: RingInterval
      type(ESMF_Time), intent(in), optional :: StopTime
      logical, intent(in), optional :: Enabled
      integer, intent(out), optional :: rc

! !DESCRIPTION:
!     Initializes an {\tt ESMF\_Alarm}
!
!     The arguments are:
!     \begin{description}
!     \item[alarm]
!          The object instance to initialize
!     \item[{[RingTime]}]
!          Optional ring time for one-shot or first repeating alarm
!     \item[{[RingInterval]}]
!          Optional ring interval for repeating alarms
!     \item[{[StopTime]}]
!          Optional stop time for repeating alarms
!     \item[Enabled]
!          Alarm enabled/disabled
!     \item[{[rc]}]
!          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
!     \end{description}
!
! !REQUIREMENTS:
!     TMG4.1, TMG4.7
!EOP
      IF ( ASSOCIATED( alarm%alarmint ) ) THEN
        alarm%alarmint%RingTimeSet = .FALSE.
        alarm%alarmint%RingIntervalSet = .FALSE.
        alarm%alarmint%StopTimeSet = .FALSE.
        IF ( PRESENT( name ) ) THEN
          alarm%alarmint%name = name
        END IF
        IF ( PRESENT( RingInterval ) ) THEN
          alarm%alarmint%RingInterval = RingInterval
          alarm%alarmint%RingIntervalSet = .TRUE.
        ENDIF
        IF ( PRESENT( RingTime ) ) THEN
          alarm%alarmint%RingTime = RingTime
          alarm%alarmint%RingTimeSet = .TRUE.
        ENDIF
        IF ( PRESENT( StopTime ) ) THEN
          alarm%alarmint%StopTime = StopTime
          alarm%alarmint%StopTimeSet = .TRUE.
        ENDIF
        alarm%alarmint%Enabled = .TRUE.
        IF ( PRESENT( Enabled ) ) THEN
          alarm%alarmint%Enabled = Enabled
        ENDIF
        IF ( PRESENT( rc ) ) THEN
          rc = ESMF_SUCCESS
        ENDIF
        alarm%alarmint%Ringing = .FALSE.
        alarm%alarmint%Enabled = .TRUE.
      ELSE
        IF ( PRESENT( rc ) ) rc = ESMF_FAILURE
      ENDIF

      end subroutine ESMF_AlarmSet



! Deallocate memory for ESMF_Alarm
! Subprogram not used       SUBROUTINE ESMF_AlarmDestroy( alarm, rc )
! Subprogram not used          TYPE(ESMF_Alarm), INTENT(INOUT) :: alarm
! Subprogram not used          INTEGER,          INTENT(  OUT), OPTIONAL :: rc
! Subprogram not used          IF ( ASSOCIATED( alarm%alarmint ) ) THEN
! Subprogram not used            DEALLOCATE( alarm%alarmint )
! Subprogram not used          ENDIF
! Subprogram not used          ! TBH:  ignore deallocate errors, for now
! Subprogram not used          IF ( PRESENT( rc ) ) rc = ESMF_SUCCESS
! Subprogram not used       END SUBROUTINE ESMF_AlarmDestroy



!------------------------------------------------------------------------------
!BOP
! !IROUTINE: ESMF_AlarmGetRingInterval - Get an alarm's ring interval
!
! !INTERFACE:
      subroutine ESMF_AlarmGetRingInterval(alarm, RingInterval, rc)

! !ARGUMENTS:
      type(ESMF_Alarm), intent(in) :: alarm
      type(ESMF_TimeInterval), intent(out) :: RingInterval
      integer, intent(out), optional :: rc

! !DESCRIPTION:
!     Get an {\tt ESMF\_Alarm}'s ring interval
!
!     The arguments are:
!     \begin{description}
!     \item[alarm]
!          The object instance to get the ring interval
!     \item[RingInterval]
!          The {\tt Alarm}'s ring interval
!     \item[{[rc]}]
!          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
!     \end{description}

! !REQUIREMENTS:
!     TMG4.7
!EOP
      IF ( ASSOCIATED( alarm%alarmint ) ) THEN
        IF ( alarm%alarmint%RingIntervalSet )THEN
           RingInterval= alarm%alarmint%RingInterval
           IF ( PRESENT( rc ) ) rc = ESMF_SUCCESS
        ELSE
           IF ( PRESENT( rc ) ) rc = ESMF_FAILURE
        END IF
      ELSE
        IF ( PRESENT( rc ) ) rc = ESMF_FAILURE
      ENDIF
      end subroutine ESMF_AlarmGetRingInterval
 
!------------------------------------------------------------------------------
!BOP
! !IROUTINE: ESMF_AlarmSetRingInterval - Set an alarm's ring interval
!
! !INTERFACE:
! Subprogram not used       subroutine ESMF_AlarmSetRingInterval(alarm, RingInterval, rc)
! Subprogram not used 
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used       type(ESMF_Alarm), intent(out) :: alarm
! Subprogram not used       type(ESMF_TimeInterval), intent(in) :: RingInterval
! Subprogram not used       integer, intent(out), optional :: rc
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !     Set an {\tt ESMF\_Alarm}'s ring interval
! Subprogram not used !
! Subprogram not used !     The arguments are:
! Subprogram not used !     \begin{description}
! Subprogram not used !     \item[alarm]
! Subprogram not used !          The object instance to set the ring interval
! Subprogram not used !     \item[RingInterval]
! Subprogram not used !          The {\tt Alarm}'s ring interval
! Subprogram not used !     \item[{[rc]}]
! Subprogram not used !          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
! Subprogram not used !     \end{description}
! Subprogram not used !
! Subprogram not used ! !REQUIREMENTS:
! Subprogram not used !     TMG4.5.2, TMG4.7
! Subprogram not used !EOP
! Subprogram not used       CALL wrf_error_fatal( 'ESMF_AlarmSetRingInterval not supported' )
! Subprogram not used       end subroutine ESMF_AlarmSetRingInterval

!------------------------------------------------------------------------------
!BOP
! !IROUTINE:  ESMF_AlarmGetRingTime - Get an alarm's time to ring
!
! !INTERFACE:
      subroutine ESMF_AlarmGetRingTime(alarm, RingTime, rc)

! !ARGUMENTS:
      type(ESMF_Alarm), intent(in) :: alarm
      type(ESMF_Time), intent(out) :: RingTime
      integer, intent(out), optional :: rc

! !DESCRIPTION:
!     Get an {\tt ESMF\_Alarm}'s time to ring
!
!     The arguments are:
!     \begin{description}
!     \item[alarm]
!          The object instance to get the ring time
!     \item[RingTime]
!          The {\tt ESMF\_Alarm}'s ring time
!     \item[{[rc]}]
!          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
!     \end{description}
!
! !REQUIREMENTS:
!     TMG4.7, TMG4.8
!EOP
      type(ESMF_Time) :: PrevRingTime
      type(ESMF_TimeInterval) :: RingInterval
      integer :: ierr

      IF ( ASSOCIATED( alarm%alarmint ) ) THEN
        IF ( alarm%alarmint%RingIntervalSet )THEN
           PrevRingTime = alarm%alarmint%PrevRingTime
           call ESMF_AlarmGetRingInterval( alarm, RingInterval, ierr)
           IF ( PRESENT( rc ) .AND. (ierr /= ESMF_SUCCESS) )THEN
              rc = ierr
              return
           END IF
           RingTime = PrevRingTime + RingInterval
           IF ( PRESENT( rc ) ) rc = ESMF_SUCCESS
        ELSE IF ( alarm%alarmint%RingTimeSet )THEN
           RingTime = alarm%alarmint%RingTime
           IF ( PRESENT( rc ) ) rc = ESMF_SUCCESS
        ELSE
           IF ( PRESENT( rc ) ) rc = ESMF_FAILURE
        END IF
      ELSE
        IF ( PRESENT( rc ) ) rc = ESMF_FAILURE
      ENDIF
      end subroutine ESMF_AlarmGetRingTime

!------------------------------------------------------------------------------
!BOP
! !IROUTINE:  ESMF_AlarmSetRingTime - Set an alarm's time to ring
!
! !INTERFACE:
! Subprogram not used       subroutine ESMF_AlarmSetRingTime(alarm, RingTime, rc)
! Subprogram not used 
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used       type(ESMF_Alarm), intent(out) :: alarm
! Subprogram not used       type(ESMF_Time), intent(in) :: RingTime
! Subprogram not used       integer, intent(out), optional :: rc
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !     Set an {\tt ESMF\_Alarm}'s time to ring
! Subprogram not used !
! Subprogram not used !     The arguments are:
! Subprogram not used !     \begin{description}
! Subprogram not used !     \item[alarm]
! Subprogram not used !          The object instance to set the ring time
! Subprogram not used !     \item[RingTime]
! Subprogram not used !          The {\tt ESMF\_Alarm}'s ring time to set
! Subprogram not used !     \item[{[rc]}]
! Subprogram not used !          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
! Subprogram not used !     \end{description}
! Subprogram not used !
! Subprogram not used ! !REQUIREMENTS:
! Subprogram not used !     TMG4.5.1, TMG4.7, TMG4.8
! Subprogram not used !EOP
! Subprogram not used       CALL wrf_error_fatal( 'ESMF_AlarmSetRingTime not supported' )
! Subprogram not used       end subroutine ESMF_AlarmSetRingTime

!------------------------------------------------------------------------------
!BOP
! !IROUTINE:  ESMF_AlarmGet - Get an alarm's parameters -- compatibility with ESMF 2.0.1
!
! !INTERFACE:
      subroutine ESMF_AlarmGet(alarm, name, RingTime, PrevRingTime, RingInterval, rc)

! !ARGUMENTS:
      type(ESMF_Alarm), intent(in) :: alarm
      character(len=*), intent(out), optional :: name
      type(ESMF_Time), intent(out), optional :: RingTime
      type(ESMF_Time), intent(out), optional :: PrevRingTime
      type(ESMF_TimeInterval), intent(out), optional :: RingInterval
      integer, intent(out), optional :: rc
      integer :: ierr

! !DESCRIPTION:
!     Get an {\tt ESMF\_Alarm}'s previous ring time
!
!     The arguments are:
!     \begin{description}
!     \item[alarm]
!          The object instance to get
!     \item[ringTime]
!          The ring time for a one-shot alarm or the next repeating alarm.
!     \item[ringInterval]
!          The ring interval for repeating (interval) alarms.
!     \item[PrevRingTime]
!          The {\tt ESMF\_Alarm}'s previous ring time
!     \item[{[rc]}]
!          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
!     \end{description}
!
! !REQUIREMENTS:
!     TMG4.7, TMG4.8
!EOP

      ierr = ESMF_SUCCESS

      IF ( PRESENT(name) ) THEN
         IF ( ASSOCIATED( alarm%alarmint ) ) THEN
           name = alarm%alarmint%name
         ELSE
           ierr = ESMF_FAILURE
         END IF
      ENDIF
      IF ( PRESENT(PrevRingTime) ) THEN
        CALL ESMF_AlarmGetPrevRingTime(alarm, PrevRingTime, rc=ierr)
      ENDIF
      IF ( PRESENT(RingTime) ) THEN
        CALL ESMF_AlarmGetRingTime(alarm, RingTime, rc=ierr)
      ENDIF
      IF ( PRESENT(RingInterval) ) THEN
        CALL ESMF_AlarmGetRingInterval(alarm, RingInterval, rc=ierr)
      ENDIF

      IF ( PRESENT(rc) ) THEN
        rc = ierr
      ENDIF

      end subroutine ESMF_AlarmGet

!------------------------------------------------------------------------------
!BOP
! !IROUTINE:  ESMF_AlarmGetPrevRingTime - Get an alarm's previous ring time
!
! !INTERFACE:
      subroutine ESMF_AlarmGetPrevRingTime(alarm, PrevRingTime, rc)

! !ARGUMENTS:
      type(ESMF_Alarm), intent(in) :: alarm
      type(ESMF_Time), intent(out) :: PrevRingTime
      integer, intent(out), optional :: rc

! !DESCRIPTION:
!     Get an {\tt ESMF\_Alarm}'s previous ring time
!
!     The arguments are:
!     \begin{description}
!     \item[alarm]
!          The object instance to get the previous ring time
!     \item[PrevRingTime]
!          The {\tt ESMF\_Alarm}'s previous ring time
!     \item[{[rc]}]
!          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
!     \end{description}
!
! !REQUIREMENTS:
!     TMG4.7, TMG4.8
!EOP
      IF ( ASSOCIATED( alarm%alarmint ) ) THEN
        PrevRingTime = alarm%alarmint%PrevRingTime
        IF ( PRESENT( rc ) ) rc = ESMF_SUCCESS
      ELSE
        IF ( PRESENT( rc ) ) rc = ESMF_FAILURE
      ENDIF
      end subroutine ESMF_AlarmGetPrevRingTime

!------------------------------------------------------------------------------
!BOP
! !IROUTINE:  ESMF_AlarmSetPrevRingTime - Set an alarm's previous ring time
!
! !INTERFACE:
! Subprogram not used       subroutine ESMF_AlarmSetPrevRingTime(alarm, PrevRingTime, rc)
! Subprogram not used 
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used       type(ESMF_Alarm), intent(out) :: alarm
! Subprogram not used       type(ESMF_Time), intent(in) :: PrevRingTime
! Subprogram not used       integer, intent(out), optional :: rc
! Subprogram not used    
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !     Set an {\tt ESMF\_Alarm}'s previous ring time
! Subprogram not used !
! Subprogram not used !     The arguments are:
! Subprogram not used !     \begin{description}
! Subprogram not used !     \item[alarm]
! Subprogram not used !          The object instance to set the previous ring time
! Subprogram not used !     \item[PrevRingTime]
! Subprogram not used !          The {\tt ESMF\_Alarm}'s previous ring time to set
! Subprogram not used !     \item[{[rc]}]
! Subprogram not used !          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
! Subprogram not used !     \end{description}
! Subprogram not used !
! Subprogram not used ! !REQUIREMENTS:
! Subprogram not used !     TMG4.7, TMG4.8
! Subprogram not used !EOP
! Subprogram not used       CALL wrf_error_fatal( 'ESMF_AlarmSetPrevRingTime not supported' )
! Subprogram not used       end subroutine ESMF_AlarmSetPrevRingTime

!------------------------------------------------------------------------------
!BOP
! !IROUTINE:  ESMF_AlarmGetStopTime - Get an alarm's stop time
!
! !INTERFACE:
! Subprogram not used       subroutine ESMF_AlarmGetStopTime(alarm, StopTime, rc)
! Subprogram not used 
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used       type(ESMF_Alarm), intent(in) :: alarm
! Subprogram not used       type(ESMF_Time), intent(out) :: StopTime
! Subprogram not used       integer, intent(out), optional :: rc
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !     Get an {\tt ESMF\_Alarm}'s stop time
! Subprogram not used !
! Subprogram not used !     The arguments are:
! Subprogram not used !     \begin{description}
! Subprogram not used !     \item[alarm]
! Subprogram not used !          The object instance to get the stop time
! Subprogram not used !     \item[StopTime]
! Subprogram not used !          The {\tt ESMF\_Alarm}'s stop time
! Subprogram not used !     \item[{[rc]}]
! Subprogram not used !          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
! Subprogram not used !     \end{description}
! Subprogram not used !
! Subprogram not used ! !REQUIREMENTS:
! Subprogram not used !     TMG4.5.2, TMG4.7
! Subprogram not used !EOP
! Subprogram not used       CALL wrf_error_fatal( 'ESMF_AlarmGetStopTime not supported' )
! Subprogram not used       end subroutine ESMF_AlarmGetStopTime

!------------------------------------------------------------------------------
!BOP
! !IROUTINE:  ESMF_AlarmSetStopTime - Set an alarm's stop time
!
! !INTERFACE:
! Subprogram not used       subroutine ESMF_AlarmSetStopTime(alarm, StopTime, rc)
! Subprogram not used 
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used       type(ESMF_Alarm), intent(out) :: alarm
! Subprogram not used       type(ESMF_Time), intent(in) :: StopTime
! Subprogram not used       integer, intent(out), optional :: rc
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !     Set an {\tt ESMF\_Alarm}'s stop time
! Subprogram not used !
! Subprogram not used !     The arguments are:
! Subprogram not used !     \begin{description}
! Subprogram not used !     \item[alarm]
! Subprogram not used !          The object instance to set the stop time
! Subprogram not used !     \item[StopTime]
! Subprogram not used !          The {\tt ESMF\_Alarm}'s stop time
! Subprogram not used !     \item[{[rc]}]
! Subprogram not used !          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
! Subprogram not used !     \end{description}
! Subprogram not used !
! Subprogram not used ! !REQUIREMENTS:
! Subprogram not used !     TMG4.5.2, TMG4.7
! Subprogram not used !EOP
! Subprogram not used       CALL wrf_error_fatal( 'ESMF_AlarmSetStopTime not supported' )
! Subprogram not used       end subroutine ESMF_AlarmSetStopTime

!------------------------------------------------------------------------------
!BOP
! !IROUTINE: ESMF_AlarmEnable - Enables an alarm

! !INTERFACE:
! Subprogram not used       subroutine ESMF_AlarmEnable(alarm, rc)
! Subprogram not used 
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used       type(ESMF_Alarm), intent(inout) :: alarm
! Subprogram not used       integer, intent(out), optional :: rc
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !     Enables an {\tt ESMF\_Alarm} to function
! Subprogram not used !
! Subprogram not used !     The arguments are:
! Subprogram not used !     \begin{description}
! Subprogram not used !     \item[alarm]
! Subprogram not used !          The object instance to enable
! Subprogram not used !     \item[{[rc]}]
! Subprogram not used !          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
! Subprogram not used !     \end{description}
! Subprogram not used 
! Subprogram not used ! !REQUIREMENTS:
! Subprogram not used !     TMG4.5.3
! Subprogram not used !EOP
! Subprogram not used       IF ( ASSOCIATED( alarm%alarmint ) ) THEN
! Subprogram not used         alarm%alarmint%Enabled = .TRUE.
! Subprogram not used         IF ( PRESENT( rc ) ) rc = ESMF_SUCCESS
! Subprogram not used       ELSE
! Subprogram not used         IF ( PRESENT( rc ) ) rc = ESMF_FAILURE
! Subprogram not used       ENDIF
! Subprogram not used       end subroutine ESMF_AlarmEnable

!------------------------------------------------------------------------------
!BOP
! !IROUTINE: ESMF_AlarmDisable - Disables an alarm

! !INTERFACE:
! Subprogram not used       subroutine ESMF_AlarmDisable(alarm, rc)
! Subprogram not used 
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used       type(ESMF_Alarm), intent(inout) :: alarm
! Subprogram not used       integer, intent(out), optional :: rc
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !     Disables an {\tt ESMF\_Alarm}
! Subprogram not used !
! Subprogram not used !     The arguments are:
! Subprogram not used !     \begin{description}
! Subprogram not used !     \item[alarm]
! Subprogram not used !          The object instance to disable
! Subprogram not used !     \item[{[rc]}]
! Subprogram not used !          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
! Subprogram not used !     \end{description}
! Subprogram not used 
! Subprogram not used ! !REQUIREMENTS:
! Subprogram not used !     TMG4.5.3
! Subprogram not used !EOP
! Subprogram not used       IF ( ASSOCIATED( alarm%alarmint ) ) THEN
! Subprogram not used         alarm%alarmint%Enabled = .FALSE.
! Subprogram not used         IF ( PRESENT( rc ) ) rc = ESMF_SUCCESS
! Subprogram not used       ELSE
! Subprogram not used         IF ( PRESENT( rc ) ) rc = ESMF_FAILURE
! Subprogram not used       ENDIF
! Subprogram not used       end subroutine ESMF_AlarmDisable

!------------------------------------------------------------------------------
!BOP
! !IROUTINE:  ESMF_AlarmRingerOn - Turn on an alarm


! !INTERFACE:
! Subprogram not used       subroutine ESMF_AlarmRingerOn(alarm, rc)
! Subprogram not used 
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used       type(ESMF_Alarm), intent(inout) :: alarm
! Subprogram not used       integer, intent(out), optional :: rc
! Subprogram not used     
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !     Turn on an {\tt ESMF\_Alarm}; sets ringing state
! Subprogram not used !
! Subprogram not used !     The arguments are:
! Subprogram not used !     \begin{description}
! Subprogram not used !     \item[alarm]
! Subprogram not used !          The object instance to turn on
! Subprogram not used !     \item[{[rc]}]
! Subprogram not used !          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
! Subprogram not used !     \end{description}
! Subprogram not used !
! Subprogram not used ! !REQUIREMENTS:
! Subprogram not used !     TMG4.6
! Subprogram not used !EOP
! Subprogram not used       IF ( ASSOCIATED( alarm%alarmint ) ) THEN
! Subprogram not used         IF ( alarm%alarmint%Enabled ) THEN
! Subprogram not used           alarm%alarmint%Ringing = .TRUE.
! Subprogram not used           IF ( PRESENT( rc ) ) rc = ESMF_SUCCESS
! Subprogram not used         ELSE
! Subprogram not used           alarm%alarmint%Ringing = .FALSE.
! Subprogram not used           IF ( PRESENT( rc ) ) rc = ESMF_FAILURE
! Subprogram not used         ENDIF
! Subprogram not used       ELSE
! Subprogram not used         IF ( PRESENT( rc ) ) rc = ESMF_FAILURE
! Subprogram not used       ENDIF
! Subprogram not used 
! Subprogram not used       end subroutine ESMF_AlarmRingerOn

!------------------------------------------------------------------------------
!BOP
! !IROUTINE:  ESMF_AlarmRingerOff - Turn off an alarm

! !INTERFACE:
      subroutine ESMF_AlarmRingerOff(alarm, rc)

! !ARGUMENTS:
      type(ESMF_Alarm), intent(inout) :: alarm
      integer, intent(out), optional :: rc
    
! !DESCRIPTION:
!     Turn off an {\tt ESMF\_Alarm}; unsets ringing state
!
!     The arguments are:
!     \begin{description}
!     \item[alarm]
!          The object instance to turn off   
!     \item[{[rc]}]
!          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
!     \end{description}

! !REQUIREMENTS:
!     TMG4.6
!EOP
      IF ( ASSOCIATED( alarm%alarmint ) ) THEN
        alarm%alarmint%Ringing = .FALSE.
        IF ( alarm%alarmint%Enabled ) THEN
          IF ( PRESENT( rc ) ) rc = ESMF_SUCCESS
        ELSE
          IF ( PRESENT( rc ) ) rc = ESMF_FAILURE
        ENDIF
      ELSE
        IF ( PRESENT( rc ) ) rc = ESMF_FAILURE
      ENDIF
      end subroutine ESMF_AlarmRingerOff

!------------------------------------------------------------------------------
!BOP
! !IROUTINE:  ESMF_AlarmIsRinging - Check if alarm is ringing

! !INTERFACE:
      function ESMF_AlarmIsRinging(alarm, rc)
!
! !RETURN VALUE:
      logical :: ESMF_AlarmIsRinging

! !ARGUMENTS:
      type(ESMF_Alarm), intent(in) :: alarm
      integer, intent(out), optional :: rc

! !DESCRIPTION:
!     Check if {\tt ESMF\_Alarm} is ringing.
!
!     The arguments are:
!     \begin{description}
!     \item[alarm]
!          The object instance to check for ringing state  
!     \item[{[rc]}]
!          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
!     \end{description}

! !REQUIREMENTS:
!     TMG4.4
!EOP
      IF ( ASSOCIATED( alarm%alarmint ) ) THEN
        IF ( alarm%alarmint%Enabled ) THEN
          ESMF_AlarmIsRinging = alarm%alarmint%Ringing
          IF ( PRESENT( rc ) ) rc = ESMF_SUCCESS
        ELSE
          ESMF_AlarmIsRinging = .FALSE.
          IF ( PRESENT( rc ) ) rc = ESMF_FAILURE
        ENDIF
      ELSE
        IF ( PRESENT( rc ) ) rc = ESMF_FAILURE
      ENDIF
      end function ESMF_AlarmIsRinging

!------------------------------------------------------------------------------
!BOP
! !IROUTINE: ESMF_AlarmCheckRingTime - Method used by a clock to check whether to trigger an alarm
!
! !INTERFACE:
! Subprogram not used       function ESMF_AlarmCheckRingTime(alarm, ClockCurrTime, positive, rc)
! Subprogram not used !
! Subprogram not used ! !RETURN VALUE:
! Subprogram not used       logical :: ESMF_AlarmCheckRingTime
! Subprogram not used !
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used       type(ESMF_Alarm), intent(inout) :: alarm
! Subprogram not used       type(ESMF_Time), intent(in) :: ClockCurrTime
! Subprogram not used       integer, intent(in) :: positive
! Subprogram not used       integer, intent(out), optional :: rc
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !     Main method used by a {\tt ESMF\_Clock} to check whether to trigger
! Subprogram not used !     the {\tt ESMF\_Alarm} 
! Subprogram not used !
! Subprogram not used !     The arguments are:
! Subprogram not used !     \begin{description}
! Subprogram not used !     \item[alarm]
! Subprogram not used !          The object instance to check if time to ring   
! Subprogram not used !     \item[ClockCurrTime]
! Subprogram not used !          The {\tt ESMF\_Clock}'s current time
! Subprogram not used !     \item[positive]
! Subprogram not used !          Whether to check ring time in the positive or negative direction
! Subprogram not used !     \item[{[rc]}]
! Subprogram not used !          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
! Subprogram not used !     \end{description}
! Subprogram not used 
! Subprogram not used ! !REQUIREMENTS:
! Subprogram not used !     TMG4.4, TMG4.6
! Subprogram not used !EOP
! Subprogram not used       CALL wrf_error_fatal( 'ESMF_AlarmCheckRingTime not supported' )
! Subprogram not used       ESMF_AlarmCheckRingTime = .FALSE.  ! keep compilers happy
! Subprogram not used       end function ESMF_AlarmCheckRingTime

!------------------------------------------------------------------------------
!BOP
! !IROUTINE:  ESMF_AlarmEQ - Compare two alarms for equality
!
! !INTERFACE:
! Subprogram not used       function ESMF_AlarmEQ(alarm1, alarm2)
! Subprogram not used !
! Subprogram not used ! !RETURN VALUE:
! Subprogram not used       logical :: ESMF_AlarmEQ
! Subprogram not used 
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used       type(ESMF_Alarm), intent(in) :: alarm1
! Subprogram not used       type(ESMF_Alarm), intent(in) :: alarm2
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !     Compare two alarms for equality; return true if equal, false otherwise
! Subprogram not used !     Maps to overloaded (==) operator interface function
! Subprogram not used !
! Subprogram not used !     The arguments are:
! Subprogram not used !     \begin{description}
! Subprogram not used !     \item[alarm1]
! Subprogram not used !          The first {\tt ESMF\_Alarm} to compare
! Subprogram not used !     \item[alarm2]
! Subprogram not used !          The second {\tt ESMF\_Alarm} to compare
! Subprogram not used !     \end{description}
! Subprogram not used !
! Subprogram not used ! !REQUIREMENTS:  
! Subprogram not used !EOP
! Subprogram not used       CALL wrf_error_fatal( 'ESMF_AlarmEQ not supported ' )
! Subprogram not used       ESMF_AlarmEQ = .FALSE.       ! keep compilers happy
! Subprogram not used       end function ESMF_AlarmEQ

!------------------------------------------------------------------------------
!
! This section defines the overridden Read, Write, Validate and Print methods
! from the ESMF_Base class
!
!------------------------------------------------------------------------------
!BOP
! !IROUTINE: ESMF_AlarmRead - restores an alarm

! !INTERFACE:
! Subprogram not used       subroutine ESMF_AlarmRead(alarm, RingInterval, RingTime, &
! Subprogram not used                            PrevRingTime, StopTime, Ringing, &
! Subprogram not used                            Enabled, ID, rc)
! Subprogram not used 
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used       type(ESMF_Alarm), intent(out) :: alarm
! Subprogram not used       type(ESMF_TimeInterval), intent(in) :: RingInterval
! Subprogram not used       type(ESMF_Time), intent(in) :: RingTime
! Subprogram not used       type(ESMF_Time), intent(in) :: PrevRingTime
! Subprogram not used       type(ESMF_Time), intent(in) :: StopTime
! Subprogram not used       logical, intent(in) :: Ringing
! Subprogram not used       logical, intent(in) :: Enabled
! Subprogram not used       integer, intent(in) :: ID
! Subprogram not used       integer, intent(out), optional :: rc
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !     Restores an {\tt ESMF\_Alarm}
! Subprogram not used !
! Subprogram not used !     The arguments are:
! Subprogram not used !     \begin{description}
! Subprogram not used !     \item[alarm]
! Subprogram not used !          The object instance to restore
! Subprogram not used !     \item[RingInterval]
! Subprogram not used !          The ring interval for repeating alarms
! Subprogram not used !     \item[RingTime]
! Subprogram not used !          Ring time for one-shot or first repeating alarm
! Subprogram not used !     \item[PrevRingTime]
! Subprogram not used !          The {\tt ESMF\_Alarm}'s previous ring time
! Subprogram not used !     \item[StopTime]
! Subprogram not used !          Stop time for repeating alarms
! Subprogram not used !     \item[Ringing]
! Subprogram not used !          The {\tt ESMF\_Alarm}'s ringing state
! Subprogram not used !     \item[Enabled]
! Subprogram not used !          {\tt ESMF\_Alarm} enabled/disabled
! Subprogram not used !     \item[ID]
! Subprogram not used !          The {\tt ESMF\_Alarm}'s ID
! Subprogram not used !     \item[{[rc]}]
! Subprogram not used !          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
! Subprogram not used !     \end{description}
! Subprogram not used !
! Subprogram not used ! !REQUIREMENTS:
! Subprogram not used !EOP
! Subprogram not used       CALL wrf_error_fatal( 'ESMF_AlarmRead not supported' )
! Subprogram not used       end subroutine ESMF_AlarmRead

!------------------------------------------------------------------------------
!BOP
! !IROUTINE: ESMF_AlarmWrite - saves an alarm

! !INTERFACE:
! Subprogram not used       subroutine ESMF_AlarmWrite(alarm, RingInterval, RingTime, &
! Subprogram not used                             PrevRingTime, StopTime, Ringing, &
! Subprogram not used                             Enabled, ID, rc)
! Subprogram not used 
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used       type(ESMF_Alarm), intent(in) :: alarm
! Subprogram not used       type(ESMF_TimeInterval), intent(out) :: RingInterval
! Subprogram not used       type(ESMF_Time), intent(out) :: RingTime
! Subprogram not used       type(ESMF_Time), intent(out) :: PrevRingTime
! Subprogram not used       type(ESMF_Time), intent(out) :: StopTime
! Subprogram not used       logical, intent(out) :: Ringing
! Subprogram not used       logical, intent(out) :: Enabled
! Subprogram not used       integer, intent(out) :: ID
! Subprogram not used       integer, intent(out), optional :: rc
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !     Saves an {\tt ESMF\_Alarm}
! Subprogram not used !
! Subprogram not used !     The arguments are:
! Subprogram not used !     \begin{description}
! Subprogram not used !     \item[alarm]
! Subprogram not used !          The object instance to save
! Subprogram not used !     \item[RingInterval]
! Subprogram not used !          Ring interval for repeating alarms
! Subprogram not used !     \item[RingTime]
! Subprogram not used !          Ring time for one-shot or first repeating alarm
! Subprogram not used !     \item[PrevRingTime]
! Subprogram not used !          The {\tt ESMF\_Alarm}'s previous ring time
! Subprogram not used !     \item[StopTime]
! Subprogram not used !          Stop time for repeating alarms
! Subprogram not used !     \item[Ringing]
! Subprogram not used !          The {\tt ESMF\_Alarm}'s ringing state
! Subprogram not used !     \item[Enabled]
! Subprogram not used !          {\tt ESMF\_Alarm} enabled/disabled
! Subprogram not used !     \item[ID]
! Subprogram not used !          The {\tt ESMF\_Alarm}'s ID
! Subprogram not used !     \item[{[rc]}]
! Subprogram not used !          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
! Subprogram not used !     \end{description}
! Subprogram not used !
! Subprogram not used ! !REQUIREMENTS:
! Subprogram not used !EOP
! Subprogram not used       CALL wrf_error_fatal( 'ESMF_AlarmWrite not supported' )
! Subprogram not used       end subroutine ESMF_AlarmWrite

!------------------------------------------------------------------------------
!BOP
! !IROUTINE:  ESMF_AlarmValidate - Validate an Alarm's properties

! !INTERFACE:
! Subprogram not used       subroutine ESMF_AlarmValidate(alarm, opts, rc)
! Subprogram not used 
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used       type(ESMF_Alarm), intent(in) :: alarm
! Subprogram not used       character (len=*), intent(in), optional :: opts
! Subprogram not used       integer, intent(out), optional :: rc
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !     Perform a validation check on a {\tt ESMF\_Alarm}'s properties
! Subprogram not used !
! Subprogram not used !     The arguments are:  
! Subprogram not used !     \begin{description}
! Subprogram not used !     \item[alarm]
! Subprogram not used !          {\tt ESMF\_Alarm} to validate
! Subprogram not used !     \item[{[opts]}]
! Subprogram not used !          Validate options
! Subprogram not used !     \item[{[rc]}]
! Subprogram not used !          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
! Subprogram not used !     \end{description} 
! Subprogram not used !
! Subprogram not used ! !REQUIREMENTS:
! Subprogram not used !     TMGn.n.n
! Subprogram not used !EOP
! Subprogram not used       CALL wrf_error_fatal( 'ESMF_AlarmValidate not supported' )
! Subprogram not used       end subroutine ESMF_AlarmValidate

!------------------------------------------------------------------------------
!BOP
! !IROUTINE:  ESMF_AlarmPrint - Print out an Alarm's properties

! !INTERFACE:
! Subprogram not used       subroutine ESMF_AlarmPrint(alarm, opts, rc)
! Subprogram not used 
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used       type(ESMF_Alarm), intent(in) :: alarm
! Subprogram not used       character (len=*), intent(in), optional :: opts
! Subprogram not used       integer, intent(out), optional :: rc
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !     To support testing/debugging, print out a {\tt ESMF\_Alarm}'s
! Subprogram not used !     properties.
! Subprogram not used ! 
! Subprogram not used !     The arguments are:
! Subprogram not used !     \begin{description}
! Subprogram not used !     \item[alarm]
! Subprogram not used !          {\tt ESMF\_Alarm} to print out
! Subprogram not used !     \item[{[opts]}]
! Subprogram not used !          Print options
! Subprogram not used !     \item[{[rc]}]
! Subprogram not used !          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
! Subprogram not used !     \end{description}
! Subprogram not used !
! Subprogram not used ! !REQUIREMENTS:
! Subprogram not used !     TMGn.n.n
! Subprogram not used !EOP
! Subprogram not used       integer :: ierr
! Subprogram not used       type(ESMF_Time) :: ringtime
! Subprogram not used       type(ESMF_Time) :: prevringtime
! Subprogram not used       type(ESMF_TimeInterval) :: ringinterval
! Subprogram not used       character(len=256) :: name
! Subprogram not used 
! Subprogram not used       IF ( ASSOCIATED( alarm%alarmint ) ) THEN
! Subprogram not used         IF ( alarm%alarmint%RingTimeSet )THEN
! Subprogram not used           call ESMF_AlarmGet( alarm, name=name, ringtime=ringtime, &
! Subprogram not used                               prevringtime=prevringtime, rc=ierr )
! Subprogram not used           IF ( PRESENT(rc) .AND. (ierr /= ESMF_SUCCESS) )THEN
! Subprogram not used              rc = ierr
! Subprogram not used           END IF 
! Subprogram not used           print *, 'Alarm name: ', trim(name)
! Subprogram not used           print *, 'Next ring time'
! Subprogram not used           call ESMF_TimePrint( ringtime )
! Subprogram not used           print *, 'Previous ring time'
! Subprogram not used           call ESMF_TimePrint( prevringtime )
! Subprogram not used         END IF
! Subprogram not used         IF ( alarm%alarmint%RingIntervalSet )THEN
! Subprogram not used           call ESMF_AlarmGet( alarm, ringinterval=ringinterval, rc=ierr )
! Subprogram not used           IF ( PRESENT(rc) .AND. (ierr /= ESMF_SUCCESS) )THEN
! Subprogram not used              rc = ierr
! Subprogram not used           END IF 
! Subprogram not used           print *, 'Ring Interval'
! Subprogram not used           call ESMF_TimeIntervalPrint( ringinterval )
! Subprogram not used         END IF
! Subprogram not used       END IF
! Subprogram not used 
! Subprogram not used       end subroutine ESMF_AlarmPrint

!------------------------------------------------------------------------------

      end module ESMF_AlarmMod
