! Define flags for compilers supporting Fortran 2008 intrinsics
! HAVE_GAMMA_INTRINSICS: gamma and log_gamma
! HAVE_ERF_INTRINSICS: erf, erfc, and erfc_scaled
! erfc_scaled(x) = (exp(x**2)*erfc(x))

! Use this flag for compilers that don't have real intrinsics, but link in
! a library for you.
! HAVE_ERF_EXTERNALS: erf and erfc

! These compilers have the intrinsics.

! PGI has external erf/derf and erfc/derfc, and will link them for you, but
! it does not consider them "intrinsics" right now.

! As of 5.3.1, NAG does not have any of these.

module shr_spfn_mod
! Module for common mathematical functions

! This #ifdef is to allow the module to be compiled with no dependencies,
! even on shr_kind_mod.
use shr_kind_mod, only: &
     r4 => shr_kind_r4, &
     r8 => shr_kind_r8
use shr_const_mod, only: &
     pi => shr_const_pi

implicit none
private
save


! Error functions
public :: shr_spfn_erf
public :: shr_spfn_erfc
public :: shr_spfn_erfc_scaled

interface shr_spfn_erf
   module procedure shr_spfn_erf_r4
   module procedure shr_spfn_erf_r8
end interface

interface shr_spfn_erfc
   module procedure shr_spfn_erfc_r4
   module procedure shr_spfn_erfc_r8
end interface

interface shr_spfn_erfc_scaled
   module procedure shr_spfn_erfc_scaled_r4
   module procedure shr_spfn_erfc_scaled_r8
end interface

! Gamma functions
! Note that we lack an implementation of log_gamma, but we do have an
! implementation of the upper incomplete gamma function, which is not in
! Fortran 2008.

! Note also that this gamma function is only for double precision. We
! haven't needed an r4 version yet.

public :: shr_spfn_gamma
public :: shr_spfn_igamma

interface shr_spfn_gamma
   module procedure shr_spfn_gamma_r8
end interface

! Mathematical constants
! sqrt(pi)
real(r8), parameter :: sqrtpi = 1.77245385090551602729_r8

! Define machine-specific constants needed in this module.
! These were used by the original gamma and calerf functions to guarantee
! safety against overflow, and precision, on many different machines.

! By defining the constants in this way, we assume that 1/xmin is
! representable (i.e. does not overflow the real type). This assumption was
! not in the original code, but is valid for IEEE single and double
! precision.

! Double precision
!---------------------------------------------------------------------
! Machine epsilon
real(r8), parameter :: epsr8 = epsilon(1._r8)
! "Huge" value is returned when actual value would be infinite.
real(r8), parameter :: xinfr8 = huge(1._r8)
! Smallest normal value.
real(r8), parameter :: xminr8 = tiny(1._r8)
! Largest number that, when added to 1., yields 1.
real(r8), parameter :: xsmallr8 = epsr8/2._r8
! Largest argument for which erfcx > 0.
real(r8), parameter :: xmaxr8 = 1._r8/(sqrtpi*xminr8)

! Single precision
!---------------------------------------------------------------------
! Machine epsilon
real(r4), parameter :: epsr4 = epsilon(1._r4)
! "Huge" value is returned when actual value would be infinite.
real(r4), parameter :: xinfr4 = huge(1._r4)
! Smallest normal value.
real(r4), parameter :: xminr4 = tiny(1._r4)
! Largest number that, when added to 1., yields 1.
real(r4), parameter :: xsmallr4 = epsr4/2._r4
! Largest argument for which erfcx > 0.
real(r4), parameter :: xmaxr4 = 1._r4/(real(sqrtpi,r4)*xminr4)


! For gamma/igamma
! Approximate value of largest acceptable argument to gamma,
! for IEEE double-precision.
real(r8), parameter :: xbig_gamma = 171.624_r8

contains

! Wrapper functions for erf
! Subprogram not used function shr_spfn_erf_r4(x) result(res)
! Subprogram not used   real(r4), intent(in) :: x
! Subprogram not used   real(r4) :: res
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   ! Call intrinsic erf.
! Subprogram not used   intrinsic erf
! Subprogram not used   res = erf(x)
! Subprogram not used 
! Subprogram not used end function shr_spfn_erf_r4

function shr_spfn_erf_r8(x) result(res)
  real(r8), intent(in) :: x
  real(r8) :: res


  ! Call intrinsic erf.
  intrinsic erf
  res = erf(x)

end function shr_spfn_erf_r8

! Wrapper functions for erfc
! Subprogram not used function shr_spfn_erfc_r4(x) result(res)
! Subprogram not used   real(r4), intent(in) :: x
! Subprogram not used   real(r4) :: res
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   ! Call intrinsic erfc.
! Subprogram not used   intrinsic erfc
! Subprogram not used   res = erfc(x)
! Subprogram not used 
! Subprogram not used end function shr_spfn_erfc_r4

function shr_spfn_erfc_r8(x) result(res)
  real(r8), intent(in) :: x
  real(r8) :: res


  ! Call intrinsic erfc.
  intrinsic erfc
  res = erfc(x)

end function shr_spfn_erfc_r8

! Wrapper functions for erfc_scaled
! Subprogram not used function shr_spfn_erfc_scaled_r4(x) result(res)
! Subprogram not used   real(r4), intent(in) :: x
! Subprogram not used   real(r4) :: res
! Subprogram not used 
! Subprogram not used   ! Call intrinsic erfc_scaled.
! Subprogram not used   intrinsic erfc_scaled
! Subprogram not used   res = erfc_scaled(x)
! Subprogram not used 
! Subprogram not used end function shr_spfn_erfc_scaled_r4

! Subprogram not used function shr_spfn_erfc_scaled_r8(x) result(res)
! Subprogram not used   real(r8), intent(in) :: x
! Subprogram not used   real(r8) :: res
! Subprogram not used 
! Subprogram not used   ! Call intrinsic erfc_scaled.
! Subprogram not used   intrinsic erfc_scaled
! Subprogram not used   res = erfc_scaled(x)
! Subprogram not used 
! Subprogram not used end function shr_spfn_erfc_scaled_r8

elemental function shr_spfn_gamma_r8(x) result(res)
  real(r8), intent(in) :: x
  real(r8) :: res

  ! Call intrinsic gamma.
  intrinsic gamma
  res = gamma(x)

end function shr_spfn_gamma_r8

!------------------------------------------------------------------
!
! 6 December 2006 -- B. Eaton
! The following comments are from the original version of CALERF.
! The only changes in implementing this module are that the function
! names previously used for the single precision versions have been
! adopted for the new generic interfaces.  To support these interfaces
! there is now both a single precision version (calerf_r4) and a
! double precision version (calerf_r8) of CALERF below.  These versions
! are hardcoded to use IEEE arithmetic.
!
!------------------------------------------------------------------
!
! This packet evaluates  erf(x),  erfc(x),  and  exp(x*x)*erfc(x)
!   for a real argument  x.  It contains three FUNCTION type
!   subprograms: ERF, ERFC, and ERFCX (or ERF_R8, ERFC_R8, and ERFCX_R8),
!   and one SUBROUTINE type subprogram, CALERF.  The calling
!   statements for the primary entries are:
!
!                   Y=ERF(X)     (or   Y=ERF_R8(X)),
!
!                   Y=ERFC(X)    (or   Y=ERFC_R8(X)),
!   and
!                   Y=ERFCX(X)   (or   Y=ERFCX_R8(X)).
!
!   The routine  CALERF  is intended for internal packet use only,
!   all computations within the packet being concentrated in this
!   routine.  The function subprograms invoke  CALERF  with the
!   statement
!
!          CALL CALERF(ARG,RESULT,JINT)
!
!   where the parameter usage is as follows
!
!      Function                     Parameters for CALERF
!       call              ARG                  Result          JINT
!
!     ERF(ARG)      ANY REAL ARGUMENT         ERF(ARG)          0
!     ERFC(ARG)     ABS(ARG) .LT. XBIG        ERFC(ARG)         1
!     ERFCX(ARG)    XNEG .LT. ARG .LT. XMAX   ERFCX(ARG)        2
!
!   The main computation evaluates near-minimax approximations
!   from "Rational Chebyshev approximations for the error function"
!   by W. J. Cody, Math. Comp., 1969, PP. 631-638.  This
!   transportable program uses rational functions that theoretically
!   approximate  erf(x)  and  erfc(x)  to at least 18 significant
!   decimal digits.  The accuracy achieved depends on the arithmetic
!   system, the compiler, the intrinsic functions, and proper
!   selection of the machine-dependent constants.
!
!*******************************************************************
!*******************************************************************
!
! Explanation of machine-dependent constants
!
!   XMIN   = the smallest positive floating-point number.
!   XINF   = the largest positive finite floating-point number.
!   XNEG   = the largest negative argument acceptable to ERFCX;
!            the negative of the solution to the equation
!            2*exp(x*x) = XINF.
!   XSMALL = argument below which erf(x) may be represented by
!            2*x/sqrt(pi)  and above which  x*x  will not underflow.
!            A conservative value is the largest machine number X
!            such that   1.0 + X = 1.0   to machine precision.
!   XBIG   = largest argument acceptable to ERFC;  solution to
!            the equation:  W(x) * (1-0.5/x**2) = XMIN,  where
!            W(x) = exp(-x*x)/[x*sqrt(pi)].
!   XHUGE  = argument above which  1.0 - 1/(2*x*x) = 1.0  to
!            machine precision.  A conservative value is
!            1/[2*sqrt(XSMALL)]
!   XMAX   = largest acceptable argument to ERFCX; the minimum
!            of XINF and 1/[sqrt(pi)*XMIN].
!
!   Approximate values for some important machines are:
!
!                          XMIN       XINF        XNEG     XSMALL
!
!  CDC 7600      (S.P.)  3.13E-294   1.26E+322   -27.220  7.11E-15
!  CRAY-1        (S.P.)  4.58E-2467  5.45E+2465  -75.345  7.11E-15
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)  1.18E-38    3.40E+38     -9.382  5.96E-8
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)  2.23D-308   1.79D+308   -26.628  1.11D-16
!  IBM 195       (D.P.)  5.40D-79    7.23E+75    -13.190  1.39D-17
!  UNIVAC 1108   (D.P.)  2.78D-309   8.98D+307   -26.615  1.73D-18
!  VAX D-Format  (D.P.)  2.94D-39    1.70D+38     -9.345  1.39D-17
!  VAX G-Format  (D.P.)  5.56D-309   8.98D+307   -26.615  1.11D-16
!
!
!                          XBIG       XHUGE       XMAX
!
!  CDC 7600      (S.P.)  25.922      8.39E+6     1.80X+293
!  CRAY-1        (S.P.)  75.326      8.39E+6     5.45E+2465
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)   9.194      2.90E+3     4.79E+37
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)  26.543      6.71D+7     2.53D+307
!  IBM 195       (D.P.)  13.306      1.90D+8     7.23E+75
!  UNIVAC 1108   (D.P.)  26.582      5.37D+8     8.98D+307
!  VAX D-Format  (D.P.)   9.269      1.90D+8     1.70D+38
!  VAX G-Format  (D.P.)  26.569      6.71D+7     8.98D+307
!
!*******************************************************************
!*******************************************************************
!
! Error returns
!
!  The program returns  ERFC = 0      for  ARG .GE. XBIG;
!
!                       ERFCX = XINF  for  ARG .LT. XNEG;
!      and
!                       ERFCX = 0     for  ARG .GE. XMAX.
!
!
! Intrinsic functions required are:
!
!     ABS, AINT, EXP
!
!
!  Author: W. J. Cody
!          Mathematics and Computer Science Division
!          Argonne National Laboratory
!          Argonne, IL 60439
!
!  Latest modification: March 19, 1990
!
!------------------------------------------------------------------

! Subprogram not used SUBROUTINE CALERF_r8(ARG, RESULT, JINT)
! Subprogram not used 
! Subprogram not used    !------------------------------------------------------------------
! Subprogram not used    !  This version uses 8-byte reals
! Subprogram not used    !------------------------------------------------------------------
! Subprogram not used    integer, parameter :: rk = r8
! Subprogram not used 
! Subprogram not used    ! arguments
! Subprogram not used    real(rk), intent(in)  :: arg
! Subprogram not used    integer,  intent(in)  :: jint
! Subprogram not used    real(rk), intent(out) :: result
! Subprogram not used 
! Subprogram not used    ! local variables
! Subprogram not used    INTEGER :: I
! Subprogram not used 
! Subprogram not used    real(rk) :: X, Y, YSQ, XNUM, XDEN, DEL
! Subprogram not used 
! Subprogram not used    !------------------------------------------------------------------
! Subprogram not used    !  Mathematical constants
! Subprogram not used    !------------------------------------------------------------------
! Subprogram not used    real(rk), parameter :: ZERO   = 0.0E0_rk
! Subprogram not used    real(rk), parameter :: FOUR   = 4.0E0_rk
! Subprogram not used    real(rk), parameter :: ONE    = 1.0E0_rk
! Subprogram not used    real(rk), parameter :: HALF   = 0.5E0_rk
! Subprogram not used    real(rk), parameter :: TWO    = 2.0E0_rk
! Subprogram not used    ! 1/sqrt(pi)
! Subprogram not used    real(rk), parameter :: SQRPI  = 5.6418958354775628695E-1_rk
! Subprogram not used    real(rk), parameter :: THRESH = 0.46875E0_rk
! Subprogram not used    real(rk), parameter :: SIXTEN = 16.0E0_rk
! Subprogram not used 
! Subprogram not used    !------------------------------------------------------------------
! Subprogram not used    !  Machine-dependent constants: IEEE double precision values
! Subprogram not used    !------------------------------------------------------------------
! Subprogram not used    real(rk), parameter :: XNEG   = -26.628E0_r8
! Subprogram not used    real(rk), parameter :: XBIG   =  26.543E0_r8
! Subprogram not used    real(rk), parameter :: XHUGE  =   6.71E7_r8
! Subprogram not used 
! Subprogram not used    !------------------------------------------------------------------
! Subprogram not used    !  Coefficients for approximation to  erf  in first interval
! Subprogram not used    !------------------------------------------------------------------
! Subprogram not used    real(rk), parameter :: A(5) = (/ 3.16112374387056560E00_rk, 1.13864154151050156E02_rk, &
! Subprogram not used                                     3.77485237685302021E02_rk, 3.20937758913846947E03_rk, &
! Subprogram not used                                     1.85777706184603153E-1_rk /)
! Subprogram not used    real(rk), parameter :: B(4) = (/ 2.36012909523441209E01_rk, 2.44024637934444173E02_rk, &
! Subprogram not used                                     1.28261652607737228E03_rk, 2.84423683343917062E03_rk /)
! Subprogram not used 
! Subprogram not used    !------------------------------------------------------------------
! Subprogram not used    !  Coefficients for approximation to  erfc  in second interval
! Subprogram not used    !------------------------------------------------------------------
! Subprogram not used    real(rk), parameter :: C(9) = (/ 5.64188496988670089E-1_rk, 8.88314979438837594E00_rk, &
! Subprogram not used                                     6.61191906371416295E01_rk, 2.98635138197400131E02_rk, &
! Subprogram not used                                     8.81952221241769090E02_rk, 1.71204761263407058E03_rk, &
! Subprogram not used                                     2.05107837782607147E03_rk, 1.23033935479799725E03_rk, &
! Subprogram not used                                     2.15311535474403846E-8_rk /)
! Subprogram not used    real(rk), parameter :: D(8) = (/ 1.57449261107098347E01_rk, 1.17693950891312499E02_rk, &
! Subprogram not used                                     5.37181101862009858E02_rk, 1.62138957456669019E03_rk, &
! Subprogram not used                                     3.29079923573345963E03_rk, 4.36261909014324716E03_rk, &
! Subprogram not used                                     3.43936767414372164E03_rk, 1.23033935480374942E03_rk /)
! Subprogram not used 
! Subprogram not used    !------------------------------------------------------------------
! Subprogram not used    !  Coefficients for approximation to  erfc  in third interval
! Subprogram not used    !------------------------------------------------------------------
! Subprogram not used    real(rk), parameter :: P(6) = (/ 3.05326634961232344E-1_rk, 3.60344899949804439E-1_rk, &
! Subprogram not used                                     1.25781726111229246E-1_rk, 1.60837851487422766E-2_rk, &
! Subprogram not used                                     6.58749161529837803E-4_rk, 1.63153871373020978E-2_rk /)
! Subprogram not used    real(rk), parameter :: Q(5) = (/ 2.56852019228982242E00_rk, 1.87295284992346047E00_rk, &
! Subprogram not used                                     5.27905102951428412E-1_rk, 6.05183413124413191E-2_rk, &
! Subprogram not used                                     2.33520497626869185E-3_rk /)
! Subprogram not used 
! Subprogram not used    !------------------------------------------------------------------
! Subprogram not used    X = ARG
! Subprogram not used    Y = ABS(X)
! Subprogram not used    IF (Y .LE. THRESH) THEN
! Subprogram not used       !------------------------------------------------------------------
! Subprogram not used       !  Evaluate  erf  for  |X| <= 0.46875
! Subprogram not used       !------------------------------------------------------------------
! Subprogram not used       YSQ = ZERO
! Subprogram not used       IF (Y .GT. XSMALLR8) YSQ = Y * Y
! Subprogram not used       XNUM = A(5)*YSQ
! Subprogram not used       XDEN = YSQ
! Subprogram not used       DO I = 1, 3
! Subprogram not used          XNUM = (XNUM + A(I)) * YSQ
! Subprogram not used          XDEN = (XDEN + B(I)) * YSQ
! Subprogram not used       end do
! Subprogram not used       RESULT = X * (XNUM + A(4)) / (XDEN + B(4))
! Subprogram not used       IF (JINT .NE. 0) RESULT = ONE - RESULT
! Subprogram not used       IF (JINT .EQ. 2) RESULT = EXP(YSQ) * RESULT
! Subprogram not used       GO TO 80
! Subprogram not used    ELSE IF (Y .LE. FOUR) THEN
! Subprogram not used       !------------------------------------------------------------------
! Subprogram not used       !  Evaluate  erfc  for 0.46875 <= |X| <= 4.0
! Subprogram not used       !------------------------------------------------------------------
! Subprogram not used       XNUM = C(9)*Y
! Subprogram not used       XDEN = Y
! Subprogram not used       DO I = 1, 7
! Subprogram not used          XNUM = (XNUM + C(I)) * Y
! Subprogram not used          XDEN = (XDEN + D(I)) * Y
! Subprogram not used       end do
! Subprogram not used       RESULT = (XNUM + C(8)) / (XDEN + D(8))
! Subprogram not used       IF (JINT .NE. 2) THEN
! Subprogram not used          YSQ = AINT(Y*SIXTEN)/SIXTEN
! Subprogram not used          DEL = (Y-YSQ)*(Y+YSQ)
! Subprogram not used          RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT
! Subprogram not used       END IF
! Subprogram not used    ELSE
! Subprogram not used       !------------------------------------------------------------------
! Subprogram not used       !  Evaluate  erfc  for |X| > 4.0
! Subprogram not used       !------------------------------------------------------------------
! Subprogram not used       RESULT = ZERO
! Subprogram not used       IF (Y .GE. XBIG) THEN
! Subprogram not used          IF ((JINT .NE. 2) .OR. (Y .GE. XMAXR8)) GO TO 30
! Subprogram not used          IF (Y .GE. XHUGE) THEN
! Subprogram not used             RESULT = SQRPI / Y
! Subprogram not used             GO TO 30
! Subprogram not used          END IF
! Subprogram not used       END IF
! Subprogram not used       YSQ = ONE / (Y * Y)
! Subprogram not used       XNUM = P(6)*YSQ
! Subprogram not used       XDEN = YSQ
! Subprogram not used       DO I = 1, 4
! Subprogram not used          XNUM = (XNUM + P(I)) * YSQ
! Subprogram not used          XDEN = (XDEN + Q(I)) * YSQ
! Subprogram not used       end do
! Subprogram not used       RESULT = YSQ *(XNUM + P(5)) / (XDEN + Q(5))
! Subprogram not used       RESULT = (SQRPI -  RESULT) / Y
! Subprogram not used       IF (JINT .NE. 2) THEN
! Subprogram not used          YSQ = AINT(Y*SIXTEN)/SIXTEN
! Subprogram not used          DEL = (Y-YSQ)*(Y+YSQ)
! Subprogram not used          RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT
! Subprogram not used       END IF
! Subprogram not used    END IF
! Subprogram not used 30 continue
! Subprogram not used    !------------------------------------------------------------------
! Subprogram not used    !  Fix up for negative argument, erf, etc.
! Subprogram not used    !------------------------------------------------------------------
! Subprogram not used    IF (JINT .EQ. 0) THEN
! Subprogram not used       RESULT = (HALF - RESULT) + HALF
! Subprogram not used       IF (X .LT. ZERO) RESULT = -RESULT
! Subprogram not used    ELSE IF (JINT .EQ. 1) THEN
! Subprogram not used       IF (X .LT. ZERO) RESULT = TWO - RESULT
! Subprogram not used    ELSE
! Subprogram not used       IF (X .LT. ZERO) THEN
! Subprogram not used          IF (X .LT. XNEG) THEN
! Subprogram not used             RESULT = XINFR8
! Subprogram not used          ELSE
! Subprogram not used             YSQ = AINT(X*SIXTEN)/SIXTEN
! Subprogram not used             DEL = (X-YSQ)*(X+YSQ)
! Subprogram not used             Y = EXP(YSQ*YSQ) * EXP(DEL)
! Subprogram not used             RESULT = (Y+Y) - RESULT
! Subprogram not used          END IF
! Subprogram not used       END IF
! Subprogram not used    END IF
! Subprogram not used 80 continue
! Subprogram not used end SUBROUTINE CALERF_r8

!------------------------------------------------------------------------------------------

! Subprogram not used SUBROUTINE CALERF_r4(ARG, RESULT, JINT)
! Subprogram not used 
! Subprogram not used    !------------------------------------------------------------------
! Subprogram not used    !  This version uses 4-byte reals
! Subprogram not used    !------------------------------------------------------------------
! Subprogram not used    integer, parameter :: rk = r4
! Subprogram not used 
! Subprogram not used    ! arguments
! Subprogram not used    real(rk), intent(in)  :: arg
! Subprogram not used    integer,  intent(in)  :: jint
! Subprogram not used    real(rk), intent(out) :: result
! Subprogram not used 
! Subprogram not used    ! local variables
! Subprogram not used    INTEGER :: I
! Subprogram not used 
! Subprogram not used    real(rk) :: X, Y, YSQ, XNUM, XDEN, DEL
! Subprogram not used 
! Subprogram not used    !------------------------------------------------------------------
! Subprogram not used    !  Mathematical constants
! Subprogram not used    !------------------------------------------------------------------
! Subprogram not used    real(rk), parameter :: ZERO   = 0.0E0_rk
! Subprogram not used    real(rk), parameter :: FOUR   = 4.0E0_rk
! Subprogram not used    real(rk), parameter :: ONE    = 1.0E0_rk
! Subprogram not used    real(rk), parameter :: HALF   = 0.5E0_rk
! Subprogram not used    real(rk), parameter :: TWO    = 2.0E0_rk
! Subprogram not used    ! 1/sqrt(pi)
! Subprogram not used    real(rk), parameter :: SQRPI  = 5.6418958354775628695E-1_rk
! Subprogram not used    real(rk), parameter :: THRESH = 0.46875E0_rk
! Subprogram not used    real(rk), parameter :: SIXTEN = 16.0E0_rk
! Subprogram not used 
! Subprogram not used    !------------------------------------------------------------------
! Subprogram not used    !  Machine-dependent constants: IEEE single precision values
! Subprogram not used    !------------------------------------------------------------------
! Subprogram not used    real(rk), parameter :: XNEG   = -9.382E0_r4
! Subprogram not used    real(rk), parameter :: XBIG   =  9.194E0_r4
! Subprogram not used    real(rk), parameter :: XHUGE  =  2.90E3_r4
! Subprogram not used 
! Subprogram not used    !------------------------------------------------------------------
! Subprogram not used    !  Coefficients for approximation to  erf  in first interval
! Subprogram not used    !------------------------------------------------------------------
! Subprogram not used    real(rk), parameter :: A(5) = (/ 3.16112374387056560E00_rk, 1.13864154151050156E02_rk, &
! Subprogram not used                                     3.77485237685302021E02_rk, 3.20937758913846947E03_rk, &
! Subprogram not used                                     1.85777706184603153E-1_rk /)
! Subprogram not used    real(rk), parameter :: B(4) = (/ 2.36012909523441209E01_rk, 2.44024637934444173E02_rk, &
! Subprogram not used                                     1.28261652607737228E03_rk, 2.84423683343917062E03_rk /)
! Subprogram not used 
! Subprogram not used    !------------------------------------------------------------------
! Subprogram not used    !  Coefficients for approximation to  erfc  in second interval
! Subprogram not used    !------------------------------------------------------------------
! Subprogram not used    real(rk), parameter :: C(9) = (/ 5.64188496988670089E-1_rk, 8.88314979438837594E00_rk, &
! Subprogram not used                                     6.61191906371416295E01_rk, 2.98635138197400131E02_rk, &
! Subprogram not used                                     8.81952221241769090E02_rk, 1.71204761263407058E03_rk, &
! Subprogram not used                                     2.05107837782607147E03_rk, 1.23033935479799725E03_rk, &
! Subprogram not used                                     2.15311535474403846E-8_rk /)
! Subprogram not used    real(rk), parameter :: D(8) = (/ 1.57449261107098347E01_rk, 1.17693950891312499E02_rk, &
! Subprogram not used                                     5.37181101862009858E02_rk, 1.62138957456669019E03_rk, &
! Subprogram not used                                     3.29079923573345963E03_rk, 4.36261909014324716E03_rk, &
! Subprogram not used                                     3.43936767414372164E03_rk, 1.23033935480374942E03_rk /)
! Subprogram not used 
! Subprogram not used    !------------------------------------------------------------------
! Subprogram not used    !  Coefficients for approximation to  erfc  in third interval
! Subprogram not used    !------------------------------------------------------------------
! Subprogram not used    real(rk), parameter :: P(6) = (/ 3.05326634961232344E-1_rk, 3.60344899949804439E-1_rk, &
! Subprogram not used                                     1.25781726111229246E-1_rk, 1.60837851487422766E-2_rk, &
! Subprogram not used                                     6.58749161529837803E-4_rk, 1.63153871373020978E-2_rk /)
! Subprogram not used    real(rk), parameter :: Q(5) = (/ 2.56852019228982242E00_rk, 1.87295284992346047E00_rk, &
! Subprogram not used                                     5.27905102951428412E-1_rk, 6.05183413124413191E-2_rk, &
! Subprogram not used                                     2.33520497626869185E-3_rk /)
! Subprogram not used 
! Subprogram not used    !------------------------------------------------------------------
! Subprogram not used    X = ARG
! Subprogram not used    Y = ABS(X)
! Subprogram not used    IF (Y .LE. THRESH) THEN
! Subprogram not used       !------------------------------------------------------------------
! Subprogram not used       !  Evaluate  erf  for  |X| <= 0.46875
! Subprogram not used       !------------------------------------------------------------------
! Subprogram not used       YSQ = ZERO
! Subprogram not used       IF (Y .GT. XSMALLR4) YSQ = Y * Y
! Subprogram not used       XNUM = A(5)*YSQ
! Subprogram not used       XDEN = YSQ
! Subprogram not used       DO I = 1, 3
! Subprogram not used          XNUM = (XNUM + A(I)) * YSQ
! Subprogram not used          XDEN = (XDEN + B(I)) * YSQ
! Subprogram not used       end do
! Subprogram not used       RESULT = X * (XNUM + A(4)) / (XDEN + B(4))
! Subprogram not used       IF (JINT .NE. 0) RESULT = ONE - RESULT
! Subprogram not used       IF (JINT .EQ. 2) RESULT = EXP(YSQ) * RESULT
! Subprogram not used       GO TO 80
! Subprogram not used    ELSE IF (Y .LE. FOUR) THEN
! Subprogram not used       !------------------------------------------------------------------
! Subprogram not used       !  Evaluate  erfc  for 0.46875 <= |X| <= 4.0
! Subprogram not used       !------------------------------------------------------------------
! Subprogram not used       XNUM = C(9)*Y
! Subprogram not used       XDEN = Y
! Subprogram not used       DO I = 1, 7
! Subprogram not used          XNUM = (XNUM + C(I)) * Y
! Subprogram not used          XDEN = (XDEN + D(I)) * Y
! Subprogram not used       end do
! Subprogram not used       RESULT = (XNUM + C(8)) / (XDEN + D(8))
! Subprogram not used       IF (JINT .NE. 2) THEN
! Subprogram not used          YSQ = AINT(Y*SIXTEN)/SIXTEN
! Subprogram not used          DEL = (Y-YSQ)*(Y+YSQ)
! Subprogram not used          RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT
! Subprogram not used       END IF
! Subprogram not used    ELSE
! Subprogram not used       !------------------------------------------------------------------
! Subprogram not used       !  Evaluate  erfc  for |X| > 4.0
! Subprogram not used       !------------------------------------------------------------------
! Subprogram not used       RESULT = ZERO
! Subprogram not used       IF (Y .GE. XBIG) THEN
! Subprogram not used          IF ((JINT .NE. 2) .OR. (Y .GE. XMAXR4)) GO TO 30
! Subprogram not used          IF (Y .GE. XHUGE) THEN
! Subprogram not used             RESULT = SQRPI / Y
! Subprogram not used             GO TO 30
! Subprogram not used          END IF
! Subprogram not used       END IF
! Subprogram not used       YSQ = ONE / (Y * Y)
! Subprogram not used       XNUM = P(6)*YSQ
! Subprogram not used       XDEN = YSQ
! Subprogram not used       DO I = 1, 4
! Subprogram not used          XNUM = (XNUM + P(I)) * YSQ
! Subprogram not used          XDEN = (XDEN + Q(I)) * YSQ
! Subprogram not used       end do
! Subprogram not used       RESULT = YSQ *(XNUM + P(5)) / (XDEN + Q(5))
! Subprogram not used       RESULT = (SQRPI -  RESULT) / Y
! Subprogram not used       IF (JINT .NE. 2) THEN
! Subprogram not used          YSQ = AINT(Y*SIXTEN)/SIXTEN
! Subprogram not used          DEL = (Y-YSQ)*(Y+YSQ)
! Subprogram not used          RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT
! Subprogram not used       END IF
! Subprogram not used    END IF
! Subprogram not used 30 continue
! Subprogram not used    !------------------------------------------------------------------
! Subprogram not used    !  Fix up for negative argument, erf, etc.
! Subprogram not used    !------------------------------------------------------------------
! Subprogram not used    IF (JINT .EQ. 0) THEN
! Subprogram not used       RESULT = (HALF - RESULT) + HALF
! Subprogram not used       IF (X .LT. ZERO) RESULT = -RESULT
! Subprogram not used    ELSE IF (JINT .EQ. 1) THEN
! Subprogram not used       IF (X .LT. ZERO) RESULT = TWO - RESULT
! Subprogram not used    ELSE
! Subprogram not used       IF (X .LT. ZERO) THEN
! Subprogram not used          IF (X .LT. XNEG) THEN
! Subprogram not used             RESULT = XINFR4
! Subprogram not used          ELSE
! Subprogram not used             YSQ = AINT(X*SIXTEN)/SIXTEN
! Subprogram not used             DEL = (X-YSQ)*(X+YSQ)
! Subprogram not used             Y = EXP(YSQ*YSQ) * EXP(DEL)
! Subprogram not used             RESULT = (Y+Y) - RESULT
! Subprogram not used          END IF
! Subprogram not used       END IF
! Subprogram not used    END IF
! Subprogram not used 80 continue
! Subprogram not used end SUBROUTINE CALERF_r4

!------------------------------------------------------------------------------------------

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! Subprogram not used pure function shr_spfn_gamma_nonintrinsic_r8(X) result(gamma)
! Subprogram not used 
! Subprogram not used !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Subprogram not used !
! Subprogram not used ! 7 Feb 2013 -- S. Santos
! Subprogram not used ! The following comments are from the original version. Changes have
! Subprogram not used ! been made to update syntax and allow inclusion into this module.
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used ! THIS ROUTINE CALCULATES THE GAMMA FUNCTION FOR A REAL ARGUMENT X.
! Subprogram not used !   COMPUTATION IS BASED ON AN ALGORITHM OUTLINED IN REFERENCE 1.
! Subprogram not used !   THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE THE GAMMA
! Subprogram not used !   FUNCTION TO AT LEAST 20 SIGNIFICANT DECIMAL DIGITS.  COEFFICIENTS
! Subprogram not used !   FOR THE APPROXIMATION OVER THE INTERVAL (1,2) ARE UNPUBLISHED.
! Subprogram not used !   THOSE FOR THE APPROXIMATION FOR X .GE. 12 ARE FROM REFERENCE 2.
! Subprogram not used !   THE ACCURACY ACHIEVED DEPENDS ON THE ARITHMETIC SYSTEM, THE
! Subprogram not used !   COMPILER, THE INTRINSIC FUNCTIONS, AND PROPER SELECTION OF THE
! Subprogram not used !   MACHINE-DEPENDENT CONSTANTS.
! Subprogram not used !
! Subprogram not used !
! Subprogram not used !*******************************************************************
! Subprogram not used !*******************************************************************
! Subprogram not used !
! Subprogram not used ! EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
! Subprogram not used !
! Subprogram not used ! BETA   - RADIX FOR THE FLOATING-POINT REPRESENTATION
! Subprogram not used ! MAXEXP - THE SMALLEST POSITIVE POWER OF BETA THAT OVERFLOWS
! Subprogram not used ! XBIG   - THE LARGEST ARGUMENT FOR WHICH GAMMA(X) IS REPRESENTABLE
! Subprogram not used !          IN THE MACHINE, I.E., THE SOLUTION TO THE EQUATION
! Subprogram not used !                  GAMMA(XBIG) = BETA**MAXEXP
! Subprogram not used ! XINF   - THE LARGEST MACHINE REPRESENTABLE FLOATING-POINT NUMBER;
! Subprogram not used !          APPROXIMATELY BETA**MAXEXP
! Subprogram not used ! EPS    - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
! Subprogram not used !          1.0+EPS .GT. 1.0
! Subprogram not used ! XMININ - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
! Subprogram not used !          1/XMININ IS MACHINE REPRESENTABLE
! Subprogram not used !
! Subprogram not used !     APPROXIMATE VALUES FOR SOME IMPORTANT MACHINES ARE:
! Subprogram not used !
! Subprogram not used !                            BETA       MAXEXP        XBIG
! Subprogram not used !
! Subprogram not used ! CRAY-1         (S.P.)        2         8191        966.961
! Subprogram not used ! CYBER 180/855
! Subprogram not used !   UNDER NOS    (S.P.)        2         1070        177.803
! Subprogram not used ! IEEE (IBM/XT,
! Subprogram not used !   SUN, ETC.)   (S.P.)        2          128        35.040
! Subprogram not used ! IEEE (IBM/XT,
! Subprogram not used !   SUN, ETC.)   (D.P.)        2         1024        171.624
! Subprogram not used ! IBM 3033       (D.P.)       16           63        57.574
! Subprogram not used ! VAX D-FORMAT   (D.P.)        2          127        34.844
! Subprogram not used ! VAX G-FORMAT   (D.P.)        2         1023        171.489
! Subprogram not used !
! Subprogram not used !                            XINF         EPS        XMININ
! Subprogram not used !
! Subprogram not used ! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
! Subprogram not used ! CYBER 180/855
! Subprogram not used !   UNDER NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
! Subprogram not used ! IEEE (IBM/XT,
! Subprogram not used !   SUN, ETC.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
! Subprogram not used ! IEEE (IBM/XT,
! Subprogram not used !   SUN, ETC.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
! Subprogram not used ! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
! Subprogram not used ! VAX D-FORMAT   (D.P.)   1.70D+38     1.39D-17    5.88D-39
! Subprogram not used ! VAX G-FORMAT   (D.P.)   8.98D+307    1.11D-16    1.12D-308
! Subprogram not used !
! Subprogram not used !*******************************************************************
! Subprogram not used !*******************************************************************
! Subprogram not used !
! Subprogram not used ! ERROR RETURNS
! Subprogram not used !
! Subprogram not used !  THE PROGRAM RETURNS THE VALUE XINF FOR SINGULARITIES OR
! Subprogram not used !     WHEN OVERFLOW WOULD OCCUR.  THE COMPUTATION IS BELIEVED
! Subprogram not used !     TO BE FREE OF UNDERFLOW AND OVERFLOW.
! Subprogram not used !
! Subprogram not used !
! Subprogram not used !  INTRINSIC FUNCTIONS REQUIRED ARE:
! Subprogram not used !
! Subprogram not used !     INT, DBLE, EXP, LOG, REAL, SIN
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! REFERENCES:  AN OVERVIEW OF SOFTWARE DEVELOPMENT FOR SPECIAL
! Subprogram not used !              FUNCTIONS   W. J. CODY, LECTURE NOTES IN MATHEMATICS,
! Subprogram not used !              506, NUMERICAL ANALYSIS DUNDEE, 1975, G. A. WATSON
! Subprogram not used !              (ED.), SPRINGER VERLAG, BERLIN, 1976.
! Subprogram not used !
! Subprogram not used !              COMPUTER APPROXIMATIONS, HART, ET. AL., WILEY AND
! Subprogram not used !              SONS, NEW YORK, 1968.
! Subprogram not used !
! Subprogram not used !  LATEST MODIFICATION: OCTOBER 12, 1989
! Subprogram not used !
! Subprogram not used !  AUTHORS: W. J. CODY AND L. STOLTZ
! Subprogram not used !           APPLIED MATHEMATICS DIVISION
! Subprogram not used !           ARGONNE NATIONAL LABORATORY
! Subprogram not used !           ARGONNE, IL 60439
! Subprogram not used !
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used   real(r8), intent(in) :: x
! Subprogram not used   real(r8) :: gamma
! Subprogram not used   real(r8) :: fact, res, sum, xden, xnum, y, y1, ysq, z
! Subprogram not used 
! Subprogram not used   integer :: i, n
! Subprogram not used   logical :: negative_odd
! Subprogram not used 
! Subprogram not used   ! log(2*pi)/2
! Subprogram not used   real(r8), parameter :: logsqrt2pi = 0.9189385332046727417803297E0_r8
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
! Subprogram not used !     APPROXIMATION OVER (1,2).
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used   real(r8), parameter :: P(8) = &
! Subprogram not used        (/-1.71618513886549492533811E+0_r8, 2.47656508055759199108314E+1_r8, &
! Subprogram not used          -3.79804256470945635097577E+2_r8, 6.29331155312818442661052E+2_r8, &
! Subprogram not used           8.66966202790413211295064E+2_r8,-3.14512729688483675254357E+4_r8, &
! Subprogram not used          -3.61444134186911729807069E+4_r8, 6.64561438202405440627855E+4_r8 /)
! Subprogram not used   real(r8), parameter :: Q(8) = &
! Subprogram not used        (/-3.08402300119738975254353E+1_r8, 3.15350626979604161529144E+2_r8, &
! Subprogram not used          -1.01515636749021914166146E+3_r8,-3.10777167157231109440444E+3_r8, &
! Subprogram not used           2.25381184209801510330112E+4_r8, 4.75584627752788110767815E+3_r8, &
! Subprogram not used          -1.34659959864969306392456E+5_r8,-1.15132259675553483497211E+5_r8 /)
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used   real(r8), parameter :: C(7) = &
! Subprogram not used        (/-1.910444077728E-03_r8,          8.4171387781295E-04_r8, &
! Subprogram not used          -5.952379913043012E-04_r8,       7.93650793500350248E-04_r8, &
! Subprogram not used          -2.777777777777681622553E-03_r8, 8.333333333333333331554247E-02_r8, &
! Subprogram not used           5.7083835261E-03_r8 /)
! Subprogram not used 
! Subprogram not used   negative_odd = .false.
! Subprogram not used   fact = 1._r8
! Subprogram not used   n = 0
! Subprogram not used   y = x
! Subprogram not used   if (y <= 0._r8) then
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !  ARGUMENT IS NEGATIVE
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used      y = -x
! Subprogram not used      y1 = aint(y)
! Subprogram not used      res = y - y1
! Subprogram not used      if (res /= 0._r8) then
! Subprogram not used         negative_odd = (y1 /= aint(y1*0.5_r8)*2._r8)
! Subprogram not used         fact = -pi/sin(pi*res)
! Subprogram not used         y = y + 1._r8
! Subprogram not used      else
! Subprogram not used         gamma = xinfr8
! Subprogram not used         return
! Subprogram not used      end if
! Subprogram not used   end if
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !  ARGUMENT IS POSITIVE
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used   if (y < epsr8) then
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !  ARGUMENT .LT. EPS
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used      if (y >= xminr8) then
! Subprogram not used         res = 1._r8/y
! Subprogram not used      else
! Subprogram not used         gamma = xinfr8
! Subprogram not used         return
! Subprogram not used      end if
! Subprogram not used   elseif (y < 12._r8) then
! Subprogram not used      y1 = y
! Subprogram not used      if (y < 1._r8) then
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !  0.0 .LT. ARGUMENT .LT. 1.0
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used         z = y
! Subprogram not used         y = y + 1._r8
! Subprogram not used      else
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used         n = int(y) - 1
! Subprogram not used         y = y - real(n, r8)
! Subprogram not used         z = y - 1._r8
! Subprogram not used      end if
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used      xnum = 0._r8
! Subprogram not used      xden = 1._r8
! Subprogram not used      do i=1,8
! Subprogram not used         xnum = (xnum+P(i))*z
! Subprogram not used         xden = xden*z + Q(i)
! Subprogram not used      end do
! Subprogram not used      res = xnum/xden + 1._r8
! Subprogram not used      if (y1 < y) then
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !  ADJUST RESULT FOR CASE  0.0 .LT. ARGUMENT .LT. 1.0
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used         res = res/y1
! Subprogram not used      elseif (y1 > y) then
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !  ADJUST RESULT FOR CASE  2.0 .LT. ARGUMENT .LT. 12.0
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used         do i = 1,n
! Subprogram not used            res = res*y
! Subprogram not used            y = y + 1._r8
! Subprogram not used         end do
! Subprogram not used      end if
! Subprogram not used   else
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !  EVALUATE FOR ARGUMENT .GE. 12.0,
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used      if (y <= xbig_gamma) then
! Subprogram not used         ysq = y*y
! Subprogram not used         sum = C(7)
! Subprogram not used         do i=1,6
! Subprogram not used            sum = sum/ysq + C(i)
! Subprogram not used         end do
! Subprogram not used         sum = sum/y - y + logsqrt2pi
! Subprogram not used         sum = sum + (y-0.5_r8)*log(y)
! Subprogram not used         res = exp(sum)
! Subprogram not used      else
! Subprogram not used         gamma = xinfr8
! Subprogram not used         return
! Subprogram not used      end if
! Subprogram not used   end if
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used !  FINAL ADJUSTMENTS AND RETURN
! Subprogram not used !----------------------------------------------------------------------
! Subprogram not used   if (negative_odd)  res = -res
! Subprogram not used   if (fact /= 1._r8) res = fact/res
! Subprogram not used   gamma = res
! Subprogram not used ! ---------- LAST LINE OF GAMMA ----------
! Subprogram not used end function shr_spfn_gamma_nonintrinsic_r8

!! Incomplete Gamma function
!!
!! @author  Tianyi Fan
!! @version August-2010
! Subprogram not used real(r8) elemental function shr_spfn_igamma(a, x)
! Subprogram not used   ! Upper incomplete gamma function.
! Subprogram not used   ! Modified for inclusion in this module and made
! Subprogram not used   ! pure elemental, September 2012
! Subprogram not used 
! Subprogram not used   real(r8), intent(in) ::      a
! Subprogram not used   real(r8), intent(in) ::      x
! Subprogram not used 
! Subprogram not used   ! local variable
! Subprogram not used   real(r8) :: xam, gin, s, r, t0
! Subprogram not used   integer  :: k
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   if (x == 0.0_r8) then
! Subprogram not used      shr_spfn_igamma = shr_spfn_gamma(a)
! Subprogram not used      return
! Subprogram not used   end if
! Subprogram not used 
! Subprogram not used   xam = -x + a * log(x)
! Subprogram not used   
! Subprogram not used   if ((xam > 700.0_r8) .or. (a > xbig_gamma)) then
! Subprogram not used      ! Out of bounds
! Subprogram not used      ! Return "huge" value.
! Subprogram not used      shr_spfn_igamma = xinfr8
! Subprogram not used      return
! Subprogram not used 
! Subprogram not used   else if (x <= (1.0_r8 + a)) then
! Subprogram not used      s = 1.0_r8 / a
! Subprogram not used      r = s
! Subprogram not used 
! Subprogram not used      do  k = 1,60
! Subprogram not used         r = r * x / (a+k)
! Subprogram not used         s = s + r
! Subprogram not used 
! Subprogram not used         if (abs(r/s) < 1.0e-15_r8) exit
! Subprogram not used      end do
! Subprogram not used         
! Subprogram not used      gin = exp(xam) * s           
! Subprogram not used      shr_spfn_igamma = shr_spfn_gamma(a) - gin
! Subprogram not used         
! Subprogram not used   else
! Subprogram not used      t0 = 0.0_r8
! Subprogram not used 
! Subprogram not used      do k = 60,1,-1
! Subprogram not used         t0 = (k - a) / (1.0_r8 + k / (x + t0))
! Subprogram not used      end do
! Subprogram not used 
! Subprogram not used      shr_spfn_igamma = exp(xam) / (x + t0)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used end function shr_spfn_igamma


end module shr_spfn_mod
