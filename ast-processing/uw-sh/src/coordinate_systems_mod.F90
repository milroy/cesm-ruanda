



module coordinate_systems_mod

! WARNING:  When using this class be sure that you know if the
! cubic coordinates are on the unit cube or the [-\pi/4,\pi/4] cube
! and if the spherical longitude is in [0,2\pi] or [-\pi,\pi]
  use kinds, only : real_kind, longdouble_kind
  implicit none
  private

  real(kind=real_kind), public, parameter :: DIST_THRESHOLD= 1.0D-9
  real(kind=real_kind), parameter :: one=1.0D0, two=2.0D0

  type, public :: cartesian2D_t
     real(real_kind) :: x             ! x coordinate
     real(real_kind) :: y             ! y coordinate
  end type cartesian2D_t

  type, public :: cartesian3D_t
     real(real_kind) :: x             ! x coordinate
     real(real_kind) :: y             ! y coordinate
     real(real_kind) :: z             ! z coordinate
  end type cartesian3D_t

  type, public :: spherical_polar_t
     real(real_kind) :: r             ! radius
     real(real_kind) :: lon           ! longitude
     real(real_kind) :: lat           ! latitude
  end type spherical_polar_t


  interface assignment ( = )
     module procedure copy_cart2d
  end interface

  interface operator( == )
     module procedure eq_cart2d
  end interface

  interface distance
     module procedure distance_cart2D
     module procedure distance_cart2D_v
     module procedure distance_cart3D
     module procedure distance_cart3D_v
  end interface

  interface change_coordinates
     module procedure spherical_to_cart_v
     module procedure spherical_to_cart
     module procedure cart_to_spherical_v
     module procedure cart_to_spherical
     module procedure aray_to_spherical
  end interface


  ! ==========================================
  ! Public Interfaces
  ! ==========================================

  public :: sphere_tri_area
  public :: surfareaxy
  public :: distance
  public :: change_coordinates
  public :: cart2cubedsphere
  public :: spherical_to_cart   !CE
  public :: projectpoint        ! should be called cubedsphere2spherical
  public :: cubedsphere2cart
  public :: sphere2cubedsphere
  public :: cube_face_number_from_cart
  public :: cube_face_number_from_sphere

! CE
  public :: cart2cubedspherexy

  public :: cart2spherical  !CE
  private :: copy_cart2d
  private :: eq_cart2d
  private :: distance_cart2D
  private :: distance_cart2D_v
  private :: distance_cart3D
  private :: distance_cart3D_v
  private :: spherical_to_cart_v
  !private :: spherical_to_cart
  private :: cart_to_spherical_v
  private :: cart_to_spherical
  private :: aray_to_spherical

contains

  ! ============================================
  ! copy_cart2d:
  !
  ! Overload assignment operator for cartesian2D_t
  ! ============================================

! Subprogram not used   subroutine copy_cart2d(cart2,cart1)
! Subprogram not used     implicit none
! Subprogram not used     type(cartesian2D_t), intent(out) :: cart2
! Subprogram not used     type(cartesian2D_t), intent(in)  :: cart1
! Subprogram not used     cart2%x=cart1%x
! Subprogram not used     cart2%y=cart1%y
! Subprogram not used   end subroutine copy_cart2d

  ! ============================================
  ! eq_cart2d:
  !
  ! Overload == operator for cartesian2D_t
  ! ============================================

! Subprogram not used   pure function eq_cart2d(cart2,cart1) result(is_same)
! Subprogram not used     implicit none
! Subprogram not used     type(cartesian2D_t), intent(in)  :: cart2
! Subprogram not used     type(cartesian2D_t), intent(in)  :: cart1
! Subprogram not used 
! Subprogram not used     logical :: is_same    
! Subprogram not used 
! Subprogram not used     if (distance(cart1,cart2)<DIST_THRESHOLD) then
! Subprogram not used        is_same=.true.
! Subprogram not used     else
! Subprogram not used        is_same=.false.
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used   end function eq_cart2d

  ! ===================================================
  ! distance_cart2D  : scalar version
  ! distance_cart2D_v: vector version
  !
  ! computes distance between cartesian 2D coordinates
  ! ===================================================

! Subprogram not used   pure function distance_cart2D(cart1,cart2) result(dist)
! Subprogram not used     implicit none
! Subprogram not used     type(cartesian2D_t), intent(in)           :: cart1
! Subprogram not used     type(cartesian2D_t), intent(in), optional :: cart2
! Subprogram not used     real(real_kind)   :: dist
! Subprogram not used 
! Subprogram not used     if (present(cart2)) then
! Subprogram not used        dist = SQRT((cart1%x-cart2%x)**2 + &
! Subprogram not used                    (cart1%y-cart2%y)**2   )
! Subprogram not used     else
! Subprogram not used        dist = SQRT(cart1%x*cart1%x + &
! Subprogram not used                    cart1%y*cart1%y   )
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used   end function distance_cart2D

! Subprogram not used   pure function distance_cart2D_v(cart1,cart2) result(dist)
! Subprogram not used     implicit none
! Subprogram not used     type(cartesian2D_t), intent(in)           :: cart1(:)
! Subprogram not used     type(cartesian2D_t), intent(in), optional :: cart2(:)
! Subprogram not used     real(real_kind)                           :: dist(SIZE(cart1))
! Subprogram not used 
! Subprogram not used     integer             :: i
! Subprogram not used 
! Subprogram not used     if (present(cart2)) then
! Subprogram not used        forall (i=1:SIZE(cart1)) dist(i) = distance_cart2D(cart1(i),cart2(i))
! Subprogram not used     else
! Subprogram not used        forall (i=1:SIZE(cart1)) dist(i) = distance_cart2D(cart1(i))
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used   end function distance_cart2D_v


  ! ===================================================
  ! distance_cart3D  : scalar version
  ! distance_cart3D_v: vector version
  ! ===================================================

! Subprogram not used   pure function distance_cart3D(cart1,cart2) result(dist)
! Subprogram not used     implicit none
! Subprogram not used     type(cartesian3D_t), intent(in)          :: cart1
! Subprogram not used     type(cartesian3D_t), intent(in),optional :: cart2
! Subprogram not used     real(real_kind)                          :: dist
! Subprogram not used 
! Subprogram not used     if (present(cart2)) then
! Subprogram not used        dist = SQRT((cart1%x-cart2%x)**2 + &
! Subprogram not used                    (cart1%y-cart2%y)**2 + &
! Subprogram not used                    (cart1%z-cart2%z)**2   )
! Subprogram not used     else
! Subprogram not used        dist = SQRT(cart1%x*cart1%x + &
! Subprogram not used                    cart1%y*cart1%y + &
! Subprogram not used                    cart1%z*cart1%z   )
! Subprogram not used     end if
! Subprogram not used   end function distance_cart3D

! Subprogram not used   pure function distance_cart3D_v(cart1,cart2) result(dist)
! Subprogram not used     implicit none
! Subprogram not used     type(cartesian3D_t), intent(in)          :: cart1(:)
! Subprogram not used     type(cartesian3D_t), intent(in),optional :: cart2(:)
! Subprogram not used     real(real_kind)                          :: dist(SIZE(cart1))
! Subprogram not used 
! Subprogram not used     integer             :: i
! Subprogram not used 
! Subprogram not used     if (present(cart2)) then
! Subprogram not used        forall (i=1:SIZE(cart1)) dist(i) = distance_cart3D(cart1(i),cart2(i))
! Subprogram not used     else
! Subprogram not used        forall (i=1:SIZE(cart1)) dist(i) = distance_cart3D(cart1(i))
! Subprogram not used     end if
! Subprogram not used   end function distance_cart3D_v

  ! ===================================================================
  ! spherical_to_cart:
  ! converts spherical polar {lon,lat}  to 3D cartesian {x,y,z}
  ! on unit sphere.  Note: spherical longitude is [0,2\pi]
  ! ===================================================================

  pure function spherical_to_cart(sphere) result (cart)
    implicit none
    type(spherical_polar_t), intent(in) :: sphere
    type(cartesian3D_t)                 :: cart

    cart%x=sphere%r*COS(sphere%lat)*COS(sphere%lon)
    cart%y=sphere%r*COS(sphere%lat)*SIN(sphere%lon)
    cart%z=sphere%r*SIN(sphere%lat)

  end function spherical_to_cart

  ! ===================================================================
  ! spherical_to_cart_v:
  ! converts spherical polar {lon,lat}  to 3D cartesian {x,y,z}
  ! on unit sphere.  Note: spherical longitude is [0,2\pi]
  ! ===================================================================

! Subprogram not used   pure function spherical_to_cart_v(sphere) result (cart)
! Subprogram not used     implicit none
! Subprogram not used     type(spherical_polar_t), intent(in) :: sphere(:)
! Subprogram not used     type(cartesian3D_t)                 :: cart(SIZE(sphere))
! Subprogram not used 
! Subprogram not used     integer                 :: i
! Subprogram not used 
! Subprogram not used     forall (i=1:SIZE(sphere)) cart(i) = spherical_to_cart(sphere(i))
! Subprogram not used   end function spherical_to_cart_v

  ! ==========================================================================
  ! cart_to_spherical:
  !
  ! converts 3D cartesian {x,y,z} to spherical polar {lon,lat} 
  ! on unit sphere. Note: spherical longitude is [0,2\pi]
  ! ==========================================================================

  ! scalar version

! Subprogram not used   pure function cart_to_spherical(cart) result (sphere)
! Subprogram not used     use physical_constants, only : dd_pi
! Subprogram not used     implicit none
! Subprogram not used     type(cartesian3D_t), intent(in) :: cart         
! Subprogram not used     type(spherical_polar_t)         :: sphere
! Subprogram not used 
! Subprogram not used     sphere%r=distance(cart)
! Subprogram not used     sphere%lat=ASIN(cart%z/sphere%r)
! Subprogram not used     sphere%lon=0
! Subprogram not used 
! Subprogram not used     ! ==========================================================
! Subprogram not used     ! enforce three facts:
! Subprogram not used     !
! Subprogram not used     ! 1) lon at poles is defined to be zero
! Subprogram not used     ! 
! Subprogram not used     ! 2) Grid points must be separated by about .01 Meter (on earth)
! Subprogram not used     !    from pole to be considered "not the pole".
! Subprogram not used     !
! Subprogram not used     ! 3) range of lon is { 0<= lon < 2*pi }
! Subprogram not used     !
! Subprogram not used     ! ==========================================================
! Subprogram not used 
! Subprogram not used !   if point is away from the POLE.  distance(cart) = distance from center of earth,
! Subprogram not used !   so this was a bug:
! Subprogram not used !    if (distance(cart) >= DIST_THRESHOLD) then 
! Subprogram not used 
! Subprogram not used     if ( abs(abs(sphere%lat)-DD_PI/2)  >= DIST_THRESHOLD ) then
! Subprogram not used        sphere%lon=ATAN2(cart%y,cart%x)
! Subprogram not used        if (sphere%lon<0) then
! Subprogram not used           sphere%lon=sphere%lon + 2*DD_PI
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used   end function cart_to_spherical

! Subprogram not used   pure function aray_to_spherical(coordinates) result (sphere)
! Subprogram not used     implicit none
! Subprogram not used     real(kind=real_kind),  intent(in)    :: coordinates(3)
! Subprogram not used     type(spherical_polar_t)              :: sphere
! Subprogram not used     type(cartesian3D_t)                  :: cart
! Subprogram not used     cart%x = coordinates(1)
! Subprogram not used     cart%y = coordinates(2)
! Subprogram not used     cart%z = coordinates(3)
! Subprogram not used     sphere = cart_to_spherical(cart) 
! Subprogram not used   end function aray_to_spherical


! Subprogram not used   pure function cart_to_spherical_v(cart) result (sphere)
! Subprogram not used     use physical_constants, only     : dd_pi
! Subprogram not used     implicit none
! Subprogram not used     type(cartesian3D_t), intent(in) :: cart(:)
! Subprogram not used     type(spherical_polar_t)         :: sphere(SIZE(cart))
! Subprogram not used 
! Subprogram not used     integer                 :: i
! Subprogram not used     forall (i=1:SIZE(cart)) sphere(i) = cart_to_spherical(cart(i))
! Subprogram not used   end function cart_to_spherical_v




  function unit_face_based_cube_to_unit_sphere(cart, face_no) result(sphere)

! Note: Output spherical longitude is [-pi,pi]

! Project from a UNIT cube to a UNIT sphere.  ie, the lenght of the cube edge is 2.
! Face 1 of the cube touches the sphere at longitude, latitude (0,0). The negative 
! x axis is negative longitude (ie. going west is negative), the positive x axis
! is increasing longitude.  Face 1 maps the Face 1 to the lat,lon on the sphere:
!    [-1,1] x [-1,1] => [-\pi/4,\pi/4] x [-\pi/4, \pi/4]

! Face 2 continues with increasing longitude (ie to the east of Face 1).  
! The left edge of Face 2 (negative x) is the right edge of Face 1 (positive x)
! The latitude is the same as Face 1, but the longitude increases:
!    [-1,1] x [-1,1] => [\pi/4, 3\pi/4] x [-\pi/4, \pi/4]

! Face 3 continues with increasing longitude (ie to the east of Face 2).  
! Face 3 is like Face 1, but the x coordinates are reversed, ie. decreasing x 
! is increasing longitude:
!    [-1,1] x [-1,1]  =    [-1,0] x [-1,1] U  [0,1] x [-1,1] =>
!            [3\pi/4,\pi] x [-\pi, -3\pi/4]

! Face 4 finally connects Face 3 to Face 1.  Like Face 2, but wtih opposite x
!    [-1,1] x [-1,1] => [-3\pi/4, -\pi/4] x [-\pi/4, \pi/4]

! Face 5 is along the bottom edges of Faces 1,2,3,and 4 so the latitude goes from
! -\pi/4 to -\pi/2.  The tricky part is lining up the longitude.  The zero longitude
! must line up with the center of Face 1. ATAN2(x,1) = 0 => x = 0.  
! So the (0,1) point on Face 5 is the zero longitude on the sphere.  The top edge of 
! Face 5 is the bottom edge of Face 1. 
! ATAN(x,0) = \pi/2 => x = 1, so the right edge of Face 5 is the bottom of Face 2.
! Continueing, the bottom edge of 5 is the bottom of 3.  Left of 5 is bottom of 4.

! Face 6 is along the top edges of Faces 1,2,3 and 4 so the latitude goes from  
! \pi/4 to \pi/2.   The zero longitude must line up with the center of Face 1.  
! This is just like Face 5, but the y axis is reversed.  So the bottom edge of Face 6
! is the top edge of Face 1.  The right edge of Face 6 is the top of Face 2.  The
! top of 6 the top of 3 and the left of 6 the top of 4.

    use physical_constants, only : DD_PI
    use parallel_mod,       only : abortmp
    implicit none
    type (cartesian2d_t), intent(in)     :: cart   ! On face_no of a unit cube
    integer,              intent(in)     :: face_no 
 
    type (spherical_polar_t)             :: sphere

    integer i,j
    real(kind=real_kind) :: r!, l_inf

! MNL: removing check that points are on the unit cube because we allow
! spherical grids to map beyond the extent of the cube (though we probably
! should still have an upper bound for how far past the edge the element lies)
!    l_inf = MAX(ABS(cart%x), ABS(cart%y)) 
!    if (1.01 < l_inf) then
!      call abortmp('unit_face_based_cube_to_unit_sphere: Input not on unit cube.')
!    end if

    sphere%r=one
    r = SQRT( one + (cart%x)**2 + (cart%y)**2)
    select case (face_no)
    case (1) 
       sphere%lat=ASIN((cart%y)/r)
       sphere%lon=ATAN2(cart%x,one)
    case (2) 
       sphere%lat=ASIN((cart%y)/r)
       sphere%lon=ATAN2(one,-cart%x)
    case (3) 
       sphere%lat=ASIN((cart%y)/r)
       sphere%lon=ATAN2(-cart%x,-one)
    case (4) 
       sphere%lat=ASIN((cart%y)/r)
       sphere%lon=ATAN2(-one,cart%x)
    case (5) 
       if (ABS(cart%y) > DIST_THRESHOLD .or. ABS(cart%x) > DIST_THRESHOLD ) then
          sphere%lon=ATAN2(cart%x, cart%y )
       else
          sphere%lon= 0.0D0     ! longitude is meaningless at south pole set to 0.0
       end if
       sphere%lat=ASIN(-one/r)
    case (6) 
       if (ABS(cart%y) > DIST_THRESHOLD .or. ABS(cart%x) > DIST_THRESHOLD ) then
          sphere%lon = ATAN2(cart%x, -cart%y)
       else
          sphere%lon= 0.0D0     ! longitude is meaningless at north pole set to 0.0
       end if
       sphere%lat=ASIN(one/r)
    case default
       call abortmp('unit_face_based_cube_to_unit_sphere: Face number not 1 to 6.')
    end select

    if (sphere%lon < 0.0D0) then
       sphere%lon=sphere%lon + two*DD_PI
    end if

  end function unit_face_based_cube_to_unit_sphere

! Subprogram not used   function cart2spherical(x,y, face_no) result(sphere)
! Subprogram not used ! IMPORTANT: INPUT ARE the REAL cartesian from the cube sphere
! Subprogram not used ! Note: Output spherical longitude is [-pi,pi]
! Subprogram not used 
! Subprogram not used ! Project from a UNIT cube to a UNIT sphere.  ie, the lenght of the cube edge is 2.
! Subprogram not used ! Face 1 of the cube touches the sphere at longitude, latitude (0,0). The negative 
! Subprogram not used ! x axis is negative longitude (ie. going west is negative), the positive x axis
! Subprogram not used ! is increasing longitude.  Face 1 maps the Face 1 to the lat,lon on the sphere:
! Subprogram not used !    [-1,1] x [-1,1] => [-\pi/4,\pi/4] x [-\pi/4, \pi/4]
! Subprogram not used 
! Subprogram not used ! Face 2 continues with increasing longitude (ie to the east of Face 1).  
! Subprogram not used ! The left edge of Face 2 (negative x) is the right edge of Face 1 (positive x)
! Subprogram not used ! The latitude is the same as Face 1, but the longitude increases:
! Subprogram not used !    [-1,1] x [-1,1] => [\pi/4, 3\pi/4] x [-\pi/4, \pi/4]
! Subprogram not used 
! Subprogram not used ! Face 3 continues with increasing longitude (ie to the east of Face 2).  
! Subprogram not used ! Face 3 is like Face 1, but the x coordinates are reversed, ie. decreasing x 
! Subprogram not used ! is increasing longitude:
! Subprogram not used !    [-1,1] x [-1,1]  =    [-1,0] x [-1,1] U  [0,1] x [-1,1] =>
! Subprogram not used !            [3\pi/4,\pi] x [-\pi, -3\pi/4]
! Subprogram not used 
! Subprogram not used ! Face 4 finally connects Face 3 to Face 1.  Like Face 2, but wtih opposite x
! Subprogram not used !    [-1,1] x [-1,1] => [-3\pi/4, -\pi/4] x [-\pi/4, \pi/4]
! Subprogram not used 
! Subprogram not used ! Face 5 is along the bottom edges of Faces 1,2,3,and 4 so the latitude goes from
! Subprogram not used ! -\pi/4 to -\pi/2.  The tricky part is lining up the longitude.  The zero longitude
! Subprogram not used ! must line up with the center of Face 1. ATAN2(x,1) = 0 => x = 0.  
! Subprogram not used ! So the (0,1) point on Face 5 is the zero longitude on the sphere.  The top edge of 
! Subprogram not used ! Face 5 is the bottom edge of Face 1. 
! Subprogram not used ! ATAN(x,0) = \pi/2 => x = 1, so the right edge of Face 5 is the bottom of Face 2.
! Subprogram not used ! Continueing, the bottom edge of 5 is the bottom of 3.  Left of 5 is bottom of 4.
! Subprogram not used 
! Subprogram not used ! Face 6 is along the top edges of Faces 1,2,3 and 4 so the latitude goes from  
! Subprogram not used ! \pi/4 to \pi/2.   The zero longitude must line up with the center of Face 1.  
! Subprogram not used ! This is just like Face 5, but the y axis is reversed.  So the bottom edge of Face 6
! Subprogram not used ! is the top edge of Face 1.  The right edge of Face 6 is the top of Face 2.  The
! Subprogram not used ! top of 6 the top of 3 and the left of 6 the top of 4.
! Subprogram not used 
! Subprogram not used     use physical_constants, only : DD_PI
! Subprogram not used     use parallel_mod,       only : abortmp
! Subprogram not used     implicit none
! Subprogram not used     real(kind=real_kind), intent(in)     :: x,y   ! On face_no of a unit cube
! Subprogram not used     integer,              intent(in)     :: face_no 
! Subprogram not used  
! Subprogram not used     type (spherical_polar_t)             :: sphere
! Subprogram not used 
! Subprogram not used     integer i,j
! Subprogram not used     real(kind=real_kind) :: r!, l_inf
! Subprogram not used 
! Subprogram not used ! MNL: removing check that points are on the unit cube because we allow
! Subprogram not used ! spherical grids to map beyond the extent of the cube (though we probably
! Subprogram not used ! should still have an upper bound for how far past the edge the element lies)
! Subprogram not used !    l_inf = MAX(ABS(cart%x), ABS(cart%y)) 
! Subprogram not used !    if (1.01 < l_inf) then
! Subprogram not used !      call abortmp('unit_face_based_cube_to_unit_sphere: Input not on unit cube.')
! Subprogram not used !    end if
! Subprogram not used 
! Subprogram not used     sphere%r=one
! Subprogram not used     r = SQRT( one + x**2 + y**2)
! Subprogram not used     select case (face_no)
! Subprogram not used     case (1) 
! Subprogram not used        sphere%lat=ASIN(y/r)
! Subprogram not used        sphere%lon=ATAN2(x,one)
! Subprogram not used     case (2) 
! Subprogram not used        sphere%lat=ASIN(y/r)
! Subprogram not used        sphere%lon=ATAN2(one,-x)
! Subprogram not used     case (3) 
! Subprogram not used        sphere%lat=ASIN(y/r)
! Subprogram not used        sphere%lon=ATAN2(-x,-one)
! Subprogram not used     case (4) 
! Subprogram not used        sphere%lat=ASIN(y/r)
! Subprogram not used        sphere%lon=ATAN2(-one,x)
! Subprogram not used     case (5) 
! Subprogram not used        if (ABS(y) > DIST_THRESHOLD .or. ABS(x) > DIST_THRESHOLD ) then
! Subprogram not used           sphere%lon=ATAN2(x, y )
! Subprogram not used        else
! Subprogram not used           sphere%lon= 0.0D0     ! longitude is meaningless at south pole set to 0.0
! Subprogram not used        end if
! Subprogram not used        sphere%lat=ASIN(-one/r)
! Subprogram not used     case (6) 
! Subprogram not used        if (ABS(y) > DIST_THRESHOLD .or. ABS(x) > DIST_THRESHOLD ) then
! Subprogram not used           sphere%lon = ATAN2(x, -y)
! Subprogram not used        else
! Subprogram not used           sphere%lon= 0.0D0     ! longitude is meaningless at north pole set to 0.0
! Subprogram not used        end if
! Subprogram not used        sphere%lat=ASIN(one/r)
! Subprogram not used     case default
! Subprogram not used        call abortmp('unit_face_based_cube_to_unit_sphere: Face number not 1 to 6.')
! Subprogram not used     end select
! Subprogram not used 
! Subprogram not used     if (sphere%lon < 0.0D0) then
! Subprogram not used        sphere%lon=sphere%lon + two*DD_PI
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used   end function cart2spherical








! Note: Output spherical longitude is [-pi,pi]
  function projectpoint(cartin, face_no) result(sphere)         

! Projection from a [-pi/4, \pi/4] sized cube.  
! This will be checked because unit_face_based_cube_to_unit_sphere checks the ranges.
! See unit_face_based_cube_to_unit_sphere for documentation.

    implicit none
    type (cartesian2d_t), intent(in)     :: cartin   
    integer,              intent(in)     :: face_no
    type (spherical_polar_t)             :: sphere
    type (cartesian2d_t)                 :: cart   

    !ASC  This is X and Y and not xhi eta ...

    cart%x = TAN(cartin%x)
    cart%y = TAN(cartin%y)

    sphere = unit_face_based_cube_to_unit_sphere(cart, face_no)

  end function projectpoint

  ! takes a 2D point on a face of the cube of size [-\pi/4, \pi/4] and projects it 
  ! onto a 3D point on a cube of size [-1,1] in R^3
  function cubedsphere2cart(cartin, face_no) result(cart)
    implicit none
    type (cartesian2d_t), intent(in)    :: cartin   ! assumed to be cartesian coordinates of cube
    integer,              intent(in)    :: face_no

    type(cartesian3D_t)                 :: cart

    cart = spherical_to_cart(projectpoint(cartin, face_no))

  end function cubedsphere2cart


  ! onto a cube of size [-\pi/2,\pi/2] in R^3
  ! the spherical longitude can be either in [0,2\pi] or [-\pi,\pi]
! Subprogram not used   pure function sphere2cubedsphere (sphere, face_no) result(cart)
! Subprogram not used     use physical_constants, only : dd_pi
! Subprogram not used     implicit none
! Subprogram not used     type(spherical_polar_t), intent(in) :: sphere
! Subprogram not used     integer,                 intent(in) :: face_no
! Subprogram not used 
! Subprogram not used     type(cartesian2d_t)                 :: cart
! Subprogram not used     real(kind=real_kind)               :: xp,yp
! Subprogram not used     real(kind=real_kind)               :: lat,lon
! Subprogram not used     real(kind=real_kind)               :: pi,twopi, pi2, pi3, pi4
! Subprogram not used 
! Subprogram not used     lat = sphere%lat
! Subprogram not used     lon = sphere%lon
! Subprogram not used 
! Subprogram not used     pi    = DD_PI 
! Subprogram not used     twopi = 2.0D0 * pi
! Subprogram not used     pi2   = pi * 0.5D0 
! Subprogram not used     pi3   = pi * 1.5D0               
! Subprogram not used     pi4   = pi * 0.25D0
! Subprogram not used 
! Subprogram not used     select case (face_no)
! Subprogram not used     case  (1) 
! Subprogram not used        xp = lon
! Subprogram not used        if (pi < lon) xp = lon - twopi !if lon in [0,2\pi]
! Subprogram not used        yp = atan(tan(lat)/cos(xp))
! Subprogram not used     case  (2) 
! Subprogram not used        xp = lon - pi2
! Subprogram not used        yp = atan(tan(lat)/cos(xp))
! Subprogram not used     case  (3) 
! Subprogram not used        xp = lon - pi
! Subprogram not used        if (lon < 0) xp = lon + pi  !if lon in [0,2\pi]
! Subprogram not used        yp = atan(tan(lat)/cos(xp))
! Subprogram not used     case  (4) 
! Subprogram not used        xp = lon - pi3
! Subprogram not used        if (lon < 0) xp = lon + pi2  !if lon in [0,2\pi]
! Subprogram not used        yp = atan(tan(lat)/cos(xp))
! Subprogram not used     case  (5) 
! Subprogram not used        xp = atan(-sin(lon)/tan(lat))
! Subprogram not used        yp = atan(-cos(lon)/tan(lat))
! Subprogram not used     case  (6) 
! Subprogram not used        xp = atan( sin(lon)/tan(lat))
! Subprogram not used        yp = atan(-cos(lon)/tan(lat))
! Subprogram not used     end select
! Subprogram not used 
! Subprogram not used     ! coordinates on the cube:
! Subprogram not used     cart%x = xp
! Subprogram not used     cart%y = yp
! Subprogram not used     
! Subprogram not used   end function sphere2cubedsphere 

! Go from an arbitrary sized cube in 3D 
! to a [-\pi/4,\pi/4] sized cube with (face,2d) coordinates.  
!
!                        Z
!                        |
!                        |
!                        |
!                        |
!                        ---------------Y
!                       /
!                      /
!                     /
!                    /
!                   X
!
! NOTE: Face 1 =>  X positive constant face of cube
!       Face 2 =>  Y positive constant face of cube
!       Face 3 =>  X negative constant face of cube
!       Face 4 =>  Y negative constant face of cube
!       Face 5 =>  Z negative constant face of cube
!       Face 6 =>  Z positive constant face of cube
! Subprogram not used   pure function cart2cubedsphere(cart3D, face_no) result(cart)
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used     type(cartesian3D_t),intent(in) :: cart3d
! Subprogram not used     integer,            intent(in) :: face_no
! Subprogram not used     type (cartesian2d_t)           :: cart   
! Subprogram not used 
! Subprogram not used     real(kind=real_kind) :: x,y
! Subprogram not used 
! Subprogram not used     select case (face_no) 
! Subprogram not used     case (1)
! Subprogram not used        x =  cart3D%y/cart3D%x
! Subprogram not used        y =  cart3D%z/cart3D%x
! Subprogram not used     case (2)
! Subprogram not used        x = -cart3D%x/cart3D%y
! Subprogram not used        y =  cart3D%z/cart3D%y
! Subprogram not used     case (3)
! Subprogram not used        x =  cart3D%y/cart3D%x
! Subprogram not used        y = -cart3D%z/cart3D%x
! Subprogram not used     case (4)
! Subprogram not used        x = -cart3D%x/cart3D%y
! Subprogram not used        y = -cart3D%z/cart3D%y
! Subprogram not used     case (5)
! Subprogram not used        x  = -cart3D%y/cart3D%z
! Subprogram not used        y  = -cart3D%x/cart3D%z
! Subprogram not used     case (6)
! Subprogram not used        x  =  cart3D%y/cart3D%z
! Subprogram not used        y  = -cart3D%x/cart3D%z
! Subprogram not used     end select
! Subprogram not used     cart%x = ATAN(x)
! Subprogram not used     cart%y = ATAN(y)
! Subprogram not used   end function cart2cubedsphere



! This function divides three dimentional space up into 
! six sectors.  These sectors are then considered as the
! faces of the cube.  It should work for any (x,y,z) coordinate
! if on a sphere or on a cube.
! Subprogram not used   pure function cube_face_number_from_cart(cart) result(face_no)
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used     type(cartesian3D_t),intent(in) :: cart  
! Subprogram not used     integer :: face_no
! Subprogram not used 
! Subprogram not used     real(real_kind) :: x,y,z
! Subprogram not used     x=cart%x
! Subprogram not used     y=cart%y
! Subprogram not used     z=cart%z
! Subprogram not used 
! Subprogram not used ! Divide the X-Y plane into for quadrants of 
! Subprogram not used ! [-\pi/2,\pi/2], [\pi/2,3\pi/2], .....
! Subprogram not used ! based on the lines X=Y and X=-Y.  This divides
! Subprogram not used ! 3D space up into four sections.  Doing the same
! Subprogram not used ! for the XZ and YZ planes divides space into six
! Subprogram not used ! sections.  Can also be thought of as conic sections
! Subprogram not used ! in the L_infinity norm.  
! Subprogram not used 
! Subprogram not used     if (y<x .and. y>-x) then      ! x>0, Face 1,5 or 6
! Subprogram not used       if (z>x) then
! Subprogram not used          face_no=6  ! north pole
! Subprogram not used       else if (z<-x) then
! Subprogram not used          face_no=5 ! south pole
! Subprogram not used       else
! Subprogram not used          face_no = 1
! Subprogram not used       endif
! Subprogram not used    else if (y>x .and. y<-x) then  ! x<0
! Subprogram not used       if (z>-x) then
! Subprogram not used          face_no=6 ! north pole
! Subprogram not used       else if (z<x) then
! Subprogram not used          face_no=5 ! south pole
! Subprogram not used       else 
! Subprogram not used          face_no=3
! Subprogram not used       endif
! Subprogram not used    else if (y>x .and. y>-x) then  ! y>0
! Subprogram not used       if (z>y) then
! Subprogram not used          face_no=6 ! north pole
! Subprogram not used       else if (z<-y) then
! Subprogram not used         face_no = 5 ! south pole
! Subprogram not used       else 
! Subprogram not used          face_no=2
! Subprogram not used       endif
! Subprogram not used    else if (y<x .and. y<-x) then  ! y<0
! Subprogram not used       if (z>-y) then
! Subprogram not used          face_no=6 ! north pole
! Subprogram not used       else if (z<y) then
! Subprogram not used          face_no=5 ! south pole
! Subprogram not used       else 
! Subprogram not used          face_no=4
! Subprogram not used       endif
! Subprogram not used    else
! Subprogram not used       ! abs(y) = abs(x).  point is on cube edge, or on face 5 or 6:
! Subprogram not used       if (abs(x)<z) then
! Subprogram not used          face_no=6  ! north pole
! Subprogram not used       else if (z<-abs(x)) then
! Subprogram not used          face_no=5 ! south pole
! Subprogram not used       else if (0 < x .and. 0 < y) then
! Subprogram not used          face_no = 1
! Subprogram not used       else if (x < 0 .and. 0 < y) then
! Subprogram not used          face_no = 2
! Subprogram not used       else if (x < 0 .and. y < 0) then
! Subprogram not used          face_no = 3
! Subprogram not used       else
! Subprogram not used          face_no = 4
! Subprogram not used       endif
! Subprogram not used    endif
! Subprogram not used    
! Subprogram not used    end function cube_face_number_from_cart

! This could be done directly by using the lon, lat coordinates,
! but call cube_face_number_from_cart just so that there is one place
! to do the conversions and they are all consistant.
! Subprogram not used   pure function cube_face_number_from_sphere(sphere) result(face_no)
! Subprogram not used     implicit none
! Subprogram not used     type (spherical_polar_t), intent(in) :: sphere
! Subprogram not used     type (cartesian3d_t)                 :: cart   
! Subprogram not used     integer                              :: face_no 
! Subprogram not used 
! Subprogram not used     cart =  spherical_to_cart(sphere)
! Subprogram not used     face_no = cube_face_number_from_cart(cart) 
! Subprogram not used   end function cube_face_number_from_sphere

! CE, need real (cartesian) xy coordinates on the cubed sphere
! Subprogram not used subroutine cart2cubedspherexy(cart3d,face_no,cartxy)
! Subprogram not used 
! Subprogram not used   type(cartesian3D_t),intent(in)      :: cart3d
! Subprogram not used   integer, intent(in)                 :: face_no
! Subprogram not used   type (cartesian2d_t), intent(out)   :: cartxy   
! Subprogram not used 
! Subprogram not used   ! a (half length of a cube side) is supposed to be 1
! Subprogram not used   select case (face_no)
! Subprogram not used   case (1)
! Subprogram not used      cartxy%x = cart3D%y/cart3D%x
! Subprogram not used      cartxy%y = cart3D%z/cart3D%x
! Subprogram not used   case (2)
! Subprogram not used      cartxy%x = -cart3D%x/cart3D%y
! Subprogram not used      cartxy%y = cart3D%z/cart3D%y
! Subprogram not used   case (3)
! Subprogram not used      cartxy%x = cart3D%y/cart3D%x
! Subprogram not used      cartxy%y = -cart3D%z/cart3D%x
! Subprogram not used   case (4)
! Subprogram not used      cartxy%x = -cart3D%x/cart3D%y
! Subprogram not used      cartxy%y = -cart3D%z/cart3D%y
! Subprogram not used   case (5)       !bottom face
! Subprogram not used      cartxy%x  = -cart3D%y/cart3D%z
! Subprogram not used      cartxy%y  = -cart3D%x/cart3D%z
! Subprogram not used   case (6)        !top face
! Subprogram not used      cartxy%x  = cart3D%y/cart3D%z
! Subprogram not used      cartxy%y  = -cart3D%x/cart3D%z
! Subprogram not used   end select
! Subprogram not used end subroutine cart2cubedspherexy
! CE END




! Subprogram not used   subroutine sphere_tri_area( v1, v2, v3, area )
! Subprogram not used   !  input: v1(3),v2(3),v3(3)  cartesian coordinates of triangle
! Subprogram not used   !  output: area
! Subprogram not used   !  based on formulas in STRI_QUAD:
! Subprogram not used   !  http://people.sc.fsu.edu/~burkardt/f_src/stri_quad/stri_quad.html
! Subprogram not used   use physical_constants, only : dd_pi
! Subprogram not used   implicit none
! Subprogram not used   real(kind=real_kind) area
! Subprogram not used   real(kind=real_kind) a,b,c,al,bl,cl,sina,sinb,sinc,sins,a1,b1,c1
! Subprogram not used   type (cartesian3D_t) v1,v2,v3
! Subprogram not used   
! Subprogram not used   ! compute great circle lengths
! Subprogram not used   al = acos( v2%x * v3%x + v2%y * v3%y + v2%z * v3%z )
! Subprogram not used   bl = acos( v3%x * v1%x + v3%y * v1%y + v3%z * v1%z )
! Subprogram not used   cl = acos( v1%x * v2%x + v1%y * v2%y + v1%z * v2%z )
! Subprogram not used 
! Subprogram not used   ! compute angles
! Subprogram not used   sina = sin( (bl+cl-al)/2 )  ! sin(sl-al)
! Subprogram not used   sinb = sin( (al+cl-bl)/2 )  ! sin(sl-bl)
! Subprogram not used   sinc = sin( (al+bl-cl)/2 ) 
! Subprogram not used   sins = sin( (al+bl+cl)/2 )
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   ! for small areas, formula above looses precision.  
! Subprogram not used   ! 2atan(x) + 2atan(1/x) = pi      
! Subprogram not used   ! 2atan(x) - pi = -2atan(1/x)
! Subprogram not used   a = sqrt( (sinb*sinc) / (sins*sina) ) 
! Subprogram not used   b = sqrt( (sina*sinc) / (sins*sinb) ) 
! Subprogram not used   c = sqrt( (sina*sinb) / (sins*sinc) ) 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   a1 = 2*atan(a)
! Subprogram not used   b1 = 2*atan(b)
! Subprogram not used   c1 = 2*atan(c)
! Subprogram not used 
! Subprogram not used   if (a.gt.b.and.a.gt.c) then
! Subprogram not used      a1 = -2*atan(1/a)
! Subprogram not used   else if (b.gt.c) then
! Subprogram not used      b1 = -2*atan(1/b)
! Subprogram not used   else 
! Subprogram not used      c1 = -2*atan(1/c)
! Subprogram not used   endif
! Subprogram not used   ! apply Girard's theorem
! Subprogram not used   area = a1+b1+c1  
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   end subroutine sphere_tri_area

!CE, 5.May 2011
!INPUT: Points in xy cubed sphere coordinates, counterclockwise
!OUTPUT: corresponding area on the sphere
! Subprogram not used function surfareaxy(x1,x2,y1,y2) result(area)
! Subprogram not used   implicit none
! Subprogram not used   real (kind=real_kind), intent(in)    :: x1, x2, y1, y2
! Subprogram not used   real (kind=real_kind)   :: area
! Subprogram not used   real (kind=real_kind)   :: a1,a2,a3,a4
! Subprogram not used 
! Subprogram not used   ! cubed-sphere cell area, from Lauritzen & Nair MWR 2008
! Subprogram not used   ! central angles:
! Subprogram not used   ! cube face: -pi/4,-pi/4 -> pi/4,pi/4
! Subprogram not used   ! this formula gives 2   so normalize by 4pi/6 / 2 = pi/3
! Subprogram not used   ! use implementation where the nodes a counterclockwise (not as in the paper)
! Subprogram not used   a1 = acos(-sin(atan(x1))*sin(atan(y1)))             
! Subprogram not used   a2 =-acos(-sin(atan(x2))*sin(atan(y1)))  
! Subprogram not used   a3 = acos(-sin(atan(x2))*sin(atan(y2)))              
! Subprogram not used   a4 =-acos(-sin(atan(x1))*sin(atan(y2)))         
! Subprogram not used   area = (a1+a2+a3+a4)
! Subprogram not used   return
! Subprogram not used end function surfareaxy



end module coordinate_systems_mod
