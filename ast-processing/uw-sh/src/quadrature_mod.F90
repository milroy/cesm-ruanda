





module quadrature_mod
  use kinds, only : longdouble_kind
  implicit none
  private

  type, public :: quadrature_t
     real (kind=longdouble_kind), dimension(:), pointer :: points
     real (kind=longdouble_kind), dimension(:), pointer :: weights
  end type quadrature_t

  public  :: gausslobatto
  public  :: test_gausslobatto
  public  :: gauss
  public  :: test_gauss
  public  :: legendre
  public  :: jacobi
  public  :: quad_norm

  public  :: trapezoid
  private :: trapN
  public  :: simpsons
  public  :: gaussian_int

  private :: gausslobatto_pts
  private :: gausslobatto_wts
  private :: gauss_pts
  private :: gauss_wts
  private :: jacobi_polynomials
  private :: jacobi_derivatives


contains

  ! ==============================================================
  ! gauss:
  !
  ! Find the Gauss collocation points and the corresponding weights.
  !
  ! ==============================================================

! Subprogram not used   function gauss(npts) result(gs)
! Subprogram not used     integer, intent(in) :: npts
! Subprogram not used     type (quadrature_t) :: gs
! Subprogram not used 
! Subprogram not used     allocate(gs%points(npts))
! Subprogram not used     allocate(gs%weights(npts))
! Subprogram not used 
! Subprogram not used     gs%points=gauss_pts(npts)
! Subprogram not used     gs%weights=gauss_wts(npts,gs%points)
! Subprogram not used 
! Subprogram not used   end function gauss


  ! ==============================================================
  ! gauss_pts:
  !
  ! Compute the Gauss Collocation points
  ! for Jacobi Polynomials
  !
  ! ==============================================================

! Subprogram not used   function gauss_pts(np1) result(pts)
! Subprogram not used     use physical_constants, only : qq_pi
! Subprogram not used 
! Subprogram not used     integer, intent(in)     :: np1        ! Number of velocity grid points
! Subprogram not used     real (kind=longdouble_kind) :: pts(np1)
! Subprogram not used 
! Subprogram not used     ! Local variables
! Subprogram not used 
! Subprogram not used     real (kind=longdouble_kind) :: alpha,beta
! Subprogram not used     real (kind=longdouble_kind) :: xjac(0:np1-1)
! Subprogram not used     real (kind=longdouble_kind) :: jac(0:np1)
! Subprogram not used     real (kind=longdouble_kind) :: djac(0:np1)
! Subprogram not used 
! Subprogram not used     integer  prec                    ! number of mantissa bits
! Subprogram not used     real (kind=longdouble_kind) eps      ! machine epsilon
! Subprogram not used     real (kind=longdouble_kind), parameter :: convthresh = 10  ! convergence threshold relative
! Subprogram not used     ! to machine epsilon
! Subprogram not used     integer, parameter :: kstop = 30 ! max iterations for polynomial deflation
! Subprogram not used 
! Subprogram not used     real (kind=longdouble_kind) :: poly
! Subprogram not used     real (kind=longdouble_kind) :: pder
! Subprogram not used     real (kind=longdouble_kind) :: recsum,thresh
! Subprogram not used     real (kind=longdouble_kind) :: dth
! Subprogram not used 
! Subprogram not used     real (kind=longdouble_kind) :: x
! Subprogram not used     real (kind=longdouble_kind) :: delx
! Subprogram not used     real (kind=longdouble_kind) :: c0,c1,c2,c10
! Subprogram not used 
! Subprogram not used     integer i,j,k
! Subprogram not used     integer n, nh
! Subprogram not used 
! Subprogram not used     n  = np1 - 1
! Subprogram not used     c0 = 0.0_longdouble_kind
! Subprogram not used     c1 = 1.0_longdouble_kind
! Subprogram not used     c2 = 2.0_longdouble_kind
! Subprogram not used     c10 = 10.0_longdouble_kind
! Subprogram not used     alpha = c0
! Subprogram not used     beta  = c0
! Subprogram not used 
! Subprogram not used     ! =========================================================
! Subprogram not used     ! compute machine precision and set the convergence
! Subprogram not used     ! threshold thresh to 10 times that level
! Subprogram not used     ! =========================================================
! Subprogram not used 
! Subprogram not used     prec   = precision(c10)
! Subprogram not used     eps    = c10**(-prec)
! Subprogram not used     thresh = convthresh*eps
! Subprogram not used 
! Subprogram not used     ! ============================================================
! Subprogram not used     ! Compute first half of the roots by "polynomial deflation".
! Subprogram not used     ! ============================================================
! Subprogram not used 
! Subprogram not used     dth = QQ_PI/(2*n+2)
! Subprogram not used 
! Subprogram not used     nh  = (n+1)/2
! Subprogram not used 
! Subprogram not used     do j=0,nh-1
! Subprogram not used        x=COS((c2*j+1)*dth)          ! first guess at root
! Subprogram not used        k=0
! Subprogram not used        delx=c1
! Subprogram not used        do while(k<kstop .and. ABS(delx) > thresh)
! Subprogram not used           call jacobi(n+1,x,alpha,beta,jac(0:n+1),djac(0:n+1))
! Subprogram not used           poly =  jac(n+1)
! Subprogram not used           pder = djac(n+1)
! Subprogram not used           recsum=c0
! Subprogram not used           do i=0,j-1
! Subprogram not used              recsum = recsum + c1/(x-xjac(i))
! Subprogram not used           end do
! Subprogram not used           delx = -poly/(pder-recsum*poly)
! Subprogram not used           x = x + delx
! Subprogram not used           k = k + 1
! Subprogram not used        end do
! Subprogram not used 
! Subprogram not used        xjac(j)=x
! Subprogram not used 
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     ! ================================================
! Subprogram not used     ! compute the second half of the roots by symmetry
! Subprogram not used     ! ================================================
! Subprogram not used 
! Subprogram not used     do j=0,nh
! Subprogram not used        xjac(n-j) = -xjac(j)
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     if (MODULO(n,2)==0) xjac(nh)=c0
! Subprogram not used 
! Subprogram not used     ! ====================================================
! Subprogram not used     ! Reverse the sign of everything so that indexing
! Subprogram not used     ! increases with position
! Subprogram not used     ! ====================================================
! Subprogram not used 
! Subprogram not used     do j=0,n
! Subprogram not used        pts(j+1) = -xjac(j)
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used   end function gauss_pts

  ! ================================================
  ! gauss_wts:
  !
  ! Gauss Legendre Weights
  ! ================================================

! Subprogram not used   function gauss_wts(np1, gpts) result(wts)
! Subprogram not used 
! Subprogram not used     integer, intent(in)                 :: np1
! Subprogram not used     real (kind=longdouble_kind), intent(in) :: gpts(np1)  ! Gauss-Legendre points
! Subprogram not used     real (kind=longdouble_kind)             :: wts(np1)   ! Gauss-Legendre weights
! Subprogram not used 
! Subprogram not used     ! Local variables
! Subprogram not used 
! Subprogram not used     real (kind=longdouble_kind) :: c0,c1,c2
! Subprogram not used     real (kind=longdouble_kind) :: alpha
! Subprogram not used     real (kind=longdouble_kind) :: beta
! Subprogram not used     real (kind=longdouble_kind) :: djac(np1)
! Subprogram not used     integer i,n
! Subprogram not used 
! Subprogram not used     c0    = 0.0_longdouble_kind
! Subprogram not used     c1    = 1.0_longdouble_kind
! Subprogram not used     c2    = 2.0_longdouble_kind
! Subprogram not used 
! Subprogram not used     alpha = c0
! Subprogram not used     beta  = c0
! Subprogram not used     n     = np1-1
! Subprogram not used 
! Subprogram not used     djac=jacobi_derivatives(np1,alpha,beta,np1,gpts)     
! Subprogram not used 
! Subprogram not used     do i=1,np1
! Subprogram not used        wts(i)=c2/((c1-gpts(i)**2)*djac(i)*djac(i))
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used   end function gauss_wts


  ! ==============================================================
  ! test_gauss:
  !
  ! Unit Tester for Gaussian Points, Weights
  ! ==============================================================

! Subprogram not used   subroutine test_gauss(npts)
! Subprogram not used     use kinds, only: real_kind
! Subprogram not used 
! Subprogram not used     integer, intent(in) :: npts
! Subprogram not used     type (quadrature_t) :: gs
! Subprogram not used 
! Subprogram not used     integer i
! Subprogram not used     real (kind=real_kind) :: gssum
! Subprogram not used     gs=gauss(npts)
! Subprogram not used 
! Subprogram not used     print *
! Subprogram not used     print *,"============================================"
! Subprogram not used     print *,"        Testing Gaussian Quadrature..."
! Subprogram not used     print *
! Subprogram not used     print *,"          points              weights"
! Subprogram not used     print *,"============================================"
! Subprogram not used     do i=1,npts
! Subprogram not used        print *,i,gs%points(i),gs%weights(i)
! Subprogram not used     end do
! Subprogram not used     print *,"============================================"
! Subprogram not used     gssum=SUM(gs%weights(:))
! Subprogram not used     print *,"sum of Gaussian weights=",gssum
! Subprogram not used     print *,"============================================"
! Subprogram not used 
! Subprogram not used     deallocate(gs%points)
! Subprogram not used     deallocate(gs%weights)
! Subprogram not used 
! Subprogram not used   end subroutine test_gauss

  ! ==============================================================
  ! gausslobatto:
  !
  ! Find the Gauss-Lobatto Legendre collocation points xgl(i) and the
  ! corresponding weights.
  !
  ! ==============================================================

  function gausslobatto(npts) result(gll)

    integer, intent(in) :: npts
    type (quadrature_t) :: gll

    allocate(gll%points(npts))
    allocate(gll%weights(npts))

    gll%points=gausslobatto_pts(npts)
    gll%weights=gausslobatto_wts(npts,gll%points)

  end function gausslobatto

  ! ==============================================================
  ! gausslobatto_pts:
  !
  ! Compute the Gauss-Lobatto Collocation points 
  ! for Jacobi Polynomials
  !
  ! ==============================================================

  function gausslobatto_pts(np1) result(pts)
    use physical_constants, only : QQ_PI

    integer, intent(in)     :: np1        ! Number of velocity grid points
    real (kind=longdouble_kind) :: pts(np1)

    ! Local variables

    real (kind=longdouble_kind) :: alpha,beta
    real (kind=longdouble_kind) :: xjac(0:np1-1)
    real (kind=longdouble_kind) :: jac(0:np1)
    real (kind=longdouble_kind) :: jacm1(0:np1)
    real (kind=longdouble_kind) :: djac(0:np1)

    integer  prec                    ! number of mantissa bits 
    real (kind=longdouble_kind) eps      ! machine epsilon
    real (kind=longdouble_kind), parameter :: convthresh = 10  ! convergence threshold relative 
    ! to machine epsilon 
    integer, parameter :: kstop = 30 ! max iterations for polynomial deflation

    real (kind=longdouble_kind) :: a,b,det
    real (kind=longdouble_kind) :: poly
    real (kind=longdouble_kind) :: pder
    real (kind=longdouble_kind) :: recsum,thresh
    real (kind=longdouble_kind) :: dth,cd,sd,cs,ss,cstmp

    real (kind=longdouble_kind) :: x
    real (kind=longdouble_kind) :: delx
    real (kind=longdouble_kind) :: c0,c1,c2,c10

    integer i,j,k
    integer n, nh

    n  = np1 - 1
    c0 = 0.0_longdouble_kind
    c1 = 1.0_longdouble_kind
    c2 = 2.0_longdouble_kind
    c10 = 10.0_longdouble_kind

    alpha = c0
    beta  = c0

    ! =========================================================
    ! compute machine precision and set the convergence
    ! threshold thresh to 10 times that level 
    ! =========================================================

    prec   = PRECISION(c10)
    eps    = c10**(-prec)
    thresh = convthresh*eps

    ! =====================================================
    ! initialize the end points
    ! =====================================================

    xjac(0) =  c1
    xjac(n) = -c1

    ! ============================================================
    ! Compute first half of the roots by "polynomial deflation".
    ! ============================================================

    ! ============================================================
    ! compute the parameters in the polynomial whose 
    ! roots are desired...
    ! ============================================================

    call jacobi(n+1, c1,alpha,beta,jac(0:n+1),djac(0:n+1))
    call jacobi(n+1,-c1,alpha,beta,jacm1(0:n+1),djac(0:n+1))

    det =   jac(n  )*jacm1(n-1)-jacm1(n  )*jac(n-1)
    a   = -(jac(n+1)*jacm1(n-1)-jacm1(n+1)*jac(n-1))/det
    b   = -(jac(n  )*jacm1(n+1)-jacm1(n  )*jac(n+1))/det

    dth = QQ_PI/(2*n+1)
    cd  = COS(c2*dth)
    sd  = SIN(c2*dth)
    cs  = COS(dth)
    ss  = SIN(dth)

    nh  = (n+1)/2

    do j=1,nh-1
       x=cs          ! first guess at root 
       k=0
       delx=c1
       do while(k<kstop .and. ABS(delx) > thresh)
          call jacobi(n+1,x,alpha,beta,jac(0:n+1),djac(0:n+1))
          poly =  jac(n+1)+a* jac(n)+b* jac(n-1)
          pder = djac(n+1)+a*djac(n)+b*djac(n-1)
          recsum=c0
          do i=0,j-1
             recsum = recsum + c1/(x-xjac(i))
          end do
          delx = -poly/(pder-recsum*poly)
          x = x + delx
          k = k + 1
       end do

       xjac(j)=x

       ! =====================================================
       !  compute the guesses for the roots
       !  for the next points, i.e :
       !
       !  ss = sn(theta) => sin(theta+2*dth)
       !  cs = cs(theta) => cs(theta+2*dth)
       ! =====================================================

       cstmp=cs*cd-ss*sd      
       ss=cs*sd+ss*cd    
       cs=cstmp          
    end do

    ! ================================================
    ! compute the second half of the roots by symmetry
    ! ================================================

    do j=1,nh 
       xjac(n-j) = -xjac(j)
    end do

    if (MODULO(n,2)==0) xjac(nh)=c0

    ! ====================================================
    ! Reverse the sign of everything so that indexing
    ! increases with position          
    ! ====================================================

    do j=0,n
       pts(j+1) = -xjac(j)
    end do

  end function gausslobatto_pts

  ! ================================================
  ! Gauss Lobatto Legendre Weights   
  ! ================================================ 

  function gausslobatto_wts(np1, glpts) result(wts)

    integer, intent(in)                 :: np1
    real (kind=longdouble_kind), intent(in) :: glpts(np1)
    real (kind=longdouble_kind)             :: wts(np1)

    ! Local variables

    real (kind=longdouble_kind) :: c0,c2
    real (kind=longdouble_kind) :: alpha
    real (kind=longdouble_kind) :: beta
    real (kind=longdouble_kind) :: jac(np1)
    integer i,n

    c0    = 0.0_longdouble_kind
    c2    = 2.0_longdouble_kind
    alpha = c0
    beta  = c0
    n     = np1-1

    jac=jacobi_polynomials(n,alpha,beta,np1,glpts)

    do i=1,np1
       wts(i)=c2/(n*(n+1)*jac(i)*jac(i))
    end do

  end function gausslobatto_wts

  ! ==============================================================
  ! test_gausslobatto:
  !
  ! Unit Tester for Gaussian Lobatto Quadrature...
  ! ==============================================================

  subroutine test_gausslobatto(npts)
    use kinds, only : real_kind
    integer, intent(in) :: npts
    type (quadrature_t) :: gll

    integer i
    real (kind=real_kind) :: gllsum
    gll=gausslobatto(npts)

    print *
    print *,"============================================"
    print *,"      Testing Gauss-Lobatto Quadrature..."
    print *
    print *,"          points              weights"
    print *,"============================================"
    do i=1,npts
       print *,i,gll%points(i),gll%weights(i)
    end do
    print *,"============================================"
    gllsum=SUM(gll%weights(:))
    print *,"sum of Gauss-Lobatto weights=",gllsum
    print *,"============================================"

    deallocate(gll%points)
    deallocate(gll%weights)

  end subroutine test_gausslobatto

  ! ================================================
  !
! Subprogram not used   ! subroutine jacobi:
! Subprogram not used   !
! Subprogram not used   !  Computes the Jacobi Polynomials (jac) and their
! Subprogram not used   !  first derivatives up to and including degree n 
! Subprogram not used   !  at point x on the interval (-1,1).
! Subprogram not used   !
! Subprogram not used   !    See for example the recurrence relations 
! Subprogram not used   !    in equation 2.5.4 (page 70) in 
! Subprogram not used   !
! Subprogram not used   !    "Spectral Methods in Fluid Dynamics",
! Subprogram not used   !    by C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A.Zang
! Subprogram not used   !    Springer-Verlag, 1988.
! Subprogram not used   ! ================================================
! Subprogram not used 
! Subprogram not used   subroutine jacobi(n, x, alpha, beta, jac, djac)
! Subprogram not used 
! Subprogram not used     integer, intent(in)                 :: n
! Subprogram not used     real (kind=longdouble_kind), intent(in) :: x
! Subprogram not used     real (kind=longdouble_kind), intent(in) :: alpha
! Subprogram not used     real (kind=longdouble_kind), intent(in) :: beta
! Subprogram not used     real (kind=longdouble_kind)             :: jac(0:n)
! Subprogram not used     real (kind=longdouble_kind)             :: djac(0:n)
! Subprogram not used 
! Subprogram not used     ! Local variables
! Subprogram not used 
! Subprogram not used     real (kind=longdouble_kind) :: a1k
! Subprogram not used     real (kind=longdouble_kind) :: a2k
! Subprogram not used     real (kind=longdouble_kind) :: a3k
! Subprogram not used     real (kind=longdouble_kind) :: da2kdx
! Subprogram not used 
! Subprogram not used     real (kind=longdouble_kind) :: c2,c1,c0
! Subprogram not used 
! Subprogram not used     integer ::  k
! Subprogram not used 
! Subprogram not used     c0 = 0.0_longdouble_kind
! Subprogram not used     c1 = 1.0_longdouble_kind
! Subprogram not used     c2 = 2.0_longdouble_kind
! Subprogram not used 
! Subprogram not used     jac(0)=c1
! Subprogram not used     jac(1)=(c1 + alpha)*x
! Subprogram not used 
! Subprogram not used     djac(0)=c0
! Subprogram not used     djac(1)=(c1 + alpha)
! Subprogram not used 
! Subprogram not used     do k=1,n-1
! Subprogram not used        a1k      = c2*( k + c1 )*( k + alpha + beta + c1 )*( c2*k + alpha + beta )
! Subprogram not used        da2kdx   = ( c2*( k + c1 ) + alpha + beta )*( c2*k + alpha + beta + c1 )*( c2*k + alpha + beta )
! Subprogram not used        a2k      = ( c2*k + alpha + beta + c1 )*( alpha*alpha - beta*beta ) + x*da2kdx
! Subprogram not used        a3k      = c2*(k + alpha)*( k + beta )*( c2*k + alpha + beta + c2 )
! Subprogram not used        jac(k+1) = ( a2k*jac(k)-a3k*jac(k-1) )/a1k
! Subprogram not used        djac(k+1)= ( a2k*djac(k) + da2kdx*jac(k) - a3k*djac(k-1) )/a1k          
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used   end subroutine jacobi


  ! ==========================================================
  ! This routine computes the Nth order Jacobi Polynomials 
  ! (jac) for a vector of positions x on the interval (-1,1),
  ! of length npoints.
  !
  !    See for example the recurrence relations 
  !    in equation 2.5.4 (page 70) in 
  !
  !     "Spectral Methods in Fluid Dynamics",
  !     by C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A.Zang
  !     Springer-Verlag, 1988.
  !
  ! ===========================================================

  function jacobi_polynomials(n, alpha, beta, npoints, x) result(jac)

    integer, intent(in)     :: n         ! order of the Jacobi Polynomial
    real (kind=longdouble_kind) :: alpha 
    real (kind=longdouble_kind) :: beta
    integer, intent(in)     :: npoints
    real (kind=longdouble_kind) :: x(npoints)
    real (kind=longdouble_kind) :: jac(npoints)

    ! Local variables

    real (kind=longdouble_kind) :: a1k
    real (kind=longdouble_kind) :: a2k
    real (kind=longdouble_kind) :: a3k
    real (kind=longdouble_kind) :: da2kdx

    real (kind=longdouble_kind) :: jacp1
    real (kind=longdouble_kind) :: jacm1
    real (kind=longdouble_kind) :: jac0
    real (kind=longdouble_kind) :: xtmp

    real (kind=longdouble_kind) :: c2,c1,c0
    integer j,k

    c0 = 0.0_longdouble_kind
    c1 = 1.0_longdouble_kind
    c2 = 2.0_longdouble_kind

    do j = 1,npoints

       xtmp=x(j)

       jacm1=c1
       jac0 =(c1+alpha)*xtmp

       do k=1,n-1
          a1k=c2*(k+c1)*(k+alpha+beta+c1)*(c2*k+alpha+beta)
          da2kdx=(c2*k+alpha+beta+c2)*(c2*k+alpha+beta+c1)*(c2*k+alpha+beta)
          a2k=(c2*k+alpha+beta+c1)*(alpha*alpha-beta*beta) + xtmp*da2kdx
          a3k=c2*(k+alpha)*(k+beta)*(c2*k+alpha+beta+c2)
          jacp1=(a2k*jac0-a3k*jacm1)/a1k
          jacm1=jac0
          jac0 =jacp1
       end do

       if (n==0)jac0=jacm1
       jac(j)=jac0
    end do

  end function jacobi_polynomials

  ! ================================================
  ! This routine computes the first derivatives of Nth
  ! order Jacobi Polynomials (djac) for a vector of 
  ! positions x on the interval (-1,1), of length npoints.
  !
  ! See for example the recurrence relations 
  ! in equation 2.5.4 (page 70) in 
  !
  ! "Spectral Methods in Fluid Dynamics",
  ! by C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A.Zang
  ! Springer-Verlag, 1988.
  !
  ! ================================================

! Subprogram not used   function jacobi_derivatives(n, alpha, beta, npoints, x) result(djac)
! Subprogram not used 
! Subprogram not used     integer                , intent(in) :: n         ! order of the Jacobi Polynomial
! Subprogram not used     real (kind=longdouble_kind), intent(in) :: alpha 
! Subprogram not used     real (kind=longdouble_kind), intent(in) :: beta
! Subprogram not used     integer                , intent(in) :: npoints
! Subprogram not used     real (kind=longdouble_kind), intent(in) :: x(npoints)
! Subprogram not used 
! Subprogram not used     real (kind=longdouble_kind)             :: djac(npoints)
! Subprogram not used 
! Subprogram not used     ! Local variables
! Subprogram not used 
! Subprogram not used     ! Local variables
! Subprogram not used 
! Subprogram not used     real (kind=longdouble_kind) :: a1k
! Subprogram not used     real (kind=longdouble_kind) :: a2k
! Subprogram not used     real (kind=longdouble_kind) :: a3k
! Subprogram not used     real (kind=longdouble_kind) :: da2kdx
! Subprogram not used 
! Subprogram not used     real (kind=longdouble_kind) :: jacp1
! Subprogram not used     real (kind=longdouble_kind) :: jacm1
! Subprogram not used     real (kind=longdouble_kind) :: jac0
! Subprogram not used     real (kind=longdouble_kind) :: djacp1
! Subprogram not used     real (kind=longdouble_kind) :: djacm1
! Subprogram not used     real (kind=longdouble_kind) :: djac0
! Subprogram not used 
! Subprogram not used     real (kind=longdouble_kind) :: xtmp
! Subprogram not used 
! Subprogram not used     real (kind=longdouble_kind) :: c2,c1,c0
! Subprogram not used     integer j,k
! Subprogram not used 
! Subprogram not used     c0 = 0.0_longdouble_kind
! Subprogram not used     c1 = 1.0_longdouble_kind
! Subprogram not used     c2 = 2.0_longdouble_kind
! Subprogram not used 
! Subprogram not used     do j = 1,npoints
! Subprogram not used 
! Subprogram not used        xtmp=x(j)
! Subprogram not used 
! Subprogram not used        jacm1=c1
! Subprogram not used        jac0 =(c1+alpha)*xtmp
! Subprogram not used 
! Subprogram not used        djacm1 = c0
! Subprogram not used        djac0  = (c1+alpha)
! Subprogram not used 
! Subprogram not used        do k=1,n-1
! Subprogram not used           a1k=c2*(k+c1)*(k+alpha+beta+c1)*(c2*k+alpha+beta)
! Subprogram not used           da2kdx=(c2*k+alpha+beta+c2)*(c2*k+alpha+beta+c1)*(c2*k+alpha+beta)
! Subprogram not used           a2k=(c2*k+alpha+beta+c1)*(alpha*alpha-beta*beta) + xtmp*da2kdx
! Subprogram not used           a3k=c2*(k+alpha)*(k+beta)*(c2*k+alpha+beta+c2)
! Subprogram not used 
! Subprogram not used           jacp1=(a2k*jac0-a3k*jacm1)/a1k
! Subprogram not used           djacp1=(a2k*djac0+da2kdx*jac0-a3k*djacm1)/a1k
! Subprogram not used 
! Subprogram not used           jacm1=jac0
! Subprogram not used           jac0=jacp1
! Subprogram not used 
! Subprogram not used           djacm1=djac0
! Subprogram not used           djac0=djacp1
! Subprogram not used 
! Subprogram not used        end do
! Subprogram not used 
! Subprogram not used        if (n==0)djac0=djacm1
! Subprogram not used        djac(j)=djac0
! Subprogram not used 
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used   end function jacobi_derivatives

  ! ===================================================
  !
  ! legendre:
  !
  ! Compute the legendre polynomials using
  ! the recurrence relationship.
  ! return leg(m+1) = P_N(x) for m=0..N
  ! p_3 = Legendre polynomial of degree N
  ! p_2 = Legendre polynomial of degree N-1 at x
  ! p_1 = Legendre polynomial of degree N-2 at x
  !
  ! ===================================================

  function legendre(x,N) result(leg)

    integer   :: N
    real (kind=longdouble_kind) :: x
    real (kind=longdouble_kind) :: leg(N+1)

    real (kind=longdouble_kind) ::  p_1, p_2, p_3
    integer   :: k

    p_3 = 1.0_longdouble_kind
    leg(1)=p_3
    if (n.ne.0) then
       p_2 = p_3
       p_3 = x
       leg(2)=p_3
       do k = 2,N
          p_1 = p_2
          p_2 = p_3
          p_3 = ( (2*k-1)*x*p_2 - (k-1)*p_1 ) / k
          leg(k+1)=p_3
       end do
    end if

  end function legendre


  ! ===========================================
  ! quad_norm:
  !
  ! compute normalization constants 
  ! for k=1,N order Legendre polynomials
  !
  ! e.g. gamma(k) in Canuto, page 58.
  !
  ! ===========================================

  function quad_norm(gquad,N) result(gamma)
    type (quadrature_t), intent(in) :: gquad
    integer            , intent(in) :: N

    real (kind=longdouble_kind) :: gamma(N)

    ! Local variables
    real (kind=longdouble_kind) :: leg(N)
    integer               :: i,k

    gamma(:)=0.0_longdouble_kind

    do i=1,N
       leg=legendre(gquad%points(i),N-1)
       do k=1,N
          gamma(k)= gamma(k)+leg(k)*leg(k)*gquad%weights(i)
       end do
    end do

  end function quad_norm

  ! =======================
  ! TrapN:
  ! Numerical recipes
  ! =======================

! Subprogram not used   subroutine trapN(f,a,b,N,it,s)
! Subprogram not used     use kinds, only : real_kind
! Subprogram not used     INTERFACE
! Subprogram not used ! Subprogram not used        FUNCTION f(x) RESULT(f_x)   ! Function to be integrated
! Subprogram not used ! Subprogram not used          use kinds, only : real_kind
! Subprogram not used ! Subprogram not used          real(kind=real_kind), INTENT(IN) :: x
! Subprogram not used ! Subprogram not used          real(kind=real_kind) :: f_x
! Subprogram not used ! Subprogram not used        END FUNCTION f
! Subprogram not used     END INTERFACE
! Subprogram not used 
! Subprogram not used     real(kind=real_kind),intent(in) :: a,b
! Subprogram not used     integer, intent(in)             :: N
! Subprogram not used     integer, intent(inout)          :: it
! Subprogram not used     real(kind=real_kind), intent(inout) :: s
! Subprogram not used 
! Subprogram not used     real(kind=real_kind) :: ssum
! Subprogram not used     real(kind=real_kind) :: del
! Subprogram not used     real(kind=real_kind) :: rtnm
! Subprogram not used     real(kind=real_kind) :: x
! Subprogram not used 
! Subprogram not used     integer :: j
! Subprogram not used 
! Subprogram not used     if (N==1) then
! Subprogram not used        s = 0.5d0*(b-a)*(f(a) + f(b))
! Subprogram not used        it =1
! Subprogram not used     else
! Subprogram not used        ssum = 0.0d0
! Subprogram not used        rtnm =1.0D0/it
! Subprogram not used        del = (b-a)*rtnm
! Subprogram not used        x=a+0.5*del
! Subprogram not used        do j=1,it
! Subprogram not used           ssum = ssum + f(x)
! Subprogram not used           x=x+del
! Subprogram not used        end do
! Subprogram not used        s=0.5d0*(s + del*ssum)
! Subprogram not used        it=2*it  
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used   end subroutine trapN

  ! ==========================================
  ! Trapezoid Rule for integrating functions 
  ! from a to b with residual error eps
  ! ==========================================

! Subprogram not used   function trapezoid(f,a,b,eps) result(Integral)
! Subprogram not used     use kinds, only : real_kind
! Subprogram not used 
! Subprogram not used     integer, parameter :: Nmax = 25  ! At most 2^Nmax + 1 points in integral
! Subprogram not used 
! Subprogram not used     INTERFACE
! Subprogram not used        FUNCTION f(x) RESULT(f_x)   ! Function to be integrated
! Subprogram not used          use kinds, only : real_kind
! Subprogram not used          real(kind=real_kind), INTENT(IN) :: x
! Subprogram not used          real(kind=real_kind) :: f_x
! Subprogram not used        END FUNCTION f
    END INTERFACE

    real(kind=real_kind), intent(in) :: a,b       ! The integral bounds
    real(kind=real_kind), intent(in) :: eps       ! relative error bound for integral
    real(kind=real_kind)             :: Integral  ! the integral result (within eps)
    real(kind=real_kind)             :: s         ! Integral approximation
    real(kind=real_kind)             :: sold      ! previous integral approx

    integer                          :: N
    integer                          :: it

    ! ==============================================================
    ! Calculate I here using trapezoid rule using f and a DO loop...
    ! ==============================================================

    s    = 1.0D30
    sold = 0.0D0
    N=1
    it=0
    do while(N<=Nmax .and. ABS(s-sold)>eps*ABS(sold))
       sold=s
       call trapN(f,a,b,N,it,s)
       N=N+1
    end do

    Integral = s

  end function trapezoid

  ! ==========================================
  ! Simpsons Rule for integrating functions 
  ! from a to b with residual error eps
  ! ==========================================

! Subprogram not used   function simpsons(f,a,b,eps) result(Integral)
! Subprogram not used     use kinds, only : real_kind
! Subprogram not used 
! Subprogram not used     integer, parameter :: Nmax = 25  ! At most 2^Nmax + 1 points in integral
! Subprogram not used 
! Subprogram not used     INTERFACE
! Subprogram not used        FUNCTION f(x) RESULT(f_x)   ! Function to be integrated
! Subprogram not used          use kinds, only : real_kind
! Subprogram not used          real(kind=real_kind), INTENT(IN) :: x
! Subprogram not used          real(kind=real_kind) :: f_x
! Subprogram not used        END FUNCTION f
    END INTERFACE

    real(kind=real_kind), intent(in) :: a,b       ! The integral bounds
    real(kind=real_kind), intent(in) :: eps       ! relative error bound for integral
    real(kind=real_kind)             :: Integral  ! the integral result (within eps)
    real(kind=real_kind)             :: s         ! Integral approximation
    real(kind=real_kind)             :: os        ! previous integral approx
    real(kind=real_kind)             :: st        ! Integral approximation
    real(kind=real_kind)             :: ost       ! previous integral approx

    integer                          :: N
    integer                          :: it

    ! ==============================================================
    ! Calculate I here using trapezoid rule using f and a DO loop...
    ! ==============================================================

    ost= 0.0D0
    s  = 1.0D30
    os = 0.0D0

    N=1
    it=0
    do while ((N<=Nmax .and. ABS(s-os)>eps*ABS(os) ) .or. N<=2)
       os = s
       call trapN(f,a,b,N,it,st)
       s=(4.0D0*st-ost)/3.0D0
       ost=st
       N=N+1
    end do

    Integral = s

  end function simpsons


  ! ==========================================
  ! gaussian_int:
  !
  ! Gaussian Quadrature Rule for integrating 
  ! function f from a to b  with gs weights and 
  ! points with precomputed gaussian quadrature 
  ! and weights.
  ! ==========================================

! Subprogram not used   function gaussian_int(f,a,b,gs) result(Integral)
! Subprogram not used     use kinds, only : real_kind
! Subprogram not used 
! Subprogram not used     integer, parameter :: Nmax = 10  ! At most 2^Nmax + 1 points in integral
! Subprogram not used 
! Subprogram not used     INTERFACE
! Subprogram not used        FUNCTION f(x) RESULT(f_x)   ! Function to be integrated
! Subprogram not used          use kinds, only : real_kind
! Subprogram not used          real(kind=real_kind), INTENT(IN) :: x
! Subprogram not used          real(kind=real_kind) :: f_x
! Subprogram not used        END FUNCTION f
    END INTERFACE

    real(kind=real_kind), intent(in) :: a,b       ! The integral bounds
    type(quadrature_t), intent(in)   :: gs        ! gaussian points/wts
    real(kind=real_kind)             :: Integral  ! the integral result (within eps)

    integer                          :: i
    real (kind=real_kind)            :: s,x
    ! ==============================================================
    ! Calculate I = S f(x)dx here using gaussian quadrature
    ! ==============================================================

    s = 0.0D0
    do i=1,SIZE(gs%points)
       x = 0.50D0*((b-a)*gs%points(i) + (b+a))
       s = s + gs%weights(i)*f(x)
    end do
    Integral = s*(0.5D0*(b-a))

  end function gaussian_int

end module quadrature_mod





