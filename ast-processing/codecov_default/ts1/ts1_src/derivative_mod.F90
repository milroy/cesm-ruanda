module derivative_mod
  use kinds, only : real_kind, longdouble_kind
  use dimensions_mod, only : np, nc, npdg, nep
  use quadrature_mod, only : quadrature_t, gauss, gausslobatto,legendre, jacobi
  use parallel_mod, only : abortmp
  ! needed for spherical differential operators:
  use physical_constants, only : rrearth 
  use element_mod, only : element_t
  use control_mod, only : hypervis_scaling, hypervis_power

implicit none
private

  type, public :: derivative_t
     real (kind=real_kind) :: Dvv(np,np)
     real (kind=real_kind) :: Dvv_diag(np,np)
     real (kind=real_kind) :: Dvv_twt(np,np)
     real (kind=real_kind) :: Mvv_twt(np,np)  ! diagonal matrix of GLL weights
     real (kind=real_kind) :: vvtemp(np,np)
     real (kind=real_kind) :: vvtempt(np,np,2)
     real (kind=real_kind) :: Mfvm(np,nc+1)
     real (kind=real_kind) :: Cfvm(np,nc)
     real (kind=real_kind) :: Sfvm(np,nep)
     real (kind=real_kind) :: legdg(np,np)
  end type derivative_t

  type, public :: derivative_stag_t
     real (kind=real_kind) :: D(np,np)
     real (kind=real_kind) :: M(np,np)
     real (kind=real_kind) :: Dpv(np,np)
     real (kind=real_kind) :: D_twt(np,np)
     real (kind=real_kind) :: M_twt(np,np)
     real (kind=real_kind) :: M_t(np,np)
     real (kind=real_kind) :: vtemp(np,np,2)
     real (kind=real_kind) :: vtempt(np,np,2)
  end type derivative_stag_t

  real (kind=real_kind), allocatable :: integration_matrix(:,:)
  private :: allocate_subcell_integration_matrix

! ======================================
! Public Interfaces
! ======================================

  public :: subcell_integration

  public :: derivinit
  public :: deriv_print

  public :: gradient
  public :: gradient_wk
  public :: vorticity
  public :: divergence

  public :: interpolate_gll2fvm_corners
  public :: interpolate_gll2fvm_points
  public :: interpolate_gll2spelt_points
  public :: remap_phys2gll


  interface divergence
      module procedure divergence_nonstag
      module procedure divergence_stag
  end interface

  interface gradient
      module procedure gradient_str_nonstag
      module procedure gradient_str_stag
  end interface

  interface gradient_wk
      module procedure gradient_wk_nonstag
      module procedure gradient_wk_stag
  end interface

  public :: v2pinit

  private :: dmatinit
  private :: dvvinit
  private :: dpvinit

! these routines compute spherical differential operators as opposed to
! the gnomonic coordinate operators above.  Vectors (input or output)
! are always expressed in lat-lon coordinates
!
! note that weak derivatives (integrated by parts form) can be defined using
! contra or co-variant test functions, so 
!
  public  :: gradient_sphere
  public  :: gradient_sphere_wk_testcov
  public  :: gradient_sphere_wk_testcontra   ! only used for debugging
  public  :: ugradv_sphere
  public  :: vorticity_sphere
  public  :: vorticity_sphere_diag
  public  :: divergence_sphere
  public  :: curl_sphere
  public  :: curl_sphere_wk_testcov
!  public  :: curl_sphere_wk_testcontra  ! not coded
  public  :: divergence_sphere_wk
  public  :: laplace_sphere_wk
  public  :: vlaplace_sphere_wk
  public  :: element_boundary_integral
  public  :: edge_flux_u_cg
  public  :: gll_to_dgmodal
  public  :: dgmodal_to_gll



contains

! ==========================================
! derivinit:
!
! Initialize the matrices for taking 
! derivatives and interpolating
! ==========================================

  subroutine derivinit(deriv,fvm_corners, fvm_points, spelt_refnep)
    type (derivative_t)      :: deriv
!    real (kind=longdouble_kind),optional :: phys_points(:)
    real (kind=longdouble_kind),optional :: fvm_corners(nc+1)
    real (kind=longdouble_kind),optional :: fvm_points(nc)
    real (kind=longdouble_kind),optional :: spelt_refnep(nep)

    ! Local variables
    type (quadrature_t) :: gp   ! Quadrature points and weights on pressure grid
    
    real (kind=longdouble_kind) :: dmat(np,np)
    real (kind=longdouble_kind) :: dpv(np,np)
    real (kind=longdouble_kind) :: v2p(np,np)
    real (kind=longdouble_kind) :: p2v(np,np)
    real (kind=longdouble_kind) :: dvv(np,np)
    real (kind=longdouble_kind) :: dvv_diag(np,np)
    real (kind=longdouble_kind) :: v2v(np,np)
    real (kind=longdouble_kind) :: xnorm
    integer i,j

    ! ============================================
    ! initialize matrices in longdouble_kind precision
    ! and transfer results into real_kind
    ! floating point precision
    ! ============================================

    gp=gausslobatto(np)

    ! Legendre polynomials of degree npdg-1, on the np GLL grid:
    if (npdg>np) call abortmp( 'FATAL ERROR: npdg>np')
    if (npdg>0 .and. npdg<np) then
       ! in this case, we use a DG basis of Legendre polynomials
       ! stored at the GLL points.  
       do i=1,np
          deriv%legdg(1:npdg,i) = legendre(gp%points(i),npdg-1)
       end do
       ! normalize
       do j=1,npdg
          xnorm=sqrt(sum(deriv%legdg(j,:)*deriv%legdg(j,:)*gp%weights(:)))
          deriv%legdg(j,:)=deriv%legdg(j,:)/xnorm
       enddo
    endif
    call dvvinit(dvv,gp)
    deriv%Dvv(:,:)   = dvv(:,:)

    do i=1,np
       do j=1,np
          if (i.eq.j) then
             deriv%dvv_diag(i,j)   = dvv(i,j)
          else
             deriv%dvv_diag(i,j) = 0.0D0
          endif 
        end do
     end do







    v2v = 0.0D0
    do i=1,np
       v2v(i,i) = gp%weights(i)
    end do

    do i=1,np
       do j=1,np
          dvv(j,i) = dvv(j,i)*gp%weights(i)
       end do
    end do

    deriv%Dvv_twt = TRANSPOSE(dvv)
    deriv%Mvv_twt = v2v
    if (present(fvm_corners)) &
         call v2pinit(deriv%Mfvm,gp%points,fvm_corners,np,nc+1)

    if (present(fvm_points)) &
         call v2pinit(deriv%Cfvm,gp%points,fvm_points,np,nc)
         
    if (present(spelt_refnep)) &     
      call v2pinit(deriv%Sfvm,gp%points,spelt_refnep,np,nep)
         
    ! notice we deallocate this memory here even though it was allocated 
    ! by the call to gausslobatto.
    deallocate(gp%points)
    deallocate(gp%weights)

  end subroutine derivinit

! Subprogram not used   subroutine deriv_print(deriv)
! Subprogram not used     type (derivative_t) :: deriv
! Subprogram not used     
! Subprogram not used     ! Local variables
! Subprogram not used 
! Subprogram not used     integer j
! Subprogram not used     print *,"Derivative Matrix Dvv"
! Subprogram not used     do j=1,np
! Subprogram not used        write(6,*)deriv%Dvv(:,j)
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     print *,"Weak Derivative Matrix Dvv_twt"
! Subprogram not used     do j=1,np
! Subprogram not used        write(6,*)deriv%Dvv_twt(:,j)
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   end subroutine deriv_print

! =======================================
! dmatinit:
!
! Compute rectangular v->p 
! derivative matrix (dmat)
! =======================================

! Subprogram not used   subroutine dmatinit(dmat)
! Subprogram not used 
! Subprogram not used     real (kind=longdouble_kind) :: dmat(np,np)
! Subprogram not used 
! Subprogram not used     ! Local variables
! Subprogram not used 
! Subprogram not used     type (quadrature_t) :: gll
! Subprogram not used     type (quadrature_t) :: gs
! Subprogram not used 
! Subprogram not used     integer i,j
! Subprogram not used     real(kind=longdouble_kind)  fact,f1,f2
! Subprogram not used     real(kind=longdouble_kind)  func0,func1
! Subprogram not used     real(kind=longdouble_kind)  dis,c0,c1
! Subprogram not used 
! Subprogram not used     real(kind=longdouble_kind)  :: leg(np,np)
! Subprogram not used     real(kind=longdouble_kind)  ::  jac(0:np-1)
! Subprogram not used     real(kind=longdouble_kind)  :: djac(0:np-1)
! Subprogram not used 
! Subprogram not used     c0 = 0.0_longdouble_kind
! Subprogram not used     c1 = 1.0_longdouble_kind
! Subprogram not used 
! Subprogram not used     gll= gausslobatto(np)
! Subprogram not used     gs = gauss(np)
! Subprogram not used 
! Subprogram not used     ! =============================================================
! Subprogram not used     ! Compute Legendre polynomials on Gauss-Lobatto grid (velocity)
! Subprogram not used     ! =============================================================
! Subprogram not used 
! Subprogram not used     do i=1,np
! Subprogram not used        leg(:,i) = legendre(gll%points(i),np-1)
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     ! ================================================================
! Subprogram not used     !  Derivatives of velocity cardinal functions on pressure grid
! Subprogram not used     !  d(i,j) = D(j,i) = D' (D-transpose) since D(i,j) = dh_j(x_i)/dx
! Subprogram not used     ! ================================================================
! Subprogram not used 
! Subprogram not used     fact = np*(np-1)
! Subprogram not used 
! Subprogram not used     do j=1,np
! Subprogram not used        call jacobi(np-1,gs%points(j),c0,c0,jac(0:np-1),djac(0:np-1))
! Subprogram not used        func0 =  jac(np-1)
! Subprogram not used        func1 = djac(np-1)
! Subprogram not used        f1 = fact*func0
! Subprogram not used        f2 = (c1 - gs%points(j))*(c1 + gs%points(j)) * func1
! Subprogram not used        do i = 1, np
! Subprogram not used           if ( gs%points(j) /= gll%points(i) ) then
! Subprogram not used              dis = gs%points(j) - gll%points(i)
! Subprogram not used              dmat(i,j) = func0 / ( leg(np,i)*dis ) + f2 / (fact*leg(np,i)*dis*dis)
! Subprogram not used !!! OTHER             dmat(i,j) = (1.0D0/(fact*leg(np,i)*dis*dis))* (func0*fact*dis + f2)
! Subprogram not used           else
! Subprogram not used              dmat(i,j) = c0
! Subprogram not used           endif
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     deallocate(gll%points)
! Subprogram not used     deallocate(gll%weights)
! Subprogram not used 
! Subprogram not used 	deallocate(gs%points)
! Subprogram not used 	deallocate(gs%weights)
! Subprogram not used 
! Subprogram not used end subroutine dmatinit

! =======================================
! dpvinit:
!
! Compute rectangular p->v
! derivative matrix (dmat) 
! for strong gradients
! =======================================

! Subprogram not used subroutine dpvinit(dmat)
! Subprogram not used 
! Subprogram not used real (kind=longdouble_kind) :: dmat(np,np)
! Subprogram not used 
! Subprogram not used ! Local variables
! Subprogram not used 
! Subprogram not used type (quadrature_t) :: gll
! Subprogram not used type (quadrature_t) :: gs
! Subprogram not used 
! Subprogram not used integer i,j
! Subprogram not used real(kind=longdouble_kind)  dis,c0,c1
! Subprogram not used 
! Subprogram not used real(kind=longdouble_kind)  :: legv(0:np,np)
! Subprogram not used real(kind=longdouble_kind)  :: dlegv(0:np,np)
! Subprogram not used 
! Subprogram not used real(kind=longdouble_kind)  :: leg(0:np)
! Subprogram not used real(kind=longdouble_kind)  :: dleg(0:np)
! Subprogram not used 
! Subprogram not used c0 = 0.0_longdouble_kind
! Subprogram not used c1 = 1.0_longdouble_kind
! Subprogram not used 
! Subprogram not used gll= gausslobatto(np)
! Subprogram not used gs = gauss(np)
! Subprogram not used 
! Subprogram not used ! =============================================================
! Subprogram not used ! Compute Legendre polynomials on Gauss-Lobatto grid (velocity)
! Subprogram not used ! =============================================================
! Subprogram not used 
! Subprogram not used do i=1,np
! Subprogram not used call jacobi(np,gll%points(i),c0,c0,legv(0:np,i),dlegv(0:np,i))
! Subprogram not used end do
! Subprogram not used 
! Subprogram not used ! ================================================================
! Subprogram not used !  Derivatives of velocity cardinal functions on pressure grid
! Subprogram not used     !  d(i,j) = D(j,i) = D' (D-transpose) since D(i,j) = dh_j(x_i)/dx
! Subprogram not used     ! ================================================================
! Subprogram not used 
! Subprogram not used     do j=1,np
! Subprogram not used        call jacobi(np,gs%points(j),c0,c0,leg(0:np),dleg(0:np))
! Subprogram not used        do i = 1, np
! Subprogram not used           if ( gs%points(j) /= gll%points(i) ) then
! Subprogram not used              dis = gll%points(i) - gs%points(j)
! Subprogram not used              dmat(j,i) = dlegv(np,i)/( dleg(np)*dis ) -  legv(np,i)/ (dleg(np)*dis*dis)
! Subprogram not used           else
! Subprogram not used              dmat(j,i) = c0
! Subprogram not used           endif
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     deallocate(gll%points)
! Subprogram not used     deallocate(gll%weights)
! Subprogram not used 
! Subprogram not used     deallocate(gs%points)
! Subprogram not used     deallocate(gs%weights)
! Subprogram not used 
! Subprogram not used   end subroutine dpvinit

! =======================================
! v2pinit:
! Compute interpolation matrix from gll(1:n1) -> gs(1:n2)
! =======================================
  subroutine v2pinit(v2p,gll,gs,n1,n2)
    integer :: n1,n2
    real(kind=real_kind)  ::  v2p(n1,n2)
    real(kind=real_kind)  ::  v2p_new(n1,n2)
    real(kind=longdouble_kind)  ::  gll(n1),gs(n2)
    ! Local variables

    integer i,j,k,m,l
    real(kind=longdouble_kind)  fact,f1, sum
    real(kind=longdouble_kind)  func0,func1

    real(kind=longdouble_kind)  :: leg(n1,n1)
    real(kind=longdouble_kind)  ::  jac(0:n1-1)
    real(kind=longdouble_kind)  :: djac(0:n1-1)
    real(kind=longdouble_kind)  :: c0,c1

    type(quadrature_t) :: gll_pts
    real(kind=longdouble_kind)  :: leg_out(n1,n2)
    real(kind=longdouble_kind)  :: gamma(n1)

    c0 = 0.0_longdouble_kind
    c1 = 1.0_longdouble_kind

    ! ==============================================================
    ! Compute Legendre polynomials on Gauss-Lobatto grid (velocity)
    ! ==============================================================

    fact = -n1*(n1-1)
    do i=1,n1
       leg(:,i) = legendre(gll(i),n1-1)
       leg(n1,i) = fact * leg(n1,i)
    end do

    ! ===================================================
    !  Velocity cardinal functions on pressure grid
    ! ===================================================

    ! NEW VERSION, with no division by (gs(j)-gll(i)):

    ! compute legendre polynomials at output points:
    gll_pts = gausslobatto(n1)

    fact = -n1*(n1-1)
    do i=1,n2
       leg_out(:,i) = legendre(gs(i),n1-1)
       leg_out(n1,i) = fact * leg_out(n1,i)
    end do


    ! compute gamma: (normalization factor for inv(leg)
    do m=1,n1
       gamma(m)=0
       do i=1,n1
          gamma(m)=gamma(m)+leg(m,i)*leg(m,i)*gll_pts%weights(i) 
       enddo
       gamma(m)=1/gamma(m)
    enddo

    ! compute product of leg_out * inv(leg):
    do j=1,n2   ! this should be fvm points
       do l=1,n1   ! this should be GLL points
          sum=0
          do k=1,n1  ! number of polynomials = number of GLL points
             sum=sum + leg_out(k,j)*gamma(k)*leg(k,l)
          enddo
          v2p_new(l,j) = gll_pts%weights(l)*sum
       enddo
    enddo
    deallocate(gll_pts%points)
    deallocate(gll_pts%weights)


    v2p=v2p_new
  end subroutine v2pinit



! =======================================
! dvvinit:
!
! Compute rectangular v->v
! derivative matrix (dvv)
! =======================================

  subroutine dvvinit(dvv,gll)

    real(kind=longdouble_kind)  ::  dvv(np,np)
    type (quadrature_t)   :: gll

    ! Local variables

    real(kind=longdouble_kind)  :: leg(np,np)
    real(kind=longdouble_kind)  :: c0,c1,c4

    integer i,j

    c0 = 0.0_longdouble_kind
    c1 = 1.0_longdouble_kind
    c4 = 4.0_longdouble_kind

    do i=1,np
       leg(:,i) = legendre(gll%points(i),np-1)
    end do

    dvv(:,:) = c0
    do j=1,np
       do i=1,j-1
          dvv(j,i) = (c1/(gll%points(i)-gll%points(j)))*leg(np,i)/leg(np,j)
       end do
       dvv(j,j) = c0
       do i=j+1,np
          dvv(j,i) = (c1/(gll%points(i)-gll%points(j)))*leg(np,i)/leg(np,j)
       end do
    end do


    dvv(np,np) = + np*(np-1)/c4
    dvv(1,1)   = - np*(np-1)/c4

  end subroutine dvvinit

!  ================================================
!  divergence_stag: 
!
!  Compute divergence (maps v grid -> p grid)
!  ================================================

! Subprogram not used   function divergence_stag(v,deriv) result(div)
! Subprogram not used 
! Subprogram not used     real(kind=real_kind), intent(in) :: v(np,np,2)
! Subprogram not used     type (derivative_stag_t)         :: deriv
! Subprogram not used 
! Subprogram not used     real(kind=real_kind) :: div(np,np)
! Subprogram not used 
! Subprogram not used     ! Local
! Subprogram not used 
! Subprogram not used     integer i
! Subprogram not used     integer j
! Subprogram not used     integer l
! Subprogram not used     logical, parameter :: UseUnroll = .TRUE.
! Subprogram not used 
! Subprogram not used     real(kind=real_kind)  sumx00,sumx01
! Subprogram not used     real(kind=real_kind)  sumy00,sumy01
! Subprogram not used     real(kind=real_kind)  sumx10,sumx11
! Subprogram not used     real(kind=real_kind)  sumy10,sumy11
! Subprogram not used 
! Subprogram not used if(MODULO(np,2) == 0 .and. UseUnroll) then 
! Subprogram not used     !JMD====================================
! Subprogram not used     !JMD  2*np*np*np Flops
! Subprogram not used     !JMD====================================
! Subprogram not used     do j=1,np,2
! Subprogram not used        do l=1,np,2
! Subprogram not used 
! Subprogram not used           sumx00=0.0d0
! Subprogram not used           sumx01=0.0d0
! Subprogram not used           sumx10=0.0d0
! Subprogram not used           sumx11=0.0d0
! Subprogram not used 
! Subprogram not used           sumy00=0.0d0
! Subprogram not used           sumy01=0.0d0
! Subprogram not used           sumy10=0.0d0
! Subprogram not used           sumy11=0.0d0
! Subprogram not used 
! Subprogram not used           do i=1,np
! Subprogram not used              sumx00 = sumx00 + deriv%D(i,l  )*v(i,j  ,1)
! Subprogram not used              sumx01 = sumx01 + deriv%D(i,l+1)*v(i,j  ,1)
! Subprogram not used              sumx10 = sumx10 + deriv%D(i,l  )*v(i,j+1,1)
! Subprogram not used              sumx11 = sumx11 + deriv%D(i,l+1)*v(i,j+1,1)
! Subprogram not used 
! Subprogram not used              sumy00 = sumy00 + deriv%M(i,l  )*v(i,j  ,2)
! Subprogram not used              sumy01 = sumy01 + deriv%M(i,l+1)*v(i,j  ,2)
! Subprogram not used              sumy10 = sumy10 + deriv%M(i,l  )*v(i,j+1,2)
! Subprogram not used              sumy11 = sumy11 + deriv%M(i,l+1)*v(i,j+1,2)
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used           deriv%vtemp(j  ,l  ,1) = sumx00
! Subprogram not used           deriv%vtemp(j  ,l+1,1) = sumx01
! Subprogram not used           deriv%vtemp(j+1,l  ,1) = sumx10
! Subprogram not used           deriv%vtemp(j+1,l+1,1) = sumx11
! Subprogram not used 
! Subprogram not used           deriv%vtemp(j  ,l  ,2) = sumy00
! Subprogram not used           deriv%vtemp(j  ,l+1,2) = sumy01
! Subprogram not used           deriv%vtemp(j+1,l  ,2) = sumy10
! Subprogram not used           deriv%vtemp(j+1,l+1,2) = sumy11
! Subprogram not used 
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     !JMD====================================
! Subprogram not used     !JMD  2*np*np*np Flops
! Subprogram not used     !JMD====================================
! Subprogram not used     do j=1,np,2
! Subprogram not used        do i=1,np,2
! Subprogram not used           sumx00=0.0d0
! Subprogram not used           sumx01=0.0d0
! Subprogram not used           sumx10=0.0d0
! Subprogram not used           sumx11=0.0d0
! Subprogram not used 
! Subprogram not used           sumy00=0.0d0
! Subprogram not used           sumy01=0.0d0
! Subprogram not used           sumy10=0.0d0
! Subprogram not used           sumy11=0.0d0
! Subprogram not used 
! Subprogram not used           do l=1,np
! Subprogram not used              sumx00 = sumx00 +  deriv%M(l,j  )*deriv%vtemp(l,i  ,1)
! Subprogram not used              sumx01 = sumx01 +  deriv%M(l,j+1)*deriv%vtemp(l,i  ,1)
! Subprogram not used              sumx10 = sumx10 +  deriv%M(l,j  )*deriv%vtemp(l,i+1,1)
! Subprogram not used              sumx11 = sumx11 +  deriv%M(l,j+1)*deriv%vtemp(l,i+1,1)
! Subprogram not used 
! Subprogram not used              sumy00 = sumy00 +  deriv%D(l,j  )*deriv%vtemp(l,i  ,2)
! Subprogram not used              sumy01 = sumy01 +  deriv%D(l,j+1)*deriv%vtemp(l,i  ,2)
! Subprogram not used              sumy10 = sumy10 +  deriv%D(l,j  )*deriv%vtemp(l,i+1,2)
! Subprogram not used              sumy11 = sumy11 +  deriv%D(l,j+1)*deriv%vtemp(l,i+1,2)
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used           div(i  ,j  ) = sumx00 + sumy00
! Subprogram not used           div(i  ,j+1) = sumx01 + sumy01
! Subprogram not used           div(i+1,j  ) = sumx10 + sumy10
! Subprogram not used           div(i+1,j+1) = sumx11 + sumy11
! Subprogram not used 
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used else
! Subprogram not used      do j=1,np
! Subprogram not used         do l=1,np
! Subprogram not used  
! Subprogram not used            sumx00=0.0d0
! Subprogram not used            sumy00=0.0d0
! Subprogram not used            do i=1,np
! Subprogram not used               sumx00 = sumx00 + deriv%D(i,l  )*v(i,j  ,1)
! Subprogram not used               sumy00 = sumy00 + deriv%M(i,l  )*v(i,j  ,2)
! Subprogram not used  	   enddo
! Subprogram not used           deriv%vtemp(j  ,l  ,1) = sumx00
! Subprogram not used           deriv%vtemp(j  ,l  ,2) = sumy00
! Subprogram not used        enddo
! Subprogram not used     enddo
! Subprogram not used     do j=1,np
! Subprogram not used        do i=1,np
! Subprogram not used           sumx00=0.0d0
! Subprogram not used 	  sumy00=0.0d0
! Subprogram not used           do l=1,np
! Subprogram not used              sumx00 = sumx00 +  deriv%M(l,j  )*deriv%vtemp(l,i  ,1)
! Subprogram not used 	     sumy00 = sumy00 +  deriv%D(l,j  )*deriv%vtemp(l,i  ,2)
! Subprogram not used 	  enddo
! Subprogram not used           div(i  ,j  ) = sumx00 + sumy00
! Subprogram not used 
! Subprogram not used 	enddo
! Subprogram not used     enddo
! Subprogram not used endif
! Subprogram not used 
! Subprogram not used   end function divergence_stag

!  ================================================
!  divergence_nonstag: 
!
!  Compute divergence (maps v->v)
!  ================================================

! Subprogram not used   function divergence_nonstag(v,deriv) result(div)
! Subprogram not used 
! Subprogram not used     real(kind=real_kind), intent(in) :: v(np,np,2)
! Subprogram not used     type (derivative_t)              :: deriv
! Subprogram not used 
! Subprogram not used     real(kind=real_kind) :: div(np,np)
! Subprogram not used 
! Subprogram not used     ! Local
! Subprogram not used 
! Subprogram not used     integer i
! Subprogram not used     integer j
! Subprogram not used     integer l
! Subprogram not used 
! Subprogram not used     logical, parameter :: UseUnroll = .TRUE.
! Subprogram not used 
! Subprogram not used     real(kind=real_kind) ::  dudx00,dudx01
! Subprogram not used     real(kind=real_kind) ::  dudx10,dudx11
! Subprogram not used 
! Subprogram not used     real(kind=real_kind) ::  dvdy00,dvdy01
! Subprogram not used     real(kind=real_kind) ::  dvdy10,dvdy11
! Subprogram not used 
! Subprogram not used if(modulo(np,2) .eq. 0 .and. UseUnroll) then
! Subprogram not used ! this is just loop unrolling - a good compiler should do it for you jpe
! Subprogram not used        do j=1,np,2
! Subprogram not used           do l=1,np,2
! Subprogram not used 
! Subprogram not used              dudx00=0.0d0
! Subprogram not used              dudx01=0.0d0
! Subprogram not used              dudx10=0.0d0
! Subprogram not used              dudx11=0.0d0
! Subprogram not used 
! Subprogram not used              dvdy00=0.0d0
! Subprogram not used              dvdy01=0.0d0
! Subprogram not used              dvdy10=0.0d0
! Subprogram not used              dvdy11=0.0d0
! Subprogram not used 
! Subprogram not used              do i=1,np
! Subprogram not used                 
! Subprogram not used                 dudx00 = dudx00 + deriv%Dvv(i,l  )*v(i,j  ,1)
! Subprogram not used                 dudx01 = dudx01 + deriv%Dvv(i,l+1)*v(i,j  ,1)
! Subprogram not used                 dudx10 = dudx10 + deriv%Dvv(i,l  )*v(i,j+1,1)
! Subprogram not used                 dudx11 = dudx11 + deriv%Dvv(i,l+1)*v(i,j+1,1)
! Subprogram not used                 
! Subprogram not used                 dvdy00 = dvdy00 + deriv%Dvv(i,l  )*v(j  ,i,2)
! Subprogram not used                 dvdy01 = dvdy01 + deriv%Dvv(i,l+1)*v(j  ,i,2)
! Subprogram not used                 dvdy10 = dvdy10 + deriv%Dvv(i,l  )*v(j+1,i,2)
! Subprogram not used                 dvdy11 = dvdy11 + deriv%Dvv(i,l+1)*v(j+1,i,2)
! Subprogram not used 
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used              div(l  ,j  ) = dudx00
! Subprogram not used              div(l+1,j  ) = dudx01
! Subprogram not used              div(l  ,j+1) = dudx10
! Subprogram not used              div(l+1,j+1) = dudx11
! Subprogram not used 
! Subprogram not used              deriv%vvtemp(j  ,l  ) = dvdy00
! Subprogram not used              deriv%vvtemp(j  ,l+1) = dvdy01
! Subprogram not used              deriv%vvtemp(j+1,l  ) = dvdy10
! Subprogram not used              deriv%vvtemp(j+1,l+1) = dvdy11
! Subprogram not used 
! Subprogram not used           end do
! Subprogram not used        end do
! Subprogram not used     else
! Subprogram not used 
! Subprogram not used        do j=1,np
! Subprogram not used           do l=1,np
! Subprogram not used              dudx00=0.0d0
! Subprogram not used              dvdy00=0.0d0
! Subprogram not used 
! Subprogram not used              do i=1,np
! Subprogram not used                 dudx00 = dudx00 + deriv%Dvv(i,l  )*v(i,j  ,1)
! Subprogram not used                 dvdy00 = dvdy00 + deriv%Dvv(i,l  )*v(j  ,i,2)
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used              div(l  ,j  ) = dudx00
! Subprogram not used              deriv%vvtemp(j  ,l  ) = dvdy00
! Subprogram not used 
! Subprogram not used 
! Subprogram not used           end do
! Subprogram not used        end do
! Subprogram not used     end if
! Subprogram not used     do j=1,np
! Subprogram not used        do i=1,np
! Subprogram not used           div(i,j)=div(i,j)+deriv%vvtemp(i,j)
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used   end function divergence_nonstag

!  ================================================
!  gradient_wk_stag:
! 
!  Compute the weak form gradient:
!  maps scalar field on the pressure grid to the
!  velocity grid
!  ================================================

! Subprogram not used   function gradient_wk_stag(p,deriv) result(dp)
! Subprogram not used 
! Subprogram not used     type (derivative_stag_t)         :: deriv
! Subprogram not used     real(kind=real_kind), intent(in) :: p(np,np)
! Subprogram not used 
! Subprogram not used     real(kind=real_kind)             :: dp(np,np,2)
! Subprogram not used 
! Subprogram not used     ! Local
! Subprogram not used       
! Subprogram not used     integer i
! Subprogram not used     integer j
! Subprogram not used     integer l
! Subprogram not used     logical, parameter :: UseUnroll = .TRUE.
! Subprogram not used 
! Subprogram not used     real(kind=real_kind)  sumx00,sumx01
! Subprogram not used     real(kind=real_kind)  sumy00,sumy01
! Subprogram not used     real(kind=real_kind)  sumx10,sumx11
! Subprogram not used     real(kind=real_kind)  sumy10,sumy11
! Subprogram not used 
! Subprogram not used     !JMD ================================
! Subprogram not used     !JMD 2*np*np*np Flops 
! Subprogram not used     !JMD ================================
! Subprogram not used 
! Subprogram not used if(MODULO(np,2) == 0 .and. UseUnroll) then 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     do j=1,np,2
! Subprogram not used        do l=1,np,2
! Subprogram not used           sumx00=0.0d0
! Subprogram not used           sumx01=0.0d0
! Subprogram not used           sumx10=0.0d0
! Subprogram not used           sumx11=0.0d0
! Subprogram not used 
! Subprogram not used           sumy00=0.0d0
! Subprogram not used           sumy01=0.0d0
! Subprogram not used           sumy10=0.0d0
! Subprogram not used           sumy11=0.0d0
! Subprogram not used 
! Subprogram not used           do i=1,np
! Subprogram not used              sumx00 = sumx00 + deriv%D_twt(i,l  )*p(i,j  )
! Subprogram not used              sumx01 = sumx01 + deriv%D_twt(i,l+1)*p(i,j  )
! Subprogram not used              sumx10 = sumx10 + deriv%D_twt(i,l  )*p(i,j+1)
! Subprogram not used              sumx11 = sumx11 + deriv%D_twt(i,l+1)*p(i,j+1)
! Subprogram not used 
! Subprogram not used              sumy00 = sumy00 + deriv%M_twt(i,l  )*p(i,j  )
! Subprogram not used              sumy01 = sumy01 + deriv%M_twt(i,l+1)*p(i,j  )
! Subprogram not used              sumy10 = sumy10 + deriv%M_twt(i,l  )*p(i,j+1)
! Subprogram not used              sumy11 = sumy11 + deriv%M_twt(i,l+1)*p(i,j+1)
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used           deriv%vtempt(j  ,l  ,1) = sumx00
! Subprogram not used           deriv%vtempt(j  ,l+1,1) = sumx01
! Subprogram not used           deriv%vtempt(j+1,l  ,1) = sumx10
! Subprogram not used           deriv%vtempt(j+1,l+1,1) = sumx11
! Subprogram not used 
! Subprogram not used           deriv%vtempt(j  ,l  ,2) = sumy00
! Subprogram not used           deriv%vtempt(j  ,l+1,2) = sumy01
! Subprogram not used           deriv%vtempt(j+1,l  ,2) = sumy10
! Subprogram not used           deriv%vtempt(j+1,l+1,2) = sumy11
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used 	
! Subprogram not used     !JMD ================================
! Subprogram not used     !JMD 2*np*np*np Flops 
! Subprogram not used     !JMD ================================
! Subprogram not used 
! Subprogram not used     do j=1,np,2
! Subprogram not used        do i=1,np,2
! Subprogram not used           sumx00=0.0d0
! Subprogram not used           sumx01=0.0d0
! Subprogram not used           sumx10=0.0d0
! Subprogram not used           sumx11=0.0d0
! Subprogram not used 
! Subprogram not used           sumy00=0.0d0
! Subprogram not used           sumy01=0.0d0
! Subprogram not used           sumy10=0.0d0
! Subprogram not used           sumy11=0.0d0
! Subprogram not used 
! Subprogram not used           do l=1,np
! Subprogram not used              sumx00 = sumx00 +  deriv%M_twt(l,j  )*deriv%vtempt(l,i  ,1)
! Subprogram not used              sumx01 = sumx01 +  deriv%M_twt(l,j+1)*deriv%vtempt(l,i  ,1)
! Subprogram not used              sumx10 = sumx10 +  deriv%M_twt(l,j  )*deriv%vtempt(l,i+1,1)
! Subprogram not used              sumx11 = sumx11 +  deriv%M_twt(l,j+1)*deriv%vtempt(l,i+1,1)
! Subprogram not used 
! Subprogram not used              sumy00 = sumy00 +  deriv%D_twt(l,j  )*deriv%vtempt(l,i  ,2)
! Subprogram not used              sumy01 = sumy01 +  deriv%D_twt(l,j+1)*deriv%vtempt(l,i  ,2)
! Subprogram not used              sumy10 = sumy10 +  deriv%D_twt(l,j  )*deriv%vtempt(l,i+1,2)
! Subprogram not used              sumy11 = sumy11 +  deriv%D_twt(l,j+1)*deriv%vtempt(l,i+1,2)
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used           dp(i  ,j  ,1) = sumx00
! Subprogram not used           dp(i  ,j+1,1) = sumx01
! Subprogram not used           dp(i+1,j  ,1) = sumx10
! Subprogram not used           dp(i+1,j+1,1) = sumx11
! Subprogram not used 
! Subprogram not used           dp(i  ,j  ,2) = sumy00
! Subprogram not used           dp(i  ,j+1,2) = sumy01
! Subprogram not used           dp(i+1,j  ,2) = sumy10
! Subprogram not used           dp(i+1,j+1,2) = sumy11
! Subprogram not used 
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used else
! Subprogram not used     do j=1,np
! Subprogram not used        do l=1,np
! Subprogram not used  	  sumx00=0.0d0
! Subprogram not used           sumy00=0.0d0
! Subprogram not used            do i=1,np
! Subprogram not used               sumx00 = sumx00 + deriv%D_twt(i,l  )*p(i,j  )
! Subprogram not used               sumy00 = sumy00 + deriv%M_twt(i,l  )*p(i,j  )
! Subprogram not used  	  enddo
! Subprogram not used            deriv%vtempt(j  ,l  ,1) = sumx00
! Subprogram not used            deriv%vtempt(j  ,l  ,2) = sumy00
! Subprogram not used         enddo
! Subprogram not used     enddo
! Subprogram not used     do j=1,np
! Subprogram not used        do i=1,np
! Subprogram not used           sumx00=0.0d0
! Subprogram not used 	  sumy00=0.0d0
! Subprogram not used           do l=1,np
! Subprogram not used              sumx00 = sumx00 +  deriv%M_twt(l,j  )*deriv%vtempt(l,i  ,1)
! Subprogram not used 	     sumy00 = sumy00 +  deriv%D_twt(l,j  )*deriv%vtempt(l,i  ,2)
! Subprogram not used 	  enddo
! Subprogram not used 	  dp(i  ,j  ,1) = sumx00
! Subprogram not used           dp(i  ,j  ,2) = sumy00
! Subprogram not used       enddo
! Subprogram not used     enddo
! Subprogram not used endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   end function gradient_wk_stag

!  ================================================
!  gradient_wk_nonstag:
! 
!  Compute the weak form gradient:
!  maps scalar field on the Gauss-Lobatto grid to the
!  weak gradient on the Gauss-Lobbatto grid
!  ================================================

! Subprogram not used   function gradient_wk_nonstag(p,deriv) result(dp)
! Subprogram not used 
! Subprogram not used     type (derivative_t)         :: deriv
! Subprogram not used     real(kind=real_kind), intent(in) :: p(np,np)
! Subprogram not used 
! Subprogram not used     real(kind=real_kind)             :: dp(np,np,2)
! Subprogram not used 
! Subprogram not used     ! Local
! Subprogram not used       
! Subprogram not used     integer i
! Subprogram not used     integer j
! Subprogram not used     integer l
! Subprogram not used     logical, parameter :: UseUnroll = .TRUE.
! Subprogram not used 
! Subprogram not used     real(kind=real_kind)  sumx00,sumx01
! Subprogram not used     real(kind=real_kind)  sumy00,sumy01
! Subprogram not used     real(kind=real_kind)  sumx10,sumx11
! Subprogram not used     real(kind=real_kind)  sumy10,sumy11
! Subprogram not used 
! Subprogram not used     !JMD ================================
! Subprogram not used     !JMD 2*np*np*np Flops 
! Subprogram not used     !JMD ================================
! Subprogram not used 
! Subprogram not used !   print *, "gradient_wk_nonstag"
! Subprogram not used     if(modulo(np,2) .eq. 0 .and. UseUnroll) then
! Subprogram not used ! this is just loop unrolling - a good compiler should do it for you jpe
! Subprogram not used 
! Subprogram not used        do j=1,np,2
! Subprogram not used           do l=1,np,2
! Subprogram not used              sumx00=0.0d0
! Subprogram not used              sumx01=0.0d0
! Subprogram not used              sumx10=0.0d0
! Subprogram not used              sumx11=0.0d0
! Subprogram not used 
! Subprogram not used              sumy00=0.0d0
! Subprogram not used              sumy01=0.0d0
! Subprogram not used              sumy10=0.0d0
! Subprogram not used              sumy11=0.0d0
! Subprogram not used 
! Subprogram not used              do i=1,np
! Subprogram not used                 sumx00 = sumx00 + deriv%Dvv_twt(i,l  )*p(i,j  )
! Subprogram not used                 sumx01 = sumx01 + deriv%Dvv_twt(i,l+1)*p(i,j  )
! Subprogram not used                 sumx10 = sumx10 + deriv%Dvv_twt(i,l  )*p(i,j+1)
! Subprogram not used                 sumx11 = sumx11 + deriv%Dvv_twt(i,l+1)*p(i,j+1)
! Subprogram not used 
! Subprogram not used                 sumy00 = sumy00 + deriv%Mvv_twt(i,l  )*p(i,j  )
! Subprogram not used                 sumy01 = sumy01 + deriv%Mvv_twt(i,l+1)*p(i,j  )
! Subprogram not used                 sumy10 = sumy10 + deriv%Mvv_twt(i,l  )*p(i,j+1)
! Subprogram not used                 sumy11 = sumy11 + deriv%Mvv_twt(i,l+1)*p(i,j+1)
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used              deriv%vvtempt(j  ,l  ,1) = sumx00
! Subprogram not used              deriv%vvtempt(j  ,l+1,1) = sumx01
! Subprogram not used              deriv%vvtempt(j+1,l  ,1) = sumx10
! Subprogram not used              deriv%vvtempt(j+1,l+1,1) = sumx11
! Subprogram not used 
! Subprogram not used              deriv%vvtempt(j  ,l  ,2) = sumy00
! Subprogram not used              deriv%vvtempt(j  ,l+1,2) = sumy01
! Subprogram not used              deriv%vvtempt(j+1,l  ,2) = sumy10
! Subprogram not used              deriv%vvtempt(j+1,l+1,2) = sumy11
! Subprogram not used 
! Subprogram not used           end do
! Subprogram not used        end do
! Subprogram not used        ! vvtempt1 = p'*Dvv_twt
! Subprogram not used        ! vvtempt2 = p'*Mvv_twt
! Subprogram not used        ! dp1 = dy*Mvv_twt*vvtempt1' = dy*Mvv_twt*(p'*Dvv_twt)' = dy*Mvv_twt*Dvv_twt'*p
! Subprogram not used        ! dp2 = dx*Dvv_twt*vvtempt2' = dx*Dvv_twt*(p'*Mvv_twt)' = dx*Dvv_twt*Mvv_twt'*p
! Subprogram not used        !     New formulation 
! Subprogram not used        ! dp1 = dy*MvvDvvt*p
! Subprogram not used        ! dp2 = dx*DvvMvvt*p
! Subprogram not used        ! MvvDvvt = Mvv_twt*Dvv_twt'
! Subprogram not used        ! DvvMvvt = Dvv_twt*Mvv_twt'
! Subprogram not used 
! Subprogram not used 
! Subprogram not used        !JMD ================================
! Subprogram not used        !JMD 2*np*np*np Flops 
! Subprogram not used        !JMD ================================
! Subprogram not used 
! Subprogram not used        do j=1,np,2
! Subprogram not used           do i=1,np,2
! Subprogram not used              sumx00=0.0d0
! Subprogram not used              sumx01=0.0d0
! Subprogram not used              sumx10=0.0d0
! Subprogram not used              sumx11=0.0d0
! Subprogram not used              
! Subprogram not used              sumy00=0.0d0
! Subprogram not used              sumy01=0.0d0
! Subprogram not used              sumy10=0.0d0
! Subprogram not used              sumy11=0.0d0
! Subprogram not used              
! Subprogram not used              do l=1,np
! Subprogram not used                 sumx00 = sumx00 +  deriv%Mvv_twt(l,j  )*deriv%vvtempt(l,i  ,1)
! Subprogram not used                 sumx01 = sumx01 +  deriv%Mvv_twt(l,j+1)*deriv%vvtempt(l,i  ,1)
! Subprogram not used                 sumx10 = sumx10 +  deriv%Mvv_twt(l,j  )*deriv%vvtempt(l,i+1,1)
! Subprogram not used                 sumx11 = sumx11 +  deriv%Mvv_twt(l,j+1)*deriv%vvtempt(l,i+1,1)
! Subprogram not used 
! Subprogram not used                 sumy00 = sumy00 +  deriv%Dvv_twt(l,j  )*deriv%vvtempt(l,i  ,2)
! Subprogram not used                 sumy01 = sumy01 +  deriv%Dvv_twt(l,j+1)*deriv%vvtempt(l,i  ,2)
! Subprogram not used                 sumy10 = sumy10 +  deriv%Dvv_twt(l,j  )*deriv%vvtempt(l,i+1,2)
! Subprogram not used                 sumy11 = sumy11 +  deriv%Dvv_twt(l,j+1)*deriv%vvtempt(l,i+1,2)
! Subprogram not used              end do
! Subprogram not used              
! Subprogram not used              dp(i  ,j  ,1) = sumx00
! Subprogram not used              dp(i  ,j+1,1) = sumx01
! Subprogram not used              dp(i+1,j  ,1) = sumx10
! Subprogram not used              dp(i+1,j+1,1) = sumx11
! Subprogram not used              
! Subprogram not used              dp(i  ,j  ,2) = sumy00
! Subprogram not used              dp(i  ,j+1,2) = sumy01
! Subprogram not used              dp(i+1,j  ,2) = sumy10
! Subprogram not used              dp(i+1,j+1,2) = sumy11
! Subprogram not used 
! Subprogram not used           end do
! Subprogram not used        end do
! Subprogram not used     else
! Subprogram not used 
! Subprogram not used        do j=1,np
! Subprogram not used           do l=1,np
! Subprogram not used              sumx00=0.0d0
! Subprogram not used 
! Subprogram not used              sumy00=0.0d0
! Subprogram not used 
! Subprogram not used              do i=1,np
! Subprogram not used                 sumx00 = sumx00 + deriv%Dvv_twt(i,l  )*p(i,j  )
! Subprogram not used 
! Subprogram not used                 sumy00 = sumy00 + deriv%Mvv_twt(i,l  )*p(i,j  )
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used              deriv%vvtempt(j  ,l  ,1) = sumx00
! Subprogram not used 
! Subprogram not used              deriv%vvtempt(j  ,l  ,2) = sumy00
! Subprogram not used 
! Subprogram not used           end do
! Subprogram not used        end do
! Subprogram not used 
! Subprogram not used        !JMD ================================
! Subprogram not used        !JMD 2*np*np*np Flops 
! Subprogram not used        !JMD ================================
! Subprogram not used 
! Subprogram not used        do j=1,np
! Subprogram not used           do i=1,np
! Subprogram not used              sumx00=0.0d0
! Subprogram not used              
! Subprogram not used              sumy00=0.0d0
! Subprogram not used              
! Subprogram not used              do l=1,np
! Subprogram not used                 sumx00 = sumx00 +  deriv%Mvv_twt(l,j  )*deriv%vvtempt(l,i  ,1)
! Subprogram not used 
! Subprogram not used                 sumy00 = sumy00 +  deriv%Dvv_twt(l,j  )*deriv%vvtempt(l,i  ,2)
! Subprogram not used              end do
! Subprogram not used              
! Subprogram not used              dp(i  ,j  ,1) = sumx00
! Subprogram not used              
! Subprogram not used              dp(i  ,j  ,2) = sumy00
! Subprogram not used 
! Subprogram not used           end do
! Subprogram not used        end do
! Subprogram not used     end if
! Subprogram not used   end function gradient_wk_nonstag

!  ================================================
!  gradient_str_stag:
! 
!  Compute the *strong* form gradient:
!  maps scalar field on the pressure grid to the
!  velocity grid
!  ================================================

! Subprogram not used   function gradient_str_stag(p,deriv) result(dp)
! Subprogram not used 
! Subprogram not used     type (derivative_stag_t)         :: deriv
! Subprogram not used     real(kind=real_kind), intent(in) :: p(np,np)
! Subprogram not used 
! Subprogram not used     real(kind=real_kind)             :: dp(np,np,2)
! Subprogram not used 
! Subprogram not used     ! Local
! Subprogram not used       
! Subprogram not used     integer i
! Subprogram not used     integer j
! Subprogram not used     integer l
! Subprogram not used 
! Subprogram not used     logical, parameter :: UseUnroll=.TRUE.
! Subprogram not used 
! Subprogram not used     real(kind=real_kind)  sumx00,sumx01
! Subprogram not used     real(kind=real_kind)  sumy00,sumy01
! Subprogram not used     real(kind=real_kind)  sumx10,sumx11
! Subprogram not used     real(kind=real_kind)  sumy10,sumy11
! Subprogram not used 
! Subprogram not used if(MODULO(np,2) == 0 .and. UseUnroll) then 
! Subprogram not used     do j=1,np,2
! Subprogram not used        do l=1,np,2
! Subprogram not used           sumx00=0.0d0
! Subprogram not used           sumx01=0.0d0
! Subprogram not used           sumx10=0.0d0
! Subprogram not used           sumx11=0.0d0
! Subprogram not used 
! Subprogram not used           sumy00=0.0d0
! Subprogram not used           sumy01=0.0d0
! Subprogram not used           sumy10=0.0d0
! Subprogram not used           sumy11=0.0d0
! Subprogram not used 
! Subprogram not used           do i=1,np
! Subprogram not used              sumx00 = sumx00 + deriv%Dpv(i,l  )*p(i,j  )
! Subprogram not used              sumx01 = sumx01 + deriv%Dpv(i,l+1)*p(i,j  )
! Subprogram not used              sumx10 = sumx10 + deriv%Dpv(i,l  )*p(i,j+1)
! Subprogram not used              sumx11 = sumx11 + deriv%Dpv(i,l+1)*p(i,j+1)
! Subprogram not used 
! Subprogram not used              sumy00 = sumy00 + deriv%M_t(i,l  )*p(i,j  )
! Subprogram not used              sumy01 = sumy01 + deriv%M_t(i,l+1)*p(i,j  )
! Subprogram not used              sumy10 = sumy10 + deriv%M_t(i,l  )*p(i,j+1)
! Subprogram not used              sumy11 = sumy11 + deriv%M_t(i,l+1)*p(i,j+1)
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used           deriv%vtempt(j  ,l  ,1) = sumx00
! Subprogram not used           deriv%vtempt(j  ,l+1,1) = sumx01
! Subprogram not used           deriv%vtempt(j+1,l  ,1) = sumx10
! Subprogram not used           deriv%vtempt(j+1,l+1,1) = sumx11
! Subprogram not used 
! Subprogram not used           deriv%vtempt(j  ,l  ,2) = sumy00
! Subprogram not used           deriv%vtempt(j  ,l+1,2) = sumy01
! Subprogram not used           deriv%vtempt(j+1,l  ,2) = sumy10
! Subprogram not used           deriv%vtempt(j+1,l+1,2) = sumy11
! Subprogram not used 
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     do j=1,np,2
! Subprogram not used        do i=1,np,2
! Subprogram not used           sumx00=0.0d0
! Subprogram not used           sumx01=0.0d0
! Subprogram not used           sumx10=0.0d0
! Subprogram not used           sumx11=0.0d0
! Subprogram not used 
! Subprogram not used           sumy00=0.0d0
! Subprogram not used           sumy01=0.0d0
! Subprogram not used           sumy10=0.0d0
! Subprogram not used           sumy11=0.0d0
! Subprogram not used 
! Subprogram not used           do l=1,np
! Subprogram not used              sumx00 = sumx00 +  deriv%M_t(l,j  )*deriv%vtempt(l,i  ,1)
! Subprogram not used              sumx01 = sumx01 +  deriv%M_t(l,j+1)*deriv%vtempt(l,i  ,1)
! Subprogram not used              sumx10 = sumx10 +  deriv%M_t(l,j  )*deriv%vtempt(l,i+1,1)
! Subprogram not used              sumx11 = sumx11 +  deriv%M_t(l,j+1)*deriv%vtempt(l,i+1,1)
! Subprogram not used 
! Subprogram not used              sumy00 = sumy00 +  deriv%Dpv(l,j  )*deriv%vtempt(l,i  ,2)
! Subprogram not used              sumy01 = sumy01 +  deriv%Dpv(l,j+1)*deriv%vtempt(l,i  ,2)
! Subprogram not used              sumy10 = sumy10 +  deriv%Dpv(l,j  )*deriv%vtempt(l,i+1,2)
! Subprogram not used              sumy11 = sumy11 +  deriv%Dpv(l,j+1)*deriv%vtempt(l,i+1,2)
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used           dp(i  ,j  ,1) = sumx00
! Subprogram not used           dp(i  ,j+1,1) = sumx01
! Subprogram not used           dp(i+1,j  ,1) = sumx10
! Subprogram not used           dp(i+1,j+1,1) = sumx11
! Subprogram not used 
! Subprogram not used           dp(i  ,j  ,2) = sumy00
! Subprogram not used           dp(i  ,j+1,2) = sumy01
! Subprogram not used           dp(i+1,j  ,2) = sumy10
! Subprogram not used           dp(i+1,j+1,2) = sumy11
! Subprogram not used 
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used else
! Subprogram not used     do j=1,np
! Subprogram not used        do l=1,np
! Subprogram not used           sumx00=0.0d0
! Subprogram not used           sumy00=0.0d0
! Subprogram not used           do i=1,np
! Subprogram not used              sumx00 = sumx00 + deriv%Dpv(i,l  )*p(i,j  )
! Subprogram not used              sumy00 = sumy00 + deriv%M_t(i,l  )*p(i,j  )
! Subprogram not used    	   enddo
! Subprogram not used 	   deriv%vtempt(j  ,l  ,1) = sumx00
! Subprogram not used    	   deriv%vtempt(j  ,l  ,2) = sumy00
! Subprogram not used 	enddo
! Subprogram not used     enddo
! Subprogram not used     do j=1,np
! Subprogram not used        do i=1,np
! Subprogram not used           sumx00=0.0d0
! Subprogram not used           sumy00=0.0d0
! Subprogram not used           do l=1,np
! Subprogram not used              sumx00 = sumx00 +  deriv%M_t(l,j  )*deriv%vtempt(l,i  ,1)
! Subprogram not used              sumy00 = sumy00 +  deriv%Dpv(l,j  )*deriv%vtempt(l,i  ,2)
! Subprogram not used 	  enddo
! Subprogram not used           dp(i  ,j  ,1) = sumx00
! Subprogram not used 	  dp(i  ,j  ,2) = sumy00
! Subprogram not used        enddo
! Subprogram not used     enddo
! Subprogram not used endif
! Subprogram not used 
! Subprogram not used   end function gradient_str_stag

!  ================================================
!  gradient_str_nonstag:
!
!  Compute the *strong* gradient on the velocity grid
!  of a scalar field on the velocity grid
!  ================================================

! Subprogram not used   function gradient_str_nonstag(s,deriv) result(ds)
! Subprogram not used 
! Subprogram not used     type (derivative_t)              :: deriv
! Subprogram not used     real(kind=real_kind), intent(in) :: s(np,np)
! Subprogram not used 
! Subprogram not used     real(kind=real_kind) :: ds(np,np,2)
! Subprogram not used 
! Subprogram not used     integer i
! Subprogram not used     integer j
! Subprogram not used     integer l
! Subprogram not used     logical, parameter :: UseUnroll = .TRUE.
! Subprogram not used 
! Subprogram not used     real(kind=real_kind) ::  dsdx00,dsdx01
! Subprogram not used     real(kind=real_kind) ::  dsdx10,dsdx11
! Subprogram not used 
! Subprogram not used     real(kind=real_kind) ::  dsdy00,dsdy01
! Subprogram not used     real(kind=real_kind) ::  dsdy10,dsdy11
! Subprogram not used     if(modulo(np,2) .eq. 0 .and. UseUnroll) then
! Subprogram not used        do j=1,np,2
! Subprogram not used           do l=1,np,2
! Subprogram not used              dsdx00=0.0d0
! Subprogram not used              dsdx01=0.0d0
! Subprogram not used              dsdx10=0.0d0
! Subprogram not used              dsdx11=0.0d0
! Subprogram not used 
! Subprogram not used              dsdy00=0.0d0
! Subprogram not used              dsdy01=0.0d0
! Subprogram not used              dsdy10=0.0d0
! Subprogram not used              dsdy11=0.0d0
! Subprogram not used 
! Subprogram not used              do i=1,np
! Subprogram not used                 dsdx00 = dsdx00 + deriv%Dvv(i,l  )*s(i,j  )
! Subprogram not used                 dsdx01 = dsdx01 + deriv%Dvv(i,l+1)*s(i,j  )
! Subprogram not used                 dsdx10 = dsdx10 + deriv%Dvv(i,l  )*s(i,j+1)
! Subprogram not used                 dsdx11 = dsdx11 + deriv%Dvv(i,l+1)*s(i,j+1)
! Subprogram not used 
! Subprogram not used                 dsdy00 = dsdy00 + deriv%Dvv(i,l  )*s(j  ,i)
! Subprogram not used                 dsdy01 = dsdy01 + deriv%Dvv(i,l+1)*s(j  ,i)
! Subprogram not used                 dsdy10 = dsdy10 + deriv%Dvv(i,l  )*s(j+1,i)
! Subprogram not used                 dsdy11 = dsdy11 + deriv%Dvv(i,l+1)*s(j+1,i)
! Subprogram not used              end do
! Subprogram not used              ds(l  ,j  ,1) = dsdx00
! Subprogram not used              ds(l+1,j  ,1) = dsdx01
! Subprogram not used              ds(l  ,j+1,1) = dsdx10
! Subprogram not used              ds(l+1,j+1,1) = dsdx11
! Subprogram not used 
! Subprogram not used              ds(j  ,l  ,2) = dsdy00
! Subprogram not used              ds(j  ,l+1,2) = dsdy01
! Subprogram not used              ds(j+1,l  ,2) = dsdy10
! Subprogram not used              ds(j+1,l+1,2) = dsdy11
! Subprogram not used 
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used        end do
! Subprogram not used     else
! Subprogram not used        do j=1,np
! Subprogram not used           do l=1,np
! Subprogram not used              dsdx00=0.0d0
! Subprogram not used 
! Subprogram not used              dsdy00=0.0d0
! Subprogram not used 
! Subprogram not used              do i=1,np
! Subprogram not used                 dsdx00 = dsdx00 + deriv%Dvv(i,l  )*s(i,j  )
! Subprogram not used 
! Subprogram not used                 dsdy00 = dsdy00 + deriv%Dvv(i,l  )*s(j  ,i)
! Subprogram not used              end do
! Subprogram not used              ds(l  ,j  ,1) = dsdx00
! Subprogram not used 
! Subprogram not used              ds(j  ,l  ,2) = dsdy00
! Subprogram not used 
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used        end do
! Subprogram not used     end if
! Subprogram not used   end function gradient_str_nonstag

!  ================================================
!  vorticity:
!
!  Compute the vorticity of the velocity field on the
!  velocity grid
!  ================================================

! Subprogram not used   function vorticity(v,deriv) result(vort)
! Subprogram not used 
! Subprogram not used     type (derivative_t)              :: deriv
! Subprogram not used     real(kind=real_kind), intent(in) :: v(np,np,2)
! Subprogram not used 
! Subprogram not used     real(kind=real_kind) :: vort(np,np)
! Subprogram not used 
! Subprogram not used     integer i
! Subprogram not used     integer j
! Subprogram not used     integer l
! Subprogram not used     
! Subprogram not used     logical, parameter :: UseUnroll = .TRUE.
! Subprogram not used 
! Subprogram not used     real(kind=real_kind) ::  dvdx00,dvdx01
! Subprogram not used     real(kind=real_kind) ::  dvdx10,dvdx11
! Subprogram not used 
! Subprogram not used     real(kind=real_kind) ::  dudy00,dudy01
! Subprogram not used     real(kind=real_kind) ::  dudy10,dudy11
! Subprogram not used 
! Subprogram not used if(MODULO(np,2) == 0 .and. UseUnroll) then 
! Subprogram not used     do j=1,np,2
! Subprogram not used        do l=1,np,2
! Subprogram not used 
! Subprogram not used           dudy00=0.0d0
! Subprogram not used           dudy01=0.0d0
! Subprogram not used           dudy10=0.0d0
! Subprogram not used           dudy11=0.0d0
! Subprogram not used 
! Subprogram not used           dvdx00=0.0d0
! Subprogram not used           dvdx01=0.0d0
! Subprogram not used           dvdx10=0.0d0
! Subprogram not used           dvdx11=0.0d0
! Subprogram not used 
! Subprogram not used           do i=1,np
! Subprogram not used 
! Subprogram not used              dvdx00 = dvdx00 + deriv%Dvv(i,l  )*v(i,j  ,2)
! Subprogram not used              dvdx01 = dvdx01 + deriv%Dvv(i,l+1)*v(i,j  ,2)
! Subprogram not used              dvdx10 = dvdx10 + deriv%Dvv(i,l  )*v(i,j+1,2)
! Subprogram not used              dvdx11 = dvdx11 + deriv%Dvv(i,l+1)*v(i,j+1,2)
! Subprogram not used 
! Subprogram not used              dudy00 = dudy00 + deriv%Dvv(i,l  )*v(j  ,i,1)
! Subprogram not used              dudy01 = dudy01 + deriv%Dvv(i,l+1)*v(j  ,i,1)
! Subprogram not used              dudy10 = dudy10 + deriv%Dvv(i,l  )*v(j+1,i,1)
! Subprogram not used              dudy11 = dudy11 + deriv%Dvv(i,l+1)*v(j+1,i,1)
! Subprogram not used 
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used           vort(l  ,j  ) = dvdx00
! Subprogram not used           vort(l+1,j  ) = dvdx01
! Subprogram not used           vort(l  ,j+1) = dvdx10
! Subprogram not used           vort(l+1,j+1) = dvdx11
! Subprogram not used 
! Subprogram not used           deriv%vvtemp(j  ,l  ) = dudy00
! Subprogram not used           deriv%vvtemp(j  ,l+1) = dudy01
! Subprogram not used           deriv%vvtemp(j+1,l  ) = dudy10
! Subprogram not used           deriv%vvtemp(j+1,l+1) = dudy11
! Subprogram not used 
! Subprogram not used         end do
! Subprogram not used     end do
! Subprogram not used else
! Subprogram not used     do j=1,np
! Subprogram not used        do l=1,np
! Subprogram not used 
! Subprogram not used           dudy00=0.0d0
! Subprogram not used 	  dvdx00=0.0d0
! Subprogram not used 
! Subprogram not used           do i=1,np
! Subprogram not used              dvdx00 = dvdx00 + deriv%Dvv(i,l  )*v(i,j  ,2)
! Subprogram not used              dudy00 = dudy00 + deriv%Dvv(i,l  )*v(j  ,i,1)
! Subprogram not used 	  enddo
! Subprogram not used  
! Subprogram not used 	  vort(l  ,j  ) = dvdx00
! Subprogram not used 	  deriv%vvtemp(j  ,l  ) = dudy00
! Subprogram not used 	enddo
! Subprogram not used      enddo
! Subprogram not used 
! Subprogram not used endif
! Subprogram not used 
! Subprogram not used     do j=1,np
! Subprogram not used        do i=1,np
! Subprogram not used           vort(i,j)=vort(i,j)-deriv%vvtemp(i,j)
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used   end function vorticity




!  ================================================
!  interpolate_gll2fvm_points:
!
!  shape funtion interpolation from data on GLL grid to cellcenters on physics grid
!  Author: Christoph Erath
!  ================================================
! Subprogram not used   function interpolate_gll2fvm_points(v,deriv) result(p)
! Subprogram not used 
! Subprogram not used     real(kind=real_kind), intent(in) :: v(np,np)
! Subprogram not used     type (derivative_t)         :: deriv
! Subprogram not used     real(kind=real_kind) :: p(nc,nc)
! Subprogram not used 
! Subprogram not used     ! Local
! Subprogram not used     integer i
! Subprogram not used     integer j
! Subprogram not used     integer l
! Subprogram not used 
! Subprogram not used     real(kind=real_kind)  sumx00,sumx01
! Subprogram not used     real(kind=real_kind)  sumx10,sumx11
! Subprogram not used     real(kind=real_kind)  vtemp(np,nc)
! Subprogram not used 
! Subprogram not used     do j=1,np
! Subprogram not used        do l=1,nc
! Subprogram not used           sumx00=0.0d0
! Subprogram not used           do i=1,np
! Subprogram not used              sumx00 = sumx00 + deriv%Cfvm(i,l  )*v(i,j  )
! Subprogram not used           enddo
! Subprogram not used           vtemp(j  ,l) = sumx00
! Subprogram not used         enddo
! Subprogram not used     enddo
! Subprogram not used     do j=1,nc
! Subprogram not used        do i=1,nc
! Subprogram not used           sumx00=0.0d0
! Subprogram not used           do l=1,np
! Subprogram not used              sumx00 = sumx00 + deriv%Cfvm(l,j  )*vtemp(l,i)
! Subprogram not used           enddo
! Subprogram not used           p(i  ,j  ) = sumx00
! Subprogram not used        enddo
! Subprogram not used     enddo
! Subprogram not used   end function interpolate_gll2fvm_points
  !  ================================================
  !  interpolate_gll2spelt_points:
  !
  !  shape function interpolation from data on GLL grid the spelt grid
  !  Author: Christoph Erath
  !  ================================================
! Subprogram not used   function interpolate_gll2spelt_points(v,deriv) result(p)
! Subprogram not used     real(kind=real_kind), intent(in) :: v(np,np)
! Subprogram not used     type (derivative_t)         :: deriv
! Subprogram not used     real(kind=real_kind) :: p(nep,nep)
! Subprogram not used 
! Subprogram not used     ! Local
! Subprogram not used     integer i,j,l
! Subprogram not used 
! Subprogram not used     real(kind=real_kind)  sumx00,sumx01
! Subprogram not used     real(kind=real_kind)  sumx10,sumx11
! Subprogram not used     real(kind=real_kind)  vtemp(np,nep)
! Subprogram not used 
! Subprogram not used     do j=1,np
! Subprogram not used        do l=1,nep
! Subprogram not used           sumx00=0.0d0
! Subprogram not used           do i=1,np
! Subprogram not used              sumx00 = sumx00 + deriv%Sfvm(i,l  )*v(i,j  )
! Subprogram not used           enddo
! Subprogram not used           vtemp(j  ,l) = sumx00
! Subprogram not used         enddo
! Subprogram not used     enddo
! Subprogram not used     do j=1,nep
! Subprogram not used        do i=1,nep
! Subprogram not used           sumx00=0.0d0
! Subprogram not used           do l=1,np
! Subprogram not used              sumx00 = sumx00 + deriv%Sfvm(l,j  )*vtemp(l,i)
! Subprogram not used           enddo
! Subprogram not used           p(i  ,j  ) = sumx00
! Subprogram not used        enddo
! Subprogram not used     enddo
! Subprogram not used   end function interpolate_gll2spelt_points

!  ================================================
!  interpolate_gll2fvm_corners:
!
!  shape funtion interpolation from data on GLL grid to physics grid
!
!  ================================================
! Subprogram not used   function interpolate_gll2fvm_corners(v,deriv) result(p)
! Subprogram not used 
! Subprogram not used     real(kind=real_kind), intent(in) :: v(np,np)
! Subprogram not used     type (derivative_t)         :: deriv
! Subprogram not used     real(kind=real_kind) :: p(nc+1,nc+1)
! Subprogram not used 
! Subprogram not used     ! Local
! Subprogram not used     integer i
! Subprogram not used     integer j
! Subprogram not used     integer l
! Subprogram not used 
! Subprogram not used     real(kind=real_kind)  sumx00,sumx01
! Subprogram not used     real(kind=real_kind)  sumx10,sumx11
! Subprogram not used     real(kind=real_kind)  vtemp(np,nc+1)
! Subprogram not used 
! Subprogram not used     do j=1,np
! Subprogram not used        do l=1,nc+1
! Subprogram not used           sumx00=0.0d0
! Subprogram not used           do i=1,np
! Subprogram not used              sumx00 = sumx00 + deriv%Mfvm(i,l  )*v(i,j  )
! Subprogram not used           enddo
! Subprogram not used           vtemp(j  ,l) = sumx00
! Subprogram not used         enddo
! Subprogram not used     enddo
! Subprogram not used     do j=1,nc+1
! Subprogram not used        do i=1,nc+1
! Subprogram not used           sumx00=0.0d0
! Subprogram not used           do l=1,np
! Subprogram not used              sumx00 = sumx00 + deriv%Mfvm(l,j  )*vtemp(l,i)
! Subprogram not used           enddo
! Subprogram not used           p(i  ,j  ) = sumx00
! Subprogram not used        enddo
! Subprogram not used     enddo
! Subprogram not used   end function interpolate_gll2fvm_corners



!  ================================================
!  remap_phys2gll:
!
!  interpolate to an equally spaced (in reference element coordinate system)
!  "physics" grid to the GLL grid
!
!  1st order, monotone, conservative
!  ================================================
! Subprogram not used   function remap_phys2gll(pin,nphys) result(pout)
! Subprogram not used     integer :: nphys
! Subprogram not used     real(kind=real_kind), intent(in) :: pin(nphys*nphys)
! Subprogram not used     real(kind=real_kind) :: pout(np,np)
! Subprogram not used     
! Subprogram not used     ! Local
! Subprogram not used     integer, save  :: nphys_init=0
! Subprogram not used     integer, save  :: nintersect
! Subprogram not used     real(kind=real_kind),save,pointer :: acell(:)  ! arrivial cell index of i'th intersection
! Subprogram not used     real(kind=real_kind),save,pointer :: dcell(:)  ! departure cell index of i'th intersection
! Subprogram not used     real(kind=real_kind),save,pointer :: delta(:)  ! length of i'th intersection
! Subprogram not used     real(kind=real_kind),save,pointer :: delta_a(:)  ! length of arrival cells
! Subprogram not used     integer in_i,in_j,ia,ja,id,jd,count,i,j
! Subprogram not used     logical :: found
! Subprogram not used 
! Subprogram not used     real(kind=real_kind) :: tol=1e-13
! Subprogram not used     real(kind=real_kind) :: weight,x1,x2,dx
! Subprogram not used     real(kind=longdouble_kind) :: gll_edges(np+1),phys_edges(nphys+1)
! Subprogram not used     type(quadrature_t) :: gll_pts
! Subprogram not used     ! setup (most be done on masterthread only) since all data is static
! Subprogram not used !OMP MASTER
! Subprogram not used     if (nphys_init/=nphys) then
! Subprogram not used        nphys_init=nphys
! Subprogram not used        ! find number of intersections
! Subprogram not used        nintersect = np+nphys-1  ! max number of possible intersections
! Subprogram not used        allocate(acell(nintersect))
! Subprogram not used        allocate(dcell(nintersect))
! Subprogram not used        allocate(delta(nintersect))
! Subprogram not used        allocate(delta_a(np))
! Subprogram not used 
! Subprogram not used        ! compute phys grid cell edges on [-1,1]
! Subprogram not used        do i=1,nphys+1
! Subprogram not used           dx = 2d0/nphys
! Subprogram not used           phys_edges(i)=-1 + (i-1)*dx
! Subprogram not used        enddo
! Subprogram not used 
! Subprogram not used        ! compute GLL cell edges on [-1,1]
! Subprogram not used        gll_pts = gausslobatto(np)
! Subprogram not used        gll_edges(1)=-1
! Subprogram not used        do i=2,np
! Subprogram not used           gll_edges(i) = gll_edges(i-1) + gll_pts%weights(i-1)
! Subprogram not used        enddo
! Subprogram not used        gll_edges(np+1)=1
! Subprogram not used        delta_a=gll_pts%weights
! Subprogram not used        deallocate(gll_pts%points)
! Subprogram not used        deallocate(gll_pts%weights)
! Subprogram not used 
! Subprogram not used        count=0
! Subprogram not used        x1=-1
! Subprogram not used        do while ( abs(x1-1) > tol )
! Subprogram not used           ! find point x2 closet to x1 and x2>x1:
! Subprogram not used           x2=1.1
! Subprogram not used           do ia=2,np+1
! Subprogram not used              if (gll_edges(ia)>x1) then
! Subprogram not used                 if ( ( gll_edges(ia)-x1) < (x2-x1) ) then
! Subprogram not used                    x2=gll_edges(ia)
! Subprogram not used                 endif
! Subprogram not used              endif
! Subprogram not used           enddo
! Subprogram not used           do id=2,nphys+1
! Subprogram not used              if (phys_edges(id)>x1) then
! Subprogram not used                 if ( ( phys_edges(id)-x1) < (x2-x1) ) then
! Subprogram not used                    x2=phys_edges(id)
! Subprogram not used                 endif
! Subprogram not used              endif
! Subprogram not used           enddo
! Subprogram not used           if (x2>1+tol) call abortmp('ERROR: did not find next intersection point')
! Subprogram not used              
! Subprogram not used           count=count+1
! Subprogram not used           delta(count)=x2-x1
! Subprogram not used           
! Subprogram not used           found=.false.
! Subprogram not used           do ia=1,np
! Subprogram not used              if (gll_edges(ia) <= x1+tol  .and.  x2-tol <= gll_edges(ia+1)) then
! Subprogram not used                 found=.true.
! Subprogram not used                 acell(count)=ia
! Subprogram not used              endif
! Subprogram not used           enddo
! Subprogram not used           if (.not. found) call abortmp('ERROR: interval search problem')
! Subprogram not used        
! Subprogram not used           found=.false.
! Subprogram not used           do id=1,nphys
! Subprogram not used              if (phys_edges(id) <= x1+tol .and.  x2-tol <= phys_edges(id+1)) then
! Subprogram not used                 found=.true.
! Subprogram not used                 dcell(count)=id
! Subprogram not used              endif
! Subprogram not used           enddo
! Subprogram not used           if (.not. found) call abortmp('ERROR: interval search problem')
! Subprogram not used           x1=x2
! Subprogram not used        enddo
! Subprogram not used        if (count>nintersect) call abortmp('ERROR: nintersect was too small')
! Subprogram not used        nintersect=count
! Subprogram not used     endif
! Subprogram not used     !OMP END MASTER
! Subprogram not used     !OMP BARRIER
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     pout=0
! Subprogram not used     do in_i = 1,nintersect
! Subprogram not used        do in_j = 1,nintersect
! Subprogram not used           ia = acell(in_i)
! Subprogram not used           ja = acell(in_j)
! Subprogram not used           id = dcell(in_i)
! Subprogram not used           jd = dcell(in_j)
! Subprogram not used           ! mass in intersection region:  value*area_intersect
! Subprogram not used           ! value_arrival = value*area_intersect/area_arrival
! Subprogram not used           weight = (  delta(in_i)*delta(in_j) ) / ( delta_a(ia)*delta_a(ja))
! Subprogram not used           ! accumulate contribution from each intersection region:
! Subprogram not used           pout(ia,ja) = pout(ia,ja) + weight*pin(id+(jd-1)*nphys)
! Subprogram not used        enddo
! Subprogram not used     enddo
! Subprogram not used     
! Subprogram not used     end function remap_phys2gll
    
!----------------------------------------------------------------


  function gradient_sphere(s,deriv,Dinv) result(ds)
!
!   input s:  scalar
!   output  ds: spherical gradient of s, lat-lon coordinates
!

    type (derivative_t)              :: deriv
    real(kind=real_kind), intent(in), dimension(2,2,np,np) :: Dinv
    real(kind=real_kind), intent(in) :: s(np,np)

    real(kind=real_kind) :: ds(np,np,2)

    integer i
    integer j
    integer l

    real(kind=real_kind) ::  dsdx00
    real(kind=real_kind) ::  dsdy00
    real(kind=real_kind) ::  v1(np,np),v2(np,np)

    do j=1,np
       do l=1,np
          dsdx00=0.0d0
          dsdy00=0.0d0
          do i=1,np
             dsdx00 = dsdx00 + deriv%Dvv(i,l  )*s(i,j  )
             dsdy00 = dsdy00 + deriv%Dvv(i,l  )*s(j  ,i)
          end do
          v1(l  ,j  ) = dsdx00*rrearth
          v2(j  ,l  ) = dsdy00*rrearth
       end do
    end do
    ! convert covarient to latlon
    do j=1,np
       do i=1,np
          ds(i,j,1)=Dinv(1,1,i,j)*v1(i,j) + Dinv(2,1,i,j)*v2(i,j)
          ds(i,j,2)=Dinv(1,2,i,j)*v1(i,j) + Dinv(2,2,i,j)*v2(i,j)
       enddo
    enddo

    end function gradient_sphere




  function curl_sphere_wk_testcov(s,deriv,elem) result(ds)
!
!   integrated-by-parts gradient, w.r.t. COVARIANT test functions
!   input s:  scalar  (assumed to be s*khat)
!   output  ds: weak curl, lat/lon coordinates
!   
! starting with: 
!   PHIcov1 = (PHI,0)  covariant vector 
!   PHIcov2 = (0,PHI)  covariant vector 
!
!   ds1 = integral[ PHIcov1 dot curl(s*khat) ] 
!   ds2 = integral[ PHIcov2 dot curl(s*khat) ] 
! integrate by parts: 
!   ds1 = integral[ vor(PHIcov1) * s ]       
!   ds2 = integral[ vor(PHIcov1) * s ]
!
!     PHIcov1 = (PHI^mn,0)   
!     PHIcov2 = (0,PHI^mn)
!  vorticity() acts on covariant vectors:
!   ds1 = sum wij g  s_ij 1/g (  (PHIcov1_2)_x  - (PHIcov1_1)_y ) 
!       = -sum wij s_ij  d/dy (PHI^mn )
! for d/dy component, only sum over i=m
!       = -sum  w_mj s_mj   d( PHI^n)(j)
!           j
!
!   ds2 = sum wij g  s_ij 1/g (  (PHIcov2_2)_x  - (PHIcov2_1)_y ) 
!       = +sum wij s_ij  d/dx (PHI^mn )
! for d/dx component, only sum over j=n
!       = +sum  w_in s_in  d( PHI^m)(i)
!           i
!
    type (derivative_t)              :: deriv
    type (element_t)                 :: elem
    real(kind=real_kind), intent(in) :: s(np,np)

    real(kind=real_kind) :: ds(np,np,2)

    integer i,j,l,m,n
    real(kind=real_kind) ::  dscontra(np,np,2)

    dscontra=0
    do n=1,np
       do m=1,np
          do j=1,np
             ! phi(n)_y  sum over second index, 1st index fixed at m
             dscontra(m,n,1)=dscontra(m,n,1)-(elem%mp(m,j)*s(m,j)*deriv%Dvv(n,j) )*rrearth
             ! phi(m)_x  sum over first index, second index fixed at n
             dscontra(m,n,2)=dscontra(m,n,2)+(elem%mp(j,n)*s(j,n)*deriv%Dvv(m,j) )*rrearth
          enddo
       enddo
    enddo

    ! convert contra -> latlon 
    do j=1,np
       do i=1,np
          ds(i,j,1)=(elem%D(1,1,i,j)*dscontra(i,j,1) + elem%D(1,2,i,j)*dscontra(i,j,2))
          ds(i,j,2)=(elem%D(2,1,i,j)*dscontra(i,j,1) + elem%D(2,2,i,j)*dscontra(i,j,2))
       enddo
    enddo
    end function curl_sphere_wk_testcov


  function gradient_sphere_wk_testcov(s,deriv,elem) result(ds)
!
!   integrated-by-parts gradient, w.r.t. COVARIANT test functions
!   input s:  scalar
!   output  ds: weak gradient, lat/lon coordinates
!   ds = - integral[ div(PHIcov) s ]
!
!     PHIcov1 = (PHI^mn,0)   
!     PHIcov2 = (0,PHI^mn)
!   div() acts on contra components, so convert test function to contra: 
!     PHIcontra1 =  metinv PHIcov1  = (a^mn,b^mn)*PHI^mn   
!                                     a = metinv(1,1)  b=metinv(2,1)
!
!   ds1 = sum wij g  s_ij 1/g ( g a PHI^mn)_x  + ( g b PHI^mn)_y ) 
!       = sum  wij s_ij  ag(m,n)  d/dx( PHI^mn ) + bg(m,n) d/dy( PHI^mn)
!          i,j 
! for d/dx component, only sum over j=n
!       = sum  w_in s_in  ag(m,n)  d( PHI^m)(i)
!          i
! for d/dy component, only sum over i=m
!       = sum  w_mj s_mj  bg(m,n)  d( PHI^n)(j)
!          j
!  
!
! This formula is identical to gradient_sphere_wk_testcontra, except that
!    g(m,n) is replaced by a(m,n)*g(m,n)   
!  and we have two terms for each componet of ds 
!
!
    type (derivative_t)              :: deriv
    type (element_t)                 :: elem
    real(kind=real_kind), intent(in) :: s(np,np)

    real(kind=real_kind) :: ds(np,np,2)

    integer i,j,l,m,n
    real(kind=real_kind) ::  dscontra(np,np,2)


    dscontra=0
    do n=1,np
       do m=1,np
          do j=1,np
             dscontra(m,n,1)=dscontra(m,n,1)-(&
                  (elem%mp(j,n)*elem%metinv(1,1,m,n)*elem%metdet(m,n)*s(j,n)*deriv%Dvv(m,j) ) +&
                  (elem%mp(m,j)*elem%metinv(2,1,m,n)*elem%metdet(m,n)*s(m,j)*deriv%Dvv(n,j) ) &
                  ) *rrearth

             dscontra(m,n,2)=dscontra(m,n,2)-(&
                  (elem%mp(j,n)*elem%metinv(1,2,m,n)*elem%metdet(m,n)*s(j,n)*deriv%Dvv(m,j) ) +&
                  (elem%mp(m,j)*elem%metinv(2,2,m,n)*elem%metdet(m,n)*s(m,j)*deriv%Dvv(n,j) ) &
                  ) *rrearth
          enddo
       enddo
    enddo
    ! convert contra -> latlon 
    do j=1,np
       do i=1,np
          ds(i,j,1)=(elem%D(1,1,i,j)*dscontra(i,j,1) + elem%D(1,2,i,j)*dscontra(i,j,2))
          ds(i,j,2)=(elem%D(2,1,i,j)*dscontra(i,j,1) + elem%D(2,2,i,j)*dscontra(i,j,2))
       enddo
    enddo

    end function gradient_sphere_wk_testcov


! Subprogram not used   function gradient_sphere_wk_testcontra(s,deriv,elem) result(ds)
! Subprogram not used !
! Subprogram not used !   integrated-by-parts gradient, w.r.t. CONTRA test functions
! Subprogram not used !   input s:  scalar
! Subprogram not used !   output  ds: weak gradient, lat/lon coordinates
! Subprogram not used !
! Subprogram not used !   integral[ div(phivec) s ] = sum  spheremp()* divergence_sphere(phivec) *s
! Subprogram not used !   ds1 = above formual with phivec=(PHI,0) in CONTRA coordinates
! Subprogram not used !   ds2 = above formual with phivec=(0,PHI) in CONTRA coordinates
! Subprogram not used !   
! Subprogram not used ! PHI = (phi,0)
! Subprogram not used !   s1 =  sum w_ij s_ij g_ij 1/g_ij ( g_ij PHI^mn )x  
! Subprogram not used !      =  sum w_ij s_ij g_mn dx(PHI^mn)_ij 
! Subprogram not used !         ij
! Subprogram not used ! because x derivative is zero for j<>n, only have to sum over j=n
! Subprogram not used !   s1(m,n)  =  sum w_i,n g_mn dx(PHI^m)_i,n s_i,n
! Subprogram not used !                i
! Subprogram not used !
! Subprogram not used     type (derivative_t)              :: deriv
! Subprogram not used     type (element_t)              :: elem
! Subprogram not used     real(kind=real_kind), intent(in) :: s(np,np)
! Subprogram not used 
! Subprogram not used     real(kind=real_kind) :: ds(np,np,2)
! Subprogram not used 
! Subprogram not used     integer i,j,l,m,n
! Subprogram not used     real(kind=real_kind) ::  dscov(np,np,2)
! Subprogram not used 
! Subprogram not used     ! debug: 
! Subprogram not used     real(kind=real_kind) ::  vcontra(np,np,2)
! Subprogram not used     real(kind=real_kind) ::  v(np,np,2)
! Subprogram not used     real(kind=real_kind) ::  div(np,np)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     dscov=0
! Subprogram not used     do n=1,np
! Subprogram not used        do m=1,np
! Subprogram not used           do j=1,np
! Subprogram not used              ! phi(m)_x  sum over first index, second index fixed at n
! Subprogram not used              dscov(m,n,1)=dscov(m,n,1)-(elem%mp(j,n)*elem%metdet(m,n)*s(j,n)*deriv%Dvv(m,j) )*rrearth
! Subprogram not used              ! phi(n)_y  sum over second index, 1st index fixed at m
! Subprogram not used              dscov(m,n,2)=dscov(m,n,2)-(elem%mp(m,j)*elem%metdet(m,n)*s(m,j)*deriv%Dvv(n,j) )*rrearth
! Subprogram not used           enddo
! Subprogram not used        enddo
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     ! convert covariant -> latlon 
! Subprogram not used     ds(:,:,1)=elem%Dinv(1,1,:,:)*dscov(:,:,1) + elem%Dinv(2,1,:,:)*dscov(:,:,2)
! Subprogram not used     ds(:,:,2)=elem%Dinv(1,2,:,:)*dscov(:,:,1) + elem%Dinv(2,2,:,:)*dscov(:,:,2)
! Subprogram not used 
! Subprogram not used     end function gradient_sphere_wk_testcontra




! Subprogram not used   function ugradv_sphere(u,v,deriv,elem) result(ugradv)
! Subprogram not used !
! Subprogram not used !   input:  vectors u and v  (latlon coordinates)
! Subprogram not used !   output: vector  [ u dot grad ] v  (latlon coordinates)
! Subprogram not used !
! Subprogram not used     type (derivative_t)              :: deriv
! Subprogram not used     type (element_t)                 :: elem
! Subprogram not used     real(kind=real_kind), intent(in) :: u(np,np,2)
! Subprogram not used     real(kind=real_kind), intent(in) :: v(np,np,2)
! Subprogram not used 
! Subprogram not used     real(kind=real_kind) :: ugradv(np,np,2)
! Subprogram not used     real(kind=real_kind) :: dum_cart(np,np,3)
! Subprogram not used 
! Subprogram not used     integer :: component
! Subprogram not used 
! Subprogram not used     ! latlon -> cartesian
! Subprogram not used     do component=1,3
! Subprogram not used        ! Summing along the third dimension is a sum over components for each point.
! Subprogram not used        ! (This is just a faster way of doing a dot product for each grid point,
! Subprogram not used        ! since reindexing the inputs to use the intrinsic effectively would be
! Subprogram not used        ! just asking for trouble.)
! Subprogram not used        dum_cart(:,:,component)=sum( elem%vec_sphere2cart(:,:,component,:)*v(:,:,:) ,3)
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     ! Do ugradv on the cartesian components.
! Subprogram not used     do component=1,3
! Subprogram not used        ! Dot u with the gradient of each component
! Subprogram not used        dum_cart(:,:,component) = sum( u(:,:,:) * &
! Subprogram not used             gradient_sphere(dum_cart(:,:,component),deriv,elem%Dinv) ,3)
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     ! cartesian -> latlon
! Subprogram not used     do component=1,2
! Subprogram not used        ! vec_sphere2cart is its own pseudoinverse.
! Subprogram not used        ugradv(:,:,component)=sum( dum_cart(:,:,:)*elem%vec_sphere2cart(:,:,:,component) ,3)
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used   end function ugradv_sphere



! Subprogram not used   function curl_sphere(s,deriv,elem) result(ds)
! Subprogram not used !
! Subprogram not used !   input s:  scalar  (assumed to be  s khat)
! Subprogram not used !   output  curl(s khat) vector in lat-lon coordinates
! Subprogram not used ! 
! Subprogram not used !   This subroutine can be used to compute divergence free velocity fields,
! Subprogram not used !   since div(ds)=0
! Subprogram not used !
! Subprogram not used !    first compute:  
! Subprogram not used !    curl(s khat) = (1/jacobian) ( ds/dy, -ds/dx ) in contra-variant coordinates
! Subprogram not used !    then map to lat-lon
! Subprogram not used !
! Subprogram not used     type (derivative_t)              :: deriv
! Subprogram not used     type (element_t)                 :: elem
! Subprogram not used     real(kind=real_kind), intent(in) :: s(np,np)
! Subprogram not used 
! Subprogram not used     real(kind=real_kind) :: ds(np,np,2)
! Subprogram not used 
! Subprogram not used     integer i
! Subprogram not used     integer j
! Subprogram not used     integer l
! Subprogram not used 
! Subprogram not used     real(kind=real_kind) ::  dsdx00
! Subprogram not used     real(kind=real_kind) ::  dsdy00
! Subprogram not used     real(kind=real_kind) ::  v1(np,np),v2(np,np)
! Subprogram not used     
! Subprogram not used     do j=1,np
! Subprogram not used        do l=1,np
! Subprogram not used           dsdx00=0.0d0
! Subprogram not used           dsdy00=0.0d0
! Subprogram not used           do i=1,np
! Subprogram not used              dsdx00 = dsdx00 + deriv%Dvv(i,l  )*s(i,j  )
! Subprogram not used              dsdy00 = dsdy00 + deriv%Dvv(i,l  )*s(j  ,i)
! Subprogram not used           end do
! Subprogram not used           v2(l  ,j  ) = -dsdx00*rrearth
! Subprogram not used           v1(j  ,l  ) =  dsdy00*rrearth
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used     ! convert contra -> latlon *and* divide by jacobian
! Subprogram not used     do j=1,np
! Subprogram not used        do i=1,np
! Subprogram not used           ds(i,j,1)=(elem%D(1,1,i,j)*v1(i,j) + elem%D(1,2,i,j)*v2(i,j))/elem%metdet(i,j)
! Subprogram not used           ds(i,j,2)= (elem%D(2,1,i,j)*v1(i,j) + elem%D(2,2,i,j)*v2(i,j))/elem%metdet(i,j)
! Subprogram not used        enddo
! Subprogram not used     enddo
! Subprogram not used  
! Subprogram not used     end function curl_sphere


!--------------------------------------------------------------------------



  function divergence_sphere_wk(v,deriv,elem) result(div)
!
!   input:  v = velocity in lat-lon coordinates
!   ouput:  div(v)  spherical divergence of v, integrated by parts
!
!   Computes  -< grad(psi) dot v > 
!   (the integrated by parts version of < psi div(v) > )
!
!   note: after DSS, divergence_sphere() and divergence_sphere_wk() 
!   are identical to roundoff, as theory predicts.
!
    real(kind=real_kind), intent(in) :: v(np,np,2)  ! in lat-lon coordinates
    type (derivative_t)              :: deriv
    type (element_t)                 :: elem
    real(kind=real_kind) :: div(np,np)

    ! Local

    integer i,j,m,n

    real(kind=real_kind) ::  vtemp(np,np,2)
    real(kind=real_kind) ::  ggtemp(np,np,2)
    real(kind=real_kind) ::  gtemp(np,np,2)
    real(kind=real_kind) ::  psi(np,np)
    real(kind=real_kind) :: xtmp

    ! latlon- > contra
    do j=1,np
       do i=1,np
          vtemp(i,j,1)=(elem%Dinv(1,1,i,j)*v(i,j,1) + elem%Dinv(1,2,i,j)*v(i,j,2))
          vtemp(i,j,2)=(elem%Dinv(2,1,i,j)*v(i,j,1) + elem%Dinv(2,2,i,j)*v(i,j,2))
       enddo
    enddo

    do n=1,np
       do m=1,np

          div(m,n)=0
          do j=1,np
             div(m,n)=div(m,n)-(elem%spheremp(j,n)*vtemp(j,n,1)*deriv%Dvv(m,j) &
                              +elem%spheremp(m,j)*vtemp(m,j,2)*deriv%Dvv(n,j)) &
                              * rrearth
          enddo

       end do
    end do
    
  end function divergence_sphere_wk



! Subprogram not used   function element_boundary_integral(v,deriv,elem) result(result)
! Subprogram not used !
! Subprogram not used !   input:  v = velocity in lat-lon coordinates
! Subprogram not used !   ouput:  result(i,j) = contour integral of PHI_ij * v dot normal
! Subprogram not used !           where PHI_ij = cardinal function at i,j GLL point 
! Subprogram not used !
! Subprogram not used !   this routine is used just to check spectral element integration by parts identities
! Subprogram not used !
! Subprogram not used     real(kind=real_kind), intent(in) :: v(np,np,2)  ! in lat-lon coordinates
! Subprogram not used     type (derivative_t)              :: deriv
! Subprogram not used     type (element_t)                 :: elem
! Subprogram not used     real(kind=real_kind) :: result(np,np)
! Subprogram not used 
! Subprogram not used     ! Local
! Subprogram not used     real(kind=real_kind) :: ucontra(np,np,2)  ! in lat-lon coordinates
! Subprogram not used     integer i,j
! Subprogram not used 
! Subprogram not used     ! latlon->contra
! Subprogram not used     do j=1,np
! Subprogram not used        do i=1,np
! Subprogram not used           ucontra(i,j,1)=(elem%Dinv(1,1,i,j)*v(i,j,1) + elem%Dinv(1,2,i,j)*v(i,j,2))
! Subprogram not used           ucontra(i,j,2)=(elem%Dinv(2,1,i,j)*v(i,j,1) + elem%Dinv(2,2,i,j)*v(i,j,2))
! Subprogram not used        enddo
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     ! note: GLL weights  weight(i) = Mvv_twt(i,i)
! Subprogram not used     result=0
! Subprogram not used     j=1
! Subprogram not used     do i=1,np
! Subprogram not used        result(i,j)=result(i,j)-deriv%Mvv_twt(i,i)*elem%metdet(i,j)*ucontra(i,j,2)*rrearth
! Subprogram not used     enddo
! Subprogram not used     
! Subprogram not used     j=np
! Subprogram not used     do i=1,np
! Subprogram not used        result(i,j)=result(i,j)+deriv%Mvv_twt(i,i)*elem%metdet(i,j)*ucontra(i,j,2)*rrearth
! Subprogram not used     enddo
! Subprogram not used     
! Subprogram not used     i=1
! Subprogram not used     do j=1,np
! Subprogram not used        result(i,j)=result(i,j)-deriv%Mvv_twt(j,j)*elem%metdet(i,j)*ucontra(i,j,1)*rrearth
! Subprogram not used     enddo
! Subprogram not used     
! Subprogram not used     i=np
! Subprogram not used     do j=1,np
! Subprogram not used        result(i,j)=result(i,j)+deriv%Mvv_twt(j,j)*elem%metdet(i,j)*ucontra(i,j,1)*rrearth
! Subprogram not used     enddo
! Subprogram not used   end function element_boundary_integral



! Subprogram not used   function edge_flux_u_cg( v,p,pedges, deriv, elem, u_is_contra) result(result)
! Subprogram not used !
! Subprogram not used !
! Subprogram not used !   input:  v = velocity in contra or lat-lon coordinates (CONTINUIOUS)
! Subprogram not used !           p      = scalar on this element
! Subprogram not used !           pedges = scalar edge data from neighbor elements
! Subprogram not used !
! Subprogram not used !   ouput:  result(i,j) = contour integral of PHI_ij * pstar * v dot normal
! Subprogram not used !           where PHI_ij = cardinal function at i,j GLL point 
! Subprogram not used !           pstar = centered or other flux
! Subprogram not used !
! Subprogram not used     real(kind=real_kind), intent(in) :: v(np,np,2) 
! Subprogram not used     real(kind=real_kind), intent(in) :: p(np,np) 
! Subprogram not used     real(kind=real_kind), intent(in) :: pedges(0:np+1,0:np+1) 
! Subprogram not used     type (derivative_t)              :: deriv
! Subprogram not used     type (element_t)                 :: elem
! Subprogram not used     real(kind=real_kind) :: result(np,np)
! Subprogram not used     logical :: u_is_contra
! Subprogram not used 
! Subprogram not used     ! Local
! Subprogram not used     real(kind=real_kind) :: ucontra(np,np,2)  ! in lat-lon coordinates
! Subprogram not used     real(kind=real_kind) :: flux,pstar
! Subprogram not used     integer i,j
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     result=0
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     if (u_is_contra) then
! Subprogram not used        ucontra=v
! Subprogram not used     else
! Subprogram not used        ! latlon->contra
! Subprogram not used        do j=1,np
! Subprogram not used           do i=1,np
! Subprogram not used              ucontra(i,j,1)=(elem%Dinv(1,1,i,j)*v(i,j,1) + elem%Dinv(1,2,i,j)*v(i,j,2))
! Subprogram not used              ucontra(i,j,2)=(elem%Dinv(2,1,i,j)*v(i,j,1) + elem%Dinv(2,2,i,j)*v(i,j,2))
! Subprogram not used           enddo
! Subprogram not used        enddo
! Subprogram not used     endif
! Subprogram not used     ! upwind
! Subprogram not used     do i=1,np
! Subprogram not used        j=1
! Subprogram not used        pstar=p(i,j)
! Subprogram not used        if (ucontra(i,j,2)>0) pstar=pedges(i,0)
! Subprogram not used        flux = -pstar*ucontra(i,j,2)*( deriv%Mvv_twt(i,i)*elem%metdet(i,j)*rrearth)
! Subprogram not used        result(i,j)=result(i,j)+flux
! Subprogram not used        
! Subprogram not used        j=np
! Subprogram not used        pstar=p(i,j)
! Subprogram not used        if (ucontra(i,j,2)<0) pstar=pedges(i,np+1)
! Subprogram not used        flux = pstar*ucontra(i,j,2)* ( deriv%Mvv_twt(i,i)*elem%metdet(i,j)*rrearth)
! Subprogram not used        result(i,j)=result(i,j)+flux
! Subprogram not used     enddo
! Subprogram not used     
! Subprogram not used     do j=1,np
! Subprogram not used        i=1
! Subprogram not used        pstar=p(i,j)
! Subprogram not used        if (ucontra(i,j,1)>0) pstar=pedges(0,j)
! Subprogram not used        flux = -pstar*ucontra(i,j,1)* ( deriv%Mvv_twt(j,j)*elem%metdet(i,j)*rrearth)
! Subprogram not used        result(i,j)=result(i,j)+flux
! Subprogram not used        
! Subprogram not used        i=np  
! Subprogram not used        pstar=p(i,j)
! Subprogram not used        if (ucontra(i,j,1)<0) pstar=pedges(np+1,j)
! Subprogram not used        flux = pstar*ucontra(i,j,1)* ( deriv%Mvv_twt(j,j)*elem%metdet(i,j)*rrearth)
! Subprogram not used        result(i,j)=result(i,j)+flux
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used   end function edge_flux_u_cg

    

  function vorticity_sphere(v,deriv,elem) result(vort)
!
!   input:  v = velocity in lat-lon coordinates
!   ouput:  spherical vorticity of v
!

    type (derivative_t)              :: deriv
    type (element_t)                 :: elem
    real(kind=real_kind), intent(in) :: v(np,np,2)

    real(kind=real_kind) :: vort(np,np)

    integer i
    integer j
    integer l
    
    real(kind=real_kind) ::  dvdx00
    real(kind=real_kind) ::  dudy00
    real(kind=real_kind) ::  vco(np,np,2)
    real(kind=real_kind) ::  vtemp(np,np)

    ! convert to covariant form
    do j=1,np
       do i=1,np
          vco(i,j,1)=(elem%D(1,1,i,j)*v(i,j,1) + elem%D(2,1,i,j)*v(i,j,2))
          vco(i,j,2)=(elem%D(1,2,i,j)*v(i,j,1) + elem%D(2,2,i,j)*v(i,j,2))
       enddo
    enddo

    do j=1,np
       do l=1,np

          dudy00=0.0d0
	  dvdx00=0.0d0

          do i=1,np
             dvdx00 = dvdx00 + deriv%Dvv(i,l  )*vco(i,j  ,2)
             dudy00 = dudy00 + deriv%Dvv(i,l  )*vco(j  ,i,1)
	  enddo
 
	  vort(l  ,j  ) = dvdx00
	  vtemp(j  ,l  ) = dudy00
	enddo
     enddo

    do j=1,np
       do i=1,np
          vort(i,j)=(vort(i,j)-vtemp(i,j))*(elem%rmetdet(i,j)*rrearth)
       end do
    end do

  end function vorticity_sphere







! Subprogram not used   function vorticity_sphere_diag(v,deriv,elem) result(vort)
! Subprogram not used   !
! Subprogram not used   !   input:  v = velocity in lat-lon coordinates
! Subprogram not used   !   ouput:  diagonal component of spherical vorticity of v
! Subprogram not used   !
! Subprogram not used 
! Subprogram not used       type (derivative_t)              :: deriv
! Subprogram not used       type (element_t)                 :: elem
! Subprogram not used       real(kind=real_kind), intent(in) :: v(np,np,2)
! Subprogram not used 
! Subprogram not used       real(kind=real_kind) :: vort(np,np)
! Subprogram not used 
! Subprogram not used       integer i
! Subprogram not used       integer j
! Subprogram not used       integer l
! Subprogram not used 
! Subprogram not used       real(kind=real_kind) ::  dvdx00
! Subprogram not used       real(kind=real_kind) ::  dudy00
! Subprogram not used       real(kind=real_kind) ::  vco(np,np,2)
! Subprogram not used       real(kind=real_kind) :: vtemp(np,np)
! Subprogram not used       real(kind=real_kind) :: rdx
! Subprogram not used       real(kind=real_kind) :: rdy
! Subprogram not used 
! Subprogram not used       ! convert to covariant form
! Subprogram not used                                                                     
! Subprogram not used       do j=1,np
! Subprogram not used          do i=1,np
! Subprogram not used             vco(i,j,1)=(elem%D(1,1,i,j)*v(i,j,1) + elem%D(2,1,i,j)*v(i,j,2))
! Subprogram not used             vco(i,j,2)=(elem%D(1,2,i,j)*v(i,j,1) + elem%D(2,2,i,j)*v(i,j,2))
! Subprogram not used 
! Subprogram not used 
! Subprogram not used          enddo
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used                                                                                                                
! Subprogram not used       do j=1,np
! Subprogram not used          do l=1,np
! Subprogram not used           
! Subprogram not used             dudy00=0.0d0
! Subprogram not used             dvdx00=0.0d0
! Subprogram not used 
! Subprogram not used             do i=1,np
! Subprogram not used                dvdx00 = dvdx00 + deriv%Dvv_diag(i,l)*vco(i,j ,2)
! Subprogram not used                dudy00 = dudy00 + deriv%Dvv_diag(i,l)*vco(j ,i,1)
! Subprogram not used             enddo 
! Subprogram not used      
! Subprogram not used             vort(l ,j) = dvdx00 
! Subprogram not used             vtemp(j ,l) = dudy00
! Subprogram not used          enddo
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       do j=1,np
! Subprogram not used          do i=1,np 
! Subprogram not used           vort(i,j)=(vort(i,j)-vtemp(i,j))*(elem%rmetdet(i,j)*rrearth)
! Subprogram not used          end do 
! Subprogram not used       end do 
! Subprogram not used      
! Subprogram not used   end function vorticity_sphere_diag



  function divergence_sphere(v,deriv,elem) result(div)
!
!   input:  v = velocity in lat-lon coordinates
!   ouput:  div(v)  spherical divergence of v
!


    real(kind=real_kind), intent(in) :: v(np,np,2)  ! in lat-lon coordinates
    type (derivative_t)              :: deriv
    type (element_t)                 :: elem
    real(kind=real_kind) :: div(np,np)

    ! Local

    integer i
    integer j
    integer l

    real(kind=real_kind) ::  dudx00
    real(kind=real_kind) ::  dvdy00
    real(kind=real_kind) ::  gv(np,np,2),vvtemp(np,np)

    ! convert to contra variant form and multiply by g
    do j=1,np
       do i=1,np
          gv(i,j,1)=elem%metdet(i,j)*(elem%Dinv(1,1,i,j)*v(i,j,1) + elem%Dinv(1,2,i,j)*v(i,j,2))
          gv(i,j,2)=elem%metdet(i,j)*(elem%Dinv(2,1,i,j)*v(i,j,1) + elem%Dinv(2,2,i,j)*v(i,j,2))
       enddo
    enddo

    ! compute d/dx and d/dy         
    do j=1,np
       do l=1,np
          dudx00=0.0d0
          dvdy00=0.0d0
          do i=1,np
             dudx00 = dudx00 + deriv%Dvv(i,l  )*gv(i,j  ,1)
             dvdy00 = dvdy00 + deriv%Dvv(i,l  )*gv(j  ,i,2)
          end do
          div(l  ,j  ) = dudx00
          vvtemp(j  ,l  ) = dvdy00
       end do
    end do

    do j=1,np
       do i=1,np
          div(i,j)=(div(i,j)+vvtemp(i,j))*(elem%rmetdet(i,j)*rrearth)
       end do
    end do
    
  end function divergence_sphere



  function laplace_sphere_wk(s,deriv,elem,var_coef) result(laplace)
!
!   input:  s = scalar
!   ouput:  -< grad(PHI), grad(s) >   = weak divergence of grad(s)
!     note: for this form of the operator, grad(s) does not need to be made C0
!            
    real(kind=real_kind), intent(in) :: s(np,np) 
    logical :: var_coef
    type (derivative_t)              :: deriv
    type (element_t)                 :: elem
    real(kind=real_kind)             :: laplace(np,np)
    real(kind=real_kind)             :: laplace2(np,np)
    integer i,j

    ! Local
    real(kind=real_kind) :: grads(np,np,2), oldgrads(np,np,2)

    grads=gradient_sphere(s,deriv,elem%Dinv)
 
    if (var_coef) then
       if (hypervis_power/=0 ) then
          ! scalar viscosity with variable coefficient
          grads(:,:,1) = grads(:,:,1)*elem%variable_hyperviscosity(:,:)
          grads(:,:,2) = grads(:,:,2)*elem%variable_hyperviscosity(:,:)
       else if (hypervis_scaling /=0 ) then
          ! tensor hv, (3)
          oldgrads=grads
          do j=1,np
             do i=1,np
                grads(i,j,1) = sum(oldgrads(i,j,:)*elem%tensorVisc(1,:,i,j))
                grads(i,j,2) = sum(oldgrads(i,j,:)*elem%tensorVisc(2,:,i,j))
             end do
          end do
       else
          ! do nothing: constant coefficient viscsoity
       endif
    endif

    ! note: divergnece_sphere and divergence_sphere_wk are identical *after* bndry_exchange
    ! if input is C_0.  Here input is not C_0, so we should use divergence_sphere_wk().  
    laplace=divergence_sphere_wk(grads,deriv,elem)

  end function laplace_sphere_wk


  function vlaplace_sphere_wk(v,deriv,elem,var_coef,nu_ratio) result(laplace)
!
!   input:  v = vector in lat-lon coordinates
!   ouput:  weak laplacian of v, in lat-lon coordinates
!
!   logic:
!      tensorHV:     requires cartesian
!      nu_div/=nu:   requires contra formulatino
!
!   One combination NOT supported:  tensorHV and nu_div/=nu then abort
!
    real(kind=real_kind), intent(in) :: v(np,np,2) 
    logical :: var_coef
    type (derivative_t)              :: deriv
    type (element_t)                 :: elem
    real(kind=real_kind), optional :: nu_ratio
    real(kind=real_kind) :: laplace(np,np,2)


    if (hypervis_scaling/=0 .and. var_coef) then
       ! tensorHV is turned on - requires cartesian formulation
       if (present(nu_ratio)) then
          if (nu_ratio /= 1) then
             call abortmp('ERROR: tensorHV can not be used with nu_div/=nu')
          endif
       endif
       laplace=vlaplace_sphere_wk_cartesian(v,deriv,elem,var_coef)
    else  
       ! all other cases, use contra formulation:
       laplace=vlaplace_sphere_wk_contra(v,deriv,elem,var_coef,nu_ratio)
    endif

  end function vlaplace_sphere_wk



! Subprogram not used   function vlaplace_sphere_wk_cartesian(v,deriv,elem,var_coef) result(laplace)
! Subprogram not used !
! Subprogram not used !   input:  v = vector in lat-lon coordinates
! Subprogram not used !   ouput:  weak laplacian of v, in lat-lon coordinates
! Subprogram not used 
! Subprogram not used     real(kind=real_kind), intent(in) :: v(np,np,2) 
! Subprogram not used     logical :: var_coef
! Subprogram not used     type (derivative_t)              :: deriv
! Subprogram not used     type (element_t)                 :: elem
! Subprogram not used     real(kind=real_kind) :: laplace(np,np,2)
! Subprogram not used     ! Local
! Subprogram not used 
! Subprogram not used     integer component
! Subprogram not used     real(kind=real_kind) :: dum_cart(np,np,3)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     ! latlon -> cartesian
! Subprogram not used     do component=1,3
! Subprogram not used        dum_cart(:,:,component)=sum( elem%vec_sphere2cart(:,:,component,:)*v(:,:,:) ,3)
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     ! Do laplace on cartesian comps
! Subprogram not used     do component=1,3
! Subprogram not used        dum_cart(:,:,component) = laplace_sphere_wk(dum_cart(:,:,component),deriv,elem,var_coef)
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     ! cartesian -> latlon
! Subprogram not used     do component=1,2
! Subprogram not used        ! vec_sphere2cart is its own pseudoinverse.
! Subprogram not used        laplace(:,:,component)=sum( dum_cart(:,:,:)*elem%vec_sphere2cart(:,:,:,component) ,3)
! Subprogram not used     end do 
! Subprogram not used 
! Subprogram not used   end function vlaplace_sphere_wk_cartesian



  function vlaplace_sphere_wk_contra(v,deriv,elem,var_coef,nu_ratio) result(laplace)
!
!   input:  v = vector in lat-lon coordinates
!   ouput:  weak laplacian of v, in lat-lon coordinates
!
!   du/dt = laplace(u) = grad(div) - curl(vor)
!   < PHI du/dt > = < PHI laplace(u) >        PHI = covariant, u = contravariant
!                 = < PHI grad(div) >  - < PHI curl(vor) >
!                 = grad_wk(div) - curl_wk(vor)               
!
    real(kind=real_kind), intent(in) :: v(np,np,2) 
    logical :: var_coef
    type (derivative_t)              :: deriv
    type (element_t)                 :: elem
    real(kind=real_kind) :: laplace(np,np,2)
    real(kind=real_kind), optional :: nu_ratio
    ! Local

    integer i,j,l,m,n
    real(kind=real_kind) :: vor(np,np),div(np,np)
    real(kind=real_kind) :: v1,v2,div1,div2,vor1,vor2,phi_x,phi_y

    div=divergence_sphere(v,deriv,elem)
    vor=vorticity_sphere(v,deriv,elem)

    if (var_coef .and. hypervis_power/=0 ) then
          ! scalar viscosity with variable coefficient
          div = div*elem%variable_hyperviscosity(:,:)
          vor = vor*elem%variable_hyperviscosity(:,:)
    endif

    if (present(nu_ratio)) div = nu_ratio*div

    laplace = gradient_sphere_wk_testcov(div,deriv,elem) - &
         curl_sphere_wk_testcov(vor,deriv,elem)

    do n=1,np
       do m=1,np
          ! add in correction so we dont damp rigid rotation
          laplace(m,n,1)=laplace(m,n,1) + 2*elem%spheremp(m,n)*v(m,n,1)*(rrearth**2)
          laplace(m,n,2)=laplace(m,n,2) + 2*elem%spheremp(m,n)*v(m,n,2)*(rrearth**2)
       enddo
    enddo
  end function vlaplace_sphere_wk_contra



!-----------------------------------------------------------------------------------


! Subprogram not used   function gll_to_dgmodal(p,deriv) result(phat)
! Subprogram not used !
! Subprogram not used !   input:  v = velocity in lat-lon coordinates
! Subprogram not used !   ouput:  phat = Legendre coefficients
! Subprogram not used !
! Subprogram not used !   Computes  < g dot p  > = SUM  g(i,j) p(i,j) w(i) w(j)
! Subprogram not used !   (the quadrature approximation on the *reference element* of the integral of p against
! Subprogram not used !    all Legendre polynomials up to degree npdg
! Subprogram not used !
! Subprogram not used !   for npdg < np, this routine gives the (exact) modal expansion of p/spheremp()
! Subprogram not used !
! Subprogram not used     real(kind=real_kind), intent(in) :: p(np,np) 
! Subprogram not used     type (derivative_t)              :: deriv
! Subprogram not used     real(kind=real_kind) :: phat(npdg,npdg)
! Subprogram not used 
! Subprogram not used     ! Local
! Subprogram not used     integer i,j,m,n
! Subprogram not used     real(kind=real_kind) :: A(np,npdg)
! Subprogram not used     A=0
! Subprogram not used     phat=0
! Subprogram not used 
! Subprogram not used     ! N^3 tensor product formulation:
! Subprogram not used     do m=1,npdg
! Subprogram not used     do j=1,np
! Subprogram not used     do i=1,np
! Subprogram not used        A(j,m)=A(j,m)+( p(i,j)*deriv%Mvv_twt(i,i)*deriv%Mvv_twt(j,j)  )*deriv%legdg(m,i)
! Subprogram not used     enddo
! Subprogram not used     enddo
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     do n=1,npdg
! Subprogram not used     do m=1,npdg
! Subprogram not used     do j=1,np
! Subprogram not used        phat(m,n)=phat(m,n)+A(j,m)*deriv%legdg(n,j)
! Subprogram not used     enddo
! Subprogram not used     enddo
! Subprogram not used     enddo
! Subprogram not used     
! Subprogram not used   end function

! Subprogram not used   function dgmodal_to_gll(phat,deriv) result(p)
! Subprogram not used !
! Subprogram not used !   input:  phat = coefficients of Legendre expansion
! Subprogram not used !   ouput:  p    = sum expansion to evaluate phat at GLL points
! Subprogram not used !
! Subprogram not used     real(kind=real_kind) :: p(np,np) 
! Subprogram not used     type (derivative_t)  :: deriv
! Subprogram not used     real(kind=real_kind) :: phat(npdg,npdg)
! Subprogram not used     ! Local
! Subprogram not used     integer i,j,m,n
! Subprogram not used     real(kind=real_kind) :: A(npdg,np)
! Subprogram not used 
! Subprogram not used     p(:,:)=0
! Subprogram not used     ! tensor product version
! Subprogram not used     A=0
! Subprogram not used     do i=1,np
! Subprogram not used     do n=1,npdg
! Subprogram not used     do m=1,npdg
! Subprogram not used        A(n,i)=A(n,i)+phat(m,n)*deriv%legdg(m,i)
! Subprogram not used     enddo
! Subprogram not used     enddo
! Subprogram not used     enddo
! Subprogram not used     do j=1,np
! Subprogram not used     do i=1,np
! Subprogram not used     do n=1,npdg
! Subprogram not used        p(i,j) = p(i,j)+A(n,i)*deriv%legdg(n,j)
! Subprogram not used     enddo
! Subprogram not used     enddo
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used   end function

  ! Given a field defined on the unit element, [-1,1]x[-1,1]
  ! sample values, sampled_val, and integration weights, metdet,
  ! at a number, np, of Gauss-Lobatto-Legendre points. Divide
  ! the square up into intervals by intervals sub-squares so that
  ! there are now intervals**2 sub-cells.  Integrate the 
  ! function defined by sampled_val and metdet over each of these
  ! sub-cells and return the integrated values as an 
  ! intervals by intervals matrix.
  !
  ! Efficiency is obtained by computing and caching the appropriate
  ! integration matrix the first time the function is called.
! Subprogram not used   function subcell_integration(sampled_val, metdet, np, intervals) result(values)
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     integer              , intent(in)  :: np
! Subprogram not used     integer              , intent(in)  :: intervals
! Subprogram not used     real (kind=real_kind), intent(in)  :: sampled_val(np,np)
! Subprogram not used     real (kind=real_kind), intent(in)  :: metdet     (np,np)
! Subprogram not used     real (kind=real_kind)              :: values(intervals,intervals)
! Subprogram not used 
! Subprogram not used     real (kind=real_kind)              :: V          (np,np)
! Subprogram not used     integer i,j
! Subprogram not used 
! Subprogram not used     V  = sampled_val * metdet
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     if (.not.ALLOCATED(integration_matrix)      .or. &
! Subprogram not used         SIZE(integration_matrix,1).ne.intervals .or. &
! Subprogram not used         SIZE(integration_matrix,2).ne.np) then
! Subprogram not used       call allocate_subcell_integration_matrix(np,intervals)
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! Multiply the sampled values by the weighted jacobians.  
! Subprogram not used     ! Symmetry allows us to write this as J^t V J
! Subprogram not used     ! where J is a vector.  
! Subprogram not used 
! Subprogram not used     values = MATMUL(integration_matrix, &
! Subprogram not used              MATMUL(V,TRANSPOSE(integration_matrix)))
! Subprogram not used 
! Subprogram not used   end function subcell_integration


  ! Helper subroutine that will fill in a matrix needed to 
  ! integrate a function defined on the GLL points of a unit
  ! square on sub-cells.  So np is the number of integration
  ! GLL points defined on the unit square (actually [-1,1]x[-1,1])
  ! and intervals is the number to cut it up into, say a 3 by 3
  ! set of uniform sub-cells.  This function will fill the 
  ! subcell_integration matrix with the correct coefficients
  ! to integrate over each subcell.  
! Subprogram not used   subroutine allocate_subcell_integration_matrix(np, intervals)
! Subprogram not used     !-----------------
! Subprogram not used     !-----------------
! Subprogram not used     use quadrature_mod, only : gausslobatto, quadrature_t
! Subprogram not used     
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     integer              , intent(in)  :: np
! Subprogram not used     integer              , intent(in)  :: intervals
! Subprogram not used     real (kind=real_kind)              :: values(intervals,intervals)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     real(kind=real_kind), parameter :: zero = 0.0D0, one=1.0D0, two=2.0D0
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     real (kind=real_kind) :: sub_gll        (intervals,np)
! Subprogram not used 
! Subprogram not used     real (kind=real_kind) :: Lagrange_interp(intervals,np,np)
! Subprogram not used     type (quadrature_t)   :: gll 
! Subprogram not used 
! Subprogram not used     real (kind=real_kind) :: legrange_div(np)
! Subprogram not used     real (kind=real_kind) :: a,b,x,y, x_j, x_i 
! Subprogram not used     real (kind=real_kind) :: r(1) 
! Subprogram not used     integer i,j,n,m
! Subprogram not used 
! Subprogram not used     if (ALLOCATED(integration_matrix)) deallocate(integration_matrix)
! Subprogram not used     allocate(integration_matrix(intervals,np))
! Subprogram not used 
! Subprogram not used     gll = gausslobatto(np)
! Subprogram not used  
! Subprogram not used     ! The GLL (Gauss-Lobatto-Legendre) points are from [-1,1], 
! Subprogram not used     ! we have a bunch of sub-intervals defined by intervals that 
! Subprogram not used     ! go from [a,b] so we need to linearly map [-1,1] -> [a,b] 
! Subprogram not used     ! all the  GLL points by  y = (a/2)(1-x) + (b/2)(1+x)
! Subprogram not used     do i=1,intervals
! Subprogram not used       a = -one + (i-one)*two/intervals   
! Subprogram not used       b = -one +  i     *two/intervals  
! Subprogram not used       sub_gll(i,:) = (a+b)/two + gll%points(:)/intervals
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     ! Now to interpolate from the values at the input GLL
! Subprogram not used     ! points to the sub-GLL points.  Do this by Lagrange
! Subprogram not used     ! interpolation.  The jth Lagrange interpolating polynomial
! Subprogram not used     ! for points x_i is 
! Subprogram not used     !              \prod_{i\ne j} (x-x_i)/(x_j-x_i)
! Subprogram not used     ! These are then multiplied by the sampled values y_i 
! Subprogram not used     ! and summed. 
! Subprogram not used     
! Subprogram not used     ! Save some time by pre-computing the denominitor. I think 
! Subprogram not used     ! this is OK since all the points are of order 1 so should
! Subprogram not used     ! be well behaved.
! Subprogram not used     do n = 1,np
! Subprogram not used       x_j = gll%points(n)
! Subprogram not used       x   = one 
! Subprogram not used       do m = 1,np 
! Subprogram not used         if (m.ne.n) then
! Subprogram not used           x_i = gll%points(m)
! Subprogram not used           x = x * (x_j-x_i)
! Subprogram not used         endif
! Subprogram not used       end do
! Subprogram not used       legrange_div(n)= x
! Subprogram not used     end do 
! Subprogram not used     do i=1,intervals
! Subprogram not used       do n=1,np
! Subprogram not used         x = sub_gll(i,n)
! Subprogram not used         do j = 1,np
! Subprogram not used           y = one
! Subprogram not used           do m = 1,np
! Subprogram not used             if (m.ne.j) then
! Subprogram not used               x_i = gll%points(m)
! Subprogram not used               y = y * (x-x_i)
! Subprogram not used             end if
! Subprogram not used           end do
! Subprogram not used           Lagrange_interp(i,n,j) = y/legrange_div(j)
! Subprogram not used         end do
! Subprogram not used       end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     ! Integration is the GLL weights times Jacobians times
! Subprogram not used     ! the interpolated values:
! Subprogram not used     !                   w^t I Y I^t w 
! Subprogram not used     ! where  
! Subprogram not used     ! w is GLL weights and Jacobians, 
! Subprogram not used     ! I is the Lagrange_interp matrix, and
! Subprogram not used     ! Y is the coefficient matrix, sampled_val.
! Subprogram not used     ! This can be written  J Y J^t where
! Subprogram not used     !                       J = w^t I
! Subprogram not used     ! J is integration_matrix
! Subprogram not used     do i=1,intervals
! Subprogram not used       integration_matrix(i,:) = MATMUL(gll%weights(:),Lagrange_interp(i,:,:))
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     ! There is still the Jacobian to consider.  We are 
! Subprogram not used     ! integrating over [a,b] x [c,d] where 
! Subprogram not used     !        |b-a| = |d-c| = 2/Intervals
! Subprogram not used     ! Multiply the weights appropriately given that 
! Subprogram not used     ! they are defined for a 2x2 square
! Subprogram not used     integration_matrix = integration_matrix/intervals
! Subprogram not used 
! Subprogram not used   end subroutine allocate_subcell_integration_matrix



end module derivative_mod








