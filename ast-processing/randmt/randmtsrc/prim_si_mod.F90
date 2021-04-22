module prim_si_mod
  implicit none
  private
  public :: preq_omegap
  public :: preq_omega_ps
  public :: preq_omega_lnps
  public :: preq_hydrostatic, geopotential_t
  public :: preq_pressure
  public :: preq_vertadv
contains
	
! ==========================================================
! Implicit system for semi-implicit primitive equations.
! ==========================================================


! Subprogram not used   subroutine preq_vertadv(T, v, eta_dot_dp_deta, rpdel, &
! Subprogram not used        T_vadv, v_vadv)
! Subprogram not used     use kinds,              only : real_kind
! Subprogram not used     use dimensions_mod,     only : nlev, np, nlevp
! Subprogram not used     implicit none
! Subprogram not used     
! Subprogram not used     real (kind=real_kind), intent(in) :: T(np,np,nlev)
! Subprogram not used     real (kind=real_kind), intent(in) :: v(np,np,2,nlev)
! Subprogram not used     real (kind=real_kind), intent(in) :: eta_dot_dp_deta(np,np,nlevp)
! Subprogram not used     real (kind=real_kind), intent(in) :: rpdel(np,np,nlev)
! Subprogram not used 
! Subprogram not used     real (kind=real_kind), intent(out) :: T_vadv(np,np,nlev)
! Subprogram not used     real (kind=real_kind), intent(out) :: v_vadv(np,np,2,nlev)
! Subprogram not used 
! Subprogram not used     ! ========================
! Subprogram not used     ! Local Variables
! Subprogram not used     ! ========================
! Subprogram not used 
! Subprogram not used     integer :: i,j,k
! Subprogram not used     real (kind=real_kind) :: facp, facm
! Subprogram not used 
! Subprogram not used     do j=1,np   !   Loop inversion (AAM)
! Subprogram not used 
! Subprogram not used     ! ===========================================================
! Subprogram not used     ! Compute vertical advection of T and v from eq. (3.b.1)
! Subprogram not used     !
! Subprogram not used     ! k = 1 case:
! Subprogram not used     ! ===========================================================
! Subprogram not used 
! Subprogram not used        k=1
! Subprogram not used        do i=1,np 
! Subprogram not used           facp            = (0.5_real_kind*rpdel(i,j,k))*eta_dot_dp_deta(i,j,k+1)
! Subprogram not used           T_vadv(i,j,k)   = facp*(T(i,j,k+1)- T(i,j,k))
! Subprogram not used           v_vadv(i,j,1,k) = facp*(v(i,j,1,k+1)- v(i,j,1,k))
! Subprogram not used           v_vadv(i,j,2,k) = facp*(v(i,j,2,k+1)- v(i,j,2,k))
! Subprogram not used        end do
! Subprogram not used 
! Subprogram not used     ! ===========================================================
! Subprogram not used     ! vertical advection
! Subprogram not used     !
! Subprogram not used     ! 1 < k < nlev case:
! Subprogram not used     ! ===========================================================
! Subprogram not used 
! Subprogram not used        do k=2,nlev-1
! Subprogram not used           do i=1,np
! Subprogram not used              facp            = (0.5_real_kind*rpdel(i,j,k))*eta_dot_dp_deta(i,j,k+1)
! Subprogram not used              facm            = (0.5_real_kind*rpdel(i,j,k))*eta_dot_dp_deta(i,j,k)
! Subprogram not used              T_vadv(i,j,k)   = facp*(T(i,j,k+1)- T(i,j,k)) + &
! Subprogram not used                   facm*(T(i,j,k)- T(i,j,k-1))
! Subprogram not used              v_vadv(i,j,1,k) = facp*(v(i,j,1,k+1)- v(i,j,1,k)) + &
! Subprogram not used                   facm*(v(i,j,1,k)- v(i,j,1,k-1))
! Subprogram not used              v_vadv(i,j,2,k) = facp*(v(i,j,2,k+1)- v(i,j,2,k)) + &
! Subprogram not used                   facm*(v(i,j,2,k)- v(i,j,2,k-1))
! Subprogram not used           end do
! Subprogram not used        end do
! Subprogram not used 
! Subprogram not used     ! ===========================================================
! Subprogram not used     ! vertical advection
! Subprogram not used     !
! Subprogram not used     ! k = nlev case:
! Subprogram not used     ! ===========================================================
! Subprogram not used 
! Subprogram not used        k=nlev
! Subprogram not used        do i=1,np
! Subprogram not used           facm            = (0.5_real_kind*rpdel(i,j,k))*eta_dot_dp_deta(i,j,k)
! Subprogram not used           T_vadv(i,j,k)   = facm*(T(i,j,k)- T(i,j,k-1))
! Subprogram not used           v_vadv(i,j,1,k) = facm*(v(i,j,1,k)- v(i,j,1,k-1))
! Subprogram not used           v_vadv(i,j,2,k) = facm*(v(i,j,2,k)- v(i,j,2,k-1))
! Subprogram not used        end do
! Subprogram not used 
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used   end subroutine preq_vertadv




!----------------------------------------------------------------------- 
! preq_omegap:

! Purpose: 
! Calculate (omega/p) needed for the Thermodynamics Equation
! 
! Method: 
! Simplified version in CAM2 for clarity
! 
! Author: Modified by Rich Loft for use in HOMME. 
! 
!-----------------------------------------------------------------------

! Subprogram not used   subroutine preq_omegap(div     ,vgrad_ps,pdel    ,rpmid, &
! Subprogram not used        hybm    ,hybd    ,omegap   )
! Subprogram not used     use kinds, only : real_kind
! Subprogram not used     use dimensions_mod, only : np, nlev
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     !------------------------------Arguments---------------------------------------------------------------
! Subprogram not used     real(kind=real_kind), intent(in) :: div(np,np,nlev)      ! divergence
! Subprogram not used     real(kind=real_kind), intent(in) :: vgrad_ps(np,np,nlev) ! v.grad(ps)
! Subprogram not used     real(kind=real_kind), intent(in) :: pdel(np,np,nlev)     ! layer thicknesses (pressure)
! Subprogram not used     real(kind=real_kind), intent(in) :: rpmid(np,np,nlev)    ! 1./pmid
! Subprogram not used     real(kind=real_kind), intent(in) :: hybm(nlev)           ! Hybrid B coefficient on mid levels
! Subprogram not used     real(kind=real_kind), intent(in) :: hybd(nlev)           ! Hybrid dB coefficient on mid levels
! Subprogram not used     real(kind=real_kind), intent(out):: omegap(np,np,nlev)   ! vertical pressure velocity
! Subprogram not used     !------------------------------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     !---------------------------Local workspace-----------------------------
! Subprogram not used     integer i,j,k                         ! longitude, level indices
! Subprogram not used     real(kind=real_kind) term             ! one half of basic term in omega/p summation 
! Subprogram not used     real(kind=real_kind) Ckk              ! diagonal term of energy conversion matrix
! Subprogram not used     real(kind=real_kind) suml(np,np)      ! partial sum over l = (1, k-1)
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     ! =========================
! Subprogram not used     ! Zero partial sum
! Subprogram not used     ! =========================
! Subprogram not used 
! Subprogram not used     do j=1,np
! Subprogram not used        do i=1,np
! Subprogram not used           suml(i,j)=0.0_real_kind
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     ! =============================
! Subprogram not used     ! Compute omegap 
! Subprogram not used     ! =============================
! Subprogram not used 
! Subprogram not used     do k=1,nlev
! Subprogram not used        do j=1,np
! Subprogram not used           do i=1,np
! Subprogram not used              Ckk       = 0.5_real_kind
! Subprogram not used              term      = Ckk*(div(i,j,k)*pdel(i,j,k) + vgrad_ps(i,j,k)*hybd(k))
! Subprogram not used              suml(i,j) = suml(i,j) + term
! Subprogram not used              omegap(i,j,k) = rpmid(i,j,k)*(hybm(k)*vgrad_ps(i,j,k) - suml(i,j))
! Subprogram not used              suml(i,j) = suml(i,j) + term
! Subprogram not used           end do
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used   end subroutine preq_omegap



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!
!  compute omega/p using ps, modeled after CCM3 formulas 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  subroutine preq_omega_ps(omega_p,hvcoord,p,vgrad_p,divdp)
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev
    use hybvcoord_mod, only : hvcoord_t
    implicit none


    !------------------------------Arguments---------------------------------------------------------------
    real(kind=real_kind), intent(in) :: divdp(np,np,nlev)      ! divergence
    real(kind=real_kind), intent(in) :: vgrad_p(np,np,nlev) ! v.grad(p)
    real(kind=real_kind), intent(in) :: p(np,np,nlev)     ! layer thicknesses (pressure)
    type (hvcoord_t),     intent(in) :: hvcoord
    real(kind=real_kind), intent(out):: omega_p(np,np,nlev)   ! vertical pressure velocity
    !------------------------------------------------------------------------------------------------------

    !---------------------------Local workspace-----------------------------
    integer i,j,k                         ! longitude, level indices
    real(kind=real_kind) term             ! one half of basic term in omega/p summation 
    real(kind=real_kind) Ckk,Ckl          ! diagonal term of energy conversion matrix
    real(kind=real_kind) suml(np,np)      ! partial sum over l = (1, k-1)
    !-----------------------------------------------------------------------

       do j=1,np   !   Loop inversion (AAM)

          do i=1,np
             ckk = 0.5d0/p(i,j,1)
             term = divdp(i,j,1)
!             omega_p(i,j,1) = hvcoord%hybm(1)*vgrad_ps(i,j,1)/p(i,j,1)
             omega_p(i,j,1) = vgrad_p(i,j,1)/p(i,j,1)
             omega_p(i,j,1) = omega_p(i,j,1) - ckk*term
             suml(i,j) = term
          end do

          do k=2,nlev-1
             do i=1,np
                ckk = 0.5d0/p(i,j,k)
                ckl = 2*ckk
                term = divdp(i,j,k)
!                omega_p(i,j,k) = hvcoord%hybm(k)*vgrad_ps(i,j,k)/p(i,j,k)
                omega_p(i,j,k) = vgrad_p(i,j,k)/p(i,j,k)
                omega_p(i,j,k) = omega_p(i,j,k) - ckl*suml(i,j) - ckk*term
                suml(i,j) = suml(i,j) + term

             end do
          end do

          do i=1,np
             ckk = 0.5d0/p(i,j,nlev)
             ckl = 2*ckk
             term = divdp(i,j,nlev)
!             omega_p(i,j,nlev) = hvcoord%hybm(nlev)*vgrad_ps(i,j,nlev)/p(i,j,nlev)
             omega_p(i,j,nlev) = vgrad_p(i,j,nlev)/p(i,j,nlev)
             omega_p(i,j,nlev) = omega_p(i,j,nlev) - ckl*suml(i,j) - ckk*term
          end do

       end do

  end subroutine preq_omega_ps




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!
!  compute omega/p using lnps 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! Subprogram not used   subroutine preq_omega_lnps(omega_p,hvcoord,ps,p,dp,vgrad_lnps,div)
! Subprogram not used     use kinds, only : real_kind
! Subprogram not used     use dimensions_mod, only : np, nlev
! Subprogram not used     use hybvcoord_mod, only : hvcoord_t
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     !------------------------------Arguments---------------------------------------------------------------
! Subprogram not used     real(kind=real_kind), intent(in) :: div(np,np,nlev)      ! divergence
! Subprogram not used     real(kind=real_kind), intent(in) :: vgrad_lnps(np,np,nlev) ! v.grad(ps)
! Subprogram not used     real(kind=real_kind), intent(in) :: p(np,np,nlev)      ! pressure
! Subprogram not used     real(kind=real_kind), intent(in) :: dp(np,np,nlev)     ! dp/dn
! Subprogram not used     real(kind=real_kind), intent(in) :: ps(np,np)
! Subprogram not used     type (hvcoord_t),     intent(in) :: hvcoord
! Subprogram not used     real(kind=real_kind), intent(out):: omega_p(np,np,nlev)   ! vertical pressure velocity
! Subprogram not used     !------------------------------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     !---------------------------Local workspace-----------------------------
! Subprogram not used     integer i,j,k                         ! longitude, level indices
! Subprogram not used     real(kind=real_kind) term             ! one half of basic term in omega/p summation 
! Subprogram not used     real(kind=real_kind) Ckk,Ckl          ! diagonal term of energy conversion matrix
! Subprogram not used     real(kind=real_kind) suml(np,np)      ! partial sum over l = (1, k-1)
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used        do j=1,np	
! Subprogram not used           do i=1,np
! Subprogram not used              ckk = 0.5d0/p(i,j,1)
! Subprogram not used              term = div(i,j,1)*dp(i,j,1) + vgrad_lnps(i,j,1)*ps(i,j)*hvcoord%hybd(1)
! Subprogram not used              omega_p(i,j,1) = hvcoord%hybm(1)*(ps(i,j)/p(i,j,1))*vgrad_lnps(i,j,1)
! Subprogram not used              omega_p(i,j,1) = omega_p(i,j,1) - ckk*term
! Subprogram not used              suml(i,j) = term
! Subprogram not used           end do
! Subprogram not used        end do
! Subprogram not used 
! Subprogram not used        do k=2,nlev-1
! Subprogram not used 
! Subprogram not used           do j=1,np
! Subprogram not used              do i=1,np
! Subprogram not used                 ckk = 0.5d0/p(i,j,k)
! Subprogram not used                 ckl = 2*ckk
! Subprogram not used                 term = div(i,j,k)*dp(i,j,k) + vgrad_lnps(i,j,k)*ps(i,j)*hvcoord%hybd(k)
! Subprogram not used                 omega_p(i,j,k) = hvcoord%hybm(k)*(ps(i,j)/p(i,j,k))*vgrad_lnps(i,j,k)
! Subprogram not used                 omega_p(i,j,k) = omega_p(i,j,k) - ckl*suml(i,j) - ckk*term
! Subprogram not used                 suml(i,j) = suml(i,j) + term
! Subprogram not used 
! Subprogram not used              end do
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used        end do
! Subprogram not used 
! Subprogram not used        do j=1,np
! Subprogram not used           do i=1,np
! Subprogram not used              ckk = 0.5d0/p(i,j,nlev)
! Subprogram not used              ckl = 2*ckk
! Subprogram not used              term = div(i,j,nlev)*dp(i,j,nlev) + vgrad_lnps(i,j,nlev)*ps(i,j)*hvcoord%hybd(nlev)
! Subprogram not used              omega_p(i,j,nlev) = hvcoord%hybm(nlev)*(ps(i,j)/p(i,j,nlev))*vgrad_lnps(i,j,nlev)
! Subprogram not used              omega_p(i,j,nlev) = omega_p(i,j,nlev) - ckl*suml(i,j) - ckk*term
! Subprogram not used           end do
! Subprogram not used        end do
! Subprogram not used 
! Subprogram not used   end subroutine preq_omega_lnps



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!
!  CCM3 hydrostatic integral
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  subroutine preq_hydrostatic(phi,phis,T_v,p,dp)
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev
    use physical_constants, only : rgas
!    use hybvcoord_mod, only : hvcoord_t
    implicit none


    !------------------------------Arguments---------------------------------------------------------------
    real(kind=real_kind), intent(out) :: phi(np,np,nlev)     
    real(kind=real_kind), intent(in) :: phis(np,np)
    real(kind=real_kind), intent(in) :: T_v(np,np,nlev)
    real(kind=real_kind), intent(in) :: p(np,np,nlev)   
    real(kind=real_kind), intent(in) :: dp(np,np,nlev)  
 !   type (hvcoord_t),     intent(in) :: hvcoord
    !------------------------------------------------------------------------------------------------------

    !---------------------------Local workspace-----------------------------
    integer i,j,k                         ! longitude, level indices
    real(kind=real_kind) Hkk,Hkl          ! diagonal term of energy conversion matrix
    real(kind=real_kind), dimension(np,np,nlev) :: phii       ! Geopotential at interfaces
    !-----------------------------------------------------------------------

       do j=1,np   !   Loop inversion (AAM)

          do i=1,np
             hkk = dp(i,j,nlev)*0.5d0/p(i,j,nlev)
             hkl = 2*hkk
             phii(i,j,nlev)  = Rgas*T_v(i,j,nlev)*hkl
             phi(i,j,nlev) = phis(i,j) + Rgas*T_v(i,j,nlev)*hkk 
          end do

          do k=nlev-1,2,-1
             do i=1,np
                ! hkk = dp*ckk
                hkk = dp(i,j,k)*0.5d0/p(i,j,k)
                hkl = 2*hkk
                phii(i,j,k) = phii(i,j,k+1) + Rgas*T_v(i,j,k)*hkl
                phi(i,j,k) = phis(i,j) + phii(i,j,k+1) + Rgas*T_v(i,j,k)*hkk
             end do
          end do

          do i=1,np
             ! hkk = dp*ckk
             hkk = 0.5d0*dp(i,j,1)/p(i,j,1)
             phi(i,j,1) = phis(i,j) + phii(i,j,2) + Rgas*T_v(i,j,1)*hkk
          end do

       end do


end subroutine preq_hydrostatic



!
!  The hydrostatic routine from 1 physics.
!  (FV stuff removed)
!  t,q input changed to take t_v
!  removed gravit, so this routine returns PHI, not zm
! Subprogram not used subroutine geopotential_t(                                 &
! Subprogram not used        pmid   , pdel   ,  tv      , rair   ,  zm)
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: 
! Subprogram not used ! Compute the geopotential height (above the surface) at the midpoints and 
! Subprogram not used ! interfaces using the input temperatures and pressures.
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used     use dimensions_mod,     only : nlev, nlevp, np
! Subprogram not used     use kinds, only : real_kind
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used !
! Subprogram not used ! Input arguments
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     real(real_kind), intent(in) :: pmid (np*np,nlev)    ! Midpoint pressures
! Subprogram not used     real(real_kind), intent(in) :: pdel (np*np,nlev)    ! layer thickness
! Subprogram not used     real(real_kind), intent(in) :: tv    (np*np,nlev)    ! temperature
! Subprogram not used     real(real_kind), intent(in) :: rair                 ! Gas constant for dry air
! Subprogram not used     ! real(real_kind), intent(in) :: gravit               ! Acceleration of gravity
! Subprogram not used     ! real(real_kind), intent(in) :: zvir                 ! rh2o/rair - 1
! Subprogram not used 
! Subprogram not used ! Output arguments
! Subprogram not used 
! Subprogram not used     real(real_kind), intent(out) :: zm(np*np,nlev)      ! Geopotential height at mid level
! Subprogram not used !
! Subprogram not used !---------------------------Local variables-----------------------------
! Subprogram not used     integer :: ncol=np*np             ! Number of longitudes
! Subprogram not used 
! Subprogram not used     integer  :: i,k                ! Lon, level indices
! Subprogram not used     real(real_kind) :: hkk(np*np)         ! diagonal element of hydrostatic matrix
! Subprogram not used     real(real_kind) :: hkl(np*np)         ! off-diagonal element
! Subprogram not used     real(real_kind) :: rog                ! Rair / gravit
! Subprogram not used     real(real_kind) :: zi(np*np,nlevp)     ! Height above surface at interfaces
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !    rog = rair/gravit
! Subprogram not used     rog = rair
! Subprogram not used 
! Subprogram not used ! The surface height is zero by definition.
! Subprogram not used     do i = 1,ncol
! Subprogram not used        zi(i,nlevp) = 0.0_real_kind
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used ! Compute zi, zm from bottom up. 
! Subprogram not used ! Note, zi(i,k) is the interface above zm(i,k)
! Subprogram not used     do k = nlev, 1, -1
! Subprogram not used ! First set hydrostatic elements consistent with dynamics
! Subprogram not used        do i = 1,ncol
! Subprogram not used           hkl(i) = pdel(i,k) / pmid(i,k)
! Subprogram not used           hkk(i) = 0.5_real_kind * hkl(i)
! Subprogram not used        end do
! Subprogram not used 
! Subprogram not used ! Now compute tv, zm, zi
! Subprogram not used        do i = 1,ncol
! Subprogram not used           ! tvfac   = 1._r8 + zvir * q(i,k)
! Subprogram not used           ! tv      = t(i,k) * tvfac
! Subprogram not used           zm(i,k) = zi(i,k+1) + rog * tv(i,k) * hkk(i)
! Subprogram not used           zi(i,k) = zi(i,k+1) + rog * tv(i,k) * hkl(i)
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     return
! Subprogram not used   end subroutine geopotential_t





!----------------------------------------------------------------------- 
! preq_pressure:
!
! Purpose: 
! Define the pressures of the interfaces and midpoints from the
! coordinate definitions and the surface pressure. Originally plevs0!
! 
! Method: 
! 
! Author: B. Boville/ Adapted for HOMME by Rich Loft
! 
!-----------------------------------------------------------------------
!
! $Id: prim_si_mod.F90,v 2.10 2005/10/14 20:17:22 jedwards Exp $
! $Author: jedwards $
!
!-----------------------------------------------------------------------

! Subprogram not used   subroutine preq_pressure (ps0,  ps,               &
! Subprogram not used        hyai, hybi, hyam, hybm, &
! Subprogram not used        pint, pmid, pdel)
! Subprogram not used     use kinds, only : real_kind
! Subprogram not used     use dimensions_mod, only : np, nlev, nlevp
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     real(kind=real_kind), intent(in)  :: ps0                ! Hybrid coordinate reference pressure (pascals)
! Subprogram not used     real(kind=real_kind), intent(in)  :: ps(np,np)          ! Surface pressure (pascals)
! Subprogram not used     real(kind=real_kind), intent(in)  :: hyai(nlevp)        ! Hybrid interface A coefficients
! Subprogram not used     real(kind=real_kind), intent(in)  :: hybi(nlevp)        ! Hybrid interface B coefficients
! Subprogram not used     real(kind=real_kind), intent(in)  :: hyam(nlev)         ! Hybrid midpoint  A coefficients
! Subprogram not used     real(kind=real_kind), intent(in)  :: hybm(nlev)         ! Hybrid midpoint  B coefficients
! Subprogram not used     real(kind=real_kind), intent(out) :: pint(np,np,nlevp)  ! Pressure at model interfaces
! Subprogram not used     real(kind=real_kind), intent(out) :: pmid(np,np,nlev)   ! Pressure at model levels
! Subprogram not used     real(kind=real_kind), intent(out) :: pdel(np,np,nlev)   ! Layer thickness (pint(k+1) - pint(k))
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     !---------------------------Local workspace-----------------------------
! Subprogram not used     integer i,j,k             ! Horizontal, level indices
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !
! Subprogram not used     ! Set interface pressures
! Subprogram not used     !
! Subprogram not used     do k=1,nlevp
! Subprogram not used        do j=1,np
! Subprogram not used           do i=1,np
! Subprogram not used              pint(i,j,k) = hyai(k)*ps0 + hybi(k)*ps(i,j)
! Subprogram not used           end do
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used     !
! Subprogram not used     ! Set midpoint pressures and layer thicknesses
! Subprogram not used     !
! Subprogram not used     do k=1,nlev
! Subprogram not used        do j=1,np
! Subprogram not used           do i=1,np
! Subprogram not used              pmid(i,j,k) = hyam(k)*ps0 + hybm(k)*ps(i,j)
! Subprogram not used              pdel(i,j,k) = pint(i,j,k+1) - pint(i,j,k)
! Subprogram not used           end do
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used   end subroutine preq_pressure




end module prim_si_mod
