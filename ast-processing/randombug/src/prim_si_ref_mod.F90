module prim_si_ref_mod
  use kinds, only: r8 => real_kind, iulog
  use dimensions_mod, only: plev => nlev, plevp => nlevp
  implicit none
  private

  type, public :: ref_state_t
     real(r8) psr                ! reference surface pressure for linearization
     real(r8) hypi(plevp)        ! reference pressures at interfaces
     real(r8) hypm(plev)         ! reference pressures at midpoints
     real(r8) hypd(plev)         ! reference pressure layer thickness

     real(r8) Tref(plev)         ! reference temperature
     real(r8) RTref(plev)        ! coefficient for ln(ps) term in velocity equation
     real(r8) Pvec(plev)         ! diagonal P matrix
     real(r8) Href(plev,plev)    ! reference hydrostatic matrix (Href)
     real(r8) Cref(plev,plev)    ! reference hydrostatic matrix (Cmat)
     real(r8) Tmat(plev,plev)    ! tau matrix (Tmat)

     real(r8) Amat(plev,plev)     ! vertical structure matrix
     real(r8) Amat_inv(plev,plev) ! inverse vertical structure matrix
     real(r8) Emat(plev,plev)     ! right eigenvector matrix
     real(r8) Emat_inv(plev,plev) ! right eigenvector matrix
     real(r8) Lambda(plev)        ! solver eigenvalues
  end type ref_state_t

  public  :: prim_si_refstate_init
  public  :: set_vert_struct_mat
  public  :: prim_set_mass
contains

  ! =====================================================
  ! prim_si_refstate_init:
  !
  ! given a hybrid vertical coordinate system, initialize
  ! the reference pressures needed by the semi-implicit.
  ! =====================================================

! Subprogram not used   function prim_si_refstate_init(lprint,masterproc,hvcoord) result(refstate)
! Subprogram not used     use hybvcoord_mod, only : hvcoord_t
! Subprogram not used 
! Subprogram not used     logical, intent(in) :: lprint
! Subprogram not used     logical, intent(in) :: masterproc
! Subprogram not used     type (hvcoord_t)    :: hvcoord
! Subprogram not used     type (ref_state_t)  :: refstate
! Subprogram not used 
! Subprogram not used     ! =============================
! Subprogram not used     ! Local variables
! Subprogram not used     ! =============================
! Subprogram not used 
! Subprogram not used     integer k                 ! Level index
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     refstate%psr    = 1000.0_r8  ! Reference surface pressure (millibars)
! Subprogram not used 
! Subprogram not used     ! ================================================================
! Subprogram not used     ! Set mean temperature
! Subprogram not used     ! NOTE: Making t0 an actual function of height ***DOES NOT WORK***
! Subprogram not used     ! ================================================================
! Subprogram not used 
! Subprogram not used     do k=1,plev
! Subprogram not used        refstate%Tref(k) = 300.0_r8
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     ! ===========================================
! Subprogram not used     ! Set reference state midpoint pressures
! Subprogram not used     ! ===========================================
! Subprogram not used 
! Subprogram not used     do k=1,plev
! Subprogram not used        refstate%hypm(k) = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*refstate%psr
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     ! ============================================
! Subprogram not used     ! Reference state interface pressures
! Subprogram not used     ! ============================================
! Subprogram not used 
! Subprogram not used     do k=1,plevp
! Subprogram not used        refstate%hypi(k) = hvcoord%hyai(k)*hvcoord%ps0 + hvcoord%hybi(k)*refstate%psr
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     ! ============================================
! Subprogram not used     ! Reference state layer thicknesses
! Subprogram not used     ! ============================================
! Subprogram not used 
! Subprogram not used     do k=1,plev
! Subprogram not used        refstate%hypd(k) = refstate%hypi(k+1) - refstate%hypi(k)
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     if (lprint.and.masterproc) then
! Subprogram not used        write(iulog,9820)
! Subprogram not used        do k=1,plev
! Subprogram not used           write(iulog,9830) k, refstate%hypi(k)
! Subprogram not used           write(iulog,9840) refstate%hypm(k), refstate%hypd(k)
! Subprogram not used        end do
! Subprogram not used        write(iulog,9830) plevp, refstate%hypi(plevp)
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     call prim_si_refmat_init(lprint,masterproc,refstate)
! Subprogram not used 
! Subprogram not used 9820 format(1x,'reference pressures (Pa)')
! Subprogram not used 9830 format(1x,i3,f15.4)
! Subprogram not used 9840 format(1x,3x,15x,2f15.4)
! Subprogram not used 
! Subprogram not used   end function prim_si_refstate_init

  ! ==============================================
  ! prim_si_refmat_init:
  ! 
  ! initialize the reference state hydrostatic and
  ! energy conversion matrices.
  ! ==============================================

! Subprogram not used   subroutine prim_si_refmat_init(lprint,masterproc,refstate)
! Subprogram not used     use physical_constants , only :  rgas
! Subprogram not used 
! Subprogram not used     logical, intent(in)  :: lprint            ! print/no print of ref matrices
! Subprogram not used     logical, intent(in)  :: masterproc        ! is the master process?
! Subprogram not used     type(ref_state_t)    :: refstate          ! reference state structure
! Subprogram not used 
! Subprogram not used     integer k,l,kk,kkk        ! level indices
! Subprogram not used     integer ik1,ik2,nkk       ! misc. integers
! Subprogram not used 
! Subprogram not used     ! =================================================================
! Subprogram not used     ! Integration matrices of hydrostatic equation(href) and conversion
! Subprogram not used     ! term(a).  href computed as in ccm0 but isothermal bottom ecref
! Subprogram not used     ! calculated to conserve energy
! Subprogram not used     ! =================================================================
! Subprogram not used 
! Subprogram not used     do k=1,plev
! Subprogram not used        do l=1,plev
! Subprogram not used 	  refstate%Href(l,k) = 0.0_r8
! Subprogram not used 	  refstate%Cref(l,k) = 0.0_r8
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     ! ======================================================================
! Subprogram not used     ! Mean atmosphere energy conversion term is consistent with continiuty
! Subprogram not used     ! equation.  In Cref, 1st index = column; 2nd index = row of matrix.
! Subprogram not used     ! Mean atmosphere energy conversion term is energy conserving
! Subprogram not used     ! ======================================================================
! Subprogram not used 
! Subprogram not used     do k=1,plev
! Subprogram not used        refstate%Cref(k,k) = 0.50_r8/refstate%hypm(k) 
! Subprogram not used        do l=1,k-1
! Subprogram not used 	  refstate%Cref(l,k) = 1.0_r8/refstate%hypm(k) 
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     ! ====================================================================
! Subprogram not used     ! Reference hydrostatic integration matrix consistent with conversion
! Subprogram not used     ! term for energy conservation.  In href, 1st index = column; 
! Subprogram not used     ! 2nd index = row of matrix.
! Subprogram not used     ! ====================================================================
! Subprogram not used 
! Subprogram not used     do k = 1,plev
! Subprogram not used        do l = k,plev
! Subprogram not used 	  refstate%Href(l,k) = refstate%Cref(k,l)*refstate%hypd(l)
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     ! ==================================
! Subprogram not used     ! Print statements
! Subprogram not used     ! ==================================
! Subprogram not used 
! Subprogram not used     if (lprint.and.masterproc) then
! Subprogram not used        nkk = plev/13
! Subprogram not used        if (mod(plev,13).ne.0) nkk = nkk + 1
! Subprogram not used        write(iulog,*)' '
! Subprogram not used        write(iulog,*)'INITCOM: Hydrostatic matrix href'
! Subprogram not used        do kk=1,nkk
! Subprogram not used 	  ik1 = 1 + (kk-1)*13
! Subprogram not used 	  ik2 = min0( ik1+12, plev )
! Subprogram not used 	  write(iulog,9920) (k,k=ik1,ik2)
! Subprogram not used 	  do kkk=1,plev
! Subprogram not used 	     write(iulog,9910) kkk,(refstate%Href(kkk,k),k=ik1,ik2)
! Subprogram not used 	  end do
! Subprogram not used        end do
! Subprogram not used        write(iulog,*)' '
! Subprogram not used        write(iulog,*)'INITCOM: Thermodynamic matrix ecref'
! Subprogram not used        do kk=1,nkk
! Subprogram not used 	  ik1 = 1 + (kk-1)*13
! Subprogram not used 	  ik2 = MIN( ik1+12, plev )
! Subprogram not used 	  write(iulog,9920) (k,k=ik1,ik2)
! Subprogram not used 	  do kkk=1,plev
! Subprogram not used 	     write(iulog,9910) kkk,(refstate%Cref(kkk,k),k=ik1,ik2)
! Subprogram not used 	  end do
! Subprogram not used        end do
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! =======================
! Subprogram not used     ! Multiply href by r
! Subprogram not used     ! =======================
! Subprogram not used 
! Subprogram not used     do k=1,plev
! Subprogram not used        do kk=1,plev
! Subprogram not used 	  refstate%Href(kk,k) = refstate%Href(kk,k)*Rgas
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 9910 format( 1x,i3,13f9.5)
! Subprogram not used 9920 format(/,      13i9)
! Subprogram not used 
! Subprogram not used   end subroutine prim_si_refmat_init

  ! =======================================================================
  ! set_vert_struct_mat:
  !
  ! Purpose: 
  ! Set time invariant hydrostatic matrices, which depend on the reference
  ! temperature and pressure in the semi-implicit time step. Note that
  ! this subroutine is actually called twice, because the effective time
  ! step changes between step 0 and step 1.
  !
  ! based on settau. 
  ! =======================================================================

! Subprogram not used   subroutine set_vert_struct_mat(dt, refstate, hvcoord, masterproc)
! Subprogram not used     use parallel_mod, only : abortmp
! Subprogram not used     use hybvcoord_mod, only : hvcoord_t
! Subprogram not used     use linear_algebra_mod, only : reigen_solver
! Subprogram not used     use physical_constants , only : g, rgas, kappa
! Subprogram not used 
! Subprogram not used     !------------------------------Arguments--------------------------------
! Subprogram not used     real(r8), intent(in)         :: dt       ! time step (or dt/2 at time 0)
! Subprogram not used     type(ref_state_t), target    :: refstate
! Subprogram not used     type(hvcoord_t) , target     :: hvcoord
! Subprogram not used     logical, intent(in)          :: masterproc ! master process
! Subprogram not used     !---------------------------Local workspace-----------------------------
! Subprogram not used 
! Subprogram not used     real(r8) zcr(plev)             ! gravity wave equivalent depth
! Subprogram not used     real(r8) zci(plev)             ! dummy, used to print phase speeds
! Subprogram not used 
! Subprogram not used     real(r8) tmp(plev,plev)
! Subprogram not used 
! Subprogram not used     real(r8) rcond
! Subprogram not used     real(r8) z(plev)
! Subprogram not used     real(r8) det(2) 
! Subprogram not used     real(r8) work(plev)
! Subprogram not used     integer ipvt(plev)
! Subprogram not used 
! Subprogram not used     real(r8) Imat(plev,plev)
! Subprogram not used 
! Subprogram not used     real(r8) factor                ! intermediate workspace
! Subprogram not used     real(r8) zdt0u                 ! vertical diff. of ref. temp (above)
! Subprogram not used     real(r8) zshu                  ! interface "sigma" (above)
! Subprogram not used     real(r8) zr2ds                 ! 1./(2.*hypd(k))
! Subprogram not used     real(r8) zdt0d                 ! vertical diff. of ref. temp (below)
! Subprogram not used     real(r8) zshd                  ! interface "sigma" (below)
! Subprogram not used     real(r8) ztd                   ! temporary accumulator
! Subprogram not used     real(r8) zcn                   ! sq(n)
! Subprogram not used 
! Subprogram not used     integer k,l,kk,kkk             ! level indices
! Subprogram not used     integer n,m                    ! n-wavenumber index
! Subprogram not used     integer nneg                   ! number of unstable mean temperatures
! Subprogram not used     integer info
! Subprogram not used     integer ik1,ik2,nkk            ! misc. integers
! Subprogram not used 
! Subprogram not used     real(r8), pointer :: Cref(:,:)
! Subprogram not used     real(r8), pointer :: Tmat(:,:)
! Subprogram not used     real(r8), pointer :: Href(:,:)
! Subprogram not used     real(r8), pointer :: Amat(:,:)
! Subprogram not used     real(r8), pointer :: Amat_inv(:,:)
! Subprogram not used     real(r8), pointer :: Emat(:,:)
! Subprogram not used     real(r8), pointer :: Emat_inv(:,:)
! Subprogram not used 
! Subprogram not used     real(r8), pointer :: RTref(:)
! Subprogram not used     real(r8), pointer :: Tref(:)
! Subprogram not used     real(r8), pointer :: Pvec(:)
! Subprogram not used     real(r8), pointer :: Lambda(:)
! Subprogram not used 
! Subprogram not used     real(r8), pointer :: hypd(:)
! Subprogram not used     real(r8), pointer :: hypi(:)
! Subprogram not used     real(r8), pointer :: hybi(:)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     if (masterproc) then
! Subprogram not used        print *,'Initializing vertical structure matrix'
! Subprogram not used     endif
! Subprogram not used     ! =====================================================
! Subprogram not used     ! Assign pointers to refstate structure for readibility
! Subprogram not used     ! =====================================================
! Subprogram not used 
! Subprogram not used     Cref     => refstate%Cref
! Subprogram not used     Tmat     => refstate%Tmat
! Subprogram not used     Href     => refstate%Href
! Subprogram not used     Amat     => refstate%Amat
! Subprogram not used     Amat_inv => refstate%Amat_inv
! Subprogram not used     Emat     => refstate%Emat
! Subprogram not used     Emat_inv => refstate%Emat_inv
! Subprogram not used     Tref     => refstate%Tref
! Subprogram not used     Pvec     => refstate%Pvec
! Subprogram not used     RTref    => refstate%RTref
! Subprogram not used     Lambda   => refstate%Lambda
! Subprogram not used 
! Subprogram not used     hybi => hvcoord%hybi
! Subprogram not used     hypd => refstate%hypd
! Subprogram not used     hypi => refstate%hypi
! Subprogram not used 
! Subprogram not used     ! =========================================
! Subprogram not used     ! Calculate hydrostatic matrix tau (Tmat)
! Subprogram not used     ! concordance with subroutine settau from 1:
! Subprogram not used     ! --------------------------------------------
! Subprogram not used     !   here     settau 
! Subprogram not used     !   ---------------
! Subprogram not used     !   Cref  == ecref
! Subprogram not used     !   Href  == href
! Subprogram not used     !   Tref  == t0
! Subprogram not used     !   RTref == bps
! Subprogram not used     !   Tmat  == tau
! Subprogram not used     !   Amat  == zb
! Subprogram not used     !
! Subprogram not used     ! =========================================
! Subprogram not used 
! Subprogram not used     ! ===========================================================================
! Subprogram not used     ! This formula for Tmat (assumes constant (in vertical) reference Temperature
! Subprogram not used     ! ===========================================================================
! Subprogram not used 
! Subprogram not used     do k=1,plev
! Subprogram not used        do l=1,k
! Subprogram not used 	  Tmat(l,k) = kappa*Tref(k)*Cref(l,k)*hypd(l)
! Subprogram not used        end do
! Subprogram not used        do l=k+1,plev
! Subprogram not used 	  Tmat(l,k) = 0.0_r8
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     ! ===============================================================
! Subprogram not used     ! Vector for linear surface pressure term in divergence
! Subprogram not used     ! Pressure gradient and diagonal term of hydrostatic components
! Subprogram not used     ! ===============================================================
! Subprogram not used     do k=1,plev
! Subprogram not used        RTref(k) = Rgas*Tref(k)
! Subprogram not used        Pvec(k)  = hypd(k)/hypi(plevp)
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     do k=1,plev
! Subprogram not used        do l=1,plev
! Subprogram not used 	  ztd = RTref(k) * Pvec(l)
! Subprogram not used 	  do kkk=1,plev
! Subprogram not used 	     ztd = ztd + Href(kkk,k)*Tmat(l,kkk)
! Subprogram not used 	  end do
! Subprogram not used 	  Amat(l,k) = ztd
! Subprogram not used 	  !       Amat(k,l) = ztd
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     tmp(:,:)           = TRANSPOSE(Amat(:,:))
! Subprogram not used     Amat(:,:) = tmp(:,:)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     ! =========================================
! Subprogram not used     ! invert the vertical structure matrix
! Subprogram not used     ! Amat_inv in Steve Thomas's notation.
! Subprogram not used     ! =========================================
! Subprogram not used 
! Subprogram not used     do k=1,plev
! Subprogram not used        do l=1,plev
! Subprogram not used 	  Amat_inv(k,l) = Amat(k,l)
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     call abortmp('not supported in cam')
! Subprogram not used 
! Subprogram not used     ! =======================================================================
! Subprogram not used     ! compute the right eigenvector matrix (Emat) and its inverse (Emat_inv)
! Subprogram not used     ! from the vertical structure matrix.
! Subprogram not used     ! =======================================================================
! Subprogram not used     Emat=Amat
! Subprogram not used     info=reigen_solver(plev,Emat,Emat_Inv, zcr)
! Subprogram not used     
! Subprogram not used     do k=1,plev
! Subprogram not used        do l=1,plev
! Subprogram not used 	  Emat_inv(k,l) = Emat(k,l)
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used     call abortmp('not supported in cam')
! Subprogram not used     tmp(:,:)      = TRANSPOSE(Emat(:,:))
! Subprogram not used     Emat(:,:)     = tmp(:,:)
! Subprogram not used     tmp(:,:)      = TRANSPOSE(Amat(:,:))
! Subprogram not used     Amat(:,:)     = tmp(:,:)
! Subprogram not used 
! Subprogram not used     tmp(:,:)      = TRANSPOSE(Emat_inv(:,:))
! Subprogram not used     Emat_inv(:,:) = tmp(:,:)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     ! =======================================================================
! Subprogram not used     ! Compute and print gravity wave equivalent depths and phase speeds
! Subprogram not used     ! co
! Subprogram not used     ! =======================================================================
! Subprogram not used 
! Subprogram not used     do k=1,plev
! Subprogram not used        zci(k) = sign(1.0_r8,zcr(k))*sqrt(abs(zcr(k)))
! Subprogram not used        Lambda(k) = zcr(k)*dt*dt               ! solver eigenvalues
! Subprogram not used        !    Lambda(k) = zcr(k)                     ! solver eigenvalues
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     if (masterproc) then
! Subprogram not used        write(iulog,910) (Tref(k),k=1,plev)
! Subprogram not used        write(iulog,920) (zci(k),k=1,plev)
! Subprogram not used        write(iulog,930) (zcr(k)/g,k=1,plev)
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! ========================================================================
! Subprogram not used     ! Test for unstable mean temperatures (negative phase speed and eqivalent
! Subprogram not used     ! depth) for at least one gravity wave.
! Subprogram not used     ! ========================================================================
! Subprogram not used 
! Subprogram not used     nneg = 0
! Subprogram not used     do k=1,plev
! Subprogram not used        if (zcr(k)<=0.0_r8) nneg = nneg + 1
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     if (nneg/=0) then
! Subprogram not used        write(iulog,*)'---------------------------------------------------'
! Subprogram not used        write(iulog,*)'SET_VERT_STRUCT_MAT: UNSTABLE MEAN TEMPERATURE. STOP',zcr
! Subprogram not used        call abortmp('prim_si_ref_mod')
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used 910 format(' REFERENCE TEMPERATURES FOR SEMI-IMPLICIT SCHEME = ',  /(1x,12f9.3))
! Subprogram not used 920 format(' GRAVITY WAVE PHASE SPEEDS (M/S) FOR MEAN STATE = '    /(1x,12f9.3))
! Subprogram not used 930 format(' GRAVITY WAVE EQUIVALENT DEPTHS (M) FOR MEAN STATE = ' /(1x,12f9.3))
! Subprogram not used 
! Subprogram not used 9910 format( 1x,i3,13f9.5)
! Subprogram not used 9920 format(/,      13i9)
! Subprogram not used 
! Subprogram not used   end subroutine set_vert_struct_mat


  subroutine prim_set_mass(elem, tl,hybrid,hvcoord,nets,nete)
  use kinds, only : real_kind
  use control_mod, only : initial_total_mass
  use physical_constants, only : g
  use element_mod, only : element_t
  use time_mod, only : timelevel_t 
  use hybvcoord_mod, only : hvcoord_t 
  use hybrid_mod, only : hybrid_t
  use dimensions_mod, only : np
  use global_norms_mod, only : global_integral 

  type (element_t), intent(inout) :: elem(:)
  type (TimeLevel_t), target, intent(in) :: tl
  type (hybrid_t),intent(in)     :: hybrid
  type (hvcoord_t), intent(in)   :: hvcoord
  integer,intent(in)             :: nets,nete
  
  ! local 
  real (kind=real_kind)  :: tmp(np,np,nets:nete)
  real (kind=real_kind)  :: scale,mass0
  integer :: n0,nm1,np1,ie

  if (initial_total_mass == 0) return;
  
  n0=tl%n0
  nm1=tl%nm1
  np1=tl%np1
  
  scale=1/g                                  ! assume code is using Pa
  if (hvcoord%ps0 <  2000 ) scale=100*scale  ! code is using mb
  ! after scaling, Energy is in J/m**2,  Mass kg/m**2
  
  do ie=nets,nete
     tmp(:,:,ie)=elem(ie)%state%ps_v(:,:,n0)
  enddo
  mass0 = global_integral(elem, tmp(:,:,nets:nete),hybrid,np,nets,nete)
  mass0 = mass0*scale;  
  
  do ie=nets,nete
     elem(ie)%state%ps_v(:,:,n0)=elem(ie)%state%ps_v(:,:,n0)*(initial_total_mass/mass0)
     elem(ie)%state%ps_v(:,:,np1)=elem(ie)%state%ps_v(:,:,n0)
     elem(ie)%state%ps_v(:,:,nm1)=elem(ie)%state%ps_v(:,:,n0)
  enddo
  if(hybrid%par%masterproc .and. hybrid%ithr==0) then 
     write (*,'(a,e24.15)') "Initializing Total Mass (kg/m^2) = ",initial_total_mass
  endif
  end subroutine prim_set_mass


end module prim_si_ref_mod
