!#define _DBG_ print *,"File:","<stdin>"," at ",5
!#define _DBG_ !DBG
!
!
module prim_advance_mod
  use edge_mod, only : EdgeBuffer_t
  use kinds, only : real_kind, iulog
  use perf_mod, only: t_startf, t_stopf, t_barrierf ! _EXTERNAL
  use parallel_mod, only : abortmp
  implicit none
  private
  save
  public :: prim_advance_exp, prim_advance_si, prim_advance_init, preq_robert3,&
       applyCAMforcing_dynamics, applyCAMforcing, smooth_phis, overwrite_SEdensity

  type (EdgeBuffer_t) :: edge1
  type (EdgeBuffer_t) :: edge2
  type (EdgeBuffer_t) :: edge3p1

  real (kind=real_kind) :: initialized_for_dt   = 0

  real (kind=real_kind), allocatable :: ur_weights(:)


contains

  subroutine prim_advance_init(integration)
    use edge_mod, only : initEdgeBuffer
    use dimensions_mod, only : nlev
    use control_mod, only : qsplit,rsplit
    character(len=*)    , intent(in) :: integration 
    integer :: i

    if (rsplit==0) then
       call initEdgeBuffer(edge3p1,3*nlev+1)
    else   
       ! need extra buffer space for dp3d
       call initEdgeBuffer(edge3p1,4*nlev+1)
    endif

    if(integration == 'semi_imp') then
       call initEdgeBuffer(edge1,nlev)
       call initEdgeBuffer(edge2,2*nlev)
    end if

    ! compute averaging weights for RK+LF (tstep_type=1) timestepping:
    allocate(ur_weights(qsplit))
    ur_weights(:)=0.0d0

    if(mod(qsplit,2).NE.0)then
       ur_weights(1)=1.0d0/qsplit
       do i=3,qsplit,2
         ur_weights(i)=2.0d0/qsplit
       enddo
    else
       do i=2,qsplit,2
         ur_weights(i)=2.0d0/qsplit
       enddo
    endif

  end subroutine prim_advance_init


  subroutine prim_advance_exp(elem, deriv, hvcoord, hybrid,&
       dt, tl,  nets, nete, compute_diagnostics)
    use bndry_mod, only : bndry_exchangev
    use control_mod, only : prescribed_wind, qsplit, tstep_type, rsplit, qsplit, moisture, integration
    use derivative_mod, only : derivative_t, vorticity, divergence, gradient, gradient_wk
    use dimensions_mod, only : np, nlev, nlevp, nvar
!    use prim_state_mod, only : prim_printstate
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack, initEdgeBuffer
    use element_mod, only : element_t
    use hybvcoord_mod, only : hvcoord_t
    use hybrid_mod, only : hybrid_t
    use reduction_mod, only : reductionbuffer_ordered_1d_t
    use time_mod, only : TimeLevel_t,  timelevel_qdp
    use diffusion_mod, only :  prim_diffusion


    implicit none

    type (element_t), intent(inout), target   :: elem(:)
    type (derivative_t), intent(in)   :: deriv
    type (hvcoord_t)                  :: hvcoord

    type (hybrid_t)    , intent(in):: hybrid

    real (kind=real_kind), intent(in) :: dt
    type (TimeLevel_t)   , intent(in) :: tl
    integer              , intent(in) :: nets
    integer              , intent(in) :: nete
    logical, intent(in)               :: compute_diagnostics

    ! =================
    ! Local
    ! =================
    real (kind=real_kind) ::  dt2, time, dt_vis, eta_ave_w
    real (kind=real_kind) ::  eta_dot_dpdn(np,np,nlevp)
    real (kind=real_kind) ::  dp(np,np)
    integer :: ie,nm1,n0,np1,nstep,method,qsplit_stage,k, qn0
    integer :: n,i,j,lx,lenx


    call t_barrierf('sync_prim_advance_exp', hybrid%par%comm)
    call t_startf('prim_advance_exp')
    nm1   = tl%nm1
    n0    = tl%n0
    np1   = tl%np1
    nstep = tl%nstep

    ! timelevel to use for accessing Qdp() to compute virtual temperature
    qn0 = -1    ! -1 = disabled (assume dry dynamics)
    if ( moisture /= "dry") then
       call TimeLevel_Qdp(tl, qsplit, qn0)  ! compute current Qdp() timelevel
    endif
     

! integration = "explicit"
!
!   tstep_type=0  pure leapfrog except for very first timestep   CFL=1
!                    typically requires qsplit=4 or 5 
!   tstep_type=1  RK2 followed by qsplit-1 leapfrog steps        CFL=close to qsplit
!                    typically requires qsplit=4 or 5 
!   tstep_type=2  RK2-SSP 3 stage (as used by tracers)           CFL=.58  
!                    optimal in terms of SSP CFL, but not        CFLSSP=2 
!                    optimal in terms of CFL
!                    typically requires qsplit=3
!                    but if windspeed > 340m/s, could use this
!                    with qsplit=1
!   tstep_type=3  classic RK3                                    CFL=1.73 (sqrt(3))
!                 
!   tstep_type=4  Kinnmark&Gray RK4 4 stage                      CFL=sqrt(8)=2.8
!                 should we replace by standard RK4 (CFL=sqrt(8))?
!                 (K&G 1st order method has CFL=3)
!   tstep_type=5  Kinnmark&Gray RK3 5 stage 3rd order            CFL=3.87  (sqrt(15))
!                 From Paul Ullrich.  3rd order for nonlinear terms also
!                 K&G method is only 3rd order for linear
!                 optimal: for windspeeds ~120m/s,gravity: 340m/2
!                 run with qsplit=1
!                 (K&G 2nd order method has CFL=4. tiny CFL improvement not worth 2nd order)
!
! integration = "full_imp"  
!
!   tstep_type=1  Backward Euler, first order  (under development)
!
!   tstep_type=2  Crank Nicolson, second order (under development)
!

! default weights for computing mean dynamics fluxes
    eta_ave_w = 1d0/qsplit   

    if(tstep_type==0)then  
       method=0                ! pure leapfrog
       if (nstep==0) method=1  ! but use RK2 on first step
    else if (tstep_type==1) then  
       method=0                           ! LF
       qsplit_stage = mod(nstep,qsplit)
       if (qsplit_stage==0) method=1      ! RK2 on first of qsplit steps
       ! RK2 + LF scheme has tricky weights:
       eta_ave_w=ur_weights(qsplit_stage+1) 
    else 
       method = tstep_type                ! other RK variants 
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    ! fix dynamical variables, skip dynamics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    if (1==prescribed_wind) then
       time=tl%nstep*dt
       do ie=nets,nete
          ! read in omega met data and store it in derived%omega_prescribed
          ! --- compute eta_dot_dpdn from omega_prescribed ...
          eta_dot_dpdn(:,:,nlev+1) = 0.0d0

          do k = nlev,2,-1
             eta_dot_dpdn(:,:,k) = eta_dot_dpdn(:,:,k+1) - 2.d0*elem(ie)%derived%omega_prescribed(:,:,k)
          enddo

          eta_dot_dpdn(:,:,1) = 0.0d0
          ! accumulate mean fluxes for advection
          elem(ie)%derived%eta_dot_dpdn(:,:,:) = &
               elem(ie)%derived%eta_dot_dpdn(:,:,:) + eta_dot_dpdn(:,:,:)*eta_ave_w

          ! subcycling code uses a mean flux to advect tracers 
          do k=1,nlev
             dp(:,:) =&
                  ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                  ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,tl%n0) 
             
             elem(ie)%derived%vn0(:,:,1,k)=elem(ie)%derived%vn0(:,:,1,k)+&
                  eta_ave_w*elem(ie)%state%v(:,:,1,k,n0)*dp(:,:)
             elem(ie)%derived%vn0(:,:,2,k)=elem(ie)%derived%vn0(:,:,2,k)+&
                  eta_ave_w*elem(ie)%state%v(:,:,2,k,n0)*dp(:,:)
          enddo
       end do
       call t_stopf('prim_advance_exp')
       return
    endif


    ! ==================================
    ! Take timestep
    ! ==================================
    dt_vis = dt                      
    if (method==0) then
       ! regular LF step
       dt2 = 2*dt
       call compute_and_apply_rhs(np1,nm1,n0,qn0,dt2,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,eta_ave_w)
       dt_vis = dt2  ! dt to use for time-split dissipation
    else if (method==1) then
       ! RK2
       ! forward euler to u(dt/2) = u(0) + (dt/2) RHS(0)  (store in u(np1))
       call compute_and_apply_rhs(np1,n0,n0,qn0,dt/2,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,0d0)
       ! leapfrog:  u(dt) = u(0) + dt RHS(dt/2)     (store in u(np1))
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,eta_ave_w)
    else if (method==2) then
       ! RK2-SSP 3 stage.  matches tracer scheme. optimal SSP CFL, but
       ! not optimal for regular CFL
       ! u1 = u0 + dt/2 RHS(u0)
       call compute_and_apply_rhs(np1,n0,n0,qn0,dt/2,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,eta_ave_w/3)
       ! u2 = u1 + dt/2 RHS(u1)
       call compute_and_apply_rhs(np1,np1,np1,qn0,dt/2,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,eta_ave_w/3)
       ! u3 = u2 + dt/2 RHS(u2)
       call compute_and_apply_rhs(np1,np1,np1,qn0,dt/2,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,eta_ave_w/3)
       ! unew = u/3 +2*u3/3  = u + 1/3 (RHS(u) + RHS(u1) + RHS(u2))
       do ie=nets,nete
          elem(ie)%state%v(:,:,:,:,np1)= elem(ie)%state%v(:,:,:,:,n0)/3 &
               + 2*elem(ie)%state%v(:,:,:,:,np1)/3
          elem(ie)%state%T(:,:,:,np1)= elem(ie)%state%T(:,:,:,n0)/3 &
               + 2*elem(ie)%state%T(:,:,:,np1)/3
          if (rsplit==0) then
             elem(ie)%state%ps_v(:,:,np1)= elem(ie)%state%ps_v(:,:,n0)/3 &
                  + 2*elem(ie)%state%ps_v(:,:,np1)/3
          else
             elem(ie)%state%dp3d(:,:,:,np1)= elem(ie)%state%dp3d(:,:,:,n0)/3 &
                  + 2*elem(ie)%state%dp3d(:,:,:,np1)/3
          endif
       enddo
    else if (method==3) then
       ! classic RK3  CFL=sqrt(3)  
       ! u1 = u0 + dt/3 RHS(u0)
       call compute_and_apply_rhs(np1,n0,n0,qn0,dt/3,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,0d0)
       ! u2 = u0 + dt/2 RHS(u1)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt/2,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0)
       ! u3 = u0 + dt RHS(u2)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,eta_ave_w)
    else if (method==4) then
       ! KG 4th order 4 stage:   CFL=sqrt(8)
       ! low storage version of classic RK4
       ! u1 = u0 + dt/4 RHS(u0)
       call compute_and_apply_rhs(np1,n0,n0,qn0,dt/4,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,0d0)
       ! u2 = u0 + dt/3 RHS(u1)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt/3,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0)
       ! u3 = u0 + dt/2 RHS(u2)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt/2,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0)
       ! u4 = u0 + dt RHS(u3)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,eta_ave_w)
    else if (method==5) then
       ! Ullrich 3nd order 5 stage:   CFL=sqrt( 4^2 -1) = 3.87
       ! u1 = u0 + dt/5 RHS(u0)  (save u1 in timelevel nm1)
       call compute_and_apply_rhs(nm1,n0,n0,qn0,dt/5,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,eta_ave_w/4)
       ! u2 = u0 + dt/5 RHS(u1)
       call compute_and_apply_rhs(np1,n0,nm1,qn0,dt/5,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0)
       ! u3 = u0 + dt/3 RHS(u2)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt/3,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0)
       ! u4 = u0 + 2dt/3 RHS(u3)
       call compute_and_apply_rhs(np1,n0,np1,qn0,2*dt/3,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0)

       ! compute (5*u1/4 - u0/4) in timelevel nm1:
       do ie=nets,nete
          elem(ie)%state%v(:,:,:,:,nm1)= (5*elem(ie)%state%v(:,:,:,:,nm1) &
               - elem(ie)%state%v(:,:,:,:,n0) ) /4
          elem(ie)%state%T(:,:,:,nm1)= (5*elem(ie)%state%T(:,:,:,nm1) &
               - elem(ie)%state%T(:,:,:,n0) )/4
          if (rsplit==0) then
             elem(ie)%state%ps_v(:,:,nm1)= ( 5*elem(ie)%state%ps_v(:,:,nm1) &
                  - elem(ie)%state%ps_v(:,:,n0) )/4
          else
             elem(ie)%state%dp3d(:,:,:,nm1)= (5*elem(ie)%state%dp3d(:,:,:,nm1) &
                  - elem(ie)%state%dp3d(:,:,:,n0) )/4
          endif
       enddo
       ! u5 = (5*u1/4 - u0/4) + 3dt/4 RHS(u4)
       call compute_and_apply_rhs(np1,nm1,np1,qn0,3*dt/4,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,3*eta_ave_w/4)
       ! final method is the same as:
       ! u5 = u0 +  dt/4 RHS(u0)) + 3dt/4 RHS(u4)

    else if (method==11) then
       ! Backward Euler fully implicit JFNK method (vertically langragian not active yet)
       if (rsplit > 0) then  
       call abortmp('ERROR: full_imp integration not yet coded for vert lagrangian adv option')
       end if
      if (hybrid%masterthread) print*, "fully implicit integration is still under development"


    else
       call abortmp('ERROR: bad choice of tstep_type')
    endif

!    call prim_printstate(elem,tl,hybrid,hvcoord,nets,nete)

    ! ==============================================
    ! Time-split Horizontal diffusion: nu.del^2 or nu.del^4
    ! U(*) = U(t+1)  + dt2 * HYPER_DIFF_TERM(t+1)
    ! ==============================================

    ! note:time step computes u(t+1)= u(t*) + RHS. 
    ! for consistency, dt_vis = t-1 - t*, so this is timestep method dependent
    if (tstep_type==0) then  
       ! leapfrog special case
       call advance_hypervis_lf(edge3p1,elem,hvcoord,hybrid,deriv,nm1,n0,np1,nets,nete,dt_vis)
    else if (tstep_type==11) then  
       ! add later once incorporated into implicit operator
!       call advance_hypervis_lf(edge3p1,elem,hvcoord,hybrid,deriv,nm1,n0,np1,nets,nete,dt_vis)
    else
       if (rsplit==0) then
          ! forward-in-time, maybe hypervis applied to PS
          call advance_hypervis(edge3p1,elem,hvcoord,hybrid,deriv,np1,nets,nete,dt_vis,eta_ave_w)
       else
          ! forward-in-time, hypervis applied to dp3d
          call advance_hypervis_dp(edge3p1,elem,hvcoord,hybrid,deriv,np1,nets,nete,dt_vis,eta_ave_w)
       endif
    endif


    call t_stopf('prim_advance_exp')
    end subroutine prim_advance_exp


! Subprogram not used subroutine prim_advance_si(elem, nets, nete, cg, blkjac, red, &
! Subprogram not used           refstate, hvcoord, deriv, flt, hybrid, tl, dt)
! Subprogram not used        use bndry_mod, only : bndry_exchangev
! Subprogram not used        use cg_mod, only : cg_t, cg_create
! Subprogram not used        use control_mod, only : filter_freq,debug_level, precon_method
! Subprogram not used        use derivative_mod, only : derivative_t, vorticity, divergence, gradient, gradient_wk
! Subprogram not used        use dimensions_mod, only : np, nlev, nlevp
! Subprogram not used        use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack, initEdgeBuffer
! Subprogram not used        use element_mod, only : element_t
! Subprogram not used        use filter_mod, only : filter_t, preq_filter
! Subprogram not used        use hybvcoord_mod, only : hvcoord_t
! Subprogram not used        use hybrid_mod, only : hybrid_t
! Subprogram not used        use prim_si_ref_mod, only : ref_state_t, set_vert_struct_mat
! Subprogram not used        use reduction_mod, only : reductionbuffer_ordered_1d_t
! Subprogram not used        use solver_mod, only : pcg_solver, blkjac_t, blkjac_init
! Subprogram not used        use time_mod, only : TimeLevel_t
! Subprogram not used        use prim_si_mod, only : preq_vertadv, preq_omegap, preq_pressure
! Subprogram not used        use diffusion_mod, only :  prim_diffusion
! Subprogram not used        use physical_constants, only : kappa, rrearth, rgas, cp, rwater_vapor
! Subprogram not used        use physics_mod, only : virtual_temperature, virtual_specific_heat
! Subprogram not used        implicit none
! Subprogram not used 
! Subprogram not used        integer, intent(in)               :: nets,nete
! Subprogram not used        type (element_t), intent(inout), target :: elem(:)
! Subprogram not used        type (blkjac_t)                   :: blkjac(nets:nete)
! Subprogram not used 
! Subprogram not used        type (cg_t)                       :: cg
! Subprogram not used 
! Subprogram not used        type (ReductionBuffer_ordered_1d_t), intent(inout) :: red
! Subprogram not used 
! Subprogram not used        type (ref_state_t), intent(in), target :: refstate
! Subprogram not used        type (hvcoord_t), intent(in)      :: hvcoord
! Subprogram not used        type (derivative_t), intent(in)   :: deriv
! Subprogram not used        type (filter_t), intent(in)       :: flt
! Subprogram not used        type (hybrid_t), intent(in)       :: hybrid
! Subprogram not used        type (TimeLevel_t), intent(in)    :: tl
! Subprogram not used        real(kind=real_kind), intent(in)  :: dt
! Subprogram not used        real(kind=real_kind)              :: time_adv
! Subprogram not used        ! ==========================
! Subprogram not used        ! Local variables...
! Subprogram not used        ! ==========================
! Subprogram not used 
! Subprogram not used        real(kind=real_kind)                           :: ps0
! Subprogram not used        real(kind=real_kind)                           :: psref
! Subprogram not used 
! Subprogram not used        real(kind=real_kind), dimension(np,np)         :: ps
! Subprogram not used        real(kind=real_kind), dimension(np,np)         :: rps
! Subprogram not used        real(kind=real_kind), dimension(np,np,nlev)    :: rpmid
! Subprogram not used        real(kind=real_kind), dimension(np,np,nlev)    :: omegap
! Subprogram not used        real(kind=real_kind), dimension(np,np,nlev)    :: rpdel
! Subprogram not used 
! Subprogram not used        real(kind=real_kind) :: pintref(nlevp)
! Subprogram not used        real(kind=real_kind) :: pdelref(nlev)
! Subprogram not used        real(kind=real_kind) :: pmidref(nlev)
! Subprogram not used        real(kind=real_kind) :: rpdelref(nlev)
! Subprogram not used        real(kind=real_kind) :: rpmidref(nlev)
! Subprogram not used 
! Subprogram not used        real(kind=real_kind) :: pint(np,np,nlevp)
! Subprogram not used        real(kind=real_kind) :: pdel(np,np,nlev)
! Subprogram not used        real(kind=real_kind) :: pmid(np,np,nlev)
! Subprogram not used 
! Subprogram not used        real(kind=real_kind), dimension(np,np,nlevp) :: eta_dot_dp_deta
! Subprogram not used        real(kind=real_kind), dimension(np,np,nlev)  :: vgrad_ps
! Subprogram not used 
! Subprogram not used        real(kind=real_kind), dimension(np,np,nlev)   :: T_vadv
! Subprogram not used        real(kind=real_kind), dimension(np,np,2,nlev) :: v_vadv
! Subprogram not used 
! Subprogram not used        real(kind=real_kind), dimension(np,np)      :: HT
! Subprogram not used        real(kind=real_kind), dimension(np,np)      :: HrefT
! Subprogram not used        real(kind=real_kind), dimension(np,np)      :: HrefTm1
! Subprogram not used 
! Subprogram not used        real(kind=real_kind), dimension(np,np)      :: Gref0
! Subprogram not used        real(kind=real_kind), dimension(np,np)      :: Grefm1
! Subprogram not used        real(kind=real_kind), dimension(np,np)      :: E
! Subprogram not used        real(kind=real_kind), dimension(np,np)      :: Phi
! Subprogram not used        real(kind=real_kind), dimension(np,np)      :: dGref
! Subprogram not used 
! Subprogram not used        real(kind=real_kind), dimension(np,np,2)    :: vco
! Subprogram not used        real(kind=real_kind), dimension(np,np,2)    :: gradT
! Subprogram not used        real(kind=real_kind), dimension(np,np,2)    :: grad_Phi
! Subprogram not used 
! Subprogram not used        real(kind=real_kind), dimension(:,:), pointer  :: Emat
! Subprogram not used        real(kind=real_kind), dimension(:,:), pointer  :: Emat_inv
! Subprogram not used        real(kind=real_kind), dimension(:,:), pointer  :: Amat
! Subprogram not used        real(kind=real_kind), dimension(:,:), pointer  :: Amat_inv
! Subprogram not used        real(kind=real_kind), dimension(:), pointer    :: Lambda
! Subprogram not used 
! Subprogram not used        real(kind=real_kind), dimension(:), pointer    :: Tref
! Subprogram not used        real(kind=real_kind), dimension(:), pointer    :: RTref
! Subprogram not used        real(kind=real_kind), dimension(:), pointer    :: Pvec
! Subprogram not used        real(kind=real_kind), dimension(:,:), pointer  :: Href
! Subprogram not used        real(kind=real_kind), dimension(:,:), pointer  :: Tmat
! Subprogram not used 
! Subprogram not used        real(kind=real_kind) :: Vscript(np,np,2,nlev,nets:nete)
! Subprogram not used        real(kind=real_kind) :: Tscript(np,np,nlev,nets:nete)
! Subprogram not used        real(kind=real_kind) :: Pscript(np,np,nets:nete)
! Subprogram not used 
! Subprogram not used        real(kind=real_kind), dimension(np,np)      :: HrefTscript
! Subprogram not used        real(kind=real_kind), dimension(np,np)      :: suml
! Subprogram not used        real(kind=real_kind), dimension(np,np,2)    :: gVscript
! Subprogram not used        real(kind=real_kind), dimension(np,np,nlev) :: div_Vscript
! Subprogram not used 
! Subprogram not used        real(kind=real_kind) :: B(np,np,nlev,nets:nete)
! Subprogram not used        real(kind=real_kind) :: C(np,np,nlev,nets:nete)
! Subprogram not used        real(kind=real_kind) :: D(np,np,nlev,nets:nete)
! Subprogram not used 
! Subprogram not used        real(kind=real_kind) :: Gamma_ref(np,np,nlev,nets:nete)
! Subprogram not used 
! Subprogram not used        real(kind=real_kind) :: Gref(np,np,nlev,nets:nete)
! Subprogram not used        real(kind=real_kind) :: grad_dGref(np,np,2,nlev)
! Subprogram not used        real(kind=real_kind) :: grad_Gref(np,np,2,nlev)
! Subprogram not used 
! Subprogram not used        real(kind=real_kind) :: div(np,np)
! Subprogram not used        real(kind=real_kind) :: gv(np,np,2)
! Subprogram not used 
! Subprogram not used        real(kind=real_kind) :: dt2
! Subprogram not used        real(kind=real_kind) :: rpsref
! Subprogram not used        real(kind=real_kind) :: rdt
! Subprogram not used        real(kind=real_kind) :: hkk, hkl
! Subprogram not used        real(kind=real_kind) :: ddiv
! Subprogram not used 
! Subprogram not used        real(kind=real_kind) :: vgradT
! Subprogram not used        real(kind=real_kind) :: hybfac
! Subprogram not used        real(kind=real_kind) :: Crkk
! Subprogram not used        real(kind=real_kind) :: v1,v2
! Subprogram not used        real(kind=real_kind) :: term
! Subprogram not used 
! Subprogram not used        real(kind=real_kind) :: Vs1,Vs2
! Subprogram not used        real(kind=real_kind) :: glnps1, glnps2
! Subprogram not used        real(kind=real_kind) :: gGr1,gGr2
! Subprogram not used 
! Subprogram not used        real (kind=real_kind),allocatable :: solver_wts(:,:)  ! solver weights array for nonstaggered grid
! Subprogram not used 
! Subprogram not used        integer              :: nm1,n0,np1,nfilt
! Subprogram not used        integer              :: nstep
! Subprogram not used        integer              :: i,j,k,l,ie,kptr
! Subprogram not used 
! Subprogram not used        call t_barrierf('sync_prim_advance_si', hybrid%par%comm)
! Subprogram not used        call t_startf('prim_advance_si')
! Subprogram not used 
! Subprogram not used        nm1   = tl%nm1
! Subprogram not used        n0    = tl%n0
! Subprogram not used        np1   = tl%np1
! Subprogram not used        nstep = tl%nstep
! Subprogram not used 
! Subprogram not used 
! Subprogram not used        if ( dt /= initialized_for_dt ) then
! Subprogram not used           if(hybrid%par%masterproc) print *,'Initializing semi-implicit matricies for dt=',dt
! Subprogram not used 
! Subprogram not used           !$OMP MASTER
! Subprogram not used           call set_vert_struct_mat(dt, refstate, hvcoord, hybrid%masterthread)
! Subprogram not used           !$OMP END MASTER
! Subprogram not used 
! Subprogram not used           allocate(solver_wts(np*np,nete-nets+1))
! Subprogram not used           do ie=nets,nete
! Subprogram not used              kptr=1
! Subprogram not used              do j=1,np
! Subprogram not used                 do i=1,np
! Subprogram not used 
! Subprogram not used                    ! so this code is BFB  with old code.  should change to simpler formula below
! Subprogram not used                    solver_wts(kptr,ie-nets+1) = 1d0/nint(1d0/(elem(ie)%mp(i,j)*elem(ie)%rmp(i,j)))
! Subprogram not used                    !solver_wts(kptr,ie-nets+1) = elem(ie)%mp(i,j)*elem(ie)%rmp(i,j)
! Subprogram not used 
! Subprogram not used                    kptr=kptr+1
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used           end do
! Subprogram not used           call cg_create(cg, np*np, nlev, nete-nets+1, hybrid, debug_level, solver_wts)
! Subprogram not used           deallocate(solver_wts)
! Subprogram not used           if (precon_method == "block_jacobi") then
! Subprogram not used              call blkjac_init(elem, deriv,refstate%Lambda,nets,nete,blkjac)
! Subprogram not used           end if
! Subprogram not used           initialized_for_dt = dt
! Subprogram not used        endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used        nfilt = tl%nm1     ! time level at which filter is applied (time level n)
! Subprogram not used        dt2   = 2.0_real_kind*dt
! Subprogram not used        rdt   = 1.0_real_kind/dt
! Subprogram not used 
! Subprogram not used        ps0      = hvcoord%ps0
! Subprogram not used        psref    = refstate%psr
! Subprogram not used 
! Subprogram not used        Emat     => refstate%Emat
! Subprogram not used        Emat_inv => refstate%Emat_inv
! Subprogram not used        Amat     => refstate%Amat
! Subprogram not used        Amat_inv => refstate%Amat_inv
! Subprogram not used        Lambda   => refstate%Lambda
! Subprogram not used 
! Subprogram not used        RTref    => refstate%RTref
! Subprogram not used        Tref     => refstate%Tref
! Subprogram not used        Href     => refstate%Href
! Subprogram not used        Tmat     => refstate%Tmat
! Subprogram not used        Pvec     => refstate%Pvec
! Subprogram not used 
! Subprogram not used        ! ============================================================
! Subprogram not used        ! If the time is right, apply a filter to the state variables
! Subprogram not used        ! ============================================================
! Subprogram not used 
! Subprogram not used        if (nstep > 0 .and. filter_freq > 0 .and. MODULO(nstep,filter_freq) == 0 ) then
! Subprogram not used           call preq_filter(elem, edge3p1, flt, cg%hybrid, nfilt, nets, nete)
! Subprogram not used        end if
! Subprogram not used 
! Subprogram not used        do ie = nets, nete
! Subprogram not used 
! Subprogram not used           elem(ie)%derived%grad_lnps(:,:,:) = gradient(elem(ie)%state%lnps(:,:,n0),deriv)*rrearth
! Subprogram not used 
! Subprogram not used        end do
! Subprogram not used 
! Subprogram not used        ! ================================================
! Subprogram not used        ! boundary exchange grad_lnps
! Subprogram not used        ! ================================================
! Subprogram not used 
! Subprogram not used        do ie = nets, nete
! Subprogram not used 
! Subprogram not used           do k=1,nlevp
! Subprogram not used              pintref(k)  = hvcoord%hyai(k)*ps0 + hvcoord%hybi(k)*psref
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used           do k=1,nlev
! Subprogram not used              pmidref(k)  = hvcoord%hyam(k)*ps0 + hvcoord%hybm(k)*psref
! Subprogram not used              pdelref(k)  = pintref(k+1) - pintref(k)
! Subprogram not used              rpmidref(k) = 1.0_real_kind/pmidref(k)
! Subprogram not used              rpdelref(k) = 1.0_real_kind/pdelref(k)
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used           rpsref   = 1.0_real_kind/psref
! Subprogram not used 
! Subprogram not used           ps(:,:) = EXP(elem(ie)%state%lnps(:,:,n0))
! Subprogram not used           rps(:,:) = 1.0_real_kind/ps(:,:)
! Subprogram not used 
! Subprogram not used           call preq_pressure(ps0,ps,hvcoord%hyai,hvcoord%hybi,hvcoord%hyam,hvcoord%hybm,pint,pmid,pdel)
! Subprogram not used 
! Subprogram not used           rpmid = 1.0_real_kind/pmid
! Subprogram not used           rpdel = 1.0_real_kind/pdel
! Subprogram not used 
! Subprogram not used           do k=1,nlev
! Subprogram not used              do j=1,np
! Subprogram not used                 do i=1,np
! Subprogram not used                    v1     = elem(ie)%state%v(i,j,1,k,n0)
! Subprogram not used                    v2     = elem(ie)%state%v(i,j,2,k,n0)
! Subprogram not used 
! Subprogram not used                    ! Contravariant velocities
! Subprogram not used 
! Subprogram not used                    vco(i,j,1) = elem(ie)%Dinv(1,1,i,j)*v1 + elem(ie)%Dinv(1,2,i,j)*v2
! Subprogram not used                    vco(i,j,2) = elem(ie)%Dinv(2,1,i,j)*v1 + elem(ie)%Dinv(2,2,i,j)*v2
! Subprogram not used 
! Subprogram not used                    vgrad_ps(i,j,k) = ps(i,j)*(vco(i,j,1)*elem(ie)%derived%grad_lnps(i,j,1) + &
! Subprogram not used                         vco(i,j,2)*elem(ie)%derived%grad_lnps(i,j,2))
! Subprogram not used 
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used           call preq_omegap(elem(ie)%derived%div(:,:,:,n0),vgrad_ps,pdel,rpmid, &
! Subprogram not used                hvcoord%hybm,hvcoord%hybd,elem(ie)%derived%omega_p)
! Subprogram not used 
! Subprogram not used           Pscript(:,:,ie)        = 0.0_real_kind
! Subprogram not used           eta_dot_dp_deta(:,:,1) = 0.0_real_kind
! Subprogram not used 
! Subprogram not used           do k=1,nlev
! Subprogram not used              do j=1,np
! Subprogram not used                 do i=1,np
! Subprogram not used                    eta_dot_dp_deta(i,j,k+1) = eta_dot_dp_deta(i,j,k) + &
! Subprogram not used                         vgrad_ps(i,j,k)*hvcoord%hybd(k) + elem(ie)%derived%div(i,j,k,n0)*pdel(i,j,k)
! Subprogram not used                    ddiv = elem(ie)%derived%div(i,j,k,n0) - 0.5_real_kind*elem(ie)%derived%div(i,j,k,nm1)
! Subprogram not used                    Pscript(i,j,ie) = Pscript(i,j,ie) + ddiv*pdelref(k)
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used           do j=1,np
! Subprogram not used              do i=1,np
! Subprogram not used                 Pscript(i,j,ie) = elem(ie)%state%lnps(i,j,nm1) + &
! Subprogram not used                      dt2*( rpsref*Pscript(i,j,ie) - rps(i,j)*eta_dot_dp_deta(i,j,nlev+1) )
! Subprogram not used              end do
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used           do k=1,nlev-1
! Subprogram not used              do j=1,np
! Subprogram not used                 do i=1,np
! Subprogram not used                    eta_dot_dp_deta(i,j,k+1) = hvcoord%hybi(k+1)*eta_dot_dp_deta(i,j,nlev+1) - &
! Subprogram not used                         eta_dot_dp_deta(i,j,k+1)
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used           eta_dot_dp_deta(:,:,nlev+1) = 0.0_real_kind
! Subprogram not used 
! Subprogram not used           call preq_vertadv(elem(ie)%state%T(:,:,:,n0),elem(ie)%state%v(:,:,:,:,n0), &
! Subprogram not used                eta_dot_dp_deta,rpdel,T_vadv,v_vadv)
! Subprogram not used 
! Subprogram not used           suml(:,:) = 0.0_real_kind
! Subprogram not used 
! Subprogram not used           do k=1,nlev
! Subprogram not used 
! Subprogram not used              gradT(:,:,:) = gradient(elem(ie)%state%T(:,:,k,n0),deriv)*rrearth
! Subprogram not used              Crkk       = 0.5_real_kind
! Subprogram not used 
! Subprogram not used              do j=1,np
! Subprogram not used                 do i=1,np
! Subprogram not used                    term = Crkk*(elem(ie)%derived%div(i,j,k,n0) - &
! Subprogram not used                         0.5_real_kind*elem(ie)%derived%div(i,j,k,nm1))*pdelref(k)
! Subprogram not used                    suml(i,j)  = suml(i,j) + term
! Subprogram not used 
! Subprogram not used                    v1     = elem(ie)%state%v(i,j,1,k,n0)
! Subprogram not used                    v2     = elem(ie)%state%v(i,j,2,k,n0)
! Subprogram not used 
! Subprogram not used                    ! Contravariant velocities
! Subprogram not used 
! Subprogram not used                    vco(i,j,1) = elem(ie)%Dinv(1,1,i,j)*v1 + elem(ie)%Dinv(1,2,i,j)*v2
! Subprogram not used                    vco(i,j,2) = elem(ie)%Dinv(2,1,i,j)*v1 + elem(ie)%Dinv(2,2,i,j)*v2
! Subprogram not used 
! Subprogram not used                    vgradT = vco(i,j,1)*gradT(i,j,1) + vco(i,j,2)*gradT(i,j,2)
! Subprogram not used 
! Subprogram not used                    Tscript(i,j,k,ie) = elem(ie)%state%T(i,j,k,nm1) &
! Subprogram not used                         + dt2*(- vgradT - T_vadv(i,j,k)           &
! Subprogram not used                         + kappa*(elem(ie)%state%T(i,j,k,n0)*elem(ie)%derived%omega_p(i,j,k) &
! Subprogram not used                         + Tref(k)*rpmidref(k)*suml(i,j)))
! Subprogram not used                    suml(i,j)  = suml(i,j) + term
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used           HrefT(:,:)   = 0.0_real_kind
! Subprogram not used           HrefTm1(:,:) = 0.0_real_kind
! Subprogram not used           HT(:,:)      = 0.0_real_kind
! Subprogram not used 
! Subprogram not used           do k=nlev,1,-1
! Subprogram not used 
! Subprogram not used              do j=1,np
! Subprogram not used                 do i=1,np
! Subprogram not used                    hkl = rpmidref(k)*pdelref(k)
! Subprogram not used                    hkk = hkl*0.5_real_kind
! Subprogram not used                    Gref0(i,j)   = HrefT(i,j) + Rgas*hkk*elem(ie)%state%T(i,j,k,n0)
! Subprogram not used                    HrefT(i,j)   = HrefT(i,j) + Rgas*hkl*elem(ie)%state%T(i,j,k,n0)
! Subprogram not used                    Grefm1(i,j)  = HrefTm1(i,j) + Rgas*hkk*elem(ie)%state%T(i,j,k,nm1)
! Subprogram not used                    HrefTm1(i,j) = HrefTm1(i,j) + Rgas*hkl*elem(ie)%state%T(i,j,k,nm1)
! Subprogram not used                    hkl = rpmid(i,j,k)*pdel(i,j,k)
! Subprogram not used                    hkk = hkl*0.5_real_kind
! Subprogram not used                    Phi(i,j) = HT(i,j) + Rgas*hkk*elem(ie)%state%T(i,j,k,n0)
! Subprogram not used                    HT(i,j)  = HT(i,j) + Rgas*hkl*elem(ie)%state%T(i,j,k,n0)
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used              do j=1,np
! Subprogram not used                 do i=1,np
! Subprogram not used                    v1     = elem(ie)%state%v(i,j,1,k,n0)
! Subprogram not used                    v2     = elem(ie)%state%v(i,j,2,k,n0)
! Subprogram not used 
! Subprogram not used                    ! covariant velocity
! Subprogram not used 
! Subprogram not used                    vco(i,j,1) = elem(ie)%D(1,1,i,j)*v1 + elem(ie)%D(2,1,i,j)*v2
! Subprogram not used                    vco(i,j,2) = elem(ie)%D(1,2,i,j)*v1 + elem(ie)%D(2,2,i,j)*v2
! Subprogram not used 
! Subprogram not used                    E(i,j) = 0.5_real_kind*( v1*v1 + v2*v2 )
! Subprogram not used 
! Subprogram not used                    Gref0(i,j)  =  Gref0(i,j)  + elem(ie)%state%phis(i,j) + RTref(k)*elem(ie)%state%lnps(i,j,n0)
! Subprogram not used                    Grefm1(i,j) =  Grefm1(i,j) + elem(ie)%state%phis(i,j) + RTref(k)*elem(ie)%state%lnps(i,j,nm1)
! Subprogram not used 
! Subprogram not used                    Phi(i,j)    =  Phi(i,j) + E(i,j) + elem(ie)%state%phis(i,j)
! Subprogram not used                    dGref(i,j)  =  -(Gref0(i,j)  - 0.5_real_kind*Grefm1(i,j))
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used              elem(ie)%derived%zeta(:,:,k) = vorticity(vco,deriv)*rrearth
! Subprogram not used              grad_Phi(:,:,:)     = gradient(Phi,deriv)*rrearth
! Subprogram not used              grad_dGref(:,:,:,k) = gradient_wk(dGref,deriv)*rrearth
! Subprogram not used 
! Subprogram not used              do j=1,np
! Subprogram not used                 do i=1,np
! Subprogram not used 
! Subprogram not used                    elem(ie)%derived%zeta(i,j,k) = elem(ie)%rmetdet(i,j)*elem(ie)%derived%zeta(i,j,k)
! Subprogram not used                    hybfac =  hvcoord%hybm(k)*(ps(i,j)*rpmid(i,j,k))
! Subprogram not used 
! Subprogram not used                    glnps1 = elem(ie)%Dinv(1,1,i,j)*elem(ie)%derived%grad_lnps(i,j,1) + &
! Subprogram not used                         elem(ie)%Dinv(2,1,i,j)*elem(ie)%derived%grad_lnps(i,j,2)
! Subprogram not used                    glnps2 = elem(ie)%Dinv(1,2,i,j)*elem(ie)%derived%grad_lnps(i,j,1) + &
! Subprogram not used                         elem(ie)%Dinv(2,2,i,j)*elem(ie)%derived%grad_lnps(i,j,2)
! Subprogram not used 
! Subprogram not used                    v1 = elem(ie)%Dinv(1,1,i,j)*grad_Phi(i,j,1) + elem(ie)%Dinv(2,1,i,j)*grad_Phi(i,j,2)
! Subprogram not used                    v2 = elem(ie)%Dinv(1,2,i,j)*grad_Phi(i,j,1) + elem(ie)%Dinv(2,2,i,j)*grad_Phi(i,j,2)
! Subprogram not used 
! Subprogram not used                    Vscript(i,j,1,k,ie) = - v_vadv(i,j,1,k) &
! Subprogram not used                         + elem(ie)%state%v(i,j,2,k,n0) * (elem(ie)%fcor(i,j) + elem(ie)%derived%zeta(i,j,k)) &
! Subprogram not used                         - v1 - Rgas*hybfac*elem(ie)%state%T(i,j,k,n0)*glnps1
! Subprogram not used 
! Subprogram not used                    Vscript(i,j,2,k,ie) = - v_vadv(i,j,2,k) &
! Subprogram not used                         - elem(ie)%state%v(i,j,1,k,n0) * (elem(ie)%fcor(i,j) + elem(ie)%derived%zeta(i,j,k)) &
! Subprogram not used                         - v2 - Rgas*hybfac*elem(ie)%state%T(i,j,k,n0)*glnps2
! Subprogram not used 
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used           do k=1,nlev
! Subprogram not used 
! Subprogram not used              do j=1,np
! Subprogram not used                 do i=1,np
! Subprogram not used                    Vs1 = elem(ie)%Dinv(1,1,i,j)*grad_dGref(i,j,1,k) + elem(ie)%Dinv(2,1,i,j)*grad_dGref(i,j,2,k)
! Subprogram not used                    Vs2 = elem(ie)%Dinv(1,2,i,j)*grad_dGref(i,j,1,k) + elem(ie)%Dinv(2,2,i,j)*grad_dGref(i,j,2,k)
! Subprogram not used 
! Subprogram not used                    Vscript(i,j,1,k,ie) = elem(ie)%mp(i,j)*Vscript(i,j,1,k,ie) + Vs1
! Subprogram not used                    Vscript(i,j,2,k,ie) = elem(ie)%mp(i,j)*Vscript(i,j,2,k,ie) + Vs2
! Subprogram not used 
! Subprogram not used                    Vscript(i,j,1,k,ie) = elem(ie)%mp(i,j)*elem(ie)%state%v(i,j,1,k,nm1) + dt2*Vscript(i,j,1,k,ie)
! Subprogram not used                    Vscript(i,j,2,k,ie) = elem(ie)%mp(i,j)*elem(ie)%state%v(i,j,2,k,nm1) + dt2*Vscript(i,j,2,k,ie)
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used           HrefTscript(:,:) = 0.0_real_kind
! Subprogram not used 
! Subprogram not used           do k=nlev,1,-1
! Subprogram not used 
! Subprogram not used              do j=1,np
! Subprogram not used                 do i=1,np
! Subprogram not used                    hkl = rpmidref(k)*pdelref(k)
! Subprogram not used                    hkk = hkl*0.5_real_kind
! Subprogram not used                    B(i,j,k,ie)      = HrefTscript(i,j) + Rgas*hkk*Tscript(i,j,k,ie)
! Subprogram not used                    B(i,j,k,ie)      = B(i,j,k,ie) +  elem(ie)%state%phis(i,j) + RTref(k)*Pscript(i,j,ie)
! Subprogram not used                    HrefTscript(i,j) = HrefTscript(i,j) + Rgas*hkl*Tscript(i,j,k,ie)
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used           kptr=0
! Subprogram not used           call edgeVpack(edge2, Vscript(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)
! Subprogram not used 
! Subprogram not used        end do
! Subprogram not used 
! Subprogram not used        call bndry_exchangeV(cg%hybrid,edge2)
! Subprogram not used 
! Subprogram not used        do ie = nets, nete
! Subprogram not used 
! Subprogram not used           kptr=0
! Subprogram not used           call edgeVunpack(edge2, Vscript(1,1,1,1,ie), 2*nlev, kptr, elem(ie)%desc)
! Subprogram not used           do k=1,nlev
! Subprogram not used 
! Subprogram not used              do j=1,np
! Subprogram not used                 do i=1,np
! Subprogram not used                    Vscript(i,j,1,k,ie) = elem(ie)%rmp(i,j)*Vscript(i,j,1,k,ie)
! Subprogram not used                    Vscript(i,j,2,k,ie) = elem(ie)%rmp(i,j)*Vscript(i,j,2,k,ie)
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used              do j=1,np
! Subprogram not used                 do i=1,np
! Subprogram not used 
! Subprogram not used                    ! Contravariant Vscript
! Subprogram not used 
! Subprogram not used                    gVscript(i,j,1) = elem(ie)%Dinv(1,1,i,j)*Vscript(i,j,1,k,ie) + &
! Subprogram not used                         elem(ie)%Dinv(1,2,i,j)*Vscript(i,j,2,k,ie)
! Subprogram not used                    gVscript(i,j,2) = elem(ie)%Dinv(2,1,i,j)*Vscript(i,j,1,k,ie) + &
! Subprogram not used                         elem(ie)%Dinv(2,2,i,j)*Vscript(i,j,2,k,ie) 
! Subprogram not used 
! Subprogram not used                    gVscript(i,j,1) = elem(ie)%metdet(i,j)*gVscript(i,j,1)
! Subprogram not used                    gVscript(i,j,2) = elem(ie)%metdet(i,j)*gVscript(i,j,2)
! Subprogram not used 
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used              div_Vscript(:,:,k) = divergence(gVscript,deriv)*rrearth
! Subprogram not used 
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used           do k=1,nlev
! Subprogram not used 
! Subprogram not used              do j=1,np
! Subprogram not used                 do i=1,np
! Subprogram not used                    C(i,j,k,ie) = elem(ie)%metdet(i,j)*B(i,j,k,ie)
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used              do l=1,nlev
! Subprogram not used                 do j=1,np
! Subprogram not used                    do i=1,np
! Subprogram not used                       C(i,j,k,ie) = C(i,j,k,ie) - dt*Amat(l,k)*div_Vscript(i,j,l)
! Subprogram not used                    end do
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used           ! ===============================================================
! Subprogram not used           !  Weight C (the RHS of the helmholtz problem) by the mass matrix
! Subprogram not used           ! ===============================================================
! Subprogram not used 
! Subprogram not used           do k=1,nlev
! Subprogram not used              do j=1,np
! Subprogram not used                 do i=1,np
! Subprogram not used                    C(i,j,k,ie) = elem(ie)%mp(i,j)*C(i,j,k,ie)
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used           ! ===================================
! Subprogram not used           ! Pack C into the edge1 buffer
! Subprogram not used           ! ===================================
! Subprogram not used 
! Subprogram not used           kptr=0
! Subprogram not used           call edgeVpack(edge1,C(1,1,1,ie),nlev,kptr,elem(ie)%desc)
! Subprogram not used 
! Subprogram not used        end do
! Subprogram not used 
! Subprogram not used        ! ==================================
! Subprogram not used        ! boundary exchange C
! Subprogram not used        ! ==================================
! Subprogram not used 
! Subprogram not used        call bndry_exchangeV(cg%hybrid,edge1)
! Subprogram not used 
! Subprogram not used        do ie=nets,nete
! Subprogram not used 
! Subprogram not used           ! ===================================
! Subprogram not used           ! Unpack C from the edge1 buffer
! Subprogram not used           ! ===================================
! Subprogram not used 
! Subprogram not used           kptr=0
! Subprogram not used           call edgeVunpack(edge1, C(1,1,1,ie), nlev, kptr, elem(ie)%desc)
! Subprogram not used 
! Subprogram not used           ! ===============================================
! Subprogram not used           ! Complete global assembly by normalizing by rmp
! Subprogram not used           ! ===============================================
! Subprogram not used 
! Subprogram not used           do k=1,nlev
! Subprogram not used 
! Subprogram not used              do j=1,np
! Subprogram not used                 do i=1,np
! Subprogram not used                    D(i,j,k,ie) = elem(ie)%rmp(i,j)*C(i,j,k,ie)
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used           do k=1,nlev
! Subprogram not used 
! Subprogram not used              do j=1,np
! Subprogram not used                 do i=1,np
! Subprogram not used                    C(i,j,k,ie) = 0.0_real_kind
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used              do l=1,nlev
! Subprogram not used                 do j=1,np
! Subprogram not used                    do i=1,np
! Subprogram not used                       C(i,j,k,ie) = C(i,j,k,ie) + Emat_inv(l,k)*D(i,j,l,ie)
! Subprogram not used                    end do
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used        end do
! Subprogram not used        ! ==========================================
! Subprogram not used        ! solve for Gamma_ref, given C as RHS input
! Subprogram not used        ! ==========================================
! Subprogram not used 
! Subprogram not used        Gamma_ref = pcg_solver(elem, &
! Subprogram not used             C,          &
! Subprogram not used             cg,         &
! Subprogram not used             red,        &
! Subprogram not used             edge1,      &
! Subprogram not used             edge2,      &   
! Subprogram not used             Lambda,     &   
! Subprogram not used             deriv,      &   
! Subprogram not used             nets,       & 
! Subprogram not used             nete,       &
! Subprogram not used             blkjac)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used        ! ================================================================
! Subprogram not used        ! Backsubstitute Gamma_ref into semi-implicit system of equations
! Subprogram not used        ! to find prognostic variables at time level n+1
! Subprogram not used        ! ================================================================
! Subprogram not used 
! Subprogram not used        do ie = nets, nete
! Subprogram not used 
! Subprogram not used           do k=1,nlev
! Subprogram not used 
! Subprogram not used              do j=1,np
! Subprogram not used                 do i=1,np
! Subprogram not used                    Gref(i,j,k,ie) = 0.0_real_kind
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used              do l=1,nlev
! Subprogram not used                 do j=1,np
! Subprogram not used                    do i=1,np
! Subprogram not used                       Gref(i,j,k,ie) = Gref(i,j,k,ie) + Emat(l,k)*Gamma_ref(i,j,l,ie)
! Subprogram not used                    end do
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used              do j=1,np  
! Subprogram not used                 do i=1,np
! Subprogram not used                    B(i,j,k,ie) = elem(ie)%mp(i,j) * dt * (B(i,j,k,ie) - Gref(i,j,k,ie)) 
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used           kptr=0
! Subprogram not used           call edgeVpack(edge1,B(:,:,:,ie),nlev,kptr,elem(ie)%desc)
! Subprogram not used 
! Subprogram not used        end do
! Subprogram not used 
! Subprogram not used        call bndry_exchangeV(cg%hybrid,edge1)
! Subprogram not used 
! Subprogram not used        do ie = nets, nete
! Subprogram not used 
! Subprogram not used           kptr=0
! Subprogram not used           call edgeVunpack(edge1, B(:,:,:,ie), nlev, kptr, elem(ie)%desc)
! Subprogram not used           do k=1,nlev
! Subprogram not used 
! Subprogram not used              do j=1,np
! Subprogram not used                 do i=1,np
! Subprogram not used                    B(i,j,k,ie) = elem(ie)%rmp(i,j)*B(i,j,k,ie)
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used           do k=1,nlev
! Subprogram not used 
! Subprogram not used              do j=1,np
! Subprogram not used                 do i=1,np
! Subprogram not used                    D(i,j,k,ie) = 0.0_real_kind
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used              do l=1,nlev
! Subprogram not used                 do j=1,np
! Subprogram not used                    do i=1,np
! Subprogram not used                       D(i,j,k,ie) = D(i,j,k,ie) + Emat_inv(l,k)*B(i,j,l,ie)
! Subprogram not used                    end do
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used           do k=1,nlev 
! Subprogram not used 
! Subprogram not used              do j=1,np
! Subprogram not used                 do i=1,np
! Subprogram not used                    elem(ie)%derived%div(i,j,k,np1) = 0.0_real_kind
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used              do l=1,nlev
! Subprogram not used                 do j=1,np
! Subprogram not used                    do i=1,np
! Subprogram not used                       elem(ie)%derived%div(i,j,k,np1) = elem(ie)%derived%div(i,j,k,np1) + Emat(l,k)*D(i,j,l,ie)/Lambda(l)
! Subprogram not used                    end do
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used           do k=1,nlev
! Subprogram not used 
! Subprogram not used              grad_Gref(:,:,:,k)=gradient_wk(Gref(:,:,k,ie),deriv)*rrearth
! Subprogram not used 
! Subprogram not used              do j=1,np
! Subprogram not used                 do i=1,np
! Subprogram not used                    gGr1 = grad_Gref(i,j,1,k)
! Subprogram not used                    gGr2 = grad_Gref(i,j,2,k)
! Subprogram not used                    elem(ie)%state%v(i,j,1,k,np1) = elem(ie)%Dinv(1,1,i,j)*gGr1 + elem(ie)%Dinv(2,1,i,j)*gGr2
! Subprogram not used                    elem(ie)%state%v(i,j,2,k,np1) = elem(ie)%Dinv(1,2,i,j)*gGr1 + elem(ie)%Dinv(2,2,i,j)*gGr2
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used              do j=1,np
! Subprogram not used                 do i=1,np
! Subprogram not used                    Pscript(i,j,ie) = Pscript(i,j,ie) - dt*Pvec(k)*elem(ie)%derived%div(i,j,k,np1)
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used 
! Subprogram not used              do l=1,nlev
! Subprogram not used                 do j=1,np
! Subprogram not used                    do i=1,np
! Subprogram not used                       Tscript(i,j,k,ie) = Tscript(i,j,k,ie) - dt*Tmat(l,k)*elem(ie)%derived%div(i,j,l,np1)
! Subprogram not used                    end do
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used           do j=1,np
! Subprogram not used              do i=1,np
! Subprogram not used                 Pscript(i,j,ie) = elem(ie)%mp(i,j)*Pscript(i,j,ie)
! Subprogram not used              end do
! Subprogram not used           end do
! Subprogram not used           do k=1,nlev
! Subprogram not used              do j=1,np
! Subprogram not used                 do i=1,np
! Subprogram not used                    Tscript(i,j,k,ie) = elem(ie)%mp(i,j)*Tscript(i,j,k,ie)
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used           ! ===============================================
! Subprogram not used           ! Pack v at time level n+1 into the edge3p1 buffer
! Subprogram not used           ! ===============================================
! Subprogram not used 
! Subprogram not used           kptr=0
! Subprogram not used           call edgeVpack(edge3p1, elem(ie)%state%v(:,:,:,:,np1),2*nlev,kptr,elem(ie)%desc)
! Subprogram not used 
! Subprogram not used           kptr=2*nlev
! Subprogram not used           call edgeVpack(edge3p1, Tscript(:,:,:,ie),nlev,kptr,elem(ie)%desc)
! Subprogram not used 
! Subprogram not used           kptr=3*nlev
! Subprogram not used           call edgeVpack(edge3p1, Pscript(:,:,ie),1,kptr,elem(ie)%desc)
! Subprogram not used 
! Subprogram not used        end do
! Subprogram not used 
! Subprogram not used        ! ======================================
! Subprogram not used        ! boundary exchange v at time level n+1
! Subprogram not used        ! ======================================
! Subprogram not used 
! Subprogram not used        call bndry_exchangeV(cg%hybrid,edge3p1)
! Subprogram not used 
! Subprogram not used        do ie=nets,nete
! Subprogram not used 
! Subprogram not used           ! ===================================
! Subprogram not used           ! Unpack v from the edge2 buffer
! Subprogram not used           ! ===================================
! Subprogram not used 
! Subprogram not used           kptr=0
! Subprogram not used           call edgeVunpack(edge3p1, elem(ie)%state%v(:,:,:,:,np1), 2*nlev, kptr, elem(ie)%desc)
! Subprogram not used 
! Subprogram not used           kptr=2*nlev
! Subprogram not used           call edgeVunpack(edge3p1, Tscript(:,:,:,ie), nlev, kptr, elem(ie)%desc)
! Subprogram not used 
! Subprogram not used           kptr=3*nlev
! Subprogram not used           call edgeVunpack(edge3p1, Pscript(:,:,ie), 1, kptr, elem(ie)%desc)
! Subprogram not used 
! Subprogram not used           ! ==========================================================
! Subprogram not used           ! Complete global assembly by normalizing velocity by rmp
! Subprogram not used           ! Vscript = Vscript - dt*grad(Gref)
! Subprogram not used           ! ==========================================================
! Subprogram not used 
! Subprogram not used           do k=1,nlev
! Subprogram not used 
! Subprogram not used              do j=1,np
! Subprogram not used                 do i=1,np
! Subprogram not used                    elem(ie)%state%v(i,j,1,k,np1) = Vscript(i,j,1,k,ie) + dt*elem(ie)%rmp(i,j)*elem(ie)%state%v(i,j,1,k,np1)
! Subprogram not used                    elem(ie)%state%v(i,j,2,k,np1) = Vscript(i,j,2,k,ie) + dt*elem(ie)%rmp(i,j)*elem(ie)%state%v(i,j,2,k,np1)
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used           do j=1,np
! Subprogram not used              do i=1,np
! Subprogram not used                 elem(ie)%state%lnps(i,j,np1) = elem(ie)%rmp(i,j)*Pscript(i,j,ie)
! Subprogram not used              end do
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used           do k=1,nlev
! Subprogram not used              do j=1,np
! Subprogram not used                 do i=1,np
! Subprogram not used                    elem(ie)%state%T(i,j,k,np1) = elem(ie)%rmp(i,j)*Tscript(i,j,k,ie)
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used        end do
! Subprogram not used        call prim_diffusion(elem, nets,nete,np1,deriv,dt2,cg%hybrid)
! Subprogram not used        call t_stopf('prim_advance_si')
! Subprogram not used   end subroutine prim_advance_si


! Subprogram not used   subroutine preq_robert3(nm1,n0,np1,elem,hvcoord,nets,nete)
! Subprogram not used   use dimensions_mod, only : np, nlev, qsize
! Subprogram not used   use hybvcoord_mod, only : hvcoord_t
! Subprogram not used   use element_mod, only : element_t
! Subprogram not used   use time_mod, only: smooth
! Subprogram not used   use control_mod, only : integration
! Subprogram not used 
! Subprogram not used   implicit none
! Subprogram not used   integer              , intent(in) :: nm1,n0,np1,nets,nete
! Subprogram not used   type (hvcoord_t), intent(in)      :: hvcoord
! Subprogram not used   type (element_t)     , intent(inout) :: elem(:)
! Subprogram not used   
! Subprogram not used   
! Subprogram not used   integer :: i,j,k,ie,q
! Subprogram not used   real (kind=real_kind) :: dp
! Subprogram not used   logical :: filter_ps = .false.
! Subprogram not used   if (integration == "explicit") filter_ps = .true.
! Subprogram not used 
! Subprogram not used   call t_startf('preq_robert')
! Subprogram not used   do ie=nets,nete
! Subprogram not used      if (filter_ps) then
! Subprogram not used         elem(ie)%state%ps_v(:,:,n0) = elem(ie)%state%ps_v(:,:,n0) + smooth*(elem(ie)%state%ps_v(:,:,nm1) &
! Subprogram not used              - 2.0D0*elem(ie)%state%ps_v(:,:,n0)   + elem(ie)%state%ps_v(:,:,np1))
! Subprogram not used         elem(ie)%state%lnps(:,:,n0) = LOG(elem(ie)%state%ps_v(:,:,n0))
! Subprogram not used      else
! Subprogram not used         elem(ie)%state%lnps(:,:,n0) = elem(ie)%state%lnps(:,:,n0) + smooth*(elem(ie)%state%lnps(:,:,nm1) &
! Subprogram not used              - 2.0D0*elem(ie)%state%lnps(:,:,n0)   + elem(ie)%state%lnps(:,:,np1))
! Subprogram not used         elem(ie)%state%ps_v(:,:,n0) = EXP(elem(ie)%state%lnps(:,:,n0))
! Subprogram not used      endif
! Subprogram not used      
! Subprogram not used      elem(ie)%state%T(:,:,:,n0) = elem(ie)%state%T(:,:,:,n0) + smooth*(elem(ie)%state%T(:,:,:,nm1) &
! Subprogram not used           - 2.0D0*elem(ie)%state%T(:,:,:,n0)   + elem(ie)%state%T(:,:,:,np1))
! Subprogram not used      elem(ie)%state%v(:,:,:,:,n0) = elem(ie)%state%v(:,:,:,:,n0) + smooth*(elem(ie)%state%v(:,:,:,:,nm1) &
! Subprogram not used           - 2.0D0*elem(ie)%state%v(:,:,:,:,n0) + elem(ie)%state%v(:,:,:,:,np1))
! Subprogram not used      
! Subprogram not used   end do
! Subprogram not used   call t_stopf('preq_robert')
! Subprogram not used   
! Subprogram not used   end subroutine preq_robert3
  
  


  subroutine applyCAMforcing(elem,hvcoord,np1,np1_qdp,dt_q,nets,nete)
  use dimensions_mod, only : np, nlev, qsize
  use element_mod, only : element_t
  use hybvcoord_mod, only : hvcoord_t
  use control_mod, only : moisture
  use physical_constants, only: Cp  
  implicit none
  type (element_t)     , intent(inout) :: elem(:)
  real (kind=real_kind), intent(in) :: dt_q
  type (hvcoord_t), intent(in)      :: hvcoord
  integer,  intent(in) :: np1,nets,nete,np1_qdp
  
  ! local
  integer :: i,j,k,ie,q
  real (kind=real_kind) :: v1,dp
  real (kind=real_kind) :: beta(np,np),E0(np,np),ED(np,np),dp0m1(np,np),dpsum(np,np)
  logical :: wet  
  
  wet = (moisture /= "dry")  
  
  do ie=nets,nete
     ! apply forcing to Qdp
     elem(ie)%derived%FQps(:,:,1)=0
     do q=1,qsize
        do k=1,nlev
           do j=1,np	
              do i=1,np
                 v1 = dt_q*elem(ie)%derived%FQ(i,j,k,q,1)
                 !if (elem(ie)%state%Qdp(i,j,k,q,np1) + v1 < 0 .and. v1<0) then
                 if (elem(ie)%state%Qdp(i,j,k,q,np1_qdp) + v1 < 0 .and. v1<0) then
                    !if (elem(ie)%state%Qdp(i,j,k,q,np1) < 0 ) then
                    if (elem(ie)%state%Qdp(i,j,k,q,np1_qdp) < 0 ) then
                       v1=0  ! Q already negative, dont make it more so
                    else
                       !v1 = -elem(ie)%state%Qdp(i,j,k,q,np1)
                       v1 = -elem(ie)%state%Qdp(i,j,k,q,np1_qdp)
                    endif
                 endif
                 !elem(ie)%state%Qdp(i,j,k,q,np1) = elem(ie)%state%Qdp(i,j,k,q,np1)+v1
                 elem(ie)%state%Qdp(i,j,k,q,np1_qdp) = elem(ie)%state%Qdp(i,j,k,q,np1_qdp)+v1
                 if (q==1) then
                    elem(ie)%derived%FQps(i,j,1)=elem(ie)%derived%FQps(i,j,1)+v1/dt_q
                 endif
              enddo
           enddo
        enddo
     enddo
     
     if (wet .and. qsize>0) then
        ! to conserve dry mass in the precese of Q1 forcing:
        elem(ie)%state%ps_v(:,:,np1) = elem(ie)%state%ps_v(:,:,np1) + &
             dt_q*elem(ie)%derived%FQps(:,:,1)
     endif





     ! Qdp(np1) and ps_v(np1) were updated by forcing - update Q(np1)
     do q=1,qsize
        do k=1,nlev
           do j=1,np	
              do i=1,np
                 dp = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                      ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,np1)
                 elem(ie)%state%Q(i,j,k,q) = elem(ie)%state%Qdp(i,j,k,q,np1_qdp)/dp
              enddo
           enddo
        enddo
     enddo
     
     elem(ie)%state%T(:,:,:,np1) = elem(ie)%state%T(:,:,:,np1) + &
          dt_q*elem(ie)%derived%FT(:,:,:,1)
     elem(ie)%state%v(:,:,:,:,np1) = elem(ie)%state%v(:,:,:,:,np1) + &
          dt_q*elem(ie)%derived%FM(:,:,:,:,1)
  enddo
  end subroutine applyCAMforcing



! Subprogram not used   subroutine applyCAMforcing_dynamics(elem,hvcoord,np1,dt_q,nets,nete)
! Subprogram not used   use dimensions_mod, only : np, nlev, qsize
! Subprogram not used   use element_mod, only : element_t
! Subprogram not used   use hybvcoord_mod, only : hvcoord_t
! Subprogram not used   
! Subprogram not used   implicit none
! Subprogram not used   type (element_t)     , intent(inout) :: elem(:)
! Subprogram not used   real (kind=real_kind), intent(in) :: dt_q
! Subprogram not used   type (hvcoord_t), intent(in)      :: hvcoord
! Subprogram not used   integer,  intent(in) :: np1,nets,nete
! Subprogram not used   
! Subprogram not used   ! local
! Subprogram not used   integer :: i,j,k,ie,q
! Subprogram not used   real (kind=real_kind) :: v1,dp
! Subprogram not used   logical :: wet  
! Subprogram not used   
! Subprogram not used   do ie=nets,nete
! Subprogram not used      elem(ie)%state%T(:,:,:,np1) = elem(ie)%state%T(:,:,:,np1) + &
! Subprogram not used           dt_q*elem(ie)%derived%FT(:,:,:,1)
! Subprogram not used      elem(ie)%state%v(:,:,:,:,np1) = elem(ie)%state%v(:,:,:,:,np1) + &
! Subprogram not used           dt_q*elem(ie)%derived%FM(:,:,:,:,1)
! Subprogram not used   enddo
! Subprogram not used   end subroutine applyCAMforcing_dynamics



! Subprogram not used   subroutine advance_hypervis(edge3,elem,hvcoord,hybrid,deriv,nt,nets,nete,dt2,eta_ave_w)
! Subprogram not used   !
! Subprogram not used   !  take one timestep of:  
! Subprogram not used   !          u(:,:,:,np) = u(:,:,:,np) +  dt2*nu*laplacian**order ( u )
! Subprogram not used   !          T(:,:,:,np) = T(:,:,:,np) +  dt2*nu_s*laplacian**order ( T )
! Subprogram not used   !
! Subprogram not used   !
! Subprogram not used   !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
! Subprogram not used   !
! Subprogram not used   !
! Subprogram not used   use dimensions_mod, only : np, np, nlev
! Subprogram not used   use control_mod, only : nu, nu_div, nu_s, hypervis_order, hypervis_subcycle, nu_p, nu_top, psurf_vis
! Subprogram not used   use hybrid_mod, only : hybrid_t
! Subprogram not used   use hybvcoord_mod, only : hvcoord_t
! Subprogram not used   use element_mod, only : element_t
! Subprogram not used   use derivative_mod, only : derivative_t, laplace_sphere_wk, vlaplace_sphere_wk
! Subprogram not used   use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack
! Subprogram not used   use bndry_mod, only : bndry_exchangev
! Subprogram not used   use viscosity_mod, only : biharmonic_wk
! Subprogram not used   use physical_constants, only: Cp
! Subprogram not used !  use time_mod, only : TimeLevel_t
! Subprogram not used   implicit none
! Subprogram not used   
! Subprogram not used   type (hybrid_t)      , intent(in) :: hybrid
! Subprogram not used   type (element_t)     , intent(inout), target :: elem(:)
! Subprogram not used   type (EdgeBuffer_t)  , intent(inout) :: edge3
! Subprogram not used   type (derivative_t)  , intent(in) :: deriv
! Subprogram not used   type (hvcoord_t), intent(in)      :: hvcoord
! Subprogram not used !  type (TimeLevel_t)   , intent(in) :: tl
! Subprogram not used   
! Subprogram not used   real (kind=real_kind) :: dt2
! Subprogram not used   integer :: nets,nete
! Subprogram not used   
! Subprogram not used   ! local
! Subprogram not used   real (kind=real_kind) :: eta_ave_w  ! weighting for mean flux terms
! Subprogram not used   real (kind=real_kind) :: nu_scale, dpdn,dpdn0, nu_scale_top
! Subprogram not used   integer :: k,kptr,i,j,ie,ic,nt
! Subprogram not used   real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)      :: vtens   
! Subprogram not used   real (kind=real_kind), dimension(np,np,nlev,nets:nete)        :: ptens
! Subprogram not used   real (kind=real_kind), dimension(np,np,nets:nete) :: pstens	
! Subprogram not used   real (kind=real_kind), dimension(np,np,nlev) :: p
! Subprogram not used   real (kind=real_kind), dimension(np,np) :: dptemp1,dptemp2
! Subprogram not used   
! Subprogram not used   
! Subprogram not used ! NOTE: PGI compiler bug: when using spheremp, rspheremp and ps as pointers to elem(ie)% members,
! Subprogram not used   !       data is incorrect (offset by a few numbers actually)
! Subprogram not used   !       removed for now.  
! Subprogram not used   !       real (kind=real_kind), dimension(:,:), pointer :: spheremp,rspheremp
! Subprogram not used   !       real (kind=real_kind), dimension(:,:,:), pointer   :: ps
! Subprogram not used   
! Subprogram not used   real (kind=real_kind), dimension(np,np) :: lap_p
! Subprogram not used   real (kind=real_kind), dimension(np,np,2) :: lap_v
! Subprogram not used   real (kind=real_kind) :: v1,v2,dt,heating,utens_tmp,vtens_tmp,ptens_tmp
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   if (nu_s == 0 .and. nu == 0 .and. nu_p==0 ) return;
! Subprogram not used   call t_barrierf('sync_advance_hypervis', hybrid%par%comm)
! Subprogram not used   call t_startf('advance_hypervis')
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   dt=dt2/hypervis_subcycle
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used   !  regular viscosity  
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used   if (hypervis_order == 1) then
! Subprogram not used      if (nu_p>0) call abortmp( 'ERROR: hypervis_order == 1 not coded for nu_p>0')
! Subprogram not used      do ic=1,hypervis_subcycle
! Subprogram not used         do ie=nets,nete
! Subprogram not used            
! Subprogram not used            do k=1,nlev
! Subprogram not used               lap_p=laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
! Subprogram not used               lap_v=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)
! Subprogram not used               ! advace in time.  (note: DSS commutes with time stepping, so we
! Subprogram not used               ! can time advance and then DSS.  this has the advantage of
! Subprogram not used               ! not letting any discontinuties accumulate in p,v via roundoff
! Subprogram not used               do j=1,np
! Subprogram not used                  do i=1,np             
! Subprogram not used                     elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt)*elem(ie)%spheremp(i,j)  +  dt*nu_s*lap_p(i,j)
! Subprogram not used                     elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%state%v(i,j,1,k,nt)*elem(ie)%spheremp(i,j) + dt*nu*lap_v(i,j,1)
! Subprogram not used                     elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%state%v(i,j,2,k,nt)*elem(ie)%spheremp(i,j) + dt*nu*lap_v(i,j,2)
! Subprogram not used                  enddo
! Subprogram not used               enddo
! Subprogram not used            enddo
! Subprogram not used            
! Subprogram not used            kptr=0
! Subprogram not used            call edgeVpack(edge3, elem(ie)%state%T(:,:,:,nt),nlev,kptr,elem(ie)%desc)
! Subprogram not used            kptr=nlev
! Subprogram not used            call edgeVpack(edge3,elem(ie)%state%v(:,:,:,:,nt),2*nlev,kptr,elem(ie)%desc)
! Subprogram not used         enddo
! Subprogram not used         
! Subprogram not used         call bndry_exchangeV(hybrid,edge3)
! Subprogram not used         
! Subprogram not used         do ie=nets,nete
! Subprogram not used            
! Subprogram not used            kptr=0
! Subprogram not used            call edgeVunpack(edge3, elem(ie)%state%T(:,:,:,nt), nlev, kptr, elem(ie)%desc)
! Subprogram not used            kptr=nlev
! Subprogram not used            call edgeVunpack(edge3, elem(ie)%state%v(:,:,:,:,nt), 2*nlev, kptr, elem(ie)%desc)
! Subprogram not used            
! Subprogram not used            ! apply inverse mass matrix
! Subprogram not used            do k=1,nlev
! Subprogram not used               do j=1,np
! Subprogram not used                  do i=1,np             
! Subprogram not used                     elem(ie)%state%T(i,j,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%T(i,j,k,nt)
! Subprogram not used                     elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,1,k,nt)
! Subprogram not used                     elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,2,k,nt)
! Subprogram not used                  enddo
! Subprogram not used               enddo
! Subprogram not used            enddo
! Subprogram not used         enddo
! Subprogram not used      enddo  ! subcycle
! Subprogram not used   endif
! Subprogram not used   
! Subprogram not used   
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used   !  hyper viscosity  
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used ! nu_p=0:
! Subprogram not used !   scale T dissipaton by dp  (conserve IE, dissipate T^2)
! Subprogram not used ! nu_p>0
! Subprogram not used !   dont scale:  T equation IE dissipation matches (to truncation error) 
! Subprogram not used !                IE dissipation from continuity equation 
! Subprogram not used !                (1 deg: to about 0.1 W/m^2)
! Subprogram not used !
! Subprogram not used   if (hypervis_order == 2) then
! Subprogram not used      do ic=1,hypervis_subcycle
! Subprogram not used         call biharmonic_wk(elem,pstens,ptens,vtens,deriv,edge3,hybrid,nt,nets,nete)
! Subprogram not used         do ie=nets,nete
! Subprogram not used 
! Subprogram not used            ! comptue mean flux
! Subprogram not used            if (nu_p>0) then
! Subprogram not used               do k=1,nlev
! Subprogram not used                  dptemp1(:,:) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
! Subprogram not used                       ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,nt)
! Subprogram not used                  elem(ie)%derived%dpdiss_ave(:,:,k)=elem(ie)%derived%dpdiss_ave(:,:,k)+eta_ave_w*dptemp1(:,:)/hypervis_subcycle
! Subprogram not used                  
! Subprogram not used                  dptemp2(:,:) = (hvcoord%hybi(k+1)-hvcoord%hybi(k))*pstens(:,:,ie)
! Subprogram not used                  elem(ie)%derived%dpdiss_biharmonic(:,:,k)=&
! Subprogram not used                       elem(ie)%derived%dpdiss_biharmonic(:,:,k)+eta_ave_w*dptemp2(:,:)/hypervis_subcycle
! Subprogram not used               enddo
! Subprogram not used            endif
! Subprogram not used            nu_scale=1
! Subprogram not used            do k=1,nlev
! Subprogram not used               ! advace in time.  
! Subprogram not used               ! note: DSS commutes with time stepping, so we can time advance and then DSS.
! Subprogram not used               ! note: weak operators alreayd have mass matrix "included"
! Subprogram not used               
! Subprogram not used               ! add regular diffusion in top 3 layers:
! Subprogram not used               if (nu_top>0 .and. k<=3) then
! Subprogram not used                  lap_p=laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
! Subprogram not used                  lap_v=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)
! Subprogram not used               endif
! Subprogram not used               nu_scale_top = 1
! Subprogram not used               if (k==1) nu_scale_top=4
! Subprogram not used               if (k==2) nu_scale_top=2
! Subprogram not used               
! Subprogram not used               do j=1,np
! Subprogram not used                  do i=1,np
! Subprogram not used                     if (nu_p==0) then
! Subprogram not used                        ! normalize so as to conserve IE  
! Subprogram not used                        ! scale by 1/rho (normalized to be O(1))
! Subprogram not used                        ! dp/dn = O(ps0)*O(delta_eta) = O(ps0)/O(nlev)
! Subprogram not used                        dpdn = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
! Subprogram not used                             ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,nt)  
! Subprogram not used                        dpdn0 = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
! Subprogram not used                             ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*hvcoord%ps0
! Subprogram not used                        nu_scale = dpdn0/dpdn
! Subprogram not used                     endif
! Subprogram not used 
! Subprogram not used                     ! biharmonic terms need a negative sign:
! Subprogram not used                     if (nu_top>0 .and. k<=3) then
! Subprogram not used                        utens_tmp=(-nu*vtens(i,j,1,k,ie) + nu_scale_top*nu_top*lap_v(i,j,1))
! Subprogram not used                        vtens_tmp=(-nu*vtens(i,j,2,k,ie) + nu_scale_top*nu_top*lap_v(i,j,2))
! Subprogram not used                        ptens_tmp=nu_scale*(-nu_s*ptens(i,j,k,ie) + nu_scale_top*nu_top*lap_p(i,j) )
! Subprogram not used                     else
! Subprogram not used                        utens_tmp=-nu*vtens(i,j,1,k,ie)
! Subprogram not used                        vtens_tmp=-nu*vtens(i,j,2,k,ie)
! Subprogram not used                        ptens_tmp=-nu_scale*nu_s*ptens(i,j,k,ie)
! Subprogram not used                     endif
! Subprogram not used                     
! Subprogram not used                     ptens(i,j,k,ie) = ptens_tmp  
! Subprogram not used                     vtens(i,j,1,k,ie)=utens_tmp
! Subprogram not used                     vtens(i,j,2,k,ie)=vtens_tmp
! Subprogram not used                  enddo
! Subprogram not used               enddo
! Subprogram not used            enddo
! Subprogram not used            
! Subprogram not used            pstens(:,:,ie)  =  -nu_p*pstens(:,:,ie)
! Subprogram not used            kptr=0
! Subprogram not used            call edgeVpack(edge3, ptens(:,:,:,ie),nlev,kptr,elem(ie)%desc)
! Subprogram not used            kptr=nlev
! Subprogram not used            call edgeVpack(edge3,vtens(:,:,:,:,ie),2*nlev,kptr,elem(ie)%desc)
! Subprogram not used            kptr=3*nlev
! Subprogram not used            call edgeVpack(edge3,pstens(:,:,ie),1,kptr,elem(ie)%desc)
! Subprogram not used         enddo
! Subprogram not used         
! Subprogram not used         
! Subprogram not used         call bndry_exchangeV(hybrid,edge3)
! Subprogram not used         
! Subprogram not used         do ie=nets,nete
! Subprogram not used            
! Subprogram not used            kptr=0
! Subprogram not used            call edgeVunpack(edge3, ptens(:,:,:,ie), nlev, kptr, elem(ie)%desc)
! Subprogram not used            kptr=nlev
! Subprogram not used            call edgeVunpack(edge3, vtens(:,:,:,:,ie), 2*nlev, kptr, elem(ie)%desc)
! Subprogram not used 
! Subprogram not used            
! Subprogram not used            ! apply inverse mass matrix, accumulate tendencies
! Subprogram not used            do k=1,nlev
! Subprogram not used               vtens(:,:,1,k,ie)=dt*vtens(:,:,1,k,ie)*elem(ie)%rspheremp(:,:)
! Subprogram not used               vtens(:,:,2,k,ie)=dt*vtens(:,:,2,k,ie)*elem(ie)%rspheremp(:,:)
! Subprogram not used               ptens(:,:,k,ie)=dt*ptens(:,:,k,ie)*elem(ie)%rspheremp(:,:)
! Subprogram not used            enddo
! Subprogram not used 
! Subprogram not used            ! apply hypervis to u -> u+utens:  
! Subprogram not used            ! E0 = dpdn * .5*u dot u + dpdn * T  + dpdn*PHIS
! Subprogram not used            ! E1 = dpdn * .5*(u+utens) dot (u+utens) + dpdn * (T-X) + dpdn*PHIS
! Subprogram not used            ! E1-E0:   dpdn (u dot utens) + dpdn .5 utens dot utens   - dpdn X
! Subprogram not used            !      X = (u dot utens) + .5 utens dot utens
! Subprogram not used            !  alt:  (u+utens) dot utens
! Subprogram not used            do k=1,nlev
! Subprogram not used               do j=1,np
! Subprogram not used                  do i=1,np
! Subprogram not used                     ! update v first (gives better results than updating v after heating)
! Subprogram not used                     elem(ie)%state%v(i,j,:,k,nt)=elem(ie)%state%v(i,j,:,k,nt) + &
! Subprogram not used                          vtens(i,j,:,k,ie)
! Subprogram not used 
! Subprogram not used                     v1=elem(ie)%state%v(i,j,1,k,nt)
! Subprogram not used                     v2=elem(ie)%state%v(i,j,2,k,nt)
! Subprogram not used                     heating = (vtens(i,j,1,k,ie)*v1  + vtens(i,j,2,k,ie)*v2 ) 
! Subprogram not used                     elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt) &
! Subprogram not used                          +ptens(i,j,k,ie)-heating/cp
! Subprogram not used 
! Subprogram not used                  enddo
! Subprogram not used               enddo
! Subprogram not used            enddo              
! Subprogram not used 
! Subprogram not used            if (nu_p>0) then
! Subprogram not used               kptr=3*nlev
! Subprogram not used               call edgeVunpack(edge3, pstens(:,:,ie), 1, kptr, elem(ie)%desc)
! Subprogram not used               pstens(:,:,ie)=dt*pstens(:,:,ie)*elem(ie)%rspheremp(:,:)
! Subprogram not used               elem(ie)%state%ps_v(:,:,nt)=elem(ie)%state%ps_v(:,:,nt) + pstens(:,:,ie)
! Subprogram not used            endif
! Subprogram not used            
! Subprogram not used         enddo
! Subprogram not used      enddo
! Subprogram not used   endif
! Subprogram not used   
! Subprogram not used   call t_stopf('advance_hypervis')
! Subprogram not used   
! Subprogram not used   end subroutine advance_hypervis
  
  


  subroutine advance_hypervis_dp(edge3,elem,hvcoord,hybrid,deriv,nt,nets,nete,dt2,eta_ave_w)
  !
  !  take one timestep of:  
  !          u(:,:,:,np) = u(:,:,:,np) +  dt2*nu*laplacian**order ( u )
  !          T(:,:,:,np) = T(:,:,:,np) +  dt2*nu_s*laplacian**order ( T )
  !
  !
  !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
  !
  !
  use dimensions_mod, only : np, np, nlev
  use control_mod, only : nu, nu_div, nu_s, hypervis_order, hypervis_subcycle, nu_p, nu_top, psurf_vis
  use hybrid_mod, only : hybrid_t
  use hybvcoord_mod, only : hvcoord_t
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t, laplace_sphere_wk, vlaplace_sphere_wk
  use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack
  use bndry_mod, only : bndry_exchangev
  use viscosity_mod, only : biharmonic_wk_dp3d
  use physical_constants, only: Cp
!  use time_mod, only : TimeLevel_t
  implicit none
  
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (EdgeBuffer_t)  , intent(inout) :: edge3
  type (derivative_t)  , intent(in) :: deriv
  type (hvcoord_t), intent(in)      :: hvcoord
!  type (TimeLevel_t)   , intent(in) :: tl
  
  real (kind=real_kind) :: dt2
  integer :: nets,nete
  
  ! local
  real (kind=real_kind) :: eta_ave_w  ! weighting for mean flux terms
  real (kind=real_kind) :: dpdn,dpdn0, nu_scale_top
  integer :: k,kptr,i,j,ie,ic,nt
  real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)      :: vtens   
  real (kind=real_kind), dimension(np,np,nlev,nets:nete)        :: ttens
  real (kind=real_kind), dimension(np,np,nlev,nets:nete)        :: dptens
  real (kind=real_kind), dimension(np,np,nlev) :: p
  real (kind=real_kind), dimension(np,np) :: dptemp1,dptemp2
  
  
! NOTE: PGI compiler bug: when using spheremp, rspheremp and ps as pointers to elem(ie)% members,
  !       data is incorrect (offset by a few numbers actually)
  !       removed for now.  
  !       real (kind=real_kind), dimension(:,:), pointer :: spheremp,rspheremp
  !       real (kind=real_kind), dimension(:,:,:), pointer   :: ps
  
  real (kind=real_kind), dimension(np,np) :: lap_t,lap_dp
  real (kind=real_kind), dimension(np,np,2) :: lap_v
  real (kind=real_kind) :: v1,v2,dt,heating,utens_tmp,vtens_tmp,ttens_tmp,dptens_tmp


  if (nu_s == 0 .and. nu == 0 .and. nu_p==0 ) return;
  call t_barrierf('sync_advance_hypervis', hybrid%par%comm)
  call t_startf('advance_hypervis_dp')


  dt=dt2/hypervis_subcycle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  regular viscosity  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (hypervis_order == 1) then
     if (nu_p>0) call abortmp( 'ERROR: hypervis_order == 1 not coded for nu_p>0')
     do ic=1,hypervis_subcycle
        do ie=nets,nete
           
           do k=1,nlev
              lap_t=laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
              lap_v=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)
              ! advace in time.  (note: DSS commutes with time stepping, so we
              ! can time advance and then DSS.  this has the advantage of
              ! not letting any discontinuties accumulate in p,v via roundoff
              do j=1,np
                 do i=1,np             
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt)*elem(ie)%spheremp(i,j)  +  dt*nu_s*lap_t(i,j)
                    elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%state%v(i,j,1,k,nt)*elem(ie)%spheremp(i,j) + dt*nu*lap_v(i,j,1)
                    elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%state%v(i,j,2,k,nt)*elem(ie)%spheremp(i,j) + dt*nu*lap_v(i,j,2)
                 enddo
              enddo
           enddo
           
           kptr=0
           call edgeVpack(edge3, elem(ie)%state%T(:,:,:,nt),nlev,kptr,elem(ie)%desc)
           kptr=nlev
           call edgeVpack(edge3,elem(ie)%state%v(:,:,:,:,nt),2*nlev,kptr,elem(ie)%desc)
        enddo
        
        call bndry_exchangeV(hybrid,edge3)
        
        do ie=nets,nete
           
           kptr=0
           call edgeVunpack(edge3, elem(ie)%state%T(:,:,:,nt), nlev, kptr, elem(ie)%desc)
           kptr=nlev
           call edgeVunpack(edge3, elem(ie)%state%v(:,:,:,:,nt), 2*nlev, kptr, elem(ie)%desc)
           
           ! apply inverse mass matrix
           do k=1,nlev
              do j=1,np
                 do i=1,np             
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%T(i,j,k,nt)
                    elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,1,k,nt)
                    elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,2,k,nt)
                 enddo
              enddo
           enddo
        enddo
     enddo  ! subcycle
  endif
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  hyper viscosity  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! nu_p=0:
!   scale T dissipaton by dp  (conserve IE, dissipate T^2)
! nu_p>0
!   dont scale:  T equation IE dissipation matches (to truncation error) 
!                IE dissipation from continuity equation 
!                (1 deg: to about 0.1 W/m^2)
!
  if (hypervis_order == 2) then
     do ic=1,hypervis_subcycle
        call biharmonic_wk_dp3d(elem,dptens,ttens,vtens,deriv,edge3,hybrid,nt,nets,nete)
        do ie=nets,nete

           ! comptue mean flux
           if (nu_p>0) then
              elem(ie)%derived%dpdiss_ave(:,:,:)=elem(ie)%derived%dpdiss_ave(:,:,:)+&
                   eta_ave_w*elem(ie)%state%dp3d(:,:,:,nt)/hypervis_subcycle
              elem(ie)%derived%dpdiss_biharmonic(:,:,:)=elem(ie)%derived%dpdiss_biharmonic(:,:,:)+&
                   eta_ave_w*dptens(:,:,:,ie)/hypervis_subcycle
           endif
           do k=1,nlev
              ! advace in time.  
              ! note: DSS commutes with time stepping, so we can time advance and then DSS.
              ! note: weak operators alreayd have mass matrix "included"
              
              ! add regular diffusion in top 3 layers:
              if (nu_top>0 .and. k<=3) then
                 lap_t=laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
                 lap_dp=laplace_sphere_wk(elem(ie)%state%dp3d(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
                 lap_v=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)
              endif
              nu_scale_top = 1
              if (k==1) nu_scale_top=4
              if (k==2) nu_scale_top=2
              
              do j=1,np
                 do i=1,np
                    ! biharmonic terms need a negative sign:
                    if (nu_top>0 .and. k<=3) then
                       utens_tmp=(-nu*vtens(i,j,1,k,ie) + nu_scale_top*nu_top*lap_v(i,j,1))
                       vtens_tmp=(-nu*vtens(i,j,2,k,ie) + nu_scale_top*nu_top*lap_v(i,j,2))
                       ttens_tmp=(-nu_s*ttens(i,j,k,ie) + nu_scale_top*nu_top*lap_t(i,j) )
                       dptens_tmp=(-nu_p*dptens(i,j,k,ie) + nu_scale_top*nu_top*lap_dp(i,j) )
                    else
                       utens_tmp=-nu*vtens(i,j,1,k,ie)
                       vtens_tmp=-nu*vtens(i,j,2,k,ie)
                       ttens_tmp=-nu_s*ttens(i,j,k,ie)
                       dptens_tmp=-nu_p*dptens(i,j,k,ie)
                    endif
                    ttens(i,j,k,ie) = ttens_tmp  
                    dptens(i,j,k,ie) =dptens_tmp  
                    vtens(i,j,1,k,ie)=utens_tmp
                    vtens(i,j,2,k,ie)=vtens_tmp
                 enddo
              enddo
           enddo

           
           kptr=0
           call edgeVpack(edge3, ttens(:,:,:,ie),nlev,kptr,elem(ie)%desc)
           kptr=nlev
           call edgeVpack(edge3,vtens(:,:,:,:,ie),2*nlev,kptr,elem(ie)%desc)
           kptr=3*nlev
           call edgeVpack(edge3,dptens(:,:,:,ie),nlev,kptr,elem(ie)%desc)
        enddo
        
        
        call bndry_exchangeV(hybrid,edge3)
        
        do ie=nets,nete
           
           kptr=0
           call edgeVunpack(edge3, ttens(:,:,:,ie), nlev, kptr, elem(ie)%desc)
           kptr=nlev
           call edgeVunpack(edge3, vtens(:,:,:,:,ie), 2*nlev, kptr, elem(ie)%desc)
           kptr=3*nlev
           call edgeVunpack(edge3, dptens(:,:,:,ie), nlev, kptr, elem(ie)%desc)

           
           ! apply inverse mass matrix, accumulate tendencies
           do k=1,nlev
              vtens(:,:,1,k,ie)=dt*vtens(:,:,1,k,ie)*elem(ie)%rspheremp(:,:)
              vtens(:,:,2,k,ie)=dt*vtens(:,:,2,k,ie)*elem(ie)%rspheremp(:,:)
              ttens(:,:,k,ie)=dt*ttens(:,:,k,ie)*elem(ie)%rspheremp(:,:)
              dptens(:,:,k,ie)=dt*dptens(:,:,k,ie)*elem(ie)%rspheremp(:,:)
           enddo

           ! apply hypervis to u -> u+utens:  
           ! E0 = dpdn * .5*u dot u + dpdn * T  + dpdn*PHIS
           ! E1 = dpdn * .5*(u+utens) dot (u+utens) + dpdn * (T-X) + dpdn*PHIS
           ! E1-E0:   dpdn (u dot utens) + dpdn .5 utens dot utens   - dpdn X
           !      X = (u dot utens) + .5 utens dot utens
           !  alt:  (u+utens) dot utens
           do k=1,nlev
              do j=1,np
                 do i=1,np
                    ! update v first (gives better results than updating v after heating)
                    elem(ie)%state%v(i,j,:,k,nt)=elem(ie)%state%v(i,j,:,k,nt) + &
                         vtens(i,j,:,k,ie)

                    v1=elem(ie)%state%v(i,j,1,k,nt)
                    v2=elem(ie)%state%v(i,j,2,k,nt)
                    heating = (vtens(i,j,1,k,ie)*v1  + vtens(i,j,2,k,ie)*v2 ) 
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt) &
                         +ttens(i,j,k,ie)-heating/cp
                    elem(ie)%state%dp3d(i,j,k,nt)=elem(ie)%state%dp3d(i,j,k,nt) + &
                         dptens(i,j,k,ie)
                 enddo
              enddo
           enddo
        enddo
     enddo
  endif
  
  call t_stopf('advance_hypervis_dp')
  
  end subroutine advance_hypervis_dp
  
  
  



! Subprogram not used   subroutine advance_hypervis_lf(edge3,elem,hvcoord,hybrid,deriv,nm1,n0,nt,nets,nete,dt2)
! Subprogram not used   !
! Subprogram not used   !  take one timestep of:  
! Subprogram not used   !          u(:,:,:,np) = u(:,:,:,np) +  dt2*nu*laplacian**order ( u )
! Subprogram not used   !          T(:,:,:,np) = T(:,:,:,np) +  dt2*nu_s*laplacian**order ( T )
! Subprogram not used   !
! Subprogram not used   !
! Subprogram not used   !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
! Subprogram not used   !
! Subprogram not used   !
! Subprogram not used   use dimensions_mod, only : np, np, nlev
! Subprogram not used   use control_mod, only : nu, nu_div, nu_s, hypervis_order, hypervis_subcycle, nu_p, nu_top, psurf_vis
! Subprogram not used   use hybrid_mod, only : hybrid_t
! Subprogram not used   use hybvcoord_mod, only : hvcoord_t
! Subprogram not used   use element_mod, only : element_t
! Subprogram not used   use derivative_mod, only : derivative_t, laplace_sphere_wk, vlaplace_sphere_wk
! Subprogram not used   use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack
! Subprogram not used   use bndry_mod, only : bndry_exchangev
! Subprogram not used   use viscosity_mod, only : biharmonic_wk
! Subprogram not used   use physical_constants, only: Cp
! Subprogram not used !  use time_mod, only : TimeLevel_t
! Subprogram not used   implicit none
! Subprogram not used   
! Subprogram not used   type (hybrid_t)      , intent(in) :: hybrid
! Subprogram not used   type (element_t)     , intent(inout), target :: elem(:)
! Subprogram not used   type (EdgeBuffer_t)  , intent(inout) :: edge3
! Subprogram not used   type (derivative_t)  , intent(in) :: deriv
! Subprogram not used   type (hvcoord_t), intent(in)      :: hvcoord
! Subprogram not used !  type (TimeLevel_t)   , intent(in) :: tl
! Subprogram not used   
! Subprogram not used   real (kind=real_kind) :: dt2
! Subprogram not used   integer :: nets,nete
! Subprogram not used   
! Subprogram not used   ! local
! Subprogram not used   real (kind=real_kind) :: nu_scale, dpdn,dpdn0, nu_scale_top
! Subprogram not used   integer :: k,kptr,i,j,ie,ic,n0,nt,nm1
! Subprogram not used   real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)      :: vtens   
! Subprogram not used   real (kind=real_kind), dimension(np,np,nlev,nets:nete)        :: ptens
! Subprogram not used   real (kind=real_kind), dimension(np,np,nets:nete) :: pstens	
! Subprogram not used   real (kind=real_kind), dimension(np,np,nlev) :: p
! Subprogram not used   real (kind=real_kind), dimension(np,np) :: dXdp
! Subprogram not used   
! Subprogram not used   
! Subprogram not used ! NOTE: PGI compiler bug: when using spheremp, rspheremp and ps as pointers to elem(ie)% members,
! Subprogram not used   !       data is incorrect (offset by a few numbers actually)
! Subprogram not used   !       removed for now.  
! Subprogram not used   !       real (kind=real_kind), dimension(:,:), pointer :: spheremp,rspheremp
! Subprogram not used   !       real (kind=real_kind), dimension(:,:,:), pointer   :: ps
! Subprogram not used   
! Subprogram not used   real (kind=real_kind), dimension(np,np) :: lap_p
! Subprogram not used   real (kind=real_kind), dimension(np,np,2) :: lap_v
! Subprogram not used   real (kind=real_kind) :: v1,v2,dt,heating,utens_tmp,vtens_tmp,ptens_tmp
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   if (nu_s == 0 .and. nu == 0 .and. nu_p==0 ) return;
! Subprogram not used   call t_barrierf('sync_advance_hypervis_lf', hybrid%par%comm)
! Subprogram not used   call t_startf('advance_hypervis_lf')
! Subprogram not used 
! Subprogram not used ! for non-leapfrog,nt=n0=nmt
! Subprogram not used !  
! Subprogram not used !  nm1 = tl%nm1   ! heating term uses U,V at average of nt and nm1 levels
! Subprogram not used !  n0 = tl%n0     ! timelevel used for ps scaling.  use n0 for leapfrog. 
! Subprogram not used !  nt = tl%np1    ! apply viscosity to this timelevel  (np1)
! Subprogram not used   
! Subprogram not used 
! Subprogram not used   dt=dt2/hypervis_subcycle
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used   !  regular viscosity  
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used   if (hypervis_order == 1) then
! Subprogram not used      if (nu_p>0) stop 'ERROR: hypervis_order == 1 not coded for nu_p>0'
! Subprogram not used      do ic=1,hypervis_subcycle
! Subprogram not used         do ie=nets,nete
! Subprogram not used            
! Subprogram not used            do k=1,nlev
! Subprogram not used               lap_p=laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
! Subprogram not used               lap_v=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)
! Subprogram not used               ! advace in time.  (note: DSS commutes with time stepping, so we
! Subprogram not used               ! can time advance and then DSS.  this has the advantage of
! Subprogram not used               ! not letting any discontinuties accumulate in p,v via roundoff
! Subprogram not used               do j=1,np
! Subprogram not used                  do i=1,np             
! Subprogram not used                     elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt)*elem(ie)%spheremp(i,j)  +  dt*nu_s*lap_p(i,j)
! Subprogram not used                     elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%state%v(i,j,1,k,nt)*elem(ie)%spheremp(i,j) + dt*nu*lap_v(i,j,1)
! Subprogram not used                     elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%state%v(i,j,2,k,nt)*elem(ie)%spheremp(i,j) + dt*nu*lap_v(i,j,2)
! Subprogram not used                  enddo
! Subprogram not used               enddo
! Subprogram not used            enddo
! Subprogram not used            
! Subprogram not used            kptr=0
! Subprogram not used            call edgeVpack(edge3, elem(ie)%state%T(:,:,:,nt),nlev,kptr,elem(ie)%desc)
! Subprogram not used            kptr=nlev
! Subprogram not used            call edgeVpack(edge3,elem(ie)%state%v(:,:,:,:,nt),2*nlev,kptr,elem(ie)%desc)
! Subprogram not used         enddo
! Subprogram not used         
! Subprogram not used         call bndry_exchangeV(hybrid,edge3)
! Subprogram not used         
! Subprogram not used         do ie=nets,nete
! Subprogram not used            
! Subprogram not used            kptr=0
! Subprogram not used            call edgeVunpack(edge3, elem(ie)%state%T(:,:,:,nt), nlev, kptr, elem(ie)%desc)
! Subprogram not used            kptr=nlev
! Subprogram not used            call edgeVunpack(edge3, elem(ie)%state%v(:,:,:,:,nt), 2*nlev, kptr, elem(ie)%desc)
! Subprogram not used            
! Subprogram not used            ! apply inverse mass matrix
! Subprogram not used            do k=1,nlev
! Subprogram not used               do j=1,np
! Subprogram not used                  do i=1,np             
! Subprogram not used                     elem(ie)%state%T(i,j,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%T(i,j,k,nt)
! Subprogram not used                     elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,1,k,nt)
! Subprogram not used                     elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,2,k,nt)
! Subprogram not used                  enddo
! Subprogram not used               enddo
! Subprogram not used            enddo
! Subprogram not used         enddo
! Subprogram not used      enddo  ! subcycle
! Subprogram not used   endif
! Subprogram not used   
! Subprogram not used   
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used   !  hyper viscosity  
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used   if (hypervis_order == 2) then
! Subprogram not used      do ic=1,hypervis_subcycle
! Subprogram not used         call biharmonic_wk(elem,pstens,ptens,vtens,deriv,edge3,hybrid,nt,nets,nete)
! Subprogram not used         do ie=nets,nete
! Subprogram not used 
! Subprogram not used            nu_scale=1
! Subprogram not used            do k=1,nlev
! Subprogram not used               ! advace in time.  
! Subprogram not used               ! note: DSS commutes with time stepping, so we can time advance and then DSS.
! Subprogram not used               ! note: weak operators alreayd have mass matrix "included"
! Subprogram not used               
! Subprogram not used               ! add regular diffusion in top 3 layers:
! Subprogram not used               if (nu_top>0 .and. k<=3) then
! Subprogram not used                  lap_p=laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
! Subprogram not used                  lap_v=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)
! Subprogram not used               endif
! Subprogram not used               nu_scale_top = 1
! Subprogram not used               if (k==1) nu_scale_top=4
! Subprogram not used               if (k==2) nu_scale_top=2
! Subprogram not used               
! Subprogram not used               do j=1,np
! Subprogram not used                  do i=1,np
! Subprogram not used                     if (psurf_vis==0) then
! Subprogram not used                        ! normalize so as to conserve IE  (not needed when using p-surface viscosity)
! Subprogram not used                        ! scale velosity by 1/rho (normalized to be O(1))
! Subprogram not used                        ! dp/dn = O(ps0)*O(delta_eta) = O(ps0)/O(nlev)
! Subprogram not used                        dpdn = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
! Subprogram not used                             ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,n0)  ! nt ?
! Subprogram not used                        dpdn0 = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
! Subprogram not used                             ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*hvcoord%ps0
! Subprogram not used                        nu_scale = dpdn0/dpdn
! Subprogram not used                     endif
! Subprogram not used 
! Subprogram not used                     ! biharmonic terms need a negative sign:
! Subprogram not used                     if (nu_top>0 .and. k<=3) then
! Subprogram not used                        utens_tmp=(-nu*vtens(i,j,1,k,ie) + nu_scale_top*nu_top*lap_v(i,j,1))
! Subprogram not used                        vtens_tmp=(-nu*vtens(i,j,2,k,ie) + nu_scale_top*nu_top*lap_v(i,j,2))
! Subprogram not used                        ptens_tmp=nu_scale*(-nu_s*ptens(i,j,k,ie) + nu_scale_top*nu_top*lap_p(i,j) )
! Subprogram not used                     else
! Subprogram not used                        utens_tmp=-nu*vtens(i,j,1,k,ie)
! Subprogram not used                        vtens_tmp=-nu*vtens(i,j,2,k,ie)
! Subprogram not used                        ptens_tmp=-nu_scale*nu_s*ptens(i,j,k,ie)
! Subprogram not used                     endif
! Subprogram not used                     
! Subprogram not used                     ptens(i,j,k,ie) = ptens_tmp  
! Subprogram not used                     vtens(i,j,1,k,ie)=utens_tmp
! Subprogram not used                     vtens(i,j,2,k,ie)=vtens_tmp
! Subprogram not used                  enddo
! Subprogram not used               enddo
! Subprogram not used            enddo
! Subprogram not used            
! Subprogram not used            pstens(:,:,ie)  =  -nu_p*pstens(:,:,ie)
! Subprogram not used            kptr=0
! Subprogram not used            call edgeVpack(edge3, ptens(:,:,:,ie),nlev,kptr,elem(ie)%desc)
! Subprogram not used            kptr=nlev
! Subprogram not used            call edgeVpack(edge3,vtens(:,:,:,:,ie),2*nlev,kptr,elem(ie)%desc)
! Subprogram not used            kptr=3*nlev
! Subprogram not used            call edgeVpack(edge3,pstens(:,:,ie),1,kptr,elem(ie)%desc)
! Subprogram not used         enddo
! Subprogram not used         
! Subprogram not used         
! Subprogram not used         call bndry_exchangeV(hybrid,edge3)
! Subprogram not used         
! Subprogram not used         do ie=nets,nete
! Subprogram not used            
! Subprogram not used            kptr=0
! Subprogram not used            call edgeVunpack(edge3, ptens(:,:,:,ie), nlev, kptr, elem(ie)%desc)
! Subprogram not used            kptr=nlev
! Subprogram not used            call edgeVunpack(edge3, vtens(:,:,:,:,ie), 2*nlev, kptr, elem(ie)%desc)
! Subprogram not used            kptr=3*nlev
! Subprogram not used            call edgeVunpack(edge3, pstens(:,:,ie), 1, kptr, elem(ie)%desc)
! Subprogram not used 
! Subprogram not used            if (psurf_vis == 1 ) then
! Subprogram not used               ! apply p-surface correction
! Subprogram not used               do k=1,nlev
! Subprogram not used                  p(:,:,k)   = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem(ie)%state%ps_v(:,:,nt)
! Subprogram not used               enddo
! Subprogram not used               do k=1,nlev
! Subprogram not used                  if (k.eq.1) then
! Subprogram not used                     ! no correction needed
! Subprogram not used                  else if (k.eq.nlev) then
! Subprogram not used                     ! one-sided difference
! Subprogram not used                     dXdp = (elem(ie)%state%T(:,:,k,nt) - elem(ie)%state%T(:,:,k-1,nt)) / &
! Subprogram not used                         (p(:,:,k)-p(:,:,k-1)) 
! Subprogram not used                     ptens(:,:,k,ie) = ptens(:,:,k,ie) - dXdp(:,:)*hvcoord%hybm(k)*pstens(:,:,ie)
! Subprogram not used                  else
! Subprogram not used                     dXdp = (elem(ie)%state%T(:,:,k+1,nt) - elem(ie)%state%T(:,:,k-1,nt)) / &
! Subprogram not used                          (p(:,:,k+1)-p(:,:,k-1)) 
! Subprogram not used                     ptens(:,:,k,ie) = ptens(:,:,k,ie) - dXdp(:,:)*hvcoord%hybm(k)*pstens(:,:,ie)
! Subprogram not used                  endif
! Subprogram not used               enddo
! Subprogram not used            endif
! Subprogram not used 
! Subprogram not used            
! Subprogram not used            ! apply inverse mass matrix, accumulate tendencies
! Subprogram not used            do k=1,nlev
! Subprogram not used               do j=1,np
! Subprogram not used                  do i=1,np
! Subprogram not used 
! Subprogram not used                     elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%state%v(i,j,1,k,nt) + &
! Subprogram not used                          dt*elem(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
! Subprogram not used                     elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%state%v(i,j,2,k,nt) +  &
! Subprogram not used                          dt*elem(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie) 
! Subprogram not used                     
! Subprogram not used                     ! better E conservation if we use v after adding in vtens:
! Subprogram not used                     v1=.5*(elem(ie)%state%v(i,j,1,k,nt)+elem(ie)%state%v(i,j,1,k,nm1))
! Subprogram not used                     v2=.5*(elem(ie)%state%v(i,j,2,k,nt)+elem(ie)%state%v(i,j,2,k,nm1))
! Subprogram not used                     heating = (vtens(i,j,1,k,ie)*v1  + vtens(i,j,2,k,ie)*v2 )
! Subprogram not used 
! Subprogram not used                     elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt)     + &
! Subprogram not used                          dt*elem(ie)%rspheremp(i,j)*(cp*ptens(i,j,k,ie) - heating)/cp
! Subprogram not used                     
! Subprogram not used                  enddo
! Subprogram not used               enddo
! Subprogram not used            enddo
! Subprogram not used            elem(ie)%state%ps_v(:,:,nt)=elem(ie)%state%ps_v(:,:,nt) + dt*elem(ie)%rspheremp(:,:)*pstens(:,:,ie)
! Subprogram not used         enddo
! Subprogram not used      enddo
! Subprogram not used   endif
! Subprogram not used   
! Subprogram not used   call t_stopf('advance_hypervis_lf')
! Subprogram not used   
! Subprogram not used   end subroutine advance_hypervis_lf
  

  subroutine compute_and_apply_rhs(np1,nm1,n0,qn0,dt2,elem,hvcoord,hybrid,&
       deriv,nets,nete,compute_diagnostics,eta_ave_w)
  ! ===================================
  ! compute the RHS, accumulate into u(np1) and apply DSS
  !
  !           u(np1) = u(nm1) + dt2*DSS[ RHS(u(n0)) ]
  !
  ! This subroutine is normally called to compute a leapfrog timestep
  ! but by adjusting np1,nm1,n0 and dt2, many other timesteps can be
  ! accomodated.  For example, setting nm1=np1=n0 this routine will
  ! take a forward euler step, overwriting the input with the output.
  !
  !    qn0 = timelevel used to access Qdp() in order to compute virtual Temperature
  !          qn0=-1 for the dry case
  !
  ! if  dt2<0, then the DSS'd RHS is returned in timelevel np1
  !
  ! Combining the RHS and DSS pack operation in one routine 
  ! allows us to fuse these two loops for more cache reuse
  !
  ! Combining the dt advance and DSS unpack operation in one routine 
  ! allows us to fuse these two loops for more cache reuse
  !
  ! note: for prescribed velocity case, velocity will be computed at
  ! "real_time", which should be the time of timelevel n0.  
  !
  !
  ! ===================================
  use kinds, only : real_kind
  use dimensions_mod, only : np, np, nlev
  use hybrid_mod, only : hybrid_t
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t, divergence_sphere, gradient_sphere, vorticity_sphere
  use edge_mod, only : edgevpack, edgevunpack
  use bndry_mod, only : bndry_exchangev
  use control_mod, only : moisture, qsplit, use_cpstar, rsplit
  use hybvcoord_mod, only : hvcoord_t

  use physical_constants, only : cp, cpwater_vapor, Rgas, kappa
  use physics_mod, only : virtual_specific_heat, virtual_temperature
  use prim_si_mod, only : preq_vertadv, preq_omega_ps, preq_hydrostatic


  implicit none
  integer :: np1,nm1,n0,qn0,nets,nete
  real*8 :: dt2
  logical  :: compute_diagnostics

  type (hvcoord_t)     , intent(in) :: hvcoord
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv
  real (kind=real_kind) :: eta_ave_w  ! weighting for eta_dot_dpdn mean flux

  ! local
  real (kind=real_kind), pointer, dimension(:,:)      :: ps         ! surface pressure for current tiime level
  real (kind=real_kind), pointer, dimension(:,:,:)   :: phi

  real (kind=real_kind), dimension(np,np,nlev)   :: omega_p       
  real (kind=real_kind), dimension(np,np,nlev)   :: T_v         
  real (kind=real_kind), dimension(np,np,nlev)   :: divdp
  real (kind=real_kind), dimension(np,np,nlev+1)   :: eta_dot_dpdn  ! half level vertical velocity on p-grid
  real (kind=real_kind), dimension(np,np)      :: sdot_sum   ! temporary field
  real (kind=real_kind), dimension(np,np,2)    :: vtemp     ! generic gradient storage
  real (kind=real_kind), dimension(np,np)      :: vgrad_T    ! v.grad(T)
  real (kind=real_kind), dimension(np,np)      :: Ephi       ! kinetic energy + PHI term
  real (kind=real_kind), dimension(np,np,2)      :: grad_ps    ! lat-lon coord version
  real (kind=real_kind), dimension(np,np,2,nlev) :: grad_p    
  real (kind=real_kind), dimension(np,np,nlev)   :: vort       ! vorticity
  real (kind=real_kind), dimension(np,np,nlev)   :: p          ! pressure
  real (kind=real_kind), dimension(np,np,nlev)   :: dp         ! delta pressure
  real (kind=real_kind), dimension(np,np,nlev)   :: rdp        ! inverse of delta pressure
  real (kind=real_kind), dimension(np,np,nlev)   :: T_vadv     ! temperature vertical advection
  real (kind=real_kind), dimension(np,np,nlev)   :: vgrad_p    ! v.grad(p)
  real (kind=real_kind), dimension(np,np,nlev+1) :: ph               ! half level pressures on p-grid
  real (kind=real_kind), dimension(np,np,2,nlev) :: v_vadv   ! velocity vertical advection
  real (kind=real_kind) ::  kappa_star(np,np,nlev)
  real (kind=real_kind) ::  vtens1(np,np,nlev)
  real (kind=real_kind) ::  vtens2(np,np,nlev)
  real (kind=real_kind) ::  ttens(np,np,nlev)

  real (kind=real_kind) ::  cp2,cp_ratio,E,de,Qt,v1,v2
  real (kind=real_kind) ::  glnps1,glnps2,gpterm
  integer :: i,j,k,kptr,ie

  call t_barrierf('sync_compute_and_apply_rhs', hybrid%par%comm)
  call t_startf('compute_and_apply_rhs')
  do ie=nets,nete
     !ps => elem(ie)%state%ps_v(:,:,n0)
     phi => elem(ie)%derived%phi(:,:,:)
     
     ! ==================================================
     ! compute pressure (p) on half levels from ps 
     ! using the hybrid coordinates relationship, i.e.
     ! e.g. equation (3.a.92) of the CCM-2 description, 
     ! (NCAR/TN-382+STR), June 1993, p. 24.
     ! ==================================================
     ! vertically eulerian only needs grad(ps)
     if (rsplit==0) &
          grad_ps = gradient_sphere(elem(ie)%state%ps_v(:,:,n0),deriv,elem(ie)%Dinv)


     ! ============================
     ! compute p and delta p
     ! ============================
!#if (defined ELEMENT_OPENMP)
!!$omp parallel do private(k,i,j)
!#endif
!     do k=1,nlev+1
!       ph(:,:,k)   = hvcoord%hyai(k)*hvcoord%ps0 + hvcoord%hybi(k)*elem(ie)%state%ps_v(:,:,n0)
!     end do

     do k=1,nlev
        if (rsplit==0) then
           dp(:,:,k) = (hvcoord%hyai(k+1)*hvcoord%ps0 + hvcoord%hybi(k+1)*elem(ie)%state%ps_v(:,:,n0)) &
                - (hvcoord%hyai(k)*hvcoord%ps0 + hvcoord%hybi(k)*elem(ie)%state%ps_v(:,:,n0))
           p(:,:,k)   = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem(ie)%state%ps_v(:,:,n0)
           grad_p(:,:,:,k) = hvcoord%hybm(k)*grad_ps(:,:,:)
        else
           ! vertically lagrangian code: we advect dp3d instead of ps_v
           ! we also need grad(p) at all levels (not just grad(ps))
           !p(k)= hyam(k)*ps0 + hybm(k)*ps
           !    = .5*(hyai(k+1)+hyai(k))*ps0 + .5*(hybi(k+1)+hybi(k))*ps
           !    = .5*(ph(k+1) + ph(k) )  = ph(k) + dp(k)/2
           !
           ! p(k+1)-p(k) = ph(k+1)-ph(k) + (dp(k+1)-dp(k))/2
           !             = dp(k) + (dp(k+1)-dp(k))/2 = (dp(k+1)+dp(k))/2
           dp(:,:,k) = elem(ie)%state%dp3d(:,:,k,n0)
           if (k==1) then
              p(:,:,k)=hvcoord%hyai(k)*hvcoord%ps0 + dp(:,:,k)/2
           else
              p(:,:,k)=p(:,:,k-1) + dp(:,:,k-1)/2 + dp(:,:,k)/2
           endif
           grad_p(:,:,:,k) = gradient_sphere(p(:,:,k),deriv,elem(ie)%Dinv)
        endif

        rdp(:,:,k) = 1.0D0/dp(:,:,k)

        ! ============================
        ! compute vgrad_lnps 
        ! ============================
        do j=1,np
           do i=1,np
              v1 = elem(ie)%state%v(i,j,1,k,n0)
              v2 = elem(ie)%state%v(i,j,2,k,n0)
!              vgrad_p(i,j,k) = &
!                   hvcoord%hybm(k)*(v1*grad_ps(i,j,1) + v2*grad_ps(i,j,2)) 
              vgrad_p(i,j,k) = (v1*grad_p(i,j,1,k) + v2*grad_p(i,j,2,k)) 
              vtemp(i,j,1) = v1*dp(i,j,k)
              vtemp(i,j,2) = v2*dp(i,j,k)
           end do
        end do


      
        ! ================================
        ! Accumulate mean Vel_rho flux in vn0
        ! ================================
        elem(ie)%derived%vn0(:,:,:,k)=elem(ie)%derived%vn0(:,:,:,k)+eta_ave_w*vtemp(:,:,:)


        ! =========================================
        !
        ! Compute relative vorticity and divergence
        !
        ! =========================================
        divdp(:,:,k)=divergence_sphere(vtemp,deriv,elem(ie))
        vort(:,:,k)=vorticity_sphere(elem(ie)%state%v(:,:,:,k,n0),deriv,elem(ie))

     enddo
     
     ! compute T_v for timelevel n0
     !if ( moisture /= "dry") then
     if (qn0 == -1 ) then
        do k=1,nlev
           do j=1,np
              do i=1,np
                 T_v(i,j,k) = elem(ie)%state%T(i,j,k,n0)
                 kappa_star(i,j,k) = kappa
              end do
           end do
        end do
     else
        do k=1,nlev
           do j=1,np
              do i=1,np
                 ! Qt = elem(ie)%state%Q(i,j,k,1) 
                 Qt = elem(ie)%state%Qdp(i,j,k,1,qn0)/dp(i,j,k)
                 T_v(i,j,k) = Virtual_Temperature(elem(ie)%state%T(i,j,k,n0),Qt)
                 if (use_cpstar==1) then
                    kappa_star(i,j,k) =  Rgas/Virtual_Specific_Heat(Qt)
                 else
                    kappa_star(i,j,k) = kappa
                 endif
              end do
           end do
        end do
     end if
     
     
     ! ====================================================
     ! Compute Hydrostatic equation, modeld after CCM-3
     ! ====================================================
     !call geopotential_t(p,dp,T_v,Rgas,phi)
     call preq_hydrostatic(phi,elem(ie)%state%phis,T_v,p,dp)
     
     ! ====================================================
     ! Compute omega_p according to CCM-3 
     ! ====================================================
     call preq_omega_ps(omega_p,hvcoord,p,vgrad_p,divdp)

     
     ! ==================================================
     ! zero partial sum for accumulating sum
     !    (div(v_k) + v_k.grad(lnps))*dsigma_k = div( v dp )
     ! used by eta_dot_dpdn and lnps tendency
     ! ==================================================
     sdot_sum=0
     

     ! ==================================================
     ! Compute eta_dot_dpdn 
     ! save sdot_sum as this is the -RHS of ps_v equation
     ! ==================================================
     if (rsplit>0) then
        ! VERTICALLY LAGRANGIAN:   no vertical motion
        eta_dot_dpdn=0    
        T_vadv=0          
        v_vadv=0
     else
        do k=1,nlev
           ! ==================================================
           ! add this term to PS equation so we exactly conserve dry mass
           ! ==================================================
           sdot_sum(:,:) = sdot_sum(:,:) + divdp(:,:,k) 
           eta_dot_dpdn(:,:,k+1) = sdot_sum(:,:)      
        end do
     
     
        ! ===========================================================
        ! at this point, eta_dot_dpdn contains integral_etatop^eta[ divdp ]
        ! compute at interfaces:
        !    eta_dot_dpdn = -dp/dt - integral_etatop^eta[ divdp ]
        ! for reference: at mid layers we have:
        !    omega = v grad p  - integral_etatop^eta[ divdp ]
        ! ===========================================================
        do k=1,nlev-1
           eta_dot_dpdn(:,:,k+1) = hvcoord%hybi(k+1)*sdot_sum(:,:) - eta_dot_dpdn(:,:,k+1)
        end do
        
        eta_dot_dpdn(:,:,1     ) = 0.0D0
        eta_dot_dpdn(:,:,nlev+1) = 0.0D0
        
        ! ===========================================================
        ! Compute vertical advection of T and v from eq. CCM2 (3.b.1)
        ! ==============================================
        call preq_vertadv(elem(ie)%state%T(:,:,:,n0),elem(ie)%state%v(:,:,:,:,n0), &
             eta_dot_dpdn,rdp,T_vadv,v_vadv)
     endif


     ! ================================
     ! accumulate mean vertical flux:
     ! ================================
     do k=1,nlev  !  Loop index added (AAM)
        elem(ie)%derived%eta_dot_dpdn(:,:,k) = &
             elem(ie)%derived%eta_dot_dpdn(:,:,k) + eta_ave_w*eta_dot_dpdn(:,:,k)
        elem(ie)%derived%omega_p(:,:,k) = &
             elem(ie)%derived%omega_p(:,:,k) + eta_ave_w*omega_p(:,:,k)
     enddo
     elem(ie)%derived%eta_dot_dpdn(:,:,nlev+1) = &
          elem(ie)%derived%eta_dot_dpdn(:,:,nlev+1) + eta_ave_w*eta_dot_dpdn(:,:,nlev+1)
     


     
     
     ! ==============================================
     ! Compute phi + kinetic energy term: 10*nv*nv Flops
     ! ==============================================
     do k=1,nlev   
        do j=1,np
           do i=1,np
              v1     = elem(ie)%state%v(i,j,1,k,n0)
              v2     = elem(ie)%state%v(i,j,2,k,n0)
              E = 0.5D0*( v1*v1 + v2*v2 )
              Ephi(i,j)=E+phi(i,j,k)+elem(ie)%derived%pecnd(i,j,k)
           end do
        end do
        ! ================================================
        ! compute gradp term (ps/p)*(dp/dps)*T
        ! ================================================
        vtemp(:,:,:)   = gradient_sphere(elem(ie)%state%T(:,:,k,n0),deriv,elem(ie)%Dinv)
        do j=1,np
           do i=1,np
              v1     = elem(ie)%state%v(i,j,1,k,n0)
              v2     = elem(ie)%state%v(i,j,2,k,n0)
              vgrad_T(i,j) =  v1*vtemp(i,j,1) + v2*vtemp(i,j,2) 
           end do
        end do
        
        
        ! vtemp = grad ( E + PHI )
        vtemp = gradient_sphere(Ephi(:,:),deriv,elem(ie)%Dinv)
        
        do j=1,np
           do i=1,np
!              gpterm = hvcoord%hybm(k)*T_v(i,j,k)/p(i,j,k)
!              glnps1 = Rgas*gpterm*grad_ps(i,j,1)
!              glnps2 = Rgas*gpterm*grad_ps(i,j,2)
              gpterm = T_v(i,j,k)/p(i,j,k)
              glnps1 = Rgas*gpterm*grad_p(i,j,1,k)
              glnps2 = Rgas*gpterm*grad_p(i,j,2,k)
              
              v1     = elem(ie)%state%v(i,j,1,k,n0)
              v2     = elem(ie)%state%v(i,j,2,k,n0)
              
              vtens1(i,j,k) =   - v_vadv(i,j,1,k)                           &
                   + v2*(elem(ie)%fcor(i,j) + vort(i,j,k))        &
                   - vtemp(i,j,1) - glnps1   
              
              vtens2(i,j,k) =   - v_vadv(i,j,2,k)                            &
                   - v1*(elem(ie)%fcor(i,j) + vort(i,j,k))        &
                   - vtemp(i,j,2) - glnps2   
              
              ttens(i,j,k)  = - T_vadv(i,j,k) - vgrad_T(i,j) + kappa_star(i,j,k)*T_v(i,j,k)*omega_p(i,j,k)

           end do
        end do

     end do

     ! =========================================================
     ! local element timestep, store in np1.  
     ! note that we allow np1=n0 or nm1
     ! apply mass matrix
     ! =========================================================
     if (dt2<0) then
        ! calling program just wanted DSS'd RHS, skip time advance
        do k=1,nlev
           elem(ie)%state%v(:,:,1,k,np1) = elem(ie)%spheremp(:,:)*vtens1(:,:,k) 
           elem(ie)%state%v(:,:,2,k,np1) = elem(ie)%spheremp(:,:)*vtens2(:,:,k) 
           elem(ie)%state%T(:,:,k,np1) = elem(ie)%spheremp(:,:)*ttens(:,:,k)
           if (rsplit>0) &
                elem(ie)%state%dp3d(:,:,k,np1) = -elem(ie)%spheremp(:,:)*&
                (divdp(:,:,k) + eta_dot_dpdn(:,:,k+1)-eta_dot_dpdn(:,:,k)) 
        enddo
        elem(ie)%state%ps_v(:,:,np1) = -elem(ie)%spheremp(:,:)*sdot_sum
     else
        do k=1,nlev
           elem(ie)%state%v(:,:,1,k,np1) = elem(ie)%spheremp(:,:)*( elem(ie)%state%v(:,:,1,k,nm1) + dt2*vtens1(:,:,k) )
           elem(ie)%state%v(:,:,2,k,np1) = elem(ie)%spheremp(:,:)*( elem(ie)%state%v(:,:,2,k,nm1) + dt2*vtens2(:,:,k) )
           elem(ie)%state%T(:,:,k,np1) = elem(ie)%spheremp(:,:)*(elem(ie)%state%T(:,:,k,nm1) + dt2*ttens(:,:,k))
           if (rsplit>0) &
                elem(ie)%state%dp3d(:,:,k,np1) = elem(ie)%spheremp(:,:)*&
                (elem(ie)%state%dp3d(:,:,k,nm1)-dt2*&
                (divdp(:,:,k) + eta_dot_dpdn(:,:,k+1)-eta_dot_dpdn(:,:,k)))
        enddo
        elem(ie)%state%ps_v(:,:,np1) = elem(ie)%spheremp(:,:)*( elem(ie)%state%ps_v(:,:,nm1) - dt2*sdot_sum )

     endif

     
     ! =========================================================
     !
     ! Pack ps(np1), T, and v tendencies into comm buffer
     !
     ! =========================================================
     kptr=0
     call edgeVpack(edge3p1, elem(ie)%state%ps_v(:,:,np1),1,kptr,elem(ie)%desc)
     
     kptr=1
     call edgeVpack(edge3p1, elem(ie)%state%T(:,:,:,np1),nlev,kptr,elem(ie)%desc)
     
     kptr=nlev+1
     call edgeVpack(edge3p1, elem(ie)%state%v(:,:,:,:,np1),2*nlev,kptr,elem(ie)%desc)

     if (rsplit>0) then
        kptr=kptr+2*nlev
        call edgeVpack(edge3p1, elem(ie)%state%dp3d(:,:,:,np1),nlev,kptr,elem(ie)%desc)
     endif
  end do
  
  ! =============================================================
    ! Insert communications here: for shared memory, just a single
  ! sync is required
  ! =============================================================
  call bndry_exchangeV(hybrid,edge3p1)
  do ie=nets,nete
     ! ===========================================================
     ! Unpack the edges for vgrad_T and v tendencies...
     ! ===========================================================
     kptr=0
     call edgeVunpack(edge3p1, elem(ie)%state%ps_v(:,:,np1), 1, kptr, elem(ie)%desc)
     
     kptr=1
     call edgeVunpack(edge3p1, elem(ie)%state%T(:,:,:,np1), nlev, kptr, elem(ie)%desc)
     
     kptr=nlev+1
     call edgeVunpack(edge3p1, elem(ie)%state%v(:,:,:,:,np1), 2*nlev, kptr, elem(ie)%desc)
     
     if (rsplit>0) then
        kptr=kptr+2*nlev
        call edgeVunpack(edge3p1, elem(ie)%state%dp3d(:,:,:,np1),nlev,kptr,elem(ie)%desc)
     endif
     
     ! ====================================================
     ! Scale tendencies by inverse mass matrix
     ! ====================================================

     do k=1,nlev
        elem(ie)%state%T(:,:,k,np1)   = elem(ie)%rspheremp(:,:)*elem(ie)%state%T(:,:,k,np1)
        elem(ie)%state%v(:,:,1,k,np1) = elem(ie)%rspheremp(:,:)*elem(ie)%state%v(:,:,1,k,np1)
        elem(ie)%state%v(:,:,2,k,np1) = elem(ie)%rspheremp(:,:)*elem(ie)%state%v(:,:,2,k,np1)
     end do

     if (rsplit>0) then
        ! vertically lagrangian: complete dp3d timestep:
        do k=1,nlev
           elem(ie)%state%dp3d(:,:,k,np1)= elem(ie)%rspheremp(:,:)*elem(ie)%state%dp3d(:,:,k,np1)
        enddo
        ! when debugging: also update ps_v
        !elem(ie)%state%ps_v(:,:,np1) = elem(ie)%rspheremp(:,:)*elem(ie)%state%ps_v(:,:,np1)
     else
        ! vertically eulerian: complete ps_v timestep:
        elem(ie)%state%ps_v(:,:,np1) = elem(ie)%rspheremp(:,:)*elem(ie)%state%ps_v(:,:,np1)
     endif

  end do


  call t_stopf('compute_and_apply_rhs')
  end subroutine compute_and_apply_rhs
  
! Subprogram not used   subroutine smooth_phis(phis,elem,hybrid,deriv,nets,nete,minf,numcycle)
! Subprogram not used   use dimensions_mod, only : np, np, nlev
! Subprogram not used   use control_mod, only : smooth_phis_nudt, hypervis_scaling
! Subprogram not used   use hybrid_mod, only : hybrid_t
! Subprogram not used   use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack, edgevunpackmax, edgevunpackmin
! Subprogram not used   use bndry_mod, only : bndry_exchangev
! Subprogram not used   use element_mod, only : element_t
! Subprogram not used   use derivative_mod, only : derivative_t , laplace_sphere_wk
! Subprogram not used   use time_mod, only : TimeLevel_t
! Subprogram not used   implicit none
! Subprogram not used   
! Subprogram not used   integer :: nets,nete
! Subprogram not used   real (kind=real_kind), dimension(np,np,nets:nete), intent(inout)   :: phis
! Subprogram not used   type (hybrid_t)      , intent(in) :: hybrid
! Subprogram not used   type (element_t)     , intent(inout), target :: elem(:)
! Subprogram not used   type (derivative_t)  , intent(in) :: deriv
! Subprogram not used   real (kind=real_kind), intent(in)   :: minf
! Subprogram not used   integer,               intent(in) :: numcycle
! Subprogram not used   
! Subprogram not used   ! local 
! Subprogram not used   real (kind=real_kind), dimension(np,np,nets:nete) :: pstens	
! Subprogram not used   real (kind=real_kind), dimension(nets:nete) :: pmin,pmax
! Subprogram not used   real (kind=real_kind) :: mx,mn
! Subprogram not used   integer :: nt,ie,ic,i,j,order,order_max, iuse
! Subprogram not used   logical :: use_var_coef
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   ! compute local element neighbor min/max
! Subprogram not used   do ie=nets,nete
! Subprogram not used      pstens(:,:,ie)=minval(phis(:,:,ie))
! Subprogram not used      call edgeVpack(edge3p1,pstens(:,:,ie),1,0,elem(ie)%desc)
! Subprogram not used   enddo
! Subprogram not used   call bndry_exchangeV(hybrid,edge3p1)
! Subprogram not used   do ie=nets,nete
! Subprogram not used      call edgeVunpackMin(edge3p1, pstens(:,:,ie), 1, 0, elem(ie)%desc)
! Subprogram not used      pmin(ie)=minval(pstens(:,:,ie))
! Subprogram not used   enddo
! Subprogram not used   do ie=nets,nete
! Subprogram not used      pstens(:,:,ie)=maxval(phis(:,:,ie))
! Subprogram not used      call edgeVpack(edge3p1,pstens(:,:,ie),1,0,elem(ie)%desc)
! Subprogram not used   enddo
! Subprogram not used   call bndry_exchangeV(hybrid,edge3p1)
! Subprogram not used   do ie=nets,nete
! Subprogram not used      call edgeVunpackMax(edge3p1, pstens(:,:,ie), 1, 0, elem(ie)%desc)
! Subprogram not used      pmax(ie)=maxval(pstens(:,:,ie))
! Subprogram not used   enddo
! Subprogram not used 
! Subprogram not used   ! order = 1   grad^2 laplacian
! Subprogram not used   ! order = 2   grad^4 (need to add a negative sign)
! Subprogram not used   ! order = 3   grad^6
! Subprogram not used   ! order = 4   grad^8 (need to add a negative sign)  
! Subprogram not used   order_max = 1
! Subprogram not used 
! Subprogram not used   
! Subprogram not used   use_var_coef=.true.
! Subprogram not used   if (hypervis_scaling/=0) then
! Subprogram not used      ! for tensorHV option, we turn off the tensor except for *last* laplace operator                             
! Subprogram not used      use_var_coef=.false.
! Subprogram not used      if (hypervis_scaling>=3) then
! Subprogram not used         ! with a 3.2 or 4 scaling, assume hyperviscosity                                                          
! Subprogram not used         order_max = 2
! Subprogram not used      endif
! Subprogram not used   endif
! Subprogram not used   
! Subprogram not used 
! Subprogram not used   do ic=1,numcycle
! Subprogram not used      pstens=phis
! Subprogram not used 
! Subprogram not used      do order=1,order_max-1
! Subprogram not used 
! Subprogram not used         do ie=nets,nete
! Subprogram not used            pstens(:,:,ie)=laplace_sphere_wk(pstens(:,:,ie),deriv,elem(ie),var_coef=use_var_coef)
! Subprogram not used            call edgeVpack(edge3p1,pstens(:,:,ie),1,0,elem(ie)%desc)
! Subprogram not used         enddo
! Subprogram not used         call bndry_exchangeV(hybrid,edge3p1)
! Subprogram not used         do ie=nets,nete
! Subprogram not used            call edgeVunpack(edge3p1, pstens(:,:,ie), 1, 0, elem(ie)%desc)
! Subprogram not used            pstens(:,:,ie)=pstens(:,:,ie)*elem(ie)%rspheremp(:,:)
! Subprogram not used         enddo
! Subprogram not used      enddo
! Subprogram not used      do ie=nets,nete
! Subprogram not used         pstens(:,:,ie)=laplace_sphere_wk(pstens(:,:,ie),deriv,elem(ie),var_coef=.true.)
! Subprogram not used      enddo
! Subprogram not used      if (mod(order_max,2)==0) pstens=-pstens
! Subprogram not used 
! Subprogram not used      do ie=nets,nete
! Subprogram not used         !  ps(t+1) = ps(t) + Minv * DSS * M * RHS
! Subprogram not used         !  ps(t+1) = Minv * DSS * M [ ps(t) +  RHS ]
! Subprogram not used         ! but output of biharminc_wk is of the form M*RHS.  rewrite as:
! Subprogram not used         !  ps(t+1) = Minv * DSS * M [ ps(t) +  M*RHS/M ]
! Subprogram not used         ! so we can apply limiter to ps(t) +  (M*RHS)/M
! Subprogram not used         mn=pmin(ie)
! Subprogram not used         mx=pmax(ie)
! Subprogram not used         iuse = numcycle+1  ! always apply min/max limiter
! Subprogram not used         phis(:,:,ie)=phis(:,:,ie) + &
! Subprogram not used            smooth_phis_nudt*pstens(:,:,ie)/elem(ie)%spheremp(:,:)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used         ! remove new extrema.  could use conservative reconstruction from advection
! Subprogram not used         ! but no reason to conserve mean PHI.
! Subprogram not used         if (ic < iuse) then
! Subprogram not used         do i=1,np
! Subprogram not used         do j=1,np
! Subprogram not used            if (phis(i,j,ie)>mx) phis(i,j,ie)=mx
! Subprogram not used            if (phis(i,j,ie)<mn) phis(i,j,ie)=mn
! Subprogram not used         enddo
! Subprogram not used         enddo
! Subprogram not used         endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used         ! user specified minimum 
! Subprogram not used         do i=1,np
! Subprogram not used         do j=1,np
! Subprogram not used            if (phis(i,j,ie)<minf) phis(i,j,ie)=minf
! Subprogram not used         enddo
! Subprogram not used         enddo
! Subprogram not used 
! Subprogram not used         phis(:,:,ie)=phis(:,:,ie)*elem(ie)%spheremp(:,:)
! Subprogram not used         call edgeVpack(edge3p1,phis(:,:,ie),1,0,elem(ie)%desc)
! Subprogram not used      enddo
! Subprogram not used      call bndry_exchangeV(hybrid,edge3p1)
! Subprogram not used      do ie=nets,nete
! Subprogram not used         call edgeVunpack(edge3p1, phis(:,:,ie), 1, 0, elem(ie)%desc)
! Subprogram not used         phis(:,:,ie)=phis(:,:,ie)*elem(ie)%rspheremp(:,:)
! Subprogram not used      enddo
! Subprogram not used   enddo
! Subprogram not used   end subroutine smooth_phis

! Subprogram not used   subroutine overwrite_SEdensity(elem, fvm, dt_q, hybrid,nets,nete, np1)
! Subprogram not used     use fvm_reconstruction_mod, only: reconstruction
! Subprogram not used     use fvm_filter_mod, only: monotonic_gradient_cart, recons_val_cart
! Subprogram not used     use dimensions_mod, only : np, nlev, nc,nhe
! Subprogram not used     use hybrid_mod, only : hybrid_t
! Subprogram not used     use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack, edgevunpackmax, edgevunpackmin
! Subprogram not used     use bndry_mod, only : bndry_exchangev
! Subprogram not used     use element_mod, only : element_t
! Subprogram not used     use derivative_mod, only : derivative_t , laplace_sphere_wk
! Subprogram not used     use time_mod, only : TimeLevel_t
! Subprogram not used     use fvm_control_volume_mod, only : fvm_struct
! Subprogram not used     use spelt_mod, only : spelt_struct
! Subprogram not used     
! Subprogram not used 
! Subprogram not used     type (element_t) , intent(inout)        :: elem(:)
! Subprogram not used     
! Subprogram not used       type(fvm_struct), intent(inout) :: fvm(:)
! Subprogram not used     type (hybrid_t), intent(in)           :: hybrid  ! distributed parallel structure (shared)
! Subprogram not used 
! Subprogram not used     integer, intent(in)                     :: nets  ! starting thread element number (private)
! Subprogram not used     integer, intent(in)                     :: nete  ! ending thread element number   (private)
! Subprogram not used     integer, intent(in)                     :: np1
! Subprogram not used     integer :: ie, k
! Subprogram not used 
! Subprogram not used     real (kind=real_kind)             :: xp,yp, tmpval, dt_q
! Subprogram not used     integer                           :: i, j,ix, jy, starti,endi,tmpi
! Subprogram not used 
! Subprogram not used     real (kind=real_kind), dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe)      :: recons
! Subprogram not used     
! Subprogram not used     if ((nc .ne. 4) .or. (np .ne. 4)) then
! Subprogram not used       if(hybrid%masterthread) then 
! Subprogram not used         print *,"You are in OVERWRITE SE AIR DENSITY MODE"
! Subprogram not used         print *,"This only works for nc=4 and np=4"
! Subprogram not used         print *,"Write a new search algorithm or pay $10000!"
! Subprogram not used       endif
! Subprogram not used       stop
! Subprogram not used     endif
! Subprogram not used   end subroutine overwrite_SEdensity


end module prim_advance_mod

