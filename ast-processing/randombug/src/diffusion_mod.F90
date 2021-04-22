module diffusion_mod
  ! =======================
  use kinds,              only : real_kind
  ! =======================
  use dimensions_mod,     only : nlev, np, qsize
  ! =======================
  use derivative_mod,     only : gradient, vorticity, derivative_t, divergence
  ! =======================
  use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack, initedgebuffer
  ! =======================
  private
  save
  public :: diffusion_init, prim_diffusion, scalar_diffusion
  type(EdgeBuffer_t)  :: edgeS1, edgeS2
  type (EdgeBuffer_t) :: edge3
  type (EdgeBuffer_t) :: edge4

contains
  subroutine diffusion_init()

       call initEdgeBuffer(edgeS1,qsize*nlev)
       call initEdgeBuffer(edgeS2,2*qsize*nlev)
       call initEdgeBuffer(edge3, 3*nlev)
       call initEdgeBuffer(edge4, 4*nlev)

  end subroutine diffusion_init

! Subprogram not used   subroutine prim_diffusion(elem, nets,nete,np1,deriv,dt2,hybrid)
! Subprogram not used     use hybrid_mod, only : hybrid_t
! Subprogram not used     use physical_constants, only : rrearth
! Subprogram not used     use element_mod, only : element_t
! Subprogram not used     use bndry_mod, only : bndry_exchangev
! Subprogram not used     use control_mod, only : nu
! Subprogram not used     use perf_mod, only: t_startf, t_stopf, t_barrierf ! _EXTERNAL
! Subprogram not used     implicit none
! Subprogram not used     type(element_t), intent(inout), target :: elem(:)
! Subprogram not used     integer, intent(in)                :: nets,nete
! Subprogram not used     integer, intent(in)                :: np1
! Subprogram not used     type(derivative_t), intent(in)     :: deriv
! Subprogram not used     real(kind=real_kind),intent(in)    :: dt2
! Subprogram not used     type(hybrid_t),intent(in)       :: hybrid
! Subprogram not used 
! Subprogram not used     ! ======================
! Subprogram not used     ! Local Variables
! Subprogram not used     ! ======================
! Subprogram not used 
! Subprogram not used     real(kind=real_kind) :: temp(np,np)
! Subprogram not used 
! Subprogram not used     real(kind=real_kind) :: grad_T_np1(np,np,2,nlev,nets:nete)
! Subprogram not used     real(kind=real_kind) :: zeta_np1(np,np,nlev,nets:nete)
! Subprogram not used     real(kind=real_kind) :: div_np1(np,np,nlev,nets:nete)
! Subprogram not used     real(kind=real_kind) :: lap_v_np1(np,np,2,nlev,nets:nete)
! Subprogram not used     real(kind=real_kind) :: lap_T_np1(np,np,nlev,nets:nete)
! Subprogram not used 
! Subprogram not used     real(kind=real_kind) :: grad_zeta(np,np,2)
! Subprogram not used     real(kind=real_kind) :: grad_div(np,np,2)
! Subprogram not used     real(kind=real_kind) :: grad_tmp(np,np,2)
! Subprogram not used     real(kind=real_kind) :: vco(np,np,2)
! Subprogram not used     real(kind=real_kind) :: gp(np,np,2)
! Subprogram not used     real(kind=real_kind) :: zeta_tmp(np,np)
! Subprogram not used     real(kind=real_kind) :: div_tmp(np,np)
! Subprogram not used 
! Subprogram not used     real(kind=real_kind) :: v1,v2
! Subprogram not used 
! Subprogram not used     real (kind=real_kind), dimension(np,np)      :: rmetdetv   
! Subprogram not used 
! Subprogram not used     real(kind=real_kind), dimension(:,:), pointer :: mp
! Subprogram not used     real(kind=real_kind), dimension(:,:), pointer :: rmp
! Subprogram not used     real(kind=real_kind), dimension(:,:), pointer :: metdet
! Subprogram not used     real(kind=real_kind), dimension(:,:,:,:), pointer :: met
! Subprogram not used     real(kind=real_kind), dimension(:,:,:,:), pointer :: metinv
! Subprogram not used     real(kind=real_kind), dimension(:,:,:,:), pointer :: D
! Subprogram not used     real(kind=real_kind), dimension(:,:,:,:), pointer :: Dinv
! Subprogram not used 
! Subprogram not used     integer :: i,j,k,ie,kptr
! Subprogram not used 
! Subprogram not used     call t_barrierf('sync_prim_diffusion', hybrid%par%comm)
! Subprogram not used     call t_startf('prim_diffusion')
! Subprogram not used 
! Subprogram not used     do ie=nets,nete
! Subprogram not used 
! Subprogram not used        mp       => elem(ie)%mp
! Subprogram not used        rmp      => elem(ie)%rmp
! Subprogram not used        met      => elem(ie)%met
! Subprogram not used        metinv   => elem(ie)%metinv
! Subprogram not used        metdet   => elem(ie)%metdet
! Subprogram not used        rmetdetv(:,:)=1.0_real_kind/elem(ie)%metdet(:,:)
! Subprogram not used 
! Subprogram not used        do k=1,nlev   
! Subprogram not used 
! Subprogram not used           grad_tmp(:,:,:) = gradient(elem(ie)%state%T(:,:,k,np1),deriv)*rrearth
! Subprogram not used 
! Subprogram not used           do j=1,np
! Subprogram not used              do i=1,np
! Subprogram not used !                 grad_T_np1(i,j,1,k,ie) = mp(i,j)*(metinv(1,1,i,j)*grad_tmp(i,j,1)+metinv(1,2,i,j)*grad_tmp(i,j,2))
! Subprogram not used !                 grad_T_np1(i,j,2,k,ie) = mp(i,j)*(metinv(2,1,i,j)*grad_tmp(i,j,1)+metinv(2,2,i,j)*grad_tmp(i,j,2))
! Subprogram not used ! Map grad_tmp to lat-lon instead of rotating it
! Subprogram not used                  grad_T_np1(i,j,1,k,ie) = mp(i,j)*(elem(ie)%Dinv(1,1,i,j)*grad_tmp(i,j,1)+elem(ie)%Dinv(2,1,i,j)*grad_tmp(i,j,2))
! Subprogram not used                  grad_T_np1(i,j,2,k,ie) = mp(i,j)*(elem(ie)%Dinv(1,2,i,j)*grad_tmp(i,j,1)+elem(ie)%Dinv(2,2,i,j)*grad_tmp(i,j,2))
! Subprogram not used                 v1     = elem(ie)%state%v(i,j,1,k,np1)
! Subprogram not used                 v2     = elem(ie)%state%v(i,j,2,k,np1)
! Subprogram not used 
! Subprogram not used                 gp(i,j,1) = metdet(i,j)*(elem(ie)%Dinv(1,1,i,j)*v1 + elem(ie)%Dinv(1,2,i,j)*v2)
! Subprogram not used                 gp(i,j,2) = metdet(i,j)*(elem(ie)%Dinv(2,1,i,j)*v1 + elem(ie)%Dinv(2,2,i,j)*v2)
! Subprogram not used                 vco(i,j,1) = elem(ie)%D(1,1,i,j)*v1 + elem(ie)%D(2,1,i,j)*v2
! Subprogram not used                 vco(i,j,2) = elem(ie)%D(1,2,i,j)*v1 + elem(ie)%D(2,2,i,j)*v2
! Subprogram not used              end do
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used           div_tmp(:,:)  = divergence(gp,deriv)*rrearth
! Subprogram not used           zeta_tmp(:,:) = vorticity(vco,deriv)*rrearth
! Subprogram not used 
! Subprogram not used           do j=1,np
! Subprogram not used              do i=1,np
! Subprogram not used                 div_np1(i,j,k,ie)  = mp(i,j)*rmetdetv(i,j)*div_tmp(i,j)
! Subprogram not used                 zeta_np1(i,j,k,ie) = mp(i,j)*rmetdetv(i,j)*zeta_tmp(i,j)
! Subprogram not used              end do
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used        end do
! Subprogram not used 
! Subprogram not used        kptr=0
! Subprogram not used        call edgeVpack(edge4,grad_T_np1(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)
! Subprogram not used 
! Subprogram not used        kptr=2*nlev
! Subprogram not used        call edgeVpack(edge4,zeta_np1(1,1,1,ie),nlev,kptr,elem(ie)%desc)
! Subprogram not used 
! Subprogram not used        kptr=3*nlev
! Subprogram not used        call edgeVpack(edge4,div_np1(1,1,1,ie),nlev,kptr,elem(ie)%desc)
! Subprogram not used 
! Subprogram not used !       kptr=0
! Subprogram not used !       call edgerotate(edge4,2*nlev,kptr,elem(ie)%desc)
! Subprogram not used 
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     call bndry_exchangeV(hybrid,edge4)
! Subprogram not used 
! Subprogram not used     do ie=nets,nete
! Subprogram not used 
! Subprogram not used        mp       => elem(ie)%mp
! Subprogram not used        rmp      => elem(ie)%rmp
! Subprogram not used        metdet   => elem(ie)%metdet
! Subprogram not used        metinv   => elem(ie)%metinv
! Subprogram not used        rmetdetv(:,:)=1.0_real_kind/elem(ie)%metdet(:,:)
! Subprogram not used        D        => elem(ie)%D
! Subprogram not used        Dinv     => elem(ie)%Dinv
! Subprogram not used 
! Subprogram not used        kptr=0
! Subprogram not used        call edgeVunpack(edge4, grad_T_np1(1,1,1,1,ie), 2*nlev, kptr, elem(ie)%desc)
! Subprogram not used 
! Subprogram not used        kptr=2*nlev
! Subprogram not used        call edgeVunpack(edge4, zeta_np1(1,1,1,ie), nlev, kptr, elem(ie)%desc)
! Subprogram not used 
! Subprogram not used        kptr=3*nlev
! Subprogram not used        call edgeVunpack(edge4, div_np1(1,1,1,ie), nlev, kptr, elem(ie)%desc)
! Subprogram not used        ! ======================================
! Subprogram not used        ! compute Laplacian of T(n+1)
! Subprogram not used        ! ======================================
! Subprogram not used 
! Subprogram not used        do k=1,nlev
! Subprogram not used           do j=1,np
! Subprogram not used              do i=1,np
! Subprogram not used 
! Subprogram not used                 ! ======================================================
! Subprogram not used                 ! complete global assembly of gradT(n+1), zeta(n+1), div(n+1)
! Subprogram not used                 ! ======================================================
! Subprogram not used 
! Subprogram not used !                v1 = rmp(i,j)*grad_T_np1(i,j,1,k,ie)
! Subprogram not used !                v2 = rmp(i,j)*grad_T_np1(i,j,2,k,ie)
! Subprogram not used 
! Subprogram not used                 zeta_np1(i,j,k,ie) = rmp(i,j)*zeta_np1(i,j,k,ie)
! Subprogram not used                 div_np1(i,j,k,ie)  = rmp(i,j)*div_np1(i,j,k,ie)
! Subprogram not used                 elem(ie)%derived%zeta(i,j,k) = zeta_np1(i,j,k,ie)
! Subprogram not used 
! Subprogram not used                 ! ==========================================
! Subprogram not used                 ! Compute contravariant gradient(T(n+1))
! Subprogram not used                 ! ==========================================
! Subprogram not used 
! Subprogram not used !                grad_T_np1(i,j,1,k,ie) = metdet(i,j)*(metinv(1,1,i,j)*v1 + metinv(1,2,i,j)*v2)
! Subprogram not used !                grad_T_np1(i,j,2,k,ie) = metdet(i,j)*(metinv(2,1,i,j)*v1 + metinv(2,2,i,j)*v2)
! Subprogram not used                 grad_tmp(i,j,:) = grad_T_np1(i,j,:,k,ie)
! Subprogram not used                 grad_T_np1(i,j,1,k,ie) = metdet(i,j)*rmp(i,j)*&
! Subprogram not used        (elem(ie)%Dinv(1,1,i,j)*grad_tmp(i,j,1)+elem(ie)%Dinv(1,2,i,j)*grad_tmp(i,j,2))
! Subprogram not used                 grad_T_np1(i,j,2,k,ie) = metdet(i,j)*rmp(i,j)*&
! Subprogram not used        (elem(ie)%Dinv(2,1,i,j)*grad_tmp(i,j,1)+elem(ie)%Dinv(2,2,i,j)*grad_tmp(i,j,2))
! Subprogram not used                 !grad_T_np1(i,j,1,k,ie) = metdet(i,j)*rmp(i,j)*grad_T_np1(i,j,1,k,ie)
! Subprogram not used                 !grad_T_np1(i,j,2,k,ie) = metdet(i,j)*rmp(i,j)*grad_T_np1(i,j,2,k,ie)
! Subprogram not used 
! Subprogram not used              end do
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used           div_tmp(:,:) = divergence(grad_T_np1(:,:,:,k,ie),deriv)*rrearth
! Subprogram not used 
! Subprogram not used           do j=1,np
! Subprogram not used              do i=1,np
! Subprogram not used                 lap_T_np1(i,j,k,ie)   = mp(i,j)*rmetdetv(i,j)*div_tmp(i,j)
! Subprogram not used              end do
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used           ! ==========================================
! Subprogram not used           ! compute covariant Laplacian of v(n+1)  
! Subprogram not used           ! ==========================================
! Subprogram not used 
! Subprogram not used           grad_div(:,:,:)  = gradient(div_np1(:,:,k,ie),deriv)*rrearth
! Subprogram not used           grad_zeta(:,:,:) = gradient(zeta_np1(:,:,k,ie),deriv)*rrearth
! Subprogram not used 
! Subprogram not used           do j=1,np
! Subprogram not used              do i=1,np
! Subprogram not used 
! Subprogram not used                 v2 = Dinv(1,1,i,j)*grad_zeta(i,j,1) + Dinv(2,1,i,j)*grad_zeta(i,j,2)
! Subprogram not used                 v1 = Dinv(1,2,i,j)*grad_zeta(i,j,1) + Dinv(2,2,i,j)*grad_zeta(i,j,2)
! Subprogram not used 
! Subprogram not used                 v2 = -v2
! Subprogram not used                 grad_zeta(i,j,1) = v1
! Subprogram not used                 grad_zeta(i,j,2) = v2
! Subprogram not used                 v1 = grad_div(i,j,1)
! Subprogram not used                 v2 = grad_div(i,j,2)
! Subprogram not used 
! Subprogram not used                 grad_div(i,j,1) = Dinv(1,1,i,j)*v1 + Dinv(2,1,i,j)*v2
! Subprogram not used                 grad_div(i,j,2) = Dinv(1,2,i,j)*v1 + Dinv(2,2,i,j)*v2
! Subprogram not used                 lap_v_np1(i,j,1,k,ie) = mp(i,j)*(grad_div(i,j,1) - grad_zeta(i,j,1))
! Subprogram not used                 lap_v_np1(i,j,2,k,ie) = mp(i,j)*(grad_div(i,j,2) - grad_zeta(i,j,2))
! Subprogram not used 
! Subprogram not used              end do
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used        end do
! Subprogram not used 
! Subprogram not used        kptr=0
! Subprogram not used        call edgeVpack(edge3, lap_v_np1(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)
! Subprogram not used 
! Subprogram not used        kptr=2*nlev
! Subprogram not used        call edgeVpack(edge3, lap_T_np1(1,1,1,ie),nlev,kptr,elem(ie)%desc)
! Subprogram not used 
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     call bndry_exchangeV(hybrid,edge3)
! Subprogram not used 
! Subprogram not used     do ie=nets,nete
! Subprogram not used 
! Subprogram not used        rmp   => elem(ie)%rmp
! Subprogram not used        metdet   => elem(ie)%metdet
! Subprogram not used        rmetdetv(:,:)=1.0_real_kind/elem(ie)%metdet(:,:)
! Subprogram not used 
! Subprogram not used        kptr=0
! Subprogram not used        call edgeVunpack(edge3, lap_v_np1(1,1,1,1,ie), 2*nlev, kptr, elem(ie)%desc)
! Subprogram not used 
! Subprogram not used        kptr=2*nlev
! Subprogram not used        call edgeVunpack(edge3, lap_T_np1(1,1,1,ie), nlev, kptr, elem(ie)%desc)
! Subprogram not used 
! Subprogram not used        do k=1,nlev   
! Subprogram not used 
! Subprogram not used           do j=1,np
! Subprogram not used              do i=1,np
! Subprogram not used                 lap_T_np1(i,j,k,ie)   = rmp(i,j)*lap_T_np1(i,j,k,ie)
! Subprogram not used                 lap_v_np1(i,j,1,k,ie) = rmp(i,j)*lap_v_np1(i,j,1,k,ie)
! Subprogram not used                 lap_v_np1(i,j,2,k,ie) = rmp(i,j)*lap_v_np1(i,j,2,k,ie)
! Subprogram not used 
! Subprogram not used                 elem(ie)%state%T(i,j,k,np1)   = elem(ie)%state%T(i,j,k,np1)   + nu*dt2*lap_T_np1(i,j,k,ie)
! Subprogram not used                 elem(ie)%state%v(i,j,1,k,np1) = elem(ie)%state%v(i,j,1,k,np1) + nu*dt2*lap_v_np1(i,j,1,k,ie)
! Subprogram not used                 elem(ie)%state%v(i,j,2,k,np1) = elem(ie)%state%v(i,j,2,k,np1) + nu*dt2*lap_v_np1(i,j,2,k,ie)
! Subprogram not used              end do
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used        end do
! Subprogram not used 
! Subprogram not used     end do
! Subprogram not used     call t_stopf('prim_diffusion')
! Subprogram not used 
! Subprogram not used   end subroutine prim_diffusion

  !only called by prim_advect_scalars_lf (no longer supported) 
  ! so the Q part is not up to date for this
! Subprogram not used   subroutine scalar_diffusion(elem, nets,nete,np1,deriv,dt2,hybrid)
! Subprogram not used     use hybrid_mod, only : hybrid_t
! Subprogram not used     use physical_constants, only : rrearth
! Subprogram not used     use element_mod, only : element_t
! Subprogram not used     use bndry_mod, only : bndry_exchangev
! Subprogram not used     use control_mod, only : nu_s
! Subprogram not used     use perf_mod, only: t_startf, t_stopf, t_barrierf ! _EXTERNAL
! Subprogram not used     implicit none
! Subprogram not used     type(element_t), intent(inout), target :: elem(:)
! Subprogram not used     integer, intent(in)                :: nets,nete
! Subprogram not used     integer, intent(in)                :: np1
! Subprogram not used     type(derivative_t), intent(in)     :: deriv
! Subprogram not used     real(kind=real_kind),intent(in)    :: dt2
! Subprogram not used     type(hybrid_t),intent(in)       :: hybrid
! Subprogram not used 
! Subprogram not used     ! ======================
! Subprogram not used     ! Local Variables
! Subprogram not used     ! ======================
! Subprogram not used     real(kind=real_kind) :: temp(np,np)
! Subprogram not used 
! Subprogram not used     real(kind=real_kind) :: grad_Q_np1(np,np,2,nlev,qsize,nets:nete)
! Subprogram not used     real(kind=real_kind) :: lap_Q_np1(np,np,nlev,qsize,nets:nete)
! Subprogram not used 
! Subprogram not used     real(kind=real_kind) :: grad_tmp(np,np,2)
! Subprogram not used     real(kind=real_kind) :: div_tmp(np,np)
! Subprogram not used 
! Subprogram not used     real(kind=real_kind) :: v1,v2
! Subprogram not used 
! Subprogram not used     real (kind=real_kind), dimension(np,np)      :: rmetdetv   
! Subprogram not used 
! Subprogram not used     real(kind=real_kind), dimension(:,:), pointer :: mp
! Subprogram not used     real(kind=real_kind), dimension(:,:), pointer :: rmp
! Subprogram not used     real(kind=real_kind), dimension(:,:), pointer :: metdet
! Subprogram not used     real(kind=real_kind), dimension(:,:,:,:), pointer :: met
! Subprogram not used     real(kind=real_kind), dimension(:,:,:,:), pointer :: metinv
! Subprogram not used     real(kind=real_kind), dimension(:,:,:,:), pointer :: D
! Subprogram not used     real(kind=real_kind), dimension(:,:,:,:), pointer :: Dinv
! Subprogram not used 
! Subprogram not used     integer :: i,j,k,ie, q
! Subprogram not used 
! Subprogram not used     call t_barrierf('sync_scalar_diffusion', hybrid%par%comm)
! Subprogram not used     call t_startf('scalar_diffusion')
! Subprogram not used     do ie=nets,nete
! Subprogram not used        do q=1,qsize
! Subprogram not used 
! Subprogram not used           mp       => elem(ie)%mp
! Subprogram not used           metinv   => elem(ie)%metinv
! Subprogram not used           metdet   => elem(ie)%metdet
! Subprogram not used 
! Subprogram not used           do k=1,nlev   
! Subprogram not used              
! Subprogram not used              grad_tmp(:,:,:) = gradient(elem(ie)%state%Q(:,:,k,q),deriv)*rrearth
! Subprogram not used 
! Subprogram not used              do j=1,np
! Subprogram not used                 do i=1,np
! Subprogram not used                    !                grad_Q_np1(i,j,1,k,ie) = mp(i,j)*grad_tmp(i,j,1)
! Subprogram not used                    !                grad_Q_np1(i,j,2,k,ie) = mp(i,j)*grad_tmp(i,j,2)
! Subprogram not used 
! Subprogram not used                    v1 = mp(i,j)*grad_tmp(i,j,1)
! Subprogram not used                    v2 = mp(i,j)*grad_tmp(i,j,2)
! Subprogram not used 
! Subprogram not used                    grad_Q_np1(i,j,1,k,q,ie) = metdet(i,j)*(metinv(1,1,i,j)*v1 + metinv(1,2,i,j)*v2)
! Subprogram not used                    grad_Q_np1(i,j,2,k,q,ie) = metdet(i,j)*(metinv(2,1,i,j)*v1 + metinv(2,2,i,j)*v2)
! Subprogram not used 
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used        end do
! Subprogram not used        
! Subprogram not used        call edgeVpack(edgeS2,grad_Q_np1(:,:,:,:,:,ie),2*nlev*qsize,0,elem(ie)%desc)
! Subprogram not used 
! Subprogram not used        call edgerotate(edgeS2,2*nlev*qsize,0,elem(ie)%desc)
! Subprogram not used 
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     call bndry_exchangeV(hybrid,edgeS2)
! Subprogram not used 
! Subprogram not used     do ie=nets,nete
! Subprogram not used 
! Subprogram not used        mp       => elem(ie)%mp
! Subprogram not used        rmp      => elem(ie)%rmp
! Subprogram not used        metdet   => elem(ie)%metdet
! Subprogram not used        metinv   => elem(ie)%metinv
! Subprogram not used        rmetdetv(:,:)=1.0_real_kind/elem(ie)%metdet(:,:)
! Subprogram not used        D        => elem(ie)%D
! Subprogram not used        Dinv     => elem(ie)%Dinv
! Subprogram not used        
! Subprogram not used        call edgeVunpack(edgeS2, grad_Q_np1(:,:,:,:,:,ie), 2*nlev*qsize, 0, elem(ie)%desc)
! Subprogram not used 
! Subprogram not used        do q=1,qsize
! Subprogram not used 
! Subprogram not used           do k=1,nlev
! Subprogram not used              do j=1,np
! Subprogram not used                 do i=1,np
! Subprogram not used 
! Subprogram not used                    !                v1 = rmp(i,j)*grad_Q_np1(i,j,1,k,ie)
! Subprogram not used                    !                v2 = rmp(i,j)*grad_Q_np1(i,j,2,k,ie)
! Subprogram not used 
! Subprogram not used                    !                grad_Q_np1(i,j,1,k,ie) = metdet(i,j)*(metinv(1,1,i,j)*v1 + metinv(1,2,i,j)*v2)
! Subprogram not used                    !                grad_Q_np1(i,j,2,k,ie) = metdet(i,j)*(metinv(2,1,i,j)*v1 + metinv(2,2,i,j)*v2)
! Subprogram not used 
! Subprogram not used                    grad_Q_np1(i,j,1,k,q,ie) = rmp(i,j)*grad_Q_np1(i,j,1,k,q,ie)
! Subprogram not used                    grad_Q_np1(i,j,2,k,q,ie) = rmp(i,j)*grad_Q_np1(i,j,2,k,q,ie)
! Subprogram not used 
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used              div_tmp(:,:) = divergence(grad_Q_np1(:,:,:,k,q,ie),deriv)*rrearth
! Subprogram not used 
! Subprogram not used              do j=1,np
! Subprogram not used                 do i=1,np
! Subprogram not used                    lap_Q_np1(i,j,k,q,ie)   = mp(i,j)*rmetdetv(i,j)*div_tmp(i,j)
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used        end do
! Subprogram not used 
! Subprogram not used        call edgeVpack(edgeS1, lap_Q_np1(1,1,1,1,ie),nlev*qsize, 0,elem(ie)%desc)
! Subprogram not used 
! Subprogram not used     end do
! Subprogram not used     call bndry_exchangeV(hybrid,edgeS1)
! Subprogram not used     do ie=nets,nete
! Subprogram not used           
! Subprogram not used        rmp   => elem(ie)%rmp
! Subprogram not used        metdet   => elem(ie)%metdet
! Subprogram not used        rmetdetv(:,:)=1.0_real_kind/elem(ie)%metdet(:,:)
! Subprogram not used 
! Subprogram not used        
! Subprogram not used        call edgeVunpack(edgeS1, lap_Q_np1(1,1,1,1,ie), nlev*qsize, 0, elem(ie)%desc)
! Subprogram not used 
! Subprogram not used        do q=1,qsize
! Subprogram not used           do k=1,nlev   
! Subprogram not used 
! Subprogram not used              do j=1,np
! Subprogram not used                 do i=1,np
! Subprogram not used                    lap_Q_np1(i,j,k,q,ie)   = rmp(i,j)*lap_Q_np1(i,j,k,q,ie)
! Subprogram not used                    elem(ie)%state%Q(i,j,k,q)   = elem(ie)%state%Q(i,j,k,q)   + nu_s*dt2*lap_Q_np1(i,j,k,q,ie)
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used     call t_stopf('scalar_diffusion')
! Subprogram not used   end subroutine scalar_diffusion

end module diffusion_mod
