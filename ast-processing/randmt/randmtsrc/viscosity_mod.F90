module viscosity_mod
!
!  This module should be renamed "global_deriv_mod.F90"
! 
!  It is a collection of derivative operators that must be applied to the field 
!  over the sphere (as opposed to derivative operators that can be applied element 
!  by element)
!
!
use kinds, only : real_kind, iulog
use dimensions_mod, only : np, nlev,qsize
use hybrid_mod, only : hybrid_t
use element_mod, only : element_t
use derivative_mod, only : derivative_t, laplace_sphere_wk, vlaplace_sphere_wk, vorticity_sphere, derivinit, divergence_sphere
use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack, edgevunpackmin, &
    edgevunpackmax, initEdgeBuffer, FreeEdgeBuffer
use bndry_mod, only : bndry_exchangev
use control_mod, only : hypervis_scaling, nu, nu_div

implicit none
save

public :: biharmonic_wk
public :: biharmonic_wk_scalar
public :: biharmonic_wk_scalar_minmax
public :: compute_zeta_C0
public :: compute_div_C0
public :: compute_zeta_C0_2d
public :: compute_div_C0_2d
public :: test_ibyp

interface compute_zeta_C0_2d
    module procedure compute_zeta_C0_2d_sphere
    module procedure compute_zeta_C0_2d_contra
end interface

interface compute_div_C0_2d
    module procedure compute_div_C0_2d_sphere
    module procedure compute_div_C0_2d_contra
end interface


type (EdgeBuffer_t)          :: edge1

contains

! Subprogram not used subroutine biharmonic_wk(elem,pstens,ptens,vtens,deriv,edge3,hybrid,nt,nets,nete)
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used ! compute weak biharmonic operator
! Subprogram not used !    input:  h,v (stored in elem()%, in lat-lon coordinates
! Subprogram not used !    output: ptens,vtens  overwritten with weak biharmonic of h,v (output in lat-lon coordinates)
! Subprogram not used !
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used type (hybrid_t)      , intent(in) :: hybrid
! Subprogram not used type (element_t)     , intent(inout), target :: elem(:)
! Subprogram not used integer :: nt,nets,nete
! Subprogram not used real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)  :: vtens
! Subprogram not used real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: ptens
! Subprogram not used type (EdgeBuffer_t)  , intent(inout) :: edge3
! Subprogram not used type (derivative_t)  , intent(in) :: deriv
! Subprogram not used real (kind=real_kind), dimension(np,np,nets:nete) :: pstens
! Subprogram not used 
! Subprogram not used ! local
! Subprogram not used integer :: k,kptr,i,j,ie,ic
! Subprogram not used real (kind=real_kind), dimension(:,:), pointer :: rspheremv
! Subprogram not used real (kind=real_kind), dimension(np,np) :: lap_ps
! Subprogram not used real (kind=real_kind), dimension(np,np,nlev) :: T
! Subprogram not used real (kind=real_kind), dimension(np,np,2) :: v
! Subprogram not used real (kind=real_kind) ::  nu_ratio1,nu_ratio2
! Subprogram not used logical var_coef1
! Subprogram not used 
! Subprogram not used    !if tensor hyperviscosity with tensor V is used, then biharmonic operator is (\grad\cdot V\grad) (\grad \cdot \grad) 
! Subprogram not used    !so tensor is only used on second call to laplace_sphere_wk
! Subprogram not used    var_coef1 = .true.
! Subprogram not used    if(hypervis_scaling > 0)  var_coef1= .false.
! Subprogram not used 
! Subprogram not used    ! note: there is a scaling bug in the treatment of nu_div
! Subprogram not used    ! nu_ratio is applied twice, once in each laplace operator
! Subprogram not used    ! so in reality:   nu_div_actual = (nu_div/nu)**2 nu
! Subprogram not used    ! We should fix this, but it requires adjusting all 1 defaults
! Subprogram not used    nu_ratio1=1
! Subprogram not used    nu_ratio2=1
! Subprogram not used    if (nu_div/=nu) then
! Subprogram not used       if(hypervis_scaling /= 0) then
! Subprogram not used          ! we have a problem with the tensor in that we cant seperate
! Subprogram not used          ! div and curl components.  So we do, with tensor V:
! Subprogram not used          ! nu * (del V del ) * ( nu_ratio * grad(div) - curl(curl))
! Subprogram not used          nu_ratio1=(nu_div/nu)**2   ! preserve buggy scaling
! Subprogram not used          nu_ratio2=1
! Subprogram not used       else
! Subprogram not used          nu_ratio1=nu_div/nu
! Subprogram not used          nu_ratio2=nu_div/nu
! Subprogram not used       endif
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    do ie=nets,nete
! Subprogram not used       
! Subprogram not used       ! should filter lnps + PHI_s/RT?
! Subprogram not used       pstens(:,:,ie)=laplace_sphere_wk(elem(ie)%state%ps_v(:,:,nt),deriv,elem(ie),var_coef=var_coef1)
! Subprogram not used       
! Subprogram not used       do k=1,nlev
! Subprogram not used          do j=1,np
! Subprogram not used             do i=1,np
! Subprogram not used                T(i,j,k)=elem(ie)%state%T(i,j,k,nt) 
! Subprogram not used             enddo
! Subprogram not used          enddo
! Subprogram not used         
! Subprogram not used          ptens(:,:,k,ie)=laplace_sphere_wk(T(:,:,k),deriv,elem(ie),var_coef=var_coef1)
! Subprogram not used          vtens(:,:,:,k,ie)=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,&
! Subprogram not used               elem(ie),var_coef=var_coef1,nu_ratio=nu_ratio1)
! Subprogram not used 
! Subprogram not used       enddo
! Subprogram not used       kptr=0
! Subprogram not used       call edgeVpack(edge3, ptens(1,1,1,ie),nlev,kptr,elem(ie)%desc)
! Subprogram not used       kptr=nlev
! Subprogram not used       call edgeVpack(edge3, vtens(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)
! Subprogram not used 
! Subprogram not used       kptr=3*nlev
! Subprogram not used       call edgeVpack(edge3, pstens(1,1,ie),1,kptr,elem(ie)%desc)
! Subprogram not used    enddo
! Subprogram not used    
! Subprogram not used    call bndry_exchangeV(hybrid,edge3)
! Subprogram not used    
! Subprogram not used    do ie=nets,nete
! Subprogram not used       rspheremv     => elem(ie)%rspheremp(:,:)
! Subprogram not used       
! Subprogram not used       kptr=0
! Subprogram not used       call edgeVunpack(edge3, ptens(1,1,1,ie), nlev, kptr, elem(ie)%desc)
! Subprogram not used       kptr=nlev
! Subprogram not used       call edgeVunpack(edge3, vtens(1,1,1,1,ie), 2*nlev, kptr, elem(ie)%desc)
! Subprogram not used       
! Subprogram not used       ! apply inverse mass matrix, then apply laplace again
! Subprogram not used       do k=1,nlev
! Subprogram not used          do j=1,np
! Subprogram not used             do i=1,np
! Subprogram not used                T(i,j,k)=rspheremv(i,j)*ptens(i,j,k,ie)
! Subprogram not used                v(i,j,1)=rspheremv(i,j)*vtens(i,j,1,k,ie)
! Subprogram not used                v(i,j,2)=rspheremv(i,j)*vtens(i,j,2,k,ie)
! Subprogram not used             enddo
! Subprogram not used          enddo
! Subprogram not used          ptens(:,:,k,ie)=laplace_sphere_wk(T(:,:,k),deriv,elem(ie),var_coef=.true.)
! Subprogram not used          vtens(:,:,:,k,ie)=vlaplace_sphere_wk(v(:,:,:),deriv,elem(ie),var_coef=.true.,&
! Subprogram not used               nu_ratio=nu_ratio2)
! Subprogram not used       enddo
! Subprogram not used          
! Subprogram not used       kptr=3*nlev
! Subprogram not used       call edgeVunpack(edge3, pstens(1,1,ie), 1, kptr, elem(ie)%desc)
! Subprogram not used       ! apply inverse mass matrix, then apply laplace again
! Subprogram not used       lap_ps(:,:)=rspheremv(:,:)*pstens(:,:,ie)
! Subprogram not used       pstens(:,:,ie)=laplace_sphere_wk(lap_ps,deriv,elem(ie),var_coef=.true.)
! Subprogram not used 
! Subprogram not used    enddo
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used end subroutine


subroutine biharmonic_wk_dp3d(elem,dptens,ptens,vtens,deriv,edge3,hybrid,nt,nets,nete)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute weak biharmonic operator
!    input:  h,v (stored in elem()%, in lat-lon coordinates
!    output: ptens,vtens  overwritten with weak biharmonic of h,v (output in lat-lon coordinates)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(inout), target :: elem(:)
integer :: nt,nets,nete
real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)  :: vtens
real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: ptens,dptens
type (EdgeBuffer_t)  , intent(inout) :: edge3
type (derivative_t)  , intent(in) :: deriv

! local
integer :: k,kptr,ie
real (kind=real_kind), dimension(:,:), pointer :: rspheremv
real (kind=real_kind), dimension(np,np) :: tmp
real (kind=real_kind), dimension(np,np,2) :: v
real (kind=real_kind) :: nu_ratio1, nu_ratio2
logical var_coef1

   !if tensor hyperviscosity with tensor V is used, then biharmonic operator is (\grad\cdot V\grad) (\grad \cdot \grad) 
   !so tensor is only used on second call to laplace_sphere_wk
   var_coef1 = .true.
   if(hypervis_scaling > 0)    var_coef1 = .false.

   ! note: there is a scaling bug in the treatment of nu_div
   ! nu_ratio is applied twice, once in each laplace operator
   ! so in reality:   nu_div_actual = (nu_div/nu)**2 nu
   ! We should fix this, but it requires adjusting all 1 defaults
   nu_ratio1=1
   nu_ratio2=1
   if (nu_div/=nu) then
      if(hypervis_scaling /= 0) then
         ! we have a problem with the tensor in that we cant seperate
         ! div and curl components.  So we do, with tensor V:
         ! nu * (del V del ) * ( nu_ratio * grad(div) - curl(curl))
         nu_ratio1=(nu_div/nu)**2   ! preserve buggy scaling
         nu_ratio2=1
      else
         nu_ratio1=nu_div/nu
         nu_ratio2=nu_div/nu
      endif
   endif



   do ie=nets,nete

      do k=1,nlev
         tmp=elem(ie)%state%T(:,:,k,nt) 
         ptens(:,:,k,ie)=laplace_sphere_wk(tmp,deriv,elem(ie),var_coef=var_coef1)
         tmp=elem(ie)%state%dp3d(:,:,k,nt) 
         dptens(:,:,k,ie)=laplace_sphere_wk(tmp,deriv,elem(ie),var_coef=var_coef1)
         vtens(:,:,:,k,ie)=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),&
              var_coef=var_coef1,nu_ratio=nu_ratio1)
      enddo
      kptr=0
      call edgeVpack(edge3, ptens(1,1,1,ie),nlev,kptr,elem(ie)%desc)
      kptr=nlev
      call edgeVpack(edge3, vtens(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)
      kptr=3*nlev
      call edgeVpack(edge3, dptens(1,1,1,ie),nlev,kptr,elem(ie)%desc)

   enddo
   
   call bndry_exchangeV(hybrid,edge3)
   
   do ie=nets,nete
      rspheremv     => elem(ie)%rspheremp(:,:)
      
      kptr=0
      call edgeVunpack(edge3, ptens(1,1,1,ie), nlev, kptr, elem(ie)%desc)
      kptr=nlev
      call edgeVunpack(edge3, vtens(1,1,1,1,ie), 2*nlev, kptr, elem(ie)%desc)
      kptr=3*nlev
      call edgeVunpack(edge3, dptens(1,1,1,ie), nlev, kptr, elem(ie)%desc)
      
      ! apply inverse mass matrix, then apply laplace again
      do k=1,nlev
         tmp(:,:)=rspheremv(:,:)*ptens(:,:,k,ie)
         ptens(:,:,k,ie)=laplace_sphere_wk(tmp,deriv,elem(ie),var_coef=.true.)
         tmp(:,:)=rspheremv(:,:)*dptens(:,:,k,ie)
         dptens(:,:,k,ie)=laplace_sphere_wk(tmp,deriv,elem(ie),var_coef=.true.)

         v(:,:,1)=rspheremv(:,:)*vtens(:,:,1,k,ie)
         v(:,:,2)=rspheremv(:,:)*vtens(:,:,2,k,ie)
         vtens(:,:,:,k,ie)=vlaplace_sphere_wk(v(:,:,:),deriv,elem(ie),&
              var_coef=.true.,nu_ratio=nu_ratio2)

      enddo
   enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine


! Subprogram not used subroutine biharmonic_wk_scalar(elem,qtens,deriv,edgeq,hybrid,nets,nete)
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used ! compute weak biharmonic operator
! Subprogram not used !    input:  qtens = Q
! Subprogram not used !    output: qtens = weak biharmonic of Q
! Subprogram not used !
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used type (hybrid_t)      , intent(in) :: hybrid
! Subprogram not used type (element_t)     , intent(inout), target :: elem(:)
! Subprogram not used integer :: nets,nete
! Subprogram not used real (kind=real_kind), dimension(np,np,nlev,qsize,nets:nete) :: qtens
! Subprogram not used type (EdgeBuffer_t)  , intent(inout) :: edgeq
! Subprogram not used type (derivative_t)  , intent(in) :: deriv
! Subprogram not used 
! Subprogram not used ! local
! Subprogram not used integer :: k,kptr,i,j,ie,ic,q
! Subprogram not used real (kind=real_kind), dimension(np,np) :: lap_p
! Subprogram not used logical var_coef1
! Subprogram not used 
! Subprogram not used    !if tensor hyperviscosity with tensor V is used, then biharmonic operator is (\grad\cdot V\grad) (\grad \cdot \grad) 
! Subprogram not used    !so tensor is only used on second call to laplace_sphere_wk
! Subprogram not used    var_coef1 = .true.
! Subprogram not used    if(hypervis_scaling > 0)    var_coef1 = .false.
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    do ie=nets,nete
! Subprogram not used       do q=1,qsize      
! Subprogram not used          do k=1,nlev    !  Potential loop inversion (AAM)
! Subprogram not used             lap_p(:,:)=qtens(:,:,k,q,ie)
! Subprogram not used ! Original use of qtens on left and right hand sides caused OpenMP errors (AAM)
! Subprogram not used            qtens(:,:,k,q,ie)=laplace_sphere_wk(lap_p,deriv,elem(ie),var_coef=var_coef1)
! Subprogram not used          enddo
! Subprogram not used       enddo
! Subprogram not used       call edgeVpack(edgeq, qtens(:,:,:,:,ie),qsize*nlev,0,elem(ie)%desc)
! Subprogram not used    enddo
! Subprogram not used 
! Subprogram not used    call bndry_exchangeV(hybrid,edgeq)
! Subprogram not used    
! Subprogram not used    do ie=nets,nete
! Subprogram not used       call edgeVunpack(edgeq, qtens(:,:,:,:,ie),qsize*nlev,0,elem(ie)%desc)
! Subprogram not used 
! Subprogram not used       ! apply inverse mass matrix, then apply laplace again
! Subprogram not used       do q=1,qsize      
! Subprogram not used       do k=1,nlev    !  Potential loop inversion (AAM)
! Subprogram not used          lap_p(:,:)=elem(ie)%rspheremp(:,:)*qtens(:,:,k,q,ie)
! Subprogram not used          qtens(:,:,k,q,ie)=laplace_sphere_wk(lap_p,deriv,elem(ie),var_coef=.true.)
! Subprogram not used       enddo
! Subprogram not used       enddo
! Subprogram not used    enddo
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used end subroutine

subroutine biharmonic_wk_scalar_minmax(elem,qtens,deriv,edgeq,hybrid,nets,nete,emin,emax)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute weak biharmonic operator
!    input:  qtens = Q
!    output: qtens = weak biharmonic of Q and Q element min/max
!
!    note: emin/emax must be initialized with Q element min/max.  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(inout), target :: elem(:)
integer :: nets,nete
real (kind=real_kind), dimension(np,np,nlev,qsize,nets:nete) :: qtens
type (EdgeBuffer_t)  , intent(inout) :: edgeq
type (derivative_t)  , intent(in) :: deriv
real (kind=real_kind), intent(out), dimension(nlev,qsize,nets:nete) :: emin,emax

! local
integer :: k,kptr,i,j,ie,ic,q
real (kind=real_kind), dimension(np,np) :: lap_p
real (kind=real_kind) :: Qmin(np,np,nlev,qsize)
real (kind=real_kind) :: Qmax(np,np,nlev,qsize)
logical var_coef1

   !if tensor hyperviscosity with tensor V is used, then biharmonic operator is (\grad\cdot V\grad) (\grad \cdot \grad) 
   !so tensor is only used on second call to laplace_sphere_wk
   var_coef1 = .true.
   if(hypervis_scaling > 0)    var_coef1 = .false.



   do ie=nets,nete
      do q=1,qsize      
      do k=1,nlev    !  Potential loop inversion (AAM)
         Qmin(:,:,k,q)=emin(k,q,ie)  ! need to set all values in element for
         Qmax(:,:,k,q)=emax(k,q,ie)  ! edgeVpack routine below
         lap_p(:,:) = qtens(:,:,k,q,ie)
! Original use of qtens on left and right hand sides caused OpenMP errors (AAM)
         qtens(:,:,k,q,ie)=laplace_sphere_wk(lap_p,deriv,elem(ie),var_coef=var_coef1)
      enddo
      enddo
      call edgeVpack(edgeq, qtens(:,:,:,:,ie),qsize*nlev,0,elem(ie)%desc)
      call edgeVpack(edgeq,Qmin,nlev*qsize,nlev*qsize,elem(ie)%desc)
      call edgeVpack(edgeq,Qmax,nlev*qsize,2*nlev*qsize,elem(ie)%desc)
   enddo
   
   call bndry_exchangeV(hybrid,edgeq)
   
   do ie=nets,nete
      do q=1,qsize      
      do k=1,nlev
         Qmin(:,:,k,q)=emin(k,q,ie)  ! restore element data.  we could avoid
         Qmax(:,:,k,q)=emax(k,q,ie)  ! this by adding a "ie" index to Qmin/max
      enddo
      enddo
! WARNING - edgeVunpackMin/Max take second argument as input/ouput
      call edgeVunpack(edgeq, qtens(:,:,:,:,ie),qsize*nlev,0,elem(ie)%desc)
      call edgeVunpackMin(edgeq, Qmin,qsize*nlev,qsize*nlev,elem(ie)%desc)
      call edgeVunpackMax(edgeq, Qmax,qsize*nlev,2*qsize*nlev,elem(ie)%desc)

      ! apply inverse mass matrix, then apply laplace again
      do q=1,qsize      
      do k=1,nlev    !  Potential loop inversion (AAM)
         lap_p(:,:)=elem(ie)%rspheremp(:,:)*qtens(:,:,k,q,ie)
         qtens(:,:,k,q,ie)=laplace_sphere_wk(lap_p,deriv,elem(ie),var_coef=.true.)
         ! note: only need to consider the corners, since the data we packed was
         ! constant within each element
         emin(k,q,ie)=min(qmin(1,1,k,q),qmin(1,np,k,q),qmin(np,1,k,q),qmin(np,np,k,q))
         emin(k,q,ie)=max(emin(k,q,ie),0d0)
         emax(k,q,ie)=max(qmax(1,1,k,q),qmax(1,np,k,q),qmax(np,1,k,q),qmax(np,np,k,q))
      enddo
      enddo
   enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine






! Subprogram not used subroutine make_C0_2d(zeta,elem,hybrid,nets,nete)
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used ! apply DSS (aka assembly procedure) to zeta.  
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used 
! Subprogram not used type (hybrid_t)      , intent(in) :: hybrid
! Subprogram not used type (element_t)     , intent(in), target :: elem(:)
! Subprogram not used integer :: nets,nete
! Subprogram not used real (kind=real_kind), dimension(np,np,nets:nete) :: zeta
! Subprogram not used 
! Subprogram not used ! local
! Subprogram not used integer :: k,i,j,ie,ic,kptr
! Subprogram not used 
! Subprogram not used 
! Subprogram not used call initEdgeBuffer(edge1,1)
! Subprogram not used 
! Subprogram not used do ie=nets,nete
! Subprogram not used    zeta(:,:,ie)=zeta(:,:,ie)*elem(ie)%spheremp(:,:)
! Subprogram not used    kptr=0
! Subprogram not used    call edgeVpack(edge1, zeta(1,1,ie),1,kptr,elem(ie)%desc)
! Subprogram not used enddo
! Subprogram not used call bndry_exchangeV(hybrid,edge1)
! Subprogram not used do ie=nets,nete
! Subprogram not used    kptr=0
! Subprogram not used    call edgeVunpack(edge1, zeta(1,1,ie),1,kptr,elem(ie)%desc)
! Subprogram not used    zeta(:,:,ie)=zeta(:,:,ie)*elem(ie)%rspheremp(:,:)
! Subprogram not used enddo
! Subprogram not used 
! Subprogram not used call FreeEdgeBuffer(edge1) 
! Subprogram not used end subroutine


! Subprogram not used subroutine make_C0(zeta,elem,hybrid,nets,nete)
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used ! apply DSS (aka assembly procedure) to zeta.  
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used 
! Subprogram not used type (hybrid_t)      , intent(in) :: hybrid
! Subprogram not used type (element_t)     , intent(in), target :: elem(:)
! Subprogram not used integer :: nets,nete
! Subprogram not used real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: zeta
! Subprogram not used 
! Subprogram not used ! local
! Subprogram not used integer :: k,i,j,ie,ic,kptr
! Subprogram not used 
! Subprogram not used 
! Subprogram not used call initEdgeBuffer(edge1,nlev)
! Subprogram not used 
! Subprogram not used do ie=nets,nete
! Subprogram not used    do k=1,nlev
! Subprogram not used       zeta(:,:,k,ie)=zeta(:,:,k,ie)*elem(ie)%spheremp(:,:)
! Subprogram not used    enddo
! Subprogram not used    kptr=0
! Subprogram not used    call edgeVpack(edge1, zeta(1,1,1,ie),nlev,kptr,elem(ie)%desc)
! Subprogram not used enddo
! Subprogram not used call bndry_exchangeV(hybrid,edge1)
! Subprogram not used do ie=nets,nete
! Subprogram not used    kptr=0
! Subprogram not used    call edgeVunpack(edge1, zeta(1,1,1,ie),nlev,kptr,elem(ie)%desc)
! Subprogram not used    do k=1,nlev
! Subprogram not used       zeta(:,:,k,ie)=zeta(:,:,k,ie)*elem(ie)%rspheremp(:,:)
! Subprogram not used    enddo
! Subprogram not used enddo
! Subprogram not used 
! Subprogram not used call FreeEdgeBuffer(edge1) 
! Subprogram not used end subroutine


! Subprogram not used subroutine make_C0_vector(v,elem,hybrid,nets,nete)
! Subprogram not used type (hybrid_t)      , intent(in) :: hybrid
! Subprogram not used type (element_t)     , intent(in), target :: elem(:)
! Subprogram not used integer :: nets,nete
! Subprogram not used real (kind=real_kind), dimension(np,np,2,nlev,nets:nete) :: v
! Subprogram not used 
! Subprogram not used ! local
! Subprogram not used integer :: k,i,j,ie,ic,kptr
! Subprogram not used type (EdgeBuffer_t)          :: edge2
! Subprogram not used 
! Subprogram not used 
! Subprogram not used call initEdgeBuffer(edge2,2*nlev)
! Subprogram not used 
! Subprogram not used do ie=nets,nete
! Subprogram not used    do k=1,nlev
! Subprogram not used       v(:,:,1,k,ie)=v(:,:,1,k,ie)*elem(ie)%spheremp(:,:)
! Subprogram not used       v(:,:,2,k,ie)=v(:,:,2,k,ie)*elem(ie)%spheremp(:,:)
! Subprogram not used    enddo
! Subprogram not used    kptr=0
! Subprogram not used    call edgeVpack(edge2, v(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)
! Subprogram not used enddo
! Subprogram not used call bndry_exchangeV(hybrid,edge2)
! Subprogram not used do ie=nets,nete
! Subprogram not used    kptr=0
! Subprogram not used    call edgeVunpack(edge2, v(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)
! Subprogram not used    do k=1,nlev
! Subprogram not used       v(:,:,1,k,ie)=v(:,:,1,k,ie)*elem(ie)%rspheremp(:,:)
! Subprogram not used       v(:,:,2,k,ie)=v(:,:,2,k,ie)*elem(ie)%rspheremp(:,:)
! Subprogram not used    enddo
! Subprogram not used enddo
! Subprogram not used 
! Subprogram not used call FreeEdgeBuffer(edge2) 
! Subprogram not used end subroutine





! Subprogram not used subroutine compute_zeta_C0_2d_sphere(zeta,elem,hybrid,nets,nete,nt,k)
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used ! compute C0 vorticity.  That is, solve:  
! Subprogram not used !     < PHI, zeta > = <PHI, curl(elem%state%v >
! Subprogram not used !
! Subprogram not used !    input:  v (stored in elem()%, in lat-lon coordinates)
! Subprogram not used !    output: zeta(:,:,:,:)   
! Subprogram not used !
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used 
! Subprogram not used type (hybrid_t)      , intent(in) :: hybrid
! Subprogram not used type (element_t)     , intent(in), target :: elem(:)
! Subprogram not used integer :: nt,nets,nete,k
! Subprogram not used real (kind=real_kind), dimension(np,np,nets:nete) :: zeta
! Subprogram not used 
! Subprogram not used ! local
! Subprogram not used integer :: i,j,ie,ic
! Subprogram not used type (derivative_t)          :: deriv
! Subprogram not used 
! Subprogram not used call derivinit(deriv)
! Subprogram not used 
! Subprogram not used do ie=nets,nete
! Subprogram not used    !    zeta(:,:,k,ie)=elem(ie)%state%zeta(:,:,k)
! Subprogram not used    zeta(:,:,ie)=vorticity_sphere(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie))
! Subprogram not used enddo
! Subprogram not used 
! Subprogram not used call make_C0_2d(zeta,elem,hybrid,nets,nete)
! Subprogram not used 
! Subprogram not used end subroutine

! Subprogram not used subroutine compute_zeta_C0_2d_contra(zeta,elem,hybrid,nets,nete,nt)
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used ! compute C0 vorticity.  That is, solve:  
! Subprogram not used !     < PHI, zeta > = <PHI, curl(elem%state%v >
! Subprogram not used !
! Subprogram not used !    input:  v (stored in elem()%, in contra-variant coordinates)
! Subprogram not used !    output: zeta(:,:,:,:)   
! Subprogram not used !
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used 
! Subprogram not used type (hybrid_t)      , intent(in) :: hybrid
! Subprogram not used type (element_t)     , intent(in), target :: elem(:)
! Subprogram not used integer :: nt,nets,nete
! Subprogram not used real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: zeta
! Subprogram not used real (kind=real_kind), dimension(np,np,2) :: ulatlon
! Subprogram not used real (kind=real_kind), dimension(np,np) :: v1,v2
! Subprogram not used 
! Subprogram not used ! local
! Subprogram not used integer :: k,ie
! Subprogram not used type (derivative_t)          :: deriv
! Subprogram not used 
! Subprogram not used call derivinit(deriv)
! Subprogram not used 
! Subprogram not used do k=1,nlev
! Subprogram not used do ie=nets,nete
! Subprogram not used     v1 = elem(ie)%state%v(:,:,1,k,nt)
! Subprogram not used     v2 = elem(ie)%state%v(:,:,2,k,nt)
! Subprogram not used     ulatlon(:,:,1) = elem(ie)%D(1,1,:,:)*v1 + elem(ie)%D(1,2,:,:)*v2
! Subprogram not used     ulatlon(:,:,2) = elem(ie)%D(2,1,:,:)*v1 + elem(ie)%D(2,2,:,:)*v2
! Subprogram not used    !    zeta(:,:,k,ie)=elem(ie)%state%zeta(:,:,k)
! Subprogram not used    zeta(:,:,k,ie)=vorticity_sphere(ulatlon,deriv,elem(ie))
! Subprogram not used enddo
! Subprogram not used enddo
! Subprogram not used 
! Subprogram not used call make_C0(zeta,elem,hybrid,nets,nete)
! Subprogram not used 
! Subprogram not used end subroutine


! Subprogram not used subroutine compute_div_C0_2d_sphere(zeta,elem,hybrid,nets,nete,nt,k)
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used ! compute C0 divergence. That is, solve:  
! Subprogram not used !     < PHI, zeta > = <PHI, div(elem%state%v >
! Subprogram not used !
! Subprogram not used !    input:  v (stored in elem()%, in lat-lon coordinates)
! Subprogram not used !    output: zeta(:,:,:,:)   
! Subprogram not used !
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used 
! Subprogram not used type (hybrid_t)      , intent(in) :: hybrid
! Subprogram not used type (element_t)     , intent(in), target :: elem(:)
! Subprogram not used integer :: nt,nets,nete,k
! Subprogram not used real (kind=real_kind), dimension(np,np,nets:nete) :: zeta
! Subprogram not used 
! Subprogram not used ! local
! Subprogram not used integer :: i,j,ie,ic
! Subprogram not used type (derivative_t)          :: deriv
! Subprogram not used 
! Subprogram not used call derivinit(deriv)
! Subprogram not used 
! Subprogram not used do ie=nets,nete
! Subprogram not used    !    zeta(:,:,k,ie)=elem(ie)%state%zeta(:,:,k)
! Subprogram not used    zeta(:,:,ie)=divergence_sphere(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie))
! Subprogram not used enddo
! Subprogram not used 
! Subprogram not used call make_C0_2d(zeta,elem,hybrid,nets,nete)
! Subprogram not used 
! Subprogram not used end subroutine

! Subprogram not used subroutine compute_div_C0_2d_contra(zeta,elem,hybrid,nets,nete,nt)
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used ! compute C0 divergence. That is, solve:  
! Subprogram not used !     < PHI, zeta > = <PHI, div(elem%state%v >
! Subprogram not used !
! Subprogram not used !    input:  v (stored in elem()%, in contra-variant coordinates)
! Subprogram not used !    output: zeta(:,:,:,:)   
! Subprogram not used !
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used 
! Subprogram not used type (hybrid_t)      , intent(in) :: hybrid
! Subprogram not used type (element_t)     , intent(in), target :: elem(:)
! Subprogram not used integer :: nt,nets,nete
! Subprogram not used real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: zeta
! Subprogram not used real (kind=real_kind), dimension(np,np,2) :: ulatlon
! Subprogram not used real (kind=real_kind), dimension(np,np) :: v1,v2
! Subprogram not used 
! Subprogram not used ! local
! Subprogram not used integer :: k,ie
! Subprogram not used type (derivative_t)          :: deriv
! Subprogram not used 
! Subprogram not used call derivinit(deriv)
! Subprogram not used 
! Subprogram not used do k=1,nlev
! Subprogram not used do ie=nets,nete
! Subprogram not used     v1 = elem(ie)%state%v(:,:,1,k,nt)
! Subprogram not used     v2 = elem(ie)%state%v(:,:,2,k,nt)
! Subprogram not used     ulatlon(:,:,1) = elem(ie)%D(1,1,:,:)*v1 + elem(ie)%D(1,2,:,:)*v2
! Subprogram not used     ulatlon(:,:,2) = elem(ie)%D(2,1,:,:)*v1 + elem(ie)%D(2,2,:,:)*v2
! Subprogram not used    !    zeta(:,:,k,ie)=elem(ie)%state%zeta(:,:,k)
! Subprogram not used    zeta(:,:,k,ie)=divergence_sphere(ulatlon,deriv,elem(ie))
! Subprogram not used enddo
! Subprogram not used enddo
! Subprogram not used 
! Subprogram not used call make_C0(zeta,elem,hybrid,nets,nete)
! Subprogram not used 
! Subprogram not used end subroutine

! Subprogram not used subroutine compute_zeta_C0(zeta,elem,hybrid,nets,nete,nt)
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used ! compute C0 vorticity.  That is, solve:  
! Subprogram not used !     < PHI, zeta > = <PHI, curl(elem%state%v >
! Subprogram not used !
! Subprogram not used !    input:  v (stored in elem()%, in lat-lon coordinates)
! Subprogram not used !    output: zeta(:,:,:,:)   
! Subprogram not used !
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used 
! Subprogram not used type (hybrid_t)      , intent(in) :: hybrid
! Subprogram not used type (element_t)     , intent(in), target :: elem(:)
! Subprogram not used integer :: nt,nets,nete
! Subprogram not used real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: zeta
! Subprogram not used 
! Subprogram not used ! local
! Subprogram not used integer :: k,i,j,ie,ic
! Subprogram not used type (derivative_t)          :: deriv
! Subprogram not used 
! Subprogram not used call derivinit(deriv)
! Subprogram not used 
! Subprogram not used do ie=nets,nete
! Subprogram not used do k=1,nlev
! Subprogram not used    !    zeta(:,:,k,ie)=elem(ie)%state%zeta(:,:,k)
! Subprogram not used    zeta(:,:,k,ie)=vorticity_sphere(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie))
! Subprogram not used enddo
! Subprogram not used enddo
! Subprogram not used 
! Subprogram not used call make_C0(zeta,elem,hybrid,nets,nete)
! Subprogram not used 
! Subprogram not used end subroutine


! Subprogram not used subroutine compute_div_C0(zeta,elem,hybrid,nets,nete,nt)
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used ! compute C0 divergence. That is, solve:  
! Subprogram not used !     < PHI, zeta > = <PHI, div(elem%state%v >
! Subprogram not used !
! Subprogram not used !    input:  v (stored in elem()%, in lat-lon coordinates)
! Subprogram not used !    output: zeta(:,:,:,:)   
! Subprogram not used !
! Subprogram not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used 
! Subprogram not used type (hybrid_t)      , intent(in) :: hybrid
! Subprogram not used type (element_t)     , intent(in), target :: elem(:)
! Subprogram not used integer :: nt,nets,nete
! Subprogram not used real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: zeta
! Subprogram not used 
! Subprogram not used ! local
! Subprogram not used integer :: k,i,j,ie,ic
! Subprogram not used type (derivative_t)          :: deriv
! Subprogram not used 
! Subprogram not used call derivinit(deriv)
! Subprogram not used 
! Subprogram not used do ie=nets,nete
! Subprogram not used do k=1,nlev
! Subprogram not used    !    zeta(:,:,k,ie)=elem(ie)%state%zeta(:,:,k)
! Subprogram not used    zeta(:,:,k,ie)=divergence_sphere(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie))
! Subprogram not used enddo
! Subprogram not used enddo
! Subprogram not used 
! Subprogram not used call make_C0(zeta,elem,hybrid,nets,nete)
! Subprogram not used 
! Subprogram not used end subroutine








subroutine neighbor_minmax(elem,hybrid,edgeMinMax,nets,nete,min_neigh,max_neigh)
!
! compute Q min&max over the element and all its neighbors
!
!
integer :: nets,nete
type (element_t)     , intent(in) :: elem(:)
type (hybrid_t)      , intent(in) :: hybrid
type (EdgeBuffer_t)  , intent(inout) :: edgeMinMax
real (kind=real_kind) :: min_neigh(nlev,qsize,nets:nete)
real (kind=real_kind) :: max_neigh(nlev,qsize,nets:nete)

! local
integer :: ie,k,q
real (kind=real_kind) :: Qmin(np,np,nlev,qsize)
real (kind=real_kind) :: Qmax(np,np,nlev,qsize)


    ! compute Qmin, Qmax
    do ie=nets,nete
       do q=1,qsize	
          do k=1,nlev
             Qmin(:,:,k,q)=min_neigh(k,q,ie)
             Qmax(:,:,k,q)=max_neigh(k,q,ie)
          enddo
       end do
       call edgeVpack(edgeMinMax,Qmin,nlev*qsize,0,elem(ie)%desc)
       call edgeVpack(edgeMinMax,Qmax,nlev*qsize,nlev*qsize,elem(ie)%desc)
    enddo

    call bndry_exchangeV(hybrid,edgeMinMax)
       
    do ie=nets,nete
       do q=1,qsize	
          do k=1,nlev         
             Qmin(:,:,k,q)=min_neigh(k,q,ie) ! restore element data.  we could avoid
             Qmax(:,:,k,q)=max_neigh(k,q,ie) ! this by adding a "ie" index to Qmin/max
          enddo
       end do
! WARNING - edgeVunpackMin/Max take second argument as input/ouput
       call edgeVunpackMin(edgeMinMax,Qmin,nlev*qsize,0,elem(ie)%desc)
       call edgeVunpackMax(edgeMinMax,Qmax,nlev*qsize,nlev*qsize,elem(ie)%desc)
       do q=1,qsize
          do k=1,nlev
             ! note: only need to consider the corners, since the data we packed was
             ! constant within each element
             min_neigh(k,q,ie)=min(qmin(1,1,k,q),qmin(1,np,k,q),qmin(np,1,k,q),qmin(np,np,k,q))
             min_neigh(k,q,ie)=max(min_neigh(k,q,ie),0d0)
             max_neigh(k,q,ie)=max(qmax(1,1,k,q),qmax(1,np,k,q),qmax(np,1,k,q),qmax(np,np,k,q))
          enddo
       end do
    end do

end subroutine







! Subprogram not used   subroutine test_ibyp(elem, hybrid,  nets,   nete)
! Subprogram not used !
! Subprogram not used ! Note: vector test functions should be co-variant since u is contra-variant
! Subprogram not used !  PHIvec = PHIcov  (test function)
! Subprogram not used !  PHIcon = DtD PHIcov
! Subprogram not used !
! Subprogram not used ! weak grad:
! Subprogram not used !  < PHIcov du/dt > = < PHIcon grad(p) >    (output of grad is covariant)
! Subprogram not used !  < PHIcov du/dt > = -< div(PHIcon) p >    (input of div is contra)
! Subprogram not used !  verify:
! Subprogram not used !    gradient_sphere_wk(p) = - <div(PHIcon) p >
! Subprogram not used !    gradient_sphere_wk(p) + MASS*grad(p) = b.c. (b.c. are covariant)
! Subprogram not used !
! Subprogram not used ! weak div:
! Subprogram not used !   < PHI div(u) > = -< grad(PHI) dot u >     u=contra, output of grad is covariant
! Subprogram not used ! verify:
! Subprogram not used !   divergence_sphere_wk(u) = -<grad(PHI) dot u>
! Subprogram not used !   divergence_sphere_wk(u) + MASS*div(u) = b.c.  (b.c. are scalars)
! Subprogram not used !
! Subprogram not used ! weak curl:
! Subprogram not used !  < PHIcov du/dt > = < PHIcov curl( a ) >    (output of curl is contra)
! Subprogram not used !  < PHIcov du/dt > = < vor(PHIcov) a >       (input to vor is covariant)
! Subprogram not used ! verify:
! Subprogram not used !    curl_sphere_wk(a) = < vor(PHIcov) a >
! Subprogram not used !    curl_sphere_wk(a) - MASS*curl(a) = b.c. (b.c. are contra)
! Subprogram not used !
! Subprogram not used     ! ---------------------
! Subprogram not used     use kinds, only : real_kind
! Subprogram not used     ! ---------------------
! Subprogram not used     use physical_constants, only : rearth 
! Subprogram not used     ! ---------------------
! Subprogram not used     use dimensions_mod, only : np, nlev
! Subprogram not used     ! ---------------------
! Subprogram not used     use element_mod, only : element_t
! Subprogram not used     ! ---------------------
! Subprogram not used     use hybrid_mod, only : hybrid_t
! Subprogram not used     ! ---------------------
! Subprogram not used     use derivative_mod, only : derivative_t, gradient_sphere, divergence_sphere,vorticity_sphere,&
! Subprogram not used                                divergence_sphere_wk, curl_sphere
! Subprogram not used     use global_norms_mod
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     type (element_t)     , intent(inout), target :: elem(:)
! Subprogram not used 
! Subprogram not used     type (hybrid_t)      , intent(in) :: hybrid
! Subprogram not used 
! Subprogram not used     integer              , intent(in) :: nets
! Subprogram not used     integer              , intent(in) :: nete
! Subprogram not used 
! Subprogram not used   end subroutine test_ibyp









! Subprogram not used   subroutine check_edge_flux(elem,deriv,nets,nete)
! Subprogram not used !
! Subprogram not used !  check local element vector dentities:
! Subprogram not used !*****
! Subprogram not used !  1. div and weak div are adjoints: (for all scalar test functions)
! Subprogram not used !     integral[  p div(u) ] + integral[ grad(p) dot u ] = boundary_integral[ p u dot n]
! Subprogram not used !       PHI div(u) spheremp - div_wk(u)(i,j) = boundary_integral[ u PHI]
! Subprogram not used !       where PHI = the delta function at (i,j)
! Subprogram not used !
! Subprogram not used !*****
! Subprogram not used !  2. grad and weak grad are adjoints: 
! Subprogram not used !     weak gradient is defined with CONTRA vector test functions
! Subprogram not used !     i.e. it returns vector:   [  integral[ p div(PHIcontra_1) ]       
! Subprogram not used !                               [  integral[ p div(PHIcontra_2) ]       
! Subprogram not used !     
! Subprogram not used !   integral[  p div(u) ] + integral[ grad(p) dot u ] = boundary_integral[ p u dot n]
! Subprogram not used ! take u = PHIcontra_1 = (1,0) vector delta funciton at (i,j):
! Subprogram not used !  -grad_wk(p)_1(i,j) + spheremp PHIcontra_1 dot grad(p) = boundary_integral[ PHIcontra_1 p]
! Subprogram not used ! and then take u = PHIcontra_2 = (0,1) vector delta function at (i,j):
! Subprogram not used !  -grad_wk(p)_2(i,j) + spheremp PHIcontra_2 dot grad(p) = boundary_integral[ PHIcontra_2 p]
! Subprogram not used !
! Subprogram not used ! which is an equation for each covariant component:
! Subprogram not used ! -grad_wk(p)_cov1 + spheremp grad(p)_cov1 = boundary_integral[ PHIcontra_1 p dot n]
! Subprogram not used ! -grad_wk(p)_cov2 + spheremp grad(p)_cov2 = boundary_integral[ PHIcontra_2 p dot n]
! Subprogram not used !
! Subprogram not used ! HOMME-SE works in latlon, so convert cov->lat/lon:
! Subprogram not used !
! Subprogram not used ! -grad_wk(p) + spheremp grad(p) = D^-t * B 
! Subprogram not used !
! Subprogram not used ! with
! Subprogram not used !    B1 = boundary_integral[ PHIcontra_1 p] 
! Subprogram not used !    B2 = boundary_integral[ PHIcontra_2 p]
! Subprogram not used !
! Subprogram not used !*****
! Subprogram not used ! 3.  weak grid with COVARIANT test functions! 
! Subprogram not used !   integral[  p div(u) ] + integral[ grad(p) dot u ] = boundary_integral[ p u dot n]
! Subprogram not used ! take u = PHIcov_1 = (1,0) vector delta funciton at (i,j):
! Subprogram not used !  -grad_wk(p)_1(i,j) + spheremp PHIcov_1 dot grad(p) = boundary_integral[ PHIcov_1 p]
! Subprogram not used ! and then take u = PHIcov_2 = (0,1) vector delta function at (i,j):
! Subprogram not used !  -grad_wk(p)_2(i,j) + spheremp PHIcov_2 dot grad(p) = boundary_integral[ PHIcov_2 p]
! Subprogram not used !
! Subprogram not used ! which is an equation for each CONTRA component:
! Subprogram not used ! -grad_wk(p)_contra1 + spheremp grad(p)_contra1 = B1
! Subprogram not used ! -grad_wk(p)_contra2 + spheremp grad(p)_contra2 = B2
! Subprogram not used !
! Subprogram not used ! HOMME-SE works in latlon, so convert contra ->lat/lon:
! Subprogram not used !
! Subprogram not used ! -grad_wk(p) + spheremp grad(p) = D * B 
! Subprogram not used !
! Subprogram not used ! with
! Subprogram not used !    B1 = boundary_integral[ PHIcov_1 p] 
! Subprogram not used !    B2 = boundary_integral[ PHIcov_2 p]
! Subprogram not used !
! Subprogram not used !*****
! Subprogram not used ! 4.  weak curl with COVARIANT test functions! 
! Subprogram not used !  integral[ u dot curl(v)] - integral[v dot curl(u)] = boundary_integral[ v cross u dot n]
! Subprogram not used !  curl(p) = curl(p*khat) = horizontal vector
! Subprogram not used !  vor(U) =  s*khat       = (which we treat as a scalar)
! Subprogram not used !   integral[ p * vor(u)  ] - integral[ u dot curl(p) ] = boundary_integral[ u cross p*khat  dot n]
! Subprogram not used !
! Subprogram not used ! take u = PHIcov_1 = (1,0) vector delta funciton at (i,j):
! Subprogram not used !   curl_wk(p)_1(i,j) - spheremp PHIcov_1 dot curl(p) = boundary_integral[ perp(PHIcov_1) p]
! Subprogram not used ! and then take u = PHIcov_2 = (0,1) vector delta function at (i,j):
! Subprogram not used !   curl_wk(p)_2(i,j) - spheremp PHIcov_2 dot curl(p) = boundary_integral[ perp(PHIcov_2) p]
! Subprogram not used !
! Subprogram not used ! which is an equation for each CONTRA component:
! Subprogram not used ! curl_wk(p)_contra1 - spheremp curl(p)_contra1 = B1
! Subprogram not used ! curl_wk(p)_contra2 - spheremp curl(p)_contra2 = B2
! Subprogram not used !
! Subprogram not used ! HOMME-SE works in latlon, so convert contra ->lat/lon:
! Subprogram not used !
! Subprogram not used ! curl_wk(p) + spheremp curl(p) = D * B 
! Subprogram not used !
! Subprogram not used ! with
! Subprogram not used !    B1 = boundary_integral[ PHIcov_1 p] 
! Subprogram not used !    B2 = boundary_integral[ PHIcov_2 p]
! Subprogram not used !
! Subprogram not used   use dimensions_mod, only : np, np, nlev
! Subprogram not used   use element_mod, only    : element_t
! Subprogram not used   use derivative_mod, only  : derivative_t, divergence_sphere, divergence_sphere_wk, &
! Subprogram not used                              element_boundary_integral, gradient_sphere, &
! Subprogram not used                              gradient_sphere_wk_testcontra,gradient_sphere_wk_testcov, &
! Subprogram not used                              curl_sphere, curl_sphere_wk_testcov
! Subprogram not used   use physical_constants, only : rrearth
! Subprogram not used 
! Subprogram not used   implicit none
! Subprogram not used   
! Subprogram not used   type (element_t)     , intent(inout), target :: elem(:)
! Subprogram not used   type (derivative_t)  , intent(in) :: deriv
! Subprogram not used   integer :: nets,nete
! Subprogram not used   ! local 
! Subprogram not used   real (kind=real_kind), dimension(np,np,2) :: ucontra,ulatlon,gradp,gradp_wk,ucov
! Subprogram not used   real (kind=real_kind), dimension(np,np) :: phidivu,ugradphi,rhs,lhs,p
! Subprogram not used   real (kind=real_kind), dimension(np,np) :: rhs2,lhs2
! Subprogram not used   integer :: i,j,ie
! Subprogram not used 
! Subprogram not used   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used   print *,'integration by parts identity: check div/weak div:'
! Subprogram not used   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used   ! test integration by parts identity for each Cardinal function PHI:
! Subprogram not used   ! div(u)*spheremp - div_wk(u) = boundary integral phi u dot n
! Subprogram not used   do ie=nets,nete
! Subprogram not used      call random_number(ucontra)
! Subprogram not used      ! contra->latlon
! Subprogram not used      ulatlon(:,:,1)=(elem(ie)%D(1,1,:,:)*ucontra(:,:,1) + elem(ie)%D(1,2,:,:)*ucontra(:,:,2))
! Subprogram not used      ulatlon(:,:,2)=(elem(ie)%D(2,1,:,:)*ucontra(:,:,1) + elem(ie)%D(2,2,:,:)*ucontra(:,:,2))
! Subprogram not used      phidivu = elem(ie)%spheremp(:,:)*divergence_sphere(ulatlon,deriv,elem(ie))
! Subprogram not used      ugradphi = divergence_sphere_wk(ulatlon,deriv,elem(ie))
! Subprogram not used      lhs = phidivu - ugradphi
! Subprogram not used      
! Subprogram not used      rhs = element_boundary_integral(ulatlon,deriv,elem(ie))
! Subprogram not used      
! Subprogram not used      
! Subprogram not used      do j=1,np
! Subprogram not used         do i=1,np
! Subprogram not used            if ( abs(lhs(i,j)-rhs(i,j)) .gt. 1d-20) then
! Subprogram not used               write(*,'(a)') 'ERROR: div/div_wk integration by parts failure!'
! Subprogram not used               write(*,'(a,2i3,a,3e12.5)') 'for test function (i,j)=',i,j,' lhs,rhs=',lhs(i,j),rhs(i,j),lhs(i,j)-rhs(i,j)
! Subprogram not used            endif
! Subprogram not used         enddo
! Subprogram not used      enddo
! Subprogram not used   enddo
! Subprogram not used 
! Subprogram not used   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used   print *,'check grad/weak grad (gradient_sphere_wk_testcontra)'
! Subprogram not used   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used   ! PHIVEC = contra cardinal function 
! Subprogram not used   !          check each contra component seperately
! Subprogram not used 
! Subprogram not used   do ie=nets,nete
! Subprogram not used      call random_number(p)
! Subprogram not used      
! Subprogram not used      ! grad(p)  (lat/lon vector)
! Subprogram not used      gradp = gradient_sphere(p,deriv,elem(ie)%Dinv)
! Subprogram not used      gradp(:,:,1)=gradp(:,:,1)*elem(ie)%spheremp(:,:)  
! Subprogram not used      gradp(:,:,2)=gradp(:,:,2)*elem(ie)%spheremp(:,:)
! Subprogram not used      gradp_wk = gradient_sphere_wk_testcontra(p,deriv,elem(ie))
! Subprogram not used      
! Subprogram not used      ucontra(:,:,1)=p  ! PHIvec_1 * p
! Subprogram not used      ucontra(:,:,2)=0
! Subprogram not used      ! contra->latlon
! Subprogram not used      ulatlon(:,:,1)=(elem(ie)%D(1,1,:,:)*ucontra(:,:,1) + elem(ie)%D(1,2,:,:)*ucontra(:,:,2))
! Subprogram not used      ulatlon(:,:,2)=(elem(ie)%D(2,1,:,:)*ucontra(:,:,1) + elem(ie)%D(2,2,:,:)*ucontra(:,:,2))
! Subprogram not used 
! Subprogram not used      rhs = element_boundary_integral(ulatlon,deriv,elem(ie))
! Subprogram not used      lhs = gradp(:,:,1)-gradp_wk(:,:,1)
! Subprogram not used 
! Subprogram not used      ucontra(:,:,1)=0  ! PHIvec_2 * p
! Subprogram not used      ucontra(:,:,2)=p
! Subprogram not used      ! contra->latlon
! Subprogram not used      ulatlon(:,:,1)=(elem(ie)%D(1,1,:,:)*ucontra(:,:,1) + elem(ie)%D(1,2,:,:)*ucontra(:,:,2))
! Subprogram not used      ulatlon(:,:,2)=(elem(ie)%D(2,1,:,:)*ucontra(:,:,1) + elem(ie)%D(2,2,:,:)*ucontra(:,:,2))
! Subprogram not used      rhs2 = element_boundary_integral(ulatlon,deriv,elem(ie))
! Subprogram not used      lhs2 = gradp(:,:,2)-gradp_wk(:,:,2)  
! Subprogram not used 
! Subprogram not used 
! Subprogram not used      ! boundary integral gives covariant components. (see above) convert to latlon:
! Subprogram not used      ! cov -> latlon
! Subprogram not used      gradp(:,:,1)=rhs
! Subprogram not used      gradp(:,:,2)=rhs2
! Subprogram not used      rhs(:,:)=elem(ie)%Dinv(1,1,:,:)*gradp(:,:,1) + elem(ie)%Dinv(2,1,:,:)*gradp(:,:,2)
! Subprogram not used      rhs2(:,:)=elem(ie)%Dinv(1,2,:,:)*gradp(:,:,1) + elem(ie)%Dinv(2,2,:,:)*gradp(:,:,2)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used      do j=1,np
! Subprogram not used         do i=1,np
! Subprogram not used            if ( abs(lhs(i,j)-rhs(i,j)) .gt. 1d-20) then
! Subprogram not used               write(*,'(a)') 'ERROR: grad/grad_wk CONTRA (1) integration by parts failure!'
! Subprogram not used               write(*,'(a,2i3,a,4e12.4)') '(i,j)=',i,j,' lhs,rhs=',lhs(i,j),rhs(i,j),&
! Subprogram not used                    lhs(i,j)-rhs(i,j),lhs(i,j)/rhs(i,j)
! Subprogram not used            endif
! Subprogram not used         enddo
! Subprogram not used      enddo
! Subprogram not used 
! Subprogram not used 
! Subprogram not used      do j=1,np
! Subprogram not used         do i=1,np
! Subprogram not used            if ( abs(lhs2(i,j)-rhs2(i,j)) .gt. 1d-20) then
! Subprogram not used               write(*,'(a)') 'ERROR: grad/grad_wk CONTRA (2) integration by parts failure!'
! Subprogram not used               write(*,'(a,2i2,a,3e12.4)') '(i,j)=',i,j,' lhs,rhs=',lhs2(i,j),rhs2(i,j),lhs2(i,j)-rhs2(i,j)
! Subprogram not used            endif
! Subprogram not used         enddo
! Subprogram not used      enddo
! Subprogram not used   enddo
! Subprogram not used 
! Subprogram not used   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used   print *,'check grad/weak grad (gradient_sphere_wk_testcov)'
! Subprogram not used   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used   do ie=nets,nete
! Subprogram not used      call random_number(p)
! Subprogram not used 
! Subprogram not used      
! Subprogram not used      ! grad(p)  (lat/lon vector)
! Subprogram not used      gradp = gradient_sphere(p,deriv,elem(ie)%Dinv)
! Subprogram not used      gradp(:,:,1)=gradp(:,:,1)*elem(ie)%spheremp(:,:)  
! Subprogram not used      gradp(:,:,2)=gradp(:,:,2)*elem(ie)%spheremp(:,:)
! Subprogram not used      gradp_wk = gradient_sphere_wk_testcov(p,deriv,elem(ie))
! Subprogram not used      lhs = gradp(:,:,1)-gradp_wk(:,:,1)
! Subprogram not used      lhs2 = gradp(:,:,2)-gradp_wk(:,:,2)  
! Subprogram not used      
! Subprogram not used      ucov(:,:,1)=p  ! PHIvec_1 * p
! Subprogram not used      ucov(:,:,2)=0
! Subprogram not used      ! cov->latlon
! Subprogram not used      ulatlon(:,:,1)=(elem(ie)%Dinv(1,1,:,:)*ucov(:,:,1) + elem(ie)%Dinv(2,1,:,:)*ucov(:,:,2))
! Subprogram not used      ulatlon(:,:,2)=(elem(ie)%Dinv(1,2,:,:)*ucov(:,:,1) + elem(ie)%Dinv(2,2,:,:)*ucov(:,:,2))
! Subprogram not used      rhs = element_boundary_integral(ulatlon,deriv,elem(ie))
! Subprogram not used 
! Subprogram not used      ucov(:,:,1)=0  ! PHIvec_2 * p
! Subprogram not used      ucov(:,:,2)=p
! Subprogram not used      ! cov->latlon
! Subprogram not used      ulatlon(:,:,1)=(elem(ie)%Dinv(1,1,:,:)*ucov(:,:,1) + elem(ie)%Dinv(2,1,:,:)*ucov(:,:,2))
! Subprogram not used      ulatlon(:,:,2)=(elem(ie)%Dinv(1,2,:,:)*ucov(:,:,1) + elem(ie)%Dinv(2,2,:,:)*ucov(:,:,2))
! Subprogram not used      rhs2 = element_boundary_integral(ulatlon,deriv,elem(ie))
! Subprogram not used 
! Subprogram not used 
! Subprogram not used      ! boundary integral gives contra components. (see above) convert to latlon:
! Subprogram not used      ! contra -> latlon
! Subprogram not used      gradp(:,:,1)=rhs
! Subprogram not used      gradp(:,:,2)=rhs2
! Subprogram not used      rhs(:,:) =elem(ie)%D(1,1,:,:)*gradp(:,:,1) + elem(ie)%D(1,2,:,:)*gradp(:,:,2)
! Subprogram not used      rhs2(:,:)=elem(ie)%D(2,1,:,:)*gradp(:,:,1) + elem(ie)%D(2,2,:,:)*gradp(:,:,2)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used      do j=1,np
! Subprogram not used         do i=1,np
! Subprogram not used            if ( abs(lhs(i,j)-rhs(i,j)) .gt. 1d-20) then
! Subprogram not used               write(*,'(a)') 'ERROR: grad/grad_wk COV (1) integration by parts failure!'
! Subprogram not used               write(*,'(a,2i2,a,4e12.4)') '(i,j)=',i,j,' lhs,rhs=',lhs(i,j),rhs(i,j),lhs(i,j)-rhs(i,j),lhs(i,j)/rhs(i,j)
! Subprogram not used            endif
! Subprogram not used         enddo
! Subprogram not used      enddo
! Subprogram not used 
! Subprogram not used      do j=1,np
! Subprogram not used         do i=1,np
! Subprogram not used            if ( abs(lhs2(i,j)-rhs2(i,j)) .gt. 1d-20) then
! Subprogram not used               write(*,'(a)') 'ERROR: grad/grad_wk COV (2) integration by parts failure!'
! Subprogram not used               write(*,'(a,2i2,a,3e12.4)') '(i,j)=',i,j,' lhs,rhs=',lhs2(i,j),rhs2(i,j),lhs2(i,j)-rhs2(i,j)
! Subprogram not used            endif
! Subprogram not used         enddo
! Subprogram not used      enddo
! Subprogram not used   enddo
! Subprogram not used 
! Subprogram not used   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used   print *,'check curl/weak curl (curl_sphere_wk_testcov)'
! Subprogram not used   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used   do ie=nets,nete
! Subprogram not used      call random_number(p)
! Subprogram not used      
! Subprogram not used      ! grad(p)  (lat/lon vector)
! Subprogram not used      gradp = curl_sphere(p,deriv,elem(ie))
! Subprogram not used      gradp(:,:,1)=gradp(:,:,1)*elem(ie)%spheremp(:,:)  
! Subprogram not used      gradp(:,:,2)=gradp(:,:,2)*elem(ie)%spheremp(:,:)
! Subprogram not used      gradp_wk = curl_sphere_wk_testcov(p,deriv,elem(ie))
! Subprogram not used      lhs =  gradp_wk(:,:,1)-gradp(:,:,1)
! Subprogram not used      lhs2 = gradp_wk(:,:,2)-gradp(:,:,2)
! Subprogram not used      
! Subprogram not used      ucov(:,:,1)=p  ! PHIvec_1 * p
! Subprogram not used      ucov(:,:,2)=0
! Subprogram not used      ! cov->latlon, and then u cross khat:
! Subprogram not used      ulatlon(:,:,2)=-(elem(ie)%Dinv(1,1,:,:)*ucov(:,:,1) + elem(ie)%Dinv(2,1,:,:)*ucov(:,:,2))
! Subprogram not used      ulatlon(:,:,1)= (elem(ie)%Dinv(1,2,:,:)*ucov(:,:,1) + elem(ie)%Dinv(2,2,:,:)*ucov(:,:,2))
! Subprogram not used      rhs = element_boundary_integral(ulatlon,deriv,elem(ie))
! Subprogram not used 
! Subprogram not used      ucov(:,:,1)=0  ! PHIvec_2 * p
! Subprogram not used      ucov(:,:,2)=p
! Subprogram not used      ! cov->latlon, and u cross khat:
! Subprogram not used      ulatlon(:,:,2)=-(elem(ie)%Dinv(1,1,:,:)*ucov(:,:,1) + elem(ie)%Dinv(2,1,:,:)*ucov(:,:,2))
! Subprogram not used      ulatlon(:,:,1)= (elem(ie)%Dinv(1,2,:,:)*ucov(:,:,1) + elem(ie)%Dinv(2,2,:,:)*ucov(:,:,2))
! Subprogram not used      rhs2 = element_boundary_integral(ulatlon,deriv,elem(ie))
! Subprogram not used 
! Subprogram not used 
! Subprogram not used      ! boundary integral gives contra components. (see above) convert to latlon:
! Subprogram not used      ! contra -> latlon
! Subprogram not used      gradp(:,:,1)=rhs
! Subprogram not used      gradp(:,:,2)=rhs2
! Subprogram not used      rhs(:,:) =elem(ie)%D(1,1,:,:)*gradp(:,:,1) + elem(ie)%D(1,2,:,:)*gradp(:,:,2)
! Subprogram not used      rhs2(:,:)=elem(ie)%D(2,1,:,:)*gradp(:,:,1) + elem(ie)%D(2,2,:,:)*gradp(:,:,2)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used      do j=1,np
! Subprogram not used         do i=1,np
! Subprogram not used            if ( abs(lhs(i,j)-rhs(i,j)) .gt. 1d-20) then
! Subprogram not used               write(*,'(a)') 'ERROR: curl/curl_wk COV (1) integration by parts failure!'
! Subprogram not used               write(*,'(a,2i2,a,4e12.4)') '(i,j)=',i,j,' lhs,rhs=',lhs(i,j),rhs(i,j),lhs(i,j)-rhs(i,j),lhs(i,j)/rhs(i,j)
! Subprogram not used            endif
! Subprogram not used         enddo
! Subprogram not used      enddo
! Subprogram not used      stop
! Subprogram not used 
! Subprogram not used      do j=1,np
! Subprogram not used         do i=1,np
! Subprogram not used            if ( abs(lhs2(i,j)-rhs2(i,j)) .gt. 1d-20) then
! Subprogram not used               write(*,'(a)') 'ERROR: curl/curl_wk COV (2) integration by parts failure!'
! Subprogram not used               write(*,'(a,2i2,a,3e12.4)') '(i,j)=',i,j,' lhs,rhs=',lhs2(i,j),rhs2(i,j),lhs2(i,j)-rhs2(i,j)
! Subprogram not used            endif
! Subprogram not used         enddo
! Subprogram not used      enddo
! Subprogram not used   enddo
! Subprogram not used 
! Subprogram not used   print *,'done. integration by parts identity check:'
! Subprogram not used   end subroutine
end module
