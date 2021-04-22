



! To remove the effect of moisture comment the following
!#define _COMPUTE_MOISTURE_
module physics_mod
  
  ! =======================
  use kinds,              only : real_kind
  ! =======================
  use physical_constants, only : rgas, Rwater_vapor, kappa, g, Rd_on_Rv, Cp, Cpd_on_Cpv, cpwater_vapor
  ! =======================
  use physical_constants, only : rearth,p0
  ! =======================
  use dimensions_mod, only : np, nlev
  ! =======================
  use element_mod, only: timelevels
  implicit none
  
  private
  
  
  public :: Saturation_Vapor_Pressure
  public :: Specific_Humidity
  public :: Saturation_Specific_Humidity
  public :: Relative_Humidity
  public :: Vapor_Pressure
  public :: Mixing_Ratio
  public :: Prim_Condense
  public :: getsurfpress
  public :: Temp2PotTemp
  public :: Virtual_Temperature
  public :: Virtual_Specific_Heat
  public :: kappastar  

 interface Virtual_Temperature
    module procedure Virtual_Temperature1d
    module procedure Virtual_Temperature3d
 end interface


contains
  
  !===========================
  !
  ! For help or information:
  ! 
  ! Amik St-Cyr
  ! 
  ! e-mail: amik@ucar.edu
  !
  !===========================
 
  !================================
  ! For reference see Emanuel 1994 
  !================================
  
  function Virtual_Temperature1d(Tin,rin) result(Tv)
    
    real (kind=real_kind),intent(in) :: Tin
    real (kind=real_kind),intent(in) :: rin
    real (kind=real_kind)            :: Tv

!    Tv = Tin*(1_real_kind + rin/Rd_on_Rv)/(1_real_kind + rin)

    Tv = Tin*(1_real_kind + (Rwater_vapor/Rgas - 1.0_real_kind)*rin)


  end function Virtual_Temperature1d

! Subprogram not used   function Virtual_Temperature3d(T,Q) result(T_v)
! Subprogram not used     real (kind=real_kind),intent(in) :: T(np,np,nlev)
! Subprogram not used     real (kind=real_kind),intent(in) :: Q(np,np,nlev)
! Subprogram not used     real (kind=real_kind) :: T_v(np,np,nlev)
! Subprogram not used     integer :: i, j, k
! Subprogram not used 
! Subprogram not used     do k=1,nlev
! Subprogram not used        do j=1,np
! Subprogram not used           do i=1,np
! Subprogram not used              T_v(i,j,k) = Virtual_Temperature1d(T(i,j,k), Q(i,j,k))
! Subprogram not used           end do
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used   end function Virtual_Temperature3d

! Subprogram not used   function kappastar(Q) result(ks)
! Subprogram not used     real(kind=real_kind), intent(in) :: Q(np,np,nlev)
! Subprogram not used     real(kind=real_kind) :: ks(np,np,nlev)
! Subprogram not used     integer i,j,k
! Subprogram not used 
! Subprogram not used     do k=1,nlev
! Subprogram not used        do j=1,np
! Subprogram not used           do i=1,np
! Subprogram not used              ks(i,j,k) =  Rgas/Virtual_Specific_Heat(Q(i,j,k))
! Subprogram not used           end do
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used   end function kappastar


! Subprogram not used   function Virtual_Specific_Heat(rin) result(Cp_star)
! Subprogram not used     
! Subprogram not used     real (kind=real_kind),intent(in) :: rin
! Subprogram not used     real (kind=real_kind)            :: Cp_star
! Subprogram not used 
! Subprogram not used !    Cp_star = Cp*(1_real_kind + rin/Cpd_on_Cpv)/(1_real_kind + rin)
! Subprogram not used  
! Subprogram not used     Cp_star = Cp*(1.0_real_kind + (Cpwater_vapor/Cp - 1.0_real_kind)*rin)
! Subprogram not used    
! Subprogram not used   end function Virtual_Specific_Heat


  !=================================================
  ! Approx. Solution to the Clausius-Clapeyron eqn.
  !=================================================
! Subprogram not used   function Saturation_Vapor_Pressure(Tin) result(estar)
! Subprogram not used 
! Subprogram not used     real (kind=real_kind),intent(in) :: Tin
! Subprogram not used     real (kind=real_kind)            :: estar
! Subprogram not used 
! Subprogram not used     ! 4.4.13 p. 116 Emanuel
! Subprogram not used     
! Subprogram not used     estar = exp(53.67957_real_kind - 6743.769_real_kind/abs(Tin) &
! Subprogram not used 	- 4.8451_real_kind * log(ABS(Tin)))
! Subprogram not used     ! convert to code units
! Subprogram not used     if (p0 <  2000 ) then
! Subprogram not used        ! code is using mb, do nothing
! Subprogram not used     else
! Subprogram not used        ! code is using Pa.  convert from mb to Pa
! Subprogram not used        estar = estar*100
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used   end function Saturation_Vapor_Pressure

! Subprogram not used   function Specific_Humidity(r) result(q)
! Subprogram not used     
! Subprogram not used     real (kind=real_kind),intent(in) :: r
! Subprogram not used     real (kind=real_kind)            :: q
! Subprogram not used     
! Subprogram not used     q = r/( r + 1_real_kind )
! Subprogram not used 
! Subprogram not used   end  function Specific_Humidity
! Subprogram not used   
! Subprogram not used   function Saturation_Specific_Humidity(p,T) result(qstar)
! Subprogram not used  
! Subprogram not used     real (kind=real_kind),intent(in) :: p,T
! Subprogram not used     real (kind=real_kind)            :: qstar,estar
! Subprogram not used 
! Subprogram not used     estar = Saturation_Vapor_Pressure(T)
! Subprogram not used     qstar = Rd_on_Rv * estar / (p - estar * (1._real_kind - Rd_on_Rv))
! Subprogram not used   end function Saturation_Specific_Humidity
  !==================================
  ! Fraction of humid air in dry air 
  !==================================
! Subprogram not used   function Relative_Humidity(e,T) result(h)
! Subprogram not used     
! Subprogram not used     real (kind=real_kind),intent(in) :: e,T
! Subprogram not used     real (kind=real_kind)            :: estar
! Subprogram not used     real (kind=real_kind)            :: h
! Subprogram not used     
! Subprogram not used     estar = Saturation_Vapor_Pressure(T)
! Subprogram not used 
! Subprogram not used     h=e/estar
! Subprogram not used 
! Subprogram not used   end function Relative_Humidity

  !=================================
  ! Partial pressure of water vapor
  !=================================
! Subprogram not used   function Vapor_Pressure(r,p) result(e_out)
! Subprogram not used 
! Subprogram not used     real (kind=real_kind),intent(in) :: r,p
! Subprogram not used     real (kind=real_kind)            :: e_out
! Subprogram not used 
! Subprogram not used     ! 4.1.2 p. 108 Emanuel
! Subprogram not used 
! Subprogram not used     ! With p -> dry pressure
! Subprogram not used     e_out = abs(r*p)/Rd_on_Rv
! Subprogram not used 
! Subprogram not used   end function Vapor_Pressure

  !==============================================
  ! Mass of water vapor per unit mass of dry air
  !==============================================
! Subprogram not used   function Mixing_Ratio(ein,p) result(r)
! Subprogram not used     
! Subprogram not used     real (kind=real_kind),intent(in) :: ein,p
! Subprogram not used         
! Subprogram not used     real (kind=real_kind)            :: r
! Subprogram not used     
! Subprogram not used     ! 4.1.2 p. 108 Emanuel
! Subprogram not used     ! p is supposed DRY
! Subprogram not used     
! Subprogram not used     r = Rd_on_Rv*abs(ein/p)
! Subprogram not used 
! Subprogram not used   end function Mixing_Ratio

  !===============================================
  ! This function creates rain.
  ! If e > e_saturation then
  ! e = e_saturation.
  ! The latent heat is not computed here for now.
  !===============================================
! Subprogram not used   subroutine Prim_Condense(r,Tin,pin)
! Subprogram not used     
! Subprogram not used     real (kind=real_kind),intent(inout)  :: r
! Subprogram not used     real (kind=real_kind),intent(in)     :: Tin,pin
! Subprogram not used 
! Subprogram not used     real (kind=real_kind)                :: estar,e_vapor
! Subprogram not used     real*8                               :: st,et
! Subprogram not used 
! Subprogram not used     ! returns pressure in mb or Pa?  
! Subprogram not used     estar   = Saturation_Vapor_Pressure(Tin)     
! Subprogram not used 
! Subprogram not used     ! returns pressure in same units as pin. 
! Subprogram not used     e_vapor = Vapor_Pressure(r,pin)
! Subprogram not used 
! Subprogram not used     if(e_vapor/estar>1_real_kind)then
! Subprogram not used        e_vapor = estar
! Subprogram not used        r       = Mixing_Ratio(e_vapor,pin)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used   end subroutine Prim_Condense

  function Temp2PotTemp(pr3d,t3d) result(pt3d)
    real (kind=real_kind), parameter  :: c_p = 1004.64D0    ! Cp 
    real (kind=real_kind), parameter  :: p_0 = 100000.0D0   ! Initial Surface pressure
    real (kind=real_kind), parameter  :: r_d = 287.04D0      ! Gas const (dry)
    real (kind=real_kind),intent(in) :: pr3d(np,np,nlev),t3d(np,np,nlev)
    real (kind=real_kind)            :: pt3d(np,np,nlev)
    integer:: i,j,k    
    real (kind=real_kind):: rdcp,pp
    
    rdcp = r_d/c_p
    do k=1,nlev
       do j=1,np
          do i=1,np
             pp = (pr3d(i,j,k) + pr3d(i,j,k+1))*0.5D0
             pt3d(i,j,k)=  t3d(i,j,k)*(p_0/pp)**rdcp 
          enddo
       enddo
    enddo
  end function Temp2PotTemp

  function getsurfpress(lnps) result (press)
    real (kind=real_kind) :: press(np,np)
    real (kind=real_kind) :: lnps(np,np)

    press(:,:) = 0.0
    
  end function getsurfpress
  
     
  end module physics_mod
