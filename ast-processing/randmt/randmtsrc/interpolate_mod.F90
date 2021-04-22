module interpolate_mod
  use kinds, only : real_kind, iulog
  use element_mod, only : element_t
  use dimensions_mod, only : np, ne, nelemd, nc, nhe, nhc
  use quadrature_mod, only : quadrature_t, legendre, quad_norm
  use coordinate_systems_mod, only : spherical_polar_t, cartesian2d_t, &
       cartesian3D_t, sphere2cubedsphere, spherical_to_cart, &
       cubedsphere2cart, distance, change_coordinates, projectpoint
  use physical_constants,     only : DD_PI
  use quadrature_mod,         only : quadrature_t, gauss, gausslobatto
  use parallel_mod,           only : abortmp, syncmp, parallel_t, MPIreal_t, MPIinteger_t
  use parallel_mod,           only : MPI_MAX, MPI_SUM, MPI_MIN
  use cube_mod,               only : convert_gbl_index, dmap, ref2sphere
  use mesh_mod,               only : MeshUseMeshFile
  use control_mod,            only : cubed_sphere_map

  implicit none
  private
  integer, parameter, public :: MAX_VECVARS=25

  logical   :: debug=.false.

  character(len=10), public :: vector_uvars(MAX_VECVARS), vector_vvars(MAX_VECVARS)

  type, public :: interpolate_t
     real (kind=real_kind), dimension(:,:), pointer :: Imat  ! P_k(xj)*wj/gamma(k)
     real (kind=real_kind), dimension(:)  , pointer :: rk    ! 1/k
     real (kind=real_kind), dimension(:)  , pointer :: vtemp ! temp results
     real (kind=real_kind), dimension(:)  , pointer :: glp   ! GLL pts (nair)
  end type interpolate_t

  type, public :: interpdata_t
     ! Output Interpolation points.  Used to output data on lat-lon (or other grid)
     ! with native element interpolation.  Each element keeps a list of points from the
     ! interpolation grid that are in this element
     type (cartesian2D_t),pointer,dimension(:):: interp_xy      ! element coordinate
     integer, pointer,dimension(:)            :: ilat,ilon   ! position of interpolation point in lat-lon grid
     integer                                  :: n_interp
     integer                                  :: nlat
     integer                                  :: nlon
     logical                                  :: first_entry = .TRUE.
  end type interpdata_t

  real (kind=real_kind), private :: delta  = 1.0D-9  ! move tiny bit off center to
  ! avoid landing on element edges
  public :: interp_init
  public :: setup_latlon_interp
  public :: interpolate_scalar
  public :: interpolate_ce
  
  public :: interpol_phys_latlon, interpol_spelt_latlon
  public :: interpolate_vector
  public :: set_interp_parameter
  public :: get_interp_parameter
  public :: get_interp_gweight
  public :: get_interp_lat
  public :: get_interp_lon
  public :: var_is_vector_uvar, var_is_vector_vvar
  public :: cube_facepoint_ne
  public :: cube_facepoint_unstructured

  public :: interpolate_tracers
  public :: minmax_tracers
  public :: interpolate_2d
  public :: interpolate_create


  interface interpolate_scalar
     module procedure interpolate_scalar2d
     module procedure interpolate_scalar3d
  end interface
  interface interpolate_vector
     module procedure interpolate_vector2d
     module procedure interpolate_vector3d
  end interface

  type (interpolate_t), target ::  interp_p

  ! store the  lat-lon grid
  ! gridtype = 1       equally spaced, including poles (FV scalars output grid)
  ! gridtype = 2       Gauss grid (1 Eulerian)
  ! gridtype = 3       equally spaced, no poles (FV staggered velocity)
  ! Seven possible history files, last one is inithist and should be native grid
  logical, public, save :: interpolate_analysis(8) = (/.true.,.false.,.false.,.false.,.false.,.false.,.false.,.false./)
  integer :: nlat,nlon
  real (kind=real_kind), pointer,dimension(:)   :: lat(:),lon(:),gweight(:)
  integer :: gridtype = 1        !
  integer :: itype = 1           ! 0 = native high order
                                 ! 1 = bilinear

  integer :: auto_grid = 0        ! 0 = interpolation grid set by namelist
                                  ! 1 = grid set via mesh resolution

contains


  subroutine set_interp_parameter(parm_name, value)
    character*(*), intent(in) :: parm_name
    character(len=80) :: msg
    integer :: value,power
    real (kind=real_kind) :: value_target

    if(parm_name .eq. 'itype') then
       itype=value
    else if(parm_name .eq. 'nlon') then
       nlon=value
    else if(parm_name .eq. 'nlat') then
       nlat=value
    else if(parm_name.eq. 'gridtype') then
       gridtype=value
    else if(parm_name.eq. 'auto') then
       auto_grid=1
       ! compute recommended nlat,nlon which has slightly higher
       ! resolution than the specifed number of points around equator given in "value"
       ! computed recommended lat-lon grid.
       ! nlon > peq   peq = points around equator cubed sphere grid
       ! take nlon power of 2, and at most 1 power of 3
       if (value.eq.0) then
           ! If reading in unstructured mesh, ne = 0
           ! This makes it hard to guess how many interpolation points to use
           ! So We'll set the default as 720 x 360
           ! BUT if you're running with an unstructured mesh, set interp_nlon and interp_nlat
           nlon = 1536
           nlat = 768
       else
           value_target=value*1.25
           power = nint(.5 +  log( value_target)/log(2d0) )
           power = max(power,7) ! min grid: 64x128
           if ( 3*2**(power-2) > value_target) then
               nlon=3*2**(power-2)   ! use 1 power of 3
           else
               nlon=2**power
           endif
       endif
       nlat=nlon/2
       if (gridtype==1) nlat=nlat+1
    else
       write(msg,*) 'Did not recognize parameter named ',parm_name,' in interpolate_mod:set_interp_parameter'
       call abortmp(msg)
    end if
  end subroutine set_interp_parameter
! Subprogram not used   function get_interp_parameter(parm_name) result(value)
! Subprogram not used     character*(*), intent(in) :: parm_name
! Subprogram not used     integer :: value
! Subprogram not used     character(len=80) :: msg
! Subprogram not used     if(parm_name .eq. 'itype') then
! Subprogram not used        value=itype
! Subprogram not used     else if(parm_name .eq. 'nlon') then
! Subprogram not used        value=nlon
! Subprogram not used     else if(parm_name .eq. 'nlat') then
! Subprogram not used        value=nlat
! Subprogram not used     else if(parm_name.eq. 'gridtype') then
! Subprogram not used        value=gridtype
! Subprogram not used     else if(parm_name.eq. 'auto_grid') then
! Subprogram not used        value=auto_grid
! Subprogram not used     else
! Subprogram not used        write(msg,*) 'Did not recognize parameter named ',parm_name,' in interpolate_mod:get_interp_parameter'
! Subprogram not used        value=-1
! Subprogram not used        call abortmp(msg)
! Subprogram not used     end if
! Subprogram not used     return
! Subprogram not used   end function get_interp_parameter
! Subprogram not used   function get_interp_gweight() result(gw)
! Subprogram not used     real(kind=real_kind) :: gw(nlat)
! Subprogram not used     gw=gweight
! Subprogram not used     return
! Subprogram not used   end function get_interp_gweight
! Subprogram not used   function get_interp_lat() result(thislat)
! Subprogram not used     use physical_constants, only : DD_PI
! Subprogram not used     real(kind=real_kind) :: thislat(nlat)
! Subprogram not used     thislat=lat*180.0D0/DD_PI
! Subprogram not used     return
! Subprogram not used   end function get_interp_lat
! Subprogram not used   function get_interp_lon() result(thislon)
! Subprogram not used     use physical_constants, only : DD_PI
! Subprogram not used     real(kind=real_kind) :: thislon(nlon)
! Subprogram not used     thislon=lon*180.0D0/DD_PI
! Subprogram not used     return
! Subprogram not used   end function get_interp_lon

! Subprogram not used   subroutine interpolate_create(gquad,interp)
! Subprogram not used     type (quadrature_t) , intent(in)   :: gquad
! Subprogram not used     type (interpolate_t), intent(out)  :: interp
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     ! Local variables
! Subprogram not used 
! Subprogram not used     integer k,j
! Subprogram not used     integer npts
! Subprogram not used     real (kind=real_kind), dimension(:), allocatable :: gamma
! Subprogram not used     real (kind=real_kind), dimension(:), allocatable :: leg
! Subprogram not used 
! Subprogram not used     npts = size(gquad%points)
! Subprogram not used 
! Subprogram not used     allocate(interp%Imat(npts,npts))
! Subprogram not used     allocate(interp%rk(npts))
! Subprogram not used     allocate(interp%vtemp(npts))
! Subprogram not used     allocate(interp%glp(npts))
! Subprogram not used     allocate(gamma(npts))
! Subprogram not used     allocate(leg(npts))
! Subprogram not used 
! Subprogram not used     gamma = quad_norm(gquad,npts)
! Subprogram not used 
! Subprogram not used     do k=1,npts
! Subprogram not used        interp%rk(k) = 1.0D0/k
! Subprogram not used        interp%glp(k) = gquad%points(k)    !nair
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     do j=1,npts
! Subprogram not used        leg=legendre(gquad%points(j),npts-1)
! Subprogram not used        do k=1,npts
! Subprogram not used           interp%Imat(j,k)=leg(k)*gquad%weights(j)/gamma(k)
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     deallocate(gamma)
! Subprogram not used     deallocate(leg)
! Subprogram not used 
! Subprogram not used   end subroutine interpolate_create


! Subprogram not used   subroutine interpolate_tracers(r, tracers, f) 
! Subprogram not used     use kinds,          only : longdouble_kind
! Subprogram not used     use dimensions_mod, only : np, qsize
! Subprogram not used     use quadrature_mod, only : quadrature_t, gausslobatto
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     type (cartesian2D_t), intent(in)  :: r
! Subprogram not used     real (kind=real_kind),intent(in)  :: tracers(np*np,qsize)
! Subprogram not used     real (kind=real_kind),intent(out) :: f(qsize)
! Subprogram not used 
! Subprogram not used     type (quadrature_t        )       :: gll        
! Subprogram not used     real (kind=real_kind      )       :: dp    (np)
! Subprogram not used     real (kind=real_kind      )       :: x     (np)
! Subprogram not used     real (kind=real_kind      )       :: y     (np)
! Subprogram not used     real (kind=real_kind      )       :: c     (np,np)
! Subprogram not used     real (kind=real_kind      )       :: xy    (np*np)
! Subprogram not used 
! Subprogram not used     integer                           :: i,j
! Subprogram not used     logical                           :: first_time=.true.
! Subprogram not used     
! Subprogram not used     save c
! Subprogram not used     save gll
! Subprogram not used    
! Subprogram not used     if (first_time) then
! Subprogram not used       first_time = .false.
! Subprogram not used       gll=gausslobatto(np)
! Subprogram not used       dp = 1
! Subprogram not used       do i=1,np
! Subprogram not used         do j=1,np
! Subprogram not used           if (i /= j) then
! Subprogram not used             dp(i) = dp(i) * (gll%points(i) - gll%points(j))
! Subprogram not used           end if
! Subprogram not used         end do
! Subprogram not used       end do 
! Subprogram not used       do i=1,np
! Subprogram not used         do j=1,np
! Subprogram not used           c(i,j) = 1/(dp(i)*dp(j))
! Subprogram not used         end do
! Subprogram not used       end do 
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     x = 1
! Subprogram not used     y = 1
! Subprogram not used     do i=1,np
! Subprogram not used       do j=1,np
! Subprogram not used         if (i /= j) then
! Subprogram not used           x(i) = x(i) * (r%x - gll%points(j))
! Subprogram not used           y(i) = y(i) * (r%y - gll%points(j))
! Subprogram not used         end if
! Subprogram not used       end do
! Subprogram not used     end do 
! Subprogram not used     do j=1,np  
! Subprogram not used       do i=1,np
! Subprogram not used         xy(i + (j-1)*np) = x(i)*y(j)*c(i,j)
! Subprogram not used       end do
! Subprogram not used     end do 
! Subprogram not used     f = MATMUL(xy,tracers)
! Subprogram not used   end subroutine interpolate_tracers

! Subprogram not used   subroutine minmax_tracers(r, tracers, mint, maxt) 
! Subprogram not used     use dimensions_mod, only : np, qsize
! Subprogram not used     use quadrature_mod, only : quadrature_t, gausslobatto
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     type (cartesian2D_t), intent(in)  :: r
! Subprogram not used     real (kind=real_kind),intent(in)  :: tracers(np,np,qsize)
! Subprogram not used     real (kind=real_kind),intent(out) :: mint(qsize)
! Subprogram not used     real (kind=real_kind),intent(out) :: maxt(qsize)
! Subprogram not used 
! Subprogram not used     type (quadrature_t        )       :: gll        
! Subprogram not used     integer                           :: i,j,k
! Subprogram not used     logical                           :: first_time=.true.
! Subprogram not used     
! Subprogram not used     save gll
! Subprogram not used    
! Subprogram not used     if (first_time) then
! Subprogram not used       first_time = .false.
! Subprogram not used       gll=gausslobatto(np)
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     do i=1,np  
! Subprogram not used       if (r%x < gll%points(i)) exit
! Subprogram not used     end do 
! Subprogram not used     do j=1,np  
! Subprogram not used       if (r%y < gll%points(j)) exit
! Subprogram not used     end do 
! Subprogram not used     if (1 < i) i = i-1
! Subprogram not used     if (1 < j) j = j-1
! Subprogram not used     mint(:) = minval(minval(tracers(i:i+1,j:j+1,:),1),1)
! Subprogram not used     maxt(:) = maxval(maxval(tracers(i:i+1,j:j+1,:),1),1)
! Subprogram not used   end subroutine minmax_tracers

! Subprogram not used   function interpolate_2d(cart, f, interp, npts, fillvalue) result(fxy)
! Subprogram not used     integer, intent(in)               :: npts
! Subprogram not used     type (cartesian2D_t), intent(in)  :: cart
! Subprogram not used     real (kind=real_kind), intent(in) :: f(npts,npts)
! Subprogram not used     type (interpolate_t)              :: interp
! Subprogram not used     real (kind=real_kind)             :: fxy     ! value of f interpolated to (x,y)
! Subprogram not used     real (kind=real_kind), intent(in), optional :: fillvalue
! Subprogram not used     ! local variables
! Subprogram not used 
! Subprogram not used     real (kind=real_kind)             :: tmp_1,tmp_2
! Subprogram not used     real (kind=real_kind)             :: fk0,fk1
! Subprogram not used     real (kind=real_kind)             :: pk
! Subprogram not used 
! Subprogram not used     integer                           :: l,j,k
! Subprogram not used 
! Subprogram not used     if(present(fillvalue)) then
! Subprogram not used        if (any(f==fillvalue)) then
! Subprogram not used           fxy = fillvalue
! Subprogram not used           return
! Subprogram not used        endif
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     do l=1,npts,2
! Subprogram not used 
! Subprogram not used        ! Compute Pk(cart%x) for Legendre order 0
! Subprogram not used 
! Subprogram not used        pk = 1.0D0
! Subprogram not used 
! Subprogram not used        fk0=0.0D0
! Subprogram not used        fk1=0.0D0
! Subprogram not used        do j=1,npts
! Subprogram not used           fk0 = fk0 + interp%Imat(j,1)*f(j,l  )
! Subprogram not used           fk1 = fk1 + interp%Imat(j,1)*f(j,l+1)
! Subprogram not used        end do
! Subprogram not used        interp%vtemp(l  ) = pk*fk0
! Subprogram not used        interp%vtemp(l+1) = pk*fk1
! Subprogram not used 
! Subprogram not used        ! Compute Pk(cart%x) for Legendre order 1
! Subprogram not used 
! Subprogram not used        tmp_2 = pk
! Subprogram not used        pk    = cart%x
! Subprogram not used 
! Subprogram not used        fk0=0.0D0
! Subprogram not used        fk1=0.0D0
! Subprogram not used        do j=1,npts
! Subprogram not used           fk0 = fk0 + interp%Imat(j,2)*f(j,l  )
! Subprogram not used           fk1 = fk1 + interp%Imat(j,2)*f(j,l+1)
! Subprogram not used        end do
! Subprogram not used        interp%vtemp(l  ) = interp%vtemp(l  ) + pk*fk0
! Subprogram not used        interp%vtemp(l+1) = interp%vtemp(l+1) + pk*fk1
! Subprogram not used 
! Subprogram not used        ! Compute Pk(cart%x) for Legendre order 2 to npts-1
! Subprogram not used 
! Subprogram not used        do k = 2,npts-1
! Subprogram not used 
! Subprogram not used           tmp_1  = tmp_2
! Subprogram not used           tmp_2  = pk
! Subprogram not used           pk = ( (2*k-1)*cart%x*tmp_2 - (k-1)*tmp_1 )*interp%rk(k)
! Subprogram not used 
! Subprogram not used           fk0=0.0D0
! Subprogram not used           fk1=0.0D0
! Subprogram not used           do j=1,npts
! Subprogram not used              fk0 = fk0 + interp%Imat(j,k+1)*f(j,l  )
! Subprogram not used              fk1 = fk1 + interp%Imat(j,k+1)*f(j,l+1)
! Subprogram not used           end do
! Subprogram not used           interp%vtemp(l  ) = interp%vtemp(l  ) + pk*fk0
! Subprogram not used           interp%vtemp(l+1) = interp%vtemp(l+1) + pk*fk1
! Subprogram not used 
! Subprogram not used        end do
! Subprogram not used 
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     ! Compute Pk(cart%y) for Legendre order 0
! Subprogram not used 
! Subprogram not used     pk = 1.0
! Subprogram not used 
! Subprogram not used     fk0 = 0.0D0
! Subprogram not used     do j=1,npts
! Subprogram not used        fk0 = fk0 + interp%Imat(j,1)*interp%vtemp(j)
! Subprogram not used     end do
! Subprogram not used     fxy = pk*fk0
! Subprogram not used 
! Subprogram not used     ! Compute Pk(cart%y) for Legendre order 1
! Subprogram not used 
! Subprogram not used     tmp_2 = pk
! Subprogram not used     pk    = cart%y
! Subprogram not used 
! Subprogram not used     fk0=0.0D0
! Subprogram not used     do j=1,npts
! Subprogram not used        fk0 = fk0 + interp%Imat(j,2)*interp%vtemp(j)
! Subprogram not used     end do
! Subprogram not used     fxy = fxy + pk*fk0
! Subprogram not used 
! Subprogram not used     ! Compute Pk(cart%y) for Legendre order 2, npts-1
! Subprogram not used 
! Subprogram not used     do k = 2,npts-1
! Subprogram not used        tmp_1  = tmp_2
! Subprogram not used        tmp_2  = pk
! Subprogram not used        pk = ( (2*k-1)*cart%y*tmp_2 - (k-1)*tmp_1 )*interp%rk(k)
! Subprogram not used 
! Subprogram not used        fk0 = 0.0D0
! Subprogram not used        do j=1,npts
! Subprogram not used           fk0 = fk0 + interp%Imat(j,k+1)*interp%vtemp(j)
! Subprogram not used        end do
! Subprogram not used 
! Subprogram not used        fxy = fxy + pk*fk0
! Subprogram not used 
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used   end function interpolate_2d

  !===============================
  !(Nair) Bilinear interpolation for every GLL grid cell
  !===============================

! Subprogram not used   function interpol_bilinear(cart, f, interp, npts, fillvalue) result(fxy)
! Subprogram not used 
! Subprogram not used     integer, intent(in)               :: npts
! Subprogram not used     type (cartesian2D_t), intent(in)  :: cart
! Subprogram not used     real (kind=real_kind), intent(in) :: f(npts,npts)
! Subprogram not used     type (interpolate_t)              :: interp
! Subprogram not used     real (kind=real_kind)             :: fxy     ! value of f interpolated to (x,y)
! Subprogram not used     real (kind=real_kind), intent(in), optional :: fillvalue
! Subprogram not used     ! local variables
! Subprogram not used 
! Subprogram not used     real (kind=real_kind)             :: xoy(npts)
! Subprogram not used     real (kind=real_kind)             :: p,q,xp,yp ,y4(4)
! Subprogram not used 
! Subprogram not used     integer                           :: l,j,k, ii, jj, na,nb,nm
! Subprogram not used 
! Subprogram not used     xp = cart%x
! Subprogram not used     yp = cart%y
! Subprogram not used 
! Subprogram not used     xoy(:) = interp%glp(:)
! Subprogram not used 
! Subprogram not used     ! Search index along "x"  (bisection method)
! Subprogram not used 
! Subprogram not used     na = 1
! Subprogram not used     nb = npts
! Subprogram not used     do
! Subprogram not used        if  ((nb-na) <=  1)  exit
! Subprogram not used        nm = (nb + na)/2
! Subprogram not used        if (xp  >  xoy(nm)) then
! Subprogram not used           na = nm
! Subprogram not used        else
! Subprogram not used           nb = nm
! Subprogram not used        endif
! Subprogram not used     enddo
! Subprogram not used     ii = na
! Subprogram not used 
! Subprogram not used     ! Search index along "y"
! Subprogram not used 
! Subprogram not used     na = 1
! Subprogram not used     nb = npts
! Subprogram not used     do
! Subprogram not used        if  ((nb-na) <=  1)  exit
! Subprogram not used        nm = (nb + na)/2
! Subprogram not used        if (yp  >  xoy(nm)) then
! Subprogram not used           na = nm
! Subprogram not used        else
! Subprogram not used           nb = nm
! Subprogram not used        endif
! Subprogram not used     enddo
! Subprogram not used     jj = na
! Subprogram not used 
! Subprogram not used     ! GLL cell containing (xp,yp)
! Subprogram not used 
! Subprogram not used     y4(1) = f(ii,jj)
! Subprogram not used     y4(2) = f(ii+1,jj)
! Subprogram not used     y4(3) = f(ii+1,jj+1)
! Subprogram not used     y4(4) = f(ii,jj+1)
! Subprogram not used 
! Subprogram not used     if(present(fillvalue)) then
! Subprogram not used        if (any(y4==fillvalue)) then
! Subprogram not used           fxy = fillvalue
! Subprogram not used           return
! Subprogram not used        endif
! Subprogram not used     endif
! Subprogram not used        p = (xp - xoy(ii))/(xoy(ii+1) - xoy(ii))
! Subprogram not used        q = (yp - xoy(jj))/(xoy(jj+1) - xoy(jj))
! Subprogram not used 
! Subprogram not used        fxy = (1.0D0 - p)*(1.0D0 - q)* y4(1) + p*(1.0D0 - q) * y4(2)   &
! Subprogram not used             + p*q* y4(3) + (1.0D0 - p)*q * y4(4)
! Subprogram not used 
! Subprogram not used   end function interpol_bilinear

! ----------------------------------------------------------------------------------!
!FUNCTION   interpol_phys_latlon----------------------------------------CE-for fvm!
! AUTHOR: CHRISTOPH ERATH, 23. May 2012                                             !
! DESCRIPTION: evaluation of the reconstruction for every physics grid cell         !
!                                                                                   !
! CALLS: 
! INPUT: 
!        
! OUTPUT: 
!-----------------------------------------------------------------------------------!
! Subprogram not used subroutine interpol_phys_latlon(interpdata,f, fvm, corners, desc, flatlon)
! Subprogram not used   use fvm_control_volume_mod, only : fvm_struct
! Subprogram not used   use fvm_reconstruction_mod, only: reconstruction
! Subprogram not used   use fvm_filter_mod, only: monotonic_gradient_cart, recons_val_cart
! Subprogram not used   use edge_mod, only : edgedescriptor_t
! Subprogram not used   
! Subprogram not used   type (interpdata_t), intent(in)     :: interpdata                        
! Subprogram not used   real (kind=real_kind), intent(in)   :: f(1-nhc:nc+nhc,1-nhc:nc+nhc)
! Subprogram not used   type (fvm_struct), intent(in)       :: fvm
! Subprogram not used   type (cartesian2d_t), intent(in)    :: corners(:)
! Subprogram not used   type (edgedescriptor_t),intent(in)  :: desc
! Subprogram not used   
! Subprogram not used                           
! Subprogram not used   real (kind=real_kind)             :: flatlon(:)
! Subprogram not used   ! local variables
! Subprogram not used   real (kind=real_kind)             :: xp,yp, tmpval
! Subprogram not used   real (kind=real_kind)             :: tmpaxp,tmpaxm, tmpayp, tmpaym
! Subprogram not used   integer                           :: i, ix, jy, starti,endi,tmpi
! Subprogram not used   real (kind=real_kind), dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe)      :: recons
! Subprogram not used   
! Subprogram not used   call reconstruction(f, fvm,recons)
! Subprogram not used   call monotonic_gradient_cart(f, fvm,recons, desc)
! Subprogram not used 
! Subprogram not used   tmpaxp=(corners(1)%x+corners(2)%x)/2
! Subprogram not used   tmpaxm=(corners(2)%x-corners(1)%x)/2
! Subprogram not used   tmpayp=(corners(1)%y+corners(4)%y)/2
! Subprogram not used   tmpaym=(corners(4)%y-corners(1)%y)/2
! Subprogram not used   do i=1,interpdata%n_interp
! Subprogram not used     ! caculation phys grid coordinate of xp point, note the interp_xy are on the reference [-1,1]x[-1,1]
! Subprogram not used     xp=tan(tmpaxp+interpdata%interp_xy(i)%x*tmpaxm)
! Subprogram not used     yp=tan(tmpayp+interpdata%interp_xy(i)%y*tmpaym)   
! Subprogram not used 
! Subprogram not used     ! Search index along "x"  (bisection method)
! Subprogram not used     starti = 1
! Subprogram not used     endi = nc+1
! Subprogram not used     do
! Subprogram not used        if  ((endi-starti) <=  1)  exit
! Subprogram not used        tmpi = (endi + starti)/2
! Subprogram not used        if (xp  >  fvm%acartx(tmpi)) then
! Subprogram not used           starti = tmpi
! Subprogram not used        else
! Subprogram not used           endi = tmpi
! Subprogram not used        endif
! Subprogram not used     enddo
! Subprogram not used     ix = starti
! Subprogram not used 
! Subprogram not used   ! Search index along "y"
! Subprogram not used     starti = 1
! Subprogram not used     endi = nc+1
! Subprogram not used     do
! Subprogram not used        if  ((endi-starti) <=  1)  exit
! Subprogram not used        tmpi = (endi + starti)/2
! Subprogram not used        if (yp  >  fvm%acarty(tmpi)) then
! Subprogram not used           starti = tmpi
! Subprogram not used        else
! Subprogram not used           endi = tmpi
! Subprogram not used        endif
! Subprogram not used     enddo
! Subprogram not used     jy = starti
! Subprogram not used     
! Subprogram not used     call recons_val_cart(f, xp,yp,fvm%spherecentroid,recons,ix,jy,tmpval)
! Subprogram not used     flatlon(i)=tmpval    
! Subprogram not used !     flatlon(i)=f(ix,jy)    
! Subprogram not used   end do
! Subprogram not used end subroutine interpol_phys_latlon


! ----------------------------------------------------------------------------------!
!FUNCTION   interpol_spelt_latlon---------------------------------------CE-for spelt!
! AUTHOR: CHRISTOPH ERATH, 24. August 2012                                          !
! DESCRIPTION: evaluation of the reconstruction for every spelt grid cell         !
!                                                                                   !
! CALLS: 
! INPUT: 
!        
! OUTPUT: 
!-----------------------------------------------------------------------------------!
! Subprogram not used subroutine interpol_spelt_latlon(interpdata,f, spelt,corners, flatlon)
! Subprogram not used   use spelt_mod, only : spelt_struct, cell_search, cip_coeff, cell_minmax, &
! Subprogram not used                         cip_interpolate, qmsl_cell_filter, metric_term
! Subprogram not used   use dimensions_mod, only: nip, nipm, nep
! Subprogram not used   use coordinate_systems_mod, only : cartesian2d_t
! Subprogram not used   
! Subprogram not used   use edge_mod, only : edgedescriptor_t
! Subprogram not used   
! Subprogram not used   type (interpdata_t), intent(in)     :: interpdata                        
! Subprogram not used   real (kind=real_kind), intent(in)   :: f(1-nipm:nep+nipm,1-nipm:nep+nipm)
! Subprogram not used   type (spelt_struct), intent(in)     :: spelt  
! Subprogram not used   type (cartesian2d_t), intent(in)    :: corners(4)
! Subprogram not used                             
! Subprogram not used   real (kind=real_kind)             :: flatlon(:)
! Subprogram not used   ! local variables
! Subprogram not used   real (kind=real_kind)             :: xp,yp,dxp,dyp, tmpval
! Subprogram not used   integer                           :: i, j, ix,jy, icell, jcell, starti,endi,tmpi
! Subprogram not used   real (kind=real_kind)             :: cf(nip,nip,1:nc,1:nc)
! Subprogram not used   real (kind=real_kind)             :: ff(nip,nip)
! Subprogram not used   real (kind=real_kind)             :: minmax(1-nhe:nc+nhe,1-nhe:nc+nhe,2)
! Subprogram not used   
! Subprogram not used   type (cartesian2d_t)              :: alphabeta   
! Subprogram not used   real(kind=real_kind)              :: pi,pj,qi,qj, sga, tmp
! Subprogram not used   
! Subprogram not used 
! Subprogram not used   do j=1,nc
! Subprogram not used     do i=1,nc
! Subprogram not used       icell=1+(i-1)*nipm
! Subprogram not used       jcell=1+(j-1)*nipm
! Subprogram not used       ff=f(icell:icell+nipm,jcell:jcell+nipm)*spelt%sga(icell:icell+nipm,jcell:jcell+nipm)
! Subprogram not used       minmax(i,j,:)=cell_minmax(ff)
! Subprogram not used       call cip_coeff(spelt%drefx(i,j),spelt%drefy(i,j),ff,ff(2,2),cf(:,:,i,j))
! Subprogram not used     enddo
! Subprogram not used   enddo
! Subprogram not used ! 
! Subprogram not used   do i=1,interpdata%n_interp
! Subprogram not used     ! caculation phys grid coordinate of xp point, note the interp_xy are on the reference [-1,1]x[-1,1] 
! Subprogram not used     xp=interpdata%interp_xy(i)%x
! Subprogram not used     yp=interpdata%interp_xy(i)%y
! Subprogram not used     ! Search index along "x"  (bisection method)
! Subprogram not used     starti = 1
! Subprogram not used     endi = nc+1
! Subprogram not used     do
! Subprogram not used        if  ((endi-starti) <=  1)  exit
! Subprogram not used        tmpi = (endi + starti)/2
! Subprogram not used        if (xp  >  spelt%pref(tmpi)) then
! Subprogram not used           starti = tmpi
! Subprogram not used        else
! Subprogram not used           endi = tmpi
! Subprogram not used        endif
! Subprogram not used     enddo
! Subprogram not used     icell = starti
! Subprogram not used 
! Subprogram not used   ! Search index along "y"
! Subprogram not used     starti = 1
! Subprogram not used     endi = nc+1
! Subprogram not used     do
! Subprogram not used        if  ((endi-starti) <=  1)  exit
! Subprogram not used        tmpi = (endi + starti)/2
! Subprogram not used        if (yp  >  spelt%pref(tmpi)) then
! Subprogram not used           starti = tmpi
! Subprogram not used        else
! Subprogram not used           endi = tmpi
! Subprogram not used        endif
! Subprogram not used     enddo
! Subprogram not used     jcell = starti
! Subprogram not used     
! Subprogram not used     if ((icell<1) .or.(icell>nc) .or. (jcell<1) .or. (jcell>nc)) then
! Subprogram not used       write(*,*) 'icell, jcell,Something is wrong in the search of interpol_spelt_latlon!'
! Subprogram not used       stop
! Subprogram not used     endif
! Subprogram not used     dxp=xp-spelt%pref(icell)
! Subprogram not used     dyp=yp-spelt%pref(jcell)
! Subprogram not used     tmp=cip_interpolate(cf(:,:,icell,jcell),dxp,dyp)      
! Subprogram not used     tmp=qmsl_cell_filter(icell,jcell,minmax,tmp) 
! Subprogram not used     
! Subprogram not used     !next lines can be deleted, once the subroutine for the metric term for ref points works
! Subprogram not used     pi = (1-xp)/2
! Subprogram not used     pj = (1-yp)/2
! Subprogram not used     qi = (1+xp)/2
! Subprogram not used     qj = (1+yp)/2
! Subprogram not used     alphabeta%x = pi*pj*corners(1)%x &
! Subprogram not used          + qi*pj*corners(2)%x &
! Subprogram not used          + qi*qj*corners(3)%x &
! Subprogram not used          + pi*qj*corners(4)%x 
! Subprogram not used     alphabeta%y = pi*pj*corners(1)%y &
! Subprogram not used          + qi*pj*corners(2)%y &
! Subprogram not used          + qi*qj*corners(3)%y &
! Subprogram not used          + pi*qj*corners(4)%y
! Subprogram not used 
! Subprogram not used     sga=metric_term(alphabeta)
! Subprogram not used     
! Subprogram not used     flatlon(i)=tmp/sga
! Subprogram not used !     flatlon(i)=f(icell*nipm,jcell*nipm)    
! Subprogram not used   end do
! Subprogram not used end subroutine interpol_spelt_latlon



! Subprogram not used   function parametric_coordinates(sphere, elem) result (ref)
! Subprogram not used     implicit none
! Subprogram not used     type (spherical_polar_t), intent(in) :: sphere
! Subprogram not used     type (element_t), intent(in) :: elem
! Subprogram not used     type (cartesian2D_t) :: ref
! Subprogram not used 
! Subprogram not used     ! local
! Subprogram not used     integer               :: face_no, i, MAX_NR_ITER=10
! Subprogram not used     real(kind=real_kind)  :: D(2,2),Dinv(2,2),detD,a,b,resa,resb,dela,delb,costh
! Subprogram not used     real(kind=real_kind)  :: tol_sq=1e-26
! Subprogram not used     type (spherical_polar_t) :: sphere1, sphere_tmp
! Subprogram not used 
! Subprogram not used     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used     ! newton iteration on: ref=ref - df^-1 (ref2sphere(ref) - sphere)
! Subprogram not used     !
! Subprogram not used     ! Generic version written in terms of HOMME's 'ref2sphere' and 'Dmap' operaters,
! Subprogram not used     ! with no assumption as to the type of map (gnomonic, equi-angular, parametric)
! Subprogram not used     !
! Subprogram not used     ! f = ref2sphere(xvec) - sphere
! Subprogram not used     ! df = d(ref2sphere) 
! Subprogram not used     !
! Subprogram not used     ! D = diag(cos(theta),1) * d(ref2sphere)       d(ref2sphere) = diag(1/cos(theta),1)*D
! Subprogram not used     ! 
! Subprogram not used     ! df = diag(1/cos(theta),1)*D 
! Subprogram not used     ! df^-1 =  D^-1 *  diag(cos(theta),1)
! Subprogram not used     ! 
! Subprogram not used     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used     costh=cos(sphere%lat)
! Subprogram not used     a=0
! Subprogram not used     b=0
! Subprogram not used     i=0
! Subprogram not used     do
! Subprogram not used        sphere1 = ref2sphere(a,b,elem)
! Subprogram not used        resa = sphere1%lon - sphere%lon
! Subprogram not used        if (resa>dd_pi) resa=resa-2*dd_pi
! Subprogram not used        if (resa<-dd_pi) resa=resa+2*dd_pi
! Subprogram not used 
! Subprogram not used        resb = sphere1%lat - sphere%lat 
! Subprogram not used 
! Subprogram not used        call Dmap(D,elem,a,b)
! Subprogram not used        detD = D(1,1)*D(2,2) - D(1,2)*D(2,1)      
! Subprogram not used        Dinv(1,1) =  D(2,2)/detD
! Subprogram not used        Dinv(1,2) = -D(1,2)/detD
! Subprogram not used        Dinv(2,1) = -D(2,1)/detD
! Subprogram not used        Dinv(2,2) =  D(1,1)/detD
! Subprogram not used        
! Subprogram not used        dela =  Dinv(1,1)*costh*resa + Dinv(1,2)*resb 
! Subprogram not used        delb =  Dinv(2,1)*costh*resa + Dinv(2,2)*resb 
! Subprogram not used        a = a - dela
! Subprogram not used        b = b - delb
! Subprogram not used        i=i+1
! Subprogram not used        if ( (costh*resa)**2 + resb**2 < tol_sq .or. MAX_NR_ITER < i) exit
! Subprogram not used     end do
! Subprogram not used     ref%x=a
! Subprogram not used     ref%y=b
! Subprogram not used   end function parametric_coordinates




!
! find element containing given point, useing HOMME's standard
! equi-angular gnomonic map.
! note that with this map, only coordinate lines are great circle arcs
!
! Subprogram not used   function point_inside_equiangular(elem, sphere, sphere_xyz) result(inside)
! Subprogram not used     implicit none
! Subprogram not used     type (spherical_polar_t), intent(in)     :: sphere
! Subprogram not used     type (cartesian3D_t),     intent(in)    :: sphere_xyz
! Subprogram not used     type (element_t)        , intent(in)     :: elem
! Subprogram not used     logical                              :: inside, inside2
! Subprogram not used     integer               :: i,j
! Subprogram not used     type (cartesian2D_t) :: corners(4),sphere_xy,cart
! Subprogram not used     type (cartesian3D_t) :: corners_xyz(4),center,a,b,cross(4)
! Subprogram not used     real (kind=real_kind) :: yp(4), y, elem_diam,dotprod
! Subprogram not used     real (kind=real_kind) :: xp(4), x, xc,yc
! Subprogram not used     real (kind=real_kind) :: tol_inside
! Subprogram not used     real (kind=real_kind) :: d1,d2
! Subprogram not used 
! Subprogram not used     type (spherical_polar_t)    :: sphere_tmp
! Subprogram not used 
! Subprogram not used     inside = .false.
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     ! first check if point is near the element:
! Subprogram not used     corners_xyz(:) = elem%corners3D(:)
! Subprogram not used     elem_diam = max( distance(corners_xyz(1),corners_xyz(3)), &
! Subprogram not used          distance(corners_xyz(2),corners_xyz(4)) )
! Subprogram not used 
! Subprogram not used     center%x = sum(corners_xyz(1:4)%x)/4
! Subprogram not used     center%y = sum(corners_xyz(1:4)%y)/4
! Subprogram not used     center%z = sum(corners_xyz(1:4)%z)/4
! Subprogram not used     if ( distance(center,sphere_xyz) > 1.0*elem_diam ) return
! Subprogram not used 
! Subprogram not used     tol_inside = 1e-10*elem_diam**2
! Subprogram not used     ! the point is close to the element, so project both to cubed sphere
! Subprogram not used     ! and perform contour integral
! Subprogram not used     sphere_xy=sphere2cubedsphere(sphere,elem%FaceNum)
! Subprogram not used     x = sphere_xy%x
! Subprogram not used     y = sphere_xy%y
! Subprogram not used     do i=1,4
! Subprogram not used       xp(i) = elem%corners(i)%x
! Subprogram not used       yp(i) = elem%corners(i)%y
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     if (debug) then
! Subprogram not used        print *,'point: ',x,y,elem%FaceNum
! Subprogram not used        print *,'element:'
! Subprogram not used        write(*,'(a,4e16.8,a)') 'x=[',xp(1:4),']'
! Subprogram not used        write(*,'(a,4e16.8,a)') 'y=[',yp(1:4),']'
! Subprogram not used 
! Subprogram not used        ! first check if centroid is in this element (sanity check)
! Subprogram not used        sphere_tmp=change_coordinates(center)
! Subprogram not used        sphere_xy=sphere2cubedsphere(sphere_tmp,elem%FaceNum)
! Subprogram not used        xc=sphere_xy%x
! Subprogram not used        yc=sphere_xy%y
! Subprogram not used        print *,'cross product with centroid: all numbers should be negative'
! Subprogram not used        j = 4
! Subprogram not used        do i=1,4
! Subprogram not used           print *,i,(xc-xp(j))*(yp(i)-yp(j))  - (yc-yp(j))*(xp(i)-xp(j))
! Subprogram not used           j = i  ! within this loopk j = i-1
! Subprogram not used        end do
! Subprogram not used 
! Subprogram not used        print *,'cross product with search point'
! Subprogram not used        j = 4
! Subprogram not used        do i=1,4
! Subprogram not used           print *,i,(x-xp(j))*(yp(i)-yp(j))  - (y-yp(j))*(xp(i)-xp(j))
! Subprogram not used           j = i  ! within this loopk j = i-1
! Subprogram not used        end do
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     
! Subprogram not used     j = 4
! Subprogram not used     do i=1,4
! Subprogram not used       ! a = x-xp(j), y-yp(j)
! Subprogram not used       ! b = xp(i)-xp(j), yp(i)-yp(j)
! Subprogram not used       ! compute a cross b:
! Subprogram not used       if ( -( (x-xp(j))*(yp(i)-yp(j))  - (y-yp(j))*(xp(i)-xp(j))) > tol_inside ) then
! Subprogram not used          return
! Subprogram not used       endif
! Subprogram not used       j = i  ! within this loopk j = i-1
! Subprogram not used     end do
! Subprogram not used     ! all cross products were negative, must be inside:
! Subprogram not used     inside=.true.
! Subprogram not used   end function point_inside_equiangular


!
! find element containing given point, with element edges assumed to be great circle arcs
! this will work with any map where straight lines are mapped to great circle arcs.
! (thus it will fail on unstructured grids using the equi-angular gnomonic map)
!
! Subprogram not used   function point_inside_gc(elem, sphere_xyz) result(inside)
! Subprogram not used     implicit none
! Subprogram not used     type (cartesian3D_t),     intent(in)    :: sphere_xyz
! Subprogram not used     type (element_t)        , intent(in)     :: elem
! Subprogram not used     logical                              :: inside, inside2
! Subprogram not used     integer               :: i,j,ii
! Subprogram not used     type (cartesian2D_t) :: corners(4),sphere_xy,cart
! Subprogram not used     type (cartesian3D_t) :: corners_xyz(4),center,a,b,cross(4)
! Subprogram not used     real (kind=real_kind) :: yp(4), y, elem_diam,dotprod
! Subprogram not used     real (kind=real_kind) :: xp(4), x
! Subprogram not used     real (kind=real_kind) :: d1,d2, tol_inside = 1e-12
! Subprogram not used 
! Subprogram not used     type (spherical_polar_t)   :: sphere  ! debug
! Subprogram not used 
! Subprogram not used     inside = .false.
! Subprogram not used 
! Subprogram not used     ! first check if point is near the element:
! Subprogram not used     corners_xyz(:) = elem%corners3D(:)
! Subprogram not used     elem_diam = max( distance(corners_xyz(1),corners_xyz(3)), &
! Subprogram not used          distance(corners_xyz(2),corners_xyz(4)) )
! Subprogram not used 
! Subprogram not used     center%x = sum(corners_xyz(1:4)%x)/4
! Subprogram not used     center%y = sum(corners_xyz(1:4)%y)/4
! Subprogram not used     center%z = sum(corners_xyz(1:4)%z)/4
! Subprogram not used     if ( distance(center,sphere_xyz) > 1.0*elem_diam ) return
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     j = 4
! Subprogram not used     do i=1,4
! Subprogram not used       ! outward normal to plane containing j->i edge:  corner(i) x corner(j)
! Subprogram not used       ! sphere dot (corner(i) x corner(j) ) = negative if inside
! Subprogram not used        cross(i)%x =  corners_xyz(i)%y*corners_xyz(j)%z - corners_xyz(i)%z*corners_xyz(j)%y
! Subprogram not used        cross(i)%y =-(corners_xyz(i)%x*corners_xyz(j)%z - corners_xyz(i)%z*corners_xyz(j)%x)
! Subprogram not used        cross(i)%z =  corners_xyz(i)%x*corners_xyz(j)%y - corners_xyz(i)%y*corners_xyz(j)%x
! Subprogram not used        dotprod = cross(i)%x*sphere_xyz%x + cross(i)%y*sphere_xyz%y +&
! Subprogram not used                cross(i)%z*sphere_xyz%z
! Subprogram not used        j = i  ! within this loopk j = i-1
! Subprogram not used 
! Subprogram not used        !if (dotprod>0 .and. dotprod/elem_diam < 1e-5) print *,dotprod/elem_diam
! Subprogram not used 
! Subprogram not used        ! dot product is proportional to elem_diam. positive means outside,
! Subprogram not used        ! but allow machine precision tolorence: 
! Subprogram not used        if (dotprod > tol_inside*elem_diam) return 
! Subprogram not used        !if (dotprod > 0) return 
! Subprogram not used     end do
! Subprogram not used     inside=.true.
! Subprogram not used     return
! Subprogram not used   end function point_inside_gc



  !================================================
  !  (Nair) Cube face index and local coordinates
  !================================================

! Subprogram not used   subroutine cube_facepoint_ne(sphere,ne,cart, number)
! Subprogram not used     use coordinate_systems_mod, only : cube_face_number_from_sphere, sphere2cubedsphere
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     type (spherical_polar_t), intent (in) :: sphere
! Subprogram not used     integer             , intent(in)      :: ne
! Subprogram not used     type (cartesian2D_t), intent(out)     :: cart
! Subprogram not used     integer             , intent(out)     :: number
! Subprogram not used 
! Subprogram not used     real (kind=real_kind) :: xp,yp
! Subprogram not used     type (cartesian2D_t)  :: cube
! Subprogram not used     integer               :: ie, je, face_no
! Subprogram not used     real (kind=real_kind) :: x1,x2
! Subprogram not used     real (kind=real_kind) :: dx
! Subprogram not used 
! Subprogram not used     face_no = cube_face_number_from_sphere(sphere)
! Subprogram not used     cube    = sphere2cubedsphere(sphere, face_no)
! Subprogram not used     xp      = cube%x
! Subprogram not used     yp      = cube%y
! Subprogram not used 
! Subprogram not used     ! MNL: for uniform grids (on cube face), analytic solution is fine
! Subprogram not used     x1 = xp + 0.25D0*DD_PI
! Subprogram not used     x2 = yp + 0.25D0*DD_PI
! Subprogram not used 
! Subprogram not used     dx = (0.5D0*DD_PI)/ne
! Subprogram not used     ie = INT(ABS(x1)/dx)
! Subprogram not used     je = INT(ABS(x2)/dx)
! Subprogram not used     ! if we are exactly on an element edge, we can put the point in
! Subprogram not used     ! either the ie or ie+1 element, EXCEPT if ie==ne.
! Subprogram not used     if ( ABS(x1) < ne*dx  ) ie=ie+1
! Subprogram not used     if ( ABS(x2) < ne*dx  ) je=je+1
! Subprogram not used     if (ie>ne .or. je>ne) then
! Subprogram not used        write(iulog,*)'ERROR: ',ie,je,ne
! Subprogram not used        write(iulog,*)'lat,lon=',sphere%lat,sphere%lon
! Subprogram not used        write(iulog,*)'face no=',face_no
! Subprogram not used        write(iulog,*)x1,x2,x1/dx,x2/dx
! Subprogram not used        call abortmp('interpolate_mod: bad argument')
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     ! bug fix MT 1/2009.  This was creating a plotting error at
! Subprogram not used     ! the row of elements in iface=2 at 50 degrees (NE=16 128x256 lat/lon grid)
! Subprogram not used     ! For point on element edge, we can have ie=2, but x1=dx
! Subprogram not used     ! but if ie>1, we must execute this statement.
! Subprogram not used     ! The only time we can skip this statement is if ie=1, but then
! Subprogram not used     ! the statement has no effect, so lets never skip it:
! Subprogram not used     !    if (x1 > dx ) then
! Subprogram not used     x1 = x1 - dble(ie-1)*dx
! Subprogram not used     !    endif
! Subprogram not used 
! Subprogram not used     x1 = 2.0D0*(x1/dx)-1.0D0
! Subprogram not used 
! Subprogram not used     !    if (x2 > dx ) then    ! removed MT 1/2009, see above
! Subprogram not used     x2 = x2 - dble(je-1)*dx
! Subprogram not used     !    endif
! Subprogram not used 
! Subprogram not used     x2 = 2.0D0*(x2/dx)-1.0D0
! Subprogram not used 
! Subprogram not used     ! coordinates within an element [-1,1]
! Subprogram not used     cart%x = x1
! Subprogram not used     cart%y = x2
! Subprogram not used     number = ie + (je-1)*ne + (face_no-1)*ne*ne
! Subprogram not used   end subroutine cube_facepoint_ne
  !================================================
  !  (Nair) Cube face index and local coordinates
  !================================================


! Subprogram not used   subroutine cube_facepoint_unstructured(sphere,cart, number, elem)
! Subprogram not used     use coordinate_systems_mod, only : cube_face_number_from_sphere, &
! Subprogram not used                                        sphere2cubedsphere,change_coordinates,cube_face_number_from_cart
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     type (element_t)     , intent(in), target :: elem(:)
! Subprogram not used     type (spherical_polar_t), intent (in) :: sphere
! Subprogram not used     type (cartesian2D_t), intent(out)     :: cart
! Subprogram not used     integer             , intent(out)     :: number
! Subprogram not used 
! Subprogram not used     integer               :: ii
! Subprogram not used     Logical               :: found
! Subprogram not used     type (cartesian3D_t)       :: sphere_xyz
! Subprogram not used     type (cartesian2D_t)  :: cube
! Subprogram not used     sphere_xyz=spherical_to_cart(sphere)
! Subprogram not used 
! Subprogram not used     number=-1
! Subprogram not used !    print *,'WARNING: using GC map'
! Subprogram not used     do ii = 1,nelemd
! Subprogram not used        ! for equiangular gnomonic map:
! Subprogram not used        ! unstructed grid element edges are NOT great circles
! Subprogram not used        if (cubed_sphere_map==0) then
! Subprogram not used           found = point_inside_equiangular(elem(ii), sphere, sphere_xyz)
! Subprogram not used        else 
! Subprogram not used           ! assume element edges are great circle arcs:
! Subprogram not used           found = point_inside_gc(elem(ii), sphere_xyz)
! Subprogram not used        endif
! Subprogram not used 
! Subprogram not used        if (found) then
! Subprogram not used           number = ii
! Subprogram not used           cart = parametric_coordinates(sphere, elem(ii))
! Subprogram not used           exit
! Subprogram not used        end if
! Subprogram not used     end do
! Subprogram not used   end subroutine cube_facepoint_unstructured


! Subprogram not used   subroutine interp_init()
! Subprogram not used     type (quadrature_t)   :: gp
! Subprogram not used 
! Subprogram not used     gp = gausslobatto(np)
! Subprogram not used     call interpolate_create(gp,interp_p)
! Subprogram not used   end subroutine interp_init


! Subprogram not used   subroutine setup_latlon_interp(elem,interpdata,par)
! Subprogram not used     !
! Subprogram not used     ! initialize interpolation data structures to interpolate to a lat-lon grid
! Subprogram not used     !
! Subprogram not used     !
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used     type (element_t)     , intent(in), target :: elem(:)
! Subprogram not used     type (parallel_t)      , intent(in)       :: par
! Subprogram not used     type (interpdata_t)  , intent(out)        :: interpdata(:)
! Subprogram not used 
! Subprogram not used     ! local
! Subprogram not used     integer i,j,ii,count_total,n_interp,count_max
! Subprogram not used     integer ngrid, number, elem_num, plat
! Subprogram not used     integer countx, missing_pts,ierr
! Subprogram not used     integer :: npts_mult_claims,max_claims
! Subprogram not used 
! Subprogram not used     real (kind=real_kind)    ::  dp,latdeg(nlat+1),clat(nlat+1),w(nlat+1),w_staggered(nlat)
! Subprogram not used     real (kind=real_kind)    ::  clat_staggered(nlat),latdeg_st(nlat),err,err2
! Subprogram not used 
! Subprogram not used     type (spherical_polar_t) :: sphere
! Subprogram not used     type (cartesian2D_t)     :: cart
! Subprogram not used     type (cartesian3D_t)     :: sphere_xyz,sphere2_xyz
! Subprogram not used 
! Subprogram not used     type (quadrature_t)       :: gp
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     logical,save                                  :: first_time=.TRUE.
! Subprogram not used 
! Subprogram not used     ! Array to make sure each interp point is on exactly one process
! Subprogram not used     type (cartesian2D_t),allocatable    :: cart_vec(:,:)
! Subprogram not used     integer :: k
! Subprogram not used     integer, allocatable :: global_elem_gid(:,:),local_elem_gid(:,:), local_elem_num(:,:)
! Subprogram not used 
! Subprogram not used     ! these arrays often are too large for stack, so lets make sure
! Subprogram not used     ! they go on the heap:
! Subprogram not used     allocate(local_elem_num(nlat,nlon))
! Subprogram not used     allocate(local_elem_gid(nlat,nlon))
! Subprogram not used     allocate(global_elem_gid(nlat,nlon))
! Subprogram not used     allocate(cart_vec(nlat,nlon))
! Subprogram not used 
! Subprogram not used     if (par%masterproc) then
! Subprogram not used        write(iulog,'(a,i4,a,i4,a)') 'Initializing ',nlat,' x ',nlon,' lat-lon interpolation grid: '
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     do ii=1,nelemd
! Subprogram not used        interpdata(ii)%n_interp=0  ! reset counter
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     if(first_time)then
! Subprogram not used        NULLIFY(lat)
! Subprogram not used        NULLIFY(gweight)
! Subprogram not used        NULLIFY(lon)
! Subprogram not used        first_time=.FALSE.
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     if (associated(lat))then
! Subprogram not used        if(size(lat)>0)deallocate(lat)
! Subprogram not used     endif
! Subprogram not used     if (associated(gweight))then
! Subprogram not used        if(size(gweight)>0)deallocate(gweight)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if (associated(lon))then
! Subprogram not used        if(size(lon)>0)deallocate(lon)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     allocate(lat(nlat))
! Subprogram not used     allocate(gweight(nlat))
! Subprogram not used     allocate(lon(nlon))
! Subprogram not used     call interp_init()
! Subprogram not used     gweight=0
! Subprogram not used     do i=1,nlon
! Subprogram not used        lon(i)=2*dd_pi*(i-1)/nlon
! Subprogram not used     enddo
! Subprogram not used     if (gridtype==1) then
! Subprogram not used        do j=1,nlat
! Subprogram not used           lat(j) = -dd_pi/2 + dd_pi*(j-1)/(nlat-1)
! Subprogram not used        end do
! Subprogram not used        plat=nlat
! Subprogram not used     endif
! Subprogram not used     if (gridtype==2) then
! Subprogram not used        gp=gauss(nlat)
! Subprogram not used        do j=1,nlat
! Subprogram not used           lat(j) = asin(gp%points(j))
! Subprogram not used           gweight(j) = gp%weights(j)
! Subprogram not used        end do
! Subprogram not used     endif
! Subprogram not used     if (gridtype==3) then
! Subprogram not used        do j=1,nlat
! Subprogram not used           lat(j) = -dd_pi/2 + dd_pi*(j-.5d0)/nlat
! Subprogram not used        end do
! Subprogram not used        plat=nlat+1
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if (gridtype==1 .or. gridtype==3) then
! Subprogram not used        ! gridtype=1    plat=nlat    gweight(1:nlat)=w(1:plat)
! Subprogram not used        ! gridtype=3    plat=nlat+1  gweight(1:nlat)=w_staggered(1:plat-1)
! Subprogram not used 
! Subprogram not used        ! L-R dynamics uses a regular latitude distribution (not gausian).
! Subprogram not used        ! The algorithm below is a bastardized version of LSM: map.F.
! Subprogram not used        dp = 180d0/(plat-1)
! Subprogram not used        do j = 1, plat
! Subprogram not used           latdeg(j) = -90d0 + (j-1)*dp
! Subprogram not used           clat(j) = latdeg(j)*dd_pi/180d0
! Subprogram not used        end do
! Subprogram not used 
! Subprogram not used        ! Calculate latitudes for the staggered grid
! Subprogram not used 
! Subprogram not used        do j = 1, plat-1
! Subprogram not used           clat_staggered(j) = (clat(j) + clat(j+1)) / 2
! Subprogram not used           latdeg_st     (j) = clat_staggered(j)*180d0/dd_pi
! Subprogram not used        end do
! Subprogram not used 
! Subprogram not used        ! Weights are defined as cos(phi)*(delta-phi)
! Subprogram not used        ! For a sanity check, the sum of w across all lats should be 2, or 1 across
! Subprogram not used        ! half of the latitudes.
! Subprogram not used 
! Subprogram not used        do j = 2, plat-1
! Subprogram not used           w(j) = sin(clat_staggered(j)) - sin(clat_staggered(j-1))
! Subprogram not used        end do
! Subprogram not used        w(1) = sin(clat_staggered(1)) + 1
! Subprogram not used        w(plat) = w(1)
! Subprogram not used 
! Subprogram not used        ! with nlat=2048, this error was 4e-16
! Subprogram not used        if (abs(sum(w(1:plat)) - 2) > 1e-8) then
! Subprogram not used           write(iulog,*) 'interpolate_mod: w weights do not sum to 2. sum=',sum(w(1:plat))
! Subprogram not used           call abortmp('interpolate_mod: weights do not sum to 2.')
! Subprogram not used        end if
! Subprogram not used 
! Subprogram not used        dp = dd_pi / (plat-1)
! Subprogram not used        do j = 1, plat-1
! Subprogram not used           w_staggered(j) = sin(clat(j+1)) - sin(clat(j))
! Subprogram not used        end do
! Subprogram not used 
! Subprogram not used 
! Subprogram not used        if (abs(sum(w_staggered(1:plat-1)) - 2) > 1e-8) then
! Subprogram not used           write(iulog,*) 'interpolate_mod: staggered weights do not sum to 2. sum=',sum(w_staggered(1:plat-1))
! Subprogram not used           call abortmp('interpolate_mod: weights do not sum to 2.')
! Subprogram not used        end if
! Subprogram not used 
! Subprogram not used        if (gridtype==1) then
! Subprogram not used           gweight(1:nlat)=w(1:plat)
! Subprogram not used        endif
! Subprogram not used        if (gridtype==3) then
! Subprogram not used           gweight(1:nlat)=w_staggered(1:plat-1)
! Subprogram not used        endif
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     ! go through once, counting the number of points on each element
! Subprogram not used     sphere%r=1
! Subprogram not used     local_elem_num  = -1
! Subprogram not used     local_elem_gid  = -1
! Subprogram not used     global_elem_gid = -1
! Subprogram not used     err=0
! Subprogram not used     do j=1,nlat
! Subprogram not used        do i=1,nlon
! Subprogram not used           sphere%lat=lat(j)
! Subprogram not used           sphere%lon=lon(i)
! Subprogram not used 
! Subprogram not used           !debug=.false.
! Subprogram not used           !if (j==18 .and. i==95) debug=.true.
! Subprogram not used           number = -1
! Subprogram not used           if ( (cubed_sphere_map /= 0) .or. MeshUseMeshFile) then
! Subprogram not used              call cube_facepoint_unstructured(sphere, cart, number, elem)
! Subprogram not used              if (number /= -1) then
! Subprogram not used                 ! If points are outside element but within tolerance, move to boundary
! Subprogram not used                 if (cart%x + 1.0d0.le.0.0d0) cart%x = -1.0d0
! Subprogram not used                 if (cart%x - 1.0d0.ge.0.0d0) cart%x = 1.0d0
! Subprogram not used                 if (cart%y + 1.0d0.le.0.0d0) cart%y = -1.0d0
! Subprogram not used                 if (cart%y - 1.0d0.ge.0.0d0) cart%y = 1.0d0
! Subprogram not used 
! Subprogram not used                 local_elem_num(j,i) = number
! Subprogram not used                 local_elem_gid(j,i) = elem(number)%vertex%number
! Subprogram not used                 cart_vec(j,i)    = cart  ! local element coordiante of interpolation point
! Subprogram not used              endif
! Subprogram not used           else
! Subprogram not used              call cube_facepoint_ne(sphere, ne, cart, number)
! Subprogram not used              ! the sphere point belongs to the element number on face = face_no.
! Subprogram not used              ! do I own this element?
! Subprogram not used              if (number /= -1) then
! Subprogram not used                 do ii=1,nelemd
! Subprogram not used                    if (number == elem(ii)%vertex%number) then
! Subprogram not used                       local_elem_gid(j,i) = number
! Subprogram not used                       local_elem_num(j,i) = ii
! Subprogram not used                       cart_vec(j,i)        = cart   ! local element coordinate found above
! Subprogram not used                       exit
! Subprogram not used                    endif
! Subprogram not used                 enddo
! Subprogram not used              endif
! Subprogram not used           endif
! Subprogram not used           ii=local_elem_num(j,i)
! Subprogram not used           if (ii /= -1) then
! Subprogram not used              ! compute error: map 'cart' back to sphere and compare with original
! Subprogram not used              ! interpolation point:
! Subprogram not used              sphere2_xyz = spherical_to_cart( ref2sphere(cart%x,cart%y,elem(ii)))
! Subprogram not used              sphere_xyz = spherical_to_cart(sphere)
! Subprogram not used              err=max(err,distance(sphere2_xyz,sphere_xyz))
! Subprogram not used           endif
! Subprogram not used        enddo
! Subprogram not used        if (par%masterproc) then
! Subprogram not used           if ((MOD(j,64).eq.1).or.(j.eq.nlat)) then
! Subprogram not used              print *,'finished latitude ',j,' of ',nlat
! Subprogram not used           endif
! Subprogram not used        endif
! Subprogram not used     enddo
! Subprogram not used     err2=err
! Subprogram not used     call MPI_Allreduce(err,err2,1,MPIreal_t,MPI_MAX,par%comm,ierr)
! Subprogram not used     if (par%masterproc) then
! Subprogram not used        write(iulog,'(a,e12.4)') 'Max interpolation point search error: ',err2
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     ! if multile elements claim a interpolation point, take the one with largest gid:
! Subprogram not used     global_elem_gid = local_elem_gid
! Subprogram not used     call MPI_Allreduce(local_elem_gid, global_elem_gid, nlat*nlon, MPIinteger_t, MPI_MAX, par%comm,ierr)
! Subprogram not used 
! Subprogram not used     missing_pts=0
! Subprogram not used     do j=1,nlat
! Subprogram not used        do i=1,nlon
! Subprogram not used           if (global_elem_gid(j,i) == -1 ) then
! Subprogram not used              missing_pts = missing_pts + 1
! Subprogram not used              if (par%masterproc) &
! Subprogram not used                   print *,'Error: point not claimed by any element j,i,lat(j),lon(i)=',j,i,lat(j),lon(i)
! Subprogram not used           else if (local_elem_gid(j,i) == global_elem_gid(j,i)  ) then
! Subprogram not used              ii = local_elem_num(j,i)
! Subprogram not used              interpdata(ii)%n_interp = interpdata(ii)%n_interp + 1
! Subprogram not used           endif
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     countx=maxval(interpdata(1:nelemd)%n_interp)
! Subprogram not used     count_max = countx
! Subprogram not used     call MPI_Allreduce(countx,count_max,1,MPIinteger_t,MPI_MAX,par%comm,ierr)
! Subprogram not used 
! Subprogram not used     if (par%masterproc) then
! Subprogram not used        write(iulog,'(a,i6)') 'Maximum number of interpolation points claimed by an element: ',count_max
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     ! allocate storage
! Subprogram not used     do ii=1,nelemd
! Subprogram not used        ngrid = interpdata(ii)%n_interp
! Subprogram not used        if(interpdata(ii)%first_entry)then
! Subprogram not used           NULLIFY(interpdata(ii)%interp_xy)
! Subprogram not used           NULLIFY(interpdata(ii)%ilat)
! Subprogram not used           NULLIFY(interpdata(ii)%ilon)
! Subprogram not used 
! Subprogram not used           interpdata(ii)%first_entry=.FALSE.
! Subprogram not used        endif
! Subprogram not used        if(associated(interpdata(ii)%interp_xy))then
! Subprogram not used           if(size(interpdata(ii)%interp_xy)>0)deallocate(interpdata(ii)%interp_xy)
! Subprogram not used        endif
! Subprogram not used        if(associated(interpdata(ii)%ilat))then
! Subprogram not used           if(size(interpdata(ii)%ilat)>0)deallocate(interpdata(ii)%ilat)
! Subprogram not used        endif
! Subprogram not used 
! Subprogram not used        if (associated(interpdata(ii)%ilon))then
! Subprogram not used           if(size(interpdata(ii)%ilon)>0)deallocate(interpdata(ii)%ilon)
! Subprogram not used        endif
! Subprogram not used        allocate(interpdata(ii)%interp_xy( ngrid ) )
! Subprogram not used        allocate(interpdata(ii)%ilat( ngrid ) )
! Subprogram not used        allocate(interpdata(ii)%ilon( ngrid ) )
! Subprogram not used        interpdata(ii)%n_interp=0  ! reset counter
! Subprogram not used     enddo
! Subprogram not used     do j=1,nlat
! Subprogram not used        do i=1,nlon
! Subprogram not used           if (local_elem_gid(j,i) == global_elem_gid(j,i) .and. &
! Subprogram not used                local_elem_gid(j,i) /= -1 ) then
! Subprogram not used              ii = local_elem_num(j,i)
! Subprogram not used              ngrid = interpdata(ii)%n_interp + 1
! Subprogram not used              interpdata(ii)%n_interp = ngrid
! Subprogram not used              interpdata(ii)%interp_xy( ngrid )   = cart_vec(j,i)
! Subprogram not used              interpdata(ii)%ilon( ngrid ) = i
! Subprogram not used              interpdata(ii)%ilat( ngrid ) = j
! Subprogram not used           endif
! Subprogram not used        enddo
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     ! now lets compute the number of points that were claimed by
! Subprogram not used     ! more than one element:
! Subprogram not used     do j=1,nlat
! Subprogram not used        do i=1,nlon
! Subprogram not used           if (local_elem_gid(j,i) == -1) then
! Subprogram not used              local_elem_gid(j,i)=0
! Subprogram not used           else
! Subprogram not used              local_elem_gid(j,i)=1
! Subprogram not used           endif
! Subprogram not used        enddo
! Subprogram not used     enddo
! Subprogram not used     global_elem_gid = local_elem_gid
! Subprogram not used     call MPI_Allreduce(local_elem_gid, global_elem_gid, nlat*nlon, MPIinteger_t, MPI_SUM, par%comm,ierr)
! Subprogram not used     if (par%masterproc) then
! Subprogram not used        countx=0
! Subprogram not used        do j=1,nlat
! Subprogram not used           do i=1,nlon
! Subprogram not used              if (global_elem_gid(j,i)>1) countx=countx+1
! Subprogram not used           enddo
! Subprogram not used        enddo
! Subprogram not used        npts_mult_claims=countx
! Subprogram not used        max_claims=maxval(global_elem_gid)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if (par%masterproc) then
! Subprogram not used        print *,'Number of interpolation points claimed by more than one element: ',npts_mult_claims
! Subprogram not used        print *,'max number of elements which claimed the same interpolation point:',max_claims
! Subprogram not used     endif
! Subprogram not used     
! Subprogram not used     deallocate(global_elem_gid)
! Subprogram not used     deallocate(local_elem_num)
! Subprogram not used     deallocate(local_elem_gid)
! Subprogram not used     deallocate(cart_vec)
! Subprogram not used 
! Subprogram not used     ! check if every point in interpolation grid was claimed by an element:
! Subprogram not used     if (missing_pts>0) then
! Subprogram not used        count_total = nlat*nlon
! Subprogram not used        if(par%masterproc) then
! Subprogram not used           write(iulog,"(3A,I4,A,I7,a,i5)")"Error:","<stdin>"," ",1393," count_total:",count_total," missing:",missing_pts
! Subprogram not used        end if
! Subprogram not used        call syncmp(par)
! Subprogram not used        call abortmp('Error: interpolation points not claimed by any element')
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   end subroutine setup_latlon_interp



! interpolate_scalar
!
! Interpolate a scalar field given in an element (fld_cube) to the points in
! interpdata%interp_xy(i), i=1 .. interpdata%n_interp.
!
! Note that it is possible the given element contains none of the interpolation points
! =======================================
! Subprogram not used subroutine interpolate_ce(cart,fld_cube,npts,fld, fillvalue)
! Subprogram not used   type (cartesian2D_t) :: cart
! Subprogram not used   integer                  ::  npts
! Subprogram not used   real (kind=real_kind)    ::  fld_cube(npts,npts) ! cube field
! Subprogram not used   real (kind=real_kind)    ::  fld          ! field at new grid lat,lon coordinates
! Subprogram not used   real (kind=real_kind), intent(in), optional :: fillvalue
! Subprogram not used   ! Local variables
! Subprogram not used   type (interpolate_t), pointer  ::  interp          ! interpolation structure
! Subprogram not used 
! Subprogram not used   integer :: ne
! Subprogram not used 
! Subprogram not used   integer :: i
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   if (npts==np) then
! Subprogram not used      interp => interp_p
! Subprogram not used   else
! Subprogram not used      call abortmp('Error in interpolate_scalar(): must be called with p or v grid data')
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   fld=interpolate_2d(cart,fld_cube,interp,npts,fillvalue)
! Subprogram not used 
! Subprogram not used end subroutine interpolate_ce



  ! =======================================
  ! interpolate_scalar
  !
  ! Interpolate a scalar field given in an element (fld_cube) to the points in
  ! interpdata%interp_xy(i), i=1 .. interpdata%n_interp.
  !
  ! Note that it is possible the given element contains none of the interpolation points
  ! =======================================
! Subprogram not used   subroutine interpolate_scalar2d(interpdata,fld_cube,npts,fld, fillvalue)
! Subprogram not used     integer                  ::  npts
! Subprogram not used     real (kind=real_kind)    ::  fld_cube(npts,npts) ! cube field
! Subprogram not used     real (kind=real_kind)    ::  fld(:)          ! field at new grid lat,lon coordinates
! Subprogram not used     type (interpdata_t)         ::  interpdata
! Subprogram not used     real (kind=real_kind), intent(in), optional :: fillvalue
! Subprogram not used     ! Local variables
! Subprogram not used     type (interpolate_t), pointer  ::  interp          ! interpolation structure
! Subprogram not used 
! Subprogram not used     integer :: ne
! Subprogram not used 
! Subprogram not used     integer :: i
! Subprogram not used 
! Subprogram not used     type (cartesian2D_t) :: cart
! Subprogram not used 
! Subprogram not used     if (npts==np) then
! Subprogram not used        interp => interp_p
! Subprogram not used     else
! Subprogram not used        call abortmp('Error in interpolate_scalar(): must be called with p or v grid data')
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used        ! Choice for Native (high-order) or Bilinear interpolations
! Subprogram not used     if(present(fillvalue)) then
! Subprogram not used        if (itype == 0) then
! Subprogram not used           do i=1,interpdata%n_interp
! Subprogram not used              fld(i)=interpolate_2d(interpdata%interp_xy(i),fld_cube,interp,npts,fillvalue)
! Subprogram not used           end do
! Subprogram not used        elseif (itype == 1) then
! Subprogram not used           do i=1,interpdata%n_interp
! Subprogram not used              fld(i)=interpol_bilinear(interpdata%interp_xy(i),fld_cube,interp,npts,fillvalue)
! Subprogram not used           end do
! Subprogram not used        end if
! Subprogram not used     else
! Subprogram not used        if (itype == 0) then
! Subprogram not used           do i=1,interpdata%n_interp
! Subprogram not used              fld(i)=interpolate_2d(interpdata%interp_xy(i),fld_cube,interp,npts)
! Subprogram not used           end do
! Subprogram not used        elseif (itype == 1) then
! Subprogram not used           do i=1,interpdata%n_interp
! Subprogram not used              fld(i)=interpol_bilinear(interpdata%interp_xy(i),fld_cube,interp,npts)
! Subprogram not used           end do
! Subprogram not used        end if
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   end subroutine interpolate_scalar2d
! Subprogram not used   subroutine interpolate_scalar3d(interpdata,fld_cube,npts,nlev,fld, fillvalue)
! Subprogram not used     integer , intent(in)                 ::  npts, nlev
! Subprogram not used     real (kind=real_kind)    ::  fld_cube(npts,npts,nlev) ! cube field
! Subprogram not used     real (kind=real_kind)    ::  fld(:,:)          ! field at new grid lat,lon coordinates
! Subprogram not used     type (interpdata_t)         ::  interpdata
! Subprogram not used     real (kind=real_kind), intent(in), optional :: fillvalue
! Subprogram not used     ! Local variables
! Subprogram not used     type (interpolate_t), pointer  ::  interp          ! interpolation structure
! Subprogram not used 
! Subprogram not used     integer :: ne
! Subprogram not used 
! Subprogram not used     integer :: i, k
! Subprogram not used 
! Subprogram not used     type (cartesian2D_t) :: cart
! Subprogram not used 
! Subprogram not used     if (npts==np) then
! Subprogram not used        interp => interp_p
! Subprogram not used     else
! Subprogram not used        call abortmp('Error in interpolate_scalar(): must be called with p or v grid data')
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     ! Choice for Native (high-order) or Bilinear interpolations
! Subprogram not used     if(present(fillvalue)) then
! Subprogram not used        if (itype == 0) then
! Subprogram not used           do k=1,nlev
! Subprogram not used              do i=1,interpdata%n_interp
! Subprogram not used                 fld(i,k)=interpolate_2d(interpdata%interp_xy(i),fld_cube(:,:,k),interp,npts,fillvalue)
! Subprogram not used              end do
! Subprogram not used           end do
! Subprogram not used        elseif (itype == 1) then
! Subprogram not used           do k=1,nlev
! Subprogram not used              do i=1,interpdata%n_interp
! Subprogram not used                 fld(i,k)=interpol_bilinear(interpdata%interp_xy(i),fld_cube(:,:,k),interp,npts,fillvalue)
! Subprogram not used              end do
! Subprogram not used           end do
! Subprogram not used        endif
! Subprogram not used     else
! Subprogram not used        if (itype == 0) then
! Subprogram not used           do k=1,nlev
! Subprogram not used              do i=1,interpdata%n_interp
! Subprogram not used                 fld(i,k)=interpolate_2d(interpdata%interp_xy(i),fld_cube(:,:,k),interp,npts)
! Subprogram not used              end do
! Subprogram not used           end do
! Subprogram not used        elseif (itype == 1) then
! Subprogram not used           do k=1,nlev
! Subprogram not used              do i=1,interpdata%n_interp
! Subprogram not used                 fld(i,k)=interpol_bilinear(interpdata%interp_xy(i),fld_cube(:,:,k),interp,npts)
! Subprogram not used              end do
! Subprogram not used           end do
! Subprogram not used        else
! Subprogram not used           write(iulog,*) itype
! Subprogram not used           call abortmp("wrong interpolation type")
! Subprogram not used        endif
! Subprogram not used     endif
! Subprogram not used   end subroutine interpolate_scalar3d



  ! =======================================
  ! interpolate_vector
  !
  ! Interpolate a vector field given in an element (fld_cube)
  ! to the points in interpdata%interp_xy(i), i=1 .. interpdata%n_interp.
  !
  ! input_coords = 0    fld_cube given in lat-lon
  ! input_coords = 1    fld_cube given in contravariant
  !
  ! Note that it is possible the given element contains none of the interpolation points
  ! =======================================
! Subprogram not used   subroutine interpolate_vector2d(interpdata,elem,fld_cube,npts,fld,input_coords, fillvalue)
! Subprogram not used     implicit none
! Subprogram not used     integer                  ::  npts
! Subprogram not used     real (kind=real_kind)    ::  fld_cube(npts,npts,2) ! vector field
! Subprogram not used     real (kind=real_kind)    ::  fld(:,:)          ! field at new grid lat,lon coordinates
! Subprogram not used     type (interpdata_t)      ::  interpdata
! Subprogram not used     type (element_t), intent(in)         :: elem
! Subprogram not used     real (kind=real_kind), intent(in), optional :: fillvalue
! Subprogram not used     integer                  ::  input_coords
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     ! Local variables
! Subprogram not used     real (kind=real_kind)    ::  fld_contra(npts,npts,2) ! vector field
! Subprogram not used     type (interpolate_t), pointer  ::  interp          ! interpolation structure
! Subprogram not used 
! Subprogram not used     real (kind=real_kind)    ::  v1,v2
! Subprogram not used     real (kind=real_kind)    ::  D(2,2)   ! derivative of gnomonic mapping
! Subprogram not used     real (kind=real_kind)    ::  JJ(2,2), tmpD(2,2)   ! derivative of gnomonic mapping
! Subprogram not used 
! Subprogram not used     integer :: i,j
! Subprogram not used 
! Subprogram not used     type (cartesian2D_t) :: cart
! Subprogram not used 
! Subprogram not used     if(present(fillvalue)) then
! Subprogram not used        if (any(fld_cube==fillvalue)) then
! Subprogram not used           fld = fillvalue
! Subprogram not used           return
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     if (input_coords==0 ) then
! Subprogram not used        ! convert to contra
! Subprogram not used        do j=1,npts
! Subprogram not used           do i=1,npts
! Subprogram not used              ! latlon->contra
! Subprogram not used              fld_contra(i,j,1) = elem%Dinv(1,1,i,j)*fld_cube(i,j,1) + elem%Dinv(1,2,i,j)*fld_cube(i,j,2)
! Subprogram not used              fld_contra(i,j,2) = elem%Dinv(2,1,i,j)*fld_cube(i,j,1) + elem%Dinv(2,2,i,j)*fld_cube(i,j,2)
! Subprogram not used           enddo
! Subprogram not used        enddo
! Subprogram not used     else
! Subprogram not used        fld_contra=fld_cube
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     if (npts==np) then
! Subprogram not used        interp => interp_p
! Subprogram not used     else if (npts==np) then
! Subprogram not used        call abortmp('Error in interpolate_vector(): input must be on velocity grid')
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used        ! Choice for Native (high-order) or Bilinear interpolations
! Subprogram not used 
! Subprogram not used     if (itype == 0) then
! Subprogram not used        do i=1,interpdata%n_interp
! Subprogram not used           fld(i,1)=interpolate_2d(interpdata%interp_xy(i),fld_contra(:,:,1),interp,npts)
! Subprogram not used           fld(i,2)=interpolate_2d(interpdata%interp_xy(i),fld_contra(:,:,2),interp,npts)
! Subprogram not used        end do
! Subprogram not used     elseif (itype == 1) then
! Subprogram not used        do i=1,interpdata%n_interp
! Subprogram not used           fld(i,1)=interpol_bilinear(interpdata%interp_xy(i),fld_contra(:,:,1),interp,npts)
! Subprogram not used           fld(i,2)=interpol_bilinear(interpdata%interp_xy(i),fld_contra(:,:,2),interp,npts)
! Subprogram not used        end do
! Subprogram not used     else
! Subprogram not used        write(iulog,*) itype
! Subprogram not used        call abortmp("wrong interpolation type")
! Subprogram not used     endif
! Subprogram not used     do i=1,interpdata%n_interp
! Subprogram not used        ! convert fld from contra->latlon
! Subprogram not used        call dmap(D,elem,interpdata%interp_xy(i)%x,interpdata%interp_xy(i)%y)
! Subprogram not used        ! convert fld from contra->latlon
! Subprogram not used        v1 = fld(i,1)
! Subprogram not used        v2 = fld(i,2)
! Subprogram not used 
! Subprogram not used        fld(i,1)=D(1,1)*v1 + D(1,2)*v2
! Subprogram not used        fld(i,2)=D(2,1)*v1 + D(2,2)*v2
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used   end subroutine interpolate_vector2d
  ! =======================================
  ! interpolate_vector
  !
  ! Interpolate a vector field given in an element (fld_cube)
  ! to the points in interpdata%interp_xy(i), i=1 .. interpdata%n_interp.
  !
  ! input_coords = 0    fld_cube given in lat-lon
  ! input_coords = 1    fld_cube given in contravariant
  !
  ! Note that it is possible the given element contains none of the interpolation points
  ! =======================================
! Subprogram not used   subroutine interpolate_vector3d(interpdata,elem,fld_cube,npts,nlev,fld,input_coords, fillvalue)
! Subprogram not used     implicit none
! Subprogram not used     type (interpdata_t),intent(in)       ::  interpdata
! Subprogram not used     type (element_t), intent(in)         :: elem
! Subprogram not used     integer, intent(in)                  ::  npts, nlev
! Subprogram not used     real (kind=real_kind), intent(in)    ::  fld_cube(npts,npts,2,nlev) ! vector field
! Subprogram not used     real (kind=real_kind), intent(out)   ::  fld(:,:,:)          ! field at new grid lat,lon coordinates
! Subprogram not used     real (kind=real_kind), intent(in),optional :: fillvalue
! Subprogram not used     integer, intent(in)                  ::  input_coords
! Subprogram not used 
! Subprogram not used     ! Local variables
! Subprogram not used     real (kind=real_kind)    ::  fld_contra(npts,npts,2,nlev) ! vector field
! Subprogram not used     type (interpolate_t), pointer  ::  interp          ! interpolation structure
! Subprogram not used 
! Subprogram not used     real (kind=real_kind)    ::  v1,v2
! Subprogram not used     real (kind=real_kind)    ::  D(2,2)   ! derivative of gnomonic mapping
! Subprogram not used     real (kind=real_kind)    ::  JJ(2,2), tmpD(2,2)   ! derivative of gnomonic mapping
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     integer :: i,j,k
! Subprogram not used 
! Subprogram not used     type (cartesian2D_t) :: cart
! Subprogram not used     if(present(fillvalue)) then
! Subprogram not used        if (any(fld_cube==fillvalue)) then
! Subprogram not used           fld = fillvalue
! Subprogram not used           return
! Subprogram not used        end if
! Subprogram not used     end if
! Subprogram not used     if (input_coords==0 ) then
! Subprogram not used        ! convert to contra
! Subprogram not used        do k=1,nlev
! Subprogram not used           do j=1,npts
! Subprogram not used              do i=1,npts
! Subprogram not used                 ! latlon->contra
! Subprogram not used                 fld_contra(i,j,1,k) = elem%Dinv(1,1,i,j)*fld_cube(i,j,1,k) + elem%Dinv(1,2,i,j)*fld_cube(i,j,2,k)
! Subprogram not used                 fld_contra(i,j,2,k) = elem%Dinv(2,1,i,j)*fld_cube(i,j,1,k) + elem%Dinv(2,2,i,j)*fld_cube(i,j,2,k)
! Subprogram not used              enddo
! Subprogram not used           enddo
! Subprogram not used        end do
! Subprogram not used     else
! Subprogram not used        fld_contra=fld_cube
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if (npts==np) then
! Subprogram not used        interp => interp_p
! Subprogram not used     else if (npts==np) then
! Subprogram not used        call abortmp('Error in interpolate_vector(): input must be on velocity grid')
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used        ! Choice for Native (high-order) or Bilinear interpolations
! Subprogram not used 
! Subprogram not used     if (itype == 0) then
! Subprogram not used        do k=1,nlev
! Subprogram not used           do i=1,interpdata%n_interp
! Subprogram not used              fld(i,k,1)=interpolate_2d(interpdata%interp_xy(i),fld_contra(:,:,1,k),interp,npts)
! Subprogram not used              fld(i,k,2)=interpolate_2d(interpdata%interp_xy(i),fld_contra(:,:,2,k),interp,npts)
! Subprogram not used           end do
! Subprogram not used        end do
! Subprogram not used     elseif (itype == 1) then
! Subprogram not used        do k=1,nlev
! Subprogram not used           do i=1,interpdata%n_interp
! Subprogram not used              fld(i,k,1)=interpol_bilinear(interpdata%interp_xy(i),fld_contra(:,:,1,k),interp,npts)
! Subprogram not used              fld(i,k,2)=interpol_bilinear(interpdata%interp_xy(i),fld_contra(:,:,2,k),interp,npts)
! Subprogram not used           end do
! Subprogram not used        end do
! Subprogram not used     else
! Subprogram not used        call abortmp("wrong interpolation type")
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     do i=1,interpdata%n_interp
! Subprogram not used        ! compute D(:,:) at the point elem%interp_cube(i)
! Subprogram not used        call dmap(D,elem,interpdata%interp_xy(i)%x,interpdata%interp_xy(i)%y)
! Subprogram not used        do k=1,nlev
! Subprogram not used           ! convert fld from contra->latlon
! Subprogram not used           v1 = fld(i,k,1)
! Subprogram not used           v2 = fld(i,k,2)
! Subprogram not used 
! Subprogram not used           fld(i,k,1)=D(1,1)*v1 + D(1,2)*v2
! Subprogram not used           fld(i,k,2)=D(2,1)*v1 + D(2,2)*v2
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used   end subroutine interpolate_vector3d
  function var_is_vector_uvar(name)
    character(len=*), intent(in) :: name
    integer :: i, var_is_vector_uvar, null_index

    var_is_vector_uvar=0
    null_index=0
    do i=1,MAX_VECVARS
       if(trim(vector_uvars(i)).eq. '') then
          null_index=i
          exit
       endif
       if(trim(vector_uvars(i)).eq.name) then
          var_is_vector_uvar=i
          exit
       end if
    end do
  end function var_is_vector_uvar


  function var_is_vector_vvar(name)
    character(len=*), intent(in) :: name
    integer :: i, var_is_vector_vvar, null_index

    var_is_vector_vvar=0
    null_index=0
    do i=1,MAX_VECVARS
       if(trim(vector_vvars(i)).eq. '') then
          null_index=i
          exit
       endif
       if(trim(vector_vvars(i)).eq.name) then
          var_is_vector_vvar=i
          exit
       end if
    end do


  end function var_is_vector_vvar



end module interpolate_mod
