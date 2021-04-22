!=======================================================================
!BOP
!
! !MODULE: ice_transport_remap - horizontal transport via incremental remapping
!
! !DESCRIPTION:
!
! Transports quantities using the second-order conservative remapping
! scheme developed by John Dukowicz and John Baumgardner (DB) and modified
! for sea ice by William Lipscomb and Elizabeth Hunke.
!
! References:
!
! Dukowicz, J. K., and J. R. Baumgardner, 2000: Incremental
!  remapping as a transport/advection algorithm, J. Comput. Phys.,
!  160, 318-335.
!
! Lipscomb, W. H., and E. C. Hunke, 2004: Modeling sea ice
!  transport using incremental remapping, Mon. Wea. Rev., 132,
!  1341-1354.
!
! !REVISION HISTORY:
!  SVN:$Id: ice_transport_remap.F 33 2006-11-13 19:51:14Z eclare $
!
! authors William H. Lipscomb, LANL
!         John Baumgardner, LANL
!
! 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb
! 2004-05: Block structure added (WHL)
! 2006: Moved remap driver to ice_transport_driver
!       Geometry changes: 
!       (1) Reconstruct fields in stretched logically rectangular coordinates
!       (2) Modify geometry so that the area flux across each edge
!           can be specified (following an idea of Mats Bentsen)
!
! !INTERFACE:
!
      module ice_transport_remap
!
! !USES:
!
      use ice_kinds_mod
      use ice_communicate, only: my_task, master_task, MPI_COMM_ICE
      use ice_domain_size
      use ice_constants
      use ice_fileunits, only: nu_diag, flush_fileunit
      use perf_mod,      only: t_startf, t_stopf, t_barrierf
!
!EOP
!
      implicit none
      save
      private
      public :: init_remap, horizontal_remap, make_masks

      logical (kind=log_kind), public :: maskhalo_remap

      integer (kind=int_kind), parameter ::                      &
         max_ntrace = 2+max_ntrcr+nilyr+nslyr  ! hice,hsno,qice,qsno,trcr

      integer (kind=int_kind), parameter ::     &
         ngroups  = 6      ,&! number of groups of triangles that
                             ! contribute transports across each edge
         nvert = 3           ! number of vertices in a triangle

      ! for triangle integral formulas
      real (kind=dbl_kind), parameter ::  & 
         p5625m = -9._dbl_kind/16._dbl_kind    ,&
         p52083 = 25._dbl_kind/48._dbl_kind

      logical (kind=log_kind), parameter :: bugcheck = .false.

!=======================================================================
! Here is some information about how the incremental remapping scheme
! works in CICE and how it can be adapted for use in other models.  
!
! The remapping routine is designed to transport a generic mass-like 
! field (in CICE, the ice fractional area) along with an arbitrary number
! of tracers in two dimensions.  The velocity components are assumed 
! to lie at grid cell corners and the transported scalars at cell centers. 
! Incremental remapping has the following desirable properties: 
! 
! (1) Tracer monotonicity is preserved.  That is, no new local 
!     extrema are produced in fields like ice thickness or internal 
!     energy. 
! (2) The reconstucted mass and tracer fields vary linearly in x and y. 
!     This means that remapping is 2nd-order accurate in space, 
!     except where horizontal gradients are limited to preserve 
!     monotonicity. 
! (3) There are economies of scale.  Transporting a single field 
!     is rather expensive, but additional fields have a relatively 
!     low marginal cost. 
! 
! The following generic conservation equations may be solved: 
! 
!            dm/dt = del*(u*m)             (0) 
!       d(m*T1)/dt = del*(u*m*T1)          (1) 
!    d(m*T1*T2)/dt = del*(u*m*T1*T2)       (2) 
! d(m*T1*T2*T3)/dt = del*(u*m*T1*T2*T3)    (3) 
!
! where d is a partial derivative, del is the 2D divergence operator,
! u is the horizontal velocity, m is the mass density field, and
! T1, T2, and T3 are tracers.
!
! In CICE, these equations have the form
! 
!               da/dt = del*(u*a)          (4)
! dv/dt =   d(a*h)/dt = del*(u*a*h)        (5)
! de/dt = d(a*h*q)/dt = del*(u*a*h*q)      (6)
!            d(aT)/dt = del*(u*a*t)        (7)
! 
! where a = fractional ice area, v = ice/snow volume, h = v/a = thickness, 
! e = ice/snow internal energy (J/m^2), q = e/v = internal energy per 
! unit volume (J/m^3), and T is a tracer.  These equations express 
! conservation of ice area, volume, internal energy, and area-weighted
! tracer, respectively. 
!
! (Note: In CICE, a, v and e are prognostic quantities from which
!  h and q are diagnosed.  The remapping routine works with tracers,
!  which means that h and q must be derived from a, v, and e before
!  calling the remapping routine.)  
!
! Earlier versions of CICE assumed fixed ice and snow density. 
! Beginning with CICE 4.0, the ice and snow density can be variable. 
! In this case, equations (5) and (6) are replaced by 
! 
! dv/dt =        d(a*h)/dt = del*(u*a*h)          (8)  
! dm/dt =    d(a*h*rho)/dt = del*(u*a*h*rho)      (9)
! de/dt = d(a*h*rho*qm)/dt = del*(u*a*h*rho*qm)   (10)
! 
! where rho = density and qm = internal energy per unit mass (J/kg). 
! Eq. (9) expresses mass conservation, which in the variable-density 
! case is no longer equivalent to volume conservation (8). 
!
! Tracers satisfying equations of the form (1) are called "type 1." 
! In CICE the paradigmatic type 1 tracers are hi and hs. 
! 
! Tracers satisfying equations of the form (2) are called "type 2". 
! The paradigmatic type 2 tracers are qi and qs (or rhoi and rhos 
!  in the variable-density case). 
! 
! Tracers satisfying equations of the form (3) are called "type 3."
! The paradigmatic type 3 tracers are qmi and qms in the variable-density
! case.  There are no such tracers in the constant-density case. 
! 
! The fields a, T1, and T2 are reconstructed in each grid cell with 
! 2nd-order accuracy.  T3 is reconstructed with 1st-order accuracy 
! (i.e., it is transported in upwind fashion) in order to avoid 
! additional mathematical complexity. 
! 
! The mass-like field lives in the array "mm" (shorthand for mean 
! mass) and the tracers fields in the array "tm" (mean tracers). 
! In order to transport tracers correctly, the remapping routine 
! needs to know the tracers types and relationships.  This is done 
! as follows: 
! 
! Each field in the "tm" array is assigned an index, 1:max_ntrace. 
! (Note: max_ntrace is not the same as max_ntrcr, the number of tracers 
! in the trcrn state variable array.  For remapping purposes we 
! have additional tracers hi, hs, qi and qs.) 
! For standard CICE with ntrcr = 1, nilyr = 4, and nslyr = 1, the 
! indexing is as follows: 
! 1   = hi 
! 2   = hs 
! 3   = Ts 
! 4-7 = qi 
! 8   = qs 
! 
! The tracer types (1,2,3) are contained in the "tracer_type" array. 
! For standard CICE: 
! 
!     tracer_type = (1 1 1 2 2 2 2 2) 
! 
! Type 2 and type 3 tracers are said to depend on type 1 tracers. 
! For instance, qi depends on hi, which is to say that 
! there is a conservation equation of the form (2) or (6). 
! Thus we define a "depend" array.  For standard CICE: 
! 
!          depend = (0 0 0 1 1 1 1 2) 
! 
! which implies that elements 1-3 (hi, hs, Ts) are type 1, 
! elements 4-7 (qi) depend on element 1 (hi), and element 8 (qs) 
! depends on element 2 (hs). 
!
! We also define a logical array "has_dependents".  In standard CICE: 
! 
!  has_dependents = (T T F F F F F F), 
! 
! which means that only elements 1 and 2 (hi and hs) have dependent 
! tracers. 
! 
! For the variable-density case, things are a bit more complicated. 
! Suppose we have 4 variable-density ice layers and one variable- 
! density snow layer.  Then the indexing is as follows: 
! 1    = hi 
! 2    = hs 
! 3    = Ts 
! 4-7  = rhoi 
! 8    = rhos 
! 9-12 = qmi 
! 13   = qms 
! 
! The key arrays are: 
! 
!    tracer_type = (1 1 1 2 2 2 2 2 3 3 3 3 3) 
! 
!         depend = (0 0 0 1 1 1 1 2 4 5 6 7 8) 
! 
! has_dependents = (T T F T T T T T F F F F F) 
! 
! which imply that hi and hs are type 1 with dependents rhoi and rhos, 
! while rhoi and rhos are type 2 with dependents qmi and qms. 
! 
! Tracers added to the ntrcr array are handled automatically 
! by the remapping with little extra coding.  It is necessary 
! only to provide the correct type and dependency information. 
!
! When using this routine in other models, most of the tracer dependency
! apparatus may be irrelevant.  In a layered ocean model, for example,
! the transported fields are the layer thickness h (the mass density
! field) and two or more tracers (T, S, and various trace species).
! Suppose there are just two tracers, T and S.  Then the tracer arrays
! have the values:
!
!    tracer_type = (1 1)
!         depend = (0 0)
! has_dependents = (F F)
!
! which is to say that all tracer transport equations are of the form (1).
!
! The tracer dependency arrays are optional input arguments for the
! main remapping subroutine.  If these arrays are not passed in, they
! take on the default values tracer_type(:) = 1, depend(:) = 0, and
! has_dependents(:) = F, which are appropriate for most purposes.
!
! Another optional argument is integral_order.  If integral_order = 1,
! then the triangle integrals are exact for linear functions of x and y.
! If integral_order = 2, these integrals are exact for both linear and
! quadratic functions.  If integral_order = 3, integrals are exact for
! cubic functions as well.  If all tracers are of type 1, then the
! integrals of mass*tracer are quadratic, and integral_order = 2 is
! sufficient.  In CICE, where there are type 2 tracers, we integrate
! functions of the form mass*tracer1*tracer2.  Thus integral_order = 3
! is required for exactness, though integral_order = 2 may be good enough
! in practice.
!
! Finally, a few words about the edgearea fields:
!
! In earlier versions of this scheme, the divergence of the velocity
! field implied by the remapping was, in general, different from the
! value of del*u computed in the dynamics.  For energetic consistency
! (in CICE as well as in layered ocean models such as HYPOP),
! these two values should agree.  This can be ensured by setting
! l_fixed_area = T and specifying the area transported across each grid
! cell edge in the arrays edgearea_e and edgearea_n.  The departure
! regions are then tweaked, following an idea by Mats Bentsen, such
! that they have the desired area.  If l_fixed_area = F, these regions
! are not tweaked, and the edgearea arrays are output variables.
!   
!=======================================================================

      contains

!=======================================================================
!
!BOP
!
! !IROUTINE: init_remap - initialize grid quantities used for remapping
!
! !INTERFACE:
!
      subroutine init_remap
!
! !DESCRIPTION:
!
! Grid quantities used by the remapping transport scheme
!
! !REVISION HISTORY:
!
! author William H. Lipscomb, LANL
!
! !USES:
!
      use ice_boundary
      use ice_domain
      use ice_blocks
      use ice_grid, only: dxt, dyt,                      &
                          xav, yav, xxav, xyav, yyav,    &
                          xxxav, xxyav, xyyav, yyyav
      use ice_exit
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) ::     &
 	 i, j, iblk     ! standard indices

      ! Compute grid cell average geometric quantities on the scaled
      ! rectangular grid with dx = 1, dy = 1.
      !
      ! Note: On a rectangular grid, the integral of any odd function
      !       of x or y = 0.

      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            xav(i,j,iblk) = c0
            yav(i,j,iblk) = c0
!!!            These formulas would be used on a rectangular grid
!!!            with dimensions (dxt, dyt):  
!!!            xxav(i,j,iblk) = dxt(i,j,iblk)**2 / c12
!!!            yyav(i,j,iblk) = dyt(i,j,iblk)**2 / c12
            xxav(i,j,iblk) = c1/c12
            yyav(i,j,iblk) = c1/c12
            xyav(i,j,iblk) = c0
            xxxav(i,j,iblk) = c0
            xxyav(i,j,iblk) = c0
            xyyav(i,j,iblk) = c0
            yyyav(i,j,iblk) = c0
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
 
      end subroutine init_remap

!=======================================================================
!BOP
!
! !IROUTINE: horizontal_remap - incremental remapping transport scheme
!
! !INTERFACE:
!
! Subprogram not used       subroutine horizontal_remap (dt,                ntrace,     &
! Subprogram not used                                    uvel,              vvel,       &
! Subprogram not used                                    mm,                tm,         &
! Subprogram not used                                    l_fixed_area,                  &
! Subprogram not used                                    edgearea_e,        edgearea_n, &
! Subprogram not used                                    tracer_type_in,    depend_in,  &
! Subprogram not used                                    has_dependents_in,             &
! Subprogram not used                                    integral_order_in,             &
! Subprogram not used                                    l_dp_midpt_in)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used 
! Subprogram not used ! Solve the transport equations for one timestep using the incremental
! Subprogram not used ! remapping scheme developed by John Dukowicz and John Baumgardner (DB)
! Subprogram not used ! and modified for sea ice by William Lipscomb and Elizabeth Hunke.
! Subprogram not used !
! Subprogram not used ! This scheme preserves monotonicity of ice area and tracers.  That is,
! Subprogram not used ! it does not produce new extrema.  It is second-order accurate in space,
! Subprogram not used ! except where gradients are limited to preserve monotonicity. 
! Subprogram not used !
! Subprogram not used ! This version of the remapping allows the user to specify the areal
! Subprogram not used ! flux across each edge, based on an idea developed by Mats Bentsen.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author William H. Lipscomb, LANL
! Subprogram not used ! 2006: Moved driver (subroutine transport_remap) into separate module. 
! Subprogram not used !       Geometry changes (logically rectangular coordinates, fixed
! Subprogram not used !        area fluxes)
! Subprogram not used !       
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use ice_boundary
! Subprogram not used       use ice_global_reductions
! Subprogram not used       use ice_domain
! Subprogram not used       use ice_blocks
! Subprogram not used       use ice_grid, only: HTE, HTN, dxt, dyt, dxu, dyu,       &
! Subprogram not used                           tarea, tarear, hm,                  &
! Subprogram not used                           xav, yav, xxav, xyav, yyav,         &
! Subprogram not used                           xxxav, xxyav, xyyav, yyyav
! Subprogram not used       use ice_exit
! Subprogram not used       use ice_calendar, only: istep1
! Subprogram not used       use ice_timers
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       real (kind=dbl_kind), intent(in) ::     &
! Subprogram not used          dt      ! time step
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), intent(in) :: &
! Subprogram not used          ntrace       ! number of tracers in use
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), intent(in),       &
! Subprogram not used                 dimension(nx_block,ny_block,max_blocks) ::           &
! Subprogram not used          uvel       ,&! x-component of velocity (m/s)
! Subprogram not used          vvel         ! y-component of velocity (m/s)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), intent(inout),     &
! Subprogram not used          dimension (nx_block,ny_block,0:ncat,max_blocks) ::          &
! Subprogram not used          mm           ! mean mass values in each grid cell
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), intent(inout),     &
! Subprogram not used          dimension (nx_block,ny_block,max_ntrace,ncat,max_blocks) ::     &
! Subprogram not used          tm           ! mean tracer values in each grid cell
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! If l_fixed_area is true, the area of each departure region is
! Subprogram not used     !  computed in advance (e.g., by taking the divergence of the 
! Subprogram not used     !  velocity field and passed to locate_triangles.  The departure 
! Subprogram not used     !  regions are adjusted to obtain the desired area.
! Subprogram not used     ! If false, edgearea is computed in locate_triangles and passed out.
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       logical, intent(in) ::    &
! Subprogram not used          l_fixed_area     ! if true, edgearea_e and edgearea_n are prescribed
! Subprogram not used                           ! if false, edgearea is computed here and passed out
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks),  &
! Subprogram not used          intent(inout) ::                                             &
! Subprogram not used          edgearea_e     ,&! area of departure regions for east edges
! Subprogram not used          edgearea_n       ! area of departure regions for north edges
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension (ntrace), intent(in),     &
! Subprogram not used          optional ::           &
! Subprogram not used          tracer_type_in       ,&! = 1, 2, or 3 (see comments above)
! Subprogram not used          depend_in              ! tracer dependencies (see above)
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind), dimension (ntrace), intent(in),     &
! Subprogram not used          optional ::     &
! Subprogram not used          has_dependents_in      ! true if a tracer has dependent tracers
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), intent(in), optional ::     &
! Subprogram not used          integral_order_in      ! polynomial order for triangle integrals
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind), intent(in), optional ::     &
! Subprogram not used          l_dp_midpt_in          ! if true, find departure points using
! Subprogram not used                                 ! corrected midpoint velocity
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       ! local variables
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension (ntrace) ::     &
! Subprogram not used          tracer_type       ,&! = 1, 2, or 3 (see comments above)
! Subprogram not used          depend              ! tracer dependencies (see above)
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind), dimension (ntrace) ::     &
! Subprogram not used          has_dependents      ! true if a tracer has dependent tracers
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind) ::     &
! Subprogram not used          integral_order      ! polynomial order for triangle integrals
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind) ::     &
! Subprogram not used          l_dp_midpt          ! if true, find departure points using
! Subprogram not used                              ! corrected midpoint velocity
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind) ::     &
! Subprogram not used          i, j           ,&! horizontal indices
! Subprogram not used          iblk           ,&! block indices
! Subprogram not used          ilo,ihi,jlo,jhi,&! beginning and end of physical domain
! Subprogram not used          n,              &! ice category index
! Subprogram not used          m                ! ice tracer index
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension(0:ncat,max_blocks) ::     &
! Subprogram not used          icellsnc         ! number of cells with ice
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind),     &
! Subprogram not used          dimension(nx_block*ny_block,0:ncat,max_blocks) ::     &
! Subprogram not used          indxinc, indxjnc   ! compressed i/j indices
! Subprogram not used 
! Subprogram not used       type (block) ::     &
! Subprogram not used          this_block       ! block information for current block
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) ::     &
! Subprogram not used          dpx            ,&! x coordinates of departure points at cell corners
! Subprogram not used          dpy              ! y coordinates of departure points at cell corners
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension(nx_block,ny_block,0:ncat,max_blocks) :: &
! Subprogram not used          mc             ,&! mass at geometric center of cell
! Subprogram not used          mx, my         ,&! limited derivative of mass wrt x and y
! Subprogram not used          mmask            ! = 1. if mass is present, = 0. otherwise
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind),      &
! Subprogram not used          dimension (nx_block,ny_block,max_ntrace,ncat,max_blocks) ::     &
! Subprogram not used          tc             ,&! tracer values at geometric center of cell
! Subprogram not used          tx, ty         ,&! limited derivative of tracer wrt x and y
! Subprogram not used          tmask            ! = 1. if tracer is present, = 0. otherwise
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,0:ncat,max_blocks) ::     &
! Subprogram not used          mflxe, mflxn     ! mass transports across E and N cell edges
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrace,ncat,max_blocks) ::     &
! Subprogram not used          mtflxe, mtflxn   ! mass*tracer transports across E and N cell edges
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,ngroups,max_blocks) ::     &
! Subprogram not used          triarea          ! area of east-edge departure triangle
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,0:nvert,ngroups,max_blocks) ::  &
! Subprogram not used          xp, yp           ! x and y coordinates of special triangle points
! Subprogram not used                           ! (need 4 points for triangle integrals)
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind),     &
! Subprogram not used          dimension (nx_block,ny_block,ngroups,max_blocks) ::     &
! Subprogram not used          iflux          ,&! i index of cell contributing transport
! Subprogram not used          jflux            ! j index of cell contributing transport
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension(ngroups,max_blocks) ::     &
! Subprogram not used          icellsng         ! number of cells with ice
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind),     &
! Subprogram not used          dimension(nx_block*ny_block,ngroups,max_blocks) ::     &
! Subprogram not used          indxing, indxjng ! compressed i/j indices
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind), dimension(max_blocks) ::     &
! Subprogram not used          l_stop           ! if true, abort the model
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension(max_blocks) ::     &
! Subprogram not used          istop, jstop     ! indices of grid cell where model aborts
! Subprogram not used 
! Subprogram not used       character (len=char_len), dimension(max_blocks) ::   &
! Subprogram not used          edge             ! 'north' or 'east'
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
! Subprogram not used          worka, &
! Subprogram not used          workb, &
! Subprogram not used          workc, &
! Subprogram not used          workd
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,2,max_blocks) ::     &
! Subprogram not used           dpwork          
! Subprogram not used  
! Subprogram not used       real (kind=dbl_kind), dimension(nx_block,ny_block,2,0:ncat,max_blocks) :: &
! Subprogram not used           mwork
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), &
! Subprogram not used          dimension(nx_block,ny_block,max_blocks) :: halomask
! Subprogram not used       type (ice_halo) :: halo_info_tracer
! Subprogram not used     !------------------------------------------------------------------- 
! Subprogram not used 
! Subprogram not used       call t_barrierf('cice_hmap_remap1_BARRIER',MPI_COMM_ICE)
! Subprogram not used       call t_startf  ('cice_hmap_remap1')
! Subprogram not used 
! Subprogram not used       l_stop = .false.
! Subprogram not used       istop = 0
! Subprogram not used       jstop = 0
! Subprogram not used 
! Subprogram not used     !------------------------------------------------------------------- 
! Subprogram not used     ! Initialize various remapping arrays and options
! Subprogram not used     ! These are either passed in as optional arguments or set to the
! Subprogram not used     ! default values.
! Subprogram not used     !------------------------------------------------------------------- 
! Subprogram not used 
! Subprogram not used       if (present(tracer_type_in)) then
! Subprogram not used          tracer_type(:) = tracer_type_in(:)
! Subprogram not used       else
! Subprogram not used          tracer_type(:) = 1
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       if (present(depend_in)) then
! Subprogram not used          depend(:) = depend_in(:)
! Subprogram not used       else
! Subprogram not used          depend(:) = 0
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       if (present(has_dependents_in)) then
! Subprogram not used          has_dependents(:) = has_dependents_in(:)
! Subprogram not used       else
! Subprogram not used          has_dependents(:) = .false.
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       if (present(integral_order_in)) then
! Subprogram not used          integral_order = integral_order_in
! Subprogram not used       else
! Subprogram not used          integral_order = 2   ! quadratic integrals
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       if (present(l_dp_midpt_in)) then
! Subprogram not used          l_dp_midpt = l_dp_midpt_in
! Subprogram not used       else
! Subprogram not used          l_dp_midpt = .false.
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       worka(:,:) = c1
! Subprogram not used       workb(:,:) = c1
! Subprogram not used       workc(:,:) = c1
! Subprogram not used       workd(:,:) = c1
! Subprogram not used 
! Subprogram not used !---!-------------------------------------------------------------------
! Subprogram not used !---! Remap the ice area and associated tracers.
! Subprogram not used !---! Remap the open water area (without tracers).
! Subprogram not used !---!-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,n)
! Subprogram not used       do iblk = 1, nblocks
! Subprogram not used 
! Subprogram not used          this_block = get_block(blocks_ice(iblk),iblk)         
! Subprogram not used          ilo = this_block%ilo
! Subprogram not used          ihi = this_block%ihi
! Subprogram not used          jlo = this_block%jlo
! Subprogram not used          jhi = this_block%jhi
! Subprogram not used 
! Subprogram not used     !------------------------------------------------------------------- 
! Subprogram not used     ! Compute masks and count ice cells.
! Subprogram not used     ! Masks are used to prevent tracer values in cells without ice from
! Subprogram not used     !  being used to compute tracer gradients.
! Subprogram not used     !------------------------------------------------------------------- 
! Subprogram not used 
! Subprogram not used          call t_startf  ('cice_hmap_remap1_masks')
! Subprogram not used          call make_masks (nx_block,           ny_block,             &
! Subprogram not used                           ilo, ihi,           jlo, jhi,             &
! Subprogram not used                           nghost,             ntrace,               &
! Subprogram not used                           has_dependents,     icellsnc (:,iblk),    &
! Subprogram not used                           indxinc(:,:,iblk),  indxjnc  (:,:,iblk),  &
! Subprogram not used                           mm   (:,:,:,iblk),  mmask  (:,:,:,iblk),  &
! Subprogram not used                           tm (:,:,:,:,iblk),  tmask(:,:,:,:,iblk))
! Subprogram not used          call t_stopf  ('cice_hmap_remap1_masks')
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Construct linear fields, limiting gradients to preserve monotonicity.
! Subprogram not used     ! Note: Pass in unit arrays instead of true distances HTE, HTN, etc.
! Subprogram not used     !       The resulting gradients are in scaled coordinates.
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          ! open water
! Subprogram not used 
! Subprogram not used          call t_startf  ('cice_hmap_remap1_cf1')
! Subprogram not used          call construct_fields(nx_block,           ny_block,           &
! Subprogram not used                                ilo, ihi,           jlo, jhi,           &
! Subprogram not used                                nghost,             ntrace,             &
! Subprogram not used                                tracer_type,        depend,             &
! Subprogram not used                                has_dependents,     icellsnc (0,iblk),  &
! Subprogram not used                                indxinc(:,0,iblk),  indxjnc(:,0,iblk),  &
! Subprogram not used !                               HTN    (:,:,iblk),  HTE    (:,:,iblk),  &
! Subprogram not used                                worka       (:,:),  workb       (:,:),  &
! Subprogram not used                                hm     (:,:,iblk),  xav    (:,:,iblk),  &
! Subprogram not used                                yav    (:,:,iblk),  xxav   (:,:,iblk),  &
! Subprogram not used                                xyav   (:,:,iblk),  yyav   (:,:,iblk),  &
! Subprogram not used                                xxxav  (:,:,iblk),  xxyav  (:,:,iblk),  &
! Subprogram not used                                xyyav  (:,:,iblk),  yyyav  (:,:,iblk),  &
! Subprogram not used !                               dxt    (:,:,iblk),  dyt    (:,:,iblk),  &
! Subprogram not used                                workc       (:,:),  workd       (:,:),  &
! Subprogram not used                                mm   (:,:,0,iblk),  mc   (:,:,0,iblk),  &
! Subprogram not used                                mx   (:,:,0,iblk),  my   (:,:,0,iblk),  &
! Subprogram not used                                mmask(:,:,0,iblk) )
! Subprogram not used          call t_stopf  ('cice_hmap_remap1_cf1')
! Subprogram not used 
! Subprogram not used          ! ice categories
! Subprogram not used 
! Subprogram not used          call t_startf  ('cice_hmap_remap1_cf2')
! Subprogram not used          do n = 1, ncat
! Subprogram not used 
! Subprogram not used             call construct_fields(nx_block,             ny_block,           &
! Subprogram not used                                   ilo, ihi,             jlo, jhi,           &
! Subprogram not used                                   nghost,               ntrace,             &
! Subprogram not used                                   tracer_type,          depend,             &
! Subprogram not used                                   has_dependents,       icellsnc (n,iblk),  &
! Subprogram not used                                   indxinc  (:,n,iblk),  indxjnc(:,n,iblk),  &
! Subprogram not used !                                  HTN      (:,:,iblk),  HTE    (:,:,iblk),  &
! Subprogram not used                                   worka         (:,:),  workb       (:,:),  &
! Subprogram not used                                   hm       (:,:,iblk),  xav    (:,:,iblk),  &
! Subprogram not used                                   yav      (:,:,iblk),  xxav   (:,:,iblk),  &
! Subprogram not used                                   xyav     (:,:,iblk),  yyav   (:,:,iblk),  &
! Subprogram not used                                   xxxav    (:,:,iblk),  xxyav  (:,:,iblk),  &
! Subprogram not used                                   xyyav    (:,:,iblk),  yyyav  (:,:,iblk),  &
! Subprogram not used !                                  dxt      (:,:,iblk),  dyt   (:,:,iblk),   &
! Subprogram not used                                   workc         (:,:),  workd       (:,:),  &
! Subprogram not used                                   mm     (:,:,n,iblk),  mc   (:,:,n,iblk),  &
! Subprogram not used                                   mx     (:,:,n,iblk),  my   (:,:,n,iblk),  &
! Subprogram not used                                   mmask  (:,:,n,iblk),                      &
! Subprogram not used                                   tm   (:,:,:,n,iblk),  tc (:,:,:,n,iblk),  &
! Subprogram not used                                   tx   (:,:,:,n,iblk),  ty (:,:,:,n,iblk),  &
! Subprogram not used                                   tmask(:,:,:,n,iblk) )
! Subprogram not used 
! Subprogram not used          enddo                  ! n
! Subprogram not used          call t_stopf  ('cice_hmap_remap1_cf2')
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Given velocity field at cell corners, compute departure points
! Subprogram not used     ! of trajectories.
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          call t_startf  ('cice_hmap_remap1_dp')
! Subprogram not used          call departure_points(nx_block,        ny_block,        &
! Subprogram not used                                ilo, ihi,        jlo, jhi,        &
! Subprogram not used                                nghost,          dt,              &
! Subprogram not used                                uvel(:,:,iblk),  vvel(:,:,iblk),  &
! Subprogram not used                                dxu (:,:,iblk),  dyu (:,:,iblk),  &
! Subprogram not used                                HTN (:,:,iblk),  HTE (:,:,iblk),  &
! Subprogram not used                                dpx (:,:,iblk),  dpy (:,:,iblk),  &
! Subprogram not used                                l_dp_midpt,      l_stop  (iblk),  &
! Subprogram not used                                istop   (iblk),  jstop   (iblk))
! Subprogram not used          call t_stopf  ('cice_hmap_remap1_dp')
! Subprogram not used 
! Subprogram not used          if (l_stop(iblk)) then
! Subprogram not used             this_block = get_block(blocks_ice(iblk),iblk)         
! Subprogram not used             write(nu_diag,*) 'istep1, my_task, iblk =',            &
! Subprogram not used                               istep1, my_task, iblk
! Subprogram not used             write (nu_diag,*) 'Global block:', this_block%block_id
! Subprogram not used             if (istop(iblk) > 0 .and. jstop(iblk) > 0)             &
! Subprogram not used                  write(nu_diag,*) 'Global i and j:',               &
! Subprogram not used                                   this_block%i_glob(istop(iblk)),  &
! Subprogram not used                                   this_block%j_glob(jstop(iblk)) 
! Subprogram not used             call abort_ice('remap transport: bad departure points')
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used       enddo                     ! iblk
! Subprogram not used       !$OMP END PARALLEL DO
! Subprogram not used 
! Subprogram not used      call t_stopf ('cice_hmap_remap1')
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Ghost cell updates
! Subprogram not used     ! If nghost >= 2, these calls are not needed
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       if (nghost==1) then
! Subprogram not used 
! Subprogram not used          call t_barrierf ('cice_hmap_halo1_BARRIER',MPI_COMM_ICE)
! Subprogram not used          call t_startf ('cice_hmap_halo1')
! Subprogram not used          call ice_timer_start(timer_bound)
! Subprogram not used 
! Subprogram not used          ! departure points
! Subprogram not used 
! Subprogram not used          ! load em up
! Subprogram not used          !$OMP PARALLEL DO PRIVATE(iblk,i,j)
! Subprogram not used          do iblk = 1, nblocks
! Subprogram not used             do j = 1, ny_block
! Subprogram not used             do i = 1, nx_block
! Subprogram not used               dpwork(i,j,1,  iblk) = dpx(i,j,iblk)
! Subprogram not used               dpwork(i,j,2,  iblk) = dpy(i,j,iblk)
! Subprogram not used               mwork (i,j,1,:,iblk) = mx (i,j,:,iblk)
! Subprogram not used               mwork (i,j,2,:,iblk) = my (i,j,:,iblk)
! Subprogram not used             enddo
! Subprogram not used             enddo
! Subprogram not used          enddo
! Subprogram not used          !$OMP END PARALLEL DO
! Subprogram not used !jw         call ice_HaloUpdate (dpx,                halo_info, &
! Subprogram not used !jw                              field_loc_NEcorner, field_type_vector)
! Subprogram not used !jw         call ice_HaloUpdate (dpy,                halo_info, &
! Subprogram not used !jw                              field_loc_NEcorner, field_type_vector)
! Subprogram not used          call ice_HaloUpdate (dpwork,             halo_info, &
! Subprogram not used                               field_loc_NEcorner, field_type_vector)
! Subprogram not used 
! Subprogram not used          ! mass field
! Subprogram not used          call ice_HaloUpdate (mc,               halo_info, &
! Subprogram not used                               field_loc_center, field_type_scalar)
! Subprogram not used !jw         call ice_HaloUpdate (mx,               halo_info, &
! Subprogram not used !jw                              field_loc_center, field_type_vector)
! Subprogram not used !jw         call ice_HaloUpdate (my,               halo_info, &
! Subprogram not used !jw                              field_loc_center, field_type_vector)
! Subprogram not used          call ice_HaloUpdate (mwork,            halo_info, &
! Subprogram not used                               field_loc_center, field_type_vector)
! Subprogram not used 
! Subprogram not used         call t_stopf ('cice_hmap_halo1')
! Subprogram not used         call t_barrierf ('cice_hmap_copy1_BARRIER',MPI_COMM_ICE)
! Subprogram not used         call t_startf ('cice_hmap_copy1')
! Subprogram not used 
! Subprogram not used          !$OMP PARALLEL DO PRIVATE(iblk,i,j)
! Subprogram not used          do iblk = 1, nblocks
! Subprogram not used             do j = 1, ny_block
! Subprogram not used             do i = 1, nx_block
! Subprogram not used               dpx(i,j,  iblk) = dpwork(i,j,1,  iblk)
! Subprogram not used               dpy(i,j,  iblk) = dpwork(i,j,2,  iblk)
! Subprogram not used               mx (i,j,:,iblk) = mwork (i,j,1,:,iblk)
! Subprogram not used               my (i,j,:,iblk) = mwork (i,j,2,:,iblk)
! Subprogram not used             enddo
! Subprogram not used             enddo
! Subprogram not used          enddo
! Subprogram not used          !$OMP END PARALLEL DO
! Subprogram not used          call t_stopf ('cice_hmap_copy1')
! Subprogram not used 
! Subprogram not used          call t_barrierf ('cice_hmap_halo2_BARRIER',MPI_COMM_ICE)
! Subprogram not used          call t_startf ('cice_hmap_halo2')
! Subprogram not used 
! Subprogram not used          if (maskhalo_remap) then
! Subprogram not used             call t_startf ('cice_hmap_halo2hm')
! Subprogram not used             halomask = 0
! Subprogram not used             !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,n,m,j,i)
! Subprogram not used             do iblk = 1,nblocks
! Subprogram not used                this_block = get_block(blocks_ice(iblk),iblk)         
! Subprogram not used                ilo = this_block%ilo
! Subprogram not used                ihi = this_block%ihi
! Subprogram not used                jlo = this_block%jlo
! Subprogram not used                jhi = this_block%jhi
! Subprogram not used                do n = 1,ncat
! Subprogram not used                do m = 1,ntrace
! Subprogram not used                do j = jlo, jhi
! Subprogram not used                do i = ilo, ihi
! Subprogram not used                   if (tc(i,j,m,n,iblk) /= c0) halomask(i,j,iblk) = 1
! Subprogram not used                   if (tx(i,j,m,n,iblk) /= c0) halomask(i,j,iblk) = 1
! Subprogram not used                   if (ty(i,j,m,n,iblk) /= c0) halomask(i,j,iblk) = 1
! Subprogram not used                enddo
! Subprogram not used                enddo
! Subprogram not used                enddo
! Subprogram not used                enddo
! Subprogram not used             enddo
! Subprogram not used             !$OMP END PARALLEL DO
! Subprogram not used             call t_stopf ('cice_hmap_halo2hm')
! Subprogram not used 
! Subprogram not used             call t_barrierf ('cice_hmap_halo2hh_BARRIER',MPI_COMM_ICE)
! Subprogram not used             call t_startf ('cice_hmap_halo2hh')
! Subprogram not used             call ice_HaloUpdate(halomask, halo_info, &
! Subprogram not used                                 field_loc_center, field_type_scalar)
! Subprogram not used             call t_stopf ('cice_hmap_halo2hh')
! Subprogram not used 
! Subprogram not used             call t_barrierf ('cice_hmap_halo2hc_BARRIER',MPI_COMM_ICE)
! Subprogram not used             call t_startf ('cice_hmap_halo2hc')
! Subprogram not used             call ice_HaloMask(halo_info_tracer, halo_info, halomask)
! Subprogram not used             call t_stopf ('cice_hmap_halo2hc')
! Subprogram not used 
! Subprogram not used             call t_barrierf ('cice_hmap_halo2hu_BARRIER',MPI_COMM_ICE)
! Subprogram not used             call t_startf ('cice_hmap_halo2hu')
! Subprogram not used 
! Subprogram not used             ! tracer fields 
! Subprogram not used             call ice_HaloUpdate (tc(:,:,1:ntrace,:,:), halo_info_tracer, &
! Subprogram not used                                  field_loc_center, field_type_scalar)
! Subprogram not used             call ice_HaloUpdate (tx(:,:,1:ntrace,:,:), halo_info_tracer, &
! Subprogram not used                                  field_loc_center, field_type_vector)
! Subprogram not used             call ice_HaloUpdate (ty(:,:,1:ntrace,:,:), halo_info_tracer, &
! Subprogram not used                                  field_loc_center, field_type_vector)
! Subprogram not used             call ice_timer_stop(timer_bound)
! Subprogram not used             call t_stopf ('cice_hmap_halo2hu')
! Subprogram not used 
! Subprogram not used             call t_barrierf ('cice_hmap_halo2hd_BARRIER',MPI_COMM_ICE)
! Subprogram not used             call t_startf ('cice_hmap_halo2hd')
! Subprogram not used             call ice_HaloDestroy(halo_info_tracer)
! Subprogram not used             call t_stopf ('cice_hmap_halo2hd')
! Subprogram not used          else
! Subprogram not used             ! tracer fields 
! Subprogram not used             call ice_HaloUpdate (tc(:,:,1:ntrace,:,:), halo_info, &
! Subprogram not used                                  field_loc_center, field_type_scalar)
! Subprogram not used             call ice_HaloUpdate (tx(:,:,1:ntrace,:,:), halo_info, &
! Subprogram not used                                  field_loc_center, field_type_vector)
! Subprogram not used             call ice_HaloUpdate (ty(:,:,1:ntrace,:,:), halo_info, &
! Subprogram not used                                  field_loc_center, field_type_vector)
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used          call ice_timer_stop(timer_bound)
! Subprogram not used          call t_stopf ('cice_hmap_halo2')
! Subprogram not used 
! Subprogram not used       endif  ! nghost
! Subprogram not used 
! Subprogram not used       call t_barrierf ('cice_hmap_remap2_BARRIER',MPI_COMM_ICE)
! Subprogram not used       call t_startf ('cice_hmap_remap2')
! Subprogram not used 
! Subprogram not used       !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,n)
! Subprogram not used       do iblk = 1, nblocks
! Subprogram not used 
! Subprogram not used          this_block = get_block(blocks_ice(iblk),iblk)         
! Subprogram not used          ilo = this_block%ilo
! Subprogram not used          ihi = this_block%ihi
! Subprogram not used          jlo = this_block%jlo
! Subprogram not used          jhi = this_block%jhi
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Transports for east cell edges.
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Compute areas and vertices of departure triangles.
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          call t_startf  ('cice_hmap_remap2e')
! Subprogram not used          edge(iblk) = 'east'
! Subprogram not used          call locate_triangles(nx_block,              ny_block,           &
! Subprogram not used                                ilo, ihi,              jlo, jhi,           &
! Subprogram not used                                nghost,                edge(iblk),         &
! Subprogram not used                                icellsng    (:,iblk),                      &
! Subprogram not used                                indxing   (:,:,iblk),  indxjng(:,:,iblk),  &
! Subprogram not used                                dpx       (:,:,iblk),  dpy    (:,:,iblk),  &
! Subprogram not used                                dxu       (:,:,iblk),  dyu    (:,:,iblk),  &
! Subprogram not used                                xp    (:,:,:,:,iblk),  yp (:,:,:,:,iblk),  &
! Subprogram not used                                iflux   (:,:,:,iblk),  jflux(:,:,:,iblk),  &
! Subprogram not used                                triarea (:,:,:,iblk),  l_fixed_area,       &
! Subprogram not used                                edgearea_e(:,:,iblk))
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Given triangle vertices, compute coordinates of triangle points
! Subprogram not used     !  needed for transport integrals.
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          call triangle_coordinates (nx_block,           ny_block,           &
! Subprogram not used                                     integral_order,     icellsng (:,iblk),  &
! Subprogram not used                                     indxing(:,:,iblk),  indxjng(:,:,iblk),  &
! Subprogram not used                                     xp (:,:,:,:,iblk),  yp (:,:,:,:,iblk))
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Compute the transport across east cell edges by summing contributions
! Subprogram not used     ! from each triangle.
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          ! open water
! Subprogram not used 
! Subprogram not used          call transport_integrals(nx_block,           ny_block,             &
! Subprogram not used                                   ntrace,             icellsng (:,iblk),    &
! Subprogram not used                                   indxing(:,:,iblk),  indxjng  (:,:,iblk),  &
! Subprogram not used                                   tracer_type,        depend,               &
! Subprogram not used                                   integral_order,     triarea(:,:,:,iblk),  &
! Subprogram not used                                   iflux(:,:,:,iblk),  jflux  (:,:,:,iblk),  &
! Subprogram not used                                   xp (:,:,:,:,iblk),  yp   (:,:,:,:,iblk),  &
! Subprogram not used                                   mc   (:,:,0,iblk),  mx     (:,:,0,iblk),  &
! Subprogram not used                                   my   (:,:,0,iblk),  mflxe  (:,:,0,iblk))
! Subprogram not used 
! Subprogram not used          ! ice categories
! Subprogram not used          do n = 1, ncat
! Subprogram not used             call transport_integrals                                       &
! Subprogram not used                                (nx_block,           ny_block,              &
! Subprogram not used                                 ntrace,             icellsng (:,iblk),     &
! Subprogram not used                                 indxing(:,:,iblk),  indxjng   (:,:,iblk),  &
! Subprogram not used                                 tracer_type,        depend,                &
! Subprogram not used                                 integral_order,     triarea (:,:,:,iblk),  &
! Subprogram not used                                 iflux(:,:,:,iblk),  jflux   (:,:,:,iblk),  &
! Subprogram not used                                 xp (:,:,:,:,iblk),  yp    (:,:,:,:,iblk),  &
! Subprogram not used                                 mc   (:,:,n,iblk),  mx      (:,:,n,iblk),  &
! Subprogram not used                                 my   (:,:,n,iblk),  mflxe   (:,:,n,iblk),  &
! Subprogram not used                                 tc (:,:,:,n,iblk),  tx    (:,:,:,n,iblk),  &
! Subprogram not used                                 ty (:,:,:,n,iblk),  mtflxe(:,:,:,n,iblk))
! Subprogram not used 
! Subprogram not used          enddo
! Subprogram not used          call t_stopf  ('cice_hmap_remap2e')
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Repeat for north edges
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          call t_startf  ('cice_hmap_remap2n')
! Subprogram not used          edge(iblk) = 'north'
! Subprogram not used          call locate_triangles(nx_block,              ny_block,           &
! Subprogram not used                                ilo, ihi,              jlo, jhi,           &
! Subprogram not used                                nghost,                edge(iblk),         &
! Subprogram not used                                icellsng    (:,iblk),                      &
! Subprogram not used                                indxing   (:,:,iblk),  indxjng(:,:,iblk),  &
! Subprogram not used                                dpx       (:,:,iblk),  dpy    (:,:,iblk),  &
! Subprogram not used                                dxu       (:,:,iblk),  dyu    (:,:,iblk),  &
! Subprogram not used                                xp    (:,:,:,:,iblk),  yp (:,:,:,:,iblk),  &
! Subprogram not used                                iflux   (:,:,:,iblk),  jflux(:,:,:,iblk),  &
! Subprogram not used                                triarea (:,:,:,iblk),  l_fixed_area,       &
! Subprogram not used                                edgearea_n(:,:,iblk))
! Subprogram not used 
! Subprogram not used          call triangle_coordinates (nx_block,           ny_block,           &
! Subprogram not used                                     integral_order,     icellsng (:,iblk),  &
! Subprogram not used                                     indxing(:,:,iblk),  indxjng(:,:,iblk),  &
! Subprogram not used                                     xp (:,:,:,:,iblk),  yp (:,:,:,:,iblk))
! Subprogram not used 
! Subprogram not used          ! open water
! Subprogram not used          call transport_integrals(nx_block,           ny_block,             &
! Subprogram not used                                   ntrace,             icellsng (:,iblk),    &
! Subprogram not used                                   indxing(:,:,iblk),  indxjng(:,:,iblk),    &
! Subprogram not used                                   tracer_type,        depend,               &
! Subprogram not used                                   integral_order,     triarea(:,:,:,iblk),  &
! Subprogram not used                                   iflux(:,:,:,iblk),  jflux  (:,:,:,iblk),  &
! Subprogram not used                                   xp (:,:,:,:,iblk),  yp   (:,:,:,:,iblk),  &
! Subprogram not used                                   mc   (:,:,0,iblk),  mx     (:,:,0,iblk),  &
! Subprogram not used                                   my   (:,:,0,iblk),  mflxn  (:,:,0,iblk))
! Subprogram not used 
! Subprogram not used          ! ice categories
! Subprogram not used          do n = 1, ncat
! Subprogram not used             call transport_integrals                                       &
! Subprogram not used                                (nx_block,           ny_block,              &
! Subprogram not used                                 ntrace,             icellsng (:,iblk),     &
! Subprogram not used                                 indxing(:,:,iblk),  indxjng   (:,:,iblk),  &
! Subprogram not used                                 tracer_type,        depend,                &
! Subprogram not used                                 integral_order,     triarea (:,:,:,iblk),  &
! Subprogram not used                                 iflux(:,:,:,iblk),  jflux   (:,:,:,iblk),  &
! Subprogram not used                                 xp (:,:,:,:,iblk),  yp    (:,:,:,:,iblk),  &
! Subprogram not used                                 mc   (:,:,n,iblk),  mx      (:,:,n,iblk),  &
! Subprogram not used                                 my   (:,:,n,iblk),  mflxn   (:,:,n,iblk),  &
! Subprogram not used                                 tc (:,:,:,n,iblk),  tx    (:,:,:,n,iblk),  &
! Subprogram not used                                 ty (:,:,:,n,iblk),  mtflxn(:,:,:,n,iblk))
! Subprogram not used 
! Subprogram not used          enddo                  ! n
! Subprogram not used          call t_stopf  ('cice_hmap_remap2n')
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Update the ice area and tracers.
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          ! open water
! Subprogram not used 
! Subprogram not used          call t_startf  ('cice_hmap_remap2_upd')
! Subprogram not used          call update_fields (nx_block,           ny_block,          &
! Subprogram not used                              ilo, ihi,           jlo, jhi,          &
! Subprogram not used                              ntrace,                                &
! Subprogram not used                              tracer_type,        depend,            &
! Subprogram not used                              tarear(:,:,iblk),   l_stop(iblk),      &
! Subprogram not used                              istop(iblk),        jstop (iblk),      &
! Subprogram not used                              mflxe(:,:,0,iblk),  mflxn(:,:,0,iblk), &
! Subprogram not used                              mm   (:,:,0,iblk))
! Subprogram not used 
! Subprogram not used          if (l_stop(iblk)) then
! Subprogram not used             this_block = get_block(blocks_ice(iblk),iblk)         
! Subprogram not used             write (nu_diag,*) 'istep1, my_task, iblk, cat =',      &
! Subprogram not used                                istep1, my_task, iblk, '0'
! Subprogram not used             write (nu_diag,*) 'Global block:', this_block%block_id
! Subprogram not used             if (istop(iblk) > 0 .and. jstop(iblk) > 0)             &
! Subprogram not used                  write(nu_diag,*) 'Global i and j:',               &
! Subprogram not used                                   this_block%i_glob(istop(iblk)),  &
! Subprogram not used                                   this_block%j_glob(jstop(iblk)) 
! Subprogram not used             call abort_ice ('ice remap_transport: negative area (open water)')
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used          ! ice categories
! Subprogram not used          do n = 1, ncat
! Subprogram not used 
! Subprogram not used             call update_fields(nx_block,              ny_block,              &
! Subprogram not used                                ilo, ihi,              jlo, jhi,              &
! Subprogram not used                                ntrace,                                       &
! Subprogram not used                                tracer_type,           depend,                &
! Subprogram not used                                tarear(:,:,iblk),      l_stop(iblk),          &
! Subprogram not used                                istop(iblk),           jstop(iblk),           &
! Subprogram not used                                mflxe (:,:,  n,iblk),  mflxn (:,:,  n,iblk),  &
! Subprogram not used                                mm    (:,:,  n,iblk),                         &
! Subprogram not used                                mtflxe(:,:,:,n,iblk),  mtflxn(:,:,:,n,iblk),  &
! Subprogram not used                                tm    (:,:,:,n,iblk))
! Subprogram not used 
! Subprogram not used             if (l_stop(iblk)) then
! Subprogram not used                this_block = get_block(blocks_ice(iblk),iblk)         
! Subprogram not used                write (nu_diag,*) 'istep1, my_task, iblk, cat =',      &
! Subprogram not used                                   istep1, my_task, iblk, n
! Subprogram not used                write (nu_diag,*) 'Global block:', this_block%block_id
! Subprogram not used                if (istop(iblk) > 0 .and. jstop(iblk) > 0)             &
! Subprogram not used                     write(nu_diag,*) 'Global i and j:',               &
! Subprogram not used                                      this_block%i_glob(istop(iblk)),  &
! Subprogram not used                                      this_block%j_glob(jstop(iblk)) 
! Subprogram not used                call abort_ice ('ice remap_transport: negative area (ice)')
! Subprogram not used             endif
! Subprogram not used          enddo                  ! n
! Subprogram not used          call t_stopf  ('cice_hmap_remap2_upd')
! Subprogram not used 
! Subprogram not used       enddo                     ! iblk
! Subprogram not used       !$OMP END PARALLEL DO
! Subprogram not used 
! Subprogram not used       call t_stopf ('cice_hmap_remap2')
! Subprogram not used 
! Subprogram not used       end subroutine horizontal_remap

!=======================================================================
!
!BOP
!
! !IROUTINE: make_masks - make area and tracer masks
!
! !INTERFACE:
!
! Subprogram not used       subroutine make_masks (nx_block, ny_block,           &
! Subprogram not used                              ilo, ihi, jlo, jhi,           &
! Subprogram not used                              nghost,   ntrace,             &
! Subprogram not used                              has_dependents,               &
! Subprogram not used                              icells,                       &
! Subprogram not used                              indxi,    indxj,              &
! Subprogram not used                              mm,       mmask,              &
! Subprogram not used                              tm,       tmask)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Make area and tracer masks.
! Subprogram not used !
! Subprogram not used ! If an area is masked out (mm < puny), then the values of tracers
! Subprogram not used !  in that grid cell are assumed to have no physical meaning.
! Subprogram not used !
! Subprogram not used ! Similarly, if a tracer with dependents is masked out
! Subprogram not used !  (abs(tm) < puny), then the values of its dependent tracers in that
! Subprogram not used !  grid cell are assumed to have no physical meaning.
! Subprogram not used ! For example, the enthalpy value has no meaning if the thickness
! Subprogram not used !  is zero.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author William H. Lipscomb, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) ::     &
! Subprogram not used            nx_block, ny_block  ,&! block dimensions
! Subprogram not used            ilo,ihi,jlo,jhi     ,&! beginning and end of physical domain
! Subprogram not used            nghost              ,&! number of ghost cells
! Subprogram not used            ntrace                ! number of tracers in use
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind), dimension (ntrace), intent(in) ::     &
! Subprogram not used            has_dependents      ! true if a tracer has dependent tracers
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension(0:ncat), intent(out) ::     &
! Subprogram not used            icells         ! number of cells with ice
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension(nx_block*ny_block,0:ncat),     &
! Subprogram not used            intent(out) ::     &
! Subprogram not used            indxi        ,&! compressed i/j indices
! Subprogram not used            indxj
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,0:ncat),     &
! Subprogram not used            intent(in) ::     &
! Subprogram not used            mm            ! mean ice area in each grid cell
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,0:ncat),     &
! Subprogram not used            intent(out) ::     &
! Subprogram not used            mmask         ! = 1. if ice is present, else = 0.
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block, ny_block, max_ntrace, ncat),  &
! Subprogram not used            intent(in), optional ::     &
! Subprogram not used            tm            ! mean tracer values in each grid cell
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block, ny_block, max_ntrace, ncat),  &
! Subprogram not used            intent(out), optional ::     &
! Subprogram not used            tmask         ! = 1. if tracer is present, else = 0.
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) ::     &
! Subprogram not used            i, j, ij       ,&! horizontal indices
! Subprogram not used            n              ,&! ice category index
! Subprogram not used            nt               ! tracer index
! Subprogram not used 
! Subprogram not used       do n = 0, ncat
! Subprogram not used          do ij = 1, nx_block*ny_block
! Subprogram not used             indxi(ij,n) = 0
! Subprogram not used             indxj(ij,n) = 0
! Subprogram not used          enddo
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! open water mask
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       icells(0) = 0
! Subprogram not used       do j = 1, ny_block
! Subprogram not used       do i = 1, nx_block
! Subprogram not used          if (mm(i,j,0) > puny) then
! Subprogram not used             mmask(i,j,0) = c1
! Subprogram not used             icells(0) = icells(0) + 1
! Subprogram not used             ij = icells(0)
! Subprogram not used             indxi(ij,0) = i
! Subprogram not used             indxj(ij,0) = j
! Subprogram not used          else
! Subprogram not used             mmask(i,j,0) = c0
! Subprogram not used          endif
! Subprogram not used       enddo
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       do n = 1, ncat
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Find grid cells where ice is present.
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          icells(n) = 0
! Subprogram not used          do j = 1, ny_block
! Subprogram not used          do i = 1, nx_block
! Subprogram not used             if (mm(i,j,n) > puny) then
! Subprogram not used                icells(n) = icells(n) + 1
! Subprogram not used                ij = icells(n)
! Subprogram not used                indxi(ij,n) = i
! Subprogram not used                indxj(ij,n) = j
! Subprogram not used             endif               ! mm > puny
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! ice area mask
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          mmask(:,:,n) = c0
! Subprogram not used          do ij = 1, icells(n)
! Subprogram not used             i = indxi(ij,n)
! Subprogram not used             j = indxj(ij,n)
! Subprogram not used             mmask(i,j,n) = c1
! Subprogram not used          enddo
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! tracer masks
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          if (present(tm)) then
! Subprogram not used 
! Subprogram not used             tmask(:,:,:,n) = c0
! Subprogram not used             do nt = 1, ntrace
! Subprogram not used                if (has_dependents(nt)) then
! Subprogram not used                   do ij = 1, icells(n)
! Subprogram not used                      i = indxi(ij,n)
! Subprogram not used                      j = indxj(ij,n)
! Subprogram not used                      if (abs(tm(i,j,nt,n)) > puny) then
! Subprogram not used                         tmask(i,j,nt,n) = c1
! Subprogram not used                      endif
! Subprogram not used                   enddo
! Subprogram not used                endif
! Subprogram not used             enddo
! Subprogram not used 
! Subprogram not used          endif                     ! present(tm)
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Redefine icells
! Subprogram not used     ! For nghost = 1, exclude ghost cells
! Subprogram not used     ! For nghost = 2, include one layer of ghost cells
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          icells(n) = 0
! Subprogram not used          do j = jlo-nghost+1, jhi+nghost-1
! Subprogram not used          do i = ilo-nghost+1, ihi+nghost-1
! Subprogram not used             if (mm(i,j,n) > puny) then
! Subprogram not used                icells(n) = icells(n) + 1
! Subprogram not used                ij = icells(n)
! Subprogram not used                indxi(ij,n) = i
! Subprogram not used                indxj(ij,n) = j
! Subprogram not used             endif               ! mm > puny
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used       
! Subprogram not used       enddo ! n
! Subprogram not used 
! Subprogram not used       end subroutine make_masks

!=======================================================================
!
!BOP
!
! !IROUTINE: construct_fields - construct fields of ice area and tracers
!
! !INTERFACE:
!
! Subprogram not used       subroutine construct_fields (nx_block,       ny_block,   &
! Subprogram not used                                    ilo, ihi,       jlo, jhi,   &
! Subprogram not used                                    nghost,         ntrace,     &
! Subprogram not used                                    tracer_type,    depend,     &
! Subprogram not used                                    has_dependents, icells,     &
! Subprogram not used                                    indxi,          indxj,      &
! Subprogram not used                                    HTN,            HTE,        &
! Subprogram not used                                    hm,             xav,        &
! Subprogram not used                                    yav,            xxav,       &
! Subprogram not used                                    xyav,           yyav,       &
! Subprogram not used                                    xxxav,          xxyav,      &
! Subprogram not used                                    xyyav,          yyyav,      &
! Subprogram not used                                    dxt,            dyt,        &
! Subprogram not used                                    mm,             mc,         &
! Subprogram not used                                    mx,             my,         &
! Subprogram not used                                    mmask,                      &
! Subprogram not used                                    tm,             tc,         &
! Subprogram not used                                    tx,             ty,         &
! Subprogram not used                                    tmask)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Construct fields of ice area and tracers.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! authors William H. Lipscomb, LANL
! Subprogram not used !         John R. Baumgardner, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) ::   &
! Subprogram not used          nx_block, ny_block  ,&! block dimensions
! Subprogram not used          ilo,ihi,jlo,jhi     ,&! beginning and end of physical domain
! Subprogram not used          nghost              ,&! number of ghost cells
! Subprogram not used          ntrace              ,&! number of tracers in use
! Subprogram not used          icells                ! number of cells with mass
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension (ntrace), intent(in) ::     &
! Subprogram not used          tracer_type       ,&! = 1, 2, or 3 (see comments above)
! Subprogram not used          depend              ! tracer dependencies (see above)
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind), dimension (ntrace), intent(in) ::     &
! Subprogram not used          has_dependents      ! true if a tracer has dependent tracers
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension(nx_block*ny_block), intent(in) :: &
! Subprogram not used          indxi          ,&! compressed i/j indices
! Subprogram not used          indxj
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block),   &
! Subprogram not used          intent(in) ::   &
! Subprogram not used          hm             ,&! land/boundary mask, thickness (T-cell)
! Subprogram not used          HTN            ,&! length of northern edge of T-cell (m)
! Subprogram not used          HTE            ,&! length of eastern edge of T-cell (m)
! Subprogram not used          xav,  yav              ,&! mean T-cell values of x, y
! Subprogram not used          xxav, xyav, yyav       ,&! mean T-cell values of xx, xy, yy
! Subprogram not used          xxxav,xxyav,xyyav,yyyav,&! mean T-cell values of , xxy, xyy, yyy
! Subprogram not used          dxt            ,&! grid cell width (m)
! Subprogram not used          dyt              ! grid cell height (m)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block),   &
! Subprogram not used          intent(in) ::   &
! Subprogram not used          mm            ,&! mean value of mass field
! Subprogram not used          mmask           ! = 1. if ice is present, = 0. otherwise
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrace),   &
! Subprogram not used          intent(in), optional ::   &
! Subprogram not used          tm             ,&! mean tracer
! Subprogram not used          tmask            ! = 1. if tracer is present, = 0. otherwise
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block),   &
! Subprogram not used          intent(out) ::   &
! Subprogram not used          mc             ,&! mass value at geometric center of cell
! Subprogram not used          mx, my           ! limited derivative of mass wrt x and y
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrace),   &
! Subprogram not used          intent(out), optional ::   &
! Subprogram not used          tc             ,&! tracer at geometric center of cell
! Subprogram not used          tx, ty           ! limited derivative of tracer wrt x and y
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) ::   &
! Subprogram not used          i, j           ,&! horizontal indices
! Subprogram not used          nt, nt1        ,&! tracer indices
! Subprogram not used          ij               ! combined i/j horizontal index
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block) ::    &
! Subprogram not used          mxav           ,&! x coordinate of center of mass
! Subprogram not used          myav             ! y coordinate of center of mass
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrace) ::  &
! Subprogram not used          mtxav          ,&! x coordinate of center of mass*tracer
! Subprogram not used          mtyav            ! y coordinate of center of mass*tracer
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind) ::   &
! Subprogram not used          w1, w2, w3, w4, w5, w6, w7   ! work variables
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Compute field values at the geometric center of each grid cell,
! Subprogram not used     ! and compute limited gradients in the x and y directions.
! Subprogram not used     !
! Subprogram not used     ! For second order accuracy, each state variable is approximated as
! Subprogram not used     ! a field varying linearly over x and y within each cell.  For each
! Subprogram not used     ! category, the integrated value of m(x,y) over the cell must
! Subprogram not used     ! equal mm(i,j,n)*tarea(i,j), where tarea(i,j) is the cell area.
! Subprogram not used     ! Similarly, the integrated value of m(x,y)*t(x,y) must equal
! Subprogram not used     ! the total mass*tracer, mm(i,j,n)*tm(i,j,n)*tarea(i,j).
! Subprogram not used     !
! Subprogram not used     ! These integral conditions are satisfied for linear fields if we
! Subprogram not used     ! stipulate the following:
! Subprogram not used     ! (1) The mean mass, mm, is equal to the mass at the cell centroid.
! Subprogram not used     ! (2) The mean value tm1 of type 1 tracers is equal to the value
! Subprogram not used     !     at the center of mass.
! Subprogram not used     ! (3) The mean value tm2 of type 2 tracers is equal to the value
! Subprogram not used     !     at the center of mass*tm1, where tm2 depends on tm1.
! Subprogram not used     !     (See comments at the top of the module.)
! Subprogram not used     !
! Subprogram not used     ! We want to find the value of each state variable at a standard
! Subprogram not used     ! reference point, which we choose to be the geometric center of
! Subprogram not used     ! the cell.  The geometric center is located at the intersection
! Subprogram not used     ! of the line joining the midpoints of the north and south edges
! Subprogram not used     ! with the line joining the midpoints of the east and west edges.
! Subprogram not used     ! To find the value at the geometric center, we must know the
! Subprogram not used     ! location of the cell centroid/center of mass, along with the
! Subprogram not used     ! mean value and the gradients with respect to x and y.
! Subprogram not used     !
! Subprogram not used     ! The cell gradients are first computed from the difference between
! Subprogram not used     ! values in the neighboring cells, then limited by requiring that
! Subprogram not used     ! no new extrema are created within the cell.
! Subprogram not used     !
! Subprogram not used     ! For rectangular coordinates the centroid and the geometric
! Subprogram not used     ! center coincide, which means that some of the equations in this
! Subprogram not used     ! subroutine could be simplified.  However, the full equations
! Subprogram not used     ! are retained for generality.
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Initialize
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       do j = 1, ny_block
! Subprogram not used       do i = 1, nx_block
! Subprogram not used          mc(i,j)  = c0
! Subprogram not used          mx(i,j)  = c0
! Subprogram not used          my(i,j)  = c0
! Subprogram not used          mxav(i,j) = c0
! Subprogram not used          myav(i,j) = c0
! Subprogram not used       enddo
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       if (present(tm)) then
! Subprogram not used          do nt = 1, ntrace
! Subprogram not used             do j = 1, ny_block
! Subprogram not used             do i = 1, nx_block
! Subprogram not used                tc(i,j,nt) = c0
! Subprogram not used                tx(i,j,nt) = c0
! Subprogram not used                ty(i,j,nt) = c0
! Subprogram not used             enddo
! Subprogram not used             enddo
! Subprogram not used          enddo
! Subprogram not used       endif
! Subprogram not used          
! Subprogram not used       ! limited gradient of mass field in each cell (except masked cells)
! Subprogram not used       ! Note: The gradient is computed in scaled coordinates with
! Subprogram not used       !       dxt = dyt = hte = htn = 1.
! Subprogram not used 
! Subprogram not used       call limited_gradient (nx_block, ny_block,   &
! Subprogram not used                              ilo, ihi, jlo, jhi,   &
! Subprogram not used                              nghost,               &
! Subprogram not used                              mm,       hm,         &
! Subprogram not used                              xav,      yav,        &
! Subprogram not used                              HTN,      HTE,        &
! Subprogram not used                              dxt,      dyt,        &
! Subprogram not used                              mx,       my)
! Subprogram not used 
! Subprogram not used       do ij = 1,icells   ! ice is present
! Subprogram not used          i = indxi(ij)
! Subprogram not used          j = indxj(ij)
! Subprogram not used 
! Subprogram not used          ! mass field at geometric center
! Subprogram not used          mc(i,j) = mm(i,j) - xav(i,j)*mx(i,j)   &
! Subprogram not used                            - yav(i,j)*my(i,j)
! Subprogram not used 
! Subprogram not used       enddo                     ! ij
! Subprogram not used 
! Subprogram not used       ! tracers
! Subprogram not used 
! Subprogram not used       if (present(tm)) then
! Subprogram not used 
! Subprogram not used        do ij = 1,icells       ! cells with mass
! Subprogram not used           i = indxi(ij)
! Subprogram not used           j = indxj(ij)
! Subprogram not used 
! Subprogram not used          ! center of mass (mxav,myav) for each cell
! Subprogram not used           mxav(i,j) = (mx(i,j)*xxav(i,j)    &
! Subprogram not used                      + my(i,j)*xyav(i,j)    &
! Subprogram not used                      + mc(i,j)*xav (i,j)) / mm(i,j)
! Subprogram not used           myav(i,j) = (mx(i,j)*xyav(i,j)    &
! Subprogram not used                      + my(i,j)*yyav(i,j)    &
! Subprogram not used                      + mc(i,j)*yav(i,j)) / mm(i,j)
! Subprogram not used        enddo
! Subprogram not used 
! Subprogram not used        do nt = 1, ntrace
! Subprogram not used 
! Subprogram not used          if (tracer_type(nt)==1) then ! independent of other tracers
! Subprogram not used 
! Subprogram not used             call limited_gradient(nx_block,     ny_block,  &
! Subprogram not used                                   ilo, ihi,     jlo, jhi,  &
! Subprogram not used                                   nghost,                  &
! Subprogram not used                                   tm(:,:,nt),   mmask,     &
! Subprogram not used                                   mxav,         myav,      &
! Subprogram not used                                   HTN,          HTE,       &
! Subprogram not used                                   dxt,          dyt,       &
! Subprogram not used                                   tx(:,:,nt),   ty(:,:,nt)) 
! Subprogram not used 
! Subprogram not used             if (has_dependents(nt)) then   ! need center of area*tracer
! Subprogram not used 
! Subprogram not used                do j = 1, ny_block
! Subprogram not used                do i = 1, nx_block
! Subprogram not used                   mtxav(i,j,nt) = c0
! Subprogram not used                   mtyav(i,j,nt) = c0
! Subprogram not used                enddo
! Subprogram not used                enddo
! Subprogram not used 
! Subprogram not used !DIR$ CONCURRENT !Cray
! Subprogram not used !cdir nodep      !NEC
! Subprogram not used !ocl novrec      !Fujitsu
! Subprogram not used                do ij = 1, icells  ! Note: no tx or ty in ghost cells
! Subprogram not used                                   ! (bound calls are later)
! Subprogram not used                   i = indxi(ij)
! Subprogram not used                   j = indxj(ij)
! Subprogram not used 
! Subprogram not used                   ! tracer value at geometric center
! Subprogram not used                   tc(i,j,nt) = tm(i,j,nt) - tx(i,j,nt)*mxav(i,j)   &
! Subprogram not used                                           - ty(i,j,nt)*myav(i,j)
! Subprogram not used 
! Subprogram not used                   if (tmask(i,j,nt) > puny) then
! Subprogram not used 
! Subprogram not used                      ! center of area*tracer
! Subprogram not used                      w1 = mc(i,j)*tc(i,j,nt)
! Subprogram not used                      w2 = mc(i,j)*tx(i,j,nt)   &
! Subprogram not used                         + mx(i,j)*tc(i,j,nt)
! Subprogram not used                      w3 = mc(i,j)*ty(i,j,nt)   &
! Subprogram not used                         + my(i,j)*tc(i,j,nt)
! Subprogram not used                      w4 = mx(i,j)*tx(i,j,nt)
! Subprogram not used                      w5 = mx(i,j)*ty(i,j,nt)   &
! Subprogram not used                         + my(i,j)*tx(i,j,nt)
! Subprogram not used                      w6 = my(i,j)*ty(i,j,nt)
! Subprogram not used                      w7 = c1 / (mm(i,j)*tm(i,j,nt))
! Subprogram not used                      mtxav(i,j,nt) = (w1*xav (i,j)  + w2*xxav (i,j)   &
! Subprogram not used                                     + w3*xyav (i,j) + w4*xxxav(i,j)   &
! Subprogram not used                                     + w5*xxyav(i,j) + w6*xyyav(i,j))   &
! Subprogram not used                                     * w7
! Subprogram not used                      mtyav(i,j,nt) = (w1*yav(i,j)   + w2*xyav (i,j)   &
! Subprogram not used                                     + w3*yyav(i,j)  + w4*xxyav(i,j)   &
! Subprogram not used                                     + w5*xyyav(i,j) + w6*yyyav(i,j))   &
! Subprogram not used                                     * w7
! Subprogram not used                   endif         ! tmask
! Subprogram not used 
! Subprogram not used                enddo            ! ij
! Subprogram not used 
! Subprogram not used             else                ! no dependents
! Subprogram not used 
! Subprogram not used                do ij = 1, icells      ! mass is present
! Subprogram not used                   i = indxi(ij)
! Subprogram not used                   j = indxj(ij)
! Subprogram not used 
! Subprogram not used                   ! tracer value at geometric center
! Subprogram not used                   tc(i,j,nt) = tm(i,j,nt) - tx(i,j,nt)*mxav(i,j)   &
! Subprogram not used                                           - ty(i,j,nt)*myav(i,j)
! Subprogram not used                enddo            ! ij
! Subprogram not used 
! Subprogram not used             endif               ! has_dependents
! Subprogram not used 
! Subprogram not used          elseif (tracer_type(nt)==2) then   ! tracer nt depends on nt1
! Subprogram not used             nt1 = depend(nt)
! Subprogram not used 
! Subprogram not used             call limited_gradient(nx_block,       ny_block,         &
! Subprogram not used                                   ilo, ihi,       jlo, jhi,         &
! Subprogram not used                                   nghost,                           &
! Subprogram not used                                   tm(:,:,nt),     tmask(:,:,nt1),   &
! Subprogram not used                                   mtxav(:,:,nt1), mtyav(:,:,nt1),   &
! Subprogram not used                                   HTN,            HTE,              &
! Subprogram not used                                   dxt,            dyt,              &
! Subprogram not used                                   tx(:,:,nt),     ty(:,:,nt))    
! Subprogram not used 
! Subprogram not used             do ij = 1, icells     ! ice is present
! Subprogram not used                i = indxi(ij)
! Subprogram not used                j = indxj(ij)
! Subprogram not used                tc(i,j,nt) = tm(i,j,nt)                    &
! Subprogram not used                           - tx(i,j,nt) * mtxav(i,j,nt1)   &
! Subprogram not used                           - ty(i,j,nt) * mtyav(i,j,nt1)
! Subprogram not used             enddo               ! ij
! Subprogram not used 
! Subprogram not used          elseif (tracer_type(nt)==3) then  ! upwind approx; gradient = 0
! Subprogram not used 
! Subprogram not used             do ij = 1, icells
! Subprogram not used                i = indxi(ij)
! Subprogram not used                j = indxj(ij)
! Subprogram not used 
! Subprogram not used                tc(i,j,nt) = tm(i,j,nt)
! Subprogram not used !               tx(i,j,nt) = c0   ! already initialized to 0.
! Subprogram not used !               ty(i,j,nt) = c0
! Subprogram not used             enddo               ! ij
! Subprogram not used 
! Subprogram not used          endif                  ! tracer_type
! Subprogram not used        enddo                    ! ntrace
! Subprogram not used 
! Subprogram not used       endif                     ! present (tm)
! Subprogram not used 
! Subprogram not used       end subroutine construct_fields

!=======================================================================
!
!BOP
!
! !IROUTINE: limited_gradient - limited gradient of a scalar field
!
! !INTERFACE:
!
! Subprogram not used       subroutine limited_gradient (nx_block, ny_block,   &
! Subprogram not used                                    ilo, ihi, jlo, jhi,   &
! Subprogram not used                                    nghost,               &
! Subprogram not used                                    phi,      phimask,    &
! Subprogram not used                                    cnx,      cny,        &
! Subprogram not used                                    HTN,      HTE,        &
! Subprogram not used                                    dxt,      dyt,        &
! Subprogram not used                                    gx,       gy)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Compute a limited gradient of the scalar field phi in scaled coordinates.
! Subprogram not used ! "Limited" means that we do not create new extrema in phi.  For
! Subprogram not used ! instance, field values at the cell corners can neither exceed the
! Subprogram not used ! maximum of phi(i,j) in the cell and its eight neighbors, nor fall
! Subprogram not used ! below the minimum.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! authors William H. Lipscomb, LANL
! Subprogram not used !         John R. Baumgardner, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) ::   &
! Subprogram not used           nx_block, ny_block,&! block dimensions
! Subprogram not used           ilo,ihi,jlo,jhi ,&! beginning and end of physical domain
! Subprogram not used           nghost              ! number of ghost cells
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block),   &
! Subprogram not used            intent (in) ::   &
! Subprogram not used           phi    ,&! input tracer field (mean values in each grid cell)
! Subprogram not used           cnx    ,&! x-coordinate of phi relative to geometric center of cell
! Subprogram not used           cny    ,&! y-coordinate of phi relative to geometric center of cell
! Subprogram not used           dxt    ,&! grid cell width (m)
! Subprogram not used           dyt    ,&! grid cell height (m)
! Subprogram not used           phimask ,&
! Subprogram not used           ! phimask(i,j) = 1 if phi(i,j) has physical meaning, = 0 otherwise.
! Subprogram not used           ! For instance, aice has no physical meaning in land cells,
! Subprogram not used           ! and hice no physical meaning where aice = 0.
! Subprogram not used           HTN    ,&! length of northern edge of T-cell (m)
! Subprogram not used           HTE      ! length of eastern edge of T-cell (m)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block),   &
! Subprogram not used           intent(out) ::   &
! Subprogram not used           gx     ,&! limited x-direction gradient
! Subprogram not used           gy       ! limited y-direction gradient
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) ::   &
! Subprogram not used           i, j, ij        ,&! standard indices
! Subprogram not used           icells            ! number of cells to limit
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension(nx_block*ny_block) ::   &
! Subprogram not used           indxi, indxj   ! combined i/j horizontal indices
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind) ::   &
! Subprogram not used           phi_nw, phi_n, phi_ne ,&! values of phi in 8 neighbor cells
! Subprogram not used           phi_w,         phi_e  ,&
! Subprogram not used           phi_sw, phi_s, phi_se ,&
! Subprogram not used           qmn, qmx     ,&! min and max value of phi within grid cell
! Subprogram not used           pmn, pmx     ,&! min and max value of phi among neighbor cells
! Subprogram not used           w1, w2, w3, w4 ! work variables
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind) ::   &
! Subprogram not used           gxtmp, gytmp   ! temporary term for x- and y- limited gradient
! Subprogram not used 
! Subprogram not used       gx(:,:) = c0
! Subprogram not used       gy(:,:) = c0
! Subprogram not used 
! Subprogram not used       ! For nghost = 1, loop over physical cells and update ghost cells later
! Subprogram not used       ! For nghost = 2, loop over a layer of ghost cells and skip the update
! Subprogram not used 
! Subprogram not used       icells = 0
! Subprogram not used       do j = jlo-nghost+1, jhi+nghost-1
! Subprogram not used       do i = ilo-nghost+1, ihi+nghost-1
! Subprogram not used          if (phimask(i,j) > puny) then
! Subprogram not used 
! Subprogram not used !jw            icells = icells + 1
! Subprogram not used !jw            indxi(icells) = i
! Subprogram not used !jw            indxj(icells) = j
! Subprogram not used !jw         endif                  ! phimask > puny
! Subprogram not used !jw      enddo
! Subprogram not used !jw      enddo
! Subprogram not used 
! Subprogram not used !jw      do ij = 1, icells
! Subprogram not used !jw         i = indxi(ij)
! Subprogram not used !jw         j = indxj(ij)
! Subprogram not used 
! Subprogram not used          ! Store values of phi in the 8 neighbor cells.
! Subprogram not used          ! Note: phimask = 1. or 0.  If phimask = 1., use the true value;
! Subprogram not used          !  if phimask = 0., use the home cell value so that non-physical
! Subprogram not used          !  values of phi do not contribute to the gradient.
! Subprogram not used          phi_nw = phimask(i-1,j+1) * phi(i-1,j+1)   &
! Subprogram not used             + (c1-phimask(i-1,j+1))* phi(i,j)
! Subprogram not used          phi_n  = phimask(i,j+1)   * phi(i,j+1)   &
! Subprogram not used             + (c1-phimask(i,j+1))  * phi(i,j)
! Subprogram not used          phi_ne = phimask(i+1,j+1) * phi(i+1,j+1)   &
! Subprogram not used             + (c1-phimask(i+1,j+1))* phi(i,j)
! Subprogram not used          phi_w  = phimask(i-1,j)   * phi(i-1,j)   &
! Subprogram not used             + (c1-phimask(i-1,j))  * phi(i,j)
! Subprogram not used          phi_e  = phimask(i+1,j)   * phi(i+1,j)   &
! Subprogram not used             + (c1-phimask(i+1,j))  * phi(i,j)
! Subprogram not used          phi_sw = phimask(i-1,j-1) * phi(i-1,j-1)   &
! Subprogram not used             + (c1-phimask(i-1,j-1))* phi(i,j)
! Subprogram not used          phi_s  = phimask(i,j-1)   * phi(i,j-1)   &
! Subprogram not used             + (c1-phimask(i,j-1))  * phi(i,j)
! Subprogram not used          phi_se = phimask(i+1,j-1) * phi(i+1,j-1)   &
! Subprogram not used             + (c1-phimask(i+1,j-1))* phi(i,j)
! Subprogram not used 
! Subprogram not used          ! unlimited gradient components
! Subprogram not used          ! (factors of two cancel out)
! Subprogram not used 
! Subprogram not used          gxtmp = (phi_e - phi(i,j)) / (dxt(i,j)   + dxt(i+1,j))   &
! Subprogram not used                + (phi(i,j) - phi_w) / (dxt(i-1,j) + dxt(i,j)  )
! Subprogram not used          gytmp = (phi_n - phi(i,j)) / (dyt(i,j)   + dyt(i,j+1))   &
! Subprogram not used                + (phi(i,j) - phi_s) / (dyt(i,j-1) + dyt(i,j)  )
! Subprogram not used 
! Subprogram not used          ! minimum and maximum among the nine local cells
! Subprogram not used          pmn = min (phi_nw, phi_n,  phi_ne, phi_w, phi(i,j),   &
! Subprogram not used                     phi_e,  phi_sw, phi_s,  phi_se)
! Subprogram not used          pmx = max (phi_nw, phi_n,  phi_ne, phi_w, phi(i,j),   &
! Subprogram not used                     phi_e,  phi_sw, phi_s,  phi_se)
! Subprogram not used 
! Subprogram not used          pmn = pmn - phi(i,j)
! Subprogram not used          pmx = pmx - phi(i,j)
! Subprogram not used 
! Subprogram not used          ! minimum and maximum deviation of phi within the cell
! Subprogram not used 
! Subprogram not used          w1  =  (p5*HTN(i,j)   - cnx(i,j)) * gxtmp   &
! Subprogram not used               + (p5*HTE(i,j)   - cny(i,j)) * gytmp
! Subprogram not used          w2  =  (p5*HTN(i,j-1) - cnx(i,j)) * gxtmp   &
! Subprogram not used               - (p5*HTE(i,j)   + cny(i,j)) * gytmp
! Subprogram not used          w3  = -(p5*HTN(i,j-1) + cnx(i,j)) * gxtmp   &
! Subprogram not used               - (p5*HTE(i-1,j) + cny(i,j)) * gytmp
! Subprogram not used          w4  =  (p5*HTE(i-1,j) - cny(i,j)) * gytmp   &
! Subprogram not used               - (p5*HTN(i,j)   + cnx(i,j)) * gxtmp
! Subprogram not used 
! Subprogram not used          qmn = min (w1, w2, w3, w4)
! Subprogram not used          qmx = max (w1, w2, w3, w4)
! Subprogram not used 
! Subprogram not used          ! Watch for underflows here
! Subprogram not used 
! Subprogram not used          ! the limiting coefficient
! Subprogram not used          if (abs(qmn) > 1.0e-300_dbl_kind) then ! 'abs(qmn) > puny' not sufficient
! Subprogram not used             w1 = max(c0, pmn/qmn)
! Subprogram not used          else
! Subprogram not used             w1 = c1
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used          if (abs(qmx) > 1.0e-300_dbl_kind) then
! Subprogram not used             w2 = max(c0, pmx/qmx)
! Subprogram not used          else
! Subprogram not used             w2 = c1
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used          w1 = min(c1, w1, w2)
! Subprogram not used 
! Subprogram not used          ! Limit the gradient components
! Subprogram not used          gx(i,j) = w1 * gxtmp
! Subprogram not used          gy(i,j) = w1 * gytmp
! Subprogram not used 
! Subprogram not used !jw      enddo                     ! ij
! Subprogram not used           endif
! Subprogram not used        enddo
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       end subroutine limited_gradient

!=======================================================================
!BOP
!
! !IROUTINE: departure_points - compute departure points of trajectories
!
! !INTERFACE:
!
! Subprogram not used       subroutine departure_points (nx_block,   ny_block,   &
! Subprogram not used                                    ilo, ihi,   jlo, jhi,   &
! Subprogram not used                                    nghost,     dt,   &
! Subprogram not used                                    uvel,       vvel,    &
! Subprogram not used                                    dxu,        dyu,     &
! Subprogram not used                                    HTN,        HTE,     &
! Subprogram not used                                    dpx,        dpy,     &
! Subprogram not used                                    l_dp_midpt, l_stop,   &
! Subprogram not used                                    istop,      jstop)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Given velocity fields on cell corners, compute departure points
! Subprogram not used ! of back trajectories in nondimensional coordinates.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author William H. Lipscomb, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) ::   &
! Subprogram not used          nx_block, ny_block,&! block dimensions
! Subprogram not used          ilo,ihi,jlo,jhi,   &! beginning and end of physical domain
! Subprogram not used          nghost              ! number of ghost cells
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), intent(in) ::   &
! Subprogram not used          dt               ! time step (s)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) ::   &
! Subprogram not used          uvel           ,&! x-component of velocity (m/s)
! Subprogram not used          vvel           ,&! y-component of velocity (m/s)
! Subprogram not used          dxu            ,&! E-W dimensions of U-cell (m)
! Subprogram not used          dyu            ,&! N-S dimensions of U-cell (m)
! Subprogram not used          HTN            ,&! length of north face of T-cell (m) 
! Subprogram not used          HTE              ! length of east face of T-cell (m) 
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out) ::   &
! Subprogram not used          dpx           ,&! coordinates of departure points (m)
! Subprogram not used          dpy             ! coordinates of departure points (m)
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind), intent(in) ::   &
! Subprogram not used          l_dp_midpt          ! if true, find departure points using
! Subprogram not used                              ! corrected midpoint velocity
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind), intent(inout) ::   &
! Subprogram not used          l_stop       ! if true, abort on return
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), intent(inout) ::   &
! Subprogram not used          istop, jstop     ! indices of grid cell where model aborts 
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) ::   &
! Subprogram not used          i, j, i2, j2     ! horizontal indices
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind) ::                  &
! Subprogram not used          mpx,  mpy      ,&! coordinates of midpoint of back trajectory,
! Subprogram not used                           ! relative to cell corner
! Subprogram not used          mpxt, mpyt     ,&! midpoint coordinates relative to cell center
! Subprogram not used          ump,  vmp        ! corrected velocity at midpoint
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Estimate departure points.
! Subprogram not used     ! This estimate is 1st-order accurate in time; improve accuracy by
! Subprogram not used     !  using midpoint approximation (to add later).
! Subprogram not used     ! For nghost = 1, loop over physical cells and update ghost cells later.
! Subprogram not used     ! For nghost = 2, loop over a layer of ghost cells and skip update.
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       dpx(:,:) = c0
! Subprogram not used       dpy(:,:) = c0
! Subprogram not used 
! Subprogram not used       do j = jlo-nghost+1, jhi+nghost-1
! Subprogram not used       do i = ilo-nghost+1, ihi+nghost-1
! Subprogram not used 
! Subprogram not used          dpx(i,j) = -dt*uvel(i,j)
! Subprogram not used          dpy(i,j) = -dt*vvel(i,j)
! Subprogram not used 
! Subprogram not used          ! Check for values out of bounds (more than one grid cell away)
! Subprogram not used          if (dpx(i,j) < -HTN(i,j) .or. dpx(i,j) > HTN(i+1,j) .or.   &
! Subprogram not used              dpy(i,j) < -HTE(i,j) .or. dpy(i,j) > HTE(i,j+1)) then
! Subprogram not used             l_stop = .true.
! Subprogram not used             istop = i
! Subprogram not used             jstop = j
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used       enddo
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       if (l_stop) then
! Subprogram not used          i = istop
! Subprogram not used          j = jstop
! Subprogram not used          write (nu_diag,*) ' '
! Subprogram not used          write (nu_diag,*)   &
! Subprogram not used                     'Warning: Departure points out of bounds in remap'
! Subprogram not used          write (nu_diag,*) 'my_task, i, j =', my_task, i, j
! Subprogram not used          write (nu_diag,*) 'dpx, dpy =', dpx(i,j), dpy(i,j)
! Subprogram not used          write (nu_diag,*) 'HTN(i,j), HTN(i+1,j) =', HTN(i,j), HTN(i+1,j)
! Subprogram not used          write (nu_diag,*) 'HTE(i,j), HTE(i,j+1) =', HTE(i,j), HTE(i,j+1)
! Subprogram not used          return
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       if (l_dp_midpt) then ! find dep pts using corrected midpt velocity 
! Subprogram not used 
! Subprogram not used        do j = jlo-nghost+1, jhi+nghost-1
! Subprogram not used        do i = ilo-nghost+1, ihi+nghost-1
! Subprogram not used          if (uvel(i,j)/=c0 .or. vvel(i,j)/=c0) then
! Subprogram not used  
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Scale departure points to coordinate system in which grid cells
! Subprogram not used     ! have sides of unit length.
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used             dpx(i,j) = dpx(i,j) / dxu(i,j)
! Subprogram not used             dpy(i,j) = dpy(i,j) / dyu(i,j)
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Estimate midpoint of backward trajectory relative to corner (i,j).
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used             mpx = p5 * dpx(i,j)
! Subprogram not used             mpy = p5 * dpy(i,j)
! Subprogram not used  
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Determine the indices (i2,j2) of the cell where the trajectory lies.
! Subprogram not used     ! Compute the coordinates of the midpoint of the backward trajectory
! Subprogram not used     !  relative to the cell center in a stretch coordinate system
! Subprogram not used     !  with vertices at (1/2, 1/2), (1/2, -1/2), etc.
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used             if (mpx >= c0 .and. mpy >= c0) then    ! cell (i+1,j+1)
! Subprogram not used                i2 = i+1
! Subprogram not used                j2 = j+1
! Subprogram not used                mpxt = mpx - p5
! Subprogram not used                mpyt = mpy - p5
! Subprogram not used             elseif (mpx < c0 .and. mpy < c0) then  ! cell (i,j)
! Subprogram not used                i2 = i
! Subprogram not used                j2 = j
! Subprogram not used                mpxt = mpx + p5
! Subprogram not used                mpyt = mpy + p5
! Subprogram not used             elseif (mpx >= c0 .and. mpy < c0) then ! cell (i+1,j)
! Subprogram not used                i2 = i+1
! Subprogram not used                j2 = j
! Subprogram not used                mpxt = mpx - p5
! Subprogram not used                mpyt = mpy + p5
! Subprogram not used             elseif (mpx < c0 .and. mpy >= c0) then ! cell (i,j+1)
! Subprogram not used                i2 = i
! Subprogram not used                j2 = j+1
! Subprogram not used                mpxt = mpx + p5
! Subprogram not used                mpyt = mpy - p5
! Subprogram not used             endif
! Subprogram not used             
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Using a bilinear approximation, estimate the velocity at the
! Subprogram not used     ! trajectory midpoint in the (i2,j2) reference frame.
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used  
! Subprogram not used             ump = uvel(i2-1,j2-1)*(mpxt-p5)*(mpyt-p5)     &
! Subprogram not used                 - uvel(i2,  j2-1)*(mpxt+p5)*(mpyt-p5)     &
! Subprogram not used                 + uvel(i2,  j2  )*(mpxt+p5)*(mpyt+p5)     &  
! Subprogram not used                 - uvel(i2-1,j2  )*(mpxt-p5)*(mpyt+p5)
! Subprogram not used  
! Subprogram not used             vmp = vvel(i2-1,j2-1)*(mpxt-p5)*(mpyt-p5)     &
! Subprogram not used                 - vvel(i2,  j2-1)*(mpxt+p5)*(mpyt-p5)     &
! Subprogram not used                 + vvel(i2,  j2  )*(mpxt+p5)*(mpyt+p5)     &
! Subprogram not used                 - vvel(i2-1,j2  )*(mpxt-p5)*(mpyt+p5)
! Subprogram not used  
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Use the midpoint velocity to estimate the coordinates of the
! Subprogram not used     !  departure point relative to corner (i,j).
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used  
! Subprogram not used             dpx(i,j) = -dt * ump
! Subprogram not used             dpy(i,j) = -dt * vmp
! Subprogram not used  
! Subprogram not used          endif               ! nonzero velocity
! Subprogram not used 
! Subprogram not used        enddo                 ! i
! Subprogram not used        enddo                 ! j
! Subprogram not used  
! Subprogram not used       endif                  ! l_dp_midpt
! Subprogram not used 
! Subprogram not used       end subroutine departure_points

!=======================================================================
!
!BOP
!
! !IROUTINE: locate_triangles - triangle info for cell edges
!
! !INTERFACE:
!
! Subprogram not used       subroutine locate_triangles (nx_block,     ny_block,   &
! Subprogram not used                                    ilo, ihi,     jlo, jhi,   &
! Subprogram not used                                    nghost,       edge,       &
! Subprogram not used                                    icells,                   &
! Subprogram not used                                    indxi,        indxj,      &
! Subprogram not used                                    dpx,          dpy,        &
! Subprogram not used                                    dxu,          dyu,        &
! Subprogram not used                                    xp,           yp,         &
! Subprogram not used                                    iflux,        jflux,      &
! Subprogram not used                                    triarea,                  &
! Subprogram not used                                    l_fixed_area, edgearea)
! Subprogram not used !
! Subprogram not used 
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Compute areas and vertices of transport triangles for north or
! Subprogram not used !  east cell edges.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! authors William H. Lipscomb, LANL
! Subprogram not used !         John R. Baumgardner, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) ::   &
! Subprogram not used          nx_block, ny_block,&! block dimensions
! Subprogram not used          ilo,ihi,jlo,jhi   ,&! beginning and end of physical domain
! Subprogram not used          nghost              ! number of ghost cells
! Subprogram not used 
! Subprogram not used       character (len=char_len), intent(in) ::   &
! Subprogram not used          edge             ! 'north' or 'east'
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) ::  &
! Subprogram not used          dpx            ,&! x coordinates of departure points at cell corners
! Subprogram not used          dpy            ,&! y coordinates of departure points at cell corners
! Subprogram not used          dxu            ,&! E-W dimension of U-cell (m)
! Subprogram not used          dyu              ! N-S dimension of U-cell (m)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,0:nvert,ngroups),   &
! Subprogram not used          intent(out) ::   &
! Subprogram not used          xp, yp           ! coordinates of triangle vertices
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,ngroups),   &
! Subprogram not used            intent(out) ::   &
! Subprogram not used          triarea          ! area of departure triangle
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension (nx_block,ny_block,ngroups),    &
! Subprogram not used          intent(out) ::   &
! Subprogram not used          iflux          ,&! i index of cell contributing transport
! Subprogram not used          jflux            ! j index of cell contributing transport
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension (ngroups), intent(out) ::   &
! Subprogram not used          icells           ! number of cells where triarea > puny
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension (nx_block*ny_block,ngroups), &
! Subprogram not used          intent(out) ::                                               &
! Subprogram not used          indxi          ,&! compressed index in i-direction
! Subprogram not used          indxj            ! compressed index in j-direction
! Subprogram not used 
! Subprogram not used       logical, intent(in) ::   &
! Subprogram not used          l_fixed_area     ! if true, the area of each departure region is
! Subprogram not used                           !  passed in as edgearea
! Subprogram not used                           ! if false, edgearea if determined internally
! Subprogram not used                           !  and is passed out
! Subprogram not used                           
! Subprogram not used       real (kind=dbl_kind), dimension(nx_block,ny_block), intent(inout) ::   &
! Subprogram not used          edgearea         ! area of departure region for each edge
! Subprogram not used                           ! edgearea > 0 for eastward/northward flow
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) ::   &
! Subprogram not used          i, j, ij, ic   ,&! horizontal indices
! Subprogram not used          ib, ie, jb, je ,&! limits for loops over edges
! Subprogram not used          ng, nv         ,&! triangle indices
! Subprogram not used          ishift, jshift   ! differences between neighbor cells
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind) ::   &
! Subprogram not used          icellsd          ! number of cells where departure area > 0.
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension (nx_block*ny_block) ::  &
! Subprogram not used          indxid         ,&! compressed index in i-direction
! Subprogram not used          indxjd           ! compressed index in j-direction
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension(nx_block,ny_block) ::   &
! Subprogram not used          dx, dy         ,&! scaled departure points
! Subprogram not used          areafac_c      ,&! area scale factor at center of edge
! Subprogram not used          areafac_l      ,&! area scale factor at left corner
! Subprogram not used          areafac_r        ! area scale factor at right corner
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind) ::   &
! Subprogram not used          xcl, ycl       ,&! coordinates of left corner point
! Subprogram not used                           ! (relative to midpoint of edge)
! Subprogram not used          xdl, ydl       ,&! left departure point
! Subprogram not used          xil, yil       ,&! left intersection point
! Subprogram not used          xcr, ycr       ,&! right corner point
! Subprogram not used          xdr, ydr       ,&! right departure point
! Subprogram not used          xir, yir       ,&! right intersection point
! Subprogram not used          xic, yic       ,&! x-axis intersection point
! Subprogram not used          xicl, yicl     ,&! left-hand x-axis intersection point
! Subprogram not used          xicr, yicr     ,&! right-hand x-axis intersection point
! Subprogram not used          xdm, ydm       ,&! midpoint of segment connecting DL and DR;
! Subprogram not used                           ! shifted if l_fixed_area = T
! Subprogram not used          dxc            ,&! xcr - xcl
! Subprogram not used          dxd            ,&! xdr - xdl
! Subprogram not used          md             ,&! slope of line connecting DL and DR
! Subprogram not used          mdl            ,&! slope of line connecting DL and DM
! Subprogram not used          mdr            ,&! slope of line connecting DR and DM
! Subprogram not used          ishift_tl, jshift_tl ,&! i,j indices of TL cell relative to edge
! Subprogram not used          ishift_bl, jshift_bl ,&! i,j indices of BL cell relative to edge
! Subprogram not used          ishift_tr, jshift_tr ,&! i,j indices of TR cell relative to edge
! Subprogram not used          ishift_br, jshift_br ,&! i,j indices of BR cell relative to edge
! Subprogram not used          ishift_tc, jshift_tc ,&! i,j indices of TC cell relative to edge
! Subprogram not used          ishift_bc, jshift_bc ,&! i,j indices of BC cell relative to edge
! Subprogram not used          area1, area2         ,&! temporary triangle areas
! Subprogram not used          area3, area4         ,&! 
! Subprogram not used          area_c               ,&! center polygon area
! Subprogram not used          w1, w2                 ! work variables
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block,ny_block,ngroups) ::   &
! Subprogram not used          areafact         ! = 1 for positive flux, -1 for negative
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension(nx_block,ny_block) ::   &
! Subprogram not used          areasum          ! sum of triangle areas for a given edge
! Subprogram not used       
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Triangle notation:
! Subprogram not used     ! For each edge, there are 20 triangles that can contribute,
! Subprogram not used     ! but many of these are mutually exclusive.  It turns out that
! Subprogram not used     ! at most 5 triangles can contribute to transport integrals at once.
! Subprogram not used     !
! Subprogram not used     ! See Figure 3 in DB for pictures of these triangles.
! Subprogram not used     ! See Table 1 in DB for logical conditions.
! Subprogram not used     !
! Subprogram not used     ! For the north edge, DB refer to these triangles as:
! Subprogram not used     ! (1) NW, NW1, W, W2
! Subprogram not used     ! (2) NE, NE1, E, E2
! Subprogram not used     ! (3) NW2, W1, NE2, E1
! Subprogram not used     ! (4) H1a, H1b, N1a, N1b
! Subprogram not used     ! (5) H2a, H2b, N2a, N2b
! Subprogram not used     !
! Subprogram not used     ! For the east edge, DB refer to these triangles as:
! Subprogram not used     ! (1) NE, NE1, N, N2
! Subprogram not used     ! (2) SE, SE1, S, S2
! Subprogram not used     ! (3) NE2, N1, SE2, S1
! Subprogram not used     ! (4) H1a, H1b, E1a, E2b
! Subprogram not used     ! (5) H2a, H2b, E2a, E2b
! Subprogram not used     !
! Subprogram not used     ! The code below works for either north or east edges.
! Subprogram not used     ! The respective triangle labels are:
! Subprogram not used     ! (1) TL,  TL1, BL,  BL2
! Subprogram not used     ! (2) TR,  TR1, BR,  BR2
! Subprogram not used     ! (3) TL2, BL1, TR2, BR1
! Subprogram not used     ! (4) BC1a, BC1b, TC1a, TC2b
! Subprogram not used     ! (5) BC2a, BC2b, TC2a, TC2b
! Subprogram not used     ! 
! Subprogram not used     ! where the cell labels are:
! Subprogram not used     ! 
! Subprogram not used     !          |        |
! Subprogram not used     !     TL   |   TC   |   TR     (top left, center, right)
! Subprogram not used     !          |        |
! Subprogram not used     !   ------------------------
! Subprogram not used     !          |        |
! Subprogram not used     !     BL   |   BC   |   BR     (bottom left, center, right)
! Subprogram not used     !          |        |
! Subprogram not used     !
! Subprogram not used     ! and the transport is across the edge between cells TC and TB.
! Subprogram not used     !
! Subprogram not used     ! Departure points are scaled to a local coordinate system
! Subprogram not used     !  whose origin is at the midpoint of the edge.
! Subprogram not used     ! In this coordinate system, the lefthand corner CL = (-0.5,0)
! Subprogram not used     !  and the righthand corner CR = (0.5, 0).
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used   
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Initialize
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       areafac_c(:,:) = c0
! Subprogram not used       areafac_l(:,:) = c0
! Subprogram not used       areafac_r(:,:) = c0
! Subprogram not used       do ng = 1, ngroups
! Subprogram not used          do j = 1, ny_block
! Subprogram not used          do i = 1, nx_block
! Subprogram not used             triarea (i,j,ng) = c0
! Subprogram not used             areafact(i,j,ng) = c0
! Subprogram not used             iflux   (i,j,ng) = i
! Subprogram not used             jflux   (i,j,ng) = j
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used          do nv = 0, nvert
! Subprogram not used             do j = 1, ny_block
! Subprogram not used             do i = 1, nx_block
! Subprogram not used                xp(i,j,nv,ng) = c0
! Subprogram not used                yp(i,j,nv,ng) = c0
! Subprogram not used             enddo
! Subprogram not used             enddo
! Subprogram not used          enddo
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       if (trim(edge) == 'north') then
! Subprogram not used 
! Subprogram not used          ! loop size
! Subprogram not used 
! Subprogram not used          ib = ilo
! Subprogram not used          ie = ihi 
! Subprogram not used          jb = jlo - nghost            ! lowest j index is a ghost cell
! Subprogram not used          je = jhi
! Subprogram not used 
! Subprogram not used          ! index shifts for neighbor cells
! Subprogram not used 
! Subprogram not used          ishift_tl = -1
! Subprogram not used          jshift_tl =  1
! Subprogram not used          ishift_bl = -1
! Subprogram not used          jshift_bl =  0
! Subprogram not used          ishift_tr =  1
! Subprogram not used          jshift_tr =  1
! Subprogram not used          ishift_br =  1
! Subprogram not used          jshift_br =  0
! Subprogram not used          ishift_tc =  0
! Subprogram not used          jshift_tc =  1
! Subprogram not used          ishift_bc =  0
! Subprogram not used          jshift_bc =  0
! Subprogram not used 
! Subprogram not used          ! area scale factor
! Subprogram not used 
! Subprogram not used          do j = jb, je
! Subprogram not used          do i = ib, ie
! Subprogram not used             areafac_l(i,j) = dxu(i-1,j)*dyu(i-1,j) 
! Subprogram not used             areafac_r(i,j) = dxu(i,j)*dyu(i,j) 
! Subprogram not used             areafac_c(i,j) = p5*(areafac_l(i,j) + areafac_r(i,j))
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used 
! Subprogram not used       else                      ! east edge
! Subprogram not used 
! Subprogram not used          ! loop size
! Subprogram not used 
! Subprogram not used          ib = ilo - nghost            ! lowest i index is a ghost cell
! Subprogram not used          ie = ihi
! Subprogram not used          jb = jlo
! Subprogram not used          je = jhi
! Subprogram not used 
! Subprogram not used          ! index shifts for neighbor cells
! Subprogram not used 
! Subprogram not used          ishift_tl =  1
! Subprogram not used          jshift_tl =  1
! Subprogram not used          ishift_bl =  0
! Subprogram not used          jshift_bl =  1
! Subprogram not used          ishift_tr =  1
! Subprogram not used          jshift_tr = -1
! Subprogram not used          ishift_br =  0
! Subprogram not used          jshift_br = -1
! Subprogram not used          ishift_tc =  1
! Subprogram not used          jshift_tc =  0
! Subprogram not used          ishift_bc =  0
! Subprogram not used          jshift_bc =  0
! Subprogram not used 
! Subprogram not used          ! area scale factors
! Subprogram not used 
! Subprogram not used          do j = jb, je
! Subprogram not used          do i = ib, ie
! Subprogram not used             areafac_l(i,j) = dxu(i,j)*dyu(i,j) 
! Subprogram not used             areafac_r(i,j) = dxu(i,j-1)*dyu(i,j-1)
! Subprogram not used             areafac_c(i,j) = p5 * (areafac_l(i,j) + areafac_r(i,j))
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used 
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Compute mask for edges with nonzero departure areas
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       if (l_fixed_area) then
! Subprogram not used          icellsd = 0
! Subprogram not used          do j = jb, je
! Subprogram not used          do i = ib, ie
! Subprogram not used             if (edgearea(i,j) /= c0) then
! Subprogram not used                icellsd = icellsd + 1
! Subprogram not used                indxid(icellsd) = i
! Subprogram not used                indxjd(icellsd) = j
! Subprogram not used             endif
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used       else
! Subprogram not used          icellsd = 0
! Subprogram not used          if (trim(edge) == 'north') then
! Subprogram not used             do j = jb, je
! Subprogram not used             do i = ib, ie
! Subprogram not used                if (dpx(i-1,j)/=c0 .or. dpy(i-1,j)/=c0   &
! Subprogram not used                                   .or.                  &
! Subprogram not used                      dpx(i,j)/=c0 .or.   dpy(i,j)/=c0) then
! Subprogram not used                   icellsd = icellsd + 1
! Subprogram not used                   indxid(icellsd) = i
! Subprogram not used                   indxjd(icellsd) = j
! Subprogram not used                endif
! Subprogram not used             enddo
! Subprogram not used             enddo
! Subprogram not used          else       ! east edge
! Subprogram not used             do j = jb, je
! Subprogram not used             do i = ib, ie
! Subprogram not used                if (dpx(i,j-1)/=c0 .or. dpy(i,j-1)/=c0   &
! Subprogram not used                                   .or.                  &
! Subprogram not used                      dpx(i,j)/=c0 .or.   dpy(i,j)/=c0) then
! Subprogram not used                   icellsd = icellsd + 1
! Subprogram not used                   indxid(icellsd) = i
! Subprogram not used                   indxjd(icellsd) = j
! Subprogram not used                endif
! Subprogram not used             enddo
! Subprogram not used             enddo
! Subprogram not used          endif       ! edge = north/east
! Subprogram not used       endif          ! l_fixed_area
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Scale the departure points
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       do j = 1, je
! Subprogram not used       do i = 1, ie
! Subprogram not used          dx(i,j) = dpx(i,j) / dxu(i,j)
! Subprogram not used          dy(i,j) = dpy(i,j) / dyu(i,j)
! Subprogram not used       enddo
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Compute departure regions, divide into triangles, and locate
! Subprogram not used     !  vertices of each triangle.
! Subprogram not used     ! Work in a nondimensional coordinate system in which lengths are
! Subprogram not used     !  scaled by the local metric coefficients (dxu and dyu).
! Subprogram not used     ! Note: The do loop includes north faces of the j = 1 ghost cells
! Subprogram not used     !       when edge = 'north'.  The loop includes east faces of i = 1
! Subprogram not used     !       ghost cells when edge = 'east'.
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       do ij = 1, icellsd
! Subprogram not used          i = indxid(ij)
! Subprogram not used          j = indxjd(ij)
! Subprogram not used   
! Subprogram not used          xcl = -p5
! Subprogram not used          ycl =  c0
! Subprogram not used 
! Subprogram not used          xcr =  p5
! Subprogram not used          ycr =  c0
! Subprogram not used 
! Subprogram not used          ! Departure points
! Subprogram not used 
! Subprogram not used          if (trim(edge) == 'north') then ! north edge
! Subprogram not used             xdl = xcl + dx(i-1,j)
! Subprogram not used             ydl = ycl + dy(i-1,j)
! Subprogram not used             xdr = xcr + dx(i,j)
! Subprogram not used             ydr = ycr + dy(i,j)
! Subprogram not used          else                   ! east edge; rotate trajectory by pi/2
! Subprogram not used             xdl = xcl - dy(i,j)
! Subprogram not used             ydl = ycl + dx(i,j)
! Subprogram not used             xdr = xcr - dy(i,j-1)
! Subprogram not used             ydr = ycr + dx(i,j-1)
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used          xdm = p5 * (xdr + xdl)
! Subprogram not used          ydm = p5 * (ydr + ydl)
! Subprogram not used 
! Subprogram not used          ! Intersection points
! Subprogram not used 
! Subprogram not used          xil = xcl
! Subprogram not used          yil = (xcl*(ydm-ydl) + xdm*ydl - xdl*ydm) / (xdm - xdl)
! Subprogram not used          
! Subprogram not used          xir = xcr
! Subprogram not used          yir = (xcr*(ydr-ydm) - xdm*ydr + xdr*ydm) / (xdr - xdm) 
! Subprogram not used          
! Subprogram not used          md = (ydr - ydl) / (xdr - xdl)
! Subprogram not used          
! Subprogram not used          if (abs(md) > puny) then
! Subprogram not used             xic = xdl - ydl/md
! Subprogram not used          else
! Subprogram not used             xic = c0
! Subprogram not used          endif
! Subprogram not used          yic = c0
! Subprogram not used 
! Subprogram not used          xicl = xic
! Subprogram not used          yicl = yic
! Subprogram not used          xicr = xic
! Subprogram not used          yicr = yic
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Locate triangles in TL cell (NW for north edge, NE for east edge)
! Subprogram not used     ! and BL cell (W for north edge, N for east edge).
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          if (yil > c0 .and. xdl < xcl .and. ydl >= c0) then
! Subprogram not used 
! Subprogram not used          ! TL (group 1)
! Subprogram not used 
! Subprogram not used             ng = 1
! Subprogram not used             xp    (i,j,1,ng) = xcl
! Subprogram not used             yp    (i,j,1,ng) = ycl
! Subprogram not used             xp    (i,j,2,ng) = xil
! Subprogram not used             yp    (i,j,2,ng) = yil
! Subprogram not used             xp    (i,j,3,ng) = xdl
! Subprogram not used             yp    (i,j,3,ng) = ydl
! Subprogram not used             iflux   (i,j,ng) = i + ishift_tl
! Subprogram not used             jflux   (i,j,ng) = j + jshift_tl
! Subprogram not used             areafact(i,j,ng) = -areafac_l(i,j)
! Subprogram not used 
! Subprogram not used          elseif (yil < c0 .and. xdl < xcl .and. ydl < c0) then
! Subprogram not used 
! Subprogram not used          ! BL (group 1)
! Subprogram not used 
! Subprogram not used             ng = 1
! Subprogram not used             xp    (i,j,1,ng) = xcl
! Subprogram not used             yp    (i,j,1,ng) = ycl
! Subprogram not used             xp    (i,j,2,ng) = xdl
! Subprogram not used             yp    (i,j,2,ng) = ydl
! Subprogram not used             xp    (i,j,3,ng) = xil
! Subprogram not used             yp    (i,j,3,ng) = yil
! Subprogram not used             iflux   (i,j,ng) = i + ishift_bl
! Subprogram not used             jflux   (i,j,ng) = j + jshift_bl
! Subprogram not used             areafact(i,j,ng) = areafac_l(i,j)
! Subprogram not used 
! Subprogram not used          elseif (yil < c0 .and. xdl < xcl .and. ydl >= c0) then
! Subprogram not used 
! Subprogram not used          ! TL1 (group 1)
! Subprogram not used 
! Subprogram not used             ng = 1
! Subprogram not used             xp    (i,j,1,ng) = xcl
! Subprogram not used             yp    (i,j,1,ng) = ycl
! Subprogram not used             xp    (i,j,2,ng) = xdl
! Subprogram not used             yp    (i,j,2,ng) = ydl
! Subprogram not used             xp    (i,j,3,ng) = xic
! Subprogram not used             yp    (i,j,3,ng) = yic
! Subprogram not used             iflux   (i,j,ng) = i + ishift_tl
! Subprogram not used             jflux   (i,j,ng) = j + jshift_tl
! Subprogram not used             areafact(i,j,ng) = areafac_l(i,j)
! Subprogram not used 
! Subprogram not used          ! BL1 (group 3)
! Subprogram not used 
! Subprogram not used             ng = 3
! Subprogram not used             xp    (i,j,1,ng) = xcl
! Subprogram not used             yp    (i,j,1,ng) = ycl
! Subprogram not used             xp    (i,j,2,ng) = xic
! Subprogram not used             yp    (i,j,2,ng) = yic
! Subprogram not used             xp    (i,j,3,ng) = xil
! Subprogram not used             yp    (i,j,3,ng) = yil
! Subprogram not used             iflux   (i,j,ng) = i + ishift_bl
! Subprogram not used             jflux   (i,j,ng) = j + jshift_bl
! Subprogram not used             areafact(i,j,ng) = areafac_l(i,j)
! Subprogram not used 
! Subprogram not used          elseif (yil > c0 .and. xdl < xcl .and. ydl < c0) then
! Subprogram not used 
! Subprogram not used          ! TL2 (group 3)
! Subprogram not used 
! Subprogram not used             ng = 3
! Subprogram not used             xp    (i,j,1,ng) = xcl
! Subprogram not used             yp    (i,j,1,ng) = ycl
! Subprogram not used             xp    (i,j,2,ng) = xil
! Subprogram not used             yp    (i,j,2,ng) = yil
! Subprogram not used             xp    (i,j,3,ng) = xic
! Subprogram not used             yp    (i,j,3,ng) = yic
! Subprogram not used             iflux   (i,j,ng) = i + ishift_tl
! Subprogram not used             jflux   (i,j,ng) = j + jshift_tl
! Subprogram not used             areafact(i,j,ng) = -areafac_l(i,j)
! Subprogram not used 
! Subprogram not used          ! BL2 (group 1)
! Subprogram not used 
! Subprogram not used             ng = 1
! Subprogram not used             xp    (i,j,1,ng) = xcl
! Subprogram not used             yp    (i,j,1,ng) = ycl
! Subprogram not used             xp    (i,j,2,ng) = xic
! Subprogram not used             yp    (i,j,2,ng) = yic
! Subprogram not used             xp    (i,j,3,ng) = xdl
! Subprogram not used             yp    (i,j,3,ng) = ydl
! Subprogram not used             iflux   (i,j,ng) = i + ishift_bl
! Subprogram not used             jflux   (i,j,ng) = j + jshift_bl
! Subprogram not used             areafact(i,j,ng) = -areafac_l(i,j)
! Subprogram not used 
! Subprogram not used          endif                  ! TL and BL triangles
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Locate triangles in TR cell (NE for north edge, SE for east edge)
! Subprogram not used     ! and in BR cell (E for north edge, S for east edge).
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          if (yir > c0 .and. xdr >= xcr .and. ydr >= c0) then
! Subprogram not used 
! Subprogram not used          ! TR (group 2)
! Subprogram not used 
! Subprogram not used             ng = 2
! Subprogram not used             xp    (i,j,1,ng) = xcr
! Subprogram not used             yp    (i,j,1,ng) = ycr
! Subprogram not used             xp    (i,j,2,ng) = xdr
! Subprogram not used             yp    (i,j,2,ng) = ydr
! Subprogram not used             xp    (i,j,3,ng) = xir
! Subprogram not used             yp    (i,j,3,ng) = yir
! Subprogram not used             iflux   (i,j,ng) = i + ishift_tr
! Subprogram not used             jflux   (i,j,ng) = j + jshift_tr
! Subprogram not used             areafact(i,j,ng) = -areafac_r(i,j)
! Subprogram not used 
! Subprogram not used          elseif (yir < c0 .and. xdr >= xcr .and. ydr < c0) then
! Subprogram not used 
! Subprogram not used          ! BR (group 2)
! Subprogram not used 
! Subprogram not used             ng = 2
! Subprogram not used             xp    (i,j,1,ng) = xcr
! Subprogram not used             yp    (i,j,1,ng) = ycr
! Subprogram not used             xp    (i,j,2,ng) = xir
! Subprogram not used             yp    (i,j,2,ng) = yir
! Subprogram not used             xp    (i,j,3,ng) = xdr
! Subprogram not used             yp    (i,j,3,ng) = ydr
! Subprogram not used             iflux   (i,j,ng) = i + ishift_br
! Subprogram not used             jflux   (i,j,ng) = j + jshift_br
! Subprogram not used             areafact(i,j,ng) = areafac_r(i,j)
! Subprogram not used 
! Subprogram not used          elseif (yir < c0 .and. xdr >= xcr  .and. ydr >= c0) then 
! Subprogram not used 
! Subprogram not used          ! TR1 (group 2)
! Subprogram not used 
! Subprogram not used             ng = 2
! Subprogram not used             xp    (i,j,1,ng) = xcr
! Subprogram not used             yp    (i,j,1,ng) = ycr
! Subprogram not used             xp    (i,j,2,ng) = xic
! Subprogram not used             yp    (i,j,2,ng) = yic
! Subprogram not used             xp    (i,j,3,ng) = xdr
! Subprogram not used             yp    (i,j,3,ng) = ydr
! Subprogram not used             iflux   (i,j,ng) = i + ishift_tr
! Subprogram not used             jflux   (i,j,ng) = j + jshift_tr
! Subprogram not used             areafact(i,j,ng) = areafac_r(i,j)
! Subprogram not used 
! Subprogram not used          ! BR1 (group 3)
! Subprogram not used 
! Subprogram not used             ng = 3
! Subprogram not used             xp    (i,j,1,ng) = xcr
! Subprogram not used             yp    (i,j,1,ng) = ycr
! Subprogram not used             xp    (i,j,2,ng) = xir
! Subprogram not used             yp    (i,j,2,ng) = yir
! Subprogram not used             xp    (i,j,3,ng) = xic
! Subprogram not used             yp    (i,j,3,ng) = yic
! Subprogram not used             iflux   (i,j,ng) = i + ishift_br
! Subprogram not used             jflux   (i,j,ng) = j + jshift_br
! Subprogram not used             areafact(i,j,ng) = areafac_r(i,j)
! Subprogram not used 
! Subprogram not used          elseif (yir > c0 .and. xdr >= xcr .and. ydr < c0) then 
! Subprogram not used 
! Subprogram not used          ! TR2 (group 3)
! Subprogram not used 
! Subprogram not used             ng = 3
! Subprogram not used             xp    (i,j,1,ng) = xcr
! Subprogram not used             yp    (i,j,1,ng) = ycr
! Subprogram not used             xp    (i,j,2,ng) = xic
! Subprogram not used             yp    (i,j,2,ng) = yic
! Subprogram not used             xp    (i,j,3,ng) = xir
! Subprogram not used             yp    (i,j,3,ng) = yir
! Subprogram not used             iflux   (i,j,ng) = i + ishift_tr
! Subprogram not used             jflux   (i,j,ng) = j + jshift_tr
! Subprogram not used             areafact(i,j,ng) = -areafac_r(i,j)
! Subprogram not used 
! Subprogram not used          ! BR2 (group 2)
! Subprogram not used 
! Subprogram not used             ng = 2                     
! Subprogram not used             xp    (i,j,1,ng) = xcr
! Subprogram not used             yp    (i,j,1,ng) = ycr
! Subprogram not used             xp    (i,j,2,ng) = xdr
! Subprogram not used             yp    (i,j,2,ng) = ydr
! Subprogram not used             xp    (i,j,3,ng) = xic
! Subprogram not used             yp    (i,j,3,ng) = yic
! Subprogram not used             iflux   (i,j,ng) = i + ishift_br
! Subprogram not used             jflux   (i,j,ng) = j + jshift_br
! Subprogram not used             areafact(i,j,ng) = -areafac_r(i,j)
! Subprogram not used 
! Subprogram not used          endif                  ! TR and BR triangles
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Redefine departure points if not located in central cells (TC or BC)
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          if (xdl < xcl) then
! Subprogram not used             xdl = xil
! Subprogram not used             ydl = yil
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used          if (xdr > xcr) then
! Subprogram not used             xdr = xir
! Subprogram not used             ydr = yir
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! For l_fixed_area = T, shift the midpoint so that the departure
! Subprogram not used     ! region has the prescribed area
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          if (l_fixed_area) then
! Subprogram not used 
! Subprogram not used             ! Sum the areas of the left and right triangles.
! Subprogram not used             ! Note that yp(i,j,1,ng) = 0 for all triangles, so we can
! Subprogram not used             !  drop those terms from the area formula.
! Subprogram not used 
! Subprogram not used             ng = 1
! Subprogram not used             area1 = p5 * ( (xp(i,j,2,ng)-xp(i,j,1,ng)) *   &
! Subprogram not used                             yp(i,j,3,ng)                   &
! Subprogram not used                          -  yp(i,j,2,ng) *                 &
! Subprogram not used                            (xp(i,j,3,ng)-xp(i,j,1,ng)) )   &
! Subprogram not used                          * areafact(i,j,ng) 
! Subprogram not used 
! Subprogram not used             ng = 2
! Subprogram not used             area2 = p5 * ( (xp(i,j,2,ng)-xp(i,j,1,ng)) *   &
! Subprogram not used                             yp(i,j,3,ng)                   &
! Subprogram not used                          -  yp(i,j,2,ng) *                 &
! Subprogram not used                            (xp(i,j,3,ng)-xp(i,j,1,ng)) )   &
! Subprogram not used                          * areafact(i,j,ng) 
! Subprogram not used 
! Subprogram not used             ng = 3
! Subprogram not used             area3 = p5 * ( (xp(i,j,2,ng)-xp(i,j,1,ng)) *   &
! Subprogram not used                             yp(i,j,3,ng)                   &
! Subprogram not used                          -  yp(i,j,2,ng) *                 &
! Subprogram not used                            (xp(i,j,3,ng)-xp(i,j,1,ng)) )   &
! Subprogram not used                          * areafact(i,j,ng) 
! Subprogram not used 
! Subprogram not used             !-----------------------------------------------------------
! Subprogram not used             ! Check whether the central triangles lie in one grid cell or two.
! Subprogram not used             ! If all are in one grid cell, then adjust the area of the central
! Subprogram not used             !  region so that the sum of all triangle areas is equal to the
! Subprogram not used             !  prescribed value.
! Subprogram not used             ! If two triangles are in one grid cell and one is in the other,
! Subprogram not used             !  then compute the area of the lone triangle using an area factor
! Subprogram not used             !  corresponding to the adjacent corner.  This is necessary to prevent
! Subprogram not used             !  negative masses in some rare cases on curved grids.  Then adjust
! Subprogram not used             !  the area of the remaining two-triangle region so that the sum of
! Subprogram not used             !  all triangle areas has the prescribed value.
! Subprogram not used             !-----------------------------------------------------------
! Subprogram not used 
! Subprogram not used             if (ydl*ydr >= c0) then   ! Both DPs lie on same side of x-axis
! Subprogram not used 
! Subprogram not used                ! compute required area of central departure region
! Subprogram not used                area_c  = edgearea(i,j) - area1 - area2 - area3
! Subprogram not used 
! Subprogram not used                ! shift midpoint so that the area of remaining triangles = area_c
! Subprogram not used                w1 = c2*area_c/areafac_c(i,j)    &
! Subprogram not used                     + (xdr-xcl)*ydl + (xcr-xdl)*ydr
! Subprogram not used                w2 = (xdr-xdl)**2 + (ydr-ydl)**2
! Subprogram not used                w1 = w1/w2
! Subprogram not used                xdm = xdm + (ydr - ydl) * w1
! Subprogram not used                ydm = ydm - (xdr - xdl) * w1
! Subprogram not used 
! Subprogram not used                ! compute left and right intersection points
! Subprogram not used                mdl = (ydm - ydl) / (xdm - xdl)
! Subprogram not used                mdr = (ydr - ydm) / (xdr - xdm)
! Subprogram not used 
! Subprogram not used                if (abs(mdl) > puny) then
! Subprogram not used                   xicl = xdl - ydl/mdl
! Subprogram not used                else
! Subprogram not used                   xicl = c0
! Subprogram not used                endif
! Subprogram not used                yicl = c0
! Subprogram not used 
! Subprogram not used                if (abs(mdr) > puny) then
! Subprogram not used                   xicr = xdr - ydr/mdr
! Subprogram not used                else
! Subprogram not used                   xicr = c0
! Subprogram not used                endif
! Subprogram not used                yicr = c0
! Subprogram not used 
! Subprogram not used             elseif (xic < c0) then  ! fix ICL = IC
! Subprogram not used 
! Subprogram not used                xicl = xic
! Subprogram not used                yicl = yic
! Subprogram not used 
! Subprogram not used                ! compute midpoint between ICL and DR
! Subprogram not used                xdm = p5 * (xdr + xicl)
! Subprogram not used                ydm = p5 *  ydr
! Subprogram not used 
! Subprogram not used                ! compute area of triangle adjacent to left corner 
! Subprogram not used                area4 = p5 * (xcl - xic) * ydl * areafac_l(i,j)
! Subprogram not used                area_c  = edgearea(i,j) - area1 - area2 - area3 - area4
! Subprogram not used 
! Subprogram not used                ! shift midpoint so that area of remaining triangles = area_c
! Subprogram not used                w1 = c2*area_c/areafac_c(i,j) + (xcr-xic)*ydr
! Subprogram not used                w2 = (xdr-xic)**2 + ydr**2
! Subprogram not used                w1 = w1/w2
! Subprogram not used                xdm = xdm + ydr*w1
! Subprogram not used                ydm = ydm - (xdr - xic) * w1
! Subprogram not used 
! Subprogram not used                ! compute ICR
! Subprogram not used                mdr = (ydr - ydm) / (xdr - xdm)
! Subprogram not used                if (abs(mdr) > puny) then
! Subprogram not used                   xicr = xdr - ydr/mdr
! Subprogram not used                else
! Subprogram not used                   xicr = c0
! Subprogram not used                endif
! Subprogram not used                yicr = c0
! Subprogram not used 
! Subprogram not used             elseif (xic >= c0) then  ! fix ICR = IR
! Subprogram not used 
! Subprogram not used                xicr = xic
! Subprogram not used                yicr = yic
! Subprogram not used 
! Subprogram not used                ! compute midpoint between ICR and DL 
! Subprogram not used                xdm = p5 * (xicr + xdl)
! Subprogram not used                ydm = p5 *  ydl
! Subprogram not used 
! Subprogram not used                area4 = p5 * (xic - xcr) * ydr * areafac_r(i,j)
! Subprogram not used                area_c  = edgearea(i,j) - area1 - area2 - area3 - area4
! Subprogram not used 
! Subprogram not used                ! shift midpoint so that area of remaining triangles = area_c
! Subprogram not used                w1 = c2*area_c/areafac_c(i,j) + (xic-xcl)*ydl
! Subprogram not used                w2 = (xic-xdl)**2 + ydl**2
! Subprogram not used                w1 = w1/w2
! Subprogram not used                xdm = xdm - ydl*w1
! Subprogram not used                ydm = ydm - (xic - xdl) * w1
! Subprogram not used 
! Subprogram not used                ! compute ICL
! Subprogram not used 
! Subprogram not used                mdl = (ydm - ydl) / (xdm - xdl)
! Subprogram not used                if (abs(mdl) > puny) then
! Subprogram not used                   xicl = xdl - ydl/mdl
! Subprogram not used                else
! Subprogram not used                   xicl = c0
! Subprogram not used                endif
! Subprogram not used                yicl = c0
! Subprogram not used 
! Subprogram not used             endif   ! ydl*ydr >= c0
! Subprogram not used 
! Subprogram not used          endif  ! l_fixed_area
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Locate triangles in BC cell (H for both north and east edges) 
! Subprogram not used     ! and TC cell (N for north edge and E for east edge).
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     ! Start with cases where both DPs lie in the same grid cell
! Subprogram not used 
! Subprogram not used          if (ydl >= c0 .and. ydr >= c0 .and. ydm >= c0) then
! Subprogram not used 
! Subprogram not used          ! TC1a (group 4)
! Subprogram not used 
! Subprogram not used             ng = 4
! Subprogram not used             xp    (i,j,1,ng) = xcl
! Subprogram not used             yp    (i,j,1,ng) = ycl
! Subprogram not used             xp    (i,j,2,ng) = xcr
! Subprogram not used             yp    (i,j,2,ng) = ycr
! Subprogram not used             xp    (i,j,3,ng) = xdl
! Subprogram not used             yp    (i,j,3,ng) = ydl
! Subprogram not used             iflux   (i,j,ng) = i + ishift_tc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_tc
! Subprogram not used             areafact(i,j,ng) = -areafac_c(i,j)
! Subprogram not used 
! Subprogram not used          ! TC2a (group 5)
! Subprogram not used 
! Subprogram not used             ng = 5
! Subprogram not used             xp    (i,j,1,ng) = xcr
! Subprogram not used             yp    (i,j,1,ng) = ycr
! Subprogram not used             xp    (i,j,2,ng) = xdr
! Subprogram not used             yp    (i,j,2,ng) = ydr
! Subprogram not used             xp    (i,j,3,ng) = xdl
! Subprogram not used             yp    (i,j,3,ng) = ydl
! Subprogram not used             iflux   (i,j,ng) = i + ishift_tc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_tc
! Subprogram not used             areafact(i,j,ng) = -areafac_c(i,j)
! Subprogram not used 
! Subprogram not used          ! TC3a (group 6)
! Subprogram not used             ng = 6
! Subprogram not used             xp    (i,j,1,ng) = xdl
! Subprogram not used             yp    (i,j,1,ng) = ydl
! Subprogram not used             xp    (i,j,2,ng) = xdr
! Subprogram not used             yp    (i,j,2,ng) = ydr
! Subprogram not used             xp    (i,j,3,ng) = xdm
! Subprogram not used             yp    (i,j,3,ng) = ydm
! Subprogram not used             iflux   (i,j,ng) = i + ishift_tc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_tc
! Subprogram not used             areafact(i,j,ng) = -areafac_c(i,j)
! Subprogram not used 
! Subprogram not used          elseif (ydl >= c0 .and. ydr >= c0 .and. ydm < c0) then  ! rare
! Subprogram not used 
! Subprogram not used          ! TC1b (group 4)
! Subprogram not used 
! Subprogram not used             ng = 4
! Subprogram not used             xp    (i,j,1,ng) = xcl
! Subprogram not used             yp    (i,j,1,ng) = ycl
! Subprogram not used             xp    (i,j,2,ng) = xicl
! Subprogram not used             yp    (i,j,2,ng) = yicl
! Subprogram not used             xp    (i,j,3,ng) = xdl
! Subprogram not used             yp    (i,j,3,ng) = ydl
! Subprogram not used             iflux   (i,j,ng) = i + ishift_tc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_tc
! Subprogram not used             areafact(i,j,ng) = -areafac_c(i,j)
! Subprogram not used 
! Subprogram not used          ! TC2b (group 5)
! Subprogram not used 
! Subprogram not used             ng = 5
! Subprogram not used             xp    (i,j,1,ng) = xcr
! Subprogram not used             yp    (i,j,1,ng) = ycr
! Subprogram not used             xp    (i,j,2,ng) = xdr
! Subprogram not used             yp    (i,j,2,ng) = ydr
! Subprogram not used             xp    (i,j,3,ng) = xicr
! Subprogram not used             yp    (i,j,3,ng) = yicr
! Subprogram not used             iflux   (i,j,ng) = i + ishift_tc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_tc
! Subprogram not used             areafact(i,j,ng) = -areafac_c(i,j)
! Subprogram not used 
! Subprogram not used          ! BC3b (group 6)
! Subprogram not used 
! Subprogram not used             ng = 6
! Subprogram not used             xp    (i,j,1,ng) = xicr
! Subprogram not used             yp    (i,j,1,ng) = yicr
! Subprogram not used             xp    (i,j,2,ng) = xicl
! Subprogram not used             yp    (i,j,2,ng) = yicl
! Subprogram not used             xp    (i,j,3,ng) = xdm
! Subprogram not used             yp    (i,j,3,ng) = ydm
! Subprogram not used             iflux   (i,j,ng) = i + ishift_bc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_bc
! Subprogram not used             areafact(i,j,ng) = areafac_c(i,j)
! Subprogram not used 
! Subprogram not used          elseif (ydl < c0 .and. ydr < c0 .and. ydm < c0) then
! Subprogram not used 
! Subprogram not used          ! BC1a (group 4)
! Subprogram not used 
! Subprogram not used             ng = 4
! Subprogram not used             xp    (i,j,1,ng) = xcl
! Subprogram not used             yp    (i,j,1,ng) = ycl
! Subprogram not used             xp    (i,j,2,ng) = xdl
! Subprogram not used             yp    (i,j,2,ng) = ydl
! Subprogram not used             xp    (i,j,3,ng) = xcr
! Subprogram not used             yp    (i,j,3,ng) = ycr
! Subprogram not used             iflux   (i,j,ng) = i + ishift_bc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_bc
! Subprogram not used             areafact(i,j,ng) = areafac_c(i,j)
! Subprogram not used 
! Subprogram not used          ! BC2a (group 5)
! Subprogram not used 
! Subprogram not used             ng = 5
! Subprogram not used             xp    (i,j,1,ng) = xcr
! Subprogram not used             yp    (i,j,1,ng) = ycr
! Subprogram not used             xp    (i,j,2,ng) = xdl
! Subprogram not used             yp    (i,j,2,ng) = ydl
! Subprogram not used             xp    (i,j,3,ng) = xdr
! Subprogram not used             yp    (i,j,3,ng) = ydr
! Subprogram not used             iflux   (i,j,ng) = i + ishift_bc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_bc
! Subprogram not used             areafact(i,j,ng) = areafac_c(i,j)
! Subprogram not used 
! Subprogram not used          ! BC3a (group 6)
! Subprogram not used 
! Subprogram not used             ng = 6
! Subprogram not used             xp    (i,j,1,ng) = xdl
! Subprogram not used             yp    (i,j,1,ng) = ydl
! Subprogram not used             xp    (i,j,2,ng) = xdm
! Subprogram not used             yp    (i,j,2,ng) = ydm
! Subprogram not used             xp    (i,j,3,ng) = xdr
! Subprogram not used             yp    (i,j,3,ng) = ydr
! Subprogram not used             iflux   (i,j,ng) = i + ishift_bc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_bc
! Subprogram not used             areafact(i,j,ng) = areafac_c(i,j)
! Subprogram not used 
! Subprogram not used          elseif (ydl < c0 .and. ydr < c0 .and. ydm >= c0) then  ! rare
! Subprogram not used 
! Subprogram not used          ! BC1b (group 4)
! Subprogram not used 
! Subprogram not used             ng = 4
! Subprogram not used             xp    (i,j,1,ng) = xcl
! Subprogram not used             yp    (i,j,1,ng) = ycl
! Subprogram not used             xp    (i,j,2,ng) = xdl
! Subprogram not used             yp    (i,j,2,ng) = ydl
! Subprogram not used             xp    (i,j,3,ng) = xicl
! Subprogram not used             yp    (i,j,3,ng) = yicl
! Subprogram not used             iflux   (i,j,ng) = i + ishift_bc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_bc
! Subprogram not used             areafact(i,j,ng) = areafac_c(i,j)
! Subprogram not used 
! Subprogram not used          ! BC2b (group 5)
! Subprogram not used 
! Subprogram not used             ng = 5
! Subprogram not used             xp    (i,j,1,ng) = xcr
! Subprogram not used             yp    (i,j,1,ng) = ycr
! Subprogram not used             xp    (i,j,2,ng) = xicr
! Subprogram not used             yp    (i,j,2,ng) = yicr
! Subprogram not used             xp    (i,j,3,ng) = xdr
! Subprogram not used             yp    (i,j,3,ng) = ydr
! Subprogram not used             iflux   (i,j,ng) = i + ishift_bc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_bc
! Subprogram not used             areafact(i,j,ng) = areafac_c(i,j)
! Subprogram not used 
! Subprogram not used          ! TC3b (group 6)
! Subprogram not used 
! Subprogram not used             ng = 6
! Subprogram not used             xp    (i,j,1,ng) = xicl
! Subprogram not used             yp    (i,j,1,ng) = yicl
! Subprogram not used             xp    (i,j,2,ng) = xicr
! Subprogram not used             yp    (i,j,2,ng) = yicr
! Subprogram not used             xp    (i,j,3,ng) = xdm
! Subprogram not used             yp    (i,j,3,ng) = ydm
! Subprogram not used             iflux   (i,j,ng) = i + ishift_tc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_tc
! Subprogram not used             areafact(i,j,ng) = -areafac_c(i,j)
! Subprogram not used 
! Subprogram not used     ! Now consider cases where the two DPs lie in different grid cells
! Subprogram not used     ! For these cases, one triangle is given the area factor associated
! Subprogram not used     !  with the adjacent corner, to avoid rare negative masses on curved grids.
! Subprogram not used 
! Subprogram not used          elseif (ydl >= c0 .and. ydr < c0 .and. xic >= c0  &
! Subprogram not used                                           .and. ydm >= c0) then
! Subprogram not used 
! Subprogram not used          ! TC1b (group 4)
! Subprogram not used 
! Subprogram not used             ng = 4
! Subprogram not used             xp    (i,j,1,ng) = xcl
! Subprogram not used             yp    (i,j,1,ng) = ycl
! Subprogram not used             xp    (i,j,2,ng) = xicr
! Subprogram not used             yp    (i,j,2,ng) = yicr
! Subprogram not used             xp    (i,j,3,ng) = xdl
! Subprogram not used             yp    (i,j,3,ng) = ydl
! Subprogram not used             iflux   (i,j,ng) = i + ishift_tc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_tc
! Subprogram not used             areafact(i,j,ng) = -areafac_c(i,j)
! Subprogram not used 
! Subprogram not used          ! BC2b (group 5)
! Subprogram not used 
! Subprogram not used             ng = 5
! Subprogram not used             xp    (i,j,1,ng) = xcr
! Subprogram not used             yp    (i,j,1,ng) = ycr
! Subprogram not used             xp    (i,j,2,ng) = xicr
! Subprogram not used             yp    (i,j,2,ng) = yicr
! Subprogram not used             xp    (i,j,3,ng) = xdr
! Subprogram not used             yp    (i,j,3,ng) = ydr
! Subprogram not used             iflux   (i,j,ng) = i + ishift_bc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_bc
! Subprogram not used             areafact(i,j,ng) = areafac_r(i,j)
! Subprogram not used 
! Subprogram not used          ! TC3b (group 6)
! Subprogram not used 
! Subprogram not used             ng = 6
! Subprogram not used             xp    (i,j,1,ng) = xdl
! Subprogram not used             yp    (i,j,1,ng) = ydl
! Subprogram not used             xp    (i,j,2,ng) = xicr
! Subprogram not used             yp    (i,j,2,ng) = yicr
! Subprogram not used             xp    (i,j,3,ng) = xdm
! Subprogram not used             yp    (i,j,3,ng) = ydm
! Subprogram not used             iflux   (i,j,ng) = i + ishift_tc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_tc
! Subprogram not used             areafact(i,j,ng) = -areafac_c(i,j)
! Subprogram not used 
! Subprogram not used          elseif (ydl >= c0 .and. ydr < c0 .and. xic >= c0  &
! Subprogram not used                                           .and. ydm < c0 ) then  ! less common
! Subprogram not used 
! Subprogram not used          ! TC1b (group 4)
! Subprogram not used 
! Subprogram not used             ng = 4
! Subprogram not used             xp    (i,j,1,ng) = xcl
! Subprogram not used             yp    (i,j,1,ng) = ycl
! Subprogram not used             xp    (i,j,2,ng) = xicl
! Subprogram not used             yp    (i,j,2,ng) = yicl
! Subprogram not used             xp    (i,j,3,ng) = xdl
! Subprogram not used             yp    (i,j,3,ng) = ydl
! Subprogram not used             iflux   (i,j,ng) = i + ishift_tc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_tc
! Subprogram not used             areafact(i,j,ng) = -areafac_c(i,j)
! Subprogram not used 
! Subprogram not used          ! BC2b (group 5)
! Subprogram not used 
! Subprogram not used             ng = 5
! Subprogram not used             xp    (i,j,1,ng) = xcr
! Subprogram not used             yp    (i,j,1,ng) = ycr
! Subprogram not used             xp    (i,j,2,ng) = xicr
! Subprogram not used             yp    (i,j,2,ng) = yicr
! Subprogram not used             xp    (i,j,3,ng) = xdr
! Subprogram not used             yp    (i,j,3,ng) = ydr
! Subprogram not used             iflux   (i,j,ng) = i + ishift_bc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_bc
! Subprogram not used             areafact(i,j,ng) = areafac_r(i,j)
! Subprogram not used 
! Subprogram not used          ! BC3b (group 6)
! Subprogram not used 
! Subprogram not used             ng = 6
! Subprogram not used             xp    (i,j,1,ng) = xicr
! Subprogram not used             yp    (i,j,1,ng) = yicr
! Subprogram not used             xp    (i,j,2,ng) = xicl
! Subprogram not used             yp    (i,j,2,ng) = yicl
! Subprogram not used             xp    (i,j,3,ng) = xdm
! Subprogram not used             yp    (i,j,3,ng) = ydm
! Subprogram not used             iflux   (i,j,ng) = i + ishift_bc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_bc
! Subprogram not used             areafact(i,j,ng) = areafac_c(i,j)
! Subprogram not used 
! Subprogram not used          elseif (ydl >= c0 .and. ydr < c0 .and. xic < c0   &
! Subprogram not used                                           .and. ydm < c0) then
! Subprogram not used 
! Subprogram not used          ! TC1b (group 4)
! Subprogram not used 
! Subprogram not used             ng = 4
! Subprogram not used             xp    (i,j,1,ng) = xcl
! Subprogram not used             yp    (i,j,1,ng) = ycl
! Subprogram not used             xp    (i,j,2,ng) = xicl
! Subprogram not used             yp    (i,j,2,ng) = yicl
! Subprogram not used             xp    (i,j,3,ng) = xdl
! Subprogram not used             yp    (i,j,3,ng) = ydl
! Subprogram not used             iflux   (i,j,ng) = i + ishift_tc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_tc
! Subprogram not used             areafact(i,j,ng) = -areafac_l(i,j)
! Subprogram not used 
! Subprogram not used          ! BC2b (group 5)
! Subprogram not used 
! Subprogram not used             ng = 5
! Subprogram not used             xp    (i,j,1,ng) = xcr
! Subprogram not used             yp    (i,j,1,ng) = ycr
! Subprogram not used             xp    (i,j,2,ng) = xicl
! Subprogram not used             yp    (i,j,2,ng) = yicl
! Subprogram not used             xp    (i,j,3,ng) = xdr
! Subprogram not used             yp    (i,j,3,ng) = ydr
! Subprogram not used             iflux   (i,j,ng) = i + ishift_bc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_bc
! Subprogram not used             areafact(i,j,ng) = areafac_c(i,j)
! Subprogram not used 
! Subprogram not used          ! BC3b (group 6)
! Subprogram not used 
! Subprogram not used             ng = 6
! Subprogram not used             xp    (i,j,1,ng) = xdr
! Subprogram not used             yp    (i,j,1,ng) = ydr
! Subprogram not used             xp    (i,j,2,ng) = xicl
! Subprogram not used             yp    (i,j,2,ng) = yicl
! Subprogram not used             xp    (i,j,3,ng) = xdm
! Subprogram not used             yp    (i,j,3,ng) = ydm
! Subprogram not used             iflux   (i,j,ng) = i + ishift_bc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_bc
! Subprogram not used             areafact(i,j,ng) = areafac_c(i,j)
! Subprogram not used 
! Subprogram not used          elseif (ydl >= c0 .and. ydr < c0 .and. xic <  c0  &
! Subprogram not used                                           .and. ydm >= c0) then  ! less common
! Subprogram not used 
! Subprogram not used          ! TC1b (group 4)
! Subprogram not used 
! Subprogram not used             ng = 4
! Subprogram not used             xp    (i,j,1,ng) = xcl
! Subprogram not used             yp    (i,j,1,ng) = ycl
! Subprogram not used             xp    (i,j,2,ng) = xicl
! Subprogram not used             yp    (i,j,2,ng) = yicl
! Subprogram not used             xp    (i,j,3,ng) = xdl
! Subprogram not used             yp    (i,j,3,ng) = ydl
! Subprogram not used             iflux   (i,j,ng) = i + ishift_tc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_tc
! Subprogram not used             areafact(i,j,ng) = -areafac_l(i,j)
! Subprogram not used 
! Subprogram not used          ! BC2b (group 5)
! Subprogram not used 
! Subprogram not used             ng = 5
! Subprogram not used             xp    (i,j,1,ng) = xcr
! Subprogram not used             yp    (i,j,1,ng) = ycr
! Subprogram not used             xp    (i,j,2,ng) = xicr
! Subprogram not used             yp    (i,j,2,ng) = yicr
! Subprogram not used             xp    (i,j,3,ng) = xdr
! Subprogram not used             yp    (i,j,3,ng) = ydr
! Subprogram not used             iflux   (i,j,ng) = i + ishift_bc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_bc
! Subprogram not used             areafact(i,j,ng) = areafac_c(i,j)
! Subprogram not used 
! Subprogram not used          ! TC3b (group 6)
! Subprogram not used 
! Subprogram not used             ng = 6
! Subprogram not used             xp    (i,j,1,ng) = xicl
! Subprogram not used             yp    (i,j,1,ng) = yicl
! Subprogram not used             xp    (i,j,2,ng) = xicr
! Subprogram not used             yp    (i,j,2,ng) = yicr
! Subprogram not used             xp    (i,j,3,ng) = xdm
! Subprogram not used             yp    (i,j,3,ng) = ydm
! Subprogram not used             iflux   (i,j,ng) = i + ishift_tc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_tc
! Subprogram not used             areafact(i,j,ng) = -areafac_c(i,j)
! Subprogram not used 
! Subprogram not used          elseif (ydl < c0 .and. ydr >= c0 .and. xic <  c0  &
! Subprogram not used                                           .and. ydm >= c0) then
! Subprogram not used 
! Subprogram not used          ! BC1b (group 4)
! Subprogram not used 
! Subprogram not used             ng = 4
! Subprogram not used             xp    (i,j,1,ng) = xcl
! Subprogram not used             yp    (i,j,1,ng) = ycl
! Subprogram not used             xp    (i,j,2,ng) = xdl
! Subprogram not used             yp    (i,j,2,ng) = ydl
! Subprogram not used             xp    (i,j,3,ng) = xicl
! Subprogram not used             yp    (i,j,3,ng) = yicl
! Subprogram not used             iflux   (i,j,ng) = i + ishift_bc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_bc
! Subprogram not used             areafact(i,j,ng) = areafac_l(i,j)
! Subprogram not used 
! Subprogram not used          ! TC2b (group 5)
! Subprogram not used 
! Subprogram not used             ng = 5
! Subprogram not used             xp    (i,j,1,ng) = xcr
! Subprogram not used             yp    (i,j,1,ng) = ycr
! Subprogram not used             xp    (i,j,2,ng) = xdr
! Subprogram not used             yp    (i,j,2,ng) = ydr
! Subprogram not used             xp    (i,j,3,ng) = xicl
! Subprogram not used             yp    (i,j,3,ng) = yicl
! Subprogram not used             iflux   (i,j,ng) = i + ishift_tc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_tc
! Subprogram not used             areafact(i,j,ng) = -areafac_c(i,j)
! Subprogram not used 
! Subprogram not used          ! TC3b (group 6)
! Subprogram not used 
! Subprogram not used             ng = 6
! Subprogram not used             xp    (i,j,1,ng) = xicl
! Subprogram not used             yp    (i,j,1,ng) = yicl
! Subprogram not used             xp    (i,j,2,ng) = xdr
! Subprogram not used             yp    (i,j,2,ng) = ydr
! Subprogram not used             xp    (i,j,3,ng) = xdm
! Subprogram not used             yp    (i,j,3,ng) = ydm
! Subprogram not used             iflux   (i,j,ng) = i + ishift_tc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_tc
! Subprogram not used             areafact(i,j,ng) = -areafac_c(i,j)
! Subprogram not used 
! Subprogram not used          elseif (ydl < c0 .and. ydr >= c0 .and. xic < c0  &
! Subprogram not used                                           .and. ydm < c0) then ! less common
! Subprogram not used 
! Subprogram not used          ! BC1b (group 4)
! Subprogram not used 
! Subprogram not used             ng = 4
! Subprogram not used             xp    (i,j,1,ng) = xcl
! Subprogram not used             yp    (i,j,1,ng) = ycl
! Subprogram not used             xp    (i,j,2,ng) = xdl
! Subprogram not used             yp    (i,j,2,ng) = ydl
! Subprogram not used             xp    (i,j,3,ng) = xicl
! Subprogram not used             yp    (i,j,3,ng) = yicl
! Subprogram not used             iflux   (i,j,ng) = i + ishift_bc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_bc
! Subprogram not used             areafact(i,j,ng) = areafac_l(i,j)
! Subprogram not used 
! Subprogram not used          ! TC2b (group 5)
! Subprogram not used 
! Subprogram not used             ng = 5
! Subprogram not used             xp    (i,j,1,ng) = xcr
! Subprogram not used             yp    (i,j,1,ng) = ycr
! Subprogram not used             xp    (i,j,2,ng) = xdr
! Subprogram not used             yp    (i,j,2,ng) = ydr
! Subprogram not used             xp    (i,j,3,ng) = xicr
! Subprogram not used             yp    (i,j,3,ng) = yicr
! Subprogram not used             iflux   (i,j,ng) = i + ishift_tc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_tc
! Subprogram not used             areafact(i,j,ng) = -areafac_c(i,j)
! Subprogram not used 
! Subprogram not used          ! BC3b (group 6)
! Subprogram not used 
! Subprogram not used             ng = 6
! Subprogram not used             xp    (i,j,1,ng) = xicr
! Subprogram not used             yp    (i,j,1,ng) = yicr
! Subprogram not used             xp    (i,j,2,ng) = xicl
! Subprogram not used             yp    (i,j,2,ng) = yicl
! Subprogram not used             xp    (i,j,3,ng) = xdm
! Subprogram not used             yp    (i,j,3,ng) = ydm
! Subprogram not used             iflux   (i,j,ng) = i + ishift_bc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_bc
! Subprogram not used             areafact(i,j,ng) = areafac_c(i,j)
! Subprogram not used 
! Subprogram not used          elseif (ydl < c0 .and. ydr >= c0 .and. xic >= c0  &
! Subprogram not used                                           .and. ydm <  c0) then
! Subprogram not used 
! Subprogram not used          ! BC1b (group 4)
! Subprogram not used 
! Subprogram not used             ng = 4
! Subprogram not used             xp    (i,j,1,ng) = xcl
! Subprogram not used             yp    (i,j,1,ng) = ycl
! Subprogram not used             xp    (i,j,2,ng) = xdl
! Subprogram not used             yp    (i,j,2,ng) = ydl
! Subprogram not used             xp    (i,j,3,ng) = xicr
! Subprogram not used             yp    (i,j,3,ng) = yicr
! Subprogram not used             iflux   (i,j,ng) = i + ishift_bc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_bc
! Subprogram not used             areafact(i,j,ng) = areafac_c(i,j)
! Subprogram not used 
! Subprogram not used          ! TC2b (group 5)
! Subprogram not used 
! Subprogram not used             ng = 5
! Subprogram not used             xp    (i,j,1,ng) = xcr
! Subprogram not used             yp    (i,j,1,ng) = ycr
! Subprogram not used             xp    (i,j,2,ng) = xdr
! Subprogram not used             yp    (i,j,2,ng) = ydr
! Subprogram not used             xp    (i,j,3,ng) = xicr
! Subprogram not used             yp    (i,j,3,ng) = yicr
! Subprogram not used             iflux   (i,j,ng) = i + ishift_tc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_tc
! Subprogram not used             areafact(i,j,ng) = -areafac_r(i,j)
! Subprogram not used 
! Subprogram not used          ! BC3b (group 6)
! Subprogram not used 
! Subprogram not used             ng = 6
! Subprogram not used             xp    (i,j,1,ng) = xicr
! Subprogram not used             yp    (i,j,1,ng) = yicr
! Subprogram not used             xp    (i,j,2,ng) = xdl
! Subprogram not used             yp    (i,j,2,ng) = ydl
! Subprogram not used             xp    (i,j,3,ng) = xdm
! Subprogram not used             yp    (i,j,3,ng) = ydm
! Subprogram not used             iflux   (i,j,ng) = i + ishift_bc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_bc
! Subprogram not used             areafact(i,j,ng) = areafac_c(i,j)
! Subprogram not used 
! Subprogram not used          elseif (ydl < c0 .and. ydr >= c0 .and. xic >= c0   &
! Subprogram not used                                           .and. ydm >= c0) then  ! less common
! Subprogram not used 
! Subprogram not used          ! BC1b (group 4)
! Subprogram not used 
! Subprogram not used             ng = 4
! Subprogram not used             xp    (i,j,1,ng) = xcl
! Subprogram not used             yp    (i,j,1,ng) = ycl
! Subprogram not used             xp    (i,j,2,ng) = xdl
! Subprogram not used             yp    (i,j,2,ng) = ydl
! Subprogram not used             xp    (i,j,3,ng) = xicl
! Subprogram not used             yp    (i,j,3,ng) = yicl
! Subprogram not used             iflux   (i,j,ng) = i + ishift_bc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_bc
! Subprogram not used             areafact(i,j,ng) = areafac_c(i,j)
! Subprogram not used 
! Subprogram not used          ! TC2b (group 5)
! Subprogram not used 
! Subprogram not used             ng = 5
! Subprogram not used             xp    (i,j,1,ng) = xcr
! Subprogram not used             yp    (i,j,1,ng) = ycr
! Subprogram not used             xp    (i,j,2,ng) = xdr
! Subprogram not used             yp    (i,j,2,ng) = ydr
! Subprogram not used             xp    (i,j,3,ng) = xicr
! Subprogram not used             yp    (i,j,3,ng) = yicr
! Subprogram not used             iflux   (i,j,ng) = i + ishift_tc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_tc
! Subprogram not used             areafact(i,j,ng) = -areafac_r(i,j)
! Subprogram not used 
! Subprogram not used          ! TC3b (group 6)
! Subprogram not used 
! Subprogram not used             ng = 6
! Subprogram not used             xp    (i,j,1,ng) = xicl
! Subprogram not used             yp    (i,j,1,ng) = yicl
! Subprogram not used             xp    (i,j,2,ng) = xicr
! Subprogram not used             yp    (i,j,2,ng) = yicr
! Subprogram not used             xp    (i,j,3,ng) = xdm
! Subprogram not used             yp    (i,j,3,ng) = ydm
! Subprogram not used             iflux   (i,j,ng) = i + ishift_tc
! Subprogram not used             jflux   (i,j,ng) = j + jshift_tc
! Subprogram not used             areafact(i,j,ng) = -areafac_c(i,j)
! Subprogram not used 
! Subprogram not used          endif                  ! TC and BC triangles
! Subprogram not used 
! Subprogram not used       enddo                     ! ij
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Compute triangle areas with appropriate sign.
! Subprogram not used     ! These are found by computing the area in scaled coordinates and
! Subprogram not used     !  multiplying by a scale factor (areafact).
! Subprogram not used     ! Note that the scale factor is positive for fluxes out of the cell 
! Subprogram not used     !  and negative for fluxes into the cell.
! Subprogram not used     !
! Subprogram not used     ! Note: The triangle area formula below gives A >=0 iff the triangle
! Subprogram not used     !        points x1, x2, and x3 are taken in counterclockwise order.
! Subprogram not used     !       These points are defined above in such a way that the
! Subprogram not used     !        order is nearly always CCW.
! Subprogram not used     !       In rare cases, we may compute A < 0.  In this case,
! Subprogram not used     !        the quadrilateral departure area is equal to the 
! Subprogram not used     !        difference of two triangle areas instead of the sum.
! Subprogram not used     !        The fluxes work out correctly in the end.
! Subprogram not used     !
! Subprogram not used     ! Also compute the cumulative area transported across each edge.
! Subprogram not used     ! If l_fixed_area = T, this area is compared to edgearea as a bug check.
! Subprogram not used     ! If l_fixed_area = F, this area is passed as an output array.
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       areasum(:,:) = c0
! Subprogram not used 
! Subprogram not used       do ng = 1, ngroups
! Subprogram not used          icells(ng) = 0
! Subprogram not used 
! Subprogram not used          do ij = 1, icellsd
! Subprogram not used             i = indxid(ij)
! Subprogram not used             j = indxjd(ij)
! Subprogram not used 
! Subprogram not used             triarea(i,j,ng) = p5 * ( (xp(i,j,2,ng)-xp(i,j,1,ng)) *   &
! Subprogram not used                                      (yp(i,j,3,ng)-yp(i,j,1,ng))   &
! Subprogram not used                                    - (yp(i,j,2,ng)-yp(i,j,1,ng)) *   &
! Subprogram not used                                      (xp(i,j,3,ng)-xp(i,j,1,ng)) )   &
! Subprogram not used                                    * areafact(i,j,ng) 
! Subprogram not used 
! Subprogram not used             if (abs(triarea(i,j,ng)) < eps16*areafac_c(i,j)) then
! Subprogram not used                triarea(i,j,ng) = c0
! Subprogram not used             else
! Subprogram not used                icells(ng) = icells(ng) + 1 
! Subprogram not used                ic = icells(ng)
! Subprogram not used                indxi(ic,ng) = i
! Subprogram not used                indxj(ic,ng) = j
! Subprogram not used             endif
! Subprogram not used 
! Subprogram not used             areasum(i,j) = areasum(i,j) + triarea(i,j,ng)
! Subprogram not used 
! Subprogram not used          enddo                  ! ij
! Subprogram not used       enddo                     ! ng
! Subprogram not used 
! Subprogram not used       if (l_fixed_area) then
! Subprogram not used        if (bugcheck) then   ! set bugcheck = F to speed up code
! Subprogram not used          do ij = 1, icellsd
! Subprogram not used             i = indxid(ij)
! Subprogram not used             j = indxjd(ij)
! Subprogram not used             if (abs(areasum(i,j) - edgearea(i,j)) > eps13*areafac_c(i,j)) then
! Subprogram not used                print*, ''
! Subprogram not used                print*, 'Areas do not add up: m, i, j, edge =',   &
! Subprogram not used                         my_task, i, j, trim(edge)
! Subprogram not used                print*, 'edgearea =', edgearea(i,j)
! Subprogram not used                print*, 'areasum =', areasum(i,j)
! Subprogram not used                print*, 'areafac_c =', areafac_c(i,j)
! Subprogram not used                print*, ''
! Subprogram not used                print*, 'Triangle areas:'
! Subprogram not used                do ng = 1, ngroups   ! not vector friendly
! Subprogram not used                   if (abs(triarea(i,j,ng)) > eps16*abs(areafact(i,j,ng))) then
! Subprogram not used                      print*, ng, triarea(i,j,ng)
! Subprogram not used                   endif
! Subprogram not used                enddo
! Subprogram not used             endif
! Subprogram not used          enddo
! Subprogram not used        endif          ! bugcheck
! Subprogram not used 
! Subprogram not used       else            ! l_fixed_area = F
! Subprogram not used          do ij = 1, icellsd
! Subprogram not used             i = indxid(ij)
! Subprogram not used             j = indxjd(ij)
! Subprogram not used             edgearea(i,j) = areasum(i,j)
! Subprogram not used          enddo
! Subprogram not used       endif     ! l_fixed_area
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Transform triangle vertices to a scaled coordinate system centered
! Subprogram not used     !  in the cell containing the triangle.
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       if (trim(edge) == 'north') then
! Subprogram not used          do ng = 1, ngroups
! Subprogram not used             do nv = 1, nvert
! Subprogram not used                do ij = 1, icells(ng)
! Subprogram not used                   i = indxi(ij,ng)
! Subprogram not used                   j = indxj(ij,ng)
! Subprogram not used                   ishift = iflux(i,j,ng) - i
! Subprogram not used                   jshift = jflux(i,j,ng) - j
! Subprogram not used                   xp(i,j,nv,ng) = xp(i,j,nv,ng) - c1*ishift
! Subprogram not used                   yp(i,j,nv,ng) = yp(i,j,nv,ng) + p5 - c1*jshift
! Subprogram not used                enddo            ! ij
! Subprogram not used             enddo               ! nv
! Subprogram not used          enddo                  ! ng
! Subprogram not used       else                      ! east edge
! Subprogram not used          do ng = 1, ngroups
! Subprogram not used             do nv = 1, nvert
! Subprogram not used                do ij = 1, icells(ng)
! Subprogram not used                   i = indxi(ij,ng)
! Subprogram not used                   j = indxj(ij,ng)
! Subprogram not used                   ishift = iflux(i,j,ng) - i
! Subprogram not used                   jshift = jflux(i,j,ng) - j
! Subprogram not used                   ! Note rotation of pi/2 here
! Subprogram not used                   w1 = xp(i,j,nv,ng)
! Subprogram not used                   xp(i,j,nv,ng) =  yp(i,j,nv,ng) + p5 - c1*ishift
! Subprogram not used                   yp(i,j,nv,ng) = -w1 - c1*jshift
! Subprogram not used                enddo            ! ij
! Subprogram not used             enddo               ! nv
! Subprogram not used          enddo                  ! ng
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       if (bugcheck) then
! Subprogram not used          do ng = 1, ngroups
! Subprogram not used          do nv = 1, nvert
! Subprogram not used             do j = jb, je
! Subprogram not used             do i = ib, ie
! Subprogram not used                if (abs(triarea(i,j,ng)) > puny) then
! Subprogram not used                   if (abs(xp(i,j,nv,ng)) > p5+puny) then
! Subprogram not used                      print*, ''
! Subprogram not used                      print*, 'WARNING: xp =', xp(i,j,nv,ng)
! Subprogram not used                      print*, 'm, i, j, ng, nv =', my_task, i, j, ng, nv
! Subprogram not used !                     print*, 'yil,xdl,xcl,ydl=',yil,xdl,xcl,ydl
! Subprogram not used !                     print*, 'yir,xdr,xcr,ydr=',yir,xdr,xcr,ydr
! Subprogram not used !                     print*, 'ydm=',ydm
! Subprogram not used !                      stop
! Subprogram not used                   endif
! Subprogram not used                   if (abs(yp(i,j,nv,ng)) > p5+puny) then
! Subprogram not used                      print*, ''
! Subprogram not used                      print*, 'WARNING: yp =', yp(i,j,nv,ng)
! Subprogram not used                      print*, 'm, i, j, ng, nv =', my_task, i, j, ng, nv
! Subprogram not used                   endif
! Subprogram not used                endif   ! triarea
! Subprogram not used             enddo
! Subprogram not used             enddo
! Subprogram not used          enddo
! Subprogram not used          enddo
! Subprogram not used       endif  ! bugcheck
! Subprogram not used 
! Subprogram not used       end subroutine locate_triangles

!=======================================================================
!
!BOP
! !IROUTINE: triangle_coordinates - find coordinates of quadrature points
!
! !INTERFACE:
!
! Subprogram not used       subroutine triangle_coordinates (nx_block,       ny_block,  &
! Subprogram not used                                        integral_order, icells,    &
! Subprogram not used                                        indxi,          indxj,     &
! Subprogram not used                                        xp,             yp)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! For each triangle, find the coordinates of the quadrature points needed
! Subprogram not used !  to compute integrals of linear, quadratic, or cubic polynomials,
! Subprogram not used !  using formulas from A.H. Stroud, Approximate Calculation of Multiple
! Subprogram not used !  Integrals, Prentice-Hall, 1971.  (Section 8.8, formula 3.1.)
! Subprogram not used ! Linear functions can be integrated exactly by evaluating the function 
! Subprogram not used !  at just one point (the midpoint).  Quadratic functions require
! Subprogram not used !  3 points, and cubics require 4 points.
! Subprogram not used ! The default is cubic, but the code can be sped up slightly using 
! Subprogram not used !  linear or quadratic integrals, usually with little loss of accuracy.
! Subprogram not used !
! Subprogram not used ! The formulas are as follows:
! Subprogram not used !
! Subprogram not used ! I1 = integral of f(x,y)*dA
! Subprogram not used !    = A * f(x0,y0)
! Subprogram not used ! where A is the traingle area and (x0,y0) is the midpoint.
! Subprogram not used !
! Subprogram not used ! I2 = A * (f(x1,y1) + f(x2,y2) + f(x3,y3))
! Subprogram not used ! where these three points are located halfway between the midpoint
! Subprogram not used ! and the three vertics of the triangle.
! Subprogram not used !
! Subprogram not used ! I3 = A * [ -9/16 *  f(x0,y0)
! Subprogram not used !           + 25/48 * (f(x1,y1) + f(x2,y2) + f(x3,y3))]
! Subprogram not used ! where (x0,y0) is the midpoint, and the other three points are
! Subprogram not used ! located 2/5 of the way from the midpoint to the three vertices.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author William H. Lipscomb, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) ::   &
! Subprogram not used            nx_block, ny_block,&! block dimensions
! Subprogram not used            integral_order      ! polynomial order for quadrature integrals 
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension (ngroups), intent(in) ::     &
! Subprogram not used            icells              ! number of cells where triarea > puny
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension (nx_block*ny_block,ngroups),     &
! Subprogram not used            intent(in) ::     &
! Subprogram not used            indxi ,&! compressed index in i-direction
! Subprogram not used            indxj   ! compressed index in j-direction
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), intent(inout),   &
! Subprogram not used            dimension (nx_block, ny_block, 0:nvert, ngroups) ::   &
! Subprogram not used            xp, yp          ! coordinates of triangle points
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) ::   &
! Subprogram not used            i, j, ij          ,&! horizontal indices
! Subprogram not used            ng                  ! triangle index
! Subprogram not used 
! Subprogram not used 
! Subprogram not used       if (integral_order == 1) then ! linear (1-point formula)
! Subprogram not used 
! Subprogram not used          do ng = 1, ngroups
! Subprogram not used          do ij = 1, icells(ng)
! Subprogram not used             i = indxi(ij,ng)
! Subprogram not used             j = indxj(ij,ng)
! Subprogram not used 
! Subprogram not used             ! coordinates of midpoint
! Subprogram not used             xp(i,j,0,ng) = p333   &
! Subprogram not used                         * (xp(i,j,1,ng) + xp(i,j,2,ng) + xp(i,j,3,ng))
! Subprogram not used             yp(i,j,0,ng) = p333   &
! Subprogram not used                         * (yp(i,j,1,ng) + yp(i,j,2,ng) + yp(i,j,3,ng))
! Subprogram not used 
! Subprogram not used          enddo                  ! ij
! Subprogram not used          enddo                  ! ng
! Subprogram not used 
! Subprogram not used       elseif (integral_order == 2) then ! quadratic (3-point formula)
! Subprogram not used 
! Subprogram not used !DIR$ CONCURRENT !Cray
! Subprogram not used !cdir nodep      !NEC
! Subprogram not used !ocl novrec      !Fujitsu
! Subprogram not used 
! Subprogram not used          do ng = 1, ngroups
! Subprogram not used          do ij = 1, icells(ng)
! Subprogram not used             i = indxi(ij,ng)
! Subprogram not used             j = indxj(ij,ng)
! Subprogram not used 
! Subprogram not used             ! coordinates of midpoint
! Subprogram not used             xp(i,j,0,ng) = p333   &
! Subprogram not used                         * (xp(i,j,1,ng) + xp(i,j,2,ng) + xp(i,j,3,ng))
! Subprogram not used             yp(i,j,0,ng) = p333   &
! Subprogram not used                         * (yp(i,j,1,ng) + yp(i,j,2,ng) + yp(i,j,3,ng))
! Subprogram not used 
! Subprogram not used             ! coordinates of the 3 points needed for integrals
! Subprogram not used 
! Subprogram not used             xp(i,j,1,ng) = p5*xp(i,j,1,ng) + p5*xp(i,j,0,ng)
! Subprogram not used             yp(i,j,1,ng) = p5*yp(i,j,1,ng) + p5*yp(i,j,0,ng)
! Subprogram not used 
! Subprogram not used             xp(i,j,2,ng) = p5*xp(i,j,2,ng) + p5*xp(i,j,0,ng)
! Subprogram not used             yp(i,j,2,ng) = p5*yp(i,j,2,ng) + p5*yp(i,j,0,ng)
! Subprogram not used 
! Subprogram not used             xp(i,j,3,ng) = p5*xp(i,j,3,ng) + p5*xp(i,j,0,ng)
! Subprogram not used             yp(i,j,3,ng) = p5*yp(i,j,3,ng) + p5*yp(i,j,0,ng)
! Subprogram not used 
! Subprogram not used          enddo                  ! ij
! Subprogram not used          enddo                  ! ng
! Subprogram not used 
! Subprogram not used       else                      ! cubic (4-point formula)
! Subprogram not used 
! Subprogram not used !DIR$ CONCURRENT !Cray
! Subprogram not used !cdir nodep      !NEC
! Subprogram not used !ocl novrec      !Fujitsu
! Subprogram not used          do ng = 1, ngroups
! Subprogram not used          do ij = 1, icells(ng)
! Subprogram not used             i = indxi(ij,ng)
! Subprogram not used             j = indxj(ij,ng)
! Subprogram not used 
! Subprogram not used             ! coordinates of midpoint
! Subprogram not used             xp(i,j,0,ng) = p333   &
! Subprogram not used                         * (xp(i,j,1,ng) + xp(i,j,2,ng) + xp(i,j,3,ng))
! Subprogram not used             yp(i,j,0,ng) = p333   &
! Subprogram not used                         * (yp(i,j,1,ng) + yp(i,j,2,ng) + yp(i,j,3,ng))
! Subprogram not used 
! Subprogram not used             ! coordinates of the other 3 points needed for integrals
! Subprogram not used 
! Subprogram not used             xp(i,j,1,ng) = p4*xp(i,j,1,ng) + p6*xp(i,j,0,ng)
! Subprogram not used             yp(i,j,1,ng) = p4*yp(i,j,1,ng) + p6*yp(i,j,0,ng)
! Subprogram not used 
! Subprogram not used             xp(i,j,2,ng) = p4*xp(i,j,2,ng) + p6*xp(i,j,0,ng)
! Subprogram not used             yp(i,j,2,ng) = p4*yp(i,j,2,ng) + p6*yp(i,j,0,ng)
! Subprogram not used                
! Subprogram not used             xp(i,j,3,ng) = p4*xp(i,j,3,ng) + p6*xp(i,j,0,ng)
! Subprogram not used             yp(i,j,3,ng) = p4*yp(i,j,3,ng) + p6*yp(i,j,0,ng)
! Subprogram not used                
! Subprogram not used          enddo                  ! ij
! Subprogram not used          enddo                  ! ng
! Subprogram not used 
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       end subroutine triangle_coordinates

!=======================================================================
!
!BOP
!
! !IROUTINE: transport_integrals - compute transports across each edge
!
! !INTERFACE:
!
! Subprogram not used       subroutine transport_integrals (nx_block,       ny_block,    &
! Subprogram not used                                       ntrace,         icells,      &
! Subprogram not used                                       indxi,          indxj,       &
! Subprogram not used                                       tracer_type,    depend,      &
! Subprogram not used                                       integral_order, triarea,     &
! Subprogram not used                                       iflux,          jflux,       &
! Subprogram not used                                       xp,             yp,          &
! Subprogram not used                                       mc,             mx,          &
! Subprogram not used                                       my,             mflx,       &
! Subprogram not used                                       tc,             tx,          &
! Subprogram not used                                       ty,             mtflx)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Compute the transports across each edge by integrating the mass
! Subprogram not used ! and tracers over each departure triangle.
! Subprogram not used ! Input variables have the same meanings as in the main subroutine.
! Subprogram not used ! Repeated use of certain sums makes the calculation more efficient.
! Subprogram not used ! Integral formulas are described in triangle_coordinates subroutine.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author William H. Lipscomb, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) ::   &
! Subprogram not used            nx_block, ny_block  ,&! block dimensions
! Subprogram not used            ntrace              ,&! number of tracers in use
! Subprogram not used            integral_order   ! polynomial order for quadrature integrals 
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension (ntrace), intent(in) ::     &
! Subprogram not used            tracer_type       ,&! = 1, 2, or 3 (see comments above)
! Subprogram not used            depend              ! tracer dependencies (see above)
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension (ngroups), intent(in) ::     &
! Subprogram not used            icells           ! number of cells where triarea > puny
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension (nx_block*ny_block,ngroups),     &
! Subprogram not used            intent(in) ::     &
! Subprogram not used            indxi ,&! compressed index in i-direction
! Subprogram not used            indxj   ! compressed index in j-direction
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), intent(in),   &
! Subprogram not used            dimension (nx_block, ny_block, 0:nvert, ngroups) ::   &
! Subprogram not used            xp, yp           ! coordinates of triangle points
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), intent(in),   &
! Subprogram not used            dimension (nx_block, ny_block, ngroups) ::   &
! Subprogram not used            triarea          ! triangle area
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), intent(in),   &
! Subprogram not used            dimension (nx_block, ny_block, ngroups) ::   &
! Subprogram not used            iflux     ,&
! Subprogram not used            jflux
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), intent(in),   &
! Subprogram not used            dimension (nx_block, ny_block) ::   &
! Subprogram not used            mc, mx, my
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), intent(out),   &
! Subprogram not used            dimension (nx_block, ny_block) ::   &
! Subprogram not used            mflx
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), intent(in),   &
! Subprogram not used            dimension (nx_block, ny_block, max_ntrace), optional ::   &
! Subprogram not used            tc, tx, ty
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), intent(out),   &
! Subprogram not used            dimension (nx_block, ny_block, max_ntrace), optional ::   &
! Subprogram not used            mtflx
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) ::   &
! Subprogram not used            i, j, ij      ,&! horizontal indices of edge
! Subprogram not used            i2, j2        ,&! horizontal indices of cell contributing transport
! Subprogram not used            ng            ,&! triangle index
! Subprogram not used            nt, nt1       ,&! tracer indices
! Subprogram not used            ilo,ihi,jlo,jhi ! beginning and end of physical domain
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind) ::   &
! Subprogram not used            m0, m1, m2, m3         ,&! mass field at internal points
! Subprogram not used            w0, w1, w2, w3           ! work variables
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block, ny_block) ::   &
! Subprogram not used            msum, mxsum, mysum     ,&! sum of mass, mass*x, and mass*y
! Subprogram not used            mxxsum, mxysum, myysum   ! sum of mass*x*x, mass*x*y, mass*y*y
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block, ny_block, max_ntrace) ::   &
! Subprogram not used            mtsum            ,&! sum of mass*tracer
! Subprogram not used            mtxsum           ,&! sum of mass*tracer*x
! Subprogram not used            mtysum             ! sum of mass*tracer*y
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Initialize
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       mflx(:,:) = c0
! Subprogram not used       if (present(mtflx)) then
! Subprogram not used          do nt = 1, ntrace
! Subprogram not used             mtflx(:,:,nt) = c0
! Subprogram not used          enddo
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Main loop
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       do ng = 1, ngroups
! Subprogram not used 
! Subprogram not used          if (integral_order == 1) then  ! linear (1-point formula)
! Subprogram not used 
! Subprogram not used !DIR$ CONCURRENT !Cray
! Subprogram not used !cdir nodep      !NEC
! Subprogram not used !ocl novrec      !Fujitsu
! Subprogram not used             do ij = 1, icells(ng)
! Subprogram not used                i = indxi(ij,ng)
! Subprogram not used                j = indxj(ij,ng)
! Subprogram not used 
! Subprogram not used                i2 = iflux(i,j,ng)
! Subprogram not used                j2 = jflux(i,j,ng)
! Subprogram not used 
! Subprogram not used                ! mass transports
! Subprogram not used 
! Subprogram not used                m0 = mc(i2,j2) + xp(i,j,0,ng)*mx(i2,j2)   &
! Subprogram not used                               + yp(i,j,0,ng)*my(i2,j2)
! Subprogram not used                msum(i,j) = m0
! Subprogram not used 
! Subprogram not used                mflx(i,j) = mflx(i,j) + triarea(i,j,ng)*msum(i,j)
! Subprogram not used 
! Subprogram not used                ! quantities needed for tracer transports
! Subprogram not used                mxsum(i,j)  =         m0*xp(i,j,0,ng) 
! Subprogram not used                mxxsum(i,j) = mxsum(i,j)*xp(i,j,0,ng) 
! Subprogram not used                mxysum(i,j) = mxsum(i,j)*yp(i,j,0,ng) 
! Subprogram not used                mysum(i,j)  =         m0*yp(i,j,0,ng) 
! Subprogram not used                myysum(i,j) = mysum(i,j)*yp(i,j,0,ng) 
! Subprogram not used             enddo               ! ij
! Subprogram not used 
! Subprogram not used          elseif (integral_order == 2) then  ! quadratic (3-point formula)
! Subprogram not used 
! Subprogram not used !DIR$ CONCURRENT !Cray
! Subprogram not used !cdir nodep      !NEC
! Subprogram not used !ocl novrec      !Fujitsu
! Subprogram not used             do ij = 1, icells(ng)
! Subprogram not used                i = indxi(ij,ng)
! Subprogram not used                j = indxj(ij,ng)
! Subprogram not used 
! Subprogram not used                i2 = iflux(i,j,ng)
! Subprogram not used                j2 = jflux(i,j,ng)
! Subprogram not used 
! Subprogram not used                ! mass transports
! Subprogram not used                ! Weighting factor of 1/3 is incorporated into the ice
! Subprogram not used                ! area terms m1, m2, and m3.
! Subprogram not used                m1 = p333 * (mc(i2,j2) + xp(i,j,1,ng)*mx(i2,j2)   &
! Subprogram not used                                       + yp(i,j,1,ng)*my(i2,j2))
! Subprogram not used                m2 = p333 * (mc(i2,j2) + xp(i,j,2,ng)*mx(i2,j2)   &
! Subprogram not used                                       + yp(i,j,2,ng)*my(i2,j2))
! Subprogram not used                m3 = p333 * (mc(i2,j2) + xp(i,j,3,ng)*mx(i2,j2)   &
! Subprogram not used                                       + yp(i,j,3,ng)*my(i2,j2))
! Subprogram not used                msum(i,j) = m1 + m2 + m3
! Subprogram not used                mflx(i,j) = mflx(i,j) + triarea(i,j,ng)*msum(i,j)
! Subprogram not used 
! Subprogram not used                ! quantities needed for mass_tracer transports
! Subprogram not used                w1 = m1 * xp(i,j,1,ng)
! Subprogram not used                w2 = m2 * xp(i,j,2,ng)
! Subprogram not used                w3 = m3 * xp(i,j,3,ng)
! Subprogram not used 
! Subprogram not used                mxsum(i,j) = w1 + w2 + w3
! Subprogram not used 
! Subprogram not used                mxxsum(i,j) = w1*xp(i,j,1,ng) + w2*xp(i,j,2,ng)   &
! Subprogram not used                            + w3*xp(i,j,3,ng) 
! Subprogram not used 
! Subprogram not used                mxysum(i,j) = w1*yp(i,j,1,ng) + w2*yp(i,j,2,ng)   &
! Subprogram not used                            + w3*yp(i,j,3,ng)
! Subprogram not used 
! Subprogram not used                w1 = m1 * yp(i,j,1,ng)
! Subprogram not used                w2 = m2 * yp(i,j,2,ng)
! Subprogram not used                w3 = m3 * yp(i,j,3,ng)
! Subprogram not used 
! Subprogram not used                mysum(i,j) = w1 + w2 + w3
! Subprogram not used 
! Subprogram not used                myysum(i,j) = w1*yp(i,j,1,ng) + w2*yp(i,j,2,ng)   &
! Subprogram not used                            + w3*yp(i,j,3,ng)
! Subprogram not used             enddo               ! ij
! Subprogram not used 
! Subprogram not used          else                   ! cubic (4-point formula)
! Subprogram not used 
! Subprogram not used !DIR$ CONCURRENT !Cray
! Subprogram not used !cdir nodep      !NEC
! Subprogram not used !ocl novrec      !Fujitsu
! Subprogram not used             do ij = 1, icells(ng)
! Subprogram not used                i = indxi(ij,ng)
! Subprogram not used                j = indxj(ij,ng)
! Subprogram not used 
! Subprogram not used                i2 = iflux(i,j,ng)
! Subprogram not used                j2 = jflux(i,j,ng)
! Subprogram not used 
! Subprogram not used                ! mass transports
! Subprogram not used 
! Subprogram not used                ! Weighting factors are incorporated into the
! Subprogram not used                ! terms m0, m1, m2, and m3.
! Subprogram not used                m0 = p5625m * (mc(i2,j2) + xp(i,j,0,ng)*mx(i2,j2)   &
! Subprogram not used                                         + yp(i,j,0,ng)*my(i2,j2))
! Subprogram not used                m1 = p52083 * (mc(i2,j2) + xp(i,j,1,ng)*mx(i2,j2)   &
! Subprogram not used                                         + yp(i,j,1,ng)*my(i2,j2))
! Subprogram not used                m2 = p52083 * (mc(i2,j2) + xp(i,j,2,ng)*mx(i2,j2)   &
! Subprogram not used                                         + yp(i,j,2,ng)*my(i2,j2))
! Subprogram not used                m3 = p52083 * (mc(i2,j2) + xp(i,j,3,ng)*mx(i2,j2)   &
! Subprogram not used                                         + yp(i,j,3,ng)*my(i2,j2))
! Subprogram not used                msum(i,j) = m0 + m1 + m2 + m3
! Subprogram not used                mflx(i,j) = mflx(i,j) + triarea(i,j,ng)*msum(i,j)
! Subprogram not used 
! Subprogram not used                ! quantities needed for tracer transports
! Subprogram not used                w0 = m0 * xp(i,j,0,ng)
! Subprogram not used                w1 = m1 * xp(i,j,1,ng)
! Subprogram not used                w2 = m2 * xp(i,j,2,ng)
! Subprogram not used                w3 = m3 * xp(i,j,3,ng)
! Subprogram not used 
! Subprogram not used                mxsum(i,j) = w0 + w1 + w2 + w3
! Subprogram not used 
! Subprogram not used                mxxsum(i,j) = w0*xp(i,j,0,ng) + w1*xp(i,j,1,ng)   &
! Subprogram not used                            + w2*xp(i,j,2,ng) + w3*xp(i,j,3,ng)
! Subprogram not used 
! Subprogram not used                mxysum(i,j) = w0*yp(i,j,0,ng) + w1*yp(i,j,1,ng)   &
! Subprogram not used                            + w2*yp(i,j,2,ng) + w3*yp(i,j,3,ng)
! Subprogram not used 
! Subprogram not used                w0 = m0 * yp(i,j,0,ng)
! Subprogram not used                w1 = m1 * yp(i,j,1,ng)
! Subprogram not used                w2 = m2 * yp(i,j,2,ng)
! Subprogram not used                w3 = m3 * yp(i,j,3,ng)
! Subprogram not used 
! Subprogram not used                mysum(i,j) = w0 + w1 + w2 + w3
! Subprogram not used 
! Subprogram not used                myysum(i,j) = w0*yp(i,j,0,ng) + w1*yp(i,j,1,ng)   &
! Subprogram not used                            + w2*yp(i,j,2,ng) + w3*yp(i,j,3,ng)
! Subprogram not used 
! Subprogram not used             enddo               ! ij
! Subprogram not used 
! Subprogram not used          endif                  ! integral_order
! Subprogram not used 
! Subprogram not used          ! mass * tracer transports
! Subprogram not used 
! Subprogram not used          if (present(mtflx)) then
! Subprogram not used 
! Subprogram not used             do nt = 1, ntrace
! Subprogram not used                if (tracer_type(nt)==1) then ! does not depend on another tracer
! Subprogram not used 
! Subprogram not used !DIR$ CONCURRENT !Cray
! Subprogram not used !cdir nodep      !NEC
! Subprogram not used !ocl novrec      !Fujitsu
! Subprogram not used                   do ij = 1, icells(ng)
! Subprogram not used                      i = indxi(ij,ng)
! Subprogram not used                      j = indxj(ij,ng)
! Subprogram not used 
! Subprogram not used                      i2 = iflux(i,j,ng)
! Subprogram not used                      j2 = jflux(i,j,ng)
! Subprogram not used 
! Subprogram not used                      mtsum(i,j,nt) =  msum(i,j) * tc(i2,j2,nt)   &
! Subprogram not used                                    + mxsum(i,j) * tx(i2,j2,nt)   &
! Subprogram not used                                    + mysum(i,j) * ty(i2,j2,nt)
! Subprogram not used 
! Subprogram not used                      mtflx(i,j,nt) = mtflx(i,j,nt)   &
! Subprogram not used                                  + triarea(i,j,ng) * mtsum(i,j,nt)
! Subprogram not used 
! Subprogram not used                      ! quantities needed for dependent tracers
! Subprogram not used 
! Subprogram not used                      mtxsum(i,j,nt) =  mxsum(i,j) * tc(i2,j2,nt)   &
! Subprogram not used                                     + mxxsum(i,j) * tx(i2,j2,nt)   &
! Subprogram not used                                     + mxysum(i,j) * ty(i2,j2,nt)
! Subprogram not used 
! Subprogram not used                      mtysum(i,j,nt) =  mysum(i,j) * tc(i2,j2,nt)   &
! Subprogram not used                                     + mxysum(i,j) * tx(i2,j2,nt)   &
! Subprogram not used                                     + myysum(i,j) * ty(i2,j2,nt)
! Subprogram not used                   enddo         ! ij
! Subprogram not used 
! Subprogram not used                elseif (tracer_type(nt)==2) then ! depends on another tracer
! Subprogram not used                   nt1 = depend(nt)
! Subprogram not used 
! Subprogram not used !DIR$ CONCURRENT !Cray
! Subprogram not used !cdir nodep      !NEC
! Subprogram not used !ocl novrec      !Fujitsu
! Subprogram not used                   do ij = 1, icells(ng)
! Subprogram not used                      i = indxi(ij,ng)
! Subprogram not used                      j = indxj(ij,ng)
! Subprogram not used 
! Subprogram not used                      i2 = iflux(i,j,ng)
! Subprogram not used                      j2 = jflux(i,j,ng)
! Subprogram not used 
! Subprogram not used                      mtsum(i,j,nt) =  mtsum(i,j,nt1) * tc(i2,j2,nt)   &
! Subprogram not used                                    + mtxsum(i,j,nt1) * tx(i2,j2,nt)   &
! Subprogram not used                                    + mtysum(i,j,nt1) * ty(i2,j2,nt)
! Subprogram not used 
! Subprogram not used                      mtflx(i,j,nt) = mtflx(i,j,nt)   &
! Subprogram not used                                    + triarea(i,j,ng) * mtsum(i,j,nt)
! Subprogram not used                   enddo         ! ij
! Subprogram not used 
! Subprogram not used 
! Subprogram not used                elseif (tracer_type(nt)==3) then ! depends on two tracers
! Subprogram not used                   nt1 = depend(nt)
! Subprogram not used 
! Subprogram not used !DIR$ CONCURRENT !Cray
! Subprogram not used !cdir nodep      !NEC
! Subprogram not used !ocl novrec      !Fujitsu
! Subprogram not used                   do ij = 1, icells(ng)
! Subprogram not used                      i = indxi(ij,ng)
! Subprogram not used                      j = indxj(ij,ng)
! Subprogram not used 
! Subprogram not used                      i2 = iflux(i,j,ng)
! Subprogram not used                      j2 = jflux(i,j,ng)
! Subprogram not used 
! Subprogram not used                      ! upwind approx (tx=ty=0) for type 3 tracers
! Subprogram not used                      mtsum(i,j,nt) =  mtsum(i,j,nt1) * tc(i2,j2,nt)
! Subprogram not used 
! Subprogram not used                      mtflx(i,j,nt) = mtflx(i,j,nt)   &
! Subprogram not used                                    + triarea(i,j,ng) * mtsum(i,j,nt)
! Subprogram not used                   enddo         ! ij
! Subprogram not used 
! Subprogram not used                endif            ! tracer type
! Subprogram not used             enddo               ! ntrace
! Subprogram not used          endif                  ! present(mtflx)
! Subprogram not used       enddo                     ! ng
! Subprogram not used 
! Subprogram not used       end subroutine transport_integrals

!=======================================================================
!
!BOP
!
! !IROUTINE: update_fields - compute new area and tracers
!
! !INTERFACE:
!
! Subprogram not used       subroutine update_fields (nx_block,    ny_block,   &
! Subprogram not used                                 ilo, ihi,    jlo, jhi,   &
! Subprogram not used                                 ntrace,                  &
! Subprogram not used                                 tracer_type, depend,     &
! Subprogram not used                                 tarear,      l_stop,     &
! Subprogram not used                                 istop,       jstop,      &
! Subprogram not used                                 mflxe,       mflxn,      &
! Subprogram not used                                 mm,                      &
! Subprogram not used                                 mtflxe,      mtflxn,     &
! Subprogram not used                                 tm)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Given transports through cell edges, compute new area and tracers.
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author William H. Lipscomb, LANL
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) ::   &
! Subprogram not used          nx_block, ny_block,&! block dimensions
! Subprogram not used          ilo,ihi,jlo,jhi   ,&! beginning and end of physical domain
! Subprogram not used          ntrace              ! number of tracers in use
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension (ntrace), intent(in) ::     &
! Subprogram not used          tracer_type       ,&! = 1, 2, or 3 (see comments above)
! Subprogram not used          depend              ! tracer dependencies (see above)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block, ny_block),   &
! Subprogram not used          intent(in) ::   &
! Subprogram not used          mflxe, mflxn   ,&! mass transport across east and north cell edges
! Subprogram not used          tarear           ! 1/tarea
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block, ny_block),   &
! Subprogram not used          intent(inout) ::   &
! Subprogram not used          mm               ! mass field (mean)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block, ny_block, max_ntrace),   &
! Subprogram not used          intent(in), optional ::   &
! Subprogram not used          mtflxe, mtflxn   ! mass*tracer transport across E and N cell edges
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension (nx_block, ny_block, max_ntrace),   &
! Subprogram not used          intent(inout), optional ::   &
! Subprogram not used          tm               ! tracer fields
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind), intent(inout) ::   &
! Subprogram not used          l_stop           ! if true, abort on return
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), intent(inout) ::   &
! Subprogram not used          istop, jstop     ! indices of grid cell where model aborts 
! Subprogram not used 
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) ::   &
! Subprogram not used          i, j           ,&! horizontal indices
! Subprogram not used          nt, nt1, nt2     ! tracer indices
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension(nx_block,ny_block,max_ntrace) ::   &
! Subprogram not used          mtold            ! old mass*tracer
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind) ::   &
! Subprogram not used          w1, w2           ! work variables
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), dimension(nx_block*ny_block) ::   &
! Subprogram not used          indxi          ,&! compressed indices in i and j directions
! Subprogram not used          indxj
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind) ::   &
! Subprogram not used          icells         ,&! number of cells with mm > 0.
! Subprogram not used          ij               ! combined i/j horizontal index
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Save starting values of mass*tracer
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       if (present(tm)) then
! Subprogram not used          do nt = 1, ntrace
! Subprogram not used             if (tracer_type(nt)==1) then ! does not depend on other tracers
! Subprogram not used                do j = jlo, jhi
! Subprogram not used                do i = ilo, ihi
! Subprogram not used                   mtold(i,j,nt) = mm(i,j) * tm(i,j,nt)
! Subprogram not used                enddo            ! i
! Subprogram not used                enddo              ! j
! Subprogram not used             elseif (tracer_type(nt)==2) then  ! depends on another tracer
! Subprogram not used                nt1 = depend(nt)
! Subprogram not used                do j = jlo, jhi
! Subprogram not used                do i = ilo, ihi
! Subprogram not used                   mtold(i,j,nt) = mm(i,j) * tm(i,j,nt1) * tm(i,j,nt)
! Subprogram not used                enddo            ! i
! Subprogram not used                enddo            ! j
! Subprogram not used             elseif (tracer_type(nt)==3) then  ! depends on two tracers
! Subprogram not used                nt1 = depend(nt)
! Subprogram not used                nt2 = depend(nt1)
! Subprogram not used                do j = jlo, jhi
! Subprogram not used                do i = ilo, ihi
! Subprogram not used                   mtold(i,j,nt) = mm(i,j)    &
! Subprogram not used                             * tm(i,j,nt2) * tm(i,j,nt1) * tm(i,j,nt)
! Subprogram not used                enddo            ! i
! Subprogram not used                enddo            ! j
! Subprogram not used 
! Subprogram not used  
! Subprogram not used             endif               ! depend(nt) = 0
! Subprogram not used          enddo                  ! nt
! Subprogram not used       endif                     ! present(tm)
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Update mass field
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       do j = jlo, jhi
! Subprogram not used       do i = ilo, ihi
! Subprogram not used 
! Subprogram not used          w1 = mflxe(i,j) - mflxe(i-1,j)   &
! Subprogram not used             + mflxn(i,j) - mflxn(i,j-1)
! Subprogram not used          mm(i,j) = mm(i,j) - w1*tarear(i,j)
! Subprogram not used 
! Subprogram not used          if (mm(i,j) < -puny) then    ! abort with negative value
! Subprogram not used             l_stop = .true.
! Subprogram not used             istop = i
! Subprogram not used             jstop = j
! Subprogram not used          elseif (mm(i,j) < c0) then   ! set to zero
! Subprogram not used             mm(i,j) = c0
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used       enddo
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       if (l_stop) then
! Subprogram not used          i = istop
! Subprogram not used          j = jstop
! Subprogram not used          w1 = mflxe(i,j) - mflxe(i-1,j)   &
! Subprogram not used             + mflxn(i,j) - mflxn(i,j-1)
! Subprogram not used          write (nu_diag,*) ' '
! Subprogram not used          write (nu_diag,*) 'New mass < 0, i, j =', i, j
! Subprogram not used          write (nu_diag,*) 'Old mass =', mm(i,j) + w1*tarear(i,j)
! Subprogram not used          write (nu_diag,*) 'New mass =', mm(i,j)
! Subprogram not used          write (nu_diag,*) 'Net transport =', -w1*tarear(i,j)
! Subprogram not used          return
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Update tracers
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       if (present(tm)) then
! Subprogram not used 
! Subprogram not used          icells = 0
! Subprogram not used          do j = jlo, jhi
! Subprogram not used          do i = ilo, ihi
! Subprogram not used             if (mm(i,j) > c0) then ! grid cells with positive areas
! Subprogram not used                icells = icells + 1
! Subprogram not used                indxi(icells) = i
! Subprogram not used                indxj(icells) = j
! Subprogram not used             endif
! Subprogram not used          enddo                  ! i
! Subprogram not used          enddo                  ! j
! Subprogram not used 
! Subprogram not used          do nt = 1, ntrace
! Subprogram not used 
! Subprogram not used             do j = jlo, jhi
! Subprogram not used             do i = ilo, ihi
! Subprogram not used                tm(i,j,nt) = c0
! Subprogram not used             enddo
! Subprogram not used             enddo
! Subprogram not used 
! Subprogram not used             if (tracer_type(nt)==1) then ! does not depend on other tracers
! Subprogram not used 
! Subprogram not used                do ij = 1, icells
! Subprogram not used                   i = indxi(ij)
! Subprogram not used                   j = indxj(ij)
! Subprogram not used 
! Subprogram not used                   w1  = mtflxe(i,j,nt) - mtflxe(i-1,j,nt)   &
! Subprogram not used                       + mtflxn(i,j,nt) - mtflxn(i,j-1,nt)
! Subprogram not used                   tm(i,j,nt) = (mtold(i,j,nt) - w1*tarear(i,j))   &
! Subprogram not used                                 / mm(i,j)
! Subprogram not used                enddo            ! ij
! Subprogram not used 
! Subprogram not used 
! Subprogram not used             elseif (tracer_type(nt)==2) then ! depends on another tracer
! Subprogram not used                nt1 = depend(nt)
! Subprogram not used 
! Subprogram not used !DIR$ CONCURRENT !Cray
! Subprogram not used !cdir nodep      !NEC
! Subprogram not used !ocl novrec      !Fujitsu
! Subprogram not used                do ij = 1, icells
! Subprogram not used                   i = indxi(ij)
! Subprogram not used                   j = indxj(ij)
! Subprogram not used 
! Subprogram not used                   if (abs(tm(i,j,nt1)) > c0) then
! Subprogram not used                      w1  = mtflxe(i,j,nt) - mtflxe(i-1,j,nt)   &
! Subprogram not used                          + mtflxn(i,j,nt) - mtflxn(i,j-1,nt)
! Subprogram not used                      tm(i,j,nt) = (mtold(i,j,nt) - w1*tarear(i,j))   &
! Subprogram not used                                  / (mm(i,j) * tm(i,j,nt1))
! Subprogram not used                   endif
! Subprogram not used                enddo            ! ij
! Subprogram not used 
! Subprogram not used             elseif (tracer_type(nt)==3) then ! depends on two tracers
! Subprogram not used                nt1 = depend(nt)
! Subprogram not used                nt2 = depend(nt1)
! Subprogram not used 
! Subprogram not used !DIR$ CONCURRENT !Cray
! Subprogram not used !cdir nodep      !NEC
! Subprogram not used !ocl novrec      !Fujitsu
! Subprogram not used                do ij = 1, icells
! Subprogram not used                   i = indxi(ij)
! Subprogram not used                   j = indxj(ij)
! Subprogram not used 
! Subprogram not used                   if (abs(tm(i,j,nt1)) > c0 .and.   &
! Subprogram not used                       abs(tm(i,j,nt2)) > c0) then
! Subprogram not used                      w1  = mtflxe(i,j,nt) - mtflxe(i-1,j,nt)   &
! Subprogram not used                          + mtflxn(i,j,nt) - mtflxn(i,j-1,nt)
! Subprogram not used                      tm(i,j,nt) = (mtold(i,j,nt) - w1*tarear(i,j))   &
! Subprogram not used                               / (mm(i,j) * tm(i,j,nt2) * tm(i,j,nt1))
! Subprogram not used                   endif
! Subprogram not used                enddo            ! ij
! Subprogram not used             endif               ! tracer_type
! Subprogram not used          enddo                  ! nt
! Subprogram not used       endif                     ! present(tm)
! Subprogram not used 
! Subprogram not used       end subroutine update_fields

!=======================================================================

      end module ice_transport_remap

!=======================================================================
