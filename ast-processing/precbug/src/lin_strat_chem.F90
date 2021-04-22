!--------------------------------------------------------------------
! linearized ozone chemistry LINOZ
! from Hsu and Prather, JGR, 2008
!
! written by Jean-Francois Lamarque (September 2008)
! modified by
!     24 Oct 2008 -- Francis Vitt
!      9 Dec 2008 -- Philip Cameron-Smith, LLNL, -- added ltrop
!--------------------------------------------------------------------
module lin_strat_chem

  use shr_kind_mod, only : r8 => shr_kind_r8
  use ppgrid,       only : begchunk, endchunk
  use physics_types,only : physics_state
  use cam_logfile,  only : iulog
  use abortutils,   only : endrun
  use spmd_utils,   only : masterproc
  !
  implicit none
  !
  private  ! all unless made public

  save
  !
  ! define public components of module
  !
  public :: lin_strat_chem_inti, lin_strat_chem_solve
  public :: do_lin_strat_chem

  integer :: index_o3
  logical :: do_lin_strat_chem


contains

!--------------------------------------------------------------------
!--------------------------------------------------------------------
  subroutine lin_strat_chem_inti(phys_state)
    !
    ! initialize linearized stratospheric chemistry by reading
    ! input parameters from netcdf file and interpolate to
    ! present model grid
    !
    use linoz_data,   only : linoz_data_init, has_linoz_data
    use ppgrid,       only : pver
    use mo_chem_utls, only : get_spc_ndx
    use cam_history,  only : addfld, phys_decomp, add_default
    use physics_buffer, only : physics_buffer_desc
    use phys_control, only : phys_getopts

    implicit none


    type(physics_state), intent(in) :: phys_state(begchunk:endchunk)

    logical :: history_chemistry

    call phys_getopts(history_chemistry_out=history_chemistry)

    if (.not.has_linoz_data) return

    !
    ! find index of ozone
    !
    index_o3 = get_spc_ndx('O3')
    do_lin_strat_chem = has_linoz_data
    if ( index_o3 <= 0 ) then
       write(iulog,*) ' No ozone in the chemical mechanism, skipping lin_strat_chem'
       do_lin_strat_chem = .false.
       return
    end if

    ! check for synoz

    if( get_spc_ndx( 'SYNOZ' ) > 0 .and. has_linoz_data) then
       call endrun('lin_strat_chem_inti: cannot have both synoz and linoz')
    endif

    ! initialize the linoz data

    call linoz_data_init()

    ! define additional output

    call addfld( 'LINOZ_DO3'    , '/s'     , pver, 'A', 'ozone vmr tendency by linearized ozone chemistry'  , phys_decomp )
    call addfld( 'LINOZ_DO3_PSC', '/s'     , pver, 'A', 'ozone vmr loss by PSCs using Carille et al. (1990)', phys_decomp )
    call addfld( 'LINOZ_SSO3'   , 'kg'     , pver, 'A', 'steady state ozone in LINOZ'                       , phys_decomp )
    call addfld( 'LINOZ_O3COL'  , 'DU'     , pver, 'A', 'ozone column above'                                , phys_decomp )
    call addfld( 'LINOZ_O3CLIM' , 'mol/mol', pver, 'A', 'climatology of ozone in LINOZ'                     , phys_decomp )
    call addfld( 'LINOZ_SZA'    , 'degrees',    1, 'A', 'solar zenith angle in LINOZ'                       , phys_decomp )

    if (history_chemistry) then
       call add_default( 'LINOZ_DO3'    , 1, ' ' )
       call add_default( 'LINOZ_DO3_PSC', 1, ' ' )
       call add_default( 'LINOZ_SSO3'   , 1, ' ' )
       call add_default( 'LINOZ_O3COL'  , 1, ' ' )
       call add_default( 'LINOZ_O3CLIM' , 1, ' ' )
       call add_default( 'LINOZ_SZA'    , 1, ' ' )
    end if

    return
  end subroutine lin_strat_chem_inti


!--------------------------------------------------------------------
!--------------------------------------------------------------------
! Subprogram not used   subroutine lin_strat_chem_solve( ncol, lchnk, o3_vmr, o3col, temp, sza, pmid, delta_t, rlats, ltrop )
! Subprogram not used  
! Subprogram not used     use chlorine_loading_data, only: chlorine_loading
! Subprogram not used 
! Subprogram not used     !
! Subprogram not used     ! this subroutine updates the ozone mixing ratio in the stratosphere
! Subprogram not used     ! using linearized chemistry 
! Subprogram not used     !
! Subprogram not used 
! Subprogram not used     use ppgrid,        only : pcols, pver
! Subprogram not used     use physconst,     only : pi, &
! Subprogram not used                               grav => gravit, &
! Subprogram not used                               mw_air => mwdry
! Subprogram not used     use cam_history,   only : outfld
! Subprogram not used     use linoz_data,    only : fields, o3_clim_ndx,t_clim_ndx,o3col_clim_ndx,PmL_clim_ndx,dPmL_dO3_ndx,&
! Subprogram not used                                       dPmL_dT_ndx,dPmL_dO3col_ndx,cariolle_pscs_ndx
! Subprogram not used     !
! Subprogram not used     ! dummy arguments
! Subprogram not used     !
! Subprogram not used     integer,  intent(in)                           :: ncol                ! number of columns in chunk
! Subprogram not used     integer,  intent(in)                           :: lchnk               ! chunk index
! Subprogram not used     real(r8), intent(inout), dimension(ncol ,pver) :: o3_vmr              ! ozone volume mixing ratio
! Subprogram not used     real(r8), intent(in)   , dimension(ncol ,pver) :: o3col               ! ozone column above box (mol/cm^2)
! Subprogram not used     real(r8), intent(in)   , dimension(pcols,pver) :: temp                ! temperature (K)
! Subprogram not used     real(r8), intent(in)   , dimension(ncol )      :: sza                 ! local solar zenith angle
! Subprogram not used     real(r8), intent(in)   , dimension(pcols,pver) :: pmid                ! midpoint pressure (Pa)
! Subprogram not used     real(r8), intent(in)                           :: delta_t             ! timestep size (secs)
! Subprogram not used     real(r8), intent(in)                           :: rlats(ncol)         ! column latitudes (radians)
! Subprogram not used     integer,  intent(in)   , dimension(pcols)      :: ltrop               ! chunk index
! Subprogram not used     !
! Subprogram not used     ! local
! Subprogram not used     !
! Subprogram not used     integer :: i,k,n !,index_lat,index_month
! Subprogram not used     real(r8) :: o3col_du,delta_temp,delta_o3col,o3_old,o3_new,delta_o3
! Subprogram not used     real(r8) :: max_sza, psc_loss
! Subprogram not used     real(r8) :: o3_clim
! Subprogram not used     real(r8), dimension(ncol) :: lats
! Subprogram not used     real(r8), dimension(ncol,pver) :: do3_linoz,do3_linoz_psc,ss_o3,o3col_du_diag,o3clim_linoz_diag
! Subprogram not used 
! Subprogram not used     real(r8), dimension(:,:), pointer :: linoz_o3_clim
! Subprogram not used     real(r8), dimension(:,:), pointer :: linoz_t_clim
! Subprogram not used     real(r8), dimension(:,:), pointer :: linoz_o3col_clim
! Subprogram not used     real(r8), dimension(:,:), pointer :: linoz_PmL_clim
! Subprogram not used     real(r8), dimension(:,:), pointer :: linoz_dPmL_dO3
! Subprogram not used     real(r8), dimension(:,:), pointer :: linoz_dPmL_dT
! Subprogram not used     real(r8), dimension(:,:), pointer :: linoz_dPmL_dO3col
! Subprogram not used     real(r8), dimension(:,:), pointer :: linoz_cariolle_psc
! Subprogram not used 
! Subprogram not used     !
! Subprogram not used     ! parameters
! Subprogram not used     !
! Subprogram not used     real(r8), parameter :: convert_to_du = 1._r8/(2.687e16_r8)      ! convert ozone column from mol/cm^2 to DU
! Subprogram not used     real(r8), parameter :: degrees_to_radians = pi/180._r8          ! conversion factors
! Subprogram not used     real(r8), parameter :: radians_to_degrees = 180._r8/pi
! Subprogram not used     real(r8), parameter :: temp_activation_cariolle = 193._r8       ! O3 loss freq when T below temp_activation_cariolle (K)
! Subprogram not used     real(r8), parameter :: chlorine_loading_1987    = 2.5977_r8     ! EESC value (ppbv)
! Subprogram not used     real(r8), parameter :: chlorine_loading_bgnd    = 0.0000_r8     ! EESC value (ppbv) for background conditions
! Subprogram not used     real(r8), parameter :: pressure_threshold       = 210.e+2_r8    ! {PJC} for diagnostics only
! Subprogram not used 
! Subprogram not used     !
! Subprogram not used     ! skip if no ozone field available
! Subprogram not used     !
! Subprogram not used     if ( .not. do_lin_strat_chem ) return
! Subprogram not used 
! Subprogram not used     ! 
! Subprogram not used     ! associate the field pointers
! Subprogram not used     !
! Subprogram not used     linoz_o3_clim      => fields(o3_clim_ndx)      %data(:,:,lchnk )
! Subprogram not used     linoz_t_clim       => fields(t_clim_ndx)       %data(:,:,lchnk )
! Subprogram not used     linoz_o3col_clim   => fields(o3col_clim_ndx)   %data(:,:,lchnk )
! Subprogram not used     linoz_PmL_clim     => fields(PmL_clim_ndx)     %data(:,:,lchnk )
! Subprogram not used     linoz_dPmL_dO3     => fields(dPmL_dO3_ndx)     %data(:,:,lchnk )
! Subprogram not used     linoz_dPmL_dT      => fields(dPmL_dT_ndx)      %data(:,:,lchnk )
! Subprogram not used     linoz_dPmL_dO3col  => fields(dPmL_dO3col_ndx)  %data(:,:,lchnk )
! Subprogram not used     linoz_cariolle_psc => fields(cariolle_pscs_ndx)%data(:,:,lchnk )
! Subprogram not used 
! Subprogram not used     !
! Subprogram not used     ! initialize output arrays
! Subprogram not used     !
! Subprogram not used     do3_linoz         = 0._r8
! Subprogram not used     do3_linoz_psc     = 0._r8
! Subprogram not used     o3col_du_diag     = 0._r8
! Subprogram not used     o3clim_linoz_diag = 0._r8
! Subprogram not used     ss_o3             = 0._r8
! Subprogram not used     !
! Subprogram not used     ! convert lats from radians to degrees
! Subprogram not used     !
! Subprogram not used     lats = rlats * radians_to_degrees
! Subprogram not used 
! Subprogram not used     LOOP_COL: do i=1,ncol
! Subprogram not used        LOOP_LEV: do k=1,ltrop(i)
! Subprogram not used           !
! Subprogram not used           ! climatological ozone
! Subprogram not used           !
! Subprogram not used           o3_clim = linoz_o3_clim(i,k)
! Subprogram not used           !
! Subprogram not used           ! skip if not in the stratosphere
! Subprogram not used           !
! Subprogram not used           if ( pmid(i,k) > pressure_threshold ) THEN   ! PJC diagnostic
! Subprogram not used              WRITE(iulog,*)'LINOZ WARNING: Exceeded PRESSURE threshold (i,k,p_threshold,pmid,o3)=',&
! Subprogram not used                i,k,nint(pressure_threshold/100._r8),'mb',nint(pmid(i,k)/100._r8),'mb',nint(o3_vmr(i,k)*1e9_r8),'ppb'   !PJC
! Subprogram not used !             cycle LOOP_LEV
! Subprogram not used           endif
! Subprogram not used           !
! Subprogram not used           ! diagnostic for output
! Subprogram not used           !
! Subprogram not used           o3clim_linoz_diag(i,k) = o3_clim
! Subprogram not used           !
! Subprogram not used           ! old ozone mixing ratio
! Subprogram not used           !
! Subprogram not used           o3_old = o3_vmr(i,k)
! Subprogram not used           !
! Subprogram not used           ! convert o3col from mol/cm2
! Subprogram not used           !
! Subprogram not used           o3col_du = o3col(i,k) * convert_to_du
! Subprogram not used           o3col_du_diag(i,k) = o3col_du
! Subprogram not used           !
! Subprogram not used           ! compute differences from climatology
! Subprogram not used           !
! Subprogram not used           delta_temp  = temp(i,k) - linoz_t_clim    (i,k)
! Subprogram not used           delta_o3col = o3col_du  - linoz_o3col_clim(i,k)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used           !
! Subprogram not used           ! steady state ozone
! Subprogram not used           !
! Subprogram not used           ss_o3(i,k) = o3_clim - (               linoz_PmL_clim   (i,k)   &
! Subprogram not used                                  + delta_o3col * linoz_dPmL_dO3col(i,k)   &
! Subprogram not used                                  + delta_temp  * linoz_dPmL_dT    (i,k)   &
! Subprogram not used                                              ) / linoz_dPmL_dO3   (i,k)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used           !
! Subprogram not used           ! ozone change
! Subprogram not used           !
! Subprogram not used           delta_o3 = (ss_o3(i,k)-o3_old) * (1._r8 - exp(linoz_dPmL_dO3(i,k)*delta_t))
! Subprogram not used           !
! Subprogram not used           ! define new ozone mixing ratio
! Subprogram not used           !
! Subprogram not used           o3_new = o3_old + delta_o3
! Subprogram not used           !
! Subprogram not used           ! output diagnostic
! Subprogram not used           !
! Subprogram not used           do3_linoz(i,k) = delta_o3/delta_t
! Subprogram not used           !
! Subprogram not used           ! PSC activation (follows Cariolle et al 1990.)
! Subprogram not used           !
! Subprogram not used           ! use only if abs(latitude) > 40.
! Subprogram not used           !
! Subprogram not used           if ( abs(lats(i)) > 40._r8 ) then   
! Subprogram not used              if ( (chlorine_loading-chlorine_loading_bgnd) > 0._r8 ) then
! Subprogram not used                 if ( temp(i,k) <= temp_activation_cariolle ) then
! Subprogram not used                    !
! Subprogram not used                    ! define maximum SZA for PSC loss (= tangent height at sunset)
! Subprogram not used                    !
! Subprogram not used                    max_sza = (90._r8 + sqrt( max( 16._r8*log10(100000._r8/pmid(i,k)),0._r8)))
! Subprogram not used                    if ( (sza(i)*radians_to_degrees) <= max_sza ) then
! Subprogram not used 
! Subprogram not used                       psc_loss = exp(-linoz_cariolle_psc(i,k) &
! Subprogram not used                            * (chlorine_loading/chlorine_loading_1987)**2 &
! Subprogram not used                            * delta_t )
! Subprogram not used 
! Subprogram not used                       o3_new = o3_old * psc_loss
! Subprogram not used                       !
! Subprogram not used                       ! output diagnostic
! Subprogram not used                       !
! Subprogram not used                       do3_linoz_psc(i,k) = (o3_new-o3_old)/delta_t
! Subprogram not used                       !
! Subprogram not used                    end if
! Subprogram not used                 end if
! Subprogram not used              end if
! Subprogram not used           end if
! Subprogram not used           !
! Subprogram not used           ! update ozone vmr
! Subprogram not used           !
! Subprogram not used           o3_vmr(i,k) = o3_new
! Subprogram not used 
! Subprogram not used        end do LOOP_LEV
! Subprogram not used     end do LOOP_COL
! Subprogram not used     !
! Subprogram not used     ! output
! Subprogram not used     !
! Subprogram not used     call outfld( 'LINOZ_DO3'    , do3_linoz              , ncol, lchnk )
! Subprogram not used     call outfld( 'LINOZ_DO3_PSC', do3_linoz_psc          , ncol, lchnk )
! Subprogram not used     call outfld( 'LINOZ_SSO3'   , ss_o3                  , ncol, lchnk )
! Subprogram not used     call outfld( 'LINOZ_O3COL'  , o3col_du_diag          , ncol, lchnk )
! Subprogram not used     call outfld( 'LINOZ_O3CLIM' , o3clim_linoz_diag      , ncol, lchnk )
! Subprogram not used     call outfld( 'LINOZ_SZA'    ,(sza*radians_to_degrees), ncol, lchnk )
! Subprogram not used 
! Subprogram not used     return
! Subprogram not used   end subroutine lin_strat_chem_solve

end module lin_strat_chem
