
      module mo_waccm_hrates

      use shr_kind_mod,      only : r8 => shr_kind_r8
      use cam_logfile,       only : iulog

      implicit none

      save

      real(r8), parameter :: secpday       = 86400._r8
      real(r8), parameter :: daypsec       = 1._r8/secpday
      real(r8), parameter :: aur_therm     = 807._r8
      real(r8), parameter :: jkcal         = 4184._r8
      real(r8), parameter :: aur_heat_eff  = .05_r8
      real(r8), parameter :: aur_hconst    = 1.e3_r8*jkcal*aur_therm*aur_heat_eff

      real(r8) :: max_zen_angle

      private
      public :: waccm_hrates, init_hrates, has_hrates

      integer :: id_co2, id_o2, id_o3, id_o2_1d, id_o2_1s, id_h2o, id_o, id_h
      logical :: has_hrates

      contains
   
      subroutine init_hrates( )
        use mo_chem_utls, only : get_spc_ndx
        use cam_history,  only : addfld, phys_decomp
        use ppgrid,       only : pver
        use ref_pres,     only : ptop_ref, psurf_ref


        implicit none

        integer :: ids(8)
        character(len=128) :: attr  ! netcdf variable attribute

        id_co2   = get_spc_ndx( 'CO2' )
        id_o2    = get_spc_ndx( 'O2' )
        id_o3    = get_spc_ndx( 'O3' )
        id_o2_1d = get_spc_ndx( 'O2_1D' )
        id_o2_1s = get_spc_ndx( 'O2_1S' )
        id_h2o   = get_spc_ndx( 'H2O' )
        id_o     = get_spc_ndx( 'O' )
        id_h     = get_spc_ndx( 'H' )

        ids = (/ id_co2, id_o2, id_o3, id_o2_1d, id_o2_1s, id_h2o, id_o, id_h /)

        has_hrates = all( ids(:) > 0 ) .and. ptop_ref < 0.0004_r8 * psurf_ref

        if (.not. has_hrates) return

        call addfld( 'CPAIR', 'J/K/kg', pver, 'I', 'specific heat cap air', phys_decomp )
        call addfld( 'QRS_AUR', 'K/s', pver, 'I', 'total auroral heating rate', phys_decomp )
        call addfld( 'QRS_CO2NIR', 'K/s', pver, 'I', 'co2 nir heating rate', phys_decomp )
        call addfld( 'QTHERMAL', 'K/s', pver, 'I', 'non-euv photolysis heating rate', phys_decomp )
        call addfld( 'QRS_MLT', 'K/s', pver, 'I', 'Total heating rate (unmerged with tropospheric RT heating)', phys_decomp )

        attr = 'O2 + hv -> O1D + O3P solar heating rate < 200nm'
        call addfld( 'QRS_SO2A', 'K/s ', pver, 'I', trim(attr), phys_decomp )
        attr = 'O2 + hv -> O3P + O3P solar heating rate < 200nm'
        call addfld( 'QRS_SO2B', 'K/s ', pver, 'I', trim(attr), phys_decomp )
        attr = 'O3 + hv -> O1D + O2_1S solar heating rate < 200nm'
        call addfld( 'QRS_SO3A', 'K/s ', pver, 'I', trim(attr), phys_decomp )
        attr = 'O3 + hv -> O3P + O2 solar heating rate < 200nm'
        call addfld( 'QRS_SO3B', 'K/s ', pver, 'I', trim(attr), phys_decomp )
        attr = 'O2 + hv -> 2*O3P solar heating rate > 200nm'
        call addfld( 'QRS_LO2B', 'K/s ', pver, 'I', trim(attr), phys_decomp )
        attr = 'O3 + hv -> O1D + O2_1S solar heating rate > 200nm'
        call addfld( 'QRS_LO3A', 'K/s ', pver, 'I', trim(attr), phys_decomp )
        attr = 'O3 + hv -> O3P + O2 solar heating rate > 200nm'
        call addfld( 'QRS_LO3B', 'K/s ', pver, 'I', trim(attr), phys_decomp )
        attr = 'Total O3 solar heating > 200nm'
        call addfld( 'QRS_LO3',  'K/s ', pver, 'I', trim(attr), phys_decomp )
        attr = 'total euv heating rate'
        call addfld( 'QRS_EUV', 'K/s', pver, 'I', trim(attr), phys_decomp )
        attr = 'total jo2 euv photolysis rate'
        call addfld( 'JO2_EUV', '/s', pver, 'I', trim(attr), phys_decomp )

      end subroutine init_hrates

! Subprogram not used       subroutine waccm_hrates(ncol, state, asdir, bot_mlt_lev, qrs_tot, pbuf )
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     ... computes the short wavelength heating rates
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       use chem_mods,         only : nabscol, nfs, gas_pcnst, rxntot, indexm
! Subprogram not used       use ppgrid,            only : pcols, pver
! Subprogram not used       use physconst,         only : rga, mbarv, cpairv
! Subprogram not used       use constituents,      only : pcnst
! Subprogram not used       use mo_gas_phase_chemdr,only: map2chm
! Subprogram not used       use mo_photo,          only : set_ub_col, setcol
! Subprogram not used       use mo_jlong,          only : jlong
! Subprogram not used       use mo_jshort,         only : jshort
! Subprogram not used       use mo_jeuv,           only : heuv
! Subprogram not used       use mo_cph,            only : cph
! Subprogram not used       use mo_heatnirco2,     only : heatnirco2
! Subprogram not used       use mo_airglow,        only : airglow
! Subprogram not used       use mo_aurora,         only : aurora
! Subprogram not used       use mo_setrxt,         only : setrxt_hrates
! Subprogram not used       use mo_adjrxt,         only : adjrxt
! Subprogram not used       use mo_usrrxt,         only : usrrxt_hrates
! Subprogram not used       use mo_setinv,         only : setinv
! Subprogram not used       use mo_mass_xforms,    only : mmr2vmr
! Subprogram not used       use physics_types,     only : physics_state
! Subprogram not used       use phys_grid,         only : get_rlat_all_p, get_rlon_all_p, &
! Subprogram not used                                     get_lat_all_p, get_lon_all_p
! Subprogram not used       use mo_mean_mass,      only : set_mean_mass
! Subprogram not used       use set_cp,            only : calc_cp
! Subprogram not used       use cam_history,       only : outfld
! Subprogram not used       use shr_orb_mod,       only : shr_orb_decl
! Subprogram not used       use time_manager,      only : get_curr_calday
! Subprogram not used       use cam_control_mod,   only : lambm0, eccen, mvelpp, obliqr
! Subprogram not used       use mo_constants,      only : r2d
! Subprogram not used       use short_lived_species,only: get_short_lived_species
! Subprogram not used       use physics_buffer,    only : physics_buffer_desc
! Subprogram not used       use phys_control,      only : waccmx_is
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !        ... dummy arguments
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       integer,             intent(in)  ::  ncol                  ! number columns in chunk
! Subprogram not used       type(physics_state), intent(in)  ::  state                 ! physics state structure
! Subprogram not used       real(r8),            intent(in)  ::  asdir(pcols)          ! shortwave, direct albedo
! Subprogram not used       integer,             intent(in)  ::  bot_mlt_lev           ! lowest model level where MLT heating is needed
! Subprogram not used       real(r8),            intent(out) ::  qrs_tot(pcols,pver)   ! total heating (K/s)
! Subprogram not used       type(physics_buffer_desc), pointer :: pbuf(:)
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !     	... local variables
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       integer             :: lchnk                 ! chunk index
! Subprogram not used       real(r8), parameter :: m2km  = 1.e-3_r8
! Subprogram not used       real(r8), parameter :: Pa2mb = 1.e-2_r8
! Subprogram not used 
! Subprogram not used       integer      ::  i, k, m, n
! Subprogram not used       integer      ::  kbot_hrates
! Subprogram not used       real(r8)     ::  esfact
! Subprogram not used       real(r8)     ::  sza                                           ! solar zenith angle (degrees)
! Subprogram not used       integer      ::  latndx(pcols)                                 ! chunk lat indicies
! Subprogram not used       integer      ::  lonndx(pcols)                                 ! chunk lon indicies
! Subprogram not used       real(r8)     ::  invariants(ncol,pver,nfs)
! Subprogram not used       real(r8)     ::  col_dens(ncol,pver,max(1,nabscol))            ! column densities (molecules/cm^2)
! Subprogram not used       real(r8)     ::  col_delta(ncol,0:pver,max(1,nabscol))         ! layer column densities (molecules/cm^2)
! Subprogram not used       real(r8)     ::  vmr(ncol,pver,gas_pcnst)                      ! xported species (vmr)
! Subprogram not used       real(r8)     ::  reaction_rates(ncol,pver,rxntot)              ! reaction rates
! Subprogram not used       real(r8)     ::  mmr(pcols,pver,gas_pcnst)                     ! chem working concentrations (kg/kg)
! Subprogram not used       real(r8)     ::  h2ovmr(ncol,pver)                             ! water vapor concentration (mol/mol)
! Subprogram not used       real(r8)     ::  mbar(ncol,pver)                               ! mean wet atmospheric mass (kg/mole)
! Subprogram not used       real(r8)     ::  zmid(ncol,pver)                               ! midpoint geopotential (km)
! Subprogram not used       real(r8)     ::  cpair(ncol,pver)                              ! specific heat capacity (J/K/kg)
! Subprogram not used       real(r8)     ::  cphrate(ncol,pver)                            ! chemical pot heat rate (K/s)
! Subprogram not used       real(r8)     ::  aghrate(ncol,pver)                            ! airglow heat rate (K/s)
! Subprogram not used       real(r8)     ::  qrs_col(pver,4)                               ! column thermal heating < 200nm
! Subprogram not used       real(r8)     ::  qrl_col(pver,4)                               ! column thermal heating > 200nm
! Subprogram not used       real(r8)     ::  qrs(ncol,pver,4)                              ! chunk thermal heating < 200nm
! Subprogram not used       real(r8)     ::  qrl(ncol,pver,4)                              ! chunk thermal heating > 200nm
! Subprogram not used       real(r8)     ::  euv_hrate_col(pver)                           ! column euv thermal heating rate
! Subprogram not used       real(r8)     ::  co2_hrate_col(pver)                           ! column co2 nir heating rate
! Subprogram not used       real(r8)     ::  euv_hrate(ncol,pver)                          ! chunk euv thermal heating rate
! Subprogram not used       real(r8)     ::  aur_hrate(ncol,pver)                          ! chunk auroral heating rate
! Subprogram not used       real(r8)     ::  co2_hrate(ncol,pver)                          ! chunk co2 nir heating rate
! Subprogram not used       real(r8)     ::  o2mmr(ncol,pver)                              ! chunk o2 concentration (kg/kg)
! Subprogram not used       real(r8)     ::  ommr(ncol,pver)                               ! chunk o concentration (kg/kg)
! Subprogram not used       real(r8)     ::  fac1(pver)                                    ! work array
! Subprogram not used       real(r8)     ::  colo3(pver)                                   ! vertical o3 column density
! Subprogram not used       real(r8)     ::  zarg(pver)                                    ! vertical height array
! Subprogram not used       real(r8)     ::  parg(pver)                                    ! vertical pressure array (hPa)
! Subprogram not used       real(r8)     ::  tline(pver)                                   ! vertical temperature array
! Subprogram not used       real(r8)     ::  eff_alb(pver)                                 ! albedo
! Subprogram not used       real(r8)     ::  mw(pver)                                      ! atms molecular weight
! Subprogram not used       real(r8)     ::  n2_line(pver)                                 ! n2 density (mol/mol)
! Subprogram not used       real(r8)     ::  o_line(pver)                                  ! o density (mol/mol)
! Subprogram not used       real(r8)     ::  o2_line(pver)                                 ! o2 density (mol/mol)
! Subprogram not used       real(r8)     ::  o3_line(pver)                                 ! o3 density (mol/mol)
! Subprogram not used       real(r8)     ::  co2_line(pver)                                ! co2 density (mol/mol)
! Subprogram not used       real(r8)     ::  scco2(pver)                                   ! co2 slant column concentration (molec/cm^2)
! Subprogram not used       real(r8)     ::  scco2i(pver)                                  ! co2 slant column concentration (molec/cm^2)
! Subprogram not used       real(r8)     ::  occ(pver)                                     ! o density (molecules/cm^3)
! Subprogram not used       real(r8)     ::  o2cc(pver)                                    ! o2 density (molecules/cm^3)
! Subprogram not used       real(r8)     ::  co2cc(pver)                                   ! co2 density (molecules/cm^3)
! Subprogram not used       real(r8)     ::  n2cc(pver)                                    ! n2 density (molecules/cm^3)
! Subprogram not used       real(r8)     ::  o3cc(pver)                                    ! o3 density (molecules/cm^3)
! Subprogram not used       real(r8)     ::  cparg(pver)                                   ! specific heat capacity
! Subprogram not used       real(r8)     ::  zen_angle(ncol)                               ! solar zenith angles (radians)
! Subprogram not used       real(r8)     ::  zsurf(ncol)                                   ! surface height (m)
! Subprogram not used       real(r8)     ::  rlats(ncol)                                   ! chunk latitudes (radians)
! Subprogram not used       real(r8)     ::  rlons(ncol)                                   ! chunk longitudes (radians)
! Subprogram not used       real(r8)     ::  calday                                        ! day of year
! Subprogram not used       real(r8)     ::  delta                                         ! solar declination (radians)
! Subprogram not used       logical      ::  do_diag
! Subprogram not used 
! Subprogram not used       qrs_tot(:ncol,:) = 0._r8
! Subprogram not used       if (.not. has_hrates) return
! Subprogram not used       
! Subprogram not used !-------------------------------------------------------------------------      
! Subprogram not used !        ... set maximum zenith angle - higher value for higher top model
! Subprogram not used !-------------------------------------------------------------------------      
! Subprogram not used       if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then 
! Subprogram not used          max_zen_angle = 116._r8
! Subprogram not used       else
! Subprogram not used          max_zen_angle = 97.01_r8 ! degrees
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used !        ... get chunk latitudes and longitudes
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used       lchnk = state%lchnk
! Subprogram not used 
! Subprogram not used       call get_lat_all_p( lchnk, ncol, latndx )
! Subprogram not used       call get_lon_all_p( lchnk, ncol, lonndx )
! Subprogram not used       call get_rlat_all_p( lchnk, ncol, rlats )
! Subprogram not used       call get_rlon_all_p( lchnk, ncol, rlons )
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used !        ... set lower limit for heating rates which is now dictated by radheat module
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used       kbot_hrates = bot_mlt_lev
! Subprogram not used       kbot_hrates = min( kbot_hrates,pver )
! Subprogram not used !     write(iulog,*) 'hrates: kbot_hrates = ',kbot_hrates
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used !        ... calculate cosine of zenith angle then cast back to angle
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used       calday = get_curr_calday()
! Subprogram not used       call zenith( calday, rlats, rlons, zen_angle, ncol )
! Subprogram not used       zen_angle(:) = acos( zen_angle(:) )
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used !        ... map incoming concentrations to working array
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used       do m = 1,pcnst
! Subprogram not used          n = map2chm(m)
! Subprogram not used          if( n > 0 ) then
! Subprogram not used             do k = 1,pver
! Subprogram not used                mmr(:ncol,k,n) = state%q(:ncol,k,m)
! Subprogram not used             end do
! Subprogram not used          end if
! Subprogram not used       end do
! Subprogram not used       call get_short_lived_species( mmr, lchnk, ncol, pbuf )
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used !        ... set atmosphere mean mass
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used       if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then 
! Subprogram not used         do k = 1,pver
! Subprogram not used           mbar(:ncol,k) = mbarv(:ncol,k,lchnk)
! Subprogram not used         enddo
! Subprogram not used       else      
! Subprogram not used         call set_mean_mass( ncol, mmr, mbar )
! Subprogram not used       endif
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used !        ... xform from mmr to vmr
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used       call mmr2vmr( mmr, vmr, mbar, ncol )
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used !        ... xform water vapor from mmr to vmr
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used       do k = 1,pver
! Subprogram not used          h2ovmr(:ncol,k) = vmr(:ncol,k,id_h2o)
! Subprogram not used       end do
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used !        ... xform geopotential height from m to km 
! Subprogram not used !            and pressure from Pa to mb
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used       zsurf(:ncol) = rga * state%phis(:ncol)
! Subprogram not used       do k = 1,pver
! Subprogram not used          zmid(:ncol,k) = m2km * (state%zm(:ncol,k) + zsurf(:ncol))
! Subprogram not used       end do
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used !        ... set the "invariants"
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used       call setinv( invariants, state%t, h2ovmr, vmr, state%pmid, ncol, lchnk, pbuf )
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used !        ... set the column densities at the upper boundary
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used       call set_ub_col( col_delta, vmr, invariants, state%pint(:,1), state%pdel, ncol, lchnk )
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used !       ...  set rates for "tabular" and user specified reactions
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used       do m = 1,rxntot
! Subprogram not used          do k = 1,pver
! Subprogram not used             reaction_rates(:,k,m) = 0._r8
! Subprogram not used          end do
! Subprogram not used       end do
! Subprogram not used       call setrxt_hrates( reaction_rates, state%t, invariants(1,1,indexm), ncol, kbot_hrates )
! Subprogram not used       call usrrxt_hrates( reaction_rates, state%t, state%t, state%t, invariants, &
! Subprogram not used                           h2ovmr, state%pmid, invariants(:,:,indexm), ncol, kbot_hrates )
! Subprogram not used       call adjrxt( reaction_rates, invariants, invariants(1,1,indexm), ncol )
! Subprogram not used       
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used !     	... set cp array
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used       if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then 
! Subprogram not used         do k = 1, pver
! Subprogram not used            cpair(:ncol,k) = cpairv(:ncol,k,lchnk)
! Subprogram not used         enddo
! Subprogram not used       else      
! Subprogram not used         call calc_cp( ncol, vmr, cpair )
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       call outfld( 'CPAIR', cpair, ncol, lchnk )
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used !     	... set the earth-sun distance factor
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used       call shr_orb_decl( calday, eccen, mvelpp, lambm0, obliqr  , &
! Subprogram not used                          delta, esfact )
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used !     	... set the column densities
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used       call setcol( col_delta, col_dens, vmr, state%pdel,  ncol )
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !        ... compute the thermal heating rates
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used       do m = 1,4
! Subprogram not used          do k = 1,pver
! Subprogram not used             qrs(:,k,m) = 0._r8
! Subprogram not used             qrl(:,k,m) = 0._r8
! Subprogram not used          end do
! Subprogram not used       end do
! Subprogram not used       do k = 1,pver
! Subprogram not used          euv_hrate(:,k) = 0._r8
! Subprogram not used          co2_hrate(:,k) = 0._r8
! Subprogram not used       end do
! Subprogram not used column_loop : &
! Subprogram not used       do i = 1,ncol
! Subprogram not used          sza = zen_angle(i)*r2d
! Subprogram not used          if( sza < max_zen_angle ) then
! Subprogram not used             zarg(:)     = zmid(i,:)
! Subprogram not used             parg(:)     = Pa2mb*state%pmid(i,:)
! Subprogram not used             colo3(:)    = col_dens(i,:,1)
! Subprogram not used             tline(:)    = state%t(i,:)
! Subprogram not used             eff_alb(:)  = asdir(i)
! Subprogram not used             o_line(:)   = vmr(i,:,id_o)
! Subprogram not used             o2_line(:)  = vmr(i,:,id_o2)
! Subprogram not used             co2_line(:) = vmr(i,:,id_co2)
! Subprogram not used             n2_line(:)  = 1._r8 - (o_line(:) + o2_line(:) + vmr(i,:,id_h))
! Subprogram not used             o3_line(:)  = vmr(i,:,id_o3)
! Subprogram not used             occ(:)      = o_line(:) * invariants(i,:,indexm)
! Subprogram not used             o2cc(:)     = o2_line(:) * invariants(i,:,indexm)
! Subprogram not used             co2cc(:)    = co2_line(:) * invariants(i,:,indexm)
! Subprogram not used             n2cc(:)     = n2_line(:) * invariants(i,:,indexm)
! Subprogram not used             o3cc(:)     = o3_line(:) * invariants(i,:,indexm)
! Subprogram not used             mw(:)       = mbar(i,:)
! Subprogram not used             cparg(:)    = cpair(i,:)
! Subprogram not used             do_diag     = .false.
! Subprogram not used             call jshort( pver, sza, o2_line, o3_line, o2cc, &
! Subprogram not used                          o3cc, tline, zarg, mw, qrs_col, &
! Subprogram not used                          cparg, lchnk, i, co2cc, scco2, do_diag )
! Subprogram not used             call jlong( pver, sza, eff_alb, parg, tline, &
! Subprogram not used                         mw, o2_line, o3_line, colo3, qrl_col, &
! Subprogram not used                         cparg, kbot_hrates )
! Subprogram not used             do m = 1,4
! Subprogram not used                qrs(i,pver:1:-1,m) = qrs_col(:,m) * esfact
! Subprogram not used             end do
! Subprogram not used             do m = 2,4
! Subprogram not used                qrl(i,:,m) = qrl_col(:,m) * esfact
! Subprogram not used             end do
! Subprogram not used             call heuv( pver, sza, occ, o2cc, n2cc, &
! Subprogram not used                        o_line, o2_line, n2_line, cparg, mw, &
! Subprogram not used                        zarg, euv_hrate_col, kbot_hrates )
! Subprogram not used             euv_hrate(i,:) = euv_hrate_col(:) * esfact
! Subprogram not used             scco2i(1:pver) = scco2(pver:1:-1)
! Subprogram not used             call heatnirco2( co2_line, scco2i, state%pmid(i,:kbot_hrates), co2_hrate_col, kbot_hrates, &
! Subprogram not used                              zarg, sza )
! Subprogram not used             co2_hrate(i,:kbot_hrates) = co2_hrate_col(:kbot_hrates) * esfact * daypsec
! Subprogram not used          end if
! Subprogram not used       end do column_loop
! Subprogram not used 
! Subprogram not used 
! Subprogram not used       call outfld( 'QRS_SO2A', qrs(:,:,1), ncol, lchnk )
! Subprogram not used       call outfld( 'QRS_SO2B', qrs(:,:,2), ncol, lchnk )
! Subprogram not used       call outfld( 'QRS_SO3A', qrs(:,:,3), ncol, lchnk )
! Subprogram not used       call outfld( 'QRS_SO3B', qrs(:,:,4), ncol, lchnk )
! Subprogram not used       call outfld( 'QRS_LO2B', qrl(:,:,2), ncol, lchnk )
! Subprogram not used       call outfld( 'QRS_LO3A', qrl(:,:,3), ncol, lchnk )
! Subprogram not used       call outfld( 'QRS_LO3B', qrl(:,:,4), ncol, lchnk )
! Subprogram not used       call outfld( 'QRS_LO3',  qrl(:,:,3)+qrl(:,:,4), ncol, lchnk )
! Subprogram not used       call outfld( 'QRS_EUV', euv_hrate(:,:), ncol, lchnk )
! Subprogram not used       call outfld( 'QRS_CO2NIR', co2_hrate(:,:), ncol, lchnk )
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used !     	... chemical pot heating rate
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used       call cph( cphrate, vmr, reaction_rates, cpair, mbar, &
! Subprogram not used                 kbot_hrates, ncol, lchnk )
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used !     	... auroral ion production
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used       do k = 1,pver
! Subprogram not used          o2mmr(:ncol,k) = mmr(:ncol,k,id_o2)
! Subprogram not used          ommr(:ncol,k)  = mmr(:ncol,k,id_o)
! Subprogram not used       end do
! Subprogram not used       call aurora( state%t, o2mmr, ommr, mbar, rlats, &
! Subprogram not used                    aur_hrate, cpair, state%pmid, lchnk, calday, &
! Subprogram not used                    ncol, rlons )
! Subprogram not used       do k = 1,pver
! Subprogram not used          aur_hrate(:,k)  = aur_hrate(:,k)/invariants(:,k,indexm)
! Subprogram not used       end do
! Subprogram not used       call outfld( 'QRS_AUR', aur_hrate(:,:), ncol, lchnk )
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used !     	... airglow heating rate
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used       call airglow( aghrate, vmr(1,1,id_o2_1s), vmr(1,1,id_o2_1d), reaction_rates, cpair, &
! Subprogram not used                     ncol, lchnk )
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used !     	... form total heating rate
! Subprogram not used !-----------------------------------------------------------------------      
! Subprogram not used       do k = 1,kbot_hrates
! Subprogram not used          qrs_tot(:ncol,k) = qrs(:,k,1) + qrs(:,k,2) + qrs(:,k,3) + qrs(:,k,4) &
! Subprogram not used                           + qrl(:,k,1) + qrl(:,k,2) + qrl(:,k,3) + qrl(:,k,4)
! Subprogram not used       end do
! Subprogram not used       call outfld( 'QTHERMAL', qrs_tot, pcols, lchnk )
! Subprogram not used       do k = 1,kbot_hrates
! Subprogram not used          qrs_tot(:ncol,k) = qrs_tot(:ncol,k) &
! Subprogram not used                           + cphrate(:,k) + euv_hrate(:,k) + aur_hrate(:,k) + co2_hrate(:,k)
! Subprogram not used       end do
! Subprogram not used       call outfld( 'QRS_MLT', qrs_tot, pcols, lchnk )
! Subprogram not used 
! Subprogram not used       end subroutine waccm_hrates

      end module mo_waccm_hrates
