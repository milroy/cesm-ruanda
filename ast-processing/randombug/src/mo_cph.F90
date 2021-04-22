!---------------------------------------------------------------------
!	... compute chemical potential heating
!---------------------------------------------------------------------

module mo_cph

  use shr_kind_mod,only : r8 => shr_kind_r8
  use chem_mods,   only : ncph=>enthalpy_cnt, exotherm=>cph_enthalpy, cph_rid

  implicit none

  save

  !==============================================================
  !... Doug Kinnison, dkin@ucar.edu
  !
  !... Enthalpy Data are taken from Atkinson et al., 
  !    Evaluated kinetic and photochemical data for atmospheric
  !    chemistry: Volume I, Atmos. Chem. Phys., 4, 1461-1738.
  !... Heats of formation at 0K.
  !... Units: kJ mol-1
  !
  !... Exception to the Atkinson et al. reference  (@0K unless noted)
  !    (4), (5), (8), (9), (10, (11), (14), (15), (27), (28) at 298K  
  !    (7)  h + o2 -> oh + o2 is multiplied by 0.6 (Mlynczak) to represent
  !         AG loss of excited OH.
  !    (25) n2d + o2 -> no + o1d taken from Roble, UMLT, Johnson and Killeen.
  !    (26) n2d + o  -> n  + o   taken from Roble, UMLT, Johnson and Killeen.
  !    (30-41) Taken from Roble, UMLT, Johnson and Killeen Ed., Geophys. Mono. 87
  !==============================================================
  !... Enthalpy Data are specified in preprocessor input mechanism file
  !... F Vitt -- 29 Oct 2013
  !==============================================================

  private
  public :: cph, init_cph

  logical :: has_cph
  character(len=24) :: fldnames(ncph)
  logical, parameter :: debug = .false.

contains

  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  subroutine init_cph

    use mo_chem_utls, only : get_rxt_ndx, get_spc_ndx
    use cam_history,  only : addfld, phys_decomp, add_default
    use ppgrid,       only : pver
    use chem_mods,    only : rxt_tag_lst, rxt_tag_map, rxt_tag_cnt
    use abortutils,   only : endrun

    implicit none

    character(len=64) :: longname
    integer :: i, n, tagndx

    has_cph = ncph > 0

    if (.not.has_cph) return

    if ( any(exotherm(:) == 0._r8) ) then
       call endrun('init_cph: Enthalpies for chemical heating must be specified in mechanism file')
    endif
    
    do i = 1,ncph

       findtagndx: do n = 1,rxt_tag_cnt
          if ( rxt_tag_map(n) == cph_rid(i) ) then
             tagndx = n
             exit findtagndx
          endif
       enddo findtagndx

       if (debug) then
          if ( i< 10 ) then
             write(fldnames(i),fmt='("CPH",i1)') i
          else if (i<100) then
             write(fldnames(i),fmt='("CPH",i2)') i
          else if (i<1000) then
             write(fldnames(i),fmt='("CPH",i3)') i
          endif
       else
          fldnames(i) = 'CPH_'//trim(rxt_tag_lst(tagndx))
       endif

       write( longname, fmt='(f12.6)') exotherm(i)
       longname = trim(adjustl(longname))//' kcal/mol chem pot heating rate for rxtn '//trim(rxt_tag_lst(tagndx))
       call addfld( fldnames(i), 'K/s ', pver, 'I', trim(longname), phys_decomp )
       if (debug) then
          call add_default( fldnames(i), 10, ' ' )
       endif
    enddo

    call addfld( 'QCP',   'K/s ', pver, 'I', 'chem pot heating rate', phys_decomp )

  end subroutine init_cph
  
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
! Subprogram not used   subroutine cph( cph_tot, vmr, rxt, cp, mbar, kbot, ncol, lchnk )
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !      	... forms the chemical potential heating rates
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     use chem_mods,     only : gas_pcnst, rxntot
! Subprogram not used     use ppgrid,        only : pver
! Subprogram not used     use cam_history,   only : outfld
! Subprogram not used     use mo_rxt_rates_conv, only : set_rates
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !     	... dummy arguments
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     integer, intent(in)   ::  ncol                                ! columns in chunck
! Subprogram not used     integer, intent(in)   ::  lchnk                               ! chunk index
! Subprogram not used     integer, intent(in)   ::  kbot                                ! bottom vert index
! Subprogram not used     real(r8), intent(in)  ::  rxt(ncol,pver,rxntot)               ! rxt rates (1/cm^3/s)
! Subprogram not used     real(r8), intent(in)  ::  cp(ncol,pver)                       ! specific heat capacity (J/K/kg)
! Subprogram not used     real(r8), intent(in)  ::  mbar(ncol,pver)                     ! atm mean mass (g/mole)
! Subprogram not used     real(r8), intent(in)  ::  vmr(ncol,pver,gas_pcnst)            ! concentrations (mol/mol)
! Subprogram not used     real(r8), intent(out) ::  cph_tot(ncol,pver)                  ! total heating (K/s)
! Subprogram not used 
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     !     	... local variables
! Subprogram not used     !-----------------------------------------------------------------------
! Subprogram not used     integer  ::  i, k
! Subprogram not used     real(r8) ::  tmp(ncol,pver)
! Subprogram not used     real(r8) ::  cph_rate(ncol,pver,ncph)
! Subprogram not used     real(r8) ::  rxt_rates(ncol,pver,rxntot)
! Subprogram not used 
! Subprogram not used     if (.not.has_cph) return
! Subprogram not used 
! Subprogram not used     ! get the reaction rates from rate constants ...
! Subprogram not used     rxt_rates(:ncol,:,:) = rxt(:ncol,:,:)
! Subprogram not used     call set_rates( rxt_rates, vmr, ncol )
! Subprogram not used 
! Subprogram not used     ! compute corresponding chem heating rates ...
! Subprogram not used     cph_rate(:,:,:) = 0._r8
! Subprogram not used     tmp(:ncol,:) =  1._r8 / (1.e-6_r8*cp(:ncol,:)*mbar(:ncol,:))
! Subprogram not used     do i = 1,ncph
! Subprogram not used        cph_rate(:ncol,:,i) = tmp(:ncol,:) * rxt_rates(:ncol,:,cph_rid(i)) * exotherm(i)
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     ! compute total heating rate ...
! Subprogram not used     cph_tot(:,:) = 0._r8
! Subprogram not used     do k = 1,kbot
! Subprogram not used        do i = 1,ncol
! Subprogram not used           cph_tot(i,k) = sum( cph_rate(i,k,:) )
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     ! output diagnostics
! Subprogram not used     do i = 1,ncph
! Subprogram not used        if ( exotherm(i)>0._r8) then
! Subprogram not used           call outfld( fldnames(i), cph_rate(:,:,i), ncol, lchnk )
! Subprogram not used        endif
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     call outfld( 'QCP', cph_tot(:,:), ncol, lchnk )
! Subprogram not used 
! Subprogram not used   end subroutine cph

end module mo_cph
