    module hirsbt

! ----- Modules -----
    use shr_kind_mod, only: r8 => shr_kind_r8
    use physconst, only: gravit, rair
    use ppgrid
    use hirsbtpar

    implicit none

! public subroutines
    public :: hirsrtm, hirsbt_init

    contains

! -------------------

! Subprogram not used     subroutine hirsrtm(lchnk, ncol, pp, tt, rmix, co2mmr, o3mix, ts, oro, &
! Subprogram not used                        tb_ir, britemp)
! Subprogram not used 
! Subprogram not used ! mji hirsrtm
! Subprogram not used ! Code includes modifications for F90 formatting and for application
! Subprogram not used ! to NCAR 1 based on the original code in CCM3.
! Subprogram not used ! M. J. Iacono, AER Inc., May 2004
! Subprogram not used ! Further structural revisions for CAM3.5
! Subprogram not used ! M. J. Iacono, AER Inc., April 2008
! Subprogram not used 
! Subprogram not used !      SUBROUTINE HIRSRTM(PP, TT, RMIX, O3MIX, TS, ORO,
! Subprogram not used !     $                   TB_IR,BRITEMP)
! Subprogram not used ! -- ------------------------------------------------------------------2
! Subprogram not used !
! Subprogram not used !                  ***       VERSION 2.0        ***
! Subprogram not used !
! Subprogram not used ! This subroutine calculates brightness temperatures at the top of the
! Subprogram not used ! atmosphere for 7 TOVS/HIRS channels (2,4,6,8,10,11,12) and 4 TOVS/MSU 
! Subprogram not used ! channels (1,2,3,4).
! Subprogram not used !
! Subprogram not used ! On Input:
! Subprogram not used !    pp      -  level pressure (hPa) at which fluxes are evaluated
! Subprogram not used !               from top of the atmosphere to surface
! Subprogram not used !    tt      -  layer temperatures (K)
! Subprogram not used !    rmix    -  layer H2O mass mixing ratio in kg/kg
! Subprogram not used !    co2mix  -  layer CO2 mass mixing ratio kg/kg 
! Subprogram not used !    o3mix   -  layer ozone mass mixing ratio kg/kg 
! Subprogram not used !    ts      -  surface temperature (K)
! Subprogram not used !    oro     - land-sea flag (sea=0, land=1)
! Subprogram not used !
! Subprogram not used ! On Ouput:
! Subprogram not used !    tb_ir   -  infrared brightness temperatures
! Subprogram not used !    britemp -  microwave brightness temperatures
! Subprogram not used !
! Subprogram not used !
! Subprogram not used !  A flag to include the 4 MSU channels can be switched on (msu_flag=1)
! Subprogram not used !  and off (msu_flag=0). To decrease the amount of computer time, the
! Subprogram not used !  microwave routine is changed to a lookup table with almost the same 
! Subprogram not used !  accuracy as the original routine.
! Subprogram not used !
! Subprogram not used ! **  last revised 3/31/97 by Richard Engelen **
! Subprogram not used !
! Subprogram not used !   This version differs from original version:
! Subprogram not used !
! Subprogram not used !     1.  New NOAA10 coefficients
! Subprogram not used !     2.  Continuum added
! Subprogram not used !     3.  Any level exceeding 100% RH is changed to 100% RH
! Subprogram not used !     4.  New channels (2,4,6,8,10,11,12)
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used !------------------------------Arguments--------------------------------
! Subprogram not used !
! Subprogram not used ! Input arguments
! Subprogram not used !
! Subprogram not used    integer, intent(in) :: lchnk                 ! chunk identifier
! Subprogram not used    integer, intent(in) :: ncol                  ! number of atmospheric columns
! Subprogram not used 
! Subprogram not used    real(r8), intent(in) :: pp(pcols,pverp)      ! Pressure for each column and layer (hPa)
! Subprogram not used    real(r8), intent(in) :: tt(pcols,pver)       ! Temperature of each column and layer (K)
! Subprogram not used    real(r8), intent(in) :: ts(pcols)            ! Surface temperature (K)
! Subprogram not used    real(r8), intent(in) :: rmix(pcols,pver)     ! Water vapor mass mixing ratio (kg/kg)
! Subprogram not used    real(r8), intent(in) :: co2mmr(pcols)        ! CO2 mass mixing ratio (kg/kg)
! Subprogram not used    real(r8), intent(in) :: o3mix(pcols,pver)    ! Ozone mass mixing ratio (kg/kg)
! Subprogram not used    real(r8), intent(in) :: oro(pcols)           ! Land surface flag, sea=0, land=1
! Subprogram not used !
! Subprogram not used ! Output arguments
! Subprogram not used !
! Subprogram not used    real(r8), intent(out) :: britemp(pcols,pnf_msu)   ! HIRS brightness temperatures
! Subprogram not used    real(r8), intent(out) :: tb_ir(pcols,pnb_hirs)    ! MSU brightness temperatures
! Subprogram not used 
! Subprogram not used !------------------------------Local variables--------------------------
! Subprogram not used 
! Subprogram not used !      real(r8), parameter :: rmwco2 = mwco2/mwdry   ! ratio of molecular weights of co2 to dry air
! Subprogram not used 
! Subprogram not used       real(r8)  uh2o(pver)
! Subprogram not used       real(r8)  uo3(pver)
! Subprogram not used       real(r8)  a_hirs(pnb_hirs)
! Subprogram not used       real(r8)  b_hirs(pnb_hirs)
! Subprogram not used       real(r8)  xh2o(pnb_hirs)
! Subprogram not used       real(r8)  yh2o(pnb_hirs)
! Subprogram not used       real(r8)  xo3(pnb_hirs)
! Subprogram not used       real(r8)  yo3(pnb_hirs)
! Subprogram not used       real(r8)  xco2(pnb_hirs)
! Subprogram not used       real(r8)  yco2(pnb_hirs)
! Subprogram not used       real(r8)  b_ir(pnb_hirs)
! Subprogram not used       real(r8)  otrans(pnb_hirs)
! Subprogram not used       real(r8)  tband(pnb_hirs)
! Subprogram not used       real(r8)  scoef10(pnb_hirs)
! Subprogram not used       real(r8)  fcoef10(pnb_hirs)
! Subprogram not used       real(r8)  cwn(pnb_hirs)
! Subprogram not used       real(r8)  dtrans(pnb_hirs)
! Subprogram not used       real(r8)  rad_lay(pnb_hirs)
! Subprogram not used       real(r8)  radir(pnb_hirs)
! Subprogram not used       real(r8)  rad(pnb_hirs)
! Subprogram not used       real(r8)  rad2(pnb_hirs)
! Subprogram not used       real(r8)  rad3(pnb_hirs)
! Subprogram not used       real(r8)  refl(pnb_hirs)
! Subprogram not used       real(r8)  otr_mw(pnf_msu)
! Subprogram not used       real(r8)  tau(pnf_msu)
! Subprogram not used       real(r8)  trans(pnf_msu)
! Subprogram not used       real(r8)  otr_mw2(pnf_msu)
! Subprogram not used       real(r8)  tau2(pnf_msu)
! Subprogram not used       real(r8)  trans2(pnf_msu)
! Subprogram not used !     real(r8)  britemp1(pcols,pnf_msu)
! Subprogram not used !     real(r8)  britemp2(pcols,pnf_msu)
! Subprogram not used !     real(r8)  britemp3(pcols,pnf_msu)
! Subprogram not used       real(r8)  freq(pnf_msu)
! Subprogram not used       real(r8)  upath_h2o
! Subprogram not used       real(r8)  upath_co2
! Subprogram not used       real(r8)  upath_o3
! Subprogram not used       real(r8)  ucont1
! Subprogram not used       real(r8)  ucont2
! Subprogram not used       real(r8)  ppath_h2o
! Subprogram not used       real(r8)  ppath_co2
! Subprogram not used       real(r8)  ppath_o3
! Subprogram not used       real(r8)  sfctemp                     ! surface temperature
! Subprogram not used       real(r8)  dp                          ! pressure depth of layer
! Subprogram not used       real(r8)  tlay                        ! layer temperature
! Subprogram not used       real(r8)  play                        ! pressure
! Subprogram not used       real(r8)  rlay                        ! water vapor mixing ratio
! Subprogram not used       real(r8)  o3lay                       ! ozone mixing ratio.
! Subprogram not used       real(r8)  rhoair                      ! layer density
! Subprogram not used       real(r8)  e                           ! saturation pressure
! Subprogram not used       real(r8)  delz                        ! height of layer
! Subprogram not used       real(r8)  dp2                         ! pressure depth of layer below
! Subprogram not used       real(r8)  tlay2
! Subprogram not used       real(r8)  play2
! Subprogram not used       real(r8)  rlay2
! Subprogram not used       real(r8)  rhoair2
! Subprogram not used       real(r8)  e2
! Subprogram not used       real(r8)  delz2
! Subprogram not used       real(r8)  b_mw
! Subprogram not used       real(r8)  b_mw2
! Subprogram not used       real(r8)  dtr_mw
! Subprogram not used       real(r8)  uco2(pver)
! Subprogram not used !      real(r8)  uco2
! Subprogram not used       real(r8)  pw
! Subprogram not used       real(r8)  psc_h2o                   ! partial pressure of h2o
! Subprogram not used       real(r8)  psc_co2                   ! partial pressure of co2
! Subprogram not used       real(r8)  psc_o3                    ! partial pressure of o3
! Subprogram not used       real(r8)  t_cont
! Subprogram not used       real(r8)  t_h2o
! Subprogram not used       real(r8)  t_co2
! Subprogram not used       real(r8)  t_o3
! Subprogram not used       real(r8)  abs(pnf_msu)              ! absorption
! Subprogram not used       real(r8)  abs2(pnf_msu)
! Subprogram not used       real(r8)  dtr_mw2
! Subprogram not used !
! Subprogram not used !      real(r8) t_malk, plnck, btemp
! Subprogram not used !      external t_malk, plnck, btemp
! Subprogram not used       integer  i, ib, if, icol    ! loop control
! Subprogram not used 
! Subprogram not used !------------------------------data statments---------------------------
! Subprogram not used       data scoef10/3.09_r8,2.64_r8,2.28_r8,1.004_r8,0.429_r8,1.945_r8,11.95_r8/
! Subprogram not used       data fcoef10/.00434_r8,.00301_r8,.0018_r8,4.33e-5_r8,0.000393_r8,0.0738_r8,1.110_r8/
! Subprogram not used       data a_hirs/.0183_r8,-.00203_r8,.0653_r8,0.21797_r8,0.29846_r8,0.04612_r8,0.06453_r8/
! Subprogram not used       data b_hirs/.99992_r8,.99994_r8,.9998_r8,.99957_r8,.9996_r8,.99963_r8,1.0006_r8/
! Subprogram not used       data cwn/680.23_r8,704.33_r8,733.13_r8,899.5_r8,1224.07_r8,1363.32_r8,1489.42_r8/
! Subprogram not used       data freq/50.31_r8,53.73_r8,54.96_r8,57.95_r8/
! Subprogram not used 
! Subprogram not used 
! Subprogram not used       data xh2o/0.41_r8,0.52_r8,0.17_r8,0.05_r8,0.68_r8,5.02_r8,17.24_r8/
! Subprogram not used       data yh2o/0.70_r8,1.41_r8,0.14_r8,0.02_r8,1.09_r8,125.9_r8,1194.5_r8/
! Subprogram not used 
! Subprogram not used       data xco2/64.32_r8,13.29_r8,5.04_r8,0.05_r8,0.03_r8,0.25_r8,0.01_r8/
! Subprogram not used       data yco2/1325.2_r8,121.01_r8,14.36_r8,0.01_r8,0.03_r8,0.02_r8,0.01_r8/
! Subprogram not used 
! Subprogram not used       data xo3/32.83_r8,28.9_r8,27.44_r8,0.01_r8,0.68_r8,0.37_r8,0.01_r8/
! Subprogram not used       data yo3/77.44_r8,66.7_r8,67.44_r8,1.6_r8,4.17_r8,0.07_r8,0.01_r8/
! Subprogram not used 
! Subprogram not used       real(r8) grav, r
! Subprogram not used !      real(r8) grav, rco2, r
! Subprogram not used !      data grav,rco2,r/9.8,0.523e-3,287./
! Subprogram not used       real(r8) eps
! Subprogram not used       data eps/0.5_r8/
! Subprogram not used 
! Subprogram not used ! use values for constants consistent with radiation code.
! Subprogram not used       grav = gravit
! Subprogram not used       r = rair
! Subprogram not used !      rco2 = co2vmr*rmwco2
! Subprogram not used 
! Subprogram not used       do icol=1,ncol
! Subprogram not used 
! Subprogram not used         sfctemp=ts(icol)
! Subprogram not used 
! Subprogram not used         do ib=1,pnb_hirs
! Subprogram not used           radir(ib)=0._r8
! Subprogram not used         enddo
! Subprogram not used 
! Subprogram not used         upath_h2o=0.0_r8
! Subprogram not used         upath_co2=0.0_r8
! Subprogram not used         upath_o3=0.0_r8
! Subprogram not used         ucont1=0._r8
! Subprogram not used         ucont2=0._r8
! Subprogram not used         ppath_h2o=0.0_r8
! Subprogram not used         ppath_co2=0.0_r8
! Subprogram not used         ppath_o3=0._r8
! Subprogram not used !       if(msu_flag.eq.1) then
! Subprogram not used           do if=1,pnf_msu
! Subprogram not used              tau(if)=0._r8
! Subprogram not used              trans(if)=0._r8
! Subprogram not used              rad(if)=0._r8
! Subprogram not used              otr_mw(if)=1._r8
! Subprogram not used              tau2(if)=0._r8
! Subprogram not used              trans2(if)=0._r8
! Subprogram not used              rad2(if)=0._r8
! Subprogram not used              otr_mw2(if)=1._r8
! Subprogram not used           end do
! Subprogram not used !       endif
! Subprogram not used         do ib=1,pnb_hirs
! Subprogram not used            tband(ib)=1.0_r8
! Subprogram not used            otrans(ib)=1._r8
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used 
! Subprogram not used         do i=1,pver
! Subprogram not used       
! Subprogram not used           dp=pp(icol,i+1)-pp(icol,i)
! Subprogram not used           tlay=tt(icol,i)
! Subprogram not used           play=sqrt(pp(icol,i)*pp(icol,i+1))
! Subprogram not used           rlay=rmix(icol,i)
! Subprogram not used           o3lay=o3mix(icol,i)
! Subprogram not used           
! Subprogram not used           rhoair=play*100._r8/(r*tt(icol,i))
! Subprogram not used           e = play*rlay/(rlay+0.6220_r8)
! Subprogram not used           delz=dp*0.1_r8/(rhoair*grav)
! Subprogram not used           
! Subprogram not used !         if(msu_flag.eq.1) then
! Subprogram not used             dp2=pp(icol,pver+2-i)-pp(icol,pver+1-i)
! Subprogram not used             tlay2=tt(icol,pver+1-i)
! Subprogram not used             play2=sqrt(pp(icol,pver+1-i)*pp(icol,pver+2-i))
! Subprogram not used             rlay2=rmix(icol,pver+1-i)
! Subprogram not used 
! Subprogram not used             rhoair2=play2*100._r8/(r*tt(icol,pver+1-i))
! Subprogram not used             e2 = play2*rlay2/(rlay2+0.6220_r8)
! Subprogram not used             delz2=dp2*0.1_r8/(rhoair2*grav)
! Subprogram not used           
! Subprogram not used !
! Subprogram not used !         microwave transfer
! Subprogram not used !
! Subprogram not used             call lookup ( play-e, tlay, abs)
! Subprogram not used             call lookup (play2-e2, tlay2, abs2)
! Subprogram not used             do if=1,pnf_msu
! Subprogram not used               tau(if)=tau(if)+abs(if)*delz
! Subprogram not used               trans(if) = exp(-tau(if))
! Subprogram not used               b_mw = 1.47445e-23_r8*freq(if)**3/(exp(0.047981_r8*freq(if) &
! Subprogram not used                      /tlay)-1)
! Subprogram not used               dtr_mw = otr_mw(if)-trans(if)
! Subprogram not used               rad(if)=rad(if) + b_mw*dtr_mw
! Subprogram not used               otr_mw(if)=trans(if)
! Subprogram not used            
! Subprogram not used               tau2(if)=tau2(if)+abs2(if)*delz2
! Subprogram not used               trans2(if) = exp(-tau2(if))
! Subprogram not used               b_mw2 = 1.47445e-23_r8*freq(if)**3/(exp(0.047981_r8*freq(if) &
! Subprogram not used                    /tlay2)-1)
! Subprogram not used               dtr_mw2 = otr_mw2(if)-trans2(if)
! Subprogram not used               rad2(if)=rad2(if) + b_mw2*dtr_mw2
! Subprogram not used               otr_mw2(if)=trans2(if)
! Subprogram not used             end do
! Subprogram not used !         endif
! Subprogram not used 
! Subprogram not used !
! Subprogram not used !                ir transfer
! Subprogram not used !
! Subprogram not used           uh2o(i)=rlay*100._r8*dp/grav
! Subprogram not used           uo3(i)=o3lay*100._r8*dp/grav
! Subprogram not used           uco2(i)=co2mmr(icol)*100._r8*dp/grav
! Subprogram not used !          uco2=rco2*100.*dp/grav
! Subprogram not used           pw=play*rlay*28.97_r8/18._r8
! Subprogram not used           ucont1=ucont1+uh2o(i)*((play-pw)/1013._r8)*(296._r8/tlay)
! Subprogram not used           ucont2=ucont2+uh2o(i)*(pw/1013._r8)*(296._r8/tlay)
! Subprogram not used 
! Subprogram not used           upath_h2o=upath_h2o+uh2o(i)
! Subprogram not used           ppath_h2o=ppath_h2o+uh2o(i)*play
! Subprogram not used           psc_h2o=ppath_h2o/upath_h2o
! Subprogram not used 
! Subprogram not used           upath_co2=upath_co2+uco2(i)
! Subprogram not used           ppath_co2=ppath_co2+uco2(i)*play
! Subprogram not used !          upath_co2=upath_co2+uco2
! Subprogram not used !          ppath_co2=ppath_co2+uco2*play
! Subprogram not used           psc_co2=ppath_co2/upath_co2
! Subprogram not used 
! Subprogram not used           upath_o3=upath_o3+uo3(i)
! Subprogram not used           ppath_o3=ppath_o3+uo3(i)*play
! Subprogram not used           psc_o3=ppath_o3/upath_o3
! Subprogram not used 
! Subprogram not used           do ib=1,pnb_hirs
! Subprogram not used             t_cont=exp(-scoef10(ib)*ucont2-fcoef10(ib)*ucont1)
! Subprogram not used 
! Subprogram not used             t_h2o=t_malk(upath_h2o,psc_h2o,xh2o(ib),yh2o(ib))
! Subprogram not used             t_co2=t_malk(upath_co2,psc_co2,xco2(ib),yco2(ib))
! Subprogram not used             t_o3=t_malk(upath_o3,psc_o3,xo3(ib),yo3(ib))
! Subprogram not used             tband(ib)=t_co2*t_o3*t_h2o*t_cont
! Subprogram not used 
! Subprogram not used 
! Subprogram not used             b_ir(ib)=plnck(cwn(ib),a_hirs(ib),b_hirs(ib),tlay)
! Subprogram not used             dtrans(ib)=otrans(ib)-tband(ib)
! Subprogram not used             rad_lay(ib)=  b_ir(ib)*dtrans(ib)
! Subprogram not used             radir(ib)= radir(ib)+ rad_lay(ib)
! Subprogram not used             otrans(ib) = tband(ib)
! Subprogram not used 
! Subprogram not used           end do
! Subprogram not used         end do     ! end of loop over vertical levels
! Subprogram not used 
! Subprogram not used 
! Subprogram not used !
! Subprogram not used !    add in the surface contribution to the radiance, including
! Subprogram not used !    reflection
! Subprogram not used !
! Subprogram not used         
! Subprogram not used !       if(msu_flag.eq.1) then
! Subprogram not used           if(oro(icol).eq.0)then
! Subprogram not used             eps=0.65_r8      ! ocean
! Subprogram not used           else
! Subprogram not used             eps=0.9_r8    ! land or sea-ice
! Subprogram not used           end if
! Subprogram not used       
! Subprogram not used           do if=1,pnf_msu
! Subprogram not used             b_mw = 1.47445e-23_r8*freq(if)**3/(exp(0.047981_r8*freq(if)/ &
! Subprogram not used                sfctemp)-1)
! Subprogram not used             rad(if) = rad(if) + eps*trans(if)*b_mw
! Subprogram not used             refl(if)=(1-eps)*rad2(if)*trans(if)
! Subprogram not used !           britemp1(icol,if) = 0.047981*freq(if)/log(1.0 
! Subprogram not used !    $                + 1.47445e-23*freq(if)**3/rad(if))
! Subprogram not used !           britemp2(icol,if) = 0.047981*freq(if)/log(1.0 
! Subprogram not used !    $                + 1.47445e-23*freq(if)**3/rad2(if))
! Subprogram not used !           britemp3(icol,if) = 0.047981*freq(if)/log(1.0 
! Subprogram not used !    $                + 1.47445e-23*freq(if)**3/refl(if))
! Subprogram not used      
! Subprogram not used             rad3(if) = rad(if) + refl(if)
! Subprogram not used             britemp(icol,if) = 0.047981_r8*freq(if)/log(1.0_r8  &
! Subprogram not used                       + 1.47445e-23_r8*freq(if)**3/rad3(if))
! Subprogram not used           end do
! Subprogram not used !       endif
! Subprogram not used       
! Subprogram not used         do ib=1,pnb_hirs
! Subprogram not used           b_ir(ib)=plnck(cwn(ib),a_hirs(ib),b_hirs(ib),sfctemp)
! Subprogram not used           radir(ib)=radir(ib)+b_ir(ib)*tband(ib)
! Subprogram not used           tb_ir(icol,ib)=btemp(cwn(ib),a_hirs(ib),b_hirs(ib),radir(ib))
! Subprogram not used         end do
! Subprogram not used 
! Subprogram not used       end do   ! end of loop over columns
! Subprogram not used 
! Subprogram not used     end subroutine hirsrtm
!
!==================================================================
!
! Subprogram not used     function btemp(wvn,a,b,chnrad)
! Subprogram not used !
! Subprogram not used !*  calculates the brightness temperature given the channel radiance.
! Subprogram not used !   uses planck function tuned to tovs frequencies
! Subprogram not used !
! Subprogram not used       real(r8) btemp
! Subprogram not used 
! Subprogram not used !------------------------------arguments--------------------------------
! Subprogram not used       real(r8) wvn
! Subprogram not used       real(r8) chnrad
! Subprogram not used       real(r8) a
! Subprogram not used       real(r8) b
! Subprogram not used !------------------------------local variables--------------------------
! Subprogram not used       real(r8) c1
! Subprogram not used       real(r8) c2
! Subprogram not used !
! Subprogram not used !   planck function parameters
! Subprogram not used !   note: parameters a and b are temperature correction factors
! Subprogram not used !   which are dependent on the channel and satellite.  these
! Subprogram not used !   parameters were extracted from the rttovs model.
! Subprogram not used !
! Subprogram not used       parameter( c1=1.191066e-08_r8 )
! Subprogram not used       parameter( c2=1.438833_r8 )
! Subprogram not used 
! Subprogram not used       btemp=(c2*wvn/log(c1*wvn**3/chnrad+1._r8)-a)/b
! Subprogram not used 
! Subprogram not used     end function btemp
!
!==================================================================
!
! Subprogram not used     function plnck(wvn,a,b,t)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! planck function
! Subprogram not used !
! Subprogram not used       real(r8) plnck
! Subprogram not used 
! Subprogram not used       real(r8) wvn,a,b,t,c1,c2
! Subprogram not used 
! Subprogram not used       parameter( c1=1.191066e-08_r8 )
! Subprogram not used       parameter( c2=1.438833_r8 )
! Subprogram not used 
! Subprogram not used       plnck=c1*(wvn)**3/(exp(c2*wvn/(a+b*t))-1._r8)
! Subprogram not used 
! Subprogram not used     end function plnck
!
!===================================================================
!
! Subprogram not used     function t_malk(u,p,x,y)
! Subprogram not used 
! Subprogram not used       real(r8) t_malk
! Subprogram not used 
! Subprogram not used       real(r8) u, p, x, y
! Subprogram not used       real(r8) p0, dnu, pi, b, bp
! Subprogram not used 
! Subprogram not used       parameter( p0=1013._r8 )
! Subprogram not used       parameter( dnu=10._r8 )
! Subprogram not used 
! Subprogram not used       pi=acos(-1._r8)
! Subprogram not used       b=4._r8*x**2/(pi*y*dnu)
! Subprogram not used       bp=pi*b*p/p0
! Subprogram not used       t_malk = exp(-0.5_r8*bp*(sqrt(1._r8+4._r8*y*u/(dnu*bp))-1._r8))
! Subprogram not used 
! Subprogram not used     end function t_malk
!
!===================================================================
!
! Subprogram not used     subroutine lookup(p,t,abs)
! Subprogram not used 
! Subprogram not used       real(r8) p, t, abs(pnf_msu)
! Subprogram not used !
! Subprogram not used       integer if, i, j, n, m
! Subprogram not used       parameter( n=17 )
! Subprogram not used       parameter( m=16 )
! Subprogram not used       real(r8) xx(n), yy(m)
! Subprogram not used       real(r8) zz1(n,m), zz2(n,m), zz3(n,m), zz4(n,m)
! Subprogram not used       real(r8) t1, t2
! Subprogram not used       data yy/5.0_r8,7.5_r8,10._r8,25._r8,50._r8,100._r8,200._r8,300._r8, &
! Subprogram not used               400._r8,500._r8,600._r8,700._r8,800._r8,900._r8,1000,1050._r8/
! Subprogram not used       data xx/160._r8,170._r8,180._r8,190._r8,200._r8,210._r8,220._r8,230._r8,240._r8,250._r8, &
! Subprogram not used               260._r8,270._r8,280._r8,290._r8,300._r8,310._r8,320._r8/
! Subprogram not used       data zz1/0.0000_r8,0.0000_r8,0.0000_r8,0.0000_r8,0.0000_r8,0.0000_r8,0.0000_r8, &
! Subprogram not used                0.0000_r8,0.0000_r8,0.0000_r8,0.0000_r8,0.0000_r8,0.0000_r8,0.0000_r8, &
! Subprogram not used                0.0000_r8,0.0000_r8,0.0000_r8,0.0000_r8,0.0000_r8,0.0000_r8,0.0000_r8, &
! Subprogram not used                0.0000_r8,0.0000_r8,0.0000_r8,0.0000_r8,0.0000_r8,0.0000_r8,0.0000_r8, &
! Subprogram not used                0.0000_r8,0.0000_r8,0.0000_r8,0.0000_r8,0.0000_r8,0.0000_r8,0.0000_r8, &
! Subprogram not used                0.0000_r8,0.0000_r8,0.0000_r8,0.0000_r8,0.0000_r8,0.0000_r8,0.0000_r8, &
! Subprogram not used                0.0000_r8,0.0000_r8,0.0000_r8,0.0000_r8,0.0000_r8,0.0000_r8,0.0000_r8, &
! Subprogram not used                0.0000_r8,0.0000_r8,0.0002_r8,0.0002_r8,0.0001_r8,0.0001_r8,0.0001_r8, &
! Subprogram not used                0.0001_r8,0.0001_r8,0.0001_r8,0.0001_r8,0.0001_r8,0.0001_r8,0.0000_r8, &
! Subprogram not used                0.0000_r8,0.0000_r8,0.0000_r8,0.0000_r8,0.0000_r8,0.0008_r8,0.0007_r8, &
! Subprogram not used                0.0006_r8,0.0005_r8,0.0004_r8,0.0004_r8,0.0003_r8,0.0003_r8,0.0003_r8, &
! Subprogram not used                0.0002_r8,0.0002_r8,0.0002_r8,0.0002_r8,0.0002_r8,0.0002_r8,0.0001_r8, &
! Subprogram not used                0.0001_r8,0.0032_r8,0.0027_r8,0.0023_r8,0.0020_r8,0.0017_r8,0.0015_r8, &
! Subprogram not used                0.0013_r8,0.0011_r8,0.0010_r8,0.0009_r8,0.0008_r8,0.0008_r8,0.0007_r8, &
! Subprogram not used                0.0006_r8,0.0006_r8,0.0006_r8,0.0005_r8,0.0130_r8,0.0108_r8,0.0091_r8, &
! Subprogram not used                0.0078_r8,0.0068_r8,0.0059_r8,0.0052_r8,0.0046_r8,0.0041_r8,0.0037_r8, &
! Subprogram not used                0.0033_r8,0.0030_r8,0.0028_r8,0.0025_r8,0.0023_r8,0.0022_r8,0.0020_r8, &
! Subprogram not used                0.0292_r8,0.0244_r8,0.0206_r8,0.0176_r8,0.0152_r8,0.0133_r8,0.0117_r8, &
! Subprogram not used                0.0103_r8,0.0092_r8,0.0083_r8,0.0074_r8,0.0067_r8,0.0061_r8,0.0056_r8, &
! Subprogram not used                0.0052_r8,0.0048_r8,0.0045_r8,0.0521_r8,0.0434_r8,0.0367_r8,0.0314_r8, &
! Subprogram not used                0.0271_r8,0.0237_r8,0.0208_r8,0.0184_r8,0.0164_r8,0.0147_r8,0.0132_r8, &
! Subprogram not used                0.0120_r8,0.0109_r8,0.0100_r8,0.0092_r8,0.0085_r8,0.0079_r8,0.0818_r8, &
! Subprogram not used                0.0681_r8,0.0576_r8,0.0493_r8,0.0426_r8,0.0371_r8,0.0326_r8,0.0289_r8, &
! Subprogram not used                0.0257_r8,0.0230_r8,0.0207_r8,0.0188_r8,0.0170_r8,0.0156_r8,0.0143_r8, &
! Subprogram not used                0.0132_r8,0.0122_r8,0.1182_r8,0.0985_r8,0.0833_r8,0.0712_r8,0.0616_r8, &
! Subprogram not used                0.0537_r8,0.0472_r8,0.0418_r8,0.0372_r8,0.0333_r8,0.0300_r8,0.0271_r8, &
! Subprogram not used                0.0246_r8,0.0225_r8,0.0206_r8,0.0189_r8,0.0175_r8,0.1617_r8,0.1348_r8, &
! Subprogram not used                0.1139_r8,0.0975_r8,0.0842_r8,0.0734_r8,0.0645_r8,0.0571_r8,0.0508_r8, &
! Subprogram not used                0.0455_r8,0.0409_r8,0.0370_r8,0.0336_r8,0.0306_r8,0.0281_r8,0.0258_r8, &
! Subprogram not used                0.0238_r8,0.2123_r8,0.1770_r8,0.1496_r8,0.1280_r8,0.1106_r8,0.0965_r8, &
! Subprogram not used                0.0848_r8,0.0750_r8,0.0667_r8,0.0597_r8,0.0537_r8,0.0486_r8,0.0441_r8, &
! Subprogram not used                0.0402_r8,0.0368_r8,0.0338_r8,0.0312_r8,0.2701_r8,0.2253_r8,0.1905_r8, &
! Subprogram not used                0.1630_r8,0.1408_r8,0.1228_r8,0.1079_r8,0.0954_r8,0.0849_r8,0.0760_r8, &
! Subprogram not used                0.0684_r8,0.0618_r8,0.0560_r8,0.0511_r8,0.0467_r8,0.0429_r8,0.0396_r8, &
! Subprogram not used                0.3353_r8,0.2798_r8,0.2366_r8,0.2024_r8,0.1749_r8,0.1525_r8,0.1340_r8, &
! Subprogram not used                0.1185_r8,0.1055_r8,0.0944_r8,0.0849_r8,0.0767_r8,0.0695_r8,0.0634_r8, &
! Subprogram not used                0.0579_r8,0.0532_r8,0.0490_r8,0.3707_r8,0.3094_r8,0.2617_r8,0.2239_r8, &
! Subprogram not used                0.1935_r8,0.1687_r8,0.1482_r8,0.1311_r8,0.1166_r8,0.1043_r8,0.0938_r8, &
! Subprogram not used                0.0848_r8,0.0769_r8,0.0700_r8,0.0640_r8,0.0588_r8,0.0541_r8/ 
! Subprogram not used       data zz2/0.0000_r8,0.0000_r8,0.0000_r8,0.0000_r8,0.0000_r8,0.0000_r8,0.0001_r8, &
! Subprogram not used                0.0001_r8,0.0001_r8,0.0001_r8,0.0001_r8,0.0001_r8,0.0001_r8,0.0001_r8, &
! Subprogram not used                0.0001_r8,0.0001_r8,0.0001_r8,0.0001_r8,0.0001_r8,0.0001_r8,0.0001_r8, &
! Subprogram not used                0.0001_r8,0.0001_r8,0.0001_r8,0.0001_r8,0.0001_r8,0.0001_r8,0.0001_r8, &
! Subprogram not used                0.0001_r8,0.0001_r8,0.0001_r8,0.0001_r8,0.0001_r8,0.0001_r8,0.0002_r8, &
! Subprogram not used                0.0002_r8,0.0002_r8,0.0002_r8,0.0002_r8,0.0002_r8,0.0002_r8,0.0002_r8, &
! Subprogram not used                0.0002_r8,0.0002_r8,0.0002_r8,0.0002_r8,0.0002_r8,0.0002_r8,0.0002_r8, &
! Subprogram not used                0.0002_r8,0.0002_r8,0.0010_r8,0.0010_r8,0.0010_r8,0.0011_r8,0.0011_r8, &
! Subprogram not used                0.0012_r8,0.0012_r8,0.0013_r8,0.0013_r8,0.0014_r8,0.0014_r8,0.0014_r8, &
! Subprogram not used                0.0015_r8,0.0015_r8,0.0015_r8,0.0015_r8,0.0015_r8,0.0040_r8,0.0038_r8, &
! Subprogram not used                0.0038_r8,0.0039_r8,0.0041_r8,0.0042_r8,0.0044_r8,0.0046_r8,0.0048_r8, &
! Subprogram not used                0.0050_r8,0.0051_r8,0.0053_r8,0.0054_r8,0.0055_r8,0.0056_r8,0.0056_r8, &
! Subprogram not used                0.0057_r8,0.0144_r8,0.0135_r8,0.0131_r8,0.0130_r8,0.0132_r8,0.0135_r8, &
! Subprogram not used                0.0140_r8,0.0145_r8,0.0150_r8,0.0155_r8,0.0161_r8,0.0166_r8,0.0170_r8, &
! Subprogram not used                0.0174_r8,0.0177_r8,0.0180_r8,0.0182_r8,0.0525_r8,0.0473_r8,0.0439_r8, &
! Subprogram not used                0.0417_r8,0.0405_r8,0.0399_r8,0.0399_r8,0.0402_r8,0.0407_r8,0.0414_r8, &
! Subprogram not used                0.0421_r8,0.0429_r8,0.0437_r8,0.0445_r8,0.0452_r8,0.0458_r8,0.0464_r8, &
! Subprogram not used                0.1135_r8,0.1005_r8,0.0913_r8,0.0848_r8,0.0802_r8,0.0772_r8,0.0752_r8, &
! Subprogram not used                0.0740_r8,0.0734_r8,0.0732_r8,0.0733_r8,0.0736_r8,0.0739_r8,0.0744_r8, &
! Subprogram not used                0.0749_r8,0.0753_r8,0.0757_r8,0.1971_r8,0.1730_r8,0.1553_r8,0.1422_r8, &
! Subprogram not used                0.1326_r8,0.1255_r8,0.1202_r8,0.1164_r8,0.1137_r8,0.1118_r8,0.1104_r8, &
! Subprogram not used                0.1095_r8,0.1088_r8,0.1084_r8,0.1081_r8,0.1079_r8,0.1077_r8,0.3029_r8, &
! Subprogram not used                0.2645_r8,0.2358_r8,0.2140_r8,0.1975_r8,0.1848_r8,0.1750_r8,0.1675_r8, &
! Subprogram not used                0.1617_r8,0.1572_r8,0.1537_r8,0.1509_r8,0.1487_r8,0.1468_r8,0.1453_r8, &
! Subprogram not used                0.1440_r8,0.1429_r8,0.4307_r8,0.3749_r8,0.3325_r8,0.3000_r8,0.2747_r8, &
! Subprogram not used                0.2550_r8,0.2395_r8,0.2272_r8,0.2174_r8,0.2095_r8,0.2030_r8,0.1978_r8, &
! Subprogram not used                0.1934_r8,0.1897_r8,0.1865_r8,0.1837_r8,0.1812_r8,0.5801_r8,0.5037_r8, &
! Subprogram not used                0.4452_r8,0.3998_r8,0.3642_r8,0.3360_r8,0.3134_r8,0.2953_r8,0.2805_r8, &
! Subprogram not used                0.2684_r8,0.2584_r8,0.2500_r8,0.2429_r8,0.2368_r8,0.2316_r8,0.2269_r8, &
! Subprogram not used                0.2228_r8,0.7503_r8,0.6505_r8,0.5734_r8,0.5131_r8,0.4655_r8,0.4273_r8, &
! Subprogram not used                0.3966_r8,0.3715_r8,0.3509_r8,0.3339_r8,0.3196_r8,0.3075_r8,0.2972_r8, &
! Subprogram not used                0.2882_r8,0.2805_r8,0.2736_r8,0.2674_r8,0.9408_r8,0.8148_r8,0.7168_r8, &
! Subprogram not used                0.6396_r8,0.5782_r8,0.5288_r8,0.4887_r8,0.4557_r8,0.4284_r8,0.4056_r8, &
! Subprogram not used                0.3864_r8,0.3700_r8,0.3560_r8,0.3437_r8,0.3331_r8,0.3236_r8,0.3151_r8, &
! Subprogram not used                1.1504_r8,0.9956_r8,0.8746_r8,0.7787_r8,0.7021_r8,0.6400_r8,0.5893_r8, &
! Subprogram not used                0.5475_r8,0.5126_r8,0.4834_r8,0.4586_r8,0.4374_r8,0.4191_r8,0.4032_r8, &
! Subprogram not used                0.3892_r8,0.3768_r8,0.3657_r8,1.2620_r8,1.0919_r8,0.9586_r8,0.8528_r8, &
! Subprogram not used                0.7679_r8,0.6991_r8,0.6427_r8,0.5961_r8,0.5572_r8,0.5245_r8,0.4967_r8, &
! Subprogram not used                0.4728_r8,0.4523_r8,0.4343_r8,0.4186_r8,0.4046_r8,0.3921_r8/
! Subprogram not used       data zz3/0.0002_r8,0.0002_r8,0.0002_r8,0.0002_r8,0.0002_r8,0.0002_r8,0.0001_r8, &
! Subprogram not used                0.0001_r8,0.0001_r8,0.0001_r8,0.0001_r8,0.0001_r8,0.0001_r8,0.0001_r8, &
! Subprogram not used                0.0001_r8,0.0001_r8,0.0001_r8,0.0004_r8,0.0004_r8,0.0003_r8,0.0003_r8, &
! Subprogram not used                0.0003_r8,0.0003_r8,0.0003_r8,0.0003_r8,0.0003_r8,0.0003_r8,0.0003_r8, &
! Subprogram not used                0.0003_r8,0.0003_r8,0.0003_r8,0.0003_r8,0.0002_r8,0.0002_r8,0.0006_r8, &
! Subprogram not used                0.0006_r8,0.0006_r8,0.0006_r8,0.0006_r8,0.0006_r8,0.0006_r8,0.0006_r8, &
! Subprogram not used                0.0006_r8,0.0005_r8,0.0005_r8,0.0005_r8,0.0005_r8,0.0005_r8,0.0005_r8, &
! Subprogram not used                0.0004_r8,0.0004_r8,0.0040_r8,0.0039_r8,0.0038_r8,0.0038_r8,0.0037_r8, &
! Subprogram not used                0.0037_r8,0.0036_r8,0.0035_r8,0.0034_r8,0.0033_r8,0.0032_r8,0.0031_r8, &
! Subprogram not used                0.0030_r8,0.0029_r8,0.0028_r8,0.0027_r8,0.0026_r8,0.0154_r8,0.0150_r8, &
! Subprogram not used                0.0147_r8,0.0145_r8,0.0143_r8,0.0141_r8,0.0139_r8,0.0136_r8,0.0133_r8, &
! Subprogram not used                0.0130_r8,0.0126_r8,0.0122_r8,0.0119_r8,0.0115_r8,0.0111_r8,0.0107_r8, &
! Subprogram not used                0.0103_r8,0.0543_r8,0.0528_r8,0.0518_r8,0.0511_r8,0.0504_r8,0.0498_r8, &
! Subprogram not used                0.0491_r8,0.0484_r8,0.0475_r8,0.0465_r8,0.0454_r8,0.0443_r8,0.0431_r8, &
! Subprogram not used                0.0419_r8,0.0406_r8,0.0393_r8,0.0380_r8,0.1702_r8,0.1617_r8,0.1560_r8, &
! Subprogram not used                0.1521_r8,0.1492_r8,0.1470_r8,0.1449_r8,0.1430_r8,0.1409_r8,0.1387_r8, &
! Subprogram not used                0.1363_r8,0.1338_r8,0.1310_r8,0.1281_r8,0.1251_r8,0.1219_r8,0.1187_r8, &
! Subprogram not used                0.3229_r8,0.3005_r8,0.2846_r8,0.2729_r8,0.2642_r8,0.2573_r8,0.2516_r8, &
! Subprogram not used                0.2467_r8,0.2421_r8,0.2377_r8,0.2334_r8,0.2290_r8,0.2246_r8,0.2200_r8, &
! Subprogram not used                0.2153_r8,0.2104_r8,0.2055_r8,0.5089_r8,0.4673_r8,0.4363_r8,0.4128_r8, &
! Subprogram not used                0.3946_r8,0.3801_r8,0.3681_r8,0.3579_r8,0.3489_r8,0.3407_r8,0.3330_r8, &
! Subprogram not used                0.3257_r8,0.3185_r8,0.3115_r8,0.3045_r8,0.2975_r8,0.2906_r8,0.7263_r8, &
! Subprogram not used                0.6605_r8,0.6104_r8,0.5716_r8,0.5410_r8,0.5162_r8,0.4957_r8,0.4783_r8, &
! Subprogram not used                0.4631_r8,0.4496_r8,0.4373_r8,0.4258_r8,0.4149_r8,0.4046_r8,0.3945_r8, &
! Subprogram not used                0.3848_r8,0.3752_r8,0.9737_r8,0.8791_r8,0.8060_r8,0.7486_r8,0.7028_r8, &
! Subprogram not used                0.6654_r8,0.6344_r8,0.6080_r8,0.5852_r8,0.5651_r8,0.5469_r8,0.5303_r8, &
! Subprogram not used                0.5149_r8,0.5005_r8,0.4867_r8,0.4735_r8,0.4609_r8,1.2496_r8,1.1217_r8, &
! Subprogram not used                1.0219_r8,0.9429_r8,0.8792_r8,0.8271_r8,0.7836_r8,0.7467_r8,0.7149_r8, &
! Subprogram not used                0.6870_r8,0.6621_r8,0.6395_r8,0.6188_r8,0.5995_r8,0.5815_r8,0.5645_r8, &
! Subprogram not used                0.5482_r8,1.5516_r8,1.3866_r8,1.2568_r8,1.1532_r8,1.0693_r8,1.0004_r8, &
! Subprogram not used                0.9428_r8,0.8939_r8,0.8518_r8,0.8150_r8,0.7824_r8,0.7531_r8,0.7263_r8, &
! Subprogram not used                0.7018_r8,0.6789_r8,0.6575_r8,0.6373_r8,1.8769_r8,1.6714_r8,1.5088_r8, &
! Subprogram not used                1.3782_r8,1.2720_r8,1.1844_r8,1.1111_r8,1.0489_r8,0.9954_r8,0.9488_r8, &
! Subprogram not used                0.9076_r8,0.8707_r8,0.8374_r8,0.8069_r8,0.7788_r8,0.7527_r8,0.7282_r8, &
! Subprogram not used                2.2222_r8,1.9737_r8,1.7758_r8,1.6162_r8,1.4859_r8,1.3780_r8,1.2876_r8, &
! Subprogram not used                1.2109_r8,1.1450_r8,1.0876_r8,1.0371_r8,0.9921_r8,0.9516_r8,0.9147_r8, &
! Subprogram not used                0.8809_r8,0.8497_r8,0.8206_r8,2.4013_r8,2.1306_r8,1.9144_r8,1.7396_r8, &
! Subprogram not used                1.5966_r8,1.4781_r8,1.3787_r8,1.2944_r8,1.2219_r8,1.1588_r8,1.1034_r8, &
! Subprogram not used                1.0541_r8,1.0098_r8,0.9696_r8,0.9328_r8,0.8989_r8,0.8673_r8/ 
! Subprogram not used       data zz4/0.0027_r8,0.0023_r8,0.0019_r8,0.0016_r8,0.0014_r8,0.0012_r8,0.0011_r8, &
! Subprogram not used                0.0009_r8,0.0008_r8,0.0007_r8,0.0006_r8,0.0006_r8,0.0005_r8,0.0005_r8, &
! Subprogram not used                0.0004_r8,0.0004_r8,0.0003_r8,0.0059_r8,0.0050_r8,0.0043_r8,0.0037_r8, &
! Subprogram not used                0.0032_r8,0.0027_r8,0.0024_r8,0.0021_r8,0.0018_r8,0.0016_r8,0.0014_r8, &
! Subprogram not used                0.0013_r8,0.0011_r8,0.0010_r8,0.0009_r8,0.0008_r8,0.0007_r8,0.0105_r8, &
! Subprogram not used                0.0089_r8,0.0076_r8,0.0065_r8,0.0056_r8,0.0049_r8,0.0042_r8,0.0037_r8, &
! Subprogram not used                0.0033_r8,0.0029_r8,0.0025_r8,0.0023_r8,0.0020_r8,0.0018_r8,0.0016_r8, &
! Subprogram not used                0.0015_r8,0.0013_r8,0.0649_r8,0.0550_r8,0.0468_r8,0.0402_r8,0.0347_r8, &
! Subprogram not used                0.0301_r8,0.0262_r8,0.0230_r8,0.0202_r8,0.0178_r8,0.0158_r8,0.0141_r8, &
! Subprogram not used                0.0126_r8,0.0113_r8,0.0101_r8,0.0091_r8,0.0082_r8,0.2465_r8,0.2097_r8, &
! Subprogram not used                0.1794_r8,0.1544_r8,0.1336_r8,0.1162_r8,0.1015_r8,0.0891_r8,0.0785_r8, &
! Subprogram not used                0.0695_r8,0.0617_r8,0.0550_r8,0.0491_r8,0.0441_r8,0.0396_r8,0.0358_r8, &
! Subprogram not used                0.0323_r8,0.8274_r8,0.7127_r8,0.6168_r8,0.5364_r8,0.4684_r8,0.4108_r8, &
! Subprogram not used                0.3616_r8,0.3196_r8,0.2834_r8,0.2522_r8,0.2251_r8,0.2015_r8,0.1810_r8, &
! Subprogram not used                0.1630_r8,0.1471_r8,0.1332_r8,0.1208_r8,2.1322_r8,1.8745_r8,1.6542_r8, &
! Subprogram not used                1.4649_r8,1.3016_r8,1.1601_r8,1.0371_r8,0.9298_r8,0.8358_r8,0.7532_r8, &
! Subprogram not used                0.6805_r8,0.6163_r8,0.5593_r8,0.5088_r8,0.4638_r8,0.4236_r8,0.3876_r8, &
! Subprogram not used                3.2645_r8,2.8943_r8,2.5762_r8,2.3012_r8,2.0625_r8,1.8542_r8,1.6717_r8, &
! Subprogram not used                1.5114_r8,1.3698_r8,1.2446_r8,1.1334_r8,1.0344_r8,0.9460_r8,0.8668_r8, &
! Subprogram not used                0.7958_r8,0.7319_r8,0.6743_r8,4.2737_r8,3.8017_r8,3.3956_r8,3.0442_r8, &
! Subprogram not used                2.7385_r8,2.4714_r8,2.2370_r8,2.0305_r8,1.8478_r8,1.6858_r8,1.5415_r8, &
! Subprogram not used                1.4126_r8,1.2972_r8,1.1936_r8,1.1004_r8,1.0162_r8,0.9400_r8,5.2092_r8, &
! Subprogram not used                4.6423_r8,4.1542_r8,3.7314_r8,3.3633_r8,3.0414_r8,2.7585_r8,2.5090_r8, &
! Subprogram not used                2.2882_r8,2.0920_r8,1.9171_r8,1.7607_r8,1.6205_r8,1.4945_r8,1.3809_r8, &
! Subprogram not used                1.2782_r8,1.1851_r8,6.0871_r8,5.4326_r8,4.8683_r8,4.3790_r8,3.9525_r8, &
! Subprogram not used                3.5791_r8,3.2507_r8,2.9608_r8,2.7039_r8,2.4754_r8,2.2716_r8,2.0892_r8, &
! Subprogram not used                1.9256_r8,1.7782_r8,1.6453_r8,1.5251_r8,1.4161_r8,6.9131_r8,6.1782_r8, &
! Subprogram not used                5.5437_r8,4.9928_r8,4.5120_r8,4.0906_r8,3.7196_r8,3.3917_r8,3.1008_r8, &
! Subprogram not used                2.8419_r8,2.6106_r8,2.4035_r8,2.2175_r8,2.0500_r8,1.8987_r8,1.7617_r8, &
! Subprogram not used                1.6374_r8,7.6903_r8,6.8822_r8,6.1833_r8,5.5755_r8,5.0445_r8,4.5784_r8, &
! Subprogram not used                4.1676_r8,3.8041_r8,3.4813_r8,3.1937_r8,2.9366_r8,2.7061_r8,2.4989_r8, &
! Subprogram not used                2.3121_r8,2.1432_r8,1.9903_r8,1.8514_r8,8.4213_r8,7.5468_r8,6.7890_r8, &
! Subprogram not used                6.1290_r8,5.5515_r8,5.0440_r8,4.5961_r8,4.1994_r8,3.8467_r8,3.5321_r8, &
! Subprogram not used                3.2507_r8,2.9981_r8,2.7708_r8,2.5657_r8,2.3801_r8,2.2119_r8,2.0590_r8, &
! Subprogram not used                9.1086_r8,8.1740_r8,7.3627_r8,6.6548_r8,6.0345_r8,5.4886_r8,5.0062_r8, &
! Subprogram not used                4.5785_r8,4.1978_r8,3.8579_r8,3.5535_r8,3.2801_r8,3.0338_r8,2.8114_r8, &
! Subprogram not used                2.6100_r8,2.4272_r8,2.2610_r8,9.4365_r8,8.4743_r8,7.6380_r8,6.9077_r8, &
! Subprogram not used                6.2673_r8,5.7033_r8,5.2047_r8,4.7622_r8,4.3683_r8,4.0163_r8,3.7009_r8, &
! Subprogram not used                3.4175_r8,3.1621_r8,2.9314_r8,2.7224_r8,2.5326_r8,2.3600_r8/ 
! Subprogram not used                     
! Subprogram not used       if(p.le.5) then
! Subprogram not used         do if = 1, pnf_msu
! Subprogram not used           abs(if)=0.0_r8
! Subprogram not used         end do
! Subprogram not used         return
! Subprogram not used       endif
! Subprogram not used       
! Subprogram not used       call locate(xx,n,t,i)
! Subprogram not used       call locate(yy,m,p,j)
! Subprogram not used             
! Subprogram not used       t1=(t-xx(i))/(xx(i+1)-xx(i))
! Subprogram not used       t2=(p-yy(j))/(yy(j+1)-yy(j))
! Subprogram not used       abs(1)=(1-t1)*(1-t2)*zz1(i,j)+t1*(1-t2)*zz1(i+1,j)+ &
! Subprogram not used                       t1*t2*zz1(i+1,j+1)+(1-t1)*t2*zz1(i,j+1)
! Subprogram not used       abs(2)=(1-t1)*(1-t2)*zz2(i,j)+t1*(1-t2)*zz2(i+1,j)+ &
! Subprogram not used                       t1*t2*zz2(i+1,j+1)+(1-t1)*t2*zz2(i,j+1)
! Subprogram not used       abs(3)=(1-t1)*(1-t2)*zz3(i,j)+t1*(1-t2)*zz3(i+1,j)+ &
! Subprogram not used                       t1*t2*zz3(i+1,j+1)+(1-t1)*t2*zz3(i,j+1)
! Subprogram not used       abs(4)=(1-t1)*(1-t2)*zz4(i,j)+t1*(1-t2)*zz4(i+1,j)+ &
! Subprogram not used                       t1*t2*zz4(i+1,j+1)+(1-t1)*t2*zz4(i,j+1)
! Subprogram not used            
! Subprogram not used     end subroutine lookup
!
!===================================================================
!
! Subprogram not used     subroutine locate(xx,n,x,j)
! Subprogram not used 
! Subprogram not used       integer n, j
! Subprogram not used       real(r8) xx(n), x
! Subprogram not used !
! Subprogram not used       integer jl, ju, jm
! Subprogram not used !
! Subprogram not used       jl=0
! Subprogram not used       ju=n+1
! Subprogram not used       do while (ju-jl.gt.1)
! Subprogram not used         jm=(ju+jl)/2
! Subprogram not used         if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
! Subprogram not used           jl=jm
! Subprogram not used         else
! Subprogram not used           ju=jm
! Subprogram not used         endif
! Subprogram not used       end do
! Subprogram not used       j=jl
! Subprogram not used 
! Subprogram not used     end subroutine locate

    subroutine hirsbt_init ()

! Initialize values used for the HIRS brightness temperature calculation.
!
    use time_manager, only: get_step_size

!------------------------------Local variables--------------------------
    integer :: dtime      ! integer timestep size

!*******************************************************************************************
!     Set constants to values used in the model
!     Leave as they are for now.  Later if this gets committed back to the main 
!     development trunk we should use the values used in the rest of the radiation code.
! mji
! This step is done in hirsrtm.f90 to use values consistent with the radiation code.
!
!     GRAV = GRAVIT
!     R = RAIR
!     RCO2 = CO2VMR*RMWCO2
!
!*******************************************************************************************

    msuname(1)  = 'MSU_1   '
    msuname(2)  = 'MSU_2   '
    msuname(3)  = 'MSU_3   '
    msuname(4)  = 'MSU_4   '
    hirsname(1) = 'HIRS_2  '
    hirsname(2) = 'HIRS_4  '
    hirsname(3) = 'HIRS_6  '
    hirsname(4) = 'HIRS_8  '
    hirsname(5) = 'HIRS_10 '
    hirsname(6) = 'HIRS_11 '
    hirsname(7) = 'HIRS_12 '

    dtime  = get_step_size()

! These should be namelist variables; 
! Set flag to do HIRS brightness temperature calculation 
    dohirs = .true.
! Set frequency of HIRS calculation
! ihirsfq is in timesteps if positive, or hours if negative; 6 hours is recommended
    ihirsfq = -6

! Convert ihirsfq from hours to timesteps if necessary
    if (ihirsfq < 0) ihirsfq = nint((-ihirsfq*3600._r8)/dtime)

    end subroutine hirsbt_init

    end module hirsbt

