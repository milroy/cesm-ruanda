!=======================================================================
!
!BOP
!
! !MODULE: ice_meltpond - Meltpond parameterization
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!  SVN:$$
!
! authors David A. Bailey (NCAR)
!         Marika M. Holland (NCAR)
!
! !INTERFACE:
!
      module ice_meltpond
!
! !USES:
!
      use ice_kinds_mod
      use ice_constants
      use ice_fileunits
      use ice_read_write
      use ice_restart, only: lenstr, restart_dir, restart_file, &
                             pointer_file, runtype
      use ice_communicate, only: my_task, master_task
      use ice_exit, only: abort_ice
!
!EOP
!
      implicit none

      logical (kind=log_kind) :: &
         restart_pond    ! if .true., read meltponds restart file

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: init_meltponds
!
! !DESCRIPTION:
!
!  Initialize melt ponds.
! 
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine init_meltponds
!
! !USES:
!
      use ice_state, only: filename_volpn
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!

      if (trim(filename_volpn) /= 'none') restart_pond = .true.

      if (restart_pond) then
         if (trim(runtype) == 'continue') then
            call read_restart_pond
         else
            call read_restart_pond(filename_volpn)
         endif
      endif

      end subroutine init_meltponds

!=======================================================================
!BOP
!
! !ROUTINE: 
!
! !INTERFACE:
!
      subroutine compute_ponds(nx_block,ny_block,          &
                               ilo, ihi, jlo, jhi,         &
                               meltt, melts,  frain,       &
                               aicen, vicen,  vsnon,       &
                               trcrn, apondn, hpondn)
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! same as module
!
! !USES:
!
      use ice_state, only: nt_Tsfc, nt_volpn
      use ice_calendar, only: dt
      use ice_domain_size, only: max_ntrcr

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ilo,ihi,jlo,jhi       ! beginning and end of physical domain

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(in) :: &
         meltt, &
         melts, &
         frain, &
         aicen, &
         vicen, &
         vsnon

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(inout) :: &
         apondn, &
         hpondn

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_ntrcr), &
         intent(inout) :: &
         trcrn

      real (kind=dbl_kind), dimension(nx_block,ny_block) :: &
         volpn, &
         Tsfcn

      integer (kind=int_kind), dimension (nx_block*ny_block) :: &
         indxi, indxj     ! compressed indices for cells with ice melting

      integer (kind=int_kind) :: i,j,ij,icells

      real (kind=dbl_kind) :: hi,hs,dTs,rfrac,asnow

      real (kind=dbl_kind), parameter :: &
         hi_min = p1

      Tsfcn(:,:) = trcrn(:,:,nt_Tsfc)
      volpn(:,:) = trcrn(:,:,nt_volpn)

      !-----------------------------------------------------------------
      ! Identify grid cells where ice can melt.
      !-----------------------------------------------------------------

      icells = 0
      do j = jlo, jhi
      do i = ilo, ihi
         if (aicen(i,j) > puny) then
            icells = icells + 1
            indxi(icells) = i
            indxj(icells) = j
         endif
      enddo                     ! i
      enddo                     ! j

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         hi = vicen(i,j)/aicen(i,j)
         hs = vsnon(i,j)/aicen(i,j)
         dTs = Timelt - Tsfcn(i,j)

! Make the runoff fraction proportional to ice area.
         rfrac = 0.85_dbl_kind - 0.7_dbl_kind*aicen(i,j)

         volpn(i,j) = volpn(i,j) + (c1-rfrac) &
            * (meltt(i,j)*(rhoi/rhofresh) + melts(i,j)*(rhos/rhofresh) &
            + frain(i,j)*dt/rhofresh)

!        Use exponential decrease in pond volume DAB
         if (Tsfcn(i,j) .lt. Timelt-c2) then
            volpn(i,j) = volpn(i,j) &
               *exp(p01*(dTs-c2)/(Timelt-c2))
         endif

         volpn(i,j) = max(volpn(i,j), c0)

         apondn(i,j) = min ( sqrt ( volpn(i,j) * 1.25 ), c1)
         hpondn(i,j) = 0.8 * apondn(i,j)

!        Limit pond depth to 90% of ice thickness
         if (hpondn(i,j) .gt. 0.9*hi) then
            hpondn(i,j) = 0.9*hi
            volpn(i,j) = hpondn(i,j)*apondn(i,j)
!           apondn(i,j) = min(volpn(i,j)/hpondn(i,j), 1.)
         endif

!        If we are freezing, do not change albedo
!        if (Tsfcn(i,j) .lt. Timelt-0.15) apondn(i,j) = c0

         ! fractional area of snow cover
         ! use linear function as in delta-Eddington
         if (hs >= hsmin) then
            asnow = min( hs/hs0, c1 )
         else
            asnow = c0
         endif

         ! fractional area of snow cover
         ! use snowpatch as in the CCSM3 shortwave
!        if (hs > puny) then
!           asnow = hs / (hs + snowpatch)
!        else
!           asnow = c0
!        endif

!        If we have snow, reduce ponds
         if ( hs .gt. puny ) then
            apondn(i,j) = (c1-asnow)*apondn(i,j)
            hpondn(i,j) = 0.8_dbl_kind * apondn(i,j)
         endif

!        remove ponds if ice becomes very thin
         if (hi .lt. hi_min) then
            apondn(i,j) = c0
            hpondn(i,j) = c0
            volpn(i,j)  = c0
         endif

         trcrn(i,j,nt_volpn) = volpn(i,j)

      enddo

      end subroutine compute_ponds

!=======================================================================
!
!BOP
!
! !IROUTINE: write_restart_pond - dumps all fields required for restart
!
! !INTERFACE:
!
! Subprogram not used       subroutine write_restart_pond(filename_spec)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Dumps all values needed for restarting
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! authors Elizabeth C. Hunke, LANL
! Subprogram not used !         David A. Bailey, NCAR
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use ice_domain_size
! Subprogram not used       use ice_calendar, only: sec, month, mday, nyr, istep1, &
! Subprogram not used                               time, time_forc, idate, year_init
! Subprogram not used       use ice_state
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       character(len=char_len_long), intent(in), optional :: filename_spec
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used           i, j, k, n, it, iblk, & ! counting indices
! Subprogram not used           iyear, imonth, iday     ! year, month, day
! Subprogram not used 
! Subprogram not used       character(len=char_len_long) :: filename
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind) :: diag
! Subprogram not used 
! Subprogram not used       ! construct path/file
! Subprogram not used       if (present(filename_spec)) then
! Subprogram not used          filename = trim(filename_spec)
! Subprogram not used       else
! Subprogram not used          iyear = nyr + year_init - 1
! Subprogram not used          imonth = month
! Subprogram not used          iday = mday
! Subprogram not used          
! Subprogram not used          write(filename,'(a,a,a,i4.4,a,i2.2,a,i2.2,a,i5.5)') &
! Subprogram not used               restart_dir(1:lenstr(restart_dir)), &
! Subprogram not used               restart_file(1:lenstr(restart_file)),'.volpn.', &
! Subprogram not used               iyear,'-',month,'-',mday,'-',sec
! Subprogram not used       end if
! Subprogram not used          
! Subprogram not used       ! begin writing restart data
! Subprogram not used       call ice_open(nu_dump_pond,filename,0)
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used         write(nu_dump_pond) istep1,time,time_forc
! Subprogram not used         write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       diag = .true.
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       do n = 1, ncat
! Subprogram not used          call ice_write(nu_dump_pond,0,trcrn(:,:,nt_volpn,n,:),'ruf8',diag)
! Subprogram not used          call ice_write(nu_dump_pond,0,apondn(:,:,n,:),'ruf8',diag)
! Subprogram not used          call ice_write(nu_dump_pond,0,hpondn(:,:,n,:),'ruf8',diag)
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) close(nu_dump_pond)
! Subprogram not used 
! Subprogram not used       end subroutine write_restart_pond

!=======================================================================
!BOP
!
! !IROUTINE: read_restart_pond - reads all fields required for restart
!
! !INTERFACE:
!
! Subprogram not used       subroutine read_restart_pond(filename_spec)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Reads all values needed for a meltpond volume restart
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! authors Elizabeth C. Hunke, LANL
! Subprogram not used !         David A. Bailey, NCAR
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use ice_domain_size
! Subprogram not used       use ice_calendar, only: sec, month, mday, nyr, istep1, &
! Subprogram not used                               time, time_forc, idate, year_init
! Subprogram not used       use ice_state
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       character(len=char_len_long), intent(in), optional :: filename_spec
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used           i, j, k, n, it, iblk, & ! counting indices
! Subprogram not used           iyear, imonth, iday     ! year, month, day
! Subprogram not used 
! Subprogram not used       character(len=char_len_long) :: &
! Subprogram not used          filename, filename0, string1, string2
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind) :: &
! Subprogram not used          diag
! Subprogram not used 
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used          ! reconstruct path/file
! Subprogram not used          if (present(filename_spec)) then
! Subprogram not used             filename = filename_spec
! Subprogram not used          else
! Subprogram not used             open(nu_rst_pointer,file=pointer_file)
! Subprogram not used             read(nu_rst_pointer,'(a)') filename0
! Subprogram not used             filename = trim(filename0)
! Subprogram not used             close(nu_rst_pointer)
! Subprogram not used 
! Subprogram not used             n = index(filename0,trim(restart_file))
! Subprogram not used             if (n == 0) call abort_ice('volpn restart: filename discrepancy')
! Subprogram not used             string1 = trim(filename0(1:n-1))
! Subprogram not used             string2 = trim(filename0(n+lenstr(restart_file):lenstr(filename0)))
! Subprogram not used             write(filename,'(a,a,a,a)') &
! Subprogram not used                string1(1:lenstr(string1)), &
! Subprogram not used                restart_file(1:lenstr(restart_file)),'.volpn', &
! Subprogram not used                string2(1:lenstr(string2))
! Subprogram not used          endif
! Subprogram not used       endif ! master_task
! Subprogram not used 
! Subprogram not used       call ice_open(nu_restart_pond,filename,0)
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used         read(nu_restart_pond) istep1,time,time_forc
! Subprogram not used         write(nu_diag,*) 'Reading ',filename(1:lenstr(filename))
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       diag = .true.
! Subprogram not used 
! Subprogram not used       !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       do n = 1, ncat
! Subprogram not used          call ice_read(nu_restart_pond,0,trcrn(:,:,nt_volpn,n,:),'ruf8',diag)
! Subprogram not used          call ice_read(nu_restart_pond,0,apondn(:,:,n,:),'ruf8',diag)
! Subprogram not used          call ice_read(nu_restart_pond,0,hpondn(:,:,n,:),'ruf8',diag)
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) close(nu_restart_pond)
! Subprogram not used 
! Subprogram not used       end subroutine read_restart_pond

!=======================================================================

      end module ice_meltpond

!=======================================================================
