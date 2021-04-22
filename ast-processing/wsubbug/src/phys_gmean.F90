module phys_gmean
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Perform mixed layer global calculations for energy conservation checks.
!
! Methods: 
! Reproducible (nonscalable): 
!    Gather to a master processor who does all the work.
! Reproducible (scalable): 
!    Convert to fixed point (integer representation) to enable
!    reproducibility when using MPI collectives. Results compared with
!    a nonreproducible (but scalable) algorithm using floating point
!    and MPI_Allreduce to verify the results are good enough.
!
! Author: Byron Boville from SOM code by Jim Rosinski/Bruce Briegleb
! Modified: P. Worley to aggregate calculations (4/04)
! Modified: J. White/P. Worley to introduce scalable algorithms;
!           B. Eaton to remove dycore-specific dependencies and to 
!           introduce gmean_mass (10/07)
! Modified: P. Worley to replace in-place implementation with call
!           to repro_sum.
! 
!-----------------------------------------------------------------------
   use shr_kind_mod,  only: r8 => shr_kind_r8
   use physconst,     only: pi
   use spmd_utils,    only: masterproc
   use ppgrid,        only: pcols, begchunk, endchunk
   use shr_reprosum_mod, only: shr_reprosum_calc, shr_reprosum_tolExceeded, &
                            shr_reprosum_reldiffmax, shr_reprosum_recompute

   use mpishorthand

   use perf_mod
   use cam_logfile,   only: iulog

   implicit none
   private
   save

   public :: &
      gmean,       &! compute global mean of 2D fields on physics decomposition
      gmean_mass    ! compute global mean mass of constituent fields on physics decomposition

   interface gmean
      module procedure gmean_arr
      module procedure gmean_scl
   endinterface

   CONTAINS

!
!========================================================================
!

   subroutine gmean_arr (arr, arr_gmean, nflds)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute the global mean of each field in "arr" in the physics 
! chunked decomposition
! 
!-----------------------------------------------------------------------
!
! Arguments
!
      integer, intent(in)  :: nflds                 ! number of fields
      real(r8), intent(in) :: arr(pcols,begchunk:endchunk,nflds) 
                                                    ! Input array, chunked
      real(r8), intent(out):: arr_gmean(nflds)      ! global means
!
! Local workspace
!
      real(r8) :: rel_diff(2,nflds)                 ! relative differences between 
                                                    !  'fast' reproducible and 
                                                    !  nonreproducible means
      integer  :: ifld                              ! field index
      logical  :: write_warning
!
!-----------------------------------------------------------------------
!
      call t_startf ('gmean_fixed_repro')
      call gmean_fixed_repro(arr, arr_gmean, rel_diff, nflds)
      call t_stopf ('gmean_fixed_repro')

      ! check that "fast" reproducible sum is accurate enough. If not, calculate
      ! using old method
      write_warning = masterproc
      if ( shr_reprosum_tolExceeded('gmean', nflds, write_warning, &
                                  iulog, rel_diff) ) then
         if ( shr_reprosum_recompute ) then
            do ifld=1,nflds
               if ( rel_diff(1,ifld) > shr_reprosum_reldiffmax ) then
                  call t_startf ('gmean_float_repro')
                  call gmean_float_repro(arr(:,:,ifld), arr_gmean(ifld), 1)
                  call t_stopf ('gmean_float_repro')
               endif
            enddo
         endif
      endif

      return
   end subroutine gmean_arr 

!
!========================================================================
!

! Subprogram not used    subroutine gmean_scl (arr, gmean)
! Subprogram not used       use phys_grid, only : get_ncols_p
! Subprogram not used 
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: 
! Subprogram not used ! Compute the global mean of each field in "arr" in the physics 
! Subprogram not used ! chunked decomposition
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used ! Arguments
! Subprogram not used !
! Subprogram not used       real(r8), intent(in) :: arr(pcols,begchunk:endchunk) 
! Subprogram not used                                                     ! Input array, chunked
! Subprogram not used       real(r8), intent(out):: gmean      ! global means
! Subprogram not used !
! Subprogram not used ! Local workspace
! Subprogram not used !
! Subprogram not used       integer, parameter :: nflds = 1
! Subprogram not used       real(r8) :: gmean_array(nflds)
! Subprogram not used       real(r8) :: array(pcols,begchunk:endchunk,nflds) 
! Subprogram not used       integer  :: ncols, lchnk
! Subprogram not used 
! Subprogram not used       do lchnk=begchunk,endchunk
! Subprogram not used          ncols = get_ncols_p(lchnk)
! Subprogram not used          array(:ncols,lchnk,1) = arr(:ncols,lchnk)
! Subprogram not used       enddo
! Subprogram not used       call gmean_arr(array,gmean_array,nflds)
! Subprogram not used       gmean = gmean_array(1)
! Subprogram not used 
! Subprogram not used    end subroutine gmean_scl 

!
!========================================================================
!

! Subprogram not used    subroutine gmean_mass(title, state)
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used ! Purpose:
! Subprogram not used ! Computes global mean mass, max and min mmr, of constituents on the 
! Subprogram not used ! physics decomposition. Prints diagnostics to log file.
! Subprogram not used !
! Subprogram not used ! Author: B. Eaton (based on gavglook)
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used       use ppgrid,         only: pver
! Subprogram not used       use physconst,      only: gravit
! Subprogram not used       use phys_grid,      only: get_ncols_p
! Subprogram not used       use physics_types,  only: physics_state
! Subprogram not used       use constituents,   only: pcnst, cnst_name
! Subprogram not used !
! Subprogram not used ! Arguments
! Subprogram not used !
! Subprogram not used       character(len=*),    intent(in) :: title    ! location of this call
! Subprogram not used       type(physics_state), intent(in) :: state(begchunk:endchunk)
! Subprogram not used !
! Subprogram not used ! Local workspace
! Subprogram not used !
! Subprogram not used       character(len=*), parameter :: sub_name='gmean_mass: '
! Subprogram not used 
! Subprogram not used       integer :: c, i, k, m
! Subprogram not used       integer :: ierr
! Subprogram not used       integer :: ncols
! Subprogram not used 
! Subprogram not used       real(r8), pointer :: mass_wet(:,:,:) ! constituent masses assuming moist mmr
! Subprogram not used       real(r8), pointer :: mass_dry(:,:,:) ! constituent masses assuming dry mmr
! Subprogram not used       real(r8) :: mass_wet_mean(pcnst)     ! global mean constituent masses assuming moist mmr
! Subprogram not used       real(r8) :: mass_dry_mean(pcnst)     ! global mean constituent masses assuming dry mmr
! Subprogram not used       real(r8) :: mmr_max(pcnst)           ! maximum constituent mmr in this process
! Subprogram not used       real(r8) :: mmr_min(pcnst)           ! minimum constituent mmr in this process
! Subprogram not used       real(r8) :: mmr_max_glob(pcnst)      ! global maximum constituent mmr
! Subprogram not used       real(r8) :: mmr_min_glob(pcnst)      ! global minimum constituent mmr
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used       allocate(mass_wet(pcols,begchunk:endchunk,pcnst), stat=ierr)
! Subprogram not used       if (ierr /= 0) write(iulog,*) sub_name // 'FAIL to allocate mass_wet'
! Subprogram not used 
! Subprogram not used       allocate(mass_dry(pcols,begchunk:endchunk,pcnst), stat=ierr)
! Subprogram not used       if (ierr /= 0) write(iulog,*) sub_name // 'FAIL to allocate mass_wet'
! Subprogram not used 
! Subprogram not used       mmr_max(:) = -1.e36_r8
! Subprogram not used       mmr_min(:) =  1.e36_r8
! Subprogram not used       do m = 1, pcnst
! Subprogram not used          do c = begchunk, endchunk
! Subprogram not used             ncols = get_ncols_p(c)
! Subprogram not used             do i = 1, ncols
! Subprogram not used 
! Subprogram not used                ! Compute column masses assuming both dry and wet mixing ratios
! Subprogram not used 
! Subprogram not used                mass_wet(i,c,m) = 0.0_r8
! Subprogram not used                do k = 1, pver
! Subprogram not used                   mass_wet(i,c,m) = mass_wet(i,c,m) + &
! Subprogram not used                                     state(c)%pdel(i,k)*state(c)%q(i,k,m)
! Subprogram not used                   mmr_max(m) = max(mmr_max(m), state(c)%q(i,k,m))
! Subprogram not used                   mmr_min(m) = min(mmr_min(m), state(c)%q(i,k,m))
! Subprogram not used                end do
! Subprogram not used                mass_wet(i,c,m) = mass_wet(i,c,m)/gravit
! Subprogram not used 
! Subprogram not used                mass_dry(i,c,m) = 0.0_r8
! Subprogram not used                do k = 1, pver
! Subprogram not used                   mass_dry(i,c,m) = mass_dry(i,c,m) + &
! Subprogram not used                                     state(c)%pdeldry(i,k)*state(c)%q(i,k,m)
! Subprogram not used                end do
! Subprogram not used                mass_dry(i,c,m) = mass_dry(i,c,m)/gravit
! Subprogram not used 
! Subprogram not used             end do
! Subprogram not used          end do
! Subprogram not used       end do
! Subprogram not used 
! Subprogram not used       ! compute global mean mass
! Subprogram not used       call gmean(mass_wet, mass_wet_mean, pcnst)
! Subprogram not used       call gmean(mass_dry, mass_dry_mean, pcnst)
! Subprogram not used 
! Subprogram not used       ! global min/max mmr
! Subprogram not used 
! Subprogram not used       call mpi_reduce(mmr_max, mmr_max_glob, pcnst, mpir8, MPI_MAX, 0, mpicom, ierr)
! Subprogram not used       call mpi_reduce(mmr_min, mmr_min_glob, pcnst, mpir8, MPI_MIN, 0, mpicom, ierr)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used       ! report to log file
! Subprogram not used       if (masterproc) then
! Subprogram not used 
! Subprogram not used          do m = 1, pcnst
! Subprogram not used                write (6,66) trim(title)//' m=',m, &
! Subprogram not used                   'name='//trim(cnst_name(m))//' gavg dry, wet, min, max ', &
! Subprogram not used                   mass_dry_mean(m), mass_wet_mean(m), mmr_min_glob(m), mmr_max_glob(m)
! Subprogram not used 66             format (a24,i2,a36,1p,4e25.13)
! Subprogram not used          end do
! Subprogram not used 
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       deallocate(mass_wet)
! Subprogram not used       deallocate(mass_dry)
! Subprogram not used 
! Subprogram not used    end subroutine gmean_mass

!
!========================================================================
!

! Subprogram not used    subroutine gmean_float_repro (arr, arr_gmean, nflds)
! Subprogram not used !----------------------------------------------------------------------- 
! Subprogram not used ! 
! Subprogram not used ! Purpose: 
! Subprogram not used ! Compute the global mean of each field in "arr" in the physics 
! Subprogram not used ! chunked decomposition - all work is done on the masterproc to avoid
! Subprogram not used ! order of operations differences and assure bfb reproducibility.
! Subprogram not used ! 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       use rgrid,        only: nlon
! Subprogram not used       use dycore,       only: dycore_is
! Subprogram not used       use phys_grid,    only: gather_chunk_to_field
! Subprogram not used       use dyn_grid,     only: get_horiz_grid_dim_d, get_horiz_grid_d, get_dyn_grid_parm_real1d
! Subprogram not used !
! Subprogram not used ! Arguments
! Subprogram not used !
! Subprogram not used       integer, intent(in)  :: nflds               ! number of fields
! Subprogram not used       real(r8), intent(in) :: &
! Subprogram not used          arr(pcols,begchunk:endchunk,nflds)       ! Input array, chunked
! Subprogram not used       real(r8), intent(out):: arr_gmean(nflds)    ! global means
! Subprogram not used !
! Subprogram not used ! Local workspace
! Subprogram not used !
! Subprogram not used       real(r8), pointer :: w(:)
! Subprogram not used       real(r8) :: zmean                       ! zonal mean value
! Subprogram not used       real(r8) :: tmean                       ! temp global mean value
! Subprogram not used       integer :: i, j, ifld, n                ! longitude, latitude, field, 
! Subprogram not used                                               !  and global column indices
! Subprogram not used       integer :: hdim1, hdim2                 ! dimensions of rectangular horizontal 
! Subprogram not used                                               !  grid data structure, If 1D data 
! Subprogram not used                                               !  structure, then hdim2_d == 1.
! Subprogram not used       integer :: ngcols                       ! global column count (all)
! Subprogram not used 
! Subprogram not used       ! rectangular version of arr
! Subprogram not used       real(r8), allocatable :: arr_field(:,:,:)   
! Subprogram not used 
! Subprogram not used       ! column integration weight (from dynamics)
! Subprogram not used       real(r8), dimension(:), allocatable :: wght_d 
! Subprogram not used 
! Subprogram not used !
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used       call get_horiz_grid_dim_d(hdim1, hdim2)
! Subprogram not used       allocate(arr_field(hdim1,hdim2,nflds))
! Subprogram not used 
! Subprogram not used       arr_field(:,:,:) = 0.0_r8
! Subprogram not used       call gather_chunk_to_field (1, 1, nflds, hdim1, arr, arr_field)
! Subprogram not used 
! Subprogram not used       if (masterproc) then
! Subprogram not used 
! Subprogram not used          if (dycore_is('UNSTRUCTURED')) then
! Subprogram not used 
! Subprogram not used             ngcols = hdim1*hdim2
! Subprogram not used             allocate ( wght_d(1:ngcols) )
! Subprogram not used 
! Subprogram not used             wght_d = 0.0_r8
! Subprogram not used             call get_horiz_grid_d(ngcols, wght_d_out=wght_d)
! Subprogram not used 
! Subprogram not used             do ifld=1,nflds
! Subprogram not used                arr_gmean(ifld) = 0._r8
! Subprogram not used                do j=1,hdim2
! Subprogram not used                   do i=1,hdim1
! Subprogram not used                      n = (j-1)*hdim1 + i
! Subprogram not used                      arr_gmean(ifld) = arr_gmean(ifld) + &
! Subprogram not used                                        arr_field(i,j,ifld)*wght_d(n)
! Subprogram not used                   end do
! Subprogram not used                end do
! Subprogram not used                arr_gmean(ifld) = arr_gmean(ifld) / (4.0_r8 * pi)
! Subprogram not used             end do
! Subprogram not used 
! Subprogram not used             deallocate ( wght_d )
! Subprogram not used 
! Subprogram not used          else
! Subprogram not used             w => get_dyn_grid_parm_real1d('w')
! Subprogram not used             do ifld=1,nflds
! Subprogram not used                tmean = 0._r8
! Subprogram not used                do j=1,hdim2
! Subprogram not used                   zmean = 0._r8
! Subprogram not used                   do i=1,hdim1
! Subprogram not used                      zmean = zmean + arr_field(i,j,ifld)
! Subprogram not used                   end do
! Subprogram not used                   tmean = tmean + zmean * 0.5_r8*w(j)/nlon(j)
! Subprogram not used                end do
! Subprogram not used                arr_gmean(ifld) = tmean
! Subprogram not used             end do
! Subprogram not used 
! Subprogram not used          end if
! Subprogram not used 
! Subprogram not used       end if
! Subprogram not used 
! Subprogram not used 
! Subprogram not used       call mpibcast (arr_gmean, nflds, mpir8, 0, mpicom)
! Subprogram not used 
! Subprogram not used       deallocate(arr_field)
! Subprogram not used 
! Subprogram not used       return
! Subprogram not used 
! Subprogram not used    end subroutine gmean_float_repro

!
!========================================================================
!
   subroutine gmean_fixed_repro (arr, arr_gmean, rel_diff, nflds)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute the global mean of each field in "arr" in the physics 
! chunked decomposition with a reproducible yet scalable implementation
! based on a fixed-point algorithm.
!
!-----------------------------------------------------------------------
      use phys_grid, only    : get_ncols_p, get_wght_all_p, ngcols_p, &
                               get_nlcols_p
!
! Arguments
!
      integer, intent(in)  :: nflds             ! number of fields
      real(r8), intent(in) :: &
         arr(pcols,begchunk:endchunk,nflds)     ! Input array, chunked
      real(r8), intent(out):: arr_gmean(nflds)  ! global means
      real(r8), intent(out):: rel_diff(2,nflds) ! relative and absolute
                                                !  differences between 
                                                !  reproducible and nonreproducible
                                                !  means
!
! Local workspace
!
      integer :: lchnk, i, ifld                 ! chunk, column, field indices
      integer :: ncols                          ! number of columns in current chunk
      integer :: count                          ! summand count
      integer :: ierr                           ! MPI error return



      
      real(r8) :: wght(pcols)                   ! column for integration weights
      real(r8), allocatable :: xfld(:,:)        ! weighted summands
      integer :: nlcols
!
!-----------------------------------------------------------------------
!
      nlcols = get_nlcols_p()
      allocate(xfld(nlcols, nflds))

! pre-weight summands
      do ifld=1,nflds
         count = 0
         do lchnk=begchunk,endchunk
            ncols = get_ncols_p(lchnk)
            call get_wght_all_p(lchnk, ncols, wght)
            do i=1,ncols
               count = count + 1
               xfld(count,ifld) = arr(i,lchnk,ifld)*wght(i)
            end do
         end do
      end do

! call fixed-point algorithm
      call shr_reprosum_calc (xfld, arr_gmean, count, nlcols, nflds, &
                      gbl_count=ngcols_p, commid=mpicom, rel_diff=rel_diff) 

      deallocate(xfld)
! final normalization
      arr_gmean(:) = arr_gmean(:) / (4.0_r8 * pi)

      return

   end subroutine gmean_fixed_repro

end module phys_gmean
