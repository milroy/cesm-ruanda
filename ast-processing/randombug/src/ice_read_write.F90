!=======================================================================
!BOP
!
! !MODULE: ice_read_write
!
! !DESCRIPTION:
!
! Routines for opening, reading and writing external files
!
! !REVISION HISTORY:
!  SVN:$Id$
!
! author: Tony Craig, NCAR
!
! 2004: Block structure added by William Lipscomb, LANL
! 2006: Converted to free source form (F90) by Elizabeth Hunke
! 2007: netcdf versions added by Alison McLaren & Ann Keen, Met Office
!
! !INTERFACE:
!
      module ice_read_write
!
! !USES:
!
      use ice_kinds_mod
      use ice_constants
      use ice_communicate, only: my_task, master_task
      use ice_broadcast
      use ice_domain_size
      use ice_blocks
      use ice_fileunits
      use netcdf      
!
!EOP

      implicit none

       public :: ice_read_global_nc

       interface ice_read_global_nc
          module procedure ice_read_global_nc_dbl, &
                           ice_read_global_nc_r4
       end interface

!=======================================================================

      contains

!=======================================================================
!
!BOP
!
! !IROUTINE: ice_open - opens an unformatted file for reading
!
! !INTERFACE:
!
! Subprogram not used       subroutine ice_open(nu, filename, nbits)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Opens an unformatted file for reading \! nbits indicates whether the file is sequential or direct access
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author: Tony Craig, NCAR
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) :: &
! Subprogram not used            nu        , & ! unit number
! Subprogram not used            nbits         ! no. of bits per variable (0 for sequential access)
! Subprogram not used 
! Subprogram not used       character (*) :: filename
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used 
! Subprogram not used          if (nbits == 0) then   ! sequential access
! Subprogram not used 
! Subprogram not used             open(nu,file=filename,form='unformatted')
! Subprogram not used 
! Subprogram not used          else                   ! direct access
! Subprogram not used             open(nu,file=filename,recl=nx_global*ny_global*nbits/8, &
! Subprogram not used                   form='unformatted',access='direct')
! Subprogram not used          endif                   ! nbits = 0
! Subprogram not used 
! Subprogram not used       endif                      ! my_task = master_task
! Subprogram not used 
! Subprogram not used       end subroutine ice_open

!=======================================================================
!BOP
!
! !IROUTINE: ice_read - read and scatter an unformatted file
!
! !INTERFACE:
!
! Subprogram not used       subroutine ice_read(nu,  nrec,  work, atype, diag, &
! Subprogram not used                           field_loc, field_type, &
! Subprogram not used                           ignore_eof, hit_eof)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Read an unformatted file and scatter to processors\! work is a real array, atype indicates the format of the data\! If the optional variables field_loc and field_type are present \! the ghost cells are filled using values from the global array.\! This prevents them from being filled with zeroes in land cells \! (subroutine ice_HaloUpdate need not be called).
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author: Tony Craig, NCAR
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use ice_domain
! Subprogram not used       use ice_gather_scatter
! Subprogram not used       use ice_work, only: work_g1, work_gr, work_gi4, work_gi8
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) :: &
! Subprogram not used            nu            , & ! unit number
! Subprogram not used            nrec              ! record number (0 for sequential access)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks), &
! Subprogram not used            intent(out) :: &
! Subprogram not used            work              ! output array (real, 8-byte)
! Subprogram not used 
! Subprogram not used       character (len=4), intent(in) :: &
! Subprogram not used            atype             ! format for input array
! Subprogram not used                              ! (real/integer, 4-byte/8-byte)
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind), intent(in) :: &
! Subprogram not used            diag              ! if true, write diagnostic output
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), optional, intent(in) :: &
! Subprogram not used            field_loc, &      ! location of field on staggered grid
! Subprogram not used            field_type        ! type of field (scalar, vector, angle)
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind), optional, intent(in)  :: ignore_eof
! Subprogram not used       logical (kind=log_kind), optional, intent(out) :: hit_eof
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) :: i, j, ios
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind) :: &
! Subprogram not used          amin, amax         ! min and max values of input array
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind) :: ignore_eof_use
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used          allocate(work_g1(nx_global,ny_global))
! Subprogram not used       else
! Subprogram not used          allocate(work_g1(1,1))   ! to save memory
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Read global array according to format atype
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used          if (present(hit_eof)) hit_eof = .false.
! Subprogram not used 
! Subprogram not used          if (atype == 'ida4') then
! Subprogram not used             allocate(work_gi4(nx_global,ny_global))
! Subprogram not used             read(nu,rec=nrec) work_gi4
! Subprogram not used             work_g1 = real(work_gi4,kind=dbl_kind)
! Subprogram not used             deallocate(work_gi4)
! Subprogram not used          elseif (atype == 'ida8') then
! Subprogram not used             allocate(work_gi8(nx_global,ny_global))
! Subprogram not used             read(nu,rec=nrec) work_gi8
! Subprogram not used             work_g1 = real(work_gi8,kind=dbl_kind)
! Subprogram not used             deallocate(work_gi8)
! Subprogram not used          elseif (atype == 'rda4') then
! Subprogram not used             allocate(work_gr(nx_global,ny_global))
! Subprogram not used             read(nu,rec=nrec) work_gr
! Subprogram not used             work_g1 = work_gr
! Subprogram not used             deallocate(work_gr)
! Subprogram not used          elseif (atype == 'rda8') then
! Subprogram not used             read(nu,rec=nrec) work_g1
! Subprogram not used          elseif (atype == 'ruf8') then
! Subprogram not used             if (present(ignore_eof)) then
! Subprogram not used                ignore_eof_use = ignore_eof
! Subprogram not used             else
! Subprogram not used                ignore_eof_use = .false.
! Subprogram not used             endif
! Subprogram not used             if (ignore_eof_use) then
! Subprogram not used              ! Read line from file, checking for end-of-file
! Subprogram not used                read(nu, iostat=ios) ((work_g1(i,j),i=1,nx_global), &
! Subprogram not used                                                    j=1,ny_global)
! Subprogram not used                if (present(hit_eof)) hit_eof = ios < 0
! Subprogram not used             else
! Subprogram not used                read(nu) ((work_g1(i,j),i=1,nx_global),j=1,ny_global)
! Subprogram not used             endif
! Subprogram not used          else
! Subprogram not used             write(nu_diag,*) ' ERROR: reading unknown atype ',atype
! Subprogram not used          endif
! Subprogram not used       endif                     ! my_task = master_task
! Subprogram not used 
! Subprogram not used       if (present(hit_eof)) then
! Subprogram not used          call broadcast_scalar(hit_eof,master_task)
! Subprogram not used          if (hit_eof) then
! Subprogram not used             deallocate(work_g1)
! Subprogram not used             return
! Subprogram not used          endif
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! optional diagnostics
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used       if (my_task==master_task .and. diag) then
! Subprogram not used          amin = minval(work_g1)
! Subprogram not used          amax = maxval(work_g1, mask = work_g1 /= spval_dbl)
! Subprogram not used          write(nu_diag,*) ' read_global ',nu, nrec, amin, amax
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Scatter data to individual processors.
! Subprogram not used     ! NOTE: Ghost cells are not updated unless field_loc is present.
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       if (present(field_loc)) then
! Subprogram not used          call scatter_global(work, work_g1, master_task, distrb_info, &
! Subprogram not used                              field_loc, field_type)
! Subprogram not used       else
! Subprogram not used          call scatter_global(work, work_g1, master_task, distrb_info, &
! Subprogram not used                              field_loc_noupdate, field_type_noupdate)
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       deallocate(work_g1)
! Subprogram not used 
! Subprogram not used       end subroutine ice_read

!=======================================================================
!BOP
!
! !IROUTINE: ice_read_global - read an unformatted file
!
! !INTERFACE:
!
! Subprogram not used       subroutine ice_read_global (nu,  nrec,  work_g, atype, diag, &
! Subprogram not used                                   ignore_eof, hit_eof)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Read an unformatted file \! Just like ice_read except that it returns a global array \! work_g is a real array, atype indicates the format of the data
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! Adapted by William Lipscomb, LANL, from ice_read
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use ice_work, only: work_gr, work_gi4, work_gi8
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) :: &
! Subprogram not used            nu            , & ! unit number
! Subprogram not used            nrec              ! record number (0 for sequential access)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension(:,:), &
! Subprogram not used            intent(out) :: &
! Subprogram not used            work_g            ! output array (real, 8-byte)
! Subprogram not used 
! Subprogram not used       character (len=4) :: &
! Subprogram not used            atype             ! format for input array
! Subprogram not used                              ! (real/integer, 4-byte/8-byte)
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind) :: &
! Subprogram not used            diag              ! if true, write diagnostic output
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind), optional, intent(in)  :: ignore_eof
! Subprogram not used       logical (kind=log_kind), optional, intent(out) :: hit_eof
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) :: i, j, ios
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind) :: &
! Subprogram not used          amin, amax         ! min and max values of input array
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind) :: ignore_eof_use
! Subprogram not used 
! Subprogram not used       work_g(:,:) = c0
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Read global array according to format atype
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used          if (present(hit_eof)) hit_eof = .false.
! Subprogram not used 
! Subprogram not used          if (atype == 'ida4') then
! Subprogram not used             allocate(work_gi4(nx_global,ny_global))
! Subprogram not used             read(nu,rec=nrec) work_gi4
! Subprogram not used             work_g = real(work_gi4,kind=dbl_kind)
! Subprogram not used             deallocate(work_gi4)
! Subprogram not used          elseif (atype == 'ida8') then
! Subprogram not used             allocate(work_gi8(nx_global,ny_global))
! Subprogram not used             read(nu,rec=nrec) work_gi8
! Subprogram not used             work_g = real(work_gi8,kind=dbl_kind)
! Subprogram not used             deallocate(work_gi8)
! Subprogram not used          elseif (atype == 'rda4') then
! Subprogram not used             allocate(work_gr(nx_global,ny_global))
! Subprogram not used             read(nu,rec=nrec) work_gr
! Subprogram not used             work_g = work_gr
! Subprogram not used             deallocate(work_gr)
! Subprogram not used          elseif (atype == 'rda8') then
! Subprogram not used             read(nu,rec=nrec) work_g
! Subprogram not used          elseif (atype == 'ruf8') then
! Subprogram not used             if (present(ignore_eof)) then
! Subprogram not used                ignore_eof_use = ignore_eof
! Subprogram not used             else
! Subprogram not used                ignore_eof_use = .false.
! Subprogram not used             endif
! Subprogram not used             if (ignore_eof_use) then
! Subprogram not used                ! Read line from file, checking for end-of-file
! Subprogram not used                read(nu, iostat=ios) ((work_g(i,j),i=1,nx_global), &
! Subprogram not used                                                   j=1,ny_global)
! Subprogram not used                if (present(hit_eof)) hit_eof = ios < 0
! Subprogram not used             else
! Subprogram not used                read(nu) ((work_g(i,j),i=1,nx_global),j=1,ny_global)
! Subprogram not used             endif
! Subprogram not used          else
! Subprogram not used             write(nu_diag,*) ' ERROR: reading unknown atype ',atype
! Subprogram not used          endif
! Subprogram not used       endif                     ! my_task = master_task
! Subprogram not used 
! Subprogram not used       if (present(hit_eof)) then
! Subprogram not used          call broadcast_scalar(hit_eof,master_task)
! Subprogram not used          if (hit_eof) return
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! optional diagnostics
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used       if (my_task == master_task .and. diag) then
! Subprogram not used          amin = minval(work_g)
! Subprogram not used          amax = maxval(work_g, mask = work_g /= spval_dbl)
! Subprogram not used          write(nu_diag,*) ' read_global ',nu, nrec, amin, amax
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       end subroutine ice_read_global

!=======================================================================
!BOP
!
! !IROUTINE: ice_write - writes an unformatted file
!
! !INTERFACE:
!
! Subprogram not used       subroutine ice_write(nu, nrec, work, atype, diag)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Writes an unformatted file \! work is a real array, atype indicates the format of the data
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author: Tony Craig, NCAR
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use ice_gather_scatter
! Subprogram not used       use ice_domain
! Subprogram not used       use ice_work, only: work_g1, work_gr, work_gi4, work_gi8
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) :: &
! Subprogram not used            nu            , & ! unit number
! Subprogram not used            nrec              ! record number (0 for sequential access)
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks), &
! Subprogram not used            intent(in) :: &
! Subprogram not used            work              ! input array (real, 8-byte)
! Subprogram not used 
! Subprogram not used       character (len=4) :: &
! Subprogram not used            atype             ! format for output array
! Subprogram not used                              ! (real/integer, 4-byte/8-byte)
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind) :: &
! Subprogram not used            diag              ! if true, write diagnostic output
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) :: i, j
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind) :: &
! Subprogram not used          amin, amax     ! min and max values of ouput array
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Gather data from individual processors
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used          allocate(work_g1(nx_global,ny_global))
! Subprogram not used       else
! Subprogram not used          allocate(work_g1(1,1)) ! to save memory
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       call gather_global(work_g1, work, master_task, distrb_info)
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Write global array according to format atype
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used          if (atype == 'ida4') then
! Subprogram not used             allocate(work_gi4(nx_global,ny_global))
! Subprogram not used             work_gi4 = nint(work_g1)
! Subprogram not used             write(nu,rec=nrec) work_gi4
! Subprogram not used             deallocate(work_gi4)
! Subprogram not used          elseif (atype == 'ida8') then
! Subprogram not used             allocate(work_gi8(nx_global,ny_global))
! Subprogram not used             work_gi8 = nint(work_g1)
! Subprogram not used             write(nu,rec=nrec) work_gi8           
! Subprogram not used             deallocate(work_gi8)
! Subprogram not used          elseif (atype == 'rda4') then
! Subprogram not used             allocate(work_gr(nx_global,ny_global))
! Subprogram not used             work_gr = work_g1
! Subprogram not used             write(nu,rec=nrec) work_gr
! Subprogram not used             deallocate(work_gr)
! Subprogram not used          elseif (atype == 'rda8') then
! Subprogram not used             write(nu,rec=nrec) work_g1
! Subprogram not used          elseif (atype == 'ruf8') then
! Subprogram not used             write(nu) ((work_g1(i,j),i=1,nx_global),j=1,ny_global)
! Subprogram not used          else
! Subprogram not used             write(nu_diag,*) ' ERROR: writing unknown atype ',atype
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! diagnostics
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used          if (diag) then
! Subprogram not used             amin = minval(work_g1)
! Subprogram not used             amax = maxval(work_g1, mask = work_g1 /= spval_dbl)
! Subprogram not used             write(nu_diag,*) ' write_global ', nu, nrec, amin, amax
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used       endif                     ! my_task = master_task
! Subprogram not used 
! Subprogram not used       deallocate(work_g1)
! Subprogram not used 
! Subprogram not used       end subroutine ice_write

!=======================================================================
!BOP
!
! !IROUTINE: ice_write_nc - writes a field to a netcdf file
!
! !INTERFACE:
!
! Subprogram not used       subroutine ice_write_nc(fid, nrec, varname, work, atype, diag)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Writes a field to a netcdf file \! work is a real array, atype indicates the format of the data
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !
! Subprogram not used ! author: David A Bailey, NCAR
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use ice_gather_scatter
! Subprogram not used       use ice_domain
! Subprogram not used       use ice_work, only: work_g1, work_gr, work_gi4, work_gi8
! Subprogram not used       use ice_exit
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) :: &
! Subprogram not used            fid          ,&   ! netcdf file id
! Subprogram not used            nrec              ! record number
! Subprogram not used 
! Subprogram not used       character (len=*), intent(in) :: varname
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks), &
! Subprogram not used            intent(in) :: &
! Subprogram not used            work              ! input array (real, 8-byte)
! Subprogram not used 
! Subprogram not used       character (len=4), intent(in) :: &
! Subprogram not used            atype             ! format for output array
! Subprogram not used                              ! (real/integer, 4-byte/8-byte)
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind), intent(in) :: &
! Subprogram not used            diag              ! if true, write diagnostic output
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind) :: i, j, varid, numDims
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind) :: &
! Subprogram not used          status        ! status variable from netCDF routine 
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind) :: &
! Subprogram not used          amin, amax     ! min and max values of ouput array
! Subprogram not used 
! Subprogram not used       integer (kind=int_kind), allocatable :: &
! Subprogram not used          start_arr(:), count_arr(:)
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Gather data from individual processors
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used          allocate(work_g1(nx_global,ny_global))
! Subprogram not used       else
! Subprogram not used          allocate(work_g1(1,1)) ! to save memory
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       call gather_global(work_g1, work, master_task, distrb_info)
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used 
! Subprogram not used          status = nf90_inq_varid(fid, trim(varname), varid)
! Subprogram not used 
! Subprogram not used          if (status /= nf90_noerr) then
! Subprogram not used            call abort_ice ( &
! Subprogram not used                'ice_write_nc: Cannot find variable '//trim(varname) )
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used          status = nf90_inquire_variable(fid, varid, ndims = numDims)
! Subprogram not used 
! Subprogram not used          if (status /= nf90_noerr) then
! Subprogram not used            call abort_ice ( &
! Subprogram not used                'ice_write_nc: Cannot find dimensions for '//trim(varname) )
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used          allocate(start_arr(numDims))
! Subprogram not used          allocate(count_arr(numDims))
! Subprogram not used 
! Subprogram not used          if (numDims > 2) then
! Subprogram not used             start_arr(1) = 1
! Subprogram not used             start_arr(2) = 1
! Subprogram not used             start_arr(3) = nrec
! Subprogram not used             count_arr(1) = nx_global
! Subprogram not used             count_arr(2) = ny_global
! Subprogram not used             count_arr(3) = 1
! Subprogram not used          else
! Subprogram not used             start_arr(1) = 1
! Subprogram not used             start_arr(2) = 1
! Subprogram not used             count_arr(1) = nx_global
! Subprogram not used             count_arr(2) = ny_global
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! Write global array according to format atype
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used          if (atype == 'ida4') then
! Subprogram not used             allocate(work_gi4(nx_global,ny_global))
! Subprogram not used             work_gi4 = nint(work_g1)
! Subprogram not used             status = nf90_put_var(fid,varid,work_gi4, &
! Subprogram not used                                   start=start_arr,  &
! Subprogram not used                                   count=count_arr)
! Subprogram not used             deallocate(work_gi4)
! Subprogram not used          elseif (atype == 'ida8') then
! Subprogram not used             allocate(work_gi8(nx_global,ny_global))
! Subprogram not used             work_gi8 = nint(work_g1)
! Subprogram not used             status = nf90_put_var(fid,varid,work_gi8, &
! Subprogram not used                                   start=start_arr,  &
! Subprogram not used                                   count=count_arr)
! Subprogram not used             deallocate(work_gi8)
! Subprogram not used          elseif (atype == 'rda4') then
! Subprogram not used             allocate(work_gr(nx_global,ny_global))
! Subprogram not used             work_gr = work_g1
! Subprogram not used             status = nf90_put_var(fid,varid,work_gr, &
! Subprogram not used                                   start=start_arr,  &
! Subprogram not used                                   count=count_arr)
! Subprogram not used             deallocate(work_gr)
! Subprogram not used          elseif (atype == 'rda8') then
! Subprogram not used             status = nf90_put_var(fid,varid,work_g1, &
! Subprogram not used                                   start=start_arr,  &
! Subprogram not used                                   count=count_arr)
! Subprogram not used          else
! Subprogram not used             write(nu_diag,*) ' ERROR: writing unknown atype ',atype
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! diagnostics
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used          if (diag) then
! Subprogram not used             amin = minval(work_g1)
! Subprogram not used             amax = maxval(work_g1, mask = work_g1 /= spval_dbl)
! Subprogram not used             write(nu_diag,*) ' write_global ', fid, varid, nrec, amin, amax
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used          deallocate(start_arr)
! Subprogram not used          deallocate(count_arr)
! Subprogram not used 
! Subprogram not used       endif                     ! my_task = master_task
! Subprogram not used 
! Subprogram not used       deallocate(work_g1)
! Subprogram not used 
! Subprogram not used       end subroutine ice_write_nc
      
!=======================================================================
!
!BOP
!
! !IROUTINE: ice_open_nc - opens a netCDF file for reading
!
! !INTERFACE:
!
      subroutine ice_open_nc(filename, fid)
!
! !DESCRIPTION:
!
! Opens a netCDF file for reading
!
! !REVISION HISTORY:
!
! Adapted by Alison McLaren, Met Office from ice_open
!
! !USES:
 
      use ice_exit
!
! !INPUT/OUTPUT PARAMETERS:
!

      character (char_len_long), intent(in) :: & 
           filename      ! netCDF filename

      integer (kind=int_kind), intent(out) :: &
           fid           ! unit number
!
!EOP
!
      integer (kind=int_kind) :: &
        status        ! status variable from netCDF routine 


      if (my_task == master_task) then

          status = nf90_open(filename, NF90_NOWRITE, fid)
          if (status /= nf90_noerr) then
             call abort_ice ( & 
                   'ice_open_nc: Cannot open '//trim(filename) )
          endif

      endif                      ! my_task = master_task

      end subroutine ice_open_nc

!=======================================================================
!BOP
!
! !IROUTINE: ice_read_nc - read and scatter one field from a netCDF file
!
! !INTERFACE:
!
      subroutine ice_read_nc(fid,  nrec,  varname, work,  diag, &
                             field_loc, field_type)
!
! !DESCRIPTION:
!
! Read a netCDF file and scatter to processors\! If the optional variables field_loc and field_type are present \! the ghost cells are filled using values from the global array.\! This prevents them from being filled with zeroes in land cells \! (subroutine ice_HaloUpdate need not be called).
!
! !REVISION HISTORY:
!
! Adapted by Alison McLaren, Met Office from ice_read
!
! !USES:
!
      use ice_domain
      use ice_gather_scatter
      use ice_work, only: work_g1
      use ice_exit
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
           fid           , & ! file id
           nrec              ! record number 

      logical (kind=log_kind), intent(in) :: &
           diag              ! if true, write diagnostic output

      character (len=*), intent(in) :: & 
           varname           ! field name in netcdf file

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks), &
           intent(out) :: &
           work              ! output array (real, 8-byte)

      integer (kind=int_kind), optional, intent(in) :: &
           field_loc, &      ! location of field on staggered grid
           field_type        ! type of field (scalar, vector, angle)
!
!EOP
!
! netCDF file diagnostics:
      integer (kind=int_kind) :: & 
         varid,           & ! netcdf id for field
         status,          & ! status output from netcdf routines
         ndim, nvar,      & ! sizes of netcdf file
         id,              & ! dimension index
         dimlen             ! size of dimension

      real (kind=dbl_kind) :: &
         amin, amax         ! min and max values of input array

      character (char_len) :: &
         dimname            ! dimension name            
!
      if (my_task == master_task) then
         allocate(work_g1(nx_global,ny_global))
      else
         allocate(work_g1(1,1))   ! to save memory
      endif


      if (my_task == master_task) then

        !-------------------------------------------------------------
        ! Find out ID of required variable
        !-------------------------------------------------------------

         status = nf90_inq_varid(fid, trim(varname), varid)
 
         if (status /= nf90_noerr) then
           call abort_ice ( & 
               'ice_read_nc: Cannot find variable '//trim(varname) )
         endif

       !--------------------------------------------------------------
       ! Read global array 
       !--------------------------------------------------------------

         status = nf90_get_var( fid, varid, work_g1, &
               start=(/1,1,nrec/), & 
               count=(/nx_global,ny_global,1/) )

      endif                     ! my_task = master_task

    !-------------------------------------------------------------------
    ! optional diagnostics
    !-------------------------------------------------------------------

      if (my_task==master_task .and. diag) then

!        write(nu_diag,*) & 
!          'ice_read_nc, fid= ',fid, ', nrec = ',nrec, & 
!          ', varname = ',trim(varname)
         status = nf90_inquire(fid, nDimensions=ndim, nVariables=nvar)
!        write(nu_diag,*) 'ndim= ',ndim,', nvar= ',nvar
         do id=1,ndim
           status = nf90_inquire_dimension(fid,id,name=dimname,len=dimlen)
!          write(nu_diag,*) 'Dim name = ',trim(dimname),', size = ',dimlen
         enddo
         amin = minval(work_g1)
         amax = maxval(work_g1, mask = work_g1 /= spval_dbl)
         write(nu_diag,*) ' read_global ',fid, varid, nrec, amin, amax

      endif

    !-------------------------------------------------------------------
    ! Scatter data to individual processors.
    ! NOTE: Ghost cells are not updated unless field_loc is present.
    !-------------------------------------------------------------------

      if (present(field_loc)) then
         call scatter_global(work, work_g1, master_task, distrb_info, &
                             field_loc, field_type)
      else
         call scatter_global(work, work_g1, master_task, distrb_info, &
                             field_loc_noupdate, field_type_noupdate)
      endif

      deallocate(work_g1)

      end subroutine ice_read_nc
!
!=======================================================================
!BOP
!
! !IROUTINE: ice_read_global_nc - read one field from a netcdf file
!
! !INTERFACE:
!
! Subprogram not used       subroutine ice_read_global_nc_dbl (fid,  nrec, varname, work_g, diag)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Read a netcdf file \! Just like ice_read_nc except that it returns a global array \! work_g is a real array
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! Adapted by William Lipscomb, LANL, from ice_read
! Subprogram not used ! Adapted by Ann Keen, Met Office, to read from a netcdf file 
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used ! 
! Subprogram not used       use ice_exit
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) :: &
! Subprogram not used            fid           , & ! file id
! Subprogram not used            nrec              ! record number 
! Subprogram not used 
! Subprogram not used      character (len=*), intent(in) :: & 
! Subprogram not used            varname           ! field name in netcdf file        
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind), dimension(:,:), &
! Subprogram not used            intent(out) :: &
! Subprogram not used            work_g            ! output array (real, 8-byte)
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind) :: &
! Subprogram not used            diag              ! if true, write diagnostic output
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used ! netCDF file diagnostics:
! Subprogram not used       integer (kind=int_kind) :: & 
! Subprogram not used          varid,           & ! netcdf id for field
! Subprogram not used          status,          & ! status output from netcdf routines
! Subprogram not used          ndim, nvar,      & ! sizes of netcdf file
! Subprogram not used          id,              & ! dimension index
! Subprogram not used          dimlen             ! size of dimension      
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind) :: &
! Subprogram not used          amin, amax         ! min and max values of input array
! Subprogram not used 
! Subprogram not used      character (char_len) :: &
! Subprogram not used          dimname            ! dimension name            
! Subprogram not used 
! Subprogram not used !
! Subprogram not used       work_g(:,:) = c0
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used 
! Subprogram not used         !-------------------------------------------------------------
! Subprogram not used         ! Find out ID of required variable
! Subprogram not used         !-------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          status = nf90_inq_varid(fid, trim(varname), varid)
! Subprogram not used 
! Subprogram not used          if (status /= nf90_noerr) then
! Subprogram not used            call abort_ice ( & 
! Subprogram not used             'ice_read_global_nc: Cannot find variable '//trim(varname) )
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used        !--------------------------------------------------------------
! Subprogram not used        ! Read global array 
! Subprogram not used        !--------------------------------------------------------------
! Subprogram not used  
! Subprogram not used          status = nf90_get_var( fid, varid, work_g, &
! Subprogram not used                start=(/1,1,nrec/), & 
! Subprogram not used                count=(/nx_global,ny_global,1/) )
! Subprogram not used 
! Subprogram not used       endif                     ! my_task = master_task
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! optional diagnostics
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       if (my_task == master_task .and. diag) then
! Subprogram not used 
! Subprogram not used !         write(nu_diag,*) & 
! Subprogram not used !           'ice_read_global_nc, fid= ',fid, ', nrec = ',nrec, & 
! Subprogram not used !           ', varname = ',trim(varname)
! Subprogram not used           status = nf90_inquire(fid, nDimensions=ndim, nVariables=nvar)
! Subprogram not used !         write(nu_diag,*) 'ndim= ',ndim,', nvar= ',nvar
! Subprogram not used           do id=1,ndim
! Subprogram not used             status = nf90_inquire_dimension(fid,id,name=dimname,len=dimlen)
! Subprogram not used !           write(nu_diag,*) 'Dim name = ',trim(dimname),', size = ',dimlen
! Subprogram not used          enddo
! Subprogram not used          amin = minval(work_g)
! Subprogram not used          amax = maxval(work_g, mask = work_g /= spval_dbl)
! Subprogram not used !        write(nu_diag,*) 'min and max = ', amin, amax
! Subprogram not used !        write(nu_diag,*) ''
! Subprogram not used          write(nu_diag,*) ' read_global ',fid, varid, nrec, amin, amax
! Subprogram not used 
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used       end subroutine ice_read_global_nc_dbl
!=======================================================================
!BOP
!
! !IROUTINE: ice_read_global_nc - read one field from a netcdf file
!
! !INTERFACE:
!
! Subprogram not used       subroutine ice_read_global_nc_r4 (fid,  nrec, varname, work_g, diag)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used !
! Subprogram not used ! Read a netcdf file \! Just like ice_read_nc except that it returns a global array \! work_g is a real array
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! Adapted by William Lipscomb, LANL, from ice_read
! Subprogram not used ! Adapted by Ann Keen, Met Office, to read from a netcdf file 
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used ! 
! Subprogram not used       use ice_exit
! Subprogram not used !
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       integer (kind=int_kind), intent(in) :: &
! Subprogram not used            fid           , & ! file id
! Subprogram not used            nrec              ! record number 
! Subprogram not used 
! Subprogram not used      character (len=*), intent(in) :: & 
! Subprogram not used            varname           ! field name in netcdf file        
! Subprogram not used 
! Subprogram not used       real (kind=real_kind), dimension(:,:), &
! Subprogram not used            intent(out) :: &
! Subprogram not used            work_g            ! output array (real, 8-byte)
! Subprogram not used 
! Subprogram not used       logical (kind=log_kind) :: &
! Subprogram not used            diag              ! if true, write diagnostic output
! Subprogram not used !
! Subprogram not used !EOP
! Subprogram not used !
! Subprogram not used ! netCDF file diagnostics:
! Subprogram not used       integer (kind=int_kind) :: & 
! Subprogram not used          varid,           & ! netcdf id for field
! Subprogram not used          status,          & ! status output from netcdf routines
! Subprogram not used          ndim, nvar,      & ! sizes of netcdf file
! Subprogram not used          id,              & ! dimension index
! Subprogram not used          dimlen             ! size of dimension      
! Subprogram not used 
! Subprogram not used       real (kind=dbl_kind) :: &
! Subprogram not used          amin, amax         ! min and max values of input array
! Subprogram not used 
! Subprogram not used      character (char_len) :: &
! Subprogram not used          dimname            ! dimension name            
! Subprogram not used 
! Subprogram not used !
! Subprogram not used       work_g(:,:) = c0
! Subprogram not used 
! Subprogram not used       if (my_task == master_task) then
! Subprogram not used 
! Subprogram not used         !-------------------------------------------------------------
! Subprogram not used         ! Find out ID of required variable
! Subprogram not used         !-------------------------------------------------------------
! Subprogram not used 
! Subprogram not used          status = nf90_inq_varid(fid, trim(varname), varid)
! Subprogram not used 
! Subprogram not used          if (status /= nf90_noerr) then
! Subprogram not used            call abort_ice ( & 
! Subprogram not used             'ice_read_global_nc: Cannot find variable '//trim(varname) )
! Subprogram not used          endif
! Subprogram not used 
! Subprogram not used        !--------------------------------------------------------------
! Subprogram not used        ! Read global array 
! Subprogram not used        !--------------------------------------------------------------
! Subprogram not used  
! Subprogram not used          status = nf90_get_var( fid, varid, work_g, &
! Subprogram not used                start=(/1,1,nrec/), & 
! Subprogram not used                count=(/nx_global,ny_global,1/) )
! Subprogram not used 
! Subprogram not used       endif                     ! my_task = master_task
! Subprogram not used 
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used     ! optional diagnostics
! Subprogram not used     !-------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used       if (my_task == master_task .and. diag) then
! Subprogram not used 
! Subprogram not used !         write(nu_diag,*) & 
! Subprogram not used !           'ice_read_global_nc, fid= ',fid, ', nrec = ',nrec, & 
! Subprogram not used !           ', varname = ',trim(varname)
! Subprogram not used           status = nf90_inquire(fid, nDimensions=ndim, nVariables=nvar)
! Subprogram not used !         write(nu_diag,*) 'ndim= ',ndim,', nvar= ',nvar
! Subprogram not used           do id=1,ndim
! Subprogram not used             status = nf90_inquire_dimension(fid,id,name=dimname,len=dimlen)
! Subprogram not used !           write(nu_diag,*) 'Dim name = ',trim(dimname),', size = ',dimlen
! Subprogram not used          enddo
! Subprogram not used          amin = minval(work_g)
! Subprogram not used          amax = maxval(work_g, mask = work_g /= spval_dbl)
! Subprogram not used !        write(nu_diag,*) 'min and max = ', amin, amax
! Subprogram not used !        write(nu_diag,*) ''
! Subprogram not used          write(nu_diag,*) ' read_global ',fid, varid, nrec, amin, amax
! Subprogram not used 
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used       end subroutine ice_read_global_nc_r4

!=======================================================================
!BOP
!
! !IROUTINE: ice_close_nc - closes a netCDF file
!
! !INTERFACE:
!
      subroutine ice_close_nc(fid)
!
! !DESCRIPTION:
!
! Closes a netCDF file
!
! !REVISION HISTORY:
!
! author: Alison McLaren, Met Office
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
           fid           ! unit number
!
!EOP
!
      integer (kind=int_kind) :: &
        status        ! status variable from netCDF routine 

      if (my_task == master_task) then

         status = nf90_close(fid)

      endif

      end subroutine ice_close_nc

!=======================================================================

      end module ice_read_write

!=======================================================================
