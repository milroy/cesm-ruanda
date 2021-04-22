module pio_utils
  use pio_types, only : file_desc_t, var_desc_t
  use pio_types, only : pio_int, pio_real, pio_double, pio_char
  use pio_types, only : iotype_netcdf, iotype_pnetcdf, PIO_internal_error
  use pio_types, only : PIO_iotype_netcdf4p, pio_iotype_netcdf4c
  use pio_types, only : PIO_bcast_error 
  use pio_kinds, only : i4, r4, r8
  use pio_support, only : checkmpireturn, piodie, Debug





  use mpi !_EXTERNAL




  implicit none
  private

  public :: check_netcdf 
  public :: bad_iotype 

  

contains

  subroutine check_netcdf(File, status, filestr, line)
    type(file_desc_t), intent(in) :: file
    integer, intent(inout) :: status
    character(len=*), intent(in) :: filestr
    integer, intent(in) :: line

    integer :: mpierr, iotype

!  Three choices for error handling:
!  1: abort on error from any task           PIO_INTERNAL_ERROR
!  2: broadcast an error from io_rank 0      PIO_BCAST_ERROR
!  3: do nothing - allow the user to handle it PIO_RETURN_ERROR
!
    iotype = file%iotype
    
    if(Debug) call mpi_barrier(file%iosystem%union_comm, mpierr)

    select case(iotype)
    case(iotype_pnetcdf)
    case(iotype_netcdf,pio_iotype_netcdf4p,pio_iotype_netcdf4c)
    end select

  end subroutine check_netcdf



!>
!! @private
!<
  subroutine bad_iotype(iotype,file,line)
    integer iotype
    character(len=*) file
    integer line

    if (iotype==iotype_pnetcdf) then
       call piodie(file,line,'PNETCDF not enabled in the build')
    endif
    if (iotype==iotype_netcdf) then
       call piodie(file,line,'NETCDF not enabled in the build')
    endif
    if (iotype==PIO_iotype_netcdf4p .or. iotype==pio_iotype_netcdf4c) then
       call piodie(file,line,'NETCDF4 not enabled in the build')
    endif
    print *,'Invalid iotype, value=',iotype
    call piodie(file,line,'Quitting')

  end subroutine bad_iotype

  

end module pio_utils
