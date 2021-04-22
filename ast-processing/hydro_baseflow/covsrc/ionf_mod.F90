
module ionf_mod



  use alloc_mod

  use pio_kinds, only: i4,r4,r8,pio_offset
  use pio_types
  use pio_utils, only: bad_iotype, check_netcdf

  use pio_support, only : Debug, DebugIO, piodie, DebugAsync   



  use pio_support, only : CheckMPIReturn



  implicit none
  private




 

   public :: create_nf
   public :: open_nf 
   public :: close_nf 
   public :: sync_nf 

contains 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! create_nf
  !

  integer function create_nf(File,fname, amode) result(ierr)

    use mpi ! _EXTERNAL



    type (File_desc_t), intent(inout) :: File
    character(len=*), intent(in)      :: fname
    integer(i4),  intent(in) :: amode
    integer(i4) :: iotype, mpierr
    integer :: nmode, tmpfh

    nmode=amode
    
    ierr=PIO_noerr
    File%fh=-1

    if(File%iosystem%ioproc) then
       iotype = File%iotype 
       select case (iotype) 
       case default
          call bad_iotype(iotype,"ionf_mod.F90",122)

       end select
       if(Debug) print *,"ionf_mod.F90",125,file%fh,ierr
    end if
    tmpfh = file%fh
    
    call mpi_bcast(tmpfh,1,mpi_integer, file%iosystem%iomaster, file%iosystem%my_comm, mpierr)
    
    if(.not. file%iosystem%ioproc) file%fh=-tmpfh

    if(Debug.or.DebugAsync) print *,"ionf_mod.F90",133,file%fh,ierr
    
    call check_netcdf(File, ierr,"ionf_mod.F90",135)

  end function create_nf



!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! open_nf
  ! 

  integer function open_nf(File,fname, mode) result(ierr)
    use mpi ! _EXTERNAL
    type (File_desc_t), intent(inout) :: File
    character(len=*), intent(in)      :: fname
    integer(i4), optional, intent(in) :: mode
    integer(i4) :: iotype, amode , mpierr, ier2
    integer :: tmpfh, format


    ierr=PIO_noerr
    File%fh=-1
    if(file%iosystem%ioproc) then
!       This subroutine seems to break pgi compiler for large files.
!       call check_file_type(File, fname)
       iotype = File%iotype 

    end if

    tmpfh = file%fh
    call mpi_bcast(tmpfh,1,mpi_integer, file%iosystem%iomaster, file%iosystem%my_comm, mpierr)

    if(.not. file%iosystem%ioproc) file%fh=-tmpfh
      
    call check_netcdf(File, ierr,"ionf_mod.F90",235)

  end function open_nf



!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! close_nf
  !


  integer function close_nf(File) result(ierr)
    type (File_desc_t), intent(inout) :: File

    ierr=PIO_noerr

    if(File%iosystem%IOproc) then
       if(Debug) print *,"ionf_mod.F90",252,'CFILE closing : ',file%fh
       select case (File%iotype) 
       case default
          call bad_iotype(File%iotype,"ionf_mod.F90",268)
       end select
    end if
    file%fh=-1
    call check_netcdf(File, ierr,"ionf_mod.F90",272)
  end function close_nf


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! sync_nf
  !


  integer function sync_nf(File) result(ierr)

    type (File_desc_t), intent(inout) :: File

    ierr=PIO_noerr

    if(File%iosystem%IOproc) then
       if(Debug) print *,"ionf_mod.F90",288,'CFILE syncing : ',file%fh
       select case (File%iotype) 
       case default
          call bad_iotype(File%iotype,"ionf_mod.F90",301)
       end select
    end if
    call check_netcdf(File, ierr,"ionf_mod.F90",304)
  end function sync_nf

  subroutine check_file_type(File, filename) 
    use mpi !_EXTERNAL

    type (File_desc_t), intent(inout) :: File
    character(len=*), intent(in) :: filename
    character(len=4) :: magic
    integer :: fh, mpierr, reclength=4, i, eof 
    logical :: UNITOK, UNITOP

!   Check format of existing files opened to read.

    inquire(file=filename, exist=UNITOK) 
    if(.not. UNITOK) return

    magic='fail'
   
    if(File%iosystem%ioproc) then      
       if(File%iosystem%io_rank==0) then
!  Find a unique unit number to open the file
          do fh=12,99
             inquire (unit=fh,exist=UNITOK,opened=UNITOP)
             if (UNITOK .and. .not. UNITOP) then
                open (unit = fh,File=filename,access='direct',recl=reclength,&
                     FORM='UNFORMATTED',STATUS='OLD',err=100)
! Read the first 4 bytes and look for the CDF or HDF stamp
                read (fh,rec=1,err=101) magic
                close(fh)
                exit
             endif
          end do
          if(magic(1:3) .eq. 'CDF') then
             ! No need to do anything here
          else if(magic(2:4).eq.'HDF') then
             call piodie("ionf_mod.F90",351,'You must link with the netcdf4 ',0,&
                  'library built with hdf5 support to read this file',0,filename)
          else 
             ! The HDF identifier could be offset further into the file.
             
             open (unit = fh,file=filename,access='direct',recl=reclength,&
                  form='UNFORMATTED',STATUS='OLD',err=100)

             i=128
             eof=0
             do while(eof>=0)
                read (fh,rec=i, iostat=eof, err=101) magic

                if(magic(2:4).eq.'HDF') then
                   if(debug) print *,'Changing file type to netcdf4p'
                   File%iotype=pio_iotype_netcdf4p
                   exit
                end if
                i=i*2
             end do
             close(fh)
             if(eof<0) call piodie("ionf_mod.F90",373,'Unrecognized file format ',0,filename)             
          end if

       end if
       
       call mpi_bcast(file%iotype,1,mpi_integer, 0, file%iosystem%io_comm, mpierr)
       call CheckMPIReturn('nf_mod',mpierr)
    end if
    return
100 call piodie("ionf_mod.F90",382,'File open error ',0,filename)
101 call piodie("ionf_mod.F90",383,'File read error ',0,filename)



  end subroutine check_file_type




end module ionf_mod
