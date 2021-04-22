!>
!! @file 
!! @brief Initialization Routines for PIO
!! 
!! $Revision: 912 $
!! $LastChangedDate: 2014-01-16 19:06:39 -0700 (Thu, 16 Jan 2014) $
!<
module piolib_mod
  !--------------
  use pio_kinds
  !--------------
  use pio_types, only : file_desc_t, iosystem_desc_t, var_desc_t, io_desc_t, &
	pio_iotype_pbinary, pio_iotype_binary, pio_iotype_direct_pbinary, &
	pio_iotype_netcdf, pio_iotype_pnetcdf, pio_iotype_netcdf4p, pio_iotype_netcdf4c, &
        pio_noerr, pio_num_ost
  !--------------
  use alloc_mod
  !--------------
  use pio_support, only : piodie, debug, debugio, debugasync, checkmpireturn
  !
  use ionf_mod, only : create_nf, open_nf,close_nf, sync_nf
  use pionfread_mod, only : read_nf
  use pionfwrite_mod, only : write_nf
  use pio_mpi_utils, only : PIO_type_to_mpi_type 
  use iompi_mod
  use rearrange
  use pio_msg_mod
  use mpi    ! _EXTERNAL
  implicit none
  private
  ! !public member functions:

  public :: PIO_init,     &
       PIO_finalize,      &
       PIO_initdecomp,    &
       PIO_openfile,      &
       PIO_syncfile,      &
       PIO_createfile,    &
       PIO_closefile,     &
       PIO_setiotype,     &
       PIO_numtoread,     &
       PIO_numtowrite,    &
       PIO_setframe,      &
       PIO_advanceframe,  &
       PIO_setdebuglevel, &
       PIO_seterrorhandling, &
       PIO_get_local_array_size, &
       PIO_freedecomp,     &
       PIO_dupiodesc,     &
       PIO_getnumiotasks, &
       PIO_set_hint,      &
       PIO_getnum_OST,    &
       PIO_setnum_OST,    &
       PIO_FILE_IS_OPEN


  !eop
  !boc
  !-----------------------------------------------------------------------
  !
  !  module variables
  !
  !-----------------------------------------------------------------------
!> 
!! @defgroup PIO_openfile PIO_openfile
!< 
  interface PIO_openfile
     module procedure PIO_openfile
  end interface

!> 
!! @defgroup PIO_syncfile PIO_syncfile
!<
  interface PIO_syncfile
     module procedure syncfile
  end interface

!> 
!! @defgroup PIO_createfile PIO_createfile
!<
  interface PIO_createfile
     module procedure createfile
  end interface

!> 
!! @defgroup PIO_setframe PIO_setframe
!! @brief sets the unlimited dimension for netcdf file for record number for binary files
!<
  interface PIO_setframe
     module procedure setframe
  end interface

!> 
!! @defgroup PIO_advanceframe PIO_advanceframe
!<
  interface PIO_advanceframe
     module procedure advanceframe
  end interface

!> 
!! @defgroup PIO_closefile PIO_closefile
!<
  interface PIO_closefile
     module procedure closefile
  end interface


!> 
!! @defgroup PIO_freedecomp PIO_freedecomp
!! free memory associated with a io descriptor
!<
  interface PIO_freedecomp
     module procedure freedecomp_ios
     module procedure freedecomp_file
  end interface

!> 
!! @defgroup PIO_init PIO_init
!! initializes the pio subsystem
!<
  interface PIO_init
     module procedure init_intracom
     module procedure init_intercom
     
  end interface

!> 
!! @defgroup PIO_finalize PIO_finalize
!! Shuts down and cleans up any memory associated with the pio library.
!<
  interface PIO_finalize
     module procedure finalize
  end interface

!>
!! @defgroup PIO_initdecomp PIO_initdecomp
!! @brief PIO_initdecomp is an overload interface the models decomposition to pio.
!<


  interface PIO_initdecomp
     module procedure PIO_initdecomp_dof_i4  ! previous name: initdecomop_1dof_nf_box
     module procedure PIO_initdecomp_dof_i8  ! previous name: initdecomop_1dof_nf_box
     module procedure PIO_initdecomp_dof_i8_vdc 
     module procedure initdecomp_1dof_nf_i4
     module procedure initdecomp_1dof_nf_i8
     module procedure initdecomp_1dof_bin_i4
     module procedure initdecomp_1dof_bin_i8
     module procedure initdecomp_2dof_nf_i4
     module procedure initdecomp_2dof_nf_i8
     module procedure initdecomp_2dof_bin_i4
     module procedure initdecomp_2dof_bin_i8
     module procedure PIO_initdecomp_bc
     module procedure PIO_initdecomp_dof_dof
  end interface

!> 
!! @defgroup PIO_dupiodesc PIO_dupiodesc
!! duplicates an eisting io descriptor
!<
  interface PIO_dupiodesc
     module procedure dupiodesc
  end interface

!> 
!! @defgroup PIO_setiotype PIO_setiotype
!!  sets the io type used by pio
!<
  interface PIO_setiotype 
     module procedure setiotype
  end interface

!> 
!! @defgroup PIO_numtoread PIO_numtoread
!! returns the total number of words to read
!<
  interface PIO_numtoread
     module procedure numtoread
  end interface

!> 
!! @defgroup PIO_numtowrite PIO_numtowrite
!! returns the total number of words to write
!<
  interface PIO_numtowrite
     module procedure numtowrite
  end interface


!> 
!! @defgroup PIO_getnumiotasks PIO_getnumiotasks
!!  returns the actual number of IO-tasks used.  PIO 
!!  will reset the total number of IO-tasks if certain 
!!  conditions are meet
!<
  interface PIO_getnumiotasks
     module procedure getnumiotasks
  end interface

!> 
!!  @defgroup PIO_setdebuglevel PIO_setdebuglevel
!!  sets the level of debug information that pio will generate.
!<
  interface PIO_setdebuglevel
     module procedure setdebuglevel
  end interface

!> 
!!  @defgroup PIO_seterrorhandling PIO_seterrorhandling
!!  sets the form of error handling for pio.
!!
!! By default pio handles errors internally by printing a string
!! describing the error and calling mpi_abort.  Application
!! developers can change this behavior for calls to the underlying netcdf
!! libraries with a call to PIO_seterrorhandling. For example if a
!! developer wanted to see if an input netcdf format file contained the variable
!! 'u' they might write the following
!! @verbinclude errorhandle
!<
  interface PIO_seterrorhandling
     module procedure seterrorhandlingf
     module procedure seterrorhandlingi
  end interface

!>
!! @defgroup PIO_get_local_array_size PIO_get_local_array_size
!<

  !eoc
  !***********************************************************************

contains
!> 
!! @public 
!! @ingroup PIO_file_is_open
!! @brief This logical function indicates if a file is open.
!! @details
!! @param File @copydoc file_desc_t
!<
  logical function PIO_FILE_IS_OPEN(File)
    type(file_desc_t), intent(in) :: file
    pio_file_is_open = file%file_is_open
  end function PIO_FILE_IS_OPEN


!> 
!! @public 
!! @ingroup PIO_get_local_array_size
!! @brief This function returns the expected local size of an array associated with iodesc
!! @details
!! @param iodesc 
!! @copydoc io_desc_t
!<
  integer function PIO_get_local_array_size(iodesc)
    type(io_desc_t), intent(in) :: iodesc   
    PIO_get_local_array_size = iodesc%compsize
  end function PIO_get_local_array_size

!> 
!! @public 
!! @ingroup PIO_advanceframe
!! @brief advances the record dimension of a variable in a netcdf format file 
!!  or the block address in a binary file
!! @details
!! @param[in,out] vardesc @copybrief var_desc_t 
!<
! Subprogram not used   subroutine advanceframe(vardesc)
! Subprogram not used     type(var_desc_t), intent(inout) :: vardesc
! Subprogram not used     vardesc%rec=vardesc%rec+1
! Subprogram not used   end subroutine advanceframe

!> 
!! @public 
!! @ingroup PIO_setframe 
!! @brief sets the record dimension of a variable in a netcdf format file 
!! or the block address in a binary file
!! @details
!! @param vardesc @copydoc var_desc_t
!! @param frame   : frame number to set
!<
  subroutine setframe(vardesc,frame)
    type(var_desc_t), intent(inout) :: vardesc
    integer(kind=PIO_offset), intent(in) :: frame
    vardesc%rec=frame
  end subroutine setframe

!>  
!! @public
!! @ingroup PIO_setdebuglevel
!! @brief sets the level of debug information output to stdout by pio 
!! @details
!! @param level : default value is 0, allowed values 0-3
!<
! Subprogram not used   subroutine setdebuglevel(level)
! Subprogram not used     integer(i4), intent(in) :: level	
! Subprogram not used     if(level.eq.0) then
! Subprogram not used        debug=.false.
! Subprogram not used        debugio=.false.
! Subprogram not used        debugasync=.false.
! Subprogram not used     else if(level.eq.1) then
! Subprogram not used        debug=.true.
! Subprogram not used        debugio=.false.
! Subprogram not used        debugasync=.false.
! Subprogram not used     else if(level.eq.2) then
! Subprogram not used        debug=.false.
! Subprogram not used        debugio=.true.
! Subprogram not used        debugasync=.false.
! Subprogram not used     else if(level.eq.3) then
! Subprogram not used        debug=.true.
! Subprogram not used        debugio=.true.
! Subprogram not used        debugasync=.false.
! Subprogram not used     else if(level.eq.4) then
! Subprogram not used        debug=.false.
! Subprogram not used        debugio=.false.
! Subprogram not used        debugasync=.true.
! Subprogram not used     else if(level.eq.5) then
! Subprogram not used        debug=.true.
! Subprogram not used        debugio=.false.
! Subprogram not used        debugasync=.true.
! Subprogram not used     else if(level.ge.6) then
! Subprogram not used        debug=.true.
! Subprogram not used        debugio=.true.
! Subprogram not used        debugasync=.true.
! Subprogram not used     
! Subprogram not used     end if
! Subprogram not used   end subroutine setdebuglevel

!>
!! @ingroup PIO_seterrorhandling
!! @public
!! @brief set the pio error handling method for a file
!!
!! @param file @copydoc file_desc_t
!! @param method :
!! @copydoc PIO_error_method
!<
  subroutine seterrorhandlingf(file, method)
    type(file_desc_t), intent(inout) :: file
    integer, intent(in) :: method

    call seterrorhandlingi(file%iosystem, method)
  end subroutine seterrorhandlingf

!>
!! @ingroup PIO_seterrorhandling 
!! @public
!! @brief set the pio error handling method for the iosystem
!! @param iosystem : a defined pio system descriptor, see PIO_types
!! @param method :
!! @copydoc PIO_error_method
!<
  subroutine seterrorhandlingi(ios, method)
    use pio_types, only : pio_internal_error, pio_return_error
    use pio_msg_mod, only : pio_msg_seterrorhandling
    type(iosystem_desc_t), intent(inout) :: ios
    integer, intent(in) :: method
    integer :: msg, ierr

    if(ios%async_interface .and. .not. ios%ioproc ) then
       msg=PIO_MSG_SETERRORHANDLING
       if(ios%comp_rank==0) call mpi_send(msg, 1, mpi_integer, ios%ioroot, 1, ios%union_comm, ierr)
       call MPI_BCAST(method,1,MPI_INTEGER,ios%CompMaster, ios%my_comm , ierr)
    end if
    if(Debugasync) print *,"piolib_mod.f90",404,method
    ios%error_handling=method

    if(method > PIO_internal_error .or. method < PIO_return_error) then
       call piodie("piolib_mod.f90",408,'invalid error handling method requested')
    end if
  end subroutine seterrorhandlingi

!> 
!! @public 
!! @ingroup PIO_initdecomp
!! @brief Implements the @ref decomp_bc for PIO_initdecomp
!! @details  This provides the ability to describe a computational 
!! decomposition in PIO that has a block-cyclic form.  That is 
!! something that can be described using start and count arrays.
!! Optional parameters for this subroutine allows for the specification
!! of io decomposition using iostart and iocount arrays.  If iostart
!! and iocount arrays are not specified by the user, and rearrangement 
!! is turned on then PIO will calculate a suitable IO decomposition
!! @param iosystem @copydoc iosystem_desc_t
!! @param basepiotype @copydoc use_PIO_kinds
!! @param dims An array of the global length of each dimesion of the variable(s)
!! @param compstart The start index into the block-cyclic computational decomposition
!! @param compcount The count for the block-cyclic computational decomposition
!! @param iodesc @copydoc iodesc_generate
!! @param iostart   The start index for the block-cyclic io decomposition
!! @param iocount   The count for the block-cyclic io decomposition
!<
! Subprogram not used   subroutine PIO_initdecomp_bc(iosystem,basepiotype,dims,compstart,compcount,iodesc,iostart,iocount)
! Subprogram not used     type (iosystem_desc_t), intent(inout) :: iosystem
! Subprogram not used     integer(i4), intent(in)               :: basepiotype
! Subprogram not used     integer(i4), intent(in)               :: dims(:)
! Subprogram not used     integer (kind=PIO_OFFSET)             :: compstart(:)  
! Subprogram not used     integer (kind=PIO_OFFSET)             :: compcount(:)    
! Subprogram not used     type (IO_desc_t), intent(out)         :: iodesc
! Subprogram not used     integer (kind=PIO_OFFSET),optional    :: iostart(:)  
! Subprogram not used     integer (kind=PIO_OFFSET),optional    :: iocount(:)    
! Subprogram not used 
! Subprogram not used !    character(len=*), parameter :: '::PIO_initdecomp_bc'
! Subprogram not used 
! Subprogram not used     call piodie("piolib_mod.f90",444,'subroutine not yet implemented')
! Subprogram not used 
! Subprogram not used   end subroutine PIO_initdecomp_bc

!>
!! @public
!! @ingroup PIO_initdecomp
!! @brief Implements the @ref decomp_dof for PIO_initdecomp
!! @details  This provides the ability to describe a computational
!! decomposition in PIO using degrees of freedom method. This is  
!! a decomposition that can not be easily described using a start  
!! and count metehod (see @ref decomp_dof).  This subroutine also 
!! requires the user to specify the IO decomposition using the 
!! degree of freedom method.  This version of the subroutine 
!! is most suitable for those who want complete control over 
!! the actions of PIO.
!! @param iosystem @copydoc iosystem_desc_t
!! @param basepiotype @copydoc use_PIO_kinds
!! @param dims An array of the global length of each dimesion of the variable(s)
!! @param compdof Mapping of the storage order for the computatinal decomposition to its memory order
!! @param iodesc @copydoc iodesc_generate
!! @param iodof Mapping of the storage order for the IO decomposition its memory order
!<
! Subprogram not used   subroutine PIO_initdecomp_dof_dof(iosystem,basepiotype,dims,compdof,iodesc,iodof)
! Subprogram not used     type (iosystem_desc_t), intent(inout)          :: iosystem
! Subprogram not used     integer(i4), intent(in)                        :: basepiotype
! Subprogram not used     integer(i4), intent(in)                        :: dims(:)
! Subprogram not used     integer(i4), intent(in)                        :: compdof(:)
! Subprogram not used     type (IO_desc_t), intent(out)                   :: iodesc
! Subprogram not used     integer(i4), intent(in)                        :: iodof(:)
! Subprogram not used 
! Subprogram not used !    character(len=*), parameter :: subName=modName//'::PIO_initdecomp_dof_dof'
! Subprogram not used 
! Subprogram not used !    call piodie(subname,477,'subroutine not yet implemented')
! Subprogram not used 
! Subprogram not used   end subroutine PIO_initdecomp_dof_dof

!> 
!! @public 
!! @ingroup PIO_initdecomp
!! @brief A deprecated interface to the PIO_initdecomp method.
!! @details
!! @deprecated
!! @param iosystem : a defined pio system descriptor, see PIO_types
!! @param basepiotype : the type of variable(s) associated with this iodesc.  
!! @copydoc PIO_kinds
!! @param dims : an array of the global length of each dimesion of the variable(s)
!! @param lenblocks :
!! @param compdof : mapping of the storage order of the variable to its memory order
!! @param iodofr :
!! @param iodofw :
!! @param iodesc @copydoc iodesc_generate
!<
! Subprogram not used   subroutine initdecomp_2dof_bin_i4(iosystem,basepiotype,dims,lenblocks,compdof,iodofr,iodofw,iodesc)
! Subprogram not used     use calcdisplace_mod, only : calcdisplace
! Subprogram not used     type (iosystem_desc_t), intent(in) :: iosystem
! Subprogram not used     integer(i4), intent(in)           :: basepiotype
! Subprogram not used     integer(i4)                       :: basetype
! Subprogram not used     integer(i4), intent(in)           :: dims(:)
! Subprogram not used     integer (i4), intent(in)          :: lenblocks
! Subprogram not used     integer (i4), intent(in)          :: compdof(:)   !> global degrees of freedom for computational decomposition
! Subprogram not used     integer (i4), intent(in)          :: iodofr(:)     !> global degrees of freedom for io decomposition 
! Subprogram not used     integer (i4), intent(in)          :: iodofw(:)     !> global degrees of freedom for io decomposition 
! Subprogram not used     type (io_desc_t), intent(inout)     :: iodesc
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     call initdecomp_2dof_bin_i8(iosystem,basepiotype,dims,lenblocks,int(compdof,kind=pio_offset),int(iodofr,kind=pio_offset), &
! Subprogram not used          int(iodofw,kind=pio_offset),iodesc)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   end subroutine initdecomp_2dof_bin_i4
! Subprogram not used   subroutine initdecomp_2dof_bin_i8(iosystem,basepiotype,dims,lenblocks,compdof,iodofr,iodofw,iodesc)
! Subprogram not used     use calcdisplace_mod, only : calcdisplace
! Subprogram not used     type (iosystem_desc_t), intent(in) :: iosystem
! Subprogram not used     integer(i4), intent(in)           :: basepiotype
! Subprogram not used     integer(i4)                       :: basetype
! Subprogram not used     integer(i4), intent(in)           :: dims(:)
! Subprogram not used     integer (i4), intent(in)          :: lenblocks
! Subprogram not used     integer (kind=pio_offset), intent(in)          :: compdof(:)   !> global degrees of freedom for computational decomposition
! Subprogram not used     integer (kind=pio_offset), intent(in)          :: iodofr(:)     !> global degrees of freedom for io decomposition 
! Subprogram not used     integer (kind=pio_offset), intent(in)          :: iodofw(:)     !> global degrees of freedom for io decomposition 
! Subprogram not used     type (io_desc_t), intent(inout)     :: iodesc
! Subprogram not used 
! Subprogram not used     integer(kind=PIO_offset) :: start(1), count(1)
! Subprogram not used 
! Subprogram not used     integer (i4) :: i,ndims,n_iotasks
! Subprogram not used     integer(kind=PIO_OFFSET) glength
! Subprogram not used     logical :: userearranger
! Subprogram not used     integer (kind=pio_offset) ::  ndispr,ndispw
! Subprogram not used     integer (kind=pio_offset) :: lengthr, lengthw
! Subprogram not used     integer (kind=pio_offset), pointer :: displacer(:),displacew(:)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     nullify(iodesc%start)
! Subprogram not used     nullify(iodesc%count)
! Subprogram not used 
! Subprogram not used     basetype=PIO_type_to_mpi_type(basepiotype)
! Subprogram not used 
! Subprogram not used     !-------------------------------------------
! Subprogram not used     ! for testing purposes set the iomap
! Subprogram not used     ! (decompmap_t) to something basic for
! Subprogram not used     ! testing.
! Subprogram not used     !-------------------------------------------
! Subprogram not used     userearranger = iosystem%userearranger
! Subprogram not used 
! Subprogram not used     !---------------------
! Subprogram not used     ! number of dimensions
! Subprogram not used     !---------------------
! Subprogram not used     ndims = size(dims)
! Subprogram not used     !---------------------
! Subprogram not used     ! total global size
! Subprogram not used     !---------------------
! Subprogram not used     glength= product(int(dims,kind=PIO_OFFSET))
! Subprogram not used     if(glength > int(huge(i),kind=pio_offset)) then
! Subprogram not used        call piodie( "piolib_mod.f90",558, &
! Subprogram not used             'requested array size too large for this interface ')       
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     lengthr = size(iodofr);
! Subprogram not used     lengthw = size(iodofw)
! Subprogram not used     if(lenblocks>0) then
! Subprogram not used        ndispw=size(iodofw)/lenblocks 
! Subprogram not used        ndispr=size(iodofr)/lenblocks
! Subprogram not used     else
! Subprogram not used        ndispw=0
! Subprogram not used        ndispr=0
! Subprogram not used     end if
! Subprogram not used     call alloc_check(displacer,int(ndispr))
! Subprogram not used     call alloc_check(displacew,int(ndispw))
! Subprogram not used 
! Subprogram not used     !--------------------------------------------
! Subprogram not used     ! calculate mpi data structure displacements
! Subprogram not used     !--------------------------------------------
! Subprogram not used     !dbg    print *,'PIO_initdecomp: before call to calcdisplace'
! Subprogram not used     if(lenblocks>0) then
! Subprogram not used        call calcdisplace(lenblocks,iodofr,displacer)
! Subprogram not used        call calcdisplace(lenblocks,iodofw,displacew)
! Subprogram not used     end if
! Subprogram not used     n_iotasks = iosystem%num_iotasks
! Subprogram not used 
! Subprogram not used     iodesc%glen = glength
! Subprogram not used 
! Subprogram not used     if(debug) print *,'iam: ',iosystem%io_rank,'initdecomp: userearranger: ',userearranger
! Subprogram not used 
! Subprogram not used     !---------------------------------------------
! Subprogram not used     !  the setup for the mpi-io type information
! Subprogram not used     !---------------------------------------------
! Subprogram not used     if(iosystem%ioproc) then
! Subprogram not used        !-----------------------------------------------
! Subprogram not used        ! setup the data structure for the read operation
! Subprogram not used        !-----------------------------------------------
! Subprogram not used        iodesc%read%n_elemtype = ndispr
! Subprogram not used        iodesc%read%n_words    = iodesc%read%n_elemtype*lenblocks
! Subprogram not used        call genindexedblock(lenblocks,basetype,iodesc%read%elemtype,iodesc%read%filetype,int(displacer))
! Subprogram not used 
! Subprogram not used        !-------------------------------------------------
! Subprogram not used        ! setup the data structure for the write operation
! Subprogram not used        !-------------------------------------------------
! Subprogram not used        iodesc%write%n_elemtype = ndispw
! Subprogram not used        iodesc%write%n_words    = iodesc%write%n_elemtype*lenblocks
! Subprogram not used 
! Subprogram not used        call genindexedblock(lenblocks,basetype,iodesc%write%elemtype,iodesc%write%filetype,int(displacew))
! Subprogram not used 
! Subprogram not used        if(debug) print *,'initdecomp: at the end of subroutine'
! Subprogram not used        !       if(iodesc%read%n_elemtype == 0 .and. iodesc%write%n_elemtype == 0) iosystem%ioproc = .false.
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     deallocate(displacer,displacew)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   end subroutine initdecomp_2dof_bin_i8


!> 
!! @public 
!! @ingroup PIO_initdecomp
!! @brief A deprecated interface to the PIO_initdecomp method.
!! @details
!! @deprecated
!! @param iosystem : a defined pio system descriptor, see PIO_types
!! @param basepiotype : the type of variable(s) associated with this iodesc.  
!! @copydoc PIO_kinds
!! @param dims : an array of the global length of each dimesion of the variable(s)
!! @param lenblocks : 
!! @param compdof : mapping of the storage order of the variable to its memory order
!! @param iodofr : 
!! @param iodesc @copydoc iodesc_generate
!<
! Subprogram not used   subroutine initdecomp_1dof_bin_i8(iosystem,basepiotype,dims,lenblocks,compdof,iodofr,iodesc)
! Subprogram not used     type (iosystem_desc_t), intent(in) :: iosystem
! Subprogram not used     integer(i4), intent(in)           :: basepiotype
! Subprogram not used     integer(i4), intent(in)           :: dims(:)
! Subprogram not used     integer(i4), intent(in)          :: lenblocks
! Subprogram not used     integer(kind=pio_offset), intent(in)          :: compdof(:)   ! global degrees of freedom for computational decomposition
! Subprogram not used     integer(kind=pio_offset), intent(in)          :: iodofr(:)     ! global degrees of freedom for io decomposition 
! Subprogram not used     type (io_desc_t), intent(inout)     :: iodesc
! Subprogram not used 
! Subprogram not used     integer(kind=PIO_offset) :: start(1), count(1)
! Subprogram not used     ! these are not used in the binary interface
! Subprogram not used 
! Subprogram not used     start(1)=-1
! Subprogram not used     count(1)=-1
! Subprogram not used     call initdecomp_1dof_nf_i8(iosystem,basepiotype,dims,lenblocks,compdof,iodofr,start, count, iodesc)
! Subprogram not used   end subroutine initdecomp_1dof_bin_i8

! Subprogram not used   subroutine initdecomp_1dof_bin_i4(iosystem,basepiotype,dims,lenblocks,compdof,iodofr,iodesc)
! Subprogram not used     type (iosystem_desc_t), intent(in) :: iosystem
! Subprogram not used     integer(i4), intent(in)           :: basepiotype
! Subprogram not used     integer(i4), intent(in)           :: dims(:)
! Subprogram not used     integer (i4), intent(in)          :: lenblocks
! Subprogram not used     integer (i4), intent(in)          :: compdof(:)   ! global degrees of freedom for computational decomposition
! Subprogram not used     integer (i4), intent(in)          :: iodofr(:)     ! global degrees of freedom for io decomposition 
! Subprogram not used     type (io_desc_t), intent(inout)     :: iodesc
! Subprogram not used 
! Subprogram not used     integer(kind=PIO_offset) :: start(1), count(1)
! Subprogram not used     ! these are not used in the binary interface
! Subprogram not used 
! Subprogram not used     start(1)=-1
! Subprogram not used     count(1)=-1
! Subprogram not used     call initdecomp_1dof_nf_i8(iosystem,basepiotype,dims,lenblocks, &
! Subprogram not used          int(compdof,kind=PIO_OFFSET),int(iodofr,kind=PIO_OFFSET),start, count, iodesc)
! Subprogram not used   end subroutine initdecomp_1dof_bin_i4

!> 
!! @public 
!! @ingroup PIO_initdecomp
!! @brief A deprecated interface to the PIO_initdecomp method.
!! @details
!! @deprecated
!! @param iosystem : a defined pio system descriptor, see PIO_types
!! @param basepiotype : the type of variable(s) associated with this iodesc.
!! @copydoc PIO_kinds
!! @param dims : an array of the global length of each dimesion of the variable(s)
!! @param lenblocks : 
!! @param compdof : mapping of the storage order of the variable to its memory order
!! @param iodofr : 
!! @param iodofw :
!! @param start : used with count to give a block description of the shape of the data
!! @param count : 
!! @param iodesc @copydoc iodesc_generate
!<
! Subprogram not used   subroutine initdecomp_2dof_nf_i4(iosystem,basepiotype,dims,lenblocks,compdof,iodofr,iodofw,start, count, iodesc)
! Subprogram not used     type (iosystem_desc_t), intent(in) :: iosystem
! Subprogram not used     integer(i4), intent(in)           :: basepiotype
! Subprogram not used     integer(i4), intent(in)           :: dims(:)
! Subprogram not used     integer (i4), intent(in)          :: lenblocks
! Subprogram not used     integer (i4), intent(in)          :: compdof(:)   ! global degrees of freedom for computational decomposition
! Subprogram not used     integer (i4), intent(in)          :: iodofr(:)     ! global degrees of freedom for io decomposition 
! Subprogram not used     integer (i4), intent(in)          :: iodofw(:)     ! global degrees of freedom for io decomposition 
! Subprogram not used 
! Subprogram not used     type (io_desc_t), intent(inout)     :: iodesc
! Subprogram not used 
! Subprogram not used     integer(kind=PIO_offset), intent(in) :: start(:), count(:)
! Subprogram not used     type (io_desc_t) :: tmp
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     call pio_initdecomp(iosystem, basepiotype,dims,lenblocks,int(compdof,kind=PIO_OFFSET),int(iodofr,kind=PIO_OFFSET), &
! Subprogram not used          int(iodofw,kind=PIO_OFFSET),start,count,iodesc)
! Subprogram not used 
! Subprogram not used   end subroutine initdecomp_2dof_nf_i4

! Subprogram not used   subroutine initdecomp_2dof_nf_i8(iosystem,basepiotype,dims,lenblocks,compdof,iodofr,iodofw,start, count, iodesc)
! Subprogram not used     type (iosystem_desc_t), intent(in) :: iosystem
! Subprogram not used     integer(i4), intent(in)           :: basepiotype
! Subprogram not used     integer(i4), intent(in)           :: dims(:)
! Subprogram not used     integer (i4), intent(in)          :: lenblocks
! Subprogram not used     integer (kind=pio_offset), intent(in)          :: compdof(:)   ! global degrees of freedom for computational decomposition
! Subprogram not used     integer (kind=pio_offset), intent(in)          :: iodofr(:)     ! global degrees of freedom for io decomposition 
! Subprogram not used     integer (kind=pio_offset), intent(in)          :: iodofw(:)     ! global degrees of freedom for io decomposition 
! Subprogram not used 
! Subprogram not used     type (io_desc_t), intent(inout)     :: iodesc
! Subprogram not used 
! Subprogram not used     integer(kind=PIO_offset), intent(in) :: start(:), count(:)
! Subprogram not used     type (io_desc_t) :: tmp
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     call initdecomp_1dof_nf_i8(iosystem, basepiotype, dims, lenblocks, compdof, iodofr, start, count, iodesc)
! Subprogram not used 
! Subprogram not used     call initdecomp_1dof_nf_i8(iosystem, basepiotype, dims, lenblocks, compdof, iodofw, start, count, tmp)
! Subprogram not used 
! Subprogram not used     call dupiodesc2(iodesc%write,tmp%write)
! Subprogram not used 
! Subprogram not used     if(debug) then
! Subprogram not used        print *, "piolib_mod.f90",729,iodesc%read%filetype,iodesc%read%elemtype,&
! Subprogram not used             iodesc%read%n_elemtype,iodesc%read%n_words   
! Subprogram not used        print *, "piolib_mod.f90",731,iodesc%write%filetype,iodesc%write%elemtype,&
! Subprogram not used             iodesc%write%n_elemtype,iodesc%write%n_words
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used   end subroutine initdecomp_2dof_nf_i8

!> 
!! @public 
!! @ingroup PIO_initdecomp
!! @brief A deprecated interface to the PIO_initdecomp method.
!! @details
!! @deprecated
!! @param iosystem : a defined PIO system descriptor, see pio_types
!! @param basepiotype : The type of variable(s) associated with this iodesc.  
!! @copydoc PIO_kinds
!! @param dims : an array of the global length of each dimesion of the variable(s)
!! @param lenblocks : 
!! @param compdof : mapping of the storage order of the variable to its memory order
!! @param iodof : 
!! @param start :
!! @param count :
!! @param iodesc @copydoc iodesc_generate
!<
! Subprogram not used   subroutine initdecomp_1dof_nf_i4(iosystem,basepiotype,dims,lenblocks,compdof,iodof,start, count, iodesc)
! Subprogram not used     use calcdisplace_mod, only : calcdisplace
! Subprogram not used     type (iosystem_desc_t), intent(in) :: iosystem
! Subprogram not used     integer(i4), intent(in)           :: basepiotype
! Subprogram not used     integer(i4), intent(in)           :: dims(:)
! Subprogram not used     integer (i4), intent(in) :: lenblocks
! Subprogram not used     integer (i4), intent(in)          :: compdof(:)   ! global degrees of freedom for computational decomposition
! Subprogram not used     integer (i4), intent(in)          :: iodof(:)     ! global degrees of freedom for io decomposition 
! Subprogram not used     type (io_desc_t), intent(inout)     :: iodesc
! Subprogram not used     integer :: piotype	
! Subprogram not used     integer(kind=PIO_offset), intent(in) :: start(:), count(:)
! Subprogram not used 
! Subprogram not used     call initdecomp_1dof_nf_i8(iosystem, basepiotype,dims,lenblocks,int(compdof,kind=pio_offset),int(iodof,kind=pio_offset),&
! Subprogram not used          start,count,iodesc)
! Subprogram not used 
! Subprogram not used   end subroutine initdecomp_1dof_nf_i4
! Subprogram not used   subroutine initdecomp_1dof_nf_i8(iosystem,basepiotype,dims,lenblocks,compdof,iodof,start, count, iodesc)
! Subprogram not used     use calcdisplace_mod, only : calcdisplace
! Subprogram not used     type (iosystem_desc_t), intent(in) :: iosystem
! Subprogram not used     integer(i4), intent(in)           :: basepiotype
! Subprogram not used     integer(i4), intent(in)           :: dims(:)
! Subprogram not used     integer (i4), intent(in) :: lenblocks
! Subprogram not used     integer (kind=pio_offset), intent(in)          :: compdof(:)   ! global degrees of freedom for computational decomposition
! Subprogram not used     integer (kind=pio_offset), intent(in)          :: iodof(:)     ! global degrees of freedom for io decomposition 
! Subprogram not used     type (io_desc_t), intent(inout)     :: iodesc
! Subprogram not used     integer :: piotype
! Subprogram not used     integer(kind=PIO_offset), intent(in) :: start(:), count(:)
! Subprogram not used 
! Subprogram not used     integer(i4) :: length,n_iotasks
! Subprogram not used     integer(i4) :: ndims
! Subprogram not used 
! Subprogram not used     integer (kind=pio_offset), pointer :: displace(:)  ! the displacements for the mpi data structure (read)
! Subprogram not used 
! Subprogram not used     integer(i4) :: prev
! Subprogram not used     integer(kind=PIO_OFFSET) :: glength    ! global length in words
! Subprogram not used     integer(i4) :: ii,i,dis,ierr
! Subprogram not used     integer(i4),pointer, dimension(:) :: blocklen,disp
! Subprogram not used     logical(log_kind) ::  userearranger
! Subprogram not used     logical, parameter :: check = .true.
! Subprogram not used     integer(kind=pio_offset) :: ndisp
! Subprogram not used     nullify(iodesc%start)
! Subprogram not used     nullify(iodesc%count)
! Subprogram not used 
! Subprogram not used     piotype=PIO_type_to_mpi_type(basepiotype)
! Subprogram not used 
! Subprogram not used     !-------------------------------------------
! Subprogram not used     ! for testing purposes set the iomap
! Subprogram not used     ! (decompmap_t) to something basic for
! Subprogram not used     ! testing.
! Subprogram not used     !-------------------------------------------
! Subprogram not used     userearranger = iosystem%userearranger
! Subprogram not used     !---------------------
! Subprogram not used     ! number of dimensions
! Subprogram not used     !---------------------
! Subprogram not used     ndims = size(dims)
! Subprogram not used     !---------------------
! Subprogram not used     ! total global size
! Subprogram not used     !---------------------
! Subprogram not used     glength= product(int(dims,kind=PIO_OFFSET))
! Subprogram not used     if(glength > huge(ndisp)) then
! Subprogram not used        print *,"<stdin>",827,dims,glength
! Subprogram not used        call piodie( "piolib_mod.f90",828, &
! Subprogram not used             'requested array size too large for this interface ')       
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if(lenblocks>0) then
! Subprogram not used        ndisp=size(iodof)/lenblocks
! Subprogram not used     else
! Subprogram not used        ndisp=0
! Subprogram not used     end if
! Subprogram not used     call alloc_check(displace,int(ndisp))
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     call alloc_check(iodesc%start,ndims)
! Subprogram not used     call alloc_check(iodesc%count,ndims)
! Subprogram not used     iodesc%start(1:size(start)) = start(:)
! Subprogram not used     iodesc%count(1:size(count)) = count(:)
! Subprogram not used     !--------------------------------------------
! Subprogram not used     ! calculate mpi data structure displacements 
! Subprogram not used     !--------------------------------------------
! Subprogram not used     if(lenblocks>0) then
! Subprogram not used        if(debug) print *,'PIO_initdecomp: calcdisplace',ndisp,size(iodof),lenblocks
! Subprogram not used        call calcdisplace(lenblocks,iodof,displace)
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     n_iotasks = iosystem%num_iotasks
! Subprogram not used     length = size(iodof)
! Subprogram not used     !
! Subprogram not used     !   this facilitates the use of seperate read and write descripters. 
! Subprogram not used     !
! Subprogram not used     iodesc%iomap%start  = iosystem%io_rank*length
! Subprogram not used     iodesc%iomap%length = length
! Subprogram not used     iodesc%glen = glength
! Subprogram not used 
! Subprogram not used     if(debug) print *,'iam: ',iosystem%io_rank,'initdecomp: userearranger: ',userearranger, glength
! Subprogram not used     if(userearranger) then 
! Subprogram not used              call piodie( "piolib_mod.f90",877, &
! Subprogram not used                   'this interface does not use rearranger')
! Subprogram not used        
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     !---------------------------------------------
! Subprogram not used     !  the setup for the mpi-io type information 
! Subprogram not used     !---------------------------------------------
! Subprogram not used     if(iosystem%ioproc) then 
! Subprogram not used        !-----------------------------------------------
! Subprogram not used        ! setup the data structure for the io operation 
! Subprogram not used        !-----------------------------------------------
! Subprogram not used        iodesc%write%n_elemtype = ndisp
! Subprogram not used        iodesc%write%n_words    = iodesc%write%n_elemtype*lenblocks
! Subprogram not used 
! Subprogram not used        call genindexedblock(lenblocks,piotype,iodesc%write%elemtype,iodesc%write%filetype,int(displace))
! Subprogram not used 
! Subprogram not used        
! Subprogram not used !       call gensubarray(dims,piotype,iodesc,iodesc%write)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used        if(debug) print *,'initdecomp: at the end of subroutine',iodesc%write%n_elemtype,iodesc%write%n_words
! Subprogram not used     endif
! Subprogram not used     call dupiodesc2(iodesc%write,iodesc%read)
! Subprogram not used     if(debug) then
! Subprogram not used        print *, "piolib_mod.f90",918,iodesc%read%filetype,iodesc%read%elemtype,&
! Subprogram not used             iodesc%read%n_elemtype,iodesc%read%n_words   
! Subprogram not used        print *, "piolib_mod.f90",920,iodesc%write%filetype,iodesc%write%elemtype,&
! Subprogram not used             iodesc%write%n_elemtype,iodesc%write%n_words
! Subprogram not used     end if
! Subprogram not used     call dealloc_check(displace)
! Subprogram not used 
! Subprogram not used   end subroutine initdecomp_1dof_nf_i8

!>
!! @public
!! @ingroup PIO_initdecomp
!! @brief Implements the @ref decomp_dof for PIO_initdecomp (previous name: \b initdecomp_1dof_nf_box)
!! @details  This provides the ability to describe a computational
!! decomposition in PIO using degrees of freedom method. This is
!! a decomposition that can not be easily described using a start
!! and count method (see @ref decomp_dof).
!! Optional parameters for this subroutine allows for the specififcation of
!! io decomposition using iostart and iocount arrays.  If iostart
!! and iocount arrays are not specified by the user, and rearrangement
!! is turned on then PIO will calculate an suitable IO decomposition.
!! Note that this subroutine was previously called \em initdecomp_1dof_nf_box
!! @param iosystem @copydoc iosystem_desc_t
!! @param basepiotype @copydoc use_PIO_kinds
!! @param dims An array of the global length of each dimesion of the variable(s)
!! @param compdof Mapping of the storage order for the computational decomposition to its memory order
!! @param iodesc @copydoc iodesc_generate
!! @param iostart   The start index for the block-cyclic io decomposition
!! @param iocount   The count for the block-cyclic io decomposition
!<
  subroutine PIO_initdecomp_dof_i4(iosystem,basepiotype,dims,compdof, iodesc, iostart, iocount, num_ts, bsize)
    type (iosystem_desc_t), intent(inout) :: iosystem
    integer(i4), intent(in)           :: basepiotype
    integer(i4), intent(in)          :: compdof(:)   ! global degrees of freedom for computational decomposition
    integer (kind=PIO_offset), optional :: iostart(:), iocount(:)
    type (io_desc_t), intent(inout)     :: iodesc
    integer(kind=PIO_OFFSET), pointer :: internal_compdof(:)
    integer(i4), intent(in)           :: dims(:)
    !vdf optionals
    integer(i4), intent(in), optional:: num_ts, bsize(3)
    allocate(internal_compdof(size(compdof)))
    internal_compdof = int(compdof,kind=pio_offset)
    
    if(present(iostart) .and. present(iocount) ) then
       call pio_initdecomp_dof_i8(iosystem, basepiotype, dims, internal_compdof, iodesc, iostart, iocount)
    else 
       call pio_initdecomp_dof_i8(iosystem, basepiotype, dims, internal_compdof, iodesc)
    endif
    deallocate(internal_compdof)

  end subroutine PIO_initdecomp_dof_i4


  subroutine PIO_initdecomp_dof_i8(iosystem,basepiotype,dims,compdof, iodesc, iostart, iocount)
    use calcdisplace_mod, only : calcdisplace_box
    use calcdecomp, only : calcstartandcount
    type (iosystem_desc_t), intent(inout) :: iosystem
    integer(i4), intent(in)           :: basepiotype
    integer(i4), intent(in)           :: dims(:)
    integer (kind=pio_offset), intent(in)          :: compdof(:)   ! global degrees of freedom for computational decomposition
    integer (kind=PIO_offset), optional :: iostart(:), iocount(:)
    type (io_desc_t), intent(inout)     :: iodesc

    integer(i4) :: length,n_iotasks
    integer(i4) :: ndims
    integer (i4)                       :: lenblocks
    integer(i4)                       ::  piotype

    integer(i4), pointer :: displace(:)  ! the displacements for the mpi data structure (read)

    integer(i4) :: prev
    integer(kind=PIO_OFFSET) :: glength    ! global length in words
    integer(i4) :: ii,i,dis,ierr
    integer(i4),pointer, dimension(:) :: blocklen,disp
    logical(log_kind) ::  userearranger
    logical, parameter :: check = .true.
    integer(kind=pio_offset) :: ndisp
    integer(i4) :: iosize               ! rml
    integer(i4) :: msg
    integer(i4), allocatable :: lstart(:),lcount(:)
    logical :: is_async=.false.
    integer ierror, dsize

    nullify(displace)

    if(iosystem%async_interface .and. .not. iosystem%ioproc) then
       msg = PIO_MSG_INITDECOMP_DOF
       is_async=.true.
       if(DebugAsync) print*,"piolib_mod.f90",1022, iosystem%ioranks
       if(iosystem%comp_rank==0) then
          call mpi_send(msg, 1, mpi_integer, iosystem%ioroot, 1, iosystem%union_comm, ierr)
       end if
       if(DebugAsync) print*,"piolib_mod.f90",1026, ierr, iosystem%ioroot, iosystem%comp_rank

       call mpi_bcast(basepiotype, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
       if(DebugAsync) print*,"piolib_mod.f90",1029
       dsize = size(dims)
       call mpi_bcast(dsize, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
       call mpi_bcast(dims, size(dims), mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)

       if(DebugAsync) print*,"piolib_mod.f90",1034
       call mpi_bcast(iodesc%async_id, 1, mpi_integer, iosystem%iomaster, iosystem%intercomm, ierr)  
       if(DebugAsync) print*,"piolib_mod.f90",1036, iodesc%async_id
    endif

    if(minval(dims)<=0) then
       print *,"piolib_mod.f90",1040,dims
       call piodie("piolib_mod.f90",1041,'bad value in dims argument')
    end if

    if (iosystem%comp_rank == 0 .and. debug) &
         print *,iosystem%comp_rank,': invoking PIO_initdecomp_dof'

    if(DebugAsync) print*,"piolib_mod.f90",1047
    piotype=PIO_type_to_mpi_type(basepiotype)

    !-------------------------------------------
    ! for testing purposes set the iomap
    ! (decompmap_t) to something basic for
    ! testing.
    !-------------------------------------------

    userearranger = iosystem%userearranger
    !---------------------
    ! number of dimensions
    !---------------------
    ndims = size(dims)
    !---------------------
    ! total global size
    !---------------------
    glength= product(int(dims,kind=PIO_OFFSET))
    if(glength > huge(int(i,kind=pio_offset))) then !not sure if this works, glength is pio_offset, if its > pio_offset range then 
       call piodie( "piolib_mod.f90",1073, & !it will simply wrap around rather than be > max_int(pio_offset)
            'requested array size too large for this interface ') !might be better to use a temp 8 byte int to store results
                                                                  !of dims product and compare to the maxint(pio_offset)       
    endif

       

    ! remember iocount() is only defined on io procs
    call alloc_check(iodesc%start,ndims)
    call alloc_check(iodesc%count,ndims)
    iodesc%basetype=piotype
       
    iodesc%compsize=size(compdof)
       
    iodesc%start=0
    iodesc%count=0

    if(debug) print*,"piolib_mod.f90",1090, 'before calcstartandcount: ', iosystem%num_tasks, iosystem%num_iotasks, &
         iosystem%io_rank, iosystem%io_comm, iosystem%ioranks

    if (iosystem%ioproc) then
       if(present(iostart) .and. present(iocount)) then
          iodesc%start = iostart
          iodesc%count = iocount
       else if(present(iostart) .or. present(iocount)) then
          call piodie( "piolib_mod.f90",1098, &
               'both optional parameters start and count must be provided')
       else	       
          call calcstartandcount(basepiotype, ndims, dims, iosystem%num_iotasks, iosystem%io_rank,&
               iodesc%start, iodesc%count,iosystem%num_aiotasks)
       endif
       iosize=1
       do i=1,ndims
          iosize=iosize*iodesc%count(i)
       end do
       call mpi_allreduce(iosize, iodesc%maxiobuflen, 1, mpi_integer, mpi_max, iosystem%io_comm, ierr)
       call checkmpireturn('mpi_allreduce in initdecomp',ierr)

       lenblocks=1
       do i=1,ndims
          if(iodesc%count(i) == dims(i)) then
             lenblocks=lenblocks*iodesc%count(i)
          else
             exit
          endif
       enddo
       if(lenblocks==1) lenblocks=iodesc%count(1)

       if(lenblocks>0) then
          ndisp=iosize/lenblocks
       else
          ndisp=0
       end if
       call alloc_check(displace,int(ndisp))

       if(debug) print *,'IAM: ',iosystem%comp_rank,' after getiostartandcount: count is: ',iodesc%count,&
            ' lenblocks =',lenblocks,' ndisp=',ndisp

       if(debug) print *,'IAM: ',iosystem%comp_rank,' after getiostartandcount, num_aiotasks is: ', iosystem%num_aiotasks       
       !--------------------------------------------
       ! calculate mpi data structure displacements 
       !--------------------------------------------
      
       if(debug) print *,'PIO_initdecomp: calcdisplace', &
            ndisp,iosize,lenblocks, iodesc%start, iodesc%count
       call calcdisplace_box(dims,lenblocks,iodesc%start,iodesc%count,ndims,displace)
          
       n_iotasks = iosystem%num_iotasks
       length = iosize                      ! rml

       !
       !   this facilitates the use of seperate read and write descripters. 
       !

       iodesc%iomap%start  = iosystem%io_rank*length
       iodesc%iomap%length = length
       iodesc%glen = glength
    endif
    if(DebugAsync) print*,"piolib_mod.f90",1151

    if(debug) print *,"piolib_mod.f90",1160,'iam: ',iosystem%io_rank, &
         'initdecomp: userearranger: ',userearranger, glength

    if(userearranger) then 
       call MPI_BCAST(iosystem%num_aiotasks,1,mpi_integer,iosystem%iomaster,&
            iosystem%my_comm,ierr)
       call rearrange_create( iosystem,compdof,dims,ndims,iodesc)
    endif

    if(DebugAsync) print*,"piolib_mod.f90",1169

    !---------------------------------------------
    !  the setup for the mpi-io type information 
    !---------------------------------------------
    if(iosystem%ioproc) then 
       !-----------------------------------------------
       ! setup the data structure for the io operation 
       !-----------------------------------------------
       call gensubarray(dims,piotype,iodesc,iodesc%write)
       
       if(debug) print *,"piolib_mod.f90",1187,iodesc%write%n_elemtype, &
        iodesc%write%n_words,iodesc%write%elemtype,iodesc%write%filetype, lenblocks

    else
       iodesc%write%n_elemtype=0
       iodesc%write%n_words=0
       iodesc%write%elemtype = mpi_datatype_null
       iodesc%write%filetype = mpi_datatype_null
    endif


    call dupiodesc2(iodesc%write,iodesc%read)
    

    if (associated(displace)) then
       call dealloc_check(displace)
    endif


  end subroutine PIO_initdecomp_dof_i8

! Subprogram not used   subroutine PIO_initdecomp_dof_i8_vdc(iosystem,dims,compdof, iodesc, num_ts, bsize)
! Subprogram not used     use calcdisplace_mod, only : calcdisplace_box
! Subprogram not used     use calcdecomp, only : calcstartandcount
! Subprogram not used     use pio_types, only : pio_real
! Subprogram not used     type (iosystem_desc_t), intent(inout) :: iosystem
! Subprogram not used     integer(i4), intent(in)           :: dims(:)
! Subprogram not used     integer (kind=pio_offset), intent(in)          :: compdof(:)   ! global degrees of freedom for computational decomposition
! Subprogram not used 
! Subprogram not used     type (io_desc_t), intent(inout)     :: iodesc
! Subprogram not used     !vdc args
! Subprogram not used     integer(i4), intent(in) :: num_ts
! Subprogram not used     integer(i4), intent(in), optional:: bsize(3)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     integer(i4) :: length,n_iotasks
! Subprogram not used     integer(i4) :: ndims
! Subprogram not used     integer (i4)                       :: lenblocks
! Subprogram not used     integer(i4)                       ::  piotype
! Subprogram not used     integer(i4), pointer :: displace(:)  ! the displacements for the mpi data structure (read)
! Subprogram not used 
! Subprogram not used     integer(i4) :: prev
! Subprogram not used     integer(kind=PIO_OFFSET) :: glength    ! global length in words
! Subprogram not used     integer(i4) :: ii,i,dis,ierr
! Subprogram not used     integer(i4),pointer, dimension(:) :: blocklen,disp
! Subprogram not used     logical(log_kind) ::  userearranger
! Subprogram not used     logical, parameter :: check = .true.
! Subprogram not used     integer(kind=pio_offset) :: ndisp
! Subprogram not used     integer(i4) :: iosize               ! rml
! Subprogram not used     integer(i4) :: msg, dsize
! Subprogram not used     logical :: is_async=.false.
! Subprogram not used 
! Subprogram not used     integer ierror
! Subprogram not used 
! Subprogram not used     nullify(iodesc%start)
! Subprogram not used     nullify(iodesc%count)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     if(iosystem%async_interface .and. .not. iosystem%ioproc) then
! Subprogram not used        msg = PIO_MSG_INITDECOMP_DOF
! Subprogram not used        is_async=.true.
! Subprogram not used        if(DebugAsync) print*,"piolib_mod.f90",1271, iosystem%ioranks
! Subprogram not used        if(iosystem%comp_rank==0) then
! Subprogram not used           call mpi_send(msg, 1, mpi_integer, iosystem%ioroot, 1, iosystem%union_comm, ierr)
! Subprogram not used        end if
! Subprogram not used        if(DebugAsync) print*,"piolib_mod.f90",1275, ierr, iosystem%ioroot, iosystem%comp_rank
! Subprogram not used 
! Subprogram not used !       call mpi_bcast(basepiotype, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
! Subprogram not used !       if(DebugAsync) print*,"piolib_mod.f90",1278
! Subprogram not used        dsize = size(dims)
! Subprogram not used        call mpi_bcast(dsize, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
! Subprogram not used        call mpi_bcast(dims, dsize, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
! Subprogram not used 
! Subprogram not used        if(DebugAsync) print*,"piolib_mod.f90",1283
! Subprogram not used        call mpi_bcast(iodesc%async_id, 1, mpi_integer, iosystem%iomaster, iosystem%intercomm, ierr)  
! Subprogram not used        if(DebugAsync) print*,"piolib_mod.f90",1285, iodesc%async_id
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if(minval(dims)<=0) then
! Subprogram not used        print *,"piolib_mod.f90",1289,dims
! Subprogram not used        call piodie("piolib_mod.f90",1290,'bad value in dims argument')
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     if (iosystem%comp_rank == 0 .and. debug) &
! Subprogram not used          print *,iosystem%comp_rank,': invoking PIO_initdecomp_dof'
! Subprogram not used 
! Subprogram not used     if(DebugAsync) print*,"piolib_mod.f90",1296
! Subprogram not used     piotype=MPI_REAL4
! Subprogram not used 
! Subprogram not used     !-------------------------------------------
! Subprogram not used     ! for testing purposes set the iomap
! Subprogram not used     ! (decompmap_t) to something basic for
! Subprogram not used     ! testing.
! Subprogram not used     !-------------------------------------------
! Subprogram not used 
! Subprogram not used     userearranger = iosystem%userearranger
! Subprogram not used     !---------------------
! Subprogram not used     ! number of dimensions
! Subprogram not used     !---------------------
! Subprogram not used     ndims = size(dims)
! Subprogram not used     !---------------------
! Subprogram not used     ! total global size
! Subprogram not used     !---------------------
! Subprogram not used     glength= product(int(dims,kind=PIO_OFFSET))
! Subprogram not used     if(glength > huge(int(i,kind=pio_offset))) then !not sure if this works, glength is pio_offset, if its > pio_offset range then 
! Subprogram not used        call piodie( "piolib_mod.f90",1322, & !it will simply wrap around rather than be > max_int(pio_offset)
! Subprogram not used             'requested array size too large for this interface ') !might be better to use a temp 8 byte int to store results
! Subprogram not used                                                                   !of dims product and compare to the maxint(pio_offset)       
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used        
! Subprogram not used 
! Subprogram not used     ! remember iocount() is only defined on io procs
! Subprogram not used     call alloc_check(iodesc%start,ndims)
! Subprogram not used     call alloc_check(iodesc%count,ndims)
! Subprogram not used     iodesc%basetype=piotype
! Subprogram not used        
! Subprogram not used     iodesc%compsize=size(compdof)
! Subprogram not used        
! Subprogram not used     iodesc%start=0
! Subprogram not used     iodesc%count=0
! Subprogram not used 
! Subprogram not used     if (iosystem%ioproc) then
! Subprogram not used 
! Subprogram not used        iosize=1
! Subprogram not used        do i=1,ndims
! Subprogram not used           iosize=iosize*iodesc%count(i)
! Subprogram not used        end do
! Subprogram not used        call mpi_allreduce(iosize, iodesc%maxiobuflen, 1, mpi_integer, mpi_max, iosystem%io_comm, ierr)
! Subprogram not used        call checkmpireturn('mpi_allreduce in initdecomp',ierr)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used        iodesc%iomap%start  = iosystem%io_rank*iosize
! Subprogram not used        iodesc%iomap%length = iosize
! Subprogram not used        iodesc%glen = glength
! Subprogram not used     endif
! Subprogram not used     if(DebugAsync) print*,"piolib_mod.f90",1373
! Subprogram not used 
! Subprogram not used     if(userearranger) then 
! Subprogram not used        call MPI_BCAST(iosystem%num_aiotasks,1,mpi_integer,iosystem%iomaster,&
! Subprogram not used             iosystem%my_comm,ierr)
! Subprogram not used        call rearrange_create( iosystem,compdof,dims,ndims,iodesc)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     iodesc%write%n_elemtype=0
! Subprogram not used     iodesc%write%n_words=0
! Subprogram not used     iodesc%write%elemtype = mpi_datatype_null
! Subprogram not used     iodesc%write%filetype = mpi_datatype_null
! Subprogram not used 
! Subprogram not used     call dupiodesc2(iodesc%write,iodesc%read)
! Subprogram not used     
! Subprogram not used 
! Subprogram not used   end subroutine PIO_initdecomp_dof_i8_vdc



  !************************************
  ! dupiodesc2
  !

  subroutine dupiodesc2(src, dest)
    use pio_types, only : io_desc2_t
    type(io_desc2_t), intent(in) :: src
    type(io_desc2_t), intent(out) :: dest

    dest%filetype = src%filetype
    dest%elemtype = src%elemtype
    dest%n_elemtype = src%n_elemtype
    dest%n_words = src%n_words
  end subroutine dupiodesc2



  !************************************
  ! genindexedblock
  !
  ! given input lenblocks, basetype, and displacement
  ! create two mpi types: 
  !   elemtype - a single block of basetype repeated lenblocks times
  !   filetype - elemtype repeated at each entry in displacement()
  !              (i.e. size(displacement) entries)
  !


! Subprogram not used   subroutine genindexedblock(lenblocks,basetype,elemtype,filetype,displace)
! Subprogram not used     use pio_types, only : pio_double, pio_int, pio_real, pio_char
! Subprogram not used     integer(i4), intent(in) :: lenblocks     ! length of blocks
! Subprogram not used     integer(i4), intent(in) :: basetype      ! base mpi type 
! Subprogram not used     integer(i4), intent(inout) :: elemtype   ! elementary mpi type
! Subprogram not used     integer(i4), intent(inout) :: filetype   ! file mpi type 
! Subprogram not used     integer(i4), intent(in) :: displace(:)   ! mpi displacement in the array
! Subprogram not used 
! Subprogram not used     integer(i4) :: numblocks,i,ierr, prev
! Subprogram not used 
! Subprogram not used     logical, parameter :: check = .true.
! Subprogram not used 
! Subprogram not used     integer:: nints, nadds, ndtypes, comb, lbasetype
! Subprogram not used 
! Subprogram not used     numblocks = size(displace)
! Subprogram not used 
! Subprogram not used     !tcx - allow empty displace array
! Subprogram not used     if (numblocks > 0) then
! Subprogram not used        prev = displace(1)
! Subprogram not used        do i=2,numblocks
! Subprogram not used           if(prev > displace(i)) then
! Subprogram not used              print *,'genindexedblock: error detected: non-monotonic increasing displace detected!'
! Subprogram not used           endif
! Subprogram not used           prev = displace(i)
! Subprogram not used        enddo
! Subprogram not used 
! Subprogram not used     endif
! Subprogram not used     select case(basetype)
! Subprogram not used     case (PIO_double)
! Subprogram not used        lbasetype=mpi_real8
! Subprogram not used     case (PIO_real  )
! Subprogram not used        lbasetype=mpi_real4
! Subprogram not used     case (PIO_int)
! Subprogram not used        lbasetype=mpi_integer
! Subprogram not used     case (PIO_char)
! Subprogram not used        lbasetype=mpi_character
! Subprogram not used     case default
! Subprogram not used        lbasetype=basetype
! Subprogram not used     end select
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     if(lenblocks<1) then
! Subprogram not used        elemtype = lbasetype
! Subprogram not used        filetype = lbasetype
! Subprogram not used     else
! Subprogram not used        elemtype = mpi_datatype_null
! Subprogram not used        filetype = mpi_datatype_null
! Subprogram not used        call mpi_type_contiguous(lenblocks,lbasetype,elemtype,ierr)
! Subprogram not used        if(check) call checkmpireturn('genindexedblock: after call to type_contiguous: ',ierr)
! Subprogram not used        call mpi_type_commit(elemtype,ierr)
! Subprogram not used        if(check) call checkmpireturn('genindexedblock: after call to type_commit: ',ierr)
! Subprogram not used        if(numblocks>0) then
! Subprogram not used          call mpi_type_create_indexed_block(numblocks,1,displace,elemtype,filetype,ierr)
! Subprogram not used          if(check) call checkmpireturn('genindexedblock: after call to type_create_indexed_block: ',ierr)
! Subprogram not used          call mpi_type_commit(filetype,ierr)
! Subprogram not used          if(check) call checkmpireturn('genindexedblock: after call to type_commit: ',ierr)
! Subprogram not used          if(debug) then
! Subprogram not used             call mpi_type_get_envelope(filetype, nints, nadds, ndtypes, comb, ierr)
! Subprogram not used             print *,"<stdin>",1498,nints,nadds,ndtypes,comb,ierr
! Subprogram not used          endif
! Subprogram not used        endif
! Subprogram not used 
! Subprogram not used     end if
! Subprogram not used     ! _MPISERIAL
! Subprogram not used 
! Subprogram not used   end subroutine genindexedblock

  subroutine gensubarray(gdims,mpidatatype, iodesc, iodesc2)
    use pio_types, only : io_desc2_t, io_desc_t
    implicit none

    integer, intent(in) :: gdims(:)
    integer, intent(in) :: mpidatatype
    type(IO_desc_t), intent(in) :: iodesc
    type(IO_desc2_t), intent(inout) :: iodesc2
    
    integer :: ndims, ierr
    integer, allocatable :: lstart(:), lcount(:)

    ndims = size(gdims)
    if(sum(iodesc%count)>0) then
       allocate(lstart(ndims),lcount(ndims))
       lstart = 0
       lcount = int(iodesc%count)
       iodesc2%n_elemtype = 1
       iodesc2%n_words = product(lcount)
       call mpi_type_contiguous(iodesc2%n_words,mpidatatype,iodesc2%elemtype,ierr)
       call checkmpireturn('mpi_type_create_subarray in initdecomp',ierr)
       call mpi_type_commit(iodesc2%elemtype,ierr)
       call checkmpireturn('mpi_type_commit in initdecomp',ierr)

       iodesc2%filetype=mpi_datatype_null
    else
       iodesc2%elemtype=mpidatatype
       iodesc2%filetype=mpidatatype          
       iodesc2%n_elemtype = 0
       iodesc2%n_words = 0
    endif
    



  end subroutine gensubarray



!> 
!! @public
!! @ingroup PIO_init
!! @brief initialize the pio subsystem. 
!! @details  This is a collective call.  Input parameters are read on comp_rank=0
!!   values on other tasks are ignored.  This variation of PIO_init locates the IO tasks on a subset 
!!   of the compute tasks.
!! @param comp_rank mpi rank of each participating task,
!! @param comp_comm the mpi communicator which defines the collective.
!! @param num_iotasks the number of iotasks to define.
!! @param num_aggregator the mpi aggregator count
!! @param stride the stride in the mpi rank between io tasks.
!! @param rearr @copydoc PIO_rearr_method
!! @param iosystem a derived type which can be used in subsequent pio operations (defined in PIO_types).
!! @param base @em optional argument can be used to offset the first io task - default base is task 1.
!<
  subroutine init_intracom(comp_rank, comp_comm, num_iotasks, num_aggregator, stride,  rearr, iosystem,base)
    use pio_types, only : pio_internal_error, pio_rearr_none
    integer(i4), intent(in) :: comp_rank
    integer(i4), intent(in) :: comp_comm
    integer(i4), intent(in) :: num_iotasks 
    integer(i4), intent(in) :: num_aggregator
    integer(i4), intent(in) :: stride
    integer(i4), intent(in) :: rearr
    type (iosystem_desc_t), intent(out)  :: iosystem  ! io descriptor to initalize

    integer(i4), intent(in),optional :: base
    
    integer(i4) :: n_iotasks
    integer(i4) :: length
    integer(i4) :: ngseg,io_rank,i,lbase, io_comm,ierr 
    integer(i4) :: lstride, itmp
    integer(i4), pointer :: iotmp(:),iotmp2(:)

    integer :: mpi_comm_io, intercomm

    character(len=5) :: cb_nodes
    logical(log_kind), parameter :: check = .true.
    logical :: async_setup = .false.

    integer(i4) :: j

    integer(i4) :: mpi_group_world, mpi_group_io, mpi_group_compute

    integer(i4) :: iotask
    integer(i4) :: rearrFlag


    iosystem%error_handling = PIO_internal_error
    iosystem%union_comm = comp_comm
    iosystem%comp_comm = comp_comm
    iosystem%comp_rank = comp_rank
    iosystem%intercomm = MPI_COMM_NULL
    iosystem%my_comm = comp_comm
    iosystem%async_interface = .false.
    iosystem%info = mpi_info_null


    if(comp_comm == MPI_COMM_NULL) then
       call piodie("piolib_mod.f90",1632,'invalid comp_comm in pio_init')
    end if

    call mpi_comm_size(comp_comm,iosystem%num_tasks,ierr)
    iosystem%num_comptasks = iosystem%num_tasks
    iosystem%union_rank = comp_rank
    iosystem%rearr = rearr

    if(check) call checkmpireturn('init: after call to comm_size: ',ierr)
    ! ---------------------------------------
    ! need some more error checking code for 
    ! setting of number of io nodes
    ! ---------------------------------------

    n_iotasks=num_iotasks

    if (n_iotasks>iosystem%num_tasks) then
       n_iotasks=iosystem%num_tasks
       if (iosystem%comp_rank==0) then
          print *,'***warning, reducing io tasks to ',n_iotasks, &
               ' because there are not enough processors'
       endif
    endif

    lbase = 0
    ! unless you are using all procs, shift off the masterproc
    if(n_iotasks<iosystem%num_tasks) then
       lbase=1
    end if
    if (present(base)) then
       if(base>=0 .and. base<iosystem%num_tasks) lbase = base
    endif

    if(debug) print *,'init: iosystem%num_tasks,n_iotasks,num_aggregator: ',iosystem%num_tasks,n_iotasks,num_aggregator

    ! --------------------------
    ! select which nodes are io
    ! nodes and set ioproc
    ! --------------------------
    lstride = stride
    ! Check sanity of input arguments

    call mpi_bcast(iosystem%rearr, 1, mpi_integer, 0, iosystem%comp_comm, ierr)
    call mpi_bcast(n_iotasks, 1, mpi_integer, 0, iosystem%comp_comm, ierr)
    call mpi_bcast(lstride, 1, mpi_integer, 0, iosystem%comp_comm, ierr)
    call mpi_bcast(lbase, 1, mpi_integer, 0, iosystem%comp_comm, ierr)

    if (lbase+(n_iotasks-1)*lstride >= iosystem%num_tasks .and. lstride > 0  .and. n_iotasks > 0) then
       print *,"piolib_mod.f90",1680,lbase,n_iotasks,lstride,iosystem%num_tasks
       call piodie("piolib_mod.f90",1681,'not enough procs for the stride')
    endif

    iosystem%ioproc = .false.


  iosystem%num_iotasks = n_iotasks
    call alloc_check(iosystem%ioranks,n_iotasks,'init:n_ioranks')

    do i=1,n_iotasks
       iosystem%ioranks(i)=(lbase + (i-1)*lstride)

       if (iosystem%ioranks(i)>=iosystem%num_tasks) then
          call piodie( "piolib_mod.f90",1758, &
               'tried to assign io processor beyond max rank ',&
               iosystem%ioranks(i), &
               ' num_tasks=',iosystem%num_tasks )
       endif

       if(comp_rank == iosystem%ioranks(i))  iosystem%ioproc = .true.
    enddo


    iosystem%iomaster = iosystem%ioranks(1)
    iosystem%ioroot = iosystem%ioranks(1)




    if(debug) print *,'init: iam: ',comp_rank,' before allocate(status): n_iotasks: ',n_iotasks

    if (iosystem%rearr == PIO_rearr_none) then
       iosystem%userearranger= .false.
    else
       iosystem%userearranger= .true.
    endif

    call mpi_info_create(iosystem%info,ierr)

    !---------------------------------
    ! initialize the rearranger system 
    !---------------------------------

    if (iosystem%userearranger) then
       call rearrange_init(iosystem)
    endif

    iosystem%io_rank=-1
    call mpi_comm_group(comp_comm,mpi_group_world,ierr)
    if(check) call checkmpireturn('init: after call to comm_group: ',ierr)

    call mpi_group_incl(mpi_group_world,n_iotasks,iosystem%ioranks,mpi_group_io,ierr)
    if(check) call checkmpireturn('init: after call to group_range_incl: ',ierr)

    if(DebugAsync) print *,"piolib_mod.f90",1802,'n: ',n_iotasks, ' r: ', &
     iosystem%ioranks, ' g: ',mpi_group_io

    !-----------------------
    ! setup io_comm and io_rank
    !-----------------------

    call mpi_comm_create(comp_comm,mpi_group_io,iosystem%io_comm,ierr)
    if(check) call checkmpireturn('init: after call to comm_create: ',ierr)

    call mpi_group_free(mpi_group_io, ierr)
    if(check) call checkmpireturn('init: after call to group_free: ',ierr)

    
    if(iosystem%ioproc) call mpi_comm_rank(iosystem%io_comm,iosystem%io_rank,ierr)
    if(check) call checkmpireturn('init: after call to comm_rank: ',ierr)
    ! turn on mpi-io aggregation 
    !DBG    print *,'PIO_init: before call to setnumagg'
    itmp = num_aggregator
    call mpi_bcast(itmp, 1, mpi_integer, 0, iosystem%comp_comm, ierr)


    if(debug) print *,1824,'init: iam: ',comp_rank,'io processor: ',iosystem%ioproc, 'io rank ',&
         iosystem%io_rank, iosystem%iomaster, iosystem%comp_comm, iosystem%io_comm


    if(itmp .gt. 0) then 
       write(cb_nodes,('(i5)')) itmp
       call PIO_set_hint(iosystem,"cb_nodes",trim(adjustl(cb_nodes)))
    endif

    iosystem%num_aiotasks = iosystem%num_iotasks
    iosystem%numost = PIO_NUM_OST
    if(debug) print *,1846,'init: iam: ',comp_rank,'io processor: ',iosystem%ioproc, 'io rank ',&
         iosystem%io_rank, iosystem%iomaster, iosystem%comp_comm, iosystem%io_comm

  end subroutine init_intracom


!> 
!! @public
!! @ingroup PIO_init
!! @brief Initialize the pio subsystem.
!! @details  This is a collective call.  Input parameters are read on comp_rank=0
!!   values on other tasks are ignored.  This variation of PIO_init sets up a distinct set of tasks
!!   to handle IO, these tasks do not return from this call.  Instead they go to an internal loop 
!!   and wait to receive further instructions from the computational tasks 
!! @param component_count The number of computational components to associate with this IO component
!! @param peer_comm  The communicator from which all other communicator arguments are derived
!! @param comp_comms The computational communicator for each of the computational components
!! @param io_comm    The io communicator 
!! @param iosystem a derived type which can be used in subsequent pio operations (defined in PIO_types).
!<
! Subprogram not used   subroutine init_intercom(component_count, peer_comm, comp_comms, io_comm, iosystem)
! Subprogram not used     use pio_types, only : pio_internal_error, pio_rearr_box
! Subprogram not used     integer, intent(in) :: component_count
! Subprogram not used     integer, intent(in) :: peer_comm
! Subprogram not used     integer, intent(in) :: comp_comms(component_count)   !  The compute communicator
! Subprogram not used     integer, intent(in) :: io_comm     !  The io communicator
! Subprogram not used 
! Subprogram not used     type (iosystem_desc_t), intent(out)  :: iosystem(component_count)  ! io descriptor to initalize
! Subprogram not used 
! Subprogram not used     integer :: ierr
! Subprogram not used     logical :: is_inter
! Subprogram not used     logical, parameter :: check=.true.
! Subprogram not used   
! Subprogram not used     integer :: i, j, iam, io_leader, comp_leader
! Subprogram not used     integer(i4), pointer :: iotmp(:)
! Subprogram not used     character(len=5) :: cb_nodes
! Subprogram not used     integer :: itmp
! Subprogram not used     
! Subprogram not used     do i=1,component_count
! Subprogram not used        iosystem(i)%error_handling = PIO_internal_error
! Subprogram not used        iosystem(i)%comp_comm = comp_comms(i)
! Subprogram not used        iosystem(i)%io_comm = io_comm
! Subprogram not used        iosystem(i)%info = mpi_info_null
! Subprogram not used        iosystem(i)%comp_rank= -1
! Subprogram not used        iosystem(i)%io_rank  = -1
! Subprogram not used        iosystem(i)%async_interface = .true.
! Subprogram not used        iosystem(i)%comproot = MPI_PROC_NULL
! Subprogram not used        iosystem(i)%ioroot = MPI_PROC_NULL
! Subprogram not used        iosystem(i)%compmaster= MPI_PROC_NULL
! Subprogram not used        iosystem(i)%iomaster = MPI_PROC_NULL 
! Subprogram not used        iosystem(i)%numOST = PIO_num_OST
! Subprogram not used 
! Subprogram not used 
! Subprogram not used        if(io_comm/=MPI_COMM_NULL) then
! Subprogram not used           ! Find the rank of the io leader in peer_comm
! Subprogram not used           call mpi_comm_rank(io_comm,iosystem(i)%io_rank, ierr)
! Subprogram not used           if(iosystem(i)%io_rank==0) then 
! Subprogram not used              call mpi_comm_rank(peer_comm, iam, ierr)
! Subprogram not used           else
! Subprogram not used              iam = -1
! Subprogram not used           end if
! Subprogram not used           call mpi_allreduce(iam, io_leader, 1, mpi_integer, MPI_MAX, peer_comm, ierr)
! Subprogram not used           call CheckMPIReturn('Call to MPI_ALLREDUCE()',ierr,"<stdin>",1918)
! Subprogram not used           ! Find the rank of the comp leader in peer_comm
! Subprogram not used           iam = -1
! Subprogram not used           call mpi_allreduce(iam, comp_leader, 1, mpi_integer, MPI_MAX, peer_comm, ierr)
! Subprogram not used           call CheckMPIReturn('Call to MPI_ALLREDUCE()',ierr,"<stdin>",1922)
! Subprogram not used           ! create the intercomm
! Subprogram not used           call mpi_intercomm_create(io_comm, 0, peer_comm, comp_leader, i, iosystem(i)%intercomm, ierr)
! Subprogram not used           ! create the union_comm
! Subprogram not used           call mpi_intercomm_merge(iosystem(i)%intercomm, .true., iosystem(i)%union_comm, ierr)
! Subprogram not used        else
! Subprogram not used           ! Find the rank of the io leader in peer_comm
! Subprogram not used           iam = -1
! Subprogram not used           call mpi_allreduce(iam, io_leader, 1, mpi_integer, MPI_MAX, peer_comm, ierr)
! Subprogram not used           call CheckMPIReturn('Call to MPI_ALLREDUCE()',ierr,"<stdin>",1931)
! Subprogram not used 
! Subprogram not used           ! Find the rank of the comp leader in peer_comm
! Subprogram not used           iosystem(i)%comp_rank = -1
! Subprogram not used           if(comp_comms(i)/=MPI_COMM_NULL) then
! Subprogram not used              call mpi_comm_rank(comp_comms(i),iosystem(i)%comp_rank, ierr)          
! Subprogram not used              if(iosystem(i)%comp_rank==0) then
! Subprogram not used                 call mpi_comm_rank(peer_comm, iam, ierr)
! Subprogram not used              else
! Subprogram not used                 iam=-1
! Subprogram not used              end if
! Subprogram not used           end if
! Subprogram not used           call mpi_allreduce(iam, comp_leader, 1, mpi_integer, MPI_MAX, peer_comm, ierr)
! Subprogram not used           call CheckMPIReturn('Call to MPI_ALLREDUCE()',ierr,"<stdin>",1944)
! Subprogram not used 
! Subprogram not used           ! create the intercomm
! Subprogram not used           call mpi_intercomm_create(comp_comms(i), 0, peer_comm, io_leader, i, iosystem(i)%intercomm, ierr)
! Subprogram not used           ! create the union comm
! Subprogram not used           call mpi_intercomm_merge(iosystem(i)%intercomm, .false., iosystem(i)%union_comm, ierr)
! Subprogram not used        end if
! Subprogram not used        if(Debugasync) print *,"piolib_mod.f90",1951,i, iosystem(i)%intercomm, iosystem(i)%union_comm
! Subprogram not used 
! Subprogram not used        if(iosystem(i)%union_comm /= MPI_COMM_NULL) then
! Subprogram not used           call mpi_comm_rank(iosystem(i)%union_comm, iosystem(i)%union_rank, ierr)
! Subprogram not used           if(check) call checkmpireturn('init: after call to comm_rank: ',ierr)
! Subprogram not used           call mpi_comm_size(iosystem(i)%union_comm, iosystem(i)%num_tasks, ierr)
! Subprogram not used           if(check) call checkmpireturn('init: after call to comm_size: ',ierr)
! Subprogram not used 
! Subprogram not used              
! Subprogram not used           if(io_comm /= MPI_COMM_NULL) then
! Subprogram not used              call mpi_comm_size(io_comm, iosystem(i)%num_iotasks, ierr)
! Subprogram not used              if(check) call checkmpireturn('init: after call to comm_size: ',ierr)
! Subprogram not used 
! Subprogram not used              if(iosystem(i)%io_rank==0) then
! Subprogram not used                 iosystem(i)%iomaster = MPI_ROOT
! Subprogram not used                 iosystem(i)%ioroot = iosystem(i)%union_rank
! Subprogram not used              end if
! Subprogram not used              iosystem(i)%ioproc = .true.
! Subprogram not used              iosystem(i)%compmaster = 0
! Subprogram not used 
! Subprogram not used              call pio_msg_handler_init(io_comm, iosystem(i)%io_rank)
! Subprogram not used           end if
! Subprogram not used 
! Subprogram not used 
! Subprogram not used           if(comp_comms(i) /= MPI_COMM_NULL) then
! Subprogram not used              call mpi_comm_size(comp_comms(i), iosystem(i)%num_comptasks, ierr)
! Subprogram not used              if(check) call checkmpireturn('init: after call to comm_size: ',ierr)
! Subprogram not used 
! Subprogram not used              iosystem(i)%iomaster = 0
! Subprogram not used              iosystem(i)%ioproc = .false.
! Subprogram not used              if(iosystem(i)%comp_rank==0) then
! Subprogram not used                 iosystem(i)%compmaster = MPI_ROOT
! Subprogram not used                 iosystem(i)%comproot = iosystem(i)%union_rank
! Subprogram not used              end if
! Subprogram not used 
! Subprogram not used           end if
! Subprogram not used 
! Subprogram not used           iosystem(i)%userearranger = .true.
! Subprogram not used           iosystem(i)%rearr = PIO_rearr_box
! Subprogram not used           
! Subprogram not used           if(Debugasync) print *,"piolib_mod.f90",1991
! Subprogram not used           
! Subprogram not used           call MPI_allreduce(iosystem(i)%comproot, j, 1, MPI_INTEGER, MPI_MAX,iosystem(i)%union_comm,ierr)
! Subprogram not used           call CheckMPIReturn('Call to MPI_ALLREDUCE()',ierr,"<stdin>",1994)
! Subprogram not used           
! Subprogram not used           iosystem%comproot=j
! Subprogram not used           call MPI_allreduce(iosystem(i)%ioroot, j, 1, MPI_INTEGER, MPI_MAX,iosystem(i)%union_comm,ierr)
! Subprogram not used           call CheckMPIReturn('Call to MPI_ALLREDUCE()',ierr,"<stdin>",1998)
! Subprogram not used 
! Subprogram not used           iosystem%ioroot=j
! Subprogram not used 
! Subprogram not used           if(Debugasync) print *,"piolib_mod.f90",2002, i, iosystem(i)%comproot, iosystem(i)%ioroot
! Subprogram not used 
! Subprogram not used           if(io_comm/=MPI_COMM_NULL) then
! Subprogram not used              call mpi_bcast(iosystem(i)%num_comptasks, 1, mpi_integer, iosystem(i)%compmaster,iosystem(i)%intercomm, ierr)
! Subprogram not used 
! Subprogram not used              call mpi_bcast(iosystem(i)%num_iotasks, 1, mpi_integer, iosystem(i)%iomaster, iosystem(i)%intercomm, ierr)
! Subprogram not used 
! Subprogram not used              call alloc_check(iotmp,iosystem(i)%num_iotasks,'init:iotmp')
! Subprogram not used              iotmp(:) = 0
! Subprogram not used              iotmp( iosystem(i)%io_rank+1)=iosystem(i)%union_rank
! Subprogram not used 
! Subprogram not used           end if
! Subprogram not used           if(comp_comms(i)/=MPI_COMM_NULL) then
! Subprogram not used              call mpi_bcast(iosystem(i)%num_comptasks, 1, mpi_integer, iosystem(i)%compmaster, iosystem(i)%intercomm, ierr)
! Subprogram not used 
! Subprogram not used              call mpi_bcast(iosystem(i)%num_iotasks, 1, mpi_integer, iosystem(i)%iomaster, iosystem(i)%intercomm, ierr)
! Subprogram not used 
! Subprogram not used              call alloc_check(iotmp,iosystem(i)%num_iotasks,'init:iotmp')
! Subprogram not used              iotmp(:)=0
! Subprogram not used 
! Subprogram not used           end if
! Subprogram not used 
! Subprogram not used           iosystem(i)%my_comm = iosystem(i)%intercomm
! Subprogram not used 
! Subprogram not used           call alloc_check(iosystem(i)%ioranks, iosystem(i)%num_iotasks,'init:n_ioranks')
! Subprogram not used           if(Debugasync) print *,"piolib_mod.f90",2027,iotmp
! Subprogram not used           call MPI_allreduce(iotmp,iosystem(i)%ioranks,iosystem(i)%num_iotasks,MPI_INTEGER,MPI_MAX,iosystem(i)%union_comm,ierr)
! Subprogram not used           call CheckMPIReturn('Call to MPI_ALLREDUCE()',ierr,"<stdin>",2029)
! Subprogram not used 
! Subprogram not used           if(Debugasync) print *,"piolib_mod.f90",2031,iosystem(i)%ioranks
! Subprogram not used           call dealloc_check(iotmp)
! Subprogram not used           
! Subprogram not used           !---------------------------------
! Subprogram not used           ! initialize the rearranger system 
! Subprogram not used           !---------------------------------
! Subprogram not used           if (iosystem(i)%userearranger) then
! Subprogram not used              call rearrange_init(iosystem(i))
! Subprogram not used           endif
! Subprogram not used        end if
! Subprogram not used     
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     if(DebugAsync) print*,"piolib_mod.f90",2069, iosystem(1)%ioranks
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     iosystem%num_aiotasks = iosystem%num_iotasks
! Subprogram not used     iosystem%numost = PIO_NUM_OST
! Subprogram not used 
! Subprogram not used     ! This routine does not return
! Subprogram not used     if(io_comm /= MPI_COMM_NULL) call pio_msg_handler(component_count,iosystem) 
! Subprogram not used     
! Subprogram not used     if(DebugAsync) print*,"piolib_mod.f90",2078, iosystem(1)%ioranks
! Subprogram not used   end subroutine init_intercom

!>
!! @public
!! @defgroup PIO_recommend_iotasks PIO_recommend_iotasks
!! @brief Recommend a subset of tasks in comm to use as IO tasks
!! @details  This subroutine will give PIO's best recommendation for the number and
!!    location of iotasks for a given system there is no requirement to follow this recommendation.
!!    Using the recommendation requires that PIO_BOX_RERRANGE be used
!! @param A communicator of mpi tasks to choose from
!! @param miniotasks \em optional The minimum number of IO tasks the caller desires
!! @param maxiotasks \em optional The maximum number of IO tasks the caller desires
!! @param iotask if true pio recommends that this task be used as an iotask
!<

! Subprogram not used   subroutine pio_recommend_iotasks(comm, ioproc, numiotasks, miniotasks, maxiotasks )
! Subprogram not used     integer, intent(in) :: comm
! Subprogram not used     logical, intent(out) :: ioproc
! Subprogram not used     integer, intent(out) :: numiotasks
! Subprogram not used     integer, optional, intent(in) :: miniotasks, maxiotasks
! Subprogram not used 
! Subprogram not used     integer :: num_tasks, ierr, iotask, iotasks, iam
! Subprogram not used 
! Subprogram not used     integer(i4), pointer :: iotmp(:),iotmp2(:)
! Subprogram not used 
! Subprogram not used     call mpi_comm_size(comm,num_tasks,ierr)
! Subprogram not used     call mpi_comm_rank(comm,iam,ierr)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   end subroutine pio_recommend_iotasks


!> 
!! @public
!! @defgroup PIO_set_hint  PIO_set_hint
!! @brief set file system hints using mpi_info_set
!! @details This is a collective call which expects the following parameters:
!! @param iosystem @copydoc io_desc_t
!! @param hint  the string name of the hint to define
!! @param hintval  the string value to set the hint to
!! @retval ierr @copydoc  error_return
!<
  subroutine PIO_set_hint(iosystem, hint, hintval)
    type (iosystem_desc_t), intent(inout)  :: iosystem  ! io descriptor to initalize
    character(len=*), intent(in) :: hint, hintval
    
    integer :: ierr
  end subroutine PIO_set_hint


!> 
!! @public
!! @ingroup PIO_finalize 
!! @brief finalizes the pio subsystem.
!! @details This is a collective call which expects the following parameters
!! @param iosystem : @copydoc io_desc_t
!! @retval ierr @copydoc  error_return
!<
  subroutine finalize(iosystem,ierr)
     type (iosystem_desc_t), intent(inout) :: iosystem 
     integer(i4), intent(out) :: ierr
     
     integer :: msg

     if(iosystem%async_interface .and. iosystem%comp_rank==0) then
        !print *,'IAM: ',iosystem%comp_rank, ' ASYNC in finalize'
        msg = PIO_MSG_EXIT
        call mpi_send(msg, 1, mpi_integer, iosystem%ioroot, 1, iosystem%union_comm, ierr)
     end if
     If (associated (iosystem%ioranks)) deallocate (iosystem%ioranks)
     if(iosystem%info .ne. mpi_info_null) then 
        call mpi_info_free(iosystem%info,ierr) 
        iosystem%info=mpi_info_null
        !print *,'IAM: ',iosystem%comp_rank, ' finalize (1) error = ', ierr
     endif
     if(iosystem%io_comm .ne. mpi_comm_null) then 
        call mpi_comm_free(iosystem%io_comm,ierr)
        iosystem%io_comm=mpi_comm_null
        !print *,'IAM: ',iosystem%comp_rank, ' finalize (2) error = ', ierr
     endif
     ierr = 0

  end subroutine finalize


!>
!! @public
!! @ingroup PIO_getnumiotasks
!! @brief This returns the number of IO-tasks that PIO is using 
!! @param iosystem : a defined pio system descriptor, see PIO_types
!! @param numiotasks : the number of IO-tasks
!<
! Subprogram not used    subroutine getnumiotasks(iosystem,numiotasks)
! Subprogram not used        type (iosystem_desc_t), intent(in) :: iosystem
! Subprogram not used        integer(i4), intent(out) :: numiotasks
! Subprogram not used        numiotasks = iosystem%num_iotasks
! Subprogram not used    end subroutine getnumiotasks




  !=============================================
  !  dupiodesc:
  !
  !   duplicate the io descriptor
  !
  !=============================================


  ! rml: possible problem here wrt dubbing the box rearranger
  ! data, as well as maybe the mct rearranger???

!> 
!! @public 
!! @ingroup PIO_dupiodesc
!! @brief duplicates an existing io descriptor
!! @details
!! @param src :  an io description handle returned from @ref PIO_initdecomp (see PIO_types)
!! @param dest : the newly created io descriptor with the same characteristcs as src.
!<
! Subprogram not used   subroutine dupiodesc(src,dest)
! Subprogram not used 
! Subprogram not used     integer :: n
! Subprogram not used     type (io_desc_t), intent(in) :: src
! Subprogram not used     type (io_desc_t), intent(inout) :: dest
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     dest%glen        =  src%glen
! Subprogram not used     if(associated(src%start)) then
! Subprogram not used        n = size(src%start)
! Subprogram not used        allocate(dest%start(n))
! Subprogram not used        dest%start(:)       =  src%start(:)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if(associated(src%count)) then
! Subprogram not used        n = size(src%count)
! Subprogram not used        allocate(dest%count(n))
! Subprogram not used        dest%count(:)       =  src%count(:)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     !dbg    print *,'before dupiodesc2'
! Subprogram not used     call dupiodesc2(src%read, dest%read)
! Subprogram not used     call dupiodesc2(src%write, dest%write)
! Subprogram not used     !dbg    print *,'after dupiodesc2'
! Subprogram not used 
! Subprogram not used     dest%basetype = src%basetype
! Subprogram not used 
! Subprogram not used     if(associated(src%dest_ioproc)) then 
! Subprogram not used        n = size(src%dest_ioproc)
! Subprogram not used        allocate(dest%dest_ioproc(n))
! Subprogram not used        dest%dest_ioproc(:) = src%dest_ioproc(:)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if(associated(src%dest_ioindex)) then 
! Subprogram not used        n = size(src%dest_ioindex)
! Subprogram not used        allocate(dest%dest_ioindex(n))
! Subprogram not used        dest%dest_ioindex(:) = src%dest_ioindex(:)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if(associated(src%rfrom)) then 
! Subprogram not used        n = size(src%rfrom)
! Subprogram not used        allocate(dest%rfrom(n))
! Subprogram not used        dest%rfrom(:) = src%rfrom(:)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if(associated(src%rtype)) then 
! Subprogram not used        n = size(src%rtype)
! Subprogram not used        allocate(dest%rtype(n))
! Subprogram not used        dest%rtype(:) = src%rtype(:)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if(associated(src%scount)) then 
! Subprogram not used        n = size(src%scount)
! Subprogram not used        allocate(dest%scount(n))
! Subprogram not used        dest%scount(:) = src%scount(:)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if(associated(src%stype)) then 
! Subprogram not used        n = size(src%stype)
! Subprogram not used        allocate(dest%stype(n))
! Subprogram not used        dest%stype(:) = src%stype(:)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     call copy_decompmap(src%iomap,dest%iomap)
! Subprogram not used     call copy_decompmap(src%compmap,dest%compmap)
! Subprogram not used 
! Subprogram not used     dest%compsize = src%compsize
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   end subroutine dupiodesc

  !=============================================
  !  copy_decompmap:
  !
  !   copy decompmap_t data structures
  !
  !=============================================

! Subprogram not used   subroutine copy_decompmap(src,dest)
! Subprogram not used     use pio_types, only : decompmap_t
! Subprogram not used     type (decompmap_t), intent(in) :: src
! Subprogram not used     type (decompmap_t), intent(inout) :: dest
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     dest%start    = src%start
! Subprogram not used     dest%length   = src%length
! Subprogram not used 
! Subprogram not used   end subroutine copy_decompmap

!> 
!! @public 
!! @ingroup PIO_setiotype
!! @brief sets the desired type of io to perform
!! @details
!! @param file @copydoc file_desc_t
!! @param iotype : @copydoc PIO_iotype
!! @param rearr : @copydoc PIO_rearr_method
!<
! Subprogram not used   subroutine setiotype(file,iotype,rearr)
! Subprogram not used 
! Subprogram not used     type (file_desc_t), intent(inout) :: file
! Subprogram not used     integer(i4), intent(in) :: iotype 
! Subprogram not used     integer(i4), intent(in) :: rearr
! Subprogram not used 
! Subprogram not used     file%iotype = iotype
! Subprogram not used     file%iosystem%rearr = rearr
! Subprogram not used 
! Subprogram not used   end subroutine setiotype

!>
!! @public
!! @ingroup PIO_numtoread
!! @brief returns the global number of words to read for this io descriptor
!! @details
!! @param iodesc : @copydoc io_desc_t
!! @retval num   :  the number of words to read 
!<
  integer function numtoread(iodesc) result(num)

    type (io_desc_t) :: iodesc

    num = iodesc%read%n_words

  end function numtoread

!>
!! @public
!! @ingroup PIO_numtowrite
!! @brief returns the global number of words to write for this io descriptor
!! @details
!! @param iodesc : @copydoc io_desc_t
!<
  integer function numtowrite(iodesc) result(num)

    type (io_desc_t) :: iodesc

    num = iodesc%write%n_words

  end function numtowrite

!> 
!! @public
!! @ingroup PIO_createfile 
!! @brief create a file using pio
!! @details  Input parameters are read on comp task 0 and ignored elsewhere
!! @param iosystem : a defined pio system descriptor created by a call to @ref PIO_init (see PIO_types)
!! @param file	:  the returned file descriptor
!! @param iotype : @copydoc PIO_iotype
!! @param fname : the name of the file to open
!! @param amode_in : the creation mode flag. the following flags are available: PIO_clobber, PIO_noclobber. 
!! @retval ierr @copydoc error_return
!<
  integer function createfile(iosystem, file,iotype, fname, amode_in) result(ierr)
    type (iosystem_desc_t), intent(inout), target :: iosystem
    type (file_desc_t), intent(out) :: file
    integer, intent(in) :: iotype
    character(len=*), intent(in)  :: fname
    integer, optional, intent(in) :: amode_in
    
    ! ===================
    !  local variables
    ! ===================
    logical :: iscallback
    integer    :: amode
    integer :: msg
    logical, parameter :: check = .true.
    character(len=9) :: rd_buffer
    character(len=4) :: stripestr
    character(len=9) :: stripestr2
    character(len=char_len)  :: myfname

    if(debug.or.debugasync) print *,'createfile: {comp,io}_rank:',iosystem%comp_rank,iosystem%io_rank, &
         'io proc: ',iosystem%ioproc,iosystem%async_interface, iotype
    ierr=PIO_noerr
    

    if(present(amode_in)) then
       amode = amode_in
    else	
       amode = 0
    end if

    file%iotype = iotype 
    myfname = fname

    if(.not. (iosystem%async_interface .and. iosystem%ioproc)) then
       call mpi_bcast(amode, 1, MPI_INTEGER, 0, iosystem%comp_comm, ierr)
       call mpi_bcast(file%iotype, 1, MPI_INTEGER, 0, iosystem%comp_comm, ierr)

       if(len(fname) > char_len) then
          print *,'Length of filename exceeds compile time max, increase char_len in pio_kinds and recompile', len(fname), char_len
          call piodie( "piolib_mod.f90",2450)
       end if

       call mpi_bcast(myfname, len(fname), mpi_character, 0, iosystem%comp_comm, ierr)
    end if

    file%iosystem => iosystem

    !--------------------------------
    ! set some iotype specific stuff
    !--------------------------------


    if(file%iotype==pio_iotype_netcdf4p .or. file%iotype==pio_iotype_netcdf4c) then
       print *, 'WARNING: PIO was not built with NETCDF 4 support changing iotype to netcdf'
       file%iotype = pio_iotype_netcdf
    end if
    if(iosystem%async_interface .and. .not. iosystem%ioproc) then
       msg = PIO_MSG_CREATE_FILE
       if(iosystem%comp_rank==0) then
          call mpi_send(msg, 1, mpi_integer, iosystem%ioroot, 1, iosystem%union_comm, ierr)
       end if

       call mpi_bcast(myfname, char_len, mpi_character, iosystem%compmaster, iosystem%intercomm, ierr)
       call mpi_bcast(iotype, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
       call mpi_bcast(amode, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)

    end if
    select case(iotype)
    case(pio_iotype_pbinary, pio_iotype_direct_pbinary)
       if(present(amode_in) .and. iosystem%io_rank==0) then
          print *, 'warning, the mode argument is currently ignored for binary file operations'
       end if
       ierr = create_mpiio(file,myfname)
    case( pio_iotype_pnetcdf, pio_iotype_netcdf, pio_iotype_netcdf4p, pio_iotype_netcdf4c)
       if(debug) print *,"piolib_mod.f90",2500,' open: ', trim(myfname), amode
       ierr = create_nf(file,trim(myfname), amode)	
       if(debug .and. iosystem%io_rank==0)print *,"piolib_mod.f90",2502,' open: ', myfname, file%fh, ierr
    case(pio_iotype_binary)
       print *,'createfile: io type not supported'
    end select
    if(ierr==0) file%file_is_open=.true.
	
    if(debug .and. file%iosystem%io_rank==0) print *,"piolib_mod.f90",2518,'open: ',file%fh, myfname

  end function createfile
!>
!! @public
!! @defgroup PIO_setnum_OST PIO_setnum_OST
!! @brief Sets the default number of Lustre Object Storage Targets (OST)
!! @details  When PIO is used on a Lustre filesystem, this subroutine sets the
!!           default number Object Storage targets (OST) to use. PIO
!!           will use min(num_aiotasks,numOST) where num_aiotasks the the
!!           actual number of active iotasks
!! @param iosystem : a defined pio system descriptor created by a call to @ref PIO_init (see PIO_types)
!! @param numOST : The number of OST to use by default
!<
  subroutine PIO_setnum_OST(iosystem,numOST)
     type (iosystem_desc_t), intent(inout), target :: iosystem
     integer(i4) :: numOST
     iosystem%numOST = numOST
  end subroutine PIO_setnum_OST
!>
!! @public
!! @defgroup PIO_getnum_OST PIO_getnum_OST
!! @brief Sets the default number of Lustre Object Storage Targets (OST)
!! @details  When PIO is used on a Lustre filesystem, this subroutine gets the
!!           default number Object Storage targets (OST) to use.
!! @param iosystem : a defined pio system descriptor created by a call to @ref PIO_init (see PIO_types)
!! @retval numOST : The number of OST to use.
!<
  integer function PIO_getnum_OST(iosystem) result(numOST)
     type (iosystem_desc_t), intent(inout), target :: iosystem
     numOST = iosystem%numOST
  end function PIO_getnum_OST
!> 
!! @public
!! @ingroup PIO_openfile 
!! @brief open an existing file using pio
!! @details  Input parameters are read on comp task 0 and ignored elsewhere.
!! @param iosystem : a defined pio system descriptor created by a call to @ref PIO_init (see PIO_types)
!! @param file	:  the returned file descriptor
!! @param iotype : @copybrief PIO_iotype
!! @param fname : the name of the file to open
!! @param mode : a zero value (or PIO_nowrite) specifies the default
!! behavior: open the dataset with read-only access, buffering and
!! caching accesses for efficiency otherwise, the creation mode is
!! PIO_write. setting the PIO_write flag opens the dataset with
!! read-write access. ("writing" means any kind of change to the dataset,
!! including appending or changing data, adding or renaming dimensions,
!! variables, and attributes, or deleting attributes.) 
!! @retval ierr @copydoc error_return
!<
  integer function PIO_openfile(iosystem, file, iotype, fname,mode, CheckMPI) result(ierr)
    type (iosystem_desc_t), intent(inout), target :: iosystem
    type (file_desc_t), intent(out) :: file
    integer, intent(in) :: iotype
    character(len=*), intent(in)  :: fname
    integer, optional, intent(in) :: mode
    logical, optional, intent(in) :: CheckMPI
    ! ===================
    !  local variables
    ! ================
    integer    :: amode, msg
    logical, parameter :: check = .true.
    character(len=9) :: rd_buffer
    character(len=char_len) :: myfname




    if(Debug .or. Debugasync) print *,'PIO_openfile: {comp,io}_rank:',iosystem%comp_rank,iosystem%io_rank,&
         'io proc: ',iosystem%ioproc
    ierr=PIO_noerr

    file%iosystem => iosystem

    if(present(mode)) then
       amode = mode
    else	
       amode = 0
    end if
    !--------------------------------
    ! set some iotype specific stuff
    !--------------------------------

    if(iosystem%num_iotasks.eq.1.and.iotype.eq.pio_iotype_pnetcdf) then	
       file%iotype = iotype 
    else
       file%iotype = iotype 
    end if


    myfname = fname

    if(file%iotype==pio_iotype_netcdf4p .or. file%iotype==pio_iotype_netcdf4c) then
       print *, 'WARNING: PIO was not built with NETCDF 4 support changing iotype to netcdf'
       file%iotype = pio_iotype_netcdf
    end if
    if(.not. (iosystem%ioproc .and. iosystem%async_interface)) then
       call mpi_bcast(amode, 1, MPI_INTEGER, 0, iosystem%comp_comm, ierr)
       call mpi_bcast(file%iotype, 1, MPI_INTEGER, 0, iosystem%comp_comm, ierr)
       if(len(fname) > char_len) then
          print *,'Length of filename exceeds compile time max, increase char_len in pio_kinds and recompile'
          call piodie( "piolib_mod.f90",2641)
       end if

       call mpi_bcast(myfname, len(fname), mpi_character, 0, iosystem%comp_comm, ierr)
    end if

    if(iosystem%async_interface .and. .not. iosystem%ioproc) then
       msg = PIO_MSG_OPEN_FILE
       if(iosystem%comp_rank==0) then
          call mpi_send(msg, 1, mpi_integer, iosystem%ioroot, 1, iosystem%union_comm, ierr)
       end if
       
       call mpi_bcast(myfname, char_len, mpi_character, iosystem%compmaster, iosystem%intercomm, ierr)
       call mpi_bcast(iotype, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
       call mpi_bcast(amode, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
    end if

    select case(iotype)
    case(pio_iotype_pbinary, pio_iotype_direct_pbinary)
       if(amode /=0) then
          print *, 'warning, the mode argument is currently ignored for binary file operations'
       end if
       if (present(CheckMPI)) then
         ierr = open_mpiio(file,myfname, CheckMPI)
       else
         ierr = open_mpiio(file,myfname)
       end if
    case( pio_iotype_pnetcdf, pio_iotype_netcdf, pio_iotype_netcdf4c, pio_iotype_netcdf4p)
       ierr = open_nf(file,myfname,amode)
       if(debug .and. iosystem%io_rank==0)print *,"piolib_mod.f90",2670,' open: ', myfname, file%fh
    case(pio_iotype_binary)   ! appears to be a no-op
    end select
    if(Debug .and. file%iosystem%io_rank==0) print *,"piolib_mod.f90",2679,'open: ',file%fh, myfname
    if(ierr==0) file%file_is_open=.true.

  end function PIO_openfile

!> 
!! @public 
!! @ingroup PIO_syncfile 
!! @brief synchronizing a file forces all writes to complete before the subroutine returns. 
!!
!! @param file @copydoc file_desc_t
!<
  subroutine syncfile(file)
    use piodarray, only : darray_write_complete
    implicit none
    type (file_desc_t), target :: file
    integer :: ierr, msg
    type(iosystem_desc_t), pointer :: ios
     
 
    ios => file%iosystem
    if(ios%async_interface .and. .not. ios%ioproc) then
       msg = PIO_MSG_SYNC_FILE
       if(ios%comp_rank==0) then
          call mpi_send(msg, 1, mpi_integer, ios%ioroot, 1, ios%union_comm, ierr)
       end if
      
       call mpi_bcast(file%fh, 1, mpi_integer, ios%compmaster, ios%intercomm, ierr)
    end if

    select case(file%iotype)
    case( pio_iotype_pnetcdf, pio_iotype_netcdf)
       call darray_write_complete(file)
       ierr = sync_nf(file)
    case(pio_iotype_pbinary, pio_iotype_direct_pbinary)
    case(pio_iotype_binary) 
    end select
  end subroutine syncfile
!> 
!! @public 
!! @ingroup PIO_freedecomp
!! @brief free all allocated storage associated with this decomposition
!! @details
!! @param ios :  a defined pio system descriptor created by call to @ref PIO_init (see PIO_types)
!! @param iodesc @copydoc io_desc_t
!<
  subroutine freedecomp_ios(ios,iodesc)
    implicit none
    type (iosystem_desc_t) :: ios
    type (io_desc_t) :: iodesc
    integer :: ierr, msg

    if(ios%async_interface .and. .not. ios%ioproc) then
       msg = PIO_MSG_FREEDECOMP
       if(ios%comp_rank==0) then
          call mpi_send(msg, 1, mpi_integer, ios%ioroot, 1, ios%union_comm, ierr)
       end if
       call MPI_Barrier(ios%comp_comm,ierr)
       call mpi_bcast(iodesc%async_id,1, mpi_integer, ios%compmaster,ios%intercomm, ierr)
    end if
    call MPI_Barrier(ios%union_comm,ierr)

    iodesc%async_id=-1
    call rearrange_free(ios,iodesc)

    if(ios%ioproc) then
!       if(debug) print *,"piolib_mod.f90",2749,iodesc%write%n_elemtype,iodesc%write%n_words, &
!       iodesc%write%elemtype,iodesc%write%filetype

       if((iodesc%read%filetype .ne. mpi_datatype_null)  &
	  .and. (iodesc%read%filetype .ne. iodesc%write%filetype) .and. &
	  iodesc%read%n_words>0) then 
          call mpi_type_free(iodesc%read%filetype,ierr)
          call checkmpireturn('freedecomp mpi_type_free: ',ierr)
          call mpi_type_free(iodesc%read%elemtype,ierr)
          call checkmpireturn('freedecomp mpi_type_free: ',ierr)
          iodesc%read%filetype=mpi_datatype_null
       endif
       if(iodesc%write%filetype .ne. mpi_datatype_null .and. &
	  iodesc%write%n_words>0) then 
          call mpi_type_free(iodesc%write%filetype,ierr)
          call checkmpireturn('freedecomp mpi_type_free: ',ierr)
          call mpi_type_free(iodesc%write%elemtype,ierr)
          call checkmpireturn('freedecomp mpi_type_free: ',ierr)
          iodesc%write%filetype=mpi_datatype_null
       endif
   
    end if

    if(associated(iodesc%start)) then
       call dealloc_check(iodesc%start,'iodesc%start')
       nullify(iodesc%start)
    end if

    if(associated(iodesc%count)) then
       call dealloc_check(iodesc%count,'iodesc%count')    
       nullify(iodesc%count)
    end if
  end subroutine freedecomp_ios
!>
!! @public 
!! @ingroup PIO_freedecomp
!! @brief free all allocated storage associated with this decomposition
!! @details
!! @param file @copydoc file_desc_t
!! @param iodesc : @copydoc io_desc_t
!! @retval ierr @copydoc error_return
!<
  subroutine freedecomp_file(file,iodesc)
    implicit none
    type (file_desc_t) :: file
    type (io_desc_t) :: iodesc

    call freedecomp_ios(file%iosystem, iodesc)

  end subroutine freedecomp_file

!> 
!! @public
!! @ingroup PIO_closefile
!! @brief close a disk file
!! @details
!! @param file @copydoc file_desc_t
!< 
  subroutine closefile(file)
    use piodarray, only : darray_write_complete
    type (file_desc_t),intent(inout)   :: file

    integer :: ierr, msg
    integer :: iotype 
    logical, parameter :: check = .true.

    if(file%iosystem%async_interface .and. .not. file%iosystem%ioproc) then
       msg = PIO_MSG_CLOSE_FILE
       if(file%iosystem%comp_rank==0) then
          call mpi_send(msg, 1, mpi_integer, file%iosystem%ioroot, 1, file%iosystem%union_comm, ierr)
       end if
       call mpi_bcast(file%fh, 1, mpi_integer, file%iosystem%compmaster, file%iosystem%intercomm, ierr)
    end if

    if(debug .and. file%iosystem%io_rank==0) &
      print *,"piolib_mod.f90",2828,'close: ',file%fh
    iotype = file%iotype 
    select case(iotype)
    case(pio_iotype_pbinary, pio_iotype_direct_pbinary)
       ierr = close_mpiio(file)
    case( pio_iotype_pnetcdf, pio_iotype_netcdf, pio_iotype_netcdf4p, pio_iotype_netcdf4c)
       call darray_write_complete(file)
       ierr = close_nf(file)
    case(pio_iotype_binary)
       print *,'closefile: io type not supported'
    end select
    if(ierr==0) file%file_is_open=.false.



  end subroutine closefile


  !******************************
  ! read_ascii
  !

  subroutine read_ascii(rank,iobuf,size)

    integer, intent(in) :: rank
    real (r8), dimension(:) :: iobuf
    integer, intent(in) :: size

    character(len=80) filename
    integer lun
    integer ios
    integer i

    lun=10+rank
    write(filename,"('fort.',i2)" ) lun
    write(6,*) 'filename is:', filename

    open(lun,file=filename,status='old',iostat=ios)
    if (ios /= 0) then
       write(6,*) rank,': could not open ascii file: ',filename
    endif

    do i=1,size
       read(unit=lun,fmt=*,iostat=ios) iobuf(i)
       if (ios /= 0) then
          write (6,*) rank,': error reading item ',i,' of ',size
          call abort
       endif

    end do

    close(lun)

  end subroutine read_ascii




end module piolib_mod

  !|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
