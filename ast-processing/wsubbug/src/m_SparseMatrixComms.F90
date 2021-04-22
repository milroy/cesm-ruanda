!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!-----------------------------------------------------------------------
! CVS $Id$
! CVS $Name$ 
!BOP -------------------------------------------------------------------
!
! !MODULE: m_SparseMatrixComms -- sparse matrix communications methods.
!
! !DESCRIPTION:
! The {\tt SparseMatrix} datatype provides sparse matrix storage for 
! the parallel matrix-vector multiplication ${\bf y} = {\bf M} {\bf x}$.
! This module provides communications services for the {\tt SparseMatrix}
! type.  These services include scattering matrix elements based on row or 
! column decompositions, gathering of matrix elements to the root, and 
! broadcasting from the root.
!
! {\bf N.B.:}  These routines will not communicate the vector portion
! of a {\tt SparseMatrix}, if it has been initialized.  A WARNING will
! be issued in most cases.  In general, do communication first,  then 
! call {\tt vecinit}.
!
! !INTERFACE:

 module m_SparseMatrixComms

      private   ! except

! !PUBLIC MEMBER FUNCTIONS:
!
      public :: ScatterByColumn
      public :: ScatterByRow
      public :: Gather
      public :: Bcast

    interface ScatterByColumn ; module procedure &
         ScatterByColumnGSMap_
    end interface

    interface ScatterByRow ; module procedure &
         ScatterByRowGSMap_
    end interface

    interface Gather ; module procedure &
	 GM_gather_, &
	 GSM_gather_
    end interface

    interface Bcast ; module procedure Bcast_ ; end interface

! !REVISION HISTORY:
! 13Apr01 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!           and API specifications.
! 10May01 - J.W. Larson <larson@mcs.anl.gov> - added GM_gather_
!           and cleaned up prologues.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='MCT::m_SparseMatrixComms'

 contains

!-------------------------------------------------------------------------
!     Math + Computer Science Division / Argonne National Laboratory     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ScatterByColumnGSMap_ - Column-based scatter for SparseMatrix.
! 
! !DESCRIPTION: This routine scatters the input {\tt SparseMatrix} 
! argument {\tt GsMat} (valid only on the root) to a distributed 
! {\tt SparseMatrix} variable {\tt LsMat} across all the processes 
! present on the communicator associated with the integer handle 
! {\tt comm}.  The decomposition defining the scatter is supplied by the 
! input {\tt GlobalSegMap} argument {\tt columnGSMap}.  The optional 
! output {\tt INTEGER} flag {\tt stat} signifies a successful (failed) 
! operation if it is returned with value zero (nonzero).
!
! {\bf N.B.:}  This routine returns an allocated {\tt SparseMatrix} 
! variable {\tt LsMat}.  The user must destroy this variable when it
! is no longer needed by invoking {\tt SparseMatrix\_Clean()}.
!
! !INTERFACE:

! Subprogram not used  subroutine ScatterByColumnGSMap_(columnGSMap, GsMat, LsMat, root, comm, stat)
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used 
! Subprogram not used    use m_die, only : MP_perr_die,die
! Subprogram not used    use m_stdio
! Subprogram not used    use m_mpif90
! Subprogram not used 
! Subprogram not used    use m_List, only: List
! Subprogram not used    use m_List, only: List_init => init
! Subprogram not used    use m_List, only: List_clean => clean
! Subprogram not used 
! Subprogram not used    use m_GlobalSegMap, only : GlobalSegMap
! Subprogram not used    use m_GlobalSegMap, only : GlobalSegMap_clean => clean
! Subprogram not used 
! Subprogram not used    use m_SparseMatrix, only : SparseMatrix
! Subprogram not used    use m_SparseMatrix, only : SparseMatrix_nRows => nRows
! Subprogram not used    use m_SparseMatrix, only : SparseMatrix_nCols => nCols
! Subprogram not used    use m_SparseMatrix, only : SparseMatrix_SortPermute => SortPermute
! Subprogram not used 
! Subprogram not used    use m_SparseMatrixDecomp, only : SparseMatrixDecompByColumn => ByColumn
! Subprogram not used 
! Subprogram not used    use m_AttrVectComms, only : AttrVect_Scatter => scatter
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS: 
! Subprogram not used !
! Subprogram not used    type(GlobalSegMap), intent(in)    :: columnGSMap
! Subprogram not used    integer,            intent(in)    :: root
! Subprogram not used    integer,            intent(in)    :: comm
! Subprogram not used 
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS: 
! Subprogram not used !
! Subprogram not used    type(SparseMatrix), intent(inout) :: GsMat
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used    type(SparseMatrix), intent(out) :: LsMat
! Subprogram not used    integer, optional,  intent(out) :: stat
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY: 
! Subprogram not used !
! Subprogram not used ! 13Apr01 - J.W. Larson <larson@mcs.anl.gov> - initial API spec.
! Subprogram not used ! 10May01 - J.W. Larson <larson@mcs.anl.gov> - cleaned up prologue.
! Subprogram not used ! 13Jun01 - J.W. Larson <larson@mcs.anl.gov> - Made status flag stat
! Subprogram not used !           optional, and ititilaze it to zero if it is present.
! Subprogram not used ! 09Jul03 - E.T. Ong <eong@mcs.anl.gov> - added sorting to distributed
! Subprogram not used !           matrix elements
! Subprogram not used !EOP
! Subprogram not used !-------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'ScatterByColumnGSMap_'
! Subprogram not used ! GlobalSegMap used to create column decomposition of GsMat
! Subprogram not used   type(GlobalSegMap) :: MatGSMap
! Subprogram not used ! Storage for the number of rows and columns in the SparseMatrix
! Subprogram not used   integer :: NumRowsColumns(2)
! Subprogram not used ! List storage for sorting keys
! Subprogram not used   type(List) :: sort_keys
! Subprogram not used ! Process ID
! Subprogram not used   integer :: myID
! Subprogram not used ! Error flag
! Subprogram not used   integer :: ierr
! Subprogram not used 
! Subprogram not used        ! Initialize stat if present
! Subprogram not used 
! Subprogram not used   if(present(stat)) stat = 0
! Subprogram not used 
! Subprogram not used        ! Which process am I?
! Subprogram not used 
! Subprogram not used   call MPI_COMM_RANK(comm, myID, ierr)
! Subprogram not used   if(ierr /= 0) then
! Subprogram not used 	call MP_perr_die(myname_,"MPI_COMM_RANK() failed",ierr)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! can't scatter vector parts.
! Subprogram not used   if((myID.eq.root) .and. GsMat%vecinit) then
! Subprogram not used       write(stderr,*) myname_,&
! Subprogram not used       "WARNING: will not scatter vector parts of GsMat"
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! Create from columnGSMap the corresponding GlobalSegMap
! Subprogram not used        ! that will decompose GsMat by column the same way.
! Subprogram not used 
! Subprogram not used   call SparseMatrixDecompByColumn(columnGSMap, GsMat, MatGSMap, root, comm)
! Subprogram not used 
! Subprogram not used        ! Broadcast the resulting GlobalSegMap across the communicator
! Subprogram not used 
! Subprogram not used        ! Scatter the matrix element data GsMat%data accordingly
! Subprogram not used 
! Subprogram not used   call AttrVect_Scatter(GsMat%data, LsMat%data, MatGSMap, root, comm, ierr)
! Subprogram not used 
! Subprogram not used   if(ierr /= 0) then
! Subprogram not used      if(present(stat)) then
! Subprogram not used 	write(stderr,*) myname_,"::  AttrVect_Scatter(GsMat%data) failed--stat=", &
! Subprogram not used 	     ierr
! Subprogram not used 	stat = ierr
! Subprogram not used 	return
! Subprogram not used      else
! Subprogram not used 	call die(myname_,"call AttrVect_Scatter(GsMat%data,..",ierr)
! Subprogram not used      endif
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! Now, distribute to all the processes the number of Rows and
! Subprogram not used        ! columns in GsMat (which are valid on the root only at this point)
! Subprogram not used 
! Subprogram not used   if(myID == root) then
! Subprogram not used      NumRowsColumns(1) = SparseMatrix_nRows(GsMat)
! Subprogram not used      NumRowsColumns(2) = SparseMatrix_nCols(GsMat)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   call MPI_Bcast(NumRowsColumns, 2, MP_INTEGER, root, comm, ierr)
! Subprogram not used 
! Subprogram not used   if(ierr /= 0) then
! Subprogram not used 	call MP_perr_die(myname_,"MPI_Bcast(NumRowsColumns...",ierr)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! Unpack NumRowsColumns
! Subprogram not used 
! Subprogram not used   LsMat%nrows = NumRowsColumns(1)
! Subprogram not used   LsMat%ncols = NumRowsColumns(2)
! Subprogram not used 
! Subprogram not used        ! Set the value of vecinit
! Subprogram not used   LsMat%vecinit = .FALSE.
! Subprogram not used 
! Subprogram not used        ! Finally, lets sort the distributed local matrix elements
! Subprogram not used 
! Subprogram not used        ! Sort the matrix entries in sMat by column, then row.  
! Subprogram not used        ! First, create the key list...
! Subprogram not used 
! Subprogram not used   call List_init(sort_keys,'gcol:grow')
! Subprogram not used 
! Subprogram not used        ! Now perform the sort/permute...
! Subprogram not used   call SparseMatrix_SortPermute(LsMat, sort_keys)
! Subprogram not used 
! Subprogram not used        ! Cleanup
! Subprogram not used 
! Subprogram not used   call List_clean(sort_keys) 
! Subprogram not used   call GlobalSegMap_clean(MatGSMap)
! Subprogram not used 
! Subprogram not used  end subroutine ScatterByColumnGSMap_

!-------------------------------------------------------------------------
!     Math + Computer Science Division / Argonne National Laboratory     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ScatterByRowGSMap_ -Row-based scatter for SparseMatrix.
! 
! !DESCRIPTION: This routine scatters the input  {\tt SparseMatrix} 
! argument {\tt GsMat} (valid only on the root) to a distributed 
! {\tt SparseMatrix} variable {\tt LsMat} across all the processes 
! present on the communicator associated with the integer handle 
! {\tt comm}.  The decomposition defining the scatter is supplied by the 
! input {\tt GlobalSegMap} argument {\tt rowGSMap}.  The output integer 
! flag {\tt stat} signifies a successful (failed) operation if it is 
! returned with value zero (nonzero).
!
! {\bf N.B.:}  This routine returns an allocated {\tt SparseMatrix} 
! variable {\tt LsMat}.  The user must destroy this variable when it
! is no longer needed by invoking {\tt SparseMatrix\_Clean()}.
!
! !INTERFACE:

 subroutine ScatterByRowGSMap_(rowGSMap, GsMat, LsMat, root, comm, stat)
!
! !USES:
!
   use m_die, only : MP_perr_die,die
   use m_stdio
   use m_mpif90

   use m_List, only: List
   use m_List, only: List_init => init
   use m_List, only: List_clean => clean

   use m_GlobalSegMap, only : GlobalSegMap
   use m_GlobalSegMap, only : GlobalSegMap_clean => clean

   use m_SparseMatrix, only : SparseMatrix
   use m_SparseMatrix, only : SparseMatrix_nRows => nRows
   use m_SparseMatrix, only : SparseMatrix_nCols => nCols
   use m_SparseMatrix, only : SparseMatrix_SortPermute => SortPermute

   use m_SparseMatrixDecomp, only : SparseMatrixDecompByRow => ByRow

   use m_AttrVectComms, only : AttrVect_Scatter => scatter

   implicit none

! !INPUT PARAMETERS: 
!
   type(GlobalSegMap), intent(in)    :: rowGSMap
   integer,            intent(in)    :: root
   integer,            intent(in)    :: comm

! !INPUT/OUTPUT PARAMETERS: 
!
   type(SparseMatrix), intent(inout) :: GsMat

! !OUTPUT PARAMETERS:
!
   type(SparseMatrix), intent(out) :: LsMat
   integer, optional,  intent(out) :: stat

! !REVISION HISTORY: 
!
! 13Apr01 - J.W. Larson <larson@mcs.anl.gov> - initial API spec.
! 26Apr01 - R.L. Jacob  <jacob@mcs.anl.gov> - fix use statement
!           from SMDecomp so it points to ByRow
! 13Jun01 - J.W. Larson <larson@mcs.anl.gov> - Made status flag stat
!           optional, and initialize it to zero if it is present.
! 09Jul03 - E.T. Ong <eong@mcs.anl.gov> - Added sorting to distributed
!           matrix elements. 
!EOP
!-------------------------------------------------------------------------

  character(len=*),parameter :: myname_=myname//'ScatterByRowGSMap_'
! GlobalSegMap used to create row decomposition of GsMat
  type(GlobalSegMap) :: MatGSMap
! Storage for the number of rows and columns in the SparseMatrix
  integer :: NumRowsColumns(2)
! List storage for sorting keys
  type(List) :: sort_keys
! Process ID
  integer :: myID
! Error flag
  integer :: ierr

       ! Initialize stat to zero (if present)

  if(present(stat)) stat = 0

       ! Which process are we?

  call MPI_COMM_RANK(comm, myID, ierr)
  if(ierr /= 0) then
     call MP_perr_die(myname_,"MPI_COMM_RANK() failed",ierr)
  endif

       ! can't scatter vector parts.
  if((myID.eq.root) .and. GsMat%vecinit) then
      write(stderr,*) myname_,&
      "WARNING: will not scatter vector parts of GsMat."
  endif

       ! Create from rowGSMap the corresponding GlobalSegMap
       ! that will decompose GsMat by row the same way.

  call SparseMatrixDecompByRow(rowGSMap, GsMat, MatGSMap, root, comm)

       ! Scatter the matrix element data GsMat%data accordingly

  call AttrVect_Scatter(GsMat%data, LsMat%data, MatGSMap, root, comm, ierr)
  if(ierr /= 0) then
     if(present(stat)) then
        write(stderr,*) myname_,"::  AttrVect_Scatter(GsMat%data) failed--stat=", &
             ierr
        stat = ierr
        return
     else
        call die(myname_,"call AttrVect_Scatter(GsMat%data,..",ierr)
     endif
  endif

       ! Now, distribute to all the processes the number of rows and
       ! columns in GsMat (which are valid on the root only at this point)

  if(myID == root) then
     NumRowsColumns(1) = SparseMatrix_nRows(GsMat)
     NumRowsColumns(2) = SparseMatrix_nCols(GsMat)
  endif

  call MPI_Bcast(NumRowsColumns, 2, MP_INTEGER, root, comm, ierr)
  if(ierr /= 0) then
     call MP_perr_die(myname_,"MPI_Bcast(NumRowsColumns...",ierr)
  endif

       ! Unpack NumRowsColumns

  LsMat%nrows = NumRowsColumns(1)
  LsMat%ncols = NumRowsColumns(2)

       ! Set the value of vecinit
  LsMat%vecinit = .FALSE.

       ! Sort the matrix entries in sMat by row, then column.  
       ! First, create the key list...

  call List_init(sort_keys,'grow:gcol')

       ! Now perform the sort/permute...
  call SparseMatrix_SortPermute(LsMat, sort_keys)

       ! Cleanup

  call List_clean(sort_keys) 
  call GlobalSegMap_clean(MatGSMap)

 end subroutine ScatterByRowGSMap_

!-------------------------------------------------------------------------
!     Math + Computer Science Division / Argonne National Laboratory     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  GM_gather_ - Gather a distributed SparseMatrix to the root.
!
! !DESCRIPTION: This routine gathers the input distributed 
! {\tt SparseMatrix} argument {\tt LsMat} to the {\tt SparseMatrix} 
! variable {\tt GsMat} on the root.  The decomposition defining the gather 
! is supplied by the input {\tt GlobalMap} argument {\tt GMap}.  The 
! status flag {\tt stat} has value zero (nonzero) if the operation has 
! succeeded (failed).
!
! {\bf N.B.:}  This routine returns an allocated {\tt SparseMatrix} 
! variable {\tt GsMat}.  The user must destroy this variable when it
! is no longer needed by invoking {\tt SparseMatrix\_Clean()}.
! 
! !INTERFACE:

! Subprogram not used  subroutine GM_gather_(LsMat, GsMat, GMap, root, comm, stat)
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used    use m_stdio
! Subprogram not used    use m_die, only : die
! Subprogram not used 
! Subprogram not used    use m_GlobalMap, only: GlobalMap
! Subprogram not used 
! Subprogram not used    use m_SparseMatrix, only: SparseMatrix
! Subprogram not used    use m_SparseMatrix, only: SparseMatrix_nRows => nRows
! Subprogram not used    use m_SparseMatrix, only: SparseMatrix_nCols => nCols
! Subprogram not used 
! Subprogram not used    use m_AttrVectComms, only : AttrVect_gather => gather
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS: 
! Subprogram not used !
! Subprogram not used    type(SparseMatrix), intent(in) :: LsMat
! Subprogram not used    type(GlobalMap),    intent(in) :: GMap
! Subprogram not used    integer,            intent(in) :: root
! Subprogram not used    integer,            intent(in) :: comm
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used    type(SparseMatrix), intent(out) :: GsMat
! Subprogram not used    integer, optional,  intent(out) :: stat
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY: 
! Subprogram not used !
! Subprogram not used ! 13Apr01 - J.W. Larson <larson@mcs.anl.gov> - initial API spec.
! Subprogram not used ! 10May01 - J.W. Larson <larson@mcs.anl.gov> - initial routine and
! Subprogram not used !           prologue
! Subprogram not used ! 13Jun01 - J.W. Larson <larson@mcs.anl.gov> - Made status flag stat
! Subprogram not used !           optional, and ititilaze it to zero if it is present.
! Subprogram not used !EOP
! Subprogram not used !-------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'GM_gather_'
! Subprogram not used   integer :: ierr
! Subprogram not used 
! Subprogram not used        ! if stat is present, initialize its value to zero (success)
! Subprogram not used 
! Subprogram not used   if(present(stat))  stat = 0
! Subprogram not used 
! Subprogram not used   if(LsMat%vecinit) then
! Subprogram not used       write(stderr,*) myname_,&
! Subprogram not used       "WARNING: will not gather vector parts of LsMat."
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   call AttrVect_gather(LsMat%data, GsMat%data, GMap, root, comm, ierr)
! Subprogram not used   if(ierr /= 0) then
! Subprogram not used      if(present(stat)) then
! Subprogram not used         write(stderr,*) myname_,"::  AttrVect_Gather(LsMat%data...) failed--stat=", &
! Subprogram not used              ierr
! Subprogram not used         stat = ierr
! Subprogram not used         return
! Subprogram not used      else
! Subprogram not used         call die(myname_,"call AttrVect_Scatter(LsMat%data...) failed",ierr)
! Subprogram not used      endif
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! For now, the GsMat inherits the number of rows and columns from
! Subprogram not used        ! the corresponding values of LsMat on the root (this should be
! Subprogram not used        ! checked in future versions).
! Subprogram not used 
! Subprogram not used   GsMat%nrows = SparseMatrix_nRows(LsMat)
! Subprogram not used   GsMat%ncols = SparseMatrix_nCols(LsMat)
! Subprogram not used 
! Subprogram not used   GsMat%vecinit = .FALSE.
! Subprogram not used 
! Subprogram not used  end subroutine GM_gather_

!-------------------------------------------------------------------------
!     Math + Computer Science Division / Argonne National Laboratory     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  GSM_gather_ - Gather a distributed SparseMatrix to the root.
!
! !DESCRIPTION: This routine gathers the input distributed 
! {\tt SparseMatrix} argument {\tt LsMat} to the {\tt SparseMatrix} 
! variable {\tt GsMat} on the root.  The decomposition defining the gather 
! is supplied by the input {\tt GlobalSegMap} argument {\tt GSMap}.  The 
! status flag {\tt stat} has value zero (nonzero) if the operation has 
! succeeded (failed).
!
! {\bf N.B.:}  This routine returns an allocated {\tt SparseMatrix} 
! variable {\tt GsMat}.  The user must destroy this variable when it
! is no longer needed by invoking {\tt SparseMatrix\_Clean()}.
! 
! !INTERFACE:

! Subprogram not used  subroutine GSM_gather_(LsMat, GsMat, GSMap, root, comm, stat)
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used    use m_stdio
! Subprogram not used    use m_die, only : die
! Subprogram not used 
! Subprogram not used    use m_GlobalSegMap, only: GlobalSegMap
! Subprogram not used 
! Subprogram not used    use m_SparseMatrix, only: SparseMatrix
! Subprogram not used    use m_SparseMatrix, only: SparseMatrix_nRows => nRows
! Subprogram not used    use m_SparseMatrix, only: SparseMatrix_nCols => nCols
! Subprogram not used 
! Subprogram not used    use m_AttrVectComms, only : AttrVect_gather => gather
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS: 
! Subprogram not used !
! Subprogram not used    type(SparseMatrix), intent(in) :: LsMat
! Subprogram not used    type(GlobalSegMap), intent(in) :: GSMap
! Subprogram not used    integer,            intent(in) :: root
! Subprogram not used    integer,            intent(in) :: comm
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used    type(SparseMatrix), intent(out) :: GsMat
! Subprogram not used    integer, optional,  intent(out) :: stat
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY: 
! Subprogram not used !
! Subprogram not used ! 13Apr01 - J.W. Larson <larson@mcs.anl.gov> - initial API spec.
! Subprogram not used ! 13Jun01 - J.W. Larson <larson@mcs.anl.gov> - Made status flag stat
! Subprogram not used !           optional, and ititilaze it to zero if it is present.
! Subprogram not used !EOP
! Subprogram not used !-------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'GSM_gather_'
! Subprogram not used   integer :: ierr
! Subprogram not used 
! Subprogram not used        ! if stat is present, initialize its value to zero (success)
! Subprogram not used 
! Subprogram not used   if(present(stat))  stat = 0
! Subprogram not used 
! Subprogram not used   if(LsMat%vecinit) then
! Subprogram not used       write(stderr,*) myname_,&
! Subprogram not used       "WARNING: will not gather vector parts of LsMat."
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! Gather the AttrVect component of LsMat to GsMat...
! Subprogram not used 
! Subprogram not used   call AttrVect_gather(LsMat%data, GsMat%data, GSMap, root, comm, ierr)
! Subprogram not used   if(ierr /= 0) then
! Subprogram not used      if(present(stat)) then
! Subprogram not used         write(stderr,*) myname_,"::  AttrVect_Gather(LsMat%data...) failed--stat=", &
! Subprogram not used              ierr
! Subprogram not used         stat = ierr
! Subprogram not used         return
! Subprogram not used      else
! Subprogram not used         call die(myname_,"call AttrVect_Gather(LsMat%data...)",ierr)
! Subprogram not used      endif
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! For now, the GsMat inherits the number of rows and columns from
! Subprogram not used        ! the corresponding values of LsMat on the root (this should be
! Subprogram not used        ! checked in future versions).
! Subprogram not used 
! Subprogram not used   GsMat%nrows = SparseMatrix_nRows(LsMat)
! Subprogram not used   GsMat%ncols = SparseMatrix_nCols(LsMat)
! Subprogram not used 
! Subprogram not used   GsMat%vecinit = .FALSE.
! Subprogram not used 
! Subprogram not used  end subroutine GSM_gather_

!-------------------------------------------------------------------------
!     Math + Computer Science Division / Argonne National Laboratory     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Bcast_ - Broadcast a SparseMatrix.
! 
! !DESCRIPTION: This routine broadcasts the {\tt SparseMatrix} argument 
! {\tt sMat} from the root to all processes on the communicator associated
! with the communicator handle {\tt comm}.  The status flag {\tt stat} 
! has value zero if the operation has succeeded.
!
! {\bf N.B.:}  This routine returns an allocated {\tt SparseMatrix} 
! variable {\tt sMat}.  The user must destroy this variable when it
! is no longer needed by invoking {\tt SparseMatrix\_Clean()}.
!
! {\bf N.B.:}  This routine will exit with an error if the vector portion
! of {\tt sMat} has been initialized prior to broadcast.
!
! !INTERFACE:

! Subprogram not used  subroutine Bcast_(sMat, root, comm, stat)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used 
! Subprogram not used    use m_die, only : MP_perr_die,die
! Subprogram not used    use m_stdio
! Subprogram not used    use m_mpif90
! Subprogram not used 
! Subprogram not used    use m_GlobalSegMap, only: GlobalSegMap
! Subprogram not used 
! Subprogram not used    use m_AttrVectComms, only : AttrVect_bcast => bcast
! Subprogram not used 
! Subprogram not used    use m_SparseMatrix, only: SparseMatrix
! Subprogram not used    use m_SparseMatrix, only: SparseMatrix_nRows => nRows
! Subprogram not used    use m_SparseMatrix, only: SparseMatrix_nCols => nCols
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS: 
! Subprogram not used !
! Subprogram not used    integer,            intent(in) :: root
! Subprogram not used    integer,            intent(in) :: comm
! Subprogram not used 
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS: 
! Subprogram not used !
! Subprogram not used    type(SparseMatrix), intent(inout) :: sMat
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used    integer, optional,  intent(out) :: stat
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY: 
! Subprogram not used !
! Subprogram not used ! 13Apr01 - J.W. Larson <larson@mcs.anl.gov> - initial API spec/code
! Subprogram not used ! 13Jun01 - J.W. Larson <larson@mcs.anl.gov> - Made status flag stat
! Subprogram not used !           optional, and ititilaze it to zero if it is present.
! Subprogram not used ! 17Jul02 - J.W. Larson <larson@mcs.anl.gov> - Bug fix--local 
! Subprogram not used !           process ID myID was uninitialized.
! Subprogram not used !EOP
! Subprogram not used !-------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'Bcast_'
! Subprogram not used 
! Subprogram not used ! Storage for the number of rows and columns in the SparseMatrix
! Subprogram not used   integer :: NumRowsColumns(2)
! Subprogram not used ! Process ID number
! Subprogram not used   integer :: myID
! Subprogram not used ! Error flag
! Subprogram not used   integer :: ierr
! Subprogram not used 
! Subprogram not used        ! Initialize stat if present
! Subprogram not used 
! Subprogram not used   if(present(stat)) stat = 0
! Subprogram not used 
! Subprogram not used        ! Determine local process ID myID:
! Subprogram not used 
! Subprogram not used   call MPI_COMM_RANK(comm, myID, ierr)
! Subprogram not used   if(ierr /= 0) then
! Subprogram not used      call MP_perr_die(myname_,"MPI_COMM_RANK() failed",ierr)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   if((myID.eq.root) .and. sMat%vecinit) then
! Subprogram not used       write(stderr,*) myname_,&
! Subprogram not used       "Cannot broadcast SparseMatrix with initialized vector parts." 
! Subprogram not used       call die(myname_,"Gather SparseMatrix with vecinit TRUE.") 
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! Broadcast sMat%data from the root
! Subprogram not used 
! Subprogram not used   call AttrVect_bcast(sMat%data, root, comm, ierr)
! Subprogram not used   if(ierr /= 0) then
! Subprogram not used      if(present(stat)) then
! Subprogram not used         write(stderr,*) myname_,"::  AttrVect_bcast(sMat%data...failed--stat=", &
! Subprogram not used              ierr
! Subprogram not used         stat = ierr
! Subprogram not used         return
! Subprogram not used      else
! Subprogram not used         call die(myname_,"call AttrVect_bcast(sMat%data...) failed",ierr)
! Subprogram not used      endif
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   if(myID == root) then
! Subprogram not used      NumRowsColumns(1) = SparseMatrix_nRows(sMat)
! Subprogram not used      NumRowsColumns(2) = SparseMatrix_nCols(sMat)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   call MPI_Bcast(NumRowsColumns, 2, MP_INTEGER, root, comm, ierr)
! Subprogram not used   if(ierr /= 0) then
! Subprogram not used      call MP_perr_die(myname_,"MPI_Bcast(NumRowsColumns...",ierr)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! Unpack NumRowsColumns on broadcast destination processes
! Subprogram not used 
! Subprogram not used   if(myID /= root) then
! Subprogram not used      sMat%nrows = NumRowsColumns(1)
! Subprogram not used      sMat%ncols = NumRowsColumns(2)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   sMat%vecinit = .FALSE.
! Subprogram not used 
! Subprogram not used  end subroutine Bcast_

 end module m_SparseMatrixComms
