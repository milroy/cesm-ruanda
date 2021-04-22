!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!-----------------------------------------------------------------------
! CVS $Id$
! CVS $Name$ 
!BOP -------------------------------------------------------------------
!
! !MODULE: m_SparseMatrix -- Sparse Matrix Object
!
! !DESCRIPTION:
! The {\tt SparseMatrix} data type is MCT's object for storing sparse 
! matrices.  In MCT, intergrid interpolation is implemented as a sparse 
! matrix-vector multiplication, with the {\tt AttrVect} type playing the 
! roles of the input and output vectors.  The interpolation matrices tend
! to be {\em extremely} sparse.  For ${\bf x} \in \Re^{N_x}$, and
! ${\bf y} \in \Re^{N_y}$, the interpolation matrix {\bf M} used to effect 
! ${\bf y} = {\bf M} {\bf x}$ will typically have ${\cal O}({N_y})$ 
! non-zero elements.  For that reason, the {\tt SparseMatrix} type 
! stores {\em only} information about non-zero matrix elements, along 
! with the number of rows and columns in the full matrix.  The nonzero  
! matrix elements are stored in {\tt AttrVect} form (see the module 
! {\tt m\_AttrVect} for more details), and the set of attributes are 
! listed below:
!
!\begin{table}[htbp]
!\begin{center}
!\begin{tabular}{|l|l|l|}
!\hline
!{\bf Attribute Name} & {\bf Significance} & {\tt Type} \!\hline

!{\tt grow} & Global Row Index & {\tt INTEGER} \!\hline

!{\tt gcol} & Global Column Index & {\tt INTEGER} \!\hline

!{\tt lrow} & Local Row Index & {\tt INTEGER} \!\hline

!{\tt lcol} & Local Column Index & {\tt INTEGER} \!\hline

!{\tt weight} & Matrix Element ${M_{ij}}$ & {\tt REAL} \!\hline

!\end{tabular}
!\end{center}
!\end{table}
! 
! The provision of both local and global column and row indices is 
! made because this datatype can be used in either shared-memory or 
! distributed-memory parallel matrix-vector products.
!
! This module contains the definition of the {\tt SparseMatrix} type, 
! creation and destruction methods, a variety of accessor methods, 
! routines for testing the suitability of the matrix for interpolation 
! (i.e. the sum of each row is either zero or unity), and methods for 
! sorting and permuting matrix entries.
!
! For better performance of the Matrix-Vector multiply on vector
! architectures, the {\tt SparseMatrix} object also contains arrays
! for holding the sparse matrix data in a more vector-friendly form.
!
!
! !INTERFACE:

 module m_SparseMatrix
!
! !USES:
!
      use m_realkinds, only : FP
      use m_AttrVect, only : AttrVect


      private   ! except

! !PUBLIC TYPES:

      public :: SparseMatrix      ! The class data structure

      Type SparseMatrix



     integer :: nrows
	 integer :: ncols
	 type(AttrVect) :: data
         logical :: vecinit       ! additional data for the vectorized sMat
         integer,dimension(:),pointer :: row_s, row_e
         integer, dimension(:,:), pointer :: tcol
         real(FP), dimension(:,:), pointer :: twgt
         integer :: row_max, row_min
         integer :: tbl_end
      End Type SparseMatrix

! !PUBLIC MEMBER FUNCTIONS:

      public :: init              ! Create a SparseMatrix
      public :: vecinit           ! Initialize the vector parts
      public :: clean             ! Destroy a SparseMatrix
      public :: lsize             ! Local number of elements
      public :: indexIA           ! Index integer attribute
      public :: indexRA           ! Index real attribute
      public :: nRows             ! Total number of rows
      public :: nCols             ! Total number of columns

      public :: exportGlobalRowIndices    ! Return global row indices 
                                          ! for matrix elements
      public :: exportGlobalColumnIndices ! Return global column indices 
                                          ! for matrix elements
      public :: exportLocalRowIndices     ! Return local row indices 
                                          ! for matrix elements
      public :: exportLocalColumnIndices  ! Return local column indices 
                                          ! for matrix elements
      public :: exportMatrixElements      ! Return matrix elements

      public :: importGlobalRowIndices    ! Set global row indices 
                                          ! using 
      public :: importGlobalColumnIndices ! Return global column indices 
                                          ! for matrix elements
      public :: importLocalRowIndices     ! Return local row indices 
                                          ! for matrix elements
      public :: importLocalColumnIndices  ! Return local column indices 
                                          ! for matrix elements
      public :: importMatrixElements      ! Return matrix elements
      public :: Copy                      ! Copy a SparseMatrix

      public :: GlobalNumElements ! Total number of nonzero elements
      public :: ComputeSparsity   ! Fraction of matrix that is nonzero
      public :: local_row_range   ! Local (on-process) row range
      public :: global_row_range  ! Local (on-process) row range
      public :: local_col_range   ! Local (on-process) column range
      public :: global_col_range  ! Local (on-process) column range
      public :: CheckBounds       ! Check row and column values
                                  ! for out-of-bounds values
      public :: row_sum           ! Return SparseMatrix row sums
      public :: row_sum_check     ! Check SparseMatrix row sums against
                                  ! input "valid" values
      public :: Sort              ! Sort matrix entries to generate an
                                  ! index permutation (to be used by
                                  ! Permute()
      public :: Permute           ! Permute matrix entries using index
                                  ! permutation gernerated by Sort()
      public :: SortPermute       ! Sort/Permute matrix entries

    interface init  ; module procedure init_  ; end interface
    interface vecinit  ; module procedure vecinit_  ; end interface
    interface clean ; module procedure clean_ ; end interface
    interface lsize ; module procedure lsize_ ; end interface
    interface indexIA ; module procedure indexIA_ ; end interface
    interface indexRA ; module procedure indexRA_ ; end interface
    interface nRows ; module procedure nRows_ ; end interface
    interface nCols ; module procedure nCols_ ; end interface

    interface exportGlobalRowIndices ; module procedure &
	 exportGlobalRowIndices_ 
    end interface

    interface exportGlobalColumnIndices ; module procedure &
	 exportGlobalColumnIndices_ 
    end interface

    interface exportLocalRowIndices ; module procedure &
	 exportLocalRowIndices_ 
    end interface

    interface exportLocalColumnIndices ; module procedure &
	 exportLocalColumnIndices_ 
    end interface

    interface exportMatrixElements ; module procedure &
	 exportMatrixElementsSP_, &
	 exportMatrixElementsDP_
    end interface

    interface importGlobalRowIndices ; module procedure &
	 importGlobalRowIndices_ 
    end interface

    interface importGlobalColumnIndices ; module procedure &
	 importGlobalColumnIndices_ 
    end interface

    interface importLocalRowIndices ; module procedure &
	 importLocalRowIndices_ 
    end interface

    interface importLocalColumnIndices ; module procedure &
	 importLocalColumnIndices_ 
    end interface

    interface importMatrixElements ; module procedure &
	 importMatrixElementsSP_, & 
	 importMatrixElementsDP_
    end interface

    interface Copy ; module procedure Copy_ ; end interface

    interface GlobalNumElements ; module procedure &
	 GlobalNumElements_ 
    end interface

    interface ComputeSparsity ; module procedure &
	 ComputeSparsitySP_,  &
	 ComputeSparsityDP_ 
    end interface

    interface local_row_range ; module procedure &
	 local_row_range_ 
    end interface

    interface global_row_range ; module procedure &
	 global_row_range_ 
    end interface

    interface local_col_range ; module procedure &
	 local_col_range_ 
    end interface

    interface global_col_range ; module procedure &
	 global_col_range_ 
    end interface

    interface CheckBounds; module procedure &
	 CheckBounds_ 
    end interface

    interface row_sum ; module procedure &
	 row_sumSP_, &
	 row_sumDP_
    end interface

    interface row_sum_check ; module procedure &
	 row_sum_checkSP_, & 
	 row_sum_checkDP_
    end interface

    interface Sort ; module procedure Sort_ ; end interface
    interface Permute ; module procedure Permute_ ; end interface
    interface SortPermute ; module procedure SortPermute_ ; end interface

! !REVISION HISTORY:
! 19Sep00 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
! 15Jan01 - J.W. Larson <larson@mcs.anl.gov> - added numerous APIs
! 25Feb01 - J.W. Larson <larson@mcs.anl.gov> - changed from row/column
!           attributes to global and local row and column attributes
! 23Apr01 - J.W. Larson <larson@mcs.anl.gov> - added number of rows
!           and columns to the SparseMatrix type.  This means the
!           SparseMatrix is no longer a straight AttrVect type.  This
!           also made necessary the addition of lsize(), indexIA(),
!           and indexRA().
! 29Oct03 - R. Jacob <jacob@mcs.anl.gov> - extend the SparseMatrix type
!           to include mods from Fujitsu for a vector-friendly MatVecMul
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='MCT::m_SparseMatrix'

! SparseMatrix_iList components:
  character(len=*),parameter :: SparseMatrix_iList='grow:gcol:lrow:lcol'
  integer,parameter :: SparseMatrix_igrow=1
  integer,parameter :: SparseMatrix_igcol=2
  integer,parameter :: SparseMatrix_ilrow=3
  integer,parameter :: SparseMatrix_ilcol=4

! SparseMatrix_rList components:
  character(len=*),parameter :: SparseMatrix_rList='weight'
  integer,parameter :: SparseMatrix_iweight=1

 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - Initialize an Empty SparseMatrix
!
! !DESCRIPTION:  This routine creates the storage space for the
! entries of a {\tt SparseMatrix}, and sets the number of rows and
! columns in it.  The input {\tt INTEGER} arguments {\tt nrows} and 
! {\tt ncols} specify the number of rows and columns respectively.
! The optional input argument {\tt lsize} specifies the number of 
! nonzero entries in the {\tt SparseMatrix}.  The initialized 
! {\tt SparseMatrix} is returned in the output argument {\tt sMat}.
!
! {\bf N.B.}:  This routine is allocating dynamical memory in the form
! of a {\tt SparseMatrix}.  The user must deallocate this space when
! the {\tt SparseMatrix} is no longer needed by invoking the routine
! {\tt clean\_()}.
!
! !INTERFACE:

 subroutine init_(sMat, nrows, ncols, lsize)
!
! !USES:
!
      use m_AttrVect, only : AttrVect_init => init
      use m_die

      implicit none

! !INPUT PARAMETERS:

      integer,            intent(in)   :: nrows
      integer,            intent(in)   :: ncols
      integer, optional,  intent(in)   :: lsize

! !OUTPUT PARAMETERS:

      type(SparseMatrix), intent(out)  :: sMat

! !REVISION HISTORY:
! 19Sep00 - Jay Larson <larson@mcs.anl.gov> - initial prototype
! 23Apr01 - Jay Larson <larson@mcs.anl.gov> - added arguments
!           nrows and ncols--number of rows and columns in the
!           SparseMatrix
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::init_'

  integer :: n

        ! if lsize is present, use it to set n; if not, set n=0

  n = 0
  if(present(lsize)) n=lsize

        ! Initialize number of rows and columns:

  sMat%nrows = nrows
  sMat%ncols = ncols

        ! Initialize sMat%data using AttrVect_init

  call AttrVect_init(sMat%data, SparseMatrix_iList, &
                     SparseMatrix_rList, n)

  ! vecinit is off by default
  sMat%vecinit = .FALSE.

 end subroutine init_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: vecinit_ - Initialize vector parts of a SparseMatrix
!
! !DESCRIPTION:  This routine creates the storage space for 
! and intializes the vector parts of a {\tt SparseMatrix}.
!
! {\bf N.B.}:  This routine assumes the locally indexed parts of a
! {\tt SparseMatrix} have been initialized.  This is
! accomplished by either importing the values directly with
! {\tt importLocalRowIndices} and {\tt importLocalColIndices} or by
! importing the Global Row and Col Indices and making two calls to 
! {\tt GlobalToLocalMatrix}.
!
! {\bf N.B.}:   The vector portion can use a large amount of
! memory so it is highly recommended that this routine only
! be called on a {\tt SparseMatrix} that has been scattered
! or otherwise sized locally.
!
! !INTERFACE:

! Subprogram not used  subroutine vecinit_(sMat)
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_die
! Subprogram not used       use m_stdio
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used       type(SparseMatrix), intent(inout)  :: sMat
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 27Oct03 - R. Jacob <jacob@mcs.anl.gov> - initial version
! Subprogram not used !           using code provided by Yoshi et. al.
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used !
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::vecinit_'
! Subprogram not used 
! Subprogram not used   integer :: irow,icol,iwgt
! Subprogram not used   integer :: num_elements
! Subprogram not used   integer :: row,col
! Subprogram not used   integer :: ier,l,n
! Subprogram not used   integer, dimension(:)  , allocatable :: nr, rn
! Subprogram not used 
! Subprogram not used   if(sMat%vecinit) then
! Subprogram not used    write(stderr,'(2a)') myname_, &
! Subprogram not used      'MCTERROR:  sMat vector parts have already been initialized...Continuing'
! Subprogram not used      RETURN
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   write(6,*) myname_,'Initializing vecMat'
! Subprogram not used   irow = indexIA_(sMat,'lrow',dieWith=myname_)
! Subprogram not used   icol = indexIA_(sMat,'lcol',dieWith=myname_)
! Subprogram not used   iwgt = indexRA_(sMat,'weight',dieWith=myname_)
! Subprogram not used 
! Subprogram not used   num_elements = lsize_(sMat)
! Subprogram not used 
! Subprogram not used   sMat%row_min = sMat%data%iAttr(irow,1)
! Subprogram not used   sMat%row_max = sMat%row_min
! Subprogram not used   do n=1,num_elements
! Subprogram not used      row = sMat%data%iAttr(irow,n)
! Subprogram not used      if ( row > sMat%row_max ) sMat%row_max = row
! Subprogram not used      if ( row < sMat%row_min ) sMat%row_min = row
! Subprogram not used   enddo
! Subprogram not used 
! Subprogram not used   allocate( nr(sMat%row_max), rn(num_elements), stat=ier)
! Subprogram not used   if(ier/=0) call die(myname_,'allocate(nr,rn)',ier)
! Subprogram not used 
! Subprogram not used   sMat%tbl_end = 0
! Subprogram not used   nr(:) = 0
! Subprogram not used   do n=1,num_elements
! Subprogram not used      row = sMat%data%iAttr(irow,n)
! Subprogram not used      nr(row) = nr(row)+1
! Subprogram not used      rn(n)   = nr(row)
! Subprogram not used   enddo
! Subprogram not used   sMat%tbl_end = maxval(rn)
! Subprogram not used 
! Subprogram not used   allocate( sMat%tcol(sMat%row_max,sMat%tbl_end),  &
! Subprogram not used             sMat%twgt(sMat%row_max,sMat%tbl_end), stat=ier )
! Subprogram not used   if(ier/=0) call die(myname_,'allocate(tcol,twgt)',ier)
! Subprogram not used 
! Subprogram not used !CDIR COLLAPSE
! Subprogram not used   sMat%tcol(:,:) = -1
! Subprogram not used   do n=1,num_elements
! Subprogram not used      row = sMat%data%iAttr(irow,n)
! Subprogram not used      sMat%tcol(row,rn(n)) = sMat%data%iAttr(icol,n)
! Subprogram not used      sMat%twgt(row,rn(n)) = sMat%data%rAttr(iwgt,n)
! Subprogram not used   enddo
! Subprogram not used 
! Subprogram not used   allocate( sMat%row_s(sMat%tbl_end) , sMat%row_e(sMat%tbl_end), &
! Subprogram not used               stat=ier )
! Subprogram not used   if(ier/=0) call die(myname_,'allocate(row_s,row_e',ier)
! Subprogram not used   sMat%row_s = sMat%row_min
! Subprogram not used   sMat%row_e = sMat%row_max
! Subprogram not used   do l=1,sMat%tbl_end
! Subprogram not used     do n=sMat%row_min,sMat%row_max
! Subprogram not used       if (nr(n) >= l) then
! Subprogram not used         sMat%row_s(l) = n
! Subprogram not used         exit
! Subprogram not used       endif
! Subprogram not used     enddo
! Subprogram not used     do n = sMat%row_max,sMat%row_min,-1
! Subprogram not used       if (nr(n) >= l) then
! Subprogram not used         sMat%row_e(l) = n
! Subprogram not used         exit
! Subprogram not used       endif
! Subprogram not used     enddo
! Subprogram not used   enddo
! Subprogram not used 
! Subprogram not used   deallocate(nr,rn, stat=ier)
! Subprogram not used   if(ier/=0) call die(myname_,'deallocate()',ier)
! Subprogram not used 
! Subprogram not used   sMat%vecinit = .TRUE.
! Subprogram not used 
! Subprogram not used  end subroutine vecinit_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - Destroy a SparseMatrix.
!
! !DESCRIPTION:  This routine deallocates dynamical memory held by the
! input {\tt SparseMatrix} argument {\tt sMat}.  It also sets the number
! of rows and columns in the {\tt SparseMatrix} to zero.
!
! !INTERFACE:

    subroutine clean_(sMat,stat)
!
! !USES:
!
      use m_AttrVect,only : AttrVect_clean => clean
      use m_die

      implicit none

! !INPUT/OUTPTU PARAMETERS:

      type(SparseMatrix), intent(inout) :: sMat

! !OUTPUT PARAMETERS:

      integer, optional,  intent(out)   :: stat

! !REVISION HISTORY:
! 19Sep00 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
! 23Apr00 - J.W. Larson <larson@mcs.anl.gov> - added changes to
!           accomodate clearing nrows and ncols.
! 01Mar02 - E.T. Ong <eong@mcs.anl.gov> Added stat argument.
! 03Oct03 - R. Jacob <jacob@mcs.anl.gov> - clean vector parts
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ier

       ! Deallocate memory held by sMat:

  if(present(stat)) then
     call AttrVect_clean(sMat%data,stat)
  else
     call AttrVect_clean(sMat%data)
  endif

       ! Set the number of rows and columns in sMat to zero:

  sMat%nrows = 0
  sMat%ncols = 0

  if(sMat%vecinit) then
    sMat%row_max = 0
    sMat%row_min = 0
    sMat%tbl_end = 0
    deallocate(sMat%row_s,sMat%row_e,stat=ier)
    if(ier/=0) then
      if(present(stat)) then
        stat=ier
      else
        call warn(myname_,'deallocate(row_s,row_e)',ier)
      endif
    endif

    deallocate(sMat%tcol,sMat%twgt,stat=ier)
    if(ier/=0) then
      if(present(stat)) then
        stat=ier
      else
        call warn(myname_,'deallocate(tcol,twgt)',ier)
      endif
    endif
    sMat%vecinit = .FALSE.
  endif
   

 end subroutine clean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: lsize_ - Local Number Non-zero Elements
!
! !DESCRIPTION:  This {\tt INTEGER} function reports on-processor storage 
! of the number of nonzero elements in the input {\tt SparseMatrix} 
! argument {\tt sMat}.  
!
! !INTERFACE:

    integer function lsize_(sMat)
!
! !USES:
!
      use m_AttrVect,only : AttrVect_lsize => lsize

      implicit none

! !INPUT PARAMETERS:

      type(SparseMatrix), intent(in) :: sMat

! !REVISION HISTORY:
! 23Apr00 - J.W. Larson <larson@mcs.anl.gov> - initial version.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::lsize_'

  lsize_ = AttrVect_lsize(sMat%data)

 end function lsize_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE:  GlobalNumElements_ - Global Number of Non-zero Elements
!
! !DESCRIPTION:  This routine computes the number of nonzero elements 
! in a distributed {\tt SparseMatrix} variable {\tt sMat}.  The input 
! {\tt SparseMatrix} argument {\tt sMat} is examined on each process 
! to determine the number of nonzero elements it holds, and this value 
! is summed across the communicator associated with the input 
! {\tt INTEGER} handle {\tt comm}, with the total returned {\em on each
! process on the communicator}.
!
! !INTERFACE:

 integer function GlobalNumElements_(sMat, comm)

!
! !USES:
!
      use m_die
      use m_mpif90

      implicit none

! !INPUT PARAMETERS: 

      type(SparseMatrix), intent(in)  :: sMat
      integer, optional,  intent(in)  :: comm

! !REVISION HISTORY:
! 24Apr01 - Jay Larson <larson@mcs.anl.gov> - New routine.
!
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//':: GlobalNumElements_'

  integer :: MyNumElements, GNumElements, ierr

       ! Determine the number of locally held nonzero elements:

  MyNumElements = lsize_(sMat)

  call MPI_ALLREDUCE(MyNumElements, GNumElements, 1, MP_INTEGER, &
                     MP_SUM, comm, ierr)
  if(ierr /= 0) then
     call MP_perr_die(myname_,"MPI_ALLREDUCE(MyNumElements...",ierr)
  endif

 GlobalNumElements_ = GNumElements

 end function GlobalNumElements_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexIA_ - Index an Integer Attribute
!
! !DESCRIPTION:  This {\tt INTEGER} function reports the row index 
! for a given {\tt INTEGER} attribute of the input {\tt SparseMatrix} 
! argument {\tt sMat}.  The attribute requested is represented by the 
! input {\tt CHARACTER} variable {\tt attribute}.  The list of integer 
! attributes one can request is defined in the description block of the 
! header of this module ({\tt m\_SparseMatrix}).
!
! Here is how {\tt indexIA\_} provides access to integer attribute data
! in a {\tt SparseMatrix} variable {\tt sMat}.  Suppose we wish to access
! global row information.  This attribute has associated with it the 
! string tag {\tt grow}.  The corresponding index returned ({\tt igrow}) 
! is determined by invoking {\tt indexIA\_}:
! \begin{verbatim}
! igrow = indexIA_(sMat, 'grow')
! \end{verbatim}
!
! Access to the global row index data in {\tt sMat} is thus obtained by 
! referencing {\tt sMat\%data\%iAttr(igrow,:)}.
!
!
! !INTERFACE:

    integer function indexIA_(sMat, item, perrWith, dieWith)
!
! !USES:
!
      use m_String, only : String
      use m_String, only : String_init => init
      use m_String, only : String_clean => clean
      use m_String, only : String_ToChar => ToChar

      use m_TraceBack, only : GenTraceBackString

      use m_AttrVect,only : AttrVect_indexIA => indexIA

      implicit none

! !INPUT PARAMETERS:

      type(SparseMatrix),         intent(in) :: sMat
      character(len=*),           intent(in) :: item
      character(len=*), optional, intent(in) :: perrWith
      character(len=*), optional, intent(in) :: dieWith

! !REVISION HISTORY:
!       23Apr00 - J.W. Larson <larson@mcs.anl.gov> - initial version.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::indexIA_'
  type(String) :: myTrace

       ! Generate a traceback String

  if(present(dieWith)) then
     call GenTraceBackString(myTrace, dieWith, myname_)
  else
     if(present(perrWith)) then
        call GenTraceBackString(myTrace, perrWith, myname_)
     else
        call GenTraceBackString(myTrace, myname_)
     endif
  endif

       ! Call AttrVect_indexIA() accordingly:

  if( present(dieWith) .or. &
     ((.not. present(dieWith)) .and. (.not. present(perrWith))) ) then
     indexIA_ = AttrVect_indexIA(sMat%data, item, &
                                 dieWith=String_ToChar(myTrace))
  else  ! perrWith but no dieWith case
     indexIA_ = AttrVect_indexIA(sMat%data, item, &
                   perrWith=String_ToChar(myTrace))
  endif

  call String_clean(myTrace)

 end function indexIA_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexRA_ - Index a Real Attribute
!
! !DESCRIPTION:  This {\tt INTEGER} function reports the row index 
! for a given {\tt REAL} attribute of the input {\tt SparseMatrix} 
! argument {\tt sMat}.  The attribute requested is represented by the 
! input {\tt CHARACTER} variable {\tt attribute}.  The list of real 
! attributes one can request is defined in the description block of the 
! header of this module ({\tt m\_SparseMatrix}).
!
! Here is how {\tt indexRA\_} provides access to integer attribute data
! in a {\tt SparseMatrix} variable {\tt sMat}.  Suppose we wish to access
! matrix element values.  This attribute has associated with it the 
! string tag {\tt weight}.  The corresponding index returned ({\tt iweight}) 
! is determined by invoking {\tt indexRA\_}:
! \begin{verbatim}
! iweight = indexRA_(sMat, 'weight')
! \end{verbatim}
!
! Access to the matrix element data in {\tt sMat} is thus obtained by 
! referencing {\tt sMat\%data\%rAttr(iweight,:)}.
!
! !INTERFACE:

    integer function indexRA_(sMat, item, perrWith, dieWith)
!
! !USES:
!
      use m_String, only : String
      use m_String, only : String_init => init
      use m_String, only : String_clean => clean
      use m_String, only : String_ToChar => ToChar

      use m_TraceBack, only : GenTraceBackString

      use m_AttrVect,only : AttrVect_indexRA => indexRA

      implicit none

! !INPUT PARAMETERS:

      type(SparseMatrix),         intent(in) :: sMat
      character(len=*),           intent(in) :: item
      character(len=*), optional, intent(in) :: perrWith
      character(len=*), optional, intent(in) :: dieWith

! !REVISION HISTORY:
! 24Apr00 - J.W. Larson <larson@mcs.anl.gov> - initial version.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::indexRA_'

  type(String) :: myTrace

       ! Generate a traceback String

  if(present(dieWith)) then ! append myname_ onto dieWith
     call GenTraceBackString(myTrace, dieWith, myname_)
  else
     if(present(perrWith)) then ! append myname_ onto perrwith
        call GenTraceBackString(myTrace, perrWith, myname_)
     else ! Start a TraceBack String
        call GenTraceBackString(myTrace, myname_)
     endif
  endif

       ! Call AttrVect_indexRA() accordingly:

  if( present(dieWith) .or. &
     ((.not. present(dieWith)) .and. (.not. present(perrWith))) ) then
     indexRA_ = AttrVect_indexRA(sMat%data, item, &
                                 dieWith=String_ToChar(myTrace))
  else  ! perrWith but no dieWith case
     indexRA_ = AttrVect_indexRA(sMat%data, item, &
                   perrWith=String_ToChar(myTrace))
  endif

  call String_clean(myTrace)

 end function indexRA_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: nRows_ - Return the Number of Rows
!
! !DESCRIPTION:  This routine returns the {\em total} number of rows
! in the input {\tt SparseMatrix} argument {\tt sMat}.  This number of
! rows is a constant, and not dependent on the decomposition of the 
! {\tt SparseMatrix}.
!
! !INTERFACE:

    integer function nRows_(sMat)
!
! !USES:
!
      implicit none

! !INPUT PARAMETERS:

      type(SparseMatrix), intent(in) :: sMat

! !REVISION HISTORY:
! 19Apr01 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::nRows_'

  nRows_ = sMat%nrows

 end function nRows_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: nCols_ - Return the Number of Columns
!
! !DESCRIPTION:  This routine returns the {\em total} number of columns
! in the input {\tt SparseMatrix} argument {\tt sMat}.  This number of
! columns is a constant, and not dependent on the decomposition of the 
! {\tt SparseMatrix}.
!
! !INTERFACE:

    integer function nCols_(sMat)
!
! !USES:
!
      implicit none

! !INPUT PARAMETERS:

      type(SparseMatrix), intent(in) :: sMat

! !REVISION HISTORY:
! 19Apr01 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::nCols_'

  nCols_ = sMat%ncols

 end function nCols_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: exportGlobalRowIndices_ - Return Global Row Indices
!
! !DESCRIPTION:
! This routine extracts from the input {\tt SparseMatrix} argument 
! {\tt sMat} its global row indices, and returns them in the {\tt INTEGER} 
! output array {\tt GlobalRows}, and its length in the output {\tt INTEGER} 
! argument {\tt length}.
!
! {\bf N.B.:}  The flexibility of this routine regarding the pointer 
! association status of the output argument {\tt GlobalRows} means the
! user must invoke this routine with care.  If the user wishes this
! routine to fill a pre-allocated array, then obviously this array
! must be allocated prior to calling this routine.  If the user wishes
! that the routine {\em create} the output argument array {\tt GlobalRows},
! then the user must ensure this pointer is not allocated (i.e. the user
! must nullify this pointer) at the time this routine is invoked.
!
! {\bf N.B.:}  If the user has relied on this routine to allocate memory
! associated with the pointer {\tt GlobalRows}, then the user is responsible 
! for deallocating this array once it is no longer needed.  Failure to 
! do so will result in a memory leak.
!
! !INTERFACE:

! Subprogram not used  subroutine exportGlobalRowIndices_(sMat, GlobalRows, length)
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_die 
! Subprogram not used       use m_stdio
! Subprogram not used 
! Subprogram not used       use m_AttrVect,      only : AttrVect_exportIAttr => exportIAttr
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS: 
! Subprogram not used 
! Subprogram not used       type(SparseMatrix),     intent(in)  :: sMat
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS: 
! Subprogram not used 
! Subprogram not used       integer,  dimension(:), pointer     :: GlobalRows
! Subprogram not used       integer,  optional,     intent(out) :: length
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  7May02 - J.W. Larson <larson@mcs.anl.gov> - initial version.
! Subprogram not used !
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::exportGlobalRowIndices_'
! Subprogram not used 
! Subprogram not used        ! Export the data (inheritance from AttrVect)
! Subprogram not used   if(present(length)) then
! Subprogram not used      call AttrVect_exportIAttr(sMat%data, 'grow', GlobalRows, length)
! Subprogram not used   else
! Subprogram not used      call AttrVect_exportIAttr(sMat%data, 'grow', GlobalRows)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used  end subroutine exportGlobalRowIndices_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: exportGlobalColumnIndices_ - Return Global Column Indices
!
! !DESCRIPTION:
! This routine extracts from the input {\tt SparseMatrix} argument 
! {\tt sMat} its global column indices, and returns them in the {\tt INTEGER} 
! output array {\tt GlobalColumns}, and its length in the output {\tt INTEGER} 
! argument {\tt length}.
!
! {\bf N.B.:}  The flexibility of this routine regarding the pointer 
! association status of the output argument {\tt GlobalColumns} means the
! user must invoke this routine with care.  If the user wishes this
! routine to fill a pre-allocated array, then obviously this array
! must be allocated prior to calling this routine.  If the user wishes
! that the routine {\em create} the output argument array {\tt GlobalColumns},
! then the user must ensure this pointer is not allocated (i.e. the user
! must nullify this pointer) at the time this routine is invoked.
!
! {\bf N.B.:}  If the user has relied on this routine to allocate memory
! associated with the pointer {\tt GlobalColumns}, then the user is responsible 
! for deallocating this array once it is no longer needed.  Failure to 
! do so will result in a memory leak.
!
! !INTERFACE:

! Subprogram not used  subroutine exportGlobalColumnIndices_(sMat, GlobalColumns, length)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_die 
! Subprogram not used       use m_stdio
! Subprogram not used 
! Subprogram not used       use m_AttrVect,      only : AttrVect_exportIAttr => exportIAttr
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS: 
! Subprogram not used 
! Subprogram not used       type(SparseMatrix),     intent(in)  :: sMat
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS: 
! Subprogram not used 
! Subprogram not used       integer,  dimension(:), pointer     :: GlobalColumns
! Subprogram not used       integer,  optional,     intent(out) :: length
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  7May02 - J.W. Larson <larson@mcs.anl.gov> - initial version.
! Subprogram not used !
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::exportGlobalColumnIndices_'
! Subprogram not used 
! Subprogram not used        ! Export the data (inheritance from AttrVect)
! Subprogram not used   if(present(length)) then
! Subprogram not used      call AttrVect_exportIAttr(sMat%data, 'gcol', GlobalColumns, length)
! Subprogram not used   else
! Subprogram not used      call AttrVect_exportIAttr(sMat%data, 'gcol', GlobalColumns)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used  end subroutine exportGlobalColumnIndices_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: exportLocalRowIndices_ - Return Local Row Indices
!
! !DESCRIPTION:
! This routine extracts from the input {\tt SparseMatrix} argument 
! {\tt sMat} its local row indices, and returns them in the {\tt INTEGER} 
! output array {\tt LocalRows}, and its length in the output {\tt INTEGER} 
! argument {\tt length}.
!
! {\bf N.B.:}  The flexibility of this routine regarding the pointer 
! association status of the output argument {\tt LocalRows} means the
! user must invoke this routine with care.  If the user wishes this
! routine to fill a pre-allocated array, then obviously this array
! must be allocated prior to calling this routine.  If the user wishes
! that the routine {\em create} the output argument array {\tt LocalRows},
! then the user must ensure this pointer is not allocated (i.e. the user
! must nullify this pointer) at the time this routine is invoked.
!
! {\bf N.B.:}  If the user has relied on this routine to allocate memory
! associated with the pointer {\tt LocalRows}, then the user is responsible 
! for deallocating this array once it is no longer needed.  Failure to 
! do so will result in a memory leak.
!
! !INTERFACE:

! Subprogram not used  subroutine exportLocalRowIndices_(sMat, LocalRows, length)
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_die 
! Subprogram not used       use m_stdio
! Subprogram not used 
! Subprogram not used       use m_AttrVect,      only : AttrVect_exportIAttr => exportIAttr
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS: 
! Subprogram not used 
! Subprogram not used       type(SparseMatrix),     intent(in)  :: sMat
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS: 
! Subprogram not used 
! Subprogram not used       integer,  dimension(:), pointer     :: LocalRows
! Subprogram not used       integer,  optional,     intent(out) :: length
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  7May02 - J.W. Larson <larson@mcs.anl.gov> - initial version.
! Subprogram not used !
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::exportLocalRowIndices_'
! Subprogram not used 
! Subprogram not used        ! Export the data (inheritance from AttrVect)
! Subprogram not used   if(present(length)) then
! Subprogram not used      call AttrVect_exportIAttr(sMat%data, 'lrow', LocalRows, length)
! Subprogram not used   else
! Subprogram not used      call AttrVect_exportIAttr(sMat%data, 'lrow', LocalRows)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used  end subroutine exportLocalRowIndices_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: exportLocalColumnIndices_ - Return Local Column Indices 
!
! !DESCRIPTION:
! This routine extracts from the input {\tt SparseMatrix} argument 
! {\tt sMat} its local column indices, and returns them in the {\tt INTEGER} 
! output array {\tt LocalColumns}, and its length in the output {\tt INTEGER} 
! argument {\tt length}.
!
! {\bf N.B.:}  The flexibility of this routine regarding the pointer 
! association status of the output argument {\tt LocalColumns} means the
! user must invoke this routine with care.  If the user wishes this
! routine to fill a pre-allocated array, then obviously this array
! must be allocated prior to calling this routine.  If the user wishes
! that the routine {\em create} the output argument array {\tt LocalColumns},
! then the user must ensure this pointer is not allocated (i.e. the user
! must nullify this pointer) at the time this routine is invoked.
!
! {\bf N.B.:}  If the user has relied on this routine to allocate memory
! associated with the pointer {\tt LocalColumns}, then the user is responsible 
! for deallocating this array once it is no longer needed.  Failure to 
! do so will result in a memory leak.
!
! !INTERFACE:

! Subprogram not used  subroutine exportLocalColumnIndices_(sMat, LocalColumns, length)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_die 
! Subprogram not used       use m_stdio
! Subprogram not used 
! Subprogram not used       use m_AttrVect,      only : AttrVect_exportIAttr => exportIAttr
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS: 
! Subprogram not used 
! Subprogram not used       type(SparseMatrix),     intent(in)  :: sMat
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS: 
! Subprogram not used 
! Subprogram not used       integer,  dimension(:), pointer     :: LocalColumns
! Subprogram not used       integer,  optional,     intent(out) :: length
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  7May02 - J.W. Larson <larson@mcs.anl.gov> - initial version.
! Subprogram not used !
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::exportLocalColumnIndices_'
! Subprogram not used 
! Subprogram not used        ! Export the data (inheritance from AttrVect)
! Subprogram not used   if(present(length)) then
! Subprogram not used      call AttrVect_exportIAttr(sMat%data, 'lcol', LocalColumns, length)
! Subprogram not used   else
! Subprogram not used      call AttrVect_exportIAttr(sMat%data, 'lcol', LocalColumns)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used  end subroutine exportLocalColumnIndices_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: exportMatrixElementsSP_ - Return Matrix Elements as Array
!
! !DESCRIPTION:
! This routine extracts the matrix elements from the input {\tt SparseMatrix} 
! argument {\tt sMat}, and returns them in the {\tt REAL} output array 
! {\tt MatrixElements}, and its length in the output {\tt INTEGER} 
! argument {\tt length}.
!
! {\bf N.B.:}  The flexibility of this routine regarding the pointer 
! association status of the output argument {\tt MatrixElements} means the
! user must invoke this routine with care.  If the user wishes this
! routine to fill a pre-allocated array, then obviously this array
! must be allocated prior to calling this routine.  If the user wishes
! that the routine {\em create} the output argument array {\tt MatrixElements},
! then the user must ensure this pointer is not allocated (i.e. the user
! must nullify this pointer) at the time this routine is invoked.
!
! {\bf N.B.:}  If the user has relied on this routine to allocate memory
! associated with the pointer {\tt MatrixElements}, then the user is responsible 
! for deallocating this array once it is no longer needed.  Failure to 
! do so will result in a memory leak.
!
! The native precision version is described here.  A double precision version
! is also available.
!
! !INTERFACE:

! Subprogram not used  subroutine exportMatrixelementsSP_(sMat, MatrixElements, length)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_die 
! Subprogram not used       use m_stdio
! Subprogram not used       use m_realkinds, only : SP
! Subprogram not used 
! Subprogram not used       use m_AttrVect,      only : AttrVect_exportRAttr => exportRAttr
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS: 
! Subprogram not used 
! Subprogram not used       type(SparseMatrix),     intent(in)  :: sMat
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS: 
! Subprogram not used 
! Subprogram not used       real(SP),  dimension(:),    pointer     :: MatrixElements
! Subprogram not used       integer,   optional,        intent(out) :: length
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  7May02 - J.W. Larson <larson@mcs.anl.gov> - initial version.
! Subprogram not used !  6Jan04 - R. Jacob <jacob@mcs.anl.gov> - SP and DP versions
! Subprogram not used !
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::exportMatrixElementsSP_'
! Subprogram not used 
! Subprogram not used        ! Export the data (inheritance from AttrVect)
! Subprogram not used   if(present(length)) then
! Subprogram not used      call AttrVect_exportRAttr(sMat%data, 'weight', MatrixElements, length)
! Subprogram not used   else
! Subprogram not used      call AttrVect_exportRAttr(sMat%data, 'weight', MatrixElements)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used  end subroutine exportMatrixElementsSP_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
! -------------------------------------------------------------------
!
! !IROUTINE: exportMatrixElementsDP_ - Return Matrix Elements as Array
!
! !DESCRIPTION:
! Double precision version of exportMatrixElementsSP_
!
! !INTERFACE:

! Subprogram not used  subroutine exportMatrixelementsDP_(sMat, MatrixElements, length)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_die 
! Subprogram not used       use m_stdio
! Subprogram not used       use m_realkinds, only : DP
! Subprogram not used 
! Subprogram not used       use m_AttrVect,      only : AttrVect_exportRAttr => exportRAttr
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS: 
! Subprogram not used 
! Subprogram not used       type(SparseMatrix),     intent(in)  :: sMat
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS: 
! Subprogram not used 
! Subprogram not used       real(DP),  dimension(:),    pointer     :: MatrixElements
! Subprogram not used       integer,   optional,        intent(out) :: length
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  7May02 - J.W. Larson <larson@mcs.anl.gov> - initial version.
! Subprogram not used !
! Subprogram not used ! ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::exportMatrixElementsDP_'
! Subprogram not used 
! Subprogram not used        ! Export the data (inheritance from AttrVect)
! Subprogram not used   if(present(length)) then
! Subprogram not used      call AttrVect_exportRAttr(sMat%data, 'weight', MatrixElements, length)
! Subprogram not used   else
! Subprogram not used      call AttrVect_exportRAttr(sMat%data, 'weight', MatrixElements)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used  end subroutine exportMatrixElementsDP_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: importGlobalRowIndices_ - Set Global Row Indices of Elements
!
! !DESCRIPTION:
! This routine imports global row index data into the {\tt SparseMatrix} 
! argument {\tt sMat}.  The user provides the index data in the input 
! {\tt INTEGER} vector {\tt inVect}.  The input {\tt INTEGER} argument 
! {\tt lsize} is used as a consistencey check to ensure the user is
! sufficient space in the {\tt SparseMatrix} to store the data.
!
! !INTERFACE:

 subroutine importGlobalRowIndices_(sMat, inVect, lsize)

!
! !USES:
!
      use m_die
      use m_stdio

      use m_AttrVect,      only : AttrVect_importIAttr => importIAttr

      implicit none

! !INPUT PARAMETERS: 

      integer,  dimension(:), pointer       :: inVect
      integer,                intent(in)    :: lsize

! !INPUT/OUTPUT PARAMETERS: 

      type(SparseMatrix),     intent(inout) :: sMat

! !REVISION HISTORY:
!  7May02 - J.W. Larson <larson@mcs.anl.gov> - initial prototype.
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::importGlobalRowIndices_'

       ! Argument Check:

  if(lsize > lsize_(sMat)) then
     write(stderr,*) myname_,':: ERROR, lsize > lsize_(sMat).', &
          'lsize = ',lsize,'lsize_(sMat) = ',lsize_(sMat)
     call die(myname_)
  endif

       ! Import the data (inheritance from AttrVect)

  call AttrVect_importIAttr(sMat%data, 'grow', inVect, lsize)

 end subroutine importGlobalRowIndices_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: importGlobalColumnIndices_ - Set Global Column Indices of Elements
!
! !DESCRIPTION:
! This routine imports global column index data into the {\tt SparseMatrix} 
! argument {\tt sMat}.  The user provides the index data in the input 
! {\tt INTEGER} vector {\tt inVect}.  The input {\tt INTEGER} argument 
! {\tt lsize} is used as a consistencey check to ensure the user is
! sufficient space in the {\tt SparseMatrix} to store the data.
!
! !INTERFACE:

 subroutine importGlobalColumnIndices_(sMat, inVect, lsize)

!
! !USES:
!
      use m_die
      use m_stdio

      use m_AttrVect,      only : AttrVect_importIAttr => importIAttr

      implicit none

! !INPUT PARAMETERS: 

      integer,  dimension(:), pointer       :: inVect
      integer,                intent(in)    :: lsize

! !INPUT/OUTPUT PARAMETERS: 

      type(SparseMatrix),     intent(inout) :: sMat

! !REVISION HISTORY:
!  7May02 - J.W. Larson <larson@mcs.anl.gov> - initial prototype.
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::importGlobalColumnIndices_'

       ! Argument Check:

  if(lsize > lsize_(sMat)) then
     write(stderr,*) myname_,':: ERROR, lsize > lsize_(sMat).', &
          'lsize = ',lsize,'lsize_(sMat) = ',lsize_(sMat)
     call die(myname_)
  endif

       ! Import the data (inheritance from AttrVect)

  call AttrVect_importIAttr(sMat%data, 'gcol', inVect, lsize)

 end subroutine importGlobalColumnIndices_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: importLocalRowIndices_ - Set Local Row Indices of Elements
!
! !DESCRIPTION:
! This routine imports local row index data into the {\tt SparseMatrix} 
! argument {\tt sMat}.  The user provides the index data in the input 
! {\tt INTEGER} vector {\tt inVect}.  The input {\tt INTEGER} argument 
! {\tt lsize} is used as a consistencey check to ensure the user is
! sufficient space in the {\tt SparseMatrix} to store the data.
!
! !INTERFACE:

! Subprogram not used  subroutine importLocalRowIndices_(sMat, inVect, lsize)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_die
! Subprogram not used       use m_stdio
! Subprogram not used 
! Subprogram not used       use m_AttrVect,      only : AttrVect_importIAttr => importIAttr
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS: 
! Subprogram not used 
! Subprogram not used       integer,  dimension(:), pointer       :: inVect
! Subprogram not used       integer,                intent(in)    :: lsize
! Subprogram not used 
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS: 
! Subprogram not used 
! Subprogram not used       type(SparseMatrix),     intent(inout) :: sMat
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  7May02 - J.W. Larson <larson@mcs.anl.gov> - initial prototype.
! Subprogram not used !
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::importLocalRowIndices_'
! Subprogram not used 
! Subprogram not used        ! Argument Check:
! Subprogram not used 
! Subprogram not used   if(lsize > lsize_(sMat)) then
! Subprogram not used      write(stderr,*) myname_,':: ERROR, lsize > lsize_(sMat).', &
! Subprogram not used           'lsize = ',lsize,'lsize_(sMat) = ',lsize_(sMat)
! Subprogram not used      call die(myname_)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! Import the data (inheritance from AttrVect)
! Subprogram not used 
! Subprogram not used   call AttrVect_importIAttr(sMat%data, 'lrow', inVect, lsize)
! Subprogram not used 
! Subprogram not used  end subroutine importLocalRowIndices_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: importLocalColumnIndices_ - Set Local Column Indices of Elements
!
! !DESCRIPTION:
! This routine imports local column index data into the {\tt SparseMatrix} 
! argument {\tt sMat}.  The user provides the index data in the input 
! {\tt INTEGER} vector {\tt inVect}.  The input {\tt INTEGER} argument 
! {\tt lsize} is used as a consistencey check to ensure the user is
! sufficient space in the {\tt SparseMatrix} to store the data.
!
! !INTERFACE:

! Subprogram not used  subroutine importLocalColumnIndices_(sMat, inVect, lsize)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_die
! Subprogram not used       use m_stdio
! Subprogram not used 
! Subprogram not used       use m_AttrVect,      only : AttrVect_importIAttr => importIAttr
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS: 
! Subprogram not used 
! Subprogram not used       integer,  dimension(:), pointer       :: inVect
! Subprogram not used       integer,                intent(in)    :: lsize
! Subprogram not used 
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS: 
! Subprogram not used 
! Subprogram not used       type(SparseMatrix),     intent(inout) :: sMat
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  7May02 - J.W. Larson <larson@mcs.anl.gov> - initial prototype.
! Subprogram not used !
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::importLocalColumnIndices_'
! Subprogram not used 
! Subprogram not used        ! Argument Check:
! Subprogram not used 
! Subprogram not used   if(lsize > lsize_(sMat)) then
! Subprogram not used      write(stderr,*) myname_,':: ERROR, lsize > lsize_(sMat).', &
! Subprogram not used           'lsize = ',lsize,'lsize_(sMat) = ',lsize_(sMat)
! Subprogram not used      call die(myname_)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! Import the data (inheritance from AttrVect)
! Subprogram not used 
! Subprogram not used   call AttrVect_importIAttr(sMat%data, 'lcol', inVect, lsize)
! Subprogram not used 
! Subprogram not used  end subroutine importLocalColumnIndices_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: importMatrixElementsSP_ - Import Non-zero Matrix Elements
!
! !DESCRIPTION:
! This routine imports matrix elements index data into the 
! {\tt SparseMatrix} argument {\tt sMat}.  The user provides the index 
! data in the input {\tt REAL} vector {\tt inVect}.  The input 
! {\tt INTEGER} argument {\tt lsize} is used as a consistencey check 
! to ensure the user is sufficient space in the {\tt SparseMatrix} 
! to store the data.
!
! !INTERFACE:

! Subprogram not used  subroutine importMatrixElementsSP_(sMat, inVect, lsize)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_die
! Subprogram not used       use m_stdio
! Subprogram not used       use m_realkinds, only : SP
! Subprogram not used 
! Subprogram not used       use m_AttrVect,      only : AttrVect_importRAttr => importRAttr
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS: 
! Subprogram not used 
! Subprogram not used       real(SP),  dimension(:),    pointer       :: inVect
! Subprogram not used       integer,                intent(in)    :: lsize
! Subprogram not used 
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS: 
! Subprogram not used 
! Subprogram not used       type(SparseMatrix),     intent(inout) :: sMat
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  7May02 - J.W. Larson <larson@mcs.anl.gov> - initial prototype.
! Subprogram not used !  6Jan04 - R. Jacob <jacob@mcs.anl.gov> - Make SP and DP versions.
! Subprogram not used !
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::importMatrixElementsSP_'
! Subprogram not used 
! Subprogram not used        ! Argument Check:
! Subprogram not used 
! Subprogram not used   if(lsize > lsize_(sMat)) then
! Subprogram not used      write(stderr,*) myname_,':: ERROR, lsize > lsize_(sMat).', &
! Subprogram not used           'lsize = ',lsize,'lsize_(sMat) = ',lsize_(sMat)
! Subprogram not used      call die(myname_)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! Import the data (inheritance from AttrVect)
! Subprogram not used 
! Subprogram not used   call AttrVect_importRAttr(sMat%data, 'weight', inVect, lsize)
! Subprogram not used 
! Subprogram not used  end subroutine importMatrixElementsSP_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
! -------------------------------------------------------------------
!
! !IROUTINE: importMatrixElementsDP_ - Import Non-zero Matrix Elements
!
! !DESCRIPTION:
! Double precision version of importMatrixElementsSP_
!
! !INTERFACE:

 subroutine importMatrixElementsDP_(sMat, inVect, lsize)

!
! !USES:
!
      use m_die
      use m_stdio
      use m_realkinds, only : DP

      use m_AttrVect,      only : AttrVect_importRAttr => importRAttr

      implicit none

! !INPUT PARAMETERS: 

      real(DP),  dimension(:),    pointer       :: inVect
      integer,                intent(in)    :: lsize

! !INPUT/OUTPUT PARAMETERS: 

      type(SparseMatrix),     intent(inout) :: sMat

! !REVISION HISTORY:
!  7May02 - J.W. Larson <larson@mcs.anl.gov> - initial prototype.
!
! ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::importMatrixElementsDP_'

       ! Argument Check:

  if(lsize > lsize_(sMat)) then
     write(stderr,*) myname_,':: ERROR, lsize > lsize_(sMat).', &
          'lsize = ',lsize,'lsize_(sMat) = ',lsize_(sMat)
     call die(myname_)
  endif

       ! Import the data (inheritance from AttrVect)

  call AttrVect_importRAttr(sMat%data, 'weight', inVect, lsize)

 end subroutine importMatrixElementsDP_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: Copy_ - Create a Copy of an Input SparseMatrix
!
! !DESCRIPTION:
! This routine creates a copy of the input {\tt SparseMatrix} argument 
! {\tt sMat}, returning it as the output {\tt SparseMatrix} argument 
! {\tt sMatCopy}.
!
! {\bf N.B.:}  The output argument {\tt sMatCopy} represents allocated 
! memory the user must deallocate when it is no longer needed.  The 
! MCT routine to use for this purpose is {\tt clean()} from this module.
!
! !INTERFACE:

 subroutine Copy_(sMat, sMatCopy)

!
! !USES:
!
      use m_die
      use m_stdio

      use m_AttrVect,      only : AttrVect
      use m_AttrVect,      only : AttrVect_init => init
      use m_AttrVect,      only : AttrVect_lsize => lsize
      use m_AttrVect,      only : AttrVect_Copy => Copy

      implicit none

! !INPUT PARAMETERS: 

      type(SparseMatrix), intent(in) :: sMat

! !OUTPUT PARAMETERS: 

      type(SparseMatrix), intent(out) :: sMatCopy

! !REVISION HISTORY:
! 27Sep02 - J.W. Larson <larson@mcs.anl.gov> - initial prototype.
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::Copy_'

       ! Step one:  copy the integer components of sMat:

  sMatCopy%nrows = sMat%nrows
  sMatCopy%ncols = sMat%ncols

  sMatCopy%vecinit = .FALSE.

       ! Step two:  Initialize the AttrVect sMatCopy%data off of sMat:

  call AttrVect_init(sMatCopy%data, sMat%data, AttrVect_lsize(sMat%data))

       ! Step three:  Copy sMat%data to sMatCopy%data:

  call AttrVect_Copy(sMat%data, aVout=sMatCopy%data)

  if(sMat%vecinit) call vecinit_(sMatCopy)

 end subroutine Copy_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: local_row_range_ - Local Row Extent of Non-zero Elements
!
! !DESCRIPTION: This routine examines the input distributed 
! {\tt SparseMatrix} variable {\tt sMat}, and returns the range of local 
! row values having nonzero elements.  The first local row with 
! nonzero elements is returned in the {\tt INTEGER} argument 
! {\tt start\_row}, the last row in {\tt end\_row}.
!
! !INTERFACE:

! Subprogram not used  subroutine local_row_range_(sMat, start_row, end_row)
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_die
! Subprogram not used 
! Subprogram not used       use m_AttrVect, only : AttrVect_lsize => lsize
! Subprogram not used       use m_AttrVect, only : AttrVect_indexIA => indexIA
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used       type(SparseMatrix), intent(in)  :: sMat
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used       integer,            intent(out) :: start_row
! Subprogram not used       integer,            intent(out) :: end_row
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 15Jan01 - Jay Larson <larson@mcs.anl.gov> - API specification.
! Subprogram not used ! 25Feb01 - Jay Larson <larson@mcs.anl.gov> - Initial prototype.
! Subprogram not used ! 23Apr01 - Jay Larson <larson@mcs.anl.gov> - Modified to accomodate
! Subprogram not used !           changes to the SparseMatrix type.
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used !
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::local_row_range_'
! Subprogram not used 
! Subprogram not used   integer :: i, ilrow, lsize
! Subprogram not used 
! Subprogram not used   ilrow = AttrVect_indexIA(sMat%data, 'lrow')
! Subprogram not used   lsize = AttrVect_lsize(sMat%data)
! Subprogram not used 
! Subprogram not used        ! Initialize start_row and end_row:
! Subprogram not used 
! Subprogram not used   start_row = sMat%data%iAttr(ilrow,1)
! Subprogram not used   end_row = sMat%data%iAttr(ilrow,1)
! Subprogram not used 
! Subprogram not used   do i=1,lsize
! Subprogram not used      start_row = min(start_row, sMat%data%iAttr(ilrow,i))
! Subprogram not used      end_row = max(end_row, sMat%data%iAttr(ilrow,i))
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used  end subroutine local_row_range_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: global_row_range_ - Global Row Extent of Non-zero Elements
!
! !DESCRIPTION:  This routine examines the input distributed 
! {\tt SparseMatrix} variable {\tt sMat}, and returns the range of 
! global row values having nonzero elements.  The first local row with 
! nonzero elements is returned in the {\tt INTEGER} argument 
! {\tt start\_row}, the last row in {\tt end\_row}. 
!
! !INTERFACE:

! Subprogram not used  subroutine global_row_range_(sMat, comm, start_row, end_row)
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_die
! Subprogram not used 
! Subprogram not used       use m_AttrVect, only : AttrVect_lsize => lsize
! Subprogram not used       use m_AttrVect, only : AttrVect_indexIA => indexIA
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used       type(SparseMatrix), intent(in)  :: sMat
! Subprogram not used       integer,            intent(in)  :: comm
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used       integer,            intent(out) :: start_row
! Subprogram not used       integer,            intent(out) :: end_row
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 15Jan01 - Jay Larson <larson@mcs.anl.gov> - API specification.
! Subprogram not used ! 25Feb01 - Jay Larson <larson@mcs.anl.gov> - Initial prototype.
! Subprogram not used ! 23Apr01 - Jay Larson <larson@mcs.anl.gov> - Modified to accomodate
! Subprogram not used !           changes to the SparseMatrix type.
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used !
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::global_row_range_'
! Subprogram not used 
! Subprogram not used   integer :: i, igrow, lsize
! Subprogram not used 
! Subprogram not used   igrow = AttrVect_indexIA(sMat%data, 'grow', dieWith=myname_)
! Subprogram not used   lsize = AttrVect_lsize(sMat%data)
! Subprogram not used 
! Subprogram not used        ! Initialize start_row and end_row:
! Subprogram not used 
! Subprogram not used   start_row = sMat%data%iAttr(igrow,1)
! Subprogram not used   end_row = sMat%data%iAttr(igrow,1)
! Subprogram not used 
! Subprogram not used   do i=1,lsize
! Subprogram not used      start_row = min(start_row, sMat%data%iAttr(igrow,i))
! Subprogram not used      end_row = max(end_row, sMat%data%iAttr(igrow,i))
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used  end subroutine global_row_range_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: local_col_range_ - Local Column Extent of Non-zero Elements
!
! !DESCRIPTION: This routine examines the input distributed 
! {\tt SparseMatrix} variable {\tt sMat}, and returns the range of
! local column values having nonzero elements.  The first local column 
! with nonzero elements is returned in the {\tt INTEGER} argument 
! {\tt start\_col}, the last column in {\tt end\_col}.
!
! !INTERFACE:

! Subprogram not used  subroutine local_col_range_(sMat, start_col, end_col)
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_die
! Subprogram not used 
! Subprogram not used       use m_AttrVect, only : AttrVect_lsize => lsize
! Subprogram not used       use m_AttrVect, only : AttrVect_indexIA => indexIA
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used       type(SparseMatrix), intent(in)  :: sMat
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used       integer,            intent(out) :: start_col
! Subprogram not used       integer,            intent(out) :: end_col
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 15Jan01 - Jay Larson <larson@mcs.anl.gov> - API specification.
! Subprogram not used ! 25Feb01 - Jay Larson <larson@mcs.anl.gov> - Initial prototype.
! Subprogram not used ! 23Apr01 - Jay Larson <larson@mcs.anl.gov> - Modified to accomodate
! Subprogram not used !           changes to the SparseMatrix type.
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used !
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::local_col_range_'
! Subprogram not used 
! Subprogram not used   integer :: i, ilcol, lsize
! Subprogram not used 
! Subprogram not used   ilcol = AttrVect_indexIA(sMat%data, 'lcol')
! Subprogram not used   lsize = AttrVect_lsize(sMat%data)
! Subprogram not used 
! Subprogram not used        ! Initialize start_col and end_col:
! Subprogram not used 
! Subprogram not used   start_col = sMat%data%iAttr(ilcol,1)
! Subprogram not used   end_col = sMat%data%iAttr(ilcol,1)
! Subprogram not used 
! Subprogram not used   do i=1,lsize
! Subprogram not used      start_col = min(start_col, sMat%data%iAttr(ilcol,i))
! Subprogram not used      end_col = max(end_col, sMat%data%iAttr(ilcol,i))
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used  end subroutine local_col_range_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: global_col_range_ - Global Column Extent of Non-zero Elements
!
! !DESCRIPTION:  This routine examines the input distributed 
! {\tt SparseMatrix} variable {\tt sMat}, and returns the range of 
! global column values having nonzero elements.  The first global 
! column with nonzero elements is returned in the {\tt INTEGER} argument 
! {\tt start\_col}, the last column in {\tt end\_col}.  
!
! !INTERFACE:

! Subprogram not used  subroutine global_col_range_(sMat, comm, start_col, end_col)
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_die
! Subprogram not used 
! Subprogram not used       use m_AttrVect, only : AttrVect_lsize => lsize
! Subprogram not used       use m_AttrVect, only : AttrVect_indexIA => indexIA
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used       type(SparseMatrix), intent(in)  :: sMat
! Subprogram not used       integer,            intent(in)  :: comm
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used       integer,            intent(out) :: start_col
! Subprogram not used       integer,            intent(out) :: end_col
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 15Jan01 - Jay Larson <larson@mcs.anl.gov> - API specification.
! Subprogram not used ! 25Feb01 - Jay Larson <larson@mcs.anl.gov> - Initial prototype.
! Subprogram not used ! 23Apr01 - Jay Larson <larson@mcs.anl.gov> - Modified to accomodate
! Subprogram not used !           changes to the SparseMatrix type.
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used !
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::global_col_range_'
! Subprogram not used 
! Subprogram not used   integer :: i, igcol, lsize
! Subprogram not used 
! Subprogram not used   igcol = AttrVect_indexIA(sMat%data, 'gcol')
! Subprogram not used   lsize = AttrVect_lsize(sMat%data)
! Subprogram not used 
! Subprogram not used        ! Initialize start_col and end_col:
! Subprogram not used 
! Subprogram not used   start_col = sMat%data%iAttr(igcol,1)
! Subprogram not used   end_col = sMat%data%iAttr(igcol,1)
! Subprogram not used 
! Subprogram not used   do i=1,lsize
! Subprogram not used      start_col = min(start_col, sMat%data%iAttr(igcol,i))
! Subprogram not used      end_col = max(end_col, sMat%data%iAttr(igcol,i))
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used  end subroutine global_col_range_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ComputeSparsitySP_ - Compute Matrix Sparsity
!
! !DESCRIPTION:  This routine computes the sparsity of a consolidated
! (all on one process) or distributed {\tt SparseMatrix}.  The input 
! {\tt SparseMatrix} argument {\tt sMat} is examined to determine the
! number of nonzero elements it holds, and this value is divided by the
! product of the number of rows and columns in {\tt sMat}.  If the 
! optional input argument {\tt comm} is given, then the distributed 
! elements are counted and the sparsity computed accordingly, and the 
! resulting value of {\tt sparsity} is returned {\em to all processes}.
!
! Given the inherent problems with multiplying and dividing large integers,
! the work in this routine is performed using floating point arithmetic on
! the logarithms of the number of rows, columns, and nonzero elements.
!
! !INTERFACE:

! Subprogram not used  subroutine ComputeSparsitySP_(sMat, sparsity, comm)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_die
! Subprogram not used       use m_mpif90
! Subprogram not used       use m_realkinds, only : SP, FP
! Subprogram not used 
! Subprogram not used       use m_AttrVect, only : AttrVect_lsize => lsize
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS: 
! Subprogram not used 
! Subprogram not used       type(SparseMatrix), intent(in)  :: sMat
! Subprogram not used       integer, optional,  intent(in)  :: comm
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used       real(SP),           intent(out) :: sparsity
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 23Apr01 - Jay Larson <larson@mcs.anl.gov> - New routine.
! Subprogram not used !
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used !
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::ComputeSparsitySP_'
! Subprogram not used 
! Subprogram not used   integer  :: num_elements, num_rows, num_cols
! Subprogram not used   real(FP) :: Lnum_elements, Lnum_rows, Lnum_cols, LMySparsity
! Subprogram not used   real(FP) :: MySparsity
! Subprogram not used   integer  :: ierr
! Subprogram not used 
! Subprogram not used        ! Extract number of nonzero elements and compute its logarithm
! Subprogram not used 
! Subprogram not used   num_elements = lsize_(sMat)
! Subprogram not used   Lnum_elements = log(REAL(num_elements,FP))
! Subprogram not used 
! Subprogram not used        ! Extract number of rows and compute its logarithm
! Subprogram not used 
! Subprogram not used   num_rows = nRows_(sMat)
! Subprogram not used   Lnum_rows = log(REAL(num_rows,FP))
! Subprogram not used 
! Subprogram not used        ! Extract number of columns and compute its logarithm
! Subprogram not used 
! Subprogram not used   num_cols = nCols_(sMat)
! Subprogram not used   Lnum_cols = log(REAL(num_cols,FP))  
! Subprogram not used 
! Subprogram not used        ! Compute logarithm of the (local) sparsity
! Subprogram not used 
! Subprogram not used   LMySparsity = Lnum_elements - Lnum_rows - Lnum_cols
! Subprogram not used 
! Subprogram not used        ! Compute the (local) sparsity from its logarithm.
! Subprogram not used 
! Subprogram not used   MySparsity = exp(LMySparsity)
! Subprogram not used 
! Subprogram not used        ! If a communicator handle is present, sum up the
! Subprogram not used        ! distributed sparsity values to all processes.  If not,
! Subprogram not used        ! return the value of MySparsity computed above.
! Subprogram not used 
! Subprogram not used   if(present(comm)) then
! Subprogram not used      call MPI_ALLREDUCE(MySparsity, sparsity, 1, MP_INTEGER, &
! Subprogram not used                         MP_SUM, comm, ierr)
! Subprogram not used      if(ierr /= 0) then
! Subprogram not used 	call MP_perr_die(myname_,"MPI_ALLREDUCE(MySparsity...",ierr)
! Subprogram not used      endif
! Subprogram not used   else
! Subprogram not used      sparsity = MySparsity
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used  end subroutine ComputeSparsitySP_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
! ----------------------------------------------------------------------
!
! !IROUTINE: ComputeSparsityDP_ - Compute Matrix Sparsity
!
! !DESCRIPTION:
! Double precision version of ComputeSparsitySP_
!
! !INTERFACE:

! Subprogram not used  subroutine ComputeSparsityDP_(sMat, sparsity, comm)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_die
! Subprogram not used       use m_mpif90
! Subprogram not used       use m_realkinds, only : DP, FP
! Subprogram not used 
! Subprogram not used       use m_AttrVect, only : AttrVect_lsize => lsize
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS: 
! Subprogram not used 
! Subprogram not used       type(SparseMatrix), intent(in)  :: sMat
! Subprogram not used       integer, optional,  intent(in)  :: comm
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used       real(DP),           intent(out) :: sparsity
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 23Apr01 - Jay Larson <larson@mcs.anl.gov> - New routine.
! Subprogram not used !
! Subprogram not used ! ______________________________________________________________________
! Subprogram not used !
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::ComputeSparsityDP_'
! Subprogram not used 
! Subprogram not used   integer  :: num_elements, num_rows, num_cols
! Subprogram not used   real(FP) :: Lnum_elements, Lnum_rows, Lnum_cols, LMySparsity
! Subprogram not used   real(FP) :: MySparsity
! Subprogram not used   integer  :: ierr
! Subprogram not used 
! Subprogram not used        ! Extract number of nonzero elements and compute its logarithm
! Subprogram not used 
! Subprogram not used   num_elements = lsize_(sMat)
! Subprogram not used   Lnum_elements = log(REAL(num_elements,FP))
! Subprogram not used 
! Subprogram not used        ! Extract number of rows and compute its logarithm
! Subprogram not used 
! Subprogram not used   num_rows = nRows_(sMat)
! Subprogram not used   Lnum_rows = log(REAL(num_rows,FP))
! Subprogram not used 
! Subprogram not used        ! Extract number of columns and compute its logarithm
! Subprogram not used 
! Subprogram not used   num_cols = nCols_(sMat)
! Subprogram not used   Lnum_cols = log(REAL(num_cols,FP))  
! Subprogram not used 
! Subprogram not used        ! Compute logarithm of the (local) sparsity
! Subprogram not used 
! Subprogram not used   LMySparsity = Lnum_elements - Lnum_rows - Lnum_cols
! Subprogram not used 
! Subprogram not used        ! Compute the (local) sparsity from its logarithm.
! Subprogram not used 
! Subprogram not used   MySparsity = exp(LMySparsity)
! Subprogram not used 
! Subprogram not used        ! If a communicator handle is present, sum up the
! Subprogram not used        ! distributed sparsity values to all processes.  If not,
! Subprogram not used        ! return the value of MySparsity computed above.
! Subprogram not used 
! Subprogram not used   if(present(comm)) then
! Subprogram not used      call MPI_ALLREDUCE(MySparsity, sparsity, 1, MP_INTEGER, &
! Subprogram not used                         MP_SUM, comm, ierr)
! Subprogram not used      if(ierr /= 0) then
! Subprogram not used 	call MP_perr_die(myname_,"MPI_ALLREDUCE(MySparsity...",ierr)
! Subprogram not used      endif
! Subprogram not used   else
! Subprogram not used      sparsity = MySparsity
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used  end subroutine ComputeSparsityDP_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: CheckBounds_ - Check for Out-of-Bounds Row/Column Values
!
! !DESCRIPTION:  This routine examines the input distributed 
! {\tt SparseMatrix} variable {\tt sMat}, and examines the global row
! and column index for each element, comparing them with the known 
! maximum values for each (as returned by the routines {\tt nRows\_()}
! and {\tt nCols\_()}, respectively).  If global row or column entries 
! are non-positive, or greater than the defined maximum values, this
! routine stops execution with an error message.  If no out-of-bounds
! values are detected, the output {\tt INTEGER} status {\tt ierror} is 
! set to zero.
!
! !INTERFACE:

! Subprogram not used  subroutine CheckBounds_(sMat, ierror)
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_die
! Subprogram not used 
! Subprogram not used       use m_AttrVect, only : AttrVect_lsize => lsize
! Subprogram not used       use m_AttrVect, only : AttrVect_indexIA => indexIA
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used       type(SparseMatrix), intent(in)  :: sMat
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used       integer,            intent(out) :: ierror
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 24Apr01 - Jay Larson <larson@mcs.anl.gov> - Initial prototype.
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used !
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::CheckBounds_'
! Subprogram not used 
! Subprogram not used   integer :: MaxRow, MaxCol, NumElements
! Subprogram not used   integer :: igrow, igcol
! Subprogram not used   integer :: i
! Subprogram not used 
! Subprogram not used        ! Initially, set ierror to zero (success):
! Subprogram not used 
! Subprogram not used   ierror = 0
! Subprogram not used 
! Subprogram not used        ! Query sMat to find the number of rows and columns:
! Subprogram not used 
! Subprogram not used   MaxRow = nRows_(sMat)
! Subprogram not used   MaxCol = nCols_(sMat)
! Subprogram not used 
! Subprogram not used        ! Query sMat for the number of nonzero elements:
! Subprogram not used 
! Subprogram not used   NumElements = lsize_(sMat)
! Subprogram not used 
! Subprogram not used        ! Query sMat to index global row and column storage indices:
! Subprogram not used 
! Subprogram not used   igrow = indexIA_(sMat=sMat,item='grow',dieWith=myname_)
! Subprogram not used   igcol = indexIA_(sMat=sMat,item='gcol',dieWith=myname_)
! Subprogram not used 
! Subprogram not used        ! Scan the entries of sMat for row or column elements that
! Subprogram not used        ! are out-of-bounds.  Here, out-of-bounds means:  1) non-
! Subprogram not used        ! positive row or column indices; 2) row or column indices
! Subprogram not used        ! exceeding the stated number of rows or columns.
! Subprogram not used 
! Subprogram not used   do i=1,NumElements
! Subprogram not used 
! Subprogram not used        ! Row index out of bounds?
! Subprogram not used 
! Subprogram not used      if((sMat%data%iAttr(igrow,i) > MaxRow) .or. &
! Subprogram not used 	  (sMat%data%iAttr(igrow,i) <= 0)) then
! Subprogram not used 	ierror = 1
! Subprogram not used 	call die(myname_,"Row index out of bounds",ierror)
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used        ! Column index out of bounds?
! Subprogram not used 
! Subprogram not used      if((sMat%data%iAttr(igcol,i) > MaxCol) .or. &
! Subprogram not used 	  (sMat%data%iAttr(igcol,i) <= 0)) then
! Subprogram not used 	ierror = 2
! Subprogram not used 	call die(myname_,"Column index out of bounds",ierror)
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used  end subroutine CheckBounds_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: row_sumSP_ - Sum Elements in Each Row
!
! !DESCRIPTION:
! Given an input {\tt SparseMatrix} argument {\tt sMat}, {\tt row\_sum\_()}
! returns the number of the rows {\tt num\_rows} in the sparse matrix and
! the sum of the elements in each row in the array {\tt sums}.  The input
! argument {\tt comm} is the Fortran 90 MPI communicator handle used to
! determine the number of rows and perform the sums.  The output arguments
! {\tt num\_rows} and {\tt sums} are valid on all processes.
!
! {\bf N.B.:  } This routine allocates an array {\tt sums}.  The user is
! responsible for deallocating this array when it is no longer needed.  
! Failure to do so will cause a memory leak.
!
! !INTERFACE:

! Subprogram not used  subroutine row_sumSP_(sMat, num_rows, sums, comm)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_die
! Subprogram not used       use m_mpif90
! Subprogram not used       use m_realkinds, only : SP, FP
! Subprogram not used 
! Subprogram not used       use m_AttrVect, only : AttrVect_lsize => lsize
! Subprogram not used       use m_AttrVect, only : AttrVect_indexIA => indexIA
! Subprogram not used       use m_AttrVect, only : AttrVect_indexRA => indexRA
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used       type(SparseMatrix), intent(in)  :: sMat
! Subprogram not used       integer,            intent(in)  :: comm
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used       integer,            intent(out) :: num_rows
! Subprogram not used       real(SP), dimension(:), pointer :: sums
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 15Jan01 - Jay Larson <larson@mcs.anl.gov> - API specification.
! Subprogram not used ! 25Jan01 - Jay Larson <larson@mcs.anl.gov> - Prototype code.
! Subprogram not used ! 23Apr01 - Jay Larson <larson@mcs.anl.gov> - Modified to accomodate
! Subprogram not used !           changes to the SparseMatrix type.
! Subprogram not used ! 18May01 - R. Jacob <jacob@mcs.anl.gov> - Use MP_TYPE function
! Subprogram not used !           to set type in the mpi_allreduce
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used !
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::row_sumSP_'
! Subprogram not used 
! Subprogram not used   integer :: i, igrow, ierr, iwgt, lsize, myID
! Subprogram not used   integer :: start_row, end_row
! Subprogram not used   integer :: mp_Type_lsums
! Subprogram not used   real(FP), dimension(:), allocatable :: lsums
! Subprogram not used   real(FP), dimension(:), allocatable :: gsums
! Subprogram not used 
! Subprogram not used        ! Determine local rank
! Subprogram not used 
! Subprogram not used   call MP_COMM_RANK(comm, myID, ierr)
! Subprogram not used 
! Subprogram not used        ! Determine on each process the row of global row indices:
! Subprogram not used 
! Subprogram not used   call global_row_range_(sMat, comm, start_row, end_row)
! Subprogram not used 
! Subprogram not used        ! Determine across the communicator the _maximum_ value of
! Subprogram not used        ! end_row, which will be assigned to num_rows on each process:
! Subprogram not used 
! Subprogram not used   call MPI_ALLREDUCE(end_row, num_rows, 1, MP_INTEGER, MP_MAX, &
! Subprogram not used                     comm, ierr)
! Subprogram not used   if(ierr /= 0) then
! Subprogram not used      call MP_perr_die(myname_,"MPI_ALLREDUCE(end_row...",ierr)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! Allocate storage for the sums on each process.
! Subprogram not used 
! Subprogram not used   allocate(lsums(num_rows), gsums(num_rows), sums(num_rows), stat=ierr)
! Subprogram not used 
! Subprogram not used   if(ierr /= 0) then
! Subprogram not used      call die(myname_,"allocate(lsums(...",ierr)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! Compute the local entries to lsum(1:num_rows) for each process:
! Subprogram not used 
! Subprogram not used   lsize = AttrVect_lsize(sMat%data)
! Subprogram not used   igrow = AttrVect_indexIA(aV=sMat%data,item='grow',dieWith=myname_)
! Subprogram not used   iwgt = AttrVect_indexRA(aV=sMat%data,item='weight',dieWith=myname_)
! Subprogram not used 
! Subprogram not used   lsums = 0._FP
! Subprogram not used   do i=1,lsize
! Subprogram not used      lsums(sMat%data%iAttr(igrow,i)) = lsums(sMat%data%iAttr(igrow,i)) + &
! Subprogram not used 	                           sMat%data%rAttr(iwgt,i)
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used        ! Compute the global sum of the entries of lsums so that all
! Subprogram not used        ! processes own the global sums.
! Subprogram not used 
! Subprogram not used   mp_Type_lsums=MP_Type(lsums)
! Subprogram not used   call MPI_ALLREDUCE(lsums, gsums, num_rows, mp_Type_lsums, MP_SUM, comm, ierr)
! Subprogram not used   if(ierr /= 0) then
! Subprogram not used      call MP_perr_die(myname_,"MPI_ALLREDUCE(lsums...",ierr)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! Copy our temporary array gsums into the output pointer sums
! Subprogram not used        ! This was done so that lsums and gsums have the same precision (FP)
! Subprogram not used        ! Precision conversion occurs here from FP to (SP or DP)
! Subprogram not used 
! Subprogram not used   sums = gsums
! Subprogram not used 
! Subprogram not used        ! Clean up...
! Subprogram not used 
! Subprogram not used   deallocate(lsums, gsums, stat=ierr)
! Subprogram not used   if(ierr /= 0) then
! Subprogram not used      call die(myname_,"deallocate(lsums...",ierr)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used  end subroutine row_sumSP_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
! ----------------------------------------------------------------------
!
! !IROUTINE: row_sumDP_ - Sum Elements in Each Row
!
! !DESCRIPTION:
! Double precision version of row_sumSP_
!
! {\bf N.B.:  } This routine allocates an array {\tt sums}.  The user is
! responsible for deallocating this array when it is no longer needed.  
! Failure to do so will cause a memory leak.
!
! !INTERFACE:

! Subprogram not used  subroutine row_sumDP_(sMat, num_rows, sums, comm)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_die
! Subprogram not used       use m_mpif90
! Subprogram not used 
! Subprogram not used       use m_realkinds, only : DP, FP
! Subprogram not used 
! Subprogram not used       use m_AttrVect, only : AttrVect_lsize => lsize
! Subprogram not used       use m_AttrVect, only : AttrVect_indexIA => indexIA
! Subprogram not used       use m_AttrVect, only : AttrVect_indexRA => indexRA
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used       type(SparseMatrix), intent(in)  :: sMat
! Subprogram not used       integer,            intent(in)  :: comm
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used       integer,            intent(out) :: num_rows
! Subprogram not used       real(DP), dimension(:), pointer :: sums
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 15Jan01 - Jay Larson <larson@mcs.anl.gov> - API specification.
! Subprogram not used ! 25Jan01 - Jay Larson <larson@mcs.anl.gov> - Prototype code.
! Subprogram not used ! 23Apr01 - Jay Larson <larson@mcs.anl.gov> - Modified to accomodate
! Subprogram not used !           changes to the SparseMatrix type.
! Subprogram not used ! 18May01 - R. Jacob <jacob@mcs.anl.gov> - Use MP_TYPE function
! Subprogram not used !           to set type in the mpi_allreduce
! Subprogram not used ! ______________________________________________________________________
! Subprogram not used !
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::row_sumDP_'
! Subprogram not used 
! Subprogram not used   integer :: i, igrow, ierr, iwgt, lsize, myID
! Subprogram not used   integer :: start_row, end_row
! Subprogram not used   integer :: mp_Type_lsums
! Subprogram not used   real(FP), dimension(:), allocatable :: lsums
! Subprogram not used   real(FP), dimension(:), allocatable :: gsums
! Subprogram not used 
! Subprogram not used        ! Determine local rank
! Subprogram not used 
! Subprogram not used   call MP_COMM_RANK(comm, myID, ierr)
! Subprogram not used 
! Subprogram not used        ! Determine on each process the row of global row indices:
! Subprogram not used 
! Subprogram not used   call global_row_range_(sMat, comm, start_row, end_row)
! Subprogram not used 
! Subprogram not used        ! Determine across the communicator the _maximum_ value of
! Subprogram not used        ! end_row, which will be assigned to num_rows on each process:
! Subprogram not used 
! Subprogram not used   call MPI_ALLREDUCE(end_row, num_rows, 1, MP_INTEGER, MP_MAX, &
! Subprogram not used                     comm, ierr)
! Subprogram not used   if(ierr /= 0) then
! Subprogram not used      call MP_perr_die(myname_,"MPI_ALLREDUCE(end_row...",ierr)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! Allocate storage for the sums on each process.
! Subprogram not used 
! Subprogram not used   allocate(lsums(num_rows), gsums(num_rows), sums(num_rows), stat=ierr)
! Subprogram not used 
! Subprogram not used   if(ierr /= 0) then
! Subprogram not used      call die(myname_,"allocate(lsums(...",ierr)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! Compute the local entries to lsum(1:num_rows) for each process:
! Subprogram not used 
! Subprogram not used   lsize = AttrVect_lsize(sMat%data)
! Subprogram not used   igrow = AttrVect_indexIA(aV=sMat%data,item='grow',dieWith=myname_)
! Subprogram not used   iwgt = AttrVect_indexRA(aV=sMat%data,item='weight',dieWith=myname_)
! Subprogram not used 
! Subprogram not used   lsums = 0._FP
! Subprogram not used   do i=1,lsize
! Subprogram not used      lsums(sMat%data%iAttr(igrow,i)) = lsums(sMat%data%iAttr(igrow,i)) + &
! Subprogram not used 	                           sMat%data%rAttr(iwgt,i)
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used        ! Compute the global sum of the entries of lsums so that all
! Subprogram not used        ! processes own the global sums.
! Subprogram not used 
! Subprogram not used   mp_Type_lsums=MP_Type(lsums)
! Subprogram not used   call MPI_ALLREDUCE(lsums, gsums, num_rows, mp_Type_lsums, MP_SUM, comm, ierr)
! Subprogram not used   if(ierr /= 0) then
! Subprogram not used      call MP_perr_die(myname_,"MPI_ALLREDUCE(lsums...",ierr)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! Copy our temporary array gsums into the output pointer sums
! Subprogram not used        ! This was done so that lsums and gsums have the same precision (FP)
! Subprogram not used        ! Precision conversion occurs here from FP to (SP or DP)
! Subprogram not used 
! Subprogram not used   sums = gsums
! Subprogram not used 
! Subprogram not used        ! Clean up...
! Subprogram not used 
! Subprogram not used   deallocate(lsums, gsums, stat=ierr)
! Subprogram not used   if(ierr /= 0) then
! Subprogram not used      call die(myname_,"deallocate(lsums...",ierr)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used  end subroutine row_sumDP_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: row_sum_checkSP_ - Check Row Sums vs. Valid Values
!
! !DESCRIPTION:  The routine {\tt row\_sum\_check()} sums the rows of 
! the input distributed (across the communicator identified by {\tt comm}) 
! {\tt SparseMatrix} variable {\tt sMat}.  It then compares these sums 
! with the {\tt num\_valid} input "valid" values stored in the array 
! {\tt valid\_sums}.  If all of the sums are within the absolute tolerence
! specified by the input argument {\tt abs\_tol} of any of the valid values,
! the output {\tt LOGICAL} flag {\tt valid} is set to {\tt .TRUE}.  
! Otherwise, this flag is returned with value {\tt .FALSE}.
!
! !INTERFACE:

! Subprogram not used  subroutine row_sum_checkSP_(sMat, comm, num_valid, valid_sums, abs_tol, valid)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_die
! Subprogram not used       use m_realkinds, only : SP, FP
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used       type(SparseMatrix), intent(in)  :: sMat
! Subprogram not used       integer,            intent(in)  :: comm
! Subprogram not used       integer,            intent(in)  :: num_valid
! Subprogram not used       real(SP),           intent(in)  :: valid_sums(num_valid)
! Subprogram not used       real(SP),           intent(in)  :: abs_tol
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used       logical,            intent(out) :: valid
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 15Jan01 - Jay Larson <larson@mcs.anl.gov> - API specification.
! Subprogram not used ! 25Feb01 - Jay Larson <larson@mcs.anl.gov> - Prototype code.
! Subprogram not used ! 06Jan03 - R. Jacob <jacob@mcs.anl.gov> - create DP and SP versions
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used !
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::row_sum_checkSP_'
! Subprogram not used 
! Subprogram not used   integer :: i, j, num_invalid, num_rows
! Subprogram not used   real(FP), dimension(:), pointer :: sums
! Subprogram not used 
! Subprogram not used        ! Compute row sums:
! Subprogram not used 
! Subprogram not used   call row_sum(sMat, num_rows, sums, comm)
! Subprogram not used 
! Subprogram not used        ! Initialize for the scanning loop (assume the matrix row
! Subprogram not used        ! sums are valid):
! Subprogram not used 
! Subprogram not used   valid = .TRUE.
! Subprogram not used   i = 1
! Subprogram not used 
! Subprogram not used   SCAN_LOOP:  do 
! Subprogram not used 
! Subprogram not used        ! Count the number of elements in valid_sums(:) that
! Subprogram not used        ! are separated from sums(i) by more than abs_tol
! Subprogram not used 
! Subprogram not used      num_invalid = 0
! Subprogram not used 
! Subprogram not used      do j=1,num_valid
! Subprogram not used 	if(abs(sums(i) - valid_sums(j)) > abs_tol) then
! Subprogram not used 	   num_invalid = num_invalid + 1
! Subprogram not used 	endif
! Subprogram not used      end do
! Subprogram not used 
! Subprogram not used        ! If num_invalid = num_valid, then we have failed to
! Subprogram not used        ! find a valid sum value within abs_tol of sums(i).  This
! Subprogram not used        ! one failure is enough to halt the process.
! Subprogram not used 
! Subprogram not used      if(num_invalid == num_valid) then
! Subprogram not used 	valid = .FALSE.
! Subprogram not used 	EXIT
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used        ! Prepare index i for the next element of sums(:)
! Subprogram not used 
! Subprogram not used      i = i + 1
! Subprogram not used      if( i > num_rows) EXIT
! Subprogram not used 
! Subprogram not used   end do SCAN_LOOP
! Subprogram not used 
! Subprogram not used  end subroutine row_sum_checkSP_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
! ----------------------------------------------------------------------
!
! !IROUTINE: row_sum_checkDP_ - Check Row Sums vs. Valid Values
!
! !DESCRIPTION:
! Double precision version of row_sum_checkSP
!
! !INTERFACE:

! Subprogram not used  subroutine row_sum_checkDP_(sMat, comm, num_valid, valid_sums, abs_tol, valid)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_die
! Subprogram not used       use m_realkinds, only : DP, FP
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used       type(SparseMatrix), intent(in)  :: sMat
! Subprogram not used       integer,            intent(in)  :: comm
! Subprogram not used       integer,            intent(in)  :: num_valid
! Subprogram not used       real(DP),           intent(in)  :: valid_sums(num_valid)
! Subprogram not used       real(DP),           intent(in)  :: abs_tol
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used       logical,            intent(out) :: valid
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 15Jan01 - Jay Larson <larson@mcs.anl.gov> - API specification.
! Subprogram not used ! 25Feb01 - Jay Larson <larson@mcs.anl.gov> - Prototype code.
! Subprogram not used ! 06Jan03 - R. Jacob <jacob@mcs.anl.gov> - create DP and SP versions
! Subprogram not used ! ______________________________________________________________________
! Subprogram not used !
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::row_sum_checkDP_'
! Subprogram not used 
! Subprogram not used   integer :: i, j, num_invalid, num_rows
! Subprogram not used   real(FP), dimension(:), pointer :: sums
! Subprogram not used 
! Subprogram not used        ! Compute row sums:
! Subprogram not used 
! Subprogram not used   call row_sum(sMat, num_rows, sums, comm)
! Subprogram not used 
! Subprogram not used        ! Initialize for the scanning loop (assume the matrix row
! Subprogram not used        ! sums are valid):
! Subprogram not used 
! Subprogram not used   valid = .TRUE.
! Subprogram not used   i = 1
! Subprogram not used 
! Subprogram not used   SCAN_LOOP:  do 
! Subprogram not used 
! Subprogram not used        ! Count the number of elements in valid_sums(:) that
! Subprogram not used        ! are separated from sums(i) by more than abs_tol
! Subprogram not used 
! Subprogram not used      num_invalid = 0
! Subprogram not used 
! Subprogram not used      do j=1,num_valid
! Subprogram not used 	if(abs(sums(i) - valid_sums(j)) > abs_tol) then
! Subprogram not used 	   num_invalid = num_invalid + 1
! Subprogram not used 	endif
! Subprogram not used      end do
! Subprogram not used 
! Subprogram not used        ! If num_invalid = num_valid, then we have failed to
! Subprogram not used        ! find a valid sum value within abs_tol of sums(i).  This
! Subprogram not used        ! one failure is enough to halt the process.
! Subprogram not used 
! Subprogram not used      if(num_invalid == num_valid) then
! Subprogram not used 	valid = .FALSE.
! Subprogram not used 	EXIT
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used        ! Prepare index i for the next element of sums(:)
! Subprogram not used 
! Subprogram not used      i = i + 1
! Subprogram not used      if( i > num_rows) EXIT
! Subprogram not used 
! Subprogram not used   end do SCAN_LOOP
! Subprogram not used 
! Subprogram not used  end subroutine row_sum_checkDP_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: Sort_ - Generate Index Permutation
!
! !DESCRIPTION:
! The subroutine {\tt Sort\_()} uses a list of sorting keys defined by 
! the input {\tt List} argument {\tt key\_list}, searches for the appropriate 
! integer or real attributes referenced by the items in {\tt key\_list} 
! ( that is, it identifies the appropriate entries in {sMat\%data\%iList} 
! and {\tt sMat\%data\%rList}), and then uses these keys to generate an index 
! permutation {\tt perm} that will put the nonzero matrix entries of stored
! in {\tt sMat\%data} in lexicographic order as defined by {\tt key\_ist} 
! (the ordering in {\tt key\_list} being from left to right.  The optional 
! {\tt LOGICAL} array input argument {\tt descend} specifies whether or
! not to sort by each key in {\em descending} order or {\em ascending} 
! order.  Entries in {\tt descend} that have value {\tt .TRUE.} correspond 
! to a sort by the corresponding key in descending order.  If the argument
! {\tt descend} is not present, the sort is performed for all keys in
! ascending order.
!
! !INTERFACE:

 subroutine Sort_(sMat, key_list, perm, descend)

!
! !USES:
!
      use m_die ,          only : die
      use m_stdio ,        only : stderr

      use m_List ,         only : List

      use m_AttrVect, only: AttrVect_Sort => Sort

      implicit none
!
! !INPUT PARAMETERS: 

      type(SparseMatrix),              intent(in) :: sMat
      type(List),                      intent(in) :: key_list
      logical, dimension(:), optional, intent(in) :: descend
!
! !OUTPUT PARAMETERS: 

      integer, dimension(:), pointer              :: perm


! !REVISION HISTORY:
! 24Apr01 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::Sort_'

  if(present(descend)) then
     call AttrVect_Sort(sMat%data, key_list, perm, descend)
  else
     call AttrVect_Sort(sMat%data, key_list, perm)
  endif

 end Subroutine Sort_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: Permute_ - Permute Matrix Elements using Supplied Index Permutation
!
! !DESCRIPTION:
! The subroutine {\tt Permute\_()} uses an input index permutation 
! {\tt perm} to re-order the entries of the {\tt SparseMatrix} argument 
! {\tt sMat}.  The index permutation {\tt perm} is generated using the 
! routine {\tt Sort\_()} (in this module).
!
! !INTERFACE:

 subroutine Permute_(sMat, perm)

!
! !USES:
!
      use m_die ,          only : die
      use m_stdio ,        only : stderr

      use m_AttrVect, only: AttrVect_Permute => Permute

      implicit none
!
! !INPUT PARAMETERS: 


      integer, dimension(:), pointer               :: perm
!
! !INPUT/OUTPUT PARAMETERS: 

      type(SparseMatrix),            intent(inout) :: sMat


! !REVISION HISTORY:
! 24Apr01 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::Permute_'

  call AttrVect_Permute(sMat%data, perm)

 end Subroutine Permute_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: SortPermute_ - Sort and Permute Matrix Elements
!
! !DESCRIPTION:
! The subroutine {\tt SortPermute\_()} uses a list of sorting keys defined 
! by the input {\tt List} argument {\tt key\_list}, searches for the 
! appropriate integer or real attributes referenced by the items in 
! {\tt key\_ist} ( that is, it identifies the appropriate entries in 
! {sMat\%data\%iList} and {\tt sMat\%data\%rList}), and then uses these 
! keys to generate an index permutation that will put the nonzero matrix 
! entries of stored in {\tt sMat\%data} in lexicographic order as defined 
! by {\tt key\_list} (the ordering in {\tt key\_list} being from left to 
! right.  The optional {\tt LOGICAL} array input argument {\tt descend} 
! specifies whether or not to sort by each key in {\em descending} order 
! or {\em ascending} order.  Entries in {\tt descend} that have value 
! {\tt .TRUE.} correspond to a sort by the corresponding key in descending 
! order.  If the argument {\tt descend} is not present, the sort is 
! performed for all keys in ascending order.
!
! Once this index permutation is created, it is applied to re-order the 
! entries of the {\tt SparseMatrix} argument {\tt sMat} accordingly.
!
! !INTERFACE:

 subroutine SortPermute_(sMat, key_list, descend)

!
! !USES:
!
      use m_die ,          only : die
      use m_stdio ,        only : stderr

      use m_List ,         only : List

      implicit none
!
! !INPUT PARAMETERS: 

      type(List),                      intent(in)    :: key_list
      logical, dimension(:), optional, intent(in)    :: descend
!
! !INPUT/OUTPUT PARAMETERS: 

      type(SparseMatrix),              intent(inout) :: sMat

! !REVISION HISTORY:
! 24Apr01 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::SortPermute_'

  integer :: ier
  integer, dimension(:), pointer :: perm

       ! Create index permutation perm(:)

  if(present(descend)) then
     call Sort_(sMat, key_list, perm, descend)
  else
     call Sort_(sMat, key_list, perm)
  endif

       ! Apply index permutation perm(:) to re-order sMat:

  call Permute_(sMat, perm)

       ! Clean up

  deallocate(perm, stat=ier)
  if(ier/=0) call die(myname_, "deallocate(perm)", ier)

 end subroutine SortPermute_

 end module m_SparseMatrix



