!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!-----------------------------------------------------------------------
! CVS $Id$
! CVS $Name$ 
!BOP -------------------------------------------------------------------
!
! !MODULE: m_GeneralGrid -- Physical Coordinate Grid Information Storage
!
! !DESCRIPTION:
! The {\tt GeneralGrid} data type is a flexible, generic structure for 
! storing physical coordinate grid information.  The {\tt GeneralGrid} 
! may be employed to store coordinate grids of arbitrary dimension, and 
! is also capable of supporting unstructured grids such as meteorological 
! observation data streams.  The grid is representated by a literal 
! listing of the gridpoint coordinates, along with other integer and real 
! {\em attributes} associated with each location.  Examples of real 
! non-coordinate attributes are grid cell length, cross-sectional area, and 
! volume elements, projections of local directional unit vectors onto 
! {\em et cetera}  A {\tt GeneralGrid} as at minimum one integer 
! attribute---{\em the global grid point number}, or {\tt GlobGridNum}, 
! which serves as a unique identifier for each physical grid location.
!
! The real attributes of of the {\tt GeneralGrid} are grouped as {\tt List} 
! components:
! \begin{itemize}
! \item {\tt GGrid\%coordinate\_list} contains the list of the physical 
! dimension names of the grid.  The user initializes a {\tt List} by 
! supplying the items in it as a string with the items delimitted by 
! colons.  For example, setting the coordinates for Euclidean 3-space 
! is accomplished by a choice of {\tt 'x:y:z'}, cylindrical coordinates 
! by {\tt 'rho:theta:z'}, spherical coordinates by {\tt 'r:theta:phi'}, 
! {\em et cetera}.
! \item {\tt GGrid\%weight\_list} contains the names of the spatial 
! cell length, area, and volume weights associated with the grid.  These 
! are also stored in {\tt List} form, and are set by the user in the same
! fashion as described above for coordinates.  For example, one might 
! wish create cell weight attributes for a cylindrical grid by defining 
! a weight list of {\tt 'drho:dphi:rhodphi:dz}.
! \item {\tt GGrid\%other\_list} is space for the user to define other 
! real attributes.  For example, one might wish to do vector calculus
! operatons in spherical coordinates.  Since the spherical coordinate 
! unit vectors ${\hat r}$, ${\hat \theta}$, and ${\hat \phi}$
! vary in space, it is sometimes useful to store their projections on 
! the fixed Euclidean unit vectors ${\bf \hat x}$, ${\bf \hat y}$, and 
! ${\bf \hat z}$.  To do this one might set up a list of attributes 
! using the string
! \begin{verbatim}
! 'rx:ry:rz:thetax:thetay:thetaz:phix:phiy:phyz'
! \end{verbatim}
! \item {\tt GGrid\%index\_list} provides space for the user to define 
! integer attributes such as alternative indexing schemes, indices for 
! defining spatial regions, {\em et cetera}.  This attribute list contains
! all the integer attributes for the {\tt GeneralGrid} save one:  the 
! with the ever-present {\em global gridpoint number attribute} 
! {\tt GlobGridNum}, which is set automatically by MCT.
! \end{itemize}
!
! This module contains the definition of the {\tt GeneralGrid} datatype, 
! various methods for creating and destroying it, query methods, and tools 
! for multiple-key sorting of gridpoints.
!
! !INTERFACE:

 module m_GeneralGrid

!
! !USES:
!
      use m_List, only : List   ! Support for List components.

      use m_AttrVect, only : AttrVect ! Support for AttrVect component.

      implicit none

      private   ! except

! !PUBLIC TYPES:

      public :: GeneralGrid      ! The class data structure

    Type GeneralGrid



      type(List)                     :: coordinate_list
      type(List)                     :: coordinate_sort_order
      logical, dimension(:), pointer :: descend
      type(List)                     :: weight_list
      type(List)                     :: other_list
      type(List)                     :: index_list
      type(AttrVect)                 :: data
    End Type GeneralGrid

! !PUBLIC MEMBER FUNCTIONS:

      public :: init             ! Create a GeneralGrid
      public :: initCartesian    !
      public :: initUnstructured !
      public :: clean            ! Destroy a GeneralGrid
      public :: zero             ! Zero data in a GeneralGrid

                             ! Query functions-----------------
      public :: dims         ! Return dimensionality of the GeneralGrid
      public :: indexIA      ! Index integer attribute (indices)
      public :: indexRA      ! Index integer attribute (coords/weights)
      public :: lsize        ! Return local number of points
      public :: exportIAttr  ! Return INTEGER attribute as a vector
      public :: exportRAttr  ! Return REAL attribute as a vector

                             ! Manipulation--------------------
      public :: importIAttr  ! Insert INTEGER vector as attribute
      public :: importRAttr  ! Insert REAL vector as attribute
      public :: Sort         ! Sort point data by coordinates -> permutation
      public :: Permute      ! Rearrange point data using input permutation
      public :: SortPermute  ! Sort and Permute point data

    interface init  ; module procedure &
        init_, &
        initl_, &
        initgg_
    end interface
    interface initCartesian ; module procedure &
        initCartesianSP_, &
	initCartesianDP_
    end interface
    interface initUnstructured ; module procedure &
        initUnstructuredSP_, &
	initUnstructuredDP_
    end interface
    interface clean ; module procedure clean_ ; end interface
    interface zero ; module procedure zero_ ; end interface

    interface dims ; module procedure dims_ ; end interface
    interface indexIA ; module procedure indexIA_ ; end interface
    interface indexRA ; module procedure indexRA_ ; end interface
    interface lsize   ; module procedure lsize_   ; end interface

    interface exportIAttr ; module procedure exportIAttr_ ; end interface
    interface exportRAttr ; module procedure &
       exportRAttrSP_, &
       exportRAttrDP_
    end interface
    interface importIAttr ; module procedure importIAttr_ ; end interface
    interface importRAttr ; module procedure &
       importRAttrSP_, &
       importRAttrDP_
    end interface

    interface Sort    ; module procedure Sort_    ; end interface
    interface Permute ; module procedure Permute_ ; end interface
    interface SortPermute ; module procedure SortPermute_ ; end interface

! !PUBLIC DATA MEMBERS:

! CHARACTER Tag for GeneralGrid Global Grid Point Identification Number

  character(len=*), parameter :: GlobGridNum='GlobGridNum'

! !SEE ALSO:
! The MCT module m_AttrVect and the mpeu module m_List.

! !REVISION HISTORY:
! 25Sep00 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
! 31Oct00 - J.W. Larson <larson@mcs.anl.gov> - modified the
!           GeneralGrid type to allow inclusion of grid cell
!           dimensions (lengths) and area/volume weights.
! 15Jan01 - J.W. Larson implemented new GeneralGrid type 
!           definition and added numerous APIs.
! 17Jan01 - J.W. Larson fixed minor bug in module header use
!           statement.
! 19Jan01 - J.W. Larson added other_list and coordinate_sort_order
!           components to the GeneralGrid type.
! 21Mar01 - J.W. Larson - deleted the initv_ API (more study
!           needed before implementation.
!  2May01 - J.W. Larson - added initgg_ API (replaces old initv_).
! 13Dec01 - J.W. Larson - added import and export methods.
! 27Mar02 - J.W. Larson <larson@mcs.anl.gov> - Corrected usage of
!           m_die routines throughout this module.
!  5Aug02 - E. Ong <eong@mcs.anl.gov> - Modified GeneralGrid usage 
!           to allow user-defined grid numbering schemes.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='MCT::m_GeneralGrid'

 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - Create an Empty GeneralGrid
!
! !DESCRIPTION:
! The routine {\tt init\_()} creates the storage space for grid point
! coordinates, area/volume weights, and other coordinate data ({\em e.g.}, 
! local cell dimensions).  These data are referenced by {\tt List}
! components that are also created by this routine (see the documentation 
! of the declaration section of this module for more details about setting 
! list information).  Each of the input {\tt CHARACTER} arguments is a 
! colon-delimited string of attribute names, each corrsponding to a 
! {\tt List} element of the output {\tt GeneralGrid} argument {\tt GGrid},
! and are summarized in the table below:
!
!\begin{table}[htbp]
!\begin{center}
!\begin{tabular}{|l|l|l|l|}
!\hline
!{\bf Argument} & {\bf Component of {\tt GGrid}} & {\bf Significance} & {\bf Required?} \!\hline

!{\tt CoordChars} & {\tt GGrid\%coordinate\_list} & Dimension Names & Yes \!\hline

!{\tt CoordSortOrder} & {\tt GGrid\%coordinate\_sort\_order} & Grid Point & No \! & & Sorting Keys & \!\hline


!{\tt WeightChars} & {\tt GGrid\%weight\_list} & Grid Cell & No \! & & Length, Area, and & \! & & Volume Weights & \!\hline



!{\tt OtherChars} & {\tt GGrid\%other\_list} & All Other & No \! & & Real Attributes & \!\hline


!{\tt IndexChars} & {\tt GGrid\%index\_list} & All Other & No \! & & Integer Attributes & \!\hline


!\end{tabular}
!\end{center}
!\end{table}
!
! The input {\tt INTEGER} argument {\tt lsize} defines the number of grid points 
! to be stored in {\tt GGrid}.
!
! If a set of sorting keys is supplied in the argument {\tt CoordSortOrder}, 
! the user can control whether the sorting by each key is in descending or 
! ascending order by supplying the input {\tt LOGICAL} array {\tt descend(:)}.  
! By default, all sorting is in {\em ascending} order for each key if the 
! argument {\tt descend} is not provided.
!
! {\bf N.B.}:  The output {\tt GeneralGrid} {\tt GGrid} is dynamically
! allocated memory.  When one no longer needs {\tt GGrid}, one should
! release this space by invoking {\tt clean()} for the {\tt GeneralGrid}.
!
! !INTERFACE:

 subroutine init_(GGrid, CoordChars, CoordSortOrder, descend, WeightChars, &
                  OtherChars, IndexChars, lsize )
!
! !USES:
!
      use m_stdio
      use m_die

      use m_List,     only : List
      use m_List,     only : List_init => init
      use m_List,     only : List_nitem => nitem
      use m_List,     only : List_shared => GetSharedListIndices
      use m_List,     only : List_append => append
      use m_List,     only : List_copy => copy
      use m_List,     only : List_nullify => nullify
      use m_List,     only : List_clean => clean

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_init => init

      implicit none

! !INPUT PARAMETERS:
!
      character(len=*),                intent(in) :: CoordChars
      character(len=*),      optional, intent(in) :: CoordSortOrder
      character(len=*),      optional, intent(in) :: WeightChars
      logical, dimension(:), optional, pointer    :: descend
      character(len=*),      optional, intent(in) :: OtherChars
      character(len=*),      optional, intent(in) :: IndexChars
      integer,               optional, intent(in) :: lsize

! !OUTPUT PARAMETERS:
!
      type(GeneralGrid), intent(out)   :: GGrid

! !REVISION HISTORY:
! 25Sep00 - Jay Larson <larson@mcs.anl.gov> - initial prototype
! 15Jan01 - Jay Larson <larson@mcs.anl.gov> - modified to fit
!           new GeneralGrid definition.  
! 19Mar01 - Jay Larson <larson@mcs.anl.gov> - added OtherChars
! 25Apr01 - Jay Larson <larson@mcs.anl.gov> - added GlobGridNum
!           as a mandatory integer attribute.
! 13Jun01 - Jay Larson <larson@mcs.anl.gov> - No longer define 
!           blank List attributes of the GeneralGrid.  Previous
!           versions of this routine had this feature, and this
!           caused problems with the GeneralGrid Send and Receive
!           operations on the AIX platform.
! 13Jun01 - R. Jacob <jacob@mcs.anl.gov> - nullify any pointers
!           for lists not declared.
! 15Feb02 - Jay Larson <larson@mcs.anl.gov> - made the input 
!           argument CoordSortOrder mandatory (rather than
!           optional).
! 18Jul02 - E. Ong <eong@mcs.anl.gov> - replaced this version of 
!           init with one that calls initl_. 
!  5Aug02 - E. Ong <eong@mcs.anl.gov> - made the input argument
!           CoordSortOrder optional to allow user-defined grid
!           numbering schemes.
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::init_'

  ! List to store real and integer attributes
  type(List) :: RAList, IAList

  ! Overlapping index storage arrays:
  integer, dimension(:), pointer :: &
       CoordListIndices, CoordSortOrderIndices

  ! Temporary vars
  integer :: NumShared, nitems, i, l, ierr

       ! Let's begin by nullifying everything:

  call List_nullify(GGrid%coordinate_list)
  call List_nullify(GGrid%coordinate_sort_order)
  call List_nullify(GGrid%weight_list)
  call List_nullify(GGrid%other_list)
  call List_nullify(GGrid%index_list)
  nullify(GGrid%descend)

        ! Convert the Character arguments to the appropriate
        ! GeneralGrid components. 

        ! Set up the integer and real attribute lists.

  call List_init(GGrid%coordinate_list,trim(CoordChars))
  call List_copy(RAList,GGrid%coordinate_list)

  if(present(CoordSortOrder)) then
     call List_init(GGrid%coordinate_sort_order,trim(CoordSortOrder))
  endif

  if(present(WeightChars)) then
     call List_init(GGrid%weight_list,trim(WeightChars))
     call List_append(RAList, GGrid%weight_list)
  endif

  if(present(OtherChars)) then
     call List_init(GGrid%other_list,trim(OtherChars))
     call List_append(RAList, GGrid%other_list)
  endif

  call List_init(IAList,GlobGridNum)

  if(present(IndexChars)) then
     call List_init(GGrid%index_list,trim(IndexChars))
     call List_append(IAList, GGrid%index_list)
  endif

        ! Check the lists that we've initialized :

  nitems = List_nitem(GGrid%coordinate_list)

        ! Check the number of coordinates

  if(nitems <= 0) then
     write(stderr,*) myname_, &
	  ':: ERROR CoordList is empty!'
     call die(myname_,'List_nitem(CoordList) <= 0',nitems)
  endif

        ! Check the items in the coordinate list and the 
        ! coordinate grid sort keys...they should contain
        ! the same items.

  if(present(CoordSortOrder)) then

     call List_shared(GGrid%coordinate_list,GGrid%coordinate_sort_order, &
	              NumShared,CoordListIndices,CoordSortOrderIndices)

     deallocate(CoordListIndices,CoordSortOrderIndices,stat=ierr)
     if(ierr/=0) call die(myname_,'deallocate(CoordListIndices..)',ierr)
     
     if(NumShared /= nitems) then
	call die(myname_,'CoordSortOrder must have the same items &
	         & as CoordList',abs(nitems-NumShared))
     endif

  endif

        ! If the LOGICAL argument descend is present, check the
        ! number of entries to ensure they match the grid dimensionality.
        ! If descend is not present, assume all coordinate grid point
        ! sortings will be in ascending order.

  if(present(descend)) then

     if( ( (.not.associated(descend)) .or. &
	   (.not.present(CoordSortOrder)) ) .or. &
	   (size(descend) /= nitems) ) then
	
	write(stderr,*) myname_, &
	     ':: ERROR using descend argument, &
	     &associated(descend) = ', associated(descend), &
	     ' present(CoordSortOrder) = ', present(CoordSortOrder), &
	     ' size(descend) = ', size(descend), &
	     ' List_nitem(CoordSortOrder) = ', &
             List_nitem(GGrid%coordinate_sort_order)
	call die(myname_, 'ERROR using -descend- argument; &
             & see stderr file for details')
     endif

  endif

       ! Finally, Initialize GGrid%descend from descend(:).
       ! If descend argument is not present, set it to the default .false.

  if(present(CoordSortOrder)) then

     allocate(GGrid%descend(nitems), stat=ierr)
     if(ierr /= 0) call die(myname_,"allocate GGrid%descend...",ierr)

     if(present(descend)) then

        do i=1,nitems
           GGrid%descend(i) = descend(i)
        enddo

     else

        do i=1,nitems
           GGrid%descend(i) = .FALSE.
        enddo

     endif
        
  endif
  
       ! Initialize GGrid%data using IAList, RAList, and lsize (if
       ! present).

  l = 0
  if(present(lsize)) l=lsize

  call AttrVect_init(GGrid%data, IAList, RAList, l)


       ! Deallocate the temporary variables

  call List_clean(IAList)
  call List_clean(RAList)
  
 end subroutine init_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initl_ - Create an Empty GeneralGrid from Lists
!
! !DESCRIPTION:
! The routine {\tt initl\_()} creates the storage space for grid point
! coordinates, area/volume weights, and other coordinate data ({\em e.g.}, 
! local cell dimensions).  These data are referenced by {\tt List}
! components that are also created by this routine (see the documentation 
! of the declaration section of this module for more details about setting 
! list information).  Each of the input {\tt List} arguments is used 
! directly to create the corresponding 
! {\tt List} element of the output {\tt GeneralGrid} argument {\tt GGrid},
! and are summarized in the table below:
!
!\begin{table}[htbp]
!\begin{center}
!\begin{tabular}{|l|l|l|l|}
!\hline
!{\bf Argument} & {\bf Component of {\tt GGrid}} & {\bf Significance} & {\bf Required?} \!\hline

!{\tt CoordList} & {\tt GGrid\%coordinate\_list} & Dimension Names & Yes \!\hline

!{\tt CoordSortOrder} & {\tt GGrid\%coordinate\_sort\_order} & Grid Point & No \! & & Sorting Keys & \!\hline


!{\tt WeightList} & {\tt GGrid\%weight\_list} & Grid Cell & No \! & & Length, Area, and & \! & & Volume Weights & \!\hline



!{\tt OtherList} & {\tt GGrid\%other\_list} & All Other & No \! & & Real Attributes & \!\hline


!{\tt IndexList} & {\tt GGrid\%index\_list} & All Other & No \! & & Integer Attributes & \!\hline


!\end{tabular}
!\end{center}
!\end{table}
!
! The input {\tt INTEGER} argument {\tt lsize} defines the number of grid points 
! to be stored in {\tt GGrid}.
!
! If a set of sorting keys is supplied in the argument {\tt CoordSortOrder}, 
! the user can control whether the sorting by each key is in descending or 
! ascending order by supplying the input {\tt LOGICAL} array {\tt descend(:)}.  
! By default, all sorting is in {\em ascending} order for each key if the 
! argument {\tt descend} is not provided.
!
! {\bf N.B.}:  The output {\tt GeneralGrid} {\tt GGrid} is dynamically
! allocated memory.  When one no longer needs {\tt GGrid}, one should
! release this space by invoking {\tt clean()} for the {\tt GeneralGrid}.
!
! !INTERFACE:

! Subprogram not used  subroutine initl_(GGrid, CoordList, CoordSortOrder, descend, WeightList, &
! Subprogram not used                    OtherList, IndexList, lsize )
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used 
! Subprogram not used       use m_stdio
! Subprogram not used       use m_die
! Subprogram not used 
! Subprogram not used       use m_List,     only : List
! Subprogram not used       use m_List,     only : List_init => init
! Subprogram not used       use m_List,     only : List_allocated => allocated
! Subprogram not used       use m_List,     only : List_nitem => nitem
! Subprogram not used       use m_List,     only : List_shared => GetSharedListIndices
! Subprogram not used       use m_List,     only : List_append => append
! Subprogram not used       use m_List,     only : List_copy => copy
! Subprogram not used       use m_List,     only : List_nullify => nullify
! Subprogram not used       use m_List,     only : List_clean => clean
! Subprogram not used 
! Subprogram not used       use m_AttrVect, only : AttrVect
! Subprogram not used       use m_AttrVect, only : AttrVect_init => init
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       Type(List),                      intent(in)  :: CoordList
! Subprogram not used       Type(List),            optional, intent(in)  :: CoordSortOrder
! Subprogram not used       Type(List),            optional, intent(in)  :: WeightList
! Subprogram not used       logical, dimension(:), optional, pointer     :: descend
! Subprogram not used       Type(List),            optional, intent(in)  :: OtherList
! Subprogram not used       Type(List),            optional, intent(in)  :: IndexList
! Subprogram not used       integer,               optional, intent(in)  :: lsize
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       type(GeneralGrid),               intent(out) :: GGrid
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 10May01 - Jay Larson <larson@mcs.anl.gov> - initial version
! Subprogram not used !  8Aug01 - E.T. Ong <eong@mcs.anl.gov> - changed list assignment(=)
! Subprogram not used !           to list copy to avoid compiler bugs with pgf90
! Subprogram not used ! 17Jul02 - E. Ong <eong@mcs.anl.gov> - general revision; 
! Subprogram not used !           added error checks
! Subprogram not used !  5Aug02 - E. Ong <eong@mcs.anl.gov> - made input argument
! Subprogram not used !           CoordSortOrder optional to allow for user-defined
! Subprogram not used !           grid numbering schemes
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used !
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::initl_'
! Subprogram not used 
! Subprogram not used   ! List to store real and integer attributes
! Subprogram not used   type(List) :: RAList, IAList
! Subprogram not used 
! Subprogram not used   ! Overlapping attribute index storage arrays:
! Subprogram not used   integer, dimension(:), pointer :: &
! Subprogram not used        CoordListIndices, CoordSortOrderIndices
! Subprogram not used 
! Subprogram not used   ! Temporary vars
! Subprogram not used   integer :: NumShared, nitems, i, l, ierr
! Subprogram not used 
! Subprogram not used        ! Let's begin by nullifying everything:
! Subprogram not used 
! Subprogram not used   call List_nullify(GGrid%coordinate_list)
! Subprogram not used   call List_nullify(GGrid%coordinate_sort_order)
! Subprogram not used   call List_nullify(GGrid%weight_list)
! Subprogram not used   call List_nullify(GGrid%other_list)
! Subprogram not used   call List_nullify(GGrid%index_list)
! Subprogram not used   nullify(GGrid%descend)
! Subprogram not used 
! Subprogram not used         ! Check the arguments:
! Subprogram not used 
! Subprogram not used   nitems = List_nitem(CoordList)
! Subprogram not used 
! Subprogram not used         ! Check the number of coordinates
! Subprogram not used 
! Subprogram not used   if(nitems <= 0) then
! Subprogram not used      write(stderr,*) myname_, &
! Subprogram not used           ':: ERROR CoordList is empty!'
! Subprogram not used      call die(myname_,'List_nitem(CoordList) <= 0',nitems)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used         ! Check the items in the coordinate list and the 
! Subprogram not used         ! coordinate grid sort keys...they should contain
! Subprogram not used         ! the same items.
! Subprogram not used 
! Subprogram not used   if(present(CoordSortOrder)) then
! Subprogram not used 
! Subprogram not used      call List_shared(CoordList,CoordSortOrder,NumShared, &
! Subprogram not used                       CoordListIndices,CoordSortOrderIndices)
! Subprogram not used 
! Subprogram not used      deallocate(CoordListIndices,CoordSortOrderIndices,stat=ierr)
! Subprogram not used      if(ierr/=0) call die(myname_,'deallocate(CoordListIndices..)',ierr)
! Subprogram not used      
! Subprogram not used      if(NumShared /= nitems) then
! Subprogram not used         call die(myname_,'CoordSortOrder must have the same items &
! Subprogram not used                  & as CoordList',abs(nitems-NumShared))
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used         ! If the LOGICAL argument descend is present, check the
! Subprogram not used         ! number of entries to ensure they match the grid dimensionality.
! Subprogram not used         ! If descend is not present, assume all coordinate grid point
! Subprogram not used         ! sortings will be in ascending order.
! Subprogram not used 
! Subprogram not used   if(present(descend)) then
! Subprogram not used 
! Subprogram not used      if( ( (.not.associated(descend)) .or. &
! Subprogram not used            (.not.present(CoordSortOrder)) ) .or. &
! Subprogram not used            (size(descend) /= nitems) ) then
! Subprogram not used         
! Subprogram not used         write(stderr,*) myname_, &
! Subprogram not used              ':: ERROR using descend argument, &
! Subprogram not used              &associated(descend) = ', associated(descend), &
! Subprogram not used              ' present(CoordSortOrder) = ', present(CoordSortOrder), &
! Subprogram not used              ' size(descend) = ', size(descend), &
! Subprogram not used              ' List_nitem(CoordSortOrder) = ', &
! Subprogram not used              List_nitem(CoordSortOrder)
! Subprogram not used         call die(myname_, 'ERROR using -descend- argument; &
! Subprogram not used              &stderr file for details')
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! Initialize GGrid%descend from descend(:), if present.  If
! Subprogram not used        ! the argument descend(:) was not passed, set GGrid%descend
! Subprogram not used        ! to the default .false.
! Subprogram not used 
! Subprogram not used   if(present(CoordSortOrder)) then
! Subprogram not used 
! Subprogram not used      allocate(GGrid%descend(nitems), stat=ierr)
! Subprogram not used      if(ierr /= 0) call die(myname_,"allocate GGrid%descend...",ierr)
! Subprogram not used 
! Subprogram not used      if(present(descend)) then
! Subprogram not used 
! Subprogram not used         do i=1,nitems
! Subprogram not used            GGrid%descend(i) = descend(i)
! Subprogram not used         enddo
! Subprogram not used 
! Subprogram not used      else
! Subprogram not used 
! Subprogram not used         do i=1,nitems
! Subprogram not used            GGrid%descend(i) = .FALSE.
! Subprogram not used         enddo
! Subprogram not used 
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used   endif
! Subprogram not used   
! Subprogram not used        ! Process input lists and create the appropriate GeneralGrid
! Subprogram not used        ! List components
! Subprogram not used 
! Subprogram not used   call List_copy(GGrid%coordinate_list,CoordList)
! Subprogram not used   call List_copy(RAList,CoordList)
! Subprogram not used 
! Subprogram not used   if(present(CoordSortOrder)) then
! Subprogram not used      if(List_allocated(CoordSortOrder)) then
! Subprogram not used         call List_copy(GGrid%coordinate_sort_order,CoordSortOrder)
! Subprogram not used      else
! Subprogram not used         call die(myname_,"Argument CoortSortOrder not allocated")
! Subprogram not used      endif
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! Concatenate present input Lists to create RAList, and
! Subprogram not used        ! at the same time assign the List components of GGrid
! Subprogram not used 
! Subprogram not used   if(present(WeightList)) then
! Subprogram not used      if(List_allocated(WeightList)) then
! Subprogram not used         call List_copy(GGrid%weight_list,WeightList)
! Subprogram not used         call List_append(RAList, WeightList)
! Subprogram not used      else
! Subprogram not used         call die(myname_,"Argument WeightList not allocated")
! Subprogram not used      endif
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   if(present(OtherList)) then
! Subprogram not used      if(List_allocated(OtherList)) then
! Subprogram not used         call List_copy(GGrid%other_list,OtherList)
! Subprogram not used         call List_append(RAList, OtherList)
! Subprogram not used      else
! Subprogram not used         call die(myname_,"Argument OtherList not allocated")
! Subprogram not used      endif
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! Concatenate present input Lists to create IAList
! Subprogram not used 
! Subprogram not used   call List_init(IAList,GlobGridNum)
! Subprogram not used 
! Subprogram not used   if(present(IndexList)) then
! Subprogram not used      call List_copy(GGrid%index_list,IndexList)
! Subprogram not used      call List_append(IAList, IndexList)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! Initialize GGrid%data using IAList, RAList, and lsize (if
! Subprogram not used        ! present).
! Subprogram not used 
! Subprogram not used   l = 0
! Subprogram not used   if(present(lsize)) l = lsize
! Subprogram not used 
! Subprogram not used   call AttrVect_init(GGrid%data, IAList, RAList, l)
! Subprogram not used 
! Subprogram not used        ! Deallocate the temporary variables
! Subprogram not used 
! Subprogram not used   call List_clean(IAList)
! Subprogram not used   call List_clean(RAList)
! Subprogram not used 
! Subprogram not used  end subroutine initl_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initgg_ - Create a GeneralGrid from Another
!
! !DESCRIPTION:
! The routine {\tt initgg\_()} creates the storage space for grid point
! coordinates, area/volume weights, and other coordinate data ({\em e.g.}, 
! nearest-neighbor coordinates).  These data are all copied from the 
! already initialized input {\tt GeneralGrid} argument {\tt iGGrid}.  This 
! routine initializes the output {\tt GeneralGrid} argument {\tt oGGrid} 
! with the same {\tt List} data as {\tt iGGrid}, but with storage space
! for {\tt lsize} gridpoints.
!
! {\bf N.B.}:  Though the attribute lists and gridpoint sorting strategy 
! of {\tt iGGrid} is copied to {\tt oGGrid}, the actual values of the 
! attributes are not.
!
! {\bf N.B.}:  It is assumed that {\tt iGGrid} has been initialized.
!
! {\bf N.B.}:  The output {\tt GeneralGrid} {\tt oGGrid} is dynamically
! allocated memory.  When one no longer needs {\tt oGGrid}, one should
! release this space by invoking {\tt GeneralGrid\_clean()}.
!
! !INTERFACE:

 subroutine initgg_(oGGrid, iGGrid, lsize)
!
! !USES:
!
      use m_stdio
      use m_die

      use m_List, only : List
      use m_List, only : List_allocated => allocated
      use m_List, only : List_copy => copy
      use m_List, only : List_nitems => nitem
      use m_List, only : List_nullify => nullify

      use m_AttrVect, only:  AttrVect
      use m_AttrVect, only:  AttrVect_init => init

      implicit none

! !INPUT PARAMETERS: 
!
      type(GeneralGrid), intent(in)  :: iGGrid
      integer, optional, intent(in)  :: lsize

! !OUTPUT PARAMETERS:
!
      type(GeneralGrid), intent(out) :: oGGrid

! !REVISION HISTORY:
!  2May01 - Jay Larson <larson@mcs.anl.gov> - Initial version.
! 13Jun01 - Jay Larson <larson@mcs.anl.gov> - Now, undefined List
!           components of the GeneralGrid iGGrid are no longer 
!           copied to oGGrid.
!  8Aug01 - E.T. Ong <eong@mcs.anl.gov> - changed list assignment(=)
!           to list copy to avoid compiler bugs with pgf90
! 24Jul02 - E.T. Ong <eong@mcs.anl.gov> - updated this init version
!           to correspond with initl_
!  5Aug02 - E. Ong <eong@mcs.anl.gov> - made input argument
!           CoordSortOrder optional to allow for user-defined
!           grid numbering schemes
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::initgg_'
! Number of grid points, number of grid dimensions
  integer :: n, ncoord, norder
! Loop index and Error Flag
  integer :: i, ierr

       ! Start by nullifying everything:

  call List_nullify(oGGrid%coordinate_list)
  call List_nullify(oGGrid%coordinate_sort_order)
  call List_nullify(oGGrid%weight_list)
  call List_nullify(oGGrid%other_list)
  call List_nullify(oGGrid%index_list)
  nullify(oGGrid%descend)

        ! Brief argument check:
        
  ncoord = dims_(iGGrid)    ! dimensionality of the GeneralGrid
  
  if(associated(iGGrid%descend)) then

     if(size(iGGrid%descend) /= ncoord) then ! size mismatch
	call die(myname_,"size(iGGrid%descend) must equal ncoord, &
	         & size(iGGrid%descend) = ", size(iGGrid%descend), &
                 "ncoord = ", ncoord )
     endif

  endif

        ! If iGGrid%descend has been allocated, copy its contents;
        ! allocate and fill oGGrid%descend

  if(associated(iGGrid%descend)) then

     allocate(oGGrid%descend(ncoord), stat=ierr)
     if(ierr /= 0) then
	call die(myname_,"allocate(oGGrid%descend...", ierr)
     endif

     do i=1,ncoord
	oGGrid%descend(i) = iGGrid%descend(i)
     end do

  endif

       ! Copy list data from iGGrid to oGGrid. 

  call List_copy(oGGrid%coordinate_list,iGGrid%coordinate_list)
  if(List_allocated(iGGrid%coordinate_sort_order)) then
     call List_copy(oGGrid%coordinate_sort_order,iGGrid%coordinate_sort_order)
  endif
  if(List_allocated(iGGrid%weight_list)) then
     call List_copy(oGGrid%weight_list,iGGrid%weight_list)
  endif
  if(List_allocated(iGGrid%other_list)) then
     call List_copy(oGGrid%other_list,iGGrid%other_list)
  endif
  if(List_allocated(iGGrid%index_list)) then
     call List_copy(oGGrid%index_list,iGGrid%index_list)
  endif

        ! if lsize is present, use it to set n; if not, set n=0

  n = 0
  if(present(lsize)) n=lsize

       ! Now, initialize oGGrid%data from iGGrid%data, but 
       ! with length n.

  call AttrVect_init(oGGrid%data, iGGrid%data, n)

 end subroutine initgg_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initCartesianSP_ - Initialize a Cartesian GeneralGrid
!
! !DESCRIPTION:
! The routine {\tt initCartesian\_()} creates the storage space for grid point
! coordinates, area and volume weights, and other coordinate data ({\em e.g.}, 
! cell area and volume weights).  The names of the Cartesian axes are supplied
! by the user as a colon-delimitted string in the input {\tt CHARACTER}
! argument {\tt CoordChars}.  For example, a Cartesian grid for Euclidian
! 3-space would have ${\tt CoordChars} = {\tt 'x:y:z'}$.  The user can 
! define named real attributes for spatial weighting data in the input 
! {\tt CHARACTER} argument {\tt WeightChars}.  For example, one could 
! define attributes for Euclidean 3-space length elements by setting 
! ${\tt WeightChars} = {\tt 'dx:dy:dz'}$.  The input {\tt CHARCTER} 
! argument {\tt OtherChars} provides space for defining other real 
! attributes (again as a colon-delimited string of attribute names).
! One can define integer attributes by supplying a colon-delimitted 
! string of names in the input {\tt CHARACTER} argument 
! {\tt IndexChars}.  For example, on could set aside storage space 
! for the {\tt x}-, {\tt y}-, and {\tt z}-indices by setting 
! ${\tt IndexChars} = {\tt 'xIndex:yIndex:zIndex'}$.
!
! Once the storage space in {\tt GGrid} is initialized, The gridpoint 
! coordinates are evaluated using the input arguments {\tt Dims} (the 
! number of points on each coordinate axis) and {\tt AxisData} (the 
! coordinate values on all of the points of all of the axes).  The user 
! presents the axes with each axis stored in a column of {\tt AxisData},
! and the axes are laid out in the same order as the ordering of the 
! axis names in {\tt CoordChars}.  The number of points on each axis 
! is defined by the entries of the input {\tt INTEGER} array 
! {\tt Dims(:)}.  Continuing with the Euclidean 3-space example given 
! above, setting ${\tt Dims(1:3)} = {\tt (256, 256, 128)}$ will result 
! in a Cartesian grid with 256 points in the {\tt x}- and {\tt y}-directions,
! and 128 points in the {\tt z}-direction.  Thus the appropriate dimensions 
! of {\tt AxisData} are 256 rows (the maximum number of axis points among
! all the axes) by 3 columns (the number of physical dimensions).  The 
! {\tt x}-axis points are stored in {\tt AxisData(1:256,1)}, the 
! {\tt y}-axis points are stored in {\tt AxisData(1:256,2)}, and the 
! {\tt z}-axis points are stored in {\tt AxisData(1:128,3)}.
!
! The sorting order of the gridpoints can be either user-defined, or 
! set automatically by MCT.  If the latter is desired, the user must 
! supply the argument {\tt CoordSortOrder}, which defines the 
! lexicographic ordering (by coordinate).  The entries optional input 
! {\tt LOGICAL} array {\tt descend(:)} stipulates whether the ordering 
! with respect to the corresponding key in {\tt CoordChars} is to be
! {\em descending}.  If {\tt CoordChars} is supplied, but {\tt descend(:)} 
! is not, the gridpoint information is placed in {\em ascending} order 
! for each key.  Returning to our Euclidian 3-space example, a choice of  
! ${\tt CoordSortOrder} = {\tt y:x:z}$ and ${\tt descend(1:3)} = 
! ({\tt .TRUE.}, {\tt .FALSE.}, {\tt .FALSE.})$ will result in the entries of 
! {\tt GGrid} being orderd lexicographically by {\tt y} (in descending 
! order), {\tt x} (in ascending order), and {\tt z} (in ascending order).
! Regardless of the gridpoint sorting strategy, MCT will number each of 
! the gridpoints in {\tt GGrid}, storing this information in the integer 
! attribute named {\tt 'GlobGridNum'}.
!
! !INTERFACE:

! Subprogram not used  subroutine initCartesianSP_(GGrid, CoordChars, CoordSortOrder, descend, &
! Subprogram not used                            WeightChars, OtherChars, IndexChars, Dims, &
! Subprogram not used                            AxisData)
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_stdio
! Subprogram not used       use m_die
! Subprogram not used       use m_realkinds,  only : SP
! Subprogram not used 
! Subprogram not used       use m_String,     only : String
! Subprogram not used       use m_String,     only : String_ToChar => ToChar
! Subprogram not used       use m_String,     only : String_clean => clean
! Subprogram not used 
! Subprogram not used       use m_List,     only : List
! Subprogram not used       use m_List,     only : List_init => init
! Subprogram not used       use m_List,     only : List_clean => clean
! Subprogram not used       use m_List,     only : List_nullify => nullify
! Subprogram not used       use m_List,     only : List_append => append
! Subprogram not used       use m_List,     only : List_nitem => nitem
! Subprogram not used       use m_List,     only : List_get => get
! Subprogram not used       use m_List,     only : List_shared => GetSharedListIndices
! Subprogram not used 
! Subprogram not used       use m_AttrVect, only : AttrVect
! Subprogram not used       use m_AttrVect, only : AttrVect_init => init
! Subprogram not used       use m_AttrVect, only : AttrVect_zero => zero
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       character(len=*),                  intent(in)  :: CoordChars
! Subprogram not used       character(len=*),        optional, intent(in)  :: CoordSortOrder
! Subprogram not used       character(len=*),        optional, intent(in)  :: WeightChars
! Subprogram not used       logical, dimension(:),   optional, pointer     :: descend
! Subprogram not used       character(len=*),        optional, intent(in)  :: OtherChars
! Subprogram not used       character(len=*),        optional, intent(in)  :: IndexChars
! Subprogram not used       integer, dimension(:),             pointer     :: Dims
! Subprogram not used       real(SP), dimension(:,:),          pointer     :: AxisData
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       type(GeneralGrid),                 intent(out) :: GGrid
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  7Jun01 - Jay Larson <larson@mcs.anl.gov> - API Specification.
! Subprogram not used ! 12Aug02 - Jay Larson <larson@mcs.anl.gov> - Implementation.
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used !
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::initCartesianSP_'
! Subprogram not used 
! Subprogram not used   type(List) :: IAList, RAList
! Subprogram not used   type(String) :: AxisName
! Subprogram not used   integer, dimension(:), pointer :: &
! Subprogram not used        CoordListIndices, CoordSortOrderIndices
! Subprogram not used   integer :: DimMax, NumDims, NumGridPoints, NumShared
! Subprogram not used   integer :: ierr, iAxis, i, j, k, n, nCycles, nRepeat
! Subprogram not used   integer :: index
! Subprogram not used 
! Subprogram not used        ! Nullify GeneralGrid components
! Subprogram not used 
! Subprogram not used   call List_nullify(GGrid%coordinate_list)
! Subprogram not used   call List_nullify(GGrid%coordinate_sort_order)
! Subprogram not used   call List_nullify(GGrid%weight_list)
! Subprogram not used   call List_nullify(GGrid%other_list)
! Subprogram not used   call List_nullify(GGrid%index_list)
! Subprogram not used   nullify(GGrid%descend)
! Subprogram not used 
! Subprogram not used        ! Sanity check on axis definition arguments:
! Subprogram not used 
! Subprogram not used        ! Ensure each axis has a positive number of points, and
! Subprogram not used        ! determine DimMax, the maximum entry in Dims(:).
! Subprogram not used 
! Subprogram not used   DimMax = 1
! Subprogram not used   do i=1,size(Dims)
! Subprogram not used      if(Dims(i) > DimMax) DimMax = Dims(i)
! Subprogram not used      if(Dims(i) <= 0) then
! Subprogram not used 	write(stderr,'(2a,i8,a,i8)') myname_, &
! Subprogram not used 	     ':: FATAL--illegal number of axis points in Dims(',i,') = ', &
! Subprogram not used 	     Dims(i)
! Subprogram not used 	call die(myname_)
! Subprogram not used      endif
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used        ! Are the definitions of Dims(:) and AxisData(:,:) compatible?
! Subprogram not used        ! The number of elements in Dims(:) should match the number of
! Subprogram not used        ! columns in AxisData(:,:), and the maximum value stored in Dims(:) 
! Subprogram not used        ! (DimMax determined above in this routine) must not exceed the 
! Subprogram not used        ! number of rows in AxisData(:,:).
! Subprogram not used 
! Subprogram not used   if(size(AxisData,2) /= size(Dims)) then
! Subprogram not used      write(stderr,'(4a,i8,a,i8)') myname_, &
! Subprogram not used 	  ':: FATAL-- The number of axes (elements) referenced in Dims(:) ', &
! Subprogram not used 	  'does not equal the number of columns in AxisData(:,:).  ', &
! Subprogram not used 	  'size(Dims) = ',size(Dims),' size(AxisData,2) = ',size(AxisData,2)
! Subprogram not used      call die(myname_)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   if(size(AxisData,1) < DimMax) then
! Subprogram not used      write(stderr,'(4a,i8,a,i8)') myname_, &
! Subprogram not used 	  ':: FATAL-- Maximum number of axis points max(Dims) is ', &
! Subprogram not used 	  'greater than the number of rows in AxisData(:,:).  ', &
! Subprogram not used 	  'max(Dims) = ',DimMax,' size(AxisData,1) = ',size(AxisData,1)
! Subprogram not used      call die(myname_)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! If the LOGICAL descend(:) flags for sorting are present, 
! Subprogram not used        ! make sure that (1) descend is associated, and 
! Subprogram not used        ! (2) CoordSortOrder is also present, and 
! Subprogram not used        ! (3) The size of descend(:) matches the size of Dims(:),
! Subprogram not used        ! both of which correspond to the number of axes on the 
! Subprogram not used        ! Cartesian Grid.
! Subprogram not used 
! Subprogram not used   if(present(descend)) then
! Subprogram not used 
! Subprogram not used      if(.not.associated(descend)) then
! Subprogram not used         call die(myname_,'descend argument must be associated')
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used      if(.not. present(CoordSortOrder)) then
! Subprogram not used         write(stderr,'(4a)') myname_, &
! Subprogram not used              ':: FATAL -- Invocation with the argument descend(:) present ', &
! Subprogram not used              'requires the presence of the argument CoordSortOrder, ', &
! Subprogram not used              'which was not provided.'
! Subprogram not used         call die(myname_, 'Argument CoordSortOrder was not provided')
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used      if(size(descend) /= size(Dims)) then
! Subprogram not used 	write(stderr,'(4a,i8,a,i8)') myname_, &
! Subprogram not used 	     ':: FATAL-- The sizes of the arrays descend(:) and Dims(:) ', &
! Subprogram not used 	     'must match (they both must equal the number of dimensions ', &
! Subprogram not used 	     'of the Cartesian Grid).  size(Dims) = ',size(Dims), &
! Subprogram not used 	     ' size(descend) = ',size(descend)
! Subprogram not used 	call die(myname_,'size of <descend> and <Dims> arguments must match')
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! Initialize GGrid%coordinate_list and use the number of items 
! Subprogram not used        ! in it to set the number of dimensions of the Cartesian 
! Subprogram not used        ! Grid (NumDims):
! Subprogram not used 
! Subprogram not used   call List_init(GGrid%coordinate_list, CoordChars)
! Subprogram not used   
! Subprogram not used   NumDims = List_nitem(GGrid%coordinate_list)
! Subprogram not used 
! Subprogram not used        ! Check the number of arguments
! Subprogram not used 
! Subprogram not used   if(NumDims <= 0) then
! Subprogram not used      write(stderr,*) myname_, &
! Subprogram not used 	  ':: ERROR CoordList is empty!'
! Subprogram not used      call die(myname_,'List_nitem(CoordList) <= 0',NumDims)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! Do the number of coordinate names specified match the number
! Subprogram not used        ! of coordinate axes (i.e., the number of columns in AxisData(:,:))?
! Subprogram not used 
! Subprogram not used   if(NumDims /= size(AxisData,2)) then
! Subprogram not used      write(stderr,'(6a,i8,a,i8)') myname_, &
! Subprogram not used 	  ':: FATAL-- Number of axes specified in argument CoordChars ', &
! Subprogram not used 	  'does not equal the number of axes stored in AxisData(:,:).  ', &
! Subprogram not used 	  'CoordChars = ', CoordChars, &
! Subprogram not used 	  'Number of axes = ',NumDims, &
! Subprogram not used 	  ' size(AxisData,2) = ',size(AxisData,2)
! Subprogram not used      call die(myname_)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! End of argument sanity checks.
! Subprogram not used 
! Subprogram not used        ! Create other List components of GGrid and build REAL
! Subprogram not used        ! and INTEGER attribute lists for the AttrVect GGrid%data
! Subprogram not used 
! Subprogram not used        ! Start off with things *guaranteed* to be in IAList and RAList.
! Subprogram not used        ! The variable GlobGridNum is a CHARACTER parameter inherited 
! Subprogram not used        ! from the declaration section of this module.
! Subprogram not used 
! Subprogram not used   call List_init(IAList, GlobGridNum)
! Subprogram not used   call List_init(RAList, CoordChars)
! Subprogram not used 
! Subprogram not used   if(present(CoordSortOrder)) then
! Subprogram not used 
! Subprogram not used      call List_init(GGrid%coordinate_sort_order, CoordSortOrder)
! Subprogram not used 
! Subprogram not used         ! Check the items in the coordinate list and the 
! Subprogram not used         ! coordinate grid sort keys...they should contain
! Subprogram not used         ! the same items.
! Subprogram not used 
! Subprogram not used      call List_shared(GGrid%coordinate_list,GGrid%coordinate_sort_order, &
! Subprogram not used 	              NumShared,CoordListIndices,CoordSortOrderIndices)
! Subprogram not used 
! Subprogram not used      deallocate(CoordListIndices,CoordSortOrderIndices,stat=ierr)
! Subprogram not used      if(ierr/=0) call die(myname_,'deallocate(CoordListIndices..)',ierr)
! Subprogram not used      
! Subprogram not used      if(NumShared /= NumDims) then
! Subprogram not used 	call die(myname_,'CoordSortOrder must have the same items &
! Subprogram not used 	         & as CoordList',abs(NumDims-NumShared))
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   if(present(WeightChars)) then
! Subprogram not used      call List_init(GGrid%weight_list, WeightChars)
! Subprogram not used      call List_append(RAList, GGrid%weight_list)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   if(present(OtherChars)) then
! Subprogram not used      call List_init(GGrid%other_list, OtherChars)
! Subprogram not used      call List_append(RAList, GGrid%other_list)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   if(present(IndexChars)) then
! Subprogram not used      call List_init(GGrid%index_list, IndexChars)
! Subprogram not used      call List_append(IAList, GGrid%index_list)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! Finally, Initialize GGrid%descend from descend(:).
! Subprogram not used        ! If descend argument is not present, set it to the default .false.
! Subprogram not used 
! Subprogram not used   if(present(CoordSortOrder)) then
! Subprogram not used 
! Subprogram not used      allocate(GGrid%descend(NumDims), stat=ierr)
! Subprogram not used      if(ierr /= 0) call die(myname_,"allocate GGrid%descend...",ierr)
! Subprogram not used 
! Subprogram not used      if(present(descend)) then
! Subprogram not used         do n=1,NumDims
! Subprogram not used            GGrid%descend(n) = descend(n)
! Subprogram not used         end do
! Subprogram not used      else
! Subprogram not used         do n=1,NumDims
! Subprogram not used            GGrid%descend(n) = .FALSE.
! Subprogram not used         end do
! Subprogram not used      endif
! Subprogram not used         
! Subprogram not used   endif ! if(present(CoordSortOrder))...
! Subprogram not used   
! Subprogram not used        ! Compute the total number of grid points in the GeneralGrid.
! Subprogram not used        ! This is merely the product of the elements of Dims(:)
! Subprogram not used 
! Subprogram not used   NumGridPoints = 1
! Subprogram not used   do i=1,NumDims
! Subprogram not used      NumGridPoints = NumGridPoints * Dims(i)
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used        ! Now we are prepared to create GGrid%data:
! Subprogram not used 
! Subprogram not used   call AttrVect_init(GGrid%data, IAList, RAList, NumGridPoints)
! Subprogram not used   call AttrVect_zero(GGrid%data)
! Subprogram not used 
! Subprogram not used        ! Now, store Cartesian gridpoint data, in the order
! Subprogram not used        ! defined by how the user laid out AxisData(:,:)
! Subprogram not used 
! Subprogram not used   do n=1,NumDims
! Subprogram not used 
! Subprogram not used        ! Retrieve first coordinate axis name from GGrid%coordinate_list 
! Subprogram not used        ! (as a String)
! Subprogram not used      call List_get(AxisName, n, GGrid%coordinate_list)
! Subprogram not used 
! Subprogram not used        ! Index this real attribute of GGrid
! Subprogram not used      iAxis = indexRA_(GGrid, String_ToChar(AxisName))
! Subprogram not used 
! Subprogram not used      if(iAxis <= 0) then
! Subprogram not used 	write(stderr,'(4a)') myname_, &
! Subprogram not used 	     ':: REAL Attribute "',String_ToChar(AxisName),'" not found.'
! Subprogram not used 	call die(myname_)
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used        ! Now, clear the String AxisName for use in the next 
! Subprogram not used        ! cycle of this loop:
! Subprogram not used 
! Subprogram not used      call String_clean(AxisName)
! Subprogram not used 
! Subprogram not used        ! Compute the number of times we cycle through the axis 
! Subprogram not used        ! values (nCycles), and the number of times each axis 
! Subprogram not used        ! value is repeated in each cycle (nRepeat)
! Subprogram not used 
! Subprogram not used      nCycles = 1
! Subprogram not used      if(n > 1) then
! Subprogram not used 	do i=1,n-1
! Subprogram not used 	   nCycles = nCycles * Dims(i)
! Subprogram not used 	end do
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used      nRepeat = 1
! Subprogram not used      if(n < NumDims) then
! Subprogram not used 	do i=n+1,NumDims
! Subprogram not used 	   nRepeat = nRepeat * Dims(i)
! Subprogram not used 	end do
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used        ! Loop over the number of cycles for which we run through
! Subprogram not used        ! all the axis points.  Within each cycle, loop over all
! Subprogram not used        ! of the axis points, repeating each value nRepeat times.
! Subprogram not used        ! This produces a set of grid entries that are in 
! Subprogram not used        ! lexicographic order with respect to how the axes are
! Subprogram not used        ! presented to this routine.
! Subprogram not used 
! Subprogram not used      index = 1
! Subprogram not used      do i=1,nCycles
! Subprogram not used 	do j=1,Dims(n)
! Subprogram not used 	   do k=1,nRepeat
! Subprogram not used 	      GGrid%data%rAttr(iAxis,index) = AxisData(j,n)
! Subprogram not used 	      index = index+1
! Subprogram not used 	   end do ! do k=1,nRepeat
! Subprogram not used 	end do ! do j=1,Dims(n)
! Subprogram not used      end do ! do i=1,nCycles
! Subprogram not used 
! Subprogram not used   end do ! do n=1,NumDims...
! Subprogram not used 
! Subprogram not used        ! If the argument CoordSortOrder was supplied, the entries
! Subprogram not used        ! of GGrid will be sorted/permuted with this lexicographic
! Subprogram not used        ! ordering, and the values of the GGrid INTEGER attribute 
! Subprogram not used        ! GlobGridNum will be numbered to reflect this new ordering
! Subprogram not used        ! scheme.
! Subprogram not used 
! Subprogram not used   index = indexIA_(GGrid, GlobGridNum)
! Subprogram not used 
! Subprogram not used   if(present(CoordSortOrder)) then ! Sort permute entries before
! Subprogram not used                                    ! numbering them
! Subprogram not used 
! Subprogram not used      call SortPermute_(GGrid) ! Sort / permute
! Subprogram not used 
! Subprogram not used   endif ! if(present(CoordSortOrder))...
! Subprogram not used 
! Subprogram not used        ! Number the gridpoints based on the AttrVect point index
! Subprogram not used        ! (i.e., the second index in GGrid%data%iAttr)
! Subprogram not used 
! Subprogram not used   do i=1, lsize_(GGrid)
! Subprogram not used      GGrid%data%iAttr(index,i) = i
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used        ! Finally, clean up intermediate Lists
! Subprogram not used 
! Subprogram not used   call List_clean(IAList)
! Subprogram not used   call List_clean(RAList)
! Subprogram not used 
! Subprogram not used  end subroutine initCartesianSP_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
! ----------------------------------------------------------------------
!
! !IROUTINE: initCartesianDP_ - Initialize a Cartesian GeneralGrid
!
! !DESCRIPTION:
! Double Precision version of initCartesianSP_
!
! !INTERFACE:

! Subprogram not used  subroutine initCartesianDP_(GGrid, CoordChars, CoordSortOrder, descend, &
! Subprogram not used                            WeightChars, OtherChars, IndexChars, Dims, &
! Subprogram not used                            AxisData)
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_stdio
! Subprogram not used       use m_die
! Subprogram not used       use m_realkinds,  only : DP
! Subprogram not used 
! Subprogram not used       use m_String,     only : String
! Subprogram not used       use m_String,     only : String_ToChar => ToChar
! Subprogram not used       use m_String,     only : String_clean => clean
! Subprogram not used 
! Subprogram not used       use m_List,     only : List
! Subprogram not used       use m_List,     only : List_init => init
! Subprogram not used       use m_List,     only : List_clean => clean
! Subprogram not used       use m_List,     only : List_nullify => nullify
! Subprogram not used       use m_List,     only : List_append => append
! Subprogram not used       use m_List,     only : List_nitem => nitem
! Subprogram not used       use m_List,     only : List_get => get
! Subprogram not used       use m_List,     only : List_shared => GetSharedListIndices
! Subprogram not used 
! Subprogram not used       use m_AttrVect, only : AttrVect
! Subprogram not used       use m_AttrVect, only : AttrVect_init => init
! Subprogram not used       use m_AttrVect, only : AttrVect_zero => zero
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       character(len=*),                  intent(in)  :: CoordChars
! Subprogram not used       character(len=*),        optional, intent(in)  :: CoordSortOrder
! Subprogram not used       character(len=*),        optional, intent(in)  :: WeightChars
! Subprogram not used       logical, dimension(:),   optional, pointer     :: descend
! Subprogram not used       character(len=*),        optional, intent(in)  :: OtherChars
! Subprogram not used       character(len=*),        optional, intent(in)  :: IndexChars
! Subprogram not used       integer, dimension(:),             pointer     :: Dims
! Subprogram not used       real(DP), dimension(:,:),          pointer     :: AxisData
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       type(GeneralGrid),                 intent(out) :: GGrid
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  7Jun01 - Jay Larson <larson@mcs.anl.gov> - API Specification.
! Subprogram not used ! 12Aug02 - Jay Larson <larson@mcs.anl.gov> - Implementation.
! Subprogram not used ! ______________________________________________________________________
! Subprogram not used !
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::initCartesianDP_'
! Subprogram not used 
! Subprogram not used   type(List) :: IAList, RAList
! Subprogram not used   type(String) :: AxisName
! Subprogram not used   integer, dimension(:), pointer :: &
! Subprogram not used        CoordListIndices, CoordSortOrderIndices
! Subprogram not used   integer :: DimMax, NumDims, NumGridPoints, NumShared
! Subprogram not used   integer :: ierr, iAxis, i, j, k, n, nCycles, nRepeat
! Subprogram not used   integer :: index
! Subprogram not used 
! Subprogram not used        ! Nullify GeneralGrid components
! Subprogram not used 
! Subprogram not used   call List_nullify(GGrid%coordinate_list)
! Subprogram not used   call List_nullify(GGrid%coordinate_sort_order)
! Subprogram not used   call List_nullify(GGrid%weight_list)
! Subprogram not used   call List_nullify(GGrid%other_list)
! Subprogram not used   call List_nullify(GGrid%index_list)
! Subprogram not used   nullify(GGrid%descend)
! Subprogram not used 
! Subprogram not used        ! Sanity check on axis definition arguments:
! Subprogram not used 
! Subprogram not used        ! Ensure each axis has a positive number of points, and
! Subprogram not used        ! determine DimMax, the maximum entry in Dims(:).
! Subprogram not used 
! Subprogram not used   DimMax = 1
! Subprogram not used   do i=1,size(Dims)
! Subprogram not used      if(Dims(i) > DimMax) DimMax = Dims(i)
! Subprogram not used      if(Dims(i) <= 0) then
! Subprogram not used 	write(stderr,'(2a,i8,a,i8)') myname_, &
! Subprogram not used 	     ':: FATAL--illegal number of axis points in Dims(',i,') = ', &
! Subprogram not used 	     Dims(i)
! Subprogram not used 	call die(myname_)
! Subprogram not used      endif
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used        ! Are the definitions of Dims(:) and AxisData(:,:) compatible?
! Subprogram not used        ! The number of elements in Dims(:) should match the number of
! Subprogram not used        ! columns in AxisData(:,:), and the maximum value stored in Dims(:) 
! Subprogram not used        ! (DimMax determined above in this routine) must not exceed the 
! Subprogram not used        ! number of rows in AxisData(:,:).
! Subprogram not used 
! Subprogram not used   if(size(AxisData,2) /= size(Dims)) then
! Subprogram not used      write(stderr,'(4a,i8,a,i8)') myname_, &
! Subprogram not used 	  ':: FATAL-- The number of axes (elements) referenced in Dims(:) ', &
! Subprogram not used 	  'does not equal the number of columns in AxisData(:,:).  ', &
! Subprogram not used 	  'size(Dims) = ',size(Dims),' size(AxisData,2) = ',size(AxisData,2)
! Subprogram not used      call die(myname_)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   if(size(AxisData,1) < DimMax) then
! Subprogram not used      write(stderr,'(4a,i8,a,i8)') myname_, &
! Subprogram not used 	  ':: FATAL-- Maximum number of axis points max(Dims) is ', &
! Subprogram not used 	  'greater than the number of rows in AxisData(:,:).  ', &
! Subprogram not used 	  'max(Dims) = ',DimMax,' size(AxisData,1) = ',size(AxisData,1)
! Subprogram not used      call die(myname_)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! If the LOGICAL descend(:) flags for sorting are present, 
! Subprogram not used        ! make sure that (1) descend is associated, and 
! Subprogram not used        ! (2) CoordSortOrder is also present, and 
! Subprogram not used        ! (3) The size of descend(:) matches the size of Dims(:),
! Subprogram not used        ! both of which correspond to the number of axes on the 
! Subprogram not used        ! Cartesian Grid.
! Subprogram not used 
! Subprogram not used   if(present(descend)) then
! Subprogram not used 
! Subprogram not used      if(.not.associated(descend)) then
! Subprogram not used         call die(myname_,'descend argument must be associated')
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used      if(.not. present(CoordSortOrder)) then
! Subprogram not used         write(stderr,'(4a)') myname_, &
! Subprogram not used              ':: FATAL -- Invocation with the argument descend(:) present ', &
! Subprogram not used              'requires the presence of the argument CoordSortOrder, ', &
! Subprogram not used              'which was not provided.'
! Subprogram not used         call die(myname_, 'Argument CoordSortOrder was not provided')
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used      if(size(descend) /= size(Dims)) then
! Subprogram not used 	write(stderr,'(4a,i8,a,i8)') myname_, &
! Subprogram not used 	     ':: FATAL-- The sizes of the arrays descend(:) and Dims(:) ', &
! Subprogram not used 	     'must match (they both must equal the number of dimensions ', &
! Subprogram not used 	     'of the Cartesian Grid).  size(Dims) = ',size(Dims), &
! Subprogram not used 	     ' size(descend) = ',size(descend)
! Subprogram not used 	call die(myname_,'size of <descend> and <Dims> arguments must match')
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! Initialize GGrid%coordinate_list and use the number of items 
! Subprogram not used        ! in it to set the number of dimensions of the Cartesian 
! Subprogram not used        ! Grid (NumDims):
! Subprogram not used 
! Subprogram not used   call List_init(GGrid%coordinate_list, CoordChars)
! Subprogram not used   
! Subprogram not used   NumDims = List_nitem(GGrid%coordinate_list)
! Subprogram not used 
! Subprogram not used        ! Check the number of arguments
! Subprogram not used 
! Subprogram not used   if(NumDims <= 0) then
! Subprogram not used      write(stderr,*) myname_, &
! Subprogram not used 	  ':: ERROR CoordList is empty!'
! Subprogram not used      call die(myname_,'List_nitem(CoordList) <= 0',NumDims)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! Do the number of coordinate names specified match the number
! Subprogram not used        ! of coordinate axes (i.e., the number of columns in AxisData(:,:))?
! Subprogram not used 
! Subprogram not used   if(NumDims /= size(AxisData,2)) then
! Subprogram not used      write(stderr,'(6a,i8,a,i8)') myname_, &
! Subprogram not used 	  ':: FATAL-- Number of axes specified in argument CoordChars ', &
! Subprogram not used 	  'does not equal the number of axes stored in AxisData(:,:).  ', &
! Subprogram not used 	  'CoordChars = ', CoordChars, &
! Subprogram not used 	  'Number of axes = ',NumDims, &
! Subprogram not used 	  ' size(AxisData,2) = ',size(AxisData,2)
! Subprogram not used      call die(myname_)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! End of argument sanity checks.
! Subprogram not used 
! Subprogram not used        ! Create other List components of GGrid and build REAL
! Subprogram not used        ! and INTEGER attribute lists for the AttrVect GGrid%data
! Subprogram not used 
! Subprogram not used        ! Start off with things *guaranteed* to be in IAList and RAList.
! Subprogram not used        ! The variable GlobGridNum is a CHARACTER parameter inherited 
! Subprogram not used        ! from the declaration section of this module.
! Subprogram not used 
! Subprogram not used   call List_init(IAList, GlobGridNum)
! Subprogram not used   call List_init(RAList, CoordChars)
! Subprogram not used 
! Subprogram not used   if(present(CoordSortOrder)) then
! Subprogram not used 
! Subprogram not used      call List_init(GGrid%coordinate_sort_order, CoordSortOrder)
! Subprogram not used 
! Subprogram not used         ! Check the items in the coordinate list and the 
! Subprogram not used         ! coordinate grid sort keys...they should contain
! Subprogram not used         ! the same items.
! Subprogram not used 
! Subprogram not used      call List_shared(GGrid%coordinate_list,GGrid%coordinate_sort_order, &
! Subprogram not used 	              NumShared,CoordListIndices,CoordSortOrderIndices)
! Subprogram not used 
! Subprogram not used      deallocate(CoordListIndices,CoordSortOrderIndices,stat=ierr)
! Subprogram not used      if(ierr/=0) call die(myname_,'deallocate(CoordListIndices..)',ierr)
! Subprogram not used      
! Subprogram not used      if(NumShared /= NumDims) then
! Subprogram not used 	call die(myname_,'CoordSortOrder must have the same items &
! Subprogram not used 	         & as CoordList',abs(NumDims-NumShared))
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   if(present(WeightChars)) then
! Subprogram not used      call List_init(GGrid%weight_list, WeightChars)
! Subprogram not used      call List_append(RAList, GGrid%weight_list)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   if(present(OtherChars)) then
! Subprogram not used      call List_init(GGrid%other_list, OtherChars)
! Subprogram not used      call List_append(RAList, GGrid%other_list)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   if(present(IndexChars)) then
! Subprogram not used      call List_init(GGrid%index_list, IndexChars)
! Subprogram not used      call List_append(IAList, GGrid%index_list)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! Finally, Initialize GGrid%descend from descend(:).
! Subprogram not used        ! If descend argument is not present, set it to the default .false.
! Subprogram not used 
! Subprogram not used   if(present(CoordSortOrder)) then
! Subprogram not used 
! Subprogram not used      allocate(GGrid%descend(NumDims), stat=ierr)
! Subprogram not used      if(ierr /= 0) call die(myname_,"allocate GGrid%descend...",ierr)
! Subprogram not used 
! Subprogram not used      if(present(descend)) then
! Subprogram not used         do n=1,NumDims
! Subprogram not used            GGrid%descend(n) = descend(n)
! Subprogram not used         end do
! Subprogram not used      else
! Subprogram not used         do n=1,NumDims
! Subprogram not used            GGrid%descend(n) = .FALSE.
! Subprogram not used         end do
! Subprogram not used      endif
! Subprogram not used         
! Subprogram not used   endif ! if(present(CoordSortOrder))...
! Subprogram not used   
! Subprogram not used        ! Compute the total number of grid points in the GeneralGrid.
! Subprogram not used        ! This is merely the product of the elements of Dims(:)
! Subprogram not used 
! Subprogram not used   NumGridPoints = 1
! Subprogram not used   do i=1,NumDims
! Subprogram not used      NumGridPoints = NumGridPoints * Dims(i)
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used        ! Now we are prepared to create GGrid%data:
! Subprogram not used 
! Subprogram not used   call AttrVect_init(GGrid%data, IAList, RAList, NumGridPoints)
! Subprogram not used   call AttrVect_zero(GGrid%data)
! Subprogram not used 
! Subprogram not used        ! Now, store Cartesian gridpoint data, in the order
! Subprogram not used        ! defined by how the user laid out AxisData(:,:)
! Subprogram not used 
! Subprogram not used   do n=1,NumDims
! Subprogram not used 
! Subprogram not used        ! Retrieve first coordinate axis name from GGrid%coordinate_list 
! Subprogram not used        ! (as a String)
! Subprogram not used      call List_get(AxisName, n, GGrid%coordinate_list)
! Subprogram not used 
! Subprogram not used        ! Index this real attribute of GGrid
! Subprogram not used      iAxis = indexRA_(GGrid, String_ToChar(AxisName))
! Subprogram not used 
! Subprogram not used      if(iAxis <= 0) then
! Subprogram not used 	write(stderr,'(4a)') myname_, &
! Subprogram not used 	     ':: REAL Attribute "',String_ToChar(AxisName),'" not found.'
! Subprogram not used 	call die(myname_)
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used        ! Now, clear the String AxisName for use in the next 
! Subprogram not used        ! cycle of this loop:
! Subprogram not used 
! Subprogram not used      call String_clean(AxisName)
! Subprogram not used 
! Subprogram not used        ! Compute the number of times we cycle through the axis 
! Subprogram not used        ! values (nCycles), and the number of times each axis 
! Subprogram not used        ! value is repeated in each cycle (nRepeat)
! Subprogram not used 
! Subprogram not used      nCycles = 1
! Subprogram not used      if(n > 1) then
! Subprogram not used 	do i=1,n-1
! Subprogram not used 	   nCycles = nCycles * Dims(i)
! Subprogram not used 	end do
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used      nRepeat = 1
! Subprogram not used      if(n < NumDims) then
! Subprogram not used 	do i=n+1,NumDims
! Subprogram not used 	   nRepeat = nRepeat * Dims(i)
! Subprogram not used 	end do
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used        ! Loop over the number of cycles for which we run through
! Subprogram not used        ! all the axis points.  Within each cycle, loop over all
! Subprogram not used        ! of the axis points, repeating each value nRepeat times.
! Subprogram not used        ! This produces a set of grid entries that are in 
! Subprogram not used        ! lexicographic order with respect to how the axes are
! Subprogram not used        ! presented to this routine.
! Subprogram not used 
! Subprogram not used      index = 1
! Subprogram not used      do i=1,nCycles
! Subprogram not used 	do j=1,Dims(n)
! Subprogram not used 	   do k=1,nRepeat
! Subprogram not used 	      GGrid%data%rAttr(iAxis,index) = AxisData(j,n)
! Subprogram not used 	      index = index+1
! Subprogram not used 	   end do ! do k=1,nRepeat
! Subprogram not used 	end do ! do j=1,Dims(n)
! Subprogram not used      end do ! do i=1,nCycles
! Subprogram not used 
! Subprogram not used   end do ! do n=1,NumDims...
! Subprogram not used 
! Subprogram not used        ! If the argument CoordSortOrder was supplied, the entries
! Subprogram not used        ! of GGrid will be sorted/permuted with this lexicographic
! Subprogram not used        ! ordering, and the values of the GGrid INTEGER attribute 
! Subprogram not used        ! GlobGridNum will be numbered to reflect this new ordering
! Subprogram not used        ! scheme.
! Subprogram not used 
! Subprogram not used   index = indexIA_(GGrid, GlobGridNum)
! Subprogram not used 
! Subprogram not used   if(present(CoordSortOrder)) then ! Sort permute entries before
! Subprogram not used                                    ! numbering them
! Subprogram not used 
! Subprogram not used      call SortPermute_(GGrid) ! Sort / permute
! Subprogram not used 
! Subprogram not used   endif ! if(present(CoordSortOrder))...
! Subprogram not used 
! Subprogram not used        ! Number the gridpoints based on the AttrVect point index
! Subprogram not used        ! (i.e., the second index in GGrid%data%iAttr)
! Subprogram not used 
! Subprogram not used   do i=1, lsize_(GGrid)
! Subprogram not used      GGrid%data%iAttr(index,i) = i
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used        ! Finally, clean up intermediate Lists
! Subprogram not used 
! Subprogram not used   call List_clean(IAList)
! Subprogram not used   call List_clean(RAList)
! Subprogram not used 
! Subprogram not used  end subroutine initCartesianDP_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initUnstructuredSP_ - Initialize an Unstructured GeneralGrid
!
! !DESCRIPTION:
! This routine creates the storage space for grid point
! coordinates, area/volume weights, and other coordinate data ({\em e.g.}, 
! local cell dimensions), and fills in user-supplied values for the grid 
! point coordinates.  These data are referenced by {\tt List}
! components that are also created by this routine (see the documentation 
! of the declaration section of this module for more details about setting 
! list information).  Each of the input {\tt CHARACTER} arguments is a 
! colon-delimited string of attribute names, each corrsponding to a 
! {\tt List} element of the output {\tt GeneralGrid} argument {\tt GGrid},
! and are summarized in the table below:
!
!\begin{table}[htbp]
!\begin{center}
!\begin{tabular}{|l|l|l|l|}
!\hline
!{\bf Argument} & {\bf Component of {\tt GGrid}} & {\bf Significance} & {\bf Required?} \!\hline

!{\tt CoordChars} & {\tt GGrid\%coordinate\_list} & Dimension Names & Yes \!\hline

!{\tt CoordSortOrder} & {\tt GGrid\%coordinate\_sort\_order} & Grid Point & No \! & & Sorting Keys & \!\hline


!{\tt WeightChars} & {\tt GGrid\%weight\_list} & Grid Cell & No \! & & Length, Area, and & \! & & Volume Weights & \!\hline



!{\tt OtherChars} & {\tt GGrid\%other\_list} & All Other & No \! & & Real Attributes & \!\hline


!{\tt IndexChars} & {\tt GGrid\%index\_list} & All Other & No \! & & Integer Attributes & \!\hline


!\end{tabular}
!\end{center}
!\end{table}
!
! The number of physical dimensions of the grid is set by the user in 
! the input {\tt INTEGER} argument {\tt nDims}, and the number of grid 
! points stored in {\tt GGrid} is set using the input {\tt INTEGER} 
! argument {\tt nPoints}.  The grid point coordinates are input via the 
! {\tt REAL} array {\tt PointData(:)}.  The number of entries in 
! {\tt PointData} must equal the product of {\tt nDims} and {\tt nPoints}.
! The grid points are grouped in {\tt nPoints} consecutive groups of 
! {\tt nDims} entries, with the coordinate values for each point set in 
! the same order as the dimensions are named in the list {\tt CoordChars}.
!
! If a set of sorting keys is supplied in the argument {\tt CoordSortOrder}, 
! the user can control whether the sorting by each key is in descending or 
! ascending order by supplying the input {\tt LOGICAL} array {\tt descend(:)}.  
! By default, all sorting is in {\em ascending} order for each key if the 
! argument {\tt descend} is not provided.
!
! {\bf N.B.}:  The output {\tt GeneralGrid} {\tt GGrid} is dynamically
! allocated memory.  When one no longer needs {\tt GGrid}, one should
! release this space by invoking {\tt clean()} for the {\tt GeneralGrid}.
!
! !INTERFACE:

! Subprogram not used  subroutine initUnstructuredSP_(GGrid, CoordChars, CoordSortOrder, descend, &
! Subprogram not used                               WeightChars, OtherChars, IndexChars, nDims, &
! Subprogram not used                               nPoints, PointData)
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_stdio
! Subprogram not used       use m_die
! Subprogram not used       use m_realkinds,only : SP
! Subprogram not used 
! Subprogram not used       use m_String,   only : String, char
! Subprogram not used       use m_List,     only : List
! Subprogram not used       use m_List,     only : List_init => init
! Subprogram not used       use m_List,     only : List_clean => clean
! Subprogram not used       use m_List,     only : List_nitem => nitem
! Subprogram not used       use m_List,     only : List_nullify => nullify
! Subprogram not used       use m_List,     only : List_copy => copy
! Subprogram not used       use m_List,     only : List_append => append
! Subprogram not used       use m_List,     only : List_shared => GetSharedListIndices
! Subprogram not used       use m_AttrVect, only : AttrVect
! Subprogram not used       use m_AttrVect, only : AttrVect_init => init
! Subprogram not used       use m_AttrVect, only : AttrVect_zero => zero
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       character(len=*),             intent(in) :: CoordChars
! Subprogram not used       character(len=*), optional,   intent(in) :: CoordSortOrder
! Subprogram not used       character(len=*), optional,   intent(in) :: WeightChars
! Subprogram not used       logical, dimension(:), optional, pointer :: descend
! Subprogram not used       character(len=*), optional,   intent(in) :: OtherChars
! Subprogram not used       character(len=*), optional,   intent(in) :: IndexChars
! Subprogram not used       integer,                      intent(in) :: nDims
! Subprogram not used       integer,                      intent(in) :: nPoints
! Subprogram not used       real(SP), dimension(:),       pointer    :: PointData
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       type(GeneralGrid), intent(out)   :: GGrid
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  7Jun01 - Jay Larson <larson@mcs.anl.gov> - API specification.
! Subprogram not used ! 22Aug02 - J. Larson <larson@mcs.anl.gov> - Implementation.
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used !
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::initUnstructuredSP_'
! Subprogram not used 
! Subprogram not used   integer :: i, ierr, index, n, nOffSet, NumShared
! Subprogram not used   integer, dimension(:), pointer :: &
! Subprogram not used        CoordListIndices, CoordSortOrderIndices
! Subprogram not used   type(List) :: IAList, RAList
! Subprogram not used 
! Subprogram not used        ! Nullify all GeneralGrid components
! Subprogram not used 
! Subprogram not used   call List_nullify(GGrid%coordinate_list)
! Subprogram not used   call List_nullify(GGrid%coordinate_sort_order)
! Subprogram not used   call List_nullify(GGrid%weight_list)
! Subprogram not used   call List_nullify(GGrid%other_list)
! Subprogram not used   call List_nullify(GGrid%index_list)
! Subprogram not used   nullify(GGrid%descend)
! Subprogram not used 
! Subprogram not used        ! Sanity checks on input arguments:
! Subprogram not used 
! Subprogram not used        ! If the LOGICAL descend(:) flags for sorting are present, 
! Subprogram not used        ! make sure that (1) it is associated, 
! Subprogram not used        ! (2) CoordSortOrder is also present, and 
! Subprogram not used        ! (3) The size of descend(:) matches the size of Dims(:),
! Subprogram not used        ! both of which correspond to the number of axes on the 
! Subprogram not used        ! Cartesian Grid.
! Subprogram not used 
! Subprogram not used   if(present(descend)) then
! Subprogram not used 
! Subprogram not used      if(.not.associated(descend)) then
! Subprogram not used         call die(myname_,'descend argument must be associated')
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used      if(.not. present(CoordSortOrder)) then
! Subprogram not used         write(stderr,'(4a)') myname_, &
! Subprogram not used              ':: FATAL -- Invocation with the argument descend(:) present ', &
! Subprogram not used              'requires the presence of the argument CoordSortOrder, ', &
! Subprogram not used              'which was not provided.'
! Subprogram not used         call die(myname_,'Argument CoordSortOrder was not provided')
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used      if(present(descend)) then
! Subprogram not used         if(size(descend) /= nDims) then
! Subprogram not used            write(stderr,'(4a,i8,a,i8)') myname_, &
! Subprogram not used                 ':: FATAL-- The size of the array descend(:) and nDims ', &
! Subprogram not used                 'must be equal (they both must equal the number of dimensions ', &
! Subprogram not used                 'of the unstructured Grid).  nDims = ',nDims, &
! Subprogram not used 	     ' size(descend) = ',size(descend)
! Subprogram not used            call die(myname_,'size(descend)/=nDims')
! Subprogram not used         endif
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! Initialize GGrid%coordinate_list and comparethe number of items 
! Subprogram not used        ! to the number of dimensions of the unstructured nDims:
! Subprogram not used 
! Subprogram not used   call List_init(GGrid%coordinate_list, CoordChars)
! Subprogram not used 
! Subprogram not used        ! Check the coordinate_list
! Subprogram not used 
! Subprogram not used   if(nDims /= List_nitem(GGrid%coordinate_list)) then
! Subprogram not used      write(stderr,'(4a,i8,3a,i8)') myname_, &
! Subprogram not used 	  ':: FATAL-- The number of coordinate names supplied in the ', &
! Subprogram not used 	  'argument CoordChars must equal the number of dimensions ', &
! Subprogram not used 	  'specified by the argument nDims.  nDims = ',nDims, &
! Subprogram not used 	  ' CoordChars = ',CoordChars, ' number of dimensions in CoordChars = ', &
! Subprogram not used 	  List_nitem(GGrid%coordinate_list) 
! Subprogram not used      call die(myname_)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   if(nDims <= 0) then
! Subprogram not used      write(stderr,*) myname_, ':: ERROR nDims=0!'
! Subprogram not used      call die(myname_,'nDims <= 0',nDims)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! PointData is a one-dimensional array containing all the gridpoint
! Subprogram not used        ! coordinates.  As such, its size must equal nDims * nPoints.  True?
! Subprogram not used 
! Subprogram not used   if(size(PointData) /= nDims * nPoints) then
! Subprogram not used      write(stderr,'(3a,3(a,i8))') myname_, &
! Subprogram not used 	  ':: FATAL-- The length of the array PointData(:) must match ', &
! Subprogram not used 	  'the product of the input arguments nDims and nPoints.  ', &
! Subprogram not used 	  'nDims = ',nDims, ' nPoints = ',nPoints,&
! Subprogram not used 	  ' size(PointData) = ',size(PointData)
! Subprogram not used      call die(myname_)
! Subprogram not used   endif
! Subprogram not used      
! Subprogram not used        ! End of input argument sanity checks.
! Subprogram not used 
! Subprogram not used        ! Create other List components of GGrid and build REAL
! Subprogram not used        ! and INTEGER attribute lists for the AttrVect GGrid%data
! Subprogram not used 
! Subprogram not used        ! Start off with things *guaranteed* to be in IAList and RAList.
! Subprogram not used        ! The variable GlobGridNum is a CHARACTER parameter inherited 
! Subprogram not used        ! from the declaration section of this module.
! Subprogram not used 
! Subprogram not used   call List_init(IAList, GlobGridNum)
! Subprogram not used   call List_init(RAList, CoordChars)
! Subprogram not used 
! Subprogram not used   if(present(CoordSortOrder)) then
! Subprogram not used 
! Subprogram not used      call List_init(GGrid%coordinate_sort_order, CoordSortOrder)
! Subprogram not used 
! Subprogram not used      call List_shared(GGrid%coordinate_list,GGrid%coordinate_sort_order, &
! Subprogram not used 	              NumShared,CoordListIndices,CoordSortOrderIndices)
! Subprogram not used 
! Subprogram not used      deallocate(CoordListIndices,CoordSortOrderIndices,stat=ierr)
! Subprogram not used      if(ierr/=0) call die(myname_,'deallocate(CoordListIndices..)',ierr)
! Subprogram not used      
! Subprogram not used      if(NumShared /= nDims) then
! Subprogram not used 	call die(myname_,'CoordSortOrder must have the same items &
! Subprogram not used 	         & as CoordList',abs(nDims-NumShared))
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   if(present(WeightChars)) then
! Subprogram not used      call List_init(GGrid%weight_list, WeightChars)
! Subprogram not used      call List_append(RAList, GGrid%weight_list)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   if(present(OtherChars)) then
! Subprogram not used      call List_init(GGrid%other_list, OtherChars)
! Subprogram not used      call List_append(RAList, GGrid%other_list)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   if(present(IndexChars)) then
! Subprogram not used      call List_init(GGrid%index_list, IndexChars)
! Subprogram not used      call List_append(IAList, GGrid%index_list)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! Initialize GGrid%descend from descend(:).
! Subprogram not used        ! If descend argument is not present, set it to the default .false.
! Subprogram not used 
! Subprogram not used   if(present(CoordSortOrder)) then
! Subprogram not used 
! Subprogram not used      allocate(GGrid%descend(nDims), stat=ierr)
! Subprogram not used      if(ierr /= 0) call die(myname_,"allocate GGrid%descend...",ierr)
! Subprogram not used 
! Subprogram not used      if(present(descend)) then
! Subprogram not used         do n=1,nDims
! Subprogram not used            GGrid%descend(n) = descend(n)
! Subprogram not used         end do
! Subprogram not used      else
! Subprogram not used         do n=1,nDims
! Subprogram not used            GGrid%descend(n) = .FALSE.
! Subprogram not used         end do
! Subprogram not used      endif
! Subprogram not used         
! Subprogram not used   endif ! if(present(CoordSortOrder))...
! Subprogram not used   
! Subprogram not used        ! Create Grid attribute data storage AttrVect GGrid%data:
! Subprogram not used 
! Subprogram not used   call AttrVect_init(GGrid%data, IAList, RAList, nPoints)
! Subprogram not used   call AttrVect_zero(GGrid%data)
! Subprogram not used 
! Subprogram not used        ! Load up gridpoint coordinate data into GGrid%data.
! Subprogram not used        ! Given how we've set up the real attributes of GGrid%data,
! Subprogram not used        ! we have guaranteed the first nDims real attributes are 
! Subprogram not used        ! the gridpoint coordinates.
! Subprogram not used 
! Subprogram not used   do n=1,nPoints
! Subprogram not used      nOffSet = (n-1) * nDims
! Subprogram not used      do i=1,nDims
! Subprogram not used 	GGrid%data%rAttr(i,n) = PointData(nOffset + i)
! Subprogram not used      end do
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used        ! If the argument CoordSortOrder was supplied, the entries
! Subprogram not used        ! of GGrid will be sorted/permuted with this lexicographic
! Subprogram not used        ! ordering, and the values of the GGrid INTEGER attribute 
! Subprogram not used        ! GlobGridNum will be numbered to reflect this new ordering
! Subprogram not used        ! scheme.
! Subprogram not used 
! Subprogram not used   index = indexIA_(GGrid, GlobGridNum)
! Subprogram not used 
! Subprogram not used   if(present(CoordSortOrder)) then ! Sort permute entries before
! Subprogram not used                                    ! numbering them
! Subprogram not used 
! Subprogram not used      call SortPermute_(GGrid) ! Sort / permute
! Subprogram not used 
! Subprogram not used   endif ! if(present(CoordSortOrder))...
! Subprogram not used 
! Subprogram not used        ! Number the gridpoints based on the AttrVect point index
! Subprogram not used        ! (i.e., the second index in GGrid%data%iAttr)
! Subprogram not used 
! Subprogram not used   do i=1, lsize_(GGrid)
! Subprogram not used      GGrid%data%iAttr(index,i) = i
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used        ! Clean up temporary allocated structures:
! Subprogram not used 
! Subprogram not used   call List_clean(IAList)
! Subprogram not used   call List_clean(RAList)
! Subprogram not used 
! Subprogram not used  end subroutine initUnstructuredSP_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
! ----------------------------------------------------------------------
!
! !IROUTINE: initUnstructuredDP_ - Initialize an Unstructured GeneralGrid
!
! !DESCRIPTION:
! Double precision version of initUnstructuredSP_
!
! !INTERFACE:

! Subprogram not used  subroutine initUnstructuredDP_(GGrid, CoordChars, CoordSortOrder, descend, &
! Subprogram not used                               WeightChars, OtherChars, IndexChars, nDims, &
! Subprogram not used                               nPoints, PointData)
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_stdio
! Subprogram not used       use m_die
! Subprogram not used       use m_realkinds,only : DP
! Subprogram not used 
! Subprogram not used       use m_String,   only : String, char
! Subprogram not used       use m_List,     only : List
! Subprogram not used       use m_List,     only : List_init => init
! Subprogram not used       use m_List,     only : List_clean => clean
! Subprogram not used       use m_List,     only : List_nitem => nitem
! Subprogram not used       use m_List,     only : List_nullify => nullify
! Subprogram not used       use m_List,     only : List_copy => copy
! Subprogram not used       use m_List,     only : List_append => append
! Subprogram not used       use m_List,     only : List_shared => GetSharedListIndices
! Subprogram not used       use m_AttrVect, only : AttrVect
! Subprogram not used       use m_AttrVect, only : AttrVect_init => init
! Subprogram not used       use m_AttrVect, only : AttrVect_zero => zero
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       character(len=*),             intent(in) :: CoordChars
! Subprogram not used       character(len=*), optional,   intent(in) :: CoordSortOrder
! Subprogram not used       character(len=*), optional,   intent(in) :: WeightChars
! Subprogram not used       logical, dimension(:), optional, pointer :: descend
! Subprogram not used       character(len=*), optional,   intent(in) :: OtherChars
! Subprogram not used       character(len=*), optional,   intent(in) :: IndexChars
! Subprogram not used       integer,                      intent(in) :: nDims
! Subprogram not used       integer,                      intent(in) :: nPoints
! Subprogram not used       real(DP), dimension(:),       pointer    :: PointData
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used       type(GeneralGrid), intent(out)   :: GGrid
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used !  7Jun01 - Jay Larson <larson@mcs.anl.gov> - API specification.
! Subprogram not used ! 22Aug02 - J. Larson <larson@mcs.anl.gov> - Implementation.
! Subprogram not used ! ______________________________________________________________________
! Subprogram not used !
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::initUnstructuredDP_'
! Subprogram not used 
! Subprogram not used   integer :: i, ierr, index, n, nOffSet, NumShared
! Subprogram not used   integer, dimension(:), pointer :: &
! Subprogram not used        CoordListIndices, CoordSortOrderIndices
! Subprogram not used   type(List) :: IAList, RAList
! Subprogram not used 
! Subprogram not used        ! Nullify all GeneralGrid components
! Subprogram not used 
! Subprogram not used   call List_nullify(GGrid%coordinate_list)
! Subprogram not used   call List_nullify(GGrid%coordinate_sort_order)
! Subprogram not used   call List_nullify(GGrid%weight_list)
! Subprogram not used   call List_nullify(GGrid%other_list)
! Subprogram not used   call List_nullify(GGrid%index_list)
! Subprogram not used   nullify(GGrid%descend)
! Subprogram not used 
! Subprogram not used        ! Sanity checks on input arguments:
! Subprogram not used 
! Subprogram not used        ! If the LOGICAL descend(:) flags for sorting are present, 
! Subprogram not used        ! make sure that (1) it is associated, 
! Subprogram not used        ! (2) CoordSortOrder is also present, and 
! Subprogram not used        ! (3) The size of descend(:) matches the size of Dims(:),
! Subprogram not used        ! both of which correspond to the number of axes on the 
! Subprogram not used        ! Cartesian Grid.
! Subprogram not used 
! Subprogram not used   if(present(descend)) then
! Subprogram not used 
! Subprogram not used      if(.not.associated(descend)) then
! Subprogram not used         call die(myname_,'descend argument must be associated')
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used      if(.not. present(CoordSortOrder)) then
! Subprogram not used         write(stderr,'(4a)') myname_, &
! Subprogram not used              ':: FATAL -- Invocation with the argument descend(:) present ', &
! Subprogram not used              'requires the presence of the argument CoordSortOrder, ', &
! Subprogram not used              'which was not provided.'
! Subprogram not used         call die(myname_,'Argument CoordSortOrder was not provided')
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used      if(present(descend)) then
! Subprogram not used         if(size(descend) /= nDims) then
! Subprogram not used            write(stderr,'(4a,i8,a,i8)') myname_, &
! Subprogram not used                 ':: FATAL-- The size of the array descend(:) and nDims ', &
! Subprogram not used                 'must be equal (they both must equal the number of dimensions ', &
! Subprogram not used                 'of the unstructured Grid).  nDims = ',nDims, &
! Subprogram not used 	     ' size(descend) = ',size(descend)
! Subprogram not used            call die(myname_,'size(descend)/=nDims')
! Subprogram not used         endif
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! Initialize GGrid%coordinate_list and comparethe number of items 
! Subprogram not used        ! to the number of dimensions of the unstructured nDims:
! Subprogram not used 
! Subprogram not used   call List_init(GGrid%coordinate_list, CoordChars)
! Subprogram not used 
! Subprogram not used        ! Check the coordinate_list
! Subprogram not used 
! Subprogram not used   if(nDims /= List_nitem(GGrid%coordinate_list)) then
! Subprogram not used      write(stderr,'(4a,i8,3a,i8)') myname_, &
! Subprogram not used 	  ':: FATAL-- The number of coordinate names supplied in the ', &
! Subprogram not used 	  'argument CoordChars must equal the number of dimensions ', &
! Subprogram not used 	  'specified by the argument nDims.  nDims = ',nDims, &
! Subprogram not used 	  ' CoordChars = ',CoordChars, ' number of dimensions in CoordChars = ', &
! Subprogram not used 	  List_nitem(GGrid%coordinate_list) 
! Subprogram not used      call die(myname_)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   if(nDims <= 0) then
! Subprogram not used      write(stderr,*) myname_, ':: ERROR nDims=0!'
! Subprogram not used      call die(myname_,'nDims <= 0',nDims)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! PointData is a one-dimensional array containing all the gridpoint
! Subprogram not used        ! coordinates.  As such, its size must equal nDims * nPoints.  True?
! Subprogram not used 
! Subprogram not used   if(size(PointData) /= nDims * nPoints) then
! Subprogram not used      write(stderr,'(3a,3(a,i8))') myname_, &
! Subprogram not used 	  ':: FATAL-- The length of the array PointData(:) must match ', &
! Subprogram not used 	  'the product of the input arguments nDims and nPoints.  ', &
! Subprogram not used 	  'nDims = ',nDims, ' nPoints = ',nPoints,&
! Subprogram not used 	  ' size(PointData) = ',size(PointData)
! Subprogram not used      call die(myname_)
! Subprogram not used   endif
! Subprogram not used      
! Subprogram not used        ! End of input argument sanity checks.
! Subprogram not used 
! Subprogram not used        ! Create other List components of GGrid and build REAL
! Subprogram not used        ! and INTEGER attribute lists for the AttrVect GGrid%data
! Subprogram not used 
! Subprogram not used        ! Start off with things *guaranteed* to be in IAList and RAList.
! Subprogram not used        ! The variable GlobGridNum is a CHARACTER parameter inherited 
! Subprogram not used        ! from the declaration section of this module.
! Subprogram not used 
! Subprogram not used   call List_init(IAList, GlobGridNum)
! Subprogram not used   call List_init(RAList, CoordChars)
! Subprogram not used 
! Subprogram not used   if(present(CoordSortOrder)) then
! Subprogram not used 
! Subprogram not used      call List_init(GGrid%coordinate_sort_order, CoordSortOrder)
! Subprogram not used 
! Subprogram not used      call List_shared(GGrid%coordinate_list,GGrid%coordinate_sort_order, &
! Subprogram not used 	              NumShared,CoordListIndices,CoordSortOrderIndices)
! Subprogram not used 
! Subprogram not used      deallocate(CoordListIndices,CoordSortOrderIndices,stat=ierr)
! Subprogram not used      if(ierr/=0) call die(myname_,'deallocate(CoordListIndices..)',ierr)
! Subprogram not used      
! Subprogram not used      if(NumShared /= nDims) then
! Subprogram not used 	call die(myname_,'CoordSortOrder must have the same items &
! Subprogram not used 	         & as CoordList',abs(nDims-NumShared))
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   if(present(WeightChars)) then
! Subprogram not used      call List_init(GGrid%weight_list, WeightChars)
! Subprogram not used      call List_append(RAList, GGrid%weight_list)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   if(present(OtherChars)) then
! Subprogram not used      call List_init(GGrid%other_list, OtherChars)
! Subprogram not used      call List_append(RAList, GGrid%other_list)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   if(present(IndexChars)) then
! Subprogram not used      call List_init(GGrid%index_list, IndexChars)
! Subprogram not used      call List_append(IAList, GGrid%index_list)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! Initialize GGrid%descend from descend(:).
! Subprogram not used        ! If descend argument is not present, set it to the default .false.
! Subprogram not used 
! Subprogram not used   if(present(CoordSortOrder)) then
! Subprogram not used 
! Subprogram not used      allocate(GGrid%descend(nDims), stat=ierr)
! Subprogram not used      if(ierr /= 0) call die(myname_,"allocate GGrid%descend...",ierr)
! Subprogram not used 
! Subprogram not used      if(present(descend)) then
! Subprogram not used         do n=1,nDims
! Subprogram not used            GGrid%descend(n) = descend(n)
! Subprogram not used         end do
! Subprogram not used      else
! Subprogram not used         do n=1,nDims
! Subprogram not used            GGrid%descend(n) = .FALSE.
! Subprogram not used         end do
! Subprogram not used      endif
! Subprogram not used         
! Subprogram not used   endif ! if(present(CoordSortOrder))...
! Subprogram not used   
! Subprogram not used        ! Create Grid attribute data storage AttrVect GGrid%data:
! Subprogram not used 
! Subprogram not used   call AttrVect_init(GGrid%data, IAList, RAList, nPoints)
! Subprogram not used   call AttrVect_zero(GGrid%data)
! Subprogram not used 
! Subprogram not used        ! Load up gridpoint coordinate data into GGrid%data.
! Subprogram not used        ! Given how we've set up the real attributes of GGrid%data,
! Subprogram not used        ! we have guaranteed the first nDims real attributes are 
! Subprogram not used        ! the gridpoint coordinates.
! Subprogram not used 
! Subprogram not used   do n=1,nPoints
! Subprogram not used      nOffSet = (n-1) * nDims
! Subprogram not used      do i=1,nDims
! Subprogram not used 	GGrid%data%rAttr(i,n) = PointData(nOffset + i)
! Subprogram not used      end do
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used        ! If the argument CoordSortOrder was supplied, the entries
! Subprogram not used        ! of GGrid will be sorted/permuted with this lexicographic
! Subprogram not used        ! ordering, and the values of the GGrid INTEGER attribute 
! Subprogram not used        ! GlobGridNum will be numbered to reflect this new ordering
! Subprogram not used        ! scheme.
! Subprogram not used 
! Subprogram not used   index = indexIA_(GGrid, GlobGridNum)
! Subprogram not used 
! Subprogram not used   if(present(CoordSortOrder)) then ! Sort permute entries before
! Subprogram not used                                    ! numbering them
! Subprogram not used 
! Subprogram not used      call SortPermute_(GGrid) ! Sort / permute
! Subprogram not used 
! Subprogram not used   endif ! if(present(CoordSortOrder))...
! Subprogram not used 
! Subprogram not used        ! Number the gridpoints based on the AttrVect point index
! Subprogram not used        ! (i.e., the second index in GGrid%data%iAttr)
! Subprogram not used 
! Subprogram not used   do i=1, lsize_(GGrid)
! Subprogram not used      GGrid%data%iAttr(index,i) = i
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used        ! Clean up temporary allocated structures:
! Subprogram not used 
! Subprogram not used   call List_clean(IAList)
! Subprogram not used   call List_clean(RAList)
! Subprogram not used 
! Subprogram not used  end subroutine initUnstructuredDP_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - Destroy a GeneralGrid
!
! !DESCRIPTION:
! This routine deallocates all attribute storage space for the input/output
! {\tt GeneralGrid} argument {\tt GGrid}, and destroys all of its {\tt List}
! components and sorting flags.  The success (failure) of this operation is
! signified by the zero (non-zero) value of the optional {\tt INTEGER} 
! output argument {\tt stat}.
!
! !INTERFACE:

    subroutine clean_(GGrid, stat)
!
! !USES:
!
      use m_stdio
      use m_die

      use m_List,     only : List_clean => clean
      use m_List,     only : List_allocated => allocated
      use m_AttrVect, only : AttrVect_clean => clean

      implicit none

! !INPUT/OUTPUT PARAMETERS: 
!
      type(GeneralGrid), intent(inout) :: GGrid
      integer, optional, intent(out)   :: stat

! !REVISION HISTORY:
! 25Sep00 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
! 20Mar01 - J.W. Larson <larson@mcs.anl.gov> - complete version.
!  1Mar01 - E.T. Ong <eong@mcs.anl.gov> - removed dies to prevent
!           crashes when cleaning uninitialized attrvects. Added
!           optional stat argument.
!  5Aug02 - E. Ong <eong@mcs.anl.gov> - a more rigorous revision
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ierr

  if(present(stat)) then

     stat=0
     call AttrVect_clean(GGrid%data,ierr)
     if(ierr/=0) stat=ierr

     call List_clean(GGrid%coordinate_list,ierr)
     if(ierr/=0) stat=ierr

     if(List_allocated(GGrid%coordinate_sort_order)) then
	call List_clean(GGrid%coordinate_sort_order,ierr)
	if(ierr/=0) stat=ierr
     endif

     if(List_allocated(GGrid%weight_list)) then
	call List_clean(GGrid%weight_list,ierr)
	if(ierr/=0) stat=ierr
     endif

     if(List_allocated(GGrid%other_list)) then 
	call List_clean(GGrid%other_list,ierr)
	if(ierr/=0) stat=ierr
     endif

     if(List_allocated(GGrid%index_list)) then 
	call List_clean(GGrid%index_list,ierr)
	if(ierr/=0) stat=ierr
     endif

     if(associated(GGrid%descend)) then
	deallocate(GGrid%descend, stat=ierr)
	if(ierr/=0) stat=ierr
     endif

  else

     call AttrVect_clean(GGrid%data)

     call List_clean(GGrid%coordinate_list)

     if(List_allocated(GGrid%coordinate_sort_order)) then
	call List_clean(GGrid%coordinate_sort_order)
     endif

     if(List_allocated(GGrid%weight_list)) then
	call List_clean(GGrid%weight_list)
     endif

     if(List_allocated(GGrid%other_list)) then 
	call List_clean(GGrid%other_list)
     endif

     if(List_allocated(GGrid%index_list)) then 
	call List_clean(GGrid%index_list)
     endif

     if(associated(GGrid%descend)) then
	deallocate(GGrid%descend, stat=ierr)
	if(ierr/=0) call die(myname_,'deallocate(GGrid%descend)',ierr) 
     endif

  endif

 end subroutine clean_

!BOP -------------------------------------------------------------------
!
! !IROUTINE: zero_ - Set GeneralGrid Data to Zero 
!
! !DESCRIPTION:
! This routine sets all of the point values of the integer and real 
! attributes of an the input/output {\tt GeneralGrid} argument {\tt GGrid}
! to zero.  The default action is to set the values of all the real and  
! integer attributes to zero.
!
! !INTERFACE:

! Subprogram not used  subroutine zero_(GGrid, zeroReals, zeroInts)
! Subprogram not used 
! Subprogram not used ! !USES:
! Subprogram not used 
! Subprogram not used 
! Subprogram not used      use m_die,only     : die
! Subprogram not used      use m_stdio,only   : stderr
! Subprogram not used 
! Subprogram not used      use m_AttrVect, only : AttrVect_zero => zero
! Subprogram not used 
! Subprogram not used      implicit none
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used !
! Subprogram not used      type(GeneralGrid),    intent(INOUT) :: GGrid
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used      logical, optional, intent(IN)    :: zeroReals
! Subprogram not used      logical, optional, intent(IN)    :: zeroInts
! Subprogram not used 
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 11May08 - R. Jacob <jacob@mcs.anl.gov> - initial prototype/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::zero_'
! Subprogram not used 
! Subprogram not used   logical myZeroReals, myZeroInts
! Subprogram not used 
! Subprogram not used   if(present(zeroReals)) then
! Subprogram not used      myZeroReals = zeroReals
! Subprogram not used   else
! Subprogram not used      myZeroReals = .TRUE.
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   if(present(zeroInts)) then
! Subprogram not used      myZeroInts = zeroInts
! Subprogram not used   else
! Subprogram not used      myZeroInts = .TRUE.
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   call AttrVect_zero(GGrid%data,myZeroReals,myZeroInts)
! Subprogram not used 
! Subprogram not used  end subroutine zero_


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: dims_ - Return the Dimensionality of a GeneralGrid
!
! !DESCRIPTION:
! This {\tt INTEGER} function returns the number of physical dimensions 
! of the input {\tt GeneralGrid} argument {\tt GGrid}.
!
! !INTERFACE:

 integer function dims_(GGrid)
!
! !USES:
!
      use m_stdio
      use m_die

      use m_List,     only : List_nitem => nitem

      implicit none

! !INPUT PARAMETERS: 
!
      type(GeneralGrid), intent(in)  :: GGrid

! !REVISION HISTORY:
! 15Jan01 - Jay Larson <larson@mcs.anl.gov> - initial version
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::dims_'


 dims_ = List_nitem(GGrid%coordinate_list)

 if(dims_<=0) then
    call die(myname_,"GGrid has zero dimensions",dims_)
 endif

 end function dims_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexIA - Index an Integer Attribute
!
! !DESCRIPTION:
! This function returns an {\tt INTEGER}, corresponding to the location 
! of an integer attribute within the input {\tt GeneralGrid} argument 
! {\tt GGrid}.  For example, every {\tt GGrid} has at least one integer 
! attribute (namely the global gridpoint index {\tt 'GlobGridNum'}).
! The array of integer values for the attribute {\tt 'GlobGridNum'} is 
! stored in 
! \begin{verbatim}
! {\tt GGrid%data%iAttr(indexIA_(GGrid,'GlobGridNum'),:)}.
! \end{verbatim}
! If {\tt indexIA\_()} is unable to match {\tt item} to any of the integer
! attributes present in {\tt GGrid}, the resulting value is zero which is 
! equivalent to an error.  The optional input {\tt CHARACTER} arguments 
! {\tt perrWith} and {\tt dieWith} control how such errors are handled.  
! Below are the rules how error handling is controlled by using 
! {\tt perrWith} and {\tt dieWith}:
! \begin{enumerate}
! \item if neither {\tt perrWith} nor {\tt dieWith} are present, 
! {\tt indexIA\_()} terminates execution with an internally generated
! error message;
! \item if {\tt perrWith} is present, but {\tt dieWith} is not, an error 
! message is written to {\tt stderr} incorporating user-supplied 
! traceback information stored in the argument {\tt perrWith};
! \item if {\tt dieWith} is present, execution terminates with an error 
! message written to {\tt stderr} that incorporates user-supplied 
! traceback information stored in the argument {\tt dieWith}; and 
! \item if both {\tt perrWith} and {\tt dieWith} are present, execution 
! terminates with an error message using {\tt dieWith}, and the argument
! {\tt perrWith} is ignored.
! \end{enumerate}
!
! !INTERFACE:

! Subprogram not used  integer function indexIA_(GGrid, item, perrWith, dieWith)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_die
! Subprogram not used       use m_stdio
! Subprogram not used 
! Subprogram not used       use m_String, only : String
! Subprogram not used       use m_String, only : String_init => init
! Subprogram not used       use m_String, only : String_clean => clean
! Subprogram not used       use m_String, only : String_ToChar => ToChar
! Subprogram not used 
! Subprogram not used       use m_TraceBack, only : GenTraceBackString
! Subprogram not used 
! Subprogram not used       use m_AttrVect,     only : AttrVect_indexIA => indexIA
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS: 
! Subprogram not used !
! Subprogram not used       type(GeneralGrid),          intent(in) :: GGrid
! Subprogram not used       character(len=*),           intent(in) :: item
! Subprogram not used       character(len=*), optional, intent(in) :: perrWith
! Subprogram not used       character(len=*), optional, intent(in) :: dieWith
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 15Jan01 - Jay Larson <larson@mcs.anl.gov> - Initial version.
! Subprogram not used ! 27Mar02 - Jay Larson <larson@mcs.anl.gov> - Cleaned up error
! Subprogram not used !           handling logic.
! Subprogram not used !  2Aug02 - Jay Larson <larson@mcs.anl.gov> - Further refinement
! Subprogram not used !           of error handling.
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used !
! Subprogram not used 
! Subprogram not used  character(len=*), parameter :: myname_=myname//'::indexIA_'
! Subprogram not used 
! Subprogram not used  type(String) :: myTrace
! Subprogram not used 
! Subprogram not used        ! Generate a traceback String
! Subprogram not used 
! Subprogram not used   if(present(dieWith)) then
! Subprogram not used      call GenTraceBackString(myTrace, dieWith, myname_)
! Subprogram not used   else
! Subprogram not used      if(present(perrWith)) then
! Subprogram not used         call GenTraceBackString(myTrace, perrWith, myname_)
! Subprogram not used      else
! Subprogram not used         call GenTraceBackString(myTrace, myname_)
! Subprogram not used      endif
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! Call AttrVect_indexIA() accordingly:
! Subprogram not used 
! Subprogram not used   if( present(dieWith) .or. &
! Subprogram not used      ((.not. present(dieWith)) .and. (.not. present(perrWith))) ) then
! Subprogram not used      indexIA_ = AttrVect_indexIA(GGrid%data, item, &
! Subprogram not used                                  dieWith=String_ToChar(myTrace))
! Subprogram not used   else  ! perrWith but no dieWith case
! Subprogram not used      indexIA_ = AttrVect_indexIA(GGrid%data, item, &
! Subprogram not used                    perrWith=String_ToChar(myTrace))
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   call String_clean(myTrace)
! Subprogram not used 
! Subprogram not used  end function indexIA_


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexRA - Index a Real Attribute
!
! !DESCRIPTION:

! This function returns an {\tt INTEGER}, corresponding to the location 
! of an integer attribute within the input {\tt GeneralGrid} argument 
! {\tt GGrid}.  For example, every {\tt GGrid} has at least one integer 
! attribute (namely the global gridpoint index {\tt 'GlobGridNum'}).
! The array of integer values for the attribute {\tt 'GlobGridNum'} is 
! stored in 
! \begin{verbatim}
! {\tt GGrid%data%iAttr(indexRA_(GGrid,'GlobGridNum'),:)}.
! \end{verbatim}
! If {\tt indexRA\_()} is unable to match {\tt item} to any of the integer
! attributes present in {\tt GGrid}, the resulting value is zero which is 
! equivalent to an error.  The optional input {\tt CHARACTER} arguments 
! {\tt perrWith} and {\tt dieWith} control how such errors are handled.  
! Below are the rules how error handling is controlled by using 
! {\tt perrWith} and {\tt dieWith}:
! \begin{enumerate}
! \item if neither {\tt perrWith} nor {\tt dieWith} are present, 
! {\tt indexRA\_()} terminates execution with an internally generated
! error message;
! \item if {\tt perrWith} is present, but {\tt dieWith} is not, an error 
! message is written to {\tt stderr} incorporating user-supplied 
! traceback information stored in the argument {\tt perrWith};
! \item if {\tt dieWith} is present, execution terminates with an error 
! message written to {\tt stderr} that incorporates user-supplied 
! traceback information stored in the argument {\tt dieWith}; and 
! \item if both {\tt perrWith} and {\tt dieWith} are present, execution 
! terminates with an error message using {\tt dieWith}, and the argument
! {\tt perrWith} is ignored.
! \end{enumerate}
!
! !INTERFACE:

 integer function indexRA_(GGrid, item, perrWith, dieWith)
!
! !USES:
!
      use m_stdio
      use m_die

      use m_String, only : String
      use m_String, only : String_init => init
      use m_String, only : String_clean => clean
      use m_String, only : String_ToChar => ToChar

      use m_TraceBack, only : GenTraceBackString

      use m_AttrVect,     only : AttrVect_indexRA => indexRA

      implicit none

! !INPUT PARAMETERS: 
!
      type(GeneralGrid),          intent(in)  :: GGrid
      character(len=*),           intent(in)  :: item
      character(len=*), optional, intent(in) :: perrWith
      character(len=*), optional, intent(in) :: dieWith

! !REVISION HISTORY:
! 15Jan01 - Jay Larson <larson@mcs.anl.gov> - Initial version.
! 27Mar02 - Jay Larson <larson@mcs.anl.gov> - Cleaned up error
!           handling logic.
!EOP ___________________________________________________________________
!
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
     indexRA_ = AttrVect_indexRA(GGrid%data, item, &
                                 dieWith=String_ToChar(myTrace))
  else  ! perrWith but no dieWith case
     indexRA_ = AttrVect_indexRA(GGrid%data, item, &
                   perrWith=String_ToChar(myTrace))
  endif

  call String_clean(myTrace)

 end function indexRA_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: lsize - Number of Grid Points
!
! !DESCRIPTION:
! This {\tt INTEGER} function returns the number of grid points stored 
! in the input {\tt GeneralGrid} argument {\tt GGrid}.  Note that the 
! value returned will be the number of points stored on a local process 
! in the case of a distributed {\tt GeneralGrid}.
!
! !INTERFACE:

 integer function lsize_(GGrid)
!
! !USES:
!
      use m_List,     only : List
      use m_List,     only : List_allocated => allocated
      use m_AttrVect, only : AttrVect_lsize => lsize
      use m_die,      only : die    
      

      implicit none

! !INPUT PARAMETERS: 
!
      type(GeneralGrid), intent(in)  :: GGrid

! !REVISION HISTORY:
! 15Jan01 - Jay Larson <larson@mcs.anl.gov> - Initial version.
! 27Mar02 - Jay Larson <larson@mcs.anl.gov> - slight logic change.
! 27Mar02 - Jay Larson <larson@mcs.anl.gov> - Bug fix and use of
!           List_allocated() function to check for existence of 
!           attributes.
!  5Aug02 - E. Ong <eong@mcs.anl.gov> - more rigorous revision
!EOP ___________________________________________________________________
!
 character(len=*),parameter :: myname_=myname//'::lsize_'

 if(List_allocated(GGrid%data%rList) .and. &
      List_allocated(GGrid%data%iList)) then

    lsize_ = AttrVect_lsize( GGrid%data )

 else

    call die(myname_,"Argument GGrid%data is not associated!")

 endif

 end function lsize_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: exportIAttr_ - Return GeneralGrid INTEGER Attribute as a Vector
!
! !DESCRIPTION:
! This routine extracts from the input {\tt GeneralGrid} argument 
! {\tt GGrid} the integer attribute corresponding to the tag defined in 
! the input {\tt CHARACTER} argument {\tt AttrTag}, and returns it in 
! the {\tt INTEGER} output array {\tt outVect}, and its length in the 
! output {\tt INTEGER} argument {\tt lsize}.
!
! {\bf N.B.:}  This routine will fail if the {\tt AttrTag} is not in 
! the {\tt GeneralGrid} {\tt List} component {\tt GGrid\%data\%iList}.
!
! {\bf N.B.:}  The flexibility of this routine regarding the pointer 
! association status of the output argument {\tt outVect} means the
! user must invoke this routine with care.  If the user wishes this
! routine to fill a pre-allocated array, then obviously this array
! must be allocated prior to calling this routine.  If the user wishes
! that the routine {\em create} the output argument array {\tt outVect},
! then the user must ensure this pointer is not allocated (i.e. the user
! must nullify this pointer) before this routine is invoked.
!
! {\bf N.B.:}  If the user has relied on this routine to allocate memory
! associated with the pointer {\tt outVect}, then the user is responsible 
! for deallocating this array once it is no longer needed.  Failure to 
! do so will result in a memory leak.
!
! !INTERFACE:

! Subprogram not used  subroutine exportIAttr_(GGrid, AttrTag, outVect, lsize)
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
! Subprogram not used       type(GeneralGrid),      intent(in)  :: GGrid
! Subprogram not used       character(len=*),       intent(in)  :: AttrTag
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS: 
! Subprogram not used 
! Subprogram not used       integer,  dimension(:), pointer     :: outVect
! Subprogram not used       integer,  optional,     intent(out) :: lsize
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 13Dec01 - J.W. Larson <larson@mcs.anl.gov> - initial prototype.
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::exportIAttr_'
! Subprogram not used 
! Subprogram not used        ! Export the data (inheritance from AttrVect)
! Subprogram not used   if(present(lsize)) then
! Subprogram not used      call AttrVect_exportIAttr(GGrid%data, AttrTag, outVect, lsize)
! Subprogram not used   else
! Subprogram not used      call AttrVect_exportIAttr(GGrid%data, AttrTag, outVect)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used  end subroutine exportIAttr_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: exportRAttrSP_ - Return GeneralGrid REAL Attribute as a Vector
!
! !DESCRIPTION:
! This routine extracts from the input {\tt GeneralGrid} argument 
! {\tt GGrid} the real attribute corresponding to the tag defined in 
! the input {\tt CHARACTER} argument {\tt AttrTag}, and returns it in 
! the {\tt REAL} output array {\tt outVect}, and its length in the 
! output {\tt INTEGER} argument {\tt lsize}.
!
! {\bf N.B.:}  This routine will fail if the {\tt AttrTag} is not in 
! the {\tt GeneralGrid} {\tt List} component {\tt GGrid\%data\%rList}.
!
! {\bf N.B.:}  The flexibility of this routine regarding the pointer 
! association status of the output argument {\tt outVect} means the
! user must invoke this routine with care.  If the user wishes this
! routine to fill a pre-allocated array, then obviously this array
! must be allocated prior to calling this routine.  If the user wishes
! that the routine {\em create} the output argument array {\tt outVect},
! then the user must ensure this pointer is not allocated (i.e. the user
! must nullify this pointer) before this routine is invoked.
!
! {\bf N.B.:}  If the user has relied on this routine to allocate memory
! associated with the pointer {\tt outVect}, then the user is responsible 
! for deallocating this array once it is no longer needed.  Failure to 
! do so will result in a memory leak.
!
! !INTERFACE:

! Subprogram not used  subroutine exportRAttrSP_(GGrid, AttrTag, outVect, lsize)
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_die
! Subprogram not used       use m_stdio
! Subprogram not used 
! Subprogram not used       use m_realkinds,  only : SP
! Subprogram not used 
! Subprogram not used       use m_AttrVect,   only : AttrVect_exportRAttr => exportRAttr
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS: 
! Subprogram not used 
! Subprogram not used       type(GeneralGrid),          intent(in)  :: GGrid
! Subprogram not used       character(len=*),           intent(in)  :: AttrTag
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS: 
! Subprogram not used 
! Subprogram not used       real(SP),  dimension(:),    pointer     :: outVect
! Subprogram not used       integer,   optional,        intent(out) :: lsize
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 13Dec01 - J.W. Larson <larson@mcs.anl.gov> - initial prototype.
! Subprogram not used !
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::exportRAttrSP_'
! Subprogram not used 
! Subprogram not used        ! Export the data (inheritance from AttrVect)
! Subprogram not used 
! Subprogram not used   if(present(lsize)) then
! Subprogram not used      call AttrVect_exportRAttr(GGrid%data, AttrTag, outVect, lsize)
! Subprogram not used   else
! Subprogram not used      call AttrVect_exportRAttr(GGrid%data, AttrTag, outVect)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used  end subroutine exportRAttrSP_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
! ---------------------------------------------------------------------
!
! !IROUTINE: exportRAttrDP_ - Return GeneralGrid REAL Attribute as a Vector
!
! !DESCRIPTION:
! double precision version of exportRAttrSP_
!
! !INTERFACE:

 subroutine exportRAttrDP_(GGrid, AttrTag, outVect, lsize)
!
! !USES:
!
      use m_die
      use m_stdio

      use m_realkinds,  only : DP

      use m_AttrVect,   only : AttrVect_exportRAttr => exportRAttr

      implicit none

! !INPUT PARAMETERS: 

      type(GeneralGrid),          intent(in)  :: GGrid
      character(len=*),           intent(in)  :: AttrTag

! !OUTPUT PARAMETERS: 

      real(DP),  dimension(:),    pointer     :: outVect
      integer,   optional,        intent(out) :: lsize

! !REVISION HISTORY:
! 13Dec01 - J.W. Larson <larson@mcs.anl.gov> - initial prototype.
!
!_______________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::exportRAttrDP_'

       ! Export the data (inheritance from AttrVect)
  if(present(lsize)) then
     call AttrVect_exportRAttr(GGrid%data, AttrTag, outVect, lsize)
  else
     call AttrVect_exportRAttr(GGrid%data, AttrTag, outVect)
  endif

 end subroutine exportRAttrDP_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: importIAttr_ - Import GeneralGrid INTEGER Attribute
!
! !DESCRIPTION:
! This routine imports data provided in the input {\tt INTEGER} vector 
! {\tt inVect} into the {\tt GeneralGrid} argument {\tt GGrid}, storing 
! it as the integer attribute corresponding to the tag defined in 
! the input {\tt CHARACTER} argument {\tt AttrTag}.  The input 
! {\tt INTEGER} argument {\tt lsize} is used to ensure there is 
! sufficient space in the {\tt GeneralGrid} to store the data.
!
! {\bf N.B.:}  This routine will fail if the {\tt AttrTag} is not in 
! the {\tt GeneralGrid} {\tt List} component {\tt GGrid\%data\%iList}.
!
! !INTERFACE:

 subroutine importIAttr_(GGrid, AttrTag, inVect, lsize)
!
! !USES:
!
      use m_die
      use m_stdio

      use m_AttrVect,      only : AttrVect_importIAttr => importIAttr

      implicit none

! !INPUT PARAMETERS: 

      character(len=*),       intent(in)    :: AttrTag
      integer,  dimension(:), pointer       :: inVect
      integer,                intent(in)    :: lsize

! !INPUT/OUTPUT PARAMETERS: 

      type(GeneralGrid),      intent(inout) :: GGrid

! !REVISION HISTORY:
! 13Dec01 - J.W. Larson <larson@mcs.anl.gov> - initial prototype.
! 27Mar02 - Jay Larson <larson@mcs.anl.gov> - improved error handling.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::importIAttr_'

       ! Argument Check:

  if(lsize > lsize_(GGrid)) then
     write(stderr,*) myname_,':: ERROR, lsize > lsize_(GGrid).', &
	  'lsize = ',lsize,'lsize_(GGrid) = ',lsize_(GGrid)
     call die(myname_)
  endif

       ! Import the data (inheritance from AttrVect)

  call AttrVect_importIAttr(GGrid%data, AttrTag, inVect, lsize)

 end subroutine importIAttr_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: importRAttrSP_ - Import GeneralGrid REAL Attribute
!
! !DESCRIPTION:
! This routine imports data provided in the input {\tt REAL} vector 
! {\tt inVect} into the {\tt GeneralGrid} argument {\tt GGrid}, storing 
! it as the real attribute corresponding to the tag defined in 
! the input {\tt CHARACTER} argument {\tt AttrTag}.  The input 
! {\tt INTEGER} argument {\tt lsize} is used to ensure there is 
! sufficient space in the {\tt GeneralGrid} to store the data.
!
! {\bf N.B.:}  This routine will fail if the {\tt AttrTag} is not in 
! the {\tt GeneralGrid} {\tt List} component {\tt GGrid\%data\%rList}.
!
! !INTERFACE:

! Subprogram not used  subroutine importRAttrSP_(GGrid, AttrTag, inVect, lsize)
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_die ,          only : die
! Subprogram not used       use m_die ,          only : MP_perr_die
! Subprogram not used       use m_stdio ,        only : stderr
! Subprogram not used 
! Subprogram not used       use m_realkinds,     only : SP
! Subprogram not used 
! Subprogram not used       use m_AttrVect,      only : AttrVect_importRAttr => importRAttr
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS: 
! Subprogram not used 
! Subprogram not used       character(len=*),           intent(in)    :: AttrTag
! Subprogram not used       real(SP), dimension(:),     pointer       :: inVect
! Subprogram not used       integer,                    intent(in)    :: lsize
! Subprogram not used 
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS: 
! Subprogram not used 
! Subprogram not used       type(GeneralGrid),          intent(inout) :: GGrid
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 13Dec01 - J.W. Larson <larson@mcs.anl.gov> - initial prototype.
! Subprogram not used ! 27Mar02 - Jay Larson <larson@mcs.anl.gov> - improved error handling.
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::importRAttrSP_'
! Subprogram not used 
! Subprogram not used        ! Argument Check:
! Subprogram not used 
! Subprogram not used   if(lsize > lsize_(GGrid)) then
! Subprogram not used      write(stderr,*) myname_,':: ERROR, lsize > lsize_(GGrid).', &
! Subprogram not used 	  'lsize = ',lsize,'lsize_(GGrid) = ',lsize_(GGrid)
! Subprogram not used      call die(myname_)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used        ! Import the data (inheritance from AttrVect)
! Subprogram not used 
! Subprogram not used   call AttrVect_importRAttr(GGrid%data, AttrTag, inVect, lsize)
! Subprogram not used 
! Subprogram not used  end subroutine importRAttrSP_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!-----------------------------------------------------------------------
!
! !IROUTINE: importRAttrDP_ - Import GeneralGrid REAL Attribute
!
! !DESCRIPTION:
! Double precision version of importRAttrSP_
!
! !INTERFACE:

 subroutine importRAttrDP_(GGrid, AttrTag, inVect, lsize)
!
! !USES:
!
      use m_die ,          only : die
      use m_die ,          only : MP_perr_die
      use m_stdio ,        only : stderr

      use m_realkinds,     only : DP

      use m_AttrVect,      only : AttrVect_importRAttr => importRAttr

      implicit none

! !INPUT PARAMETERS: 

      character(len=*),           intent(in)    :: AttrTag
      real(DP), dimension(:),     pointer       :: inVect
      integer,                    intent(in)    :: lsize

! !INPUT/OUTPUT PARAMETERS: 

      type(GeneralGrid),          intent(inout) :: GGrid

! !REVISION HISTORY:
! 13Dec01 - J.W. Larson <larson@mcs.anl.gov> - initial prototype.
! 27Mar02 - Jay Larson <larson@mcs.anl.gov> - improved error handling.
!_______________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::importRAttrDP_'

       ! Argument Check:

  if(lsize > lsize_(GGrid)) then
     write(stderr,*) myname_,':: ERROR, lsize > lsize_(GGrid).', &
	  'lsize = ',lsize,'lsize_(GGrid) = ',lsize_(GGrid)
     call die(myname_)
  endif

       ! Import the data (inheritance from AttrVect)

  call AttrVect_importRAttr(GGrid%data, AttrTag, inVect, lsize)

 end subroutine importRAttrDP_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: Sort_ - Generate Sort Permutation Defined by Arbitrary Keys.
!
! !DESCRIPTION:
! The subroutine {\tt Sort\_()} uses the list of keys present in the
! input  {\tt List} variable {\tt key\_List}.  This list of keys is 
! checked to ensure that {\em only} coordinate attributes are present 
! in the sorting keys, and that there are no redundant keys.  Once 
! checked, this list is used to find the appropriate real attributes
! referenced by the items in {\tt key\_list} ( that is, it identifies the 
! appropriate entries in {\tt GGrid\%data\%rList}), and then uses these 
! keys to generate a an output permutation {\tt perm} that will put
! the entries of the attribute vector {\tt GGrid\%data} in lexicographic 
! order as defined by {\tt key\_list} (the ordering in {\tt key\_list} 
! being from left to right.
!
! !INTERFACE:

! Subprogram not used  subroutine Sort_(GGrid, key_List, perm, descend)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_stdio
! Subprogram not used       use m_die
! Subprogram not used 
! Subprogram not used       use m_AttrVect,     only : AttrVect_Sort => Sort
! Subprogram not used       use m_List,        only : List_nitem => nitem
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS: 
! Subprogram not used !
! Subprogram not used       type(GeneralGrid),               intent(in) :: GGrid
! Subprogram not used       type(List),                      intent(in) :: key_list
! Subprogram not used       logical, dimension(:), optional, intent(in) :: descend
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS: 
! Subprogram not used !
! Subprogram not used       integer, dimension(:), pointer              :: perm
! Subprogram not used 
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 15Jan01 - Jay Larson <larson@mcs.anl.gov> - Initial version.
! Subprogram not used ! 20Mar01 - Jay Larson <larson@mcs.anl.gov> - Final working version.
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used !
! Subprogram not used  character(len=*),parameter :: myname_=myname//'::Sort_'
! Subprogram not used  logical, dimension(:), allocatable :: descending
! Subprogram not used  integer :: n, ierr
! Subprogram not used 
! Subprogram not used         ! Here is how we transmit the sort order keys stored
! Subprogram not used         ! in descending (if present):
! Subprogram not used 
! Subprogram not used   n = List_nitem(key_list)
! Subprogram not used   allocate(descending(n), stat=ierr)
! Subprogram not used   if(ierr /= 0) then
! Subprogram not used      call die(myname_,"allocate(descending...",ierr)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   if(present(descend)) then
! Subprogram not used      descending = descend
! Subprogram not used   else
! Subprogram not used      descending = .false.
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used         ! This is a straightforward call to AttrVect_Sort().
! Subprogram not used   
! Subprogram not used   call AttrVect_Sort(GGrid%data, key_list, perm, descending) 
! Subprogram not used 
! Subprogram not used         ! Clean up...
! Subprogram not used 
! Subprogram not used   deallocate(descending, stat=ierr)
! Subprogram not used   if(ierr /= 0) then
! Subprogram not used      call die(myname_,"deallocate(descending...",ierr)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used  end subroutine Sort_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: Sortg_ - Generate Sort Permutation Based on GeneralGrid Keys.
!
! !DESCRIPTION:
! The subroutine {\tt Sortg\_()} uses the list of sorting keys present in 
! the input  {\tt GeneralGrid} variable {\tt GGrid\%coordinate\_sort\_order}
! to create a sort permutation {\tt perm(:)}.  Sorting is either in ascending
! or descending order based on the entries of {\tt GGrid\%descend(:)}.
! The output index permutation is stored in the array {\tt perm(:)} that 
! will put the entries of the attribute vector {\tt GGrid\%data} in 
! lexicographic order as defined by {\tt GGrid\%coordinate\_sort\_order}. The 
! ordering in {\tt GGrid\%coordinate\_sort\_order} being from left to right.
!
! {\bf N.B.:}  This routine returnss an allocatable array perm(:).  This 
! allocated array must be deallocated when the user no longer needs it.  
! Failure to do so will cause a memory leak.
!
! {\bf N.B.:}  This routine will fail if {\tt GGrid} has not been initialized 
! with sort keys in the {\tt List} component {\tt GGrid\%coordinate\_sort\_order}.
!
! !INTERFACE:

! Subprogram not used  subroutine Sortg_(GGrid, perm)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_List, only : List_allocated => allocated
! Subprogram not used       use m_die,  only : die
! Subprogram not used    
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS: 
! Subprogram not used !
! Subprogram not used       type(GeneralGrid),     intent(in) :: GGrid
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS: 
! Subprogram not used !
! Subprogram not used       integer, dimension(:), pointer    :: perm
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 22Mar01 - Jay Larson <larson@mcs.anl.gov> - Initial version.
! Subprogram not used !  5Aug02 - E. Ong <eong@mcs.anl.gov> - revise with more error checking.
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used !
! Subprogram not used  character(len=*),parameter :: myname_=myname//'::Sortg_'
! Subprogram not used 
! Subprogram not used      if(.not.List_allocated(GGrid%coordinate_sort_order)) then
! Subprogram not used 	call die(myname_, "GGrid%coordinate_aort_order must be &
! Subprogram not used 	         &allocated for use in any sort function")
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used      if(associated(GGrid%descend)) then
! Subprogram not used 	call Sort_(GGrid, GGrid%coordinate_sort_order, &
! Subprogram not used 	               perm, GGrid%descend)
! Subprogram not used      else
! Subprogram not used 	call Sort_(GGrid=GGrid, key_list=GGrid%coordinate_sort_order, &
! Subprogram not used 	               perm=perm)
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used  end subroutine Sortg_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: Permute_ - Permute GeneralGrid Attributes Using Supplied Index Permutation
!
! !DESCRIPTION:
! The subroutine {\tt Permute\_()} uses an input index permutation {\tt perm} 
! to re-order the coordinate data stored in the {\tt GeneralGrid} argument 
! {\tt GGrid}.  This permutation can be generated by either of the routines
! {\tt Sort\_()} or {\tt Sortg\_()} contained in this module.
!
! !INTERFACE:

! Subprogram not used  subroutine Permute_(GGrid, perm)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used 
! Subprogram not used       use m_stdio
! Subprogram not used       use m_die
! Subprogram not used 
! Subprogram not used       use m_AttrVect,     only : AttrVect
! Subprogram not used       use m_AttrVect, only : AttrVect_Permute => Permute
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS: 
! Subprogram not used !
! Subprogram not used       integer, dimension(:), intent(in)    :: perm
! Subprogram not used 
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS: 
! Subprogram not used !
! Subprogram not used       type(GeneralGrid),     intent(inout) :: GGrid
! Subprogram not used 
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 15Jan01 - Jay Larson <larson@mcs.anl.gov> - API specification.
! Subprogram not used ! 10Apr01 - Jay Larson <larson@mcs.anl.gov> - API modified, working
! Subprogram not used !           code.
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used !
! Subprogram not used  character(len=*),parameter :: myname_=myname//'::Permute_'
! Subprogram not used 
! Subprogram not used         ! This is a straightforward call to AttrVect_Permute:
! Subprogram not used 
! Subprogram not used   call AttrVect_Permute(GGrid%data, perm)
! Subprogram not used 
! Subprogram not used  end subroutine Permute_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: SortPermute_ - Sort and Permute GeneralGrid Attributes
!
! !DESCRIPTION:
! The subroutine {\tt SortPermute\_()} uses the list of keys defined in 
! {\tt GGrid\%coordinate\_sort\_order} to create an index permutation 
! {\tt perm}, which is then applied to re-order the coordinate data stored 
! in the {\tt GeneralGrid} argument {\tt GGrid} (more specifically, the 
! gridpoint data stored in {\tt GGrid\%data}.  This permutation is generated  
! by the routine {\tt Sortg\_()} contained in this module.  The permutation 
! is carried out by the routine {\tt Permute\_()} contained in this module.
!
! {\bf N.B.:}  This routine will fail if {\tt GGrid} has not been initialized 
! with sort keys in the {\tt List} component {\tt GGrid\%coordinate\_sort\_order}.
!
! !INTERFACE:

! Subprogram not used  subroutine SortPermute_(GGrid)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used       use m_stdio
! Subprogram not used       use m_die
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS: 
! Subprogram not used !
! Subprogram not used       type(GeneralGrid),     intent(inout)   :: GGrid
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 15Jan01 - Jay Larson <larson@mcs.anl.gov> - API specification.
! Subprogram not used ! 10Apr01 - Jay Larson <larson@mcs.anl.gov> - API modified, working
! Subprogram not used !           code.
! Subprogram not used ! 13Apr01 - Jay Larson <larson@mcs.anl.gov> - Simplified API and
! Subprogram not used !           code (Thanks to Tony Craig of NCAR for detecting the
! Subprogram not used !           bug that inspired these changes).
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used !
! Subprogram not used  character(len=*),parameter :: myname_=myname//'::SortPermute_'
! Subprogram not used 
! Subprogram not used  integer, dimension(:), pointer :: perm
! Subprogram not used  integer :: ierr
! Subprogram not used 
! Subprogram not used   call Sortg_(GGrid, perm)
! Subprogram not used 
! Subprogram not used   call Permute_(GGrid, perm)
! Subprogram not used 
! Subprogram not used ! Clean up--deallocate temporary permutation array:
! Subprogram not used 
! Subprogram not used   deallocate(perm, stat=ierr)
! Subprogram not used   if(ierr /= 0) then
! Subprogram not used      call die(myname_,"deallocate(perm)",ierr)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used  end subroutine SortPermute_

 end module m_GeneralGrid
















































