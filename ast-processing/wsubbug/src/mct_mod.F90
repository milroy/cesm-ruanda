!===============================================================================
! SVN $Id: mct_mod.F90 61510 2014-06-26 21:58:56Z tcraig $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_140626/shr/mct_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: mct_mod -- provides a standard API naming convention for MCT code
!
! !DESCRIPTION:
!    This module should be used instead of accessing mct modules directly.  
!    This module:
!    \begin{itemize}
!    \item Uses Fortran {\sf use} renaming of MCT routines and data types so that they
!          all have an mct\_ prefix and related data types and routines have related names. 
!    \item Provides easy and uniform access to 
!          all MCT routines and data types that must be accessed.
!    \item Provides a convienient list of 
!          all MCT routines and data types that can be accessed.
!    \item Blocks access to MCT routines that are not used in cpl6.
!    \end{itemize}
!    This module also includes some MCT-only functions to augment
!    the MCT library.
!
! !REVISION HISTORY:
!     2001-Aug-14 - B. Kauffman - first prototype
!     2006-Apr-13 - M. Vertenstein - modified for sequential mode
!     2007-Mar-01 - R. Jacob - moved to shr
!
! !INTERFACE: ------------------------------------------------------------------
module mct_mod

! !USES:

   use shr_kind_mod         ! shared kinds
   use shr_sys_mod          ! share system routines
   use shr_mpi_mod          ! mpi layer
   use shr_const_mod        ! constants
   use shr_string_mod       ! string functions

   use shr_log_mod          ,only: s_loglev               => shr_log_Level
   use shr_log_mod          ,only: s_logunit              => shr_log_Unit

   use m_MCTWorld           ,only: mct_world_init         => init

   use m_AttrVect           ,only: mct_aVect              => AttrVect
   use m_AttrVect           ,only: mct_aVect_init         => init
   use m_AttrVect           ,only: mct_aVect_clean        => clean
   use m_AttrVect           ,only: mct_aVect_zero         => zero
   use m_AttrVect           ,only: mct_aVect_lsize        => lsize
   use m_AttrVect           ,only: mct_aVect_indexIA      => indexIA
   use m_AttrVect           ,only: mct_aVect_indexRA      => indexRA
   use m_AttrVect           ,only: mct_aVect_importIattr  => importIattr
   use m_AttrVect           ,only: mct_aVect_exportIattr  => exportIattr
   use m_AttrVect           ,only: mct_aVect_importRattr  => importRattr
   use m_AttrVect           ,only: mct_aVect_exportRattr  => exportRattr
   use m_AttrVect           ,only: mct_aVect_getIList     => getIList
   use m_AttrVect           ,only: mct_aVect_getRList     => getRList
   use m_AttrVect           ,only: mct_aVect_getIList2c   => getIListToChar
   use m_AttrVect           ,only: mct_aVect_getRList2c   => getRListToChar
   use m_AttrVect           ,only: mct_aVect_exportIList2c=> exportIListToChar
   use m_AttrVect           ,only: mct_aVect_exportRList2c=> exportRListToChar
   use m_AttrVect           ,only: mct_aVect_nIAttr       => nIAttr
   use m_AttrVect           ,only: mct_aVect_nRAttr       => nRAttr
   use m_AttrVect           ,only: mct_aVect_copy         => Copy
   use m_AttrVect           ,only: mct_aVect_permute      => Permute
   use m_AttrVect           ,only: mct_aVect_unpermute    => Unpermute
   use m_AttrVect           ,only: mct_aVect_SharedIndices=> AVSharedIndices
   use m_AttrVect           ,only: mct_aVect_setSharedIndices=> SharedIndices
   use m_AttrVectComms      ,only: mct_aVect_scatter      => scatter
   use m_AttrVectComms      ,only: mct_aVect_gather       => gather 
   use m_AttrVectComms      ,only: mct_aVect_bcast        => bcast  

   use m_GeneralGrid        ,only: mct_gGrid              => GeneralGrid
   use m_GeneralGrid        ,only: mct_gGrid_init         => init
   use m_GeneralGrid        ,only: mct_gGrid_clean        => clean
   use m_GeneralGrid        ,only: mct_gGrid_dims         => dims
   use m_GeneralGrid        ,only: mct_gGrid_lsize        => lsize
   use m_GeneralGrid        ,only: mct_ggrid_indexIA      => indexIA
   use m_GeneralGrid        ,only: mct_gGrid_indexRA      => indexRA
   use m_GeneralGrid        ,only: mct_gGrid_exportRattr  => exportRattr
   use m_GeneralGrid        ,only: mct_gGrid_importRattr  => importRattr
   use m_GeneralGrid        ,only: mct_gGrid_exportIattr  => exportIattr
   use m_GeneralGrid        ,only: mct_gGrid_importIattr  => importIattr
   use m_GeneralGrid        ,only: mct_gGrid_permute      => permute
   use m_GeneralGridComms   ,only: mct_gGrid_scatter      => scatter
   use m_GeneralGridComms   ,only: mct_gGrid_gather       => gather 
   use m_GeneralGridComms   ,only: mct_gGrid_bcast        => bcast  

   use m_Transfer           ,only: mct_send               => Send
   use m_Transfer           ,only: mct_recv               => Recv
      
   use m_GlobalSegMap       ,only: mct_gsMap              => GlobalSegMap
   use m_GlobalSegMap       ,only: mct_gsMap_init         => init
   use m_GlobalSegMap       ,only: mct_gsMap_clean        => clean
   use m_GlobalSegMap       ,only: mct_gsMap_lsize        => lsize
   use m_GlobalSegMap       ,only: mct_gsMap_gsize        => gsize
   use m_GlobalSegMap       ,only: mct_gsMap_gstorage     => GlobalStorage
   use m_GlobalSegMap       ,only: mct_gsMap_ngseg        => ngseg
   use m_GlobalSegMap       ,only: mct_gsMap_nlseg        => nlseg
   use m_GlobalSegMap       ,only: mct_gsMap_OP           => OrderedPoints
   use m_GlobalSegMap       ,only: mct_gsMap_maxnlseg     => max_nlseg
   use m_GlobalSegMap       ,only: mct_gsMap_activepes    => active_pes
   use m_GlobalSegMap       ,only: mct_gsMap_copy         => copy
   use m_GlobalSegMap       ,only: mct_gsMap_increasing   => increasing
   use m_GlobalSegMap       ,only: mct_gsMap_orderedPoints=> OrderedPoints
   use m_GlobalSegMapComms  ,only: mct_gsMap_bcast        => bcast 

   use m_Rearranger         ,only: mct_rearr              => Rearranger
   use m_Rearranger         ,only: mct_rearr_init         => init
   use m_Rearranger         ,only: mct_rearr_clean        => clean
   use m_Rearranger         ,only: mct_rearr_print        => print
   use m_Rearranger         ,only: mct_rearr_rearrange    => rearrange

   use m_Router             ,only: mct_router             => Router
   use m_Router             ,only: mct_router_init        => init

   use m_SparseMatrixToMaps ,only: mct_sMat_2XgsMap       => SparseMatrixToXGlobalSegMap
   use m_SparseMatrixToMaps ,only: mct_sMat_2YgsMap       => SparseMatrixToYGlobalSegMap
   use m_SparseMatrix       ,only: mct_sMat               => SparseMatrix
   use m_SparseMatrix       ,only: mct_sMat_Init          => init
   use m_SparseMatrix       ,only: mct_sMat_Vecinit       => vecinit
   use m_SparseMatrix       ,only: mct_sMat_Clean         => clean
   use m_SparseMatrix       ,only: mct_sMat_indexIA       => indexIA
   use m_SparseMatrix       ,only: mct_sMat_indexRA       => indexRA
   use m_SparseMatrix       ,only: mct_sMat_lsize         => lsize
   use m_SparseMatrix       ,only: mct_sMat_nrows         => nRows
   use m_SparseMatrix       ,only: mct_sMat_ncols         => nCols
   use m_SparseMatrix       ,only: mct_sMat_SortPermute   => SortPermute
   use m_SparseMatrix       ,only: mct_sMat_GNumEl        => GlobalNumElements
   use m_SparseMatrix       ,only: mct_sMat_ImpGRowI      => ImportGlobalRowIndices
   use m_SparseMatrix       ,only: mct_sMat_ImpGColI      => ImportGlobalColumnIndices
   use m_SparseMatrix       ,only: mct_sMat_ImpLRowI      => ImportLocalRowIndices
   use m_SparseMatrix       ,only: mct_sMat_ImpLColI      => ImportLocalColumnIndices
   use m_SparseMatrix       ,only: mct_sMat_ImpMatrix     => ImportMatrixElements
   use m_SparseMatrix       ,only: mct_sMat_ExpGRowI      => ExportGlobalRowIndices
   use m_SparseMatrix       ,only: mct_sMat_ExpGColI      => ExportGlobalColumnIndices
   use m_SparseMatrix       ,only: mct_sMat_ExpLRowI      => ExportLocalRowIndices
   use m_SparseMatrix       ,only: mct_sMat_ExpLColI      => ExportLocalColumnIndices
   use m_SparseMatrix       ,only: mct_sMat_ExpMatrix     => ExportMatrixElements
   use m_SparseMatrixComms  ,only: mct_sMat_ScatterByRow  => ScatterByRow
   use m_SparseMatrixComms  ,only: mct_sMat_ScatterByCol  => ScatterByColumn
   use m_SparseMatrixPlus   ,only: mct_sMatP              => SparseMatrixPlus
   use m_SparseMatrixPlus   ,only: mct_sMatP_Init         => init
   use m_SparseMatrixPlus   ,only: mct_sMatP_Vecinit      => vecinit
   use m_SparseMatrixPlus   ,only: mct_sMatP_clean        => clean
   use m_MatAttrVectMul     ,only: mct_sMat_avMult        => sMatAvMult
   use m_GlobalToLocal      ,only: mct_sMat_g2lMat        => GlobalToLocalMatrix

   use m_List               ,only: mct_list               => list     
   use m_List               ,only: mct_list_init          => init
   use m_List               ,only: mct_list_get           => get 
   use m_List               ,only: mct_list_nitem         => nitem 
   use m_List               ,only: mct_list_clean         => clean
   use m_string             ,only: mct_string             => string 
   use m_string             ,only: mct_string_clean       => clean
   use m_string             ,only: mct_string_toChar      => toChar 
   use m_die                ,only: mct_perr_die           => mp_perr_die
   use m_die                ,only: mct_die                => die
   use m_inpak90

   use m_Permuter           ,only: mct_permute            => Permute

   use m_MergeSorts         ,only: mct_indexset           => IndexSet
   use m_MergeSorts         ,only: mct_indexsort          => IndexSort

   implicit none

   public :: mct_aVect_info
   public :: mct_aVect_fldIndex
   public :: mct_aVect_sharedFields
   public :: mct_aVect_initSharedFields
   public :: mct_aVect_getRAttr
   public :: mct_aVect_putRAttr
   public :: mct_aVect_accum
   public :: mct_aVect_avg
   public :: mct_avect_mult
   public :: mct_avect_vecmult
   public :: mct_rearr_rearrange_fldlist
   public :: mct_gsmap_identical

   logical,public :: mct_usealltoall = .false.
   logical,public :: mct_usevector = .false.

!EOP

   !--- local kinds ---
   integer,parameter,private :: R8 = SHR_KIND_R8
   integer,parameter,private :: IN = SHR_KIND_IN
   integer,parameter,private :: CL = SHR_KIND_CL
   integer,parameter,private :: CX = SHR_KIND_CX
   integer,parameter,private :: CXX = SHR_KIND_CXX

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: mct_aVect_info - print out aVect info for debugging
!
! !DESCRIPTION:
!     Print out information about the input MCT {\it AttributeVector}
!     {\tt aVect} to stdout. {\tt flag} sets the level of information:
!     \begin{enumerate}
!     \item  print out names of attributes in {\tt aVect}.
!     \item  also print out local max and min of data in {\tt aVect}.
!     \item  also print out global max and min of data in {\tt aVect}.
!     \item  Same as 3 but include name of this routine.
!     \end{enumerate}
!     If {\tt flag} is 3 or higher, then optional argument {\tt comm}
!     must be provided.
!     If optional argument {\tt fld} is present, only information for
!     that field will be printed.
!     If optional argument {\tt istr} is present, it will be output
!     before any of the information.
!
!
! !REVISION HISTORY:
!     2003 Jul 01 - B. Kauffman, T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used subroutine mct_aVect_info(flag,aVect,comm,pe,fld,istr)
! Subprogram not used 
! Subprogram not used ! !USES:
! Subprogram not used   
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    integer(IN)    ,intent(in)           :: flag  ! info level flag
! Subprogram not used    type(mct_aVect),intent(in)           :: aVect ! Attribute vector
! Subprogram not used    integer(IN)    ,intent(in),optional  :: comm  ! MPI communicator
! Subprogram not used    integer(IN)    ,intent(in),optional  :: pe    ! processor number
! Subprogram not used    character(*)   ,intent(in),optional  :: fld   ! fld
! Subprogram not used    character(*)   ,intent(in),optional  :: istr  ! string for print
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used 
! Subprogram not used    !--- local ---
! Subprogram not used    integer(IN)          :: i,j,k,n      ! generic indicies
! Subprogram not used    integer(IN)          :: ks,ke        ! start and stop k indices
! Subprogram not used    integer(IN)          :: nflds        ! number of flds in AV to diagnose
! Subprogram not used    integer(IN)          :: nsize        ! grid point size of AV
! Subprogram not used    type(mct_string)     :: item         ! mct string
! Subprogram not used    character(CL)        :: itemc        ! item converted to char
! Subprogram not used    integer(IN)          :: comm_loc     ! local variable for comm
! Subprogram not used    integer(IN)          :: pe_loc       ! local variable for pe
! Subprogram not used    logical              :: commOK       ! is comm available
! Subprogram not used    logical              :: peOK         ! is pe available
! Subprogram not used    real(R8),allocatable :: minl(:)      ! local  min
! Subprogram not used    real(R8),allocatable :: ming(:)      ! global min
! Subprogram not used    real(R8),allocatable :: maxl(:)      ! local  max
! Subprogram not used    real(R8),allocatable :: maxg(:)      ! global max
! Subprogram not used 
! Subprogram not used    !--- formats ---
! Subprogram not used    character(*),parameter :: subName = '(mct_aVect_info) '
! Subprogram not used    character(*),parameter :: F00 = "('(mct_aVect_info) ',8a)"
! Subprogram not used    character(*),parameter :: F01 = "('(mct_aVect_info) ',a,i9)"
! Subprogram not used    character(*),parameter :: F02 = "('(mct_aVect_info) ',240a)"
! Subprogram not used    character(*),parameter :: F03 = "('(mct_aVect_info) ',a,2es11.3,i4,2x,a)"
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used ! NOTE: has hard-coded knowledge/assumptions about mct aVect data type internals
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    commOK = .false.
! Subprogram not used    peOK   = .false.
! Subprogram not used 
! Subprogram not used    if (present(pe)) then
! Subprogram not used      peOK = .true.
! Subprogram not used      pe_loc = pe
! Subprogram not used    endif
! Subprogram not used    if (present(comm)) then
! Subprogram not used      commOK = .true.
! Subprogram not used      comm_loc = comm
! Subprogram not used      if (.not.PEOK) then
! Subprogram not used        call shr_mpi_commrank(comm,pe_loc,subName)
! Subprogram not used        peOK = .true.
! Subprogram not used      endif
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    nsize = mct_aVect_lsize(aVect)
! Subprogram not used 
! Subprogram not used    if (present(fld)) then
! Subprogram not used      nflds = 1
! Subprogram not used      ks = mct_aVect_indexRA(aVect,fld,perrWith=subName)
! Subprogram not used      ke = ks
! Subprogram not used    else
! Subprogram not used      nflds = mct_aVect_nRAttr(aVect)
! Subprogram not used      ks = 1
! Subprogram not used      ke = nflds
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    if (flag >= 1) then
! Subprogram not used      if (present(istr)) then
! Subprogram not used         if (s_loglev > 0) write(s_logunit,*) trim(istr)
! Subprogram not used      endif
! Subprogram not used      if (s_loglev > 0) write(s_logunit,F01) "local size =",nsize
! Subprogram not used      if (associated(aVect%iList%bf)) then
! Subprogram not used         if (s_loglev > 0) write(s_logunit,F02) "iList = ",aVect%iList%bf
! Subprogram not used      endif
! Subprogram not used      if (associated(aVect%rList%bf)) then
! Subprogram not used         if (s_loglev > 0) write(s_logunit,F02) "rList = ",aVect%rList%bf
! Subprogram not used      endif
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    if (flag >= 2) then
! Subprogram not used 
! Subprogram not used      allocate(minl(nflds))
! Subprogram not used      allocate(maxl(nflds))
! Subprogram not used 
! Subprogram not used      do k=ks,ke
! Subprogram not used        minl(k) = minval(aVect%rAttr(k,:))
! Subprogram not used        maxl(k) = maxval(aVect%rAttr(k,:))
! Subprogram not used      enddo
! Subprogram not used 
! Subprogram not used      if (flag >= 4 .and. commOK) then
! Subprogram not used        allocate(ming(nflds))
! Subprogram not used        allocate(maxg(nflds))
! Subprogram not used        ming = 0._R8
! Subprogram not used        maxg = 0._R8
! Subprogram not used        call shr_mpi_min(minl,ming,comm,subName)
! Subprogram not used        call shr_mpi_max(maxl,maxg,comm,subName)
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used      do k=ks,ke
! Subprogram not used        call mct_aVect_getRList(item,k,aVect)
! Subprogram not used        itemc = mct_string_toChar(item)
! Subprogram not used        call mct_string_clean(item)
! Subprogram not used        if (s_loglev > 0) write(s_logunit,F03) 'l min/max ',minl(k),maxl(k),k,trim(itemc)
! Subprogram not used        if (flag >= 3 .and. commOK) then
! Subprogram not used          if ((peOK .and. pe_loc == 0) .or. .not.peOK) then
! Subprogram not used            if (s_loglev > 0) write(s_logunit,F03) 'g min/max ',ming(k),maxg(k),k,trim(itemc)
! Subprogram not used          endif
! Subprogram not used        endif
! Subprogram not used        if (flag >= 4 .and. commOK) then
! Subprogram not used          if ((peOK .and. pe_loc == 0) .or. .not.peOK) then
! Subprogram not used            if (s_loglev > 0) write(s_logunit,*) trim(subName),'g min/max ',ming(k),maxg(k),k,trim(itemc)
! Subprogram not used          endif
! Subprogram not used        endif
! Subprogram not used      enddo
! Subprogram not used 
! Subprogram not used       deallocate(minl)
! Subprogram not used       deallocate(maxl)
! Subprogram not used       if (flag >= 4 .and. commOK) then
! Subprogram not used          deallocate(ming)
! Subprogram not used          deallocate(maxg)
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    call shr_sys_flush(s_logunit)
! Subprogram not used 
! Subprogram not used end subroutine mct_aVect_info

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: mct_aVect_fldIndex - get a real fld index from an AVect
!
! !DESCRIPTION:
!     Get the field index for a real field in an attribute vector.
!     This is like mct_aVect_indexRA but with a calling interface
!     that returns the index without any error messages.
!
! !REMARKS:
!   This is like the MCT routine indexRA
!
! !REVISION HISTORY:
!     2010 Oct 27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used integer function mct_aVect_fldIndex(aVect,fld)
! Subprogram not used 
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    type(mct_aVect),intent(in)  :: aVect    ! an Attribute vector
! Subprogram not used    character(*)   ,intent(in)  :: fld      ! field name string
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used 
! Subprogram not used    !--- local ---
! Subprogram not used 
! Subprogram not used    !--- formats ---
! Subprogram not used    character(*),parameter :: subName = "(mct_aVect_fldIndex) "
! Subprogram not used    character(*),parameter :: F00 = "('(mct_aVect_fldIndex) ',8a)"
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    mct_aVect_fldIndex = mct_aVect_indexRA(aVect,trim(fld),perrWith='quiet')
! Subprogram not used 
! Subprogram not used end function mct_aVect_fldIndex

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: mct_aVect_sharedFields - get a shared real fld index from two AVects
!
! !DESCRIPTION:
!     Get the shared field index for a real field in two attribute vectors.
!
! !REMARKS:
!
! !REVISION HISTORY:
!     2013 Jul 17 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine mct_aVect_sharedFields(aVect1, aVect2, rlistout, ilistout)

! !INPUT/OUTPUT PARAMETERS:

   type(mct_aVect),intent(in)  :: aVect1   ! an Attribute vector
   type(mct_aVect),intent(in)  :: aVect2   ! an Attribute vector
   character(*)   ,intent(inout),optional  :: rlistout      ! field name string
   character(*)   ,intent(inout),optional  :: ilistout      ! field name string

!EOP

   !--- local ---
   integer(IN) :: nflds1,nflds2
   character(len=CXX) :: list1,list2

   !--- formats ---
   character(*),parameter :: subName = "(mct_aVect_sharedFields) "
   character(*),parameter :: F00 = "('(mct_aVect_sharedFields) ',8a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (present(rlistout)) then
      nflds1 = mct_aVect_nRAttr(aVect1)
      nflds2 = mct_aVect_nRAttr(aVect2)
      rlistout = ''
      list1 = ''
      list2 = ''
      if (nflds1 > 0 .and. nflds2 > 0) then
         list1 = mct_aVect_exportRList2c(aVect1)
         list2 = mct_aVect_exportRlist2c(aVect2)
         call shr_string_listIntersect(list1,list2,rlistout)
      endif
   endif

   if (present(ilistout)) then
      nflds1 = mct_aVect_nIAttr(aVect1)
      nflds2 = mct_aVect_nIAttr(aVect2)
      ilistout = ''
      list1 = ''
      list2 = ''
      if (nflds1 > 0 .and. nflds2 > 0) then
         list1 = mct_aVect_exportIList2c(aVect1)
         list2 = mct_aVect_exportIlist2c(aVect2)
         call shr_string_listIntersect(list1,list2,ilistout)
      endif
   endif

end subroutine mct_aVect_sharedFields

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: mct_aVect_initSharedFields - init new AVect based on shared fields 
!     from two input aVects
!
! !DESCRIPTION:
!     Init new AVect based on shared fields of two input AVects
!
! !REMARKS:
!
! !REVISION HISTORY:
!     2013 Jul 17 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine mct_aVect_initSharedFields(aVect1, aVect2, aVect3, lsize)

! !INPUT/OUTPUT PARAMETERS:

   type(mct_aVect),intent(in)  :: aVect1   ! an Attribute vector
   type(mct_aVect),intent(in)  :: aVect2   ! an Attribute vector
   type(mct_aVect),intent(inout)  :: aVect3   ! new Attribute vector
   integer(IN)    ,intent(in)  :: lsize    ! aVect3 size

!EOP

   !--- local ---
   character(len=CXX) :: rlist,ilist

   !--- formats ---
   character(*),parameter :: subName = "(mct_aVect_initSharedFields) "
   character(*),parameter :: F00 = "('(mct_aVect_initSharedFields) ',8a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   call mct_aVect_sharedFields(aVect1,aVect2,rlist,ilist)
   call mct_aVect_init(aVect3,ilist,rlist,lsize)

end subroutine mct_aVect_initSharedFields

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: mct_aVect_getRAttr - get real F90 array data out of an aVect
!
! !DESCRIPTION:
!     Get the data associated with attribute {\tt str} in 
!     {\it AttributeVector} {\tt aVect} and return in the
!     real F90 array data {\tt data}.
!     {\tt rcode} will be 0 if succesful, 1 if size of {\tt data}
!     does not match size  of {\tt aVect} and 2 if {\tt str} is
!     not found.
!
! !REMARKS:
!   This is like the MCT routine exportRAttr except the output argument
!   is not a pointer.
!
! !REVISION HISTORY:
!     2002 Apr xx - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine mct_aVect_getRAttr(aVect,str,data,rcode)

! !INPUT/OUTPUT PARAMETERS:

   type(mct_aVect)    ,intent(in)  :: aVect    ! an Attribute vector
   character(*)       ,intent(in)  :: str      ! field name string
   real(R8)           ,intent(out) :: data(:)  ! an F90 array
   integer(IN)        ,intent(out) :: rcode    ! return code

!EOP

   !--- local ---
   integer(IN) :: k,n,m
   integer(IN) :: aVsize

   !--- formats ---
   character(*),parameter :: subName = "(mct_aVect_getRAttr) "
   character(*),parameter :: F00 = "('(mct_aVect_getRAttr) ',8a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   rcode = 0

   n = mct_aVect_lsize(aVect)
   m = size(data)
   if (n /= m) then
      if (s_loglev > 0) write(s_logunit,*) subName,"ERROR: size aV,data,attr = ",n,m,trim(str)
      data = SHR_CONST_SPVAL
      rcode = 1
      return
   end if
   
   k = mct_aVect_indexRA(aVect,trim(str) ,perrWith=subName)
   if ( k < 1) then
      if (s_loglev > 0) write(s_logunit,*) subName,"ERROR: attribute not found, var = ",trim(str),", k=",k
      data = SHR_CONST_SPVAL
      rcode = 2
      return
   end if

   data(:) = aVect%rAttr(k,:)

end subroutine mct_aVect_getRAttr

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: mct_aVect_putRAttr - put real F90 array data into an aVect
!
! !DESCRIPTION:
!     Put the data in array {\tt data} into the  {\it AttributeVector}
!     {\tt aVect} under the attribute {\tt str}.
!     {\tt rcode} will be 0 if succesful, 1 if size of {\tt data}
!     does not match size  of {\tt aVect} and 2 if {\tt str} is not
!     found.
!
! !REMARKS:
!   This is like the MCT routine importRAttr except the output argument
!   is not a pointer.

! !REVISION HISTORY:
!     2002 Apr xx - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used subroutine mct_aVect_putRAttr(aVect,str,data,rcode)
! Subprogram not used 
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    type(mct_aVect),intent(inout) :: aVect ! Attribute vector
! Subprogram not used    character(*)   ,intent(in)  :: str
! Subprogram not used    real(R8)       ,intent(in)  :: data(:)
! Subprogram not used    integer(IN)    ,intent(out) :: rcode
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used 
! Subprogram not used    !--- local ---
! Subprogram not used    integer(IN) :: k,n,m
! Subprogram not used    integer(IN) :: aVsize
! Subprogram not used 
! Subprogram not used    !--- formats ---
! Subprogram not used    character(*),parameter :: subName = "(mct_aVect_putRAttr) "
! Subprogram not used    character(*),parameter :: F00 = "('(mct_aVect_putRAttr) ',8a)"
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    rcode = 0
! Subprogram not used 
! Subprogram not used    n = mct_aVect_lsize(aVect)
! Subprogram not used    m = size(data)
! Subprogram not used    if (n /= m) then
! Subprogram not used       if (s_loglev > 0) write(s_logunit,*) subName,"ERROR: size aV,data,attr = ",n,m,trim(str)
! Subprogram not used       rcode = 1
! Subprogram not used       return
! Subprogram not used    end if
! Subprogram not used    
! Subprogram not used    k = mct_aVect_indexRA(aVect,trim(str) ,perrWith=subName)
! Subprogram not used    if ( k < 1) then
! Subprogram not used       if (s_loglev > 0) write(s_logunit,*) subName,"ERROR: attribute not found, var = ",trim(str),", k=",k
! Subprogram not used       rcode = 2
! Subprogram not used       return
! Subprogram not used    end if
! Subprogram not used 
! Subprogram not used    aVect%rAttr(k,:) = data(:) 
! Subprogram not used 
! Subprogram not used end subroutine mct_aVect_putRAttr

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: mct_aVect_accum - accumulate attributes from one aVect to another
!
! !DESCRIPTION:
! This routine accumulates from input argment {\tt aVin} into the output 
! {\it AttrVect} argument {\tt aVout} the real and integer attributes specified in 
! input {\tt CHARACTER} argument {\tt iList} and {\tt rList}. The attributes can
! be listed in any order.  If neither {\tt iList} nor {\tt rList} are provided, 
! all attributes shared between {\tt aVin} and {\tt aVout} will be copied.
!
! If any attributes in {\tt aVout} have different names but represent the
! the same quantity and should still be copied, you must provide a translation
! argument {\tt TrList} and/or {\tt TiList}.  The translation arguments should
! be identical to the {\tt rList} or {\tt iList} but with the correct {\tt aVout}
! name subsititued at the appropriate place.
!
! This routine leverages the mct copy routines directly
!
! {\bf N.B.:}  This routine will fail if the {\tt aVout} is not initialized or
! if any of the specified attributes are not present in either {\tt aVout} or {\tt aVin}.
!
! !REVISION HISTORY:
!    2002 Sep 15 - ? - initial version.
!     2013-Jul-20 - T. Craig -- updated
!
! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used  subroutine mct_avect_accum(aVin, aVout, rList, TrList, iList, TiList, vector, sharedIndices,counter)
! Subprogram not used 
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !INPUT PARAMETERS: 
! Subprogram not used 
! Subprogram not used       type(mct_avect),            intent(in)    :: aVin
! Subprogram not used       character(len=*), optional, intent(in)    :: iList
! Subprogram not used       character(len=*), optional, intent(in)    :: rList
! Subprogram not used       character(len=*), optional, intent(in)    :: TiList
! Subprogram not used       character(len=*), optional, intent(in)    :: TrList
! Subprogram not used       logical, optional,          intent(in)    :: vector 
! Subprogram not used       type(mct_avect_SharedIndices), optional, intent(in) :: sharedIndices
! Subprogram not used 
! Subprogram not used ! !OUTPUT PARAMETERS: 
! Subprogram not used 
! Subprogram not used       type(mct_avect),         intent(inout) :: aVout
! Subprogram not used       integer, optional,       intent(inout) :: counter
! Subprogram not used 
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used 
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used    !--- local ---
! Subprogram not used   logical :: usevector
! Subprogram not used   integer(IN) :: lsize,nflds,npts,i,j
! Subprogram not used   type(mct_avect) :: avotmp  ! temporary aVout copy
! Subprogram not used   character(*),parameter :: subName = '(mct_aVect_accum) '
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------
! Subprogram not used 
! Subprogram not used   usevector = .false.
! Subprogram not used   if (present(vector)) then
! Subprogram not used      usevector = vector
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   if (present(counter)) then
! Subprogram not used      counter = counter + 1
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   ! --- allocate avotmp, a duplciate of aVout
! Subprogram not used 
! Subprogram not used   lsize = mct_aVect_lsize(aVout)
! Subprogram not used   call mct_avect_init(avotmp,aVout,lsize)
! Subprogram not used   call mct_avect_zero(avotmp)
! Subprogram not used 
! Subprogram not used   ! --- copy aVin fields into avotmp
! Subprogram not used 
! Subprogram not used   if (present(sharedIndices)) then
! Subprogram not used 
! Subprogram not used      if (present(rList) .and. present(iList)) then
! Subprogram not used         if (present(trList) .and. present(tilist)) then
! Subprogram not used            call mct_avect_copy(aVin, avotmp, rList, TrList, iList, tiList, vector = usevector, sharedIndices=sharedIndices)
! Subprogram not used         elseif (present(trList)) then
! Subprogram not used            call mct_avect_copy(aVin, avotmp, rList, TrList, iList, vector = usevector, sharedIndices=sharedIndices)
! Subprogram not used         elseif (present(tiList)) then
! Subprogram not used            call mct_avect_copy(aVin, avotmp, rList, iList=iList, tiList=tiList, vector = usevector, sharedIndices=sharedIndices)
! Subprogram not used         else
! Subprogram not used            call mct_avect_copy(aVin, avotmp, rList=rList, iList=iList, vector = usevector, sharedIndices=sharedIndices)
! Subprogram not used         endif
! Subprogram not used      else if (present(rList)) then
! Subprogram not used         if (present(trList)) then
! Subprogram not used            call mct_avect_copy(aVin, avotmp, rList, TrList, vector = usevector, sharedIndices=sharedIndices)
! Subprogram not used         else
! Subprogram not used            call mct_avect_copy(aVin, avotmp, rList, vector = usevector, sharedIndices=sharedIndices)
! Subprogram not used         endif
! Subprogram not used 
! Subprogram not used      else if (present(iList)) then
! Subprogram not used         if (present(tiList)) then
! Subprogram not used            call mct_avect_copy(aVin, avotmp, ilist=iList, tiList=tiList, vector = usevector, sharedIndices=sharedIndices)
! Subprogram not used         else
! Subprogram not used            call mct_avect_copy(aVin, avotmp, ilist=iList, vector = usevector, sharedIndices=sharedIndices)
! Subprogram not used         endif
! Subprogram not used 
! Subprogram not used      else
! Subprogram not used         call mct_avect_copy(aVin, avotmp, vector=usevector, sharedIndices=sharedIndices)
! Subprogram not used 
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used   else   ! sharedIndices
! Subprogram not used 
! Subprogram not used      if (present(rList) .and. present(iList)) then
! Subprogram not used         if (present(trList) .and. present(tilist)) then
! Subprogram not used            call mct_avect_copy(aVin, avotmp, rList, TrList, iList, tiList, vector = usevector)
! Subprogram not used         elseif (present(trList)) then
! Subprogram not used            call mct_avect_copy(aVin, avotmp, rList, TrList, iList, vector = usevector)
! Subprogram not used         elseif (present(tiList)) then
! Subprogram not used            call mct_avect_copy(aVin, avotmp, rList, iList=iList, tiList=tiList, vector = usevector)
! Subprogram not used         else
! Subprogram not used            call mct_avect_copy(aVin, avotmp, rList=rList, iList=iList, vector = usevector)
! Subprogram not used         endif
! Subprogram not used      else if (present(rList)) then
! Subprogram not used         if (present(trList)) then
! Subprogram not used            call mct_avect_copy(aVin, avotmp, rList, TrList, vector = usevector)
! Subprogram not used         else
! Subprogram not used            call mct_avect_copy(aVin, avotmp, rList, vector = usevector)
! Subprogram not used         endif
! Subprogram not used 
! Subprogram not used      else if (present(iList)) then
! Subprogram not used         if (present(tiList)) then
! Subprogram not used            call mct_avect_copy(aVin, avotmp, ilist=iList, tiList=tiList, vector = usevector)
! Subprogram not used         else
! Subprogram not used            call mct_avect_copy(aVin, avotmp, ilist=iList, vector = usevector)
! Subprogram not used         endif
! Subprogram not used 
! Subprogram not used      else
! Subprogram not used         call mct_avect_copy(aVin, avotmp, vector=usevector)
! Subprogram not used 
! Subprogram not used      endif
! Subprogram not used 
! Subprogram not used   endif ! shared indices
! Subprogram not used 
! Subprogram not used   ! --- accumulate avotmp into avout
! Subprogram not used 
! Subprogram not used   nflds = mct_aVect_nRAttr(aVout)
! Subprogram not used   npts  = mct_aVect_lsize (aVout)
! Subprogram not used !DIR$ CONCURRENT
! Subprogram not used !DIR$ PREFERVECTOR
! Subprogram not used   do i=1,npts 
! Subprogram not used   do j=1,nflds
! Subprogram not used      aVout%rattr(j,i) = aVout%rattr(j,i) + avotmp%rattr(j,i)
! Subprogram not used   enddo
! Subprogram not used   enddo
! Subprogram not used 
! Subprogram not used   ! --- clean avotmp
! Subprogram not used 
! Subprogram not used   call mct_avect_clean(avotmp)
! Subprogram not used 
! Subprogram not used  end subroutine mct_avect_accum

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: mct_aVect_avg - averages an accumulated attribute vector
!
! !DESCRIPTION:
!     Average the data in attribute vector {\tt aVect}.  Divides all fields in 
!     the attribute vector {\tt aVect} by the value of the input counter.
!
! !REVISION HISTORY:
!     2002-Sep-15 - T. Craig -- initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine mct_aVect_avg(aVect, counter)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(mct_aVect),intent(inout) :: aVect   ! bundle to read
   integer        ,intent(in)    :: counter ! counter 

!EOP

   !--- local ---
   integer(IN) :: i,j    ! generic indicies
   integer(IN) :: npts   ! number of points (local) in an aVect field
   integer(IN) :: nflds  ! number of aVect fields (real)
   real(R8)    :: ravg   ! accumulation count

   !--- formats ---
   character(*),parameter :: subName = '(mct_aVect_avg) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (counter == 0 .or. counter == 1) return

   ravg = 1.0_R8/real(counter,R8)

   nflds = mct_aVect_nRAttr(aVect)
   npts  = mct_aVect_lsize (aVect)
!DIR$ CONCURRENT
!DIR$ PREFERVECTOR
   do i=1,npts 
   do j=1,nflds
      aVect%rattr(j,i) = aVect%rattr(j,i)*ravg
   enddo
   enddo

end subroutine mct_aVect_avg

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: mct_avect_mult - multiply an attribute vector by a field.
!
! !DESCRIPTION:
!     Replace each field in {\tt av} by the product of that field and the
!     field {\tt fld1} from input argument {\tt av1}.
!
!     If optional argument {\tt bunlist} is present, only those attributes 
!     in {\tt bun} will be replaced.
!
!     If optional argument {\tt initav} is present, then the data in {\tt av}
!     is replaced by the product of the data in {\tt initav} and {\tt fld1}
!     from {\tt av1}. NOTE:  this assume {\tt initav} has the exact same
!     attributes in the same order as {\tt av}.
!
!
! !REVISION HISTORY:
!     2007-Jun-11 - M. Vertenstein -- initial version
!
! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used subroutine mct_avect_mult(av,av1,fld1,avlist)
! Subprogram not used 
! Subprogram not used ! !USES:
! Subprogram not used 
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    type(mct_aVect)      ,intent(inout) :: av       ! attribute vector output
! Subprogram not used    type(mct_aVect)      ,intent(in)    :: av1      ! attribute vector input
! Subprogram not used    character(*)         ,intent(in)    :: fld1     ! av1 field name
! Subprogram not used    character(*),optional,intent(in)    :: avlist   ! sublist of field in av
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used 
! Subprogram not used    !--- local ---
! Subprogram not used    integer(IN) :: n,m            ! generic indicies
! Subprogram not used    integer(IN) :: npts           ! number of points (local) in an aVect field
! Subprogram not used    integer(IN) :: nfld           ! number of fields (local) in an aVect field
! Subprogram not used    integer(IN) :: nfldi          ! number of fields (local) in an aVect field
! Subprogram not used    integer(IN) :: nptsx          ! number of points (local) in an aVect field
! Subprogram not used    integer(IN) :: nptsi          ! number of points (local) in an aVect field
! Subprogram not used    integer(IN) :: kfld           ! field number of fld1 in av1
! Subprogram not used    integer(IN),dimension(:),allocatable :: kfldin   ! field numbers of avlist in av
! Subprogram not used    type(mct_list)   :: blist     ! avlist as a List
! Subprogram not used    type(mct_string) :: tattr     ! an attribute
! Subprogram not used 
! Subprogram not used    !--- formats ---
! Subprogram not used    character(*),parameter :: subName = '(mct_aVect_mult) '
! Subprogram not used 
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used !
! Subprogram not used !-------------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used    nptsx = mct_aVect_lsize(av1)
! Subprogram not used    npts  = mct_aVect_lsize(av)
! Subprogram not used    if (nptsx /= npts .and. s_loglev > 0) write(s_logunit,*) subName,' ERROR: npts error1 ',npts,nptsx
! Subprogram not used 
! Subprogram not used    kfld  = mct_aVect_indexRA(av1,fld1,perrWith=subName)
! Subprogram not used 
! Subprogram not used    if (present(avlist)) then
! Subprogram not used 
! Subprogram not used      call mct_list_init(blist,avlist)
! Subprogram not used 
! Subprogram not used      nfld=mct_list_nitem(blist)
! Subprogram not used 
! Subprogram not used      allocate(kfldin(nfld))
! Subprogram not used      do m=1,nfld
! Subprogram not used        call mct_list_get(tattr,m,blist)
! Subprogram not used        kfldin(m) = mct_aVect_indexRA(av,mct_string_toChar(tattr))
! Subprogram not used        call mct_string_clean(tattr)
! Subprogram not used      enddo
! Subprogram not used      call mct_list_clean(blist)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used      do n=1,npts
! Subprogram not used      do m=1,nfld
! Subprogram not used 
! Subprogram not used         av%rAttr(kfldin(m),n) = av%rAttr(kfldin(m),n)*av1%rAttr(kfld,n)
! Subprogram not used      enddo
! Subprogram not used      enddo
! Subprogram not used 
! Subprogram not used      deallocate(kfldin)
! Subprogram not used 
! Subprogram not used    else
! Subprogram not used 
! Subprogram not used      nfld  = mct_aVect_nRAttr(av)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used      do n=1,npts
! Subprogram not used      do m=1,nfld
! Subprogram not used 
! Subprogram not used         av%rAttr(m,n) = av%rAttr(m,n)*av1%rAttr(kfld,n)
! Subprogram not used      enddo
! Subprogram not used      enddo
! Subprogram not used 
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used end subroutine mct_aVect_mult

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: mct_avect_vecmult - multiply an attribute vector by a field.
!
! !DESCRIPTION:
!     Replace each field in {\tt av} by the product of that field and the
!     field {\tt fld1} from input argument {\tt av1}.
!
!     If optional argument {\tt bunlist} is present, only those attributes 
!     in {\tt bun} will be replaced.
!
!     If optional argument {\tt initav} is present, then the data in {\tt av}
!     is replaced by the product of the data in {\tt initav} and {\tt fld1}
!     from {\tt av1}. NOTE:  this assume {\tt initav} has the exact same
!     attributes in the same order as {\tt av}.
!
!
! !REVISION HISTORY:
!     2007-Jun-11 - M. Vertenstein -- initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine mct_avect_vecmult(av,vec,avlist,mask_spval)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(mct_aVect)      ,intent(inout) :: av       ! attribute vector output
   real(R8)             ,intent(in)    :: vec(:)
   character(*),optional,intent(in)    :: avlist   ! sublist of field in av
   logical, optional    ,intent(in)    :: mask_spval

!EOP

   !--- local ---
   integer(IN) :: n,m            ! generic indicies
   integer(IN) :: npts           ! number of points (local) in an aVect field
   integer(IN) :: nfld           ! number of fields (local) in an aVect field
   integer(IN) :: nfldi          ! number of fields (local) in an aVect field
   integer(IN) :: nptsx          ! number of points (local) in an aVect field
   integer(IN) :: nptsi          ! number of points (local) in an aVect field
   logical     :: lmspval        ! local mask spval
   integer(IN),dimension(:),allocatable :: kfldin   ! field numbers of avlist in av
   type(mct_list)   :: blist     ! avlist as a List
   type(mct_string) :: tattr     ! an attribute

   !--- formats ---
   character(*),parameter :: subName = '(mct_aVect_vecmult) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   lmspval = .false.
   if (present(mask_spval)) then
      lmspval = mask_spval
   endif

   nptsx = size(vec,1)
   npts  = mct_aVect_lsize(av)
   if (nptsx /= npts .and. s_loglev > 0) write(s_logunit,*) subName,' ERROR: npts error1 ',npts,nptsx


   if (present(avlist)) then

     call mct_list_init(blist,avlist)

     nfld=mct_list_nitem(blist)

     allocate(kfldin(nfld))
     do m=1,nfld
       call mct_list_get(tattr,m,blist)
       kfldin(m) = mct_aVect_indexRA(av,mct_string_toChar(tattr))
       call mct_string_clean(tattr)
     enddo
     call mct_list_clean(blist)

     if (lmspval) then







        do n=1,npts
        do m=1,nfld

           if (.not. shr_const_isspval(av%rAttr(kfldin(m),n))) then
              av%rAttr(kfldin(m),n) = av%rAttr(kfldin(m),n)*vec(n)
           endif
        enddo
        enddo

     else  ! lmspval







        do n=1,npts
        do m=1,nfld

           av%rAttr(kfldin(m),n) = av%rAttr(kfldin(m),n)*vec(n)
        enddo
        enddo

     endif  ! lmspval

     deallocate(kfldin)

   else  ! avlist

     nfld  = mct_aVect_nRAttr(av)

     if (lmspval) then







        do n=1,npts
        do m=1,nfld

           if (.not. shr_const_isspval(av%rAttr(m,n))) then
              av%rAttr(m,n) = av%rAttr(m,n)*vec(n)
           endif
        enddo
        enddo

     else  ! lmspval







        do n=1,npts
        do m=1,nfld

           av%rAttr(m,n) = av%rAttr(m,n)*vec(n)
        enddo
        enddo

     endif   ! lmspval

   endif   ! avlist

end subroutine mct_aVect_vecmult

!===============================================================================
! !BOP ===========================================================================
!
! !IROUTINE:  subroutine mct_rearr_rearrange_fldlst - rearrange on a fieldlist
!
! !DESCRIPTION: 
!     Perform regarranger between two attribute vectors only on the fieldlist
!     that is provided
!
!
! !REVISION HISTORY: 
!    2007-Jun-22 - M. Vertenstein - first version
! 
! !INTERFACE:  -----------------------------------------------------------------

subroutine mct_rearr_rearrange_fldlist(avi, avo, Rearr, vector, alltoall, fldlist, tag)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(mct_aVect) , intent(in)  :: avi
   type(mct_aVect) , intent(inout):: avo
   type(mct_rearr) , intent(in)  :: Rearr
   logical         , intent(in)  :: vector
   logical         , intent(in)  :: alltoall
   character(len=*), intent(in)  :: fldlist   
   integer(IN)     , intent(in),optional :: tag
! !EOP

   !---local ---
   type(mct_aVect) :: avi_fl
   type(mct_aVect) :: avo_fl
   integer(IN)     :: lsize
   integer(IN)     :: ltag

   !--- formats ---
   character(*),parameter :: subName = '(mct_rearr_rearrange_fldlist) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (present(tag)) then
      ltag = tag
   else
      ltag = 3000
   endif

   lsize = mct_aVect_lsize(avi)
   call mct_aVect_init (avi_fl, rlist=fldlist, lsize=lsize)
   call mct_aVect_zero (avi_fl)

   lsize = mct_aVect_lsize(avo)
   call mct_aVect_init (avo_fl, rlist=fldlist, lsize=lsize)
   call mct_aVect_zero (avo_fl)
   
   call mct_aVect_copy (aVin=avi, aVout=avi_fl)
   call mct_rearr_rearrange(avi_fl, avo_fl, Rearr, VECTOR=vector, ALLTOALL=alltoall, tag=ltag)
   call mct_aVect_copy (aVin=avo_fl, aVout=avo, vector=vector)

   call mct_aVect_clean(avi_fl)
   call mct_aVect_clean(avo_fl)

end subroutine mct_rearr_rearrange_fldlist

!=======================================================================
logical function mct_gsmap_Identical(gsmap1,gsmap2)

  implicit none
  type(mct_gsMap), intent(IN):: gsmap1
  type(mct_gsMap), intent(IN):: gsmap2

  ! Local variables

  character(len=*),parameter :: subname = "(mct_gsmap_Identical) "
  integer :: n
  logical :: identical

  !-----------------------

  identical = .true.

  ! --- continue compare ---
  if (identical) then
     if (mct_gsMap_gsize(gsmap1) /= mct_gsMap_gsize(gsmap2)) identical = .false.
     if (mct_gsMap_ngseg(gsmap1) /= mct_gsMap_ngseg(gsmap2)) identical = .false.
  endif

  ! --- continue compare ---
  if (identical) then
     do n = 1,mct_gsMap_ngseg(gsmap1)
        if (gsmap1%start(n)  /= gsmap2%start(n) ) identical = .false.
        if (gsmap1%length(n) /= gsmap2%length(n)) identical = .false.
        if (gsmap1%pe_loc(n) /= gsmap2%pe_loc(n)) identical = .false.
     enddo
  endif

  mct_gsmap_Identical = identical

end function mct_gsmap_Identical
    
!===============================================================================
! !BOP ===========================================================================
!
! !IROUTINE:  mct_myindex - binary search for index in list
!
! !DESCRIPTION: 
!     Do a binary search to see if a value is contained in a list of
!     values.  return true or false.  starti must be monotonically
!     increasing, function does NOT check this.
!
!
! !REVISION HISTORY: 
!    2007-Jan-17 - T. Craig -- first version
!    2007-Mar-20 - R. Jacob - move to mct_mod
! 
! !INTERFACE:  -----------------------------------------------------------------

logical function mct_myindex(index,starti,counti)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   integer(IN) :: index       ! is this index in start/count list
   integer(IN) :: starti(:)   ! start list
   integer(IN) :: counti(:)   ! count list

! !EOP

   !--- local ---
   integer(IN)    :: nl,nc,nr,ncprev 
   integer(IN)    :: lsize
   logical        :: stopnow

   !--- formats ---
   character(*),parameter :: subName = '(mct_myindex) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   mct_myindex = .false.

   lsize = size(starti)
   if (lsize < 1) return

   nl = 0
   nr = lsize + 1
   nc = (nl+nr)/2
   stopnow = .false.
   do while (.not.stopnow)
      if (index < starti(nc)) then
         nr = nc
      elseif (index > (starti(nc) + counti(nc) - 1)) then
         nl = nc
      else
         mct_myindex = .true.
         return
      endif
      ncprev = nc
      nc = (nl + nr)/2
      if (nc == ncprev .or. nc < 1 .or. nc > lsize) stopnow = .true.
   enddo

   mct_myindex = .false.
   return

end function mct_myindex
!===============================================================================

end module mct_mod

