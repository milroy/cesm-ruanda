!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
! CVS $Id$
! CVS $Name$  
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: m_mpif90 - a Fortran 90 style MPI module interface.
!
! !DESCRIPTION:
!
!   By wrapping \verb'include "mpif.h"' into a module, \verb"m_mpif()"
!   provides an easy way to
!\begin{itemize}
!  \item avoid the problem with {\sl fixed} or {\sl free} formatted
!	Fortran 90 files;
!  \item provide protections with only a limited set of \verb"PUBLIC"
!	variables; and
!  \item be extended to a MPI Fortran 90 binding.
!\end{itemize}
!
! !INTERFACE:

    module m_mpif90
      use m_mpif, only : MP_INTEGER	=> MPI_INTEGER
      use m_mpif, only : MP_REAL	=> MPI_REAL
      use m_mpif, only : MP_DOUBLE_PRECISION	&
					=> MPI_DOUBLE_PRECISION
      use m_mpif, only : MP_LOGICAL	=> MPI_LOGICAL
      use m_mpif, only : MP_CHARACTER	=> MPI_CHARACTER

      use m_mpif, only : MP_REAL4	=> MPI_REAL4
      use m_mpif, only : MP_REAL8	=> MPI_REAL8

      use m_mpif, only : MP_COMM_WORLD	=> MPI_COMM_WORLD
      use m_mpif, only : MP_COMM_NULL	=> MPI_COMM_NULL
      use m_mpif, only : MP_SUM		=> MPI_SUM
      use m_mpif, only : MP_PROD	=> MPI_PROD
      use m_mpif, only : MP_MIN 	=> MPI_MIN
      use m_mpif, only : MP_MAX 	=> MPI_MAX
      use m_mpif, only : MP_MAX_ERROR_STRING	&
					=> MPI_MAX_ERROR_STRING
      use m_mpif, only : MP_STATUS_SIZE => MPI_STATUS_SIZE
      use m_mpif, only : MP_ANY_SOURCE	=> MPI_ANY_SOURCE

      implicit none
      private

      public :: MP_type

      public :: MP_INTEGER
      public :: MP_REAL
      public :: MP_DOUBLE_PRECISION
      public :: MP_LOGICAL
      public :: MP_CHARACTER

      public :: MP_REAL4
      public :: MP_REAL8

      public :: MP_COMM_WORLD
      public :: MP_COMM_NULL

      public :: MP_SUM
      public :: MP_PROD
      public :: MP_MIN
      public :: MP_MAX

      public :: MP_ANY_SOURCE

      public :: MP_MAX_ERROR_STRING

      public :: MP_init
      public :: MP_initialized
      public :: MP_finalize
      public :: MP_abort

      public :: MP_wtime
      public :: MP_wtick

      public :: MP_comm_size
      public :: MP_comm_rank
      public :: MP_comm_dup
      public :: MP_comm_free

      public :: MP_cart_create
      public :: MP_dims_create
      public :: MP_cart_coords
      public :: MP_cart_rank

      public :: MP_error_string

      public :: MP_perr

      public :: MP_STATUS_SIZE
      public :: MP_status

      public :: MP_log2

! !REVISION HISTORY:
! 	09Dec97 - Jing Guo <guo@thunder> - initial prototyping/coding.
!		. started with everything public, without any interface
!		  declaration.
!		. Then limited to only variables current expected to
!		  be used.
!	
!EOP
!_______________________________________________________________________

integer,dimension(MP_STATUS_SIZE) :: MP_status

	!----------------------------------------

interface MP_init
  subroutine MPI_init(ier)
    integer :: ier
  end subroutine MPI_init
end interface

interface MP_initialized
  subroutine MPI_initialized(flag,ier)
    logical :: flag
    integer :: ier
  end subroutine MPI_initialized
end interface

interface MP_finalize
  subroutine MPI_finalize(ier)
    integer :: ier
  end subroutine MPI_finalize
end interface

interface MP_error_string
  subroutine MPI_error_string(ierror,cerror,ln,ier)
    integer :: ierror
    character(len=*) :: cerror
    integer :: ln
    integer :: ier
  end subroutine MPI_error_string
end interface

interface MP_type; module procedure	&
  typeI_,	& ! MPI_INTEGER
  typeL_,	& ! MPI_LOGICAL
  typeC_,	& ! MPI_CHARACTER
  typeSP_,	& ! MPI_REAL
  typeDP_,	& ! MPI_DOUBLE_PRECISION
  typeI1_,	& ! MPI_INTEGER
  typeL1_,	& ! MPI_LOGICAL
  typeC1_,	& ! MPI_CHARACTER
  typeSP1_,	& ! MPI_REAL
  typeDP1_,	& ! MPI_DOUBLE_PRECISION
  typeI2_,	& ! MPI_INTEGER
  typeL2_,	& ! MPI_LOGICAL
  typeC2_,	& ! MPI_CHARACTER
  typeSP2_,	& ! MPI_REAL
  typeDP2_	  ! MPI_DOUBLE_PRECISION
end interface

interface MP_perr; module procedure perr_; end interface

interface MP_abort
  subroutine MPI_abort(comm,errorcode,ier)
    integer :: comm
    integer :: errorcode
    integer :: ier
  end subroutine MPI_abort
end interface

	!----------------------------------------
interface MP_wtime
  function MPI_wtime()
    double precision :: MPI_wtime
  end function MPI_wtime
end interface

interface MP_wtick
  function MPI_wtick()
    double precision :: MPI_wtick
  end function MPI_wtick
end interface

	!----------------------------------------
interface MP_comm_size
  subroutine MPI_comm_size(comm,size,ier)
    integer :: comm
    integer :: size
    integer :: ier
  end subroutine MPI_comm_size
end interface

interface MP_comm_rank
  subroutine MPI_comm_rank(comm,rank,ier)
    integer :: comm
    integer :: rank
    integer :: ier
  end subroutine MPI_comm_rank
end interface

interface MP_comm_dup
  subroutine MPI_comm_dup(comm,newcomm,ier)
    integer :: comm
    integer :: newcomm
    integer :: ier
  end subroutine MPI_comm_dup
end interface

interface MP_comm_free
  subroutine MPI_comm_free(comm,ier)
    integer :: comm
    integer :: ier
  end subroutine MPI_comm_free
end interface

	!----------------------------------------
interface MP_cart_create
  subroutine MPI_cart_create(comm_old,ndims,dims,periods,	&
  	reorder,comm_cart,ier)
    integer :: comm_old
    integer :: ndims
    integer,dimension(*) :: dims
    logical,dimension(*) :: periods
    logical :: reorder
    integer :: comm_cart
    integer :: ier
  end subroutine MPI_cart_create
end interface

interface MP_dims_create
  subroutine MPI_dims_create(nnodes,ndims,dims,ier)
    integer :: nnodes
    integer :: ndims
    integer,dimension(*) :: dims
    integer :: ier
  end subroutine MPI_dims_create
end interface

interface MP_cart_coords
  subroutine MPI_cart_coords(comm,rank,maxdims,coords,ier)
    integer :: comm
    integer :: rank
    integer :: maxdims
    integer,dimension(*) :: coords
    integer :: ier
  end subroutine MPI_cart_coords
end interface

interface MP_cart_rank
  subroutine MPI_cart_rank(comm,coords,rank,ier)
    integer :: comm
    integer,dimension(*) :: coords
    integer :: rank
    integer :: ier
  end subroutine MPI_cart_rank
end interface
	!----------------------------------------

  character(len=*),parameter :: myname='m_mpif90'
contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: typeI_ - return MPI datatype of INTEGER
!
! !DESCRIPTION:
!
! !INTERFACE:

    function typeI_(ival)
      implicit none
      integer,intent(in) :: ival
      integer :: typeI_

! !REVISION HISTORY:
! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::typeI_'

  typeI_=MP_INTEGER

end function typeI_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: typeL_ - return MPI datatype of LOGICAL
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     function typeL_(lval)
! Subprogram not used       implicit none
! Subprogram not used       logical,intent(in) :: lval
! Subprogram not used       integer :: typeL_
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::typeL_'
! Subprogram not used 
! Subprogram not used   typeL_=MP_LOGICAL
! Subprogram not used 
! Subprogram not used end function typeL_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: typeC_ - return MPI datatype of CHARACTER
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     function typeC_(cval)
! Subprogram not used       implicit none
! Subprogram not used       character(len=*),intent(in) :: cval
! Subprogram not used       integer :: typeC_
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::typeC_'
! Subprogram not used 
! Subprogram not used   typeC_=MP_CHARACTER
! Subprogram not used 
! Subprogram not used end function typeC_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: typeSP_ - return MPI datatype of single precision REAL
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     function typeSP_(rval)
! Subprogram not used       use m_realkinds,only : SP
! Subprogram not used       implicit none
! Subprogram not used       real(SP),intent(in) :: rval
! Subprogram not used       integer :: typeSP_
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::typeSP_'
! Subprogram not used 
! Subprogram not used   typeSP_=MP_REAL
! Subprogram not used 
! Subprogram not used end function typeSP_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: typeDP_ - return MPI datatype of double precision REAL
!
! !DESCRIPTION:
!
! !INTERFACE:

    function typeDP_(rval)
      use m_realkinds,only : DP
      implicit none
      real(DP),intent(in) :: rval
      integer :: typeDP_

! !REVISION HISTORY:
! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::typeDP_'

  typeDP_=MP_DOUBLE_PRECISION

end function typeDP_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: typeI1_ - return MPI datatype of INTEGER
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     function typeI1_(ival)
! Subprogram not used       implicit none
! Subprogram not used       integer,dimension(:),intent(in) :: ival
! Subprogram not used       integer :: typeI1_
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::typeI1_'
! Subprogram not used 
! Subprogram not used   typeI1_=MP_INTEGER
! Subprogram not used 
! Subprogram not used end function typeI1_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: typeL1_ - return MPI datatype of LOGICAL
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     function typeL1_(lval)
! Subprogram not used       implicit none
! Subprogram not used       logical,dimension(:),intent(in) :: lval
! Subprogram not used       integer :: typeL1_
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::typeL1_'
! Subprogram not used 
! Subprogram not used   typeL1_=MP_LOGICAL
! Subprogram not used 
! Subprogram not used end function typeL1_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: typeC1_ - return MPI datatype of CHARACTER
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     function typeC1_(cval)
! Subprogram not used       implicit none
! Subprogram not used       character(len=*),dimension(:),intent(in) :: cval
! Subprogram not used       integer :: typeC1_
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::typeC1_'
! Subprogram not used 
! Subprogram not used   typeC1_=MP_CHARACTER
! Subprogram not used 
! Subprogram not used end function typeC1_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: typeSP1_ - return MPI datatype of single precision REAL
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     function typeSP1_(rval)
! Subprogram not used       use m_realkinds,only : SP
! Subprogram not used       implicit none
! Subprogram not used       real(SP),dimension(:),intent(in) :: rval
! Subprogram not used       integer :: typeSP1_
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::typeSP1_'
! Subprogram not used 
! Subprogram not used   typeSP1_=MP_REAL
! Subprogram not used 
! Subprogram not used end function typeSP1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: typeDP1_ - return MPI datatype of double precision REAL
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     function typeDP1_(rval)
! Subprogram not used       use m_realkinds,only : DP
! Subprogram not used       implicit none
! Subprogram not used       real(DP),dimension(:),intent(in) :: rval
! Subprogram not used       integer :: typeDP1_
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::typeDP1_'
! Subprogram not used 
! Subprogram not used   typeDP1_=MP_DOUBLE_PRECISION
! Subprogram not used 
! Subprogram not used end function typeDP1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: typeI2_ - return MPI datatype of INTEGER
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     function typeI2_(ival)
! Subprogram not used       implicit none
! Subprogram not used       integer,dimension(:,:),intent(in) :: ival
! Subprogram not used       integer :: typeI2_
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::typeI2_'
! Subprogram not used 
! Subprogram not used   typeI2_=MP_INTEGER
! Subprogram not used 
! Subprogram not used end function typeI2_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: typeL2_ - return MPI datatype of LOGICAL
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     function typeL2_(lval)
! Subprogram not used       implicit none
! Subprogram not used       logical,dimension(:,:),intent(in) :: lval
! Subprogram not used       integer :: typeL2_
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::typeL2_'
! Subprogram not used 
! Subprogram not used   typeL2_=MP_LOGICAL
! Subprogram not used 
! Subprogram not used end function typeL2_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: typeC2_ - return MPI datatype of CHARACTER
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     function typeC2_(cval)
! Subprogram not used       implicit none
! Subprogram not used       character(len=*),dimension(:,:),intent(in) :: cval
! Subprogram not used       integer :: typeC2_
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::typeC2_'
! Subprogram not used 
! Subprogram not used   typeC2_=MP_CHARACTER
! Subprogram not used 
! Subprogram not used end function typeC2_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: typeSP2_ - return MPI datatype of single precision REAL
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     function typeSP2_(rval)
! Subprogram not used       use m_realkinds,only : SP
! Subprogram not used       implicit none
! Subprogram not used       real(SP),dimension(:,:),intent(in) :: rval
! Subprogram not used       integer :: typeSP2_
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::typeSP2_'
! Subprogram not used 
! Subprogram not used   typeSP2_=MP_REAL
! Subprogram not used 
! Subprogram not used end function typeSP2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: typeDP2_ - return MPI datatype of double precision REAL
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     function typeDP2_(rval)
! Subprogram not used       use m_realkinds,only : DP
! Subprogram not used       implicit none
! Subprogram not used       real(DP),dimension(:,:),intent(in) :: rval
! Subprogram not used       integer :: typeDP2_
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	28Sep99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::typeDP2_'
! Subprogram not used 
! Subprogram not used   typeDP2_=MP_DOUBLE_PRECISION
! Subprogram not used 
! Subprogram not used end function typeDP2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: perr_ - MPI error information hanlder
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine perr_(proc,MP_proc,ierror)
! Subprogram not used       use m_stdio, only : stderr
! Subprogram not used       implicit none
! Subprogram not used       character(len=*),intent(in) :: proc
! Subprogram not used       character(len=*),intent(in) :: MP_proc
! Subprogram not used       integer,intent(in) :: ierror
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::perr_'
! Subprogram not used 
! Subprogram not used   character(len=MP_MAX_ERROR_STRING) :: estr
! Subprogram not used   integer :: ln,ier
! Subprogram not used 
! Subprogram not used   call MP_error_string(ierror,estr,ln,ier)
! Subprogram not used   if(ier /= 0 .or. ln<=0) then
! Subprogram not used     write(stderr,'(4a,i4)') proc,': ',	&
! Subprogram not used 	MP_proc,' error, ierror =',ierror
! Subprogram not used   else
! Subprogram not used     write(stderr,'(6a)') proc,': ',	&
! Subprogram not used 	MP_proc,' error, "',estr(1:ln),'"'
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used end subroutine perr_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: MP_log2 - The smallest integer its power of 2 is >= nPE
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     function MP_log2(nPE)
! Subprogram not used       implicit none
! Subprogram not used       integer,intent(in) :: nPE
! Subprogram not used       integer :: MP_log2
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	01Feb00	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::MP_log2'
! Subprogram not used 
! Subprogram not used   integer :: n2
! Subprogram not used 
! Subprogram not used   MP_log2=0
! Subprogram not used   n2=1
! Subprogram not used   do while(n2<nPE)
! Subprogram not used     MP_log2 = MP_log2+1
! Subprogram not used     n2 = n2+n2
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used end function MP_log2

end module m_mpif90
!.
