!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
! CVS $Id$
! CVS $Name$  
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: m_mall - A bookkeeper of user allocated memories
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_mall
      implicit none
      private	! except

      public :: mall_ci
      public :: mall_co
      public :: mall_mci
      public :: mall_mco
      public :: mall_flush
      public :: mall_reset

		! mall_ activity controls

      public :: mall_ison
      public :: mall_set

      interface mall_ci;    module procedure ci_; end interface
      interface mall_co;    module procedure co_; end interface

      interface mall_mci;    module procedure	&
	ciI0_,	&
	ciI1_,	&
	ciI2_,	&
	ciI3_,	&
	ciR0_,	&
	ciR1_,	&
	ciR2_,	&
	ciR3_,	&
	ciD0_,	&
	ciD1_,	&
	ciD2_,	&
	ciD3_,	&
	ciL0_,	&
	ciL1_,	&
	ciL2_,	&
	ciL3_,	&
	ciC0_,	&
	ciC1_,	&
	ciC2_,	&
	ciC3_
      end interface

      interface mall_mco;    module procedure	&
	coI0_,	&
	coI1_,	&
	coI2_,	&
	coI3_,	&
	coR0_,	&
	coR1_,	&
	coR2_,	&
	coR3_,	&
	coD0_,	&
	coD1_,	&
	coD2_,	&
	coD3_,	&
	coL0_,	&
	coL1_,	&
	coL2_,	&
	coL3_,	&
	coC0_,	&
	coC1_,	&
	coC2_,	&
	coC3_
      end interface

      interface mall_flush; module procedure flush_; end interface
      interface mall_reset; module procedure reset_; end interface

      interface mall_ison; module procedure ison_; end interface
      interface mall_set;  module procedure set_;  end interface

! !REVISION HISTORY:
! 	13Feb98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP
!_______________________________________________________________________
  character(len=*),parameter :: myname='MCT(MPEU)::m_mall'

  integer,parameter :: NBYTE_PER_WORD = 4

  integer,parameter :: NSZ= 32
  integer,parameter :: MXL=250

  integer, save :: nreset = 0		! number of reset_() calls
  logical, save :: started = .false.	! the module is in use

  integer, save :: n_ =0		! number of accouting bins.
  character(len=NSZ),dimension(MXL),save :: name_

  ! integer, dimension(1) :: mall
					! names of the accouting bins

  logical,save :: mall_on=.false.	! mall activity switch

  integer,save :: mci
  integer,dimension(MXL),save :: mci_	! maximum ci_() calls
  integer,save :: nci
  integer,dimension(MXL),save :: nci_	! net ci_() calls
  integer,save :: hwm
  integer,dimension(MXL),save :: hwm_	! high-water-mark of allocate()
  integer,save :: nwm
  integer,dimension(MXL),save :: nwm_	! net-water-mark of allocate()

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ison_ -
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ison_()
      implicit none
      logical :: ison_

! !REVISION HISTORY:
! 	25Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ison_'

  ison_=mall_on

end function ison_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: set_ - set the switch on
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine set_(on)
! Subprogram not used       implicit none
! Subprogram not used       logical,optional,intent(in) :: on
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	25Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::set_'
! Subprogram not used 
! Subprogram not used   mall_on=.true.
! Subprogram not used   if(present(on)) mall_on=on
! Subprogram not used 
! Subprogram not used end subroutine set_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ciI0_ - check in as an integer scalar
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine ciI0_(marg,thread)
! Subprogram not used       implicit none
! Subprogram not used       integer,intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::ciI0_'
! Subprogram not used 
! Subprogram not used   if(mall_on) call ci_(1,thread)
! Subprogram not used 
! Subprogram not used end subroutine ciI0_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ciI1_ - check in as an integer rank 1 array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine ciI1_(marg,thread)
! Subprogram not used       implicit none
! Subprogram not used       integer,dimension(:),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::ciI1_'
! Subprogram not used 
! Subprogram not used   if(mall_on) call ci_(size(marg),thread)
! Subprogram not used 
! Subprogram not used end subroutine ciI1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ciI2_ - check in as an integer rank 2 array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine ciI2_(marg,thread)
! Subprogram not used       implicit none
! Subprogram not used       integer,dimension(:,:),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::ciI2_'
! Subprogram not used 
! Subprogram not used   if(mall_on) call ci_(size(marg),thread)
! Subprogram not used 
! Subprogram not used end subroutine ciI2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ciI3_ - check in as an integer rank 3 array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine ciI3_(marg,thread)
! Subprogram not used       implicit none
! Subprogram not used       integer,dimension(:,:,:),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::ciI3_'
! Subprogram not used 
! Subprogram not used   if(mall_on) call ci_(size(marg),thread)
! Subprogram not used 
! Subprogram not used end subroutine ciI3_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ciR0_ - check in as a real(SP) scalar
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine ciR0_(marg,thread)
! Subprogram not used       use m_realkinds, only : SP
! Subprogram not used       implicit none
! Subprogram not used       real(SP),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::ciR0_'
! Subprogram not used 
! Subprogram not used   if(mall_on) call ci_(1,thread)
! Subprogram not used 
! Subprogram not used end subroutine ciR0_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ciR1_ - check in as a real(SP) rank 1 array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine ciR1_(marg,thread)
! Subprogram not used       use m_realkinds, only : SP
! Subprogram not used       implicit none
! Subprogram not used       real(SP),dimension(:),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::ciR1_'
! Subprogram not used 
! Subprogram not used   if(mall_on) call ci_(size(marg),thread)
! Subprogram not used 
! Subprogram not used end subroutine ciR1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ciR2_ - check in as a real(SP) rank 2 array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine ciR2_(marg,thread)
! Subprogram not used       use m_realkinds, only : SP
! Subprogram not used       implicit none
! Subprogram not used       real(SP),dimension(:,:),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::ciR2_'
! Subprogram not used 
! Subprogram not used   if(mall_on) call ci_(size(marg),thread)
! Subprogram not used 
! Subprogram not used end subroutine ciR2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ciR3_ - check in as a real(SP) rank 3 array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine ciR3_(marg,thread)
! Subprogram not used       use m_realkinds, only : SP
! Subprogram not used       implicit none
! Subprogram not used       real(SP),dimension(:,:,:),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::ciR3_'
! Subprogram not used 
! Subprogram not used   if(mall_on) call ci_(size(marg),thread)
! Subprogram not used 
! Subprogram not used end subroutine ciR3_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ciD0_ - check in as a real(DP) scalar
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine ciD0_(marg,thread)
! Subprogram not used       use m_realkinds, only : DP
! Subprogram not used       implicit none
! Subprogram not used       real(DP),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::ciD0_'
! Subprogram not used 
! Subprogram not used   if(mall_on) call ci_(2,thread)
! Subprogram not used 
! Subprogram not used end subroutine ciD0_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ciD1_ - check in as a real(DP) rank 1 array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine ciD1_(marg,thread)
! Subprogram not used       use m_realkinds, only : DP
! Subprogram not used       implicit none
! Subprogram not used       real(DP),dimension(:),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::ciD1_'
! Subprogram not used 
! Subprogram not used   if(mall_on) call ci_(2*size(marg),thread)
! Subprogram not used 
! Subprogram not used end subroutine ciD1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ciD2_ - check in as a real(DP) rank 2 array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine ciD2_(marg,thread)
! Subprogram not used       use m_realkinds, only : DP
! Subprogram not used       implicit none
! Subprogram not used       real(DP),dimension(:,:),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::ciD2_'
! Subprogram not used 
! Subprogram not used   if(mall_on) call ci_(2*size(marg),thread)
! Subprogram not used 
! Subprogram not used end subroutine ciD2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ciD3_ - check in as a real(DP) rank 3 array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine ciD3_(marg,thread)
! Subprogram not used       use m_realkinds, only : DP
! Subprogram not used       implicit none
! Subprogram not used       real(DP),dimension(:,:,:),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::ciD3_'
! Subprogram not used 
! Subprogram not used   if(mall_on) call ci_(2*size(marg),thread)
! Subprogram not used 
! Subprogram not used end subroutine ciD3_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ciL0_ - check in as a logical scalar
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine ciL0_(marg,thread)
! Subprogram not used       implicit none
! Subprogram not used       logical,intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::ciL0_'
! Subprogram not used 
! Subprogram not used   if(mall_on) call ci_(1,thread)
! Subprogram not used 
! Subprogram not used end subroutine ciL0_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ciL1_ - check in as a logical rank 1 array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine ciL1_(marg,thread)
! Subprogram not used       implicit none
! Subprogram not used       logical,dimension(:),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::ciL1_'
! Subprogram not used 
! Subprogram not used   if(mall_on) call ci_(size(marg),thread)
! Subprogram not used 
! Subprogram not used end subroutine ciL1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ciL2_ - check in as a logical rank 2 array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine ciL2_(marg,thread)
! Subprogram not used       implicit none
! Subprogram not used       logical,dimension(:,:),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::ciL2_'
! Subprogram not used 
! Subprogram not used   if(mall_on) call ci_(size(marg),thread)
! Subprogram not used 
! Subprogram not used end subroutine ciL2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ciL3_ - check in as a logical rank 3 array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine ciL3_(marg,thread)
! Subprogram not used       implicit none
! Subprogram not used       logical,dimension(:,:,:),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::ciL3_'
! Subprogram not used 
! Subprogram not used   if(mall_on) call ci_(size(marg),thread)
! Subprogram not used 
! Subprogram not used end subroutine ciL3_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ciC0_ - check in as a character scalar
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine ciC0_(marg,thread)
! Subprogram not used       implicit none
! Subprogram not used       character(len=*),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::ciC0_'
! Subprogram not used   integer :: nw
! Subprogram not used 
! Subprogram not used   if(.not.mall_on) return
! Subprogram not used   nw=(len(marg)+NBYTE_PER_WORD-1)/NBYTE_PER_WORD
! Subprogram not used   call ci_(nw,thread)
! Subprogram not used 
! Subprogram not used end subroutine ciC0_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ciC1_ - check in as a character rank 1 array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine ciC1_(marg,thread)
! Subprogram not used       implicit none
! Subprogram not used       character(len=*),dimension(:),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::ciC1_'
! Subprogram not used   integer :: nw
! Subprogram not used 
! Subprogram not used   if(.not.mall_on) return
! Subprogram not used   nw=(len(marg(1))+NBYTE_PER_WORD-1)/NBYTE_PER_WORD
! Subprogram not used   call ci_(size(marg)*nw,thread)
! Subprogram not used 
! Subprogram not used end subroutine ciC1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ciC2_ - check in as a character rank 2 array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine ciC2_(marg,thread)
! Subprogram not used       implicit none
! Subprogram not used       character(len=*),dimension(:,:),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::ciC2_'
! Subprogram not used   integer :: nw
! Subprogram not used 
! Subprogram not used   if(.not.mall_on) return
! Subprogram not used   nw=(len(marg(1,1))+NBYTE_PER_WORD-1)/NBYTE_PER_WORD
! Subprogram not used   call ci_(size(marg)*nw,thread)
! Subprogram not used 
! Subprogram not used end subroutine ciC2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ciC3_ - check in as a character rank 3 array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine ciC3_(marg,thread)
! Subprogram not used       implicit none
! Subprogram not used       character(len=*),dimension(:,:,:),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::ciC3_'
! Subprogram not used   integer :: nw
! Subprogram not used 
! Subprogram not used   if(.not.mall_on) return
! Subprogram not used   nw=(len(marg(1,1,1))+NBYTE_PER_WORD-1)/NBYTE_PER_WORD
! Subprogram not used   call ci_(size(marg)*nw,thread)
! Subprogram not used 
! Subprogram not used end subroutine ciC3_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ci_ - check-in allocate activity
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine ci_(nword,thread)
! Subprogram not used       use m_stdio, only : stderr
! Subprogram not used       use m_die, only : die
! Subprogram not used       implicit none
! Subprogram not used       integer,intent(in) :: nword
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	13Feb98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! Subprogram not used !EOP
! Subprogram not used !_______________________________________________________________________
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::ci_'
! Subprogram not used   integer :: ith
! Subprogram not used 
! Subprogram not used   if(.not.mall_on) return
! Subprogram not used 
! Subprogram not used   if(nword < 0) then
! Subprogram not used     write(stderr,'(2a,i4)') myname_,	&
! Subprogram not used 	': invalide argument, nword = ',nword
! Subprogram not used     call die(myname_)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   ith=lookup_(thread)
! Subprogram not used 
! Subprogram not used 	! update the account
! Subprogram not used 
! Subprogram not used   nci_(ith)=nci_(ith)+1
! Subprogram not used   mci_(ith)=mci_(ith)+1
! Subprogram not used   nwm_(ith)=nwm_(ith)+nword
! Subprogram not used   if(hwm_(ith).lt.nwm_(ith)) hwm_(ith)=nwm_(ith)
! Subprogram not used 
! Subprogram not used 	! update the total budget
! Subprogram not used 
! Subprogram not used   nci=nci+1
! Subprogram not used   mci=mci+1
! Subprogram not used   nwm=nwm+nword
! Subprogram not used   if(hwm.lt.nwm) hwm=nwm
! Subprogram not used 
! Subprogram not used end subroutine ci_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: coI0_ - check in as an integer scalar
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine coI0_(marg,thread)
! Subprogram not used       implicit none
! Subprogram not used       integer,intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::coI0_'
! Subprogram not used 
! Subprogram not used   if(mall_on) call co_(1,thread)
! Subprogram not used 
! Subprogram not used end subroutine coI0_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: coI1_ - check in as an integer rank 1 array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine coI1_(marg,thread)
! Subprogram not used       implicit none
! Subprogram not used       integer,dimension(:),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::coI1_'
! Subprogram not used 
! Subprogram not used   if(mall_on) call co_(size(marg),thread)
! Subprogram not used 
! Subprogram not used end subroutine coI1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: coI2_ - check in as an integer rank 2 array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine coI2_(marg,thread)
! Subprogram not used       implicit none
! Subprogram not used       integer,dimension(:,:),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::coI2_'
! Subprogram not used 
! Subprogram not used   if(mall_on) call co_(size(marg),thread)
! Subprogram not used 
! Subprogram not used end subroutine coI2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: coI3_ - check in as an integer rank 3 array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine coI3_(marg,thread)
! Subprogram not used       implicit none
! Subprogram not used       integer,dimension(:,:,:),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::coI3_'
! Subprogram not used 
! Subprogram not used   if(mall_on) call co_(size(marg),thread)
! Subprogram not used 
! Subprogram not used end subroutine coI3_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: coR0_ - check in as a real(SP) scalar
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine coR0_(marg,thread)
! Subprogram not used       use m_realkinds, only : SP
! Subprogram not used       implicit none
! Subprogram not used       real(SP),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::coR0_'
! Subprogram not used 
! Subprogram not used   if(mall_on) call co_(1,thread)
! Subprogram not used 
! Subprogram not used end subroutine coR0_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: coR1_ - check in as a real(SP) rank 1 array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine coR1_(marg,thread)
! Subprogram not used       use m_realkinds, only : SP
! Subprogram not used       implicit none
! Subprogram not used       real(SP),dimension(:),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::coR1_'
! Subprogram not used 
! Subprogram not used   if(mall_on) call co_(size(marg),thread)
! Subprogram not used 
! Subprogram not used end subroutine coR1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: coR2_ - check in as a real(SP) rank 2 array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine coR2_(marg,thread)
! Subprogram not used       use m_realkinds, only : SP
! Subprogram not used       implicit none
! Subprogram not used       real(SP),dimension(:,:),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::coR2_'
! Subprogram not used 
! Subprogram not used   if(mall_on) call co_(size(marg),thread)
! Subprogram not used 
! Subprogram not used end subroutine coR2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: coR3_ - check in as a real(SP) rank 3 array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine coR3_(marg,thread)
! Subprogram not used       use m_realkinds, only : SP
! Subprogram not used       implicit none
! Subprogram not used       real(SP),dimension(:,:,:),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::coR3_'
! Subprogram not used 
! Subprogram not used   if(mall_on) call co_(size(marg),thread)
! Subprogram not used 
! Subprogram not used end subroutine coR3_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: coD0_ - check in as a real(DP) scalar
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine coD0_(marg,thread)
! Subprogram not used       use m_realkinds, only : DP
! Subprogram not used       implicit none
! Subprogram not used       real(DP),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::coD0_'
! Subprogram not used 
! Subprogram not used   if(mall_on) call co_(2,thread)
! Subprogram not used 
! Subprogram not used end subroutine coD0_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: coD1_ - check in as a real(DP) rank 1 array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine coD1_(marg,thread)
! Subprogram not used       use m_realkinds, only : DP
! Subprogram not used       implicit none
! Subprogram not used       real(DP),dimension(:),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::coD1_'
! Subprogram not used 
! Subprogram not used   if(mall_on) call co_(2*size(marg),thread)
! Subprogram not used 
! Subprogram not used end subroutine coD1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: coD2_ - check in as a real(DP) rank 2 array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine coD2_(marg,thread)
! Subprogram not used       use m_realkinds, only : DP
! Subprogram not used       implicit none
! Subprogram not used       real(DP),dimension(:,:),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::coD2_'
! Subprogram not used 
! Subprogram not used   if(mall_on) call co_(2*size(marg),thread)
! Subprogram not used 
! Subprogram not used end subroutine coD2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: coD3_ - check in as a real(DP) rank 3 array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine coD3_(marg,thread)
! Subprogram not used       use m_realkinds, only : DP
! Subprogram not used       implicit none
! Subprogram not used       real(DP),dimension(:,:,:),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::coD3_'
! Subprogram not used 
! Subprogram not used   if(mall_on) call co_(2*size(marg),thread)
! Subprogram not used 
! Subprogram not used end subroutine coD3_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: coL0_ - check in as a logical scalar
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine coL0_(marg,thread)
! Subprogram not used       implicit none
! Subprogram not used       logical,intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::coL0_'
! Subprogram not used 
! Subprogram not used   if(mall_on) call co_(1,thread)
! Subprogram not used 
! Subprogram not used end subroutine coL0_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: coL1_ - check in as a logical rank 1 array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine coL1_(marg,thread)
! Subprogram not used       implicit none
! Subprogram not used       logical,dimension(:),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::coL1_'
! Subprogram not used 
! Subprogram not used   if(mall_on) call co_(size(marg),thread)
! Subprogram not used 
! Subprogram not used end subroutine coL1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: coL2_ - check in as a logical rank 2 array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine coL2_(marg,thread)
! Subprogram not used       implicit none
! Subprogram not used       logical,dimension(:,:),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::coL2_'
! Subprogram not used 
! Subprogram not used   if(mall_on) call co_(size(marg),thread)
! Subprogram not used 
! Subprogram not used end subroutine coL2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: coL3_ - check in as a logical rank 3 array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine coL3_(marg,thread)
! Subprogram not used       implicit none
! Subprogram not used       logical,dimension(:,:,:),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::coL3_'
! Subprogram not used 
! Subprogram not used   if(mall_on) call co_(size(marg),thread)
! Subprogram not used 
! Subprogram not used end subroutine coL3_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: coC0_ - check in as a character scalar
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine coC0_(marg,thread)
! Subprogram not used       implicit none
! Subprogram not used       character(len=*),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::coC0_'
! Subprogram not used   integer :: nw
! Subprogram not used 
! Subprogram not used   if(.not.mall_on) return
! Subprogram not used   nw=(len(marg)+NBYTE_PER_WORD-1)/NBYTE_PER_WORD
! Subprogram not used   call co_(nw,thread)
! Subprogram not used 
! Subprogram not used end subroutine coC0_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: coC1_ - check in as a character rank 1 array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine coC1_(marg,thread)
! Subprogram not used       implicit none
! Subprogram not used       character(len=*),dimension(:),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::coC1_'
! Subprogram not used   integer :: nw
! Subprogram not used 
! Subprogram not used   if(.not.mall_on) return
! Subprogram not used   nw=(len(marg(1))+NBYTE_PER_WORD-1)/NBYTE_PER_WORD
! Subprogram not used   call co_(size(marg)*nw,thread)
! Subprogram not used 
! Subprogram not used end subroutine coC1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: coC2_ - check in as a character rank 2 array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine coC2_(marg,thread)
! Subprogram not used       implicit none
! Subprogram not used       character(len=*),dimension(:,:),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::coC2_'
! Subprogram not used   integer :: nw
! Subprogram not used 
! Subprogram not used   if(.not.mall_on) return
! Subprogram not used   nw=(len(marg(1,1))+NBYTE_PER_WORD-1)/NBYTE_PER_WORD
! Subprogram not used   call co_(size(marg)*nw,thread)
! Subprogram not used 
! Subprogram not used end subroutine coC2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: coC3_ - check in as a character rank 3 array
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine coC3_(marg,thread)
! Subprogram not used       implicit none
! Subprogram not used       character(len=*),dimension(:,:,:),intent(in) :: marg
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	21Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
! Subprogram not used !		- initial prototype/prolog/code
! Subprogram not used !EOP ___________________________________________________________________
! Subprogram not used 
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::coC3_'
! Subprogram not used   integer :: nw
! Subprogram not used 
! Subprogram not used   if(.not.mall_on) return
! Subprogram not used   nw=(len(marg(1,1,1))+NBYTE_PER_WORD-1)/NBYTE_PER_WORD
! Subprogram not used   call co_(size(marg)*nw,thread)
! Subprogram not used 
! Subprogram not used end subroutine coC3_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: co_ - check-out allocate activity
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine co_(nword,thread)
! Subprogram not used       use m_stdio, only : stderr
! Subprogram not used       use m_die, only : die
! Subprogram not used       implicit none
! Subprogram not used       integer,intent(in) :: nword
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	13Feb98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! Subprogram not used !EOP
! Subprogram not used !_______________________________________________________________________
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::co_'
! Subprogram not used   integer :: ith
! Subprogram not used 
! Subprogram not used   if(.not.mall_on) return
! Subprogram not used 
! Subprogram not used   if(nword < 0) then
! Subprogram not used     write(stderr,'(2a,i4)') myname_,	&
! Subprogram not used 	': invalide argument, nword = ',nword
! Subprogram not used     call die(myname_)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used 	! if the thread is "unknown", it would be treated as a
! Subprogram not used 	! new thread with net negative memory activity.
! Subprogram not used 
! Subprogram not used   ith=lookup_(thread)
! Subprogram not used 
! Subprogram not used 	! update the account
! Subprogram not used 
! Subprogram not used   nci_(ith)=nci_(ith)-1
! Subprogram not used   nwm_(ith)=nwm_(ith)-nword
! Subprogram not used 
! Subprogram not used 	! update the total budget
! Subprogram not used 
! Subprogram not used   nci=nci-1
! Subprogram not used   nwm=nwm-nword
! Subprogram not used 
! Subprogram not used end subroutine co_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: cix_ - handling macro ALLOC_() error
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine cix_(thread,stat,fnam,line)
! Subprogram not used       use m_stdio, only : stderr
! Subprogram not used       use m_die, only : die
! Subprogram not used       implicit none
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used       integer,intent(in) :: stat
! Subprogram not used       character(len=*),intent(in) :: fnam
! Subprogram not used       integer,intent(in) :: line
! Subprogram not used 
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	13Feb98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! Subprogram not used !EOP
! Subprogram not used !_______________________________________________________________________
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::cix_'
! Subprogram not used 
! Subprogram not used   write(stderr,'(2a,i4)') trim(thread),	&
! Subprogram not used 	': ALLOC_() error, stat =',stat
! Subprogram not used   call die('ALLOC_',fnam,line)
! Subprogram not used 
! Subprogram not used end subroutine cix_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: cox_ - handling macro DEALLOC_() error
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine cox_(thread,stat,fnam,line)
! Subprogram not used       use m_stdio, only : stderr
! Subprogram not used       use m_die, only : die
! Subprogram not used       implicit none
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used       integer,intent(in) :: stat
! Subprogram not used       character(len=*),intent(in) :: fnam
! Subprogram not used       integer,intent(in) :: line
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	13Feb98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! Subprogram not used !EOP
! Subprogram not used !_______________________________________________________________________
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::cox_'
! Subprogram not used 
! Subprogram not used   write(stderr,'(2a,i4)') trim(thread),	&
! Subprogram not used 	': DEALLOC_() error, stat =',stat
! Subprogram not used   call die('DEALLOC_',fnam,line)
! Subprogram not used 
! Subprogram not used end subroutine cox_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: flush_ - balancing the up-to-date ci/co calls
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine flush_(lu)
! Subprogram not used       use m_stdio, only : stderr
! Subprogram not used       use m_ioutil, only : luflush
! Subprogram not used       use m_die, only : die
! Subprogram not used       implicit none
! Subprogram not used       integer,intent(in) :: lu
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	17Feb98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! Subprogram not used !EOP
! Subprogram not used !_______________________________________________________________________
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::flush_'
! Subprogram not used 
! Subprogram not used   integer,parameter :: lnmax=38
! Subprogram not used   character(len=max(lnmax,NSZ)) :: name
! Subprogram not used 
! Subprogram not used   character(len=6) :: hwm_wd,nwm_wd
! Subprogram not used   character(len=1) :: flag_ci,flag_wm
! Subprogram not used   integer :: i,ier,ln
! Subprogram not used 
! Subprogram not used   if(.not.mall_on) return
! Subprogram not used 
! Subprogram not used   if(.not.started) call reset_()
! Subprogram not used 
! Subprogram not used   write(lu,'(72a/)',iostat=ier) ('_',i=1,72)
! Subprogram not used   if(ier /= 0) then
! Subprogram not used     write(stderr,'(2a,i3)') myname_,': can not write(), unit =',lu
! Subprogram not used     call die(myname_)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   write(lu,'(a,t39,4(2x,a))',iostat=ier) '[MALL]',	&
! Subprogram not used   		'max-ci','net-ci ','max-wm','net-wm'
! Subprogram not used   if(ier /= 0) then
! Subprogram not used     write(stderr,'(2a,i4)') myname_,': can not write(), unit =',lu
! Subprogram not used     call die(myname_)
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   call luflush(lu)
! Subprogram not used 
! Subprogram not used !23.|....1....|....2....|....3....|....4....|....5....|....6....|....7..
! Subprogram not used !_______________________________________________________________________
! Subprogram not used !
! Subprogram not used ![MALL]                                 max_ci  net-ci   max-wm  net-wm
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used !total.                                 ...333  ...333*  ..333M  ..333i*
! Subprogram not used !_______________________________________________________________________
! Subprogram not used 
! Subprogram not used   write(lu,'(72a)') ('-',i=1,72)
! Subprogram not used 
! Subprogram not used   do i=1,min(n_,MXL)
! Subprogram not used     call wcount_(hwm_(i),hwm_wd)
! Subprogram not used     call wcount_(nwm_(i),nwm_wd)
! Subprogram not used       
! Subprogram not used     flag_ci=' '
! Subprogram not used     if(nci_(i) /= 0) flag_ci='*'
! Subprogram not used 
! Subprogram not used     flag_wm=' '
! Subprogram not used     if(nwm_(i) /= 0) flag_wm='*'
! Subprogram not used 
! Subprogram not used     name=name_(i)
! Subprogram not used     ln=max(len_trim(name),lnmax)
! Subprogram not used     write(lu,'(a,2(2x,i6),a,2(2x,a6),a)') name(1:ln),	&
! Subprogram not used 	mci_(i),nci_(i),flag_ci,hwm_wd,nwm_wd,flag_wm
! Subprogram not used   end do
! Subprogram not used 
! Subprogram not used   call wcount_(hwm,hwm_wd)
! Subprogram not used   call wcount_(nwm,nwm_wd)
! Subprogram not used       
! Subprogram not used   flag_ci=' '
! Subprogram not used   if(nci /= 0) flag_ci='*'
! Subprogram not used   flag_wm=' '
! Subprogram not used   if(nwm /= 0) flag_wm='*'
! Subprogram not used 
! Subprogram not used   name='.total.'
! Subprogram not used   ln=max(len_trim(name),lnmax)
! Subprogram not used   write(lu,'(a,2(2x,i6),a,2(2x,a6),a)') name(1:ln),	&
! Subprogram not used 	mci,nci,flag_ci,hwm_wd,nwm_wd,flag_wm
! Subprogram not used 
! Subprogram not used   write(lu,'(72a/)') ('_',i=1,72)
! Subprogram not used 
! Subprogram not used   if(nreset /= 1) write(lu,'(2a,i3,a)') myname_,	&
! Subprogram not used 	': reset_ ',nreset,' times'
! Subprogram not used 
! Subprogram not used   call luflush(lu)
! Subprogram not used end subroutine flush_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: wcount_ - generate word count output with unit
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine wcount_(wknt,cknt)
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used       integer,         intent(in)  :: wknt ! given an integer value
! Subprogram not used       character(len=6),intent(out) :: cknt ! return a string value
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	17Feb98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! Subprogram not used !EOP
! Subprogram not used !_______________________________________________________________________
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::wcount_'
! Subprogram not used 
! Subprogram not used character(len=1) :: cwd
! Subprogram not used integer,parameter :: KWD=1024
! Subprogram not used integer,parameter :: MWD=1024*1024
! Subprogram not used integer,parameter :: GWD=1024*1024*1024
! Subprogram not used 
! Subprogram not used integer :: iwd
! Subprogram not used 
! Subprogram not used if(wknt < 0) then
! Subprogram not used   cknt='------'
! Subprogram not used else
! Subprogram not used   cwd='i'
! Subprogram not used   iwd=wknt
! Subprogram not used   if(iwd > 9999) then
! Subprogram not used     cwd='K'
! Subprogram not used     iwd=(wknt+KWD-1)/KWD
! Subprogram not used   endif
! Subprogram not used   if(iwd > 9999) then
! Subprogram not used     cwd='M'
! Subprogram not used     iwd=(wknt+MWD-1)/MWD
! Subprogram not used   endif
! Subprogram not used   if(iwd > 9999) then
! Subprogram not used     cwd='G'
! Subprogram not used     iwd=(wknt+GWD-1)/GWD
! Subprogram not used   endif
! Subprogram not used   write(cknt,'(i5,a)') iwd,cwd
! Subprogram not used endif
! Subprogram not used end subroutine wcount_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: lookup_ - search/insert a name in a list
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     function lookup_(thread)
! Subprogram not used       use m_chars, only : uppercase
! Subprogram not used       implicit none
! Subprogram not used       character(len=*),intent(in) :: thread
! Subprogram not used       integer :: lookup_
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	17Feb98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! Subprogram not used !EOP
! Subprogram not used !_______________________________________________________________________
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::lookup_'
! Subprogram not used 
! Subprogram not used   logical :: found
! Subprogram not used   integer :: ith
! Subprogram not used 
! Subprogram not used   if(.not.started) call reset_()
! Subprogram not used 
! Subprogram not used !----------------------------------------
! Subprogram not used ith=0
! Subprogram not used found=.false.
! Subprogram not used do while(.not.found .and. ith < min(n_,MXL))
! Subprogram not used   ith=ith+1
! Subprogram not used   found= uppercase(thread) == uppercase(name_(ith))
! Subprogram not used end do
! Subprogram not used 
! Subprogram not used if(.not.found) then
! Subprogram not used   if(n_==0) then
! Subprogram not used     nci=0
! Subprogram not used     mci=0
! Subprogram not used     nwm=0
! Subprogram not used     hwm=0
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   n_=n_+1
! Subprogram not used   if(n_ == MXL) then
! Subprogram not used     ith=MXL
! Subprogram not used     name_(ith)='.overflow.'
! Subprogram not used   else
! Subprogram not used     ith=n_
! Subprogram not used     name_(ith)=thread
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used   nci_(ith)=0
! Subprogram not used   mci_(ith)=0
! Subprogram not used   nwm_(ith)=0
! Subprogram not used   hwm_(ith)=0
! Subprogram not used endif
! Subprogram not used 
! Subprogram not used lookup_=ith
! Subprogram not used 
! Subprogram not used end function lookup_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: reset_ - initialize the module data structure
!
! !DESCRIPTION:
!
! !INTERFACE:

! Subprogram not used     subroutine reset_()
! Subprogram not used       implicit none
! Subprogram not used 
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! 	16Mar98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! Subprogram not used !EOP
! Subprogram not used !_______________________________________________________________________
! Subprogram not used   character(len=*),parameter :: myname_=myname//'::reset_'
! Subprogram not used 
! Subprogram not used   if(.not.mall_on) return
! Subprogram not used 
! Subprogram not used   nreset=nreset+1
! Subprogram not used   started=.true.
! Subprogram not used 
! Subprogram not used   name_(1:n_)=' '
! Subprogram not used 
! Subprogram not used   mci_(1:n_)=0
! Subprogram not used   nci_(1:n_)=0
! Subprogram not used   hwm_(1:n_)=0
! Subprogram not used   nwm_(1:n_)=0
! Subprogram not used 
! Subprogram not used   n_ =0
! Subprogram not used 
! Subprogram not used   mci=0
! Subprogram not used   nci=0
! Subprogram not used   hwm=0
! Subprogram not used   nwm=0
! Subprogram not used 
! Subprogram not used end subroutine reset_
!=======================================================================
end module m_mall
