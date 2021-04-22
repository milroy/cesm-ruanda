
module pspect
!----------------------------------------------------------------------- 
! 
! Purpose: Parameters related to spectral domain
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
   integer, parameter :: ptrm = 1                ! M truncation parameter
   integer, parameter :: ptrn = 1                ! N truncation parameter
   integer, parameter :: ptrk = 1                ! K truncation parameter
                                                   
   integer, parameter :: pmax = ptrn+1              ! number of diagonals
   integer, parameter :: pmaxp = pmax+1             ! Number of diagonals plus 1
   integer, parameter :: pnmax = ptrk+1             ! Number of values of N
   integer, parameter :: pmmax = ptrm+1             ! Number of values of M
   integer, parameter :: par0 = ptrm+ptrn-ptrk      ! intermediate parameter
   integer, parameter :: par2 = par0*(par0+1)/2     ! intermediate parameter
   integer, parameter :: pspt = (ptrn+1)*pmmax-par2 ! Total num complex spectral coeffs retained
   integer, parameter :: psp = 2*pspt               ! 2*pspt (real) size of coeff array per level

end module pspect
