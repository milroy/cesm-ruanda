module subgridAveMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: subgridAveMod
!
! !DESCRIPTION:
! Utilities to perfrom subgrid averaging
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clmtype
  use clm_varcon, only : spval, isturb,  icol_roof, icol_sunwall, icol_shadewall, &
                         icol_road_perv, icol_road_imperv
  use clm_varctl, only : iulog
  use abortutils, only : endrun

! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: p2c   ! Perfrom an average from pfts to columns
  public :: p2l   ! Perfrom an average from pfts to landunits
  public :: p2g   ! Perfrom an average from pfts to gridcells
  public :: c2l   ! Perfrom an average from columns to landunits
  public :: c2g   ! Perfrom an average from columns to gridcells
  public :: l2g   ! Perfrom an average from landunits to gridcells

  interface p2c
     module procedure p2c_1d
     module procedure p2c_2d
     module procedure p2c_1d_filter
     module procedure p2c_2d_filter
  end interface
  interface p2l
     module procedure p2l_1d
     module procedure p2l_2d
  end interface
  interface p2g
     module procedure p2g_1d
     module procedure p2g_2d
  end interface
  interface c2l
     module procedure c2l_1d
     module procedure c2l_2d
  end interface
  interface c2g
     module procedure c2g_1d
     module procedure c2g_2d
  end interface
  interface l2g
     module procedure l2g_1d
     module procedure l2g_2d
  end interface
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: build_scale_l2g
  private :: create_scale_l2g_lookup
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/03
!
!EOP

! WJS (10-14-11): TODO:
! 
! - I believe that scale_p2c, scale_c2l and scale_l2g should be included in the sumwt
! accumulations (e.g., sumwt = sumwt + wtgcell * scale_p2c * scale_c2l * scale_l2g), but
! that requires some more thought to (1) make sure that is correct, and (2) make sure it
! doesn't break the urban scaling. (See also my notes in create_scale_l2g_lookup.)
!   - Once that is done, you could use a scale of 0, avoiding the need for the use of
!   spval and the special checks that requires.
!
! - Currently, there is a lot of repeated code to calculate scale_c2l. This should be
! cleaned up.
!   - At a minimum, should collect the repeated code into a subroutine to eliminate this
!   repitition
!   - The best thing might be to use a lookup array, as is done for scale_l2g
! -----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: p2c_1d
!
! !INTERFACE:
! Subprogram not used   subroutine p2c_1d (lbp, ubp, lbc, ubc, parr, carr, p2c_scale_type)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used ! Perfrom subgrid-average from pfts to columns.
! Subprogram not used ! Averaging is only done for points that are not equal to "spval".
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used     use clm_varpar, only : max_pft_per_col
! Subprogram not used !
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     integer , intent(in)  :: lbp, ubp              ! beginning and ending pft
! Subprogram not used     integer , intent(in)  :: lbc, ubc              ! beginning and ending column
! Subprogram not used     real(r8), intent(in)  :: parr(lbp:ubp)         ! pft array
! Subprogram not used     real(r8), intent(out) :: carr(lbc:ubc)         ! column array
! Subprogram not used     character(len=*), intent(in) :: p2c_scale_type ! scale type
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! Created by Mariana Vertenstein 12/03
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! !LOCAL VARIABLES:
! Subprogram not used !EOP
! Subprogram not used     integer  :: pi,p,c,index           ! indices
! Subprogram not used     real(r8) :: scale_p2c(lbp:ubp)     ! scale factor for column->landunit mapping
! Subprogram not used     logical  :: found                  ! temporary for error check
! Subprogram not used     real(r8) :: sumwt(lbc:ubc)         ! sum of weights
! Subprogram not used     real(r8), pointer :: wtcol(:)      ! weight of pft relative to column
! Subprogram not used     integer , pointer :: pcolumn(:)    ! column index of corresponding pft
! Subprogram not used     integer , pointer :: npfts(:)      ! number of pfts in column
! Subprogram not used     integer , pointer :: pfti(:)       ! initial pft index in column
! Subprogram not used !------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     wtcol    => pft%wtcol
! Subprogram not used     pcolumn  => pft%column
! Subprogram not used     npfts    => col%npfts
! Subprogram not used     pfti     => col%pfti
! Subprogram not used 
! Subprogram not used     if (p2c_scale_type == 'unity') then
! Subprogram not used        do p = lbp,ubp
! Subprogram not used           scale_p2c(p) = 1.0_r8
! Subprogram not used        end do
! Subprogram not used     else
! Subprogram not used        write(iulog,*)'p2c_1d error: scale type ',p2c_scale_type,' not supported'
! Subprogram not used        call endrun()
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     carr(lbc:ubc) = spval
! Subprogram not used     sumwt(lbc:ubc) = 0._r8
! Subprogram not used     do p = lbp,ubp
! Subprogram not used        if (wtcol(p) /= 0._r8) then
! Subprogram not used           if (parr(p) /= spval) then
! Subprogram not used              c = pcolumn(p)
! Subprogram not used              if (sumwt(c) == 0._r8) carr(c) = 0._r8
! Subprogram not used              carr(c) = carr(c) + parr(p) * scale_p2c(p) * wtcol(p)
! Subprogram not used              sumwt(c) = sumwt(c) + wtcol(p)
! Subprogram not used           end if
! Subprogram not used        end if
! Subprogram not used     end do
! Subprogram not used     found = .false.
! Subprogram not used     do c = lbc,ubc
! Subprogram not used        if (sumwt(c) > 1.0_r8 + 1.e-6_r8) then
! Subprogram not used           found = .true.
! Subprogram not used           index = c
! Subprogram not used        else if (sumwt(c) /= 0._r8) then
! Subprogram not used           carr(c) = carr(c)/sumwt(c)
! Subprogram not used        end if
! Subprogram not used     end do
! Subprogram not used     if (found) then
! Subprogram not used        write(iulog,*)'p2c error: sumwt is greater than 1.0 at c= ',index
! Subprogram not used        call endrun()
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used   end subroutine p2c_1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: p2c_2d
!
! !INTERFACE:
! Subprogram not used   subroutine p2c_2d (lbp, ubp, lbc, ubc, num2d, parr, carr, p2c_scale_type)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used ! Perfrom subgrid-average from landunits to gridcells.
! Subprogram not used ! Averaging is only done for points that are not equal to "spval".
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used     use clm_varpar, only : max_pft_per_col
! Subprogram not used !
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     integer , intent(in)  :: lbp, ubp              ! beginning and ending pft
! Subprogram not used     integer , intent(in)  :: lbc, ubc              ! beginning and ending column
! Subprogram not used     integer , intent(in)  :: num2d                 ! size of second dimension
! Subprogram not used     real(r8), intent(in)  :: parr(lbp:ubp,num2d)   ! pft array
! Subprogram not used     real(r8), intent(out) :: carr(lbc:ubc,num2d)   ! column array
! Subprogram not used     character(len=*), intent(in) :: p2c_scale_type ! scale type
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! Created by Mariana Vertenstein 12/03
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! !LOCAL VARIABLES:
! Subprogram not used !EOP
! Subprogram not used     integer  :: j,pi,p,c,index         ! indices
! Subprogram not used     real(r8) :: scale_p2c(lbp:ubp)     ! scale factor for column->landunit mapping
! Subprogram not used     logical  :: found                  ! temporary for error check
! Subprogram not used     real(r8) :: sumwt(lbc:ubc)         ! sum of weights
! Subprogram not used     real(r8), pointer :: wtcol(:)      ! weight of pft relative to column
! Subprogram not used     integer , pointer :: pcolumn(:)    ! column index of corresponding pft
! Subprogram not used     integer , pointer :: npfts(:)      ! number of pfts in column
! Subprogram not used     integer , pointer :: pfti(:)       ! initial pft index in column
! Subprogram not used !------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     wtcol    => pft%wtcol
! Subprogram not used     pcolumn  => pft%column
! Subprogram not used     npfts    => col%npfts
! Subprogram not used     pfti     => col%pfti
! Subprogram not used 
! Subprogram not used     if (p2c_scale_type == 'unity') then
! Subprogram not used        do p = lbp,ubp
! Subprogram not used           scale_p2c(p) = 1.0_r8
! Subprogram not used        end do
! Subprogram not used     else
! Subprogram not used        write(iulog,*)'p2c_2d error: scale type ',p2c_scale_type,' not supported'
! Subprogram not used        call endrun()
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     carr(:,:) = spval
! Subprogram not used     do j = 1,num2d
! Subprogram not used        sumwt(:) = 0._r8
! Subprogram not used        do p = lbp,ubp
! Subprogram not used           if (wtcol(p) /= 0._r8) then
! Subprogram not used              if (parr(p,j) /= spval) then
! Subprogram not used                 c = pcolumn(p)
! Subprogram not used                 if (sumwt(c) == 0._r8) carr(c,j) = 0._r8
! Subprogram not used                 carr(c,j) = carr(c,j) + parr(p,j) * scale_p2c(p) * wtcol(p)
! Subprogram not used                 sumwt(c) = sumwt(c) + wtcol(p)
! Subprogram not used              end if
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used        found = .false.
! Subprogram not used        do c = lbc,ubc
! Subprogram not used           if (sumwt(c) > 1.0_r8 + 1.e-6_r8) then
! Subprogram not used              found = .true.
! Subprogram not used              index = c
! Subprogram not used           else if (sumwt(c) /= 0._r8) then
! Subprogram not used              carr(c,j) = carr(c,j)/sumwt(c)
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used        if (found) then
! Subprogram not used           write(iulog,*)'p2c_2d error: sumwt is greater than 1.0 at c= ',index,' lev= ',j
! Subprogram not used           call endrun()
! Subprogram not used        end if
! Subprogram not used     end do 
! Subprogram not used   end subroutine p2c_2d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: p2c_1d_filter
!
! !INTERFACE:
  subroutine p2c_1d_filter (numfc, filterc, pftarr, colarr)
!
! !DESCRIPTION:
! perform pft to column averaging for single level pft arrays
!
! !USES:
    use clm_varpar, only : max_pft_per_col
    use clm_varcon, only : istice_mec
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: numfc
    integer , intent(in)  :: filterc(numfc)
    real(r8), pointer     :: pftarr(:)
    real(r8), pointer     :: colarr(:)
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/03
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: fc,c,pi,p,l         ! indices
    integer , pointer :: npfts(:)
    integer , pointer :: pfti(:)
    integer , pointer :: pftf(:)
    integer , pointer :: clandunit(:)
    integer , pointer :: ltype(:)
    real(r8), pointer :: wtcol(:)
    real(r8), pointer :: wtgcell(:)
!-----------------------------------------------------------------------

    npfts     => col%npfts
    pfti      => col%pfti
    pftf      => col%pftf
    clandunit => col%landunit
    ltype     => lun%itype
    wtcol     => pft%wtcol
    wtgcell   => pft%wtgcell

    do fc = 1,numfc
       c = filterc(fc)
       l = clandunit(c)
       colarr(c) = 0._r8
       do p = pfti(c), pftf(c)
          ! Note: some glacier_mec pfts may have zero weight
          if (wtgcell(p) > 0._r8 .or. ltype(l)==istice_mec) colarr(c) = colarr(c) + pftarr(p) * wtcol(p)
       end do
    end do

  end subroutine p2c_1d_filter

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: p2c_2d_filter
!
! !INTERFACE:
! Subprogram not used   subroutine p2c_2d_filter (lev, numfc, filterc, pftarr, colarr)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used ! perform pft to column averaging for multi level pft arrays
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used     use clm_varpar, only : max_pft_per_col
! Subprogram not used 
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     integer , intent(in)  :: lev
! Subprogram not used     integer , intent(in)  :: numfc
! Subprogram not used     integer , intent(in)  :: filterc(numfc)
! Subprogram not used     real(r8), pointer     :: pftarr(:,:)
! Subprogram not used     real(r8), pointer     :: colarr(:,:)
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! Created by Mariana Vertenstein 12/03
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! !LOCAL VARIABLES:
! Subprogram not used !EOP
! Subprogram not used     integer :: fc,c,pi,p,j    ! indices
! Subprogram not used     integer , pointer :: npfts(:)
! Subprogram not used     integer , pointer :: pfti(:)
! Subprogram not used     integer , pointer :: pftf(:)
! Subprogram not used     real(r8), pointer :: wtcol(:)
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     npfts => col%npfts
! Subprogram not used     pfti  => col%pfti
! Subprogram not used     pftf  => col%pftf
! Subprogram not used     wtcol => pft%wtcol
! Subprogram not used 
! Subprogram not used     do j = 1,lev
! Subprogram not used        do fc = 1,numfc
! Subprogram not used           c = filterc(fc)
! Subprogram not used           colarr(c,j) = 0._r8
! Subprogram not used           do p = pfti(c), pftf(c)
! Subprogram not used              colarr(c,j) = colarr(c,j) + pftarr(p,j) * wtcol(p)
! Subprogram not used           end do
! Subprogram not used        end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used   end subroutine p2c_2d_filter

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: p2l_1d
!
! !INTERFACE:
! Subprogram not used   subroutine p2l_1d (lbp, ubp, lbc, ubc, lbl, ubl, parr, larr, &
! Subprogram not used        p2c_scale_type, c2l_scale_type)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used ! Perfrom subgrid-average from pfts to landunits
! Subprogram not used ! Averaging is only done for points that are not equal to "spval".
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used     use clm_varpar, only : max_pft_per_lu
! Subprogram not used !
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     integer , intent(in)  :: lbp, ubp              ! beginning and ending pft indices
! Subprogram not used     integer , intent(in)  :: lbc, ubc              ! beginning and ending column indices
! Subprogram not used     integer , intent(in)  :: lbl, ubl              ! beginning and ending landunit indices
! Subprogram not used     real(r8), intent(in)  :: parr(lbp:ubp)         ! input column array
! Subprogram not used     real(r8), intent(out) :: larr(lbl:ubl)         ! output landunit array
! Subprogram not used     character(len=*), intent(in) :: p2c_scale_type ! scale factor type for averaging
! Subprogram not used     character(len=*), intent(in) :: c2l_scale_type ! scale factor type for averaging
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! Created by Mariana Vertenstein 12/03
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! !LOCAL VARIABLES:
! Subprogram not used !EOP
! Subprogram not used     integer  :: pi,p,c,l,index         ! indices
! Subprogram not used     logical  :: found                  ! temporary for error check
! Subprogram not used     real(r8) :: sumwt(lbl:ubl)         ! sum of weights
! Subprogram not used     real(r8) :: scale_p2c(lbc:ubc)     ! scale factor for pft->column mapping
! Subprogram not used     real(r8) :: scale_c2l(lbc:ubc)     ! scale factor for column->landunit mapping
! Subprogram not used     real(r8), pointer :: wtlunit(:)    ! weight of pft relative to landunit
! Subprogram not used     integer , pointer :: pcolumn(:)    ! column of corresponding pft
! Subprogram not used     integer , pointer :: plandunit(:)  ! landunit of corresponding pft
! Subprogram not used     integer , pointer :: npfts(:)      ! number of pfts in landunit
! Subprogram not used     integer , pointer :: pfti(:)       ! initial pft index in landunit
! Subprogram not used     integer , pointer :: clandunit(:)  ! landunit of corresponding column
! Subprogram not used     integer , pointer :: ctype(:)      ! column type
! Subprogram not used     integer , pointer :: ltype(:)      ! landunit type
! Subprogram not used     real(r8), pointer :: canyon_hwr(:) ! urban canyon height to width ratio
! Subprogram not used !------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     canyon_hwr => lun%canyon_hwr
! Subprogram not used     ltype      => lun%itype
! Subprogram not used     ctype      => col%itype
! Subprogram not used     clandunit  => col%landunit
! Subprogram not used     wtlunit    => pft%wtlunit
! Subprogram not used     pcolumn    => pft%column
! Subprogram not used     plandunit  => pft%landunit
! Subprogram not used     npfts      => lun%npfts
! Subprogram not used     pfti       => lun%pfti
! Subprogram not used 
! Subprogram not used     if (c2l_scale_type == 'unity') then
! Subprogram not used        do c = lbc,ubc
! Subprogram not used           scale_c2l(c) = 1.0_r8
! Subprogram not used        end do
! Subprogram not used     else if (c2l_scale_type == 'urbanf') then
! Subprogram not used        do c = lbc,ubc
! Subprogram not used           l = clandunit(c) 
! Subprogram not used           if (ltype(l) == isturb) then
! Subprogram not used              if (ctype(c) == icol_sunwall) then
! Subprogram not used                 scale_c2l(c) = 3.0 * canyon_hwr(l) 
! Subprogram not used              else if (ctype(c) == icol_shadewall) then
! Subprogram not used                 scale_c2l(c) = 3.0 * canyon_hwr(l) 
! Subprogram not used              else if (ctype(c) == icol_road_perv .or. ctype(c) == icol_road_imperv) then
! Subprogram not used                 scale_c2l(c) = 3.0_r8
! Subprogram not used              else if (ctype(c) == icol_roof) then
! Subprogram not used                 scale_c2l(c) = 1.0_r8
! Subprogram not used              end if
! Subprogram not used           else
! Subprogram not used              scale_c2l(c) = 1.0_r8
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used     else if (c2l_scale_type == 'urbans') then
! Subprogram not used        do c = lbc,ubc
! Subprogram not used           l = clandunit(c) 
! Subprogram not used           if (ltype(l) == isturb) then
! Subprogram not used              if (ctype(c) == icol_sunwall) then
! Subprogram not used                 scale_c2l(c) = (3.0 * canyon_hwr(l)) / (2.*canyon_hwr(l) + 1.)
! Subprogram not used              else if (ctype(c) == icol_shadewall) then
! Subprogram not used                 scale_c2l(c) = (3.0 * canyon_hwr(l)) / (2.*canyon_hwr(l) + 1.)
! Subprogram not used              else if (ctype(c) == icol_road_perv .or. ctype(c) == icol_road_imperv) then
! Subprogram not used                 scale_c2l(c) = 3.0 / (2.*canyon_hwr(l) + 1.)
! Subprogram not used              else if (ctype(c) == icol_roof) then
! Subprogram not used                 scale_c2l(c) = 1.0_r8
! Subprogram not used              end if
! Subprogram not used           else
! Subprogram not used              scale_c2l(c) = 1.0_r8
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used     else
! Subprogram not used        write(iulog,*)'p2l_1d error: scale type ',c2l_scale_type,' not supported'
! Subprogram not used        call endrun()
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     if (p2c_scale_type == 'unity') then
! Subprogram not used        do p = lbp,ubp
! Subprogram not used           scale_p2c(p) = 1.0_r8
! Subprogram not used        end do
! Subprogram not used     else
! Subprogram not used        write(iulog,*)'p2l_1d error: scale type ',p2c_scale_type,' not supported'
! Subprogram not used        call endrun()
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     larr(:) = spval
! Subprogram not used     sumwt(:) = 0._r8
! Subprogram not used     do p = lbp,ubp
! Subprogram not used        if (wtlunit(p) /= 0._r8) then
! Subprogram not used           c = pcolumn(p)
! Subprogram not used           if (parr(p) /= spval .and. scale_c2l(c) /= spval) then
! Subprogram not used              l = plandunit(p)
! Subprogram not used              if (sumwt(l) == 0._r8) larr(l) = 0._r8
! Subprogram not used              larr(l) = larr(l) + parr(p) * scale_p2c(p) * scale_c2l(c) * wtlunit(p)
! Subprogram not used              sumwt(l) = sumwt(l) + wtlunit(p)
! Subprogram not used           end if
! Subprogram not used        end if
! Subprogram not used     end do
! Subprogram not used     found = .false.
! Subprogram not used     do l = lbl,ubl
! Subprogram not used        if (sumwt(l) > 1.0_r8 + 1.e-6_r8) then
! Subprogram not used           found = .true.
! Subprogram not used           index = l
! Subprogram not used        else if (sumwt(l) /= 0._r8) then
! Subprogram not used           larr(l) = larr(l)/sumwt(l)
! Subprogram not used        end if
! Subprogram not used     end do
! Subprogram not used     if (found) then
! Subprogram not used        write(iulog,*)'p2l_1d error: sumwt is greater than 1.0 at l= ',index
! Subprogram not used        call endrun()
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used   end subroutine p2l_1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: p2l_2d
!
! !INTERFACE:
! Subprogram not used   subroutine p2l_2d(lbp, ubp, lbc, ubc, lbl, ubl, num2d, parr, larr, &
! Subprogram not used        p2c_scale_type, c2l_scale_type)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used ! Perfrom subgrid-average from pfts to landunits
! Subprogram not used ! Averaging is only done for points that are not equal to "spval".
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used     use clm_varpar, only : max_pft_per_lu
! Subprogram not used !
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     integer , intent(in)  :: lbp, ubp              ! beginning and ending pft indices
! Subprogram not used     integer , intent(in)  :: lbc, ubc              ! beginning and ending column indices
! Subprogram not used     integer , intent(in)  :: lbl, ubl              ! beginning and ending landunit indices
! Subprogram not used     integer , intent(in)  :: num2d                 ! size of second dimension
! Subprogram not used     real(r8), intent(in)  :: parr(lbp:ubp,num2d)   ! input pft array
! Subprogram not used     real(r8), intent(out) :: larr(lbl:ubl,num2d)   ! output gridcell array
! Subprogram not used     character(len=*), intent(in) :: p2c_scale_type ! scale factor type for averaging
! Subprogram not used     character(len=*), intent(in) :: c2l_scale_type ! scale factor type for averaging
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! Created by Mariana Vertenstein 12/03
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! !LOCAL VARIABLES:
! Subprogram not used !EOP
! Subprogram not used     integer  :: j,pi,p,c,l,index       ! indices
! Subprogram not used     logical  :: found                  ! temporary for error check
! Subprogram not used     real(r8) :: sumwt(lbl:ubl)         ! sum of weights
! Subprogram not used     real(r8) :: scale_p2c(lbc:ubc)     ! scale factor for pft->column mapping
! Subprogram not used     real(r8) :: scale_c2l(lbc:ubc)     ! scale factor for column->landunit mapping
! Subprogram not used     real(r8), pointer :: wtlunit(:)    ! weight of pft relative to landunit
! Subprogram not used     integer , pointer :: pcolumn(:)    ! column of corresponding pft
! Subprogram not used     integer , pointer :: plandunit(:)  ! landunit of corresponding pft
! Subprogram not used     integer , pointer :: npfts(:)      ! number of pfts in landunit
! Subprogram not used     integer , pointer :: pfti(:)       ! initial pft index in landunit
! Subprogram not used     integer , pointer :: clandunit(:)  ! landunit of corresponding column
! Subprogram not used     integer , pointer :: ctype(:)      ! column type
! Subprogram not used     integer , pointer :: ltype(:)      ! landunit type
! Subprogram not used     real(r8), pointer :: canyon_hwr(:) ! urban canyon height to width ratio
! Subprogram not used !------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     canyon_hwr => lun%canyon_hwr
! Subprogram not used     ltype      => lun%itype
! Subprogram not used     clandunit  => col%landunit
! Subprogram not used     ctype      => col%itype
! Subprogram not used     wtlunit   => pft%wtlunit
! Subprogram not used     pcolumn   => pft%column
! Subprogram not used     plandunit => pft%landunit
! Subprogram not used     npfts     => lun%npfts
! Subprogram not used     pfti      => lun%pfti
! Subprogram not used 
! Subprogram not used     if (c2l_scale_type == 'unity') then
! Subprogram not used        do c = lbc,ubc
! Subprogram not used           scale_c2l(c) = 1.0_r8
! Subprogram not used        end do
! Subprogram not used     else if (c2l_scale_type == 'urbanf') then
! Subprogram not used        do c = lbc,ubc
! Subprogram not used           l = clandunit(c) 
! Subprogram not used           if (ltype(l) == isturb) then
! Subprogram not used              if (ctype(c) == icol_sunwall) then
! Subprogram not used                 scale_c2l(c) = 3.0 * canyon_hwr(l) 
! Subprogram not used              else if (ctype(c) == icol_shadewall) then
! Subprogram not used                 scale_c2l(c) = 3.0 * canyon_hwr(l) 
! Subprogram not used              else if (ctype(c) == icol_road_perv .or. ctype(c) == icol_road_imperv) then
! Subprogram not used                 scale_c2l(c) = 3.0_r8
! Subprogram not used              else if (ctype(c) == icol_roof) then
! Subprogram not used                 scale_c2l(c) = 1.0_r8
! Subprogram not used              end if
! Subprogram not used           else
! Subprogram not used              scale_c2l(c) = 1.0_r8
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used     else if (c2l_scale_type == 'urbans') then
! Subprogram not used        do c = lbc,ubc
! Subprogram not used           l = clandunit(c) 
! Subprogram not used           if (ltype(l) == isturb) then
! Subprogram not used              if (ctype(c) == icol_sunwall) then
! Subprogram not used                 scale_c2l(c) = (3.0 * canyon_hwr(l)) / (2.*canyon_hwr(l) + 1.)
! Subprogram not used              else if (ctype(c) == icol_shadewall) then
! Subprogram not used                 scale_c2l(c) = (3.0 * canyon_hwr(l)) / (2.*canyon_hwr(l) + 1.)
! Subprogram not used              else if (ctype(c) == icol_road_perv .or. ctype(c) == icol_road_imperv) then
! Subprogram not used                 scale_c2l(c) = 3.0 / (2.*canyon_hwr(l) + 1.)
! Subprogram not used              else if (ctype(c) == icol_roof) then
! Subprogram not used                 scale_c2l(c) = 1.0_r8
! Subprogram not used              end if
! Subprogram not used           else
! Subprogram not used              scale_c2l(c) = 1.0_r8
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used     else
! Subprogram not used        write(iulog,*)'p2l_2d error: scale type ',c2l_scale_type,' not supported'
! Subprogram not used        call endrun()
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     if (p2c_scale_type == 'unity') then
! Subprogram not used        do p = lbp,ubp
! Subprogram not used           scale_p2c(p) = 1.0_r8
! Subprogram not used        end do
! Subprogram not used     else
! Subprogram not used        write(iulog,*)'p2l_2d error: scale type ',p2c_scale_type,' not supported'
! Subprogram not used        call endrun()
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     larr(:,:) = spval
! Subprogram not used     do j = 1,num2d
! Subprogram not used        sumwt(:) = 0._r8
! Subprogram not used        do p = lbp,ubp
! Subprogram not used           if (wtlunit(p) /= 0._r8) then
! Subprogram not used              c = pcolumn(p)
! Subprogram not used              if (parr(p,j) /= spval .and. scale_c2l(c) /= spval) then
! Subprogram not used                 l = plandunit(p)
! Subprogram not used                 if (sumwt(l) == 0._r8) larr(l,j) = 0._r8
! Subprogram not used                 larr(l,j) = larr(l,j) + parr(p,j) * scale_p2c(p) * scale_c2l(c) * wtlunit(p)
! Subprogram not used                 sumwt(l) = sumwt(l) + wtlunit(p)
! Subprogram not used              end if
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used        found = .false.
! Subprogram not used        do l = lbl,ubl
! Subprogram not used           if (sumwt(l) > 1.0_r8 + 1.e-6_r8) then
! Subprogram not used              found = .true.
! Subprogram not used              index = l
! Subprogram not used           else if (sumwt(l) /= 0._r8) then
! Subprogram not used              larr(l,j) = larr(l,j)/sumwt(l)
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used        if (found) then
! Subprogram not used           write(iulog,*)'p2l_2d error: sumwt is greater than 1.0 at l= ',index,' j= ',j
! Subprogram not used           call endrun()
! Subprogram not used        end if
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used   end subroutine p2l_2d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: p2g_1d
!
! !INTERFACE:
  subroutine p2g_1d(lbp, ubp, lbc, ubc, lbl, ubl, lbg, ubg, parr, garr, &
       p2c_scale_type, c2l_scale_type, l2g_scale_type)
!
! !DESCRIPTION:
! Perfrom subgrid-average from pfts to gridcells.
! Averaging is only done for points that are not equal to "spval".
!
! !USES:
    use clm_varpar, only : max_pft_per_gcell
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbp, ubp            ! beginning and ending pft indices
    integer , intent(in)  :: lbc, ubc            ! beginning and ending column indices
    integer , intent(in)  :: lbl, ubl            ! beginning and ending landunit indices
    integer , intent(in)  :: lbg, ubg            ! beginning and ending gridcell indices
    real(r8), intent(in)  :: parr(lbp:ubp)       ! input pft array
    real(r8), intent(out) :: garr(lbg:ubg)       ! output gridcell array
    character(len=*), intent(in) :: p2c_scale_type ! scale factor type for averaging
    character(len=*), intent(in) :: c2l_scale_type ! scale factor type for averaging
    character(len=*), intent(in) :: l2g_scale_type ! scale factor type for averaging
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/03
!
!  !LOCAL VARIABLES:
!EOP
    integer  :: pi,p,c,l,g,index       ! indices
    logical  :: found                  ! temporary for error check
    real(r8) :: scale_p2c(lbp:ubp)     ! scale factor
    real(r8) :: scale_c2l(lbc:ubc)     ! scale factor
    real(r8) :: scale_l2g(lbl:ubl)     ! scale factor
    real(r8) :: sumwt(lbg:ubg)         ! sum of weights
    real(r8), pointer :: wtgcell(:)    ! weight of pfts relative to gridcells
    integer , pointer :: pcolumn(:)    ! column of corresponding pft
    integer , pointer :: plandunit(:)  ! landunit of corresponding pft
    integer , pointer :: pgridcell(:)  ! gridcell of corresponding pft
    integer , pointer :: npfts(:)      ! number of pfts in gridcell
    integer , pointer :: pfti(:)       ! initial pft index in gridcell
    integer , pointer :: ctype(:)      ! column type
    integer , pointer :: clandunit(:)  ! landunit of corresponding column
    integer , pointer :: ltype(:)      ! landunit type
    real(r8), pointer :: canyon_hwr(:) ! urban canyon height to width ratio
!------------------------------------------------------------------------

    canyon_hwr => lun%canyon_hwr
    ltype      => lun%itype
    clandunit  => col%landunit
    ctype      => col%itype
    wtgcell   => pft%wtgcell
    pcolumn   => pft%column
    pgridcell => pft%gridcell
    plandunit => pft%landunit
    npfts     => grc%npfts
    pfti      => grc%pfti

    call build_scale_l2g(l2g_scale_type, lbl, ubl, scale_l2g)

    if (c2l_scale_type == 'unity') then
       do c = lbc,ubc
          scale_c2l(c) = 1.0_r8
       end do
    else if (c2l_scale_type == 'urbanf') then
       do c = lbc,ubc
          l = clandunit(c) 
          if (ltype(l) == isturb) then
             if (ctype(c) == icol_sunwall) then
                scale_c2l(c) = 3.0 * canyon_hwr(l) 
             else if (ctype(c) == icol_shadewall) then
                scale_c2l(c) = 3.0 * canyon_hwr(l) 
             else if (ctype(c) == icol_road_perv .or. ctype(c) == icol_road_imperv) then
                scale_c2l(c) = 3.0_r8
             else if (ctype(c) == icol_roof) then
                scale_c2l(c) = 1.0_r8
             end if
          else
             scale_c2l(c) = 1.0_r8
          end if
       end do
    else if (c2l_scale_type == 'urbans') then
       do c = lbc,ubc
          l = clandunit(c) 
          if (ltype(l) == isturb) then
             if (ctype(c) == icol_sunwall) then
                scale_c2l(c) = (3.0 * canyon_hwr(l)) / (2.*canyon_hwr(l) + 1.)
             else if (ctype(c) == icol_shadewall) then
                scale_c2l(c) = (3.0 * canyon_hwr(l)) / (2.*canyon_hwr(l) + 1.)
             else if (ctype(c) == icol_road_perv .or. ctype(c) == icol_road_imperv) then
                scale_c2l(c) = 3.0 / (2.*canyon_hwr(l) + 1.)
             else if (ctype(c) == icol_roof) then
                scale_c2l(c) = 1.0_r8
             end if
          else
             scale_c2l(c) = 1.0_r8
          end if
       end do
    else
       write(iulog,*)'p2g_1d error: scale type ',c2l_scale_type,' not supported'
       call endrun()
    end if

    if (p2c_scale_type == 'unity') then
       do p = lbp,ubp
          scale_p2c(p) = 1.0_r8
       end do
    else
       write(iulog,*)'p2g_1d error: scale type ',c2l_scale_type,' not supported'
       call endrun()
    end if

    garr(:) = spval
    sumwt(:) = 0._r8
    do p = lbp,ubp
       if (wtgcell(p) /= 0._r8) then
          c = pcolumn(p)
          l = plandunit(p)
          if (parr(p) /= spval .and. scale_c2l(c) /= spval .and. scale_l2g(l) /= spval) then
             g = pgridcell(p)
             if (sumwt(g) == 0._r8) garr(g) = 0._r8
             garr(g) = garr(g) + parr(p) * scale_p2c(p) * scale_c2l(c) * scale_l2g(l) * wtgcell(p)
             sumwt(g) = sumwt(g) + wtgcell(p)
          end if
       end if
    end do
    found = .false.
    do g = lbg, ubg
       if (sumwt(g) > 1.0_r8 + 1.e-6_r8) then
          found = .true.
          index = g
       else if (sumwt(g) /= 0._r8) then
          garr(g) = garr(g)/sumwt(g)
       end if
    end do
    if (found) then
       write(iulog,*)'p2g_1d error: sumwt is greater than 1.0 at g= ',index
       call endrun()
    end if

  end subroutine p2g_1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: p2g_2d
!
! !INTERFACE:
  subroutine p2g_2d(lbp, ubp, lbc, ubc, lbl, ubl, lbg, ubg, num2d, &
       parr, garr, p2c_scale_type, c2l_scale_type, l2g_scale_type)
!
! !DESCRIPTION:
! Perfrom subgrid-average from pfts to gridcells.
! Averaging is only done for points that are not equal to "spval".
!
! !USES:
    use clm_varpar, only : max_pft_per_gcell
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbp, ubp              ! beginning and ending pft indices
    integer , intent(in)  :: lbc, ubc              ! beginning and ending column indices
    integer , intent(in)  :: lbl, ubl              ! beginning and ending landunit indices
    integer , intent(in)  :: lbg, ubg              ! beginning and ending gridcell indices
    integer , intent(in)  :: num2d                 ! size of second dimension
    real(r8), intent(in)  :: parr(lbp:ubp,num2d)   ! input pft array
    real(r8), intent(out) :: garr(lbg:ubg,num2d)   ! output gridcell array
    character(len=*), intent(in) :: p2c_scale_type ! scale factor type for averaging
    character(len=*), intent(in) :: c2l_scale_type ! scale factor type for averaging
    character(len=*), intent(in) :: l2g_scale_type ! scale factor type for averaging
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/03
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: j,pi,p,c,l,g,index     ! indices
    logical  :: found                  ! temporary for error check
    real(r8) :: scale_p2c(lbp:ubp)     ! scale factor
    real(r8) :: scale_c2l(lbc:ubc)     ! scale factor
    real(r8) :: scale_l2g(lbl:ubl)     ! scale factor
    real(r8) :: sumwt(lbg:ubg)         ! sum of weights
    real(r8), pointer :: wtgcell(:)    ! weight of pfts relative to gridcells
    integer , pointer :: pcolumn(:)    ! column of corresponding pft
    integer , pointer :: plandunit(:)  ! landunit of corresponding pft
    integer , pointer :: pgridcell(:)  ! gridcell of corresponding pft
    integer , pointer :: npfts(:)      ! number of pfts in gridcell
    integer , pointer :: pfti(:)       ! initial pft index in gridcell
    integer , pointer :: clandunit(:)  ! landunit of corresponding column
    integer , pointer :: ctype(:)      ! column type
    integer , pointer :: ltype(:)      ! landunit type
    real(r8), pointer :: canyon_hwr(:) ! urban canyon height to width ratio
!------------------------------------------------------------------------

    canyon_hwr   => lun%canyon_hwr
    ltype        => lun%itype
    clandunit    => col%landunit
    ctype        => col%itype
    wtgcell      => pft%wtgcell
    pcolumn      => pft%column
    pgridcell    => pft%gridcell
    plandunit    => pft%landunit
    npfts        => grc%npfts
    pfti         => grc%pfti

    call build_scale_l2g(l2g_scale_type, lbl, ubl, scale_l2g)

    if (c2l_scale_type == 'unity') then
       do c = lbc,ubc
          scale_c2l(c) = 1.0_r8
       end do
    else if (c2l_scale_type == 'urbanf') then
       do c = lbc,ubc
          l = clandunit(c) 
          if (ltype(l) == isturb) then
             if (ctype(c) == icol_sunwall) then
                scale_c2l(c) = 3.0 * canyon_hwr(l) 
             else if (ctype(c) == icol_shadewall) then
                scale_c2l(c) = 3.0 * canyon_hwr(l) 
             else if (ctype(c) == icol_road_perv .or. ctype(c) == icol_road_imperv) then
                scale_c2l(c) = 3.0_r8
             else if (ctype(c) == icol_roof) then
                scale_c2l(c) = 1.0_r8
             end if
          else
             scale_c2l(c) = 1.0_r8
          end if
       end do
    else if (c2l_scale_type == 'urbans') then
       do c = lbc,ubc
          l = clandunit(c) 
          if (ltype(l) == isturb) then
             if (ctype(c) == icol_sunwall) then
                scale_c2l(c) = (3.0 * canyon_hwr(l)) / (2.*canyon_hwr(l) + 1.)
             else if (ctype(c) == icol_shadewall) then
                scale_c2l(c) = (3.0 * canyon_hwr(l)) / (2.*canyon_hwr(l) + 1.)
             else if (ctype(c) == icol_road_perv .or. ctype(c) == icol_road_imperv) then
                scale_c2l(c) = 3.0 / (2.*canyon_hwr(l) + 1.)
             else if (ctype(c) == icol_roof) then
                scale_c2l(c) = 1.0_r8
             end if
          else
             scale_c2l(c) = 1.0_r8
          end if
       end do
    else
       write(iulog,*)'p2g_2d error: scale type ',c2l_scale_type,' not supported'
       call endrun()
    end if

    if (p2c_scale_type == 'unity') then
       do p = lbp,ubp
          scale_p2c(p) = 1.0_r8
       end do
    else
       write(iulog,*)'p2g_2d error: scale type ',c2l_scale_type,' not supported'
       call endrun()
    end if

    garr(:,:) = spval
    do j = 1,num2d
       sumwt(:) = 0._r8
       do p = lbp,ubp
          if (wtgcell(p) /= 0._r8) then
             c = pcolumn(p)
             l = plandunit(p)
             if (parr(p,j) /= spval .and. scale_c2l(c) /= spval .and. scale_l2g(l) /= spval) then
                g = pgridcell(p)
                if (sumwt(g) == 0._r8) garr(g,j) = 0._r8
                garr(g,j) = garr(g,j) + parr(p,j) * scale_p2c(p) * scale_c2l(c) * scale_l2g(l) * wtgcell(p)
                sumwt(g) = sumwt(g) + wtgcell(p)
             end if
          end if
       end do
       found = .false.
       do g = lbg, ubg
          if (sumwt(g) > 1.0_r8 + 1.e-6_r8) then
             found = .true.
             index = g
          else if (sumwt(g) /= 0._r8) then
             garr(g,j) = garr(g,j)/sumwt(g)
          end if
       end do
       if (found) then
          write(iulog,*)'p2g_2d error: sumwt gt 1.0 at g/sumwt = ',index,sumwt(index)
          call endrun()
       end if
    end do

  end subroutine p2g_2d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: c2l_1d
!
! !INTERFACE:
! Subprogram not used   subroutine c2l_1d (lbc, ubc, lbl, ubl, carr, larr, c2l_scale_type)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used ! Perfrom subgrid-average from columns to landunits
! Subprogram not used ! Averaging is only done for points that are not equal to "spval".
! Subprogram not used !
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     integer , intent(in)  :: lbc, ubc      ! beginning and ending column indices
! Subprogram not used     integer , intent(in)  :: lbl, ubl      ! beginning and ending landunit indices
! Subprogram not used     real(r8), intent(in)  :: carr(lbc:ubc) ! input column array
! Subprogram not used     real(r8), intent(out) :: larr(lbl:ubl) ! output landunit array
! Subprogram not used     character(len=*), intent(in) :: c2l_scale_type ! scale factor type for averaging
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! Created by Mariana Vertenstein 12/03
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! !LOCAL VARIABLES:
! Subprogram not used !EOP
! Subprogram not used     integer  :: ci,c,l,index           ! indices
! Subprogram not used     integer  :: max_col_per_lu         ! max columns per landunit; on the fly
! Subprogram not used     logical  :: found                  ! temporary for error check
! Subprogram not used     real(r8) :: scale_c2l(lbc:ubc)     ! scale factor for column->landunit mapping
! Subprogram not used     real(r8) :: sumwt(lbl:ubl)         ! sum of weights
! Subprogram not used     real(r8), pointer :: wtlunit(:)    ! weight of landunits relative to gridcells
! Subprogram not used     integer , pointer :: clandunit(:)  ! gridcell of corresponding column
! Subprogram not used     integer , pointer :: ncolumns(:)   ! number of columns in landunit
! Subprogram not used     integer , pointer :: coli(:)       ! initial column index in landunit
! Subprogram not used     integer , pointer :: ctype(:)      ! column type
! Subprogram not used     integer , pointer :: ltype(:)      ! landunit type
! Subprogram not used     real(r8), pointer :: canyon_hwr(:) ! urban canyon height to width ratio
! Subprogram not used !------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     ctype      => col%itype
! Subprogram not used     ltype      => lun%itype
! Subprogram not used     canyon_hwr => lun%canyon_hwr
! Subprogram not used     wtlunit    => col%wtlunit
! Subprogram not used     clandunit  => col%landunit
! Subprogram not used     ncolumns   => lun%ncolumns
! Subprogram not used     coli       => lun%coli
! Subprogram not used 
! Subprogram not used     if (c2l_scale_type == 'unity') then
! Subprogram not used        do c = lbc,ubc
! Subprogram not used           scale_c2l(c) = 1.0_r8
! Subprogram not used        end do
! Subprogram not used     else if (c2l_scale_type == 'urbanf') then
! Subprogram not used        do c = lbc,ubc
! Subprogram not used           l = clandunit(c) 
! Subprogram not used           if (ltype(l) == isturb) then
! Subprogram not used              if (ctype(c) == icol_sunwall) then
! Subprogram not used                 scale_c2l(c) = 3.0 * canyon_hwr(l) 
! Subprogram not used              else if (ctype(c) == icol_shadewall) then
! Subprogram not used                 scale_c2l(c) = 3.0 * canyon_hwr(l) 
! Subprogram not used              else if (ctype(c) == icol_road_perv .or. ctype(c) == icol_road_imperv) then
! Subprogram not used                 scale_c2l(c) = 3.0_r8
! Subprogram not used              else if (ctype(c) == icol_roof) then
! Subprogram not used                 scale_c2l(c) = 1.0_r8
! Subprogram not used              end if
! Subprogram not used           else
! Subprogram not used              scale_c2l(c) = 1.0_r8
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used     else if (c2l_scale_type == 'urbans') then
! Subprogram not used        do c = lbc,ubc
! Subprogram not used           l = clandunit(c) 
! Subprogram not used           if (ltype(l) == isturb) then
! Subprogram not used              if (ctype(c) == icol_sunwall) then
! Subprogram not used                 scale_c2l(c) = (3.0 * canyon_hwr(l)) / (2.*canyon_hwr(l) + 1.)
! Subprogram not used              else if (ctype(c) == icol_shadewall) then
! Subprogram not used                 scale_c2l(c) = (3.0 * canyon_hwr(l)) / (2.*canyon_hwr(l) + 1.)
! Subprogram not used              else if (ctype(c) == icol_road_perv .or. ctype(c) == icol_road_imperv) then
! Subprogram not used                 scale_c2l(c) = 3.0 / (2.*canyon_hwr(l) + 1.)
! Subprogram not used              else if (ctype(c) == icol_roof) then
! Subprogram not used                 scale_c2l(c) = 1.0_r8
! Subprogram not used              end if
! Subprogram not used           else
! Subprogram not used              scale_c2l(c) = 1.0_r8
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used     else
! Subprogram not used        write(iulog,*)'c2l_1d error: scale type ',c2l_scale_type,' not supported'
! Subprogram not used        call endrun()
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     larr(:) = spval
! Subprogram not used     sumwt(:) = 0._r8
! Subprogram not used     do c = lbc,ubc
! Subprogram not used        if (wtlunit(c) /= 0._r8) then
! Subprogram not used           if (carr(c) /= spval .and. scale_c2l(c) /= spval) then
! Subprogram not used              l = clandunit(c)
! Subprogram not used              if (sumwt(l) == 0._r8) larr(l) = 0._r8
! Subprogram not used              larr(l) = larr(l) + carr(c) * scale_c2l(c) * wtlunit(c)
! Subprogram not used              sumwt(l) = sumwt(l) + wtlunit(c)
! Subprogram not used           end if
! Subprogram not used        end if
! Subprogram not used     end do
! Subprogram not used     found = .false.
! Subprogram not used     do l = lbl,ubl
! Subprogram not used        if (sumwt(l) > 1.0_r8 + 1.e-6_r8) then
! Subprogram not used           found = .true.
! Subprogram not used           index = l
! Subprogram not used        else if (sumwt(l) /= 0._r8) then
! Subprogram not used           larr(l) = larr(l)/sumwt(l)
! Subprogram not used        end if
! Subprogram not used     end do
! Subprogram not used     if (found) then
! Subprogram not used        write(iulog,*)'c2l_1d error: sumwt is greater than 1.0 at l= ',index
! Subprogram not used        call endrun()
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used   end subroutine c2l_1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: c2l_2d
!
! !INTERFACE:
! Subprogram not used   subroutine c2l_2d (lbc, ubc, lbl, ubl, num2d, carr, larr, c2l_scale_type)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used ! Perfrom subgrid-average from columns to landunits
! Subprogram not used ! Averaging is only done for points that are not equal to "spval".
! Subprogram not used !
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     integer , intent(in)  :: lbc, ubc            ! beginning and ending column indices
! Subprogram not used     integer , intent(in)  :: lbl, ubl            ! beginning and ending landunit indices
! Subprogram not used     integer , intent(in)  :: num2d               ! size of second dimension
! Subprogram not used     real(r8), intent(in)  :: carr(lbc:ubc,num2d) ! input column array
! Subprogram not used     real(r8), intent(out) :: larr(lbl:ubl,num2d) ! output landunit array
! Subprogram not used     character(len=*), intent(in) :: c2l_scale_type ! scale factor type for averaging
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! Created by Mariana Vertenstein 12/03
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! !LOCAL VARIABLES:
! Subprogram not used !EOP
! Subprogram not used     integer  :: j,l,ci,c,index         ! indices
! Subprogram not used     integer  :: max_col_per_lu         ! max columns per landunit; on the fly
! Subprogram not used     logical  :: found                  ! temporary for error check
! Subprogram not used     real(r8) :: scale_c2l(lbc:ubc)        ! scale factor for column->landunit mapping
! Subprogram not used     real(r8) :: sumwt(lbl:ubl)         ! sum of weights
! Subprogram not used     real(r8), pointer :: wtlunit(:)    ! weight of column relative to landunit
! Subprogram not used     integer , pointer :: clandunit(:)  ! landunit of corresponding column
! Subprogram not used     integer , pointer :: ncolumns(:)   ! number of columns in landunit
! Subprogram not used     integer , pointer :: coli(:)       ! initial column index in landunit
! Subprogram not used     integer , pointer :: ctype(:)      ! column type
! Subprogram not used     integer , pointer :: ltype(:)      ! landunit type
! Subprogram not used     real(r8), pointer :: canyon_hwr(:) ! urban canyon height to width ratio
! Subprogram not used !------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     ctype      => col%itype
! Subprogram not used     ltype      => lun%itype
! Subprogram not used     canyon_hwr => lun%canyon_hwr
! Subprogram not used     wtlunit    => col%wtlunit
! Subprogram not used     clandunit  => col%landunit
! Subprogram not used     ncolumns   => lun%ncolumns
! Subprogram not used     coli       => lun%coli
! Subprogram not used 
! Subprogram not used     if (c2l_scale_type == 'unity') then
! Subprogram not used        do c = lbc,ubc
! Subprogram not used           scale_c2l(c) = 1.0_r8
! Subprogram not used        end do
! Subprogram not used     else if (c2l_scale_type == 'urbanf') then
! Subprogram not used        do c = lbc,ubc
! Subprogram not used           l = clandunit(c) 
! Subprogram not used           if (ltype(l) == isturb) then
! Subprogram not used              if (ctype(c) == icol_sunwall) then
! Subprogram not used                 scale_c2l(c) = 3.0 * canyon_hwr(l) 
! Subprogram not used              else if (ctype(c) == icol_shadewall) then
! Subprogram not used                 scale_c2l(c) = 3.0 * canyon_hwr(l) 
! Subprogram not used              else if (ctype(c) == icol_road_perv .or. ctype(c) == icol_road_imperv) then
! Subprogram not used                 scale_c2l(c) = 3.0_r8
! Subprogram not used              else if (ctype(c) == icol_roof) then
! Subprogram not used                 scale_c2l(c) = 1.0_r8
! Subprogram not used              end if
! Subprogram not used           else
! Subprogram not used              scale_c2l(c) = 1.0_r8
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used     else if (c2l_scale_type == 'urbans') then
! Subprogram not used        do c = lbc,ubc
! Subprogram not used           l = clandunit(c) 
! Subprogram not used           if (ltype(l) == isturb) then
! Subprogram not used              if (ctype(c) == icol_sunwall) then
! Subprogram not used                 scale_c2l(c) = (3.0 * canyon_hwr(l)) / (2.*canyon_hwr(l) + 1.)
! Subprogram not used              else if (ctype(c) == icol_shadewall) then
! Subprogram not used                 scale_c2l(c) = (3.0 * canyon_hwr(l)) / (2.*canyon_hwr(l) + 1.)
! Subprogram not used              else if (ctype(c) == icol_road_perv .or. ctype(c) == icol_road_imperv) then
! Subprogram not used                 scale_c2l(c) = 3.0 / (2.*canyon_hwr(l) + 1.)
! Subprogram not used              else if (ctype(c) == icol_roof) then
! Subprogram not used                 scale_c2l(c) = 1.0_r8
! Subprogram not used              end if
! Subprogram not used           else
! Subprogram not used              scale_c2l(c) = 1.0_r8
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used     else
! Subprogram not used        write(iulog,*)'c2l_2d error: scale type ',c2l_scale_type,' not supported'
! Subprogram not used        call endrun()
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     larr(:,:) = spval
! Subprogram not used     do j = 1,num2d
! Subprogram not used        sumwt(:) = 0._r8
! Subprogram not used        do c = lbc,ubc
! Subprogram not used           if (wtlunit(c) /= 0._r8) then
! Subprogram not used              if (carr(c,j) /= spval .and. scale_c2l(c) /= spval) then
! Subprogram not used                 l = clandunit(c)
! Subprogram not used                 if (sumwt(l) == 0._r8) larr(l,j) = 0._r8
! Subprogram not used                 larr(l,j) = larr(l,j) + carr(c,j) * scale_c2l(c) * wtlunit(c)
! Subprogram not used                 sumwt(l) = sumwt(l) + wtlunit(c)
! Subprogram not used              end if
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used        found = .false.
! Subprogram not used        do l = lbl,ubl
! Subprogram not used           if (sumwt(l) > 1.0_r8 + 1.e-6_r8) then
! Subprogram not used              found = .true.
! Subprogram not used              index = l
! Subprogram not used           else if (sumwt(l) /= 0._r8) then
! Subprogram not used              larr(l,j) = larr(l,j)/sumwt(l)
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used        if (found) then
! Subprogram not used           write(iulog,*)'c2l_2d error: sumwt is greater than 1.0 at l= ',index,' lev= ',j
! Subprogram not used           call endrun()
! Subprogram not used        end if
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used   end subroutine c2l_2d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: c2g_1d
!
! !INTERFACE:
  subroutine c2g_1d(lbc, ubc, lbl, ubl, lbg, ubg, carr, garr, &
       c2l_scale_type, l2g_scale_type)
!
! !DESCRIPTION:
! Perfrom subgrid-average from columns to gridcells.
! Averaging is only done for points that are not equal to "spval".
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbc, ubc              ! beginning and ending column indices
    integer , intent(in)  :: lbl, ubl              ! beginning and ending landunit indices
    integer , intent(in)  :: lbg, ubg              ! beginning and ending landunit indices
    real(r8), intent(in)  :: carr(lbc:ubc)         ! input column array
    real(r8), intent(out) :: garr(lbg:ubg)         ! output gridcell array
    character(len=*), intent(in) :: c2l_scale_type ! scale factor type for averaging
    character(len=*), intent(in) :: l2g_scale_type ! scale factor type for averaging
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/03
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: ci,c,l,g,index         ! indices
    integer  :: max_col_per_gcell      ! max columns per gridcell; on the fly
    logical  :: found                  ! temporary for error check
    real(r8) :: scale_c2l(lbc:ubc)     ! scale factor
    real(r8) :: scale_l2g(lbl:ubl)     ! scale factor
    real(r8) :: sumwt(lbg:ubg)         ! sum of weights
    real(r8), pointer :: wtgcell(:)    ! weight of columns relative to gridcells
    integer , pointer :: clandunit(:)  ! landunit of corresponding column
    integer , pointer :: cgridcell(:)  ! gridcell of corresponding column
    integer , pointer :: ncolumns(:)   ! number of columns in gridcell
    integer , pointer :: coli(:)       ! initial column index in gridcell
    integer , pointer :: ctype(:)      ! column type
    integer , pointer :: ltype(:)      ! landunit type
    real(r8), pointer :: canyon_hwr(:) ! urban canyon height to width ratio
!------------------------------------------------------------------------

    ctype      => col%itype
    ltype      => lun%itype
    canyon_hwr => lun%canyon_hwr
    wtgcell    => col%wtgcell
    clandunit  => col%landunit
    cgridcell  => col%gridcell
    ncolumns   => grc%ncolumns
    coli       => grc%coli

    call build_scale_l2g(l2g_scale_type, lbl, ubl, scale_l2g)

    if (c2l_scale_type == 'unity') then
       do c = lbc,ubc
          scale_c2l(c) = 1.0_r8
       end do
    else if (c2l_scale_type == 'urbanf') then
       do c = lbc,ubc
          l = clandunit(c) 
          if (ltype(l) == isturb) then
             if (ctype(c) == icol_sunwall) then
                scale_c2l(c) = 3.0 * canyon_hwr(l) 
             else if (ctype(c) == icol_shadewall) then
                scale_c2l(c) = 3.0 * canyon_hwr(l) 
             else if (ctype(c) == icol_road_perv .or. ctype(c) == icol_road_imperv) then
                scale_c2l(c) = 3.0_r8
             else if (ctype(c) == icol_roof) then
                scale_c2l(c) = 1.0_r8
             end if
          else
             scale_c2l(c) = 1.0_r8
          end if
       end do
    else if (c2l_scale_type == 'urbans') then
       do c = lbc,ubc
          l = clandunit(c) 
          if (ltype(l) == isturb) then
             if (ctype(c) == icol_sunwall) then
                scale_c2l(c) = (3.0 * canyon_hwr(l)) / (2.*canyon_hwr(l) + 1.)
             else if (ctype(c) == icol_shadewall) then
                scale_c2l(c) = (3.0 * canyon_hwr(l)) / (2.*canyon_hwr(l) + 1.)
             else if (ctype(c) == icol_road_perv .or. ctype(c) == icol_road_imperv) then
                scale_c2l(c) = 3.0 / (2.*canyon_hwr(l) + 1.)
             else if (ctype(c) == icol_roof) then
                scale_c2l(c) = 1.0_r8
             end if
          else
             scale_c2l(c) = 1.0_r8
          end if
       end do
    else
       write(iulog,*)'c2l_1d error: scale type ',c2l_scale_type,' not supported'
       call endrun()
    end if

    garr(:) = spval
    sumwt(:) = 0._r8
    do c = lbc,ubc
       if ( wtgcell(c) /= 0._r8) then
          l = clandunit(c)
          if (carr(c) /= spval .and. scale_c2l(c) /= spval .and. scale_l2g(l) /= spval) then
             g = cgridcell(c)
             if (sumwt(g) == 0._r8) garr(g) = 0._r8
             garr(g) = garr(g) + carr(c) * scale_c2l(c) * scale_l2g(l) * wtgcell(c)
             sumwt(g) = sumwt(g) + wtgcell(c)
          end if
       end if
    end do
    found = .false.
    do g = lbg, ubg
       if (sumwt(g) > 1.0_r8 + 1.e-6_r8) then
          found = .true.
          index = g
       else if (sumwt(g) /= 0._r8) then
          garr(g) = garr(g)/sumwt(g)
       end if
    end do
    if (found) then
       write(iulog,*)'c2g_1d error: sumwt is greater than 1.0 at g= ',index
       call endrun()
    end if

  end subroutine c2g_1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: c2g_2d
!
! !INTERFACE:
  subroutine c2g_2d(lbc, ubc, lbl, ubl, lbg, ubg, num2d, carr, garr, &
       c2l_scale_type, l2g_scale_type)
!
! !DESCRIPTION:
! Perfrom subgrid-average from columns to gridcells.
! Averaging is only done for points that are not equal to "spval".
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbc, ubc              ! beginning and ending column indices
    integer , intent(in)  :: lbl, ubl              ! beginning and ending landunit indices
    integer , intent(in)  :: lbg, ubg              ! beginning and ending gridcell indices
    integer , intent(in)  :: num2d                 ! size of second dimension
    real(r8), intent(in)  :: carr(lbc:ubc,num2d)   ! input column array
    real(r8), intent(out) :: garr(lbg:ubg,num2d)   ! output gridcell array
    character(len=*), intent(in) :: c2l_scale_type ! scale factor type for averaging
    character(len=*), intent(in) :: l2g_scale_type ! scale factor type for averaging
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/03
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: j,ci,c,g,l,index       ! indices
    integer  :: max_col_per_gcell      ! max columns per gridcell; on the fly
    logical  :: found                  ! temporary for error check
    real(r8) :: scale_c2l(lbc:ubc)     ! scale factor
    real(r8) :: scale_l2g(lbl:ubl)     ! scale factor
    real(r8) :: sumwt(lbg:ubg)         ! sum of weights
    real(r8), pointer :: wtgcell(:)    ! weight of columns relative to gridcells
    integer , pointer :: clandunit(:)  ! landunit of corresponding column
    integer , pointer :: cgridcell(:)  ! gridcell of corresponding column
    integer , pointer :: ncolumns(:)   ! number of columns in gridcell
    integer , pointer :: coli(:)       ! initial column index in gridcell
    integer , pointer :: ctype(:)      ! column type
    integer , pointer :: ltype(:)      ! landunit type
    real(r8), pointer :: canyon_hwr(:) ! urban canyon height to width ratio
!------------------------------------------------------------------------

    ctype      => col%itype
    ltype      => lun%itype
    canyon_hwr => lun%canyon_hwr
    wtgcell    => col%wtgcell
    clandunit  => col%landunit
    cgridcell  => col%gridcell
    ncolumns   => grc%ncolumns
    coli       => grc%coli

    call build_scale_l2g(l2g_scale_type, lbl, ubl, scale_l2g)

    if (c2l_scale_type == 'unity') then
       do c = lbc,ubc
          scale_c2l(c) = 1.0_r8
       end do
    else if (c2l_scale_type == 'urbanf') then
       do c = lbc,ubc
          l = clandunit(c) 
          if (ltype(l) == isturb) then
             if (ctype(c) == icol_sunwall) then
                scale_c2l(c) = 3.0 * canyon_hwr(l) 
             else if (ctype(c) == icol_shadewall) then
                scale_c2l(c) = 3.0 * canyon_hwr(l) 
             else if (ctype(c) == icol_road_perv .or. ctype(c) == icol_road_imperv) then
                scale_c2l(c) = 3.0_r8
             else if (ctype(c) == icol_roof) then
                scale_c2l(c) = 1.0_r8
             end if
          else
             scale_c2l(c) = 1.0_r8
          end if
       end do
    else if (c2l_scale_type == 'urbans') then
       do c = lbc,ubc
          l = clandunit(c) 
          if (ltype(l) == isturb) then
             if (ctype(c) == icol_sunwall) then
                scale_c2l(c) = (3.0 * canyon_hwr(l)) / (2.*canyon_hwr(l) + 1.)
             else if (ctype(c) == icol_shadewall) then
                scale_c2l(c) = (3.0 * canyon_hwr(l)) / (2.*canyon_hwr(l) + 1.)
             else if (ctype(c) == icol_road_perv .or. ctype(c) == icol_road_imperv) then
                scale_c2l(c) = 3.0 / (2.*canyon_hwr(l) + 1.)
             else if (ctype(c) == icol_roof) then
                scale_c2l(c) = 1.0_r8
             end if
          else
             scale_c2l(c) = 1.0_r8
          end if
       end do
    else
       write(iulog,*)'c2g_2d error: scale type ',c2l_scale_type,' not supported'
       call endrun()
    end if

    garr(:,:) = spval
    do j = 1,num2d
       sumwt(:) = 0._r8
       do c = lbc,ubc
          if (wtgcell(c) /= 0._r8) then
             l = clandunit(c)
             if (carr(c,j) /= spval .and. scale_c2l(c) /= spval .and. scale_l2g(l) /= spval) then
                g = cgridcell(c)
                if (sumwt(g) == 0._r8) garr(g,j) = 0._r8
                garr(g,j) = garr(g,j) + carr(c,j) * scale_c2l(c) * scale_l2g(l) * wtgcell(c)
                sumwt(g) = sumwt(g) + wtgcell(c)
             end if
          end if
       end do
       found = .false.
       do g = lbg, ubg
          if (sumwt(g) > 1.0_r8 + 1.e-6_r8) then
             found = .true.
             index = g
          else if (sumwt(g) /= 0._r8) then
             garr(g,j) = garr(g,j)/sumwt(g)
          end if
       end do
       if (found) then
          write(iulog,*)'c2g_2d error: sumwt is greater than 1.0 at g= ',index
          call endrun()
       end if
    end do

  end subroutine c2g_2d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: l2g_1d
!
! !INTERFACE:
  subroutine l2g_1d(lbl, ubl, lbg, ubg, larr, garr, l2g_scale_type)
!
! !DESCRIPTION:
! Perfrom subgrid-average from landunits to gridcells.
! Averaging is only done for points that are not equal to "spval".
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbl, ubl       ! beginning and ending sub landunit indices
    integer , intent(in)  :: lbg, ubg       ! beginning and ending gridcell indices
    real(r8), intent(in)  :: larr(lbl:ubl)  ! input landunit array
    real(r8), intent(out) :: garr(lbg:ubg)  ! output gridcell array
    character(len=*), intent(in) :: l2g_scale_type ! scale factor type for averaging
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/03
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: li,l,g,index           ! indices
    integer  :: max_lu_per_gcell       ! max landunits per gridcell; on the fly
    logical  :: found                  ! temporary for error check
    real(r8) :: scale_l2g(lbl:ubl)     ! scale factor
    real(r8) :: sumwt(lbg:ubg)         ! sum of weights
    real(r8), pointer :: wtgcell(:)    ! weight of landunits relative to gridcells
    integer , pointer :: lgridcell(:)  ! gridcell of corresponding landunit
    integer , pointer :: nlandunits(:) ! number of landunits in gridcell
    integer , pointer :: luni(:)       ! initial landunit index in gridcell
!------------------------------------------------------------------------

    wtgcell    => lun%wtgcell
    lgridcell  => lun%gridcell
    nlandunits => grc%nlandunits
    luni       => grc%luni

    call build_scale_l2g(l2g_scale_type, lbl, ubl, scale_l2g)

    garr(:) = spval
    sumwt(:) = 0._r8
    do l = lbl,ubl
       if (wtgcell(l) /= 0._r8) then
          if (larr(l) /= spval .and. scale_l2g(l) /= spval) then
             g = lgridcell(l)
             if (sumwt(g) == 0._r8) garr(g) = 0._r8
             garr(g) = garr(g) + larr(l) * scale_l2g(l) * wtgcell(l)
             sumwt(g) = sumwt(g) + wtgcell(l)
          end if
       end if
    end do
    found = .false.
    do g = lbg, ubg
       if (sumwt(g) > 1.0_r8 + 1.e-6_r8) then
          found = .true.
          index = g
       else if (sumwt(g) /= 0._r8) then
          garr(g) = garr(g)/sumwt(g)
       end if
    end do
    if (found) then
       write(iulog,*)'l2g_1d error: sumwt is greater than 1.0 at g= ',index
       call endrun()
    end if

  end subroutine l2g_1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: l2g_2d
!
! !INTERFACE:
! Subprogram not used   subroutine l2g_2d(lbl, ubl, lbg, ubg, num2d, larr, garr, l2g_scale_type)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used ! Perfrom subgrid-average from landunits to gridcells.
! Subprogram not used ! Averaging is only done for points that are not equal to "spval".
! Subprogram not used !
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     integer , intent(in)  :: lbl, ubl             ! beginning and ending column indices
! Subprogram not used     integer , intent(in)  :: lbg, ubg             ! beginning and ending gridcell indices
! Subprogram not used     integer , intent(in)  :: num2d                ! size of second dimension
! Subprogram not used     real(r8), intent(in)  :: larr(lbl:ubl,num2d)  ! input landunit array
! Subprogram not used     real(r8), intent(out) :: garr(lbg:ubg,num2d)  ! output gridcell array
! Subprogram not used     character(len=*), intent(in) :: l2g_scale_type ! scale factor type for averaging
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! Created by Mariana Vertenstein 12/03
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! !LOCAL VARIABLES:
! Subprogram not used !EOP
! Subprogram not used     integer  :: j,g,li,l,index         ! indices
! Subprogram not used     integer  :: max_lu_per_gcell       ! max landunits per gridcell; on the fly
! Subprogram not used     logical  :: found                  ! temporary for error check
! Subprogram not used     real(r8) :: scale_l2g(lbl:ubl)     ! scale factor
! Subprogram not used     real(r8) :: sumwt(lbg:ubg)         ! sum of weights
! Subprogram not used     real(r8), pointer :: wtgcell(:)    ! weight of landunits relative to gridcells
! Subprogram not used     integer , pointer :: lgridcell(:)  ! gridcell of corresponding landunit
! Subprogram not used     integer , pointer :: nlandunits(:) ! number of landunits in gridcell
! Subprogram not used     integer , pointer :: luni(:)       ! initial landunit index in gridcell
! Subprogram not used !------------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     wtgcell   => lun%wtgcell
! Subprogram not used     lgridcell => lun%gridcell
! Subprogram not used     nlandunits => grc%nlandunits
! Subprogram not used     luni       => grc%luni
! Subprogram not used 
! Subprogram not used     call build_scale_l2g(l2g_scale_type, lbl, ubl, scale_l2g)
! Subprogram not used 
! Subprogram not used     garr(:,:) = spval
! Subprogram not used     do j = 1,num2d
! Subprogram not used        sumwt(:) = 0._r8
! Subprogram not used        do l = lbl,ubl
! Subprogram not used           if (wtgcell(l) /= 0._r8) then
! Subprogram not used              if (larr(l,j) /= spval .and. scale_l2g(l) /= spval) then
! Subprogram not used                 g = lgridcell(l)
! Subprogram not used                 if (sumwt(g) == 0._r8) garr(g,j) = 0._r8
! Subprogram not used                 garr(g,j) = garr(g,j) + larr(l,j) * scale_l2g(l) * wtgcell(l)
! Subprogram not used                 sumwt(g) = sumwt(g) + wtgcell(l)
! Subprogram not used              end if
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used        found = .false.
! Subprogram not used        do g = lbg,ubg
! Subprogram not used           if (sumwt(g) > 1.0_r8 + 1.e-6_r8) then
! Subprogram not used              found = .true.
! Subprogram not used              index= g
! Subprogram not used           else if (sumwt(g) /= 0._r8) then
! Subprogram not used              garr(g,j) = garr(g,j)/sumwt(g)
! Subprogram not used           end if
! Subprogram not used        end do
! Subprogram not used        if (found) then
! Subprogram not used           write(iulog,*)'l2g_2d error: sumwt is greater than 1.0 at g= ',index,' lev= ',j
! Subprogram not used           call endrun()
! Subprogram not used        end if
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used   end subroutine l2g_2d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: build_scale_l2g
!
! !INTERFACE:
  subroutine build_scale_l2g(l2g_scale_type, lbl, ubl, scale_l2g)
!
! !DESCRIPTION:
! Fill the scale_l2g(lbl:ubl) array with appropriate values for the given l2g_scale_type.
! This array can later be used to scale each landunit in forming grid cell averages.
!
! !USES:
     use clm_varcon, only : max_lunit
!
! !ARGUMENTS:
     implicit none
     character(len=*), intent(in)  :: l2g_scale_type     ! scale factor type for averaging
     integer         , intent(in)  :: lbl, ubl           ! beginning and ending column indices
     real(r8)        , intent(out) :: scale_l2g(lbl:ubl) ! scale factor 
!
! !REVISION HISTORY:
! Created by Bill Sacks 10/11
!
!
! !LOCAL VARIABLES:
!EOP
     real(r8) :: scale_lookup(max_lunit) ! scale factor for each landunit type
     integer  :: l                       ! index
     integer , pointer :: ltype(:)       ! landunit type
!-----------------------------------------------------------------------

     ltype      => lun%itype
     
     call create_scale_l2g_lookup(l2g_scale_type, scale_lookup)

     do l = lbl,ubl
        scale_l2g(l) = scale_lookup(ltype(l))
     end do

  end subroutine build_scale_l2g

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: create_scale_l2g_lookup
!
! !INTERFACE:
  subroutine create_scale_l2g_lookup(l2g_scale_type, scale_lookup)
! 
! DESCRIPTION:
! Create a lookup array, scale_lookup(1..max_lunit), which gives the scale factor for
! each landunit type depending on l2g_scale_type
!
! !USES:
     use clm_varcon, only : istsoil, istice, istdlak, istslak, istwet, isturb, istice_mec,&
                            istcrop, max_lunit, spval
!
! !ARGUMENTS:
     implicit none
     character(len=*), intent(in)  :: l2g_scale_type           ! scale factor type for averaging
     real(r8)        , intent(out) :: scale_lookup(max_lunit)  ! scale factor for each landunit type
!
! !REVISION HISTORY:
! Created by Bill Sacks 10/11
!
!EOP
!-----------------------------------------------------------------------

     ! ------------ WJS (10-14-11): IMPORTANT GENERAL NOTES ------------
     !
     ! Since scale_l2g is not currently included in the sumwt accumulations, you need to
     ! be careful about the scale values you use. Values of 1 and spval are safe
     ! (including having multiple landunits with value 1), but only use other values if
     ! you know what you are doing! For example, using a value of 0 is NOT the correct way
     ! to exclude a landunit from the average, because the normalization will be done
     ! incorrectly in this case: instead, use spval to exclude a landunit from the
     ! average. Similarly, using a value of 2 is NOT the correct way to give a landunit
     ! double relative weight in general, because the normalization won't be done
     ! correctly in this case, either.
     !
     ! In the longer-term, I believe that the correct solution to this problem is to
     ! include scale_l2g (and the other scale factors) in the sumwt accumulations
     ! (e.g., sumwt = sumwt + wtgcell * scale_p2c * scale_c2l * scale_l2g), but that
     ! requires some more thought to (1) make sure that is correct, and (2) make sure it
     ! doesn't break the urban scaling.
     !
     ! -----------------------------------------------------------------


     ! Initialize scale_lookup to spval for all landunits. Thus, any landunit that keeps
     ! the default value will be excluded from grid cell averages.
     scale_lookup(:) = spval

     if (l2g_scale_type == 'unity') then
        scale_lookup(:) = 1.0_r8
     else if (l2g_scale_type == 'veg') then
        scale_lookup(istsoil) = 1.0_r8
        scale_lookup(istcrop) = 1.0_r8
     else if (l2g_scale_type == 'ice') then
        scale_lookup(istice) = 1.0_r8
        scale_lookup(istice_mec) = 1.0_r8
     else if (l2g_scale_type == 'nonurb') then
        scale_lookup(:) = 1.0_r8
        scale_lookup(isturb) = spval
     else
        write(iulog,*)'scale_l2g_lookup_array error: scale type ',l2g_scale_type,' not supported'
        call endrun()
     end if

  end subroutine create_scale_l2g_lookup

end module subgridAveMod
