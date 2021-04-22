module ice_probability

    use ice_kinds_mod 
    use ice_communicate, only: my_task, master_task, get_num_procs
    use ice_constants, only: pi, puny, c1
    use ice_fileunits, only: nu_timing
    use ice_domain_size, only: nx_global, ny_global, max_blocks
    use ice_domain, only: nblocks, blocks_ice, distrb_info, distribution_wght_file
    use ice_blocks, only: nx_block, ny_block,  nblocks_tot, block, get_block
    use ice_read_write, only: ice_open_nc, ice_read_global_nc, ice_close_nc   
    use ice_broadcast, only: broadcast_array
    use ice_global_reductions, only: sum_vector_dbl
    use ice_timers, only: 
    use ice_work, only: work_gr
    use ice_probability_tools

    implicit none 

   integer (int_kind), public, parameter ::  &    ! types of blocks:
         lndType     = 0,                    &    !     Land
         icefreeType = 1,                    &    !     ice free (ocean only)
         iceType     = 2                          !     sea ice 

!    include "netcdf.inc"
   integer (int_kind), public              :: dynCnt          ! number of calls to step_dynamic
   real (dbl_kind), allocatable :: lnumIceCells(:)  ! number of active ice cells

   public :: ReadProbabilityFile

   public :: CalcWorkPerBlock

!   public :: WriteProbabilityStats

   public :: init_numIceCells,   &
             accum_numIceCells,  &
             accum_numIceCells2, &
             print_numIceCells

   public :: set_numIceCells,write_numIceCells 
   

contains 
   
   subroutine init_numIceCells()

       dynCnt = 0

!       print *,'init_numIceCells: nblocks_tot is: ',nblocks_tot
       allocate(lnumIceCells(nblocks_tot))
       lnumIceCells = 0.0

   end subroutine init_numIceCells

! Subprogram not used    subroutine accum_numIceCells(iblk,icells)
! Subprogram not used 
! Subprogram not used      integer (int_kind) :: iblk,icells
! Subprogram not used 
! Subprogram not used      lnumIceCells(iblk) = lnumIceCells(iblk) + real(icells,kind=dbl_kind)
! Subprogram not used 
! Subprogram not used    end subroutine accum_numIceCells

   subroutine accum_numIceCells2(aice)

     real (dbl_kind) :: aice(nx_block,ny_block,max_blocks)
 
     integer (int_kind) :: igblk,iblk

     real (dbl_kind) :: tmp
     type (block) :: this_block
     integer (int_kind) :: i,j,ihi,ilo,jhi,jlo
  
  
     do iblk = 1,nblocks     
        igblk = blocks_ice(iblk)
        this_block = get_block(igblk,iblk)
        ilo = this_block%ilo
        ihi = this_block%ihi
        jlo = this_block%jlo
        jhi = this_block%jhi
        tmp = 0.0
        do j = jlo,jhi
        do i = ilo,ihi
           if(aice(i,j,iblk) > puny) then 
              tmp = tmp + c1
           endif
        enddo
        enddo 
        lnumIceCells(igblk) = lnumIceCells(igblk) + tmp
     enddo

   end subroutine accum_numIceCells2



! Subprogram not used    subroutine set_numIceCells(iblk,ncells)
! Subprogram not used       integer (int_kind) :: iblk
! Subprogram not used       real (dbl_kind)  :: ncells
! Subprogram not used     
! Subprogram not used      lnumIceCells(iblk) = ncells 
! Subprogram not used   
! Subprogram not used    end subroutine set_numIceCells

! Subprogram not used    subroutine write_numIceCells
! Subprogram not used     real (dbl_kind), allocatable :: gnumIceCells(:)
! Subprogram not used     integer :: n
! Subprogram not used 
! Subprogram not used     allocate(gnumIceCells(nblocks_tot))
! Subprogram not used     
! Subprogram not used     call sum_vector_dbl(lnumIceCells,gnumIceCells,distrb_info)
! Subprogram not used 
! Subprogram not used     gnumIceCells=gnumIceCells/real(dynCnt,kind=dbl_kind)
! Subprogram not used     if(my_task == master_task) then
! Subprogram not used        open(nu_timing,file='numCells2.bin',recl=8*nblocks_tot, &
! Subprogram not used           form = 'unformatted', access = 'direct', status = 'unknown')
! Subprogram not used         write(nu_timing,rec=1) gnumIceCells
! Subprogram not used         close(nu_timing)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used  
! Subprogram not used 
! Subprogram not used    end subroutine write_numIceCells

! Subprogram not used    subroutine print_numIceCells
! Subprogram not used 
! Subprogram not used       real (dbl_kind), allocatable :: gnumIceCells(:)
! Subprogram not used       real (dbl_kind), allocatable :: gnumIceCells2(:)
! Subprogram not used       integer :: ii,n
! Subprogram not used 
! Subprogram not used       allocate(gnumIceCells(nblocks_tot))
! Subprogram not used       allocate(gnumIceCells2(nblocks_tot))
! Subprogram not used 
! Subprogram not used       call sum_vector_dbl(lnumIceCells, gnumIceCells, distrb_info)
! Subprogram not used 
! Subprogram not used       gnumIceCells = gnumIceCells/real(dynCnt,kind=dbl_kind)
! Subprogram not used       if(my_task == master_task) then 
! Subprogram not used          ! compress out land blocks 
! Subprogram not used          ii =0
! Subprogram not used          do n=1,nblocks_tot
! Subprogram not used             if(nocn(n) > 0) then 
! Subprogram not used                ii = ii+1
! Subprogram not used                gnumIceCells2(ii) = gnumIceCells(n)
! Subprogram not used             endif
! Subprogram not used          enddo
! Subprogram not used           open(nu_timing,file='numCells.bin',recl=8*ii, &
! Subprogram not used             form = 'unformatted', access = 'direct', status = 'unknown')
! Subprogram not used           write(nu_timing,rec=1) gnumIceCells2(1:ii)
! Subprogram not used           close(nu_timing)
! Subprogram not used           print *,'numCells: ',gnumIceCells2(1:ii)
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       deallocate(gnumIceCells,gnumIceCells2)
! Subprogram not used     
! Subprogram not used 
! Subprogram not used    end subroutine print_numIceCells

! Subprogram not used    subroutine ReadProbabilityFile(distribution_wght_file,Prob)
! Subprogram not used    
! Subprogram not used       character(char_len_long), intent(in) :: distribution_wght_file
! Subprogram not used       real(real_kind), intent(inout)           :: Prob(:,:)
! Subprogram not used       
! Subprogram not used       type(block) :: this_block
! Subprogram not used       integer(int_kind) :: ilo,ihi
! Subprogram not used       integer(int_kind) :: jlo,jhi
! Subprogram not used 
! Subprogram not used       integer(int_kind) :: fid_prob
! Subprogram not used       integer(int_kind) :: amode,ncid,iostat,varid
! Subprogram not used       integer(int_kind) :: i,j,n,ncnt,ierr,ig,jg
! Subprogram not used       real(dbl_kind) :: val,sum,avg
! Subprogram not used 
! Subprogram not used       character (char_len) :: &
! Subprogram not used          fieldname       ! field name in netCDF file
! Subprogram not used 
! Subprogram not used 
! Subprogram not used       call ice_open_nc(distribution_wght_file,fid_prob)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used       fieldname = 'ice_present'
! Subprogram not used       call ice_read_global_nc(fid_prob,1,fieldname,Prob,.true.)
! Subprogram not used 
! Subprogram not used       if(my_task == master_task) then 
! Subprogram not used         call ice_close_nc(fid_prob)
! Subprogram not used         print *,'MAXVAL(Prob): ',MAXVAL(Prob) 
! Subprogram not used         print *,'MINVAL(Prob): ',MINVAL(Prob) 
! Subprogram not used         print *,'COUNT(Prob > 0.5): ',COUNT(Prob > 0.5)
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       ! normalize probability [0,1] 
! Subprogram not used !      prob=prob/100.0
! Subprogram not used       
! Subprogram not used    end subroutine ReadProbabilityFile

!***********************************************************************

!***********************************************************************
 subroutine CalcWorkPerBlock(distribution_wght, KMTG,ULATG,work_per_block, prob_per_block,blockType,bStats) 

   character (char_len), intent(in) :: distribution_wght

   real (dbl_kind), dimension(nx_global,ny_global), intent(in) :: &
      KMTG           ,&! global topography
      ULATG            ! global latitude field (radians)

   integer (int_kind), intent(inout)   :: &
      work_per_block(nblocks_tot)         ! number of work units per block

   real (dbl_kind), intent(inout) :: &
      prob_per_block(nblocks_tot)         ! probability of sea-ice in block

   integer (int_kind), intent(inout) ::  &
      blockType(nblocks_tot)              ! Type of block: 
                                          ! one of the following
                                          !    lnd,ice,icefree

   real (dbl_kind), intent(inout) :: bStats(:,:) 


   type (block) :: &
      this_block           ! block information for current block

   integer (int_kind), parameter :: &
      max_work_unit=10     ! quantize the work into values from 1,max

   integer(int_kind) :: jg,j,n,ig,i,work_unit
   integer(int_kind) :: maxkmt,mkmt,ilo,jlo,ihi,jhi
   integer(int_kind) :: activePts

   integer(int_kind) :: bid(max_blocks)
   integer(int_kind) :: bsize_x,bsize_y

   integer(int_kind), allocatable :: iLocation(:)

   real(dbl_kind), allocatable, dimension(:)   :: prob,work

   real(dbl_kind) :: lat
   real(dbl_kind), parameter :: w0 = 1., &  ! Constants for work blocks:
                          w1 = 10.   ! w0 for all blocks
                                      ! w1 for sea-ice blocks
   real(dbl_kind) :: maxlat,maxalat
   real(dbl_kind) :: shlatT,nhlatT
   real(real_kind) :: tmp

   integer(int_kind) :: ActiveBlocks,totWork

   logical, parameter :: Debug = .FALSE.

   integer(int_kind) :: numLnd,numIcefree,numIce

   integer(int_kind) :: ierr
   
!----------------------------------------------------------------------
!
!  estimate the amount of work per processor using the topography
!  and latitude
!
!----------------------------------------------------------------------



   allocate(work_gr(nx_global,ny_global)) 
   allocate(prob(nblocks_tot),work(nblocks_tot))
   allocate(nocn(nblocks_tot))
   allocate(nice005(nblocks_tot),nice010(nblocks_tot),nice050(nblocks_tot), &
	    nice100(nblocks_tot),nice250(nblocks_tot),nice500(nblocks_tot))

   nocn=0
   nice005=0
   nice010=0
   nice050=0
   nice100=0
   nice250=0
   nice500=0

   ActiveBlocks = 0
   select case (distribution_wght)
     case('file')
        call ReadProbabilityFile(distribution_wght_file,work_gr)
     case('erfc')
        work_gr = ErfcProbability(ULATG)
     case default
        work_gr = ErfcProbability(ULATG)
   end select
   !--------------------------------------------------
   ! It would be nice if we can not find the 
   ! the distribution_wgt_file,
   ! fall back to erfc weight function.
   ! However currently the error trapping does not prevent 
   ! such a fault recovery 
   !-------------------------------------------------------
   !if(ierr < 0) then 
   !    print *,'CalcWorkPerBlock: Could not open file: ',trim(distribution_wght_file)
   !    print *,'CalcWorkPerBlock:  Using ERFC probability function instead'
   !    work_gr = ErfcProbability(ULATG)
   !endif
       
     
!   print *,'IAM: ',my_task,'prod: ',prob

if(my_task == master_task) then 
   work_per_block=0
   do n=1,nblocks_tot

      this_block=get_block(n,n)
      ilo = this_block%ilo;ihi = this_block%ihi
      jlo = this_block%jlo;jhi = this_block%jhi
      !----------------------------------------------------
      ! calculate the probability of sea-ice in this block 
      !----------------------------------------------------
      tmp = 0.0
      do j=jlo,jhi
         jg=this_block%j_glob(j)
         if(jg>0) then  
            do i=ilo,ihi
              ig = this_block%i_glob(i)
              if(ig>0) then 
                 if(KMTG(ig,jg)>puny) then 
                    nocn(n) = nocn(n) + 1
                    if(work_gr(ig,jg) > 0.005) nice005(n) = nice005(n) + 1  
                    if(work_gr(ig,jg) > 0.010) nice010(n) = nice010(n) + 1  
                    if(work_gr(ig,jg) > 0.050) nice050(n) = nice050(n) + 1  
                    if(work_gr(ig,jg) > 0.100) nice100(n) = nice100(n) + 1  
                    if(work_gr(ig,jg) > 0.250) nice250(n) = nice250(n) + 1  
                    if(work_gr(ig,jg) > 0.500) nice500(n) = nice500(n) + 1  
                    tmp = tmp + work_gr(ig,jg)
                 endif
              endif
           enddo
         endif
      enddo
      if(nocn(n) > 0) then 
!         print *,'n:',n,' tmp:',tmp,' nocn(n): ',nocn(n)
         prob_per_block(n) = real(tmp,kind=dbl_kind)/real(nocn(n),kind=dbl_kind)
!         print *,'prob_per_block(n):' ,prob_per_block(n)
!         print *,' ihi,ilo,(ihi-ilo+1): ',ihi,ilo,(ihi-ilo+1)
!         print *,' jhi,jlo,(jhi-jlo+1): ',jhi,jlo,(jhi-jlo+1)
!         print *,' w0,w1: ',w0,w1
         work_per_block(n) =  &
           CEILING(w0 + (real(nocn(n),kind=dbl_kind)/real((ihi-ilo+1)*(jhi-jlo+1),kind=dbl_kind)) & 
                            *prob_per_block(n)*w1,kind=int_kind)
      else
         prob_per_block(n) = 0.0d0
      endif

      !--------------------------------------------------
      ! set the type of block (used latter for partition) 
      !--------------------------------------------------
      if(prob_per_block(n) >0.005 .and. nocn(n) > 0) then
        blockType(n) = iceType
      elseif (nocn(n) > 0) then 
        blockType(n) = icefreeType
      elseif (nocn(n) == 0) then
        blockType(n) = lndType
      endif

   enddo
 
   numIceFree = COUNT(blockType .eq. iceFreeType) 
   numLnd     = COUNT(blockType .eq. lndType)
   numIce     = COUNT(blockTYpe .eq. iceType)
    
!   print *,'Total blocks:',nblocks_tot,' land blocks: ',numLnd,' Ice blocks: ', &
!		numIce,' IceFree blocks: ',numIceFree
   write(*,23) nblocks_tot,numIce,numIceFree,numLnd 

23   format('CalcWorkPerBlock: Total blocks: ',i5,' Ice blocks: ',i5,' IceFree blocks: ',i5,' Land blocks: ',i5)

endif
   !-----------------------------------------------------------------
   ! broadcast info
   !-----------------------------------------------------------------
   call broadcast_array(work_per_block,master_task)
   call broadcast_array(prob_per_block,master_task)
   call broadcast_array(blockType,master_task)
   call broadcast_array(nocn,master_task)
   call broadcast_array(nice005,master_task)
   call broadcast_array(nice010,master_task)
   call broadcast_array(nice050,master_task)
   call broadcast_array(nice100,master_task)
   call broadcast_array(nice250,master_task)
   call broadcast_array(nice500,master_task)

   allocate(iLocation(nblocks_tot))
   do i=1,nblocks_tot
       iLocation(i) = i
   enddo

   call BuildProbabilityStats2(iLocation,bStats)

   deallocate(work_gr)
   deallocate(iLocation) 
!DBG   print *,'CalcWorkPerBlock: at the end of the subroutine' 
      
 end subroutine CalcWorkPerBlock

 function ErfcProbability(lat) result(prob)

      real (kind=dbl_kind), intent(in), &
         dimension (nx_global,ny_global) :: lat

      real (kind=real_kind),  &
	 dimension (nx_global,ny_global) :: prob


      real (kind=dbl_kind) :: ltmp,thetai,sigmai,arg
      integer :: i,j

      
      do j=1,ny_global
      do i=1,nx_global
         ltmp = lat(i,j)         
         if(ltmp > 0.0d0) then 
                   !---------------------
           ! northern latitude
           !---------------------
!JMD           thetai = (70.0_dbl_kind/180.0_dbl_kind)*pi
!JMD  For startup lots of random ice at low lattitudes
           thetai = (55.0_dbl_kind/180.0_dbl_kind)*pi
           sigmai = (5.0_dbl_kind/180.0_dbl_kind)*pi
           arg = (thetai-ltmp)/sigmai
        else 
           !---------------------
           ! southern hemisphere
           !---------------------
!JMD           thetai = (60.0_dbl_kind/180.0_dbl_kind)*pi
!JMD  For startup lots of random ice at low lattitudes
           thetai = (60.0_dbl_kind/180.0_dbl_kind)*pi
           sigmai = (5.0_dbl_kind/180.0_dbl_kind)*pi
           arg = (thetai+ltmp)/sigmai
	 endif
         prob(i,j) = real(erfc06(arg),kind=real_kind)/2.0_real_kind
      enddo
      enddo

 end function ErfcProbability 

 function erfc06(x)  result(erfc)

      implicit none
      real(dbl_kind), intent(in)     :: x
      real(dbl_kind)                 :: erfc
      real(dbl_kind)                 :: t,u
      real(dbl_kind), parameter      :: pa  =  3.97886080735226000e+00_dbl_kind
      real(dbl_kind), parameter      :: p00 =  2.75374741597376782e-01_dbl_kind
      real(dbl_kind), parameter      :: p01 =  4.90165080585318424e-01_dbl_kind
      real(dbl_kind), parameter      :: p02 =  7.74368199119538609e-01_dbl_kind
      real(dbl_kind), parameter      :: p03 =  1.07925515155856677e+00_dbl_kind
      real(dbl_kind), parameter      :: p04 =  1.31314653831023098e+00_dbl_kind
      real(dbl_kind), parameter      :: p05 =  1.37040217682338167e+00_dbl_kind
      real(dbl_kind), parameter      :: p06 =  1.18902982909273333e+00_dbl_kind
      real(dbl_kind), parameter      :: p07 =  8.05276408752910567e-01_dbl_kind
      real(dbl_kind), parameter      :: p08 =  3.57524274449531043e-01_dbl_kind
      real(dbl_kind), parameter      :: p09 =  1.66207924969367356e-02_dbl_kind
      real(dbl_kind), parameter      :: p10 = -1.19463959964325415e-01_dbl_kind
      real(dbl_kind), parameter      :: p11 = -8.38864557023001992e-02_dbl_kind
      real(dbl_kind), parameter      :: p12 =  2.49367200053503304e-03_dbl_kind
      real(dbl_kind), parameter      :: p13 =  3.90976845588484035e-02_dbl_kind
      real(dbl_kind), parameter      :: p14 =  1.61315329733252248e-02_dbl_kind
      real(dbl_kind), parameter      :: p15 = -1.33823644533460069e-02_dbl_kind
      real(dbl_kind), parameter      :: p16 = -1.27223813782122755e-02_dbl_kind
      real(dbl_kind), parameter      :: p17 =  3.83335126264887303e-03_dbl_kind
      real(dbl_kind), parameter      :: p18 =  7.73672528313526668e-03_dbl_kind
      real(dbl_kind), parameter      :: p19 = -8.70779635317295828e-04_dbl_kind
      real(dbl_kind), parameter      :: p20 = -3.96385097360513500e-03_dbl_kind
      real(dbl_kind), parameter      :: p21 =  1.19314022838340944e-04_dbl_kind
      real(dbl_kind), parameter      :: p22 =  1.27109764952614092e-03_dbl_kind




      t = pa/(pa + abs(x))
      u = t - 0.5_dbl_kind
      erfc = ((((((((((((((((((((((p22 * u + p21) * u + p20)   &
                                       * u + p19) * u + p18)   &
                                       * u + p17) * u + p16)   &
                                       * u + p15) * u + p14)   &
                                       * u + p13) * u + p12)   &
                                       * u + p11) * u + p10)   &
                                       * u + p09) * u + p08)   &
                                       * u + p07) * u + p06)   &
                                       * u + p05) * u + p04)   &
                                       * u + p03) * u + p02)   &
                                       * u + p01) * u + p00) * t * exp(-x**2)
      if (x < 0.0_dbl_kind) erfc = 2.0_dbl_kind - erfc

      end function erfc06

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

end module ice_probability
