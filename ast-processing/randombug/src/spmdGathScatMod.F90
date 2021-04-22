module spmdGathScatMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: spmdGathScatMod
!
! !DESCRIPTION:
! Perform SPMD gather and scatter operations.
!
! !USES:
  use clm_varcon, only: spval, ispval
  use decompMod, only : get_clmlevel_gsmap
  use shr_kind_mod, only: r8 => shr_kind_r8
  use spmdMod
  use mct_mod
  use abortutils, only : endrun
  use clm_varctl, only : iulog
  use perf_mod
!
! !PUBLIC TYPES:
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
  public  scatter_data_from_master, gather_data_to_master

  interface scatter_data_from_master
     module procedure scatter_1darray_int
     module procedure scatter_1darray_real
  end interface

  interface gather_data_to_master
     module procedure gather_1darray_int
     module procedure gather_1darray_real
  end interface
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
  integer,private,parameter :: debug = 0

!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: scatter_1darray_int
!
! !INTERFACE:
  subroutine scatter_1darray_int (alocal, aglobal, clmlevel)
!
! !DESCRIPTION:
! Wrapper routine to scatter int 1d array
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    integer , pointer            :: alocal(:)       ! local  data (output)
    integer , pointer            :: aglobal(:)      ! global data (input)
    character(len=*) ,intent(in) :: clmlevel    ! type of input grid
!
! !REVISION HISTORY:
! Author: T Craig
!
!
! !LOCAL VARIABLES:
!EOP
    integer            :: n1,n2,lb1,ub1,lb2,ub2 ! indices
    integer            :: lsize      ! size of local array
    type(mct_aVect)    :: AVi, AVo   ! attribute vectors
    integer ,pointer   :: adata(:)   ! local data array
    character(len=256) :: rstring    ! real field list string
    character(len=256) :: istring    ! int field list string
    character(len=8)   :: fname      ! arbitrary field name
    type(mct_gsMap),pointer       :: gsmap   ! global seg map
    character(len=*),parameter :: subname = 'scatter_1darray_int'

!-----------------------------------------------------------------------

    call t_startf(trim(subname)//'_total')
    call get_clmlevel_gsmap(clmlevel,gsmap)

    lb1 = lbound(alocal,dim=1)
    ub1 = ubound(alocal,dim=1)
    lb2 = 1
    ub2 = 1

    rstring = ""
    istring = ""

    do n2 = lb2,ub2
       write(fname,'(a1,i3.3)') 'f',n2-lb2+1
       if (len_trim(istring) == 0) then
          istring = trim(fname)
       else
          istring = trim(istring)//":"//trim(fname)
       endif
    enddo

    if (masterproc .and. debug > 2) then
       write(iulog,*) trim(subname),' strings:',trim(rstring),' ',trim(istring)
    endif

    if (debug > 1) call t_startf(trim(subname)//'_pack')

    if (masterproc) then
       lsize = size(aglobal,dim=1)
       call mct_aVect_init(AVi,rList=trim(rstring),iList=trim(istring),lsize=lsize)
       allocate(adata(lsize))
       do n2 = lb2,ub2
          adata(1:lsize) = aglobal(1:lsize)
          write(fname,'(a1,i3.3)') 'f',n2-lb2+1
          call mct_aVect_importIattr(AVi,trim(fname),adata,lsize)
       enddo
       deallocate(adata)
    endif

    if (debug > 1) call t_stopf(trim(subname)//'_pack')
    if (debug > 1) call t_startf(trim(subname)//'_scat')

    call mct_aVect_scatter(AVi, AVo, gsmap, 0, mpicom)

    if (debug > 1) call t_stopf(trim(subname)//'_scat')
    if (debug > 1) call t_startf(trim(subname)//'_upck')

    lsize = size(alocal,dim=1)
    allocate(adata(lsize))
    do n2 = lb2,ub2
       write(fname,'(a1,i3.3)') 'f',n2-lb2+1
       call mct_aVect_exportIattr(AVo,trim(fname),adata,lsize)
       do n1 = lb1,ub1
          alocal(n1) = adata(n1-lb1+1)
       enddo
    enddo
    deallocate(adata)

    if (debug > 1) call t_stopf(trim(subname)//'_upck')

    if (masterproc) then
       call mct_aVect_clean(AVi)
    endif
    call mct_aVect_clean(AVo)

    call t_stopf(trim(subname)//'_total')

  end subroutine scatter_1darray_int

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: gather_1darray_int
!
! !INTERFACE:
  subroutine gather_1darray_int (alocal, aglobal, clmlevel, missing)
!
! !DESCRIPTION:
! Wrapper routine to gather int 1d array
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    integer , pointer            :: alocal(:)       ! local  data (output)
    integer , pointer            :: aglobal(:)      ! global data (input)
    character(len=*) ,intent(in) :: clmlevel    ! type of input grid
    integer ,optional,intent(in) :: missing     ! missing value
!
! !REVISION HISTORY:
! Author: T Craig
!
!
! !LOCAL VARIABLES:
!EOP
    integer            :: n1,n2,lb1,ub1,lb2,ub2 ! indices
    integer            :: lsize      ! size of local array
    type(mct_aVect)    :: AVi, AVo   ! attribute vectors
    integer ,pointer   :: adata(:)   ! temporary data array
    integer ,pointer   :: mvect(:)   ! local array for mask
    character(len=256) :: rstring    ! real field list string
    character(len=256) :: istring    ! int field list string
    character(len=8)   :: fname      ! arbitrary field name
    type(mct_gsMap),pointer       :: gsmap   ! global seg map
    character(len=*),parameter :: subname = 'gather_1darray_int'

!-----------------------------------------------------------------------

    call t_startf(trim(subname)//'_total')
    call get_clmlevel_gsmap(clmlevel,gsmap)

    lsize = size(alocal,dim=1)
    lb1 = lbound(alocal,dim=1)
    ub1 = ubound(alocal,dim=1)
    lb2 = 1
    ub2 = 1
   
    rstring = ""
    istring = ""

    if (present(missing)) then
       istring = "mask"
    endif

    do n2 = lb2,ub2
       write(fname,'(a1,i3.3)') 'f',n2-lb2+1
       if (len_trim(istring) == 0) then
          istring = trim(fname)
       else
          istring = trim(istring)//":"//trim(fname)
       endif
    enddo

    if (masterproc .and. debug > 2) then
       write(iulog,*) trim(subname),' strings:',trim(rstring),' ',trim(istring)
    endif

    call mct_aVect_init(AVi,rList=trim(rstring),iList=trim(istring),lsize=lsize)

    if (debug > 1) call t_startf(trim(subname)//'_pack')
    allocate(adata(lsize))
    do n2 = lb2,ub2
       do n1 = lb1,ub1
          adata(n1-lb1+1) = alocal(n1)
       enddo
       write(fname,'(a1,i3.3)') 'f',n2-lb2+1
       call mct_aVect_importIattr(AVi,trim(fname),adata,lsize)
    enddo
    deallocate(adata)

    if (present(missing)) then
       allocate(mvect(lsize))
       do n1 = lb1,ub1
          mvect(n1-lb1+1) = 1
       enddo
       call mct_aVect_importIattr(AVi,"mask",mvect,lsize)
       deallocate(mvect)
    endif

    if (debug > 1) call t_stopf(trim(subname)//'_pack')
    if (debug > 1) call t_startf(trim(subname)//'_gath')

    if (present(missing)) then
! tcx wait for update in mct, then get rid of "mask"
!       call mct_aVect_gather(AVi, AVo, gsmap, 0, mpicom, missing = missing)
       call mct_aVect_gather(AVi, AVo, gsmap, 0, mpicom)
    else
       call mct_aVect_gather(AVi, AVo, gsmap, 0, mpicom)
    endif

    if (debug > 1) call t_stopf(trim(subname)//'_gath')
    if (debug > 1) call t_startf(trim(subname)//'_upck')

    if (masterproc) then
       lsize = size(aglobal,dim=1)
       allocate(adata(lsize))
       do n2 = lb2,ub2
          write(fname,'(a1,i3.3)') 'f',n2-lb2+1
          call mct_aVect_exportIattr(AVo,trim(fname),adata,lsize)
          aglobal(1:lsize) = adata(1:lsize)
       enddo
       deallocate(adata)
       if (present(missing)) then
          allocate(mvect(lsize))
          call mct_aVect_exportIattr(AVo,"mask",mvect,lsize)
          do n1 = 1,lsize
             if (mvect(n1) == 0) then
                do n2 = lb2,ub2
                   aglobal(n1) = missing
                enddo
             endif
          enddo
          deallocate(mvect)
       endif
    endif

    if (debug > 1) call t_stopf(trim(subname)//'_upck')

    if (masterproc) then
       call mct_aVect_clean(AVo)
    endif

    call mct_aVect_clean(AVi)

    call t_stopf(trim(subname)//'_total')

  end subroutine gather_1darray_int

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: scatter_1darray_real
!
! !INTERFACE:
! Subprogram not used   subroutine scatter_1darray_real (alocal, aglobal, clmlevel)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used ! Wrapper routine to scatter real 1d array
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     real(r8), pointer            :: alocal(:)       ! local  data (output)
! Subprogram not used     real(r8), pointer            :: aglobal(:)      ! global data (input)
! Subprogram not used     character(len=*) ,intent(in) :: clmlevel    ! type of input grid
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! Author: T Craig
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! !LOCAL VARIABLES:
! Subprogram not used !EOP
! Subprogram not used     integer            :: n1,n2,lb1,ub1,lb2,ub2 ! indices
! Subprogram not used     integer            :: lsize      ! size of local array
! Subprogram not used     type(mct_aVect)    :: AVi, AVo   ! attribute vectors
! Subprogram not used     real(r8),pointer   :: adata(:)   ! local data array
! Subprogram not used     character(len=256) :: rstring    ! real field list string
! Subprogram not used     character(len=256) :: istring    ! int field list string
! Subprogram not used     character(len=8)   :: fname      ! arbitrary field name
! Subprogram not used     type(mct_gsMap),pointer       :: gsmap   ! global seg map
! Subprogram not used     character(len=*),parameter :: subname = 'scatter_1darray_real'
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     call t_startf(trim(subname)//'_total')
! Subprogram not used     call get_clmlevel_gsmap(clmlevel,gsmap)
! Subprogram not used 
! Subprogram not used     lb1 = lbound(alocal,dim=1)
! Subprogram not used     ub1 = ubound(alocal,dim=1)
! Subprogram not used     lb2 = 1
! Subprogram not used     ub2 = 1
! Subprogram not used 
! Subprogram not used     rstring = ""
! Subprogram not used     istring = ""
! Subprogram not used 
! Subprogram not used     do n2 = lb2,ub2
! Subprogram not used        write(fname,'(a1,i3.3)') 'f',n2-lb2+1
! Subprogram not used        if (len_trim(rstring) == 0) then
! Subprogram not used           rstring = trim(fname)
! Subprogram not used        else
! Subprogram not used           rstring = trim(rstring)//":"//trim(fname)
! Subprogram not used        endif
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     if (masterproc .and. debug > 2) then
! Subprogram not used        write(iulog,*) trim(subname),' strings:',trim(rstring),' ',trim(istring)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if (debug > 1) call t_startf(trim(subname)//'_pack')
! Subprogram not used 
! Subprogram not used     if (masterproc) then
! Subprogram not used        lsize = size(aglobal,dim=1)
! Subprogram not used        call mct_aVect_init(AVi,rList=trim(rstring),iList=trim(istring),lsize=lsize)
! Subprogram not used        allocate(adata(lsize))
! Subprogram not used        do n2 = lb2,ub2
! Subprogram not used           adata(1:lsize) = aglobal(1:lsize)
! Subprogram not used           write(fname,'(a1,i3.3)') 'f',n2-lb2+1
! Subprogram not used           call mct_aVect_importRattr(AVi,trim(fname),adata,lsize)
! Subprogram not used        enddo
! Subprogram not used        deallocate(adata)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if (debug > 1) call t_stopf(trim(subname)//'_pack')
! Subprogram not used     if (debug > 1) call t_startf(trim(subname)//'_scat')
! Subprogram not used 
! Subprogram not used     call mct_aVect_scatter(AVi, AVo, gsmap, 0, mpicom)
! Subprogram not used 
! Subprogram not used     if (debug > 1) call t_stopf(trim(subname)//'_scat')
! Subprogram not used     if (debug > 1) call t_startf(trim(subname)//'_upck')
! Subprogram not used 
! Subprogram not used     lsize = size(alocal,dim=1)
! Subprogram not used     allocate(adata(lsize))
! Subprogram not used     do n2 = lb2,ub2
! Subprogram not used        write(fname,'(a1,i3.3)') 'f',n2-lb2+1
! Subprogram not used        call mct_aVect_exportRattr(AVo,trim(fname),adata,lsize)
! Subprogram not used        do n1 = lb1,ub1
! Subprogram not used           alocal(n1) = adata(n1-lb1+1)
! Subprogram not used        enddo
! Subprogram not used     enddo
! Subprogram not used     deallocate(adata)
! Subprogram not used 
! Subprogram not used     if (debug > 1) call t_stopf(trim(subname)//'_upck')
! Subprogram not used 
! Subprogram not used     if (masterproc) then
! Subprogram not used        call mct_aVect_clean(AVi)
! Subprogram not used     endif
! Subprogram not used     call mct_aVect_clean(AVo)
! Subprogram not used 
! Subprogram not used     call t_stopf(trim(subname)//'_total')
! Subprogram not used 
! Subprogram not used   end subroutine scatter_1darray_real

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: gather_1darray_real
!
! !INTERFACE:
! Subprogram not used   subroutine gather_1darray_real (alocal, aglobal, clmlevel, missing)
! Subprogram not used !
! Subprogram not used ! !DESCRIPTION:
! Subprogram not used ! Wrapper routine to gather real 1d array
! Subprogram not used !
! Subprogram not used ! !USES:
! Subprogram not used !
! Subprogram not used ! !ARGUMENTS:
! Subprogram not used     implicit none
! Subprogram not used     real(r8), pointer            :: alocal(:)       ! local  data (output)
! Subprogram not used     real(r8), pointer            :: aglobal(:)      ! global data (input)
! Subprogram not used     character(len=*) ,intent(in) :: clmlevel    ! type of input grid
! Subprogram not used     real(r8),optional,intent(in) :: missing     ! missing value
! Subprogram not used !
! Subprogram not used ! !REVISION HISTORY:
! Subprogram not used ! Author: T Craig
! Subprogram not used !
! Subprogram not used !
! Subprogram not used ! !LOCAL VARIABLES:
! Subprogram not used !EOP
! Subprogram not used     integer            :: n1,n2,lb1,ub1,lb2,ub2 ! indices
! Subprogram not used     integer            :: lsize      ! size of local array
! Subprogram not used     type(mct_aVect)    :: AVi, AVo   ! attribute vectors
! Subprogram not used     real(r8),pointer   :: adata(:)   ! temporary data array
! Subprogram not used     integer ,pointer   :: mvect(:)   ! local array for mask
! Subprogram not used     character(len=256) :: rstring    ! real field list string
! Subprogram not used     character(len=256) :: istring    ! int field list string
! Subprogram not used     character(len=8)   :: fname      ! arbitrary field name
! Subprogram not used     type(mct_gsMap),pointer       :: gsmap   ! global seg map
! Subprogram not used     character(len=*),parameter :: subname = 'gather_1darray_real'
! Subprogram not used 
! Subprogram not used !-----------------------------------------------------------------------
! Subprogram not used 
! Subprogram not used     call t_startf(trim(subname)//'_total')
! Subprogram not used     call get_clmlevel_gsmap(clmlevel,gsmap)
! Subprogram not used 
! Subprogram not used     lsize = size(alocal,dim=1)
! Subprogram not used     lb1 = lbound(alocal,dim=1)
! Subprogram not used     ub1 = ubound(alocal,dim=1)
! Subprogram not used     lb2 = 1
! Subprogram not used     ub2 = 1
! Subprogram not used    
! Subprogram not used     rstring = ""
! Subprogram not used     istring = ""
! Subprogram not used 
! Subprogram not used     if (present(missing)) then
! Subprogram not used        istring = "mask"
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     do n2 = lb2,ub2
! Subprogram not used        write(fname,'(a1,i3.3)') 'f',n2-lb2+1
! Subprogram not used        if (len_trim(rstring) == 0) then
! Subprogram not used           rstring = trim(fname)
! Subprogram not used        else
! Subprogram not used           rstring = trim(rstring)//":"//trim(fname)
! Subprogram not used        endif
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     if (masterproc .and. debug > 2) then
! Subprogram not used        write(iulog,*) trim(subname),' strings:',trim(rstring),' ',trim(istring)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     call mct_aVect_init(AVi,rList=trim(rstring),iList=trim(istring),lsize=lsize)
! Subprogram not used 
! Subprogram not used     if (debug > 1) call t_startf(trim(subname)//'_pack')
! Subprogram not used     allocate(adata(lsize))
! Subprogram not used     do n2 = lb2,ub2
! Subprogram not used        do n1 = lb1,ub1
! Subprogram not used           adata(n1-lb1+1) = alocal(n1)
! Subprogram not used        enddo
! Subprogram not used        write(fname,'(a1,i3.3)') 'f',n2-lb2+1
! Subprogram not used        call mct_aVect_importRattr(AVi,trim(fname),adata,lsize)
! Subprogram not used     enddo
! Subprogram not used     deallocate(adata)
! Subprogram not used 
! Subprogram not used     if (present(missing)) then
! Subprogram not used        allocate(mvect(lsize))
! Subprogram not used        do n1 = lb1,ub1
! Subprogram not used           mvect(n1-lb1+1) = 1
! Subprogram not used        enddo
! Subprogram not used        call mct_aVect_importIattr(AVi,"mask",mvect,lsize)
! Subprogram not used        deallocate(mvect)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if (debug > 1) call t_stopf(trim(subname)//'_pack')
! Subprogram not used     if (debug > 1) call t_startf(trim(subname)//'_gath')
! Subprogram not used 
! Subprogram not used     if (present(missing)) then
! Subprogram not used ! tcx wait for update in mct, then get rid of "mask"
! Subprogram not used !       call mct_aVect_gather(AVi, AVo, gsmap, 0, mpicom, missing = missing)
! Subprogram not used        call mct_aVect_gather(AVi, AVo, gsmap, 0, mpicom)
! Subprogram not used     else
! Subprogram not used        call mct_aVect_gather(AVi, AVo, gsmap, 0, mpicom)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if (debug > 1) call t_stopf(trim(subname)//'_gath')
! Subprogram not used     if (debug > 1) call t_startf(trim(subname)//'_upck')
! Subprogram not used 
! Subprogram not used     if (masterproc) then
! Subprogram not used        lsize = size(aglobal,dim=1)
! Subprogram not used        allocate(adata(lsize))
! Subprogram not used        do n2 = lb2,ub2
! Subprogram not used           write(fname,'(a1,i3.3)') 'f',n2-lb2+1
! Subprogram not used           call mct_aVect_exportRattr(AVo,trim(fname),adata,lsize)
! Subprogram not used           aglobal(1:lsize) = adata(1:lsize)
! Subprogram not used        enddo
! Subprogram not used        deallocate(adata)
! Subprogram not used        if (present(missing)) then
! Subprogram not used           allocate(mvect(lsize))
! Subprogram not used           call mct_aVect_exportIattr(AVo,"mask",mvect,lsize)
! Subprogram not used           do n1 = 1,lsize
! Subprogram not used              if (mvect(n1) == 0) then
! Subprogram not used                 do n2 = lb2,ub2
! Subprogram not used                    aglobal(n1) = missing
! Subprogram not used                 enddo
! Subprogram not used              endif
! Subprogram not used           enddo
! Subprogram not used           deallocate(mvect)
! Subprogram not used        endif
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     if (debug > 1) call t_stopf(trim(subname)//'_upck')
! Subprogram not used 
! Subprogram not used     if (masterproc) then
! Subprogram not used        call mct_aVect_clean(AVo)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     call mct_aVect_clean(AVi)
! Subprogram not used 
! Subprogram not used     call t_stopf(trim(subname)//'_total')
! Subprogram not used 
! Subprogram not used   end subroutine gather_1darray_real

end module spmdGathScatMod
