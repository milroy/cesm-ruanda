module seq_map_mod

!---------------------------------------------------------------------
!
! Purpose:
!
! General mapping routines
! including self normalizing mapping routine with optional fraction
!       
! Author: T. Craig, Jan-28-2011
!
!---------------------------------------------------------------------

  use shr_kind_mod      ,only: R8 => SHR_KIND_R8, IN=>SHR_KIND_IN
  use shr_kind_mod      ,only: CL => SHR_KIND_CL, CX => SHR_KIND_CX
  use shr_sys_mod
  use shr_const_mod
  use shr_mct_mod, only: shr_mct_sMatPInitnc, shr_mct_queryConfigFile
  use mct_mod
  use seq_comm_mct





  use component_type_mod
  use seq_map_type_mod

  implicit none
  save
  private  ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: seq_map_init_rcfile     ! cpl pes
  public :: seq_map_init_rearrolap  ! cpl pes
  public :: seq_map_initvect        ! cpl pes
  public :: seq_map_map             ! cpl pes 
  public :: seq_map_mapvect         ! cpl pes
  public :: seq_map_readdata        ! cpl pes




  interface seq_map_avNorm 
     module procedure seq_map_avNormArr
     module procedure seq_map_avNormAvF
  end interface

!--------------------------------------------------------------------------
! Public data
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

  character(*),parameter :: seq_map_stroff = 'variable_unset'
  character(*),parameter :: seq_map_stron  = 'StrinG_is_ON'
  real(R8),parameter,private :: deg2rad = shr_const_pi/180.0_R8  ! deg to rads

!=======================================================================
contains
!=======================================================================


 !===============================================================================

  subroutine seq_map_init_rcfile( mapper, comp_s, comp_d, &
       maprcfile, maprcname, maprctype, samegrid, string, esmf_map)

    implicit none
    !-----------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_map)        ,intent(inout),pointer :: mapper
    type(component_type) ,intent(inout)         :: comp_s 
    type(component_type) ,intent(inout)         :: comp_d 
    character(len=*)     ,intent(in)            :: maprcfile
    character(len=*)     ,intent(in)            :: maprcname
    character(len=*)     ,intent(in)            :: maprctype
    logical              ,intent(in)            :: samegrid
    character(len=*)     ,intent(in),optional   :: string
    logical              ,intent(in),optional   :: esmf_map
    !
    ! Local Variables
    !
    type(mct_gsmap), pointer    :: gsmap_s ! temporary pointers
    type(mct_gsmap), pointer    :: gsmap_d ! temporary pointers
    integer(IN)                 :: mpicom
    character(CX)               :: mapfile
    character(CL)               :: maptype
    integer(IN)                 :: mapid
    integer(IN)                 :: ssize,dsize
    character(len=*),parameter  :: subname = "(seq_map_init_rcfile) "
    !-----------------------------------------------------

    if (seq_comm_iamroot(CPLID) .and. present(string)) then
       write(logunit,'(A)') subname//' called for '//trim(string)
    endif

    call seq_comm_setptrs(CPLID, mpicom=mpicom)

    gsmap_s => component_get_gsmap_cx(comp_s)
    gsmap_d => component_get_gsmap_cx(comp_d)

    if (mct_gsmap_Identical(gsmap_s,gsmap_d)) then
       call seq_map_mapmatch(mapid,gsmap_s=gsmap_s,gsmap_d=gsmap_d,strategy="copy")

       if (mapid > 0) then
          call seq_map_mappoint(mapid,mapper)
       else
          call seq_map_mapinit(mapper,mpicom)
          mapper%copy_only = .true.
          mapper%strategy = "copy"
          mapper%gsmap_s => component_get_gsmap_cx(comp_s)
          mapper%gsmap_d => component_get_gsmap_cx(comp_d)
       endif

    elseif (samegrid) then
       call seq_map_mapmatch(mapid,gsmap_s=gsmap_s,gsmap_d=gsmap_d,strategy="rearrange")

       if (mapid > 0) then
          call seq_map_mappoint(mapid,mapper)
       else
          ! --- Initialize rearranger
          call seq_map_mapinit(mapper,mpicom)
          mapper%rearrange_only = .true.
          mapper%strategy = "rearrange"
          mapper%gsmap_s => component_get_gsmap_cx(comp_s)
          mapper%gsmap_d => component_get_gsmap_cx(comp_d)
          call seq_map_gsmapcheck(gsmap_s, gsmap_d)
          call mct_rearr_init(gsmap_s, gsmap_d, mpicom, mapper%rearr)
       endif

    else

       ! --- Initialize Smatp
       call shr_mct_queryConfigFile(mpicom,maprcfile,maprcname,mapfile,maprctype,maptype)

       call seq_map_mapmatch(mapid,gsMap_s=gsMap_s,gsMap_d=gsMap_d,mapfile=mapfile,strategy=maptype)

       if (mapid > 0) then
          call seq_map_mappoint(mapid,mapper)
       else
          call seq_map_mapinit(mapper,mpicom)
          mapper%mapfile = trim(mapfile)
          mapper%strategy= trim(maptype)
          mapper%gsmap_s => component_get_gsmap_cx(comp_s)
          mapper%gsmap_d => component_get_gsmap_cx(comp_d)

          call shr_mct_sMatPInitnc(mapper%sMatp, mapper%gsMap_s, mapper%gsMap_d, trim(mapfile),trim(maptype),mpicom)
          if (present(esmf_map)) mapper%esmf_map = esmf_map

          if (mapper%esmf_map) then
             call shr_sys_abort(subname//' ERROR: esmf SMM not allowed without USE_ESMF_LIB')
          endif  ! esmf_map

       endif  ! mapid >= 0
    endif

    if (seq_comm_iamroot(CPLID)) then
       write(logunit,'(2A,I6,4A)') subname,' mapper counter, strategy, mapfile = ', &
          mapper%counter,' ',trim(mapper%strategy),' ',trim(mapper%mapfile)
       call shr_sys_flush(logunit)
    endif

  end subroutine seq_map_init_rcfile

  !=======================================================================

  subroutine seq_map_init_rearrolap(mapper, comp_s, comp_d, string)

    implicit none
    !-----------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_map)        ,intent(inout),pointer :: mapper
    type(component_type) ,intent(inout)         :: comp_s 
    type(component_type) ,intent(inout)         :: comp_d 
    character(len=*)     ,intent(in),optional   :: string
    !
    ! Local Variables
    !
    integer(IN)                :: mapid
    type(mct_gsmap), pointer   :: gsmap_s
    type(mct_gsmap), pointer   :: gsmap_d
    integer(IN)                :: mpicom
    character(len=*),parameter :: subname = "(seq_map_init_rearrolap) "
    !-----------------------------------------------------

    if (seq_comm_iamroot(CPLID) .and. present(string)) then
       write(logunit,'(A)') subname//' called for '//trim(string)
    endif

    call seq_comm_setptrs(CPLID, mpicom=mpicom)

    gsmap_s => component_get_gsmap_cx(comp_s)
    gsmap_d => component_get_gsmap_cx(comp_d)

    if (mct_gsmap_Identical(gsmap_s,gsmap_d)) then
       call seq_map_mapmatch(mapid,gsmap_s=gsmap_s,gsmap_d=gsmap_d,strategy="copy")

       if (mapid > 0) then
          call seq_map_mappoint(mapid,mapper)
       else
          call seq_map_mapinit(mapper,mpicom)
          mapper%copy_only = .true.
          mapper%strategy = "copy"
          mapper%gsmap_s => component_get_gsmap_cx(comp_s)
          mapper%gsmap_d => component_get_gsmap_cx(comp_d)
       endif

    else
       call seq_map_mapmatch(mapid,gsmap_s=gsmap_s,gsmap_d=gsmap_d,strategy="rearrange")

       if (mapid > 0) then
          call seq_map_mappoint(mapid,mapper)
       else
          ! --- Initialize rearranger
          call seq_map_mapinit(mapper, mpicom)
          mapper%rearrange_only = .true.
          mapper%strategy = "rearrange"
          mapper%gsmap_s => component_get_gsmap_cx(comp_s)
          mapper%gsmap_d => component_get_gsmap_cx(comp_d)
          call seq_map_gsmapcheck(gsmap_s, gsmap_d)
          call mct_rearr_init(gsmap_s, gsmap_d, mpicom, mapper%rearr)
       endif

    endif

    if (seq_comm_iamroot(CPLID)) then
       write(logunit,'(2A,I6,4A)') subname,' mapper counter, strategy, mapfile = ', &
          mapper%counter,' ',trim(mapper%strategy),' ',trim(mapper%mapfile)
       call shr_sys_flush(logunit)
    endif

  end subroutine seq_map_init_rearrolap

  !=======================================================================

  subroutine seq_map_map( mapper, av_s, av_d, fldlist, norm, avwts_s, avwtsfld_s, &
                          string, msgtag )

    implicit none
    !-----------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_map)   ,intent(inout)       :: mapper
    type(mct_aVect) ,intent(in)          :: av_s
    type(mct_aVect) ,intent(inout)       :: av_d
    character(len=*),intent(in),optional :: fldlist
    logical         ,intent(in),optional :: norm
    type(mct_aVect) ,intent(in),optional :: avwts_s
    character(len=*),intent(in),optional :: avwtsfld_s
    character(len=*),intent(in),optional :: string
    integer(IN)     ,intent(in),optional :: msgtag
    !
    ! Local Variables
    !
    logical :: lnorm
    integer(IN),save :: ltag    ! message tag for rearrange
    character(len=*),parameter :: subname = "(seq_map_map) "
    !-----------------------------------------------------

    if (seq_comm_iamroot(CPLID) .and. present(string)) then
       write(logunit,'(A)') subname//' called for '//trim(string)
    endif

    lnorm = .true.
    if (present(norm)) then
       lnorm = norm
    endif

    if (present(msgtag)) then
       ltag = msgtag
    else
       ltag = 2000
    endif

    if (present(avwts_s) .and. .not. present(avwtsfld_s)) then
       write(logunit,*) subname,' ERROR: avwts_s present but avwtsfld_s not'
       call shr_sys_abort(subname//' ERROR: avwts present')
    endif
    if (.not. present(avwts_s) .and. present(avwtsfld_s)) then
       write(logunit,*) subname,' ERROR: avwtsfld_s present but avwts_s not'
       call shr_sys_abort(subname//' ERROR: avwtsfld present')
    endif

    if (mapper%copy_only) then
       !-------------------------------------------
       ! COPY data
       !-------------------------------------------
       if (present(fldlist)) then
          call mct_aVect_copy(aVin=av_s,aVout=av_d,rList=fldlist,vector=mct_usevector)
       else
          call mct_aVect_copy(aVin=av_s,aVout=av_d,vector=mct_usevector)
       endif

    else if (mapper%rearrange_only) then
       !-------------------------------------------
       ! REARRANGE data
       !-------------------------------------------
       if (present(fldlist)) then
          call mct_rearr_rearrange_fldlist(av_s, av_d, mapper%rearr, tag=ltag, VECTOR=mct_usevector, &
               ALLTOALL=mct_usealltoall, fldlist=fldlist)
       else
          call mct_rearr_rearrange(av_s, av_d, mapper%rearr, tag=ltag, VECTOR=mct_usevector, &
               ALLTOALL=mct_usealltoall)
       endif

    else
       !-------------------------------------------
       ! MAP data
       !-------------------------------------------
       if (present(avwts_s)) then
          if (present(fldlist)) then
             call seq_map_avNorm(mapper, av_s, av_d, avwts_s, trim(avwtsfld_s), &
                  rList=fldlist, norm=lnorm)
          else
             call seq_map_avNorm(mapper, av_s, av_d, avwts_s, trim(avwtsfld_s), &
                  norm=lnorm)
          endif
       else
          if (present(fldlist)) then
             call seq_map_avNorm(mapper, av_s, av_d, rList=fldlist, norm=lnorm)
          else
             call seq_map_avNorm(mapper, av_s, av_d, norm=lnorm)
          endif
       endif
    end if

  end subroutine seq_map_map

  !=======================================================================

  subroutine seq_map_initvect(mapper, type, comp_s, comp_d, string)

    !-----------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_map)        ,intent(inout)       :: mapper
    character(len=*)     ,intent(in)          :: type
    type(component_type) ,intent(inout)       :: comp_s 
    type(component_type) ,intent(inout)       :: comp_d 
    character(len=*)     ,intent(in),optional :: string
    !
    ! Local Variables
    !
    type(mct_gGrid), pointer   :: dom_s
    type(mct_gGrid), pointer   :: dom_d
    integer(IN)                :: klon, klat, lsize, n
    logical                    :: lnorm
    character(len=CL)          :: lstring
    character(len=*),parameter :: subname = "(seq_map_initvect) "
    !-----------------------------------------------------

    lstring = ' '
    if (present(string)) then
       if (seq_comm_iamroot(CPLID)) write(logunit,'(A)') subname//' called for '//trim(string)
       lstring = trim(string)
    endif

    dom_s => component_get_dom_cx(comp_s)
    dom_d => component_get_dom_cx(comp_d)

    if (trim(type(1:6)) == 'cart3d') then
       if (mapper%cart3d_init == trim(seq_map_stron)) return

       !--- compute these up front for vector mapping ---
       lsize = mct_aVect_lsize(dom_s%data)
       allocate(mapper%slon_s(lsize),mapper%clon_s(lsize), &
                mapper%slat_s(lsize),mapper%clat_s(lsize))
       klon = mct_aVect_indexRa(dom_s%data, "lon" )
       klat = mct_aVect_indexRa(dom_s%data, "lat" )
       do n = 1,lsize
          mapper%slon_s(n) = sin(dom_s%data%rAttr(klon,n)*deg2rad)
          mapper%clon_s(n) = cos(dom_s%data%rAttr(klon,n)*deg2rad)
          mapper%slat_s(n) = sin(dom_s%data%rAttr(klat,n)*deg2rad)
          mapper%clat_s(n) = cos(dom_s%data%rAttr(klat,n)*deg2rad)
       enddo
       
       lsize = mct_aVect_lsize(dom_d%data)
       allocate(mapper%slon_d(lsize),mapper%clon_d(lsize), &
                mapper%slat_d(lsize),mapper%clat_d(lsize))
       klon = mct_aVect_indexRa(dom_d%data, "lon" )
       klat = mct_aVect_indexRa(dom_d%data, "lat" )
       do n = 1,lsize
          mapper%slon_d(n) = sin(dom_d%data%rAttr(klon,n)*deg2rad)
          mapper%clon_d(n) = cos(dom_d%data%rAttr(klon,n)*deg2rad)
          mapper%slat_d(n) = sin(dom_d%data%rAttr(klat,n)*deg2rad)
          mapper%clat_d(n) = cos(dom_d%data%rAttr(klat,n)*deg2rad)
       enddo
       mapper%cart3d_init = trim(seq_map_stron)
    endif

  end subroutine seq_map_initvect

  !=======================================================================

  subroutine seq_map_mapvect( mapper, type, av_s, av_d, fldu, fldv, norm, string )

    implicit none
    !-----------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_map)   ,intent(inout)       :: mapper
    character(len=*),intent(in)          :: type
    type(mct_aVect) ,intent(in)          :: av_s
    type(mct_aVect) ,intent(inout)       :: av_d
    character(len=*),intent(in)          :: fldu
    character(len=*),intent(in)          :: fldv
    logical         ,intent(in),optional :: norm
    character(len=*),intent(in),optional :: string
    !
    ! Local Variables
    !
    logical :: lnorm
    character(len=CL) :: lstring
    character(len=*),parameter :: subname = "(seq_map_mapvect) "
    !-----------------------------------------------------

    lstring = ' '
    if (present(string)) then
       if (seq_comm_iamroot(CPLID)) write(logunit,'(A)') subname//' called for '//trim(string)
       lstring = trim(string)
    endif

    if (mapper%copy_only .or. mapper%rearrange_only) then
       return
    endif

    lnorm = .true.
    if (present(norm)) then
       lnorm = norm
    endif

    if (trim(type(1:6)) == 'cart3d') then
       if (mapper%cart3d_init /= trim(seq_map_stron)) then
          call shr_sys_abort(trim(subname)//' ERROR: cart3d not initialized '//trim(lstring))
       endif
       call seq_map_cart3d(mapper, type, av_s, av_d, fldu, fldv, norm=lnorm, string=string)
    elseif (trim(type) == 'none') then
       call seq_map_map(mapper, av_s, av_d, fldlist=trim(fldu)//':'//trim(fldv), norm=lnorm)
    else
       write(logunit,*) subname,' ERROR: type unsupported ',trim(type)
       call shr_sys_abort(trim(subname)//' ERROR in type='//trim(type))
    end if

  end subroutine seq_map_mapvect

  !=======================================================================

! Subprogram not used   subroutine seq_map_cart3d( mapper, type, av_s, av_d, fldu, fldv, norm, string)
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used     !-----------------------------------------------------
! Subprogram not used     ! 
! Subprogram not used     ! Arguments
! Subprogram not used     !
! Subprogram not used     type(seq_map)   ,intent(inout)       :: mapper
! Subprogram not used     character(len=*),intent(in)          :: type
! Subprogram not used     type(mct_aVect) ,intent(in)          :: av_s
! Subprogram not used     type(mct_aVect) ,intent(inout)       :: av_d
! Subprogram not used     character(len=*),intent(in)          :: fldu
! Subprogram not used     character(len=*),intent(in)          :: fldv
! Subprogram not used     logical         ,intent(in),optional :: norm
! Subprogram not used     character(len=*),intent(in),optional :: string
! Subprogram not used     !
! Subprogram not used     ! Local Variables
! Subprogram not used     !
! Subprogram not used     integer           :: lsize
! Subprogram not used     logical           :: lnorm
! Subprogram not used     integer           :: ku,kv,kux,kuy,kuz,n
! Subprogram not used     real(r8)          :: ue,un,ur,ux,uy,uz,speed
! Subprogram not used     real(r8)          :: urmaxl,urmax,uravgl,uravg,spavgl,spavg
! Subprogram not used     type(mct_aVect)   :: av3_s, av3_d
! Subprogram not used     integer(in)       :: mpicom,my_task,ierr,urcnt,urcntl
! Subprogram not used     character(len=*),parameter :: subname = "(seq_map_cart3d) "
! Subprogram not used 
! Subprogram not used     lnorm = .true.
! Subprogram not used     if (present(norm)) then
! Subprogram not used        lnorm=norm
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     mpicom = mapper%mpicom
! Subprogram not used 
! Subprogram not used     ku = mct_aVect_indexRA(av_s, trim(fldu), perrwith='quiet')
! Subprogram not used     kv = mct_aVect_indexRA(av_s, trim(fldv), perrwith='quiet')
! Subprogram not used 
! Subprogram not used     if (ku /= 0 .and. kv /= 0) then
! Subprogram not used        lsize = mct_aVect_lsize(av_s)
! Subprogram not used        call mct_avect_init(av3_s,rList='ux:uy:uz',lsize=lsize)
! Subprogram not used 
! Subprogram not used        lsize = mct_aVect_lsize(av_d)
! Subprogram not used        call mct_avect_init(av3_d,rList='ux:uy:uz',lsize=lsize)
! Subprogram not used 
! Subprogram not used        kux = mct_aVect_indexRA(av3_s,'ux')
! Subprogram not used        kuy = mct_aVect_indexRA(av3_s,'uy')
! Subprogram not used        kuz = mct_aVect_indexRA(av3_s,'uz')
! Subprogram not used        lsize = mct_aVect_lsize(av_s)
! Subprogram not used        do n = 1,lsize
! Subprogram not used           ur = 0.0_r8
! Subprogram not used           ue = av_s%rAttr(ku,n)
! Subprogram not used           un = av_s%rAttr(kv,n)
! Subprogram not used           ux = mapper%clon_s(n)*mapper%clat_s(n)*ur - &
! Subprogram not used                mapper%clon_s(n)*mapper%slat_s(n)*un - &
! Subprogram not used                mapper%slon_s(n)*ue
! Subprogram not used           uy = mapper%slon_s(n)*mapper%clon_s(n)*ur - &
! Subprogram not used                mapper%slon_s(n)*mapper%slat_s(n)*un + &
! Subprogram not used                mapper%clon_s(n)*ue
! Subprogram not used           uz = mapper%slat_s(n)*ur + &
! Subprogram not used                mapper%clat_s(n)*un
! Subprogram not used           av3_s%rAttr(kux,n) = ux
! Subprogram not used           av3_s%rAttr(kuy,n) = uy
! Subprogram not used           av3_s%rAttr(kuz,n) = uz
! Subprogram not used        enddo
! Subprogram not used 
! Subprogram not used        call seq_map_map(mapper, av3_s, av3_d, norm=lnorm)
! Subprogram not used 
! Subprogram not used        kux = mct_aVect_indexRA(av3_d,'ux')
! Subprogram not used        kuy = mct_aVect_indexRA(av3_d,'uy')
! Subprogram not used        kuz = mct_aVect_indexRA(av3_d,'uz')
! Subprogram not used        lsize = mct_aVect_lsize(av_d)
! Subprogram not used        urmaxl = -1.0_r8
! Subprogram not used        uravgl = 0.0_r8
! Subprogram not used        urcntl = 0
! Subprogram not used        spavgl = 0.0_r8
! Subprogram not used        do n = 1,lsize
! Subprogram not used           ux = av3_d%rAttr(kux,n)
! Subprogram not used           uy = av3_d%rAttr(kuy,n)
! Subprogram not used           uz = av3_d%rAttr(kuz,n)
! Subprogram not used           ue = -mapper%slon_d(n)          *ux + &
! Subprogram not used                 mapper%clon_d(n)          *uy
! Subprogram not used           un = -mapper%clon_d(n)*mapper%slat_d(n)*ux - &
! Subprogram not used                 mapper%slon_d(n)*mapper%slat_d(n)*uy + &
! Subprogram not used                 mapper%clat_d(n)*uz
! Subprogram not used           ur =  mapper%clon_d(n)*mapper%clat_d(n)*ux + &
! Subprogram not used                 mapper%slon_d(n)*mapper%clat_d(n)*uy - &
! Subprogram not used                 mapper%slat_d(n)*uz
! Subprogram not used           speed = sqrt(ur*ur + ue*ue + un*un)
! Subprogram not used           if (trim(type) == 'cart3d_diag' .or. trim(type) == 'cart3d_uvw_diag') then
! Subprogram not used              if (speed /= 0.0_r8) then
! Subprogram not used                 urmaxl = max(urmaxl,abs(ur))
! Subprogram not used                 uravgl = uravgl + abs(ur)
! Subprogram not used                 spavgl = spavgl + speed
! Subprogram not used                 urcntl = urcntl + 1
! Subprogram not used              endif
! Subprogram not used           endif
! Subprogram not used           if (type(1:10) == 'cart3d_uvw') then
! Subprogram not used              !--- this adds ur to ue and un, while preserving u/v angle and total speed ---
! Subprogram not used              if (un == 0.0_R8) then
! Subprogram not used                 !--- if ue is also 0.0 then just give speed to ue, this is arbitrary ---
! Subprogram not used                 av_d%rAttr(ku,n) = sign(speed,ue)
! Subprogram not used                 av_d%rAttr(kv,n) = 0.0_r8
! Subprogram not used              else if (ue == 0.0_R8) then
! Subprogram not used                 av_d%rAttr(ku,n) = 0.0_r8
! Subprogram not used                 av_d%rAttr(kv,n) = sign(speed,un)
! Subprogram not used              else
! Subprogram not used                 av_d%rAttr(ku,n) = sign(speed/sqrt(1.0_r8 + ((un*un)/(ue*ue))),ue)
! Subprogram not used                 av_d%rAttr(kv,n) = sign(speed/sqrt(1.0_r8 + ((ue*ue)/(un*un))),un)
! Subprogram not used              endif
! Subprogram not used           else
! Subprogram not used              !--- this ignores ur ---
! Subprogram not used              av_d%rAttr(ku,n) = ue
! Subprogram not used              av_d%rAttr(kv,n) = un
! Subprogram not used           endif
! Subprogram not used        enddo
! Subprogram not used        if (trim(type) == 'cart3d_diag' .or. trim(type) == 'cart3d_uvw_diag') then
! Subprogram not used           call mpi_comm_rank(mpicom,my_task,ierr)
! Subprogram not used           call shr_mpi_max(urmaxl,urmax,mpicom,'urmax')
! Subprogram not used           call shr_mpi_sum(uravgl,uravg,mpicom,'uravg')
! Subprogram not used           call shr_mpi_sum(spavgl,spavg,mpicom,'spavg')
! Subprogram not used           call shr_mpi_sum(urcntl,urcnt,mpicom,'urcnt')
! Subprogram not used           if (my_task == 0 .and. urcnt > 0) then
! Subprogram not used              uravg = uravg / urcnt
! Subprogram not used              spavg = spavg / urcnt
! Subprogram not used              write(logunit,*) trim(subname),' cart3d uravg,urmax,spavg = ',uravg,urmax,spavg
! Subprogram not used           endif
! Subprogram not used        endif
! Subprogram not used 
! Subprogram not used       call mct_avect_clean(av3_s)
! Subprogram not used       call mct_avect_clean(av3_d)
! Subprogram not used 
! Subprogram not used    endif  ! ku,kv
! Subprogram not used 
! Subprogram not used   end subroutine seq_map_cart3d

  !=======================================================================

! Subprogram not used   subroutine seq_map_readdata(maprcfile, maprcname, mpicom, ID, &
! Subprogram not used          ni_s, nj_s, av_s, gsmap_s, avfld_s, filefld_s, &
! Subprogram not used          ni_d, nj_d, av_d, gsmap_d, avfld_d, filefld_d, string)
! Subprogram not used 
! Subprogram not used     !--- lifted from work by J Edwards, April 2011
! Subprogram not used 
! Subprogram not used     use shr_pio_mod, only : shr_pio_getiosys, shr_pio_getiotype
! Subprogram not used     use pio, only : pio_openfile, pio_closefile, pio_read_darray, pio_inq_dimid, &
! Subprogram not used        pio_inq_dimlen, pio_inq_varid, file_desc_t, io_desc_t, iosystem_desc_t, &
! Subprogram not used        var_desc_t, pio_int, pio_get_var, pio_double, pio_initdecomp, pio_freedecomp
! Subprogram not used     implicit none
! Subprogram not used     !-----------------------------------------------------
! Subprogram not used     ! 
! Subprogram not used     ! Arguments
! Subprogram not used     !
! Subprogram not used     character(len=*),intent(in)             :: maprcfile
! Subprogram not used     character(len=*),intent(in)             :: maprcname
! Subprogram not used     integer(IN)     ,intent(in)             :: mpicom
! Subprogram not used     integer(IN)     ,intent(in)             :: ID
! Subprogram not used     integer(IN)     ,intent(out)  ,optional :: ni_s
! Subprogram not used     integer(IN)     ,intent(out)  ,optional :: nj_s
! Subprogram not used     type(mct_avect) ,intent(inout),optional :: av_s
! Subprogram not used     type(mct_gsmap) ,intent(in)   ,optional :: gsmap_s
! Subprogram not used     character(len=*),intent(in)   ,optional :: avfld_s
! Subprogram not used     character(len=*),intent(in)   ,optional :: filefld_s
! Subprogram not used     integer(IN)     ,intent(out)  ,optional :: ni_d
! Subprogram not used     integer(IN)     ,intent(out)  ,optional :: nj_d
! Subprogram not used     type(mct_avect) ,intent(inout),optional :: av_d
! Subprogram not used     type(mct_gsmap) ,intent(in)   ,optional :: gsmap_d
! Subprogram not used     character(len=*),intent(in)   ,optional :: avfld_d
! Subprogram not used     character(len=*),intent(in)   ,optional :: filefld_d
! Subprogram not used     character(len=*),intent(in)   ,optional :: string
! Subprogram not used     !
! Subprogram not used     ! Local Variables
! Subprogram not used     !
! Subprogram not used     type(iosystem_desc_t), pointer :: pio_subsystem
! Subprogram not used     integer(IN)       :: pio_iotype
! Subprogram not used     type(file_desc_t) :: File    ! PIO file pointer
! Subprogram not used     type(io_desc_t)   :: iodesc  ! PIO parallel io descriptor
! Subprogram not used     integer(IN)       :: rcode   ! pio routine return code
! Subprogram not used     type(var_desc_t)  :: vid     ! pio variable  ID
! Subprogram not used     integer(IN)       :: did     ! pio dimension ID
! Subprogram not used     integer(IN)       :: na      ! size of source domain
! Subprogram not used     integer(IN)       :: nb      ! size of destination domain
! Subprogram not used     integer(IN)       :: i       ! index
! Subprogram not used     integer(IN)       :: mytask  ! my task
! Subprogram not used     integer(IN), pointer :: dof(:)    ! DOF pointers for parallel read
! Subprogram not used     character(len=256):: fileName
! Subprogram not used     character(len=64) :: lfld_s, lfld_d, lfile_s, lfile_d
! Subprogram not used     character(*),parameter :: areaAV_field = 'aream'
! Subprogram not used     character(*),parameter :: areafile_s   = 'area_a'
! Subprogram not used     character(*),parameter :: areafile_d   = 'area_b'
! Subprogram not used     character(len=*),parameter :: subname  = "(seq_map_readdata) "
! Subprogram not used     !-----------------------------------------------------
! Subprogram not used 
! Subprogram not used     if (seq_comm_iamroot(CPLID) .and. present(string)) then
! Subprogram not used        write(logunit,'(A)') subname//' called for '//trim(string)
! Subprogram not used        call shr_sys_flush(logunit)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     call MPI_COMM_RANK(mpicom,mytask,rcode)
! Subprogram not used 
! Subprogram not used     lfld_s = trim(areaAV_field)
! Subprogram not used     if (present(avfld_s)) then
! Subprogram not used        lfld_s = trim(avfld_s)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     lfld_d = trim(areaAV_field)
! Subprogram not used     if (present(avfld_d)) then
! Subprogram not used        lfld_s = trim(avfld_d)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     lfile_s = trim(areafile_s)
! Subprogram not used     if (present(filefld_s)) then
! Subprogram not used        lfile_s = trim(filefld_s)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     lfile_d = trim(areafile_d)
! Subprogram not used     if (present(filefld_d)) then
! Subprogram not used        lfile_d = trim(filefld_d)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     call I90_allLoadF(trim(maprcfile),0,mpicom,rcode)
! Subprogram not used     if(rcode /= 0) then
! Subprogram not used        write(logunit,*)"Cant find maprcfile file ",trim(maprcfile)
! Subprogram not used        call shr_sys_abort(trim(subname)//"i90_allLoadF File Not Found")
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     call i90_label(trim(maprcname),rcode)
! Subprogram not used     if(rcode /= 0) then
! Subprogram not used        write(logunit,*)"Cant find label ",maprcname
! Subprogram not used        call shr_sys_abort(trim(subname)//"i90_label Not Found")
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     call i90_gtoken(filename,rcode)
! Subprogram not used     if(rcode /= 0) then
! Subprogram not used        write(logunit,*)"Error reading token ",filename
! Subprogram not used        call shr_sys_abort(trim(subname)//"i90_gtoken Error on filename read")
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     pio_subsystem => shr_pio_getiosys(ID)
! Subprogram not used     pio_iotype = shr_pio_getiotype(ID)
! Subprogram not used 
! Subprogram not used     rcode = pio_openfile(pio_subsystem, File, pio_iotype, filename)
! Subprogram not used 
! Subprogram not used     if (present(ni_s)) then 
! Subprogram not used        rcode = pio_inq_dimid (File, 'ni_a', did)  ! number of lons in input grid
! Subprogram not used        rcode = pio_inq_dimlen(File, did  , ni_s)
! Subprogram not used     end if
! Subprogram not used     if(present(nj_s)) then
! Subprogram not used        rcode = pio_inq_dimid (File, 'nj_a', did)  ! number of lats in input grid
! Subprogram not used        rcode = pio_inq_dimlen(File, did  , nj_s)
! Subprogram not used     end if
! Subprogram not used     if(present(ni_d)) then
! Subprogram not used        rcode = pio_inq_dimid (File, 'ni_b', did)  ! number of lons in output grid
! Subprogram not used        rcode = pio_inq_dimlen(File, did  , ni_d)
! Subprogram not used     end if
! Subprogram not used     if(present(nj_d)) then
! Subprogram not used        rcode = pio_inq_dimid (File, 'nj_b', did)  ! number of lats in output grid
! Subprogram not used        rcode = pio_inq_dimlen(File, did  , nj_d)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     !--- read and load area_a ---
! Subprogram not used     if (present(av_s)) then
! Subprogram not used        if (.not.present(gsmap_s)) then
! Subprogram not used           call shr_sys_abort(trim(subname)//' ERROR av_s must have gsmap_s')
! Subprogram not used        endif
! Subprogram not used        rcode = pio_inq_dimid (File, 'n_a', did)  ! size of  input vector
! Subprogram not used        rcode = pio_inq_dimlen(File, did  , na)
! Subprogram not used        i = mct_avect_indexra(av_s, trim(lfld_s))
! Subprogram not used        call mct_gsmap_OrderedPoints(gsMap_s, mytask, dof)
! Subprogram not used        call pio_initdecomp(pio_subsystem, pio_double, (/na/), dof, iodesc)
! Subprogram not used        deallocate(dof)
! Subprogram not used        rcode = pio_inq_varid(File,trim(lfile_s),vid)
! Subprogram not used        call pio_read_darray(File, vid, iodesc, av_s%rattr(i,:), rcode)
! Subprogram not used        call pio_freedecomp(File,iodesc)
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     !--- read and load area_b ---
! Subprogram not used     if (present(av_d)) then
! Subprogram not used        if (.not.present(gsmap_d)) then
! Subprogram not used           call shr_sys_abort(trim(subname)//' ERROR av_d must have gsmap_d')
! Subprogram not used        endif
! Subprogram not used        rcode = pio_inq_dimid (File, 'n_b', did)  ! size of output vector
! Subprogram not used        rcode = pio_inq_dimlen(File, did  , nb)
! Subprogram not used        i = mct_avect_indexra(av_d, trim(lfld_d))
! Subprogram not used        call mct_gsmap_OrderedPoints(gsMap_d, mytask, dof)
! Subprogram not used        call pio_initdecomp(pio_subsystem, pio_double, (/nb/), dof, iodesc)
! Subprogram not used        deallocate(dof)
! Subprogram not used        rcode = pio_inq_varid(File,trim(lfile_d),vid)
! Subprogram not used        call pio_read_darray(File, vid, iodesc, av_d%rattr(i,:), rcode)
! Subprogram not used        call pio_freedecomp(File,iodesc)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     call pio_closefile(File)
! Subprogram not used 
! Subprogram not used   end subroutine seq_map_readdata

!=======================================================================

  subroutine seq_map_avNormAvF(mapper, av_i, av_o, avf_i, avfifld, rList, norm)

    implicit none
    !-----------------------------------------------------
    !
    ! Arguments
    !
    type(seq_map)   , intent(inout)       :: mapper  ! mapper
    type(mct_aVect) , intent(in)          :: av_i    ! input 
    type(mct_aVect) , intent(inout)       :: av_o    ! output
    type(mct_aVect) , intent(in)          :: avf_i   ! extra src "weight"
    character(len=*), intent(in)          :: avfifld ! field name in avf_i
    character(len=*), intent(in),optional :: rList   ! fields list
    logical         , intent(in),optional :: norm    ! normalize at end
    !
    integer(IN) :: lsize_i, lsize_f, lsize_o, kf, j
    real(r8),allocatable :: frac_i(:),frac_o(:)
    logical :: lnorm
    character(*),parameter :: subName = '(seq_map_avNormAvF) '
    !-----------------------------------------------------

    lnorm = .true.
    if (present(norm)) then
       lnorm = norm
    endif

    lsize_i = mct_aVect_lsize(av_i)
    lsize_f = mct_aVect_lsize(avf_i)

    if (lsize_i /= lsize_f) then
       write(logunit,*) subname,' ERROR: lsize_i ne lsize_f ',lsize_i,lsize_f
       call shr_sys_abort(subname//' ERROR size_i ne lsize_f')
    endif

    !--- extract frac_i field from avf_i to pass to seq_map_avNormArr ---
    allocate(frac_i(lsize_i))
    do j = 1,lsize_i
       kf = mct_aVect_indexRA(avf_i,trim(avfifld))
       frac_i(j) = avf_i%rAttr(kf,j)
    enddo

    if (present(rList)) then
       call seq_map_avNormArr(mapper, av_i, av_o, frac_i, rList=rList, norm=lnorm)
    else
       call seq_map_avNormArr(mapper, av_i, av_o, frac_i, norm=lnorm)
    endif

    deallocate(frac_i)

  end subroutine seq_map_avNormAvF

!=======================================================================

  subroutine seq_map_avNormArr(mapper, av_i, av_o, norm_i, rList, norm)

    implicit none
    !-----------------------------------------------------
    !
    ! Arguments
    !
    type(seq_map)   , intent(inout) :: mapper! mapper
    type(mct_aVect) , intent(in)    :: av_i  ! input 
    type(mct_aVect) , intent(inout) :: av_o  ! output
    real(r8)        , intent(in), optional :: norm_i(:)  ! source "weight"
    character(len=*), intent(in), optional :: rList ! fields list
    logical         , intent(in), optional :: norm  ! normalize at end
    !
    ! Local variables
    !
    type(mct_sMatp)        :: sMatp ! sMat
    type(mct_aVect)        :: avp_i , avp_o
    integer(IN)            :: i,j,ier,kf
    integer(IN)            :: lsize_i,lsize_o
    real(r8)               :: normval
    character(CX)          :: lrList
    logical                :: lnorm
    character(*),parameter :: subName = '(seq_map_avNormArr) '
    character(len=*),parameter :: ffld = 'norm8wt'  ! want something unique
    !-----------------------------------------------------

    sMatp   = mapper%sMatp
    lsize_i = mct_aVect_lsize(av_i)
    lsize_o = mct_aVect_lsize(av_o)

    lnorm = .true.
    if (present(norm)) then
       lnorm = norm
    endif

    if (present(norm_i) .and..not.lnorm) then
       write(logunit,*) subname,' ERROR: norm_i and norm = false'
       call shr_sys_abort(subname//' ERROR norm_i and norm = false')
    endif

    if (present(norm_i)) then
       if (size(norm_i) /= lsize_i) then
          write(logunit,*) subname,' ERROR: size(norm_i) ne lsize_i ',size(norm_i),lsize_i
          call shr_sys_abort(subname//' ERROR size(norm_i) ne lsize_i')
       endif
    endif

    !--- create temporary avs for mapping ---

    if (present(rList)) then
       call mct_aVect_init(avp_i, rList=trim( rList)//':'//ffld, lsize=lsize_i)
       call mct_aVect_init(avp_o, rList=trim( rList)//':'//ffld, lsize=lsize_o)
    else
       lrList = trim(mct_aVect_exportRList2c(av_i))
       call mct_aVect_init(avp_i, rList=trim(lrList)//':'//ffld, lsize=lsize_i)
       lrList = trim(mct_aVect_exportRList2c(av_o))
       call mct_aVect_init(avp_o, rList=trim(lrList)//':'//ffld, lsize=lsize_o)
    endif

    !--- copy av_i to avp_i and set ffld value to 1.0
    !--- then multiply all fields by norm_i if norm_i exists 
    !--- this will do the right thing for the norm_i normalization 

    call mct_aVect_copy(aVin=av_i, aVout=avp_i, VECTOR=mct_usevector)
    kf = mct_aVect_indexRA(avp_i,ffld)
    do j = 1,lsize_i
       avp_i%rAttr(kf,j) = 1.0_r8
    enddo

    if (present(norm_i)) then
       do j = 1,lsize_i
          avp_i%rAttr(:,j) = avp_i%rAttr(:,j)*norm_i(j)
       enddo
    endif

    !--- map ---

    if (mapper%esmf_map) then
       call shr_sys_abort(subname//' ERROR: esmf SMM not allowed without USE_ESMF_LIB')

    else
       ! MCT based SMM
       call mct_sMat_avMult(avp_i, sMatp, avp_o, VECTOR=mct_usevector)
    endif


    !--- renormalize avp_o by mapped norm_i  ---

    if (lnorm) then
       do j = 1,lsize_o
          kf = mct_aVect_indexRA(avp_o,ffld)
          normval = avp_o%rAttr(kf,j)
          if (normval /= 0.0_r8) then
             normval = 1.0_r8/normval
          endif
          avp_o%rAttr(:,j) = avp_o%rAttr(:,j)*normval
       enddo
    endif

    !--- copy back into av_o and we are done ---

    call mct_aVect_copy(aVin=avp_o, aVout=av_o, VECTOR=mct_usevector)

    call mct_aVect_clean(avp_i)
    call mct_aVect_clean(avp_o)

  end subroutine seq_map_avNormArr

end module seq_map_mod
