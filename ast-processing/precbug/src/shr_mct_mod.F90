!===============================================================================
! SVN $Id: shr_mct_mod.F90 18548 2009-09-26 23:55:51Z tcraig $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_091114/shr/shr_mct_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: shr_mct_mod -- higher level mct type routines
!     needed to prevent some circular dependencies
!
! !REVISION HISTORY:
!     2009-Dec-16 - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------
module shr_mct_mod

! !USES:

   use shr_kind_mod         ! shared kinds
   use shr_sys_mod          ! share system routines
   use shr_mpi_mod          ! mpi layer
   use shr_const_mod        ! constants
   use mct_mod

   use shr_log_mod          ,only: s_loglev               => shr_log_Level
   use shr_log_mod          ,only: s_logunit              => shr_log_Unit

   implicit none
   private

! PUBLIC: Public interfaces

   public :: shr_mct_sMatReadnc
   interface shr_mct_sMatPInitnc
      module procedure shr_mct_sMatPInitnc_mapfile
   end interface
   public :: shr_mct_sMatPInitnc
   public :: shr_mct_sMatReaddnc
   public :: shr_mct_sMatWritednc
   public :: shr_mct_queryConfigFile

!EOP

   !--- local kinds ---
   integer,parameter,private :: R8 = SHR_KIND_R8
   integer,parameter,private :: IN = SHR_KIND_IN
   integer,parameter,private :: CL = SHR_KIND_CL

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_mct_sMatReadnc - read all mapping data from a NetCDF SCRIP file
!                              in to a full SparseMatrix
!
! !DESCRIPTION:
!   Read in mapping matrix data from a SCRIP netCDF data file so a sMat.
!
! !REMARKS:
!   Based on cpl_map_read
!
! !REVISION HISTORY:
!     2006 Nov 27: R. Jacob
!
! !INTERFACE: ------------------------------------------------------------------

! Subprogram not used subroutine shr_mct_sMatReadnc(sMat,fileName)
! Subprogram not used 
! Subprogram not used !     NetCDF-3.
! Subprogram not used !
! Subprogram not used ! netcdf version 3 fortran interface:
! Subprogram not used !
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! external netcdf data types:
! Subprogram not used !
! Subprogram not used       integer nf_byte
! Subprogram not used       integer nf_int1
! Subprogram not used       integer nf_char
! Subprogram not used       integer nf_short
! Subprogram not used       integer nf_int2
! Subprogram not used       integer nf_int
! Subprogram not used       integer nf_float
! Subprogram not used       integer nf_real
! Subprogram not used       integer nf_double
! Subprogram not used 
! Subprogram not used       parameter (nf_byte = 1)
! Subprogram not used       parameter (nf_int1 = nf_byte)
! Subprogram not used       parameter (nf_char = 2)
! Subprogram not used       parameter (nf_short = 3)
! Subprogram not used       parameter (nf_int2 = nf_short)
! Subprogram not used       parameter (nf_int = 4)
! Subprogram not used       parameter (nf_float = 5)
! Subprogram not used       parameter (nf_real = nf_float)
! Subprogram not used       parameter (nf_double = 6)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! default fill values:
! Subprogram not used !
! Subprogram not used       integer           nf_fill_byte
! Subprogram not used       integer           nf_fill_int1
! Subprogram not used       integer           nf_fill_char
! Subprogram not used       integer           nf_fill_short
! Subprogram not used       integer           nf_fill_int2
! Subprogram not used       integer           nf_fill_int
! Subprogram not used       real              nf_fill_float
! Subprogram not used       real              nf_fill_real
! Subprogram not used       doubleprecision   nf_fill_double
! Subprogram not used 
! Subprogram not used       parameter (nf_fill_byte = -127)
! Subprogram not used       parameter (nf_fill_int1 = nf_fill_byte)
! Subprogram not used       parameter (nf_fill_char = 0)
! Subprogram not used       parameter (nf_fill_short = -32767)
! Subprogram not used       parameter (nf_fill_int2 = nf_fill_short)
! Subprogram not used       parameter (nf_fill_int = -2147483647)
! Subprogram not used       parameter (nf_fill_float = 9.9692099683868690e+36)
! Subprogram not used       parameter (nf_fill_real = nf_fill_float)
! Subprogram not used       parameter (nf_fill_double = 9.9692099683868690d+36)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! mode flags for opening and creating a netcdf dataset:
! Subprogram not used !
! Subprogram not used       integer nf_nowrite
! Subprogram not used       integer nf_write
! Subprogram not used       integer nf_clobber
! Subprogram not used       integer nf_noclobber
! Subprogram not used       integer nf_fill
! Subprogram not used       integer nf_nofill
! Subprogram not used       integer nf_lock
! Subprogram not used       integer nf_share
! Subprogram not used       integer nf_64bit_offset
! Subprogram not used       integer nf_sizehint_default
! Subprogram not used       integer nf_align_chunk
! Subprogram not used       integer nf_format_classic
! Subprogram not used       integer nf_format_64bit
! Subprogram not used       integer nf_diskless
! Subprogram not used       integer nf_mmap
! Subprogram not used 
! Subprogram not used       parameter (nf_nowrite = 0)
! Subprogram not used       parameter (nf_write = 1)
! Subprogram not used       parameter (nf_clobber = 0)
! Subprogram not used       parameter (nf_noclobber = 4)
! Subprogram not used       parameter (nf_fill = 0)
! Subprogram not used       parameter (nf_nofill = 256)
! Subprogram not used       parameter (nf_lock = 1024)
! Subprogram not used       parameter (nf_share = 2048)
! Subprogram not used       parameter (nf_64bit_offset = 512)
! Subprogram not used       parameter (nf_sizehint_default = 0)
! Subprogram not used       parameter (nf_align_chunk = -1)
! Subprogram not used       parameter (nf_format_classic = 1)
! Subprogram not used       parameter (nf_format_64bit = 2)
! Subprogram not used       parameter (nf_diskless = 8)
! Subprogram not used       parameter (nf_mmap = 16)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! size argument for defining an unlimited dimension:
! Subprogram not used !
! Subprogram not used       integer nf_unlimited
! Subprogram not used       parameter (nf_unlimited = 0)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! global attribute id:
! Subprogram not used !
! Subprogram not used       integer nf_global
! Subprogram not used       parameter (nf_global = 0)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! implementation limits:
! Subprogram not used !
! Subprogram not used       integer nf_max_dims
! Subprogram not used       integer nf_max_attrs
! Subprogram not used       integer nf_max_vars
! Subprogram not used       integer nf_max_name
! Subprogram not used       integer nf_max_var_dims
! Subprogram not used 
! Subprogram not used       parameter (nf_max_dims = 1024)
! Subprogram not used       parameter (nf_max_attrs = 8192)
! Subprogram not used       parameter (nf_max_vars = 8192)
! Subprogram not used       parameter (nf_max_name = 256)
! Subprogram not used       parameter (nf_max_var_dims = nf_max_dims)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! error codes:
! Subprogram not used !
! Subprogram not used       integer nf_noerr
! Subprogram not used       integer nf_ebadid
! Subprogram not used       integer nf_eexist
! Subprogram not used       integer nf_einval
! Subprogram not used       integer nf_eperm
! Subprogram not used       integer nf_enotindefine
! Subprogram not used       integer nf_eindefine
! Subprogram not used       integer nf_einvalcoords
! Subprogram not used       integer nf_emaxdims
! Subprogram not used       integer nf_enameinuse
! Subprogram not used       integer nf_enotatt
! Subprogram not used       integer nf_emaxatts
! Subprogram not used       integer nf_ebadtype
! Subprogram not used       integer nf_ebaddim
! Subprogram not used       integer nf_eunlimpos
! Subprogram not used       integer nf_emaxvars
! Subprogram not used       integer nf_enotvar
! Subprogram not used       integer nf_eglobal
! Subprogram not used       integer nf_enotnc
! Subprogram not used       integer nf_ests
! Subprogram not used       integer nf_emaxname
! Subprogram not used       integer nf_eunlimit
! Subprogram not used       integer nf_enorecvars
! Subprogram not used       integer nf_echar
! Subprogram not used       integer nf_eedge
! Subprogram not used       integer nf_estride
! Subprogram not used       integer nf_ebadname
! Subprogram not used       integer nf_erange
! Subprogram not used       integer nf_enomem
! Subprogram not used       integer nf_evarsize
! Subprogram not used       integer nf_edimsize
! Subprogram not used       integer nf_etrunc
! Subprogram not used 
! Subprogram not used       parameter (nf_noerr = 0)
! Subprogram not used       parameter (nf_ebadid = -33)
! Subprogram not used       parameter (nf_eexist = -35)
! Subprogram not used       parameter (nf_einval = -36)
! Subprogram not used       parameter (nf_eperm = -37)
! Subprogram not used       parameter (nf_enotindefine = -38)
! Subprogram not used       parameter (nf_eindefine = -39)
! Subprogram not used       parameter (nf_einvalcoords = -40)
! Subprogram not used       parameter (nf_emaxdims = -41)
! Subprogram not used       parameter (nf_enameinuse = -42)
! Subprogram not used       parameter (nf_enotatt = -43)
! Subprogram not used       parameter (nf_emaxatts = -44)
! Subprogram not used       parameter (nf_ebadtype = -45)
! Subprogram not used       parameter (nf_ebaddim = -46)
! Subprogram not used       parameter (nf_eunlimpos = -47)
! Subprogram not used       parameter (nf_emaxvars = -48)
! Subprogram not used       parameter (nf_enotvar = -49)
! Subprogram not used       parameter (nf_eglobal = -50)
! Subprogram not used       parameter (nf_enotnc = -51)
! Subprogram not used       parameter (nf_ests = -52)
! Subprogram not used       parameter (nf_emaxname = -53)
! Subprogram not used       parameter (nf_eunlimit = -54)
! Subprogram not used       parameter (nf_enorecvars = -55)
! Subprogram not used       parameter (nf_echar = -56)
! Subprogram not used       parameter (nf_eedge = -57)
! Subprogram not used       parameter (nf_estride = -58)
! Subprogram not used       parameter (nf_ebadname = -59)
! Subprogram not used       parameter (nf_erange = -60)
! Subprogram not used       parameter (nf_enomem = -61)
! Subprogram not used       parameter (nf_evarsize = -62)
! Subprogram not used       parameter (nf_edimsize = -63)
! Subprogram not used       parameter (nf_etrunc = -64)
! Subprogram not used !
! Subprogram not used ! error handling modes:
! Subprogram not used !
! Subprogram not used       integer  nf_fatal
! Subprogram not used       integer nf_verbose
! Subprogram not used 
! Subprogram not used       parameter (nf_fatal = 1)
! Subprogram not used       parameter (nf_verbose = 2)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! miscellaneous routines:
! Subprogram not used !
! Subprogram not used       character*80   nf_inq_libvers
! Subprogram not used       external       nf_inq_libvers
! Subprogram not used 
! Subprogram not used       character*80   nf_strerror
! Subprogram not used !                         (integer             ncerr)
! Subprogram not used       external       nf_strerror
! Subprogram not used 
! Subprogram not used       logical        nf_issyserr
! Subprogram not used !                         (integer             ncerr)
! Subprogram not used       external       nf_issyserr
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! control routines:
! Subprogram not used !
! Subprogram not used       integer         nf_inq_base_pe
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             pe)
! Subprogram not used       external        nf_inq_base_pe
! Subprogram not used 
! Subprogram not used       integer         nf_set_base_pe
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             pe)
! Subprogram not used       external        nf_set_base_pe
! Subprogram not used 
! Subprogram not used       integer         nf_create
! Subprogram not used !                         (character*(*)       path,
! Subprogram not used !                          integer             cmode,
! Subprogram not used !                          integer             ncid)
! Subprogram not used       external        nf_create
! Subprogram not used 
! Subprogram not used       integer         nf__create
! Subprogram not used !                         (character*(*)       path,
! Subprogram not used !                          integer             cmode,
! Subprogram not used !                          integer             initialsz,
! Subprogram not used !                          integer             chunksizehint,
! Subprogram not used !                          integer             ncid)
! Subprogram not used       external        nf__create
! Subprogram not used 
! Subprogram not used       integer         nf__create_mp
! Subprogram not used !                         (character*(*)       path,
! Subprogram not used !                          integer             cmode,
! Subprogram not used !                          integer             initialsz,
! Subprogram not used !                          integer             basepe,
! Subprogram not used !                          integer             chunksizehint,
! Subprogram not used !                          integer             ncid)
! Subprogram not used       external        nf__create_mp
! Subprogram not used 
! Subprogram not used       integer         nf_open
! Subprogram not used !                         (character*(*)       path,
! Subprogram not used !                          integer             mode,
! Subprogram not used !                          integer             ncid)
! Subprogram not used       external        nf_open
! Subprogram not used 
! Subprogram not used       integer         nf__open
! Subprogram not used !                         (character*(*)       path,
! Subprogram not used !                          integer             mode,
! Subprogram not used !                          integer             chunksizehint,
! Subprogram not used !                          integer             ncid)
! Subprogram not used       external        nf__open
! Subprogram not used 
! Subprogram not used       integer         nf__open_mp
! Subprogram not used !                         (character*(*)       path,
! Subprogram not used !                          integer             mode,
! Subprogram not used !                          integer             basepe,
! Subprogram not used !                          integer             chunksizehint,
! Subprogram not used !                          integer             ncid)
! Subprogram not used       external        nf__open_mp
! Subprogram not used 
! Subprogram not used       integer         nf_set_fill
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             fillmode,
! Subprogram not used !                          integer             old_mode)
! Subprogram not used       external        nf_set_fill
! Subprogram not used 
! Subprogram not used       integer         nf_set_default_format
! Subprogram not used !                          (integer             format,
! Subprogram not used !                          integer             old_format)
! Subprogram not used       external        nf_set_default_format
! Subprogram not used 
! Subprogram not used       integer         nf_redef
! Subprogram not used !                         (integer             ncid)
! Subprogram not used       external        nf_redef
! Subprogram not used 
! Subprogram not used       integer         nf_enddef
! Subprogram not used !                         (integer             ncid)
! Subprogram not used       external        nf_enddef
! Subprogram not used 
! Subprogram not used       integer         nf__enddef
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             h_minfree,
! Subprogram not used !                          integer             v_align,
! Subprogram not used !                          integer             v_minfree,
! Subprogram not used !                          integer             r_align)
! Subprogram not used       external        nf__enddef
! Subprogram not used 
! Subprogram not used       integer         nf_sync
! Subprogram not used !                         (integer             ncid)
! Subprogram not used       external        nf_sync
! Subprogram not used 
! Subprogram not used       integer         nf_abort
! Subprogram not used !                         (integer             ncid)
! Subprogram not used       external        nf_abort
! Subprogram not used 
! Subprogram not used       integer         nf_close
! Subprogram not used !                         (integer             ncid)
! Subprogram not used       external        nf_close
! Subprogram not used 
! Subprogram not used       integer         nf_delete
! Subprogram not used !                         (character*(*)       ncid)
! Subprogram not used       external        nf_delete
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! general inquiry routines:
! Subprogram not used !
! Subprogram not used 
! Subprogram not used       integer         nf_inq
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             ndims,
! Subprogram not used !                          integer             nvars,
! Subprogram not used !                          integer             ngatts,
! Subprogram not used !                          integer             unlimdimid)
! Subprogram not used       external        nf_inq
! Subprogram not used 
! Subprogram not used ! new inquire path
! Subprogram not used 
! Subprogram not used       integer nf_inq_path
! Subprogram not used       external nf_inq_path
! Subprogram not used 
! Subprogram not used       integer         nf_inq_ndims
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             ndims)
! Subprogram not used       external        nf_inq_ndims
! Subprogram not used 
! Subprogram not used       integer         nf_inq_nvars
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             nvars)
! Subprogram not used       external        nf_inq_nvars
! Subprogram not used 
! Subprogram not used       integer         nf_inq_natts
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             ngatts)
! Subprogram not used       external        nf_inq_natts
! Subprogram not used 
! Subprogram not used       integer         nf_inq_unlimdim
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             unlimdimid)
! Subprogram not used       external        nf_inq_unlimdim
! Subprogram not used 
! Subprogram not used       integer         nf_inq_format
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             format)
! Subprogram not used       external        nf_inq_format
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! dimension routines:
! Subprogram not used !
! Subprogram not used 
! Subprogram not used       integer         nf_def_dim
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             len,
! Subprogram not used !                          integer             dimid)
! Subprogram not used       external        nf_def_dim
! Subprogram not used 
! Subprogram not used       integer         nf_inq_dimid
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             dimid)
! Subprogram not used       external        nf_inq_dimid
! Subprogram not used 
! Subprogram not used       integer         nf_inq_dim
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             dimid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             len)
! Subprogram not used       external        nf_inq_dim
! Subprogram not used 
! Subprogram not used       integer         nf_inq_dimname
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             dimid,
! Subprogram not used !                          character(*)        name)
! Subprogram not used       external        nf_inq_dimname
! Subprogram not used 
! Subprogram not used       integer         nf_inq_dimlen
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             dimid,
! Subprogram not used !                          integer             len)
! Subprogram not used       external        nf_inq_dimlen
! Subprogram not used 
! Subprogram not used       integer         nf_rename_dim
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             dimid,
! Subprogram not used !                          character(*)        name)
! Subprogram not used       external        nf_rename_dim
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! general attribute routines:
! Subprogram not used !
! Subprogram not used 
! Subprogram not used       integer         nf_inq_att
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             xtype,
! Subprogram not used !                          integer             len)
! Subprogram not used       external        nf_inq_att
! Subprogram not used 
! Subprogram not used       integer         nf_inq_attid
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             attnum)
! Subprogram not used       external        nf_inq_attid
! Subprogram not used 
! Subprogram not used       integer         nf_inq_atttype
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             xtype)
! Subprogram not used       external        nf_inq_atttype
! Subprogram not used 
! Subprogram not used       integer         nf_inq_attlen
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             len)
! Subprogram not used       external        nf_inq_attlen
! Subprogram not used 
! Subprogram not used       integer         nf_inq_attname
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             attnum,
! Subprogram not used !                          character(*)        name)
! Subprogram not used       external        nf_inq_attname
! Subprogram not used 
! Subprogram not used       integer         nf_copy_att
! Subprogram not used !                         (integer             ncid_in,
! Subprogram not used !                          integer             varid_in,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             ncid_out,
! Subprogram not used !                          integer             varid_out)
! Subprogram not used       external        nf_copy_att
! Subprogram not used 
! Subprogram not used       integer         nf_rename_att
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        curname,
! Subprogram not used !                          character(*)        newname)
! Subprogram not used       external        nf_rename_att
! Subprogram not used 
! Subprogram not used       integer         nf_del_att
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name)
! Subprogram not used       external        nf_del_att
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! attribute put/get routines:
! Subprogram not used !
! Subprogram not used 
! Subprogram not used       integer         nf_put_att_text
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             len,
! Subprogram not used !                          character(*)        text)
! Subprogram not used       external        nf_put_att_text
! Subprogram not used 
! Subprogram not used       integer         nf_get_att_text
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          character(*)        text)
! Subprogram not used       external        nf_get_att_text
! Subprogram not used 
! Subprogram not used       integer         nf_put_att_int1
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             xtype,
! Subprogram not used !                          integer             len,
! Subprogram not used !                          nf_int1_t           i1vals(1))
! Subprogram not used       external        nf_put_att_int1
! Subprogram not used 
! Subprogram not used       integer         nf_get_att_int1
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          nf_int1_t           i1vals(1))
! Subprogram not used       external        nf_get_att_int1
! Subprogram not used 
! Subprogram not used       integer         nf_put_att_int2
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             xtype,
! Subprogram not used !                          integer             len,
! Subprogram not used !                          nf_int2_t           i2vals(1))
! Subprogram not used       external        nf_put_att_int2
! Subprogram not used 
! Subprogram not used       integer         nf_get_att_int2
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          nf_int2_t           i2vals(1))
! Subprogram not used       external        nf_get_att_int2
! Subprogram not used 
! Subprogram not used       integer         nf_put_att_int
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             xtype,
! Subprogram not used !                          integer             len,
! Subprogram not used !                          integer             ivals(1))
! Subprogram not used       external        nf_put_att_int
! Subprogram not used 
! Subprogram not used       integer         nf_get_att_int
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             ivals(1))
! Subprogram not used       external        nf_get_att_int
! Subprogram not used 
! Subprogram not used       integer         nf_put_att_real
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             xtype,
! Subprogram not used !                          integer             len,
! Subprogram not used !                          real                rvals(1))
! Subprogram not used       external        nf_put_att_real
! Subprogram not used 
! Subprogram not used       integer         nf_get_att_real
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          real                rvals(1))
! Subprogram not used       external        nf_get_att_real
! Subprogram not used 
! Subprogram not used       integer         nf_put_att_double
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             xtype,
! Subprogram not used !                          integer             len,
! Subprogram not used !                          double              dvals(1))
! Subprogram not used       external        nf_put_att_double
! Subprogram not used 
! Subprogram not used       integer         nf_get_att_double
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          double              dvals(1))
! Subprogram not used       external        nf_get_att_double
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! general variable routines:
! Subprogram not used !
! Subprogram not used 
! Subprogram not used       integer         nf_def_var
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             datatype,
! Subprogram not used !                          integer             ndims,
! Subprogram not used !                          integer             dimids(1),
! Subprogram not used !                          integer             varid)
! Subprogram not used       external        nf_def_var
! Subprogram not used 
! Subprogram not used       integer         nf_inq_var
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             datatype,
! Subprogram not used !                          integer             ndims,
! Subprogram not used !                          integer             dimids(1),
! Subprogram not used !                          integer             natts)
! Subprogram not used       external        nf_inq_var
! Subprogram not used 
! Subprogram not used       integer         nf_inq_varid
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             varid)
! Subprogram not used       external        nf_inq_varid
! Subprogram not used 
! Subprogram not used       integer         nf_inq_varname
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name)
! Subprogram not used       external        nf_inq_varname
! Subprogram not used 
! Subprogram not used       integer         nf_inq_vartype
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             xtype)
! Subprogram not used       external        nf_inq_vartype
! Subprogram not used 
! Subprogram not used       integer         nf_inq_varndims
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             ndims)
! Subprogram not used       external        nf_inq_varndims
! Subprogram not used 
! Subprogram not used       integer         nf_inq_vardimid
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             dimids(1))
! Subprogram not used       external        nf_inq_vardimid
! Subprogram not used 
! Subprogram not used       integer         nf_inq_varnatts
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             natts)
! Subprogram not used       external        nf_inq_varnatts
! Subprogram not used 
! Subprogram not used       integer         nf_rename_var
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name)
! Subprogram not used       external        nf_rename_var
! Subprogram not used 
! Subprogram not used       integer         nf_copy_var
! Subprogram not used !                         (integer             ncid_in,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             ncid_out)
! Subprogram not used       external        nf_copy_var
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! entire variable put/get routines:
! Subprogram not used !
! Subprogram not used 
! Subprogram not used       integer         nf_put_var_text
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        text)
! Subprogram not used       external        nf_put_var_text
! Subprogram not used 
! Subprogram not used       integer         nf_get_var_text
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        text)
! Subprogram not used       external        nf_get_var_text
! Subprogram not used 
! Subprogram not used       integer         nf_put_var_int1
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          nf_int1_t           i1vals(1))
! Subprogram not used       external        nf_put_var_int1
! Subprogram not used 
! Subprogram not used       integer         nf_get_var_int1
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          nf_int1_t           i1vals(1))
! Subprogram not used       external        nf_get_var_int1
! Subprogram not used 
! Subprogram not used       integer         nf_put_var_int2
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          nf_int2_t           i2vals(1))
! Subprogram not used       external        nf_put_var_int2
! Subprogram not used 
! Subprogram not used       integer         nf_get_var_int2
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          nf_int2_t           i2vals(1))
! Subprogram not used       external        nf_get_var_int2
! Subprogram not used 
! Subprogram not used       integer         nf_put_var_int
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             ivals(1))
! Subprogram not used       external        nf_put_var_int
! Subprogram not used 
! Subprogram not used       integer         nf_get_var_int
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             ivals(1))
! Subprogram not used       external        nf_get_var_int
! Subprogram not used 
! Subprogram not used       integer         nf_put_var_real
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          real                rvals(1))
! Subprogram not used       external        nf_put_var_real
! Subprogram not used 
! Subprogram not used       integer         nf_get_var_real
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          real                rvals(1))
! Subprogram not used       external        nf_get_var_real
! Subprogram not used 
! Subprogram not used       integer         nf_put_var_double
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          doubleprecision     dvals(1))
! Subprogram not used       external        nf_put_var_double
! Subprogram not used 
! Subprogram not used       integer         nf_get_var_double
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          doubleprecision     dvals(1))
! Subprogram not used       external        nf_get_var_double
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! single variable put/get routines:
! Subprogram not used !
! Subprogram not used 
! Subprogram not used       integer         nf_put_var1_text
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             index(1),
! Subprogram not used !                          character*1         text)
! Subprogram not used       external        nf_put_var1_text
! Subprogram not used 
! Subprogram not used       integer         nf_get_var1_text
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             index(1),
! Subprogram not used !                          character*1         text)
! Subprogram not used       external        nf_get_var1_text
! Subprogram not used 
! Subprogram not used       integer         nf_put_var1_int1
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             index(1),
! Subprogram not used !                          nf_int1_t           i1val)
! Subprogram not used       external        nf_put_var1_int1
! Subprogram not used 
! Subprogram not used       integer         nf_get_var1_int1
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             index(1),
! Subprogram not used !                          nf_int1_t           i1val)
! Subprogram not used       external        nf_get_var1_int1
! Subprogram not used 
! Subprogram not used       integer         nf_put_var1_int2
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             index(1),
! Subprogram not used !                          nf_int2_t           i2val)
! Subprogram not used       external        nf_put_var1_int2
! Subprogram not used 
! Subprogram not used       integer         nf_get_var1_int2
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             index(1),
! Subprogram not used !                          nf_int2_t           i2val)
! Subprogram not used       external        nf_get_var1_int2
! Subprogram not used 
! Subprogram not used       integer         nf_put_var1_int
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             index(1),
! Subprogram not used !                          integer             ival)
! Subprogram not used       external        nf_put_var1_int
! Subprogram not used 
! Subprogram not used       integer         nf_get_var1_int
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             index(1),
! Subprogram not used !                          integer             ival)
! Subprogram not used       external        nf_get_var1_int
! Subprogram not used 
! Subprogram not used       integer         nf_put_var1_real
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             index(1),
! Subprogram not used !                          real                rval)
! Subprogram not used       external        nf_put_var1_real
! Subprogram not used 
! Subprogram not used       integer         nf_get_var1_real
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             index(1),
! Subprogram not used !                          real                rval)
! Subprogram not used       external        nf_get_var1_real
! Subprogram not used 
! Subprogram not used       integer         nf_put_var1_double
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             index(1),
! Subprogram not used !                          doubleprecision     dval)
! Subprogram not used       external        nf_put_var1_double
! Subprogram not used 
! Subprogram not used       integer         nf_get_var1_double
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             index(1),
! Subprogram not used !                          doubleprecision     dval)
! Subprogram not used       external        nf_get_var1_double
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! variable array put/get routines:
! Subprogram not used !
! Subprogram not used 
! Subprogram not used       integer         nf_put_vara_text
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          character(*)        text)
! Subprogram not used       external        nf_put_vara_text
! Subprogram not used 
! Subprogram not used       integer         nf_get_vara_text
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          character(*)        text)
! Subprogram not used       external        nf_get_vara_text
! Subprogram not used 
! Subprogram not used       integer         nf_put_vara_int1
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          nf_int1_t           i1vals(1))
! Subprogram not used       external        nf_put_vara_int1
! Subprogram not used 
! Subprogram not used       integer         nf_get_vara_int1
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          nf_int1_t           i1vals(1))
! Subprogram not used       external        nf_get_vara_int1
! Subprogram not used 
! Subprogram not used       integer         nf_put_vara_int2
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          nf_int2_t           i2vals(1))
! Subprogram not used       external        nf_put_vara_int2
! Subprogram not used 
! Subprogram not used       integer         nf_get_vara_int2
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          nf_int2_t           i2vals(1))
! Subprogram not used       external        nf_get_vara_int2
! Subprogram not used 
! Subprogram not used       integer         nf_put_vara_int
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             ivals(1))
! Subprogram not used       external        nf_put_vara_int
! Subprogram not used 
! Subprogram not used       integer         nf_get_vara_int
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             ivals(1))
! Subprogram not used       external        nf_get_vara_int
! Subprogram not used 
! Subprogram not used       integer         nf_put_vara_real
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          real                rvals(1))
! Subprogram not used       external        nf_put_vara_real
! Subprogram not used 
! Subprogram not used       integer         nf_get_vara_real
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          real                rvals(1))
! Subprogram not used       external        nf_get_vara_real
! Subprogram not used 
! Subprogram not used       integer         nf_put_vara_double
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          doubleprecision     dvals(1))
! Subprogram not used       external        nf_put_vara_double
! Subprogram not used 
! Subprogram not used       integer         nf_get_vara_double
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          doubleprecision     dvals(1))
! Subprogram not used       external        nf_get_vara_double
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! strided variable put/get routines:
! Subprogram not used !
! Subprogram not used 
! Subprogram not used       integer         nf_put_vars_text
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          character(*)        text)
! Subprogram not used       external        nf_put_vars_text
! Subprogram not used 
! Subprogram not used       integer         nf_get_vars_text
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          character(*)        text)
! Subprogram not used       external        nf_get_vars_text
! Subprogram not used 
! Subprogram not used       integer         nf_put_vars_int1
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          nf_int1_t           i1vals(1))
! Subprogram not used       external        nf_put_vars_int1
! Subprogram not used 
! Subprogram not used       integer         nf_get_vars_int1
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          nf_int1_t           i1vals(1))
! Subprogram not used       external        nf_get_vars_int1
! Subprogram not used 
! Subprogram not used       integer         nf_put_vars_int2
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          nf_int2_t           i2vals(1))
! Subprogram not used       external        nf_put_vars_int2
! Subprogram not used 
! Subprogram not used       integer         nf_get_vars_int2
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          nf_int2_t           i2vals(1))
! Subprogram not used       external        nf_get_vars_int2
! Subprogram not used 
! Subprogram not used       integer         nf_put_vars_int
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          integer             ivals(1))
! Subprogram not used       external        nf_put_vars_int
! Subprogram not used 
! Subprogram not used       integer         nf_get_vars_int
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          integer             ivals(1))
! Subprogram not used       external        nf_get_vars_int
! Subprogram not used 
! Subprogram not used       integer         nf_put_vars_real
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          real                rvals(1))
! Subprogram not used       external        nf_put_vars_real
! Subprogram not used 
! Subprogram not used       integer         nf_get_vars_real
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          real                rvals(1))
! Subprogram not used       external        nf_get_vars_real
! Subprogram not used 
! Subprogram not used       integer         nf_put_vars_double
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          doubleprecision     dvals(1))
! Subprogram not used       external        nf_put_vars_double
! Subprogram not used 
! Subprogram not used       integer         nf_get_vars_double
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          doubleprecision     dvals(1))
! Subprogram not used       external        nf_get_vars_double
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! mapped variable put/get routines:
! Subprogram not used !
! Subprogram not used 
! Subprogram not used       integer         nf_put_varm_text
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          integer             imap(1),
! Subprogram not used !                          character(*)        text)
! Subprogram not used       external        nf_put_varm_text
! Subprogram not used 
! Subprogram not used       integer         nf_get_varm_text
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          integer             imap(1),
! Subprogram not used !                          character(*)        text)
! Subprogram not used       external        nf_get_varm_text
! Subprogram not used 
! Subprogram not used       integer         nf_put_varm_int1
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          integer             imap(1),
! Subprogram not used !                          nf_int1_t           i1vals(1))
! Subprogram not used       external        nf_put_varm_int1
! Subprogram not used 
! Subprogram not used       integer         nf_get_varm_int1
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          integer             imap(1),
! Subprogram not used !                          nf_int1_t           i1vals(1))
! Subprogram not used       external        nf_get_varm_int1
! Subprogram not used 
! Subprogram not used       integer         nf_put_varm_int2
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          integer             imap(1),
! Subprogram not used !                          nf_int2_t           i2vals(1))
! Subprogram not used       external        nf_put_varm_int2
! Subprogram not used 
! Subprogram not used       integer         nf_get_varm_int2
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          integer             imap(1),
! Subprogram not used !                          nf_int2_t           i2vals(1))
! Subprogram not used       external        nf_get_varm_int2
! Subprogram not used 
! Subprogram not used       integer         nf_put_varm_int
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          integer             imap(1),
! Subprogram not used !                          integer             ivals(1))
! Subprogram not used       external        nf_put_varm_int
! Subprogram not used 
! Subprogram not used       integer         nf_get_varm_int
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          integer             imap(1),
! Subprogram not used !                          integer             ivals(1))
! Subprogram not used       external        nf_get_varm_int
! Subprogram not used 
! Subprogram not used       integer         nf_put_varm_real
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          integer             imap(1),
! Subprogram not used !                          real                rvals(1))
! Subprogram not used       external        nf_put_varm_real
! Subprogram not used 
! Subprogram not used       integer         nf_get_varm_real
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          integer             imap(1),
! Subprogram not used !                          real                rvals(1))
! Subprogram not used       external        nf_get_varm_real
! Subprogram not used 
! Subprogram not used       integer         nf_put_varm_double
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          integer             imap(1),
! Subprogram not used !                          doubleprecision     dvals(1))
! Subprogram not used       external        nf_put_varm_double
! Subprogram not used 
! Subprogram not used       integer         nf_get_varm_double
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          integer             imap(1),
! Subprogram not used !                          doubleprecision     dvals(1))
! Subprogram not used       external        nf_get_varm_double
! Subprogram not used 
! Subprogram not used 
! Subprogram not used !     NetCDF-4.
! Subprogram not used !     This is part of netCDF-4. Copyright 2006, UCAR, See COPYRIGHT
! Subprogram not used !     file for distribution information.
! Subprogram not used 
! Subprogram not used !     Netcdf version 4 fortran interface.
! Subprogram not used 
! Subprogram not used !     $Id: netcdf4.inc,v 1.28 2010/05/25 13:53:02 ed Exp $
! Subprogram not used 
! Subprogram not used !     New netCDF-4 types.
! Subprogram not used       integer nf_ubyte
! Subprogram not used       integer nf_ushort
! Subprogram not used       integer nf_uint
! Subprogram not used       integer nf_int64
! Subprogram not used       integer nf_uint64
! Subprogram not used       integer nf_string
! Subprogram not used       integer nf_vlen
! Subprogram not used       integer nf_opaque
! Subprogram not used       integer nf_enum
! Subprogram not used       integer nf_compound
! Subprogram not used 
! Subprogram not used       parameter (nf_ubyte = 7)
! Subprogram not used       parameter (nf_ushort = 8)
! Subprogram not used       parameter (nf_uint = 9)
! Subprogram not used       parameter (nf_int64 = 10)
! Subprogram not used       parameter (nf_uint64 = 11)
! Subprogram not used       parameter (nf_string = 12)
! Subprogram not used       parameter (nf_vlen = 13)
! Subprogram not used       parameter (nf_opaque = 14)
! Subprogram not used       parameter (nf_enum = 15)
! Subprogram not used       parameter (nf_compound = 16)
! Subprogram not used 
! Subprogram not used !     New netCDF-4 fill values.
! Subprogram not used       integer           nf_fill_ubyte
! Subprogram not used       integer           nf_fill_ushort
! Subprogram not used !      real              nf_fill_uint
! Subprogram not used !      real              nf_fill_int64
! Subprogram not used !      real              nf_fill_uint64
! Subprogram not used       parameter (nf_fill_ubyte = 255)
! Subprogram not used       parameter (nf_fill_ushort = 65535)
! Subprogram not used 
! Subprogram not used !     New constants.
! Subprogram not used       integer nf_format_netcdf4
! Subprogram not used       parameter (nf_format_netcdf4 = 3)
! Subprogram not used 
! Subprogram not used       integer nf_format_netcdf4_classic
! Subprogram not used       parameter (nf_format_netcdf4_classic = 4)
! Subprogram not used 
! Subprogram not used       integer nf_netcdf4
! Subprogram not used       parameter (nf_netcdf4 = 4096)
! Subprogram not used 
! Subprogram not used       integer nf_classic_model
! Subprogram not used       parameter (nf_classic_model = 256)
! Subprogram not used 
! Subprogram not used       integer nf_chunk_seq
! Subprogram not used       parameter (nf_chunk_seq = 0)
! Subprogram not used       integer nf_chunk_sub
! Subprogram not used       parameter (nf_chunk_sub = 1)
! Subprogram not used       integer nf_chunk_sizes
! Subprogram not used       parameter (nf_chunk_sizes = 2)
! Subprogram not used 
! Subprogram not used       integer nf_endian_native
! Subprogram not used       parameter (nf_endian_native = 0)
! Subprogram not used       integer nf_endian_little
! Subprogram not used       parameter (nf_endian_little = 1)
! Subprogram not used       integer nf_endian_big
! Subprogram not used       parameter (nf_endian_big = 2)
! Subprogram not used 
! Subprogram not used !     For NF_DEF_VAR_CHUNKING
! Subprogram not used       integer nf_chunked
! Subprogram not used       parameter (nf_chunked = 0)
! Subprogram not used       integer nf_contiguous
! Subprogram not used       parameter (nf_contiguous = 1)
! Subprogram not used 
! Subprogram not used !     For NF_DEF_VAR_FLETCHER32
! Subprogram not used       integer nf_nochecksum
! Subprogram not used       parameter (nf_nochecksum = 0)
! Subprogram not used       integer nf_fletcher32
! Subprogram not used       parameter (nf_fletcher32 = 1)
! Subprogram not used 
! Subprogram not used !     For NF_DEF_VAR_DEFLATE
! Subprogram not used       integer nf_noshuffle
! Subprogram not used       parameter (nf_noshuffle = 0)
! Subprogram not used       integer nf_shuffle
! Subprogram not used       parameter (nf_shuffle = 1)
! Subprogram not used 
! Subprogram not used !     For NF_DEF_VAR_SZIP
! Subprogram not used       integer nf_szip_ec_option_mask
! Subprogram not used       parameter (nf_szip_ec_option_mask = 4)
! Subprogram not used       integer nf_szip_nn_option_mask
! Subprogram not used       parameter (nf_szip_nn_option_mask = 32)
! Subprogram not used 
! Subprogram not used !     For parallel I/O.
! Subprogram not used       integer nf_mpiio      
! Subprogram not used       parameter (nf_mpiio = 8192)
! Subprogram not used       integer nf_mpiposix
! Subprogram not used       parameter (nf_mpiposix = 16384)
! Subprogram not used       integer nf_pnetcdf
! Subprogram not used       parameter (nf_pnetcdf = 32768)
! Subprogram not used 
! Subprogram not used !     For NF_VAR_PAR_ACCESS.
! Subprogram not used       integer nf_independent
! Subprogram not used       parameter (nf_independent = 0)
! Subprogram not used       integer nf_collective
! Subprogram not used       parameter (nf_collective = 1)
! Subprogram not used 
! Subprogram not used !     New error codes.
! Subprogram not used       integer nf_ehdferr        ! Error at HDF5 layer. 
! Subprogram not used       parameter (nf_ehdferr = -101)
! Subprogram not used       integer nf_ecantread      ! Can't read. 
! Subprogram not used       parameter (nf_ecantread = -102)
! Subprogram not used       integer nf_ecantwrite     ! Can't write. 
! Subprogram not used       parameter (nf_ecantwrite = -103)
! Subprogram not used       integer nf_ecantcreate    ! Can't create. 
! Subprogram not used       parameter (nf_ecantcreate = -104)
! Subprogram not used       integer nf_efilemeta      ! Problem with file metadata. 
! Subprogram not used       parameter (nf_efilemeta = -105)
! Subprogram not used       integer nf_edimmeta       ! Problem with dimension metadata. 
! Subprogram not used       parameter (nf_edimmeta = -106)
! Subprogram not used       integer nf_eattmeta       ! Problem with attribute metadata. 
! Subprogram not used       parameter (nf_eattmeta = -107)
! Subprogram not used       integer nf_evarmeta       ! Problem with variable metadata. 
! Subprogram not used       parameter (nf_evarmeta = -108)
! Subprogram not used       integer nf_enocompound    ! Not a compound type. 
! Subprogram not used       parameter (nf_enocompound = -109)
! Subprogram not used       integer nf_eattexists     ! Attribute already exists. 
! Subprogram not used       parameter (nf_eattexists = -110)
! Subprogram not used       integer nf_enotnc4        ! Attempting netcdf-4 operation on netcdf-3 file.   
! Subprogram not used       parameter (nf_enotnc4 = -111)
! Subprogram not used       integer nf_estrictnc3     ! Attempting netcdf-4 operation on strict nc3 netcdf-4 file.   
! Subprogram not used       parameter (nf_estrictnc3 = -112)
! Subprogram not used       integer nf_enotnc3        ! Attempting netcdf-3 operation on netcdf-4 file.   
! Subprogram not used       parameter (nf_enotnc3 = -113)
! Subprogram not used       integer nf_enopar         ! Parallel operation on file opened for non-parallel access.   
! Subprogram not used       parameter (nf_enopar = -114)
! Subprogram not used       integer nf_eparinit       ! Error initializing for parallel access.   
! Subprogram not used       parameter (nf_eparinit = -115)
! Subprogram not used       integer nf_ebadgrpid      ! Bad group ID.   
! Subprogram not used       parameter (nf_ebadgrpid = -116)
! Subprogram not used       integer nf_ebadtypid      ! Bad type ID.   
! Subprogram not used       parameter (nf_ebadtypid = -117)
! Subprogram not used       integer nf_etypdefined    ! Type has already been defined and may not be edited. 
! Subprogram not used       parameter (nf_etypdefined = -118)
! Subprogram not used       integer nf_ebadfield      ! Bad field ID.   
! Subprogram not used       parameter (nf_ebadfield = -119)
! Subprogram not used       integer nf_ebadclass      ! Bad class.   
! Subprogram not used       parameter (nf_ebadclass = -120)
! Subprogram not used       integer nf_emaptype       ! Mapped access for atomic types only.   
! Subprogram not used       parameter (nf_emaptype = -121)
! Subprogram not used       integer nf_elatefill      ! Attempt to define fill value when data already exists. 
! Subprogram not used       parameter (nf_elatefill = -122)
! Subprogram not used       integer nf_elatedef       ! Attempt to define var properties, like deflate, after enddef. 
! Subprogram not used       parameter (nf_elatedef = -123)
! Subprogram not used       integer nf_edimscale      ! Probem with HDF5 dimscales. 
! Subprogram not used       parameter (nf_edimscale = -124)
! Subprogram not used       integer nf_enogrp       ! No group found.
! Subprogram not used       parameter (nf_enogrp = -125)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used !     New functions.
! Subprogram not used 
! Subprogram not used !     Parallel I/O.
! Subprogram not used       integer nf_create_par
! Subprogram not used       external nf_create_par
! Subprogram not used 
! Subprogram not used       integer nf_open_par
! Subprogram not used       external nf_open_par
! Subprogram not used 
! Subprogram not used       integer nf_var_par_access
! Subprogram not used       external nf_var_par_access
! Subprogram not used 
! Subprogram not used !     Functions to handle groups.
! Subprogram not used       integer nf_inq_ncid
! Subprogram not used       external nf_inq_ncid
! Subprogram not used 
! Subprogram not used       integer nf_inq_grps
! Subprogram not used       external nf_inq_grps
! Subprogram not used 
! Subprogram not used       integer nf_inq_grpname
! Subprogram not used       external nf_inq_grpname
! Subprogram not used 
! Subprogram not used       integer nf_inq_grpname_full
! Subprogram not used       external nf_inq_grpname_full
! Subprogram not used 
! Subprogram not used       integer nf_inq_grpname_len
! Subprogram not used       external nf_inq_grpname_len
! Subprogram not used 
! Subprogram not used       integer nf_inq_grp_parent
! Subprogram not used       external nf_inq_grp_parent
! Subprogram not used 
! Subprogram not used       integer nf_inq_grp_ncid
! Subprogram not used       external nf_inq_grp_ncid
! Subprogram not used 
! Subprogram not used       integer nf_inq_grp_full_ncid
! Subprogram not used       external nf_inq_grp_full_ncid
! Subprogram not used 
! Subprogram not used       integer nf_inq_varids
! Subprogram not used       external nf_inq_varids
! Subprogram not used 
! Subprogram not used       integer nf_inq_dimids
! Subprogram not used       external nf_inq_dimids
! Subprogram not used 
! Subprogram not used       integer nf_def_grp
! Subprogram not used       external nf_def_grp
! Subprogram not used 
! Subprogram not used !     New rename grp function
! Subprogram not used 
! Subprogram not used       integer nf_rename_grp
! Subprogram not used       external nf_rename_grp
! Subprogram not used 
! Subprogram not used !     New options for netCDF variables.
! Subprogram not used       integer nf_def_var_deflate
! Subprogram not used       external nf_def_var_deflate
! Subprogram not used 
! Subprogram not used       integer nf_inq_var_deflate
! Subprogram not used       external nf_inq_var_deflate
! Subprogram not used 
! Subprogram not used       integer nf_def_var_fletcher32
! Subprogram not used       external nf_def_var_fletcher32
! Subprogram not used 
! Subprogram not used       integer nf_inq_var_fletcher32
! Subprogram not used       external nf_inq_var_fletcher32
! Subprogram not used 
! Subprogram not used       integer nf_def_var_chunking
! Subprogram not used       external nf_def_var_chunking
! Subprogram not used 
! Subprogram not used       integer nf_inq_var_chunking
! Subprogram not used       external nf_inq_var_chunking
! Subprogram not used 
! Subprogram not used       integer nf_def_var_fill
! Subprogram not used       external nf_def_var_fill
! Subprogram not used 
! Subprogram not used       integer nf_inq_var_fill
! Subprogram not used       external nf_inq_var_fill
! Subprogram not used 
! Subprogram not used       integer nf_def_var_endian
! Subprogram not used       external nf_def_var_endian
! Subprogram not used 
! Subprogram not used       integer nf_inq_var_endian
! Subprogram not used       external nf_inq_var_endian
! Subprogram not used 
! Subprogram not used !     User defined types.
! Subprogram not used       integer nf_inq_typeids
! Subprogram not used       external nf_inq_typeids
! Subprogram not used 
! Subprogram not used       integer nf_inq_typeid
! Subprogram not used       external nf_inq_typeid
! Subprogram not used 
! Subprogram not used       integer nf_inq_type
! Subprogram not used       external nf_inq_type
! Subprogram not used 
! Subprogram not used       integer nf_inq_user_type
! Subprogram not used       external nf_inq_user_type
! Subprogram not used 
! Subprogram not used !     User defined types - compound types.
! Subprogram not used       integer nf_def_compound
! Subprogram not used       external nf_def_compound
! Subprogram not used 
! Subprogram not used       integer nf_insert_compound
! Subprogram not used       external nf_insert_compound
! Subprogram not used 
! Subprogram not used       integer nf_insert_array_compound
! Subprogram not used       external nf_insert_array_compound
! Subprogram not used 
! Subprogram not used       integer nf_inq_compound
! Subprogram not used       external nf_inq_compound
! Subprogram not used 
! Subprogram not used       integer nf_inq_compound_name
! Subprogram not used       external nf_inq_compound_name
! Subprogram not used 
! Subprogram not used       integer nf_inq_compound_size
! Subprogram not used       external nf_inq_compound_size
! Subprogram not used 
! Subprogram not used       integer nf_inq_compound_nfields
! Subprogram not used       external nf_inq_compound_nfields
! Subprogram not used 
! Subprogram not used       integer nf_inq_compound_field
! Subprogram not used       external nf_inq_compound_field
! Subprogram not used 
! Subprogram not used       integer nf_inq_compound_fieldname
! Subprogram not used       external nf_inq_compound_fieldname
! Subprogram not used 
! Subprogram not used       integer nf_inq_compound_fieldindex
! Subprogram not used       external nf_inq_compound_fieldindex
! Subprogram not used 
! Subprogram not used       integer nf_inq_compound_fieldoffset
! Subprogram not used       external nf_inq_compound_fieldoffset
! Subprogram not used 
! Subprogram not used       integer nf_inq_compound_fieldtype
! Subprogram not used       external nf_inq_compound_fieldtype
! Subprogram not used 
! Subprogram not used       integer nf_inq_compound_fieldndims
! Subprogram not used       external nf_inq_compound_fieldndims
! Subprogram not used 
! Subprogram not used       integer nf_inq_compound_fielddim_sizes
! Subprogram not used       external nf_inq_compound_fielddim_sizes
! Subprogram not used 
! Subprogram not used !     User defined types - variable length arrays.
! Subprogram not used       integer nf_def_vlen
! Subprogram not used       external nf_def_vlen
! Subprogram not used 
! Subprogram not used       integer nf_inq_vlen
! Subprogram not used       external nf_inq_vlen
! Subprogram not used 
! Subprogram not used       integer nf_free_vlen
! Subprogram not used       external nf_free_vlen
! Subprogram not used 
! Subprogram not used !     User defined types - enums.
! Subprogram not used       integer nf_def_enum
! Subprogram not used       external nf_def_enum
! Subprogram not used 
! Subprogram not used       integer nf_insert_enum
! Subprogram not used       external nf_insert_enum
! Subprogram not used 
! Subprogram not used       integer nf_inq_enum
! Subprogram not used       external nf_inq_enum
! Subprogram not used 
! Subprogram not used       integer nf_inq_enum_member
! Subprogram not used       external nf_inq_enum_member
! Subprogram not used 
! Subprogram not used       integer nf_inq_enum_ident
! Subprogram not used       external nf_inq_enum_ident
! Subprogram not used 
! Subprogram not used !     User defined types - opaque.
! Subprogram not used       integer nf_def_opaque
! Subprogram not used       external nf_def_opaque
! Subprogram not used 
! Subprogram not used       integer nf_inq_opaque
! Subprogram not used       external nf_inq_opaque
! Subprogram not used 
! Subprogram not used !     Write and read attributes of any type, including user defined
! Subprogram not used !     types.
! Subprogram not used       integer nf_put_att
! Subprogram not used       external nf_put_att
! Subprogram not used       integer nf_get_att
! Subprogram not used       external nf_get_att
! Subprogram not used 
! Subprogram not used !     Write and read variables of any type, including user defined
! Subprogram not used !     types.
! Subprogram not used       integer nf_put_var
! Subprogram not used       external nf_put_var
! Subprogram not used       integer nf_put_var1
! Subprogram not used       external nf_put_var1
! Subprogram not used       integer nf_put_vara
! Subprogram not used       external nf_put_vara
! Subprogram not used       integer nf_put_vars
! Subprogram not used       external nf_put_vars
! Subprogram not used       integer nf_get_var
! Subprogram not used       external nf_get_var
! Subprogram not used       integer nf_get_var1
! Subprogram not used       external nf_get_var1
! Subprogram not used       integer nf_get_vara
! Subprogram not used       external nf_get_vara
! Subprogram not used       integer nf_get_vars
! Subprogram not used       external nf_get_vars
! Subprogram not used 
! Subprogram not used !     64-bit int functions.
! Subprogram not used       integer nf_put_var1_int64
! Subprogram not used       external nf_put_var1_int64
! Subprogram not used       integer nf_put_vara_int64
! Subprogram not used       external nf_put_vara_int64
! Subprogram not used       integer nf_put_vars_int64
! Subprogram not used       external nf_put_vars_int64
! Subprogram not used       integer nf_put_varm_int64
! Subprogram not used       external nf_put_varm_int64
! Subprogram not used       integer nf_put_var_int64
! Subprogram not used       external nf_put_var_int64
! Subprogram not used       integer nf_get_var1_int64
! Subprogram not used       external nf_get_var1_int64
! Subprogram not used       integer nf_get_vara_int64
! Subprogram not used       external nf_get_vara_int64
! Subprogram not used       integer nf_get_vars_int64
! Subprogram not used       external nf_get_vars_int64
! Subprogram not used       integer nf_get_varm_int64
! Subprogram not used       external nf_get_varm_int64
! Subprogram not used       integer nf_get_var_int64
! Subprogram not used       external nf_get_var_int64
! Subprogram not used 
! Subprogram not used !     For helping F77 users with VLENs.
! Subprogram not used       integer nf_get_vlen_element
! Subprogram not used       external nf_get_vlen_element
! Subprogram not used       integer nf_put_vlen_element
! Subprogram not used       external nf_put_vlen_element
! Subprogram not used 
! Subprogram not used !     For dealing with file level chunk cache.
! Subprogram not used       integer nf_set_chunk_cache
! Subprogram not used       external nf_set_chunk_cache
! Subprogram not used       integer nf_get_chunk_cache
! Subprogram not used       external nf_get_chunk_cache
! Subprogram not used 
! Subprogram not used !     For dealing with per variable chunk cache.
! Subprogram not used       integer nf_set_var_chunk_cache
! Subprogram not used       external nf_set_var_chunk_cache
! Subprogram not used       integer nf_get_var_chunk_cache
! Subprogram not used       external nf_get_var_chunk_cache
! Subprogram not used 
! Subprogram not used !     NetCDF-2.
! Subprogram not used !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Subprogram not used ! begin netcdf 2.4 backward compatibility:
! Subprogram not used !
! Subprogram not used 
! Subprogram not used !      
! Subprogram not used ! functions in the fortran interface
! Subprogram not used !
! Subprogram not used       integer nccre
! Subprogram not used       integer ncopn
! Subprogram not used       integer ncddef
! Subprogram not used       integer ncdid
! Subprogram not used       integer ncvdef
! Subprogram not used       integer ncvid
! Subprogram not used       integer nctlen
! Subprogram not used       integer ncsfil
! Subprogram not used 
! Subprogram not used       external nccre
! Subprogram not used       external ncopn
! Subprogram not used       external ncddef
! Subprogram not used       external ncdid
! Subprogram not used       external ncvdef
! Subprogram not used       external ncvid
! Subprogram not used       external nctlen
! Subprogram not used       external ncsfil
! Subprogram not used 
! Subprogram not used 
! Subprogram not used       integer ncrdwr
! Subprogram not used       integer nccreat
! Subprogram not used       integer ncexcl
! Subprogram not used       integer ncindef
! Subprogram not used       integer ncnsync
! Subprogram not used       integer nchsync
! Subprogram not used       integer ncndirty
! Subprogram not used       integer nchdirty
! Subprogram not used       integer nclink
! Subprogram not used       integer ncnowrit
! Subprogram not used       integer ncwrite
! Subprogram not used       integer ncclob
! Subprogram not used       integer ncnoclob
! Subprogram not used       integer ncglobal
! Subprogram not used       integer ncfill
! Subprogram not used       integer ncnofill
! Subprogram not used       integer maxncop
! Subprogram not used       integer maxncdim
! Subprogram not used       integer maxncatt
! Subprogram not used       integer maxncvar
! Subprogram not used       integer maxncnam
! Subprogram not used       integer maxvdims
! Subprogram not used       integer ncnoerr
! Subprogram not used       integer ncebadid
! Subprogram not used       integer ncenfile
! Subprogram not used       integer nceexist
! Subprogram not used       integer nceinval
! Subprogram not used       integer nceperm
! Subprogram not used       integer ncenotin
! Subprogram not used       integer nceindef
! Subprogram not used       integer ncecoord
! Subprogram not used       integer ncemaxds
! Subprogram not used       integer ncename
! Subprogram not used       integer ncenoatt
! Subprogram not used       integer ncemaxat
! Subprogram not used       integer ncebadty
! Subprogram not used       integer ncebadd
! Subprogram not used       integer ncests
! Subprogram not used       integer nceunlim
! Subprogram not used       integer ncemaxvs
! Subprogram not used       integer ncenotvr
! Subprogram not used       integer nceglob
! Subprogram not used       integer ncenotnc
! Subprogram not used       integer ncfoobar
! Subprogram not used       integer ncsyserr
! Subprogram not used       integer ncfatal
! Subprogram not used       integer ncverbos
! Subprogram not used       integer ncentool
! Subprogram not used 
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! netcdf data types:
! Subprogram not used !
! Subprogram not used       integer ncbyte
! Subprogram not used       integer ncchar
! Subprogram not used       integer ncshort
! Subprogram not used       integer nclong
! Subprogram not used       integer ncfloat
! Subprogram not used       integer ncdouble
! Subprogram not used 
! Subprogram not used       parameter(ncbyte = 1)
! Subprogram not used       parameter(ncchar = 2)
! Subprogram not used       parameter(ncshort = 3)
! Subprogram not used       parameter(nclong = 4)
! Subprogram not used       parameter(ncfloat = 5)
! Subprogram not used       parameter(ncdouble = 6)
! Subprogram not used 
! Subprogram not used !     
! Subprogram not used !     masks for the struct nc flag field; passed in as 'mode' arg to
! Subprogram not used !     nccreate and ncopen.
! Subprogram not used !     
! Subprogram not used 
! Subprogram not used !     read/write, 0 => readonly 
! Subprogram not used       parameter(ncrdwr = 1)
! Subprogram not used !     in create phase, cleared by ncendef 
! Subprogram not used       parameter(nccreat = 2)
! Subprogram not used !     on create destroy existing file 
! Subprogram not used       parameter(ncexcl = 4)
! Subprogram not used !     in define mode, cleared by ncendef 
! Subprogram not used       parameter(ncindef = 8)
! Subprogram not used !     synchronise numrecs on change (x'10')
! Subprogram not used       parameter(ncnsync = 16)
! Subprogram not used !     synchronise whole header on change (x'20')
! Subprogram not used       parameter(nchsync = 32)
! Subprogram not used !     numrecs has changed (x'40')
! Subprogram not used       parameter(ncndirty = 64)  
! Subprogram not used !     header info has changed (x'80')
! Subprogram not used       parameter(nchdirty = 128)
! Subprogram not used !     prefill vars on endef and increase of record, the default behavior
! Subprogram not used       parameter(ncfill = 0)
! Subprogram not used !     do not fill vars on endef and increase of record (x'100')
! Subprogram not used       parameter(ncnofill = 256)
! Subprogram not used !     isa link (x'8000')
! Subprogram not used       parameter(nclink = 32768)
! Subprogram not used 
! Subprogram not used !     
! Subprogram not used !     'mode' arguments for nccreate and ncopen
! Subprogram not used !     
! Subprogram not used       parameter(ncnowrit = 0)
! Subprogram not used       parameter(ncwrite = ncrdwr)
! Subprogram not used       parameter(ncclob = nf_clobber)
! Subprogram not used       parameter(ncnoclob = nf_noclobber)
! Subprogram not used 
! Subprogram not used !     
! Subprogram not used !     'size' argument to ncdimdef for an unlimited dimension
! Subprogram not used !     
! Subprogram not used       integer ncunlim
! Subprogram not used       parameter(ncunlim = 0)
! Subprogram not used 
! Subprogram not used !     
! Subprogram not used !     attribute id to put/get a global attribute
! Subprogram not used !     
! Subprogram not used       parameter(ncglobal  = 0)
! Subprogram not used 
! Subprogram not used !     
! Subprogram not used !     advisory maximums:
! Subprogram not used !     
! Subprogram not used       parameter(maxncop = 64)
! Subprogram not used       parameter(maxncdim = 1024)
! Subprogram not used       parameter(maxncatt = 8192)
! Subprogram not used       parameter(maxncvar = 8192)
! Subprogram not used !     not enforced 
! Subprogram not used       parameter(maxncnam = 256)
! Subprogram not used       parameter(maxvdims = maxncdim)
! Subprogram not used 
! Subprogram not used !     
! Subprogram not used !     global netcdf error status variable
! Subprogram not used !     initialized in error.c
! Subprogram not used !     
! Subprogram not used 
! Subprogram not used !     no error 
! Subprogram not used       parameter(ncnoerr = nf_noerr)
! Subprogram not used !     not a netcdf id 
! Subprogram not used       parameter(ncebadid = nf_ebadid)
! Subprogram not used !     too many netcdfs open 
! Subprogram not used       parameter(ncenfile = -31)   ! nc_syserr
! Subprogram not used !     netcdf file exists && ncnoclob
! Subprogram not used       parameter(nceexist = nf_eexist)
! Subprogram not used !     invalid argument 
! Subprogram not used       parameter(nceinval = nf_einval)
! Subprogram not used !     write to read only 
! Subprogram not used       parameter(nceperm = nf_eperm)
! Subprogram not used !     operation not allowed in data mode 
! Subprogram not used       parameter(ncenotin = nf_enotindefine )   
! Subprogram not used !     operation not allowed in define mode 
! Subprogram not used       parameter(nceindef = nf_eindefine)   
! Subprogram not used !     coordinates out of domain 
! Subprogram not used       parameter(ncecoord = nf_einvalcoords)
! Subprogram not used !     maxncdims exceeded 
! Subprogram not used       parameter(ncemaxds = nf_emaxdims)
! Subprogram not used !     string match to name in use 
! Subprogram not used       parameter(ncename = nf_enameinuse)   
! Subprogram not used !     attribute not found 
! Subprogram not used       parameter(ncenoatt = nf_enotatt)
! Subprogram not used !     maxncattrs exceeded 
! Subprogram not used       parameter(ncemaxat = nf_emaxatts)
! Subprogram not used !     not a netcdf data type 
! Subprogram not used       parameter(ncebadty = nf_ebadtype)
! Subprogram not used !     invalid dimension id 
! Subprogram not used       parameter(ncebadd = nf_ebaddim)
! Subprogram not used !     ncunlimited in the wrong index 
! Subprogram not used       parameter(nceunlim = nf_eunlimpos)
! Subprogram not used !     maxncvars exceeded 
! Subprogram not used       parameter(ncemaxvs = nf_emaxvars)
! Subprogram not used !     variable not found 
! Subprogram not used       parameter(ncenotvr = nf_enotvar)
! Subprogram not used !     action prohibited on ncglobal varid 
! Subprogram not used       parameter(nceglob = nf_eglobal)
! Subprogram not used !     not a netcdf file 
! Subprogram not used       parameter(ncenotnc = nf_enotnc)
! Subprogram not used       parameter(ncests = nf_ests)
! Subprogram not used       parameter (ncentool = nf_emaxname) 
! Subprogram not used       parameter(ncfoobar = 32)
! Subprogram not used       parameter(ncsyserr = -31)
! Subprogram not used 
! Subprogram not used !     
! Subprogram not used !     global options variable. used to determine behavior of error handler.
! Subprogram not used !     initialized in lerror.c
! Subprogram not used !     
! Subprogram not used       parameter(ncfatal = 1)
! Subprogram not used       parameter(ncverbos = 2)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used !     default fill values.  these must be the same as in the c interface.
! Subprogram not used !
! Subprogram not used       integer filbyte
! Subprogram not used       integer filchar
! Subprogram not used       integer filshort
! Subprogram not used       integer fillong
! Subprogram not used       real filfloat
! Subprogram not used       doubleprecision fildoub
! Subprogram not used 
! Subprogram not used       parameter (filbyte = -127)
! Subprogram not used       parameter (filchar = 0)
! Subprogram not used       parameter (filshort = -32767)
! Subprogram not used       parameter (fillong = -2147483647)
! Subprogram not used       parameter (filfloat = 9.9692099683868690e+36)
! Subprogram not used       parameter (fildoub = 9.9692099683868690e+36)
! Subprogram not used 
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    type(mct_sMat),intent(inout)  :: sMat
! Subprogram not used    character(*),intent(in)  :: filename  ! netCDF file to read
! Subprogram not used 
! Subprogram not used !EOP
! Subprogram not used 
! Subprogram not used    !--- local ---
! Subprogram not used    integer(IN)           :: n       ! generic loop indicies
! Subprogram not used    integer(IN)           :: na      ! size of source domain
! Subprogram not used    integer(IN)           :: nb      ! size of destination domain
! Subprogram not used    integer(IN)           :: ns      ! number of non-zero elements in matrix
! Subprogram not used    integer(IN)           :: ni,nj   ! number of row and col in the matrix
! Subprogram not used    integer(IN)           :: igrow   ! aVect index for matrix row
! Subprogram not used    integer(IN)           :: igcol   ! aVect index for matrix column
! Subprogram not used    integer(IN)           :: iwgt    ! aVect index for matrix element
! Subprogram not used 
! Subprogram not used    real(R8)   ,allocatable :: rtemp(:)  ! reals
! Subprogram not used    integer(IN),allocatable :: itemp(:)  ! ints
! Subprogram not used 
! Subprogram not used    integer(IN)           :: rcode   ! netCDF routine return code
! Subprogram not used    integer(IN)           :: fid     ! netCDF file      ID
! Subprogram not used    integer(IN)           :: vid     ! netCDF variable  ID
! Subprogram not used    integer(IN)           :: did     ! netCDF dimension ID
! Subprogram not used 
! Subprogram not used    character(*),parameter :: subName = '(shr_mct_sMatReadnc) '
! Subprogram not used    character(*),parameter :: F00 = "('(shr_mct_sMatReadnc) ',4a)"
! Subprogram not used    character(*),parameter :: F01 = '("(shr_mct_sMatReadnc) ",2(a,i9))'
! Subprogram not used 
! Subprogram not used    if (s_loglev > 0) write(s_logunit,F00) "reading mapping matrix data..."
! Subprogram not used 
! Subprogram not used    !----------------------------------------------------------------------------
! Subprogram not used    ! open & read the file
! Subprogram not used    !----------------------------------------------------------------------------
! Subprogram not used    if (s_loglev > 0) write(s_logunit,F00) "* file name                  : ",trim(fileName)
! Subprogram not used    rcode = nf_open(filename,NF_NOWRITE,fid)
! Subprogram not used    if (rcode /= NF_NOERR) then
! Subprogram not used       write(s_logunit,F00) nf_strerror(rcode)
! Subprogram not used       call mct_die(subName,"error opening Netcdf file")
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    !--- allocate memory & get matrix data ----------
! Subprogram not used    rcode = nf_inq_dimid (fid, 'n_s', did)  ! size of sparse matrix
! Subprogram not used    rcode = nf_inq_dimlen(fid, did  , ns)
! Subprogram not used    rcode = nf_inq_dimid (fid, 'n_a', did)  ! size of  input vector
! Subprogram not used    rcode = nf_inq_dimlen(fid, did  , na)
! Subprogram not used    rcode = nf_inq_dimid (fid, 'n_b', did)  ! size of output vector
! Subprogram not used    rcode = nf_inq_dimlen(fid, did  , nb)
! Subprogram not used 
! Subprogram not used    if (s_loglev > 0) write(s_logunit,F01) "* matrix dimensions src x dst: ",na,' x',nb
! Subprogram not used    if (s_loglev > 0) write(s_logunit,F01) "* number of non-zero elements: ",ns
! Subprogram not used 
! Subprogram not used    !----------------------------------------------------------------------------
! Subprogram not used    ! init the mct sMat data type
! Subprogram not used    !----------------------------------------------------------------------------
! Subprogram not used    ! mct_sMat_init must be given the number of rows and columns that
! Subprogram not used    ! would be in the full matrix.  Nrows= size of output vector=nb.
! Subprogram not used    ! Ncols = size of input vector = na.
! Subprogram not used    call mct_sMat_init(sMat, nb, na, ns)
! Subprogram not used 
! Subprogram not used    igrow = mct_sMat_indexIA(sMat,'grow')
! Subprogram not used    igcol = mct_sMat_indexIA(sMat,'gcol')
! Subprogram not used    iwgt  = mct_sMat_indexRA(sMat,'weight')
! Subprogram not used 
! Subprogram not used    !!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used    ! read and load matrix weights
! Subprogram not used    allocate(rtemp(ns),stat=rcode)
! Subprogram not used    if (rcode /= 0) &
! Subprogram not used      call mct_die(subName,':: allocate weights',rcode)
! Subprogram not used 
! Subprogram not used    rcode = nf_inq_varid     (fid,'S'  ,vid)
! Subprogram not used    rcode = nf_get_var_double(fid,vid  ,rtemp  )
! Subprogram not used    if (rcode /= NF_NOERR .and. s_loglev > 0) write(s_logunit,F00) nf_strerror(rcode)
! Subprogram not used 
! Subprogram not used    sMat%data%rAttr(iwgt ,:) =   rtemp(:)
! Subprogram not used 
! Subprogram not used    deallocate(rtemp, stat=rcode)
! Subprogram not used    if (rcode /= 0) call mct_perr_die(subName,':: deallocate weights',rcode)
! Subprogram not used 
! Subprogram not used    !!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used    ! read and load rows
! Subprogram not used    allocate(itemp(ns),stat=rcode)
! Subprogram not used    if (rcode /= 0) call mct_perr_die(subName,':: allocate rows',rcode)
! Subprogram not used 
! Subprogram not used    rcode = nf_inq_varid     (fid,'row',vid)
! Subprogram not used    rcode = nf_get_var_int   (fid,vid  ,itemp)
! Subprogram not used    if (rcode /= NF_NOERR .and. s_loglev > 0) write(s_logunit,F00) nf_strerror(rcode)
! Subprogram not used 
! Subprogram not used    sMat%data%iAttr(igrow,:) = itemp(:)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used    !!!!!!!!!!!!!!!!!!!!!!!!!!
! Subprogram not used    ! read and load columns
! Subprogram not used    itemp(:) = 0
! Subprogram not used 
! Subprogram not used    rcode = nf_inq_varid     (fid,'col',vid)
! Subprogram not used    rcode = nf_get_var_int   (fid,vid  ,itemp)
! Subprogram not used    if (rcode /= NF_NOERR .and. s_loglev > 0) write(s_logunit,F00) nf_strerror(rcode)
! Subprogram not used 
! Subprogram not used    sMat%data%iAttr(igcol,:) = itemp(:)
! Subprogram not used 
! Subprogram not used    deallocate(itemp, stat=rcode)
! Subprogram not used    if (rcode /= 0) call mct_perr_die(subName,':: deallocate cols',rcode)
! Subprogram not used 
! Subprogram not used    rcode = nf_close(fid)
! Subprogram not used 
! Subprogram not used    if (s_loglev > 0) write(s_logunit,F00) "... done reading file"
! Subprogram not used    call shr_sys_flush(s_logunit)
! Subprogram not used 
! Subprogram not used end subroutine shr_mct_sMatReadnc

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_mct_queryConfigFile - get mct config file info
!
! !DESCRIPTION:
!   Query MCT config file variables
!
! !REMARKS:
!
! !REVISION HISTORY:
!     2013 Aug 17: T. Craig
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_mct_queryConfigFile(mpicom, ConfigFileName, &
           Label1,Value1,Label2,Value2,Label3,Value3)

! !INPUT/OUTPUT PARAMETERS:
   integer          ,intent(in)  :: mpicom
   character(len=*), intent(in)  :: ConfigFileName
   character(len=*), intent(in)  :: Label1
   character(len=*), intent(out) :: Value1
   character(len=*), intent(in) ,optional :: Label2
   character(len=*), intent(out),optional :: Value2
   character(len=*), intent(in) ,optional :: Label3
   character(len=*), intent(out),optional :: Value3

!EOP
   integer :: iret
   character(*),parameter :: subName = '(shr_mct_queryConfigFile) '

   call I90_allLoadF(ConfigFileName,0,mpicom,iret)
   if(iret /= 0) then
      write(s_logunit,*) trim(subname),"Cant find config file ",ConfigFileName
      call shr_sys_abort(trim(subname)//' File Not Found')
   endif

   call i90_label(trim(Label1),iret)
   if(iret /= 0) then
      write(s_logunit,*) trim(subname),"Cant find label ",Label1
      call shr_sys_abort(trim(subname)//' Label1 Not Found')
   endif

   call i90_gtoken(Value1,iret)
   if(iret /= 0) then
      write(s_logunit,*) trim(subname),"Error reading token ",Value1
      call shr_sys_abort(trim(subname)//' Error on read value1')
   endif

   if (present(Label2) .and. present(Value2)) then

      call i90_label(trim(Label2),iret)
      if(iret /= 0) then
         write(s_logunit,*) trim(subname),"Cant find label ",Label2
         call shr_sys_abort(trim(subname)//' Label2 Not Found')
      endif

      call i90_gtoken(Value2,iret)
      if(iret /= 0) then
         write(s_logunit,*)"Error reading token ",Value2
         call shr_sys_abort(trim(subname)//' Error on read value2')
      endif

   endif

   if (present(Label3) .and. present(Value3)) then

      call i90_label(trim(Label3),iret)
      if(iret /= 0) then
         write(s_logunit,*) trim(subname),"Cant find label ",Label3
         call shr_sys_abort(trim(subname)//' Label3 Not Found')
      endif

      call i90_gtoken(Value3,iret)
      if(iret /= 0) then
         write(s_logunit,*)"Error reading token ",Value3
         call shr_sys_abort(trim(subname)//' Error on read value3')
      endif

   endif

   call I90_Release(iret)

end subroutine shr_mct_queryConfigFile

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_mct_sMatPInitnc_mapfile - initialize a SparseMatrixPlus.
!
! !DESCRIPTION:
!   Read in mapping matrix data from a SCRIP netCDF data file in first an
!   Smat and then an SMatPlus
!
! !REMARKS:
!
! !REVISION HISTORY:
!     2012 Feb 27: M. Vertenstein
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_mct_sMatPInitnc_mapfile(sMatP, gsMapX, gsMapY, &
                                       filename, maptype, mpicom, &
                                       ni_i, nj_i, ni_o, nj_o, &
                                       areasrc, areadst)

! !INPUT/OUTPUT PARAMETERS:

   type(mct_sMatP),intent(inout)         :: sMatP
   type(mct_gsMap),intent(in)            :: gsMapX
   type(mct_gsMap),intent(in)            :: gsMapY
   character(*)   ,intent(in)            :: filename        ! scrip map file to read
   character(*)   ,intent(in)            :: maptype         ! map type
   integer        ,intent(in)            :: mpicom
   integer        ,intent(out), optional :: ni_i            ! number of longitudes on input grid
   integer        ,intent(out), optional :: nj_i            ! number of latitudes  on input grid
   integer        ,intent(out), optional :: ni_o            ! number of longitudes on output grid
   integer        ,intent(out), optional :: nj_o            ! number of latitudes  on output grid
   type(mct_Avect),intent(out), optional :: areasrc         ! area of src grid from mapping file
   type(mct_Avect),intent(out), optional :: areadst         ! area of src grid from mapping file

!EOP
   type(mct_sMat ) :: sMati    ! initial sMat from read (either root or decomp)
   type(mct_Avect) :: areasrc_map ! area of src grid from mapping file
   type(mct_Avect) :: areadst_map ! area of dst grid from mapping file

   integer          :: lsize
   integer          :: iret
   integer          :: pe_loc
   logical          :: usevector 
   character(len=3) :: Smaptype
   character(*),parameter :: areaAV_field = 'aream'
   character(*),parameter :: F00 = "('(shr_mct_sMatPInitnc) ',4a)"
   character(*),parameter :: F01 = "('(shr_mct_sMatPInitnc) ',a,i10)"

   call shr_mpi_commrank(mpicom,pe_loc)

   if (s_loglev > 0) write(s_logunit,*) " "
   if (s_loglev > 0) write(s_logunit,F00) "Initializing SparseMatrixPlus"
   if (s_loglev > 0) write(s_logunit,F00) "SmatP mapname ",trim(filename)
   if (s_loglev > 0) write(s_logunit,F00) "SmatP maptype ",trim(maptype)

   if (maptype == "X") then
      Smaptype = "src"
   else if(maptype == "Y") then
      Smaptype = "dst"
   end if

   call shr_mpi_commrank(mpicom, pe_loc)

   lsize = mct_gsMap_lsize(gsMapX, mpicom)
   call mct_aVect_init(areasrc_map, rList=areaAV_field, lsize=lsize)

   lsize = mct_gsMap_lsize(gsMapY, mpicom)
   call mct_aVect_init(areadst_map, rList=areaAV_field, lsize=lsize)

   if (present(ni_i) .and. present(nj_i) .and. present(ni_o) .and. present(nj_o)) then
      call shr_mct_sMatReaddnc(sMati, gsMapX, gsMapY, Smaptype, areasrc_map, areadst_map, &
           fileName, pe_loc, mpicom, ni_i, nj_i, ni_o, nj_o)
   else
      call shr_mct_sMatReaddnc(sMati, gsMapX, gsMapY, Smaptype, areasrc_map, areadst_map, &
           fileName, pe_loc, mpicom)
   end if
   call mct_sMatP_Init(sMatP, sMati, gsMapX, gsMapY, 0, mpicom, gsMapX%comp_id)


   lsize = mct_smat_gNumEl(sMatP%Matrix,mpicom)
   if (s_loglev > 0) write(s_logunit,F01) "Done initializing SmatP, nElements = ",lsize

   usevector = .false.
   if (present(areasrc)) then
      call mct_aVect_copy(aVin=areasrc_map, aVout=areasrc, vector=usevector)
   end if
   if (present(areadst)) then
      call mct_aVect_copy(aVin=areadst_map, aVout=areadst, vector=usevector)
   end if

   call mct_aVect_clean(areasrc_map)
   call mct_aVect_clean(areadst_map)

   call mct_sMat_Clean(sMati)

end subroutine shr_mct_sMatPInitnc_mapfile

!BOP ===========================================================================
!
! !IROUTINE:  shr_mct_sMatReaddnc - Do a distributed read of a NetCDF SCRIP file and
!                                return weights in a distributed SparseMatrix
!
! !DESCRIPTION: 
!     Read in mapping matrix data from a SCRIP netCDF data file using
!     a low memory method and then scatter to all pes.
!
! !REMARKS:
!   This routine leverages gsmaps to determine scatter pattern
!   The scatter is implemented as a bcast of all weights then a local
!     computation on each pe to determine with weights to keep based
!     on gsmap information.
!   The algorithm to determine whether a weight belongs on a pe involves
!     creating a couple local arrays (lsstart and lscount) which are
!     the local values of start and length from the gsmap.  these are
!     sorted via a bubble sort and then searched via a binary search
!     to check whether a global index is on the local pe.
!   The local buffer sizes are estimated up front based on ngridcell/npes
!     plus 20% (see 1.2 below).  If the local buffer size fills up, then
!     the buffer is reallocated 50% large (see 1.5 below) and the fill
!     continues.  The idea is to trade off memory reallocation and copy
!     with memory usage.  1.2 and 1.5 are arbitary, other values may
!     result in better performance.
!   Once all the matrix weights have been read, the sMat is initialized,
!     the values from the buffers are copied in, and everything is deallocated.

! !SEE ALSO: 
!    mct/m_SparseMatrix.F90 (MCT source code)
!
! !REVISION HISTORY: 
!    2007-Jan-18 - T. Craig -- first version
!    2007-Mar-20 - R. Jacob -- rename to shr_mct_sMatReaddnc.  Remove use of cpl_
!                  variables and move to shr_mct_mod
! 
! !INTERFACE:  -----------------------------------------------------------------

subroutine shr_mct_sMatReaddnc(sMat,SgsMap,DgsMap,newdom,areasrc,areadst, &
                            fileName,mytask, mpicom, ni_i,nj_i,ni_o,nj_o )

! !USES:

!     NetCDF-3.
!
! netcdf version 3 fortran interface:
!

!
! external netcdf data types:
!
      integer nf_byte
      integer nf_int1
      integer nf_char
      integer nf_short
      integer nf_int2
      integer nf_int
      integer nf_float
      integer nf_real
      integer nf_double

      parameter (nf_byte = 1)
      parameter (nf_int1 = nf_byte)
      parameter (nf_char = 2)
      parameter (nf_short = 3)
      parameter (nf_int2 = nf_short)
      parameter (nf_int = 4)
      parameter (nf_float = 5)
      parameter (nf_real = nf_float)
      parameter (nf_double = 6)

!
! default fill values:
!
      integer           nf_fill_byte
      integer           nf_fill_int1
      integer           nf_fill_char
      integer           nf_fill_short
      integer           nf_fill_int2
      integer           nf_fill_int
      real              nf_fill_float
      real              nf_fill_real
      doubleprecision   nf_fill_double

      parameter (nf_fill_byte = -127)
      parameter (nf_fill_int1 = nf_fill_byte)
      parameter (nf_fill_char = 0)
      parameter (nf_fill_short = -32767)
      parameter (nf_fill_int2 = nf_fill_short)
      parameter (nf_fill_int = -2147483647)
      parameter (nf_fill_float = 9.9692099683868690e+36)
      parameter (nf_fill_real = nf_fill_float)
      parameter (nf_fill_double = 9.9692099683868690d+36)

!
! mode flags for opening and creating a netcdf dataset:
!
      integer nf_nowrite
      integer nf_write
      integer nf_clobber
      integer nf_noclobber
      integer nf_fill
      integer nf_nofill
      integer nf_lock
      integer nf_share
      integer nf_64bit_offset
      integer nf_sizehint_default
      integer nf_align_chunk
      integer nf_format_classic
      integer nf_format_64bit
      integer nf_diskless
      integer nf_mmap

      parameter (nf_nowrite = 0)
      parameter (nf_write = 1)
      parameter (nf_clobber = 0)
      parameter (nf_noclobber = 4)
      parameter (nf_fill = 0)
      parameter (nf_nofill = 256)
      parameter (nf_lock = 1024)
      parameter (nf_share = 2048)
      parameter (nf_64bit_offset = 512)
      parameter (nf_sizehint_default = 0)
      parameter (nf_align_chunk = -1)
      parameter (nf_format_classic = 1)
      parameter (nf_format_64bit = 2)
      parameter (nf_diskless = 8)
      parameter (nf_mmap = 16)

!
! size argument for defining an unlimited dimension:
!
      integer nf_unlimited
      parameter (nf_unlimited = 0)

!
! global attribute id:
!
      integer nf_global
      parameter (nf_global = 0)

!
! implementation limits:
!
      integer nf_max_dims
      integer nf_max_attrs
      integer nf_max_vars
      integer nf_max_name
      integer nf_max_var_dims

      parameter (nf_max_dims = 1024)
      parameter (nf_max_attrs = 8192)
      parameter (nf_max_vars = 8192)
      parameter (nf_max_name = 256)
      parameter (nf_max_var_dims = nf_max_dims)

!
! error codes:
!
      integer nf_noerr
      integer nf_ebadid
      integer nf_eexist
      integer nf_einval
      integer nf_eperm
      integer nf_enotindefine
      integer nf_eindefine
      integer nf_einvalcoords
      integer nf_emaxdims
      integer nf_enameinuse
      integer nf_enotatt
      integer nf_emaxatts
      integer nf_ebadtype
      integer nf_ebaddim
      integer nf_eunlimpos
      integer nf_emaxvars
      integer nf_enotvar
      integer nf_eglobal
      integer nf_enotnc
      integer nf_ests
      integer nf_emaxname
      integer nf_eunlimit
      integer nf_enorecvars
      integer nf_echar
      integer nf_eedge
      integer nf_estride
      integer nf_ebadname
      integer nf_erange
      integer nf_enomem
      integer nf_evarsize
      integer nf_edimsize
      integer nf_etrunc

      parameter (nf_noerr = 0)
      parameter (nf_ebadid = -33)
      parameter (nf_eexist = -35)
      parameter (nf_einval = -36)
      parameter (nf_eperm = -37)
      parameter (nf_enotindefine = -38)
      parameter (nf_eindefine = -39)
      parameter (nf_einvalcoords = -40)
      parameter (nf_emaxdims = -41)
      parameter (nf_enameinuse = -42)
      parameter (nf_enotatt = -43)
      parameter (nf_emaxatts = -44)
      parameter (nf_ebadtype = -45)
      parameter (nf_ebaddim = -46)
      parameter (nf_eunlimpos = -47)
      parameter (nf_emaxvars = -48)
      parameter (nf_enotvar = -49)
      parameter (nf_eglobal = -50)
      parameter (nf_enotnc = -51)
      parameter (nf_ests = -52)
      parameter (nf_emaxname = -53)
      parameter (nf_eunlimit = -54)
      parameter (nf_enorecvars = -55)
      parameter (nf_echar = -56)
      parameter (nf_eedge = -57)
      parameter (nf_estride = -58)
      parameter (nf_ebadname = -59)
      parameter (nf_erange = -60)
      parameter (nf_enomem = -61)
      parameter (nf_evarsize = -62)
      parameter (nf_edimsize = -63)
      parameter (nf_etrunc = -64)
!
! error handling modes:
!
      integer  nf_fatal
      integer nf_verbose

      parameter (nf_fatal = 1)
      parameter (nf_verbose = 2)

!
! miscellaneous routines:
!
      character*80   nf_inq_libvers
      external       nf_inq_libvers

      character*80   nf_strerror
!                         (integer             ncerr)
      external       nf_strerror

      logical        nf_issyserr
!                         (integer             ncerr)
      external       nf_issyserr

!
! control routines:
!
      integer         nf_inq_base_pe
!                         (integer             ncid,
!                          integer             pe)
      external        nf_inq_base_pe

      integer         nf_set_base_pe
!                         (integer             ncid,
!                          integer             pe)
      external        nf_set_base_pe

      integer         nf_create
!                         (character*(*)       path,
!                          integer             cmode,
!                          integer             ncid)
      external        nf_create

      integer         nf__create
!                         (character*(*)       path,
!                          integer             cmode,
!                          integer             initialsz,
!                          integer             chunksizehint,
!                          integer             ncid)
      external        nf__create

      integer         nf__create_mp
!                         (character*(*)       path,
!                          integer             cmode,
!                          integer             initialsz,
!                          integer             basepe,
!                          integer             chunksizehint,
!                          integer             ncid)
      external        nf__create_mp

      integer         nf_open
!                         (character*(*)       path,
!                          integer             mode,
!                          integer             ncid)
      external        nf_open

      integer         nf__open
!                         (character*(*)       path,
!                          integer             mode,
!                          integer             chunksizehint,
!                          integer             ncid)
      external        nf__open

      integer         nf__open_mp
!                         (character*(*)       path,
!                          integer             mode,
!                          integer             basepe,
!                          integer             chunksizehint,
!                          integer             ncid)
      external        nf__open_mp

      integer         nf_set_fill
!                         (integer             ncid,
!                          integer             fillmode,
!                          integer             old_mode)
      external        nf_set_fill

      integer         nf_set_default_format
!                          (integer             format,
!                          integer             old_format)
      external        nf_set_default_format

      integer         nf_redef
!                         (integer             ncid)
      external        nf_redef

      integer         nf_enddef
!                         (integer             ncid)
      external        nf_enddef

      integer         nf__enddef
!                         (integer             ncid,
!                          integer             h_minfree,
!                          integer             v_align,
!                          integer             v_minfree,
!                          integer             r_align)
      external        nf__enddef

      integer         nf_sync
!                         (integer             ncid)
      external        nf_sync

      integer         nf_abort
!                         (integer             ncid)
      external        nf_abort

      integer         nf_close
!                         (integer             ncid)
      external        nf_close

      integer         nf_delete
!                         (character*(*)       ncid)
      external        nf_delete

!
! general inquiry routines:
!

      integer         nf_inq
!                         (integer             ncid,
!                          integer             ndims,
!                          integer             nvars,
!                          integer             ngatts,
!                          integer             unlimdimid)
      external        nf_inq

! new inquire path

      integer nf_inq_path
      external nf_inq_path

      integer         nf_inq_ndims
!                         (integer             ncid,
!                          integer             ndims)
      external        nf_inq_ndims

      integer         nf_inq_nvars
!                         (integer             ncid,
!                          integer             nvars)
      external        nf_inq_nvars

      integer         nf_inq_natts
!                         (integer             ncid,
!                          integer             ngatts)
      external        nf_inq_natts

      integer         nf_inq_unlimdim
!                         (integer             ncid,
!                          integer             unlimdimid)
      external        nf_inq_unlimdim

      integer         nf_inq_format
!                         (integer             ncid,
!                          integer             format)
      external        nf_inq_format

!
! dimension routines:
!

      integer         nf_def_dim
!                         (integer             ncid,
!                          character(*)        name,
!                          integer             len,
!                          integer             dimid)
      external        nf_def_dim

      integer         nf_inq_dimid
!                         (integer             ncid,
!                          character(*)        name,
!                          integer             dimid)
      external        nf_inq_dimid

      integer         nf_inq_dim
!                         (integer             ncid,
!                          integer             dimid,
!                          character(*)        name,
!                          integer             len)
      external        nf_inq_dim

      integer         nf_inq_dimname
!                         (integer             ncid,
!                          integer             dimid,
!                          character(*)        name)
      external        nf_inq_dimname

      integer         nf_inq_dimlen
!                         (integer             ncid,
!                          integer             dimid,
!                          integer             len)
      external        nf_inq_dimlen

      integer         nf_rename_dim
!                         (integer             ncid,
!                          integer             dimid,
!                          character(*)        name)
      external        nf_rename_dim

!
! general attribute routines:
!

      integer         nf_inq_att
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len)
      external        nf_inq_att

      integer         nf_inq_attid
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             attnum)
      external        nf_inq_attid

      integer         nf_inq_atttype
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype)
      external        nf_inq_atttype

      integer         nf_inq_attlen
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             len)
      external        nf_inq_attlen

      integer         nf_inq_attname
!                         (integer             ncid,
!                          integer             varid,
!                          integer             attnum,
!                          character(*)        name)
      external        nf_inq_attname

      integer         nf_copy_att
!                         (integer             ncid_in,
!                          integer             varid_in,
!                          character(*)        name,
!                          integer             ncid_out,
!                          integer             varid_out)
      external        nf_copy_att

      integer         nf_rename_att
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        curname,
!                          character(*)        newname)
      external        nf_rename_att

      integer         nf_del_att
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name)
      external        nf_del_att

!
! attribute put/get routines:
!

      integer         nf_put_att_text
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             len,
!                          character(*)        text)
      external        nf_put_att_text

      integer         nf_get_att_text
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          character(*)        text)
      external        nf_get_att_text

      integer         nf_put_att_int1
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          nf_int1_t           i1vals(1))
      external        nf_put_att_int1

      integer         nf_get_att_int1
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          nf_int1_t           i1vals(1))
      external        nf_get_att_int1

      integer         nf_put_att_int2
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          nf_int2_t           i2vals(1))
      external        nf_put_att_int2

      integer         nf_get_att_int2
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          nf_int2_t           i2vals(1))
      external        nf_get_att_int2

      integer         nf_put_att_int
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          integer             ivals(1))
      external        nf_put_att_int

      integer         nf_get_att_int
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             ivals(1))
      external        nf_get_att_int

      integer         nf_put_att_real
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          real                rvals(1))
      external        nf_put_att_real

      integer         nf_get_att_real
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          real                rvals(1))
      external        nf_get_att_real

      integer         nf_put_att_double
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          double              dvals(1))
      external        nf_put_att_double

      integer         nf_get_att_double
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          double              dvals(1))
      external        nf_get_att_double

!
! general variable routines:
!

      integer         nf_def_var
!                         (integer             ncid,
!                          character(*)        name,
!                          integer             datatype,
!                          integer             ndims,
!                          integer             dimids(1),
!                          integer             varid)
      external        nf_def_var

      integer         nf_inq_var
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             datatype,
!                          integer             ndims,
!                          integer             dimids(1),
!                          integer             natts)
      external        nf_inq_var

      integer         nf_inq_varid
!                         (integer             ncid,
!                          character(*)        name,
!                          integer             varid)
      external        nf_inq_varid

      integer         nf_inq_varname
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name)
      external        nf_inq_varname

      integer         nf_inq_vartype
!                         (integer             ncid,
!                          integer             varid,
!                          integer             xtype)
      external        nf_inq_vartype

      integer         nf_inq_varndims
!                         (integer             ncid,
!                          integer             varid,
!                          integer             ndims)
      external        nf_inq_varndims

      integer         nf_inq_vardimid
!                         (integer             ncid,
!                          integer             varid,
!                          integer             dimids(1))
      external        nf_inq_vardimid

      integer         nf_inq_varnatts
!                         (integer             ncid,
!                          integer             varid,
!                          integer             natts)
      external        nf_inq_varnatts

      integer         nf_rename_var
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name)
      external        nf_rename_var

      integer         nf_copy_var
!                         (integer             ncid_in,
!                          integer             varid,
!                          integer             ncid_out)
      external        nf_copy_var

!
! entire variable put/get routines:
!

      integer         nf_put_var_text
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        text)
      external        nf_put_var_text

      integer         nf_get_var_text
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        text)
      external        nf_get_var_text

      integer         nf_put_var_int1
!                         (integer             ncid,
!                          integer             varid,
!                          nf_int1_t           i1vals(1))
      external        nf_put_var_int1

      integer         nf_get_var_int1
!                         (integer             ncid,
!                          integer             varid,
!                          nf_int1_t           i1vals(1))
      external        nf_get_var_int1

      integer         nf_put_var_int2
!                         (integer             ncid,
!                          integer             varid,
!                          nf_int2_t           i2vals(1))
      external        nf_put_var_int2

      integer         nf_get_var_int2
!                         (integer             ncid,
!                          integer             varid,
!                          nf_int2_t           i2vals(1))
      external        nf_get_var_int2

      integer         nf_put_var_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             ivals(1))
      external        nf_put_var_int

      integer         nf_get_var_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             ivals(1))
      external        nf_get_var_int

      integer         nf_put_var_real
!                         (integer             ncid,
!                          integer             varid,
!                          real                rvals(1))
      external        nf_put_var_real

      integer         nf_get_var_real
!                         (integer             ncid,
!                          integer             varid,
!                          real                rvals(1))
      external        nf_get_var_real

      integer         nf_put_var_double
!                         (integer             ncid,
!                          integer             varid,
!                          doubleprecision     dvals(1))
      external        nf_put_var_double

      integer         nf_get_var_double
!                         (integer             ncid,
!                          integer             varid,
!                          doubleprecision     dvals(1))
      external        nf_get_var_double

!
! single variable put/get routines:
!

      integer         nf_put_var1_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          character*1         text)
      external        nf_put_var1_text

      integer         nf_get_var1_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          character*1         text)
      external        nf_get_var1_text

      integer         nf_put_var1_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          nf_int1_t           i1val)
      external        nf_put_var1_int1

      integer         nf_get_var1_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          nf_int1_t           i1val)
      external        nf_get_var1_int1

      integer         nf_put_var1_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          nf_int2_t           i2val)
      external        nf_put_var1_int2

      integer         nf_get_var1_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          nf_int2_t           i2val)
      external        nf_get_var1_int2

      integer         nf_put_var1_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          integer             ival)
      external        nf_put_var1_int

      integer         nf_get_var1_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          integer             ival)
      external        nf_get_var1_int

      integer         nf_put_var1_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          real                rval)
      external        nf_put_var1_real

      integer         nf_get_var1_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          real                rval)
      external        nf_get_var1_real

      integer         nf_put_var1_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          doubleprecision     dval)
      external        nf_put_var1_double

      integer         nf_get_var1_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          doubleprecision     dval)
      external        nf_get_var1_double

!
! variable array put/get routines:
!

      integer         nf_put_vara_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          character(*)        text)
      external        nf_put_vara_text

      integer         nf_get_vara_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          character(*)        text)
      external        nf_get_vara_text

      integer         nf_put_vara_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          nf_int1_t           i1vals(1))
      external        nf_put_vara_int1

      integer         nf_get_vara_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          nf_int1_t           i1vals(1))
      external        nf_get_vara_int1

      integer         nf_put_vara_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          nf_int2_t           i2vals(1))
      external        nf_put_vara_int2

      integer         nf_get_vara_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          nf_int2_t           i2vals(1))
      external        nf_get_vara_int2

      integer         nf_put_vara_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             ivals(1))
      external        nf_put_vara_int

      integer         nf_get_vara_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             ivals(1))
      external        nf_get_vara_int

      integer         nf_put_vara_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          real                rvals(1))
      external        nf_put_vara_real

      integer         nf_get_vara_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          real                rvals(1))
      external        nf_get_vara_real

      integer         nf_put_vara_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          doubleprecision     dvals(1))
      external        nf_put_vara_double

      integer         nf_get_vara_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          doubleprecision     dvals(1))
      external        nf_get_vara_double

!
! strided variable put/get routines:
!

      integer         nf_put_vars_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          character(*)        text)
      external        nf_put_vars_text

      integer         nf_get_vars_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          character(*)        text)
      external        nf_get_vars_text

      integer         nf_put_vars_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          nf_int1_t           i1vals(1))
      external        nf_put_vars_int1

      integer         nf_get_vars_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          nf_int1_t           i1vals(1))
      external        nf_get_vars_int1

      integer         nf_put_vars_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          nf_int2_t           i2vals(1))
      external        nf_put_vars_int2

      integer         nf_get_vars_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          nf_int2_t           i2vals(1))
      external        nf_get_vars_int2

      integer         nf_put_vars_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             ivals(1))
      external        nf_put_vars_int

      integer         nf_get_vars_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             ivals(1))
      external        nf_get_vars_int

      integer         nf_put_vars_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          real                rvals(1))
      external        nf_put_vars_real

      integer         nf_get_vars_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          real                rvals(1))
      external        nf_get_vars_real

      integer         nf_put_vars_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          doubleprecision     dvals(1))
      external        nf_put_vars_double

      integer         nf_get_vars_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          doubleprecision     dvals(1))
      external        nf_get_vars_double

!
! mapped variable put/get routines:
!

      integer         nf_put_varm_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          character(*)        text)
      external        nf_put_varm_text

      integer         nf_get_varm_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          character(*)        text)
      external        nf_get_varm_text

      integer         nf_put_varm_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          nf_int1_t           i1vals(1))
      external        nf_put_varm_int1

      integer         nf_get_varm_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          nf_int1_t           i1vals(1))
      external        nf_get_varm_int1

      integer         nf_put_varm_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          nf_int2_t           i2vals(1))
      external        nf_put_varm_int2

      integer         nf_get_varm_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          nf_int2_t           i2vals(1))
      external        nf_get_varm_int2

      integer         nf_put_varm_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          integer             ivals(1))
      external        nf_put_varm_int

      integer         nf_get_varm_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          integer             ivals(1))
      external        nf_get_varm_int

      integer         nf_put_varm_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          real                rvals(1))
      external        nf_put_varm_real

      integer         nf_get_varm_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          real                rvals(1))
      external        nf_get_varm_real

      integer         nf_put_varm_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          doubleprecision     dvals(1))
      external        nf_put_varm_double

      integer         nf_get_varm_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          doubleprecision     dvals(1))
      external        nf_get_varm_double


!     NetCDF-4.
!     This is part of netCDF-4. Copyright 2006, UCAR, See COPYRIGHT
!     file for distribution information.

!     Netcdf version 4 fortran interface.

!     $Id: netcdf4.inc,v 1.28 2010/05/25 13:53:02 ed Exp $

!     New netCDF-4 types.
      integer nf_ubyte
      integer nf_ushort
      integer nf_uint
      integer nf_int64
      integer nf_uint64
      integer nf_string
      integer nf_vlen
      integer nf_opaque
      integer nf_enum
      integer nf_compound

      parameter (nf_ubyte = 7)
      parameter (nf_ushort = 8)
      parameter (nf_uint = 9)
      parameter (nf_int64 = 10)
      parameter (nf_uint64 = 11)
      parameter (nf_string = 12)
      parameter (nf_vlen = 13)
      parameter (nf_opaque = 14)
      parameter (nf_enum = 15)
      parameter (nf_compound = 16)

!     New netCDF-4 fill values.
      integer           nf_fill_ubyte
      integer           nf_fill_ushort
!      real              nf_fill_uint
!      real              nf_fill_int64
!      real              nf_fill_uint64
      parameter (nf_fill_ubyte = 255)
      parameter (nf_fill_ushort = 65535)

!     New constants.
      integer nf_format_netcdf4
      parameter (nf_format_netcdf4 = 3)

      integer nf_format_netcdf4_classic
      parameter (nf_format_netcdf4_classic = 4)

      integer nf_netcdf4
      parameter (nf_netcdf4 = 4096)

      integer nf_classic_model
      parameter (nf_classic_model = 256)

      integer nf_chunk_seq
      parameter (nf_chunk_seq = 0)
      integer nf_chunk_sub
      parameter (nf_chunk_sub = 1)
      integer nf_chunk_sizes
      parameter (nf_chunk_sizes = 2)

      integer nf_endian_native
      parameter (nf_endian_native = 0)
      integer nf_endian_little
      parameter (nf_endian_little = 1)
      integer nf_endian_big
      parameter (nf_endian_big = 2)

!     For NF_DEF_VAR_CHUNKING
      integer nf_chunked
      parameter (nf_chunked = 0)
      integer nf_contiguous
      parameter (nf_contiguous = 1)

!     For NF_DEF_VAR_FLETCHER32
      integer nf_nochecksum
      parameter (nf_nochecksum = 0)
      integer nf_fletcher32
      parameter (nf_fletcher32 = 1)

!     For NF_DEF_VAR_DEFLATE
      integer nf_noshuffle
      parameter (nf_noshuffle = 0)
      integer nf_shuffle
      parameter (nf_shuffle = 1)

!     For NF_DEF_VAR_SZIP
      integer nf_szip_ec_option_mask
      parameter (nf_szip_ec_option_mask = 4)
      integer nf_szip_nn_option_mask
      parameter (nf_szip_nn_option_mask = 32)

!     For parallel I/O.
      integer nf_mpiio      
      parameter (nf_mpiio = 8192)
      integer nf_mpiposix
      parameter (nf_mpiposix = 16384)
      integer nf_pnetcdf
      parameter (nf_pnetcdf = 32768)

!     For NF_VAR_PAR_ACCESS.
      integer nf_independent
      parameter (nf_independent = 0)
      integer nf_collective
      parameter (nf_collective = 1)

!     New error codes.
      integer nf_ehdferr        ! Error at HDF5 layer. 
      parameter (nf_ehdferr = -101)
      integer nf_ecantread      ! Can't read. 
      parameter (nf_ecantread = -102)
      integer nf_ecantwrite     ! Can't write. 
      parameter (nf_ecantwrite = -103)
      integer nf_ecantcreate    ! Can't create. 
      parameter (nf_ecantcreate = -104)
      integer nf_efilemeta      ! Problem with file metadata. 
      parameter (nf_efilemeta = -105)
      integer nf_edimmeta       ! Problem with dimension metadata. 
      parameter (nf_edimmeta = -106)
      integer nf_eattmeta       ! Problem with attribute metadata. 
      parameter (nf_eattmeta = -107)
      integer nf_evarmeta       ! Problem with variable metadata. 
      parameter (nf_evarmeta = -108)
      integer nf_enocompound    ! Not a compound type. 
      parameter (nf_enocompound = -109)
      integer nf_eattexists     ! Attribute already exists. 
      parameter (nf_eattexists = -110)
      integer nf_enotnc4        ! Attempting netcdf-4 operation on netcdf-3 file.   
      parameter (nf_enotnc4 = -111)
      integer nf_estrictnc3     ! Attempting netcdf-4 operation on strict nc3 netcdf-4 file.   
      parameter (nf_estrictnc3 = -112)
      integer nf_enotnc3        ! Attempting netcdf-3 operation on netcdf-4 file.   
      parameter (nf_enotnc3 = -113)
      integer nf_enopar         ! Parallel operation on file opened for non-parallel access.   
      parameter (nf_enopar = -114)
      integer nf_eparinit       ! Error initializing for parallel access.   
      parameter (nf_eparinit = -115)
      integer nf_ebadgrpid      ! Bad group ID.   
      parameter (nf_ebadgrpid = -116)
      integer nf_ebadtypid      ! Bad type ID.   
      parameter (nf_ebadtypid = -117)
      integer nf_etypdefined    ! Type has already been defined and may not be edited. 
      parameter (nf_etypdefined = -118)
      integer nf_ebadfield      ! Bad field ID.   
      parameter (nf_ebadfield = -119)
      integer nf_ebadclass      ! Bad class.   
      parameter (nf_ebadclass = -120)
      integer nf_emaptype       ! Mapped access for atomic types only.   
      parameter (nf_emaptype = -121)
      integer nf_elatefill      ! Attempt to define fill value when data already exists. 
      parameter (nf_elatefill = -122)
      integer nf_elatedef       ! Attempt to define var properties, like deflate, after enddef. 
      parameter (nf_elatedef = -123)
      integer nf_edimscale      ! Probem with HDF5 dimscales. 
      parameter (nf_edimscale = -124)
      integer nf_enogrp       ! No group found.
      parameter (nf_enogrp = -125)


!     New functions.

!     Parallel I/O.
      integer nf_create_par
      external nf_create_par

      integer nf_open_par
      external nf_open_par

      integer nf_var_par_access
      external nf_var_par_access

!     Functions to handle groups.
      integer nf_inq_ncid
      external nf_inq_ncid

      integer nf_inq_grps
      external nf_inq_grps

      integer nf_inq_grpname
      external nf_inq_grpname

      integer nf_inq_grpname_full
      external nf_inq_grpname_full

      integer nf_inq_grpname_len
      external nf_inq_grpname_len

      integer nf_inq_grp_parent
      external nf_inq_grp_parent

      integer nf_inq_grp_ncid
      external nf_inq_grp_ncid

      integer nf_inq_grp_full_ncid
      external nf_inq_grp_full_ncid

      integer nf_inq_varids
      external nf_inq_varids

      integer nf_inq_dimids
      external nf_inq_dimids

      integer nf_def_grp
      external nf_def_grp

!     New rename grp function

      integer nf_rename_grp
      external nf_rename_grp

!     New options for netCDF variables.
      integer nf_def_var_deflate
      external nf_def_var_deflate

      integer nf_inq_var_deflate
      external nf_inq_var_deflate

      integer nf_def_var_fletcher32
      external nf_def_var_fletcher32

      integer nf_inq_var_fletcher32
      external nf_inq_var_fletcher32

      integer nf_def_var_chunking
      external nf_def_var_chunking

      integer nf_inq_var_chunking
      external nf_inq_var_chunking

      integer nf_def_var_fill
      external nf_def_var_fill

      integer nf_inq_var_fill
      external nf_inq_var_fill

      integer nf_def_var_endian
      external nf_def_var_endian

      integer nf_inq_var_endian
      external nf_inq_var_endian

!     User defined types.
      integer nf_inq_typeids
      external nf_inq_typeids

      integer nf_inq_typeid
      external nf_inq_typeid

      integer nf_inq_type
      external nf_inq_type

      integer nf_inq_user_type
      external nf_inq_user_type

!     User defined types - compound types.
      integer nf_def_compound
      external nf_def_compound

      integer nf_insert_compound
      external nf_insert_compound

      integer nf_insert_array_compound
      external nf_insert_array_compound

      integer nf_inq_compound
      external nf_inq_compound

      integer nf_inq_compound_name
      external nf_inq_compound_name

      integer nf_inq_compound_size
      external nf_inq_compound_size

      integer nf_inq_compound_nfields
      external nf_inq_compound_nfields

      integer nf_inq_compound_field
      external nf_inq_compound_field

      integer nf_inq_compound_fieldname
      external nf_inq_compound_fieldname

      integer nf_inq_compound_fieldindex
      external nf_inq_compound_fieldindex

      integer nf_inq_compound_fieldoffset
      external nf_inq_compound_fieldoffset

      integer nf_inq_compound_fieldtype
      external nf_inq_compound_fieldtype

      integer nf_inq_compound_fieldndims
      external nf_inq_compound_fieldndims

      integer nf_inq_compound_fielddim_sizes
      external nf_inq_compound_fielddim_sizes

!     User defined types - variable length arrays.
      integer nf_def_vlen
      external nf_def_vlen

      integer nf_inq_vlen
      external nf_inq_vlen

      integer nf_free_vlen
      external nf_free_vlen

!     User defined types - enums.
      integer nf_def_enum
      external nf_def_enum

      integer nf_insert_enum
      external nf_insert_enum

      integer nf_inq_enum
      external nf_inq_enum

      integer nf_inq_enum_member
      external nf_inq_enum_member

      integer nf_inq_enum_ident
      external nf_inq_enum_ident

!     User defined types - opaque.
      integer nf_def_opaque
      external nf_def_opaque

      integer nf_inq_opaque
      external nf_inq_opaque

!     Write and read attributes of any type, including user defined
!     types.
      integer nf_put_att
      external nf_put_att
      integer nf_get_att
      external nf_get_att

!     Write and read variables of any type, including user defined
!     types.
      integer nf_put_var
      external nf_put_var
      integer nf_put_var1
      external nf_put_var1
      integer nf_put_vara
      external nf_put_vara
      integer nf_put_vars
      external nf_put_vars
      integer nf_get_var
      external nf_get_var
      integer nf_get_var1
      external nf_get_var1
      integer nf_get_vara
      external nf_get_vara
      integer nf_get_vars
      external nf_get_vars

!     64-bit int functions.
      integer nf_put_var1_int64
      external nf_put_var1_int64
      integer nf_put_vara_int64
      external nf_put_vara_int64
      integer nf_put_vars_int64
      external nf_put_vars_int64
      integer nf_put_varm_int64
      external nf_put_varm_int64
      integer nf_put_var_int64
      external nf_put_var_int64
      integer nf_get_var1_int64
      external nf_get_var1_int64
      integer nf_get_vara_int64
      external nf_get_vara_int64
      integer nf_get_vars_int64
      external nf_get_vars_int64
      integer nf_get_varm_int64
      external nf_get_varm_int64
      integer nf_get_var_int64
      external nf_get_var_int64

!     For helping F77 users with VLENs.
      integer nf_get_vlen_element
      external nf_get_vlen_element
      integer nf_put_vlen_element
      external nf_put_vlen_element

!     For dealing with file level chunk cache.
      integer nf_set_chunk_cache
      external nf_set_chunk_cache
      integer nf_get_chunk_cache
      external nf_get_chunk_cache

!     For dealing with per variable chunk cache.
      integer nf_set_var_chunk_cache
      external nf_set_var_chunk_cache
      integer nf_get_var_chunk_cache
      external nf_get_var_chunk_cache

!     NetCDF-2.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! begin netcdf 2.4 backward compatibility:
!

!      
! functions in the fortran interface
!
      integer nccre
      integer ncopn
      integer ncddef
      integer ncdid
      integer ncvdef
      integer ncvid
      integer nctlen
      integer ncsfil

      external nccre
      external ncopn
      external ncddef
      external ncdid
      external ncvdef
      external ncvid
      external nctlen
      external ncsfil


      integer ncrdwr
      integer nccreat
      integer ncexcl
      integer ncindef
      integer ncnsync
      integer nchsync
      integer ncndirty
      integer nchdirty
      integer nclink
      integer ncnowrit
      integer ncwrite
      integer ncclob
      integer ncnoclob
      integer ncglobal
      integer ncfill
      integer ncnofill
      integer maxncop
      integer maxncdim
      integer maxncatt
      integer maxncvar
      integer maxncnam
      integer maxvdims
      integer ncnoerr
      integer ncebadid
      integer ncenfile
      integer nceexist
      integer nceinval
      integer nceperm
      integer ncenotin
      integer nceindef
      integer ncecoord
      integer ncemaxds
      integer ncename
      integer ncenoatt
      integer ncemaxat
      integer ncebadty
      integer ncebadd
      integer ncests
      integer nceunlim
      integer ncemaxvs
      integer ncenotvr
      integer nceglob
      integer ncenotnc
      integer ncfoobar
      integer ncsyserr
      integer ncfatal
      integer ncverbos
      integer ncentool


!
! netcdf data types:
!
      integer ncbyte
      integer ncchar
      integer ncshort
      integer nclong
      integer ncfloat
      integer ncdouble

      parameter(ncbyte = 1)
      parameter(ncchar = 2)
      parameter(ncshort = 3)
      parameter(nclong = 4)
      parameter(ncfloat = 5)
      parameter(ncdouble = 6)

!     
!     masks for the struct nc flag field; passed in as 'mode' arg to
!     nccreate and ncopen.
!     

!     read/write, 0 => readonly 
      parameter(ncrdwr = 1)
!     in create phase, cleared by ncendef 
      parameter(nccreat = 2)
!     on create destroy existing file 
      parameter(ncexcl = 4)
!     in define mode, cleared by ncendef 
      parameter(ncindef = 8)
!     synchronise numrecs on change (x'10')
      parameter(ncnsync = 16)
!     synchronise whole header on change (x'20')
      parameter(nchsync = 32)
!     numrecs has changed (x'40')
      parameter(ncndirty = 64)  
!     header info has changed (x'80')
      parameter(nchdirty = 128)
!     prefill vars on endef and increase of record, the default behavior
      parameter(ncfill = 0)
!     do not fill vars on endef and increase of record (x'100')
      parameter(ncnofill = 256)
!     isa link (x'8000')
      parameter(nclink = 32768)

!     
!     'mode' arguments for nccreate and ncopen
!     
      parameter(ncnowrit = 0)
      parameter(ncwrite = ncrdwr)
      parameter(ncclob = nf_clobber)
      parameter(ncnoclob = nf_noclobber)

!     
!     'size' argument to ncdimdef for an unlimited dimension
!     
      integer ncunlim
      parameter(ncunlim = 0)

!     
!     attribute id to put/get a global attribute
!     
      parameter(ncglobal  = 0)

!     
!     advisory maximums:
!     
      parameter(maxncop = 64)
      parameter(maxncdim = 1024)
      parameter(maxncatt = 8192)
      parameter(maxncvar = 8192)
!     not enforced 
      parameter(maxncnam = 256)
      parameter(maxvdims = maxncdim)

!     
!     global netcdf error status variable
!     initialized in error.c
!     

!     no error 
      parameter(ncnoerr = nf_noerr)
!     not a netcdf id 
      parameter(ncebadid = nf_ebadid)
!     too many netcdfs open 
      parameter(ncenfile = -31)   ! nc_syserr
!     netcdf file exists && ncnoclob
      parameter(nceexist = nf_eexist)
!     invalid argument 
      parameter(nceinval = nf_einval)
!     write to read only 
      parameter(nceperm = nf_eperm)
!     operation not allowed in data mode 
      parameter(ncenotin = nf_enotindefine )   
!     operation not allowed in define mode 
      parameter(nceindef = nf_eindefine)   
!     coordinates out of domain 
      parameter(ncecoord = nf_einvalcoords)
!     maxncdims exceeded 
      parameter(ncemaxds = nf_emaxdims)
!     string match to name in use 
      parameter(ncename = nf_enameinuse)   
!     attribute not found 
      parameter(ncenoatt = nf_enotatt)
!     maxncattrs exceeded 
      parameter(ncemaxat = nf_emaxatts)
!     not a netcdf data type 
      parameter(ncebadty = nf_ebadtype)
!     invalid dimension id 
      parameter(ncebadd = nf_ebaddim)
!     ncunlimited in the wrong index 
      parameter(nceunlim = nf_eunlimpos)
!     maxncvars exceeded 
      parameter(ncemaxvs = nf_emaxvars)
!     variable not found 
      parameter(ncenotvr = nf_enotvar)
!     action prohibited on ncglobal varid 
      parameter(nceglob = nf_eglobal)
!     not a netcdf file 
      parameter(ncenotnc = nf_enotnc)
      parameter(ncests = nf_ests)
      parameter (ncentool = nf_emaxname) 
      parameter(ncfoobar = 32)
      parameter(ncsyserr = -31)

!     
!     global options variable. used to determine behavior of error handler.
!     initialized in lerror.c
!     
      parameter(ncfatal = 1)
      parameter(ncverbos = 2)

!
!     default fill values.  these must be the same as in the c interface.
!
      integer filbyte
      integer filchar
      integer filshort
      integer fillong
      real filfloat
      doubleprecision fildoub

      parameter (filbyte = -127)
      parameter (filchar = 0)
      parameter (filshort = -32767)
      parameter (fillong = -2147483647)
      parameter (filfloat = 9.9692099683868690e+36)
      parameter (fildoub = 9.9692099683868690e+36)

! !INPUT/OUTPUT PARAMETERS:

   type(mct_sMat)  ,intent(out)           :: sMat    ! mapping data
   type(mct_gsMap) ,intent(in) ,target    :: SgsMap  ! src gsmap
   type(mct_gSMap) ,intent(in) ,target    :: DgsMap  ! dst gsmap
   character(*)    ,intent(in)            :: newdom  ! type of sMat (src or dst)
   type(mct_Avect) ,intent(out), optional :: areasrc ! area of src grid from mapping file
   type(mct_Avect) ,intent(out), optional :: areadst ! area of dst grid from mapping file
   character(*)    ,intent(in)            :: filename! netCDF file to read
   integer(IN)     ,intent(in)            :: mytask   ! processor id
   integer(IN)     ,intent(in)            :: mpicom  ! communicator
   integer(IN)     ,intent(out), optional :: ni_i    ! number of lons on input grid   
   integer(IN)     ,intent(out), optional :: nj_i    ! number of lats on input grid   
   integer(IN)     ,intent(out), optional :: ni_o    ! number of lons on output grid   
   integer(IN)     ,intent(out), optional :: nj_o    ! number of lats on output grid   

! !EOP

   !--- local ---
   integer(IN)           :: n,m     ! generic loop indicies
   integer(IN)           :: na      ! size of source domain
   integer(IN)           :: nb      ! size of destination domain
   integer(IN)           :: ns      ! number of non-zero elements in matrix
   integer(IN)           :: ni,nj   ! number of row and col in the matrix
   integer(IN)           :: igrow   ! aVect index for matrix row
   integer(IN)           :: igcol   ! aVect index for matrix column
   integer(IN)           :: iwgt    ! aVect index for matrix element
   integer(IN)           :: iarea   ! aVect index for area
   integer(IN)           :: rsize   ! size of read buffer
   integer(IN)           :: cnt     ! local num of wgts
   integer(IN)           :: cntold  ! local num of wgts, previous read
   integer(IN)           :: start(1)! netcdf read
   integer(IN)           :: count(1)! netcdf read
   integer(IN)           :: bsize   ! buffer size
   integer(IN)           :: nread   ! number of reads 
   logical               :: mywt    ! does this weight belong on my pe

   !--- buffers for i/o ---
   real(R8)   ,allocatable :: rtemp(:) ! real temporary
   real(R8)   ,allocatable :: Sbuf(:)  ! real weights
   integer(IN),allocatable :: Rbuf(:)  ! ints rows
   integer(IN),allocatable :: Cbuf(:)  ! ints cols

   !--- variables associated with local computation of global indices
   integer(IN)             :: lsize     ! size of local seg map
   integer(IN)             :: commsize  ! size of local communicator
   integer(IN),allocatable :: lsstart(:) ! local seg map info
   integer(IN),allocatable :: lscount(:) ! local seg map info
   type(mct_gsMap),pointer :: mygsmap ! pointer to one of the gsmaps
   integer(IN)             :: l1,l2     ! generice indices for sort
   logical                 :: found     ! for sort

   !--- variable assocaited with local data buffers and reallocation
   real(R8)   ,allocatable :: Snew(:),Sold(:)  ! reals
   integer(IN),allocatable :: Rnew(:),Rold(:)  ! ints
   integer(IN),allocatable :: Cnew(:),Cold(:)  ! ints

   character,allocatable :: str(:)  ! variable length char string
   character(CL)         :: attstr  ! netCDF attribute name string
   integer(IN)           :: rcode   ! netCDF routine return code
   integer(IN)           :: fid     ! netCDF file      ID
   integer(IN)           :: vid     ! netCDF variable  ID
   integer(IN)           :: did     ! netCDF dimension ID
   !--- arbitrary size of read buffer, this is the chunk size weights reading
   integer(IN),parameter :: rbuf_size = 100000

   !--- global source and destination areas ---
   type(mct_Avect) :: areasrc0   ! area of src grid from mapping file
   type(mct_Avect) :: areadst0   ! area of src grid from mapping file

   character(*),parameter :: areaAV_field = 'aream'

   !--- formats ---
   character(*),parameter :: subName = '(shr_mct_sMatReaddnc) '
   character(*),parameter :: F00 = '("(shr_mct_sMatReaddnc) ",4a)'
   character(*),parameter :: F01 = '("(shr_mct_sMatReaddnc) ",2(a,i10))'

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

 call shr_mpi_commsize(mpicom,commsize)
 if (mytask == 0) then
   if (s_loglev > 0) write(s_logunit,F00) "reading mapping matrix data decomposed..."

   !----------------------------------------------------------------------------
   ! open & read the file
   !----------------------------------------------------------------------------
   if (s_loglev > 0) write(s_logunit,F00) "* file name                  : ",trim(fileName)
   rcode = nf_open(filename,NF_NOWRITE,fid)
   if (rcode /= NF_NOERR) then 
      print *,'Failed to open file ',trim(filename)
      call shr_sys_abort(trim(subName)//nf_strerror(rcode))
   end if
   

   !--- get matrix dimensions ----------
   rcode = nf_inq_dimid (fid, 'n_s', did)  ! size of sparse matrix
   rcode = nf_inq_dimlen(fid, did  , ns)
   rcode = nf_inq_dimid (fid, 'n_a', did)  ! size of  input vector
   rcode = nf_inq_dimlen(fid, did  , na)
   rcode = nf_inq_dimid (fid, 'n_b', did)  ! size of output vector
   rcode = nf_inq_dimlen(fid, did  , nb)
   
   if (present(ni_i) .and. present(nj_i) .and. present(ni_o) .and. present(nj_o)) then
      rcode = nf_inq_dimid (fid, 'ni_a', did)  ! number of lons in input grid
      rcode = nf_inq_dimlen(fid, did  , ni_i)
      rcode = nf_inq_dimid (fid, 'nj_a', did)  ! number of lats in input grid
      rcode = nf_inq_dimlen(fid, did  , nj_i)
      rcode = nf_inq_dimid (fid, 'ni_b', did)  ! number of lons in output grid
      rcode = nf_inq_dimlen(fid, did  , ni_o)
      rcode = nf_inq_dimid (fid, 'nj_b', did)  ! number of lats in output grid
      rcode = nf_inq_dimlen(fid, did  , nj_o)
   end if

   if (s_loglev > 0) write(s_logunit,F01) "* matrix dims src x dst      : ",na,' x',nb
   if (s_loglev > 0) write(s_logunit,F01) "* number of non-zero elements: ",ns

 endif
 
   !--- read and load area_a ---
   if (present(areasrc)) then
   if (mytask == 0) then
      call mct_aVect_init(areasrc0,' ',areaAV_field,na)
      rcode = nf_inq_varid     (fid,'area_a',vid)
      if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)
      rcode = nf_get_var_double(fid, vid, areasrc0%rAttr)
      if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)
   endif
   call mct_aVect_scatter(areasrc0, areasrc, SgsMap, 0, mpicom, rcode)
   if (rcode /= 0) call mct_die("shr_mct_sMatReaddnc","Error on scatter of areasrc0")
   if (mytask == 0) then
!      if (present(dbug)) then
!         if (dbug > 2) then
!            write(6,*) subName,'Size of src ',mct_aVect_lSize(areasrc0)
!            write(6,*) subName,'min/max src ',minval(areasrc0%rAttr(1,:)),maxval(areasrc0%rAttr(1,:))
!         endif
!      end if
      call mct_aVect_clean(areasrc0)
   end if
   end if

   !--- read and load area_b ---
   if (present(areadst)) then
   if (mytask == 0) then
      call mct_aVect_init(areadst0,' ',areaAV_field,nb)
      rcode = nf_inq_varid     (fid,'area_b',vid)
      if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)
      rcode = nf_get_var_double(fid, vid, areadst0%rAttr)
      if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)
   endif
   call mct_aVect_scatter(areadst0, areadst, DgsMap, 0, mpicom, rcode)
   if (rcode /= 0) call mct_die("shr_mct_sMatReaddnc","Error on scatter of areadst0")
   if (mytask == 0) then
!      if (present(dbug)) then
!         if (dbug > 2) then
!            write(6,*) subName,'Size of dst ',mct_aVect_lSize(areadst0)
!            write(6,*) subName,'min/max dst ',minval(areadst0%rAttr(1,:)),maxval(areadst0%rAttr(1,:))
!         endif
!      end if
      call mct_aVect_clean(areadst0)
   endif
   endif

   if (present(ni_i) .and. present(nj_i) .and. present(ni_o) .and. present(nj_o)) then
      call shr_mpi_bcast(ni_i,mpicom,subName//" MPI in ni_i bcast")
      call shr_mpi_bcast(nj_i,mpicom,subName//" MPI in nj_i bcast")
      call shr_mpi_bcast(ni_o,mpicom,subName//" MPI in ni_o bcast")
      call shr_mpi_bcast(nj_o,mpicom,subName//" MPI in nj_o bcast")
   end if

   call shr_mpi_bcast(ns,mpicom,subName//" MPI in ns bcast")
   call shr_mpi_bcast(na,mpicom,subName//" MPI in na bcast")
   call shr_mpi_bcast(nb,mpicom,subName//" MPI in nb bcast")

   !--- setup local seg map, sorted
   if (newdom == 'src') then
      mygsmap => DgsMap
   elseif (newdom == 'dst') then
      mygsmap => SgsMap
   else
      write(s_logunit,F00) 'ERROR: invalid newdom value = ',newdom
      call shr_sys_abort(trim(subName)//" invalid newdom value")
   endif
   lsize = 0
   do n = 1,size(mygsmap%start)
      if (mygsmap%pe_loc(n) == mytask) then
         lsize=lsize+1
      endif
   enddo
   allocate(lsstart(lsize),lscount(lsize),stat=rcode)
   if (rcode /= 0) call mct_perr_die(subName,':: allocate Lsstart',rcode)

   lsize = 0
   do n = 1,size(mygsmap%start)
      if (mygsmap%pe_loc(n) == mytask) then  ! on my pe
         lsize=lsize+1
         found = .false.
         l1 = 1
         do while (.not.found .and. l1 < lsize)         ! bubble sort copy
            if (mygsmap%start(n) < lsstart(l1)) then
               do l2 = lsize, l1+1, -1
                  lsstart(l2) = lsstart(l2-1)
                  lscount(l2) = lscount(l2-1)
               enddo
               found = .true.
            else
               l1 = l1 + 1
            endif
         enddo
         lsstart(l1) = mygsmap%start(n)
         lscount(l1) = mygsmap%length(n)
      endif
   enddo
   do n = 1,lsize-1
      if (lsstart(n) > lsstart(n+1)) then
         write(s_logunit,F00) ' ERROR: lsstart not properly sorted'
         call shr_sys_abort()
      endif
   enddo

   rsize = min(rbuf_size,ns)                     ! size of i/o chunks
   bsize = ((ns/commsize) + 1 ) * 1.2   ! local temporary buffer size
   if (ns == 0) then
      nread = 0
   else
      nread = (ns-1)/rsize + 1                      ! num of reads to do
   endif

   allocate(Sbuf(rsize),Rbuf(rsize),Cbuf(rsize),stat=rcode)
   if (rcode /= 0) call mct_perr_die(subName,':: allocate Sbuf',rcode)
   allocate(Snew(bsize),Cnew(bsize),Rnew(bsize),stat=rcode)
   if (rcode /= 0) call mct_perr_die(subName,':: allocate Snew1',rcode)

   cnt = 0
   do n = 1,nread
      start(1) = (n-1)*rsize + 1
      count(1) = min(rsize,ns-start(1)+1)

      !--- read data on root pe
      if (mytask== 0) then
         rcode = nf_inq_varid      (fid,'S'  ,vid)
         rcode = nf_get_vara_double(fid,vid,start,count,Sbuf)
         if (rcode /= NF_NOERR .and. s_loglev > 0) write(s_logunit,F00) nf_strerror(rcode)

         rcode = nf_inq_varid      (fid,'row',vid)
         rcode = nf_get_vara_int   (fid,vid,start,count,Rbuf)
         if (rcode /= NF_NOERR .and. s_loglev > 0) write(s_logunit,F00) nf_strerror(rcode)

         rcode = nf_inq_varid      (fid,'col',vid)
         rcode = nf_get_vara_int   (fid,vid,start,count,Cbuf)
         if (rcode /= NF_NOERR .and. s_loglev > 0) write(s_logunit,F00) nf_strerror(rcode)
      endif

      !--- send S, row, col to all pes
      call shr_mpi_bcast(Sbuf,mpicom,subName//" MPI in Sbuf bcast")
      call shr_mpi_bcast(Rbuf,mpicom,subName//" MPI in Rbuf bcast")
      call shr_mpi_bcast(Cbuf,mpicom,subName//" MPI in Cbuf bcast")

      !--- now each pe keeps what it should
      do m = 1,count(1)
         !--- should this weight be on my pe
         if (newdom == 'src') then
            mywt = mct_myindex(Rbuf(m),lsstart,lscount)
         elseif (newdom == 'dst') then
            mywt = mct_myindex(Cbuf(m),lsstart,lscount)
         endif

         if (mywt) then
            cntold = cnt
            cnt = cnt + 1

            !--- new arrays need to be bigger
            if (cnt > bsize) then
               !--- allocate old arrays and copy new into old
               allocate(Sold(cntold),Rold(cntold),Cold(cntold),stat=rcode)
               if (rcode /= 0) call mct_perr_die(subName,':: allocate old',rcode)
               Sold(1:cntold) = Snew(1:cntold)
               Rold(1:cntold) = Rnew(1:cntold)
               Cold(1:cntold) = Cnew(1:cntold)

               !--- reallocate new to bigger size, increase buffer by 50% (arbitrary)
               deallocate(Snew,Rnew,Cnew,stat=rcode)
               if (rcode /= 0) call mct_perr_die(subName,':: allocate new',rcode)
               bsize = 1.5 * bsize
               if (s_loglev > 1) write(s_logunit,F01) ' reallocate bsize to ',bsize
               allocate(Snew(bsize),Rnew(bsize),Cnew(bsize),stat=rcode)
               if (rcode /= 0) call mct_perr_die(subName,':: allocate old',rcode)

               !--- copy data back into new
               Snew(1:cntold) = Sold(1:cntold)
               Rnew(1:cntold) = Rold(1:cntold)
               Cnew(1:cntold) = Cold(1:cntold)
               deallocate(Sold,Rold,Cold,stat=rcode)
               if (rcode /= 0) call mct_perr_die(subName,':: deallocate old',rcode)
            endif

            Snew(cnt) = Sbuf(m)
            Rnew(cnt) = Rbuf(m)
            Cnew(cnt) = Cbuf(m)
         endif
      enddo  ! count
   enddo   ! nread

   deallocate(Sbuf,Rbuf,Cbuf, stat=rcode)
   if (rcode /= 0) call mct_perr_die(subName,':: deallocate Sbuf',rcode)

   !----------------------------------------------------------------------------
   ! init the mct sMat data type
   !----------------------------------------------------------------------------
   ! mct_sMat_init must be given the number of rows and columns that
   ! would be in the full matrix.  Nrows= size of output vector=nb.
   ! Ncols = size of input vector = na.
   call mct_sMat_init(sMat, nb, na, cnt)

   igrow = mct_sMat_indexIA(sMat,'grow')
   igcol = mct_sMat_indexIA(sMat,'gcol')
   iwgt  = mct_sMat_indexRA(sMat,'weight')

   if (cnt /= 0) then
      sMat%data%rAttr(iwgt ,1:cnt) = Snew(1:cnt)
      sMat%data%iAttr(igrow,1:cnt) = Rnew(1:cnt)
      sMat%data%iAttr(igcol,1:cnt) = Cnew(1:cnt)
   endif
   deallocate(Snew,Rnew,Cnew, stat=rcode)
   deallocate(lsstart,lscount,stat=rcode)
   if (rcode /= 0) call mct_perr_die(subName,':: deallocate new',rcode)

   if (mytask == 0) then
      rcode = nf_close(fid)
      if (s_loglev > 0) write(s_logunit,F00) "... done reading file"
      call shr_sys_flush(s_logunit)
   endif

end subroutine shr_mct_sMatReaddnc

!BOP ===========================================================================
!
! !IROUTINE:  shr_mct_sMatWritednc - Do a distributed write of a NetCDF SCRIP file
!                                 based on a distributed SparseMatrix
!
! !DESCRIPTION: 
!     Write out mapping matrix data from a SCRIP netCDF data file using
!     a low memory method.
!
! !SEE ALSO: 
!    mct/m_SparseMatrix.F90 (MCT source code)
!
! !REVISION HISTORY: 
!    2009-Dec-15 - T. Craig -- first version
! 
! !INTERFACE:  -----------------------------------------------------------------

! Subprogram not used subroutine shr_mct_sMatWritednc(sMat,iosystem, io_type, fileName,compid, mpicom)
! Subprogram not used 
! Subprogram not used ! !USES:
! Subprogram not used   use pio, only : iosystem_desc_t
! Subprogram not used    use shr_pcdf_mod, only : shr_pcdf_readwrite
! Subprogram not used    implicit none
! Subprogram not used !     NetCDF-3.
! Subprogram not used !
! Subprogram not used ! netcdf version 3 fortran interface:
! Subprogram not used !
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! external netcdf data types:
! Subprogram not used !
! Subprogram not used       integer nf_byte
! Subprogram not used       integer nf_int1
! Subprogram not used       integer nf_char
! Subprogram not used       integer nf_short
! Subprogram not used       integer nf_int2
! Subprogram not used       integer nf_int
! Subprogram not used       integer nf_float
! Subprogram not used       integer nf_real
! Subprogram not used       integer nf_double
! Subprogram not used 
! Subprogram not used       parameter (nf_byte = 1)
! Subprogram not used       parameter (nf_int1 = nf_byte)
! Subprogram not used       parameter (nf_char = 2)
! Subprogram not used       parameter (nf_short = 3)
! Subprogram not used       parameter (nf_int2 = nf_short)
! Subprogram not used       parameter (nf_int = 4)
! Subprogram not used       parameter (nf_float = 5)
! Subprogram not used       parameter (nf_real = nf_float)
! Subprogram not used       parameter (nf_double = 6)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! default fill values:
! Subprogram not used !
! Subprogram not used       integer           nf_fill_byte
! Subprogram not used       integer           nf_fill_int1
! Subprogram not used       integer           nf_fill_char
! Subprogram not used       integer           nf_fill_short
! Subprogram not used       integer           nf_fill_int2
! Subprogram not used       integer           nf_fill_int
! Subprogram not used       real              nf_fill_float
! Subprogram not used       real              nf_fill_real
! Subprogram not used       doubleprecision   nf_fill_double
! Subprogram not used 
! Subprogram not used       parameter (nf_fill_byte = -127)
! Subprogram not used       parameter (nf_fill_int1 = nf_fill_byte)
! Subprogram not used       parameter (nf_fill_char = 0)
! Subprogram not used       parameter (nf_fill_short = -32767)
! Subprogram not used       parameter (nf_fill_int2 = nf_fill_short)
! Subprogram not used       parameter (nf_fill_int = -2147483647)
! Subprogram not used       parameter (nf_fill_float = 9.9692099683868690e+36)
! Subprogram not used       parameter (nf_fill_real = nf_fill_float)
! Subprogram not used       parameter (nf_fill_double = 9.9692099683868690d+36)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! mode flags for opening and creating a netcdf dataset:
! Subprogram not used !
! Subprogram not used       integer nf_nowrite
! Subprogram not used       integer nf_write
! Subprogram not used       integer nf_clobber
! Subprogram not used       integer nf_noclobber
! Subprogram not used       integer nf_fill
! Subprogram not used       integer nf_nofill
! Subprogram not used       integer nf_lock
! Subprogram not used       integer nf_share
! Subprogram not used       integer nf_64bit_offset
! Subprogram not used       integer nf_sizehint_default
! Subprogram not used       integer nf_align_chunk
! Subprogram not used       integer nf_format_classic
! Subprogram not used       integer nf_format_64bit
! Subprogram not used       integer nf_diskless
! Subprogram not used       integer nf_mmap
! Subprogram not used 
! Subprogram not used       parameter (nf_nowrite = 0)
! Subprogram not used       parameter (nf_write = 1)
! Subprogram not used       parameter (nf_clobber = 0)
! Subprogram not used       parameter (nf_noclobber = 4)
! Subprogram not used       parameter (nf_fill = 0)
! Subprogram not used       parameter (nf_nofill = 256)
! Subprogram not used       parameter (nf_lock = 1024)
! Subprogram not used       parameter (nf_share = 2048)
! Subprogram not used       parameter (nf_64bit_offset = 512)
! Subprogram not used       parameter (nf_sizehint_default = 0)
! Subprogram not used       parameter (nf_align_chunk = -1)
! Subprogram not used       parameter (nf_format_classic = 1)
! Subprogram not used       parameter (nf_format_64bit = 2)
! Subprogram not used       parameter (nf_diskless = 8)
! Subprogram not used       parameter (nf_mmap = 16)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! size argument for defining an unlimited dimension:
! Subprogram not used !
! Subprogram not used       integer nf_unlimited
! Subprogram not used       parameter (nf_unlimited = 0)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! global attribute id:
! Subprogram not used !
! Subprogram not used       integer nf_global
! Subprogram not used       parameter (nf_global = 0)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! implementation limits:
! Subprogram not used !
! Subprogram not used       integer nf_max_dims
! Subprogram not used       integer nf_max_attrs
! Subprogram not used       integer nf_max_vars
! Subprogram not used       integer nf_max_name
! Subprogram not used       integer nf_max_var_dims
! Subprogram not used 
! Subprogram not used       parameter (nf_max_dims = 1024)
! Subprogram not used       parameter (nf_max_attrs = 8192)
! Subprogram not used       parameter (nf_max_vars = 8192)
! Subprogram not used       parameter (nf_max_name = 256)
! Subprogram not used       parameter (nf_max_var_dims = nf_max_dims)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! error codes:
! Subprogram not used !
! Subprogram not used       integer nf_noerr
! Subprogram not used       integer nf_ebadid
! Subprogram not used       integer nf_eexist
! Subprogram not used       integer nf_einval
! Subprogram not used       integer nf_eperm
! Subprogram not used       integer nf_enotindefine
! Subprogram not used       integer nf_eindefine
! Subprogram not used       integer nf_einvalcoords
! Subprogram not used       integer nf_emaxdims
! Subprogram not used       integer nf_enameinuse
! Subprogram not used       integer nf_enotatt
! Subprogram not used       integer nf_emaxatts
! Subprogram not used       integer nf_ebadtype
! Subprogram not used       integer nf_ebaddim
! Subprogram not used       integer nf_eunlimpos
! Subprogram not used       integer nf_emaxvars
! Subprogram not used       integer nf_enotvar
! Subprogram not used       integer nf_eglobal
! Subprogram not used       integer nf_enotnc
! Subprogram not used       integer nf_ests
! Subprogram not used       integer nf_emaxname
! Subprogram not used       integer nf_eunlimit
! Subprogram not used       integer nf_enorecvars
! Subprogram not used       integer nf_echar
! Subprogram not used       integer nf_eedge
! Subprogram not used       integer nf_estride
! Subprogram not used       integer nf_ebadname
! Subprogram not used       integer nf_erange
! Subprogram not used       integer nf_enomem
! Subprogram not used       integer nf_evarsize
! Subprogram not used       integer nf_edimsize
! Subprogram not used       integer nf_etrunc
! Subprogram not used 
! Subprogram not used       parameter (nf_noerr = 0)
! Subprogram not used       parameter (nf_ebadid = -33)
! Subprogram not used       parameter (nf_eexist = -35)
! Subprogram not used       parameter (nf_einval = -36)
! Subprogram not used       parameter (nf_eperm = -37)
! Subprogram not used       parameter (nf_enotindefine = -38)
! Subprogram not used       parameter (nf_eindefine = -39)
! Subprogram not used       parameter (nf_einvalcoords = -40)
! Subprogram not used       parameter (nf_emaxdims = -41)
! Subprogram not used       parameter (nf_enameinuse = -42)
! Subprogram not used       parameter (nf_enotatt = -43)
! Subprogram not used       parameter (nf_emaxatts = -44)
! Subprogram not used       parameter (nf_ebadtype = -45)
! Subprogram not used       parameter (nf_ebaddim = -46)
! Subprogram not used       parameter (nf_eunlimpos = -47)
! Subprogram not used       parameter (nf_emaxvars = -48)
! Subprogram not used       parameter (nf_enotvar = -49)
! Subprogram not used       parameter (nf_eglobal = -50)
! Subprogram not used       parameter (nf_enotnc = -51)
! Subprogram not used       parameter (nf_ests = -52)
! Subprogram not used       parameter (nf_emaxname = -53)
! Subprogram not used       parameter (nf_eunlimit = -54)
! Subprogram not used       parameter (nf_enorecvars = -55)
! Subprogram not used       parameter (nf_echar = -56)
! Subprogram not used       parameter (nf_eedge = -57)
! Subprogram not used       parameter (nf_estride = -58)
! Subprogram not used       parameter (nf_ebadname = -59)
! Subprogram not used       parameter (nf_erange = -60)
! Subprogram not used       parameter (nf_enomem = -61)
! Subprogram not used       parameter (nf_evarsize = -62)
! Subprogram not used       parameter (nf_edimsize = -63)
! Subprogram not used       parameter (nf_etrunc = -64)
! Subprogram not used !
! Subprogram not used ! error handling modes:
! Subprogram not used !
! Subprogram not used       integer  nf_fatal
! Subprogram not used       integer nf_verbose
! Subprogram not used 
! Subprogram not used       parameter (nf_fatal = 1)
! Subprogram not used       parameter (nf_verbose = 2)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! miscellaneous routines:
! Subprogram not used !
! Subprogram not used       character*80   nf_inq_libvers
! Subprogram not used       external       nf_inq_libvers
! Subprogram not used 
! Subprogram not used       character*80   nf_strerror
! Subprogram not used !                         (integer             ncerr)
! Subprogram not used       external       nf_strerror
! Subprogram not used 
! Subprogram not used       logical        nf_issyserr
! Subprogram not used !                         (integer             ncerr)
! Subprogram not used       external       nf_issyserr
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! control routines:
! Subprogram not used !
! Subprogram not used       integer         nf_inq_base_pe
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             pe)
! Subprogram not used       external        nf_inq_base_pe
! Subprogram not used 
! Subprogram not used       integer         nf_set_base_pe
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             pe)
! Subprogram not used       external        nf_set_base_pe
! Subprogram not used 
! Subprogram not used       integer         nf_create
! Subprogram not used !                         (character*(*)       path,
! Subprogram not used !                          integer             cmode,
! Subprogram not used !                          integer             ncid)
! Subprogram not used       external        nf_create
! Subprogram not used 
! Subprogram not used       integer         nf__create
! Subprogram not used !                         (character*(*)       path,
! Subprogram not used !                          integer             cmode,
! Subprogram not used !                          integer             initialsz,
! Subprogram not used !                          integer             chunksizehint,
! Subprogram not used !                          integer             ncid)
! Subprogram not used       external        nf__create
! Subprogram not used 
! Subprogram not used       integer         nf__create_mp
! Subprogram not used !                         (character*(*)       path,
! Subprogram not used !                          integer             cmode,
! Subprogram not used !                          integer             initialsz,
! Subprogram not used !                          integer             basepe,
! Subprogram not used !                          integer             chunksizehint,
! Subprogram not used !                          integer             ncid)
! Subprogram not used       external        nf__create_mp
! Subprogram not used 
! Subprogram not used       integer         nf_open
! Subprogram not used !                         (character*(*)       path,
! Subprogram not used !                          integer             mode,
! Subprogram not used !                          integer             ncid)
! Subprogram not used       external        nf_open
! Subprogram not used 
! Subprogram not used       integer         nf__open
! Subprogram not used !                         (character*(*)       path,
! Subprogram not used !                          integer             mode,
! Subprogram not used !                          integer             chunksizehint,
! Subprogram not used !                          integer             ncid)
! Subprogram not used       external        nf__open
! Subprogram not used 
! Subprogram not used       integer         nf__open_mp
! Subprogram not used !                         (character*(*)       path,
! Subprogram not used !                          integer             mode,
! Subprogram not used !                          integer             basepe,
! Subprogram not used !                          integer             chunksizehint,
! Subprogram not used !                          integer             ncid)
! Subprogram not used       external        nf__open_mp
! Subprogram not used 
! Subprogram not used       integer         nf_set_fill
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             fillmode,
! Subprogram not used !                          integer             old_mode)
! Subprogram not used       external        nf_set_fill
! Subprogram not used 
! Subprogram not used       integer         nf_set_default_format
! Subprogram not used !                          (integer             format,
! Subprogram not used !                          integer             old_format)
! Subprogram not used       external        nf_set_default_format
! Subprogram not used 
! Subprogram not used       integer         nf_redef
! Subprogram not used !                         (integer             ncid)
! Subprogram not used       external        nf_redef
! Subprogram not used 
! Subprogram not used       integer         nf_enddef
! Subprogram not used !                         (integer             ncid)
! Subprogram not used       external        nf_enddef
! Subprogram not used 
! Subprogram not used       integer         nf__enddef
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             h_minfree,
! Subprogram not used !                          integer             v_align,
! Subprogram not used !                          integer             v_minfree,
! Subprogram not used !                          integer             r_align)
! Subprogram not used       external        nf__enddef
! Subprogram not used 
! Subprogram not used       integer         nf_sync
! Subprogram not used !                         (integer             ncid)
! Subprogram not used       external        nf_sync
! Subprogram not used 
! Subprogram not used       integer         nf_abort
! Subprogram not used !                         (integer             ncid)
! Subprogram not used       external        nf_abort
! Subprogram not used 
! Subprogram not used       integer         nf_close
! Subprogram not used !                         (integer             ncid)
! Subprogram not used       external        nf_close
! Subprogram not used 
! Subprogram not used       integer         nf_delete
! Subprogram not used !                         (character*(*)       ncid)
! Subprogram not used       external        nf_delete
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! general inquiry routines:
! Subprogram not used !
! Subprogram not used 
! Subprogram not used       integer         nf_inq
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             ndims,
! Subprogram not used !                          integer             nvars,
! Subprogram not used !                          integer             ngatts,
! Subprogram not used !                          integer             unlimdimid)
! Subprogram not used       external        nf_inq
! Subprogram not used 
! Subprogram not used ! new inquire path
! Subprogram not used 
! Subprogram not used       integer nf_inq_path
! Subprogram not used       external nf_inq_path
! Subprogram not used 
! Subprogram not used       integer         nf_inq_ndims
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             ndims)
! Subprogram not used       external        nf_inq_ndims
! Subprogram not used 
! Subprogram not used       integer         nf_inq_nvars
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             nvars)
! Subprogram not used       external        nf_inq_nvars
! Subprogram not used 
! Subprogram not used       integer         nf_inq_natts
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             ngatts)
! Subprogram not used       external        nf_inq_natts
! Subprogram not used 
! Subprogram not used       integer         nf_inq_unlimdim
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             unlimdimid)
! Subprogram not used       external        nf_inq_unlimdim
! Subprogram not used 
! Subprogram not used       integer         nf_inq_format
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             format)
! Subprogram not used       external        nf_inq_format
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! dimension routines:
! Subprogram not used !
! Subprogram not used 
! Subprogram not used       integer         nf_def_dim
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             len,
! Subprogram not used !                          integer             dimid)
! Subprogram not used       external        nf_def_dim
! Subprogram not used 
! Subprogram not used       integer         nf_inq_dimid
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             dimid)
! Subprogram not used       external        nf_inq_dimid
! Subprogram not used 
! Subprogram not used       integer         nf_inq_dim
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             dimid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             len)
! Subprogram not used       external        nf_inq_dim
! Subprogram not used 
! Subprogram not used       integer         nf_inq_dimname
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             dimid,
! Subprogram not used !                          character(*)        name)
! Subprogram not used       external        nf_inq_dimname
! Subprogram not used 
! Subprogram not used       integer         nf_inq_dimlen
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             dimid,
! Subprogram not used !                          integer             len)
! Subprogram not used       external        nf_inq_dimlen
! Subprogram not used 
! Subprogram not used       integer         nf_rename_dim
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             dimid,
! Subprogram not used !                          character(*)        name)
! Subprogram not used       external        nf_rename_dim
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! general attribute routines:
! Subprogram not used !
! Subprogram not used 
! Subprogram not used       integer         nf_inq_att
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             xtype,
! Subprogram not used !                          integer             len)
! Subprogram not used       external        nf_inq_att
! Subprogram not used 
! Subprogram not used       integer         nf_inq_attid
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             attnum)
! Subprogram not used       external        nf_inq_attid
! Subprogram not used 
! Subprogram not used       integer         nf_inq_atttype
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             xtype)
! Subprogram not used       external        nf_inq_atttype
! Subprogram not used 
! Subprogram not used       integer         nf_inq_attlen
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             len)
! Subprogram not used       external        nf_inq_attlen
! Subprogram not used 
! Subprogram not used       integer         nf_inq_attname
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             attnum,
! Subprogram not used !                          character(*)        name)
! Subprogram not used       external        nf_inq_attname
! Subprogram not used 
! Subprogram not used       integer         nf_copy_att
! Subprogram not used !                         (integer             ncid_in,
! Subprogram not used !                          integer             varid_in,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             ncid_out,
! Subprogram not used !                          integer             varid_out)
! Subprogram not used       external        nf_copy_att
! Subprogram not used 
! Subprogram not used       integer         nf_rename_att
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        curname,
! Subprogram not used !                          character(*)        newname)
! Subprogram not used       external        nf_rename_att
! Subprogram not used 
! Subprogram not used       integer         nf_del_att
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name)
! Subprogram not used       external        nf_del_att
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! attribute put/get routines:
! Subprogram not used !
! Subprogram not used 
! Subprogram not used       integer         nf_put_att_text
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             len,
! Subprogram not used !                          character(*)        text)
! Subprogram not used       external        nf_put_att_text
! Subprogram not used 
! Subprogram not used       integer         nf_get_att_text
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          character(*)        text)
! Subprogram not used       external        nf_get_att_text
! Subprogram not used 
! Subprogram not used       integer         nf_put_att_int1
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             xtype,
! Subprogram not used !                          integer             len,
! Subprogram not used !                          nf_int1_t           i1vals(1))
! Subprogram not used       external        nf_put_att_int1
! Subprogram not used 
! Subprogram not used       integer         nf_get_att_int1
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          nf_int1_t           i1vals(1))
! Subprogram not used       external        nf_get_att_int1
! Subprogram not used 
! Subprogram not used       integer         nf_put_att_int2
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             xtype,
! Subprogram not used !                          integer             len,
! Subprogram not used !                          nf_int2_t           i2vals(1))
! Subprogram not used       external        nf_put_att_int2
! Subprogram not used 
! Subprogram not used       integer         nf_get_att_int2
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          nf_int2_t           i2vals(1))
! Subprogram not used       external        nf_get_att_int2
! Subprogram not used 
! Subprogram not used       integer         nf_put_att_int
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             xtype,
! Subprogram not used !                          integer             len,
! Subprogram not used !                          integer             ivals(1))
! Subprogram not used       external        nf_put_att_int
! Subprogram not used 
! Subprogram not used       integer         nf_get_att_int
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             ivals(1))
! Subprogram not used       external        nf_get_att_int
! Subprogram not used 
! Subprogram not used       integer         nf_put_att_real
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             xtype,
! Subprogram not used !                          integer             len,
! Subprogram not used !                          real                rvals(1))
! Subprogram not used       external        nf_put_att_real
! Subprogram not used 
! Subprogram not used       integer         nf_get_att_real
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          real                rvals(1))
! Subprogram not used       external        nf_get_att_real
! Subprogram not used 
! Subprogram not used       integer         nf_put_att_double
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             xtype,
! Subprogram not used !                          integer             len,
! Subprogram not used !                          double              dvals(1))
! Subprogram not used       external        nf_put_att_double
! Subprogram not used 
! Subprogram not used       integer         nf_get_att_double
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          double              dvals(1))
! Subprogram not used       external        nf_get_att_double
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! general variable routines:
! Subprogram not used !
! Subprogram not used 
! Subprogram not used       integer         nf_def_var
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             datatype,
! Subprogram not used !                          integer             ndims,
! Subprogram not used !                          integer             dimids(1),
! Subprogram not used !                          integer             varid)
! Subprogram not used       external        nf_def_var
! Subprogram not used 
! Subprogram not used       integer         nf_inq_var
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             datatype,
! Subprogram not used !                          integer             ndims,
! Subprogram not used !                          integer             dimids(1),
! Subprogram not used !                          integer             natts)
! Subprogram not used       external        nf_inq_var
! Subprogram not used 
! Subprogram not used       integer         nf_inq_varid
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          character(*)        name,
! Subprogram not used !                          integer             varid)
! Subprogram not used       external        nf_inq_varid
! Subprogram not used 
! Subprogram not used       integer         nf_inq_varname
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name)
! Subprogram not used       external        nf_inq_varname
! Subprogram not used 
! Subprogram not used       integer         nf_inq_vartype
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             xtype)
! Subprogram not used       external        nf_inq_vartype
! Subprogram not used 
! Subprogram not used       integer         nf_inq_varndims
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             ndims)
! Subprogram not used       external        nf_inq_varndims
! Subprogram not used 
! Subprogram not used       integer         nf_inq_vardimid
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             dimids(1))
! Subprogram not used       external        nf_inq_vardimid
! Subprogram not used 
! Subprogram not used       integer         nf_inq_varnatts
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             natts)
! Subprogram not used       external        nf_inq_varnatts
! Subprogram not used 
! Subprogram not used       integer         nf_rename_var
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        name)
! Subprogram not used       external        nf_rename_var
! Subprogram not used 
! Subprogram not used       integer         nf_copy_var
! Subprogram not used !                         (integer             ncid_in,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             ncid_out)
! Subprogram not used       external        nf_copy_var
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! entire variable put/get routines:
! Subprogram not used !
! Subprogram not used 
! Subprogram not used       integer         nf_put_var_text
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        text)
! Subprogram not used       external        nf_put_var_text
! Subprogram not used 
! Subprogram not used       integer         nf_get_var_text
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          character(*)        text)
! Subprogram not used       external        nf_get_var_text
! Subprogram not used 
! Subprogram not used       integer         nf_put_var_int1
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          nf_int1_t           i1vals(1))
! Subprogram not used       external        nf_put_var_int1
! Subprogram not used 
! Subprogram not used       integer         nf_get_var_int1
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          nf_int1_t           i1vals(1))
! Subprogram not used       external        nf_get_var_int1
! Subprogram not used 
! Subprogram not used       integer         nf_put_var_int2
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          nf_int2_t           i2vals(1))
! Subprogram not used       external        nf_put_var_int2
! Subprogram not used 
! Subprogram not used       integer         nf_get_var_int2
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          nf_int2_t           i2vals(1))
! Subprogram not used       external        nf_get_var_int2
! Subprogram not used 
! Subprogram not used       integer         nf_put_var_int
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             ivals(1))
! Subprogram not used       external        nf_put_var_int
! Subprogram not used 
! Subprogram not used       integer         nf_get_var_int
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             ivals(1))
! Subprogram not used       external        nf_get_var_int
! Subprogram not used 
! Subprogram not used       integer         nf_put_var_real
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          real                rvals(1))
! Subprogram not used       external        nf_put_var_real
! Subprogram not used 
! Subprogram not used       integer         nf_get_var_real
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          real                rvals(1))
! Subprogram not used       external        nf_get_var_real
! Subprogram not used 
! Subprogram not used       integer         nf_put_var_double
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          doubleprecision     dvals(1))
! Subprogram not used       external        nf_put_var_double
! Subprogram not used 
! Subprogram not used       integer         nf_get_var_double
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          doubleprecision     dvals(1))
! Subprogram not used       external        nf_get_var_double
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! single variable put/get routines:
! Subprogram not used !
! Subprogram not used 
! Subprogram not used       integer         nf_put_var1_text
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             index(1),
! Subprogram not used !                          character*1         text)
! Subprogram not used       external        nf_put_var1_text
! Subprogram not used 
! Subprogram not used       integer         nf_get_var1_text
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             index(1),
! Subprogram not used !                          character*1         text)
! Subprogram not used       external        nf_get_var1_text
! Subprogram not used 
! Subprogram not used       integer         nf_put_var1_int1
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             index(1),
! Subprogram not used !                          nf_int1_t           i1val)
! Subprogram not used       external        nf_put_var1_int1
! Subprogram not used 
! Subprogram not used       integer         nf_get_var1_int1
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             index(1),
! Subprogram not used !                          nf_int1_t           i1val)
! Subprogram not used       external        nf_get_var1_int1
! Subprogram not used 
! Subprogram not used       integer         nf_put_var1_int2
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             index(1),
! Subprogram not used !                          nf_int2_t           i2val)
! Subprogram not used       external        nf_put_var1_int2
! Subprogram not used 
! Subprogram not used       integer         nf_get_var1_int2
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             index(1),
! Subprogram not used !                          nf_int2_t           i2val)
! Subprogram not used       external        nf_get_var1_int2
! Subprogram not used 
! Subprogram not used       integer         nf_put_var1_int
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             index(1),
! Subprogram not used !                          integer             ival)
! Subprogram not used       external        nf_put_var1_int
! Subprogram not used 
! Subprogram not used       integer         nf_get_var1_int
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             index(1),
! Subprogram not used !                          integer             ival)
! Subprogram not used       external        nf_get_var1_int
! Subprogram not used 
! Subprogram not used       integer         nf_put_var1_real
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             index(1),
! Subprogram not used !                          real                rval)
! Subprogram not used       external        nf_put_var1_real
! Subprogram not used 
! Subprogram not used       integer         nf_get_var1_real
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             index(1),
! Subprogram not used !                          real                rval)
! Subprogram not used       external        nf_get_var1_real
! Subprogram not used 
! Subprogram not used       integer         nf_put_var1_double
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             index(1),
! Subprogram not used !                          doubleprecision     dval)
! Subprogram not used       external        nf_put_var1_double
! Subprogram not used 
! Subprogram not used       integer         nf_get_var1_double
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             index(1),
! Subprogram not used !                          doubleprecision     dval)
! Subprogram not used       external        nf_get_var1_double
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! variable array put/get routines:
! Subprogram not used !
! Subprogram not used 
! Subprogram not used       integer         nf_put_vara_text
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          character(*)        text)
! Subprogram not used       external        nf_put_vara_text
! Subprogram not used 
! Subprogram not used       integer         nf_get_vara_text
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          character(*)        text)
! Subprogram not used       external        nf_get_vara_text
! Subprogram not used 
! Subprogram not used       integer         nf_put_vara_int1
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          nf_int1_t           i1vals(1))
! Subprogram not used       external        nf_put_vara_int1
! Subprogram not used 
! Subprogram not used       integer         nf_get_vara_int1
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          nf_int1_t           i1vals(1))
! Subprogram not used       external        nf_get_vara_int1
! Subprogram not used 
! Subprogram not used       integer         nf_put_vara_int2
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          nf_int2_t           i2vals(1))
! Subprogram not used       external        nf_put_vara_int2
! Subprogram not used 
! Subprogram not used       integer         nf_get_vara_int2
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          nf_int2_t           i2vals(1))
! Subprogram not used       external        nf_get_vara_int2
! Subprogram not used 
! Subprogram not used       integer         nf_put_vara_int
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             ivals(1))
! Subprogram not used       external        nf_put_vara_int
! Subprogram not used 
! Subprogram not used       integer         nf_get_vara_int
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             ivals(1))
! Subprogram not used       external        nf_get_vara_int
! Subprogram not used 
! Subprogram not used       integer         nf_put_vara_real
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          real                rvals(1))
! Subprogram not used       external        nf_put_vara_real
! Subprogram not used 
! Subprogram not used       integer         nf_get_vara_real
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          real                rvals(1))
! Subprogram not used       external        nf_get_vara_real
! Subprogram not used 
! Subprogram not used       integer         nf_put_vara_double
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          doubleprecision     dvals(1))
! Subprogram not used       external        nf_put_vara_double
! Subprogram not used 
! Subprogram not used       integer         nf_get_vara_double
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          doubleprecision     dvals(1))
! Subprogram not used       external        nf_get_vara_double
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! strided variable put/get routines:
! Subprogram not used !
! Subprogram not used 
! Subprogram not used       integer         nf_put_vars_text
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          character(*)        text)
! Subprogram not used       external        nf_put_vars_text
! Subprogram not used 
! Subprogram not used       integer         nf_get_vars_text
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          character(*)        text)
! Subprogram not used       external        nf_get_vars_text
! Subprogram not used 
! Subprogram not used       integer         nf_put_vars_int1
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          nf_int1_t           i1vals(1))
! Subprogram not used       external        nf_put_vars_int1
! Subprogram not used 
! Subprogram not used       integer         nf_get_vars_int1
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          nf_int1_t           i1vals(1))
! Subprogram not used       external        nf_get_vars_int1
! Subprogram not used 
! Subprogram not used       integer         nf_put_vars_int2
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          nf_int2_t           i2vals(1))
! Subprogram not used       external        nf_put_vars_int2
! Subprogram not used 
! Subprogram not used       integer         nf_get_vars_int2
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          nf_int2_t           i2vals(1))
! Subprogram not used       external        nf_get_vars_int2
! Subprogram not used 
! Subprogram not used       integer         nf_put_vars_int
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          integer             ivals(1))
! Subprogram not used       external        nf_put_vars_int
! Subprogram not used 
! Subprogram not used       integer         nf_get_vars_int
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          integer             ivals(1))
! Subprogram not used       external        nf_get_vars_int
! Subprogram not used 
! Subprogram not used       integer         nf_put_vars_real
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          real                rvals(1))
! Subprogram not used       external        nf_put_vars_real
! Subprogram not used 
! Subprogram not used       integer         nf_get_vars_real
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          real                rvals(1))
! Subprogram not used       external        nf_get_vars_real
! Subprogram not used 
! Subprogram not used       integer         nf_put_vars_double
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          doubleprecision     dvals(1))
! Subprogram not used       external        nf_put_vars_double
! Subprogram not used 
! Subprogram not used       integer         nf_get_vars_double
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          doubleprecision     dvals(1))
! Subprogram not used       external        nf_get_vars_double
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! mapped variable put/get routines:
! Subprogram not used !
! Subprogram not used 
! Subprogram not used       integer         nf_put_varm_text
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          integer             imap(1),
! Subprogram not used !                          character(*)        text)
! Subprogram not used       external        nf_put_varm_text
! Subprogram not used 
! Subprogram not used       integer         nf_get_varm_text
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          integer             imap(1),
! Subprogram not used !                          character(*)        text)
! Subprogram not used       external        nf_get_varm_text
! Subprogram not used 
! Subprogram not used       integer         nf_put_varm_int1
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          integer             imap(1),
! Subprogram not used !                          nf_int1_t           i1vals(1))
! Subprogram not used       external        nf_put_varm_int1
! Subprogram not used 
! Subprogram not used       integer         nf_get_varm_int1
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          integer             imap(1),
! Subprogram not used !                          nf_int1_t           i1vals(1))
! Subprogram not used       external        nf_get_varm_int1
! Subprogram not used 
! Subprogram not used       integer         nf_put_varm_int2
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          integer             imap(1),
! Subprogram not used !                          nf_int2_t           i2vals(1))
! Subprogram not used       external        nf_put_varm_int2
! Subprogram not used 
! Subprogram not used       integer         nf_get_varm_int2
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          integer             imap(1),
! Subprogram not used !                          nf_int2_t           i2vals(1))
! Subprogram not used       external        nf_get_varm_int2
! Subprogram not used 
! Subprogram not used       integer         nf_put_varm_int
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          integer             imap(1),
! Subprogram not used !                          integer             ivals(1))
! Subprogram not used       external        nf_put_varm_int
! Subprogram not used 
! Subprogram not used       integer         nf_get_varm_int
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          integer             imap(1),
! Subprogram not used !                          integer             ivals(1))
! Subprogram not used       external        nf_get_varm_int
! Subprogram not used 
! Subprogram not used       integer         nf_put_varm_real
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          integer             imap(1),
! Subprogram not used !                          real                rvals(1))
! Subprogram not used       external        nf_put_varm_real
! Subprogram not used 
! Subprogram not used       integer         nf_get_varm_real
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          integer             imap(1),
! Subprogram not used !                          real                rvals(1))
! Subprogram not used       external        nf_get_varm_real
! Subprogram not used 
! Subprogram not used       integer         nf_put_varm_double
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          integer             imap(1),
! Subprogram not used !                          doubleprecision     dvals(1))
! Subprogram not used       external        nf_put_varm_double
! Subprogram not used 
! Subprogram not used       integer         nf_get_varm_double
! Subprogram not used !                         (integer             ncid,
! Subprogram not used !                          integer             varid,
! Subprogram not used !                          integer             start(1),
! Subprogram not used !                          integer             count(1),
! Subprogram not used !                          integer             stride(1),
! Subprogram not used !                          integer             imap(1),
! Subprogram not used !                          doubleprecision     dvals(1))
! Subprogram not used       external        nf_get_varm_double
! Subprogram not used 
! Subprogram not used 
! Subprogram not used !     NetCDF-4.
! Subprogram not used !     This is part of netCDF-4. Copyright 2006, UCAR, See COPYRIGHT
! Subprogram not used !     file for distribution information.
! Subprogram not used 
! Subprogram not used !     Netcdf version 4 fortran interface.
! Subprogram not used 
! Subprogram not used !     $Id: netcdf4.inc,v 1.28 2010/05/25 13:53:02 ed Exp $
! Subprogram not used 
! Subprogram not used !     New netCDF-4 types.
! Subprogram not used       integer nf_ubyte
! Subprogram not used       integer nf_ushort
! Subprogram not used       integer nf_uint
! Subprogram not used       integer nf_int64
! Subprogram not used       integer nf_uint64
! Subprogram not used       integer nf_string
! Subprogram not used       integer nf_vlen
! Subprogram not used       integer nf_opaque
! Subprogram not used       integer nf_enum
! Subprogram not used       integer nf_compound
! Subprogram not used 
! Subprogram not used       parameter (nf_ubyte = 7)
! Subprogram not used       parameter (nf_ushort = 8)
! Subprogram not used       parameter (nf_uint = 9)
! Subprogram not used       parameter (nf_int64 = 10)
! Subprogram not used       parameter (nf_uint64 = 11)
! Subprogram not used       parameter (nf_string = 12)
! Subprogram not used       parameter (nf_vlen = 13)
! Subprogram not used       parameter (nf_opaque = 14)
! Subprogram not used       parameter (nf_enum = 15)
! Subprogram not used       parameter (nf_compound = 16)
! Subprogram not used 
! Subprogram not used !     New netCDF-4 fill values.
! Subprogram not used       integer           nf_fill_ubyte
! Subprogram not used       integer           nf_fill_ushort
! Subprogram not used !      real              nf_fill_uint
! Subprogram not used !      real              nf_fill_int64
! Subprogram not used !      real              nf_fill_uint64
! Subprogram not used       parameter (nf_fill_ubyte = 255)
! Subprogram not used       parameter (nf_fill_ushort = 65535)
! Subprogram not used 
! Subprogram not used !     New constants.
! Subprogram not used       integer nf_format_netcdf4
! Subprogram not used       parameter (nf_format_netcdf4 = 3)
! Subprogram not used 
! Subprogram not used       integer nf_format_netcdf4_classic
! Subprogram not used       parameter (nf_format_netcdf4_classic = 4)
! Subprogram not used 
! Subprogram not used       integer nf_netcdf4
! Subprogram not used       parameter (nf_netcdf4 = 4096)
! Subprogram not used 
! Subprogram not used       integer nf_classic_model
! Subprogram not used       parameter (nf_classic_model = 256)
! Subprogram not used 
! Subprogram not used       integer nf_chunk_seq
! Subprogram not used       parameter (nf_chunk_seq = 0)
! Subprogram not used       integer nf_chunk_sub
! Subprogram not used       parameter (nf_chunk_sub = 1)
! Subprogram not used       integer nf_chunk_sizes
! Subprogram not used       parameter (nf_chunk_sizes = 2)
! Subprogram not used 
! Subprogram not used       integer nf_endian_native
! Subprogram not used       parameter (nf_endian_native = 0)
! Subprogram not used       integer nf_endian_little
! Subprogram not used       parameter (nf_endian_little = 1)
! Subprogram not used       integer nf_endian_big
! Subprogram not used       parameter (nf_endian_big = 2)
! Subprogram not used 
! Subprogram not used !     For NF_DEF_VAR_CHUNKING
! Subprogram not used       integer nf_chunked
! Subprogram not used       parameter (nf_chunked = 0)
! Subprogram not used       integer nf_contiguous
! Subprogram not used       parameter (nf_contiguous = 1)
! Subprogram not used 
! Subprogram not used !     For NF_DEF_VAR_FLETCHER32
! Subprogram not used       integer nf_nochecksum
! Subprogram not used       parameter (nf_nochecksum = 0)
! Subprogram not used       integer nf_fletcher32
! Subprogram not used       parameter (nf_fletcher32 = 1)
! Subprogram not used 
! Subprogram not used !     For NF_DEF_VAR_DEFLATE
! Subprogram not used       integer nf_noshuffle
! Subprogram not used       parameter (nf_noshuffle = 0)
! Subprogram not used       integer nf_shuffle
! Subprogram not used       parameter (nf_shuffle = 1)
! Subprogram not used 
! Subprogram not used !     For NF_DEF_VAR_SZIP
! Subprogram not used       integer nf_szip_ec_option_mask
! Subprogram not used       parameter (nf_szip_ec_option_mask = 4)
! Subprogram not used       integer nf_szip_nn_option_mask
! Subprogram not used       parameter (nf_szip_nn_option_mask = 32)
! Subprogram not used 
! Subprogram not used !     For parallel I/O.
! Subprogram not used       integer nf_mpiio      
! Subprogram not used       parameter (nf_mpiio = 8192)
! Subprogram not used       integer nf_mpiposix
! Subprogram not used       parameter (nf_mpiposix = 16384)
! Subprogram not used       integer nf_pnetcdf
! Subprogram not used       parameter (nf_pnetcdf = 32768)
! Subprogram not used 
! Subprogram not used !     For NF_VAR_PAR_ACCESS.
! Subprogram not used       integer nf_independent
! Subprogram not used       parameter (nf_independent = 0)
! Subprogram not used       integer nf_collective
! Subprogram not used       parameter (nf_collective = 1)
! Subprogram not used 
! Subprogram not used !     New error codes.
! Subprogram not used       integer nf_ehdferr        ! Error at HDF5 layer. 
! Subprogram not used       parameter (nf_ehdferr = -101)
! Subprogram not used       integer nf_ecantread      ! Can't read. 
! Subprogram not used       parameter (nf_ecantread = -102)
! Subprogram not used       integer nf_ecantwrite     ! Can't write. 
! Subprogram not used       parameter (nf_ecantwrite = -103)
! Subprogram not used       integer nf_ecantcreate    ! Can't create. 
! Subprogram not used       parameter (nf_ecantcreate = -104)
! Subprogram not used       integer nf_efilemeta      ! Problem with file metadata. 
! Subprogram not used       parameter (nf_efilemeta = -105)
! Subprogram not used       integer nf_edimmeta       ! Problem with dimension metadata. 
! Subprogram not used       parameter (nf_edimmeta = -106)
! Subprogram not used       integer nf_eattmeta       ! Problem with attribute metadata. 
! Subprogram not used       parameter (nf_eattmeta = -107)
! Subprogram not used       integer nf_evarmeta       ! Problem with variable metadata. 
! Subprogram not used       parameter (nf_evarmeta = -108)
! Subprogram not used       integer nf_enocompound    ! Not a compound type. 
! Subprogram not used       parameter (nf_enocompound = -109)
! Subprogram not used       integer nf_eattexists     ! Attribute already exists. 
! Subprogram not used       parameter (nf_eattexists = -110)
! Subprogram not used       integer nf_enotnc4        ! Attempting netcdf-4 operation on netcdf-3 file.   
! Subprogram not used       parameter (nf_enotnc4 = -111)
! Subprogram not used       integer nf_estrictnc3     ! Attempting netcdf-4 operation on strict nc3 netcdf-4 file.   
! Subprogram not used       parameter (nf_estrictnc3 = -112)
! Subprogram not used       integer nf_enotnc3        ! Attempting netcdf-3 operation on netcdf-4 file.   
! Subprogram not used       parameter (nf_enotnc3 = -113)
! Subprogram not used       integer nf_enopar         ! Parallel operation on file opened for non-parallel access.   
! Subprogram not used       parameter (nf_enopar = -114)
! Subprogram not used       integer nf_eparinit       ! Error initializing for parallel access.   
! Subprogram not used       parameter (nf_eparinit = -115)
! Subprogram not used       integer nf_ebadgrpid      ! Bad group ID.   
! Subprogram not used       parameter (nf_ebadgrpid = -116)
! Subprogram not used       integer nf_ebadtypid      ! Bad type ID.   
! Subprogram not used       parameter (nf_ebadtypid = -117)
! Subprogram not used       integer nf_etypdefined    ! Type has already been defined and may not be edited. 
! Subprogram not used       parameter (nf_etypdefined = -118)
! Subprogram not used       integer nf_ebadfield      ! Bad field ID.   
! Subprogram not used       parameter (nf_ebadfield = -119)
! Subprogram not used       integer nf_ebadclass      ! Bad class.   
! Subprogram not used       parameter (nf_ebadclass = -120)
! Subprogram not used       integer nf_emaptype       ! Mapped access for atomic types only.   
! Subprogram not used       parameter (nf_emaptype = -121)
! Subprogram not used       integer nf_elatefill      ! Attempt to define fill value when data already exists. 
! Subprogram not used       parameter (nf_elatefill = -122)
! Subprogram not used       integer nf_elatedef       ! Attempt to define var properties, like deflate, after enddef. 
! Subprogram not used       parameter (nf_elatedef = -123)
! Subprogram not used       integer nf_edimscale      ! Probem with HDF5 dimscales. 
! Subprogram not used       parameter (nf_edimscale = -124)
! Subprogram not used       integer nf_enogrp       ! No group found.
! Subprogram not used       parameter (nf_enogrp = -125)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used !     New functions.
! Subprogram not used 
! Subprogram not used !     Parallel I/O.
! Subprogram not used       integer nf_create_par
! Subprogram not used       external nf_create_par
! Subprogram not used 
! Subprogram not used       integer nf_open_par
! Subprogram not used       external nf_open_par
! Subprogram not used 
! Subprogram not used       integer nf_var_par_access
! Subprogram not used       external nf_var_par_access
! Subprogram not used 
! Subprogram not used !     Functions to handle groups.
! Subprogram not used       integer nf_inq_ncid
! Subprogram not used       external nf_inq_ncid
! Subprogram not used 
! Subprogram not used       integer nf_inq_grps
! Subprogram not used       external nf_inq_grps
! Subprogram not used 
! Subprogram not used       integer nf_inq_grpname
! Subprogram not used       external nf_inq_grpname
! Subprogram not used 
! Subprogram not used       integer nf_inq_grpname_full
! Subprogram not used       external nf_inq_grpname_full
! Subprogram not used 
! Subprogram not used       integer nf_inq_grpname_len
! Subprogram not used       external nf_inq_grpname_len
! Subprogram not used 
! Subprogram not used       integer nf_inq_grp_parent
! Subprogram not used       external nf_inq_grp_parent
! Subprogram not used 
! Subprogram not used       integer nf_inq_grp_ncid
! Subprogram not used       external nf_inq_grp_ncid
! Subprogram not used 
! Subprogram not used       integer nf_inq_grp_full_ncid
! Subprogram not used       external nf_inq_grp_full_ncid
! Subprogram not used 
! Subprogram not used       integer nf_inq_varids
! Subprogram not used       external nf_inq_varids
! Subprogram not used 
! Subprogram not used       integer nf_inq_dimids
! Subprogram not used       external nf_inq_dimids
! Subprogram not used 
! Subprogram not used       integer nf_def_grp
! Subprogram not used       external nf_def_grp
! Subprogram not used 
! Subprogram not used !     New rename grp function
! Subprogram not used 
! Subprogram not used       integer nf_rename_grp
! Subprogram not used       external nf_rename_grp
! Subprogram not used 
! Subprogram not used !     New options for netCDF variables.
! Subprogram not used       integer nf_def_var_deflate
! Subprogram not used       external nf_def_var_deflate
! Subprogram not used 
! Subprogram not used       integer nf_inq_var_deflate
! Subprogram not used       external nf_inq_var_deflate
! Subprogram not used 
! Subprogram not used       integer nf_def_var_fletcher32
! Subprogram not used       external nf_def_var_fletcher32
! Subprogram not used 
! Subprogram not used       integer nf_inq_var_fletcher32
! Subprogram not used       external nf_inq_var_fletcher32
! Subprogram not used 
! Subprogram not used       integer nf_def_var_chunking
! Subprogram not used       external nf_def_var_chunking
! Subprogram not used 
! Subprogram not used       integer nf_inq_var_chunking
! Subprogram not used       external nf_inq_var_chunking
! Subprogram not used 
! Subprogram not used       integer nf_def_var_fill
! Subprogram not used       external nf_def_var_fill
! Subprogram not used 
! Subprogram not used       integer nf_inq_var_fill
! Subprogram not used       external nf_inq_var_fill
! Subprogram not used 
! Subprogram not used       integer nf_def_var_endian
! Subprogram not used       external nf_def_var_endian
! Subprogram not used 
! Subprogram not used       integer nf_inq_var_endian
! Subprogram not used       external nf_inq_var_endian
! Subprogram not used 
! Subprogram not used !     User defined types.
! Subprogram not used       integer nf_inq_typeids
! Subprogram not used       external nf_inq_typeids
! Subprogram not used 
! Subprogram not used       integer nf_inq_typeid
! Subprogram not used       external nf_inq_typeid
! Subprogram not used 
! Subprogram not used       integer nf_inq_type
! Subprogram not used       external nf_inq_type
! Subprogram not used 
! Subprogram not used       integer nf_inq_user_type
! Subprogram not used       external nf_inq_user_type
! Subprogram not used 
! Subprogram not used !     User defined types - compound types.
! Subprogram not used       integer nf_def_compound
! Subprogram not used       external nf_def_compound
! Subprogram not used 
! Subprogram not used       integer nf_insert_compound
! Subprogram not used       external nf_insert_compound
! Subprogram not used 
! Subprogram not used       integer nf_insert_array_compound
! Subprogram not used       external nf_insert_array_compound
! Subprogram not used 
! Subprogram not used       integer nf_inq_compound
! Subprogram not used       external nf_inq_compound
! Subprogram not used 
! Subprogram not used       integer nf_inq_compound_name
! Subprogram not used       external nf_inq_compound_name
! Subprogram not used 
! Subprogram not used       integer nf_inq_compound_size
! Subprogram not used       external nf_inq_compound_size
! Subprogram not used 
! Subprogram not used       integer nf_inq_compound_nfields
! Subprogram not used       external nf_inq_compound_nfields
! Subprogram not used 
! Subprogram not used       integer nf_inq_compound_field
! Subprogram not used       external nf_inq_compound_field
! Subprogram not used 
! Subprogram not used       integer nf_inq_compound_fieldname
! Subprogram not used       external nf_inq_compound_fieldname
! Subprogram not used 
! Subprogram not used       integer nf_inq_compound_fieldindex
! Subprogram not used       external nf_inq_compound_fieldindex
! Subprogram not used 
! Subprogram not used       integer nf_inq_compound_fieldoffset
! Subprogram not used       external nf_inq_compound_fieldoffset
! Subprogram not used 
! Subprogram not used       integer nf_inq_compound_fieldtype
! Subprogram not used       external nf_inq_compound_fieldtype
! Subprogram not used 
! Subprogram not used       integer nf_inq_compound_fieldndims
! Subprogram not used       external nf_inq_compound_fieldndims
! Subprogram not used 
! Subprogram not used       integer nf_inq_compound_fielddim_sizes
! Subprogram not used       external nf_inq_compound_fielddim_sizes
! Subprogram not used 
! Subprogram not used !     User defined types - variable length arrays.
! Subprogram not used       integer nf_def_vlen
! Subprogram not used       external nf_def_vlen
! Subprogram not used 
! Subprogram not used       integer nf_inq_vlen
! Subprogram not used       external nf_inq_vlen
! Subprogram not used 
! Subprogram not used       integer nf_free_vlen
! Subprogram not used       external nf_free_vlen
! Subprogram not used 
! Subprogram not used !     User defined types - enums.
! Subprogram not used       integer nf_def_enum
! Subprogram not used       external nf_def_enum
! Subprogram not used 
! Subprogram not used       integer nf_insert_enum
! Subprogram not used       external nf_insert_enum
! Subprogram not used 
! Subprogram not used       integer nf_inq_enum
! Subprogram not used       external nf_inq_enum
! Subprogram not used 
! Subprogram not used       integer nf_inq_enum_member
! Subprogram not used       external nf_inq_enum_member
! Subprogram not used 
! Subprogram not used       integer nf_inq_enum_ident
! Subprogram not used       external nf_inq_enum_ident
! Subprogram not used 
! Subprogram not used !     User defined types - opaque.
! Subprogram not used       integer nf_def_opaque
! Subprogram not used       external nf_def_opaque
! Subprogram not used 
! Subprogram not used       integer nf_inq_opaque
! Subprogram not used       external nf_inq_opaque
! Subprogram not used 
! Subprogram not used !     Write and read attributes of any type, including user defined
! Subprogram not used !     types.
! Subprogram not used       integer nf_put_att
! Subprogram not used       external nf_put_att
! Subprogram not used       integer nf_get_att
! Subprogram not used       external nf_get_att
! Subprogram not used 
! Subprogram not used !     Write and read variables of any type, including user defined
! Subprogram not used !     types.
! Subprogram not used       integer nf_put_var
! Subprogram not used       external nf_put_var
! Subprogram not used       integer nf_put_var1
! Subprogram not used       external nf_put_var1
! Subprogram not used       integer nf_put_vara
! Subprogram not used       external nf_put_vara
! Subprogram not used       integer nf_put_vars
! Subprogram not used       external nf_put_vars
! Subprogram not used       integer nf_get_var
! Subprogram not used       external nf_get_var
! Subprogram not used       integer nf_get_var1
! Subprogram not used       external nf_get_var1
! Subprogram not used       integer nf_get_vara
! Subprogram not used       external nf_get_vara
! Subprogram not used       integer nf_get_vars
! Subprogram not used       external nf_get_vars
! Subprogram not used 
! Subprogram not used !     64-bit int functions.
! Subprogram not used       integer nf_put_var1_int64
! Subprogram not used       external nf_put_var1_int64
! Subprogram not used       integer nf_put_vara_int64
! Subprogram not used       external nf_put_vara_int64
! Subprogram not used       integer nf_put_vars_int64
! Subprogram not used       external nf_put_vars_int64
! Subprogram not used       integer nf_put_varm_int64
! Subprogram not used       external nf_put_varm_int64
! Subprogram not used       integer nf_put_var_int64
! Subprogram not used       external nf_put_var_int64
! Subprogram not used       integer nf_get_var1_int64
! Subprogram not used       external nf_get_var1_int64
! Subprogram not used       integer nf_get_vara_int64
! Subprogram not used       external nf_get_vara_int64
! Subprogram not used       integer nf_get_vars_int64
! Subprogram not used       external nf_get_vars_int64
! Subprogram not used       integer nf_get_varm_int64
! Subprogram not used       external nf_get_varm_int64
! Subprogram not used       integer nf_get_var_int64
! Subprogram not used       external nf_get_var_int64
! Subprogram not used 
! Subprogram not used !     For helping F77 users with VLENs.
! Subprogram not used       integer nf_get_vlen_element
! Subprogram not used       external nf_get_vlen_element
! Subprogram not used       integer nf_put_vlen_element
! Subprogram not used       external nf_put_vlen_element
! Subprogram not used 
! Subprogram not used !     For dealing with file level chunk cache.
! Subprogram not used       integer nf_set_chunk_cache
! Subprogram not used       external nf_set_chunk_cache
! Subprogram not used       integer nf_get_chunk_cache
! Subprogram not used       external nf_get_chunk_cache
! Subprogram not used 
! Subprogram not used !     For dealing with per variable chunk cache.
! Subprogram not used       integer nf_set_var_chunk_cache
! Subprogram not used       external nf_set_var_chunk_cache
! Subprogram not used       integer nf_get_var_chunk_cache
! Subprogram not used       external nf_get_var_chunk_cache
! Subprogram not used 
! Subprogram not used !     NetCDF-2.
! Subprogram not used !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Subprogram not used ! begin netcdf 2.4 backward compatibility:
! Subprogram not used !
! Subprogram not used 
! Subprogram not used !      
! Subprogram not used ! functions in the fortran interface
! Subprogram not used !
! Subprogram not used       integer nccre
! Subprogram not used       integer ncopn
! Subprogram not used       integer ncddef
! Subprogram not used       integer ncdid
! Subprogram not used       integer ncvdef
! Subprogram not used       integer ncvid
! Subprogram not used       integer nctlen
! Subprogram not used       integer ncsfil
! Subprogram not used 
! Subprogram not used       external nccre
! Subprogram not used       external ncopn
! Subprogram not used       external ncddef
! Subprogram not used       external ncdid
! Subprogram not used       external ncvdef
! Subprogram not used       external ncvid
! Subprogram not used       external nctlen
! Subprogram not used       external ncsfil
! Subprogram not used 
! Subprogram not used 
! Subprogram not used       integer ncrdwr
! Subprogram not used       integer nccreat
! Subprogram not used       integer ncexcl
! Subprogram not used       integer ncindef
! Subprogram not used       integer ncnsync
! Subprogram not used       integer nchsync
! Subprogram not used       integer ncndirty
! Subprogram not used       integer nchdirty
! Subprogram not used       integer nclink
! Subprogram not used       integer ncnowrit
! Subprogram not used       integer ncwrite
! Subprogram not used       integer ncclob
! Subprogram not used       integer ncnoclob
! Subprogram not used       integer ncglobal
! Subprogram not used       integer ncfill
! Subprogram not used       integer ncnofill
! Subprogram not used       integer maxncop
! Subprogram not used       integer maxncdim
! Subprogram not used       integer maxncatt
! Subprogram not used       integer maxncvar
! Subprogram not used       integer maxncnam
! Subprogram not used       integer maxvdims
! Subprogram not used       integer ncnoerr
! Subprogram not used       integer ncebadid
! Subprogram not used       integer ncenfile
! Subprogram not used       integer nceexist
! Subprogram not used       integer nceinval
! Subprogram not used       integer nceperm
! Subprogram not used       integer ncenotin
! Subprogram not used       integer nceindef
! Subprogram not used       integer ncecoord
! Subprogram not used       integer ncemaxds
! Subprogram not used       integer ncename
! Subprogram not used       integer ncenoatt
! Subprogram not used       integer ncemaxat
! Subprogram not used       integer ncebadty
! Subprogram not used       integer ncebadd
! Subprogram not used       integer ncests
! Subprogram not used       integer nceunlim
! Subprogram not used       integer ncemaxvs
! Subprogram not used       integer ncenotvr
! Subprogram not used       integer nceglob
! Subprogram not used       integer ncenotnc
! Subprogram not used       integer ncfoobar
! Subprogram not used       integer ncsyserr
! Subprogram not used       integer ncfatal
! Subprogram not used       integer ncverbos
! Subprogram not used       integer ncentool
! Subprogram not used 
! Subprogram not used 
! Subprogram not used !
! Subprogram not used ! netcdf data types:
! Subprogram not used !
! Subprogram not used       integer ncbyte
! Subprogram not used       integer ncchar
! Subprogram not used       integer ncshort
! Subprogram not used       integer nclong
! Subprogram not used       integer ncfloat
! Subprogram not used       integer ncdouble
! Subprogram not used 
! Subprogram not used       parameter(ncbyte = 1)
! Subprogram not used       parameter(ncchar = 2)
! Subprogram not used       parameter(ncshort = 3)
! Subprogram not used       parameter(nclong = 4)
! Subprogram not used       parameter(ncfloat = 5)
! Subprogram not used       parameter(ncdouble = 6)
! Subprogram not used 
! Subprogram not used !     
! Subprogram not used !     masks for the struct nc flag field; passed in as 'mode' arg to
! Subprogram not used !     nccreate and ncopen.
! Subprogram not used !     
! Subprogram not used 
! Subprogram not used !     read/write, 0 => readonly 
! Subprogram not used       parameter(ncrdwr = 1)
! Subprogram not used !     in create phase, cleared by ncendef 
! Subprogram not used       parameter(nccreat = 2)
! Subprogram not used !     on create destroy existing file 
! Subprogram not used       parameter(ncexcl = 4)
! Subprogram not used !     in define mode, cleared by ncendef 
! Subprogram not used       parameter(ncindef = 8)
! Subprogram not used !     synchronise numrecs on change (x'10')
! Subprogram not used       parameter(ncnsync = 16)
! Subprogram not used !     synchronise whole header on change (x'20')
! Subprogram not used       parameter(nchsync = 32)
! Subprogram not used !     numrecs has changed (x'40')
! Subprogram not used       parameter(ncndirty = 64)  
! Subprogram not used !     header info has changed (x'80')
! Subprogram not used       parameter(nchdirty = 128)
! Subprogram not used !     prefill vars on endef and increase of record, the default behavior
! Subprogram not used       parameter(ncfill = 0)
! Subprogram not used !     do not fill vars on endef and increase of record (x'100')
! Subprogram not used       parameter(ncnofill = 256)
! Subprogram not used !     isa link (x'8000')
! Subprogram not used       parameter(nclink = 32768)
! Subprogram not used 
! Subprogram not used !     
! Subprogram not used !     'mode' arguments for nccreate and ncopen
! Subprogram not used !     
! Subprogram not used       parameter(ncnowrit = 0)
! Subprogram not used       parameter(ncwrite = ncrdwr)
! Subprogram not used       parameter(ncclob = nf_clobber)
! Subprogram not used       parameter(ncnoclob = nf_noclobber)
! Subprogram not used 
! Subprogram not used !     
! Subprogram not used !     'size' argument to ncdimdef for an unlimited dimension
! Subprogram not used !     
! Subprogram not used       integer ncunlim
! Subprogram not used       parameter(ncunlim = 0)
! Subprogram not used 
! Subprogram not used !     
! Subprogram not used !     attribute id to put/get a global attribute
! Subprogram not used !     
! Subprogram not used       parameter(ncglobal  = 0)
! Subprogram not used 
! Subprogram not used !     
! Subprogram not used !     advisory maximums:
! Subprogram not used !     
! Subprogram not used       parameter(maxncop = 64)
! Subprogram not used       parameter(maxncdim = 1024)
! Subprogram not used       parameter(maxncatt = 8192)
! Subprogram not used       parameter(maxncvar = 8192)
! Subprogram not used !     not enforced 
! Subprogram not used       parameter(maxncnam = 256)
! Subprogram not used       parameter(maxvdims = maxncdim)
! Subprogram not used 
! Subprogram not used !     
! Subprogram not used !     global netcdf error status variable
! Subprogram not used !     initialized in error.c
! Subprogram not used !     
! Subprogram not used 
! Subprogram not used !     no error 
! Subprogram not used       parameter(ncnoerr = nf_noerr)
! Subprogram not used !     not a netcdf id 
! Subprogram not used       parameter(ncebadid = nf_ebadid)
! Subprogram not used !     too many netcdfs open 
! Subprogram not used       parameter(ncenfile = -31)   ! nc_syserr
! Subprogram not used !     netcdf file exists && ncnoclob
! Subprogram not used       parameter(nceexist = nf_eexist)
! Subprogram not used !     invalid argument 
! Subprogram not used       parameter(nceinval = nf_einval)
! Subprogram not used !     write to read only 
! Subprogram not used       parameter(nceperm = nf_eperm)
! Subprogram not used !     operation not allowed in data mode 
! Subprogram not used       parameter(ncenotin = nf_enotindefine )   
! Subprogram not used !     operation not allowed in define mode 
! Subprogram not used       parameter(nceindef = nf_eindefine)   
! Subprogram not used !     coordinates out of domain 
! Subprogram not used       parameter(ncecoord = nf_einvalcoords)
! Subprogram not used !     maxncdims exceeded 
! Subprogram not used       parameter(ncemaxds = nf_emaxdims)
! Subprogram not used !     string match to name in use 
! Subprogram not used       parameter(ncename = nf_enameinuse)   
! Subprogram not used !     attribute not found 
! Subprogram not used       parameter(ncenoatt = nf_enotatt)
! Subprogram not used !     maxncattrs exceeded 
! Subprogram not used       parameter(ncemaxat = nf_emaxatts)
! Subprogram not used !     not a netcdf data type 
! Subprogram not used       parameter(ncebadty = nf_ebadtype)
! Subprogram not used !     invalid dimension id 
! Subprogram not used       parameter(ncebadd = nf_ebaddim)
! Subprogram not used !     ncunlimited in the wrong index 
! Subprogram not used       parameter(nceunlim = nf_eunlimpos)
! Subprogram not used !     maxncvars exceeded 
! Subprogram not used       parameter(ncemaxvs = nf_emaxvars)
! Subprogram not used !     variable not found 
! Subprogram not used       parameter(ncenotvr = nf_enotvar)
! Subprogram not used !     action prohibited on ncglobal varid 
! Subprogram not used       parameter(nceglob = nf_eglobal)
! Subprogram not used !     not a netcdf file 
! Subprogram not used       parameter(ncenotnc = nf_enotnc)
! Subprogram not used       parameter(ncests = nf_ests)
! Subprogram not used       parameter (ncentool = nf_emaxname) 
! Subprogram not used       parameter(ncfoobar = 32)
! Subprogram not used       parameter(ncsyserr = -31)
! Subprogram not used 
! Subprogram not used !     
! Subprogram not used !     global options variable. used to determine behavior of error handler.
! Subprogram not used !     initialized in lerror.c
! Subprogram not used !     
! Subprogram not used       parameter(ncfatal = 1)
! Subprogram not used       parameter(ncverbos = 2)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used !     default fill values.  these must be the same as in the c interface.
! Subprogram not used !
! Subprogram not used       integer filbyte
! Subprogram not used       integer filchar
! Subprogram not used       integer filshort
! Subprogram not used       integer fillong
! Subprogram not used       real filfloat
! Subprogram not used       doubleprecision fildoub
! Subprogram not used 
! Subprogram not used       parameter (filbyte = -127)
! Subprogram not used       parameter (filchar = 0)
! Subprogram not used       parameter (filshort = -32767)
! Subprogram not used       parameter (fillong = -2147483647)
! Subprogram not used       parameter (filfloat = 9.9692099683868690e+36)
! Subprogram not used       parameter (fildoub = 9.9692099683868690e+36)
! Subprogram not used 
! Subprogram not used !
! Subprogram not used !    Copyright (C) Silicon Graphics International Corp.
! Subprogram not used !    All rights reserved.
! Subprogram not used !
! Subprogram not used !    SGI RESERVES THE RIGHT TO WITHDRAW, MODIFY, OR REPLACE THIS SOFTWARE AT
! Subprogram not used !    ANY TIME, WITHOUT NOTICE. THE SOFTWARE IS "AS IS." IN CONNECTION WITH OR
! Subprogram not used !    ARISING IN RELATION TO THE SOFTWARE AND/OR THIS NOTICE, (1) IN NO EVENT
! Subprogram not used !    SHALL SGI OR ITS SUPPLIERS BE LIABLE FOR ANY SPECIAL, CONSEQUENTIAL,
! Subprogram not used !    INCIDENTAL OR INDIRECT DAMAGES, EVEN IF PRE-ADVISED OF THEIR PROSPECT,
! Subprogram not used !    HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY; AND, (2) SGI AND
! Subprogram not used !    ITS SUPPLIERS DISCLAIM ANY AND ALL LIABILITY FOR: (a) WARRANTIES AND
! Subprogram not used !    CONDITIONS, WHETHER EXPRESSED, IMPLIED, OR STATUTORY,  ARISING IN RELATION
! Subprogram not used !    TO THE SOFTWARE AND/OR THIS NOTICE, INCLUDING WITHOUT LIMITATION ANY
! Subprogram not used !    WARRANTY AND/OR CONDITION OF ERROR-FREE AND/OR UNINTERRUPTED OPERATION,
! Subprogram not used !    MERCHANTABILITY, SATISFACTORY QUALITY, FITNESS FOR A PARTICULAR PURPOSE,
! Subprogram not used !    AND NON-INFRINGEMENT; AND (b) INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
! Subprogram not used !    OR CONSEQUENTIAL DAMAGES RELATING TO THE SOFTWARE, ITS USE, MISUSE,
! Subprogram not used !    AND/OR FAILURE OF USE. ALL OF THE FOREGOING APPLY NOTWITHSTANDING THE
! Subprogram not used !    FAILURE OF ESSENTIAL PURPOSE OF ANY CONTRACTUAL REMEDY.
! Subprogram not used !
! Subprogram not used ! Copyright Notice
! Subprogram not used !  + 1993 University of Chicago
! Subprogram not used !  + 1993 Mississippi State University
! Subprogram not used 
! Subprogram not used !
! Subprogram not used !    Copyright (C) Silicon Graphics International Corp.
! Subprogram not used !    All rights reserved.
! Subprogram not used !
! Subprogram not used !    SGI RESERVES THE RIGHT TO WITHDRAW, MODIFY, OR REPLACE THIS SOFTWARE AT
! Subprogram not used !    ANY TIME, WITHOUT NOTICE. THE SOFTWARE IS "AS IS." IN CONNECTION WITH OR
! Subprogram not used !    ARISING IN RELATION TO THE SOFTWARE AND/OR THIS NOTICE, (1) IN NO EVENT
! Subprogram not used !    SHALL SGI OR ITS SUPPLIERS BE LIABLE FOR ANY SPECIAL, CONSEQUENTIAL,
! Subprogram not used !    INCIDENTAL OR INDIRECT DAMAGES, EVEN IF PRE-ADVISED OF THEIR PROSPECT,
! Subprogram not used !    HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY; AND, (2) SGI AND
! Subprogram not used !    ITS SUPPLIERS DISCLAIM ANY AND ALL LIABILITY FOR: (a) WARRANTIES AND
! Subprogram not used !    CONDITIONS, WHETHER EXPRESSED, IMPLIED, OR STATUTORY,  ARISING IN RELATION
! Subprogram not used !    TO THE SOFTWARE AND/OR THIS NOTICE, INCLUDING WITHOUT LIMITATION ANY
! Subprogram not used !    WARRANTY AND/OR CONDITION OF ERROR-FREE AND/OR UNINTERRUPTED OPERATION,
! Subprogram not used !    MERCHANTABILITY, SATISFACTORY QUALITY, FITNESS FOR A PARTICULAR PURPOSE,
! Subprogram not used !    AND NON-INFRINGEMENT; AND (b) INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
! Subprogram not used !    OR CONSEQUENTIAL DAMAGES RELATING TO THE SOFTWARE, ITS USE, MISUSE,
! Subprogram not used !    AND/OR FAILURE OF USE. ALL OF THE FOREGOING APPLY NOTWITHSTANDING THE
! Subprogram not used !    FAILURE OF ESSENTIAL PURPOSE OF ANY CONTRACTUAL REMEDY.
! Subprogram not used !
! Subprogram not used 
! Subprogram not used         integer MPI_VERSION
! Subprogram not used         integer MPI_SUBVERSION
! Subprogram not used         logical MPI_SUBARRAYS_SUPPORTED
! Subprogram not used         logical MPI_ASYNC_PROTECTS_NONBLOCKING
! Subprogram not used 
! Subprogram not used         parameter (MPI_VERSION                    = 3)
! Subprogram not used         parameter (MPI_SUBVERSION                 = 1)
! Subprogram not used         parameter (MPI_SUBARRAYS_SUPPORTED        = .FALSE.)
! Subprogram not used         parameter (MPI_ASYNC_PROTECTS_NONBLOCKING = .FALSE.)
! Subprogram not used 
! Subprogram not used ! MPI_Status
! Subprogram not used 
! Subprogram not used         integer MPI_STATUS_SIZE
! Subprogram not used         parameter (MPI_STATUS_SIZE      = 6)
! Subprogram not used 
! Subprogram not used ! Misc Fortran declarations
! Subprogram not used 
! Subprogram not used         integer MPI_BOTTOM
! Subprogram not used         common /MPI_SGI_PRIVATE/ MPI_BOTTOM
! Subprogram not used 
! Subprogram not used         external MPI_NULL_COPY_FN, MPI_NULL_DELETE_FN, MPI_DUP_FN
! Subprogram not used 
! Subprogram not used         integer MPI_IN_PLACE(1)
! Subprogram not used         common /MPI_SGI_PRIVATE_INPLACE/ MPI_IN_PLACE
! Subprogram not used 
! Subprogram not used ! MPI-2 Section 5.3
! Subprogram not used 
! Subprogram not used         CHARACTER*(1) MPI_ARGV_NULL(1)
! Subprogram not used         common /MPI_SGI_PRIVATE_CHAR/ MPI_ARGV_NULL
! Subprogram not used 
! Subprogram not used         CHARACTER*(1) MPI_ARGVS_NULL(1)
! Subprogram not used         EQUIVALENCE(MPI_ARGV_NULL,MPI_ARGVS_NULL(1))
! Subprogram not used 
! Subprogram not used         integer MPI_ERRCODES_IGNORE(1)
! Subprogram not used         EQUIVALENCE(MPI_BOTTOM,MPI_ERRCODES_IGNORE(1))
! Subprogram not used 
! Subprogram not used ! MPI-1 error codes and classes
! Subprogram not used 
! Subprogram not used         integer MPI_SUCCESS
! Subprogram not used         integer MPI_ERR_BUFFER
! Subprogram not used         integer MPI_ERR_COUNT
! Subprogram not used         integer MPI_ERR_TYPE
! Subprogram not used         integer MPI_ERR_TAG
! Subprogram not used         integer MPI_ERR_COMM
! Subprogram not used         integer MPI_ERR_RANK
! Subprogram not used         integer MPI_ERR_REQUEST
! Subprogram not used         integer MPI_ERR_ROOT
! Subprogram not used         integer MPI_ERR_GROUP
! Subprogram not used         integer MPI_ERR_OP
! Subprogram not used         integer MPI_ERR_TOPOLOGY
! Subprogram not used         integer MPI_ERR_DIMS
! Subprogram not used         integer MPI_ERR_ARG
! Subprogram not used         integer MPI_ERR_UNKNOWN
! Subprogram not used         integer MPI_ERR_TRUNCATE
! Subprogram not used         integer MPI_ERR_OTHER
! Subprogram not used         integer MPI_ERR_INTERN
! Subprogram not used         integer MPI_ERR_IN_STATUS
! Subprogram not used         integer MPI_ERR_PENDING
! Subprogram not used 
! Subprogram not used         parameter (MPI_SUCCESS          = 0)
! Subprogram not used         parameter (MPI_ERR_BUFFER       = 1)
! Subprogram not used         parameter (MPI_ERR_COUNT        = 2)
! Subprogram not used         parameter (MPI_ERR_TYPE         = 3)
! Subprogram not used         parameter (MPI_ERR_TAG          = 4)
! Subprogram not used         parameter (MPI_ERR_COMM         = 5)
! Subprogram not used         parameter (MPI_ERR_RANK         = 6)
! Subprogram not used         parameter (MPI_ERR_REQUEST      = 7)
! Subprogram not used         parameter (MPI_ERR_ROOT         = 8)
! Subprogram not used         parameter (MPI_ERR_GROUP        = 9)
! Subprogram not used         parameter (MPI_ERR_OP           = 10)
! Subprogram not used         parameter (MPI_ERR_TOPOLOGY     = 11)
! Subprogram not used         parameter (MPI_ERR_DIMS         = 12)
! Subprogram not used         parameter (MPI_ERR_ARG          = 13)
! Subprogram not used         parameter (MPI_ERR_UNKNOWN      = 14)
! Subprogram not used         parameter (MPI_ERR_TRUNCATE     = 15)
! Subprogram not used         parameter (MPI_ERR_OTHER        = 16)
! Subprogram not used         parameter (MPI_ERR_INTERN       = 17)
! Subprogram not used         parameter (MPI_ERR_IN_STATUS    = 18)
! Subprogram not used         parameter (MPI_ERR_PENDING      = 19)
! Subprogram not used 
! Subprogram not used ! MPI-2 error codes and classes
! Subprogram not used 
! Subprogram not used         integer MPI_ERR_ACCESS
! Subprogram not used         integer MPI_ERR_AMODE
! Subprogram not used         integer MPI_ERR_ASSERT
! Subprogram not used         integer MPI_ERR_BAD_FILE
! Subprogram not used         integer MPI_ERR_BASE
! Subprogram not used         integer MPI_ERR_CONVERSION
! Subprogram not used         integer MPI_ERR_DISP
! Subprogram not used         integer MPI_ERR_DUP_DATAREP
! Subprogram not used         integer MPI_ERR_FILE_EXISTS
! Subprogram not used         integer MPI_ERR_FILE_IN_USE
! Subprogram not used         integer MPI_ERR_FILE
! Subprogram not used         integer MPI_ERR_INFO_KEY
! Subprogram not used         integer MPI_ERR_INFO_NOKEY
! Subprogram not used         integer MPI_ERR_INFO_VALUE
! Subprogram not used         integer MPI_ERR_INFO
! Subprogram not used         integer MPI_ERR_IO
! Subprogram not used         integer MPI_ERR_KEYVAL
! Subprogram not used         integer MPI_ERR_LOCKTYPE
! Subprogram not used         integer MPI_ERR_NAME
! Subprogram not used         integer MPI_ERR_NO_MEM
! Subprogram not used         integer MPI_ERR_NOT_SAME
! Subprogram not used         integer MPI_ERR_NO_SPACE
! Subprogram not used         integer MPI_ERR_NO_SUCH_FILE
! Subprogram not used         integer MPI_ERR_PORT
! Subprogram not used         integer MPI_ERR_QUOTA
! Subprogram not used         integer MPI_ERR_READ_ONLY
! Subprogram not used         integer MPI_ERR_RMA_CONFLICT
! Subprogram not used         integer MPI_ERR_RMA_SYNC
! Subprogram not used         integer MPI_ERR_SERVICE
! Subprogram not used         integer MPI_ERR_SIZE
! Subprogram not used         integer MPI_ERR_SPAWN
! Subprogram not used         integer MPI_ERR_UNSUPPORTED_DATAREP
! Subprogram not used         integer MPI_ERR_UNSUPPORTED_OPERATION
! Subprogram not used         integer MPI_ERR_WIN
! Subprogram not used         integer MPI_ERR_RMA_RANGE
! Subprogram not used         integer MPI_ERR_RMA_ATTACH
! Subprogram not used         integer MPI_ERR_RMA_SHARED
! Subprogram not used         integer MPI_ERR_RMA_FLAVOR
! Subprogram not used         integer MPI_T_ERR_CANNOT_INIT
! Subprogram not used         integer MPI_T_ERR_NOT_INITIALIZED
! Subprogram not used         integer MPI_T_ERR_MEMORY
! Subprogram not used         integer MPI_T_ERR_INVALID_INDEX
! Subprogram not used         integer MPI_T_ERR_INVALID_ITEM
! Subprogram not used         integer MPI_T_ERR_INVALID_SESSION
! Subprogram not used         integer MPI_T_ERR_INVALID_HANDLE
! Subprogram not used         integer MPI_T_ERR_OUT_OF_HANDLES
! Subprogram not used         integer MPI_T_ERR_OUT_OF_SESSIONS
! Subprogram not used         integer MPI_T_ERR_CVAR_SET_NOT_NOW
! Subprogram not used         integer MPI_T_ERR_CVAR_SET_NEVER
! Subprogram not used         integer MPI_T_ERR_PVAR_NO_WRITE
! Subprogram not used         integer MPI_T_ERR_PVAR_NO_STARTSTOP
! Subprogram not used         integer MPI_T_ERR_PVAR_NO_ATOMIC
! Subprogram not used         integer MPI_T_ERR_INVALID_NAME
! Subprogram not used         integer MPI_T_ERR_INVALID
! Subprogram not used 
! Subprogram not used         parameter (MPI_ERR_ACCESS               = 28)
! Subprogram not used         parameter (MPI_ERR_AMODE                = 29)
! Subprogram not used         parameter (MPI_ERR_ASSERT               = 30)
! Subprogram not used         parameter (MPI_ERR_BAD_FILE             = 31)
! Subprogram not used         parameter (MPI_ERR_BASE                 = 32)
! Subprogram not used         parameter (MPI_ERR_CONVERSION           = 33)
! Subprogram not used         parameter (MPI_ERR_DISP                 = 34)
! Subprogram not used         parameter (MPI_ERR_DUP_DATAREP          = 35)
! Subprogram not used         parameter (MPI_ERR_FILE_EXISTS          = 36)
! Subprogram not used         parameter (MPI_ERR_FILE_IN_USE          = 37)
! Subprogram not used         parameter (MPI_ERR_FILE                 = 38)
! Subprogram not used         parameter (MPI_ERR_INFO_KEY             = 39)
! Subprogram not used         parameter (MPI_ERR_INFO_NOKEY           = 40)
! Subprogram not used         parameter (MPI_ERR_INFO_VALUE           = 41)
! Subprogram not used         parameter (MPI_ERR_INFO                 = 42)
! Subprogram not used         parameter (MPI_ERR_IO                   = 43)
! Subprogram not used         parameter (MPI_ERR_KEYVAL               = 44)
! Subprogram not used         parameter (MPI_ERR_LOCKTYPE             = 45)
! Subprogram not used         parameter (MPI_ERR_NAME                 = 46)
! Subprogram not used         parameter (MPI_ERR_NO_MEM               = 47)
! Subprogram not used         parameter (MPI_ERR_NOT_SAME             = 48)
! Subprogram not used         parameter (MPI_ERR_NO_SPACE             = 49)
! Subprogram not used         parameter (MPI_ERR_NO_SUCH_FILE         = 50)
! Subprogram not used         parameter (MPI_ERR_PORT                 = 51)
! Subprogram not used         parameter (MPI_ERR_QUOTA                = 52)
! Subprogram not used         parameter (MPI_ERR_READ_ONLY            = 53)
! Subprogram not used         parameter (MPI_ERR_RMA_CONFLICT         = 54)
! Subprogram not used         parameter (MPI_ERR_RMA_SYNC             = 55)
! Subprogram not used         parameter (MPI_ERR_SERVICE              = 56)
! Subprogram not used         parameter (MPI_ERR_SIZE                 = 57)
! Subprogram not used         parameter (MPI_ERR_SPAWN                = 58)
! Subprogram not used         parameter (MPI_ERR_UNSUPPORTED_DATAREP  = 59)
! Subprogram not used         parameter (MPI_ERR_UNSUPPORTED_OPERATION= 60)
! Subprogram not used         parameter (MPI_ERR_WIN                  = 61)
! Subprogram not used         parameter (MPI_ERR_RMA_RANGE            = 62)
! Subprogram not used         parameter (MPI_ERR_RMA_ATTACH           = 63)
! Subprogram not used         parameter (MPI_ERR_RMA_SHARED           = 64)
! Subprogram not used         parameter (MPI_ERR_RMA_FLAVOR           = 65)
! Subprogram not used         parameter (MPI_T_ERR_CANNOT_INIT        = 66)
! Subprogram not used         parameter (MPI_T_ERR_NOT_INITIALIZED    = 67)
! Subprogram not used         parameter (MPI_T_ERR_MEMORY             = 68)
! Subprogram not used         parameter (MPI_T_ERR_INVALID_INDEX      = 69)
! Subprogram not used         parameter (MPI_T_ERR_INVALID_ITEM       = 70)
! Subprogram not used         parameter (MPI_T_ERR_INVALID_SESSION    = 71)
! Subprogram not used         parameter (MPI_T_ERR_INVALID_HANDLE     = 72)
! Subprogram not used         parameter (MPI_T_ERR_OUT_OF_HANDLES     = 73)
! Subprogram not used         parameter (MPI_T_ERR_OUT_OF_SESSIONS    = 74)
! Subprogram not used         parameter (MPI_T_ERR_CVAR_SET_NOT_NOW   = 75)
! Subprogram not used         parameter (MPI_T_ERR_CVAR_SET_NEVER     = 76)
! Subprogram not used         parameter (MPI_T_ERR_PVAR_NO_WRITE      = 77)
! Subprogram not used         parameter (MPI_T_ERR_PVAR_NO_STARTSTOP  = 78)
! Subprogram not used         parameter (MPI_T_ERR_PVAR_NO_ATOMIC     = 79)
! Subprogram not used         parameter (MPI_T_ERR_INVALID_NAME       = 80)
! Subprogram not used         parameter (MPI_T_ERR_INVALID            = 81)
! Subprogram not used 
! Subprogram not used         integer MPI_ERR_LASTCODE
! Subprogram not used         parameter (MPI_ERR_LASTCODE             = 100)
! Subprogram not used 
! Subprogram not used ! Permanent keyvals
! Subprogram not used 
! Subprogram not used         integer MPI_KEYVAL_INVALID
! Subprogram not used         integer MPI_TAG_UB
! Subprogram not used         integer MPI_HOST
! Subprogram not used         integer MPI_IO
! Subprogram not used         integer MPI_WTIME_IS_GLOBAL
! Subprogram not used         integer MPI_UNIVERSE_SIZE
! Subprogram not used         integer MPI_APPNUM
! Subprogram not used         integer MPI_LASTUSEDCODE
! Subprogram not used 
! Subprogram not used         parameter (MPI_KEYVAL_INVALID   = 0)
! Subprogram not used         parameter (MPI_TAG_UB           = 5)
! Subprogram not used         parameter (MPI_HOST             = 6)
! Subprogram not used         parameter (MPI_IO               = 7)
! Subprogram not used         parameter (MPI_WTIME_IS_GLOBAL  = 8)
! Subprogram not used         parameter (MPI_UNIVERSE_SIZE   = 10)
! Subprogram not used         parameter (MPI_APPNUM          = 12)
! Subprogram not used         parameter (MPI_LASTUSEDCODE    = 13)
! Subprogram not used 
! Subprogram not used ! Results of the compare operations
! Subprogram not used 
! Subprogram not used         integer MPI_IDENT
! Subprogram not used         integer MPI_CONGRUENT
! Subprogram not used         integer MPI_SIMILAR
! Subprogram not used         integer MPI_UNEQUAL
! Subprogram not used 
! Subprogram not used         parameter (MPI_IDENT            = 0)
! Subprogram not used         parameter (MPI_CONGRUENT        = 1)
! Subprogram not used         parameter (MPI_SIMILAR          = 2)
! Subprogram not used         parameter (MPI_UNEQUAL          = 3)
! Subprogram not used 
! Subprogram not used ! Topology types
! Subprogram not used 
! Subprogram not used         integer MPI_GRAPH
! Subprogram not used         integer MPI_CART
! Subprogram not used         integer MPI_DIST_GRAPH
! Subprogram not used 
! Subprogram not used         parameter (MPI_GRAPH    = 1)
! Subprogram not used         parameter (MPI_CART     = 2)
! Subprogram not used         parameter (MPI_DIST_GRAPH = 3)
! Subprogram not used 
! Subprogram not used         integer MPI_UNWEIGHTED
! Subprogram not used         common /MPI_SGI_PRIVATE_UNWEIGHTED/ MPI_UNWEIGHTED
! Subprogram not used         integer MPI_WEIGHTS_EMPTY
! Subprogram not used         common /MPI_SGI_PRIVATE_WEIGHTS_EMPTY/ MPI_WEIGHTS_EMPTY
! Subprogram not used 
! Subprogram not used ! Misc constants
! Subprogram not used 
! Subprogram not used         integer MPI_MAX_PROCESSOR_NAME
! Subprogram not used         parameter (MPI_MAX_PROCESSOR_NAME = 255)
! Subprogram not used 
! Subprogram not used         integer MPI_MAX_ERROR_STRING
! Subprogram not used         parameter (MPI_MAX_ERROR_STRING = 255)
! Subprogram not used 
! Subprogram not used         integer MPI_MAX_LIBRARY_VERSION_STRING
! Subprogram not used         parameter (MPI_MAX_LIBRARY_VERSION_STRING = 255)
! Subprogram not used 
! Subprogram not used         integer MPI_MAX_OBJECT_NAME
! Subprogram not used         parameter (MPI_MAX_OBJECT_NAME = 127)
! Subprogram not used 
! Subprogram not used         integer MPI_BSEND_OVERHEAD
! Subprogram not used         parameter (MPI_BSEND_OVERHEAD = 32)
! Subprogram not used 
! Subprogram not used         integer MPI_ROOT
! Subprogram not used         parameter (MPI_ROOT = -4)
! Subprogram not used 
! Subprogram not used         integer MPI_UNDEFINED
! Subprogram not used         parameter (MPI_UNDEFINED = -3)
! Subprogram not used 
! Subprogram not used         integer MPI_ANY_SOURCE
! Subprogram not used         parameter (MPI_ANY_SOURCE = -2)
! Subprogram not used 
! Subprogram not used         integer MPI_PROC_NULL
! Subprogram not used         parameter (MPI_PROC_NULL = -1)
! Subprogram not used 
! Subprogram not used         integer MPI_ANY_TAG
! Subprogram not used         parameter (MPI_ANY_TAG = -1)
! Subprogram not used 
! Subprogram not used ! The following 2 lines are included in the main mpif.h
! Subprogram not used !       double precision MPI_WTIME, MPI_WTICK
! Subprogram not used !       external MPI_WTIME, MPI_WTICK
! Subprogram not used 
! Subprogram not used ! MPI-2 Section 4.10
! Subprogram not used 
! Subprogram not used         integer MPI_MAX_INFO_KEY
! Subprogram not used         parameter (MPI_MAX_INFO_KEY = 254)
! Subprogram not used 
! Subprogram not used         integer MPI_MAX_INFO_VAL
! Subprogram not used         parameter (MPI_MAX_INFO_VAL = 1023)
! Subprogram not used 
! Subprogram not used ! MPI-2 Section 5.4
! Subprogram not used 
! Subprogram not used         integer MPI_MAX_PORT_NAME
! Subprogram not used         parameter (MPI_MAX_PORT_NAME = 255)
! Subprogram not used 
! Subprogram not used ! Kind values for MPI-2
! Subprogram not used 
! Subprogram not used         integer MPI_INTEGER_KIND
! Subprogram not used         parameter (MPI_INTEGER_KIND = 4)
! Subprogram not used 
! Subprogram not used         integer MPI_OFFSET_KIND
! Subprogram not used         parameter (MPI_OFFSET_KIND = 8)
! Subprogram not used 
! Subprogram not used         integer MPI_ADDRESS_KIND
! Subprogram not used         parameter (MPI_ADDRESS_KIND = 8)
! Subprogram not used 
! Subprogram not used         integer MPI_COUNT_KIND
! Subprogram not used         parameter (MPI_COUNT_KIND = 8)
! Subprogram not used 
! Subprogram not used ! Section 6.4 bindings for one-sided communication
! Subprogram not used 
! Subprogram not used         integer MPI_MODE_NOCHECK
! Subprogram not used         integer MPI_MODE_NOSTORE
! Subprogram not used         integer MPI_MODE_NOPUT
! Subprogram not used         integer MPI_MODE_NOPRECEDE
! Subprogram not used         integer MPI_MODE_NOSUCCEED
! Subprogram not used 
! Subprogram not used         parameter (MPI_MODE_NOCHECK             = 1)
! Subprogram not used         parameter (MPI_MODE_NOSTORE             = 2)
! Subprogram not used         parameter (MPI_MODE_NOPUT               = 4)
! Subprogram not used         parameter (MPI_MODE_NOPRECEDE           = 8)
! Subprogram not used         parameter (MPI_MODE_NOSUCCEED           = 16)
! Subprogram not used 
! Subprogram not used         integer MPI_LOCK_SHARED
! Subprogram not used         parameter (MPI_LOCK_SHARED              = 1)
! Subprogram not used         integer MPI_LOCK_EXCLUSIVE
! Subprogram not used         parameter (MPI_LOCK_EXCLUSIVE           = 2)
! Subprogram not used 
! Subprogram not used ! Thread-safety support levels
! Subprogram not used 
! Subprogram not used         integer MPI_THREAD_SINGLE
! Subprogram not used         integer MPI_THREAD_FUNNELED
! Subprogram not used         integer MPI_THREAD_SERIALIZED
! Subprogram not used         integer MPI_THREAD_MULTIPLE
! Subprogram not used 
! Subprogram not used         parameter (MPI_THREAD_SINGLE            = 0)
! Subprogram not used         parameter (MPI_THREAD_FUNNELED          = 1)
! Subprogram not used         parameter (MPI_THREAD_SERIALIZED        = 2)
! Subprogram not used         parameter (MPI_THREAD_MULTIPLE          = 3)
! Subprogram not used 
! Subprogram not used ! Datatype Decoding constants
! Subprogram not used 
! Subprogram not used         integer MPI_COMBINER_NAMED
! Subprogram not used         integer MPI_COMBINER_CONTIGUOUS
! Subprogram not used         integer MPI_COMBINER_VECTOR
! Subprogram not used         integer MPI_COMBINER_HVECTOR
! Subprogram not used         integer MPI_COMBINER_INDEXED
! Subprogram not used         integer MPI_COMBINER_HINDEXED
! Subprogram not used         integer MPI_COMBINER_STRUCT
! Subprogram not used         integer MPI_COMBINER_DARRAY
! Subprogram not used         integer MPI_COMBINER_DUP
! Subprogram not used         integer MPI_COMBINER_F90_COMPLEX
! Subprogram not used         integer MPI_COMBINER_F90_INTEGER
! Subprogram not used         integer MPI_COMBINER_F90_REAL
! Subprogram not used         integer MPI_COMBINER_HINDEXED_INTEGER
! Subprogram not used         integer MPI_COMBINER_HVECTOR_INTEGER
! Subprogram not used         integer MPI_COMBINER_INDEXED_BLOCK
! Subprogram not used         integer MPI_COMBINER_RESIZED
! Subprogram not used         integer MPI_COMBINER_STRUCT_INTEGER
! Subprogram not used         integer MPI_COMBINER_SUBARRAY
! Subprogram not used         integer MPI_COMBINER_HINDEXED_BLOCK
! Subprogram not used 
! Subprogram not used         parameter (MPI_COMBINER_NAMED            = (-1))
! Subprogram not used         parameter (MPI_COMBINER_CONTIGUOUS       = 0)
! Subprogram not used         parameter (MPI_COMBINER_VECTOR           = 1)
! Subprogram not used         parameter (MPI_COMBINER_HVECTOR          = 2)
! Subprogram not used         parameter (MPI_COMBINER_INDEXED          = 3)
! Subprogram not used         parameter (MPI_COMBINER_HINDEXED         = 4)
! Subprogram not used         parameter (MPI_COMBINER_STRUCT           = 5)
! Subprogram not used         parameter (MPI_COMBINER_DARRAY           = 6)
! Subprogram not used         parameter (MPI_COMBINER_DUP              = 7)
! Subprogram not used         parameter (MPI_COMBINER_F90_COMPLEX      = 8)
! Subprogram not used         parameter (MPI_COMBINER_F90_INTEGER      = 9)
! Subprogram not used         parameter (MPI_COMBINER_F90_REAL         = 10)
! Subprogram not used         parameter (MPI_COMBINER_HINDEXED_INTEGER = 11)
! Subprogram not used         parameter (MPI_COMBINER_HVECTOR_INTEGER  = 12)
! Subprogram not used         parameter (MPI_COMBINER_INDEXED_BLOCK    = 13)
! Subprogram not used         parameter (MPI_COMBINER_RESIZED          = 14)
! Subprogram not used         parameter (MPI_COMBINER_STRUCT_INTEGER   = 15)
! Subprogram not used         parameter (MPI_COMBINER_SUBARRAY         = 16)
! Subprogram not used         parameter (MPI_COMBINER_HINDEXED_BLOCK   = 17)
! Subprogram not used 
! Subprogram not used ! Permanent window keyvals
! Subprogram not used 
! Subprogram not used         integer MPI_WIN_BASE
! Subprogram not used         integer MPI_WIN_SIZE
! Subprogram not used         integer MPI_WIN_DISP_UNIT
! Subprogram not used         integer MPI_WIN_CREATE_FLAVOR
! Subprogram not used         integer MPI_WIN_MODEL
! Subprogram not used 
! Subprogram not used         parameter (MPI_WIN_BASE         = 4)
! Subprogram not used         parameter (MPI_WIN_SIZE         = 5)
! Subprogram not used         parameter (MPI_WIN_DISP_UNIT    = 6)
! Subprogram not used         parameter (MPI_WIN_CREATE_FLAVOR = 8)
! Subprogram not used         parameter (MPI_WIN_MODEL        = 10)
! Subprogram not used 
! Subprogram not used ! Typeclasses
! Subprogram not used 
! Subprogram not used         integer MPI_TYPECLASS_INTEGER
! Subprogram not used         integer MPI_TYPECLASS_REAL
! Subprogram not used         integer MPI_TYPECLASS_COMPLEX
! Subprogram not used 
! Subprogram not used         parameter (MPI_TYPECLASS_INTEGER = 1)
! Subprogram not used         parameter (MPI_TYPECLASS_REAL    = 2)
! Subprogram not used         parameter (MPI_TYPECLASS_COMPLEX = 3)
! Subprogram not used 
! Subprogram not used ! Communicator types
! Subprogram not used 
! Subprogram not used         integer MPI_COMM_TYPE_SHARED
! Subprogram not used 
! Subprogram not used         parameter (MPI_COMM_TYPE_SHARED = 1)
! Subprogram not used 
! Subprogram not used ! Window flavors
! Subprogram not used 
! Subprogram not used         integer MPI_WIN_FLAVOR_CREATE
! Subprogram not used         integer MPI_WIN_FLAVOR_ALLOCATE
! Subprogram not used         integer MPI_WIN_FLAVOR_DYNAMIC
! Subprogram not used         integer MPI_WIN_FLAVOR_SHARED
! Subprogram not used 
! Subprogram not used         parameter (MPI_WIN_FLAVOR_CREATE    = 1)
! Subprogram not used         parameter (MPI_WIN_FLAVOR_ALLOCATE  = 2)
! Subprogram not used         parameter (MPI_WIN_FLAVOR_DYNAMIC   = 3)
! Subprogram not used         parameter (MPI_WIN_FLAVOR_SHARED    = 4)
! Subprogram not used 
! Subprogram not used         integer MPI_WIN_SEPARATE
! Subprogram not used         integer MPI_WIN_UNIFIED
! Subprogram not used 
! Subprogram not used         parameter (MPI_WIN_SEPARATE = 1)
! Subprogram not used         parameter (MPI_WIN_UNIFIED  = 2)
! Subprogram not used 
! Subprogram not used ! MPI-2 I/O definitions
! Subprogram not used 
! Subprogram not used !
! Subprogram not used !    Copyright (C) Silicon Graphics International Corp.
! Subprogram not used !    All rights reserved.
! Subprogram not used !
! Subprogram not used !    SGI RESERVES THE RIGHT TO WITHDRAW, MODIFY, OR REPLACE THIS SOFTWARE AT
! Subprogram not used !    ANY TIME, WITHOUT NOTICE. THE SOFTWARE IS "AS IS." IN CONNECTION WITH OR
! Subprogram not used !    ARISING IN RELATION TO THE SOFTWARE AND/OR THIS NOTICE, (1) IN NO EVENT
! Subprogram not used !    SHALL SGI OR ITS SUPPLIERS BE LIABLE FOR ANY SPECIAL, CONSEQUENTIAL,
! Subprogram not used !    INCIDENTAL OR INDIRECT DAMAGES, EVEN IF PRE-ADVISED OF THEIR PROSPECT,
! Subprogram not used !    HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY; AND, (2) SGI AND
! Subprogram not used !    ITS SUPPLIERS DISCLAIM ANY AND ALL LIABILITY FOR: (a) WARRANTIES AND
! Subprogram not used !    CONDITIONS, WHETHER EXPRESSED, IMPLIED, OR STATUTORY,  ARISING IN RELATION
! Subprogram not used !    TO THE SOFTWARE AND/OR THIS NOTICE, INCLUDING WITHOUT LIMITATION ANY
! Subprogram not used !    WARRANTY AND/OR CONDITION OF ERROR-FREE AND/OR UNINTERRUPTED OPERATION,
! Subprogram not used !    MERCHANTABILITY, SATISFACTORY QUALITY, FITNESS FOR A PARTICULAR PURPOSE,
! Subprogram not used !    AND NON-INFRINGEMENT; AND (b) INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
! Subprogram not used !    OR CONSEQUENTIAL DAMAGES RELATING TO THE SOFTWARE, ITS USE, MISUSE,
! Subprogram not used !    AND/OR FAILURE OF USE. ALL OF THE FOREGOING APPLY NOTWITHSTANDING THE
! Subprogram not used !    FAILURE OF ESSENTIAL PURPOSE OF ANY CONTRACTUAL REMEDY.
! Subprogram not used !
! Subprogram not used !    Fortran MPI-IO programs 
! Subprogram not used !    Copyright (C) 1997 University of Chicago. 
! Subprogram not used !
! Subprogram not used       INTEGER MPI_MODE_RDONLY, MPI_MODE_RDWR, MPI_MODE_WRONLY
! Subprogram not used       INTEGER MPI_MODE_DELETE_ON_CLOSE, MPI_MODE_UNIQUE_OPEN
! Subprogram not used       INTEGER MPI_MODE_CREATE, MPI_MODE_EXCL
! Subprogram not used       INTEGER MPI_MODE_APPEND, MPI_MODE_SEQUENTIAL
! Subprogram not used       PARAMETER (MPI_MODE_RDONLY=2, MPI_MODE_RDWR=8, MPI_MODE_WRONLY=4)
! Subprogram not used       PARAMETER (MPI_MODE_CREATE=1, MPI_MODE_DELETE_ON_CLOSE=16)
! Subprogram not used       PARAMETER (MPI_MODE_UNIQUE_OPEN=32, MPI_MODE_EXCL=64)
! Subprogram not used       PARAMETER (MPI_MODE_APPEND=128, MPI_MODE_SEQUENTIAL=256)
! Subprogram not used !
! Subprogram not used       INTEGER MPI_FILE_NULL
! Subprogram not used       PARAMETER (MPI_FILE_NULL=0)
! Subprogram not used !
! Subprogram not used       INTEGER MPI_MAX_DATAREP_STRING
! Subprogram not used       PARAMETER (MPI_MAX_DATAREP_STRING=128)
! Subprogram not used !
! Subprogram not used       INTEGER MPI_SEEK_SET, MPI_SEEK_CUR, MPI_SEEK_END
! Subprogram not used       PARAMETER (MPI_SEEK_SET=600, MPI_SEEK_CUR=602, MPI_SEEK_END=604)
! Subprogram not used !
! Subprogram not used       INTEGER MPIO_REQUEST_NULL
! Subprogram not used       PARAMETER (MPIO_REQUEST_NULL=0)
! Subprogram not used !
! Subprogram not used !      INTEGER MPI_OFFSET_KIND
! Subprogram not used !      PARAMETER (MPI_OFFSET_KIND=8)
! Subprogram not used !
! Subprogram not used       integer(kind=8) MPI_DISPLACEMENT_CURRENT
! Subprogram not used       PARAMETER (MPI_DISPLACEMENT_CURRENT=-54278278)
! Subprogram not used 
! Subprogram not used       INTEGER MPI_ORDER_C, MPI_ORDER_FORTRAN
! Subprogram not used       PARAMETER (MPI_ORDER_C=56, MPI_ORDER_FORTRAN=57)
! Subprogram not used       INTEGER MPI_DISTRIBUTE_BLOCK, MPI_DISTRIBUTE_CYCLIC
! Subprogram not used       INTEGER MPI_DISTRIBUTE_NONE, MPI_DISTRIBUTE_DFLT_DARG
! Subprogram not used       PARAMETER (MPI_DISTRIBUTE_BLOCK=121, MPI_DISTRIBUTE_CYCLIC=122)
! Subprogram not used       PARAMETER (MPI_DISTRIBUTE_NONE=123)
! Subprogram not used       PARAMETER (MPI_DISTRIBUTE_DFLT_DARG=-49767)
! Subprogram not used 
! Subprogram not used       EXTERNAL MPI_CONVERSION_FN_NULL
! Subprogram not used !
! Subprogram not used !   End Fortran MPI-IO
! Subprogram not used 
! Subprogram not used !
! Subprogram not used !    Copyright (C) Silicon Graphics International Corp.
! Subprogram not used !    All rights reserved.
! Subprogram not used !
! Subprogram not used !    SGI RESERVES THE RIGHT TO WITHDRAW, MODIFY, OR REPLACE THIS SOFTWARE AT
! Subprogram not used !    ANY TIME, WITHOUT NOTICE. THE SOFTWARE IS "AS IS." IN CONNECTION WITH OR
! Subprogram not used !    ARISING IN RELATION TO THE SOFTWARE AND/OR THIS NOTICE, (1) IN NO EVENT
! Subprogram not used !    SHALL SGI OR ITS SUPPLIERS BE LIABLE FOR ANY SPECIAL, CONSEQUENTIAL,
! Subprogram not used !    INCIDENTAL OR INDIRECT DAMAGES, EVEN IF PRE-ADVISED OF THEIR PROSPECT,
! Subprogram not used !    HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY; AND, (2) SGI AND
! Subprogram not used !    ITS SUPPLIERS DISCLAIM ANY AND ALL LIABILITY FOR: (a) WARRANTIES AND
! Subprogram not used !    CONDITIONS, WHETHER EXPRESSED, IMPLIED, OR STATUTORY,  ARISING IN RELATION
! Subprogram not used !    TO THE SOFTWARE AND/OR THIS NOTICE, INCLUDING WITHOUT LIMITATION ANY
! Subprogram not used !    WARRANTY AND/OR CONDITION OF ERROR-FREE AND/OR UNINTERRUPTED OPERATION,
! Subprogram not used !    MERCHANTABILITY, SATISFACTORY QUALITY, FITNESS FOR A PARTICULAR PURPOSE,
! Subprogram not used !    AND NON-INFRINGEMENT; AND (b) INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
! Subprogram not used !    OR CONSEQUENTIAL DAMAGES RELATING TO THE SOFTWARE, ITS USE, MISUSE,
! Subprogram not used !    AND/OR FAILURE OF USE. ALL OF THE FOREGOING APPLY NOTWITHSTANDING THE
! Subprogram not used !    FAILURE OF ESSENTIAL PURPOSE OF ANY CONTRACTUAL REMEDY.
! Subprogram not used !
! Subprogram not used 
! Subprogram not used ! MPI_Status
! Subprogram not used 
! Subprogram not used         integer MPI_STATUS_IGNORE(MPI_STATUS_SIZE)
! Subprogram not used         integer MPI_STATUSES_IGNORE(MPI_STATUS_SIZE,1)
! Subprogram not used         equivalence (MPI_STATUS_IGNORE, MPI_STATUSES_IGNORE)
! Subprogram not used 
! Subprogram not used         common /MPI_SGI_PRIVATE_STATUS/ MPI_STATUS_IGNORE
! Subprogram not used 
! Subprogram not used ! Permanent window keyvals
! Subprogram not used 
! Subprogram not used         external MPI_COMM_NULL_COPY_FN, MPI_COMM_NULL_DELETE_FN
! Subprogram not used         external MPI_COMM_DUP_FN
! Subprogram not used         external MPI_WIN_NULL_COPY_FN, MPI_WIN_NULL_DELETE_FN
! Subprogram not used         external MPI_WIN_DUP_FN
! Subprogram not used 
! Subprogram not used !
! Subprogram not used         integer MPI_FLOAT_INT
! Subprogram not used         integer MPI_DOUBLE_INT
! Subprogram not used         integer MPI_LONG_INT
! Subprogram not used         integer MPI_2INT
! Subprogram not used         integer MPI_SHORT_INT
! Subprogram not used         integer MPI_LONG_DOUBLE_INT
! Subprogram not used 
! Subprogram not used         parameter (MPI_FLOAT_INT        = 31)
! Subprogram not used         parameter (MPI_DOUBLE_INT       = 32)
! Subprogram not used         parameter (MPI_LONG_INT         = 33)
! Subprogram not used         parameter (MPI_2INT             = 34)
! Subprogram not used         parameter (MPI_SHORT_INT        = 35)
! Subprogram not used         parameter (MPI_LONG_DOUBLE_INT  = 36)
! Subprogram not used 
! Subprogram not used         integer MPI_BYTE
! Subprogram not used         integer MPI_PACKED
! Subprogram not used         integer MPI_UB
! Subprogram not used         integer MPI_LB
! Subprogram not used 
! Subprogram not used         parameter (MPI_BYTE             = 27)
! Subprogram not used         parameter (MPI_PACKED           = 28)
! Subprogram not used         parameter (MPI_UB               = 29)
! Subprogram not used         parameter (MPI_LB               = 30)
! Subprogram not used 
! Subprogram not used         integer MPI_2REAL
! Subprogram not used         integer MPI_2INTEGER
! Subprogram not used 
! Subprogram not used         parameter (MPI_2REAL            = 37)
! Subprogram not used         parameter (MPI_2INTEGER         = 39)
! Subprogram not used 
! Subprogram not used         integer MPI_AINT
! Subprogram not used         integer MPI_OFFSET
! Subprogram not used 
! Subprogram not used         parameter (MPI_AINT             = 55)
! Subprogram not used         parameter (MPI_OFFSET           = 56)
! Subprogram not used 
! Subprogram not used ! MPI_Op
! Subprogram not used 
! Subprogram not used         integer MPI_OP_NULL
! Subprogram not used         integer MPI_MAX
! Subprogram not used         integer MPI_MIN
! Subprogram not used         integer MPI_SUM
! Subprogram not used         integer MPI_PROD
! Subprogram not used         integer MPI_LAND
! Subprogram not used         integer MPI_BAND
! Subprogram not used         integer MPI_LOR
! Subprogram not used         integer MPI_BOR
! Subprogram not used         integer MPI_LXOR
! Subprogram not used         integer MPI_BXOR
! Subprogram not used         integer MPI_MAXLOC
! Subprogram not used         integer MPI_MINLOC
! Subprogram not used         integer MPI_REPLACE
! Subprogram not used         integer MPI_NO_OP
! Subprogram not used 
! Subprogram not used         parameter (MPI_OP_NULL  = 0)
! Subprogram not used         parameter (MPI_MAX      = 1)
! Subprogram not used         parameter (MPI_MIN      = 2)
! Subprogram not used         parameter (MPI_SUM      = 3)
! Subprogram not used         parameter (MPI_PROD     = 4)
! Subprogram not used         parameter (MPI_LAND     = 5)
! Subprogram not used         parameter (MPI_BAND     = 6)
! Subprogram not used         parameter (MPI_LOR      = 7)
! Subprogram not used         parameter (MPI_BOR      = 8)
! Subprogram not used         parameter (MPI_LXOR     = 9)
! Subprogram not used         parameter (MPI_BXOR     = 10)
! Subprogram not used         parameter (MPI_MAXLOC   = 11)
! Subprogram not used         parameter (MPI_MINLOC   = 12)
! Subprogram not used         parameter (MPI_REPLACE  = 13)
! Subprogram not used         parameter (MPI_NO_OP    = 14)
! Subprogram not used 
! Subprogram not used ! MPI_Datatype
! Subprogram not used 
! Subprogram not used         integer MPI_DATATYPE_NULL
! Subprogram not used 
! Subprogram not used         integer MPI_CHAR
! Subprogram not used         integer MPI_SHORT
! Subprogram not used         integer MPI_INT
! Subprogram not used         integer MPI_LONG
! Subprogram not used         integer MPI_UNSIGNED_CHAR
! Subprogram not used         integer MPI_UNSIGNED_SHORT
! Subprogram not used         integer MPI_UNSIGNED
! Subprogram not used         integer MPI_UNSIGNED_LONG
! Subprogram not used         integer MPI_FLOAT
! Subprogram not used         integer MPI_DOUBLE
! Subprogram not used         integer MPI_LONG_DOUBLE
! Subprogram not used         integer MPI_LONG_LONG
! Subprogram not used         integer MPI_LONG_LONG_INT
! Subprogram not used 
! Subprogram not used         integer MPI_INTEGER
! Subprogram not used         integer MPI_REAL
! Subprogram not used         integer MPI_DOUBLE_PRECISION
! Subprogram not used         integer MPI_COMPLEX
! Subprogram not used         integer MPI_DOUBLE_COMPLEX
! Subprogram not used         integer MPI_LOGICAL
! Subprogram not used         integer MPI_CHARACTER
! Subprogram not used         integer MPI_INTEGER1
! Subprogram not used         integer MPI_INTEGER2
! Subprogram not used         integer MPI_INTEGER4
! Subprogram not used         integer MPI_INTEGER8
! Subprogram not used         integer MPI_REAL4
! Subprogram not used         integer MPI_REAL8
! Subprogram not used         integer MPI_REAL16
! Subprogram not used 
! Subprogram not used         integer MPI_2DOUBLE_PRECISION
! Subprogram not used 
! Subprogram not used         integer MPI_WCHAR
! Subprogram not used         integer MPI_SIGNED_CHAR
! Subprogram not used         integer MPI_UNSIGNED_LONG_LONG
! Subprogram not used 
! Subprogram not used         integer MPI_INTEGER16
! Subprogram not used         integer MPI_COMPLEX8
! Subprogram not used         integer MPI_COMPLEX16
! Subprogram not used         integer MPI_COMPLEX32
! Subprogram not used 
! Subprogram not used         integer MPI_INT8_T
! Subprogram not used         integer MPI_INT16_T
! Subprogram not used         integer MPI_INT32_T
! Subprogram not used         integer MPI_INT64_T
! Subprogram not used         integer MPI_UINT8_T
! Subprogram not used         integer MPI_UINT16_T
! Subprogram not used         integer MPI_UINT32_T
! Subprogram not used         integer MPI_UINT64_T
! Subprogram not used 
! Subprogram not used         integer MPI_C_BOOL
! Subprogram not used         integer MPI_C_FLOAT_COMPLEX
! Subprogram not used         integer MPI_C_COMPLEX
! Subprogram not used         integer MPI_C_DOUBLE_COMPLEX
! Subprogram not used         integer MPI_C_LONG_DOUBLE_COMPLEX
! Subprogram not used 
! Subprogram not used         integer MPI_COUNT
! Subprogram not used 
! Subprogram not used         integer MPI_CXX_BOOL
! Subprogram not used         integer MPI_CXX_FLOAT_COMPLEX
! Subprogram not used         integer MPI_CXX_DOUBLE_COMPLEX
! Subprogram not used         integer MPI_CXX_LONG_DOUBLE_COMPLEX
! Subprogram not used 
! Subprogram not used         parameter (MPI_DATATYPE_NULL    = 0)
! Subprogram not used 
! Subprogram not used         parameter (MPI_CHAR             = 1)
! Subprogram not used         parameter (MPI_SHORT            = 2)
! Subprogram not used         parameter (MPI_INT              = 3)
! Subprogram not used         parameter (MPI_LONG             = 4)
! Subprogram not used         parameter (MPI_UNSIGNED_CHAR    = 5)
! Subprogram not used         parameter (MPI_UNSIGNED_SHORT   = 6)
! Subprogram not used         parameter (MPI_UNSIGNED         = 7)
! Subprogram not used         parameter (MPI_UNSIGNED_LONG    = 8)
! Subprogram not used         parameter (MPI_FLOAT            = 9)
! Subprogram not used         parameter (MPI_DOUBLE           = 10)
! Subprogram not used         parameter (MPI_LONG_DOUBLE      = 11)
! Subprogram not used         parameter (MPI_LONG_LONG        = 12)
! Subprogram not used         parameter (MPI_LONG_LONG_INT    = 12)
! Subprogram not used 
! Subprogram not used         parameter (MPI_INTEGER          = 13)
! Subprogram not used         parameter (MPI_REAL             = 14)
! Subprogram not used         parameter (MPI_DOUBLE_PRECISION = 15)
! Subprogram not used         parameter (MPI_COMPLEX          = 16)
! Subprogram not used         parameter (MPI_DOUBLE_COMPLEX   = 17)
! Subprogram not used         parameter (MPI_LOGICAL          = 18)
! Subprogram not used         parameter (MPI_CHARACTER        = 19)
! Subprogram not used         parameter (MPI_INTEGER1         = 20)
! Subprogram not used         parameter (MPI_INTEGER2         = 21)
! Subprogram not used         parameter (MPI_INTEGER4         = 22)
! Subprogram not used         parameter (MPI_INTEGER8         = 23)
! Subprogram not used         parameter (MPI_REAL4            = 24)
! Subprogram not used         parameter (MPI_REAL8            = 25)
! Subprogram not used         parameter (MPI_REAL16           = 26)
! Subprogram not used 
! Subprogram not used         parameter (MPI_2DOUBLE_PRECISION= 38)
! Subprogram not used 
! Subprogram not used         parameter (MPI_WCHAR            = 40)
! Subprogram not used         parameter (MPI_SIGNED_CHAR      = 41)
! Subprogram not used         parameter (MPI_UNSIGNED_LONG_LONG = 42)
! Subprogram not used 
! Subprogram not used         parameter (MPI_INTEGER16        = 43)
! Subprogram not used         parameter (MPI_COMPLEX8         = 44)
! Subprogram not used         parameter (MPI_COMPLEX16        = 45)
! Subprogram not used         parameter (MPI_COMPLEX32        = 46)
! Subprogram not used 
! Subprogram not used         parameter (MPI_INT8_T           = 47)
! Subprogram not used         parameter (MPI_INT16_T          = 48)
! Subprogram not used         parameter (MPI_INT32_T          = 49)
! Subprogram not used         parameter (MPI_INT64_T          = 50)
! Subprogram not used         parameter (MPI_UINT8_T          = 51)
! Subprogram not used         parameter (MPI_UINT16_T         = 52)
! Subprogram not used         parameter (MPI_UINT32_T         = 53)
! Subprogram not used         parameter (MPI_UINT64_T         = 54)
! Subprogram not used 
! Subprogram not used         parameter (MPI_C_BOOL           = 57)
! Subprogram not used         parameter (MPI_C_FLOAT_COMPLEX  = 58)
! Subprogram not used         parameter (MPI_C_COMPLEX        = 58)
! Subprogram not used         parameter (MPI_C_DOUBLE_COMPLEX = 59)
! Subprogram not used         parameter (MPI_C_LONG_DOUBLE_COMPLEX = 60)
! Subprogram not used 
! Subprogram not used         parameter (MPI_COUNT            = 61)
! Subprogram not used 
! Subprogram not used         parameter (MPI_CXX_BOOL                = 62)
! Subprogram not used         parameter (MPI_CXX_FLOAT_COMPLEX       = 63)
! Subprogram not used         parameter (MPI_CXX_DOUBLE_COMPLEX      = 64)
! Subprogram not used         parameter (MPI_CXX_LONG_DOUBLE_COMPLEX = 65)
! Subprogram not used 
! Subprogram not used ! MPI_Comm
! Subprogram not used 
! Subprogram not used         integer MPI_COMM_NULL
! Subprogram not used         integer MPI_COMM_WORLD
! Subprogram not used         integer MPI_COMM_SELF
! Subprogram not used 
! Subprogram not used         parameter (MPI_COMM_NULL        = 0)
! Subprogram not used         parameter (MPI_COMM_WORLD       = 1)
! Subprogram not used         parameter (MPI_COMM_SELF        = 2)
! Subprogram not used 
! Subprogram not used ! MPI_Errhandler
! Subprogram not used 
! Subprogram not used         integer MPI_ERRHANDLER_NULL
! Subprogram not used         integer MPI_ERRORS_ARE_FATAL
! Subprogram not used         integer MPI_ERRORS_RETURN
! Subprogram not used 
! Subprogram not used         parameter (MPI_ERRHANDLER_NULL  = 0)
! Subprogram not used         parameter (MPI_ERRORS_ARE_FATAL = 1)
! Subprogram not used         parameter (MPI_ERRORS_RETURN    = 2)
! Subprogram not used 
! Subprogram not used         integer MPI_SOURCE
! Subprogram not used         integer MPI_TAG
! Subprogram not used         integer MPI_ERROR
! Subprogram not used 
! Subprogram not used         parameter (MPI_SOURCE           = 1)
! Subprogram not used         parameter (MPI_TAG              = 2)
! Subprogram not used         parameter (MPI_ERROR            = 3)
! Subprogram not used 
! Subprogram not used ! MPI_Group
! Subprogram not used 
! Subprogram not used         integer MPI_GROUP_NULL
! Subprogram not used         integer MPI_GROUP_EMPTY
! Subprogram not used 
! Subprogram not used         parameter (MPI_GROUP_NULL  = 0)
! Subprogram not used         parameter (MPI_GROUP_EMPTY = 1)
! Subprogram not used 
! Subprogram not used         integer MPI_INFO_NULL
! Subprogram not used         parameter (MPI_INFO_NULL = 0)
! Subprogram not used 
! Subprogram not used         integer MPI_INFO_ENV
! Subprogram not used         parameter (MPI_INFO_ENV = 1)
! Subprogram not used 
! Subprogram not used ! MPI_Message
! Subprogram not used 
! Subprogram not used         integer MPI_MESSAGE_NO_PROC
! Subprogram not used         integer MPI_MESSAGE_NULL
! Subprogram not used 
! Subprogram not used         parameter (MPI_MESSAGE_NO_PROC = -1)
! Subprogram not used         parameter (MPI_MESSAGE_NULL = 0)
! Subprogram not used 
! Subprogram not used ! MPI_Request
! Subprogram not used 
! Subprogram not used         integer MPI_REQUEST_NULL
! Subprogram not used         parameter (MPI_REQUEST_NULL = 0)
! Subprogram not used 
! Subprogram not used         integer MPI_WIN_NULL
! Subprogram not used         parameter (MPI_WIN_NULL = 0)
! Subprogram not used 
! Subprogram not used         double precision MPI_WTIME, MPI_WTICK, PMPI_WTIME, PMPI_WTICK
! Subprogram not used         external MPI_WTIME, MPI_WTICK, PMPI_WTIME, PMPI_WTICK
! Subprogram not used 
! Subprogram not used ! !INPUT/OUTPUT PARAMETERS:
! Subprogram not used 
! Subprogram not used    type(mct_sMat)  ,intent(in)   :: sMat     ! mapping data
! Subprogram not used    type(iosystem_desc_t)         :: iosystem ! PIO subsystem description
! Subprogram not used    integer(IN)     ,intent(in)   :: io_type  ! type of io interface for this file
! Subprogram not used    character(*)    ,intent(in)   :: filename ! netCDF file to read
! Subprogram not used    integer(IN)     ,intent(in)   :: compid   ! processor id
! Subprogram not used    integer(IN)     ,intent(in)   :: mpicom   ! communicator
! Subprogram not used 
! Subprogram not used  ! !local
! Subprogram not used    integer(IN) :: na,nb,ns,lsize,npes,ierr,my_task,n
! Subprogram not used    integer(IN), pointer :: start(:),count(:),ssize(:),pe_loc(:)
! Subprogram not used    integer(IN), pointer :: expvari(:)
! Subprogram not used    real(R8)   , pointer :: expvarr(:)
! Subprogram not used    type(mct_gsmap) :: gsmap
! Subprogram not used    type(mct_avect) :: AV
! Subprogram not used    character(*),parameter :: subName = '(shr_mct_sMatWritednc) '
! Subprogram not used 
! Subprogram not used !----------------------------------------
! Subprogram not used 
! Subprogram not used    call MPI_COMM_SIZE(mpicom,npes,ierr)
! Subprogram not used    call MPI_COMM_RANK(mpicom,my_task,ierr)
! Subprogram not used    allocate(start(npes),count(npes),ssize(npes),pe_loc(npes))
! Subprogram not used 
! Subprogram not used    na = mct_sMat_ncols(sMat)
! Subprogram not used    nb = mct_sMat_nrows(sMat)
! Subprogram not used    ns = mct_sMat_gNumEl(sMat,mpicom)
! Subprogram not used    lsize = mct_sMat_lsize(sMat)
! Subprogram not used 
! Subprogram not used    count(:) = -999
! Subprogram not used    pe_loc(:) = -999
! Subprogram not used    ssize(:) = 1
! Subprogram not used    call MPI_GATHER(lsize,1,MPI_INTEGER,count,ssize,MPI_INTEGER,0,mpicom,ierr)
! Subprogram not used 
! Subprogram not used    if (my_task == 0) then
! Subprogram not used       if (minval(count) < 0) then
! Subprogram not used          call shr_sys_abort(subname//' ERROR: count invalid')
! Subprogram not used       endif
! Subprogram not used 
! Subprogram not used       start(1) = 1
! Subprogram not used       pe_loc(1) = 0
! Subprogram not used       do n = 2,npes
! Subprogram not used          start(n) = start(n-1)+count(n-1)
! Subprogram not used          pe_loc(n) = n-1
! Subprogram not used       enddo
! Subprogram not used 
! Subprogram not used    endif
! Subprogram not used 
! Subprogram not used    call mct_gsmap_init(gsmap,npes,start,count,pe_loc,0,mpicom,compid,ns)
! Subprogram not used    deallocate(start,count,ssize,pe_loc)
! Subprogram not used 
! Subprogram not used    call mct_aVect_init(AV,iList='row:col',rList='S',lsize=lsize)
! Subprogram not used    allocate(expvari(lsize))
! Subprogram not used    call mct_sMat_ExpGRowI(sMat,expvari)
! Subprogram not used    AV%iAttr(1,:) = expvari(:)
! Subprogram not used    call mct_sMat_ExpGColI(sMat,expvari)
! Subprogram not used    AV%iAttr(2,:) = expvari(:)
! Subprogram not used    deallocate(expvari)
! Subprogram not used    allocate(expvarr(lsize))
! Subprogram not used    call mct_sMat_ExpMatrix(sMat,expvarr)
! Subprogram not used    AV%rAttr(1,:) = expvarr(:)
! Subprogram not used    deallocate(expvarr)
! Subprogram not used 
! Subprogram not used    call shr_pcdf_readwrite('write',iosystem,io_type, trim(filename),mpicom,gsmap,clobber=.false.,cdf64=.true., &
! Subprogram not used       id1=na,id1n='n_a',id2=nb,id2n='n_b',id3=ns,id3n='n_s',av1=AV,av1n='')
! Subprogram not used 
! Subprogram not used    call mct_gsmap_clean(gsmap)
! Subprogram not used    call mct_avect_clean(AV)
! Subprogram not used 
! Subprogram not used end subroutine shr_mct_sMatWritednc
!===============================================================================

end module shr_mct_mod

