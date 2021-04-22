module interp_mod
  use cam_logfile, only : iulog
  use shr_kind_mod, only : r8 => shr_kind_r8
  use dimensions_mod, only : nelemd, np
  use interpolate_mod, only : interpolate_scalar, setup_latlon_interp, set_interp_parameter, get_interp_lat, get_interp_lon, &
       var_is_vector_uvar, var_is_vector_vvar, interpolate_vector, interpdata_t, get_interp_gweight
  use dyn_grid,       only : elem, w
  use spmd_utils,       only : masterproc, iam
  use cam_pio_utils,  only: phys_decomp, fillvalue
  use hybrid_mod,     only : hybrid_t, hybrid_create
  use abortutils, only: endrun

  implicit none
  private
  type(interpdata_t), pointer :: cam_interpolate(:)

  public get_interp_lat, get_interp_lon, setup_history_interpolation, write_interpolated
  public var_is_vector_uvar, var_is_vector_vvar, latlon_interpolation, add_interp_attributes

  interface write_interpolated
     module procedure write_interpolated_scalar
     module procedure write_interpolated_vector
  end interface
  type(hybrid_t) :: hybrid

contains

! Subprogram not used   subroutine add_interp_attributes(file)
! Subprogram not used     use pio, only : file_desc_t, pio_put_att, pio_global
! Subprogram not used     use interpolate_mod, only : get_interp_parameter
! Subprogram not used     type(file_desc_t) :: file
! Subprogram not used 
! Subprogram not used     integer :: ierr
! Subprogram not used     integer :: itmp
! Subprogram not used     
! Subprogram not used     itmp = get_interp_parameter('itype')
! Subprogram not used     if(itmp == 0) then
! Subprogram not used        ierr = pio_put_att(file, PIO_GLOBAL, 'interp_type', 'se basis functions')
! Subprogram not used     else if(itmp == 1) then
! Subprogram not used        ierr = pio_put_att(file, PIO_GLOBAL, 'interp_type', 'bilinear')
! Subprogram not used     else
! Subprogram not used        ierr = pio_put_att(file, PIO_GLOBAL, 'interp_type', itmp)
! Subprogram not used     end if
! Subprogram not used     
! Subprogram not used     itmp = get_interp_parameter('gridtype')
! Subprogram not used     select case(itmp)
! Subprogram not used     case(1)
! Subprogram not used        ierr = pio_put_att(file, PIO_GLOBAL, 'interp_outputgridtype', 'equally spaced with poles')
! Subprogram not used     case(2)        
! Subprogram not used        ierr = pio_put_att(file, PIO_GLOBAL, 'interp_outputgridtype', 'Gauss')
! Subprogram not used     case(3)
! Subprogram not used        ierr = pio_put_att(file, PIO_GLOBAL, 'interp_outputgridtype', 'equally spaced no poles')
! Subprogram not used     case default
! Subprogram not used        ierr = pio_put_att(file, PIO_GLOBAL, 'interp_outputgridtype', itmp)
! Subprogram not used     end select
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   end subroutine add_interp_attributes

  subroutine setup_history_interpolation(mtapes)

    use dyn_comp, only : dom_mt
    use parallel_mod,   only: par
    use thread_mod,     only: omp_get_thread_num
    use interpolate_mod, only : interpolate_analysis, get_interp_parameter
    implicit none
  
    integer, intent(in) :: mtapes
    integer :: ithr, nthreads

    if(iam>= par%nprocs) return

    ithr=omp_get_thread_num()
    hybrid = hybrid_create(par,ithr,1)
       
    if(any(interpolate_analysis)) then
       allocate(cam_interpolate(nelemd))
       call setup_latlon_interp(elem, cam_interpolate, par)
       allocate(w(get_interp_parameter('nlat')))
       w = get_interp_gweight()
    end if

  end subroutine setup_history_interpolation

  function latlon_interpolation(t)
    use interpolate_mod, only : interpolate_analysis
    integer, intent(in) :: t

    logical :: latlon_interpolation

    if (t<=size(interpolate_analysis)) then
       latlon_interpolation = interpolate_analysis(t)
    else
       latlon_interpolation = .false.
    endif

  end function latlon_interpolation



! Subprogram not used   subroutine write_interpolated_scalar(File, varid, fld, numlev, data_type, decomp_type) 
! Subprogram not used     use pio, only : file_desc_t, io_desc_t, var_desc_t, pio_write_darray, iosystem_desc_t, &
! Subprogram not used          pio_initdecomp, pio_freedecomp, pio_setdebuglevel
! Subprogram not used     use cam_instance, only: atm_id
! Subprogram not used     use spmd_dyn,       only: local_dp_map, block_buf_nrecs, chunk_buf_nrecs
! Subprogram not used     use ppgrid, only : begchunk, endchunk, pcols, pver
! Subprogram not used     use phys_grid, only : get_gcol_all_p, get_ncols_p, chunk_to_block_send_pters, chunk_to_block_recv_pters, &
! Subprogram not used          transpose_chunk_to_block
! Subprogram not used     use dyn_grid,       only: get_gcol_block_d
! Subprogram not used     use dimensions_mod, only: npsq
! Subprogram not used     use element_mod, only : element_t
! Subprogram not used     use dof_mod, only : PutUniquePoints
! Subprogram not used     use interpolate_mod, only : get_interp_parameter
! Subprogram not used     use shr_pio_mod, only : shr_pio_getiosys
! Subprogram not used     use edge_mod, only : edgebuffer_t, edgevpack, edgevunpack, initedgebuffer, freeedgebuffer
! Subprogram not used     use bndry_mod, only : bndry_exchangeV
! Subprogram not used     use parallel_mod,   only: par
! Subprogram not used     use abortutils, only : endrun
! Subprogram not used     
! Subprogram not used     implicit none
! Subprogram not used     type(file_desc_t), intent(inout) :: File
! Subprogram not used     type(var_desc_t), intent(inout) :: varid
! Subprogram not used     real(r8), intent(in) :: fld(:,:,:)
! Subprogram not used     integer, intent(in) :: numlev, data_type, decomp_type
! Subprogram not used 
! Subprogram not used     type(io_desc_t) :: iodesc
! Subprogram not used 
! Subprogram not used     integer :: lchnk, i, j, m, icol, ncols, pgcols(pcols), ierr
! Subprogram not used     integer :: idmb1(1), idmb2(1), idmb3(1)
! Subprogram not used     integer :: bpter(npsq,0:pver)    ! offsets into block buffer for packing data
! Subprogram not used     integer :: cpter(pcols,0:pver)   ! offsets into chunk buffer for unpacking data
! Subprogram not used 
! Subprogram not used     real(r8), pointer :: dest(:,:,:,:) 
! Subprogram not used     real(r8), pointer :: bbuffer(:), cbuffer(:), fldout(:,:)
! Subprogram not used     real(r8) :: fld_dyn(npsq,numlev,nelemd)
! Subprogram not used     integer :: st, en, ie, ioff, ncnt_out, k
! Subprogram not used     integer, pointer :: idof(:)
! Subprogram not used     integer :: nlon, nlat, ncol
! Subprogram not used     logical :: usefillvalues=.false.
! Subprogram not used     type(iosystem_desc_t), pointer :: pio_subsystem
! Subprogram not used     type (EdgeBuffer_t) :: edgebuf              ! edge buffer
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     nlon=get_interp_parameter('nlon')
! Subprogram not used     nlat=get_interp_parameter('nlat')
! Subprogram not used     pio_subsystem => shr_pio_getiosys(atm_id)
! Subprogram not used 
! Subprogram not used     if(decomp_type==phys_decomp) then
! Subprogram not used        fld_dyn = -999_R8
! Subprogram not used        if(local_dp_map) then
! Subprogram not used           !$omp parallel do private (lchnk, ncols, pgcols, icol, idmb1, idmb2, idmb3, ie, ioff, m)
! Subprogram not used           do lchnk=begchunk,endchunk
! Subprogram not used              ncols=get_ncols_p(lchnk)
! Subprogram not used              call get_gcol_all_p(lchnk,pcols,pgcols)
! Subprogram not used              do icol=1,ncols
! Subprogram not used                 call get_gcol_block_d(pgcols(icol),1,idmb1,idmb2,idmb3)
! Subprogram not used                 ie = idmb3(1)
! Subprogram not used                 ioff=idmb2(1)
! Subprogram not used                 do k=1,numlev
! Subprogram not used                    fld_dyn(ioff,k,ie)      = fld(icol, k, lchnk-begchunk+1)
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used              
! Subprogram not used           end do
! Subprogram not used        else
! Subprogram not used 
! Subprogram not used           allocate( bbuffer(block_buf_nrecs*numlev) )
! Subprogram not used           allocate( cbuffer(chunk_buf_nrecs*numlev) )
! Subprogram not used 
! Subprogram not used           !$omp parallel do private (lchnk, ncols, cpter, i, icol)
! Subprogram not used           do lchnk = begchunk,endchunk
! Subprogram not used              ncols = get_ncols_p(lchnk)
! Subprogram not used              
! Subprogram not used              call chunk_to_block_send_pters(lchnk,pcols,pver+1,1,cpter)
! Subprogram not used              
! Subprogram not used              do i=1,ncols
! Subprogram not used                 cbuffer(cpter(i,1):cpter(i,1)) = 0.0_r8
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used              do icol=1,ncols
! Subprogram not used                 
! Subprogram not used                 cbuffer   (cpter(icol,:))     = fld(icol,:,lchnk-begchunk+1)
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used           call transpose_chunk_to_block(1, cbuffer, bbuffer)
! Subprogram not used           if(iam < par%nprocs) then
! Subprogram not used !$omp parallel do private (ie, bpter, icol)
! Subprogram not used              do ie=1,nelemd
! Subprogram not used           
! Subprogram not used                 call chunk_to_block_recv_pters(elem(ie)%GlobalID,npsq,pver+1,1,bpter)
! Subprogram not used                 ncols = elem(ie)%idxp%NumUniquePts
! Subprogram not used                 do icol=1,ncols
! Subprogram not used                    fld_dyn   (icol,:,ie)   = bbuffer(bpter(icol,:))
! Subprogram not used                 end do
! Subprogram not used 
! Subprogram not used              end do
! Subprogram not used           end if
! Subprogram not used           deallocate( bbuffer )
! Subprogram not used           deallocate( cbuffer )
! Subprogram not used 
! Subprogram not used        end if
! Subprogram not used        allocate(dest(np,np,numlev,nelemd))
! Subprogram not used        call initEdgeBuffer(edgebuf, numlev)
! Subprogram not used 
! Subprogram not used        do ie=1,nelemd
! Subprogram not used           ncols = elem(ie)%idxp%NumUniquePts
! Subprogram not used           call putUniquePoints(elem(ie)%idxP, numlev, fld_dyn(1:ncols,:,ie), dest(:,:,:,ie))
! Subprogram not used           call edgeVpack(edgebuf, dest(:,:,:,ie), numlev, 0, elem(ie)%desc)
! Subprogram not used        enddo
! Subprogram not used        if(iam < par%nprocs) then
! Subprogram not used           call bndry_exchangeV(par, edgebuf)
! Subprogram not used        end if
! Subprogram not used        do ie=1,nelemd
! Subprogram not used           call edgeVunpack(edgebuf, dest(:,:,:,ie), numlev, 0, elem(ie)%desc)
! Subprogram not used        end do
! Subprogram not used        call freeEdgeBuffer(edgebuf)
! Subprogram not used        usefillvalues = any(dest == fillvalue)
! Subprogram not used     else
! Subprogram not used         usefillvalues=any(fld==fillvalue)
! Subprogram not used        allocate(dest(np,np,numlev,1))
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ncnt_out = sum(cam_interpolate(1:nelemd)%n_interp)
! Subprogram not used     allocate(fldout(ncnt_out,numlev))
! Subprogram not used     allocate(idof(ncnt_out*numlev))
! Subprogram not used     fldout = -999_r8
! Subprogram not used     idof = 0
! Subprogram not used     st = 1
! Subprogram not used     
! Subprogram not used     
! Subprogram not used 
! Subprogram not used     do ie=1,nelemd
! Subprogram not used        ncol = cam_interpolate(ie)%n_interp
! Subprogram not used        do k=0,numlev-1
! Subprogram not used           do i=1,ncol
! Subprogram not used              idof(st+i-1+k*ncnt_out)=cam_interpolate(ie)%ilon(i)+nlon*(cam_interpolate(ie)%ilat(i)-1)+nlon*nlat*k
! Subprogram not used           enddo
! Subprogram not used        enddo
! Subprogram not used        
! Subprogram not used        ! Now that we have the field on the dyn grid we need to interpolate
! Subprogram not used        en = st+cam_interpolate(ie)%n_interp-1
! Subprogram not used        if(decomp_type==phys_decomp) then
! Subprogram not used           if(usefillvalues) then
! Subprogram not used              call interpolate_scalar(cam_interpolate(ie),dest(:,:,:,ie), np, numlev, fldout(st:en,:), fillvalue) 
! Subprogram not used           else
! Subprogram not used              call interpolate_scalar(cam_interpolate(ie),dest(:,:,:,ie), np, numlev, fldout(st:en,:)) 
! Subprogram not used           end if
! Subprogram not used        else
! Subprogram not used           do j=1,np
! Subprogram not used              do i=1,np
! Subprogram not used                 dest(i,j,:,1) = fld(i+(j-1)*np,:,ie)
! Subprogram not used              end do
! Subprogram not used           end do
! Subprogram not used           if(usefillvalues) then
! Subprogram not used              call interpolate_scalar(cam_interpolate(ie),dest(:,:,:,1), &
! Subprogram not used                   np, numlev, fldout(st:en,:), fillvalue) 
! Subprogram not used           else
! Subprogram not used              call interpolate_scalar(cam_interpolate(ie),dest(:,:,:,1), &
! Subprogram not used                   np, numlev, fldout(st:en,:)) 
! Subprogram not used           end if
! Subprogram not used        end if
! Subprogram not used 
! Subprogram not used 
! Subprogram not used        st = en+1
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     if(numlev==1) then
! Subprogram not used        call pio_initdecomp(pio_subsystem, data_type, (/nlon,nlat/), idof, iodesc)
! Subprogram not used     else
! Subprogram not used        call pio_initdecomp(pio_subsystem, data_type, (/nlon,nlat,numlev/), idof, iodesc)
! Subprogram not used     end if
! Subprogram not used     call pio_write_darray(File, varid, iodesc, fldout, ierr)
! Subprogram not used 
! Subprogram not used     deallocate(dest)
! Subprogram not used 
! Subprogram not used     deallocate(fldout)
! Subprogram not used     deallocate(idof)
! Subprogram not used     call pio_freedecomp(file,iodesc)
! Subprogram not used 
! Subprogram not used   end subroutine write_interpolated_scalar




! Subprogram not used   subroutine write_interpolated_vector(File, varidu, varidv, fldu, fldv, numlev, data_type, decomp_type)
! Subprogram not used     use pio, only : file_desc_t, io_desc_t, var_desc_t, pio_write_darray, iosystem_desc_t, &
! Subprogram not used          pio_initdecomp, pio_freedecomp, pio_setdebuglevel
! Subprogram not used     use spmd_dyn,       only: local_dp_map, block_buf_nrecs, chunk_buf_nrecs
! Subprogram not used     use ppgrid, only : begchunk, endchunk, pcols, pver
! Subprogram not used     use phys_grid, only : get_gcol_all_p, get_ncols_p, chunk_to_block_send_pters, chunk_to_block_recv_pters, &
! Subprogram not used          transpose_chunk_to_block
! Subprogram not used     use dyn_grid,       only: get_gcol_block_d
! Subprogram not used     use dimensions_mod, only: npsq
! Subprogram not used     use element_mod, only : element_t
! Subprogram not used     use dof_mod, only : PutUniquePoints
! Subprogram not used     use interpolate_mod, only : get_interp_parameter
! Subprogram not used     use shr_pio_mod, only : shr_pio_getiosys
! Subprogram not used     use edge_mod, only : edgebuffer_t, edgevpack, edgevunpack, initedgebuffer, freeedgebuffer
! Subprogram not used     use bndry_mod, only : bndry_exchangeV
! Subprogram not used     use parallel_mod,   only: par
! Subprogram not used     implicit none
! Subprogram not used     type(file_desc_t), intent(inout) :: File
! Subprogram not used     type(var_desc_t), intent(inout) :: varidu, varidv
! Subprogram not used     real(r8), intent(in) :: fldu(:,:,:), fldv(:,:,:)
! Subprogram not used     integer, intent(in) :: numlev, data_type, decomp_type
! Subprogram not used 
! Subprogram not used     type(io_desc_t) :: iodesc
! Subprogram not used 
! Subprogram not used     integer :: lchnk, i, j, m, icol, ncols, pgcols(pcols), ierr
! Subprogram not used     integer :: idmb1(1), idmb2(1), idmb3(1)
! Subprogram not used     integer :: bpter(npsq,0:pver)    ! offsets into block buffer for packing data
! Subprogram not used     integer :: cpter(pcols,0:pver)   ! offsets into chunk buffer for unpacking data
! Subprogram not used 
! Subprogram not used     real(r8), allocatable :: dest(:,:,:,:,:)
! Subprogram not used     real(r8), pointer :: bbuffer(:), cbuffer(:), fldout(:,:,:)
! Subprogram not used     real(r8) :: fld_dyn(npsq,2,numlev,nelemd)
! Subprogram not used     integer :: st, en, ie, ioff, ncnt_out, k
! Subprogram not used     integer, pointer :: idof(:)
! Subprogram not used     integer :: nlon, nlat, ncol
! Subprogram not used     logical :: usefillvalues=.false.
! Subprogram not used 
! Subprogram not used     type(iosystem_desc_t), pointer :: pio_subsystem
! Subprogram not used     type (EdgeBuffer_t) :: edgebuf              ! edge buffer
! Subprogram not used 
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     nlon=get_interp_parameter('nlon')
! Subprogram not used     nlat=get_interp_parameter('nlat')
! Subprogram not used     pio_subsystem => shr_pio_getiosys('ATM')
! Subprogram not used     fld_dyn = -999_R8
! Subprogram not used     if(decomp_type==phys_decomp) then
! Subprogram not used        allocate(dest(np,np,2,numlev,nelemd))
! Subprogram not used        if(local_dp_map) then
! Subprogram not used           !$omp parallel do private (lchnk, ncols, pgcols, icol, idmb1, idmb2, idmb3, ie, ioff, m)
! Subprogram not used           do lchnk=begchunk,endchunk
! Subprogram not used              ncols=get_ncols_p(lchnk)
! Subprogram not used              call get_gcol_all_p(lchnk,pcols,pgcols)
! Subprogram not used              
! Subprogram not used              do icol=1,ncols
! Subprogram not used                 call get_gcol_block_d(pgcols(icol),1,idmb1,idmb2,idmb3)
! Subprogram not used                 ie = idmb3(1)
! Subprogram not used                 ioff=idmb2(1)
! Subprogram not used                 do k=1,numlev
! Subprogram not used                    fld_dyn(ioff,1,k,ie)      = fldu(icol, k, lchnk-begchunk+1)
! Subprogram not used                    fld_dyn(ioff,2,k,ie)      = fldv(icol, k, lchnk-begchunk+1)
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used              
! Subprogram not used           end do
! Subprogram not used        else
! Subprogram not used 
! Subprogram not used           allocate( bbuffer(2*block_buf_nrecs*numlev) )
! Subprogram not used           allocate( cbuffer(2*chunk_buf_nrecs*numlev) )
! Subprogram not used 
! Subprogram not used           !$omp parallel do private (lchnk, ncols, cpter, i, icol)
! Subprogram not used           do lchnk = begchunk,endchunk
! Subprogram not used              ncols = get_ncols_p(lchnk)
! Subprogram not used              
! Subprogram not used              call chunk_to_block_send_pters(lchnk,pcols,pver+1,2,cpter)
! Subprogram not used              
! Subprogram not used              do i=1,ncols
! Subprogram not used                 cbuffer(cpter(i,1):cpter(i,1)) = 0.0_r8
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used              do icol=1,ncols
! Subprogram not used                 do k=1,numlev
! Subprogram not used                    cbuffer   (cpter(icol,k))     = fldu(icol,k,lchnk-begchunk+1)
! Subprogram not used                    cbuffer   (cpter(icol,k)+1)   = fldv(icol,k,lchnk-begchunk+1)
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used 
! Subprogram not used           end do
! Subprogram not used 
! Subprogram not used           call transpose_chunk_to_block(2, cbuffer, bbuffer)
! Subprogram not used           if(iam < par%nprocs) then
! Subprogram not used              !$omp parallel do private (ie, bpter, icol)
! Subprogram not used              do ie=1,nelemd
! Subprogram not used                 
! Subprogram not used                 call chunk_to_block_recv_pters(elem(ie)%GlobalID,npsq,pver+1,2,bpter)
! Subprogram not used                 ncols = elem(ie)%idxp%NumUniquePts
! Subprogram not used                 do icol=1,ncols
! Subprogram not used                    do k=1,numlev
! Subprogram not used                       fld_dyn   (icol,1,k,ie)   = bbuffer(bpter(icol,k))
! Subprogram not used                       fld_dyn   (icol,2,k,ie)   = bbuffer(bpter(icol,k)+1)
! Subprogram not used                    enddo
! Subprogram not used                 end do
! Subprogram not used                 
! Subprogram not used              end do
! Subprogram not used           end if
! Subprogram not used           deallocate( bbuffer )
! Subprogram not used           deallocate( cbuffer )
! Subprogram not used 
! Subprogram not used        end if
! Subprogram not used        call initEdgeBuffer(edgebuf, 2*numlev)
! Subprogram not used 
! Subprogram not used        do ie=1,nelemd
! Subprogram not used           ncols = elem(ie)%idxp%NumUniquePts
! Subprogram not used           call putUniquePoints(elem(ie)%idxP, 2, numlev, fld_dyn(1:ncols,:,:,ie), dest(:,:,:,:,ie))
! Subprogram not used           
! Subprogram not used           call edgeVpack(edgebuf, dest(:,:,:,:,ie), 2*numlev, 0, elem(ie)%desc)
! Subprogram not used        enddo
! Subprogram not used        if(iam < par%nprocs) then
! Subprogram not used           call bndry_exchangeV(par, edgebuf)
! Subprogram not used        end if
! Subprogram not used 
! Subprogram not used        do ie=1,nelemd
! Subprogram not used           call edgeVunpack(edgebuf, dest(:,:,:,:,ie), 2*numlev, 0, elem(ie)%desc)
! Subprogram not used        enddo
! Subprogram not used        call freeEdgeBuffer(edgebuf)
! Subprogram not used        usefillvalues = any(dest==fillvalue)
! Subprogram not used     else
! Subprogram not used        usefillvalues = (any(fldu==fillvalue) .or. any(fldv==fillvalue))
! Subprogram not used        allocate(dest(np,np,2,numlev,1))
! Subprogram not used     endif
! Subprogram not used     ncnt_out = sum(cam_interpolate(1:nelemd)%n_interp)
! Subprogram not used     allocate(fldout(ncnt_out,numlev,2))
! Subprogram not used     allocate(idof(ncnt_out*numlev))
! Subprogram not used     
! Subprogram not used     fldout = -999_r8
! Subprogram not used     idof = 0
! Subprogram not used     st = 1
! Subprogram not used     do ie=1,nelemd
! Subprogram not used        ncol = cam_interpolate(ie)%n_interp
! Subprogram not used        do k=0,numlev-1
! Subprogram not used           do i=1,ncol
! Subprogram not used              idof(st+i-1+k*ncnt_out)=cam_interpolate(ie)%ilon(i)+nlon*(cam_interpolate(ie)%ilat(i)-1)+nlon*nlat*k
! Subprogram not used           enddo
! Subprogram not used        enddo
! Subprogram not used        
! Subprogram not used 
! Subprogram not used     ! Now that we have the field on the dyn grid we need to interpolate
! Subprogram not used        en = st+cam_interpolate(ie)%n_interp-1
! Subprogram not used        if(decomp_type==phys_decomp) then
! Subprogram not used           if(usefillvalues) then
! Subprogram not used              call interpolate_vector(cam_interpolate(ie),elem(ie), &
! Subprogram not used                   dest(:,:,:,:,ie), np, numlev, fldout(st:en,:,:), 0, fillvalue) 
! Subprogram not used           else
! Subprogram not used              call interpolate_vector(cam_interpolate(ie),elem(ie),&
! Subprogram not used                   dest(:,:,:,:,ie), np, numlev, fldout(st:en,:,:), 0) 
! Subprogram not used           endif
! Subprogram not used        else
! Subprogram not used           do k=1,numlev
! Subprogram not used              do j=1,np
! Subprogram not used                 do i=1,np
! Subprogram not used                    dest(i,j,1,k,1) = fldu(i+(j-1)*np,k,ie)
! Subprogram not used                    dest(i,j,1,k,1) = fldv(i+(j-1)*np,k,ie)
! Subprogram not used                 end do
! Subprogram not used              end do
! Subprogram not used           end do
! Subprogram not used           if(usefillvalues) then
! Subprogram not used              call interpolate_vector(cam_interpolate(ie),elem(ie),&
! Subprogram not used                   dest(:,:,:,:,1), np, numlev, fldout(st:en,:,:), 0, fillvalue) 
! Subprogram not used           else
! Subprogram not used              call interpolate_vector(cam_interpolate(ie),elem(ie),&
! Subprogram not used                   dest(:,:,:,:,1), np, numlev, fldout(st:en,:,:), 0) 
! Subprogram not used           end if
! Subprogram not used        end if
! Subprogram not used 
! Subprogram not used        st = en+1
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     if(numlev==1) then
! Subprogram not used        call pio_initdecomp(pio_subsystem, data_type, (/nlon,nlat/), idof, iodesc)
! Subprogram not used     else
! Subprogram not used        call pio_initdecomp(pio_subsystem, data_type, (/nlon,nlat,numlev/), idof, iodesc)
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     call pio_write_darray(File, varidu, iodesc, fldout(:,:,1), ierr)
! Subprogram not used 
! Subprogram not used     call pio_write_darray(File, varidv, iodesc, fldout(:,:,2), ierr)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     deallocate(fldout)
! Subprogram not used     deallocate(idof)
! Subprogram not used     deallocate(dest)
! Subprogram not used     call pio_freedecomp(file,iodesc)
! Subprogram not used 
! Subprogram not used   end subroutine write_interpolated_vector












end module interp_mod

