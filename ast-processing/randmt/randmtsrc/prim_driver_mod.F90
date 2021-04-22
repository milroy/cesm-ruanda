!#define _DBG_ print *,"file: ","<stdin>"," line: ",5," ithr: ",hybrid%ithr
module prim_driver_mod
  use kinds, only : real_kind, iulog, longdouble_kind
  use dimensions_mod, only : np, nlev, nlevp, nelem, nelemd, nelemdmax, GlobalUniqueCols, ntrac, qsize, nc,nhc, nep, nipm
  use cg_mod, only : cg_t
  use hybrid_mod, only : hybrid_t
  use quadrature_mod, only : quadrature_t, test_gauss, test_gausslobatto, gausslobatto
  use prim_si_ref_mod, only : ref_state_t
  use solver_mod, only : blkjac_t
  use filter_mod, only : filter_t
  use derivative_mod, only : derivative_t
  use reduction_mod, only : reductionbuffer_ordered_1d_t, red_min, red_max, &
         red_sum, red_sum_int, red_flops, initreductionbuffer

  use fvm_mod, only : fvm_init1,fvm_init2, fvm_init3  
  use fvm_control_volume_mod, only : fvm_struct
  use spelt_mod, only : spelt_struct, spelt_init1,spelt_init2, spelt_init3
  
  use element_mod, only : element_t, timelevels,  allocate_element_desc

  implicit none
  private
  public :: prim_init1, prim_init2 , prim_run, prim_run_subcycle, prim_finalize, leapfrog_bootstrap
  public :: smooth_topo_datasets

  type (cg_t), allocatable  :: cg(:)              ! conjugate gradient struct (nthreads)
  type (quadrature_t)   :: gp                     ! element GLL points
  real(kind=longdouble_kind)  :: fvm_corners(nc+1)     ! fvm cell corners on reference element
  real(kind=longdouble_kind)  :: fvm_points(nc)     ! fvm cell centers on reference element
  real (kind=longdouble_kind) :: spelt_refnep(1:nep)


  type (ref_state_t)    :: refstate        ! semi-implicit vertical reference state
  type (blkjac_t),allocatable  :: blkjac(:)  ! (nets:nete)
  type (filter_t)       :: flt             ! Filter struct for v and p grid
  type (filter_t)       :: flt_advection   ! Filter struct for v grid for advection only
  type (derivative_t), allocatable   :: deriv(:) ! derivative struct (nthreads)
  real*8  :: tot_iter
  type (ReductionBuffer_ordered_1d_t), save :: red   ! reduction buffer               (shared)

contains

  subroutine prim_init1(elem, fvm, par, dom_mt, Tl)

    ! --------------------------------
    use thread_mod, only : nthreads, omp_get_thread_num, omp_set_num_threads
    ! --------------------------------
    use control_mod, only : runtype, restartfreq, filter_counter, integration, topology, &
         partmethod, while_iter
    ! --------------------------------
    use prim_state_mod, only : prim_printstate_init
    ! --------------------------------
    use namelist_mod, only : readnl
    ! --------------------------------
    use mesh_mod, only : MeshUseMeshFile
    ! -------------------------------- 
    use time_mod, only : nmax, time_at, timelevel_init, timelevel_t
    ! --------------------------------
    ! --------------------------------
    use mass_matrix_mod, only : mass_matrix
    ! --------------------------------
    use cube_mod,  only : cubeedgecount , cubeelemcount, cubetopology
    ! --------------------------------
    use mesh_mod, only : MeshSetCoordinates, MeshUseMeshFile, MeshCubeTopology, &
         MeshCubeElemCount, MeshCubeEdgeCount
    use cube_mod, only : cube_init_atomic, rotation_init_atomic, set_corner_coordinates, assign_node_numbers_to_elem
    ! --------------------------------
    use metagraph_mod, only : metavertex_t, metaedge_t, localelemcount, initmetagraph
    ! --------------------------------
    use gridgraph_mod, only : gridvertex_t, gridedge_t, allocate_gridvertex_nbrs, deallocate_gridvertex_nbrs
    ! --------------------------------
    use schedule_mod, only : schedule, genEdgeSched,  PrintSchedule
    ! --------------------------------
    use prim_advection_mod, only: prim_advec_init
    ! --------------------------------    
    use prim_advance_mod, only: prim_advance_init
    ! --------------------------------    
    use diffusion_mod, only      : diffusion_init
    ! --------------------------------    
    use parallel_mod, only : iam, parallel_t, syncmp, abortmp, global_shared_buf, nrepro_vars
    use parallel_mod, only : mpiinteger_t, mpireal_t, mpi_max, mpi_sum, haltmp
    ! --------------------------------
    use metis_mod, only : genmetispart
    ! --------------------------------
    use spacecurve_mod, only : genspacepart
    ! --------------------------------
    use dof_mod, only : global_dof, CreateUniqueIndex, SetElemOffset
    ! --------------------------------
    use params_mod, only : SFCURVE
    ! --------------------------------
    use domain_mod, only : domain1d_t, decompose
    ! --------------------------------
    use physical_constants, only : dd_pi
    ! --------------------------------
    use shr_reprosum_mod, only: repro_sum => shr_reprosum_calc
    implicit none
    type (element_t), pointer :: elem(:)
     type (fvm_struct), pointer   :: fvm(:)    
    type (parallel_t), intent(in) :: par
    type (domain1d_t), pointer :: dom_mt(:)      
    type (timelevel_t), intent(out) :: Tl
    ! Local Variables

    type (GridVertex_t), target,allocatable :: GridVertex(:)
    type (GridEdge_t),   target,allocatable :: Gridedge(:)
    type (MetaVertex_t), target,allocatable :: MetaVertex(:)

    integer :: ii,ie, ith
    integer :: nets, nete
    integer :: nelem_edge,nedge
    integer :: nstep
    integer :: nlyr
    integer :: iMv
    integer :: err, ierr, l, j

    real(kind=real_kind), allocatable :: aratio(:,:)
    real(kind=real_kind) :: area(1),xtmp
    character(len=80) rot_type   ! cube edge rotation type

    integer  :: i
    integer,allocatable :: TailPartition(:)
    integer,allocatable :: HeadPartition(:)

    integer total_nelem
    real(kind=real_kind) :: approx_elements_per_task
    integer :: n_domains


    ! =====================================
    ! Read in model control information
    ! =====================================
    ! cam readnl is called in spmd_dyn (needed prior to mpi_init)
    ! ====================================
    ! Set cube edge rotation type for model
    ! unnecessary complication here: all should
    ! be on the same footing. RDL
    ! =====================================
    rot_type="contravariant"



    ! ===============================================================
    ! Allocate and initialize the graph (array of GridVertex_t types)
    ! ===============================================================

    if (topology=="cube") then

       if (par%masterproc) then
          write(iulog,*)"creating cube topology..."
       end if

       if (MeshUseMeshFile) then
           nelem = MeshCubeElemCount()
           nelem_edge = MeshCubeEdgeCount()
       else
           nelem      = CubeElemCount()
           nelem_edge = CubeEdgeCount()
       end if

       allocate(GridVertex(nelem))
       allocate(GridEdge(nelem_edge))

       do j =1,nelem
          call allocate_gridvertex_nbrs(GridVertex(j))
       end do

       if (MeshUseMeshFile) then
           if (par%masterproc) then
               write(iulog,*) "Set up grid vertex from mesh..."
           end if
           call MeshCubeTopology(GridEdge, GridVertex)
       else
           call CubeTopology(GridEdge,GridVertex)
        end if
       
       if(par%masterproc)       write(iulog,*)"...done."
    end if


    !debug  call PrintGridVertex(GridVertex)


    if(partmethod .eq. SFCURVE) then 
       if(par%masterproc) write(iulog,*)"partitioning graph using SF Curve..."
       call genspacepart(GridEdge,GridVertex)
    else
        if(par%masterproc) write(iulog,*)"partitioning graph using Metis..."
       call genmetispart(GridEdge,GridVertex)
    endif

    ! ===========================================================
    ! given partition, count number of local element descriptors
    ! ===========================================================
    allocate(MetaVertex(1))
    allocate(Schedule(1))


    nelem_edge=SIZE(GridEdge)

    allocate(TailPartition(nelem_edge))
    allocate(HeadPartition(nelem_edge))
    do i=1,nelem_edge
       TailPartition(i)=GridEdge(i)%tail%processor_number
       HeadPartition(i)=GridEdge(i)%head%processor_number
    enddo

    ! ====================================================
    !  Generate the communication graph
    ! ====================================================
    call initMetaGraph(iam,MetaVertex(1),GridVertex,GridEdge)


    nelemd = LocalElemCount(MetaVertex(1))

    if(nelemd .le. 0) then 
       call abortmp('Not yet ready to handle nelemd = 0 yet' )
       stop
    endif
    call mpi_allreduce(nelemd,nelemdmax,1,MPIinteger_t,MPI_MAX,par%comm,ierr)


    if (nelemd>0) then
       allocate(elem(nelemd))
       call allocate_element_desc(elem)

    endif

    if (ntrac>0) then
       allocate(fvm(nelemd))
    else
       ! Even if fvm not needed, still desirable to allocate it as empty
       ! so it can be passed as a (size zero) array rather than pointer.
       allocate(fvm(0))
    end if

    ! ====================================================
    !  Generate the communication schedule
    ! ====================================================

    call genEdgeSched(elem,iam,Schedule(1),MetaVertex(1))


    allocate(global_shared_buf(nelemd,nrepro_vars))
    !  nlyr=edge3p1%nlyr
    !  call MessageStats(nlyr)
    !  call testchecksum(par,GridEdge)

    ! ========================================================
    ! load graph information into local element descriptors
    ! ========================================================

    !  do ii=1,nelemd
    !     elem(ii)%vertex = MetaVertex(iam)%members(ii)
    !  enddo

    call syncmp(par)

    ! =================================================================
    ! Set number of domains (for 'decompose') equal to number of threads
    !  for OpenMP across elements, equal to 1 for OpenMP within element
    ! =================================================================
    n_domains = min(Nthreads,nelemd)


    ! =================================================================
    ! Initialize shared boundary_exchange and reduction buffers
    ! =================================================================
    if(par%masterproc) write(iulog,*) 'init shared boundary_exchange buffers'
    call InitReductionBuffer(red,3*nlev,n_domains)
    call InitReductionBuffer(red_sum,5)
    call InitReductionBuffer(red_sum_int,1)
    call InitReductionBuffer(red_max,1)
    call InitReductionBuffer(red_min,1)
    call initReductionBuffer(red_flops,1)


    gp=gausslobatto(np)  ! GLL points

    ! fvm nodes are equally spaced in alpha/beta
    ! HOMME with equ-angular gnomonic projection maps alpha/beta space
    ! to the reference element via simple scale + translation
    ! thus, fvm nodes in reference element [-1,1] are a tensor product of
    ! array 'fvm_corners(:)' computed below:
    xtmp=nc 
    do i=1,nc+1
       fvm_corners(i)= 2*(i-1)/xtmp - 1  ! [-1,1] including end points
    end do
    do i=1,nc
       fvm_points(i)= ( fvm_corners(i)+fvm_corners(i+1) ) /2
    end do
     
    xtmp=nep-1
    do i=1,nep
      spelt_refnep(i)= 2*(i-1)/xtmp - 1
    end do
     
    if (topology=="cube") then
       if(par%masterproc) write(iulog,*) "initializing cube elements..."
       if (MeshUseMeshFile) then
           call MeshSetCoordinates(elem)
       else
           do ie=1,nelemd
               call set_corner_coordinates(elem(ie))
           end do
           call assign_node_numbers_to_elem(elem, GridVertex)
       end if
       do ie=1,nelemd
          call cube_init_atomic(elem(ie),gp%points)
       enddo
    end if

    ! =================================================================
    ! Initialize mass_matrix
    ! =================================================================
    if(par%masterproc) write(iulog,*) 'running mass_matrix'
    call mass_matrix(par,elem)
    allocate(aratio(nelemd,1))

    if (topology=="cube") then
       area = 0
       do ie=1,nelemd
          aratio(ie,1) = sum(elem(ie)%mp(:,:)*elem(ie)%metdet(:,:))
       enddo
       call repro_sum(aratio, area, nelemd, nelemd, 1, commid=par%comm)
       area(1) = 4*dd_pi/area(1)  ! ratio correction
       deallocate(aratio)
       if (par%masterproc) &
            write(iulog,'(a,f20.17)') " re-initializing cube elements: area correction=",area(1)

       do ie=1,nelemd
          call cube_init_atomic(elem(ie),gp%points,area(1))
          call rotation_init_atomic(elem(ie),rot_type)
       enddo
    end if


    if(par%masterproc) write(iulog,*) 're-running mass_matrix'
    call mass_matrix(par,elem)


    ! =================================================================
    ! Determine the global degree of freedome for each gridpoint
    ! =================================================================
    if(par%masterproc) write(iulog,*) 'running global_dof'
    call global_dof(par,elem)

    ! =================================================================
    ! Create Unique Indices
    ! =================================================================

    do ie=1,nelemd
       call CreateUniqueIndex(elem(ie)%GlobalId,elem(ie)%gdofP,elem(ie)%idxP)
    enddo

    call SetElemOffset(par,elem, GlobalUniqueCols)

    do ie=1,nelemd
       elem(ie)%idxV=>elem(ie)%idxP
    end do

    !JMD call PrintDofP(elem)
    !JMD call PrintDofV(elem)



    call prim_printstate_init(par)
    ! Initialize output fields for plotting...


    while_iter = 0
    filter_counter = 0

    ! initialize flux terms to 0

    do ie=1,nelemd
       elem(ie)%derived%FM=0.0
       elem(ie)%derived%FQ=0.0
       elem(ie)%derived%FQps=0.0
       elem(ie)%derived%FT=0.0
       elem(ie)%derived%pecnd=0.0

       elem(ie)%accum%Qvar=0
       elem(ie)%accum%Qmass=0
       elem(ie)%accum%Q1mass=0

       elem(ie)%derived%Omega_p=0
       elem(ie)%state%dp3d=0
    enddo


    ! ==========================================================
    !  This routines initalizes a Restart file.  This involves:
    !      I)  Setting up the MPI datastructures
    ! ==========================================================
    !DBG  write(iulog,*) 'prim_init: after call to initRestartFile'

    deallocate(GridEdge)
    do j =1,nelem
       call deallocate_gridvertex_nbrs(GridVertex(j))
    end do
    deallocate(GridVertex)

    do j = 1, MetaVertex(1)%nmembers
       call deallocate_gridvertex_nbrs(MetaVertex(1)%members(j))
    end do
    deallocate(MetaVertex)
    deallocate(TailPartition)
    deallocate(HeadPartition)

    n_domains = min(Nthreads,nelemd)
    call omp_set_num_threads(n_domains)
    ! =====================================
    ! Set number of threads...
    ! =====================================
    if(par%masterproc) then
       write(iulog,*) "Main:NThreads=",NThreads
       write(iulog,*) "Main:n_domains = ",n_domains
    endif

    allocate(dom_mt(0:n_domains-1))
    do ith=0,n_domains-1
       dom_mt(ith)=decompose(1,nelemd,n_domains,ith)
    end do
    ith=0
    nets=1
    nete=nelemd
    allocate(deriv(0:n_domains-1))
    allocate(cg(0:n_domains-1))
    call prim_advance_init(integration)
    call Prim_Advec_Init()
    call diffusion_init()
    if (ntrac>0) then
      call fvm_init1(par)    
    endif
    call TimeLevel_init(tl)
    if(par%masterproc) write(iulog,*) 'end of prim_init'
  end subroutine prim_init1
!=======================================================================================================! 

  subroutine prim_init2(elem, fvm, hybrid, nets, nete, tl, hvcoord)

    use parallel_mod, only : parallel_t, haltmp, syncmp, abortmp
    use time_mod, only : timelevel_t, tstep, phys_tscale, timelevel_init, nendstep, smooth, nsplit, TimeLevel_Qdp
    use prim_state_mod, only : prim_printstate, prim_diag_scalars
    use filter_mod, only : filter_t, fm_filter_create, taylor_filter_create, &
         fm_transfer, bv_transfer
    use control_mod, only : runtype, integration, filter_mu, filter_mu_advection, test_case, &
         debug_level, vfile_int, filter_freq, filter_freq_advection, &
         transfer_type, vform, vfile_mid, filter_type, kcut_fm, wght_fm, p_bv, &
         s_bv, topology,columnpackage, moisture, precon_method, rsplit, qsplit, rk_stage_user,&
         sub_case, &
         limiter_option, nu, nu_q, nu_div, tstep_type, hypervis_subcycle, &
         hypervis_subcycle_q
    use prim_si_ref_mod, only: prim_si_refstate_init, prim_set_mass
    use thread_mod, only : nthreads
    use derivative_mod, only : derivinit, interpolate_gll2fvm_points, interpolate_gll2spelt_points, v2pinit
    use global_norms_mod, only : test_global_integral, print_cfl
    use hybvcoord_mod, only : hvcoord_t

    type (element_t), intent(inout) :: elem(:)
     type (fvm_struct), intent(inout)    :: fvm(:)    
    type (hybrid_t), intent(in) :: hybrid

    type (TimeLevel_t), intent(inout)    :: tl              ! time level struct
    type (hvcoord_t), intent(inout)      :: hvcoord         ! hybrid vertical coordinate struct

     integer, intent(in)                     :: nets  ! starting thread element number (private)
    integer, intent(in)                     :: nete  ! ending thread element number   (private)


    ! ==================================
    ! Local variables
    ! ==================================

    real (kind=real_kind) :: dt              ! "timestep dependent" timestep
!   variables used to calculate CFL
    real (kind=real_kind) :: dtnu            ! timestep*viscosity parameter
    real (kind=real_kind) :: dt_dyn_vis      ! viscosity timestep used in dynamics
    real (kind=real_kind) :: dt_tracer_vis      ! viscosity timestep used in tracers

    real (kind=real_kind) :: dp        


    real (kind=real_kind) :: ps(np,np)       ! surface pressure

    character(len=80)     :: fname
    character(len=8)      :: njusn
    character(len=4)      :: charnum

    real (kind=real_kind) :: Tp(np)     ! transfer function 

    integer :: simday
    integer :: i,j,k,ie,iptr,t,q
    integer :: ierr
    integer :: nfrc
    integer :: n0_qdp    


    ! ==========================
    ! begin executable code
    ! ==========================
    if (topology == "cube") then
       call test_global_integral(elem, hybrid,nets,nete)  
    end if


    ! compute most restrictive dt*nu for use by variable res viscosity: 
    if (tstep_type == 0) then
       ! LF case: no tracers, timestep seen by viscosity is 2*tstep
       dt_tracer_vis = 0  
       dt_dyn_vis = 2*tstep
       dtnu = 2.0d0*tstep*max(nu,nu_div)
    else
       ! compute timestep seen by viscosity operator:
       dt_dyn_vis = tstep
       if (qsplit>1 .and. tstep_type == 1) then
          ! tstep_type==1: RK2 followed by LF.  internal LF stages apply viscosity at 2*dt
          dt_dyn_vis = 2*tstep  
       endif
       dt_tracer_vis=tstep*qsplit

       ! compute most restrictive condition: 
       ! note: dtnu ignores subcycling 
       dtnu=max(dt_dyn_vis*max(nu,nu_div), dt_tracer_vis*nu_q)
       ! compute actual viscosity timesteps with subcycling
       dt_tracer_vis = dt_tracer_vis/hypervis_subcycle_q
       dt_dyn_vis = dt_dyn_vis/hypervis_subcycle
    endif


    ! ==================================
    ! Initialize derivative structure
    ! ==================================
    call derivinit(deriv(hybrid%ithr),fvm_corners, fvm_points, spelt_refnep)
    ! ================================================
    ! fvm initialization
    ! ================================================
    if (ntrac>0) then
      call fvm_init2(elem,fvm,hybrid,nets,nete,tl)    
    endif
    ! ====================================
    ! In the semi-implicit case:
    ! initialize vertical structure and 
    ! related matrices..
    ! ====================================
!$OMP MASTER
    if (integration == "semi_imp") then
       refstate = prim_si_refstate_init(.false.,hybrid%masterthread,hvcoord)
       if (precon_method == "block_jacobi") then
          allocate(blkjac(nets:nete))
       endif
    endif
!$OMP END MASTER
    ! ==========================================
    ! Initialize pressure and velocity grid 
    ! filter matrix...
    ! ==========================================
    if (transfer_type == "bv") then
       Tp    = bv_transfer(p_bv,s_bv,np)
    else if (transfer_type == "fm") then
       Tp    = fm_transfer(kcut_fm,wght_fm,np)
    end if
    if (filter_type == "taylor") then
       flt           = taylor_filter_create(Tp, filter_mu,gp)
       flt_advection = taylor_filter_create(Tp, filter_mu_advection,gp)
    else if (filter_type == "fischer") then
       flt           = fm_filter_create(Tp, filter_mu, gp)
       flt_advection = fm_filter_create(Tp, filter_mu_advection, gp)
    end if



    if (hybrid%masterthread) then
       if (filter_freq>0 .or. filter_freq_advection>0) then
          write(iulog,*) "transfer function type in preq=",transfer_type
          write(iulog,*) "filter type            in preq=",filter_type
          write(*,'(a,99f10.6)') "dynamics: I-mu + mu*Tp(:) = ",&
               (1-filter_mu)+filter_mu*Tp(:)
          write(*,'(a,99f10.6)') "advection: I-mu + mu*Tp(:) = ",&
               (1-filter_mu_advection)+filter_mu_advection*Tp(:)
       endif
    endif

    !$OMP BARRIER
    if (hybrid%ithr==0) then
       call syncmp(hybrid%par)
    end if
    !$OMP BARRIER

    if (topology /= "cube") then
       call abortmp('Error: only cube topology supported for primaitve equations') 
    endif


    ! For new runs, and branch runs, convert state variable to (Qdp)
    ! because initial conditon reads in Q, not Qdp
    ! restart runs will read dpQ from restart file
    ! need to check what 1 does on a branch run
    if (runtype==0 .or. runtype==2) then
       do ie=nets,nete
          elem(ie)%derived%omega_p(:,:,:) = 0D0
       end do

       call TimeLevel_Qdp( tl, qsplit, n0_qdp)
       do ie=nets,nete
          do k=1,nlev    !  Loop inversion (AAM)
             do t=1,3
                do q=1,qsize       
                   do i=1,np
                      do j=1,np          
                         dp = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                              ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,t)
                         
                         elem(ie)%state%Qdp(i,j,k,q,n0_qdp)=elem(ie)%state%Q(i,j,k,q)*dp  
                         
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    endif
    
 ! do it only for SPELT/FVM tracers, FIRST TRACER will be the AIR DENSITY   
 ! should be optimize and combined with the above caculation 
    if (ntrac>0) then 
      ! do it only for FVM tracers, FIRST TRACER will be the AIR DENSITY   
      ! should be optimize and combined with the above caculation 
      do ie=nets,nete 
        do k=1,nlev
     	    do i=1,np
     	      do j=1,np      
         		  elem(ie)%derived%dp(i,j,k)=( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
    		       ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,tl%n0)
     	      enddo
     	    enddo
          !write air density in tracer 1 of FVM
          fvm(ie)%c(1:nc,1:nc,k,1,tl%n0)=interpolate_gll2fvm_points(elem(ie)%derived%dp(:,:,k),deriv(hybrid%ithr))
!            fvm(ie)%c(:,:,k,1,tl%n0)=1.0D0
        enddo
      enddo
      call fvm_init3(elem,fvm,hybrid,nets,nete,tl%n0)
      do ie=nets,nete 
   	    do i=1-nhc,nc+nhc
   	      do j=1-nhc,nc+nhc  
   	        fvm(ie)%psc(i,j) = sum(fvm(ie)%c(i,j,:,1,tl%n0)) +  hvcoord%hyai(1)*hvcoord%ps0
   	      enddo
   	    enddo
      enddo
      if (hybrid%masterthread) then
         write(iulog,*) 'FVM tracers (incl. in halo zone) initialized. FIRST tracer has air density!'	       
      end if
    endif  

    ! for restart runs, we read in Qdp for exact restart, and rederive Q
    if (runtype==1) then
       call TimeLevel_Qdp( tl, qsplit, n0_qdp)
       do ie=nets,nete
          do k=1,nlev    !  Loop inversion (AAM)
             do t=tl%n0,tl%n0
                do q=1,qsize       
                   do i=1,np
                      do j=1,np          
                         dp = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                              ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,t)
                         elem(ie)%state%Q(i,j,k,q)=elem(ie)%state%Qdp(i,j,k,q, n0_qdp)/dp 
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    endif


    ! timesteps to use for advective stability:  tstep*qsplit and tstep
    call print_cfl(elem,hybrid,nets,nete,dtnu)

    if (hybrid%masterthread) then 
       ! 1 has set tstep based on dtime before calling prim_init2(), 
       ! so only now does HOMME learn the timstep.  print them out:
       write(iulog,'(a,2f9.2)') "dt_remap: (0=disabled)   ",tstep*qsplit*rsplit
       write(iulog,'(a,2f9.2)') "dt_tracer, per RK stage: ",tstep*qsplit,(tstep*qsplit)/(rk_stage_user-1)
       write(iulog,'(a,2f9.2)') "dt_dyn:                  ",tstep
       write(iulog,'(a,2f9.2)') "dt_dyn (viscosity):      ",dt_dyn_vis
       write(iulog,'(a,2f9.2)') "dt_tracer (viscosity):   ",dt_tracer_vis

 
       if (phys_tscale/=0) then
          write(iulog,'(a,2f9.2)') "CAM physics timescale:       ",phys_tscale
       endif
       write(iulog,'(a,2f9.2)') "CAM dtime (dt_phys):         ",tstep*nsplit*qsplit*max(rsplit,1)
    end if


    if (hybrid%masterthread) write(iulog,*) "initial state:"
    call prim_printstate(elem, tl, hybrid,hvcoord,nets,nete, fvm)
  end subroutine prim_init2

!=======================================================================================================! 



! Subprogram not used   subroutine leapfrog_bootstrap(elem, hybrid,nets,nete,tstep,tl,hvcoord)
! Subprogram not used 
! Subprogram not used   !
! Subprogram not used   ! leapfrog bootstrap code.  
! Subprogram not used   !
! Subprogram not used   ! take the equivilent of one timestep, but do it with a 
! Subprogram not used   ! dt/2 euler and a dt/2 leapfrog step  
! Subprogram not used   !
! Subprogram not used   use hybvcoord_mod, only : hvcoord_t
! Subprogram not used   use time_mod, only : TimeLevel_t
! Subprogram not used 
! Subprogram not used   type (element_t) , intent(inout)        :: elem(:)
! Subprogram not used   type (hybrid_t), intent(in)           :: hybrid  ! distributed parallel structure (shared)
! Subprogram not used   type (hvcoord_t), intent(in)      :: hvcoord         ! hybrid vertical coordinate struct
! Subprogram not used   integer, intent(in)                     :: nets  ! starting thread element number (private)
! Subprogram not used   integer, intent(in)                     :: nete  ! ending thread element number   (private)
! Subprogram not used   real(kind=real_kind), intent(in)        :: tstep          ! "timestep dependent" timestep
! Subprogram not used   type (TimeLevel_t), intent(inout)       :: tl
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   ! local
! Subprogram not used   real(kind=real_kind) :: tstep_tmp,tstep_dyn
! Subprogram not used   integer :: i,ie
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   tstep_dyn = tstep
! Subprogram not used   ! forward euler to get to tstep_dyn/2 (keep t=0 in nm1 timelevel) 
! Subprogram not used   ! (note: leapfrog tstep_dyn/4 with nm1=n0 is Euler with tstep_dyn/2 )
! Subprogram not used   tstep_tmp=tstep_dyn/4        
! Subprogram not used 
! Subprogram not used   call prim_run(elem, hybrid,nets,nete, tstep_tmp, tl, hvcoord, "forward")
! Subprogram not used   
! Subprogram not used   ! leapfrog with tstep_dyn/2 to get to tstep_dyn (keep t=0 in nm1 timelevel)
! Subprogram not used   tstep_tmp=tstep_dyn/2
! Subprogram not used   call prim_run(elem, hybrid,nets,nete, tstep_tmp, tl, hvcoord, "forward")
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   tl%nstep=tl%nstep-1        ! count all of that as 1 timestep
! Subprogram not used 
! Subprogram not used   end subroutine leapfrog_bootstrap


!=======================================================================================================! 


! Subprogram not used   subroutine prim_run(elem, hybrid,nets,nete, dt, tl, hvcoord, advance_name)
! Subprogram not used     use hybvcoord_mod, only : hvcoord_t
! Subprogram not used     use time_mod, only : TimeLevel_t, timelevel_update, smooth
! Subprogram not used     use control_mod, only: statefreq, integration, ftype, qsplit
! Subprogram not used     use prim_advance_mod, only : prim_advance_exp, prim_advance_si, preq_robert3
! Subprogram not used     use prim_state_mod, only : prim_printstate, prim_diag_scalars, prim_energy_halftimes
! Subprogram not used     use parallel_mod, only : abortmp
! Subprogram not used 
! Subprogram not used     type (element_t) , intent(inout)        :: elem(:)
! Subprogram not used     type (hybrid_t), intent(in)           :: hybrid  ! distributed parallel structure (shared)
! Subprogram not used 
! Subprogram not used     type (hvcoord_t), intent(in)      :: hvcoord         ! hybrid vertical coordinate struct
! Subprogram not used 
! Subprogram not used     integer, intent(in)                     :: nets  ! starting thread element number (private)
! Subprogram not used     integer, intent(in)                     :: nete  ! ending thread element number   (private)
! Subprogram not used     real(kind=real_kind), intent(in)        :: dt              ! "timestep dependent" timestep
! Subprogram not used     type (TimeLevel_t), intent(inout)       :: tl
! Subprogram not used     character(len=*), intent(in) :: advance_name
! Subprogram not used     real(kind=real_kind) :: st, st1, dp
! Subprogram not used     integer :: ie, t, q,k,i,j
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     logical :: compute_diagnostics
! Subprogram not used 
! Subprogram not used     ! ===================================
! Subprogram not used     ! Main timestepping loop
! Subprogram not used     ! ===================================
! Subprogram not used 
! Subprogram not used     ! compute diagnostics and energy for STDOUT 
! Subprogram not used     ! compute energy if we are using an energy fixer
! Subprogram not used 
! Subprogram not used     compute_diagnostics=.false.
! Subprogram not used     if (MODULO(tl%nstep+1,statefreq)==0 .or. tl%nstep+1==tl%nstep0) then
! Subprogram not used        compute_diagnostics=.true.  
! Subprogram not used     endif
! Subprogram not used     tot_iter=0.0       
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     ! Forcing options for testing 1-HOMME energy balance:
! Subprogram not used     if (ftype == -1) then
! Subprogram not used        ! disable all forcing, but allow moisture:
! Subprogram not used        do ie=nets,nete
! Subprogram not used           elem(ie)%derived%FQ = 0
! Subprogram not used           elem(ie)%derived%FM = 0
! Subprogram not used           elem(ie)%derived%FT = 0
! Subprogram not used        enddo
! Subprogram not used     endif
! Subprogram not used     if (ftype == -2) then
! Subprogram not used        ! disable moisture, but allow dynamics forcing
! Subprogram not used        do ie=nets,nete
! Subprogram not used           elem(ie)%state%Q = 0
! Subprogram not used           elem(ie)%state%Qdp = 0
! Subprogram not used           elem(ie)%derived%FQ = 0
! Subprogram not used        enddo
! Subprogram not used     endif
! Subprogram not used     if (ftype == -3) then
! Subprogram not used        ! disable forcing & moisture
! Subprogram not used        do ie=nets,nete
! Subprogram not used           elem(ie)%state%Q = 0
! Subprogram not used           elem(ie)%state%Qdp = 0
! Subprogram not used           elem(ie)%derived%FQ = 0
! Subprogram not used           elem(ie)%derived%FM = 0
! Subprogram not used           elem(ie)%derived%FT = 0
! Subprogram not used        enddo
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     ! =================================
! Subprogram not used     ! energy, dissipation rate diagnostics.  Uses data at t-1,t 
! Subprogram not used     ! to compute diagnostics at t - 0.5.  
! Subprogram not used     ! small error in the t+.5 terms because at this
! Subprogram not used     ! point only state variables at t-1 has been Robert filtered.  
! Subprogram not used     ! =================================
! Subprogram not used     if (compute_diagnostics) then
! Subprogram not used        call prim_energy_halftimes(elem,hvcoord,tl,1,.true.,nets,nete)
! Subprogram not used        call prim_diag_scalars(elem,hvcoord,tl,1,.true.,nets,nete)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     ! ===============
! Subprogram not used     ! initialize mean flux accumulation variables
! Subprogram not used     ! ===============
! Subprogram not used     do ie=nets,nete
! Subprogram not used        elem(ie)%derived%eta_dot_dpdn=0
! Subprogram not used        elem(ie)%derived%omega_p=0
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     ! ===============
! Subprogram not used     ! Dynamical Step  uses Q at tl%n0
! Subprogram not used     ! ===============
! Subprogram not used !$OMP BARRIER
! Subprogram not used     if (integration == "semi_imp") then
! Subprogram not used        call prim_advance_si(elem, nets, nete, cg(hybrid%ithr), blkjac, red, &
! Subprogram not used             refstate, hvcoord, deriv(hybrid%ithr), flt, hybrid, tl, dt)
! Subprogram not used        tot_iter=tot_iter+cg(hybrid%ithr)%iter
! Subprogram not used     else if (integration == "full_imp") then
! Subprogram not used        call abortmp('full_imp integration requires tstep_type > 0')
! Subprogram not used     else 
! Subprogram not used        call prim_advance_exp(elem, deriv(hybrid%ithr), hvcoord,   &
! Subprogram not used             hybrid, dt, tl, nets, nete, compute_diagnostics)
! Subprogram not used 
! Subprogram not used        ! keep lnps up to date (we should get rid of this requirement)
! Subprogram not used        do ie=nets,nete
! Subprogram not used           elem(ie)%state%lnps(:,:,tl%np1)= LOG(elem(ie)%state%ps_v(:,:,tl%np1))
! Subprogram not used        enddo
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     ! =================================
! Subprogram not used     ! energy, dissipation rate diagnostics.  Uses data at t and t+1
! Subprogram not used     ! to compute diagnostics at t + 0.5.
! Subprogram not used     ! =================================
! Subprogram not used     if (compute_diagnostics) then
! Subprogram not used        call prim_energy_halftimes(elem,hvcoord,tl,2,.false.,nets,nete)
! Subprogram not used        call prim_diag_scalars(elem,hvcoord,tl,2,.false.,nets,nete)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used     ! =================================== 
! Subprogram not used     ! Compute Forcing Tendencies from nm1 data (for PROCESS SPLIT)
! Subprogram not used     ! or np1 data (for TIMESPLIT) and add tendencies into soluiton at timelevel np1
! Subprogram not used     ! ===================================       
! Subprogram not used     call abortmp('CAM-HOMME-SE requires RK timestepping option turned on')
! Subprogram not used     ! measure the effects of forcing
! Subprogram not used     if (compute_diagnostics) then
! Subprogram not used        call prim_energy_halftimes(elem,hvcoord,tl,3,.false.,nets,nete)
! Subprogram not used        call prim_diag_scalars(elem,hvcoord,tl,3,.false.,nets,nete)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     ! =================================
! Subprogram not used     ! timestep is complete.  
! Subprogram not used     ! =================================
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     !Now apply robert filter to all prognostic variables
! Subprogram not used     if (smooth/=0) &
! Subprogram not used        call preq_robert3(tl%nm1,tl%n0,tl%np1,elem,hvcoord,nets,nete)
! Subprogram not used     ! measure the effects of Robert filter
! Subprogram not used     if (compute_diagnostics) then
! Subprogram not used        call prim_energy_halftimes(elem,hvcoord,tl,4,.false.,nets,nete)
! Subprogram not used        call prim_diag_scalars(elem,hvcoord,tl,4,.false.,nets,nete)
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     ! =================================
! Subprogram not used     ! update dynamics time level pointers 
! Subprogram not used     ! =================================
! Subprogram not used     call TimeLevel_update(tl,advance_name)
! Subprogram not used 
! Subprogram not used   ! ============================================================
! Subprogram not used     ! Print some diagnostic information 
! Subprogram not used     ! ============================================================
! Subprogram not used 
! Subprogram not used     if (compute_diagnostics) then
! Subprogram not used        if (hybrid%masterthread) then 
! Subprogram not used           if (integration == "semi_imp") write(iulog,*) "cg its=",cg(0)%iter
! Subprogram not used        end if
! Subprogram not used        call prim_printstate(elem, tl, hybrid,hvcoord,nets,nete)
! Subprogram not used     end if
! Subprogram not used   end subroutine prim_run

!=======================================================================================================! 


  subroutine prim_run_subcycle(elem, fvm, hybrid,nets,nete, dt, tl, hvcoord,nsubstep)
!
!   advance all variables (u,v,T,ps,Q,C) from time t to t + dt_q
!     
!   input: 
!       tl%nm1   not used
!       tl%n0    data at time t
!       tl%np1   new values at t+dt_q
!   
!   then we update timelevel pointers:
!       tl%nm1 = tl%n0
!       tl%n0  = tl%np1
!   so that:
!       tl%nm1   tracers:  t    dynamics:  t+(qsplit-1)*dt
!       tl%n0    time t + dt_q
!   
!
    use hybvcoord_mod, only : hvcoord_t
    use time_mod, only : TimeLevel_t, timelevel_update, timelevel_qdp, nsplit
    use control_mod, only: statefreq,&
           energy_fixer, ftype, qsplit, rsplit, test_cfldep
    use prim_advance_mod, only : applycamforcing, &
                                 applycamforcing_dynamics
    use prim_state_mod, only : prim_printstate, prim_diag_scalars, prim_energy_halftimes
    use parallel_mod, only : abortmp
    use reduction_mod, only : parallelmax
    use prim_advection_mod, only : vertical_remap
    

    type (element_t) , intent(inout)        :: elem(:)
    
      type(fvm_struct), intent(inout) :: fvm(:)
    type (hybrid_t), intent(in)           :: hybrid  ! distributed parallel structure (shared)

    type (hvcoord_t), intent(in)      :: hvcoord         ! hybrid vertical coordinate struct

    integer, intent(in)                     :: nets  ! starting thread element number (private)
    integer, intent(in)                     :: nete  ! ending thread element number   (private)
    real(kind=real_kind), intent(in)        :: dt  ! "timestep dependent" timestep
    type (TimeLevel_t), intent(inout)       :: tl
    integer, intent(in)                     :: nsubstep  ! nsubstep = 1 .. nsplit
    real(kind=real_kind) :: st, st1, dp, dt_q, dt_remap
    integer :: ie, t, q,k,i,j,n, n_Q
    integer :: n0_qdp,np1_qdp,r, nstep_end

    real (kind=real_kind)                          :: maxcflx, maxcfly  
    real (kind=real_kind) :: dp_np1(np,np)
    logical :: compute_diagnostics, compute_energy


    ! ===================================
    ! Main timestepping loop
    ! ===================================
    dt_q = dt*qsplit
    dt_remap = dt_q
    nstep_end = tl%nstep + qsplit
    if (rsplit>0) then
       dt_remap=dt_q*rsplit   ! rsplit=0 means use eulerian code, not vert. lagrange
       nstep_end = tl%nstep + qsplit*rsplit  ! nstep at end of this routine
    endif



    ! compute diagnostics and energy for STDOUT 
    ! compute energy if we are using an energy fixer
    compute_diagnostics=.false.
    compute_energy=energy_fixer > 0

    if (MODULO(nstep_end,statefreq)==0 .or. nstep_end==tl%nstep0) then
       compute_diagnostics=.true.  
       compute_energy = .true.
    endif
    if (compute_diagnostics) &
       call prim_diag_scalars(elem,hvcoord,tl,4,.true.,nets,nete)



    ! ftype=2  Q was adjusted by physics, but apply u,T forcing here
    ! ftype=1  forcing was applied time-split in 1 coupling layer
    ! ftype=0 means forcing apply here
    ! ftype=-1 do not apply forcing
    call TimeLevel_Qdp( tl, qsplit, n0_qdp)
    if (ftype==0) call ApplyCAMForcing(elem, hvcoord,tl%n0,n0_qdp, dt_remap,nets,nete)
    if (ftype==2) call ApplyCAMForcing_dynamics(elem, hvcoord,tl%n0,dt_remap,nets,nete)

    ! E(1) Energy after 1 forcing
    if (compute_energy) call prim_energy_halftimes(elem,hvcoord,tl,1,.true.,nets,nete)

    ! qmass and variance, using Q(n0),Qdp(n0)
    if (compute_diagnostics) call prim_diag_scalars(elem,hvcoord,tl,1,.true.,nets,nete)

    ! initialize dp3d from ps
    if (rsplit>0) then
    do ie=nets,nete
       do k=1,nlev
          elem(ie)%state%dp3d(:,:,k,tl%n0)=&
               ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
               ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,tl%n0)
       enddo
       ! DEBUGDP step: ps_v should not be used for rsplit>0 code during prim_step
       ! vertical_remap.  so to this for debugging:
       elem(ie)%state%ps_v(:,:,tl%n0)=-9e9
    enddo
    endif


    ! loop over rsplit vertically lagrangian timesteps
    call prim_step(elem, fvm, hybrid,nets,nete, dt, tl, hvcoord,compute_diagnostics)
    do r=2,rsplit   
       call TimeLevel_update(tl,"leapfrog")
       call prim_step(elem, fvm, hybrid,nets,nete, dt, tl, hvcoord,.false.)
    enddo
    ! defer final timelevel update until after remap and diagnostics


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !  apply vertical remap
    !  always for tracers  
    !  if rsplit>0:  also remap dynamics and compute reference level ps_v
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !compute timelevels for tracers (no longer the same as dynamics)
    call TimeLevel_Qdp( tl, qsplit, n0_qdp, np1_qdp) 
    call vertical_remap(elem,fvm,hvcoord,dt_remap,tl%np1,np1_qdp,nets,nete)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! time step is complete.  update some diagnostic variables:
    ! lnps (we should get rid of this)
    ! Q    (mixing ratio)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    do ie=nets,nete
       elem(ie)%state%lnps(:,:,tl%np1)= LOG(elem(ie)%state%ps_v(:,:,tl%np1))
       do k=1,nlev    !  Loop inversion (AAM)
          dp_np1(:,:) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
               ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,tl%np1)
          do q=1,qsize
             elem(ie)%state%Q(:,:,k,q)=elem(ie)%state%Qdp(:,:,k,q,np1_qdp)/dp_np1(:,:) 
          enddo
       enddo
    enddo
    


   

    ! now we have:
    !   u(nm1)   dynamics at  t+dt_remap - 2*dt 
    !   u(n0)    dynamics at  t+dt_remap - dt    
    !   u(np1)   dynamics at  t+dt_remap
    !
    !   Q(1)   Q at t+dt_remap
    if (compute_diagnostics) call prim_diag_scalars(elem,hvcoord,tl,2,.false.,nets,nete)
    if (compute_energy) call prim_energy_halftimes(elem,hvcoord,tl,2,.false.,nets,nete)

    if (energy_fixer > 0) then
       call prim_energy_fixer(elem,hvcoord,hybrid,tl,nets,nete,nsubstep)
    endif

    if (compute_diagnostics) then
       call prim_diag_scalars(elem,hvcoord,tl,3,.false.,nets,nete)
       call prim_energy_halftimes(elem,hvcoord,tl,3,.false.,nets,nete)
     endif

    ! =================================
    ! update dynamics time level pointers 
    ! =================================
    call TimeLevel_update(tl,"leapfrog")

    ! now we have:
    !   u(nm1)   dynamics at  t+dt_remap - dt       (Robert-filtered)
    !   u(n0)    dynamics at  t+dt_remap 
    !   u(np1)   undefined
 

    ! ============================================================
    ! Print some diagnostic information 
    ! ============================================================
    if (compute_diagnostics) then
       call prim_printstate(elem, tl, hybrid,hvcoord,nets,nete, fvm)
    end if
  end subroutine prim_run_subcycle






  subroutine prim_step(elem, fvm, hybrid,nets,nete, dt, tl, hvcoord, compute_diagnostics)
!
!   Take qsplit dynamics steps and one tracer step
!   for vertically lagrangian option, this subroutine does only the horizontal step
!     
!   input: 
!       tl%nm1   not used
!       tl%n0    data at time t
!       tl%np1   new values at t+dt_q
!   
!   then we update timelevel pointers:
!       tl%nm1 = tl%n0
!       tl%n0  = tl%np1
!   so that:
!       tl%nm1   tracers:  t    dynamics:  t+(qsplit-1)*dt
!       tl%n0    time t + dt_q
!   
!
    use hybvcoord_mod, only : hvcoord_t
    use time_mod, only : TimeLevel_t, timelevel_update, nsplit
    use control_mod, only: statefreq, integration, ftype, qsplit, nu_p, test_cfldep, rsplit
    use prim_advance_mod, only : prim_advance_exp, overwrite_SEdensity
    use prim_advection_mod, only : prim_advec_tracers_remap_rk2, prim_advec_tracers_fvm, &
         prim_advec_tracers_spelt
    use parallel_mod, only : abortmp
    use reduction_mod, only : parallelmax
    use derivative_mod, only : interpolate_gll2spelt_points

    type (element_t) , intent(inout)        :: elem(:)
    
      type(fvm_struct), intent(inout) :: fvm(:)
    type (hybrid_t), intent(in)           :: hybrid  ! distributed parallel structure (shared)

    type (hvcoord_t), intent(in)      :: hvcoord         ! hybrid vertical coordinate struct

    integer, intent(in)                     :: nets  ! starting thread element number (private)
    integer, intent(in)                     :: nete  ! ending thread element number   (private)
    real(kind=real_kind), intent(in)        :: dt  ! "timestep dependent" timestep
    type (TimeLevel_t), intent(inout)       :: tl
    real(kind=real_kind) :: st, st1, dp, dt_q
    integer :: ie, t, q,k,i,j,n, n_Q

    real (kind=real_kind)                          :: maxcflx, maxcfly 
     
    real (kind=real_kind) :: dp_np1(np,np)
    logical :: compute_diagnostics

    dt_q = dt*qsplit

    ! ===============
    ! initialize mean flux accumulation variables and save some variables at n0 
    ! for use by advection
    ! ===============
    do ie=nets,nete
      elem(ie)%derived%eta_dot_dpdn=0     ! mean vertical mass flux
      elem(ie)%derived%vn0=0              ! mean horizontal mass flux
      elem(ie)%derived%omega_p=0          
      if (nu_p>0) then
         elem(ie)%derived%dpdiss_ave=0
         elem(ie)%derived%dpdiss_biharmonic=0
      endif
      if (ntrac>0) then
        ! save velocity at time t for fvm
        fvm(ie)%vn0=elem(ie)%state%v(:,:,:,:,tl%n0)
      end if

      if (rsplit==0) then
      ! save dp at time t for use in tracers
         do k=1,nlev
            elem(ie)%derived%dp(:,:,k)=&
                 ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                 ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,tl%n0)
         enddo
      else
         ! dp at time t:  use floating lagrangian levels:
         elem(ie)%derived%dp(:,:,:)=elem(ie)%state%dp3d(:,:,:,tl%n0)
      endif
    enddo

    ! ===============
    ! Dynamical Step 
    ! ===============
    n_Q = tl%n0  ! n_Q = timelevel of FV tracers at time t.  need to save this
                 ! FV tracers still carry 3 timelevels 
                 ! SE tracers only carry 2 timelevels 
    call prim_advance_exp(elem, deriv(hybrid%ithr), hvcoord,   &
         hybrid, dt, tl, nets, nete, compute_diagnostics)
    do n=2,qsplit
       call TimeLevel_update(tl,"leapfrog")
       call prim_advance_exp(elem, deriv(hybrid%ithr), hvcoord,   &
            hybrid, dt, tl, nets, nete, .false.)
       ! defer final timelevel update until after Q update.
    enddo

    ! ===============
    ! Tracer Advection.  SE advection uses mean flux variables:
    !        derived%dp              =  dp at start of timestep
    !        derived%vn0             =  mean horiz. flux:   U*dp
    ! in addition, this routine will apply the DSS to:
    !        derived%eta_dot_dpdn    =  mean vertical velocity (used for remap below)
    !        derived%omega           =  
    ! ===============
    if (qsize>0) call Prim_Advec_Tracers_remap_rk2(elem, deriv(hybrid%ithr),hvcoord,flt_advection,hybrid,&
         dt_q,tl,nets,nete)


    if (ntrac>0) then
      if ( n_Q /= tl%n0 ) then
        ! make sure tl%n0 contains tracers at start of timestep
        do ie=nets,nete
          fvm(ie)%c(:,:,:,1:ntrac,tl%n0)  = fvm(ie)%c(:,:,:,1:ntrac,n_Q)
        enddo
      endif 
      call Prim_Advec_Tracers_fvm(elem, fvm, deriv(hybrid%ithr),hvcoord,hybrid,&
           dt_q,tl,nets,nete)
           ! values in the halo zone are only in np1 at this time
       do ie=nets,nete 
         do i=1-nhc,nc+nhc
           do j=1-nhc,nc+nhc  
             fvm(ie)%psc(i,j) = sum(fvm(ie)%c(i,j,:,1,tl%np1)) +  hvcoord%hyai(1)*hvcoord%ps0
           enddo
         enddo
       enddo

       if(test_cfldep) then
         maxcflx=0.0D0
         maxcfly=0.0D0
         do k=1, nlev

!            maxcflx = parallelmax(fvm(:)%maxcfl(1,k),hybrid)
!            maxcfly = parallelmax(fvm(:)%maxcfl(2,k),hybrid) 
           maxcflx = max(maxcflx,parallelmax(fvm(:)%maxcfl(1,k),hybrid))
           maxcfly = max(maxcfly,parallelmax(fvm(:)%maxcfl(2,k),hybrid))
          end do
           
           if(hybrid%masterthread) then 
             write(*,*) "nstep",tl%nstep,"dt_q=", dt_q, "maximum over all Level"
             write(*,*) "CFL: maxcflx=", maxcflx, "maxcfly=", maxcfly 
             print *
           endif 
       endif   
       !overwrite SE density by fvm(ie)%psc
!        call overwrite_SEdensity(elem,fvm,dt_q,hybrid,nets,nete,tl%np1) 
    endif

  end subroutine prim_step


!=======================================================================================================! 


! Subprogram not used   subroutine prim_finalize(hybrid)
! Subprogram not used     type (hybrid_t), intent(in)           :: hybrid  ! distributed parallel structure (shared)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     ! ==========================
! Subprogram not used     ! end of the hybrid program
! Subprogram not used     ! ==========================
! Subprogram not used   end subroutine prim_finalize



!=======================================================================================================! 
! Subprogram not used   subroutine prim_energy_fixer(elem,hvcoord,hybrid,tl,nets,nete,nsubstep)
! Subprogram not used ! 
! Subprogram not used ! non-subcycle code:
! Subprogram not used !  Solution is given at times u(t-1),u(t),u(t+1)
! Subprogram not used !  E(n=1) = energy before dynamics
! Subprogram not used !  E(n=2) = energy after dynamics
! Subprogram not used !
! Subprogram not used !  fixer will add a constant to the temperature so E(n=2) = E(n=1)
! Subprogram not used !
! Subprogram not used     use parallel_mod, only: global_shared_buf, global_shared_sum
! Subprogram not used     use kinds, only : real_kind
! Subprogram not used     use hybvcoord_mod, only : hvcoord_t
! Subprogram not used     use physical_constants, only : Cp 
! Subprogram not used     use time_mod, only : timelevel_t
! Subprogram not used     use control_mod, only : use_cpstar, energy_fixer
! Subprogram not used     use hybvcoord_mod, only : hvcoord_t
! Subprogram not used     use global_norms_mod, only: wrap_repro_sum
! Subprogram not used     use parallel_mod, only : abortmp
! Subprogram not used     type (hybrid_t), intent(in)           :: hybrid  ! distributed parallel structure (shared)
! Subprogram not used     integer :: t2,n,nets,nete
! Subprogram not used     type (element_t)     , intent(inout), target :: elem(:)
! Subprogram not used     type (hvcoord_t)                  :: hvcoord
! Subprogram not used     type (TimeLevel_t), intent(inout)       :: tl
! Subprogram not used     integer, intent(in)                    :: nsubstep
! Subprogram not used 
! Subprogram not used     integer :: ie,k,i,j,nmax
! Subprogram not used     real (kind=real_kind), dimension(np,np,nlev)  :: dp   ! delta pressure
! Subprogram not used     real (kind=real_kind), dimension(np,np,nlev)  :: sumlk
! Subprogram not used     real (kind=real_kind), pointer  :: PEner(:,:,:)
! Subprogram not used     real (kind=real_kind), dimension(np,np)  :: suml
! Subprogram not used     real (kind=real_kind) :: psum(nets:nete,4),psum_g(4),beta
! Subprogram not used 
! Subprogram not used     ! when forcing is applied during dynamics timstep, actual forcing is
! Subprogram not used     ! slightly different (about 0.1 W/m^2) then expected by the physics
! Subprogram not used     ! since u & T are changing while FU and FT are held constant.
! Subprogram not used     ! to correct for this, save compute de_from_forcing at step 1
! Subprogram not used     ! and then adjust by:  de_from_forcing_step1 - de_from_forcing_stepN
! Subprogram not used     real (kind=real_kind),save :: de_from_forcing_step1
! Subprogram not used     real (kind=real_kind)      :: de_from_forcing
! Subprogram not used 
! Subprogram not used     t2=tl%np1    ! timelevel for T
! Subprogram not used     if (use_cpstar /= 0 ) then
! Subprogram not used        call abortmp('Energy fixer requires use_cpstar=0')
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     psum = 0
! Subprogram not used     do ie=nets,nete
! Subprogram not used 
! Subprogram not used        do k=1,nlev
! Subprogram not used           dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
! Subprogram not used                ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,t2)
! Subprogram not used        enddo
! Subprogram not used        suml=0
! Subprogram not used        do k=1,nlev
! Subprogram not used           do i=1,np
! Subprogram not used           do j=1,np
! Subprogram not used                 sumlk(i,j,k) = cp*dp(i,j,k) 
! Subprogram not used           enddo
! Subprogram not used           enddo
! Subprogram not used        enddo
! Subprogram not used        suml=0
! Subprogram not used        do k=1,nlev
! Subprogram not used           do i=1,np
! Subprogram not used           do j=1,np
! Subprogram not used              suml(i,j) = suml(i,j) + sumlk(i,j,k)
! Subprogram not used           enddo
! Subprogram not used           enddo
! Subprogram not used        enddo
! Subprogram not used        PEner => elem(ie)%accum%PEner(:,:,:)
! Subprogram not used 
! Subprogram not used        ! psum(:,4) = energy before forcing
! Subprogram not used        ! psum(:,1) = energy after forcing, before dynamics
! Subprogram not used        ! psum(:,2) = energy after dynamics
! Subprogram not used        ! psum(:,3) = cp*dp (internal energy added is beta*psum(:,3))
! Subprogram not used        psum(ie,3) = psum(ie,3) + SUM(suml(:,:)*elem(ie)%spheremp(:,:))
! Subprogram not used        do n=1,2
! Subprogram not used           psum(ie,n) = psum(ie,n) + SUM(  elem(ie)%spheremp(:,:)*&
! Subprogram not used                (PEner(:,:,n) + &
! Subprogram not used                elem(ie)%accum%IEner(:,:,n) + &
! Subprogram not used                elem(ie)%accum%KEner(:,:,n) ) )
! Subprogram not used        enddo
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     nmax=3
! Subprogram not used 
! Subprogram not used     do ie=nets,nete
! Subprogram not used        do n=1,nmax
! Subprogram not used           global_shared_buf(ie,n) = psum(ie,n)
! Subprogram not used        enddo
! Subprogram not used     enddo
! Subprogram not used     call wrap_repro_sum(nvars=nmax, comm=hybrid%par%comm)
! Subprogram not used     do n=1,nmax
! Subprogram not used        psum_g(n) = global_shared_sum(n)
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     beta = ( psum_g(1)-psum_g(2) )/psum_g(3)
! Subprogram not used 
! Subprogram not used     ! apply fixer
! Subprogram not used     do ie=nets,nete
! Subprogram not used        elem(ie)%state%T(:,:,:,t2) =  elem(ie)%state%T(:,:,:,t2) + beta
! Subprogram not used     enddo
! Subprogram not used     end subroutine prim_energy_fixer
!=======================================================================================================! 



! Subprogram not used     subroutine smooth_topo_datasets(phis,sghdyn,sgh30dyn,elem,hybrid,nets,nete)
! Subprogram not used     use control_mod, only : smooth_phis_numcycle,smooth_sgh_numcycle
! Subprogram not used     use hybrid_mod, only : hybrid_t
! Subprogram not used     use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack
! Subprogram not used     use bndry_mod, only : bndry_exchangev
! Subprogram not used     use derivative_mod, only : derivative_t , laplace_sphere_wk
! Subprogram not used     use viscosity_mod, only : biharmonic_wk
! Subprogram not used     use prim_advance_mod, only : smooth_phis
! Subprogram not used     implicit none
! Subprogram not used     
! Subprogram not used     integer , intent(in) :: nets,nete
! Subprogram not used     real (kind=real_kind), intent(inout)   :: phis(np,np,nets:nete)
! Subprogram not used     real (kind=real_kind), intent(inout)   :: sghdyn(np,np,nets:nete)
! Subprogram not used     real (kind=real_kind), intent(inout)   :: sgh30dyn(np,np,nets:nete)
! Subprogram not used     type (hybrid_t)      , intent(in) :: hybrid
! Subprogram not used     type (element_t)     , intent(inout), target :: elem(:)
! Subprogram not used     ! local
! Subprogram not used     integer :: ie
! Subprogram not used     real (kind=real_kind) :: minf 
! Subprogram not used 
! Subprogram not used     minf=-9e9
! Subprogram not used     if (hybrid%masterthread) &
! Subprogram not used        write(iulog,*) "Applying hyperviscosity smoother to PHIS"
! Subprogram not used     call smooth_phis(phis,elem,hybrid,deriv(hybrid%ithr),nets,nete,minf,smooth_phis_numcycle)
! Subprogram not used 
! Subprogram not used     minf=0
! Subprogram not used     if (hybrid%masterthread) &
! Subprogram not used        write(iulog,*) "Applying hyperviscosity smoother to SGH"
! Subprogram not used     call smooth_phis(sghdyn,elem,hybrid,deriv(hybrid%ithr),nets,nete,minf,smooth_sgh_numcycle)
! Subprogram not used     if (hybrid%masterthread) &
! Subprogram not used        write(iulog,*) "Applying hyperviscosity smoother to SGH30"
! Subprogram not used     call smooth_phis(sgh30dyn,elem,hybrid,deriv(hybrid%ithr),nets,nete,minf,smooth_sgh_numcycle)
! Subprogram not used 
! Subprogram not used     end subroutine smooth_topo_datasets

end module prim_driver_mod



