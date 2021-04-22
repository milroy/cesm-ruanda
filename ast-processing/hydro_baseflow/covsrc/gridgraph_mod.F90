



module GridGraph_mod
  !-------------------------
  use kinds, only : real_kind, iulog
  !-------------------------------
  use dimensions_mod, only : max_neigh_edges
  !-------------------------
  use control_mod, only : north, south, east, west, neast, nwest, seast, swest
  !-----
  implicit none


  private

  integer, public, parameter :: num_neighbors=8 ! for north, south, east, west, neast, nwest, seast, swest


  type, public :: GridVertex_t

      integer, pointer          :: nbrs(:) => null()           ! The numbers of the neighbor elements
      integer, pointer          :: nbrs_face(:) => null()      ! The cube face number of the neighbor element (nbrs array)
      integer, pointer          :: nbrs_wgt(:) => null()       ! The weights for edges defined by nbrs array
      integer, pointer          :: nbrs_wgt_ghost(:) => null() ! The weights for edges defined by nbrs array
      integer                   :: nbrs_ptr(num_neighbors + 1) !index into the nbrs array for each neighbor direction

      integer                   :: face_number           ! which face of the cube this vertex is on
      integer                   :: number                ! element number
      integer                   :: processor_number      ! processor number
      integer                   :: SpaceCurve  ! index in Space-Filling curve
  end type GridVertex_t

  type, public :: EdgeIndex_t
      integer, pointer            :: ixP(:) => null()
      integer, pointer            :: iyP(:) => null()
  end type EdgeIndex_t

  type, public :: GridEdge_t
      integer                      :: head_face  ! needed if head vertex has shape (i.e. square)
      integer                      :: tail_face  ! needed if tail vertex has shape (i.e. square)
      integer                      :: head_dir   !which of 8 neighbor directions is the head
      integer                      :: tail_dir   !which of 8 neighbor directions is the tail
      type (GridVertex_t),pointer  :: head => null()  ! edge head vertex
      type (GridVertex_t),pointer  :: tail => null()  ! edge tail vertex
      logical                      :: reverse

  end type GridEdge_t
  
! ==========================================
! Public Interfaces
! ==========================================

  public :: set_GridVertex_number
  public :: PrintGridVertex
 
  public :: allocate_gridvertex_nbrs
  public :: deallocate_gridvertex_nbrs
  public :: initgridedge
  public :: gridedge_search
  public :: gridedge_type
  public :: grid_edge_uses_vertex
  public :: PrintGridEdge
  public :: CheckGridNeighbors
  public :: PrintChecksum

  public :: CreateSubGridGraph
  public :: FreeGraph

  public :: assignment ( = ) 

  interface assignment ( = )
      module procedure copy_gridedge
      module procedure copy_edgeindex
      module procedure copy_gridvertex
  end interface

contains

!======================================================================

  subroutine allocate_gridvertex_nbrs(vertex, dim)

    type (GridVertex_t), intent(inout)   :: vertex
    integer, optional, intent(in)        :: dim
    integer                              :: num

    if (present(dim)) then
       num = dim
    else
       num = max_neigh_edges
    end if

    allocate(vertex%nbrs(num))
    allocate(vertex%nbrs_face(num))
    allocate(vertex%nbrs_wgt(num))
    allocate(vertex%nbrs_wgt_ghost(num))
 

  end subroutine allocate_gridvertex_nbrs
!======================================================================

  subroutine deallocate_gridvertex_nbrs(vertex)

    type (GridVertex_t), intent(inout)   :: vertex

    deallocate(vertex%nbrs)
    deallocate(vertex%nbrs_face)
    deallocate(vertex%nbrs_wgt)
    deallocate(vertex%nbrs_wgt_ghost)
 
  end subroutine deallocate_gridvertex_nbrs

!======================================================================

! =====================================
! copy edge:
! copy device for overloading = sign.
! =====================================


  recursive subroutine copy_gridedge(edge2, edge1)

    type (GridEdge_t), intent(out) :: edge2
    type (GridEdge_t), intent(in)  :: edge1

    edge2%tail_face = edge1%tail_face
    edge2%head_face = edge1%head_face
    edge2%tail_dir = edge1%tail_dir
    edge2%head_dir = edge1%head_dir
    edge2%reverse   = edge1%reverse

    if (associated(edge1%tail)) then
       edge2%tail=>edge1%tail
    end if
    if (associated(edge1%head)) then
       edge2%head=>edge1%head
    end if

  end subroutine copy_gridedge

!======================================================================

  recursive subroutine copy_gridvertex(vertex2, vertex1)
        
    implicit none 

    type (GridVertex_t), intent(out)   :: vertex2
    type (GridVertex_t), intent(in)    :: vertex1

    integer                            :: i,j,n
   
     n = SIZE(vertex1%nbrs)

     if (associated(vertex2%nbrs)) then
        nullify(vertex2%nbrs)
     end if
     if (associated(vertex2%nbrs_face)) then
        nullify(vertex2%nbrs_face)
     end if
     if (associated(vertex2%nbrs_wgt)) then
        nullify(vertex2%nbrs_wgt)
     end if
     if (associated(vertex2%nbrs_wgt_ghost)) then
        nullify(vertex2%nbrs_wgt_ghost)
     end if

     call allocate_gridvertex_nbrs(vertex2)

     do i=1,n
        vertex2%nbrs(i) = vertex1%nbrs(i)
        vertex2%nbrs_face(i) = vertex1%nbrs_face(i)
        vertex2%nbrs_wgt(i)  = vertex1%nbrs_wgt(i)
        vertex2%nbrs_wgt_ghost(i)  = vertex1%nbrs_wgt_ghost(i)
     enddo

     do i=1, num_neighbors+1
        vertex2%nbrs_ptr(i) = vertex1%nbrs_ptr(i)
     enddo

     vertex2%face_number     = vertex1%face_number
     vertex2%number     = vertex1%number
     vertex2%processor_number  = vertex1%processor_number 
     vertex2%SpaceCurve = vertex1%SpaceCurve 

  end subroutine copy_gridvertex

!======================================================================
! Subprogram not used   recursive subroutine copy_edgeindex(index2,index1)
! Subprogram not used   
! Subprogram not used   type (EdgeIndex_t), intent(out) :: index2
! Subprogram not used   type (EdgeIndex_t), intent(in)  :: index1
! Subprogram not used 
! Subprogram not used   if(associated(index1%iyP)) then 
! Subprogram not used     index2%iyP => index1%iyP
! Subprogram not used   endif
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   if(associated(index1%ixP)) then 
! Subprogram not used      index2%ixP => index1%ixP
! Subprogram not used   endif
! Subprogram not used  
! Subprogram not used   end subroutine copy_edgeindex

!======================================================================
! Subprogram not used   subroutine FreeGraph(Vertex)
! Subprogram not used 
! Subprogram not used      implicit none
! Subprogram not used      type (GridVertex_t)           :: Vertex(:)
! Subprogram not used      integer                       :: i,nelem
! Subprogram not used 
! Subprogram not used      nelem = SIZE(Vertex)
! Subprogram not used 
! Subprogram not used !JMD     do i=1,nelem
! Subprogram not used !JMD        deallocate(Vertex(i)%wgtV)
! Subprogram not used !JMD        deallocate(Vertex(i)%wgtG)
! Subprogram not used !JMD        deallocate(Vertex(i)%nbrs)
! Subprogram not used !JMD     enddo
! Subprogram not used 
! Subprogram not used   end subroutine FreeGraph

!======================================================================

!===========================
! search edge list for match
!===========================

! Subprogram not used   function gridedge_search(nvert1, nvert2, edge) result(number)
! Subprogram not used 
! Subprogram not used     integer, intent(in) :: nvert1
! Subprogram not used     integer, intent(in) :: nvert2
! Subprogram not used     type(GridEdge_t), intent(in) :: edge(:)
! Subprogram not used     integer :: number
! Subprogram not used 
! Subprogram not used     integer :: tmp
! Subprogram not used     integer :: head
! Subprogram not used     integer :: tail
! Subprogram not used 
! Subprogram not used     integer :: nedge
! Subprogram not used     integer :: i
! Subprogram not used 
! Subprogram not used     nedge=SIZE(edge)
! Subprogram not used 
! Subprogram not used     tail=nvert1
! Subprogram not used     head=nvert2
! Subprogram not used 
! Subprogram not used     if (tail > head) then
! Subprogram not used        tmp  = tail
! Subprogram not used        tail = head
! Subprogram not used        head = tmp
! Subprogram not used     end if
! Subprogram not used 
! Subprogram not used     do i=1,nedge
! Subprogram not used        if (edge(i)%tail%number==tail .and. edge(i)%head%number==head)then
! Subprogram not used           number=i
! Subprogram not used        end if
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used   end function gridedge_search

!======================================================================

  function gridedge_type(edge) result(type)

    use params_mod, only : INTERNAL_EDGE, EXTERNAL_EDGE
    type (GridEdge_t), intent(in)  :: edge
    integer                        :: type

    if (edge%head%processor_number==edge%tail%processor_number) then
        type=INTERNAL_EDGE
    else
        type=EXTERNAL_EDGE
    endif

  end function gridedge_type

!======================================================================



! Subprogram not used   function grid_edge_uses_vertex(Vertex,Edge) result(log)
! Subprogram not used 
! Subprogram not used     type(GridVertex_t), intent(in) :: Vertex
! Subprogram not used     type(GridEdge_t),   intent(in) :: Edge
! Subprogram not used     logical :: log
! Subprogram not used     integer  :: number
! Subprogram not used 
! Subprogram not used     number = Vertex%number
! Subprogram not used     if(number == Edge%head%number .or. number == Edge%tail%number) then
! Subprogram not used         log = .TRUE.
! Subprogram not used     else
! Subprogram not used         log = .FALSE.
! Subprogram not used     endif
! Subprogram not used 
! Subprogram not used   end function grid_edge_uses_vertex

!======================================================================

! Subprogram not used   subroutine PrintChecksum(TestPattern,Checksum)
! Subprogram not used 
! Subprogram not used    use dimensions_mod, only : nlev, nelemd, np
! Subprogram not used 
! Subprogram not used    implicit none
! Subprogram not used 
! Subprogram not used    real(kind=real_kind), target,intent(in)   :: TestPattern(:,:,:,:)
! Subprogram not used    real(kind=real_kind), target,intent(in)   :: Checksum(:,:,:,:)
! Subprogram not used 
! Subprogram not used    integer                                  :: i,k,ix,iy
! Subprogram not used 
! Subprogram not used    print *
! Subprogram not used    write (iulog,*) 'checksums:'
! Subprogram not used    do i=1,nelemd
! Subprogram not used      !  Lets start out only looking at the first element
! Subprogram not used         write(iulog,*)
! Subprogram not used         do k=1,nlev
! Subprogram not used         do iy=1,np
! Subprogram not used         do ix=1,np
! Subprogram not used            write(iulog,*)INT(TestPattern(ix,iy,k,i))," checksum = ",INT(Checksum(ix,iy,k,i))
! Subprogram not used         enddo
! Subprogram not used         enddo
! Subprogram not used         enddo
! Subprogram not used    enddo
! Subprogram not used 
! Subprogram not used 
! Subprogram not used   end subroutine PrintChecksum

!======================================================================

! Subprogram not used   subroutine CreateSubGridGraph(Vertex, SVertex, local2global)
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used         
! Subprogram not used     type (GridVertex_t),intent(in)         :: Vertex(:)
! Subprogram not used     type (GridVertex_t),intent(inout)      :: SVertex(:)
! Subprogram not used     integer,intent(in)                     :: local2global(:)
! Subprogram not used 
! Subprogram not used     integer                                :: nelem,nelem_s,n,ncount,cnt,pos, orig_start
! Subprogram not used     integer                                :: inbr,i,ig,j,k, new_pos
! Subprogram not used     
! Subprogram not used     integer,allocatable                    :: global2local(:)
! Subprogram not used     logical, parameter    :: Debug = .FALSE.
! Subprogram not used 
! Subprogram not used     nelem   = SIZE(Vertex)
! Subprogram not used     nelem_s = SiZE(SVertex) 
! Subprogram not used 
! Subprogram not used     if(Debug) write(iulog,*)'CreateSubGridGraph: point #1'
! Subprogram not used 
! Subprogram not used     allocate(global2local(nelem))
! Subprogram not used 
! Subprogram not used     if(Debug) write(iulog,*)'CreateSubGridGraph: point #2'
! Subprogram not used 
! Subprogram not used     global2local(:) = 0
! Subprogram not used     do i=1,nelem_s
! Subprogram not used         ig = local2global(i)
! Subprogram not used         global2local(ig) = i
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used     if(Debug) write(iulog,*)'CreateSubGridGraph: point #3'
! Subprogram not used 
! Subprogram not used     do i=1,nelem_s
! Subprogram not used         ig = local2global(i)
! Subprogram not used 
! Subprogram not used         if(Debug) write(iulog,*)'CreateSubGridGraph: point #4'
! Subprogram not used         call copy_gridvertex(SVertex(i),Vertex(ig))  !svertex(i) = vertex(ig)
! Subprogram not used 
! Subprogram not used         n = SIZE(SVertex(i)%nbrs(:))
! Subprogram not used         ! ==============================================
! Subprogram not used         ! Apply the correction to the neighbors list to 
! Subprogram not used         ! reflect new subgraph numbers 
! Subprogram not used         ! ==============================================
! Subprogram not used 
! Subprogram not used 
! Subprogram not used         if(Debug) write(iulog,*)'CreateSubGridGraph: point #5'
! Subprogram not used         
! Subprogram not used         orig_start = 1
! Subprogram not used 
! Subprogram not used         do j=1,num_neighbors
! Subprogram not used 
! Subprogram not used            cnt = Svertex(i)%nbrs_ptr(j+1) - orig_start  !number of neighbors for this direction
! Subprogram not used            ncount = 0
! Subprogram not used            do k = 1, cnt
! Subprogram not used               pos = orig_start + k-1
! Subprogram not used               if(Debug) write(iulog,*)'CreateSubGridGraph: point #5.1 size(global2local) Svertex(i)%nbrs(j) ', &
! Subprogram not used                    size(global2local), Svertex(i)%nbrs(pos)
! Subprogram not used 
! Subprogram not used                  inbr = global2local(Svertex(i)%nbrs(pos))
! Subprogram not used 
! Subprogram not used                  if(Debug) write(iulog,*)'CreateSubGridGraph: point #5.2'
! Subprogram not used                  if(inbr .gt. 0) then 
! Subprogram not used                     new_pos = Svertex(i)%nbrs_ptr(j) + ncount
! Subprogram not used 
! Subprogram not used                     Svertex(i)%nbrs(new_pos) = inbr
! Subprogram not used                     Svertex(i)%nbrs_face(new_pos) = Svertex(i)%nbrs_face(pos)
! Subprogram not used                     Svertex(i)%nbrs_wgt(new_pos) = Svertex(i)%nbrs_wgt(pos)
! Subprogram not used                     Svertex(i)%nbrs_wgt_ghost(new_pos) = Svertex(i)%nbrs_wgt_ghost(pos)
! Subprogram not used                     ncount = ncount+1
! Subprogram not used                  endif
! Subprogram not used            enddo
! Subprogram not used            !set neighbors ptr
! Subprogram not used            orig_start =  Svertex(i)%nbrs_ptr(j+1);
! Subprogram not used            Svertex(i)%nbrs_ptr(j+1) =  Svertex(i)%nbrs_ptr(j) + ncount
! Subprogram not used            
! Subprogram not used 
! Subprogram not used            if(Debug) write(iulog,*)'CreateSubGridGraph: point #5.3'
! Subprogram not used         enddo !num_neighbors loop
! Subprogram not used 
! Subprogram not used 
! Subprogram not used         if(Debug) write(iulog,*)'CreateSubGridGraph: point #6'
! Subprogram not used         Svertex(i)%number = i
! Subprogram not used      enddo !nelem_s loop
! Subprogram not used      if(Debug) write(iulog,*)'CreateSubGridGraph: point #7'
! Subprogram not used      deallocate(global2local)
! Subprogram not used      if(Debug) write(iulog,*)'CreateSubGridGraph: point #8'
! Subprogram not used 
! Subprogram not used   end subroutine CreateSubGridGraph

!======================================================================

! Subprogram not used   subroutine PrintGridEdge(Edge)
! Subprogram not used 
! Subprogram not used     implicit none
! Subprogram not used     type (GridEdge_t), intent(in) :: Edge(:)
! Subprogram not used 
! Subprogram not used     integer           :: i,nedge,ii,wgtP
! Subprogram not used 
! Subprogram not used     nedge = SIZE(Edge)
! Subprogram not used 
! Subprogram not used     write(iulog,95)
! Subprogram not used     do i=1,nedge
! Subprogram not used           ii=Edge(i)%tail_face
! Subprogram not used 
! Subprogram not used           !map to correct location - for now all on same nbr side have same wgt, so take the first one
! Subprogram not used           ii = Edge(i)%tail%nbrs_ptr(ii)
! Subprogram not used 
! Subprogram not used           wgtP=Edge(i)%tail%nbrs_wgt(ii)       
! Subprogram not used           write(iulog,100) i, &
! Subprogram not used                Edge(i)%tail%number,Edge(i)%tail_face, wgtP, &
! Subprogram not used                Edge(i)%head%number,Edge(i)%head_face, gridedge_type(Edge(i))
! Subprogram not used     enddo
! Subprogram not used   95 format(5x,'GRIDEDGE #',3x,'Tail (face)',5x,'Head (face)',3x,'Type')
! Subprogram not used  100 format(10x,I6,8x,I4,1x,'(',I1,')  --',I2,'--> ',I6,1x,'(',I1,')',5x,'[',I1,']')
! Subprogram not used 
! Subprogram not used   end subroutine PrintGridEdge

!======================================================================
! ==========================================
! set_GridVertex_neighbors:
!
! Set global element number for element elem
! ==========================================

! Subprogram not used   subroutine set_GridVertex_number(elem,number)
! Subprogram not used 
! Subprogram not used     type(GridVertex_t)         :: elem
! Subprogram not used     integer                 :: number
! Subprogram not used 
! Subprogram not used     elem%number=number
! Subprogram not used 
! Subprogram not used   end subroutine set_GridVertex_number

!======================================================================
! Subprogram not used   subroutine PrintGridVertex(Vertex)
! Subprogram not used 
! Subprogram not used     implicit none 
! Subprogram not used     type (GridVertex_t), intent(in),target :: Vertex(:)
! Subprogram not used   
! Subprogram not used     integer        :: i,nvert
! Subprogram not used     integer ::n_west, n_east, n_south, n_north, n_swest, n_seast, n_nwest, n_neast
! Subprogram not used     integer ::w_west, w_east, w_south, w_north, w_swest, w_seast, w_nwest, w_neast
! Subprogram not used     integer ::n, print_buf(90), nbr(8), j, k, start, cnt, nbrs_cnt(8)
! Subprogram not used 
! Subprogram not used     integer, pointer :: np(:)
! Subprogram not used 
! Subprogram not used     nbr = (/ west, east, south, north, swest, seast, nwest, neast/)
! Subprogram not used 
! Subprogram not used     nvert = SIZE(Vertex)
! Subprogram not used         
! Subprogram not used     write(iulog,98)
! Subprogram not used     do i=1,nvert
! Subprogram not used 
! Subprogram not used        print_buf(:) = 0
! Subprogram not used        nbrs_cnt(:) = 0
! Subprogram not used        np =>  Vertex(i)%nbrs_ptr !alias
! Subprogram not used        cnt = 1
! Subprogram not used        do j = 1,num_neighbors
! Subprogram not used           n = np(nbr(j)+1) - np(nbr(j)) !num neigbors in that directions
! Subprogram not used           start =  np(nbr(j)) !start in array
! Subprogram not used           nbrs_cnt(j) = n
! Subprogram not used           do k = 1, n
! Subprogram not used              print_buf(cnt) = Vertex(i)%nbrs(start+k-1)
! Subprogram not used              print_buf(cnt+1) = Vertex(i)%nbrs_wgt(start+k-1)
! Subprogram not used              print_buf(cnt+2) = Vertex(i)%nbrs_face(start+k-1)
! Subprogram not used              cnt = cnt + 3
! Subprogram not used           end do
! Subprogram not used        enddo
! Subprogram not used 
! Subprogram not used        write(iulog,991) Vertex(i)%number, Vertex(i)%processor_number, &
! Subprogram not used             Vertex(i)%face_number, &
! Subprogram not used             print_buf(1:cnt-1)
! Subprogram not used 
! Subprogram not used        write(iulog,992) nbrs_cnt(1:8)
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     enddo
! Subprogram not used   98  format(5x,'GRIDVERTEX #',2x,'PART',2x,'DEG',4x,'W',9x,'E',9x, &
! Subprogram not used                 'S',9x,'N',9x,'SW',9x,'SE',9x,'NW',9x,'NE')
! Subprogram not used 
! Subprogram not used   991  format(10x,I3,8x,I4,8x,I4,2x,30(1x,I4,1x,'(',I2,I2,')'))
! Subprogram not used   992 format(30x,'nbrs_cnt:', 2x,8(1x,I4))
! Subprogram not used 
! Subprogram not used   end subroutine PrintGridVertex


!======================================================================

! Subprogram not used   subroutine CheckGridNeighbors(Vertex)
! Subprogram not used   
! Subprogram not used   implicit none
! Subprogram not used   type (GridVertex_t), intent(in) :: Vertex(:)
! Subprogram not used 
! Subprogram not used   integer :: i,j,k,l,m,nnbrs,inbrs,nvert
! Subprogram not used   nvert = SIZE(Vertex)
! Subprogram not used 
! Subprogram not used   do i=1,nvert
! Subprogram not used         nnbrs = SIZE(Vertex(i)%nbrs)
! Subprogram not used         do j=1,nnbrs
! Subprogram not used            inbrs = Vertex(i)%nbrs(j)
! Subprogram not used            if(inbrs > 0) then
! Subprogram not used               do k=1,nnbrs
! Subprogram not used                  if( inbrs .eq. Vertex(i)%nbrs(k) .and. (j/=k) ) &
! Subprogram not used                       write(iulog,*)'CheckGridNeighbors: ERROR identical neighbors detected  for Vertex ',i
! Subprogram not used                     
! Subprogram not used               enddo
! Subprogram not used            endif
! Subprogram not used         enddo
! Subprogram not used      enddo
! Subprogram not used 
! Subprogram not used   end subroutine CheckGridNeighbors

!======================================================================
  subroutine initgridedge(GridEdge,GridVertex)
  use parallel_mod, only : abortmp
  use dimensions_mod, only : max_corner_elem
  implicit none

  type (GridEdge_t), intent(inout)       :: GridEdge(:)
  type (GridVertex_t), intent(in),target :: GridVertex(:)

  integer                                :: i,j,k,iptr,m,n,wgtV,wgtP
  integer                                :: nelem,nelem_edge,inbr
  logical                                :: Verbose=.FALSE.
  integer                                :: mynbr_cnt, cnt, mystart, start

  nelem      = SIZE(GridVertex)
  nelem_edge = SIZE(GridEdge)

  GridEdge(:)%reverse=.FALSE.

  iptr=1
  do j=1,nelem
     do i=1,num_neighbors    
        mynbr_cnt = GridVertex(j)%nbrs_ptr(i+1) - GridVertex(j)%nbrs_ptr(i) !length of neighbor location  
        mystart = GridVertex(j)%nbrs_ptr(i) 
        do m=0,mynbr_cnt-1
           if((GridVertex(j)%nbrs_wgt(mystart + m) .gt. 0)) then    ! Do this only if has a non-zero weight
              if (nelem_edge<iptr) call abortmp('Error in initgridedge: Number of edges greater than expected.')
              GridEdge(iptr)%tail      => GridVertex(j)
              GridEdge(iptr)%tail_face =  mystart + m ! needs to be mystart + m (location in array)
              GridEdge(iptr)%tail_dir = i*max_corner_elem + m !conversion needed for setcycle
              inbr                     =  GridVertex(j)%nbrs(mystart+m)
              GridEdge(iptr)%head      => GridVertex(inbr)

              ! ===========================================
              ! Need this awful piece of code to determine
              ! which "face" of the neighbor element the
              ! edge links (i.e. the "head_face")
              ! ===========================================
              do k=1,num_neighbors
                 cnt = GridVertex(inbr)%nbrs_ptr(k+1) -GridVertex(inbr)%nbrs_ptr(k)                     
                 start = GridVertex(inbr)%nbrs_ptr(k)
                 do  n = 0, cnt-1
                    if(GridVertex(inbr)%nbrs(start+n) == GridVertex(j)%number) then
                       GridEdge(iptr)%head_face=start+n !needs to be start + n (location in array)
                       GridEdge(iptr)%head_dir=k*max_corner_elem+n !conversion (un-done in setcycle)
                    endif
                 enddo
              enddo
              iptr=iptr+1
           end if
        end do ! m loop
     end do !end i loop
  end do !end j loop
  if (nelem_edge+1 /= iptr) call abortmp('Error in initgridedge: Number of edges less than expected.')
  if (Verbose) then

     print *
     write(iulog,*)"element edge tail,head list: (TEST)"
     do i=1,nelem_edge
        write(iulog,*)GridEdge(i)%tail%number,GridEdge(i)%head%number
     end do

     print *
     write(iulog,*)"element edge tail_face, head_face list: (TEST)"
     do i=1,nelem_edge
        write(iulog,*)GridEdge(i)%tail_face,GridEdge(i)%head_face
     end do
  end if

  end subroutine initgridedge
!======================================================================

end module GridGraph_mod
