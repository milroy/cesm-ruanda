



module ll_mod
  implicit none
  private
  type :: node_t
     integer :: id
     integer :: Src,Dest
     logical :: valid
     type(node_t), pointer :: prev,next
  end type node_t

  type :: root_t
     integer     :: number
     type(node_t), pointer :: first
  end type root_t
  public :: node_t, root_t
  integer, public :: NumEdges 

  public :: PrintEdgeList
  public :: LLAddEdge,LLFindEdge, LLInsertEdge
  public :: LLSetEdgeCount,LLGetEdgeCount
  public :: LLFree

  interface assignment(=) 
     module procedure copy_node
  end interface

contains 

  subroutine LLSetEdgeCount(value)
    implicit none
    integer,intent(in)   :: value
    NumEdges=value
  end subroutine LLSetEdgeCount

  subroutine LLGetEdgeCount(value)
    implicit none
    integer,intent(out)  :: value
    value=NumEdges
  end subroutine LLGetEdgeCount

! Subprogram not used   recursive subroutine copy_node(node2,node1)
! Subprogram not used 
! Subprogram not used     type (node_t), intent(out) :: node2
! Subprogram not used     type (node_t), intent(in)  :: node1
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     node2%Src    = node1%Src
! Subprogram not used     node2%Dest    = node1%Dest
! Subprogram not used     node2%id      = node1%id
! Subprogram not used     node2%valid    = node1%valid
! Subprogram not used     node2%prev       = node1%prev
! Subprogram not used     node2%next       = node1%next
! Subprogram not used 
! Subprogram not used   end subroutine copy_node

! Subprogram not used   subroutine PrintEdgeList(EdgeList)
! Subprogram not used 
! Subprogram not used     type(root_t) :: EdgeList(:)
! Subprogram not used     type(node_t), pointer :: temp_node
! Subprogram not used     integer :: nlist, i
! Subprogram not used     nlist = SIZE(EdgeList)
! Subprogram not used 
! Subprogram not used     do i=1,nlist
! Subprogram not used        temp_node => EdgeList(i)%first
! Subprogram not used        do while(associated(temp_node)) 
! Subprogram not used           print *,'Vertex: ',EdgeList(i)%number ,temp_node%Src,'->' ,temp_node%dest, '(',temp_node%id,')'
! Subprogram not used           temp_node => temp_node%next
! Subprogram not used        enddo
! Subprogram not used     enddo
! Subprogram not used 
! Subprogram not used   end subroutine PrintEdgeList

  subroutine LLFree(List)

    implicit none
    type(root_t) :: List
    type(node_t), pointer :: temp_node
    integer :: nlist,i


    temp_node => List%first
    ! Find the end of the list
    do while(associated(temp_node%next))
       temp_node => temp_node%next
    enddo

    temp_node => temp_node%prev
    !Now step back and deallocate all entries  
    do while(associated(temp_node))
       deallocate(temp_node%next)
       temp_node => temp_node%prev
    enddo

  end subroutine LLFree

  subroutine LLInsertEdge(EdgeList,src,dest,eNum)
    type (root_t), intent(inout) :: EdgeList
    integer, intent(in) :: src,dest
    integer, intent(out) :: eNum
    logical :: found

    call LLFindEdge(EdgeList,src,dest,eNum,found) 
    if(.not. found) then 
       call LLAddEdge(EdgeList,src,dest,eNum) 
    endif

  end subroutine LLInsertEdge

  subroutine LLFindEdge(Edge,src,dest,id,found)

    type (root_t), intent(in) :: Edge
    integer, intent(in)  :: src,dest
    integer, intent(out) :: id
    logical, intent(out) :: found

    type (node_t), pointer :: temp_node

    found =.FALSE.

    temp_node => Edge%first
    do while(associated(temp_node) .and. (.not. found))
       if((dest .eq. temp_node%dest) .and. (src .eq. temp_node%Src) ) then 
          found = .TRUE. 
          id=temp_node%id
       else
          temp_node => temp_node%next
       endif
    enddo
  end subroutine LLFindEdge

  subroutine LLAddEdge(EdgeList,src,dest,id)
    type (root_t), intent(inout) :: EdgeList
    integer, intent(in)  :: src
    integer, intent(in)  :: dest
    integer, intent(out)  :: id

    type(node_t), pointer :: temp_node
    type(node_t), pointer  :: new_node
    type(node_t), pointer :: parent

    temp_node => EdgeList%first
    parent    => EdgeList%first

    do while(associated(temp_node))
       parent => temp_node
       temp_node => parent%next
    enddo
    allocate(new_node)
    NumEdges = NumEdges + 1

    new_node%src=src
    new_node%dest=dest
    new_node%id=NumEdges
    NULLIFY(new_node%next)
    new_node%prev => parent

    if(associated(EdgeList%first)) then
       parent%next => new_node 
    else
       EdgeList%first => new_node 
    endif
    id = NumEdges

  end subroutine LLAddEdge

end module ll_mod
