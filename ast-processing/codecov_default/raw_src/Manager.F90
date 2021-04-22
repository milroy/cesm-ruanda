!------------------------------------------------------------------------------
! High Order Method Modeling Environment (HOMME)
!------------------------------------------------------------------------------
!
! MODULE: Manager
!
!
! DESCRIPTION:
!> @brief The Manager controls global shared resources in applications using HOMME.
!!
!! This module is designed after the <b>singleton</b> pattern for Object Oriented programming.
!! See \a Design \a Patterns by Gamma, Helm, Johnson and Vlissides.
!! The idea is to avoid passing resources from high level routines in HOMME to lower
!! level routines via subroutines parameters. \n
!! For example, if one has an array of 'cell'
!! which is utilized everywhere in HOMME, rather than forcing the arguments of the routines
!! to take a 'cell' argument and passed forward like this
!! \code{.F90}
!! subroutine higher(..., cell my_cell, ...)
!!    ...some code  doing something...
!!          \! Pass along the argument cell
!!          call lower(.., my_cell, ...)
!! end subroutine higher
!! \endcode
!! the function at the higher level can avoid the parameter my_cell all together and instead the function
!! at the lower level, the one that really needs the resource obtains the resource like this
!! \code{.F90}
!! subroutine lower(...)
!!   use Manager, only: cell_get
!!    ...some code  doing something...
!!          \! get the cell resource
!!          call cell_get(cell)
!! end subroutine lower
!! \endcode
!! Notice the concept of a resource is an entity in the program that is unique and that can be obtained with
!! a handle. In Fortran 90 we implement this idea via pointers to the resource, so a caller must create
!! a pointer outside, nullify it, and call the appropiatte subroutine to get its
!! own pointer assigned to the resource. For example,
!! \code{.F90}
!! subroutine InitColumnModel(elem, cm,hvcoord,hybrid,tl,nets,nete,runtype)
!!
!!    use Manager
!!    ...
!!    nullify(elem_physics)
!!    call element_physics_get(elem_physics)
!!    ...
!! end subroutine InitColumnModel
!! \endcode
!! Important
!! 1. Clearly this is not thread safe by default unless you use a mutex to control the resource. \n
!! 2. The naming convention for the subroutines to "get" the resources is \<resource name\>\_get. \n
!! 3. If a pointer is not null when using a get routine, the program will abort because we can't override
!!    a valid address.
!! \n
!! \n
!> @author Jose Garcia (jgarcia@ucar.edu). NCAR
!------------------------------------------------------------------------------
module Manager

! In the absence of a clear alternative, I have removed this module
! and all uses of it from 1 compilations. -Sean Santos (santos@ucar.edu)
end module Manager
