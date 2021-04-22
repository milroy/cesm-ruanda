module cmparray_mod

  use shr_kind_mod, only : r8 => shr_kind_r8
  
  implicit none
  private
  save
  
  public expdaynite, cmpdaynite

  interface CmpDayNite
    module procedure CmpDayNite_1d_R
    module procedure CmpDayNite_2d_R
    module procedure CmpDayNite_3d_R
    module procedure CmpDayNite_1d_R_Copy
    module procedure CmpDayNite_2d_R_Copy
    module procedure CmpDayNite_3d_R_Copy
    module procedure CmpDayNite_1d_I
    module procedure CmpDayNite_2d_I
    module procedure CmpDayNite_3d_I
  end interface ! CmpDayNite

  interface ExpDayNite
    module procedure ExpDayNite_1d_R
    module procedure ExpDayNite_2d_R
    module procedure ExpDayNite_3d_R
    module procedure ExpDayNite_1d_I
    module procedure ExpDayNite_2d_I
    module procedure ExpDayNite_3d_I
  end interface ! ExpDayNite

  interface cmparray
    module procedure cmparray_1d_R
    module procedure cmparray_2d_R
    module procedure cmparray_3d_R
  end interface ! cmparray

  interface chksum
    module procedure chksum_1d_R
    module procedure chksum_2d_R
    module procedure chksum_3d_R
    module procedure chksum_1d_I
    module procedure chksum_2d_I
    module procedure chksum_3d_I
  end interface ! chksum

  contains

! Subprogram not used   subroutine CmpDayNite_1d_R(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1)
! Subprogram not used     integer, intent(in) :: Nday, Nnite
! Subprogram not used     integer, intent(in) :: il1, iu1
! Subprogram not used     integer, intent(in), dimension(Nday) :: IdxDay
! Subprogram not used     integer, intent(in), dimension(Nnite) :: IdxNite
! Subprogram not used     real(r8), intent(inout), dimension(il1:iu1) :: Array
! Subprogram not used 
! Subprogram not used     call CmpDayNite_3d_R(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, 1, 1, 1, 1)
! Subprogram not used 
! Subprogram not used     return
! Subprogram not used   end subroutine CmpDayNite_1d_R

! Subprogram not used   subroutine CmpDayNite_2d_R(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2, iu2)
! Subprogram not used     integer, intent(in) :: Nday, Nnite
! Subprogram not used     integer, intent(in) :: il1, iu1
! Subprogram not used     integer, intent(in) :: il2, iu2
! Subprogram not used     integer, intent(in), dimension(Nday) :: IdxDay
! Subprogram not used     integer, intent(in), dimension(Nnite) :: IdxNite
! Subprogram not used     real(r8), intent(inout), dimension(il1:iu1,il2:iu2) :: Array
! Subprogram not used 
! Subprogram not used     call CmpDayNite_3d_R(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2, iu2, 1, 1)
! Subprogram not used 
! Subprogram not used     return
! Subprogram not used   end subroutine CmpDayNite_2d_R

! Subprogram not used   subroutine CmpDayNite_3d_R(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2,iu2, il3, iu3)
! Subprogram not used     integer, intent(in) :: Nday, Nnite
! Subprogram not used     integer, intent(in) :: il1, iu1
! Subprogram not used     integer, intent(in) :: il2, iu2
! Subprogram not used     integer, intent(in) :: il3, iu3
! Subprogram not used     integer, intent(in), dimension(Nday) :: IdxDay
! Subprogram not used     integer, intent(in), dimension(Nnite) :: IdxNite
! Subprogram not used     real(r8), intent(inout), dimension(il1:iu1,il2:iu2,il3:iu3) :: Array
! Subprogram not used 
! Subprogram not used     real(r8), dimension(il1:iu1) :: tmp
! Subprogram not used     integer :: i, j, k
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     do k = il3, iu3
! Subprogram not used       do j = il2, iu2
! Subprogram not used 
! Subprogram not used         tmp(1:Nnite) = Array(IdxNite(1:Nnite),j,k)
! Subprogram not used         Array(il1:il1+Nday-1,j,k) = Array(IdxDay(1:Nday),j,k)
! Subprogram not used         Array(il1+Nday:il1+Nday+Nnite-1,j,k) = tmp(1:Nnite)
! Subprogram not used 
! Subprogram not used       end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     return
! Subprogram not used   end subroutine CmpDayNite_3d_R

  subroutine CmpDayNite_1d_R_Copy(InArray, OutArray, Nday, IdxDay, Nnite, IdxNite, il1, iu1)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    real(r8), intent(in), dimension(il1:iu1) :: InArray
    real(r8), intent(out), dimension(il1:iu1) :: OutArray

    call CmpDayNite_3d_R_Copy(InArray, OutArray, Nday, IdxDay, Nnite, IdxNite, il1, iu1, 1, 1, 1, 1)

    return
  end subroutine CmpDayNite_1d_R_Copy

  subroutine CmpDayNite_2d_R_Copy(InArray, OutArray, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2, iu2)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in) :: il2, iu2
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    real(r8), intent(in), dimension(il1:iu1,il2:iu2) :: InArray
    real(r8), intent(out), dimension(il1:iu1,il2:iu2) :: OutArray

    call CmpDayNite_3d_R_Copy(InArray, OutArray, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2, iu2, 1, 1)

    return
  end subroutine CmpDayNite_2d_R_Copy

  subroutine CmpDayNite_3d_R_Copy(InArray, OutArray, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2,iu2, il3, iu3)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in) :: il2, iu2
    integer, intent(in) :: il3, iu3
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    real(r8), intent(in), dimension(il1:iu1,il2:iu2,il3:iu3) :: InArray
    real(r8), intent(out), dimension(il1:iu1,il2:iu2,il3:iu3) :: OutArray

    integer :: i, j, k


    do k = il3, iu3
      do j = il2, iu2

         do i=il1,il1+Nday-1
            OutArray(i,j,k) = InArray(IdxDay(i-il1+1),j,k)
         enddo
         do i=il1+Nday,il1+Nday+Nnite-1
            OutArray(i,j,k) = InArray(IdxNite(i-(il1+Nday)+1),j,k)
         enddo
        

      end do
    end do

    return
  end subroutine CmpDayNite_3d_R_Copy

! Subprogram not used   subroutine CmpDayNite_1d_I(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1)
! Subprogram not used     integer, intent(in) :: Nday, Nnite
! Subprogram not used     integer, intent(in) :: il1, iu1
! Subprogram not used     integer, intent(in), dimension(Nday) :: IdxDay
! Subprogram not used     integer, intent(in), dimension(Nnite) :: IdxNite
! Subprogram not used     integer, intent(inout), dimension(il1:iu1) :: Array
! Subprogram not used 
! Subprogram not used     call CmpDayNite_3d_I(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, 1, 1, 1, 1)
! Subprogram not used 
! Subprogram not used     return
! Subprogram not used   end subroutine CmpDayNite_1d_I

! Subprogram not used   subroutine CmpDayNite_2d_I(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2, iu2)
! Subprogram not used     integer, intent(in) :: Nday, Nnite
! Subprogram not used     integer, intent(in) :: il1, iu1
! Subprogram not used     integer, intent(in) :: il2, iu2
! Subprogram not used     integer, intent(in), dimension(Nday) :: IdxDay
! Subprogram not used     integer, intent(in), dimension(Nnite) :: IdxNite
! Subprogram not used     integer, intent(inout), dimension(il1:iu1,il2:iu2) :: Array
! Subprogram not used 
! Subprogram not used     call CmpDayNite_3d_I(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2, iu2, 1, 1)
! Subprogram not used 
! Subprogram not used     return
! Subprogram not used   end subroutine CmpDayNite_2d_I

! Subprogram not used   subroutine CmpDayNite_3d_I(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2,iu2, il3, iu3)
! Subprogram not used     integer, intent(in) :: Nday, Nnite
! Subprogram not used     integer, intent(in) :: il1, iu1
! Subprogram not used     integer, intent(in) :: il2, iu2
! Subprogram not used     integer, intent(in) :: il3, iu3
! Subprogram not used     integer, intent(in), dimension(Nday) :: IdxDay
! Subprogram not used     integer, intent(in), dimension(Nnite) :: IdxNite
! Subprogram not used     integer, intent(inout), dimension(il1:iu1,il2:iu2,il3:iu3) :: Array
! Subprogram not used 
! Subprogram not used     integer, dimension(il1:iu1) :: tmp
! Subprogram not used     integer :: i, j, k
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     do k = il3, iu3
! Subprogram not used       do j = il2, iu2
! Subprogram not used 
! Subprogram not used         tmp(1:Nnite) = Array(IdxNite(1:Nnite),j,k)
! Subprogram not used         Array(il1:il1+Nday-1,j,k) = Array(IdxDay(1:Nday),j,k)
! Subprogram not used         Array(il1+Nday:il1+Nday+Nnite-1,j,k) = tmp(1:Nnite)
! Subprogram not used 
! Subprogram not used       end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     return
! Subprogram not used   end subroutine CmpDayNite_3d_I

  subroutine ExpDayNite_1d_R(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    real(r8), intent(inout), dimension(il1:iu1) :: Array

    call ExpDayNite_3d_R(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, 1, 1, 1, 1)

    return
  end subroutine ExpDayNite_1d_R

  subroutine ExpDayNite_2d_R(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2, iu2)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in) :: il2, iu2
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    real(r8), intent(inout), dimension(il1:iu1,il2:iu2) :: Array

    call ExpDayNite_3d_R(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2, iu2, 1, 1)

    return
  end subroutine ExpDayNite_2d_R

  subroutine ExpDayNite_3d_R(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2,iu2, il3, iu3)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in) :: il2, iu2
    integer, intent(in) :: il3, iu3
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    real(r8), intent(inout), dimension(il1:iu1,il2:iu2,il3:iu3) :: Array

    real(r8), dimension(il1:iu1) :: tmp
    integer :: i, j, k


    do k = il3, iu3
      do j = il2, iu2

        tmp(1:Nday) = Array(1:Nday,j,k)
        Array(IdxNite(1:Nnite),j,k) = Array(il1+Nday:il1+Nday+Nnite-1,j,k)
        Array(IdxDay(1:Nday),j,k) = tmp(1:Nday)

      end do
    end do

    return
  end subroutine ExpDayNite_3d_R

! Subprogram not used   subroutine ExpDayNite_1d_I(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1)
! Subprogram not used     integer, intent(in) :: Nday, Nnite
! Subprogram not used     integer, intent(in) :: il1, iu1
! Subprogram not used     integer, intent(in), dimension(Nday) :: IdxDay
! Subprogram not used     integer, intent(in), dimension(Nnite) :: IdxNite
! Subprogram not used     integer, intent(inout), dimension(il1:iu1) :: Array
! Subprogram not used 
! Subprogram not used     call ExpDayNite_3d_I(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, 1, 1, 1, 1)
! Subprogram not used 
! Subprogram not used     return
! Subprogram not used   end subroutine ExpDayNite_1d_I

! Subprogram not used   subroutine ExpDayNite_2d_I(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2, iu2)
! Subprogram not used     integer, intent(in) :: Nday, Nnite
! Subprogram not used     integer, intent(in) :: il1, iu1
! Subprogram not used     integer, intent(in) :: il2, iu2
! Subprogram not used     integer, intent(in), dimension(Nday) :: IdxDay
! Subprogram not used     integer, intent(in), dimension(Nnite) :: IdxNite
! Subprogram not used     integer, intent(inout), dimension(il1:iu1,il2:iu2) :: Array
! Subprogram not used 
! Subprogram not used     call ExpDayNite_3d_I(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2, iu2, 1, 1)
! Subprogram not used 
! Subprogram not used     return
! Subprogram not used   end subroutine ExpDayNite_2d_I

! Subprogram not used   subroutine ExpDayNite_3d_I(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2,iu2, il3, iu3)
! Subprogram not used     integer, intent(in) :: Nday, Nnite
! Subprogram not used     integer, intent(in) :: il1, iu1
! Subprogram not used     integer, intent(in) :: il2, iu2
! Subprogram not used     integer, intent(in) :: il3, iu3
! Subprogram not used     integer, intent(in), dimension(Nday) :: IdxDay
! Subprogram not used     integer, intent(in), dimension(Nnite) :: IdxNite
! Subprogram not used     integer, intent(inout), dimension(il1:iu1,il2:iu2,il3:iu3) :: Array
! Subprogram not used 
! Subprogram not used     integer, dimension(il1:iu1) :: tmp
! Subprogram not used     integer :: i, j, k
! Subprogram not used 
! Subprogram not used 
! Subprogram not used     do k = il3, iu3
! Subprogram not used       do j = il2, iu2
! Subprogram not used 
! Subprogram not used         tmp(1:Nday) = Array(1:Nday,j,k)
! Subprogram not used         Array(IdxNite(1:Nnite),j,k) = Array(il1+Nday:il1+Nday+Nnite-1,j,k)
! Subprogram not used         Array(IdxDay(1:Nday),j,k) = tmp(1:Nday)
! Subprogram not used 
! Subprogram not used       end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     return
! Subprogram not used   end subroutine ExpDayNite_3d_I

!******************************************************************************!
!                                                                              !
!                                 DEBUG                                        !
!                                                                              !
!******************************************************************************!

! Subprogram not used   subroutine cmparray_1d_R(name, Ref, New, id1, is1, ie1)
! Subprogram not used     character(*), intent(in) :: name
! Subprogram not used     integer,  intent(in) :: id1, is1, ie1
! Subprogram not used     real(r8), intent(in), dimension(id1) :: Ref
! Subprogram not used     real(r8), intent(in), dimension(id1) :: New
! Subprogram not used 
! Subprogram not used     call cmparray_3d_R(name, Ref, New, id1, is1, ie1, 1, 1, 1, 1, 1, 1)
! Subprogram not used   end subroutine cmparray_1d_R

! Subprogram not used   subroutine cmparray_2d_R(name, Ref, New, id1, is1, ie1, id2, is2, ie2)
! Subprogram not used     character(*), intent(in) :: name
! Subprogram not used     integer,  intent(in) :: id1, is1, ie1
! Subprogram not used     integer,  intent(in) :: id2, is2, ie2
! Subprogram not used     real(r8), intent(in), dimension(id1, id2) :: Ref
! Subprogram not used     real(r8), intent(in), dimension(id1, id2) :: New
! Subprogram not used 
! Subprogram not used     call cmparray_3d_R(name, Ref, New, id1, is1, ie1, id2, is2, ie2, 1, 1, 1)
! Subprogram not used   end subroutine cmparray_2d_R

! Subprogram not used   subroutine cmparray_3d_R(name, Ref, New, id1, is1, ie1, id2, is2, ie2, id3, is3, ie3)
! Subprogram not used     character(*), intent(in) :: name
! Subprogram not used     integer,  intent(in) :: id1, is1, ie1
! Subprogram not used     integer,  intent(in) :: id2, is2, ie2
! Subprogram not used     integer,  intent(in) :: id3, is3, ie3
! Subprogram not used     real(r8), intent(in), dimension(id1, id2, id3) :: Ref
! Subprogram not used     real(r8), intent(in), dimension(id1, id2, id3) :: New
! Subprogram not used 
! Subprogram not used     integer :: i, j, k
! Subprogram not used     integer :: nerr
! Subprogram not used     logical :: found
! Subprogram not used     real(r8):: rdiff
! Subprogram not used     real(r8), parameter :: rtol = 1.0e-13_r8
! Subprogram not used 
! Subprogram not used     nerr = 0
! Subprogram not used 
! Subprogram not used     do k = is3, ie3
! Subprogram not used       do j = is2, ie2
! Subprogram not used 
! Subprogram not used         found = .false.
! Subprogram not used         do i = is1, ie1
! Subprogram not used           rdiff = abs(New(i,j,k)-Ref(i,j,k))
! Subprogram not used           rdiff = rdiff / merge(abs(Ref(i,j,k)), 1.0_r8, Ref(i,j,k) /= 0.0_r8)
! Subprogram not used           if ( rdiff > rtol ) then
! Subprogram not used             found = .true.
! Subprogram not used             exit
! Subprogram not used           end if
! Subprogram not used         end do
! Subprogram not used 
! Subprogram not used         if ( found ) then
! Subprogram not used           do i = is1, ie1
! Subprogram not used             rdiff = abs(New(i,j,k)-Ref(i,j,k))
! Subprogram not used             rdiff = rdiff / merge(abs(Ref(i,j,k)), 1.0_r8, Ref(i,j,k) /= 0.0_r8)
! Subprogram not used             if ( rdiff > rtol ) then
! Subprogram not used               print 666, name, i, j, k, Ref(i, j, k), New(i, j, k), rdiff
! Subprogram not used               nerr = nerr + 1
! Subprogram not used               if ( nerr > 10 ) stop
! Subprogram not used             end if
! Subprogram not used           end do
! Subprogram not used         end if
! Subprogram not used 
! Subprogram not used       end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     return
! Subprogram not used 666 format('cmp3d: ', a10, 3(1x, i4), 3(1x, e20.14))
! Subprogram not used 
! Subprogram not used   end subroutine cmparray_3d_R

! Subprogram not used   subroutine chksum_1d_R(name, Ref, id1, is1, ie1)
! Subprogram not used     character(*), intent(in) :: name
! Subprogram not used     integer,  intent(in) :: id1, is1, ie1
! Subprogram not used     real(r8), intent(in), dimension(id1) :: Ref
! Subprogram not used 
! Subprogram not used     call chksum_3d_R(name, Ref, id1, is1, ie1, 1, 1, 1, 1, 1, 1)
! Subprogram not used   end subroutine chksum_1d_R

! Subprogram not used   subroutine chksum_1d_I(name, Ref, id1, is1, ie1)
! Subprogram not used     character(*), intent(in) :: name
! Subprogram not used     integer,  intent(in) :: id1, is1, ie1
! Subprogram not used     integer,  intent(in), dimension(id1) :: Ref
! Subprogram not used 
! Subprogram not used     call chksum_3d_I(name, Ref, id1, is1, ie1, 1, 1, 1, 1, 1, 1)
! Subprogram not used   end subroutine chksum_1d_I

! Subprogram not used   subroutine chksum_2d_R(name, Ref, id1, is1, ie1, id2, is2, ie2)
! Subprogram not used     character(*), intent(in) :: name
! Subprogram not used     integer,  intent(in) :: id1, is1, ie1
! Subprogram not used     integer,  intent(in) :: id2, is2, ie2
! Subprogram not used     real(r8), intent(in), dimension(id1, id2) :: Ref
! Subprogram not used 
! Subprogram not used     call chksum_3d_R(name, Ref, id1, is1, ie1, id2, is2, ie2, 1, 1, 1)
! Subprogram not used   end subroutine chksum_2d_R

! Subprogram not used   subroutine chksum_2d_I(name, Ref, id1, is1, ie1, id2, is2, ie2)
! Subprogram not used     character(*), intent(in) :: name
! Subprogram not used     integer,  intent(in) :: id1, is1, ie1
! Subprogram not used     integer,  intent(in) :: id2, is2, ie2
! Subprogram not used     integer,  intent(in), dimension(id1, id2) :: Ref
! Subprogram not used 
! Subprogram not used     call chksum_3d_I(name, Ref, id1, is1, ie1, id2, is2, ie2, 1, 1, 1)
! Subprogram not used   end subroutine chksum_2d_I

! Subprogram not used   subroutine chksum_3d_R(name, Ref, id1, is1, ie1, id2, is2, ie2, id3, is3, ie3)
! Subprogram not used     character(*), intent(in) :: name
! Subprogram not used     integer,  intent(in) :: id1, is1, ie1
! Subprogram not used     integer,  intent(in) :: id2, is2, ie2
! Subprogram not used     integer,  intent(in) :: id3, is3, ie3
! Subprogram not used !orig    real(r8), intent(in), dimension(id1, id2, id3) :: Ref
! Subprogram not used     real(r8), intent(in), dimension(is1:ie1, is2:ie2, is3:ie3) :: Ref
! Subprogram not used 
! Subprogram not used     real(r8) :: chksum
! Subprogram not used     real(r8) :: rmin, rmax
! Subprogram not used     integer :: i, j, k
! Subprogram not used     integer :: imin, jmin, kmin
! Subprogram not used     integer :: imax, jmax, kmax
! Subprogram not used 
! Subprogram not used     imin = is1 ; jmin = is2 ; kmin = is3
! Subprogram not used     imax = is1 ; jmax = is2 ; kmax = is3
! Subprogram not used     rmin = Ref(is1, is2, is3) ; rmax = rmin
! Subprogram not used 
! Subprogram not used     chksum = 0.0_r8
! Subprogram not used 
! Subprogram not used     do k = is3, ie3
! Subprogram not used       do j = is2, ie2
! Subprogram not used         do i = is1, ie1
! Subprogram not used           chksum = chksum + abs(Ref(i,j,k))
! Subprogram not used           if ( Ref(i,j,k) < rmin ) then
! Subprogram not used             rmin = Ref(i,j,k)
! Subprogram not used             imin = i ; jmin = j ; kmin = k
! Subprogram not used           end if
! Subprogram not used           if ( Ref(i,j,k) > rmax ) then
! Subprogram not used             rmax = Ref(i,j,k)
! Subprogram not used             imax = i ; jmax = j ; kmax = k
! Subprogram not used           end if
! Subprogram not used         end do
! Subprogram not used       end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     print 666, name, chksum, imin, jmin, kmin, imax, jmax, kmax
! Subprogram not used 666 format('chksum: ', a8, 1x, e20.14, 6(1x, i4))
! Subprogram not used 
! Subprogram not used   end subroutine chksum_3d_R

! Subprogram not used   subroutine chksum_3d_I(name, Ref, id1, is1, ie1, id2, is2, ie2, id3, is3, ie3)
! Subprogram not used     character(*), intent(in) :: name
! Subprogram not used     integer,  intent(in) :: id1, is1, ie1
! Subprogram not used     integer,  intent(in) :: id2, is2, ie2
! Subprogram not used     integer,  intent(in) :: id3, is3, ie3
! Subprogram not used     integer,  intent(in), dimension(id1, id2, id3) :: Ref
! Subprogram not used 
! Subprogram not used     integer :: i, j, k
! Subprogram not used     integer :: chksum
! Subprogram not used     chksum = 0
! Subprogram not used 
! Subprogram not used     do k = is3, ie3
! Subprogram not used       do j = is2, ie2
! Subprogram not used         do i = is1, ie1
! Subprogram not used           chksum = chksum + abs(Ref(i,j,k))
! Subprogram not used         end do
! Subprogram not used       end do
! Subprogram not used     end do
! Subprogram not used 
! Subprogram not used     print 666, name, chksum
! Subprogram not used 666 format('chksum: ', a8, 1x, i8)
! Subprogram not used 
! Subprogram not used   end subroutine chksum_3d_I

end module cmparray_mod
