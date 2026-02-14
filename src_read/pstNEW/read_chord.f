!***********************************************************************
      subroutine read_chord(unitno, capacity, xst, yst, xen, yen,
     >                            nlines, error)
!***********************************************************************
      use cunit
      implicit none

      integer, intent(in)  :: unitno
      integer, intent(in)  :: capacity
      real*8, intent(out)  :: xst(capacity)
      real*8, intent(out)  :: yst(capacity)
      real*8, intent(out)  :: xen(capacity)
      real*8, intent(out)  :: yen(capacity)
      integer, intent(out) :: nlines 
      integer, intent(out), optional :: error

      integer :: stat ! io status
      integer :: nl
      integer :: no   ! line No. read,  dummy
      integer :: i
      real*8 :: sline(4,100)

      namelist /sight_lines/ nlines, sline

!     write(n6,*) '****  read chord (sight lines) ****'

      if (present(error)) error=0
      nlines=0

      read(unit=unitno, nml=sight_lines,
     >         end=100, err=200, iostat=stat)

      do i=1, nlines
         xst(i) = sline(1,i)
         yst(i) = sline(2,i)
         xen(i) = sline(3,i)
         yen(i) = sline(4,i)
!         write(6,*) i, xst(i), yst(i), xen(i), yen(i)
      end do

      return

!
  100     continue
      write(n6, '("read_chord: EOF (nlines); stat=")') stat
      if (present(error)) error=stat
      return
!
  200     continue
      write(n6, '("read_chord: ERROR (nlines): stat=")') stat
      if (present(error)) error=stat
      return
!
  300     continue
      write(n6, '("read_chord: EOF (points) stat=")') stat
      if (present(error)) error=stat
      return
!
  400     continue
      write(n6, '("read_chord: ERROR (points) stat=")') stat
      if (present(error)) error=stat
      return
      end subroutine read_chord

