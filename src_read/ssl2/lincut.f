!*********************************************************************
!xxx      program test_lincut
      subroutine test_lincut
!*********************************************************************
      implicit none
      character :: clin*80
      integer :: msta

      clin = "* formatI   (2i5)  "
      call lincut(clin,"formatIX",msta)
      write(6,'(2x,"[",a,"]")') trim(clin)
      if( msta > 0 ) then
      write(6,'(2x,"[",a,"]")') trim(clin(msta:))
      endif

      stop
      end

!*********************************************************************
      subroutine lincut(clin,ched,msta)
!*********************************************************************
      implicit none

!::arguments
      character,intent(in) :: clin*(*), ched*(*)
      integer,  intent(out):: msta

!::local varaibales
      integer :: mj, j

      mj = index(clin,ched)

      if(mj == 0 ) then
        msta = 0
        return
      endif

      mj = mj + len(ched)
      do j = mj, len(clin)
        if( clin(j:j) /= " " ) then
          msta = j
          goto 200
        endif
      enddo

      msta = 0
      return

 200  continue
      return
      end

