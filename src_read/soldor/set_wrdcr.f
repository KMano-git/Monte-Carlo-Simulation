!xx      call set_wrdcr
!xx      stop
!xx      end

!**********************************************************************
      subroutine set_wrdcr
!**********************************************************************
!
!  wrdcr  wrd : radiation  cr : corona model
!
!
! &uinpt
!  wcr_nty = 3,
!  wcr_typ(1)="W2", wcr_cnc(1:6,1)= 1.0d-5,1.0d-5,1.0d-5,1.0d-5,1.0d-5,1.0d-5,
!  wcr_typ(2)="Ne", wcr_cnc(1:6,2)= 1.0d-2,0.5d-2,1.0d-2,0.2d-2,0.4d-2,2.0d-2,
!  wcr_typ(3)="Ar", wcr_cnc(1:6,3)= 1.0d-4,0.5d-4,1.0d-4,0.2d-4,0.4d-4,2.0d-4,
! &end
!
!---------------------------------------------------------------------
      use cplcom, only : cimp
      use cplwrd, only : ndprg, prg_typ, wcr_cnc, wcr_ity, wcr_nty
     >    , wcr_typ
      use cunit,  only : n6
      implicit none

!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer :: ii, i, j, k, ir
      integer :: i, j, k, ir
      character(4) :: typ

      write(n6,'(/2x,"*** set_wrdcr ***")')

!::impurity speces to be programed   2015/06/12
!xx  prg_typ = (/"C", "W", "W2", "Ar", "Ne"/)

      write(n6,'(2x,"programing for imputity  ",10(a,4x))')
     >  (trim(prg_typ(i)),i=1,ndprg)

!::define ity
      do i = 1, wcr_nty
        typ = wcr_typ(i)
        do j = 1, ndprg
          if( trim(typ) == trim(prg_typ(j)) ) then
            wcr_ity(i) = j
            goto 120
          endif
        enddo
        wcr_ity(i) = 0
 120    continue
      enddo

!::debug write
      write(n6,'(2x,"wcr_nty =",i2)') wcr_nty
      do i = 1, wcr_nty
        k = wcr_ity(i)
        if( k > 0 ) then
        write(n6,'(2x,i2,2x,a,1p6e12.3,2x,":",i2,2x,a)')
     >    i, wcr_typ(i), (wcr_cnc(ir,i),ir =1,6), k, prg_typ(k)
       else
       write(n6,'(2x,i2,2x,a,72x,2x,":",i2)') i, wcr_typ(i), k
       endif
      enddo

!::cimp  max value of wcr_cnc
      cimp = -1.0d20
      do i = 1, wcr_nty
        do ir = 1, 6
          cimp = dmax1( cimp, wcr_cnc(ir,i) )
        enddo
      enddo
      write(n6,'(2x,"cimp =",1pe12.3)') cimp

      return
      end
