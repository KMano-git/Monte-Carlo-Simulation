!***********************************************************************
      subroutine selnty(wmc_nty,nty)
!***********************************************************************
      implicit none

      integer wmc_nty, nty,mji
      character*80 cmsg,cinp

      nty = 1
      if(wmc_nty .gt. 1) then

      write(cmsg,'(a,i2,a)')
     >     ">> select impurity spicies (1<rtn>:",wmc_nty,") ==> "
      call gdget(cmsg,cinp,mji)

      if(mji.gt.0)then
        read(cinp,*)nty
        if(nty.gt.wmc_nty) call wexit("select_nty","imvalid nty")
      endif
      call gdprint(nty)

      endif

      return
      end



