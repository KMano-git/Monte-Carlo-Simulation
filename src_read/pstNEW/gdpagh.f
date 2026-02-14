!***********************************************************************
      subroutine gdpagh(kk)
!***********************************************************************
      use com_gdhedr, only : hd_fno, hd_day, hd_dir
      implicit none

      integer  kk

!::local variables
      integer lp; data lp/0/; save lp
      real*4  x0, y0, y1, cfhh
      character ctx1*256, ctx2*256, lbnum*6
      integer ifig
!
!
      if( kk.eq.0 ) goto 100
      lp = lp + 1
!
      x0 = 15.0
      y0 = 28.0
      y1 = 27.2
      cfhh = 0.4
!
      call gdtbpo(14.0,18.2)
!
!::date
      ifig = hd_fno + lp
      ctx1 = trim(hd_day)//"   "//"Fig."
      call gdtbv1(trim(ctx1),real(ifig),'(f3.0)')
!
!::dir
      ctx2 = trim(hd_dir)
      call gdtbv1(trim(ctx2),99.0,'  ')
!
!::debug
      write(6,'(2x,"gdpagh ",a)') trim(ctx1)
      write(6,'(2x,"gdpagh ",a)') trim(ctx2)
!
 100  continue
      call gdpag(kk)
!
      return
      end

