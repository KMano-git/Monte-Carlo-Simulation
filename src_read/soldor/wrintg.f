!***********************************************************************
      subroutine wrintg(wrd,wrsm,wrgn)
!***********************************************************************
!
!      wrd(ndmc)  In
!      wrsm = Integ(div+sol+prv+edge) wrgn(irg)
!      wrgn = region
!
!-----------------------------------------------------------------------
      use cntcom, only : mrgn, ncmax, volm
      use csize,  only : ndmc
      implicit none
!
!::argument
! modified 3/3 lines organize local variables and include files by kamata 2021/05/31
!ik   real*8, dimension(ndmc) :: wrd
!ik   real*8 :: wrsm
!ik   real*8, dimension(10)   :: wrgn
      real(8), dimension(ndmc), intent(in)  :: wrd
      real(8),                  intent(out) :: wrsm
      real(8), dimension(10),   intent(out) :: wrgn
!
!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer :: ic, j, i, ir
      integer :: ic, ir
!
      wrgn(1:10) = 0.0d0
      wrsm = 0.0d0
!
      do ic = 1, ncmax
! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik   j = iplx(ic)
!ik   i = iply(ic)
      ir  = mrgn(ic)
!xx   ir2 = mrgnp(ic)  ! <== Not defined
      if( ir.gt.0 ) then
      wrgn(ir) = wrgn(ir) + wrd(ic)*volm(ic)
      endif
      enddo
!
      wrsm  = 0.0d0
      do ir = 1, 6  !  6:core edge 7:hot core
      wrsm  = wrsm  + wrgn(ir)
      enddo
!
!xx      write(n6,'(2x,"wrdisk wrsm =",1pe12.3)') wrsm
      return
      end
!
!***********************************************************************
      subroutine wrdisk
!***********************************************************************
      use cntcom, only : iplx, iply, ncmax
      use cplcom, only : wrad
      use cplwrd, only : wcre, wime, wimi, wmce, wmci
      use csize,  only : ndmc
      use csonic, only : limp
      use cunit,  only : n6
      implicit none
!
!::local variables
      integer   ift
      character dsn*80, dsn2*80
      integer   ir
      logical   lex
      integer :: ii, ic, j, i
!
      real*8, dimension(10) :: hwrgn
      real*8 ::  hwrsm
!
      ift = 221
      dsn = "wrad.dat"
!
!::KH111108  non-coronal model
      if( limp.eq.0 ) then
      wcre(1:ndmc) = 0.0d0
      wmci(1:ndmc) = 0.0d0
      wmce(1:ndmc) = 0.0d0
      wimi(1:ndmc) = 0.0d0
      wime(1:ndmc) = 0.0d0
      ii = 0
      do ic = 1, ncmax
        j = iplx(ic)
        i = iply(ic)
        if( j.le.0 .or. i.le.0 ) cycle
        ii = ii + 1
        wcre(ic) = wrad(j,i)
        wime(ic) = wrad(j,i)
      enddo
      write(n6,'(2x,"sub. wrdisk set wcre and wime in limp=0 case")')
      write(n6,'(2x,"ii =",i6,"  ncmax =",i6)') ii, ncmax
      endif
!
!::wrad-data       CR-e  MC-i  MC-e  Tot-i Tot-e
      call nopen(ift,dsn,"binary write",dsn2,lex)
      write(ift)  wcre, wmci, wmce, wimi, wime
      close(ift)
!
!::check
      call wrintg(wime,hwrsm,hwrgn)
!
      write(n6,'(1x,a,2x,1pe14.6,1x,1p12e11.3)')
     >   "wrdisk: cradTT",hwrsm,(hwrgn(ir),ir=1,7)
!
      return
      end
