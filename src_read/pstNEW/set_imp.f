!**********************************************************************
      subroutine set_imp
!***********************************************************************
!
!    read(trim(cinp),*) wfinp   error at iferc   by K. Hoshino
!     ==>  read(cinp,*) wfinp       KH120912
!
!      after call set_imp   tdnz, wrad  absolute value  2015/04/07
!                           wmce, wime  absolute value  2015/05/12
!---------------------------------------------------------------------
      use csize
      use cimden
      use cntcom
      use cplimp
      use cplimp_plimp
      use cplwrd
      use csonic
      use cunit
      use cimset
      implicit none

!
!::local variables
      integer :: nskp, ir, iz, ic, mji, nty
      real(8) :: hwrsm
      real(8), dimension(10) :: hwrgn
      character(80) :: cinp
      character(160) :: cmsg
!
      nskp = 500
!
      write(n6,'(/2x,"*** set_imp ***")')
!
      do nty = 1, wmc_nty
        write(n6,'(2x,"nty ="i3)') nty
        write(n6,'(2x,"  ncmax,ncmax2 =",2i6,"  nzmx =",i3)')
     >      ncmax, ncmax2, nzmxL(nty)
!
        write(n6,'(2x)')
        write(n6,'(6x,"ic",2x,"icx",1x,"icy",1x,"ir",2x,"twrd",8x,
     >      6("tdnz_",i1,6x))') (iz,iz=0,5)
        do ic = nskp, ncmax, nskp
        ir = mrgn(ic)
        write(n6,'(2x,i6,i5,i4,i3,1p10e12.4)')
     >    ic, iplx(ic), iply(ic), ir, twrdL(ic,nty), 
     >    (tdnzL(iz,ic,nty),iz=0,5)
        enddo !ic
      enddo !nty
!
      cmsg = "  "
      write(cmsg,'("wfac =",f8.4)') wfac
      write(n6,'(2x)')
      write(n6,'(2x,"wfac =",f8.4)') wfac
      call gdput(trim(cmsg))
c     call gdget(">> enter wfac (2.71/def:<rtn>) ==> ",cinp,mji)
c     if( mji.eq.0 ) cinp = "0.0"
c     read(cinp,*) wfinp
      if( wfinp.le.0.0 ) wfinp = wfac
c     write(cmsg,'("enhanced twrd & tdnz by wfinp =",1pe12.4)') wfinp
c     call gdput(trim(cmsg))
c     write(n6,'(2x,a)') trim(cmsg)
!
!::debug
c     call wrintg(twrd,hwrsm,hwrgn)
c     write(n6,602)   " twrd",hwrsm, (hwrgn(ir),ir=1,7)
c     write(cmsg,602) " twrd",hwrsm, (hwrgn(ir),ir=1,7)
c     call gdput(trim(cmsg))
 602  format(1x,a,2x,1pe14.6,1x,1p12e11.3)
!
!::absolute value Wrad and density
      do ic = 1, ncmax2
        twrdL(ic,1) = wfinp*twrdL(ic,1)
        do iz = 0, nzmxL(1)
          tdnzL(iz,ic,1) = wfinp*tdnzL(iz,ic,1)
          tionZL(iz,ic,1) = wfinp*tionZL(iz,ic,1)
          trecZL(iz,ic,1) = wfinp*trecZL(iz,ic,1)
! delete 3 lines unused
!          tradiZL(iz,ic,1) = wfinp*tradiZL(iz,ic,1)
!          tradliZL(iz,ic,1) = wfinp*tradliZL(iz,ic,1)
!          tradrZL(iz,ic,1) = wfinp*tradrZL(iz,ic,1)
        enddo
      enddo
!
!::IMPMC rad (include wfac)
c     call wrintg(twrd,hwrsm,hwrgn)
c     write(n6,602)   "Ftwrd",hwrsm, (hwrgn(ir),ir=1,7)
c     write(cmsg,602) "Ftwrd",hwrsm, (hwrgn(ir),ir=1,7)
c     call gdput(trim(cmsg))
!
!::non corona
c     call wrintg(wcre,hwrsm,hwrgn)
c     write(n6,602)   "-wcre",-hwrsm, (-hwrgn(ir),ir=1,7)
c     write(cmsg,602) "-wcre",-hwrsm, (-hwrgn(ir),ir=1,7)
c     call gdput(trim(cmsg))
!
!::IMPMC rad
c     do ic = 1, ncmax2
c       wmce(ic) = -twrd(ic)
c     enddo
c     call wrintg(wmce,hwrsm,hwrgn)
c     write(n6,602)   "-wmce",-hwrsm, (-hwrgn(ir),ir=1,7)
c     write(cmsg,602) "-wmce",-hwrsm, (-hwrgn(ir),ir=1,7)
c     call gdput(trim(cmsg))
!
!::non-corona + IMPMC rad
c     do ic = 1, ncmax2
c     wime(ic) = wcre(ic) + wmce(ic)
c     enddo
c     call wrintg(wime,hwrsm,hwrgn)
c     write(n6,602)   "-wime",-hwrsm, (-hwrgn(ir),ir=1,7)
c     write(cmsg,602) "-wime",-hwrsm, (-hwrgn(ir),ir=1,7)
c     call gdput(trim(cmsg))
!
c     write(n6,'(2x)')
c     write(n6,'(6x,"ic",2x,"icx",1x,"icy",1x,"ir",2x,"twrd",8x,
c    >    6("tdnz_",i1,6x))') (iz,iz=0,5)
c     do ic = nskp, ncmax, nskp
c     ir = mrgn(ic)
c     write(n6,'(2x,i6,i5,i4,i3,1p10e12.4)')
c    >  ic, iplx(ic), iply(ic), ir, twrd(ic), (tdnz(iz,ic),iz=0,5)
c     enddo
!
      return
      end
