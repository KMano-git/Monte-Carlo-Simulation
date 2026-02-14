!**********************************************************************
      subroutine ntdbsp(nft)
!**********************************************************************
      use cntcom, only : lbrc, nrcmx, rcfc, scfc
      use cntmnt, only : dotn, ssn, ssp, swe, swi
      use cplcom, only : vte
      use cplmet, only : hvol, icel, itsle, jcel, jtmax, kreg
      use csonic, only : itim, time
      implicit none
!
!::argument
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer nft
      integer, intent(in) :: nft
!
!::local variables
      real*8   zsn(10),zsp(10),zwe(10),zwi(10),zvl(10),zte(10)
      integer  i, it, jt, jmax, icx, icy, ia, n
!
      do i = 1, 10
      zsn(i) = 0.0d0
      zsp(i) = 0.0d0
      zwe(i) = 0.0d0
      zwi(i) = 0.0d0
      zvl(i) = 0.0d0
      zte(i) = 0.0d0
      enddo
!
      zte(1) = 1.0d20
      zte(2) = 10.0d0
      zte(3) = 5.0d0
      zte(4) = 2.0d0
      zte(5) = 1.0d0
!
      do it = 2, itsle
      jmax = jtmax(it)
      do jt = 2, jmax-1
      icx = jcel(jt,it)
      icy = icel(jt,it)
      if( kreg(icx,icy).ne.3 ) cycle
!
!::outer divertor
      ia = 1
      n  = 1
      zsn(n) = zsn(n) + ssn(icx,icy,ia)*hvol(icx,icy)
      zsp(n) = zsp(n) + ssp(icx,icy,ia)*hvol(icx,icy)
      zwi(n) = zwi(n) + swi(icx,icy)*hvol(icx,icy)
      zwe(n) = zwe(n) + swe(icx,icy)*hvol(icx,icy)
      zvl(n) = zvl(n) + hvol(icx,icy)
!::Te < 10.0 eV
      if( vte(icx,icy).le.zte(2) ) then
      n  = 2
      zsn(n) = zsn(n) + ssn(icx,icy,ia)*hvol(icx,icy)
      zsp(n) = zsp(n) + ssp(icx,icy,ia)*hvol(icx,icy)
      zwi(n) = zwi(n) + swi(icx,icy)*hvol(icx,icy)
      zwe(n) = zwe(n) + swe(icx,icy)*hvol(icx,icy)
      zvl(n) = zvl(n) + hvol(icx,icy)
      endif
!::Te < 5.0 eV
      if( vte(icx,icy).le.zte(3) ) then
      n  = 3
      zsn(n) = zsn(n) + ssn(icx,icy,ia)*hvol(icx,icy)
      zsp(n) = zsp(n) + ssp(icx,icy,ia)*hvol(icx,icy)
      zwi(n) = zwi(n) + swi(icx,icy)*hvol(icx,icy)
      zwe(n) = zwe(n) + swe(icx,icy)*hvol(icx,icy)
      zvl(n) = zvl(n) + hvol(icx,icy)
      endif
!::Te < 2.0 eV
      if( vte(icx,icy).le.zte(4) ) then
      n  = 4
      zsn(n) = zsn(n) + ssn(icx,icy,ia)*hvol(icx,icy)
      zsp(n) = zsp(n) + ssp(icx,icy,ia)*hvol(icx,icy)
      zwi(n) = zwi(n) + swi(icx,icy)*hvol(icx,icy)
      zwe(n) = zwe(n) + swe(icx,icy)*hvol(icx,icy)
      zvl(n) = zvl(n) + hvol(icx,icy)
      endif
!::Te < 1.0 eV
      if( vte(icx,icy).le.zte(5) ) then
      n  = 5
      zsn(n) = zsn(n) + ssn(icx,icy,ia)*hvol(icx,icy)
      zsp(n) = zsp(n) + ssp(icx,icy,ia)*hvol(icx,icy)
      zwi(n) = zwi(n) + swi(icx,icy)*hvol(icx,icy)
      zwe(n) = zwe(n) + swe(icx,icy)*hvol(icx,icy)
      zvl(n) = zvl(n) + hvol(icx,icy)
      endif
!
      enddo
      enddo
!
!::debug write
      write(nft,'(/2x,"*** ntdbsp ***   itim =",i5,"  time =",1pe12.3
     >   ,"  dotn =",1pe12.3)') itim,time,dotn
      write(nft,'(2x,a,2x,15(a5,2x))')   "lbrc",(lbrc(i),i=1,nrcmx)
      write(nft,'(2x,a,2x,15(f5.1,2x))') "rcfc",(rcfc(i),i=1,nrcmx)
      write(nft,'(2x,a,2x,15(f5.1,2x))') "scfc",(scfc(i),i=1,nrcmx)
      write(nft,'(2x)')
      write(nft,'(2x,a3,2x,1p5e12.3)') "Te ",(zte(i),i=1,5)
      write(nft,'(2x,a3,2x,1p5e12.3)') "vol",(zvl(i),i=1,5)
      write(nft,'(2x,a3,2x,1p5e12.3)') "sn ",(zsn(i),i=1,5)
      write(nft,'(2x,a3,2x,1p5e12.3)') "sp ",(zsp(i),i=1,5)
      write(nft,'(2x,a3,2x,1p5e12.3)') "wi ",(zwi(i),i=1,5)
      write(nft,'(2x,a3,2x,1p5e12.3)') "we ",(zwe(i),i=1,5)
!
      return
      end
