!**********************************************************************
      subroutine figdat( cpath, ctime )
!**********************************************************************
!
!  N0-data  DTMCP  deni(ic),teme,temi,  tden0(ic),teng0
!                  IMPMC                soldor (call plsorc)
!                                       /cntsrc/
!
!        a005128/PRJ/sonicV3/src/soldor     2013/06/05
!
!  group in charge of the code   2-GRP type
!    [PLS]: 1 / NTL]: 2 /[IMP]: 2 /  <== cdgr(i) i=1,2,3(ncode)
!
!----------------------------------------------------------------------
      use cimctl, only : nprg
      use cimden, only : csput, fsput, ncmx, nsput, nzmx, tdnz, tionz
     >    , trecz, tthz, tvlz, twrd
      use cntcom, only : iplx, iply, mrgn, ncmax
      use csize,  only : ndy
      use cunit,  only : lmspe, lmype, mype, n6
      use dbg_mod, only : kimpsy
      implicit none
! arguments
      character, intent(in) :: cpath*(*), ctime*(*)
! cpath : path name
! ctime : time

! local variables
      integer :: i60
      integer :: nskp = 100
      integer :: nft, ia, ic, iz, ir, i
      character :: cdsn*200, ched*20, ctyp*10
      character(10) :: cact
      real*8  ::  hwrsm
      real*8, dimension(10) :: hwrgn
      real(8) ::  twr
      real(8), dimension(10) ::  twr_rg
      real(8), dimension(ndy) :: twr_it
!
      if( lmype.ne.lmspe ) return
!
!::impmc
      write(n6,'(/2x,"*** figdat ***",2x,a)') "DTIMP  w"
!
      i60 = 180000 + mype
!
      cact = "w"
!
!---------------------------------------------------------------------
!::DTIMP/DTIPF : impurity data
!---------------------------------------------------------------------
!  nsizp, nsizs, nsizc, nsput, nzmx, ncmx, fsput, twrd, tdnz, csput
!
!  IMP_2: nsput, nzmx, ncmx, fsput, twrd, tdnz, csput
!  nsizp, nsizs, nsizc
!
!      ctyp = "DTIMP"
      write(ctyp,'(i5)') nprg
      ctyp = adjustl( ctyp )
      if( ctime == ' ' ) then
        cdsn = './DTIMP' // trim( ctyp )
      else
        cdsn =  trim( cpath ) // 'DTIMP' // trim( ctyp )
     >          // '_' // trim( ctime )
      endif

      ched = "*** imdisk ***"
      call imdisk(cdsn,cact)
!
      write(i60,'(/2x,a,2x,a,2x,a"  twrd,tdnz")')
     >     trim(ched), trim(cdsn),cact
      call totrgn(twrd, twr, twr_rg, twr_it)
      write(i60,'(2x,"nsput =",i2,"  nzmx =",i2,"  ncmx =",i5)')
     >   nsput, nzmx, ncmx
      write(i60,'(2x,5(a," =",1pe12.4,2x))')
     >  (csput(i), fsput(i), i = 1, nsput)
      write(i60,'(2x,"totwr=",1pe14.6,2x,a)')   twr,
     >   " = Intg_twrd(1:6) + wrd_core"
      call wrintg(twrd,hwrsm,hwrgn)
      write(i60,'(2x,"Intg-twrd =",1pe14.6,2x,1p10e12.4)')
     >   hwrsm,(hwrgn(ir),ir=1,7)
      write(i60,'(6x,"ic",2x,"icx",1x,"icy",1x,"ir",2x,"twrd",8x,
     >    6("tdnz_",i1,6x))') (iz,iz=0,5)
      do ic = nskp, ncmax, nskp
        ir = mrgn(ic)
        write(i60,'(2x,i6,i5,i4,i3,1p10e12.4)')
     >    ic, iplx(ic), iply(ic), ir, twrd(ic), (tdnz(iz,ic),iz=0,5)
      enddo

! Additional impurity informations: YAMOTO
      cact = "y"
      write(ctyp,'(i5)') nprg
      ctyp = adjustl( ctyp )
      if( ctime == ' ' ) then
        cdsn = './IMPSY' // trim( ctyp )
      else
        cdsn =  trim( cpath ) // 'IMPSY' // trim( ctyp )
     >          // '_' // trim( ctime )
      endif

      if( ctime == ' ' .or. kimpsy == 1 ) then
        call imdisk(cdsn,cact)

        write(i60,'(6x,"ic",2x,"icx",1x,"icy",1x,"ir",8x,
     >      6("tthz_",i1,6x))') (iz,iz=1,6)
        do ic = nskp, ncmax, nskp
          ir = mrgn(ic)
          write(i60,'(2x,i6,i5,i4,i3,1p10e12.4)')
     >      ic, iplx(ic), iply(ic), ir, (tthz(iz,ic),iz=1,6)
        enddo

        write(i60,'(6x,"ic",2x,"icx",1x,"icy",1x,"ir",8x,
     >      6("tvlz_",i1,6x))') (iz,iz=1,6)
        do ic = nskp, ncmax, nskp
          ir = mrgn(ic)
          write(i60,'(2x,i6,i5,i4,i3,1p10e12.4)')
     >     ic, iplx(ic), iply(ic), ir, (tvlz(iz,ic),iz=1,6)
        enddo

        write(i60,'(6x,"ic",2x,"icx",1x,"icy",1x,"ir",8x,
     >      6("tionz_",i1,6x))') (iz,iz=1,6)
        do ic = nskp, ncmax, nskp
          ir = mrgn(ic)
          write(i60,'(2x,i6,i5,i4,i3,1p10e12.4)')
     >      ic, iplx(ic), iply(ic), ir, (tionZ(iz,ic),iz=1,6)
        enddo

        write(i60,'(6x,"ic",2x,"icx",1x,"icy",1x,"ir",8x,
     >      6("trecz_",i1,6x))') (iz,iz=1,6)
        do ic = nskp, ncmax, nskp
          ir = mrgn(ic)
          write(i60,'(2x,i6,i5,i4,i3,1p10e12.4)')
     >      ic, iplx(ic), iply(ic), ir, (trecZ(iz,ic),iz=1,6)
        enddo
      endif
      return
      end
