!***********************************************************************
      subroutine impprof
!***********************************************************************
!
!    profile along poloidal length
!    avaraged over 2~3 flux tube   
!               under development for plasma parameter
!               omit this function
!
!----------------------------------------------------------------------
      use csize
      use cimcom
      use cimden
      use cntcom
      use cplcom
      use cplimp
      use cplmet
      use cplwrd
      use cplwrd2
      implicit none

!
!::local variables
      integer :: n6 = 6
      integer :: nft, n, n1, n2
      character(80) :: drimp, dsn
!
      integer :: i, j, it, nty
      logical :: lex
!
      integer :: np
      real(8), allocatable :: pls(:), pne(:), pni(:), pte(:), pti(:)
     > ,  gne(:), gni(:), gte(:), gti(:)

      integer    istat
      character  cmsg*80

! allocate pls (ndx is not parameter) 23/1/25
      if( .not. allocated( pls ) ) then
        allocate(pls(ndx), pne(ndx), pni(ndx), pte(ndx), pti(ndx)
     > ,  gne(ndx), gni(ndx), gte(ndx), gti(ndx)
     > ,  stat = istat)
        if( istat /= 0 ) then
          write(cmsg,'(a,i6)')
     >     'pls allocate error in impprof,
     >      istat = ', istat
          call wexit( 'impprof', trim( cmsg ) )
        endif
        pls(1:ndx) = 0.0_8
        pne(1:ndx) = 0.0_8
        pni(1:ndx) = 0.0_8
        pte(1:ndx) = 0.0_8
        pti(1:ndx) = 0.0_8
        gne(1:ndx) = 0.0_8
        gni(1:ndx) = 0.0_8
        gte(1:ndx) = 0.0_8
        gti(1:ndx) = 0.0_8
      endif

      call getenv("IMP1D",drimp)

!::tube number
      n1  = itsls + 1
      n2  = itmpe
!
!::output file  (all tube)
      nft = 810
      do n = n1, n2
        it = n
        if( it.eq.itpve ) cycle
        do nty = 1, wmc_nty
          call imprf1d(drimp,nty,it,ndx)
        enddo
        call prfvgrd(it,np,pls,pne,pni,pte,pti,gne,gni,gte,gti,ndx)
        call NTLtxt(drimp,it)
      enddo
!
      return
! force for multi-imp case is under development
!
!::thermal force & frcition force
      call set_forc    !  set gdte, gdti & clear cfgte, cfgti
      call set_gcoef   !  set cfgte, cfgti
!
      do n = n1, n2
        it = n
        if( it.eq.itpve ) cycle
        do nty = 1, wmc_nty
          call lst_gdfrc2(nty,it)
        enddo
      enddo
!
      return

      contains
      subroutine NTLtxt(drimp,it)
      use cntsrc, only : tden0, teng0, tdeng, tengg
      use cimset, only : wfinp
      use cntpls, only : deni, temi
      implicit none

!::argument
      character(*), intent(in) :: drimp
      integer, intent(in) :: it
!
!::local variables
      integer nft, i, j, jw, ic
      integer jwst, jwen
      real(8) zlmx
      character(80) :: dsn
      character(3) :: cno
      real(8), allocatable :: alnx(:)
      integer    istat
      character  cmsg*80
      real(8) p0,pi,pg ! neutral particle pressure
! function
      integer imox, imoy

! allocate alnx (ndx is not parameter) 23/1/25
      if( .not. allocated( alnx ) ) then
        allocate(alnx(ndx), stat = istat)
        if( istat /= 0 ) then
          write(cmsg,'(a,i6)')
     >     'alnx allocate error in imprf1d,
     >      istat = ', istat
          call wexit( 'imprf1d', trim( cmsg ) )
        endif
        alnx(1:ndx) = 0.0_8
      endif

      write(cno,'(i3)') it
      nft = 103
      dsn = trim(drimp) // "/it" // trim(cno) // "NTL.txt"
      call delspc(dsn)

      open(unit=nft,file=dsn)
      write(nft,'(2x,a,2x,i4)') trim(dsn), it
      write(nft,'(3x,"it",3x,"ix",3x,"iy",2x,"lp",9x,"lp2",8x,
     >  "N0",9x,"E0",9x,"P0",9x,
     >  "Ng",9x,"Eg",9x,"Pg",9x,
     >  "Ni",9x,"Ei",9x,"Pi")')

      jwst = jtmin(it)
      jwen = jtmax(it)
      if( it.gt.itmps ) then
        jwst = jwst + 1
        jwen = jwen - 1
      endif
      call plenc(it,alnx)
      zlmx = alnx(jwen)
!
      do jw = jwst, jwen
        j = jcel(jw,it)
        i = icel(jw,it)
        ic = mcel(j,i)
        if( ic.le.0 ) cycle
!
        p0 = tden0(ic,1)*teng0(ic,1)*1.60210d-19
        pg = tdeng(ic,1)*tengg(ic,1)*1.60210d-19
        pi = deni(ic,1)*temi(ic)*1.60210d-19
        write(nft,'(2x,i3,2i5,1p2e11.3,1p9e11.3)') 
     >   it, imox(ic), imoy(ic), alnx(jw), zlmx-alnx(jw),
     >   tden0(ic,1), teng0(ic,1), p0,
     >   tdeng(ic,1), tengg(ic,1), pg,
     >   deni(ic,1),  temi(ic),    pi
      enddo  ! loop(jw)
      close(nft)

      return
      end subroutine NTLtxt
!
      end subroutine impprof
!
!***********************************************************************
      subroutine imprf1d(drimp,nty,nt,ndm)
!***********************************************************************
      use csize
      use cimcom
      use cimden
      use cntcom
      use cntpls, only: dene, deni, teme, temi
      use cplcom
      use cplimp
      use cplmet
      use cplwrd
      use cplwrd2
      use cunit
      use cimset
      use cplimp_plimp
      implicit none

!::argument
      character(*) :: drimp
      integer :: nt, nty
!
!::local variables
      character(3) :: cno
      character(1) :: cnty
      integer :: ndm, ic, iz, it, izmx, jwst, jwen, j, i, ia
      real(8) :: pls(ndm)
      real(8) :: dls(ndm), plsh2(ndm)
      real(8) :: plsh(ndm), pneh(ndm), pnih(ndm), pteh(ndm), ptih(ndm)
      integer :: its, jw, k, nnty
      integer :: imox, imoy
      integer :: nft1, nft2, nft_y1, nft_y2, nft_y3, nft_y4
      integer :: nft_y5, nft_y6, nft_y7, nft_y8, nft_y9
      real(8) :: z0, vwrt, vwrc, vwrm, wrdm, vwrm0
      real(8) :: zlmx, zef, tnz, eff_counter
      real(8),dimension(1:nzsmx) :: zav, znzi, znzt
      real(8),dimension(1:nzsmx) :: zeff_ech
      real(8) :: zeff_local
! modified 4/4 lines treat 4 or more impurities with IMPMC
      ! real(8) :: hvl, hne, hni(ndsp), hnz(0:ndis2L)
      ! real(8) :: hfrz(0:ndis2L), hthz(0:ndis2L), hvlz(0:ndis2L)
      ! real(8) :: hionz(0:ndis2L), hrecz(0:ndis2L)
      ! real(8) :: hradiz(0:ndis2L), hradliz(0:ndis2L), hradrz(0:ndis2L)
      real(8) :: hvl, hne, hni(ndsp), hnz(0:ndmis)
      real(8) :: hfrz(0:ndmis), hthz(0:ndmis), hvlz(0:ndmis)
      real(8) :: hionz(0:ndmis), hrecz(0:ndmis)
      real(8) :: hradiz(0:ndmis), hradliz(0:ndmis), hradrz(0:ndmis)
      real(8) :: nfrz, nthz, nvlz
      real(8),dimension(1:nzsmx) :: tnfrz, tnthz, tnvlz
      real(8) :: hflw(ndsp), flw
      character(80) :: dsn1, dsn2, dsn_y1, dsn_y2
      character(80) :: dsn_y3, dsn_y4, dsn_y5, dsn_y6
      character(80) :: dsn_y7, dsn_y8, dsn_y9
      real(8) :: zfac = 0.0_8
      real(8), allocatable :: alnx(:)

      integer    istat
      character  cmsg*80
      character(len=100) :: nft2_label

! allocate alnx (ndx is not parameter) 23/1/25
      if( .not. allocated( alnx ) ) then
        allocate(alnx(ndx), stat = istat)
        if( istat /= 0 ) then
          write(cmsg,'(a,i6)')
     >     'alnx allocate error in imprf1d,
     >      istat = ', istat
          call wexit( 'imprf1d', trim( cmsg ) )
        endif
        alnx(1:ndx) = 0.0_8
      endif

      it = nt
      izmx = min0(ismaxL(nty),74)  !  W 74
      if(nty.eq.1) then
        zfac = wfac
      else
        zfac = 1.0d0
      endif
!
      write(n6,'(2x,"*** imprf1d ***   it =",i3)') it
      write(cno,'(i3)') it
      write(cnty,'(i1)')nty
!
      nft1 = 101
      nft2 = 102
      ! Additional parameters : Yamoto
      nft_y1 = 1101
      nft_y2 = 1102
      nft_y3 = 1103
      nft_y4 = 1104
      nft_y5 = 1105
      nft_y6 = 1106
      nft_y7 = 1107
      nft_y8 = 1108
      nft_y9 = 1109

      dsn1 = trim(drimp) // "/it" // trim(cno) // "WRD"//cnty//".txt"
      dsn2 = trim(drimp) // "/it" // trim(cno) // "IMP"//cnty//".txt"
      dsn_y1 = trim(drimp) // "/it" // trim(cno) // "YDT"//cnty//".txt"
      dsn_y2 = trim(drimp) // "/it" // trim(cno) // "FRI"//cnty//".txt"
      dsn_y3 = trim(drimp) // "/it" // trim(cno) // "THF"//cnty//".txt"
      dsn_y4 = trim(drimp) // "/it" // trim(cno) // "VLZ"//cnty//".txt"
      dsn_y5 = trim(drimp) // "/it" // trim(cno) // "ION"//cnty//".txt"
      dsn_y6 = trim(drimp) // "/it" // trim(cno) // "REC"//cnty//".txt"
      dsn_y7 = trim(drimp) // "/it" // trim(cno) // "RDI"//cnty//".txt"
      dsn_y8 = trim(drimp) // "/it" // trim(cno) // "RDL"//cnty//".txt"
      dsn_y9 = trim(drimp) // "/it" // trim(cno) // "RDR"//cnty//".txt"
      call delspc(dsn1)
      call delspc(dsn2)
      call delspc(dsn_y1)
      call delspc(dsn_y2)
      call delspc(dsn_y3)
      call delspc(dsn_y4)
      call delspc(dsn_y5)
      call delspc(dsn_y6)
      call delspc(dsn_y7)
      call delspc(dsn_y8)
      call delspc(dsn_y9)

!::position and distance between cell boundaries
      call plenc(it,alnx)  ! Lp(j)
      call slenb(it,plsh2) ! Lsh(j+1/2)
      call gdlen(it,dls)   ! dLs(j)
!
      jwst = jtmin(it)
      jwen = jtmax(it)

!xx      if( it.gt.itmps ) then
!xx        jwst = jwst + 1
!xx        jwen = jwen - 1
!xx      endif
!
!::Lsh(j+1/2),Ls(j)
      plsh(jwst) = 0.0d0
      do jw = jwst+1, jwen-1
        plsh(jw) = plsh(jw-1) + dls(jw)
      enddo
      plsh(jwen) = 0.0d0

      pls(jwst) = plsh(jwst) - 0.5d0*dls(jwst)
      do jw = jwst+1, jwen
        pls(jw) = plsh(jw-1) + 0.5d0*dls(jw)
      enddo

!
      open(unit=nft1,file=dsn1)
      open(unit=nft2,file=dsn2)
      open(unit=nft_y1,file=dsn_y1)
      open(unit=nft_y2,file=dsn_y2)
      open(unit=nft_y3,file=dsn_y3)
      open(unit=nft_y4,file=dsn_y4)
      open(unit=nft_y5,file=dsn_y5)
      open(unit=nft_y6,file=dsn_y6)
      open(unit=nft_y7,file=dsn_y7)
      open(unit=nft_y8,file=dsn_y8)
      open(unit=nft_y9,file=dsn_y9)

      write(nft1,'(2x,a,2x,i4,5x,"Z = ",i2)') trim(dsn1), it,ismaxL(nty)
      write(nft2,'(2x,a,2x,i4,5x,"Z = ",i2)') trim(dsn2), it,ismaxL(nty)
      write(nft_y1,'(2x,a,2x,i4,5x,"Z = ",i2)')
     >  trim(dsn2), it,ismaxL(nty)
      write(nft_y2,'(2x,a,2x,i4,5x,"Z = ",i2)')
     >  trim(dsn2), it,ismaxL(nty)
      write(nft_y3,'(2x,a,2x,i4,5x,"Z = ",i2)')
     >  trim(dsn2), it,ismaxL(nty)
      write(nft_y4,'(2x,a,2x,i4,5x,"Z = ",i2)')
     >  trim(dsn2), it,ismaxL(nty)
      write(nft_y5,'(2x,a,2x,i4,5x,"Z = ",i2)')
     >  trim(dsn2), it,ismaxL(nty)
      write(nft_y6,'(2x,a,2x,i4,5x,"Z = ",i2)')
     >  trim(dsn2), it,ismaxL(nty)
!
      write(nft1,'(2x,"wfac =",f9.5,"  wfinp =",f9.5,2x,a)') zfac,wfinp,
     >  "vwrt = vwrc(cr) + vwrm(IMPMC)  vwrm = vwrm1 (twrd)+vwrm2"
      write(nft2,'(2x,"wfac =",f9.5,"  wfinp =",f9.5,2x,a)') zfac,wfinp,
     >  "vwrt = vwrc(cr) + vwrm(IMPMC)  vwrm = vwrm1 (twrd)+vwrm2"
      write(nft_y1,'(2x,"wfac =",f9.5,"  wfinp =",f9.5,2x,a)') 
     >  zfac,wfinp,
     >  "vwrt = vwrc(cr) + vwrm(IMPMC)  vwrm = vwrm1 (twrd)+vwrm2"
      write(nft_y2,'(2x,"wfac =",f9.5,"  wfinp =",f9.5,2x,a)') 
     >  zfac,wfinp,
     >  "vwrt = vwrc(cr) + vwrm(IMPMC)  vwrm = vwrm1 (twrd)+vwrm2"
      write(nft_y3,'(2x,"wfac =",f9.5,"  wfinp =",f9.5,2x,a)') 
     >  zfac,wfinp,
     >  "vwrt = vwrc(cr) + vwrm(IMPMC)  vwrm = vwrm1 (twrd)+vwrm2"
      write(nft_y4,'(2x,"wfac =",f9.5,"  wfinp =",f9.5,2x,a)') 
     >  zfac,wfinp,
     >  "vwrt = vwrc(cr) + vwrm(IMPMC)  vwrm = vwrm1 (twrd)+vwrm2"
      write(nft_y5,'(2x,"wfac =",f9.5,"  wfinp =",f9.5,2x,a)') 
     >  zfac,wfinp,
     >  "vwrt = vwrc(cr) + vwrm(IMPMC)  vwrm = vwrm1 (twrd)+vwrm2"
      write(nft_y6,'(2x,"wfac =",f9.5,"  wfinp =",f9.5,2x,a)') 
     >  zfac,wfinp,
     >  "vwrt = vwrc(cr) + vwrm(IMPMC)  vwrm = vwrm1 (twrd)+vwrm2"
      write(nft_y7,'(2x,"wfac =",f9.5,"  wfinp =",f9.5,2x,a)') 
     >  zfac,wfinp,
     >  "vwrt = vwrc(cr) + vwrm(IMPMC)  vwrm = vwrm1 (twrd)+vwrm2"
      write(nft_y8,'(2x,"wfac =",f9.5,"  wfinp =",f9.5,2x,a)') 
     >  zfac,wfinp,
     >  "vwrt = vwrc(cr) + vwrm(IMPMC)  vwrm = vwrm1 (twrd)+vwrm2"
      write(nft_y9,'(2x,"wfac =",f9.5,"  wfinp =",f9.5,2x,a)') 
     >  zfac,wfinp,
     >  "vwrt = vwrc(cr) + vwrm(IMPMC)  vwrm = vwrm1 (twrd)+vwrm2"
!
      write(nft1,'(3x,"it",3x,"ix",3x,"iy",2x,"lp",9x,"lp2",8x,
     >  "dene",7x,"vnezef",5x,"hne",8x,"deni",7x,"teme",7x,
     >  "temi",7x,"vwrt",7x,"vwrc",7x,"vwrm",7x,"vwrm",i1,6x
     >  "Nz/Ni",6x,"Nz/Ne",6x,"<Z>",8x,"Zeff")') nty
!!!
      write(nft2_label,"(i0)") izmx+1
      nft2_label = 
     > "(3x,a2,3x,a2,3x,a2,2x,a2,9x,a3,8x,a3,8x,a3,8x," 
     > // trim(nft2_label)
     > // "(i2.2,9x),a3,8x,a3)" 
      write(nft2,nft2_label)
     >   "it","ix","iy","lp","lp2","Nzi","Nzt",
     >   (iz,iz=0,izmx),
     >   "E00","P00"
!!!
      write(nft_y1,'(3x,"it",2x,"Lp",9x,"Ls",9x,"Vf",9x,"nD+"
     >     8x,"Ffrt",7x,"Ftht",7x,"Vzt")')
      write(nft_y2,'(3x,"it",2x,"Lp",9x,"Ls",9x,41(i2.2,9x))')
     >    (iz,iz=0,izmx)
      write(nft_y3,'(3x,"it",2x,"Lp",9x,"Ls",9x,41(i2.2,9x))')
     >    (iz,iz=0,izmx)
      write(nft_y4,'(3x,"it",2x,"Lp",9x,"Ls",9x,41(i2.2,9x))')
     >    (iz,iz=0,izmx)
      write(nft_y5,'(3x,"it",2x,"Lp",9x,"Ls",9x,41(i2.2,9x))')
     >    (iz,iz=0,izmx)
      write(nft_y6,'(3x,"it",2x,"Lp",9x,"Ls",9x,41(i2.2,9x))')
     >    (iz,iz=0,izmx)
      write(nft_y7,'(3x,"it",2x,"Lp",9x,"Ls",9x,41(i2.2,9x))')
     >    (iz,iz=0,izmx)
      write(nft_y8,'(3x,"it",2x,"Lp",9x,"Ls",9x,41(i2.2,9x))')
     >    (iz,iz=0,izmx)
      write(nft_y9,'(3x,"it",2x,"Lp",9x,"Ls",9x,41(i2.2,9x))')
     >    (iz,iz=0,izmx)
!
      jwst = jtmin(it)
      jwen = jtmax(it)
      if( it.gt.itmps ) then
        jwst = jwst + 1
        jwen = jwen - 1
      endif
      call plenc(it,alnx)
      zlmx = alnx(jwen)
!
      do jw = jwst, jwen
      j = jcel(jw,it)
      i = icel(jw,it)
      ic = mcel(j,i)
      if( ic.le.0 ) cycle
!
      hni(1:nion) = deni(ic,1:nion)
      hflw(1:nion) = vva(j,i,1:nion)

      do nnty=1,wmc_nty
      hnz(0:ismaxL(nnty)) = tdnzL(0:ismaxL(nnty),ic,nnty) ! *zfac 
                             ! wfac(wfinp) already multiplyed in set_imp
      hfrz(0:ismaxL(nnty)) = tfrzL(0:ismaxL(nnty),ic,nnty) ! *zfac 
                             ! wfac(wfinp) already multiplyed in set_imp
      hthz(0:ismaxL(nnty)) = tthzL(0:ismaxL(nnty),ic,nnty) ! *zfac 
                             ! wfac(wfinp) already multiplyed in set_im
      hvlz(0:ismaxL(nnty)) = tvlzL(0:ismaxL(nnty),ic,nnty) ! *zfac 
                             ! wfac(wfinp) already multiplyed in set_imp
      hionz(0:ismaxL(nnty)) = tionZL(0:ismaxL(nnty),ic,nnty) ! *zfac 
                             ! wfac(wfinp) already multiplyed in set_imp
      hrecz(0:ismaxL(nnty)) = trecZL(0:ismaxL(nnty),ic,nnty) ! *zfac 
                             ! wfac(wfinp) already multiplyed in set_imp
      hradiz(0:ismaxL(nnty)) = sradi(0:ismaxL(nnty),ic,nnty) ! *zfac 
                             ! wfac(wfinp) already multiplyed in set_imp
      hradliz(0:ismaxL(nnty)) = sradli(0:ismaxL(nnty),ic,nnty) ! *zfac 
                             ! wfac(wfinp) already multiplyed in set_imp
      hradrz(0:ismaxL(nnty)) = sradr(0:ismaxL(nnty),ic,nnty) ! *zfac 
                             ! wfac(wfinp) already multiplyed in set_imp

!
!::Ne, Nz, <Z>, Zef
      hne = 0.0d0
      zef = 0.0d0
      flw = 0.0d0
      do ia = 1, nion
      hne = hne + aza(ia)*hni(ia)
      zef = zef + aza(ia)**2*hni(ia)
      flw = flw + hflw(ia)
      enddo
! Yamoto modified//
      tnz = 0.0d0
      nfrz = 0.d0
      nthz = 0.d0
      nvlz = 0.d0
      zav = 0.0d0
      zeff_local = 0.0d0
      eff_counter = 0.0d0
      do iz = 1, ismaxL(nnty)
      tnz = tnz + hnz(iz)
      hne = hne + dfloat(iz)*hnz(iz)
      nfrz = nfrz + hfrz(iz)/(dfloat(iz)**2)
      nthz = nthz + hthz(iz)/(dfloat(iz)**2)
      nvlz = nvlz + hvlz(iz)
      zav = zav + dfloat(iz)*hnz(iz)
      zef = zef + dfloat(iz)**2*hnz(iz)
      if(hvlz(iz).ne.0.0d0) eff_counter=eff_counter+1.0d0
      enddo
      zeff_ech(nnty) = zef
      znzi(nnty) = tnz          ! ion
      znzt(nnty) = tnz + hnz(0) ! total
      if(eff_counter.ne.0.0d0)then
      tnfrz(nnty) = nfrz/eff_counter
      tnthz(nnty) = nthz/eff_counter
      tnvlz(nnty) = nvlz/eff_counter
      else
      tnfrz(nnty) = 0.0d0
      tnthz(nnty) = 0.0d0
      tnvlz(nnty) = 0.0d0
      endif

      if( tnz.ne.0.0d0 ) zav(nnty) = zav(nnty)/tnz

! Filtering SYamoto
! For noisy Ffrz, and Fthz with low nZ
      do iz= 1, ismaxL(nnty)
         if(hnz(iz).lt.1.0d16)then
            hfrz(iz) = 0.0d0
            hthz(iz) = 0.0d0
         end if
      enddo ! iz

      enddo ! nnty

      if( hne.ne.0.0d0 ) zef = zef/hne
!
!::radiation
!  v : soldor variables  t:total m:monte c:corona   vwrm = wrdm
      vwrt = -wime(ic)  ! total radiation in soldor
      vwrc = -wcre(ic)  ! Corona Model in soldor
      vwrm = -wmce(ic)  ! Total MC radiation in soldor
      vwrm0= -wmc_zwe(ic,nty)   ! Each MC radiation in soldor
!c$$$      vwrm0= sum(sradi(:,ic,nty))
!c$$$     > +sum(sradli(:,ic,nty))
!c$$$     > +sum(sradr(:,ic,nty)) - wmc_zwe(ic,nty)     ! =0 check
!      wrdm = twrd(ic)   ! MC radiation in IMPMC 
! 
      ia = 1
      z0 = 1.0d-30

      write(nft1,'(2x,i3,2i5,1p15e11.3,1p50e11.3)')
     > it, imox(ic), imoy(ic), alnx(jw), zlmx-alnx(jw),
     > dmax1(z0,dene(ic)), dmax1(z0,vnezef(j,i)), dmax1(z0,hne), 
     > dmax1(z0,deni(ic,ia)), dmax1(z0,teme(ic)),dmax1(z0,temi(ic)),
     > dmax1(z0,vwrt),dmax1(z0,vwrc),dmax1(z0,vwrm),
     > dmax1(z0,vwrm0), dmax1(z0,znzt(nty)/deni(ic,ia)), 
     > dmax1(z0,znzt(nty)/dene(ic)), dmax1(z0,zav(nty)), dmax1(z0,zef)
!
      write(nft2,'(2x,i3,2i5,1p2e11.3,1p50e11.3,1p2e11.3)') 
     > it, imox(ic), imoy(ic), alnx(jw), zlmx-alnx(jw),
     > znzi(nty), znzt(nty),
     >  (dmax1(z0,tdnzL(iz,ic,nty)),iz=0,izmx),
     > tengzL(0,ic,nty),tprzL(0,ic,nty)
!
      write(nft_y1,'(2x,i3,1p50e11.3)') it, alnx(jw), pls(jw), flw,
     >  hni(1), tnfrz(nty), tnthz(nty), tnvlz(nty)
      write(nft_y2,'(2x,i3,1p50e11.3)') it, alnx(jw), pls(jw),
     >  ((tfrzL(iz,ic,nty)+z0),iz=0,izmx)
      write(nft_y3,'(2x,i3,1p50e11.3)') it, alnx(jw), pls(jw),
     >  ((tthzL(iz,ic,nty)+z0),iz=0,izmx)
      write(nft_y4,'(2x,i3,1p50e11.3)') it, alnx(jw), pls(jw),
     >  ((tvlzL(iz,ic,nty)+z0),iz=0,izmx)
      write(nft_y5,'(2x,i3,1p50e11.3)') it, alnx(jw), pls(jw),
     >  ((tionZL(iz,ic,nty)+z0),iz=0,izmx)
      write(nft_y6,'(2x,i3,1p50e11.3)') it, alnx(jw), pls(jw),
     >  ((trecZL(iz,ic,nty)+z0),iz=0,izmx)
      write(nft_y7,'(2x,i3,1p50e11.3)') it, alnx(jw), pls(jw),
     >  ((sradi(iz,ic,nty)+z0),iz=0,izmx)
      write(nft_y8,'(2x,i3,1p50e11.3)') it, alnx(jw), pls(jw),
     >  ((sradli(iz,ic,nty)+z0),iz=0,izmx)
      write(nft_y9,'(2x,i3,1p50e11.3)') it, alnx(jw), pls(jw),
     >  ((sradr(iz,ic,nty)+z0),iz=0,izmx)
      enddo  ! loop(jw)
!
      close(nft1)
      close(nft2)
      close(nft_y1)
      close(nft_y2)
      close(nft_y3)
      close(nft_y4)
      close(nft_y5)
      close(nft_y6)
      close(nft_y7)
      close(nft_y8)
      close(nft_y9)

      return
!
!::error
 910  continue
      write(n6,'(/2x,"index error at sub. impprof ")')
      call wexit("impprof","index error")
      end
