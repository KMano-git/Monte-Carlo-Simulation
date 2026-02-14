!***********************************************************************
      subroutine plhist(kk)
!***********************************************************************
      use cplcom, only : cfth, nfth
      use cplhst, only : cvnm, hstt, hstv, iphs, jphs, ndhp, ndht, ndhv
     >    , nhpn, nhtm, nhvl
      use cplmet, only : icel, itmpe, itmps, itpve, itsle, jcdp1, jcdp2
     >    , jcxp1, jcxp2
      use cpmpls, only : imd1, imd2, jmd1, jmd2
      use csonic, only : itend, itim, nhsav
      use cunit,  only : lmspe, lmype, n6
      implicit none
!
!::argument
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  kk
      integer, intent(in) :: kk
!
!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  ii, it, ic, jmax, n, jc, nft, irec, ll, itmh, lenx, i, j
      integer  ii, it, ic, n, jc, nft, irec, i, j
      integer  md
!
      character cv*20, cft*80, clin*120, cpos*10
      character cfrm*80, cdsn*80
      logical lex
!
      save  nft, cft, irec
!
      if( nfth.le.0 ) return
!
!::KSFUJI
      if( lmype /= lmspe ) return
!
!----------------------------------------------------------------------
!::initial set
!----------------------------------------------------------------------
      if( kk.eq.1 ) then
!
!::clear   KSFUJI
      cvnm(1:ndhp) = "   "
      iphs(1:ndhp) = 0
      jphs(1:ndhp) = 0
      nhpn = 0
      nhtm = 0
      nhvl = 0
      hstt(1:ndht) = 0.0_4
      hstv(1:ndht,1:ndhv) = 0.0_4
!
      ii = 0
!
!::(1) time
      ii = ii + 1
      cvnm(ii) = "time"
      jphs(ii) = 99
      iphs(ii) = 99
!
!::(2) max value of residuals
      ii = ii + 1
      cvnm(ii) = "max-dt"
      jphs(ii) = +99
      iphs(ii) = +99
!
!::separatrix
      it = itsle
      ic = icel(1,it)
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   jmax = jtmax(it)
!::(3)
      ii = ii + 1
      cvnm(ii) = "odiv spx"
      jphs(ii) = jcdp1
      iphs(ii) = ic
!::(4)
      ii = ii + 1
      cvnm(ii) = "idiv spx"
      jphs(ii) = jcdp2
      iphs(ii) = ic
!::(5)
      ii = ii + 1
      cvnm(ii) = "odiv-1/2 spx"
      jphs(ii) = (jcdp1+jcxp1)/2
      iphs(ii) = ic
!::(6)
      ii = ii + 1
      cvnm(ii) = "idiv-1/2 spx"
      jphs(ii) = (jcdp2+jcxp2)/2
      iphs(ii) = ic
!::(7)
      ii = ii + 1
      cvnm(ii) = "sol-1/2 spx"
      jphs(ii) = (jcxp1+jcxp2)/2
      iphs(ii) = ic
!
!::sol-wall
      it = 1
      ic = icel(1,it)
!::(8)
      ii = ii + 1
      cvnm(ii) = "odiv s-wall"
      jphs(ii) = jcdp1
      iphs(ii) = ic
!::(9)
      ii = ii + 1
      cvnm(ii) = "idiv s-wall"
      jphs(ii) = jcdp2
      iphs(ii) = ic
!::(10)
      ii = ii + 1
      cvnm(ii) = "sol-1/2 s-wall"
      jphs(ii) = (jcxp1+jcxp2)/2
      iphs(ii) = ic
!
!::private wall
      it = itpve
      ic = icel(1,it)
!::(11)
      ii = ii + 1
      cvnm(ii) = "odiv p-wall"
      jphs(ii) = jcdp1
      iphs(ii) = ic
!::(12)
      ii = ii + 1
      cvnm(ii) = "idiv p-wall"
      jphs(ii) = jcdp2
      iphs(ii) = ic
!
!::main
!::(13)
      it = itmps
      ic = icel(1,it)
      ii = ii + 1
      cvnm(ii) = "main itmps"
      jphs(ii) = (jcxp1+jcxp2)/2
      iphs(ii) = ic
!::(14)
      ii = ii + 1
      it = itmpe
      ic = icel(1,it)
      cvnm(ii) = "main itmpe"
      jphs(ii) = (jcxp1+jcxp2)/2
      iphs(ii) = ic
!
!::(15) Flux
      ii = ii + 1
      cvnm(ii) = "Flux & Ndot"
      jphs(ii) = 99
      iphs(ii) = 99
!
!::(16) Wrad
      ii = ii + 1
      cvnm(ii) = "Wrad"
      jphs(ii) = 99
      iphs(ii) = 99
!
!::(17) O-Sol
      ii = ii + 1
      cvnm(ii) = "osol spx"
      jphs(ii) = jmd1
      iphs(ii) = imd1
!
!::(18) I-Sol
      ii = ii + 1
      cvnm(ii) = "isol spx"
      jphs(ii) = jmd2
      iphs(ii) = imd2
!
!::last
      nhpn = ii
      nhtm = 0
!
      write(n6,'(/2x,"*** plhist ***  nhpn =",i3,"  ndhp =",i3)')
     >    nhpn,ndhp
      if( nhpn.gt.ndhp ) then
      write(clin,'("dimension error  nhpn.gt.ndhp  ",2i5)') nhpn,ndhp
      call wexit( "plhist", clin )
      endif
!
!::debug write
      write(n6,'(5x,"n",4x,"jc",3x,"ic",5x,"cv")')
      do 110 n = 1, nhpn
      jc = jphs(n)
      ic = iphs(n)
      cv = cvnm(n)
      write(n6,'(2x,3i5,2x,a)') n, jc,ic,cv
 110  continue
!
!::initialization of file
      nft = nfth
      cft = cfth
      irec = 0
      nhtm = 0
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   itmh = -1
      return
!
!----------------------------------------------------------------------
!::write file during loop        kk.eq.2
!----------------------------------------------------------------------
      elseif( kk.eq.2 ) then
      nhtm = nhtm + 1
      call plhsav
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   itmh = itim
!
      if( nhtm.eq.nhsav ) then
      irec = irec + 1
      cpos = "append"
      if( irec.eq.1 ) cpos = "rewind"
!xx   open(unit=nft,file=cft,form='unformatted',position=cpos)
      cfrm = "write,binary,"//trim(cpos)
      call nopen(nft,cft,cfrm,cdsn,lex)
!-----
      write(nft)
     >   cvnm,iphs,jphs,nhpn,hstt,hstv,nhtm,nhvl
      close(nft)
!
!::KSFUJI
      nhtm = 0
      hstt(1:ndht) = 0.0_4
      hstv(1:ndht,1:ndhv) = 0.0_4
      endif
!
!----------------------------------------------------------------------
!::write file after calculation  kk.eq.3
!----------------------------------------------------------------------
      elseif( kk.eq.3 ) then
      if( nhtm.gt.0 ) then
      irec = irec + 1
      cpos = "append"
      if( irec.eq.1 ) cpos = "rewind"
!xx   open(unit=nft,file=cft,form='unformatted',position=cpos)
      cfrm = "write,binary,"//trim(cpos)
      call nopen(nft,cft,cfrm,cdsn,lex)
!-----
      write(nft)
     >   cvnm,iphs,jphs,nhpn,hstt,hstv,nhtm,nhvl
      close(nft)
      write(n6,'(2x,"#plhist  write irec =",i6,"  itim =",i6,
     >  "  nhtm =",i5)') irec, itim, nhtm
      endif
!
      write(n6,'(2x,"*** plhist ***   kk= 3   irec =",i6,"  nhtm =",i4,
     >   "  file =",a)') irec,nhtm,cft(1:20)
!
!xx   open(unit=nft,file=cft,form='unformatted',position='rewind')
      cfrm = "read,binary,rewind"
      call nopen(nft,cft,cfrm,cdsn,lex)
!-----
      md = itend/5
      if( md.eq.0 ) md = 1
      j = 0
 310  continue
      read(nft,end=320,err=330)
     >   cvnm,iphs,jphs,nhpn,hstt,hstv,nhtm,nhvl
      do i = 1, nhtm
      j = j + 1
      if( mod(j,md).ne.0 ) cycle
      write(n6,'(2x,"plhist",2x,i6,1p4e12.3)')
     >   j, hstv(i,1), hstv(i,9), hstv(i,11), hstv(i,12)
      enddo
      goto 310
 320  continue
      write(n6,'(2x,"plhist  find end  loop =",i7)') j
      goto 350
 330  continue
      write(n6,'(2x,"plhist  find err")')
 350  continue
      close(nft)
      return
      endif
!
      end
!
!***********************************************************************
      subroutine plhsav
!***********************************************************************
      use cntmnt, only : pfl_ion, psm_abs, psm_man, psm_pmp, tnpuf
      use cplcom, only : edtq1, edtq2, edtq3, edtq4, trad, vna, vne, vte
     >    , vti
      use cplhst, only : hstt, hstv, iphs, jphs, nhtm
      use cplqcn, only : pcn_pdt
      use csonic, only : dtim, itend, itim, time
      use cunit,  only : lmspe, lmype
      implicit none
!
!::local variables
      integer  ih, ii, k, ia, nst, nen, n, j, i
!
!::KSFUJI
!KH      if( mype.ne.mspe ) return
      if( lmype /= lmspe ) return
      if( itim.gt.itend )  return
!
      ih = nhtm
!
      hstt(ih) = real(time,kind=4)
!
      ii = 0
!
!::time
      ii = ii + 1
      hstv(ih,ii) = real(time,kind=4)
      ii = ii + 1
      hstv(ih,ii) = real(dtim,kind=4)
      ii = ii + 1
      hstv(ih,ii) = 99.0
      ii = ii + 1
      hstv(ih,ii) = 99.0
!
!::residual
      do 110 k = 1, 4
      ii = ii + 1
      if( k.eq.1 ) hstv(ih,ii) = real(edtq1,kind=4)
      if( k.eq.2 ) hstv(ih,ii) = real(edtq2,kind=4)
      if( k.eq.3 ) hstv(ih,ii) = real(edtq3,kind=4)
      if( k.eq.4 ) hstv(ih,ii) = real(edtq4,kind=4)
 110  continue
!
!::plasma parameter
      ia  = 1
      nst = 3
      nen = 14
      do 210 n = nst, nen
      j = jphs(n)
      i = iphs(n)
      do 220 k = 1, 4
      ii = ii + 1
      if( k.eq.1 ) hstv(ih,ii) = real(vne(j,i),kind=4)
      if( k.eq.2 ) hstv(ih,ii) = real(vna(j,i,ia),kind=4) ! vva => vna
      if( k.eq.3 ) hstv(ih,ii) = real(vti(j,i),kind=4)
      if( k.eq.4 ) hstv(ih,ii) = real(vte(j,i),kind=4)
 220  continue
 210  continue
!
!::Flux
      ia = 1
      ii = ii + 1; hstv(ih,ii) = real(-pfl_ion(5)+tnpuf, kind=4)
      ii = ii + 1; hstv(ih,ii) = real(pcn_pdt, kind=4)
      ii = ii + 1; hstv(ih,ii) = real(psm_abs+psm_man, kind=4)
      ii = ii + 1; hstv(ih,ii) = real(psm_pmp, kind=4)
!
!::Wrad   ! KSFUJI
      ii = ii + 1
      hstv(ih,ii) = real(-trad/1.0d6, kind=4)
!
!::plasma parameter at spx (Ne,Na,Ti,Te)
      ia = 1
      nst = 17
      nen = 18
      do n = nst, nen
        j = jphs(n)
        i = iphs(n)
        do k = 1, 4
          ii = ii + 1
          if( k.eq.1 ) hstv(ih,ii) = real(vne(j,i), kind=4)
          if( k.eq.2 ) hstv(ih,ii) = real(vna(j,i,ia), kind=4)
          if( k.eq.3 ) hstv(ih,ii) = real(vti(j,i), kind=4)
          if( k.eq.4 ) hstv(ih,ii) = real(vte(j,i), kind=4)
        enddo
      enddo
!
!::debug write
!x      write(n6,'(2x,"***  plhsav ***  ii =",i3)') ii
!x      write(n6,'(2x,1p10e12.3)') (hstv(ih,i),i=1,ii)
!
      return
      end
