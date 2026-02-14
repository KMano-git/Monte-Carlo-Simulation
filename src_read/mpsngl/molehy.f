!----------------------------------------------------------------------
!::molecular data of degas2-code
!
!   original-data   nanasvr   ~j5253/degas2.2.6/degas/data/h2dis.nc
!   cdf-data        nanasvr   ~j3859/dtdegas2/cdfdt/h2dis.cdf
!                   comet    shimizu/dtdegas2/cdfdt/h2dis.cdf
!
!   compact-data    comet    shimizu/dtdegas2/cdfdt/molhy.dat
!                            "/home/shimizu/dat/molhy.dat"
!
!  1    h2dis      [e + H2 -> e + H + H]
!  2    h2ion      [e + H2 -> e + H2+ + e]
!  3    h2dision   [e + H2 -> e + H+ + H + e]
!  4    h2pdis     [e + H2+ -> e + H+ + H]
!  5    h2pdision  [e + H2+ -> e + H+ + H+ + e]
!  6    h2pdisrec  [e + H2+ -> H + H *]
!                  123456789012345678901234567890
!
!    sgvm(ite,irc) :  <sig*v>(Te)  [cm3/s]  reaction rate
!    stem(ite)     :  Te           [eV]     electron temperature
!    elsm(irc)     :  Eloss(Te)    [eV]     energy loss
!    engm(irc)     :  E0           [eV]     dissociation energy
!
!**********************************************************************
      subroutine test_molehy
!**********************************************************************
      use cmolhy, only : cpid, nrc
      use cunit,  only : lmspe, lmype, n6
      implicit none
!
      real*8  ate1, ate2, dlt1, dlt2, dlt, ate, ane
      integer nmax, i, irc
!
      integer n61,n62,n63,n64
      data n61,n62,n63,n64/41,42,43,44/
!
      real*8 asgv(6),aels(6),aeng(6)
!
!::mpi variables for nopen routine
      lmype = 0
      lmspe = 0
      n6 = 6
!
      call molehy
!
      ane = 1.0d16
      ate1 = 0.1d0
      ate2 = 1.0d4
      dlt1 = dlog(ate1)
      dlt2 = dlog(ate2)
      nmax = 501
!
      open(unit=n61,file="moldt_sgv.dat",form='formatted')
      open(unit=n62,file="moldt_els.dat",form='formatted')
      open(unit=n63,file="moldt_eng.dat",form='formatted')
      open(unit=n64,file="moldt_wls.dat",form='formatted')
!
      write(n61,'(2x,"sgv_j",2x,"Te",10x,6(a10,2x))')
     >  (cpid(i),i=1,nrc)
      write(n62,'(2x,"els_j",2x,"Te",10x,6(a10,2x))')
     >  (cpid(i),i=1,nrc)
      write(n63,'(2x,"eng_j",2x,"Te",10x,6(a10,2x))')
     >  (cpid(i),i=1,nrc)
      write(n64,'(2x,"wls_j",2x,"Te",10x,6(a10,2x))')
     >  (cpid(i),i=1,nrc)
!
      do 310 i = 1, nmax
      dlt = dlt1 + (dlt2-dlt1)/dfloat(nmax-1)*dfloat(i-1)
      ate = dexp(dlt)
      call sgvmhy( ate, ane, asgv, aels, aeng )
!
      write(n61,'(2x,i4,1p7e12.3)') i,ate,(asgv(irc)*1.0d6,irc=1,6)
      write(n62,'(2x,i4,1p7e12.3)') i,ate,(aels(irc),irc=1,6)
      write(n63,'(2x,i4,1p7e12.3)') i,ate,(aeng(irc),irc=1,6)
      write(n64,'(2x,i4,1p7e12.3)') i,ate,(aels(irc)*asgv(irc)*1.0d6,
     >    irc=1,6)
 310  continue
!
      stop
      end
!
!**********************************************************************
      subroutine molehy
!**********************************************************************
      use cmolhy,      only : cpid, cract, els3, elsm, eng3, engm, jt3e
     >    , jt3l, jtmx, lmolhcr, ndtmm, nrc, sgvm, stedlt, stem, stemax
     >    , stemin, te3e, te3l
      use cunit,       only : lmspe, lmype, n6
      use mod_shexe,   only: use_namelist, DEGAS_MOLHY
      implicit none
!
      integer nft, i, mj, j, irc, mje, mjs, num
      real*8  eps, zelos, zdeng
      logical lex
      character dsn*80, clin*120
      character ckey1*7, ckey2*14
      real(8)  zte(ndtmm)

! function
      integer  lenx

      call trmark("molehy","start")
!
!--MASTER PE--start-------------------------------------------
      if( lmype.eq.lmspe ) then
!
!::file
      nft = 21
      if(use_namelist) then
!..   use variable in inppls
         call nopen(nft,DEGAS_MOLHY,"text",dsn,lex)
      else
!..   use environment variable
         call nopen(nft,"DEGAS_MOLHY","text",dsn,lex)
      end if
!
      write(n6,'(/2x,"*** molehy ***")')
      write(n6,'(4x,"file =",a)') dsn(1:lenx(dsn))
!
!::ste
      do 110 i = 1, 9
      read(nft,'(a)') clin
 110  continue
      mj = index(clin,"Te")
      read(clin(mj+3:),*) jtmx
      read(nft,*) (zte(j),j=1,jtmx)
      do 120 j = 1, jtmx
      stem(j) = 10.0d0**(0.1d0*dfloat(j-1))
 120  continue
      eps = 1.0d-8
      stemin = stem(1) + eps*(stem(2)-stem(1))
      stemax = stem(jtmx) - eps*(stem(jtmx)-stem(jtmx-1))
      stedlt = 10.0d0**0.1d0
!
      write(n6,'(4x,"jtmx =",i3,"  temin =",1pe12.3,"  temax =",1pe12.3
     >  )') jtmx,zte(1),zte(jtmx)
      write(n6,'(4x,2x,3x,"reaction",23x,"elsm",4x,1x,"  engm",4x,
     >  "  stem =",1p15e11.3)')  (stem(j),j=1,jtmx,10),stem(jtmx)
!
!::file
      irc = 0
 200  continue
      read(nft,'(a)',end=290) clin
      if( index(clin,"file").eq.0 ) goto 200
      mje = index(clin,".cdf") - 1
      do 205 j = mje,1,-1
      mjs = j
      if( index(clin(j:j),"/").gt.0 ) goto 207
 205  continue
 207  continue
      irc  = irc + 1
      cpid(irc) = clin(mjs+1:mje)
!
!::reaction
      read(nft,'(a)') clin
      if( index(clin,"reaction").gt.0 ) then
      mj = index(clin,"[")
      cract(irc) = clin(mj:)
      endif
!
!::elos
      read(nft,'(a)') clin
      if( index(clin,"elos").gt.0 ) then
      read(clin,'(a7,2x,i3,2x,a14)') ckey1,num,ckey2
      if( num.eq.1 ) then
        read(ckey2,*) zelos
      else
        if( irc.ne.3 ) goto 910
        read(nft,*) (te3l(j),j=1,num)
        read(nft,*) (els3(j),j=1,num)
        jt3l = num
        zelos = els3(1)
      endif
      elsm(irc) = zelos
      endif
!
!::deng
      read(nft,'(a)') clin
      if( index(clin,"deng").gt.0 ) then
      read(clin,'(a7,2x,i3,2x,a14)') ckey1,num,ckey2
      if( num.eq.1 ) then
        read(ckey2,*) zdeng
      else
        if( irc.ne.3 ) goto 920
        read(nft,*) (te3e(j),j=1,num)
        read(nft,*) (eng3(j),j=1,num)
        zdeng = eng3(1)
        jt3e = num
      endif
      engm(irc) = zdeng
      endif
!
!::rate
      read(nft,'(a)') clin
      if( index(clin,"rate").gt.0 ) then
      read(clin,'(a7,2x,i3,2x,a14)') ckey1,num,ckey2
      read(nft,*) (sgvm(j,irc),j=1,num)
      endif
!
!::debug write
      write(n6,'(4x,i2,2x,a28,2x,1pe10.3,1x,1pe10.3"  sgvm ="
     > ,1p15e11.3)') irc,cract(irc)(1:lenx(cract(irc)))
     > ,elsm(irc),engm(irc),(sgvm(j,irc),j=1,jtmx,10),sgvm(jtmx,irc)
!
!::loop
      goto 200
!
 290  continue
      nrc = irc
!
!::log
      stemin = dlog(stemin)
      stemax = dlog(stemax)
      stedlt = dlog(stedlt)
      do 410 j = 1, jtmx
      stem(j) = dlog(stem(j))
      do 420 irc = 1, nrc
      sgvm(j,irc) = dlog(sgvm(j,irc))
 420  continue
 410  continue
!
!--MASTER PE--end----------------------------------------
      close( nft )
      endif

!----------------------------------------------------------------------
!::molecular reaction rate by amjuel
!    overwrite in sgvmhy
!----------------------------------------------------------------------
      if(lmolhcr==1) call set_molhcr

      call trmark("molehy","return")
      return
!
 910  continue
      write(clin,'("elos  irc,num =",2i4)') irc,num
      call wexit("molehy",clin)
!
 920  continue
      write(clin,'("deng  irc,num =",2i4)') irc,num
      call wexit("mokehy",clin)
!
      end
!
!**********************************************************************
      subroutine sgvmhy( ate, ane, asgv, aels, aeng )
!**********************************************************************
      use cmolhy, only : els3, elsm, eng3, engm, jt3e, jt3l, lmolhcr
     >    , sgvm, stedlt, stemax, stemin, te3e, te3l
      implicit none
!
      real(8), intent(in)  :: ate, ane
      real(8), intent(out) :: asgv(6), aels(6), aeng(6)
!
      real(8)  zte, zne, zite, zlte, zdte, zfn
      integer  itea, iteb, irc
! external function
      real(8)  dintp, sv_molhcr
!
!::MKS => CGS
      zte = ate
      zne = ane
!
!::index of Te
      zlte = dlog(zte)
      zlte = dmax1( zlte, stemin )
      zlte = dmin1( zlte, stemax )
      zite = (zlte-stemin)/stedlt + 1.0d0
      itea = zite
      iteb = itea+1
      zdte = zite - dfloat(itea)
!
!::linear interpolation
!::CGS => MKS [m3/s]
      do 110 irc = 1, 6
      zfn  = sgvm(itea,irc)*(1.0d0-zdte) + sgvm(iteb,irc)*zdte
      zfn  = dexp(zfn)
      asgv(irc) = zfn*1.0d-6
 110  continue
!
      if(lmolhcr==1)then
      do irc = 1, 6
        asgv(irc) = sv_molhcr(zte,zne,irc)
      enddo
      endif
!
!::els,eng
      do 120 irc = 1, 6
      aels(irc) = elsm(irc)
      aeng(irc) = engm(irc)
 120  continue
!
!::h2dision
      irc = 3
      aels(irc) = dintp(ate,jt3l,te3l,els3)
      aeng(irc) = dintp(ate,jt3e,te3e,eng3)
!
!::h2pdisrec
      irc = 6
      aels(irc) = aels(irc)*ate
!
      return
      end
