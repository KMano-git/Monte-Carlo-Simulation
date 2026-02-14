!
!----------------------------------------------------------------------
!::atomic data of degas2-code    nanasvr/j5253/degas2
!
!      H0 + e -> H+ + e + e   ionization
!      H+ + e + e -> H0 + e   recombination
!
!  sghio :  <sig*v>(Te,ne) [cm3/w]
!  elhio :   Eloss(Te,ne)  [eV]
!  emhio :   ndot(H_alpha) [cm3/s]    2009/4/14 shimizu-san, 160603KH
!
!  electron energy source
!     We = - elhio*cev*Ne*N0*sghio - elhrc*cev*NH*Ne*sghrc
!          +  13.6*cev*NH*Ne*sghrc
!
!----------------------------------------------------------------------
!x  inc/soldori/catmhy
!
!x      parameter (ndtm=60,nddn=15)
!x      common /catmhy/ stemin,stemax,stedlt,snemin,snemax,snedlt
!x     >  ,ste(ndtm),sne(nddn)
!x     >  ,sghio(ndtm,nddn),elhio(ndtm,nddn),emhio(ndtm,nddn)
!x     >  ,sghrc(ndtm,nddn),elhrc(ndtm,nddn),emhrc(ndtm,nddn)
!x
!x  inc/soldori/csonic
!x      integer  lmype, lmspe, lnope, mywld, n6
!
!
!**********************************************************************
      subroutine atomhy
!**********************************************************************
      use catmhy,      only : elhio, elhrc, emhio, emhrc, nddn, ndtm
     >    , sghio, sghrc, sne, snedlt, snemax, snemin, ste, stedlt
     >    , stemax, stemin
      use cunit,       only : lmspe, lmype, lnope, mywld
      use mod_dtypdef, only : nddt, typ_dnam, typ_itag
      use mod_shexe,   only : use_namelist, DEGAS_ATMHY
      use mpi!,         only : mpi_bcast
      use cfatmhy, only : latomhcr
      implicit none
!::
      real*8    dat(ndtm,nddn)
      real*8    cev_mk, eion_ev, cev_cg, eion_cg, eps
      real*8    zne, zwe, e32, a32, zemi
      integer   nft, jtmx, jt, jnmx, jn, nout
!****** use, jtmx => ndtm, jnmx => nddn, scheduled to be fixed *****
      character dsn*80, ctag*80, ckey*10
      logical lex
      integer  ldbg
!::
      integer  ierr, itag, ityp

!
      ldbg = 0   ! atomic data
!
!::const
      cev_mk  = 1.60217733d-19         ! [eV] in MKS-unit
      eion_ev = 13.595d0               ! [eV]
      cev_cg  = cev_mk*1.0d7           ! [eV] in CGS-unit
      eion_cg = eion_ev*cev_cg         ! [erg]
      eps = 1.0d-8
!
!::file
      nft = 21
      if(use_namelist) then
!.. use variable in inppls
         call nopen(nft,DEGAS_ATMHY,"text",dsn,lex)
      else
!.. use environment variable
         call nopen(nft,"DEGAS_ATMHY","text",dsn,lex)
      end if
!
!::ste
      jtmx = ndtm
      do 110 jt = 1, jtmx
      ste(jt) = 10.0d0**(-1.2d0+0.1d0*dfloat(jt-1))
 110  continue
      stemin = ste(1) + eps*(ste(2)-ste(1))
      stemax = ste(jtmx) - eps*(ste(jtmx)-ste(jtmx-1))
      stedlt = 10.0d0**0.1d0
!
!::sne
      jnmx = nddn
      do 120 jn = 1, jnmx
      sne(jn) = 10.0d0**(10.0d0+0.5d0*(jn-1))
 120  continue
      snemin = sne(1) + eps*(sne(2)-sne(1))
      snemax = sne(jnmx) - eps*(sne(jnmx)-sne(jnmx-1))
      snedlt = 10.0d0**0.5d0
!
!--MASTER PE--start-------------------------------------------
      if( lmype.eq.lmspe ) then
!
!----------------------------------------------------------------------
!::Ionization  reaction rate
!----------------------------------------------------------------------
      ctag = "Ionization Rate"
      call flread( nft, ctag, sghio )
!
!----------------------------------------------------------------------
!::Ionization  background energy loss rate
!----------------------------------------------------------------------
      ctag = "Neutral"//" Electron Losses"
      call flread( nft, ctag, dat )
!
! [erg/s] ==> [erg*cm3/s]
! Add ionization losses to exciation losses
      do 220 jn = 1, jnmx
      zne = sne(jn)
      do 225 jt = 1, jtmx
      zwe = dat(jt,jn)/zne + eion_cg*sghio(jt,jn)
      elhio(jt,jn) = zwe/sghio(jt,jn)/cev_cg
 225  continue
 220  continue
!
!----------------------------------------------------------------------
!::Ionization  emission rate
!----------------------------------------------------------------------
      ctag = "Neutral"//" n=3 / n=1"
      call flread( nft, ctag, dat )
!
! |e32| = transition energy from n=3 to n=2 (H_alpha)
! |a32| = radiative decay rate for this process ("Einstein coefficient")
! ionization/emission [event/photon]
!    Halpha 1 photon = e32 = 3.025d-12 [erg]
!
      e32 = eion_cg*(1.0d0/2.0d0**2-1.0d0/3.0d0**2)
      a32 = 4.41d7
!
      do 230 jn = 1, jnmx
      zne = sne(jn)
      do 235 jt = 1, jtmx
      zemi = dat(jt,jn)*e32*a32/zne
      emhio(jt,jn) = zemi/e32
 235  continue
 230  continue
!
!::file
      rewind nft
!
!----------------------------------------------------------------------
!::Recombination  reaction rate
!----------------------------------------------------------------------
      ctag = "Recombination Rate"
      call flread( nft, ctag, sghrc )
!
!----------------------------------------------------------------------
!::Recombination  background energy loss rate
!----------------------------------------------------------------------
      ctag = "Continuum"//" Electron Losses"
      call flread( nft, ctag, dat )
!
! [erg/s] ==> [erg*cm3/s]
! [erg*cm3/s] ==> [eV]
      do 320 jn = 1, jnmx
      zne = sne(jn)
      do 325 jt = 1, jtmx
      elhrc(jt,jn) = dat(jt,jn)/zne/sghrc(jt,jn)/cev_cg
 325  continue
 320  continue
!
!----------------------------------------------------------------------
!::Recombination  emission rate
!----------------------------------------------------------------------
      ctag = "Continuum"//" n=3 / n=1"
      call flread( nft, ctag, dat )
!
! |e32| = transition energy from n=3 to n=2 (H_alpha)
! |a32| = radiative decay rate for this process ("Einstein coefficient")
! ionization/emission [event/photon]
!    Halpha 1 photon = e32 = 3.025d-12 [erg]
!
      e32 = eion_cg*(1.0d0/2.0d0**2-1.0d0/3.0d0**2)
      a32 = 4.41d7
!
      do 330 jn = 1, jnmx
      zne =sne(jn)
      do 335 jt = 1, jtmx
      zemi = dat(jt,jn)*e32*a32/zne
      emhrc(jt,jn) = zemi/e32
 335  continue
 330  continue
!
!----------------------------------------------------------------------
!::debug
!----------------------------------------------------------------------
      if( ldbg.eq.1 ) then
      nout = 23
!--ionH_sigv
      ckey = "ionH_sigv"
      open( unit=nout, file='atm_'//ckey(1:9), form ='formatted' )
      call lstrdt( nout, ckey, sghio, ste, sne, jtmx, jnmx )
      close( nout )
!--ionH_Eels
      ckey = "ionH_Eels"
      open( unit=nout, file='atm_'//ckey(1:9), form ='formatted' )
      call lstrdt( nout, ckey, elhio, ste, sne, jtmx, jnmx )
      close( nout )
!--ionH_fcem
      ckey = "ionH_fcem"
      open( unit=nout, file='atm_'//ckey(1:9), form ='formatted' )
      call lstrdt( nout, ckey, emhio, ste, sne, jtmx, jnmx )
      close( nout )
!--recH_sigv
      ckey = "recH_sigv"
      open( unit=nout, file='atm_'//ckey(1:9), form ='formatted' )
      call lstrdt( nout, ckey, sghrc, ste, sne, jtmx, jnmx )
      close( nout )
!--recH_Eels
      ckey = "recH_Eels"
      open( unit=nout, file='atm_'//ckey(1:9), form ='formatted' )
      call lstrdt( nout, ckey, elhrc, ste, sne, jtmx, jnmx )
      close( nout )
!--recH_fcem
      ckey = "recH_fcem"
      open( unit=nout, file='atm_'//ckey(1:9), form ='formatted' )
      call lstrdt( nout, ckey, emhrc, ste, sne, jtmx, jnmx )
      close( nout )
      endif
!
!--MASTER PE--end----------------------------------------
      close( nft )
      endif
!
!----------------------------------------------------------------------
!::send data
!----------------------------------------------------------------------
      if( lnope.gt.1 ) then
!::[MPI_Bcast in atomhy]    catmhy  (stemin-emhrc)  10/04/21
      call tbfind( nddt, typ_dnam, ityp, 'CATMHY' )
      itag = typ_itag(ityp)
      call MPI_Bcast( stemin, 1, itag, lmspe, mywld, ierr )
      endif
!
!----------------------------------------------------------------------
!::molecular reaction rate by amjuel
!    overwrite in svhion, elhion 20211112 SY
!----------------------------------------------------------------------
      if(latomhcr==1) call set_atomhcr

!----------------------------------------------------------------------
!::log
!----------------------------------------------------------------------
      stemin = dlog(stemin)
      stemax = dlog(stemax)
      stedlt = dlog(stedlt)
      snemin = dlog(snemin)
      snemax = dlog(snemax)
      snedlt = dlog(snedlt)
!
      do 410 jt = 1, jtmx
      ste(jt) = dlog(ste(jt))
 410  continue
      do 420 jn = 1, jnmx
      sne(jn) = dlog(sne(jn))
 420  continue
!
      do 430 jn = 1, jnmx
      do 435 jt = 1, jtmx
      sghio(jt,jn) = dlog(sghio(jt,jn))
      sghrc(jt,jn) = dlog(sghrc(jt,jn))
      elhio(jt,jn) = dlog(elhio(jt,jn))
      elhrc(jt,jn) = dlog(elhrc(jt,jn))
      emhio(jt,jn) = dlog(emhio(jt,jn))
      emhrc(jt,jn) = dlog(emhrc(jt,jn))
 435  continue
 430  continue
!
!----------------------------------------------------------------------
!::return
!----------------------------------------------------------------------
 900  continue
      return
      end
!
!**********************************************************************
      subroutine flread( nft, tag, dat )
!**********************************************************************
      use catmhy, only : nddn, ndtm
      use cunit,  only : n6
      implicit none
!
      integer,   intent(in)  :: nft
      character, intent(in)  :: tag*(*)
      real(8),   intent(out) :: dat(ndtm,nddn)
!
      integer  i, mj, jtmx, ii, m, jt, jn
      character clin*80, cmsg*80
      integer   ldbg
!
      ldbg = 0
!
      do 10 i = len(tag), 1, -1
      mj = i
      if( tag(i:i).ne." " ) goto 20
 10   continue
      mj = 1
      goto 190
 20   continue
      if( ldbg.eq.1 ) then
      write(n6,'(/1x,"flread   key word ",i3,2x,a)') mj,tag(1:mj)
      endif
!
 110  continue
      read( nft, '(a)', end=190 ) clin
      if( index( clin, tag(1:mj) ).eq.0 ) goto 110
!
!::find tag
      cmsg = clin
      jtmx = 60
      ii = 0
 210  continue
      read( nft, '(a)', end=192 ) clin
      m = index( clin, "jn = " )
      if( m.eq.0 ) goto 250
      ii = ii + 1
!
!::get data
      read( nft, * ) (dat(jt,ii),jt=1,jtmx)
      read( nft, '(a)' ) clin
      goto 210
!
 250  continue
!
!::debug
      if( ldbg.eq.1 ) then
      write(n6,'(a)') cmsg
      do 310 jn = 1, ii
      write(n6,'(1x,"jn = ",i3)') jn
      write(n6,'(1p6e13.5)') (dat(jt,jn),jt=1,jtmx)
 310  continue
      endif
!
      return
!
 190  continue
      call wexit("flread","no found tag : "//tag(1:mj))
!
 192  continue
      call wexit("flread","no found key word jn")
!
      end
!
!**********************************************************************
      subroutine lstrdt( nout, ckey, dat, ate, ane, jtmx, jnmx )
!**********************************************************************
      use catmhy, only : nddn, ndtm
      implicit none
!
      integer,   intent(in) :: nout
      character, intent(in) :: ckey*(*)
      real(8),   intent(in) :: dat(ndtm,nddn), ate(ndtm), ane(nddn)
      integer,   intent(in) :: jtmx, jnmx
!
      real*8  zdt(nddn), zmin
      integer jn, jt
!
      zmin = 1.0d-50
      write(nout,'(2x,a,2x,a,8x,1p15e12.4)')
     >   ckey,"Te",(ane(jn),jn=1,jnmx)
      do 210 jt = 1, jtmx
      do 220 jn = 1, jnmx
      zdt(jn) = dat(jt,jn)
      if( dabs(zdt(jn)).lt.zmin ) zdt(jn) = zmin
 220  continue
      write(nout,'(1x,i4,6x,1p16e12.4)')
     >   jt,ate(jt),(zdt(jn),jn=1,jnmx)
 210  continue
!
      return
      end
!
!**********************************************************************
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   function svhion( ate, ane )
      real(8) function svhion( ate, ane )
!**********************************************************************
      use catmhy, only : sghio, snedlt, snemax, snemin, stedlt, stemax
     >    , stemin
      use cfatmhy, only : latomhcr
      implicit none
!
      real(8), intent(in) :: ate, ane
!
      real*8   zte, zne, zlte, zite, zdte, zlne, zine, zdne, zfn
      integer  itea, iteb, inea, ineb
      real*8   sv_atomhcr
!
!::MKS => CGS
      zte = ate
      zne = ane*1.0d-6
!
!::index of Te
      zlte = dlog(zte)
      zlte = dmax1( zlte, stemin )
      zlte = dmin1( zlte, stemax )
      zite = (zlte-stemin)/stedlt + 1.0d0
      itea = int(zite)
      iteb = itea+1
      zdte = zite - dfloat(itea)
!
!::index of Ne
      zlne = dlog(zne)
      zlne = dmax1( zlne, snemin )
      zlne = dmin1( zlne, snemax )
      zine = (zlne-snemin)/snedlt + 1.0d0
      inea = int(zine)
      ineb = inea+1
      zdne = zine - dfloat(inea)
!
!::bi-linear interpolation
      zfn  = sghio(itea,inea)*(1.0d0-zdte)*(1.0d0-zdne) +
     >       sghio(iteb,inea)*zdte*(1.0d0-zdne) +
     >       sghio(iteb,ineb)*zdte*zdne +
     >       sghio(itea,ineb)*(1.0d0-zdte)*zdne
      zfn  = dexp(zfn)
!
!::CGS => MKS
      svhion = zfn*1.0d-6

      if(latomhcr==1) svhion = sv_atomhcr(ate,ane,1) !20211112 SY

      return
      end
!
!**********************************************************************
      real(8) function elhion( ate, ane )
!**********************************************************************
      use catmhy, only : elhio, snedlt, snemax, snemin, stedlt, stemax
     >    , stemin
      use cfatmhy, only : latomhcr
      implicit none
!
      real(8), intent(in) :: ate, ane
!
      real*8  zte, zne, zlte, zite, zdte, zlne, zine, zdne, zfn
      integer itea, iteb, inea, ineb
      real*8   sv_atomhcr
!
!::MKS => CGS
      zte = ate
      zne = ane*1.0d-6
!
!::index of Te
      zlte = dlog(zte)
      zlte = dmax1( zlte, stemin )
      zlte = dmin1( zlte, stemax )
      zite = (zlte-stemin)/stedlt + 1.0d0
      itea = int(zite)
      iteb = itea+1
      zdte = zite - dfloat(itea)
!
!::index of Ne
      zlne = dlog(zne)
      zlne = dmax1( zlne, snemin )
      zlne = dmin1( zlne, snemax )
      zine = (zlne-snemin)/snedlt + 1.0d0
      inea = int(zine)
      ineb = inea+1
      zdne = zine - dfloat(inea)
!
!::bi-linear interpolation
      zfn  = elhio(itea,inea)*(1.0d0-zdte)*(1.0d0-zdne) +
     >       elhio(iteb,inea)*zdte*(1.0d0-zdne) +
     >       elhio(iteb,ineb)*zdte*zdne +
     >       elhio(itea,ineb)*(1.0d0-zdte)*zdne
      zfn  = dexp(zfn)
!
!::CGS => MKS  [eV]
      elhion = zfn
!      
      if(latomhcr==1)then
!         elhion = sv_atomhcr(ate,ane,2)/sv_atomhcr(ate,ane,1)
      elhion = 13.6d0
          !20211112 SY
      endif

      return
      end
!
!**********************************************************************
      real(8) function svhrcm( ate, ane )
!**********************************************************************
      use catmhy, only : sghrc, snedlt, snemax, snemin, stedlt, stemax
     >    , stemin
      implicit none
!
      real(8), intent(in) :: ate, ane

      real*8  zte, zne, zlte, zite, zdte, zlne, zine, zdne, zfn
      integer itea, iteb, inea, ineb
!
!::MKS => CGS
      zte = ate
      zne = ane*1.0d-6
!
!::index of Te
      zlte = dlog(zte)
      zlte = dmax1( zlte, stemin )
      zlte = dmin1( zlte, stemax )
      zite = (zlte-stemin)/stedlt + 1.0d0
      itea = int(zite)
      iteb = itea+1
      zdte = zite - dfloat(itea)
!
!::index of Ne
      zlne = dlog(zne)
      zlne = dmax1( zlne, snemin )
      zlne = dmin1( zlne, snemax )
      zine = (zlne-snemin)/snedlt + 1.0d0
      inea = int(zine)
      ineb = inea+1
      zdne = zine - dfloat(inea)
!
!::bi-linear interpolation
      zfn  = sghrc(itea,inea)*(1.0d0-zdte)*(1.0d0-zdne) +
     >       sghrc(iteb,inea)*zdte*(1.0d0-zdne) +
     >       sghrc(iteb,ineb)*zdte*zdne +
     >       sghrc(itea,ineb)*(1.0d0-zdte)*zdne
      zfn  = dexp(zfn)
!
!::CGS => MKS  [cm3/s]
      svhrcm = zfn*1.0d-6
!
      return
      end
!
!**********************************************************************
      real(8) function elhrcm( ate, ane )
!**********************************************************************
      use catmhy, only : elhrc, snedlt, snemax, snemin, stedlt, stemax
     >    , stemin
      implicit none
!
      real(8), intent(in) :: ate, ane
!
      real*8  zte, zne, zlte, zite, zdte, zlne, zine, zdne, zfn
      integer itea, iteb, inea, ineb
!
!::MKS => CGS
      zte = ate
      zne = ane*1.0d-6
!
!::index of Te
      zlte = dlog(zte)
      zlte = dmax1( zlte, stemin )
      zlte = dmin1( zlte, stemax )
      zite = (zlte-stemin)/stedlt + 1.0d0
      itea = int(zite)
      iteb = itea+1
      zdte = zite - dfloat(itea)
!
!::index of Ne
      zlne = dlog(zne)
      zlne = dmax1( zlne, snemin )
      zlne = dmin1( zlne, snemax )
      zine = (zlne-snemin)/snedlt + 1.0d0
      inea = int(zine)
      ineb = inea+1
      zdne = zine - dfloat(inea)
!
!::bi-linear interpolation
      zfn  = elhrc(itea,inea)*(1.0d0-zdte)*(1.0d0-zdne) +
     >       elhrc(iteb,inea)*zdte*(1.0d0-zdne) +
     >       elhrc(iteb,ineb)*zdte*zdne +
     >       elhrc(itea,ineb)*(1.0d0-zdte)*zdne
      zfn  = dexp(zfn)
!
!::CGS => MKS  [eV]
      elhrcm = zfn
!
      return
      end
!
!**********************************************************************
      real(8) function emhion( ate, ane )
!**********************************************************************
      use catmhy, only : emhio, snedlt, snemax, snemin, stedlt, stemax
     >    , stemin
      implicit none
!
      real(8), intent(in) :: ate, ane
!
      real*8  zte, zne, zlte, zite, zdte, zlne, zine, zdne, zfn
      integer itea, iteb, inea, ineb
!
!::MKS => CGS
      zte = ate
      zne = ane*1.0d-6
!
!::index of Te
      zlte = dlog(zte)
      zlte = dmax1( zlte, stemin )
      zlte = dmin1( zlte, stemax )
      zite = (zlte-stemin)/stedlt + 1.0d0
      itea = int(zite)
      iteb = itea+1
      zdte = zite - dfloat(itea)
!
!::index of Ne
      zlne = dlog(zne)
      zlne = dmax1( zlne, snemin )
      zlne = dmin1( zlne, snemax )
      zine = (zlne-snemin)/snedlt + 1.0d0
      inea = int(zine)
      ineb = inea+1
      zdne = zine - dfloat(inea)
!
!::bi-linear interpolation
      zfn  = emhio(itea,inea)*(1.0d0-zdte)*(1.0d0-zdne) +
     >       emhio(iteb,inea)*zdte*(1.0d0-zdne) +
     >       emhio(iteb,ineb)*zdte*zdne +
     >       emhio(itea,ineb)*(1.0d0-zdte)*zdne
      zfn  = dexp(zfn)
!
!::CGS => MKS  [m3/s]
      emhion = zfn*1.0d-6
!
      return
      end
!
!**********************************************************************
      real(8) function emhrcm( ate, ane )
!**********************************************************************
      use catmhy, only : emhrc, snedlt, snemax, snemin, stedlt, stemax
     >    , stemin
      implicit none
!
      real(8), intent(in) :: ate, ane
!
      real*8  zte, zne, zlte, zite, zdte, zlne, zine, zdne, zfn
      integer itea, iteb, inea, ineb
!
!::MKS => CGS
      zte = ate
      zne = ane*1.0d-6
!
!::index of Te
      zlte = dlog(zte)
      zlte = dmax1( zlte, stemin )
      zlte = dmin1( zlte, stemax )
      zite = (zlte-stemin)/stedlt + 1.0d0
      itea = int(zite)
      iteb = itea+1
      zdte = zite - dfloat(itea)
!
!::index of Ne
      zlne = dlog(zne)
      zlne = dmax1( zlne, snemin )
      zlne = dmin1( zlne, snemax )
      zine = (zlne-snemin)/snedlt + 1.0d0
      inea = int(zine)
      ineb = inea+1
      zdne = zine - dfloat(inea)
!
!::bi-linear interpolation
      zfn  = emhrc(itea,inea)*(1.0d0-zdte)*(1.0d0-zdne) +
     >       emhrc(iteb,inea)*zdte*(1.0d0-zdne) +
     >       emhrc(iteb,ineb)*zdte*zdne +
     >       emhrc(itea,ineb)*(1.0d0-zdte)*zdne
      zfn  = dexp(zfn)
!
!::CGS => MKS  [m3/s]
      emhrcm = zfn*1.0d-6
!
      return
      end
