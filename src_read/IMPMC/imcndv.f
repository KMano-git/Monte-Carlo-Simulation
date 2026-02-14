!**********************************************************************
      subroutine imcndv
!**********************************************************************
!
!   imcndv  im: impurity  cnd: conductance  v: vacume region
!
!     wvacn(jgi,jgo)
!       jgi : gate number when enter vac-region
!       jgo : out of vac-region   jgo4 = nogt + 4
!
!       wvacn(jgi,0) : weght to enter from gate(jgi)
!       wvacn(jgi,1) : exit from gate(1)
!       wvacn(jgi,2) : exit from gate(2)       2 = nogt
!       wvacn(jgi,3) : pump   ien = 7          3 = nogt + 1
!       wvacn(jgi,4) : absorption at vac-wall  4 = nogt + 2
!       wvacn(jgi,5) : survive                 5 = nogt + 3
!       wvacn(jgi,6) : stick  ien = 6          6 = nogt + 4
!
!  Note These data are scoring variables
!      wvacn(gi,1) wvacn(gi,2)  set back after sumup
!
!    /cimcom_12c/
!       integer, parameter :: ndext = ndgt + 5
!       real*8  wvacn(0:ndgt,0:ndext), wvaci(0:ndgt,0:ndext)
!       parameter (nwkmp_im4 = (ndgt+1)*(ndext+1)*2 )
!
!----------------------------------------------------------------------
      use cimcom,      only : bkflux, bkfrac, htout, ien, igtno, ir, is
     >    , itimz, lgstk, lhst, lpcm, lprf, npmax, ns4 => nwkmp_im4
     >    , sptyc, timez, tout, wght, wstp, wvaci, wvacn
      use cimctl,      only : icalz, lstdy
      use cntcom,      only : mrgn, nogt
      use csize,       only : ndgt
      use csonic,      only : itim
      use cunit,       only : lmspe, lmype, mype, mywld, n6
      use mod_dtypdef, only : nddt, typ_dnam, typ_itag
      use mod_shexe,   only : impmc_model
      use mpi!,         only : mpi_bcast, mpi_real8, mpi_reduce, mpi_sum
      implicit none
!
!::mpi variables
      real(8), dimension(ns4) :: wwork4, vwork4, twork4
!
!::local variables
      integer :: jp, jc, jrg, jgo3, jpm
      real(8) :: smgt0n, smgt0i, zwght
      real(8), dimension(1:ndgt) :: smcaln, smcali, smstkn, smstki
      integer :: ih60
      integer :: ierr, itag, ityp
      integer :: jgo, jgo4, jgi, jex, jab, i
      real(8), dimension(0:ndgt) :: zwern, zweri, zwton, zwtoi, zcndv
      real(8) :: zwghn, zwghi, zwoun, zwoui
      real(8) :: zsum, zfout, zero
!
!
!::classify test particles  for conductance  wvacn
!    flux onto wall in vac region  wvacn(0:4)  cal in sub. imwrfn
!    total weight of survive particles(5) and stick particles(6)
!
!::vac-conductance (see imwrfn)
!   wvacn(0,0) = Er2   a:absoprtin/in:exist in Vac/ stk:stick
!   wvacn(1,0) = g1-in, g1/g2-ou(1,1/2), a(+1) aV(+2) inV(+3) stk(+4)
!   wvacn(2,0) = g2-in, g1/g2-ou(2,1/2), a(+1) aV(+2) inV(+3) stk(+4)
!     stick : ien(ip) = 6 wght > 0 ,  cal : ien(ip) = 0
!
!   Note.  when stick at gate, ir = 4/5  not. 8
!          when exit from gate, igtno = 0
!
!
      if( impmc_model == 0 ) then
        write(n6,'(2x,"*** imcndv ***  lpcm =",i10,"  lprf =",i2,
     >   "  lhst =",i2)') lpcm, lprf, lhst
      else
!::GTst  lgbkf= 1 (cal. bk_mag)  lgstk= 0 (no stick)
!::GTab       = 0 (no use)            = 1 (stick at gate)
!::GTbk       = 0 (use  bk_mag)       = 1 (stick at gate)

        write(n6,'(/2x,a,2x,"lgbkf/lgstk =",2i3,2x,
     >    "itim/icalZ =",i6,i5)') "imcndv",0 ,lgstk,itim,icalZ
        if( lgstk == 1 ) return
      endif
!
      smgt0n = 0.0d0          ! no found gate number 0
      smgt0i = 0.0d0
      smcaln(1:nogt) = 0.0d0  ! calculation  nogt + 3
      smcali(1:nogt) = 0.0d0
      smstkn(1:nogt) = 0.0d0  ! stick  nogt + 4
      smstki(1:nogt) = 0.0d0
!
      do jp = 1, npmax
        jc  = ir(jp)
        if(jc.le.0) cycle ! escape error(imcndv should be call after imp_pack)
        jgi = igtno(jp)       ! gate number when enter into Vac region
        jrg = mrgn(jc)
        if( jrg.ne.8 ) cycle  ! Vac region
!
!::error: cannot enter Vac region without passing gate
        if( jgi.eq.0 ) then
          if( is(jp).eq.0 ) then
            smgt0n = smgt0n + wght(jp)
          else
            smgt0i = smgt0i + wght(jp)
          endif
          goto 120
        endif
!::cal
        if( ien(jp).eq.0 ) then
          if( is(jp).eq.0 ) then
            smcaln(jgi) = smcaln(jgi) + wght(jp)
          else
            smcali(jgi) = smcali(jgi) + wght(jp)
          endif
          goto 120
        endif
!::stick in Vac region
!:: wght =0.95 when stick => wght = 0.0  wstp = 0.95 for IMPMC_TD SY VV5.2
        if( ien(jp).eq.6 ) then
          zwght = (1-impmc_model)*wght(jp)+impmc_model*wstp(jp) ! VV_5.2
          if( is(jp).eq.0 ) then
            smstkn(jgi) = smstkn(jgi) + zwght
          else
            smstki(jgi) = smstki(jgi) + zwght
          endif
          goto 120
        endif
!
 120    continue
      enddo
!
!::common variables  wvacn, wvaci
      wvacn(0,0) = smgt0n  ! igtnp(ip) = 0
      wvaci(0,0) = smgt0i
      jgo3 = nogt + 3  ! in Vc(cal)
      jgo4 = nogt + 4  ! in Vc(stk)
      do jgi = 1, nogt
        wvacn(jgi,jgo3) = smcaln(jgi)
        wvaci(jgi,jgo3) = smcali(jgi)
        wvacn(jgi,jgo4) = smstkn(jgi)
        wvaci(jgi,jgo4) = smstki(jgi)
      enddo
!
!::sumup
      call setvwork4( 2, wwork4 )
      twork4(1:ns4) = wwork4(1:ns4)
!
!::[MPI_Reduce in imwcon]  cimcom_12c  (wvacn,wvaci) 12/05/25
      call MPI_Reduce( wwork4, vwork4, ns4, MPI_REAL8, MPI_SUM,
     >                 lmspe, mywld, ierr )
      if( lmype.ne.lmspe ) goto 200
!
      call setvwork4( 1, vwork4 )
!
!::vac-conductance (see imwrfn)
!   wvacn(1,0) = g1-in, g1/g2-ou(1,1/2), a(nogt+1) V(+2) in imwrfn
!    cal in Vac  inV(+3)   stk in Vac  stk(+4) in this routine
!
      jgo  = nogt + 3
      jgo4 = nogt + 4
      do jgi = 1, nogt
        zwern(jgi) = wvacn(jgi,0)
        zweri(jgi) = wvaci(jgi,0)
        do jex = 1, jgo4
          zwern(jgi) = zwern(jgi) - wvacn(jgi,jex)
          zweri(jgi) = zweri(jgi) - wvaci(jgi,jex)
        enddo
      enddo
      zwghn = wvacn(0,0) ! ghost : NO found gate
      zwghi = wvacn(0,0) ! ghost : NO found gate
!
!::out flux
      zwton(0:jgo) = 0.0d0
      zwtoi(0:jgo) = 0.0d0
      do jex = 0, jgo
        do jgi = 1, nogt
          zwton(jex) = zwton(jex) + wvacn(jgi,jex)
          zwtoi(jex) = zwtoi(jex) + wvaci(jgi,jex)
        enddo
      enddo
!
!::total out flux
      zwoun = 0.0d0
      zwoui = 0.0d0
      do jex = 1, jgo
        zwoun = zwoun + zwton(jex)
        zwoui = zwoui + zwtoi(jex)
      enddo
!
!::conductance
      jpm = nogt + 1     ! pump
      zfout = 0.0d0
      do jex = 1, jpm
        zfout = zfout + zwton(jex)
      enddo
!
      zcndv(0:nogt) = 0.0d0
      if( zwton(jpm).gt.0.0d0 ) then
        zcndv(0) = zfout/zwton(jpm)
        do jex = 1, nogt
          zcndv(jex) = zwton(jex)/zwton(jpm)
        enddo
      endif
!
      if( lprf.eq.1 ) then
        write(n6,'(/2x,"*** wvacn *** : neutral flux ")')
        write(n6,'(4x,"in",2x,"g0",10x,"g1",10x,"g2",10x,"a1",10x,"V1",
     >    10x,"inV",9x,"Err",9x,"Egt",9x,"stk")')
        do jgi = 1, nogt
          write(n6,'(4x,i2,1p10e12.4)') jgi
     >     ,(wvacn(jgi,i), i=0,6)
     >     ,zwern(jgi)
     >     ,zwghn, wvacn(jgi,jgo4)
        enddo
        write(n6,'(4x,a2,1p7e12.4)') " T"
     >     ,(zwton(i), i=0,ndgt), zwoun
!
        if( zwtoi(0).gt.0 ) then
          write(n6,'(/2x,"*** wvaci *** : ion flux")')
          write(n6,'(4x,"in",2x,"g0",10x,"g1",10x,"g2",10x,"a1",10x,
     >      "V1",10x,"inV",9x,"Err",9x,"Egt",9x,"stk")')
          do jgi = 1, nogt
            write(n6,'(4x,i2,1p10e12.4)') jgi
     >       ,(wvaci(jgi,i), i=0,6)
     >       ,zweri(jgi)
     >       ,zwghi, wvaci(jgi,jgo4)
          enddo
          write(n6,'(4x,a2,1p7e12.4)') " T"
     >     ,(zwtoi(i),i=0,5), zwoui
        endif
      endif
!
!::conductance  wvacn(jgi,jex) => bkflux(jex)
!  jgi = 1, 2(nogt), jex = 0(in), 1(bk1), 2(bk2), 3(pmp&exV)
!
      if( ( impmc_model == 0 .and. lgstk == 0 )
     >    .or.  impmc_model == 1 ) then
        do jex = 0, nogt
          zsum = 0.0d0
          do jgi = 1, nogt
            zsum = zsum + wvacn(jgi,jex)
          enddo
          bkflux(jex) = zsum
        enddo
!
        jab = nogt + 1           ! ab: absorption (pmp/Vex)
        zsum = 0.0d0
        do jex = nogt+1, nogt+2  ! (pmp & Vex )
          do jgi = 1, nogt
            zsum = zsum + wvacn(jgi,jex)
          enddo
        enddo
        bkflux(jab) = zsum
      ! bugescape for zero division
        if(bkflux(jab).ne.0.0_8) then
          do jex = 0, nogt
            bkfrac(jex) = bkflux(jex)/bkflux(jab)
          enddo
        else
          bkfrac = 0.0_8
          write(n6,'(4x,"imcndv: Warning, sum of wvacn is 0")')
        endif
!
      endif
!
!::set back after sumup
      call setvwork4( 1, twork4 )
!
 200  continue
      if( lgstk.eq.1 ) return
!
      call tbfind( nddt, typ_dnam, ityp, 'CIMCOM_12D' )
      itag = typ_itag(ityp)
      call MPI_Bcast( bkfrac, 1, itag, lmspe, mywld, ierr )
!
      if( lstdy.ge.1 .and. lmype.eq.lmspe .and. lhst.eq.1 ) then
        ih60 = 380000 + 2000 + mype
        open(unit=ih60,file="Hist_imcndv.txt")
        write(ih60,'(2x,"spty",7x,"lpcm",3x,"itimZ",2x,
     >    "timeZ",7x,"tout",8x,"cndI",8x,"cndg1",7x,"cndg2",7x,
     >    "flxIn",7x,"flxg1",7x,"flxg2",7x,"flxpm",7x,"flxVs",7x,
     >    "flxcal",7x,"flxou")')
        zero = 1.0d-30
        write(ih60,'(2x,a4,2x,i10,i7,1p2e12.4,1p15e12.4)')
     >    sptyc, lpcm, itimZ, timeZ, tout,
     >    dmax1(zcndv(0),zero),        ! in flux normalized by pump-flux
     >    dmax1(zcndv(1:nogt),zero),   ! conductance
     >    dmax1(zwton(0),zero),        ! in-flux
     >   (dmax1(zwton(jex),zero),jex=1,nogt),   ! out-gate
     >    dmax1(zwton(nogt+1),zero),   ! pump
     >    dmax1(zwton(nogt+2),zero),   ! abs Vessel
     >    dmax1(zwton(nogt+3),zero),   ! cal in Vac
     >    dmax1(zwoun,zero)            ! total for check
!
        write(n6,'(1x,a,2x,a4,i8,1pe12.4,1pe12.4,2x,1p20e12.4)')
     >    "imcndv",sptyc, itim, htout, (bkfrac(jex),jex=0,nogt),
     >    (bkflux(jex),jex=0,jab)
      endif
!
      return
      end
