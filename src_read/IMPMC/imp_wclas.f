!***********************************************************************
      subroutine imp_wclas(grpp,grpw,grpn,regp,regw,regn)
!***********************************************************************
!)
!)    imp_wclas  imp(urity) w(eight) classified
!)
!)   myptl = sum(1)                : particle
!)   mywgt = sum(wght(ip))         : weight
!)   mypcn = sum(pemt(ip)*wght(ip) : (particle) number
!)
!)   grp  : (emt/cal/sgt/sce/min/err/exh/pmp/nul/run/tst/ten/tch/tav
!)   reg  : (odv/sol/idv/opv/ipv/edg/man/vac)
!)
!)  igo : treatment of particle
!)    1: enter  2: survive  3: terminate
!)    4: exhaust (continue to trace)
!)
!)     ityp           igo                        ien(ip)
!)        1  igemt   :  1  emitted                -  ==> imp_lnch
!)        2  igcal   :  2  reach tmend            0
!)        3  igsgt   :  3  stick at gate/wall     6
!)        4  igsce   :  3  stick at core edge     8
!)        5  igmin   :  3  wght < wmin            7
!)        6  igerr   :  3  error                  9
!)        7  igexh   :  4  exhaust  at wall       -  ==> imwrfi/imwrfn
!)        8  igpmp   :  4  pump                   -  ==> imwrfi/imwrfn
!)
!)        9  ignul   ????  (idea of generation)
!)       10  igrun   :  2  loop in imtrcI/N      -1  Serious error
!)
!)       11  igtst   tst: at start time
!)       12  igten   ten: at end time
!)       13  igtbl   tbl: check conservation
!)       14  igtav   tav: Integ(<Nz>*dvl)
!)
!::conservation check
!)     Wsum = Wcal + Wsgt + Wsce + Wmin + Werr + ... [0,10]
!)     Wtst = Wsrv   from difinition  srv: suvive particle at Tst
!)     Wten = Wsum   from difinition  sum: sumup
!)     Wtbl = Wtst + Wemt - Wexh      bal: balance
!)
!)  O  Wten = Wbal
!)  O  Wtav ~ (Wtst + Wen)/2   (~ :include MC process)
!)
!)  how to scatter  10:region for survive particle (igo = 2 )
!)     1                 ==> regp(irg)
!)     weit(ip)          ==> regw(irg)
!)     pemt(ip)*weit(ip) ==> regn(irg)
!)
!)Stedy grpw emt  = cal  + exh + stk +  min + err + ...
!)
!)IMPTC grpw srv + emt = sum + wxh  sum = ten
!)
!)    regw(irg)   cal = odv + sol + idv + opv + ipv + ceg + cre + vac
!)    Error in case of run > 0.0
!)
!)::refer common variables
!)      integer :: npmax, ien(1), ihwl(1), il(1), ir(1), mrgn(1)
!)      character(1) :: chwl(1)
!)      real(8) :: wght(1)
!)
!)    difference of krfl = 5 and 6
!)      krfl = 5  lwstk(ih) = 1                     ien(ip) = 6
!)      krfl = 6  lwstk(ih) = 0 and wght < wrfmin   ien(ip) = 7
!)
!---------------------------------------------------------------------
!)::important
!)     ien(ip) = -1(run) 0(tmend) 1(puf) 2(rec) 3(ion) 4(wln) 5(wli)
!)               6(stkGT) 7(absMN) 8(stkCE) 9(err)
!)         ==>  10(excl.)  ! ien /= 0 of already emitted particle
!)
!)   terminate trace
!)      in case of  ien = 6(stkGT) 7(absMN) 8(stkCE) 9(err)
!)
!)   OLd version : Do not set wght(ip) = 0         ! Important
!)   New version : when ien = 6-10  wstp= wght wght=0  ! VV_5.2
!)   2 type     wght = 0.0  death  wght > 0  life
!)
!----------------------------------------------------------------------
      use cimcom, only : ien, igcal, igemt, igerr, igexh, igmin, igpmp
     >    , igrun, igsce, igsgt, igtav, igtbl, igten, igtst, ir, mypcn
     >    , myptl, mywgt, npmax, ntags, pemt, stimz, wght, wstp
      use cimden, only : nsrc_spt
      use cntcom, only : mrgn
      use cunit,  only : n6
      implicit none
!
!::argument  grp: group  reg: region
      real(8), intent(out) :: grpp(15), grpw(15), grpn(15)
      real(8), intent(out) :: regp(10), regw(10), regn(10)

!::local variables
      integer :: ip, icnd, ic, irg, i
      real(8) :: zwgt
      real(8) :: tot_regp, tot_regw, tot_regn
      character :: cmsg*(80)
      integer :: ignum(13)

      ignum(1:13) =
     >   (/igtst,igten,igtbl,igtav,
     >     igemt,igexh,igpmp,igcal,igsgt,igsce,igmin,igerr,igrun/)
!)     ignul : no effec
!)
!----------------------------------------------------------------------
!::classify test particles  for conservation
!----------------------------------------------------------------------
!::classfy group and region sorvive prtciles ien = 0
!::p: particles  w: weight  n: p-number
      grpp(1:15) = 0.0d0
      grpw(1:15) = 0.0d0
      grpn(1:15) = 0.0d0
      regp(1:10) = 0.0d0
      regw(1:10) = 0.0d0
      regn(1:10) = 0.0d0
!
!::sample partcile 1 PE
      grpp(igemt) = myptl(igemt)
      grpp(igpmp) = myptl(igpmp)
      grpp(igexh) = myptl(igexh)
      grpw(igemt) = mywgt(igemt)
      grpw(igpmp) = mywgt(igpmp)
      grpw(igexh) = mywgt(igexh)
      grpn(igemt) = mypcn(igemt)
      grpn(igpmp) = mypcn(igpmp)
      grpn(igexh) = mypcn(igexh)

      grpp(igtst) = myptl(igtst)
      grpp(igten) = myptl(igten)
      grpp(igtav) = myptl(igtav)
      grpw(igtst) = mywgt(igtst)
      grpw(igten) = mywgt(igten)
      grpw(igtav) = mywgt(igtav)
      grpn(igtst) = mypcn(igtst)
      grpn(igten) = mypcn(igten)
      grpn(igtav) = mypcn(igtav)

      do ip = 1, npmax
        if( ien(ip) == 10 ) cycle  ! important
        if( ntags(ip).ne.nsrc_spt) cycle !important
        icnd = ien(ip)

!::ien = -1, 0     life particle    cal   wght = 0.92 wstp = 0.0
!::ien =  6,7,8,9  death particle   stop  wght = 0.0  wstp = 0.92

        select case(icnd)
!
!::run
        case(-1)  ! life partcile
          zwgt = wght(ip)
          grpp(igrun) = grpp(igrun) + 1
          grpw(igrun) = grpw(igrun) + zwgt
          grpn(igrun) = grpn(igrun) + pemt(ip)*zwgt
!::cal
        case(0)   ! life particle
          zwgt = wght(ip)
          grpp(igcal) = grpp(igcal) + 1
          grpw(igcal) = grpw(igcal) + zwgt
          grpn(igcal) = grpn(igcal) + pemt(ip)*zwgt
!::stkGT
        case(6)  ! death particle
          zwgt = wstp(ip)
          grpp(igsgt) = grpp(igsgt) + 1
          grpw(igsgt) = grpw(igsgt) + zwgt
          grpn(igsgt) = grpn(igsgt) + pemt(ip)*zwgt
!::stkCE
        case(8)  ! death particle
          zwgt = wstp(ip)
          grpp(igsce) = grpp(igsce) + 1
          grpw(igsce) = grpw(igsce) + zwgt
          grpn(igsce) = grpn(igsce) + pemt(ip)*zwgt
!::wmin
        case(7)  ! death particle (not error)
          zwgt = wstp(ip)
          grpp(igmin) = grpp(igmin) + 1
          grpw(igmin) = grpw(igmin) + zwgt
          grpn(igmin) = grpn(igmin) + pemt(ip)*zwgt
!::err
        case(9)  ! death particle
          zwgt = wstp(ip)
          grpp(igerr) = grpp(igerr) + 1
          grpw(igerr) = grpw(igerr) + zwgt
          grpn(igerr) = grpn(igerr) + pemt(ip)*zwgt

!::nul  cal 0 ==> nul 10 for CErs
!       No use

!::others
        case default
           write(cmsg,'("ien(ip) =",i3,2x,a)')
     >       icnd,"0(cal),6(stkGT),8(stkCE),7(wmin),9(err)"
           call wexit("imwclas",trim(cmsg))
        end select

        if( ien(ip) == 0 ) then  ! life particle
          ic  = ir(ip)
          irg = mrgn(ic)
          regp(irg) = regp(irg) + 1.0d0
          regw(irg) = regw(irg) + wght(ip)
          regn(irg) = regn(irg) + pemt(ip)*wght(ip)
        endif
      enddo

!::total region
      tot_regp = 0.0d0
      tot_regw = 0.0d0
      tot_regn = 0.0d0
      do irg = 1, 8
        tot_regp = tot_regp + regp(irg)
        tot_regw = tot_regw + regw(irg)
        tot_regn = tot_regn + regn(irg)
      enddo

      grpp(igtbl) = grpp(igtst) + grpp(igemt) - grpp(igexh)
      grpw(igtbl) = grpw(igtst) + grpw(igemt) - grpw(igexh)
      grpn(igtbl) = grpn(igtst) + grpn(igemt) - grpn(igexh)

!::debug write
      write(n6,'(2x,a,2x,a,2x,"stimz =",1pe11.4)') "imp_wclas",
     >  "tst/ten/tbl/tav  grp(emt/exh/pmp)"//
     >  "  grp(cal/sgt/sce/min/err/run)", stimz

      write(n6,602) (grpp(ignum(i)),i=1,4), "grpp",
     >  (grpp(ignum(i)),i=5,13)
      write(n6,602) (grpw(ignum(i)),i=1,4), "grpw",
     >  (grpw(ignum(i)),i=5,13)
      write(n6,602) (grpn(ignum(i)),i=1,4), "grpn",
     >  (grpn(ignum(i)),i=5,13)
      write(n6,'(2x,1p2e12.4,2x,a)')  (grpn(igtst)+grpn(igten))/2.0d0,
     >  grpn(igtav),  "(tst+ten)/2~tav[Nz]"

      write(n6,604) grpp(igcal), tot_regp, "regp", regp(1:8)
      write(n6,604) grpw(igcal), tot_regw, "regw", regw(1:8)
      write(n6,604) grpn(igcal), tot_regn, "regn", regn(1:8)

 602  format(2x,1p4e12.4,3x,a,5x,1p9e12.4)
 604  format(2x,1p2e12.4,3x,a,5x,1p8e12.4)

      return
      end
