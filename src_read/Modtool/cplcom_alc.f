! added dynamic allocation of arrays by kamata 2022/05/29
! allocate variables in cplcom
      subroutine cplcom_alc( kk )
      use cplcom, only : alnr, anamp, anapv, anasl, anemp, anepv, anesl
     >    , animp, anipv, anisl, atemp, atepv, atesl, atimp, atipv
     >    , atisl, cc, daan, dacfE, dacfW, dacl, dd, dfl1a, dfl2a, dfl3
     >    , dfl4, dmid, dq1a, dq2a, dq3, dq4, dtrateq, ee, ehcfE, ehcfW
     >    , etan, etcfE, etcfW, etcl, feqp, ff, fl1a, fl2a, fl3, fl4
     >    , fl4m, fl4p, flmet, flmxe, flmxi, flps, flvl, fvpmp, gcse
     >    , gcsi, gg, hea, heat_flux_para_by_kappa_para
     >    , heat_flux_para_elec, hia, ibcyl, iupwd, kappa_helander_test
     >    , q1a, q2a, q3, q4, qb1a, qb2a, qb3, qb4, r1ae, r1ai, r1av
     >    , r1aw, r2ai, r2av, r2aw, r3i, r4e, ss, ssnc, ssnv, sspc, sspv
     >    , swec, swev, swic, swiv, tN0, tNg, tT0, tTg, tV0, temndp, vcs
     >    , vdda, vdet, vdxe, vdxi, vea, vna, vnag, vne, vnezef, vni
     >    , vte, vteg, vti, vtig, vva, vvag, vve, vzf, wdnz, wq1a, wq2a
     >    , wq3, wq4, wrad, xean, xecfE, xecfW, xecl, xfdfa, xian, xicfE
     >    , xicfW, xicl, xicl_test, xmid, xste, xsti, xsva, xwtam, xwtap
     >    , ymid
      use csize,  only : ndeq, ndgs, ndmc, ndsp, ndwp, ndx, ndxy, ndy
      implicit none
! arguments
      integer, intent(in) :: kk
! kk : flag = 1: allocate,    = 2: deallocate

! local variables
      integer    istat
      character  cmsg*80

      select case( kk )
      case( 1 )
! allocate
        if( .not. allocated( alnr ) ) then
          allocate( alnr(ndxy), anamp(ndy,ndsp), anapv(ndy,ndsp)
     >      , anasl(ndy,ndsp), anemp(ndy), anepv(ndy), anesl(ndy)
     >      , animp(ndy), anipv(ndy), anisl(ndy), atemp(ndy), atepv(ndy)
     >      , atesl(ndy), atimp(ndy), atipv(ndy), atisl(ndy)
     >      , cc(ndeq,ndeq,ndxy), daan(ndxy,ndsp), dacfE(ndxy,ndsp)
     >      , dacfW(ndxy,ndsp), dacl(ndxy,ndsp), dd(ndeq,ndeq,ndxy)
     >      , dfl1a(ndxy,ndsp), dfl2a(ndxy,ndsp), dfl3(ndxy), dfl4(ndxy)
     >      , dmid(ndy), dq1a(ndx,ndy,ndsp), dq2a(ndx,ndy,ndsp)
     >      , dq3(ndx,ndy), dq4(ndx,ndy), dtrateq(ndx,ndy,ndeq)
     >      , ee(ndeq,ndeq,ndxy), ehcfE(ndxy,ndsp), ehcfW(ndxy,ndsp)
     >      , etan(ndxy,ndsp), etcfE(ndxy,ndsp), etcfW(ndxy,ndsp)
     >      , etcl(ndxy,ndsp), feqp(ndxy), ff(ndeq,ndxy), fl1a(ndx,ndsp)
     >      , fl2a(ndx,ndsp), fl3(ndx), fl4(ndx), fl4m(ndx), fl4p(ndx)
     >      , flmet(ndx,ndy), flmxe(ndx,ndy), flmxi(ndx,ndy)
     >      , flps(ndwp,ndsp), flvl(ndx,ndy,ndsp), fvpmp(ndy), gcse(ndy)
     >      , gcsi(ndy), gg(ndeq,ndxy), hea(ndxy,ndsp)
     >      , heat_flux_para_by_kappa_para(ndx,ndy)
     >      , heat_flux_para_elec(ndx,ndy), hia(ndxy,ndsp), ibcyl(ndy)
     >      , iupwd(ndx,ndy), kappa_helander_test(ndx,ndy)
     >      , q1a(ndx,ndy,ndsp), q2a(ndx,ndy,ndsp), q3(ndx,ndy)
     >      , q4(ndx,ndy), qb1a(ndxy,ndsp,4), qb2a(ndxy,ndsp,4)
     >      , qb3(ndxy,4), qb4(ndxy,4), r1ae(ndxy,ndsp), r1ai(ndxy,ndsp)
     >      , r1av(ndxy,ndsp), r1aw(ndxy,ndsp), r2ai(ndxy,ndsp)
     >      , r2av(ndxy,ndsp), r2aw(ndxy,ndsp), r3i(ndxy), r4e(ndxy)
     >      , ss(ndeq,ndeq,ndxy), ssnc(ndx,ndy,ndsp), ssnv(ndx,ndy,ndsp)
     >      , sspc(ndx,ndy,ndsp), sspv(ndx,ndy,ndsp), swec(ndx,ndy)
     >      , swev(ndx,ndy), swic(ndx,ndy), swiv(ndx,ndy)
     >      , tN0(0:ndmc,ndgs), tNg(0:ndmc,2), tT0(0:ndmc,ndgs)
     >      , tTg(0:ndmc,2), tV0(0:ndmc,ndgs), temndp(ndy), vcs(ndx,ndy)
     >      , vdda(ndx,ndy,ndsp), vdet(ndx,ndy,ndsp), vdxe(ndx,ndy)
     >      , vdxi(ndx,ndy), vea(ndx,ndy,ndsp), vna(ndx,ndy,ndsp)
     >      , vnag(ndx,ndy,ndsp), vne(ndx,ndy), vnezef(ndx,ndy)
     >      , vni(ndx,ndy), vte(ndx,ndy), vteg(ndx,ndy), vti(ndx,ndy)
     >      , vtig(ndx,ndy), vva(ndx,ndy,ndsp), vvag(ndx,ndy,ndsp)
     >      , vve(ndx,ndy), vzf(ndx,ndy), wdnz(ndx,ndy)
     >      , wq1a(ndx,ndy,ndsp), wq2a(ndx,ndy,ndsp), wq3(ndx,ndy)
     >      , wq4(ndx,ndy), wrad(ndx,ndy), xean(ndxy), xecfE(ndxy)
     >      , xecfW(ndxy), xecl(ndxy), xfdfa(ndx,ndy,ndsp), xian(ndxy)
     >      , xicfE(ndxy), xicfW(ndxy), xicl(ndxy), xicl_test(ndx,ndy)
     >      , xmid(ndy), xste(ndxy,ndsp), xsti(ndxy,ndsp)
     >      , xsva(ndxy,ndsp), xwtam(ndx,ndy,ndsp), xwtap(ndx,ndy,ndsp)
     >      , ymid(ndy)
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'alnr allocate error in cplcom_alc, istat = ', istat
            call wexit( 'cplcom', trim( cmsg ) )
          endif

! initial set
          alnr(1:ndxy) = 0.0_8 
          anamp(1:ndy,1:ndsp) = 0.0_8
          anapv(1:ndy,1:ndsp) = 0.0_8
          anasl(1:ndy,1:ndsp) = 0.0_8
          anemp(1:ndy) = 0.0_8
          anepv(1:ndy) = 0.0_8
          anesl(1:ndy) = 0.0_8
          animp(1:ndy) = 0.0_8
          anipv(1:ndy) = 0.0_8
          anisl(1:ndy) = 0.0_8
          atemp(1:ndy) = 0.0_8
          atepv(1:ndy) = 0.0_8
          atesl(1:ndy) = 0.0_8
          atimp(1:ndy) = 0.0_8
          atipv(1:ndy) = 0.0_8
          atisl(1:ndy) = 0.0_8
          cc(1:ndeq,1:ndeq,1:ndxy) = 0.0_8
          daan(1:ndxy,1:ndsp) = 0.0_8
          dacfE(1:ndxy,1:ndsp) = 0.0_8
          dacfW(1:ndxy,1:ndsp) = 0.0_8
          dacl(1:ndxy,1:ndsp) = 0.0_8
          dd(1:ndeq,1:ndeq,1:ndxy) = 0.0_8
          dfl1a(1:ndxy,1:ndsp) = 0.0_8
          dfl2a(1:ndxy,1:ndsp) = 0.0_8
          dfl3(1:ndxy) = 0.0_8
          dfl4(1:ndxy) = 0.0_8
          dmid(1:ndy) = 0.0_8
          dq1a(1:ndx,1:ndy,1:ndsp) = 0.0_8
          dq2a(1:ndx,1:ndy,1:ndsp) = 0.0_8
          dq3(1:ndx,1:ndy) = 0.0_8
          dq4(1:ndx,1:ndy) = 0.0_8
!         dtrat(1:ndx,1:ndy) = 0.0_8
          dtrateq(1:ndx,1:ndy,1:ndeq) = 0.0_8
          ee(1:ndeq,1:ndeq,1:ndxy) = 0.0_8
          ehcfE(1:ndxy,1:ndsp) = 0.0_8
          ehcfW(1:ndxy,1:ndsp) = 0.0_8
          etan(1:ndxy,1:ndsp) = 0.0_8
          etcfE(1:ndxy,1:ndsp) = 0.0_8
          etcfW(1:ndxy,1:ndsp) = 0.0_8
!         etch(1:ndxy,1:ndsp) = 0.0_8
          etcl(1:ndxy,1:ndsp) = 0.0_8
          feqp(1:ndxy) = 0.0_8
          ff(1:ndeq,1:ndxy) = 0.0_8
          fl1a(1:ndx,1:ndsp) = 0.0_8
!         fl1am(1:ndx,1:ndsp) = 0.0_8
!         fl1ap(1:ndx,1:ndsp) = 0.0_8
          fl2a(1:ndx,1:ndsp) = 0.0_8
!         fl2am(1:ndx,1:ndsp) = 0.0_8
!         fl2ap(1:ndx,1:ndsp) = 0.0_8
          fl3(1:ndx) = 0.0_8
!         fl3m(1:ndx) = 0.0_8
!         fl3p(1:ndx) = 0.0_8
          fl4(1:ndx) = 0.0_8
          fl4m(1:ndx) = 0.0_8
          fl4p(1:ndx) = 0.0_8
          flmet(1:ndx,1:ndy) = 0.0_8
          flmxe(1:ndx,1:ndy) = 0.0_8
          flmxi(1:ndx,1:ndy) = 0.0_8
          flps(1:ndwp,1:ndsp) = 0.0_8
          flvl(1:ndx,1:ndy,1:ndsp) = 0.0_8
          fvpmp(1:ndy) = 0.0_8
          gcse(1:ndy) = 0.0_8
          gcsi(1:ndy) = 0.0_8
          gg(1:ndeq,1:ndxy) = 0.0_8
          hea(1:ndxy,1:ndsp) = 0.0_8
          heat_flux_para_by_kappa_para(1:ndx,1:ndy) = 0.0_8
          heat_flux_para_elec(1:ndx,1:ndy) = 0.0_8
          hia(1:ndxy,1:ndsp) = 0.0_8
!         ibcdpi(1:ndy) = 0.0_8
!         ibcdpo(1:ndy) = 0.0_8
          ibcyl(1:ndy) = 0.0_8
          iupwd(1:ndx,1:ndy) = 0.0_8
          kappa_helander_test(1:ndx,1:ndy) = 0.0_8
          q1a(1:ndx,1:ndy,1:ndsp) = 0.0_8
          q2a(1:ndx,1:ndy,1:ndsp) = 0.0_8
          q3(1:ndx,1:ndy) = 0.0_8
          q4(1:ndx,1:ndy) = 0.0_8
          qb1a(1:ndxy,1:ndsp,1:4) = 0.0_8
          qb2a(1:ndxy,1:ndsp,1:4) = 0.0_8
          qb3(1:ndxy,1:4) = 0.0_8
          qb4(1:ndxy,1:4) = 0.0_8
          r1ae(1:ndxy,1:ndsp) = 0.0_8
          r1ai(1:ndxy,1:ndsp) = 0.0_8
          r1av(1:ndxy,1:ndsp) = 0.0_8
          r1aw(1:ndxy,1:ndsp) = 0.0_8
!         r2ae(1:ndxy,1:ndsp) = 0.0_8
          r2ai(1:ndxy,1:ndsp) = 0.0_8
          r2av(1:ndxy,1:ndsp) = 0.0_8
          r2aw(1:ndxy,1:ndsp) = 0.0_8
          r3i(1:ndxy) = 0.0_8
          r4e(1:ndxy) = 0.0_8
!         seqp(1:ndxy) = 0.0_8
!         serec(1:ndx,1:ndy) = 0.0_8
!         sfrec(1:ndx,1:ndy) = 0.0_8
          ss(1:ndeq,1:ndeq,1:ndxy) = 0.0_8
          ssnc(1:ndx,1:ndy,1:ndsp) = 0.0_8
          ssnv(1:ndx,1:ndy,1:ndsp) = 0.0_8
          sspc(1:ndx,1:ndy,1:ndsp) = 0.0_8
          sspv(1:ndx,1:ndy,1:ndsp) = 0.0_8
          swec(1:ndx,1:ndy) = 0.0_8
          swev(1:ndx,1:ndy) = 0.0_8
          swic(1:ndx,1:ndy) = 0.0_8
          swiv(1:ndx,1:ndy) = 0.0_8
          tN0(0:ndmc,1:ndgs) = 0.0_8
          tNg(0:ndmc,1:2) = 0.0_8
          tT0(0:ndmc,1:ndgs) = 0.0_8
          tTg(0:ndmc,1:2) = 0.0_8
          tV0(0:ndmc,1:ndgs) = 0.0_8
          temndp(1:ndy) = 0.0_8
          vcs(1:ndx,1:ndy) = 0.0_8
          vdda(1:ndx,1:ndy,1:ndsp) = 0.0_8
          vdet(1:ndx,1:ndy,1:ndsp) = 0.0_8
          vdxe(1:ndx,1:ndy) = 0.0_8
          vdxi(1:ndx,1:ndy) = 0.0_8
          vea(1:ndx,1:ndy,1:ndsp) = 0.0_8
          vna(1:ndx,1:ndy,1:ndsp) = 0.0_8
          vnag(1:ndx,1:ndy,1:ndsp) = 0.0_8
          vne(1:ndx,1:ndy) = 0.0_8
          vnezef(1:ndx,1:ndy) = 0.0_8
          vni(1:ndx,1:ndy) = 0.0_8
          vte(1:ndx,1:ndy) = 0.0_8
          vteg(1:ndx,1:ndy) = 0.0_8
          vti(1:ndx,1:ndy) = 0.0_8
          vtig(1:ndx,1:ndy) = 0.0_8
          vva(1:ndx,1:ndy,1:ndsp) = 0.0_8
          vvag(1:ndx,1:ndy,1:ndsp) = 0.0_8
          vve(1:ndx,1:ndy) = 0.0_8
          vzf(1:ndx,1:ndy) = 0.0_8
          wdnz(1:ndx,1:ndy) = 0.0_8
!         wlos(1:ndx,1:ndy) = 0.0_8
          wq1a(1:ndx,1:ndy,1:ndsp) = 0.0_8
          wq2a(1:ndx,1:ndy,1:ndsp) = 0.0_8
          wq3(1:ndx,1:ndy) = 0.0_8
          wq4(1:ndx,1:ndy) = 0.0_8
          wrad(1:ndx,1:ndy) = 0.0_8
          xean(1:ndxy) = 0.0_8
          xecfE(1:ndxy) = 0.0_8
          xecfW(1:ndxy) = 0.0_8
          xecl(1:ndxy) = 0.0_8
          xfdfa(1:ndx,1:ndy,1:ndsp) = 0.0_8
          xian(1:ndxy) = 0.0_8
          xicfE(1:ndxy) = 0.0_8
          xicfW(1:ndxy) = 0.0_8
          xicl(1:ndxy) = 0.0_8
          xicl_test(1:ndx,1:ndy) = 0.0_8
          xmid(1:ndy) = 0.0_8
          xste(1:ndxy,1:ndsp) = 0.0_8
          xsti(1:ndxy,1:ndsp) = 0.0_8
          xsva(1:ndxy,1:ndsp) = 0.0_8
          xwtam(1:ndx,1:ndy,1:ndsp) = 0.0_8
          xwtap(1:ndx,1:ndy,1:ndsp) = 0.0_8
          ymid(1:ndy) = 0.0_8
        endif
      case( 2 )
! deallocate
        if( allocated( alnr ) ) then
          deallocate( alnr, anamp, anapv, anasl, anemp, anepv, anesl
     >      , animp, anipv, anisl, atemp, atepv, atesl, atimp, atipv
     >      , atisl, cc, daan, dacfE, dacfW, dacl, dd, dfl1a, dfl2a
     >      , dfl3, dfl4, dmid, dq1a, dq2a, dq3, dq4, dtrateq, ee, ehcfE
     >      , ehcfW, etan, etcfE, etcfW, etcl, feqp, ff, fl1a, fl2a, fl3
     >      , fl4, fl4m, fl4p, flmet, flmxe, flmxi, flps, flvl, fvpmp
     >      , gcse, gcsi, gg, hea, heat_flux_para_by_kappa_para
     >      , heat_flux_para_elec, hia, ibcyl, iupwd
     >      , kappa_helander_test, q1a, q2a, q3, q4, qb1a, qb2a, qb3
     >      , qb4, r1ae, r1ai, r1av, r1aw, r2ai, r2av, r2aw, r3i, r4e
     >      , ss, ssnc, ssnv, sspc, sspv, swec, swev, swic, swiv, tN0
     >      , tNg, tT0, tTg, tV0, temndp, vcs, vdda, vdet, vdxe, vdxi
     >      , vea, vna, vnag, vne, vnezef, vni, vte, vteg, vti, vtig
     >      , vva, vvag, vve, vzf, wdnz, wq1a, wq2a, wq3, wq4, wrad
     >      , xean, xecfE, xecfW, xecl, xfdfa, xian, xicfE, xicfW, xicl
     >      , xicl_test, xmid, xste, xsti, xsva, xwtam, xwtap, ymid
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'alnr deallocate error in cplcom_alc, istat = ', istat
            call wexit( 'cplcom', trim( cmsg ) )
          endif
        endif
      end select

      return
      end
