! added dynamic allocation of arrays by kamata 2022/05/29
! allocate variables in cntcom
      subroutine cntcom_alc( kk )
      use cntcom, only : bvx, bvy, bvz, csbr, csps, cswl, den0, deng
     >    , e0br, e0ps, e0wl, eibr, eips, eng0, engg, flbr, icbr
     >    , icps, icwl, iebr, igtw, ihwl, ikbr, ikps, ikwl, imps, imwl
     >    , inwl, ipfiw, ipfmp, iplx, iply, ipmp, ipps, ipwl, isbr, isxp
     >    , isxw, isyp, isyw, iwbr, iwpbr, ixyw, mcel, mclx, mcly, mcly2
     >    , mgrd, migx, migy, mknd, mrgn, mseg, mxjw, mxjw2, myit, next
     >    , nocl, nogd
     >    , nogdv, pfflx, prbr, rcps, rcwl, snbr, snps, snwl, tbde_ds
     >    , tbel_ds, tbel_ei, tbsg_ds, tbsg_ei, tywl, vion, vlp0, volm
     >    , xpnt, ypnt
     >    , dsgt, iwgt, ipgt, icgt, nfl0x, nfl0y, nfl0z
     >    , iwpbr_wall
      use csize, only : ndbr, ndgs, ndmc, ndmg, ndms, ndwp, ndx, ndy
     >    , nvxd, nvyd, ndgtp
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
        if( .not. allocated( bvx ) ) then
          allocate( bvx(0:ndmc), bvy(0:ndmc), bvz(0:ndmc), csbr(ndbr)
     >      , csps(ndwp), cswl(ndwp), den0(0:ndmc,ndgs), deng(0:ndmc,2)
     >      , e0br(ndbr), e0ps(ndwp), e0wl(ndwp), eibr(ndbr), eips(ndwp)
     >      , eng0(0:ndmc,ndgs), engg(0:ndmc,2), flbr(ndbr)
     >      , icbr(ndbr), icps(ndwp), icwl(ndwp), iebr(ndbr), igtw(ndwp)
     >      , ihwl(ndwp), ikbr(ndbr), ikps(ndwp), ikwl(ndwp), imps(ndwp)
     >      , imwl(ndwp), inwl(ndwp), ipfiw(ndwp), ipfmp(ndwp)
     >      , iplx(0:ndmc), iply(0:ndmc), ipmp(ndwp), ipps(ndwp)
     >      , ipwl(ndwp), isbr(ndbr), isxp(ndwp), isxw(ndwp), isyp(ndwp)
     >      , isyw(ndwp), iwbr(ndbr), iwpbr(ndbr), ixyw(ndwp)
     >      , mcel(0:ndx,0:ndy), mclx(ndx), mcly(ndy), mcly2(ndy)
     >      , mgrd(ndmc,ndms+1), migx(ndmc), migy(ndmc), mknd(ndmc,ndms)
     >      , mrgn(ndmc), mseg(ndmc), mxjw(ndx), mxjw2(ndx), myit(ndy)
     >      , next(ndmc,ndms)
     >      , nocl(0:ndx,0:ndy)
     >      , nogd(0:ndx,0:ndy), nogdv(0:nvxd,0:nvyd), pfflx(ndwp)
     >      , prbr(ndbr), rcps(ndwp), rcwl(ndwp), snbr(ndbr), snps(ndwp)
     >      , snwl(ndwp), tbde_ds(0:ndmc,6), tbel_ds(0:ndmc,6)
     >      , tbel_ei(0:ndmc), tbsg_ds(0:ndmc,6), tbsg_ei(0:ndmc)
     >      , tywl(ndwp), vion(0:ndmc,ndgs), vlp0(0:ndmc,ndgs)
     >      , volm(ndmc), xpnt(ndmg), ypnt(ndmg)
     >      , dsgt(ndgtp), iwgt(ndgtp), ipgt(ndgtp), icgt(ndgtp)
     >      , nfl0x(0:ndmc,ndgs), nfl0y(0:ndmc,ndgs), nfl0z(0:ndmc,ndgs)
     >      , iwpbr_wall(ndbr)
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'bvx allocate error in cntcom_alc, istat = ', istat
            call wexit( 'cntcom', trim( cmsg ) )
          endif

! initial set
          bvx(0:ndmc) = 0.0_8
          bvy(0:ndmc) = 0.0_8 
          bvz(0:ndmc) = 0.0_8 
          csbr(1:ndbr) = 0.0_8 
          csps(1:ndwp) = 0.0_8 
          cswl(1:ndwp) = 0.0_8 
          den0(0:ndmc,1:ndgs) = 0.0_8 
          deng(0:ndmc,1:2) = 0.0_8
          e0br(1:ndbr) = 0.0_8
          e0ps(1:ndwp) = 0.0_8
          e0wl(1:ndwp) = 0.0_8  
          eibr(1:ndbr) = 0.0_8 
          eips(1:ndwp) = 0.0_8 
          eng0(0:ndmc,1:ndgs) = 0.0_8
          engg(0:ndmc,1:2) = 0.0_8
          flbr(1:ndbr) = 0.0_8 
          icbr(1:ndbr) = 0
          icps(1:ndwp) = 0 
          icwl(1:ndwp) = 0
          iebr(1:ndbr) = 0
          igtw(1:ndwp) = 0
          ihwl(1:ndwp) = 0
          ikbr(1:ndbr) = 0
          ikps(1:ndwp) = 0
          ikwl(1:ndwp) = 0
          imps(1:ndwp) = 0
          imwl(1:ndwp) = 0
          inwl(1:ndwp) = 0
          ipfiw(1:ndwp) = 0
          ipfmp(1:ndwp) = 0
          iplx(0:ndmc) = 0
          iply(0:ndmc) = 0
          ipmp(1:ndwp) = 0
          ipps(1:ndwp) = 0
          ipwl(1:ndwp) = 0
          isbr(1:ndbr) = 0
          isxp(1:ndwp) = 0
          isxw(1:ndwp) = 0
          isyp(1:ndwp) = 0
          isyw(1:ndwp) = 0 
          iwbr(1:ndbr) = 0
          iwpbr(1:ndbr) = 0
          ixyw(1:ndwp) = 0
          mcel(0:ndx,0:ndy) = 0
          mclx(1:ndx) = 0
          mcly(1:ndy) = 0
          mcly2(1:ndy) = 0
          mgrd(1:ndmc,1:ndms+1) = 0
          migx(1:ndmc) = 0
          migy(1:ndmc) = 0
          mknd(1:ndmc,1:ndms) = 0
          mrgn(1:ndmc) = 0
          mseg(1:ndmc) = 0
          mxjw(1:ndx) = 0
          mxjw2(1:ndx) = 0
          myit(1:ndy) = 0
          next(1:ndmc,1:ndms) = 0
          nocl(0:ndx,0:ndy) = 0
          nogd(0:ndx,0:ndy) = 0 
          nogdv(0:nvxd,0:nvyd) = 0
          pfflx(1:ndwp) = 0.0_8
          prbr(1:ndbr) = 0.0_8
          rcps(1:ndwp) = 0.0_8
          rcwl(1:ndwp) = 0.0_8
          snbr(1:ndbr) = 0.0_8
          snps(1:ndwp) = 0.0_8
          snwl(1:ndwp) = 0.0_8
          tbde_ds(0:ndmc,1:6) = 0.0_8
          tbel_ds(0:ndmc,1:6) = 0.0_8
          tbel_ei(0:ndmc) = 0.0_8
          tbsg_ds(0:ndmc,1:6) = 0.0_8
          tbsg_ei(0:ndmc) = 0.0_8
          tywl(1:ndwp) = ' '
          vion(0:ndmc,1:ndgs) = 0.0_8
          vlp0(0:ndmc,1:ndgs) = 0.0_8
          volm(1:ndmc) = 0.0_8
          xpnt(1:ndmg) = 0.0_8
          ypnt(1:ndmg) = 0.0_8
          dsgt(1:ndgtp) = 0.0_8
          iwgt(1:ndgtp) = 0
          ipgt(1:ndgtp) = 0
          icgt(1:ndgtp) = 0
          nfl0x(0:ndmc,1:ndgs) = 0.0_8
          nfl0y(0:ndmc,1:ndgs) = 0.0_8
          nfl0z(0:ndmc,1:ndgs) = 0.0_8
          iwpbr_wall(1:ndbr) = 0
        endif
      case( 2 )
! deallocate
        if( allocated( bvx ) ) then
          deallocate( bvx, bvy, bvz, csbr, csps, cswl, den0, deng, e0br
     >      , e0ps, e0wl, eibr, eips, eng0, engg, flbr, icbr, icps
     >      , icwl, iebr, igtw, ihwl, ikbr, ikps, ikwl, imps, imwl, inwl
     >      , ipfiw, ipfmp, iplx, iply, ipmp, ipps, ipwl, isbr, isxp
     >      , isxw, isyp, isyw, iwbr, iwpbr, ixyw, mcel, mclx, mcly
     >      , mcly2, mgrd, migx, migy, mknd, mrgn, mseg, mxjw, mxjw2
     >      , myit, next, nocl
     >      , nogd, nogdv, pfflx, prbr, rcps, rcwl, snbr, snps
     >      , snwl, tbde_ds, tbel_ds, tbel_ei, tbsg_ds, tbsg_ei, tywl
     >      , vion, vlp0, volm, xpnt, ypnt
     >      , dsgt, iwgt, ipgt, icgt, nfl0x, nfl0y, nfl0z
     >      , iwpbr_wall
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'bvx deallocate error in cntcom_alc, istat = ', istat
            call wexit( 'cntcom', trim( cmsg ) )
          endif
        endif
      end select

      return
      end
