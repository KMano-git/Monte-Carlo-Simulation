! added integrated calculation with TOPICS by kamata 20/11/16
! remade addition and reconsider of passed data by kamata 2022/02/23
! set of send / receive data for SONIC
      subroutine set_top_sr( ktyp )
      use cgdcom,     only : mpman, psman
      use cntcom,     only : mcel, volm
      use cntmnt,     only : ssn, swe, swi
      use cplcom,     only : anamp, anemp, animp, atemp, atimp, vdda
     >    , vdxe, vdxi, vna, vte, vti
      use cplimp,     only : wmc_nty
      use cplmet,     only : icaxs, icmpe, icmps, jcxp1, jcxp2
      use cpmpls,     only : adamp, ahpmp, avpmp, axemp, aximp, flxni
     >    , flxqe, flxqi
!s* for debug by kamata 2022/02/23 ********
     >    , jmd1
!e* for debug by kamata 2022/02/23 ********
      use cplwrd,     only : wime
      use csize,      only : ndmis
      use csonic,     only : time
      use ctopics,    only : ktop_di, ktop_fni, ktop_fqe, ktop_fqi
     >    , ktop_hp, ktop_ne, ktop_ni, ktop_pp, ktop_te, ktop_ti
     >    , ktop_xe, ktop_xi
      use cunit,      only : n6
      use topics_mod, only : gdn, gmion, gnro, gro, groh, gtim, gvr
     >    , iybf, kden, khdf, khfl, khsn, khsr, khvp, kpdf, kpfl, kpsn
     >    , kpvp, ktem, nddnz, ndmk
      implicit none
! from sonicV1/src/mpexdt/exd_tokrcv.f
!        icmps                      icmpe                     icaxs
!         V                          V                         V
!        |31|32|33|34|35|36|37|38|39|40|41|42|43|....   |78|79|80|
!
!        nro-1                      nrbc (iybc)
!         V                          V
!        |50|49|48|47|46|45|44|43|42|41|40|39|38|....   |03|02|01|
!                                   A  A
!                                  42  41 
!                                 (iybf)
!
!        iycn = nro-1 + icmps = 81
!        iybc = iycn  - icmpe = 81 - 40 = 41
!        iybf = iybc + 1 = 41 + 1  = 42

! arguments
      integer, intent(in) :: ktyp
! ktyp : data type
!        = 1 : TOPICS mesh data,    = 2 : SONIC mesh data
!        = 3 : TOPICS profile data, = 4 : SONIC edge data

! local variables
      integer, parameter :: kion = 1
! modified 1/1 lines TOPICS calculation using arbitrary SONIC mesh data by kamata 2022/08/14
!ik   integer :: i, ia = 1, ic, is, iy, iyp, j, jcxen, jcxst, nroblm
      integer :: i, ia = 1, ic, is, iy, iyp, j, jcxen, jcxst, nrom
! added 2 lines density of SONIC impurities to TOPICS by kamata 2022/03/12
      integer :: icnd, nst = 0
! modified 3/3 lines TOPICS calculation using arbitrary SONIC mesh data by kamata 2022/08/14
!ik   real(8)    dgvl(nro), gvl(nro), sumd(10), wdfh(nro,0:kion)
!ik  >         , wdnh(nro,0:kion), whnh(nro,0:kion), whrh(nro)
!ik  >         , wkih(nro,0:kion), wpnh(nro,0:kion), wtmh(nro,0:kion)
      real(8)    dgvl(gnro), sumd(10), wdfh(gnro,0:kion)
     >         , wdnh(gnro,0:kion), whnh(gnro,0:kion), whrh(gnro)
     >         , wkih(gnro,0:kion), wpnh(gnro,0:kion), wtmh(gnro,0:kion)
     >         , zdvl, ztvl
! added 1 line bug for SONIC + IMPACT by kamata 2022/09/05
! deleted 1 line changed base of sent/received data from ρ_vol to ρ_psi by kamata 2022/11/24
!ik   real(8)    gvlh(gnro-1)

! added 7 lines bug of real number operation exception by kamata 2022/07/04
! modified 7/7 lines TOPICS calculation using arbitrary SONIC mesh data by kamata 2022/08/14
!ik   wdfh(1:nro,0:kion) = 0.0_8
!ik   wdnh(1:nro,0:kion) = 0.0_8
!ik   whnh(1:nro,0:kion) = 0.0_8
!ik   whrh(1:nro) = 0.0_8
!ik   wkih(1:nro,0:kion) = 0.0_8
!ik   wpnh(1:nro,0:kion) = 0.0_8
!ik   wtmh(1:nro,0:kion) = 0.0_8
      wdfh(1:gnro,0:kion) = 0.0_8
      wdnh(1:gnro,0:kion) = 0.0_8
      whnh(1:gnro,0:kion) = 0.0_8
      whrh(1:gnro) = 0.0_8
      wkih(1:gnro,0:kion) = 0.0_8
      wpnh(1:gnro,0:kion) = 0.0_8
      wtmh(1:gnro,0:kion) = 0.0_8

      select case( ktyp )
      case( 1 ) ! recieve
      case( 2 ) ! send
! added 14 lines TOPICS calculation using arbitrary SONIC mesh data by kamata 2022/08/14
! modified 14/3 lines changed base of sent/received data from ρ_vol to ρ_psi by kamata 2022/11/24
!ik     ztvl    = 0.0_8
!ik     iy      = 1
!ik     gvl(iy) = 0.0_8
!ik     do i = icaxs, icmps, -1
!ik       iy  = iy + 1
!ik       if( iy > gnro )
!ik  >       call wexit( 'set_top_sr', 'mismatched number of meshes.' )
!ik       zdvl = 0.0_8
!ik       do j = jcxp1+1, jcxp2-1
!ik         zdvl = zdvl + volm(mcel(j,i))
!ik       enddo
!ik       ztvl = ztvl + zdvl
!ik       gvl(iy) = ztvl
!ik     enddo
        gro(1:mpman) = sqrt( ( psman(1) - psman(1:mpman) )
     >                     / ( psman(1) - psman(mpman ) ) )
        groh(1:mpman-1) = ( gro(1:mpman-1) + gro(2:mpman) ) * 0.5_8
      case( 3 ) ! recieve
        time = gtim

        iyp  = iybf - 1
        if( ktop_fni == 1 )
     >    flxni = ( gvr(iyp,ia,kpfl) + gvr(iybf,ia,kpfl) ) * 0.5_8 ! ion particle flux
        if( ktop_fqe == 1 )
     >    flxqe = ( gvr(iyp, 0,khfl) + gvr(iybf, 0,khfl) ) * 0.5_8 ! ele. heat flux
        if( ktop_fqi == 1 )
     >    flxqi = ( gvr(iyp,ia,khfl) + gvr(iybf,ia,khfl) ) * 0.5_8 ! ion heat flux

! ni, Ti, Te in core plasma (cell center)
! Da, Xi, Xe in core plasma (cell center)
        iy = 0
        do i = icaxs, icmps, -1
          iy  = iy + 1
          iyp = iy + 1

          if( ktop_ne == 1 )
     >      anemp(i)    = ( gvr(iy, 0,kden) + gvr(iyp, 0,kden) ) * 0.5_8
          if( ktop_ni == 1 ) then
            anamp(i,ia) = ( gvr(iy,ia,kden) + gvr(iyp,ia,kden) ) * 0.5_8
            animp(i)    = anamp(i,ia)
          endif
          if( ktop_te == 1 )
     >      atemp(i)    = ( gvr(iy, 0,ktem) + gvr(iyp, 0,ktem) ) * 0.5_8
          if( ktop_ti == 1 )
     >      atimp(i)    = ( gvr(iy,ia,ktem) + gvr(iyp,ia,ktem) ) * 0.5_8

          if( ktop_di == 1 ) adamp(i,ia) = gvr(iy,ia,kpdf)
          if( ktop_xe == 1 ) axemp(i)    = gvr(iy, 0,khdf)
          if( ktop_xi == 1 ) aximp(i)    = gvr(iy,ia,khdf)
          if( ktop_pp == 1 ) avpmp(i,ia) = gvr(iy,ia,kpvp)
          if( ktop_hp == 1 ) ahpmp(i,ia) = gvr(iy,ia,khvp)
        enddo

! 2D data  ( exclude dummy cells in the main plasma )
        call pldcof

        write(n6,'(2x,a,i3,2x,a,3es12.3)') 'Flux iybf = ', iybf
     >      , 'flxni,flxqi,flxqe = ', flxni, flxqi, flxqe
      case( 4 ) ! send
! modified 1/1 lines TOPICS calculation using arbitrary SONIC mesh data by kamata 2022/08/14
!ik     nroblm = nroblk - 1
        nrom = gnro - 1
! clear
! modified 1/1 lines TOPICS calculation using arbitrary SONIC mesh data by kamata 2022/08/14
!ik     gvr(1:nro,0:gmion,1:ndmk) = 0.0_8
        gvr(1:gnro,0:gmion,1:ndmk) = 0.0_8

! plasma parameter at cell center
        jcxst = jcxp1 + 1
        jcxen = jcxp2 - 1

        ztvl   = 0.0_8
! deleted 1 line changed base of sent/received data from ρ_vol to ρ_psi by kamata 2022/11/24
!ik     gvl(1) = 0.0_8
        iy = 0
        do i = icaxs, icmps, -1  ! gnro-1
          zdvl = 0.0_8
          do j = jcxst, jcxen
            ic = mcel(j,i)
            zdvl = zdvl + volm(ic)
          enddo
          ztvl = ztvl + zdvl
          iy  = iy + 1
          dgvl(iy)  = zdvl
! deleted 1 line changed base of sent/received data from ρ_vol to ρ_psi by kamata 2022/11/24
!ik       gvl(iy+1) = ztvl

! hot core
          if( i > icmpe ) then
            wdnh(iy,1) = anamp(i,ia)
            wtmh(iy,0) = atemp(i)
            wtmh(iy,1) = atimp(i)
! edge
          else
            sumd(1:10) = 0.0_8
            do j = jcxst, jcxen
              ic = mcel(j,i)
              sumd( 1) = sumd( 1) + vna(j,i,ia)  * volm(ic)
              sumd( 2) = sumd( 2) + vte(j,i)     * volm(ic)
              sumd( 3) = sumd( 3) + vti(j,i)     * volm(ic)
              sumd( 4) = sumd( 4) + vdda(j,i,ia) * volm(ic)
              sumd( 5) = sumd( 5) + vdxe(j,i)    * volm(ic)
              sumd( 6) = sumd( 6) + vdxi(j,i)    * volm(ic)
              sumd( 7) = sumd( 7) + ssn(j,i,ia)  * volm(ic)
              sumd( 8) = sumd( 8) + swe(j,i)     * volm(ic)
              sumd( 9) = sumd( 9) + swi(j,i)     * volm(ic)
              sumd(10) = sumd(10) + wime(ic)     * volm(ic)
            enddo
            wdnh(iy,1) = sumd( 1) / dgvl(iy) ! ni
            wtmh(iy,0) = sumd( 2) / dgvl(iy) ! Te
            wtmh(iy,1) = sumd( 3) / dgvl(iy) ! Ti
            wdfh(iy,1) = sumd( 4) / dgvl(iy) ! Di
            wkih(iy,0) = sumd( 5) / dgvl(iy) ! Xe
            wkih(iy,1) = sumd( 6) / dgvl(iy) ! Xi
            wpnh(iy,1) = sumd( 7) / dgvl(iy) ! particle source by neutral i
            whnh(iy,0) = sumd( 8) / dgvl(iy) ! energy   source by neutral e
            whnh(iy,1) = sumd( 9) / dgvl(iy) ! energy   source by neutral i
            whrh(iy)   = sumd(10) / dgvl(iy) ! radiation loss due to impurities
          endif
        enddo

! separatrix 
        i  = icmps - 1 ! gnro
        zdvl = 0.0_8
        do j = jcxst, jcxen
          ic = mcel(j,i)
          zdvl = zdvl + volm(ic)
        enddo
        iy = iy + 1
        dgvl(iy) = zdvl
        sumd(1:3) = 0.0_8
        do j = jcxst, jcxen
          ic = mcel(j,i)
          sumd(1) = sumd(1) + vna(j,i,ia) * volm(ic)
          sumd(2) = sumd(2) + vte(j,i)    * volm(ic)
          sumd(3) = sumd(3) + vti(j,i)    * volm(ic)
        enddo
        wdnh(iy,1) = sumd(1) / dgvl(iy) ! ni
        wtmh(iy,0) = sumd(2) / dgvl(iy) ! Te
        wtmh(iy,1) = sumd(3) / dgvl(iy) ! Ti

! set gtim, gro
        gtim = time
! modified 1/1 lines TOPICS calculation using arbitrary SONIC mesh data by kamata 2022/08/14
!ik     gro(1:nroblk) = sqrt( gvl(1:nroblk) / gvl(nroblk) )
! deleted 1 line changed base of sent/received data from ρ_vol to ρ_psi by kamata 2022/11/24
!ik     gro(1:gnro) = sqrt( gvl(1:gnro) / gvl(gnro) )
! added 5 lines bug for SONIC + IMPACT by kamata 2022/09/05
! deleted 5 lines changed base of sent/received data from ρ_vol to ρ_psi by kamata 2022/11/24
!ik     gvlh(1) = gvl(2) * 0.5_8
!ik     do i = 2, gnro-1
!ik       gvlh(i) = ( gvl(i) + gvl(i+1) ) * 0.5_8
!ik     enddo
!ik     groh(1:gnro-1) = sqrt( gvlh(1:gnro-1) / gvl(gnro) )

! set gvr
! modified 1/1 lines TOPICS calculation using arbitrary SONIC mesh data by kamata 2022/08/14
!ik     do iy = 1, nroblm
        do iy = 1, nrom
          gvr(iy+1,1,kden) = ( wdnh(iy,1) + wdnh(iy+1,1) ) * 0.5_8 ! ni
          gvr(iy+1,0,ktem) = ( wtmh(iy,0) + wtmh(iy+1,0) ) * 0.5_8 ! Te
          gvr(iy+1,1,ktem) = ( wtmh(iy,1) + wtmh(iy+1,1) ) * 0.5_8 ! Ti
          gvr(iy,1,kpdf)   = wdfh(iy,1)                            ! Di
          gvr(iy,0,khdf)   = wkih(iy,0)                            ! Xe
          gvr(iy,1,khdf)   = wkih(iy,1)                            ! Xi
          gvr(iy,1,kpsn)   = wpnh(iy,1)                            ! Sn
          gvr(iy,0,khsn)   = whnh(iy,0)                            ! We
          gvr(iy,1,khsn)   = whnh(iy,1)                            ! Wi
          gvr(iy,0,khsr)   = whrh(iy)                              ! CR
        enddo
!- axis
        gvr(1,1,kden) = wdnh(1,1)
        gvr(1,0,ktem) = wtmh(1,0)
        gvr(1,1,ktem) = wtmh(1,1)

        write(n6,'(10(5x,a,8x))')
     >      'ni', 'Te', 'Ti', 'Di', 'Xe', 'Xi', 'Sn', 'We', 'Wi', 'CR'
        write(n6,'(10es15.7)')
! modified 4/4 lines TOPICS calculation using arbitrary SONIC mesh data by kamata 2022/08/14
!ik  >      gvr(nroblk,1,kden), gvr(nroblk,0,ktem), gvr(nroblk,1,ktem)
!ik  >    , gvr(nroblm,1,kpdf), gvr(nroblm,0,khdf), gvr(nroblm,1,khdf)
!ik  >    , gvr(nroblm,1,kpsn), gvr(nroblm,0,khsn), gvr(nroblm,1,khsn)
!ik  >    , gvr(nroblm,0,khsr)
     >      gvr(gnro,1,kden), gvr(gnro,0,ktem), gvr(gnro,1,ktem)
     >    , gvr(nrom,1,kpdf), gvr(nrom,0,khdf), gvr(nrom,1,khdf)
     >    , gvr(nrom,1,kpsn), gvr(nrom,0,khsn), gvr(nrom,1,khsn)
     >    , gvr(nrom,0,khsr)
!s* for debug by kamata 2022/02/23 ********
        write(n6,'(a,es15.7)') 'ni_mid = '
     >    , (vna(jmd1,icmps,ia)+vna(jmd1,icmps-1,ia))*0.5_8
!e* for debug by kamata 2022/02/23 ********

!- impurities
        do is = 2, gmion
! modified 3/3 lines TOPICS calculation using arbitrary SONIC mesh data by kamata 2022/08/14
!ik       gvr(1:nroblk,is,ktem) = gvr(1:nroblk,1,ktem)
!ik       gvr(1:nroblm,is,kpdf) = gvr(1:nroblm,1,kpdf)
!ik       gvr(1:nroblm,is,khdf) = gvr(1:nroblm,1,khdf)
          gvr(1:gnro,is,ktem) = gvr(1:gnro,1,ktem)
          gvr(1:nrom,is,kpdf) = gvr(1:nrom,1,kpdf)
          gvr(1:nrom,is,khdf) = gvr(1:nrom,1,khdf)
        enddo
!- electron
! modified 1/1 lines TOPICS calculation using arbitrary SONIC mesh data by kamata 2022/08/14
!ik     gvr(1:nroblm,0,kpdf) = gvr(1:nroblm,1,kpdf)
        gvr(1:nrom,0,kpdf) = gvr(1:nrom,1,kpdf)

! added 10 lines density of SONIC impurities to TOPICS by kamata 2022/03/12
! calculation of impurity density for each valence
        if( .not. allocated ( gdn ) ) then
! modified 1/1 lines bug for SONIC + IMPACT by kamata 2022/09/05
!ik       nddnz(1) = icmpe - icmps + 2 !  on mesh ( add gnro )
          nddnz(1) = icmpe - icmps + 1 !  half mesh
          nddnz(2) = ndmis * wmc_nty
! modified 1/1 lines TOPICS calculation using arbitrary SONIC mesh data by kamata 2022/08/14
!ik       nst = nroblk - nddnz(1) + 1
! modified 1/1 lines bug for SONIC + IMPACT by kamata 2022/09/05
!ik       nst = gnro - nddnz(1) + 1
          nst = gnro - nddnz(1)
          allocate( gdn(nddnz(1),nddnz(2)), stat = icnd )
          if( icnd /= 0 )
     >     call wexit( 'set_top_sr', 'gdn allocate error.' )
        endif
! modified 1/1 lines bug for SONIC + IMPACT by kamata 2022/09/05
!ik     call cal_denz( gro(nst), gdn, nddnz(1), nddnz(2) )
        call cal_denz( groh(nst), gdn, nddnz(1), nddnz(2) )
      end select

      return
      end
