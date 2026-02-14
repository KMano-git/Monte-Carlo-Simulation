! added dynamic allocation of arrays by kamata 2022/05/29
! allocate variables in cimcom
      subroutine cimcom_alc( kk )
      use cimcom, only : cdLz_cx, cdLz_i, cdLz_r, cdi, cdi0, cdr
     >    , denFlxZtoWall, denZ, dfcf, eneFlxZtoWall, flmxe_at_impmc
     >    , flmxi_at_impmc, friZ, gdez, gdte, gdti, glte, glti
     >    , heat_flux_para_at_ic_center, heat_flux_para_e_ic_center
     >    , ionZ, phdnz, phfrz, phionz, phrecz, phthz, phvlz, phwci
     >    , phwrd, recZ, slnv, slw0, slw1, temZ, thfZ, vzpZ, weflx
     >    , wpflx, wsct, v_imp_conv
      use csize,  only : ndmc, ndmis, ndwp
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
        if( .not. allocated( cdLz_cx ) ) then
          allocate( cdLz_cx(0:ndmis,0:ndmc), cdLz_i(0:ndmis,0:ndmc)
     >      , cdLz_r(0:ndmis,0:ndmc), cdi(0:ndmis,0:ndmc), cdi0(0:ndmc)
     >      , cdr(0:ndmis,0:ndmc), denFlxZtoWall(0:ndmis,ndwp)
     >      , denZ(0:ndmis,0:ndmc), dfcf(0:ndmis,ndmc)
     >      , eneFlxZtoWall(0:ndmis,ndwp), flmxe_at_impmc(0:ndmc)
     >      , flmxi_at_impmc(0:ndmc), friZ(0:ndmis,0:ndmc), gdez(0:ndmc)
     >      , gdte(0:ndmc), gdti(0:ndmc), glte(0:ndmc), glti(0:ndmc)
     >      , heat_flux_para_at_ic_center(0:ndmc)
     >      , heat_flux_para_e_ic_center(0:ndmc), ionZ(0:ndmis,0:ndmc)
     >      , phdnz(0:ndmis,ndmc), phfrz(0:ndmis,ndmc)
     >      , phionz(0:ndmis,ndmc), phrecz(0:ndmis,ndmc)
     >      , phthz(0:ndmis,ndmc), phvlz(0:ndmis,ndmc), phwci(ndmc)
     >      , phwrd(ndmc), recZ(0:ndmis,0:ndmc), slnv(0:ndmc,2)
     >      , slw0(0:ndmc,2), slw1(0:ndmc,2), temZ(0:ndmis,0:ndmc)
     >      , thfZ(0:ndmis,0:ndmc), vzpZ(0:ndmis,0:ndmc)
     >      , weflx(0:ndmis,ndwp,5), wpflx(0:ndmis,ndwp,5), wsct(0:ndmc)
     >      , v_imp_conv(0:ndmis,ndmc)
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'cdLz_cx allocate error in cimcom_alc, istat = ', istat
            call wexit( 'cimcom', trim( cmsg ) )
          endif

! initial set
          cdLz_cx(0:ndmis,0:ndmc) = 0.0_8
          cdLz_i(0:ndmis,0:ndmc) = 0.0_8
          cdLz_r(0:ndmis,0:ndmc) = 0.0_8
          cdi(0:ndmis,0:ndmc) = 0.0_8
          cdi0(0:ndmc) = 0.0_8
          cdr(0:ndmis,0:ndmc) = 0.0_8
          denFlxZtoWall(0:ndmis,ndwp) = 0.0_8
          denZ(0:ndmis,0:ndmc) = 0.0_8
          dfcf(0:ndmis,1:ndmc) = 0.0_8
          eneFlxZtoWall(0:ndmis,ndwp) = 0.0_8
          flmxe_at_impmc(0:ndmc) = 0.0_8
          flmxi_at_impmc(0:ndmc) = 0.0_8
          friZ(0:ndmis,0:ndmc) = 0.0_8
          gdez(0:ndmc) = 0.0_8
          gdte(0:ndmc) = 0.0_8
          gdti(0:ndmc) = 0.0_8
          glte(0:ndmc) = 0.0_8
          glti(0:ndmc) = 0.0_8
          heat_flux_para_at_ic_center(0:ndmc) = 0.0_8
          heat_flux_para_e_ic_center(0:ndmc) = 0.0_8
          ionZ(0:ndmis,0:ndmc) = 0.0_8
          phdnz(0:ndmis,ndmc) = 0.0_8
          phfrz(0:ndmis,ndmc) = 0.0_8
          phionz(0:ndmis,ndmc) = 0.0_8
          phrecz(0:ndmis,ndmc) = 0.0_8
          phthz(0:ndmis,ndmc) = 0.0_8
          phvlz(0:ndmis,ndmc) = 0.0_8
          phwci(ndmc) = 0.0_8
          phwrd(ndmc) = 0.0_8
          recZ(0:ndmis,0:ndmc) = 0.0_8
          slnv(0:ndmc,2) = 0.0_8
          slw0(0:ndmc,2) = 0.0_8
          slw1(0:ndmc,2) = 0.0_8
          temZ(0:ndmis,0:ndmc) = 0.0_8
          thfZ(0:ndmis,0:ndmc) = 0.0_8
          vzpZ(0:ndmis,0:ndmc) = 0.0_8
          weflx(0:ndmis,ndwp,5) = 0.0_8
          wpflx(0:ndmis,ndwp,5) = 0.0_8
          wsct(0:ndmc) = 0.0_8
          v_imp_conv(0:ndmis,1:ndmc) = 0.0_8
        endif
      case( 2 )
! deallocate
        if( allocated( cdLz_cx ) ) then
          deallocate( cdLz_cx, cdLz_i, cdLz_r, cdi, cdi0, cdr
     >      , denFlxZtoWall, denZ, dfcf, eneFlxZtoWall, flmxe_at_impmc
     >      , flmxi_at_impmc, friZ, gdez, gdte, gdti, glte, glti
     >      , heat_flux_para_at_ic_center, heat_flux_para_e_ic_center
     >      , ionZ, phdnz, phfrz, phionz, phrecz, phthz, phvlz, phwci
     >      , phwrd, recZ, slnv, slw0, slw1, temZ, thfZ, vzpZ, weflx
     >      , wpflx, wsct, v_imp_conv
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'cdLz_cx deallocate error in cimcom_alc, istat = ', istat
            call wexit( 'cimcom', trim( cmsg ) )
          endif
        endif
      end select

      return
      end
