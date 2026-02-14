! added replace all include files with module files by kamata 2021/08/18
! data copy of vwork ( IMPMC/inc/cimcom /cimcom_11/ )
      subroutine setvwork( ktd, kk, work, md )
      use cimcom, only : denz, scpum, sdmy1, sdmy2, sdmy3, sdtmz, sflux
     >    , sitmp, sitmz, sptcl, stimp, stimz, swtot, temz, vzpz, wsct
      use csize,  only : ndmc, ndis => ndmis
      implicit none
! arguments
      integer, intent(in)    :: kk, ktd, md
      real(8), intent(inout) :: work(md)
! kk   : flag to set data, = 1 : work to module variables
!                          = 2 : module variables to work
! ktd  : calculation mode of IMPMC ( = 0 : steady state, = 1 : time dependent )
! md   : dimension size
! work : work area

! local variables
      integer    i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, ia, ib, ic, id
     >         , ie, ig, ih, nd, ndp

! initial set
      ndp = ndmc + 1
      nd  = ( ndis + 1 ) * ndp

      select case( ktd )
! steady state
      case( 0 )
! set start posiotion of each variable
        i0 = 1        ! sflux
        i1 = i0 + 1   ! swtot
        i2 = i1 + 1   ! sptcl
        i3 = i2 + 1   ! stimp
        i4 = i3 + 1   ! sitmp
        i5 = i4 + 1   ! stimz
        i6 = i5 + 1   ! sdtmz
        i7 = i6 + 1   ! sdmy1
        i8 = i7 + 1   ! sdmy2
        i9 = i8 + 1   ! sdmy3
        ia = i9 + 1   ! denz
        ib = ia + nd  ! temz
        ic = ib + nd  ! wsct
        id = ic + ndp ! (last)+1
      
        select case ( kk )
        case( 1 )
! work to module variables
          sflux = work(i0)
          swtot = work(i1)
          sptcl = work(i2)
          stimp = work(i3)
          sitmp = work(i4)
          stimz = work(i5)
          sdtmz = work(i6)
          sdmy1 = work(i7)
          sdmy2 = work(i8)
          sdmy3 = work(i9)

          call datcp( work(ia), denz, nd )
          call datcp( work(ib), temz, nd )

          wsct(0:ndmc) = work(ic:id-1)
        case( 2 )
! module variables to work
          work(i0) = sflux
          work(i1) = swtot
          work(i2) = sptcl
          work(i3) = stimp
          work(i4) = sitmp
          work(i5) = stimz
          work(i6) = sdtmz
          work(i7) = sdmy1
          work(i8) = sdmy2
          work(i9) = sdmy3

          call datcp( denz, work(ia), nd )
          call datcp( temz, work(ib), nd )

          work(ic:id-1) = wsct(0:ndmc)
        end select
! time dependent
      case( 1 )
! set start posiotion of each variable
        i0 = 1        ! sflux
        i1 = i0 + 1   ! swtot
        i2 = i1 + 1   ! sptcl
        i3 = i2 + 1   ! stimp
        i4 = i3 + 1   ! sitmp
        i5 = i4 + 1   ! stimz
        i6 = i5 + 1   ! sitmz
        i7 = i6 + 1   ! sdtmz
        i8 = i7 + 1   ! scpum
        i9 = i8 + 1   ! sdmy1
        ia = i9 + 1   ! sdmy2
        ib = ia + 1   ! sdmy3
        ic = ib + 1   ! denz
        id = ic + nd  ! temz
        ie = id + nd  ! vzpz
        ig = ie + nd  ! wsct
        ih = ig + ndp ! (last)+1

        select case ( kk )
        case( 1 )
! work to module variables
          sflux = work(i0)
          swtot = work(i1)
          sptcl = work(i2)
          stimp = work(i3)
          sitmp = work(i4)
          stimz = work(i5)
          sitmz = work(i6)
          sdtmz = work(i7)
          scpum = work(i8)
          sdmy1 = work(i9)
          sdmy2 = work(ia)
          sdmy3 = work(ib)

          call datcp( work(ic), denz, nd )
          call datcp( work(id), temz, nd )
          call datcp( work(ie), vzpz, nd )

          wsct(0:ndmc) = work(ig:ih-1)
        case( 2 )
! module variables to work
          work(i0) = sflux
          work(i1) = swtot
          work(i2) = sptcl
          work(i3) = stimp
          work(i4) = sitmp
          work(i5) = stimz
          work(i6) = sitmz
          work(i7) = sdtmz
          work(i8) = scpum
          work(i9) = sdmy1
          work(ia) = sdmy2
          work(ib) = sdmy3

          call datcp( denz, work(ic), nd )
          call datcp( temz, work(id), nd )
          call datcp( vzpz, work(ie), nd )

          work(ig:ih-1) = wsct(0:ndmc)
        end select
      end select

      return
      end
