! added replace all include files with module files by kamata 2021/08/18
! data copy of vwork3
!  IMPMC/inc/cimcom /cimcom_12b/ or IMPMC_TD/inc/cimcom /cimcom_12b/ 
!  or ( IMPMC/imhist.f /com_imhist/ --> clocal.f )
      subroutine setvwork3( ktd, kk, work3, md )
      use cimcom, only : hcnpt, hcnwt, hcpum, hgrpn, hgrpp, hgrpw, hmype
     >    , hpmax, hregn, hregp, hregw, hrgpa, hrgpc, hstep, htout
      use clocal, only : xcnpt, xcnwt, xcpum, xmype, xpmax, xrgpa, xrgpc
     >    , xstep, xtout
      implicit none
! arguments
      integer, intent(in)    :: kk, ktd, md
      real(8), intent(inout) :: work3(md)
! kk    : flag to set data, = 1 : work3 to module variables
!                           = 2 : module variables to work3
! ktd   : calculation mode of IMPMC or local common
!         =  0 : steady state cimcom.f, =  1 : time dependent cimcom.f
!         = -1 : local common in IMPMC/imhist.f /com_imhist/ --> clocal.f
! md    : dimension size
! work3 : work area

! local variables
      integer    i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, ia
      integer :: nd = 11, nd1 = 15, nd2 = 10

      select case( ktd )
! local common in steady state
      case( -1 )
! set start posiotion of each variable
        i0 = 1       ! xmype
        i1 = i0 + 1  ! xtout
        i2 = i1 + 1  ! xstep
        i3 = i2 + 1  ! xcpum
        i4 = i3 + 1  ! xpmax
        i5 = i4 + 1  ! xrgpa
        i6 = i5 + nd ! xrgpc
        i7 = i6 + nd ! xcnpt
        i8 = i7 + nd ! xcnwt
      
        select case ( kk )
        case( 1 )
! work3 to module variables
          xmype = work3(i0) 
          xtout = work3(i1)
          xstep = work3(i2)
          xcpum = work3(i3)
          xpmax = work3(i4)

          call datcp( work3(i5), xrgpa, nd )
          call datcp( work3(i6), xrgpc, nd )
          call datcp( work3(i7), xcnpt, nd )
          call datcp( work3(i8), xcnwt, nd )
        case( 2 )
! module variables to work3
          work3(i0) = xmype 
          work3(i1) = xtout
          work3(i2) = xstep
          work3(i3) = xcpum
          work3(i4) = xpmax

          call datcp( xrgpa, work3(i5), nd )
          call datcp( xrgpc, work3(i6), nd )
          call datcp( xcnpt, work3(i7), nd )
          call datcp( xcnwt, work3(i8), nd )
        end select
! steady state cimcom
      case( 0 )
! set start posiotion of each variable
        i0 = 1       ! hmype
        i1 = i0 + 1  ! htout
        i2 = i1 + 1  ! hstep
        i3 = i2 + 1  ! hcpum
        i4 = i3 + 1  ! hpmax
        i5 = i4 + 1  ! hrgpa
        i6 = i5 + nd ! hrgpc
        i7 = i6 + nd ! hcnpt
        i8 = i7 + nd ! hcnwt

        select case ( kk )
        case( 1 )
! work3 to module variables
          hmype = work3(i0) 
          htout = work3(i1)
          hstep = work3(i2)
          hcpum = work3(i3)
          hpmax = work3(i4)

          call datcp( work3(i5), hrgpa, nd )
          call datcp( work3(i6), hrgpc, nd )
          call datcp( work3(i7), hcnpt, nd )
          call datcp( work3(i8), hcnwt, nd )
        case( 2 )
! module variables to work3
          work3(i0) = hmype 
          work3(i1) = htout
          work3(i2) = hstep
          work3(i3) = hcpum
          work3(i4) = hpmax

          call datcp( hrgpa, work3(i5), nd )
          call datcp( hrgpc, work3(i6), nd )
          call datcp( hcnpt, work3(i7), nd )
          call datcp( hcnwt, work3(i8), nd )
        end select
! time dependent cimcom
      case( 1 )
! set start posiotion of each variable
        i0 = 1        ! hmype
        i1 = i0 + 1   ! htout
        i2 = i1 + 1   ! hstep
        i3 = i2 + 1   ! hcpum
        i4 = i3 + 1   ! hpmax
        i5 = i4 + 1   ! hgrpp
        i6 = i5 + nd1 ! hgrpw
        i7 = i6 + nd1 ! hgrpn
        i8 = i7 + nd1 ! hregp
        i9 = i8 + nd2 ! hregw
        ia = i9 + nd2 ! hregn

        select case ( kk )
        case( 1 )
! work3 to module variables
          hmype = work3(i0) 
          htout = work3(i1)
          hstep = work3(i2)
          hcpum = work3(i3)
          hpmax = work3(i4)

          call datcp( work3(i5), hgrpp, nd1 )
          call datcp( work3(i6), hgrpw, nd1 )
          call datcp( work3(i7), hgrpn, nd1 )
          call datcp( work3(i8), hregp, nd2 )
          call datcp( work3(i9), hregw, nd2 )
          call datcp( work3(ia), hregn, nd2 )
        case( 2 )
! module variables to work3
          work3(i0) = hmype 
          work3(i1) = htout
          work3(i2) = hstep
          work3(i3) = hcpum
          work3(i4) = hpmax

          call datcp( hgrpp, work3(i5), nd1 )
          call datcp( hgrpw, work3(i6), nd1 )
          call datcp( hgrpn, work3(i7), nd1 )
          call datcp( hregp, work3(i8), nd2 )
          call datcp( hregw, work3(i9), nd2 )
          call datcp( hregn, work3(ia), nd2 )
        end select
      end select

      return
      end
