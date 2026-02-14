! added replace all include files with module files by kamata 2021/08/18
! data copy of vworky ( IMPMC/inc/cimcom /cimcom_11y/ )
      subroutine setvworky( ktd, kk, worky, md )
      use cimcom, only : friz, ionz, recz, thfz, vzpz
      use csize,  only : ndmc, ndis => ndmis
      implicit none
! arguments
      integer, intent(in)    :: kk, ktd, md
      real(8), intent(inout) :: worky(md)
! kk    : flag to set data, = 1 : worky to module variables
!                           = 2 : module variables to worky
! ktd   : calculation mode of IMPMC ( = 0 : steady state, = 1 : time dependent )
! md    : dimension size
! worky : work area

! local variables
      integer    i0, i1, i2, i3, i4, nd

! initial set
      nd  = ( ndis + 1 ) * ( ndmc + 1 )

      select case( ktd )
! steady state
      case( 0 )
! set start posiotion of each variable
        i0 = 1       ! friz
        i1 = i0 + nd ! thfz
        i2 = i1 + nd ! vzpz
        i3 = i2 + nd ! ionz
        i4 = i3 + nd ! recz
      
        select case ( kk )
        case( 1 )
! worky to module variables
          call datcp( worky(i0), friz, nd )
          call datcp( worky(i1), thfz, nd )
          call datcp( worky(i2), vzpz, nd )
          call datcp( worky(i3), ionz, nd )
          call datcp( worky(i4), recz, nd )
        case( 2 )
! module variables to worky
          call datcp( friz, worky(i0), nd )
          call datcp( thfz, worky(i1), nd )
          call datcp( vzpz, worky(i2), nd )
          call datcp( ionz, worky(i3), nd )
          call datcp( recz, worky(i4), nd )
        end select
! time dependent
      case( 1 )
! set start posiotion of each variable
        i0 = 1       ! friz
        i1 = i0 + nd ! thfz
        i2 = i1 + nd ! ionz
        i3 = i2 + nd ! recz

        select case ( kk )
        case( 1 )
          call datcp( worky(i0), friz, nd )
          call datcp( worky(i1), thfz, nd )
          call datcp( worky(i2), ionz, nd )
          call datcp( worky(i3), recz, nd )
! worky to module variables
        case( 2 )
! module variables to worky
          call datcp( friz, worky(i0), nd )
          call datcp( thfz, worky(i1), nd )
          call datcp( ionz, worky(i2), nd )
          call datcp( recz, worky(i3), nd )
        end select
      end select

      return
      end
