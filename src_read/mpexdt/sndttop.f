! added integrated calculation with TOPICS by kamata 20/11/16
! send data to TOPICS in SONIC
      subroutine sndttop( kflg )
      use cplmet,      only : icaxs, icmpe, icmps
      use cunit,       only : n6
      use mod_dtypdef, only : nddt, typ_dnam, typ_itag
      use mod_mpicomm, only : nwld_all
      use mpi
      use topics_mod,  only : gdn, gmion, gnro, gro, gtim, gvr, kden
     >    , ktem, ktoprnk, lhist, nddnz, ndmk
      implicit none
! arguments
      integer, intent(in) :: kflg
! kflg : = 1 : mesh data
!        = 2 : core profile data

! local variables
! modified 1/1 lines TOPICS calculation using arbitrary SONIC mesh data by kamata 2022/08/14
!ik   integer   :: ierr, iroed, itag = 1, ityp, lnm, nedgs, nedgt, nmans
      integer   :: ierr, itag = 1, ityp, lnm
! modified 1/1 lines addition and reconsider of passed data by kamata 2022/02/23
!ik  >           , nroblm
! modified 1/1 lines TOPICS calculation using arbitrary SONIC mesh data by kamata 2022/08/14
!ik  >           , nroblm, nsiz
     >           , nsiz
      logical   :: lcre(2) = .false.
      character    dnam*10

      select case( kflg )
      case( 1 )
! mesh data
!- set type struct
        dnam = 'MESH'
        lnm = len_trim( dnam )
        if( .not. lcre(1) ) then
          call strct_sonic( dnam(1:lnm) )
          lcre(1) = .true.
        endif

!- send data
        call tbfind( nddt, typ_dnam, ityp, dnam(1:lnm) )
        ityp = typ_itag(ityp)
        call mpi_send( icaxs, 1, ityp, ktoprnk, itag, nwld_all, ierr )
! added 2 lines TOPICS calculation using arbitrary SONIC mesh data by kamata 2022/08/14
! modified 1/1 lines changed base of sent/received data from ρ_vol to ρ_psi by kamata 2022/11/24
!ik     call mpi_send( gvl, gnro, MPI_REAL8, ktoprnk, itag, nwld_all
        call mpi_send( gro, gnro, MPI_REAL8, ktoprnk, itag, nwld_all
     >               , ierr )

! save
! deleted 2 lines TOPICS calculation using arbitrary SONIC mesh data by kamata 2022/08/14
!ik     nroblm = nroblk - 1
!ik     iroed  = nroblm + icmps - icmpe

! check
! deleted 3 lines TOPICS calculation using arbitrary SONIC mesh data by kamata 2022/08/14
!ik     nmans = icaxs  - icmps + 1
!ik     nedgs = icmpe  - icmps + 1
!ik     nedgt = nroblm - iroed + 1
  
        write(n6,900)
     >      'icaxs = ', icaxs, 'icmpe = ', icmpe, 'icmps = ', icmps
! deleted 7 lines TOPICS calculation using arbitrary SONIC mesh data by kamata 2022/08/14
!ik     write(n6,910)
!ik  >      'nman = ', nroblm, nmans, 'nedg = ', nedgt, nedgs
!ik     write(n6,900)
!ik  >      'iroed = ', iroed

!ik     if( nroblm /= nmans .or. nedgt /= nedgs )
!ik  >    call wexit( 'sndttop', 'wrong mesh number' )
      case( 2 )
! edge data
!- set type struct
        dnam = 'EDGE'
        lnm = len_trim( dnam )
        if( .not. lcre(2) ) then
          call strct_sonic( dnam(1:lnm) )
          lcre(2) = .true.
        endif

!-  send data
        call tbfind( nddt, typ_dnam, ityp, dnam(1:lnm) )
        ityp = typ_itag(ityp)
! modified 1/6 lines addition and reconsider of passed data by kamata 2022/02/23
!ik     call mpi_send( gvr, 1, ityp, ktoprnk, itag, nwld_all, ierr )
        call mpi_send( gtim, 1, ityp, ktoprnk, itag, nwld_all, ierr )
! modified 1/1 lines TOPICS calculation using arbitrary SONIC mesh data by kamata 2022/08/14
!ik     nsiz = nro * ( gmion + 1 ) * ndmk
        nsiz = gnro * ( gmion + 1 ) * ndmk
        call mpi_send( gvr, nsiz, MPI_REAL8, ktoprnk, itag, nwld_all
     >               , ierr )
! modified 1/1 lines TOPICS calculation using arbitrary SONIC mesh data by kamata 2022/08/14
!ik     call mpi_send( gro, nro, MPI_REAL8, ktoprnk, itag, nwld_all
! deleted 2 lines changed base of sent/received data from ρ_vol to ρ_psi by kamata 2022/11/24
!ik     call mpi_send( gro, gnro, MPI_REAL8, ktoprnk, itag, nwld_all
!ik  >               , ierr )
! added 6 lines density of SONIC impurities to TOPICS by kamata 2022/03/12
        call mpi_send( nddnz, 3, MPI_INTEGER, ktoprnk, itag, nwld_all
     >               , ierr )
! modified 1/1 lines bug of real number operation exception by kamata 2022/07/04
! returned 1/1 lines TOPICS calculation using arbitrary SONIC mesh data by kamata 2022/08/14
        nsiz = nddnz(1) * nddnz(3)
!ik     nsiz = nddnz(1) * nddnz(2)
        if( nsiz > 0 )
     >    call mpi_send( gdn, nsiz, MPI_REAL8, ktoprnk, itag, nwld_all
     >                 , ierr )

        write(n6,'(2x,"sndttop  gtim = ", es15.7, "  ierr = ", i3 )')
     >      gtim, ierr
        write(n6,'(2x,"ni,ne = ", 2es15.7, 2x, "Ti,Te = ", 2es15.7)')
! modified 3/3 lines TOPICS calculation using arbitrary SONIC mesh data by kamata 2022/08/14
!ik  >      gvr(nroblk,1,kden), gvr(nroblk,0,kden), gvr(nroblk,1,ktem)
!ik  >    , gvr(nroblk,0,ktem)
!ik     if( lhist == 1 ) call dbg_top_sr( 2, 2, gtim, n6 )
     >      gvr(gnro,1,kden), gvr(gnro,0,kden), gvr(gnro,1,ktem)
     >    , gvr(gnro,0,ktem)
! modified 1/1 lines changed base of sent/received data from ρ_vol to ρ_psi by kamata 2022/11/24
!ik     if( lhist == 1 ) call dbg_snc_sr( 2, 2, gtim, gmion, n6 )
        if( lhist == 1 ) call dbg_top_sr( 2, 2, gtim, gmion, n6 )
      end select

! error message
      if( ierr /= 0 )
     >  call wf_exit( 'sndttop', 'data ' // dnam(1:lnm) //
     >      ' not sent.' )
      return
  900 format( 3(2x,a,i5) )
! deleted 1 line TOPICS calculation using arbitrary SONIC mesh data by kamata 2022/08/14
!ik  910 format( 2(2x,a,2i5) )
      end
