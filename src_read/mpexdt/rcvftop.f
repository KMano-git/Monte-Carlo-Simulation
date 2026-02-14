! added integrated calculation with TOPICS by kamata 20/11/16
! recieve data form TOPICS in SONIC
      subroutine rcvftop( kflg )
      use cgdcom,      only : mpman
      use cplmet,      only : icaxs, icmpe, icmps
      use cunit,       only : n6
      use mod_dtypdef, only : nddt, typ_dnam, typ_itag
      use mod_mpicomm, only : nwld_all
      use mpi!,         only : mpi_recv, mpi_status_size
      use topics_mod,  only : gmion, gnro, gro, groh, gtim, gvr, iybf
     >    , ktoprnk, lhist, ndmk, nro
      implicit none
! arguments
      integer, intent(in) :: kflg
! kflg : = 1 : mesh data
!        = 2 : core profile data

! local variables
! modified 2/2 lines TOPICS calculation using arbitrary SONIC mesh data by kamata 2022/08/14
!ik   integer   :: ierr, istat(MPI_STATUS_SIZE), itag = 1, ityp, iybc
!ik  >           , iycn, lnm, nroblm, nsedgs, nsedgt, nsmans, nsmant
      integer   :: ierr, istat(MPI_STATUS_SIZE), itag = 1, ityp
     >           , lnm
! added 1 line addition and reconsider of passed data by kamata 2022/02/23
      integer      icond, nsiz
      logical   :: lcre(2) = .false.
      character    dnam*10

      select case( kflg )
      case( 1 )
! mesh data
!- set type struct
        dnam = 'NRO'
        lnm = len_trim( dnam )
        if( .not. lcre(1) ) then
          call strct_sonic( dnam(1:lnm) )
          lcre(1) = .true.
        endif

!- recv data
        call tbfind( nddt, typ_dnam, ityp, dnam(1:lnm) )
        ityp = typ_itag(ityp )
        call mpi_recv( nro, 1, ityp, ktoprnk, itag, nwld_all, istat
     >               , ierr )

!- check
!
!        icmps                      icmpe
!         V                          V
!        |31|32|33|34|35|36|37|38|39|40|41|42|43|....   |78|79|80|
!
!        nroblk-1                   nrbc (iybc)
!         V                          V
!        |50|49|48|47|46|45|44|43|42|41|40|39|38|....   |03|02|01|
!                                   A  A
!                                  42  41 
!                                 (iybf)
!
!        iycn = nroblk-1 + icmps = 50 + 31 = 81
!        iybc = iycn  - icmpe = 81 - 40 = 41
!        iybf = iybc + 1 = 41 + 1  = 42
!
!        nsedg = icmpe    - icmps + 1 = 40 - 31 + 1 = 10
!              = nroblk-1 - nrbc  + 1 = 50 - 41 + 1 = 10
!
!        nsman = icaxs    - icmps + 1 = 80 - 31 + 1 = 50
!              = nroiblk - 1 = 50

! modified 17/3 lines TOPICS calculation using arbitrary SONIC mesh data by kamata 2022/08/14
!ik     nroblm = nroblk - 1
!ik     iycn = nroblm + icmps ! on mesh
!ik     iybc = iycn   - icmpe ! half mesh
!ik     iybf = iybc   + 1     ! on mesh

!ik     nsmant = nroblm
!ik     nsmans = icaxs  - icmps + 1
!ik     nsedgt = nroblm - iybc  + 1
!ik     nsedgs = icmpe  - icmps + 1

!ik     if( nsedgt /= nsedgs .or. nsmant /= nsmans ) then
!ik       write(n6,900) 'icmps, icmpe,    icaxs = ', icmps, icmpe, icaxs
!ik       write(n6,900) 'iycn,  iybc/nrbc, iybf = ', iycn, iybc, iybf
!ik       write(n6,910) 'nsedg = ', nsedgt, nsedgs
!ik  >                , 'nsman = ', nsmant, nsmans
!ik       call wexit( 'rcvftop', 'wrong mesh number' )
!ik     endif
        gnro   = icaxs - icmps + 2 ! SONIC on mesh size for TOPICS
        iybf   = icaxs - icmpe + 2 ! on mesh
! modified 1/3 lines changed base of sent/received data from ρ_vol to ρ_psi by kamata 2022/11/24
!ik     write(n6,900) 'gnro, iybf = ', gnro, iybf
        write(n6,900) 'gnro, mpman, iybf = ', gnro, mpman, iybf
        if( gnro /= mpman )
     >    call wf_exit( 'rcvftop', "gnro doesn't match mpman" )

! added 7 lines addition and reconsider of passed data by kamata 2022/02/23
! modified 1/2 lines TOPICS calculation using arbitrary SONIC mesh data by kamata 2022/08/14
!ik     allocate( gro(nro), gvr(nro,0:gmion,ndmk), stat=icond )
! modified 2/2 lines bug for SONIC + IMPACT by kamata 2022/09/05
!ik     allocate( gro(gnro), gvl(gnro), gvr(gnro,0:gmion,ndmk)
!ik  >          , stat=icond )
! modified 2/2 lines changed base of sent/received data from ρ_vol to ρ_psi by kamata 2022/11/24
!ik     allocate( gro(gnro), groh(gnro-1), gvl(gnro)
!ik  >    , gvr(gnro,0:gmion,ndmk), stat=icond )
        allocate( gro(gnro), groh(gnro), gvr(gnro,0:gmion,ndmk)
     >          , stat=icond )
        if( icond /= 0 ) then
! modified 1/1 lines TOPICS calculation using arbitrary SONIC mesh data by kamata 2022/08/14
!ik       write(n6,900) 'nro, gmion, ndmk = ', nro, gmion, ndmk
          write(n6,900) 'gnro, gmion, ndmk = ', gnro, gmion, ndmk
          call wexit( 'rcvftop', 'gro allocate error' )
        endif 
! modified 2/2 lines TOPICS calculation using arbitrary SONIC mesh data by kamata 2022/08/14
!ik     gro(1:nro) = 0.0_8
!ik     gvr(1:nro,0:gmion,1:ndmk) = 0.0_8
        gro(1:gnro) = 0.0_8
        gvr(1:gnro,0:gmion,1:ndmk) = 0.0_8
! added 1 line bug for SONIC + IMPACT by kamata 2022/09/05
! modified 1/1 lines changed base of sent/received data from ρ_vol to ρ_psi by kamata 2022/11/24
!ik     groh(1:gnro-1) = 0.0_8
        groh(1:gnro) = 0.0_8
      case( 2 )
! core profile data
!- set type struct
        dnam = 'PROF'
        lnm = len_trim( dnam )
        if( .not. lcre(2) ) then
          call strct_sonic( dnam(1:lnm) )
          lcre(2) = .true.
        endif

!- recv data
        call tbfind( nddt, typ_dnam, ityp, dnam(1:lnm) )
        ityp = typ_itag(ityp )
! modified 2/7 lines addition and reconsider of passed data by kamata 2022/02/23
!ik     call mpi_recv( gvr, 1, ityp, ktoprnk, itag, nwld_all, istat
!ik  >               , ierr )
        call mpi_recv( gtim, 1, ityp, ktoprnk, itag, nwld_all, istat
     >               , ierr )
! modified 1/1 lines TOPICS calculation using arbitrary SONIC mesh data by kamata 2022/08/14
!ik     nsiz = nro * ( gmion + 1 ) * ndmk
        nsiz = gnro * ( gmion + 1 ) * ndmk
        call mpi_recv( gvr, nsiz, MPI_REAL8, ktoprnk, itag, nwld_all
     >               , istat, ierr )
! modified 1/1 lines TOPICS calculation using arbitrary SONIC mesh data by kamata 2022/08/14
!ik     call mpi_recv( gro, nro, MPI_REAL8, ktoprnk, itag, nwld_all
! deleted 2 lines changed base of sent/received data from ρ_vol to ρ_psi by kamata 2022/11/24
!ik     call mpi_recv( gro, gnro, MPI_REAL8, ktoprnk, itag, nwld_all
!ik  >               , istat, ierr )
! modified 1/1 lines TOPICS calculation using arbitrary SONIC mesh data by kamata 2022/08/14
!ik     if( lhist == 1 ) call dbg_top_sr( 1, 1, gtim, n6 )
! modified 1/1 lines changed base of sent/received data from ρ_vol to ρ_psi by kamata 2022/11/24
!ik     if( lhist == 1 ) call dbg_snc_sr( 1, 1, gtim, gmion, n6 )
        if( lhist == 1 ) call dbg_top_sr( 1, 1, gtim, gmion, n6 )
      end select
! error message
      if( ierr /= 0 )
     >  call wf_exit( 'rcvftop', 'data ' // dnam(1:lnm) //
     >      ' not recieved.' )

      return
  900 format( 2x, a, 3i5 )
! deleted 1 line TOPICS calculation using arbitrary SONIC mesh data by kamata 2022/08/14
!ik  910 format( 2( 2x, a, 2i5 ) )
      end
