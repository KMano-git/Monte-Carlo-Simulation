! added get the code number from the run-time parameters by kamata 2021/08/18
! remaked treat 4 or more impurities with IMPMC by kamata 2022/04/21
      subroutine  set_cdno( cprg, clev, nprg )
      use mod_mpicomm, only : m6
      implicit none
! arguments
      character, intent(inout) :: clev*(*), cprg*(*)
      integer,   intent(out)   :: nprg
! clev : code name for log (in:'#= IMPMC',out:'#= IMPMC','#= IMPMC2',...)
! cprg : code name         (in:'IMPMC',   out:'IMPMC','IMPMC2',...)
! nprg : affiliation group number

! local variables
      integer      i            ! counter
      integer      lp           ! last position
      integer      narg         ! number of parameters
      character :: codnum(2)*10 ! code number ( character )

      call get_par( codnum, 2, narg )
      if( narg == 0 ) then
        nprg = 1
      elseif( narg == 1 ) then
! set code name
        if( codnum(narg) /= '1' ) then
          lp = len_trim( cprg ) + 1
          cprg(lp:) = trim( codnum(narg) )
          lp = len_trim( clev ) + 1
          clev(lp:) = trim( codnum(narg) )
        endif
        read(codnum(narg),*) nprg
      else
        do i = 1, narg
          write(m6,'(i3,a)') i, ' codnum = ' // trim( codnum(i) )
        enddo 
        call wexit( 'set_cdno', 'run-time parameter error. ' )
      endif

      ! m6 is not defined as %IMPMC_XXX at that time, it is just "fort.17".
      ! Do not use write-stete here for avoid file write error in plasma-simulator.
      ! write(m6,'(a,i5)') 'cprg = ' // trim( cprg ) // ' nprg = ', nprg

      return
      end
