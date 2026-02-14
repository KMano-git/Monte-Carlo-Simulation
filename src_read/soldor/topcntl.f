! added integrated calculation with TOPICS by kamata 20/11/16
! control for integrated calculation with TOPICS
      subroutine topcntl( kflg )
      use cplcom,      only : dtmin
      use csonic,      only : dtim, mxcpu, time
      use cunit,       only : n6
      use mod_keylist, only : knend, kstop
      use mod_loc_tim, only : msts
      use topics_mod,  only : dtcal, dtduc, dtim_s, kcon, lexdt, tduc
      implicit none
! arguments
      integer, intent(in) :: kflg
! kflg : = 1 : recieve data from TOPICS
!        = 2 : send    data to   TOPICS

! local variables
      real(8) :: dtminh = 0.0_8
 
      dtminh = dtmin * 0.5_8

      select case( kflg )
      case( 1 )
        if( dtduc < dtminh ) then
! recieve core profile data from TOPICS
          call rcvftop( 2 )
          call set_top_sr( 3 )
          if( lexdt >= 0 ) then
            tduc  = time + dtcal
            dtduc = dtcal
          endif
          write(n6,*) 'topcntl rcvftop 2 kflg,lexdt,time,dtcal = '
     >        , kflg, lexdt, time, dtcal
        else
          if( dtim + dtminh >= dtduc ) then
            dtim_s = dtim
            dtim   = dtduc
          endif
        endif
      case( 2 )
! dtduc update
        dtduc = dtduc - dtim
        if( msts == kstop ) then
! send kcon to TOPICS exiting due to error
          kcon = 1
          call sndttop( 2 )
        elseif( msts == knend .and. mxcpu == 999999 ) then
! send kcon to TOPICS exiting due to time limit
          kcon = -1
          call sndttop( 2 )
        endif
        if( lexdt > 0 ) then
          if( abs( dtduc ) < dtminh ) then
! send edge data to TOPICS
            kcon = 0
            call set_top_sr( 4 )
            call sndttop( 2 )
! reset dtim
            dtim  = dtim_s
! near equal 0.0 to 0.0
            dtduc = 0.0_8
          elseif( dtduc < -dtminh ) then
! time control error
            kcon = 1
            call sndttop( 2 )
          endif
        endif
      end select

      return
      end
