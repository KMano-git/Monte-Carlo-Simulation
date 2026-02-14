! added rperiodic invocation of IMPMC,monte,soldor/figdat by kamata 2023/08/04
      module dbg_mod
      implicit none

! common use
      integer, parameter :: mtim = 10
      integer   :: iwflg = -6
! iwflg : decimal places (iwflg = -6 : XXX.YYYYYY )
! mtim  : dimension size of timdbg*, dtmdbg*

! IMPMC
      character, target :: cpathi*10 = './OUT_IMP/'
      logical,   target :: kseti = .false.
      integer,   target :: ntimi = 0
      integer           :: kimpsy = 0
      real(8),   target :: timbi = 0.0_8
      real(8),   target :: dtmdbgimp(mtim) = 1.0d10
     >                   , timdbgimp(mtim) = 0.0_8
! cpathi    : data storage directory
! dtmdbgimp : output time increment [s]
!            ( dtmdbgimp(i) between timdbgimp(i) and timdbgimp(i+1) )
! kimpsy    : IMPSY output flag ( = 1 : output, = 0 : no output )
! kseti     : flag for initialization
! ntimi     : number of data points set in timdbgimp
! timbi     : last output time [s]
! timdbgimp : time to change the output time step [s]

! NEUT2D
      character, target :: cpathn*10 = './OUT_NEU/'
      logical,   target :: ksetn = .false.
      integer,   target :: ntimn = 0
      real(8),   target :: timbn = 0.0_8
      real(8),   target :: dtmdbgneu(mtim) = 1.0d10
     >                   , timdbgneu(mtim) = 0.0_8
! cpathn    : data storage directory
! dtmdbgneu : output time increment [s]
!            ( dtmdbgneu(i) between timdbgneu(i) and timdbgneu(i+1) )
! ksetn     : flag for initialization
! ntimn     : number of data points set in timdbgneu
! timbn     : last output time [s]
! timdbgneu : time to change the output time step [s]

! SOLDOR
      character, target :: cpaths*10 = './OUT_SOL/'
      logical,   target :: ksets = .false.
      integer,   target :: ntims = 0
      real(8),   target :: timbs = 0.0_8
      real(8),   target :: dtmdbgsol(mtim) = 1.0d10
     >                   , timdbgsol(mtim) = 0.0_8
! cpaths    : data storage directory
! dtmdbgsol : output time increment [s]
!            ( dtmdbgsol(i) between timdbgsol(i) and timdbgsol(i+1) )
! ksets     : flag for initialization
! ntims     : number of data points set in timdbgsol
! timbs     : last output time [s]
! timdbgsol : time to change the output time step [s]

      contains
! time control of file output for debugging
      subroutine dbgcntl( kc )
      use csonic, only : time
      use cunit,  only : lmspe, lmype
      implicit none
! arguments
      integer, intent(in) :: kc
! kc        : number to code
!             = 1 : IMPMC, = 2 : MEUT2D, = 3 : SOLDOR

! local variables
      character, pointer :: cpath*10
      logical,   pointer :: kset
      integer,   pointer :: ntim
      real(8),   pointer :: dtmtbl(:), timb, timtbl(:)
 
      character :: ctime*20 = ' '
      integer   :: i, judgement

! execute only for master PE
      if( lmype.ne.lmspe ) return

! initial set
      select case( kc )
      case( 1 )
        cpath  => cpathi
        dtmtbl => dtmdbgimp
        kset   => kseti
        ntim   => ntimi
        timb   => timbi
        timtbl => timdbgimp
      case( 2 )
        cpath  => cpathn
        dtmtbl => dtmdbgneu
        kset   => ksetn
        ntim   => ntimn
        timb   => timbn
        timtbl => timdbgneu
      case( 3 )
        cpath  => cpaths
        dtmtbl => dtmdbgsol
        kset   => ksets
        ntim   => ntims
        timb   => timbs
        timtbl => timdbgsol
      end select

      if( .not. kset ) then
! set ntim
        ntim = mtim
        do i = 1, mtim-1
          if( timtbl(i) >= timtbl(i+1) ) then
            if( i > 1 ) then
              ntim = i
            else
              ntim = 0
            endif
            exit
          endif
        enddo

! create output directory
        if( ntim > 0 ) call system( 'mkdir -p ' // cpath )

        kset = .true.
      endif

! no output
      if( ntim <= 0 ) return

! output judgement
      call timcntl( timtbl, dtmtbl, ntim, time, timb, judgement )

! output
      if( judgement == 1 ) then
! set ctime
        call rn2a( time, -iwflg, len( ctime ), ctime )
! output data
        call figdat( cpath, ctime )
        if(kc==3) call pst_qdpl( cpath, ctime ) ! for soldor only
! update timb
        timb = time
      endif

      return
      end
      end module dbg_mod
