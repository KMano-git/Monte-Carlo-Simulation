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
      end module dbg_mod
