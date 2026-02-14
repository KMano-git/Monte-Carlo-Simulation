! added replace cplcom from include file to module file by kamata 2021/08/18
! removed ndwf because it is unused by kamata 2022/05/29
      module csize
      implicit none
!----------------------------------------------------------------------
! size    size_68/csize   ! ITER  NCT
!----------------------------------------------------------------------
! plasma
      integer, parameter :: ndsp = 5, ndeq = ndsp*2 + 2
      integer :: ndx = 160, ndxy, ndy = 80

! neutral/impmc
      integer, parameter :: ndms = 4, ndgs = 1, ndwl = 30, ndwh = 40
      integer, parameter :: nvxd = 100, ndgt = 5
     >    , ndad =   4, ndpm =  5
      integer :: ndbr, ndmc = 6500, ndmg = 6800, ndwp, nvyd = 20
      integer :: ndgtp = 100

! No use in gmesh
! neutral scoreing variables (see cntwcn)
      integer nwkmp_nt

! neutral (flux & source)
      integer, parameter :: ndmfl = 8, ndmsr = 6, ndptl = 2000
  
! nzsmx : maximum kind of impurity (C + Ar + etc)
! ndmis : impurity atomic number
!  C,  Ne,  Ar,  Kr,  Xe,   W
!  6,  10,  18,  36,  54,  74
! mdsp : collision frequency including impurity
      integer, parameter :: ndmis = 74
      integer :: mdsp, nzsmx = 3

      end module csize
  ! Note :
  ! csize
  ! vnam    68    79    79B   120    400
  ! ndx    160   160   160    160    320
  ! ndy     80    80    80    150    120
  ! ndxy   160   160   160    160    320
  ! ndmg  6800  7900  7900  12000  40000
  ! ndmc  6500  7600  7600  11000  38400
  ! nvyd    20    20    50     20     20
  ! ndwp  2*(ndx+ndy),600,600,2*(ndx+ndy),2*(ndx+ndy)

  ! cgdcom
  ! ndq    160   160   160    160    320
  ! ndp     80    80    80    150    120
  
  ! cplcom_alc
  ! ndeq, ndgs, ndmc, ndmfl, ndsp, ndwl, ndwp, ndx, ndy
  ! ndq = ndx, ndp = ndy
  