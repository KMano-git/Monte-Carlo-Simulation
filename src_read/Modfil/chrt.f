! added replace all include files with module files by kamata 2021/08/18
      module chrt
      implicit none

! hrt common
      real(8), dimension(30,30,30,30,30):: sgion = 0.0_8, sgalf = 0.0_8
     >    , sgrec = 0.0_8
      real(8), dimension(30):: ste = 0.0_8, sne = 0.0_8, st0 = 0.0_8
     >    , sn0 = 0.0_8, sng = 0.0_8
      real(8) :: stemin_hrt = 0.0_8, snemin_hrt = 0.0_8
     >    , st0min_hrt = 0.0_8, sn0min_hrt = 0.0_8, sngmin_hrt = 0.0_8
      real(8) :: stemax_hrt = 0.0_8, snemax_hrt = 0.0_8
     >    , st0max_hrt = 0.0_8, sn0max_hrt = 0.0_8, sngmax_hrt = 0.0_8
      real(8) :: facte = 0.0_8, facne = 0.0_8, fact0 = 0.0_8
     >    , facn0 = 0.0_8, facng = 0.0_8
      integer :: mxte = 0, mxne = 0, mxt0 = 0, mxn0 = 0, mxng = 0
      integer :: lfte = 0, lfne = 0, lft0 = 0, lfn0 = 0, lfng = 0

      character(80) :: cdhrt = ' '

      end module chrt
