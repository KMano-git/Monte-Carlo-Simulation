! added replace all include files with module files by kamata 2021/08/18
      module cplwrd
      implicit none

      real(8) :: wtime = 0.0_8, witim = 0.0_8, wnlp = 0.0_8
     >    , wfac = 0.0_8

      real(8), allocatable :: wcre(:)  ! CR radiation
     >             , wmci(:), wmce(:)  ! MC radiation
     >             , wimi(:), wime(:)  ! total radiation
     >             , wcr_wrd(:,:)

!::2015/06/10  K.Shimizu  multi non-corona model
      integer, parameter :: ndzty = 3
      integer      :: wcr_nty = 0               ! number of impurity kind  
      integer      :: wcr_ity(ndzty) = 0        !  3 (W)  2 (Ar)
      character(4) :: wcr_typ(ndzty) = ' '      ! "C", "Ar", "W", "W2"
      real(8)      :: wcr_cnc(10,ndzty) = 0.0_8 ! concentration in each region

      integer, parameter :: ndprg = 8
      character(4) :: prg_typ(ndprg) = 
     >           (/"C ", "W ", "W2", "Be", "Ar", "Ne", "Kr", "Xe"/)

      end module cplwrd
