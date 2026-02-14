      module com_msdat
      implicit none

      integer,parameter :: ndq = 320, ndp = 120
      integer,parameter :: ndfx = 129, ndfy = 129
      real*4 :: fx(ndfx) = 0.0_4, fy(ndfy) = 0.0_4
     > , fz(ndfx*ndfy) = 0.0_4

      integer :: mqd1 = 0, mqx1 = 0, mqs1 = 0, mqs2 = 0
     > , mqh1 = 0, mqh2 = 0, mqx2 = 0, mqd2 = 0, nqmx = 0

      integer :: mpw1 = 0, mpsp = 0, mpw2 = 0, mpax = 0, npmx = 0
      real*8 :: pssol(ndp) = 0.0_8, psprv(ndp) = 0.0_8
     > , psman(ndp) = 0.0_8
      real*8 :: grdx(ndq,ndp) = 0.0_8, grdy(ndq,ndp) = 0.0_8

      real*8 :: hbr(ndq,ndp) = 0.0_8
     > , hbz(ndq,ndp) = 0.0_8, hbt(ndq,ndp) = 0.0_8

      end module com_msdat