! added replace all include files with module files by kamata 2021/08/18
      module cmolhy
      implicit none

      integer, parameter :: ndtmm = 50, ndrcm = 6

!::[MPI_Bcast in molehy]   cmolhy (stemin,jt3e)  10/04/21
      real(8) :: stemin = 0.0_8, stemax = 0.0_8, stedlt = 0.0_8
     >    , stem(ndtmm) = 0.0_8, sgvm(ndtmm,ndrcm) = 0.0_8
     >    , elsm(ndrcm) = 0.0_8, engm(ndrcm) = 0.0_8
     >    , te3l(ndtmm) = 0.0_8, els3(ndtmm) = 0.0_8
     >    , te3e(ndtmm) = 0.0_8, eng3(ndtmm) = 0.0_8
      integer :: jtmx = 0, nrc = 0, jt3l = 0, jt3e = 0

!::[MPI_Bcast in molehy]  cmolhy (cract,cpid)  10/04/21
      character :: cract(ndrcm)*40 = ' ', cpid(ndrcm)*10 = ' '

      integer :: lmolhcr = 0

! added 2 lines replace all include files with module files by kamata 2021/08/18
! move from sv_molhcr ( molhcr.f )
      real(8) :: bij(0:8,0:8,1:6) = 0.0_8

      end module cmolhy
