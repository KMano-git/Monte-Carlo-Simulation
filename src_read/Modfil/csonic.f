! added replace all include files with module files by kamata 2021/08/18
      module csonic
      implicit none
!---------------------------------------------------------------------
!::modelling & flag
!---------------------------------------------------------------------
!
!   move (itim,itend,lstop)     in csonic-/csonic/
!   move (time,tend,dtim,lcstg) in cplcom-/cmtime/
!                                            ==> csonic-/cscntl/
!
!---------------------------------------------------------------------
      integer :: lstep = 0, lcode = 0
      integer :: ndsk = 0, nhst = 0, nhsav = 0, nprf = 0, mxcpu = 0
     >    , mxdsk = 0, mxhst = 0, mxprf = 0, lfopt(100) = 0
     >    , lfdbg(100) = 0, lpst_bprf(20) = 0, nfopt = 0, nfdbg = 0
     >    , lpls = 0, npls = 0, kpls = 0, lpcn = 0, npcn = 0, kpcn = 0
     >    , lpchk = 0, kdsk = 0, khst = 0, lntl = 0, limp = 0, lmstd = 0
     >    , lqick = 0, lpost = 0, lrand = 0

!::[MPI_Bcast in lnkall_ini]  csonic/csonic/ (lstep,lrand)  10/04/21
!::[MPI_Bcast in opinpt]  csonic/csonic/ (lstep,lrand)  10/04/21

      logical :: is_rsta = .false.

!::[MPI_Bcast in lnkall_ini]  csonic/cscntl/ (time,emrk)    10/04/21
!::[MPI_Bcast in lnkall_tim]  csonic/cscntl/ (time,emrk)    10/04/21

! modified 1/1 lines undefined use by kamata 2021/08/18
!ik   real(8)    time, tend, dtim
      real(8) :: time = 0.0d0, tend = 0.0_8, dtim = 0.0_8
      integer :: itim = 0, itend = 0, mstop = 0, lstop = 0, lcstg = 0

      end module csonic
