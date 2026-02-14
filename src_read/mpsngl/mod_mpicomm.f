      module mod_mpicomm     ! mpi communication
      use mod_sizeDef
      implicit none

!::program name
      character(lnnam) :: grpprg(0:ndgrp)
!
!::groups
      integer :: ngrp, mygrp
!

!::
!::whole world  _all
!::group world  _grp
!::commnucation world _cmm
!
      integer :: nope_all, nrnk_all, nmst_all, nwld_all
      integer :: nope_grp, nrnk_grp, nmst_grp, nwld_grp
      integer :: nope_cmm, nrnk_cmm, nmst_cmm, nwld_cmm
!
!:: program name & IO
      character(16):: myfile
      character(6) :: fmtpe
      integer      :: n5=15, n6=16, m6=17, iw6=18, port_inppls=111

!
!:: cpu time
      real*8:: cpu0,cpu1
!

      interface

      subroutine wf_exit(csub,cmsg)
      character(len=*),intent(in):: csub,cmsg
      end subroutine wf_exit

      subroutine wf_comm(pnam)
      character,intent(in):: pnam*(*)
      end subroutine wf_comm

      subroutine wf_term(cmsg)
      character(*),intent(in) :: cmsg
      end subroutine wf_term

      end interface

      end module mod_mpicomm
