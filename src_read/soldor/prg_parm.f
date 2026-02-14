!**********************************************************************
      subroutine prg_parm(kitr,kcnd)
!**********************************************************************
!)
!)  include csonic
!)  <  time, tend, dtim, itim, itend
!)
!)  us mod_loc_tim
!)     time, tend, dtim, itim, iend
!)
!)   parm 1 : soldor   =>  PLS_1                     !  kitr = 1
!)   parm 2 : neut2d  MST_1  PLS_1         => NTL_1
!)   parm 3 : IMPMC   MST_1  PLS_1  NTL_1  => IMP_1
!)   parm 4 : soldor  NTL_1  IMP_1                   !  kitr = 2
!)--------------------------------------------------------------------
      use catcom,      only : catmz
      use cplimp,      only : wmc_nty
      use csonic,      only : dtim, itend, itim, limp, lrand, lstop
     >    , tend, time
      use cunit,       only : n6
      use mod_dtflw,   only : mlpmax
      use mod_keylist, only : knorm
      use mod_loc_tim, only : dtimm, iend, itimm, tendm, timem
      use topics_mod,  only : lgtpc
      implicit none

!::argument
! modified 1/2 lines organize local variables and include files by kamata 2021/05/31
!ik   integer :: kcnd, kitr
      integer, intent(in)  :: kitr
      integer, intent(out) :: kcnd

!::local variables
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   character(12) :: chnum
      character(80) :: cmsg
! added 1 line treat 4 or more impurities with IMPMC by kamata 2022/04/15
      integer    i

      select case(kitr)
!--------------------------------------------------------------------
      case(1)
!--------------------------------------------------------------------
      call ranini(lrand)
! added 1 line integrated calculation with TOPICS by kamata 2020/11/16
      if( lgtpc == 0 ) then
        call pls_ini
! added 6 lines integrated calculation with TOPICS by kamata 2020/11/16
      else
        call pls_ini_t
        tend   = 1.0d38
        itend  = 2**31-1
        mlpmax = itend
      endif

      lstop = 0
      timeM = time
      tendM = tend
      dtimM = dtim
      itimM = itim
      iend  = itend

! modified 2/2 lines with TOPICS by kamata 2021/12/22
!ik   write(n6,'(4x,"Time in mod_loc_tim from pls_ini",2x,
!ik  >  "time/tend/dtim =",1p3e12.4,2x,"itim/itend =",i3,i6)')
      write(n6,'(4x,"Time in mod_loc_tim from pls_ini",2x,
     >  "time/tend/dtim =",1p3e12.4,2x,"itim/itend =",2i11)')
     >  timeM, tendM, dtimM, itimM, iend

!--------------------------------------------------------------------
      case(2)
!--------------------------------------------------------------------
      call plzset           !  mass, charge of ion and impurity
      if(limp.eq.3) then
        call chk_imp
!
! modified 13/4 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik     write(n6,'(/2x, "ADAS data set for ",i3,a)'), 1, catmzL(1)
!ik     call set_atcom_B1
!ik     call atm_ini
!ik     if(wmc_nty>1)then
!ik        write(n6,'(/2x, "ADAS data set for ",i3,a)'), 2, catmzL(2)
!ik        call set_atcom_B2
!ik        call atm_ini_2
!ik     endif
!ik     if(wmc_nty>2)then
!ik        write(n6,'(/2x, "ADAS data set for ",i3,a)'), 3, catmzL(3)
!ik        call set_atcom_B3
!ik        call atm_ini_3
!ik     endif
        do i = 1, wmc_nty
          write(n6,'(/2x, "ADAS data set for ",i3,a)') i, catmz(i)
          call atm_ini( i )
        enddo
      endif

!--------------------------------------------------------------------
      case default
!--------------------------------------------------------------------
      write(cmsg,'("kitr =",i2," > 2")') kitr
      call wexit("prg_init(soldor)",trim(cmsg))
      end select

      kcnd = knorm
      return
      end
