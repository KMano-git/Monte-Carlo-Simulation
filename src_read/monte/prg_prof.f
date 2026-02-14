!**********************************************************************
      subroutine prg_prof(kitr,kcnd)
!**********************************************************************
!)
!)      sequential type    2015/07/23
!)
!)   Mnkntl_pls (vna,vti), (flps)  => NEUT2D
!)   Mnkimp_pls (vna,vti), (flps), (q1a), (N0, T0) => IMPMC
!)
!)
!)    order  code     IN-dtyp      => OUT-dtyp
!)   ----------------------------------------------------------
!)      1    soldor                => PLS_2B, PLS_2
!)      2    neut2d   PLS_2B       => NTL_2
!)      3    IMPMC    PLS_2        => IMP_2
!)      4    soldor   NTL_2, IMP_2 => PLS_2B, PLS_2
!)
!::SOLDOR
!)    PLS_2, PLS_2B contain flps calculated in sub. plpflx
!)    plmrgn Not depend on PLS_2(B)
!)    plwtfy depend on q1a
!)    plpflx depend on plwtfy and send to neut2d IMPMC
!)    plmflx tflne tflqi tflqe for output
!)
!----------------------------------------------------------------------
      use cimctl,      only : kimp, lstdy
      use cntctl,      only : kntl, nini_cal
      use cputim,      only : cputot, eltm0, eltm00, eltm5
      use csonic,      only : kpcn, limp
      use cunit,       only : n6
      use mod_keylist, only : knorm
      use mod_mpicomm, only : m6
      use mpi!,         only : mpi_wtime
      implicit none

!::argument
      integer, intent(in)  :: kitr
      integer, intent(out) :: kcnd

!::local variables
      character(80) :: cmsg
      integer :: k, iend

      iend = 0

      select case(kitr)
!--------------------------------------------------------------------
      case(1,2)
!--------------------------------------------------------------------
!::option => simple
      kntl = 1
      kimp = 1
      kpcn = 1

        write(m6,'(2x,"prg_prof neut2d STA")')
        call flush(m6)
        call ntl_pls
        do k = 1, 10
          call ntl_cal
          if( kitr == 1 ) exit
          if( k == nini_cal ) exit
        enddo
        write(m6,'(2x,"prg_prof neut2d END")')
        call flush(m6)

!::stedy cal.
      if( lstdy.ne.0 ) then
        eltm5 = MPI_WTIME()
        cputot = eltm5-eltm00
        write(n6,'(/2x," total cpu = ", f10.3, " sec.",f6.2," hour")')
     >    cputot,  cputot/3600.0
      call wexit("prg_prof","lstdy == 1")
      stop
      endif


!--------------------------------------------------------------------
      case default
!--------------------------------------------------------------------
      write(cmsg,'("kitr =",i2," > 2")') kitr
      call wexit("prg_prof(monte)",trim(cmsg))
      end select

      if( limp.ne.3 .and. kitr ==1 ) iend = 1
      if( kitr == 2 ) iend = 1
      if( iend == 1 ) then
        eltm0 = MPI_WTIME()
        write(n6,'(2x," initial cpu ",f10.3,"  sec.")') eltm0-eltm00
      endif

      kcnd = knorm
      return
      end
