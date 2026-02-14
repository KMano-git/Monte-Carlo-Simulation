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
!)    order  code          IN-dtyp       => OUT-dtyp
!)   ----------------------------------------------------------
!)      1    soldor  (1)                 =>  PLS_2B, PLS_2
!)      2    neut2d  (1)   PLS_2B        =>  NTL_2
!)      3    soldor  (2)   NTL_2         =>  PLS_2B, PLS_2
!)      4    IMPMC   (1)   PLS_2         =>  IMP_2
!)      5    soldor  (3)   NTL_2, IMP_2  =>  PLS_2B, PLS_2
!)
!::SOLDOR
!)    PLS_2, PLS_2B contain flps calculated in sub. plpflx
!)    plmrgn Not depend on PLS_2(B)
!)    plwtfy depend on q1a
!)    plpflx depend on plwtfy and send to neut2d IMPMC
!)    plmflx tflne tflqi tflqe for output
!)
!)
!) explanation of call pldens in soldor(2)
!)   neut2d   Ni, Ti, Te  (PLS)      => (Si^, Wcx^)  (NTL_2)
!)   soldor   cal. snflx, (Si, Wcx)  =>  N0, T0      (PLS_2)
!)   IMPMC    N0, T0
!)
!----------------------------------------------------------------------
      use cntctl,      only : kntl
      use cntmnt,      only : dotn, sdotn
      use cplcom,      only : dlpq, dtimb, dtmax, dttb, edtb, edtmx
     >    , edtq, elpmx, eltb, itcn, ittb, lttb, snflx
      use cplwrd,      only : wfac
      use cputim,      only : eltm0, eltm00
      use csonic,      only : dtim, itim, kpcn, limp
      use cunit,       only : n6
      use mod_dtflw,   only : rst_drr
      use mod_keylist, only : knorm
      use mod_mpicomm, only : m6
      use mpi!,         only : mpi_wtime
      implicit none

!::argument
! modified 1/2 lines organize local variables and include files by kamata 2021/05/31
!ik   integer :: kitr, kcnd
      integer, intent(in)  :: kitr
      integer, intent(out) :: kcnd

!::local variables
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   character(12) :: chnum
      character(80) :: cmsg
      real(8) :: tfac
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer :: ierr, k, iend
      integer :: iend

      select case(kitr)
!--------------------------------------------------------------------
      case(1,2,3,4,5)
!--------------------------------------------------------------------
!::option => simple
      kntl = 1
!      kimp = 1
      kpcn = 1

      write(m6,'(2x,"prg_prof soldor STA",i4)')kitr
      write(n6,'(2x,"prg_prof soldor STA",i4)')kitr

      if( limp .ne. 3 .and. kitr > 2) then
        write(cmsg,'("kitr > 2 in no IMPMC case")')
        call wexit("prg_prof(soldor)",trim(cmsg))
      endif

        call flush(m6)
        call plmrgn
        call plwtfy  !  weit for y-flux (before plpflx)
        call plpflx  !  flux for neut2d
        call plmflx  !  flux at the plasma edge
      write(m6,'(2x,"prg_prof soldor END")')
      call flush(m6)

      write(n6,'(2x,"sonic_old   soldor finished")')
      call flush(n6)

!::stedy cal.
!kh      if( lstdy.ne.0 ) then
!kh        eltm5 = MPI_WTIME()
!kh        cputot = eltm5-eltm00
!kh        write(n6,'(/2x," total cpu = ", f10.3, " sec.",f6.2," hour")')
!kh     >    cputot,  cputot/3600.0
!kh      call wexit("prg_prof","lstdy == 1")
!kh      endif

        call plnflx(snflx)      ! Note. snflx defined in pwstep
        call pldens("prg_mcini")   ! (snflx, N^0, T^0) => N0, T0
        call dbg_snflx("prg_prof-soldor(2)")

        dotn = 0.0d0
        sdotn = 0.0d0
        call plsorc(0)
        call plcrad_mc
        call plauxv
!-------
        wfac = 0.0d0
        write(n6,'(2x,"prg_prof  wfac = 0.0")')
!-------
        call plfimp
        call plbman
        call plmflx

!::KH170202 particle flux to wall for C-genaration in IMPMC
      if( kitr.eq.1 ) call ntwdeg
      if( limp.eq.3 .and. kitr.gt.1 )then
          call ntwflx
          call plwflx
      endif
!
!--------------------------------------------------------------------
      case default
!--------------------------------------------------------------------
      write(cmsg,'("kitr =",i2," > 5")') kitr
      call wexit("prg_prof(soldor)",trim(cmsg))
!--------------------------------------------------------------------
      end select

!::time control
      iend=0
      if( limp .ne. 3 .and. kitr == 2) iend=1
      if( kitr == 5 ) iend=1
      if( iend == 1 ) then
        itim = 0
        dlpq = 1.0d0
        edtq = 1.0d0
        ittb  = 1
        itcn  = lttb(ittb)
        dtim  = dttb(ittb)
        dtmax = dttb(ittb)
        elpmx = eltb(ittb)
        edtmx = edtb(ittb)
        dtimb = dtim
        tfac  = dtim/dtimb
        write(n6,'(2x,"<dtcntl-so>  set  dtim =",1pe11.3,"  tfac ="
     >  ,0pf8.3,"  elpmx,edtmx =",1p2e11.3/)') dtim,tfac,elpmx,edtmx
        if( dtim.le.0.0 .or. elpmx.le.0.0 .or. edtmx.le.0.0 ) goto 290

      eltm0 = MPI_WTIME()
      write(n6,'(2x," initial cpu ",f10.3,"  sec.")') eltm0-eltm00
      endif

!::rstfile
      if( kitr == 1 .and. rst_drR(1:1) == ".") then
        call rstfile(2)
      endif

      kcnd = knorm
      return

!-----------------------------------------------------------------------
!::error
!-----------------------------------------------------------------------
 290  continue
      write(n6,'("invalid time step control  dtim,elpmx,edtmx =",
     >  1p3e12.3)') dtim,elpmx,edtmx
      call wexit( "aasonic", "invalid time step control" )
      end
