!**********************************************************************
      subroutine rstfile(kact)
!**********************************************************************
!)
!)  kact = 0: set restart dir name
!)       = 1: write restart file     rst_drW = "../wxdr_24"
!)       = 2: read  restart file     rst_drR = "../wxdr_07"
!)
!)  file name                     No use
!)     "../wxdr_07/RSTF-PLS_2" ==>   "  ../wxdr_07/RSTF-PLS_2
!)
!)RSTF    soldor/prg_main.f:  call rstfile(0)   set drW/R mlp
!)RSTF    WFL/wf_maindrv.f:   call rstfile(1)   write
!)RSTF    soldor/prg_prof.f:  call rstfile(2)   read
!)
!)  dtyp-data in system  : GL <==> GL
!)
!---------------------------------------------------------------------
      use csonic,      only : is_rsta, itim, time
      use cunit,       only : n6, wh_dir
      use mod_keylist, only : kexec
      use mod_sizedef, only : lnnam
      use mod_mpicomm, only : m6, mygrpM => mygrp, nmst_cmm, nrnk_grp
     >    , nwld_grp
      use mod_dtflw,   only : mlp, ntyp, rst_drr, rst_drw, rst_mlp
     >    , rst_nxt, typ_dnam, typ_igrp, typ_irst
      use mod_shexe,   only : use_namelist, INP_PLS
      use mpi!,         only : mpi_bcast, mpi_character
      use cplcom, only : cftr
      use cimcom, only : mrgn_max
      implicit none

!::argument
      integer, intent(in) :: kact

!::local variables
      integer :: nft
      character(80) :: cdsn
      logical :: lex
      integer :: nch, ierr

      character(80) :: cmsg
      integer :: mjs(20), mje(20), nk
      integer :: i, igrp
      character(lnnam) :: dnam
      character(80) :: dsn, rNumber_st
      character(len=9), parameter :: dummy_char = "NoReadDir"
      integer error, mj
      integer         nc
      character(2)    np
      logical is_imprst
!
! for dummy value
      integer, parameter :: rgsize=10, dpsize=2, xplmt_len=80
      real(8), dimension(rgsize) :: dtq1rg, dtq2rg, dtq3rg, dtq4rg
      real(8), dimension(dpsize)  :: dtq1dp, dtq2dp, dtq3dp, dtq4dp
      character xplmt(xplmt_len)
      integer timeNum, interNum
      real(8) :: vlda_dummy(mrgn_max)
!
      namelist /urest/ rst_drW, rst_drR, rst_mlp ! rst_drW, rst_drR in inppls is dummy.
!
      is_imprst = .true.
!
      select case(kact)

!::set dir name
      case(0)
        write(n6,'(/a,2x,"*** rstfile(0) ***  mlp =",i6)') "RSTF",
     >    mlp
        write(m6,*)
        write(m6,'(/a,2x,"*** rstfile(0) ***  mlp =",i6)') "RSTF",
     >    mlp
        write(n6,'("** rst_drw in inppls is ingnored**")')
        write(m6,'("** rst_drw in inppls is ingnored**")')
        write(n6,'("** rst_drR in inppls is ingnored**")')
        write(m6,'("** rst_drR in inppls is ingnored**")')

        if( nrnk_grp == nmst_cmm ) then
          nft = 21
          if(use_namelist) then
            call nopen(nft,INP_PLS,"text",cdsn,lex)
          else
            call nopen(nft,"INP_PLS","text",cdsn,lex)
          end if
          read(nft,urest,iostat=error) ! read rst_mlp
          call read_uplinp(nft,dtq1rg, dtq2rg, dtq3rg, dtq4rg, rgsize
     >     , dtq1dp, dtq2dp, dtq3dp, dtq4dp, dpsize
     >     , xplmt, xplmt_len, timeNum, interNum, vlda_dummy)
          close(nft)

          ! set rst_drW as current directory
          rst_drW = "../"//trim(wh_dir)

          ! set rst_drR as consistent with soldor(cftr=../dtprf_XX)
          inquire(file=trim(cftr), exist=lex)
          if(lex) then
            mj = index(cftr,"_")
            rNumber_st = trim(cftr(mj+1:))
            rst_drR = "../wxdr_"//trim(rNumber_st)
          else
            rst_drR = dummy_char
          endif
          !
          write(m6,'("  [",a,"]",2x,"[",a,"]",2x,"[",i6,"]")')
     >     trim(rst_drW), trim(rst_drR), rst_mlp
          write(cmsg,'("  [",a,"]",2x,"[",a,"]",2x,"[",i6,"]")')
     >     trim(rst_drW), trim(rst_drR), rst_mlp
        endif

        if( nrnk_grp >= 0 ) then
          nch = 80
          write(m6,*) 'RSTFILE BEFORE BCAST'
          call MPI_Bcast(cmsg,nch,MPI_CHARACTER,0,nwld_grp,ierr)
          call MPI_Bcast(cftr,nch,MPI_CHARACTER,0,nwld_grp,ierr)
          write(m6,'(a)') cmsg
          call linsep(cmsg,"[]",nk,mjs,mje,20)
          rst_drW = cmsg(mjs(1):mje(1))
          rst_drR = cmsg(mjs(2):mje(2))
          if( rst_drR == dummy_char ) rst_drR = ' '
          read(cmsg(mjs(3):mje(3)),*) rst_mlp

          rst_nxt = kexec + rst_mlp
          write(n6,'(2x,"rst_drW = [",a,"]  rst_drR = [",a,"]",
     >     "  rst_mlp =",i6,"  rst_nxt =",i4)')
     >     trim(rst_drW), trim(rst_drR), rst_mlp, rst_nxt
        else
          write(n6,'(2x,"nrnk_cmm < 0  No rstfile data")')
        endif

!::write restart data
      case(1)
        write(n6,'(/a,2x,"*** rstfile(1) ***  mlp =",i6)') "RSTF",
     >    mlp
        do i = 1, ntyp
          if( typ_irst(i) == 0 ) cycle
          igrp = typ_igrp(i)
          dnam = typ_dnam(i)
          if( mygrpM /= igrp ) cycle
          dsn = trim(rst_drW)//"/RSTF-"//trim(dnam)
          write(n6,'(2x,"RSTF",2x,"mygrp =",i2,2x,a,2x,"mlp/itim =",
     >     2i6,"  time =",1pe14.6)') mygrpM,trim(dsn),mlp,itim,time
          call dtypfile(dnam,dsn,1,is_imprst)
        enddo

!::read restart data
      case(2)
        write(n6,'(/a,2x,"*** rstfile(2) ***  mlp =",i6)') "RSTF",
     >    mlp
!
        ! only for IMPMC, read IMP(nc)_2 at first.
        if(mygrpM .ge. 3) then
          do i = 1, ntyp
            ! check group number
            igrp = typ_igrp(i)
            if( mygrpM /= igrp ) cycle
            ! check IMP...
            dnam = trim(typ_dnam(i))
            if( dnam(1:3)/="IMP" ) cycle
            ! check IMP(nc)_2
            call getncnp( dnam, nc, np )
            if(np == '2 ') then
              if( typ_irst(i) == 0 ) then
                write(n6,'(2x,"No IMP(nc)_2 in dtflw, skip RSTF")')
                is_rsta = .false.
                is_imprst = .false.
                return
              endif
              dsn = trim(rst_drR)//"/RSTF-"//trim(dnam)
              inquire(file=trim(dsn), exist=is_rsta)
              if(is_rsta) then
                ! read IMP(nc)_2 
                call dtypfile(dnam,dsn,2,is_imprst)
                if(is_imprst) then
                  exit
                else
                  write(n6,'(2x,"Invalid IMP(nc)_2, skip RSTF")')
                  is_rsta = .false.
                  return
                endif
              else ! file not found
                write(n6,'(2x,"Not found IMP(nc)_2, skip RSTF")')
                is_imprst = .false.
                return
              endif
            endif
          enddo
        endif
!
        do i = 1, ntyp
          if( typ_irst(i) == 0 )then
          write(n6,*) 'cycle due to typ_irst(i) eq 0.', i
            cycle
          end if
          igrp = typ_igrp(i)
          dnam = typ_dnam(i)
          if( mygrpM /= igrp )then
            write(n6,*) 'cycle mygrpM =', mygrpM, 'igrp =', igrp
            cycle
          end if
          dsn = trim(rst_drR)//"/RSTF-"//trim(dnam)
          inquire(file=trim(dsn), exist=is_rsta)
          if(.not. is_rsta) then
            write(n6,'(2x,"RSTF not used or not found: ",a)') trim(dsn)
            cycle
          endif
          write(n6,'(2x,"RSTF",2x,"mygrp =",i2,2x,a,2x,"mlp/itim =",
     >     2i6,"  time =",1pe14.6)') mygrpM,trim(dsn),mlp,itim,time
          call dtypfile(dnam,dsn,2,is_imprst)
        enddo

      end select

      return
      end
