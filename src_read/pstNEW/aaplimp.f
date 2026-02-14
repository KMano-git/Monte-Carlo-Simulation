!***********************************************************************
!
!   New impmc code
!
!                      2005/02/14  K. Shimizu
!                      2012/09/19  K. Shimizu
!                      2016/05/30  K. Hohshino 
!                               (merge work by Advance Soft)
!
!***********************************************************************
      program post_impmc
!***********************************************************************
      use mod_shexe, only: use_namelist, INP_PLS, set_shexe_variables
      use csize
      use cimcom
      use cimden
      use cntcom
      use cntsrc
      use cplcom
      use cplimp
      use cplmet
      use cplwrd
      use cimctl
      use cimi60
      use csonic ! limp
      use cunit ! cprg
      use catcom, only : catmz
      implicit none
!
!::local variables
      integer  ntyp, i, nv, nft, ift, it, iskp
      integer  kcolr, ic, lenx, mji
      character dsn*80, cinp*80, dsn2*80, drwk*20
      character(120) :: cdsn
      logical lex
      integer   ir2
      real*8, dimension(10) :: hwrgn
      real*8 ::  hwrsm
      integer :: mj
      character(80) :: drimp, comd
      integer :: i6 = 180000 + 21
!         plot_mode: valid value is one of 'plot_2d' and 'sight_line'
      integer :: plot_mode
      integer, parameter :: plot_2d=1,  plot_sight_line=2

      namelist /uplflg/ lfopt,lfdbg,lpst_bprf

!::directory  ../wxdr_17
      call stjob
      drwk = trim(wh_dir)
      write(n6,'(4x,"drexe  = ",a)') trim(drwk)

      call set_shexe_variables
! allocate module valuales
      cprg = 'plimp'
! set the size required for allocate
      call alc_size
      nzsmx = 10 ! expect IMPMC speces less than 10
! allocate module variables
      call gvar_alc ( cprg, 1 )

!::mpifile
      stcpy = "delete"
      n5 = 10
      n6 = 6
      ufl6 = "@@@p"
      open( unit=n6, file=ufl6 )
      call git_mesg(n6)
!
!::dbg_absimp
      open(unit=i6,file="@impabs")
!
!::start
      nope = 0
      mype = 0
      mspe = 0
      lnope = 1
      lmype = 0
      lmspe = 0
      mjpe  = 1
      cmype = '0'
      mygrp = 1
      cdgr(1) = 1
      cdgr(2) = 1
      cdgr(3) = 1
!
      kpcn = 1
!
!::define i60  edit inc/cimi60
      i60_imcndv = n6
      i60_imhist = n6
      i60_tmwgtn = n6
!
!::header
      call gdhedr
!
!::initial set
      call opinpt
      call opnfrt   !  set cdrout
      call phcnst
      lpost = 1     !  post process
!
!::plasma
      call pls_ini
      cdsn = trim(drwk) // "/DTPLV"      
      call plfgdt(cdsn,"r")
!
!::neutral
      call ntl_ini
      call mcplas
      call ntplas
      nftnr = 21
      cdsn = trim(drwk) // "/DTNTL"
      call plfgdt(cdsn,"r")
! add 1 line
      call plzset           !  mass, charge of ion and impurity
      call set_flx(1)
      call set_flx(2)
      call set_ntl
!
!
!::impmc
      wfac = 0.0d0
      cdirZ = "."
      nprg = 1 ! for avoid error in imp_pls. dfcf is not used in plimp, OK just nprg=1.
      call imp_ini_x(nprg) ! OK for kis=1. We call iminpt for another IMPMC later.
      call ntl_pls
      call imp_pls
! old file name
      cdsn = trim(drwk) // "/DTIMP"
      inquire(file=cdsn,exist=lex)
      if( lex ) then
        call plfgdt(cdsn,"r")
        call set_plimp_1(1)
        call set_plimp_2(1)
        wmc_nty = 1
        limp = 3
      endif

! new file name
        do i = 1, nzsmx
          write (cdsn,'(i0)') i
          cdsn = trim(drwk) // "/DTIMP" // trim(cdsn)
          inquire(file=cdsn,exist=lex)
          if( lex ) then
            limp = 3
            if(i>1) then
              nft = 21
              close(nft)
            endif
            call plfgdt(cdsn,"r")
            if(i>1) ismax = nzmx
            call set_plimp_1(i)
            call set_plimp_2(i)
            wmc_nty = i
          endif
        enddo

! new file name: Yamoto
        do i = 1, nzsmx
          write (cdsn,'(i0)') i
          cdsn = trim(drwk) // "/IMPSY" // trim(cdsn)
          inquire(file=cdsn,exist=lex)
          if( lex ) then
            limp = 3
            if(i>1) then
              nft = 21
              close(nft)
            endif
            call plfgdt(cdsn,"z")
            ismax = nzmx
            call set_plimp_1y(i)
            call set_plimp_2y(i)
            wmc_nty = i
          endif
        enddo
      nzsmx = wmc_nty
!
!
!::radiation
      cdsn = trim(drwk) // "/DTWRD"
      call plfgdt(cdsn,"r")
!
!::graphic mode
      call gdhedr
      call gdtrm("X")
      call gdinit(1)
!
!::normalization   (use gdget)
      call set_imp
!
!::text data
!     ::fig size
      call listmc
      nft = 21
      if(use_namelist) then
         call nopen(nft,INP_PLS,"text",cdsn,lex)
      else
         call nopen(nft,"INP_PLS","text",cdsn,lex)
      end if
      read(nft,uplflg)
      close(nft)
      call pst_bprf
!
!::profile along poloidal length & radial profile
      open(unit=9999,file="inppls",form="formatted")
      if( wmc_nty > 0) then
        do i = 1, wmc_nty
          call iminpt_x(9999,i) !catcom
        enddo
      endif
      close(9999)

!::IMPMC results extracted from plcrad_mcN SYamoto
      if(limp.eq.3) then
        do i = 1, wmc_nty
          call set_adas_plimp(i) !catcom
          write(n6,'(/2x, "ADAS data set for ",i3,a)')
     >       i, catmz(i)
          call atm_ini(i)
        enddo
      endif
      call plwrdMC ! nty=1 wmc_zwi, wmc_zwe

      call impprof
      call imprfrd
!
!
      plot_mode=plot_2d

      call gdget(">> next step (2:plot2d<rtn>/s:sight line/q:quit) ==> "
     >                                                     ,cinp,mji)
      if( mji.gt.0 .and. cinp(1:1).eq."q" ) goto 910
      if( mji.gt.0 .and. cinp(1:1).eq."s" ) plot_mode=plot_sight_line
!
      if(plot_mode==plot_2d) then
!::data
          kcolr = 1   !  colour
          call gdget(">> enter kcolr (0:bw/1:color<rtn>) ==> ",cinp,mji)
          if( mji.gt.0 .and. cinp(1:1).eq."0" ) then      
          kcolr = 0   !  black/white
          endif
          call gdprint(kcolr)
      endif
!
      if(plot_mode==plot_2d) then
        do 
!         ::fig size
          nft = 21
          call nopen(nft,"inpfig","text",dsn,lex)
          call gpsta(nft)
          close(nft)
          call gtfgmx
!
          call plotp2d(plot_mode)
        enddo
      else
        call plotp2d(plot_mode)
      endif
!
!::graphic mode
 910  continue
      call gdmod(0)
!
      stop
      end
