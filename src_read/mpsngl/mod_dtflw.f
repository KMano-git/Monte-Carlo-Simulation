      module mod_dtflw
      use mpi
      use mod_sizeDef
      implicit none

! deleted 1 line organize local variables and include files by kamata 2021/05/31, yamamoto 22/10/07 
!ik   character(lnnam) :: ctg_DTFLW
      integer :: itg_DTFLW

      integer :: hd_DTFLW
      integer :: tmgrp
      integer :: mlpmax

!::derived type data
      integer :: ntyp
! deleted 2 lines generalization of derived data type by kamata 2020/10/31
!ik   character(lnnam), dimension(ndtyp) :: typ_dnam
!ik   character(lnnam), dimension(ndtyp) :: typ_prog
      integer, dimension(ndtyp) :: typ_igrp
! deleted 1 line generalization of derived data type by kamata 2020/10/31
!ik   integer, dimension(ndtyp) :: typ_itag
      integer, dimension(ndtyp) :: typ_irst     ! = 1  restart data
! deleted 2 lines generalization of derived data type by kamata 2020/10/31
!ik   integer, dimension(ndtyp) :: typ_ierr
!ik   integer(kind=MPI_Address_Kind), dimension(ndtyp) :: typ_adr0

      real(8), dimension(0:ndgrp,ndtyp) :: typ_idgr !mlp + 0.01*mitr 3.02 25.01
      integer, dimension(0:ndgrp,ndtyp) :: typ_rcvg !(0:no use/1:recv/2:send)
      integer, dimension(0:ndnpe,ndtyp) :: typ_itag_pe

      integer :: nprc
      integer, dimension(ndprc)       :: prc_istp, prc_igrp, prc_iord
      integer, dimension(ndprc)       :: prc_itrn
      integer, dimension(ndtyp,ndprc) :: prc_dtrw
      integer, dimension(5)           :: prc_lvst, prc_lven
      character(120), dimension(ndprc) :: prc_scmd

      integer :: mlp        ! loop number in wf_maindrv (kmlp)
      integer :: mitr       ! iteration number (1/2/...)
      real(8) :: midno      ! ID number of typ-data  mlp + 0.01*mitr
      integer :: mstp       ! step number
!      integer :: idnum      ! ID number of typ-data
! deleted 2 lines organize local variables and include files by kamata 2021/05/31, yamamoto 22/10/07
!ik   integer :: ldemo = 0  ! demonstration of data flow
!ik   integer :: dbglv = 2  ! 0 (off) 1 (summary) 2 (all inf.)

!::control of restart  set value in rstfile(0)...merged from Shimizusan branch
!:: SYamoto
      integer       :: rst_nxt  ! when mlp = rst_nxt write rst-data
      integer       :: rst_mlp  ! set flag = 1 at every mlp
      character(30) :: rst_drW  ! write rest-data in  "../wxdr_09/PLS_2
      character(30) :: rst_drR  ! read  rest-data in  "../wxdr_05/PLS_2

!::data type (ityp <= ntyp <= ndtyp)
!    ierr: =0 normally defined dtype
!    adr0: top address of dtype data

!  ityp  typ_dnam   typ_igrp   typ_itag  typ_ierr  typ_adr0  typ_idgr(0,it)   typ_rcvg(ig,it)
!    1    "TOK"        1         21         0/1     324309      621             2,1,1
!
!::ID number of derived data
! ==> typ_idno(it) => typ_idgr(0,it)
!     typ_idgr(0,it)  : Lp number of data when it was created.
!     typ_idgr(ig,it) : Lp number of data when it was sent to ig-group.
!
!::process to make data   "exec  soldor   NTL IMP => TOK"
!   prc_istp  prc_igrp    prc_iord
!    4 (exec)  1 (soldor)  2 (ntl:1,imp:1) => TOK
!
!                                 TOK   NTL   IMP   FLX
!  "11 exec neut2d TOK => IMP"
!   prc_dtrw(:,11)             =   1     2     0     0
!  "12 exec impmc  TOK NTL => IMP FLX
!   prc_dtrw(:,12)             =   1     1     2     2
!                   1(read)  2(write)  0(No use)
!
!     del     prc_dtty(i,iprc),  prc_dtno(iprc)
!     change  prc_dtrw(i,iprc)
!
!               init  parm  prof  exec  term
!    prc_lvst =  1,    4,    7,    10,   13
!    prc_lven =  3,    6,    9,    12,   15    sol/ntl/imp
!
!    prc_itrn : iteration number in a step
!    prc_scmd : single command
!
!---------------------------------------------------------------------------

      contains

!----------------------------------------------------------------
      subroutine dtypdef_DTFLW(kerr)
!----------------------------------------------------------------
      use mpi
      use mod_mpicomm

      implicit none
      integer, intent(out) :: kerr
      integer :: ierr

      call dtyp_init( "DTFLW", "master" )

      call dtyp_addv( hd_DTFLW )
      call dtyp_addv( tmgrp )
      call dtyp_addv( mlpmax )

      call dtyp_addv( ntyp )
      call dtyp_addv( typ_dnam )
      call dtyp_addv( typ_igrp )
      call dtyp_addv( typ_itag )
      call dtyp_addv( typ_irst )
      call dtyp_addv( typ_idgr )
      call dtyp_addv( typ_rcvg )

      call dtyp_addv( nprc )
      call dtyp_addv( prc_istp )
      call dtyp_addv( prc_igrp )
      call dtyp_addv( prc_iord )
      call dtyp_addv( prc_itrn )
      call dtyp_addv( prc_dtrw )
      call dtyp_addv( prc_lvst )
      call dtyp_addv( prc_lven )
      call dtyp_addv( prc_scmd )

      call dtyp_term( ierr )

      kerr = ierr
      return
      end subroutine dtypdef_DTFLW

!----------------------------------------------------------------
      subroutine dtyplst_DTFLW
!----------------------------------------------------------------
      use mod_keylist, only : cstep
      use mod_mpicomm, only : grpprg, m6, ngrp
      use mpi

      implicit none
      character(120) :: cmsg
      integer :: ityp, i, iprc, igrp
      integer :: idtr(10), idtw(10), ndtr, ndtw
      character(80) :: cdtr, cdtw

      write(m6,'(/2x,
     >  "*** dtyplst_DTFLW ***  WFL/mod_dyflw.f")')
!       master make no derived data(igrp=0)
      write(m6,'(2x,"   tmgrp = ",i1,2x,a8)') tmgrp, grpprg(tmgrp)
      write(m6,'(2x,"   mlpmax = ",i6)') mlpmax
   
      write(m6,'(2x,"   ntyp/ndtyp = ",2i5)') ntyp, ndtyp
      write(m6,'(2x,"   typ_dnam = ",7(a8,2x)/(16x,7(a8,2x)))')
     >    typ_dnam(1:ntyp)
      write(m6,'(2x,"   prog(tx) = ",7(a8,2x)/(16x,7(a8,2x)))')
     >              (grpprg(typ_igrp(ityp)),ityp=1,ntyp)
      write(m6,'(2x,"   typ_irst = ",7(i3,7x)/(16x,7(i3,7x)))')
     >              (typ_irst(ityp),ityp=1,ntyp)
      write(m6,'(2x)')
      write(m6,'(2x,"   typ_rcvg   ",7(a8,2x))') (grpprg(igrp),
     >  igrp=1,ngrp-1)
! modified 2/1 lines organize local variables and include files by kamata 2021/05/31, yamamoto 22/10/07
!ik   do i = 1, ntyp
!ik   ityp = i
      do ityp = 1, ntyp
      write(m6,'(2x,"     ",a,10(i4,6x))') typ_dnam(ityp),
     >   (typ_rcvg(igrp,ityp),igrp=1,ngrp-1)
      enddo

      write(m6,'(/2x,"dtflw   nprc/ndprc =",2i5)') nprc, ndprc
! modified 2/1 lines organize local variables and include files by kamata 2021/05/31, yamamoto 22/10/07
!ik   do i = 1, nprc
!ik   iprc = i
      do iprc = 1, nprc
      
      call dtrwprc(iprc,idtr,idtw,ndtr,ndtw,cdtr,cdtw)
      cmsg = "  "
      if( len_trim(cdtr)+len_trim(cdtw) > 0 ) then
      cmsg = trim(cdtr)//"  =>  "//trim(cdtw)
      endif
      if( cmsg(1:2) == "  " ) cmsg(1:) = cmsg(3:)
      write(m6,'(2x,i3,2x,a4,3x,i1,2x,a,2x,i3,"  :  ",a)') iprc,
     >   cstep(prc_istp(iprc)), prc_iord(iprc),
     >   grpprg(prc_igrp(iprc)), prc_itrn(iprc), trim(cmsg)
      enddo

      write(m6,'(/2x,"process : ",5(a,1x,2i3,2x))') (trim(cstep(i)),
     >  prc_lvst(i),prc_lven(i),i=1,5)

      write(m6,'(/2x,"prc_scmd")')
      do i = 1, nprc
        write(m6,'(2x,i3,2x,"[",a,"]")') i, trim(prc_scmd(i))
      enddo

      return

      call wf_exit("mod_dtflw","test")
      end subroutine dtyplst_DTFLW

      end module mod_dtflw
