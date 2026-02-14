!**********************************************************************
      subroutine dtrwprc(iprc,idtr,idtw,ndtr,ndtw,cdtr,cdtw)
!**********************************************************************
!    dtrwprc  dt(data) rw(read/write) prc(for a prcoess)
!
!    11 exec  soldor    NTL IMP => TOK
!
!                     TOK  NTL  IMP  FLX
!     prc_dtrw(:11)    2    1    1    0
!
!       iprc : (i) 11
!       idtr : (o) 2, 3       idtw : (o) 1
!       ndtr : (o) 2          ndtw : (o) 1
!       cdtr : (o) "NTL IMP"  cdtw : (o) "TOK"
!
!----------------------------------------------------------------------
      use mod_dtflw,   only : nprc, ntyp, prc_dtrw, typ_dnam
      use mod_mpicomm, only : m6
      use mpi

      implicit none
! modified 2/3 lines organize local variables and include files by kamata 2021/05/31, yamamoto 22/10/07
!ik   integer :: iprc, idtr(10), idtw(10), ndtr, ndtw
!ik   character(80) :: cdtr, cdtw
      integer, intent(in)        :: iprc
      integer, intent(out)       :: idtr(10), idtw(10), ndtr, ndtw
      character(80), intent(out) :: cdtr, cdtw

! modified 2/1 lines organize local variables and include files by kamata 2021/05/31, yamamoto 22/10/07
!ik   integer :: iir, iiw, i, ity, irw, mj
!ik   character(80) :: ctmp
      integer :: iir, iiw, ity, irw

      if( iprc < 1 .or. iprc > nprc ) then
      write(m6,'(/2x,"*** Stop at sub. dtrwprc  wrong iprc <= nprc",
     > 2x,2i4)') iprc, nprc
      call wf_exit("dtrwprc","wrong iprc")
      endif

!::clear
      ndtr = 0
      ndtw = 0
      idtr(1:10) = 0
      idtw(1:10) = 0
      cdtr = " "
      cdtw = " "

      iir = 0
      iiw = 0
! modified 2/1 lines organize local variables and include files by kamata 2021/05/31, yamamoto 22/10/07
!ik   do i = 1, ntyp
!ik   ity = i
      do ity = 1, ntyp

      irw = prc_dtrw(ity,iprc)
      if( irw == 1 ) then
      iir = iir + 1
      if( iir > 10 ) goto 920
      idtr(iir) = ity
      cdtr = trim(cdtr)//"  "// trim(typ_dnam(ity))
      endif
      if( irw == 2 ) then
      iiw = iiw + 1
      if( iiw > 10 ) goto 940
      idtw(iiw) = ity
      cdtw = trim(cdtw)//"  "//trim(typ_dnam(ity))
      endif
      enddo

      cdtr(1:) = cdtr(3:)
      cdtw(1:) = cdtw(3:)
      ndtr = iir
      ndtw = iiw
    
      return

!::error
  920 continue
      write(m6,'(/2x,"*** Stop at sub dtrwprc  iir > 10",
     >  3x,2i4)') iir, 10
      call wf_exit("dtrwprc","iir > 10")

  940 continue
      write(m6,'(/2x,"*** Stop at sub dtrwprc  iiw > 10",
     >  3x,2i4)') iiw, 10
      call wf_exit("dtrwprc","iiw > 10")

      end subroutine dtrwprc
