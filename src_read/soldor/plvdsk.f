!***********************************************************************
      subroutine plvdsk(nft,kk,cft)
!***********************************************************************
      use cplcom, only : vcs, vea, vna, vne, vnezef, vni, vva, vve, vte
     >    , vti, vzf 
      use csize,  only : ndsp, ndx, ndy
      use cunit,  only : lmspe, lmype, n6
      implicit none
!
!::argument
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   integer :: nft, kk
!ik   character(*) :: cft
      integer,      intent(in) :: nft, kk
      character(*), intent(in) :: cft
!
!::local variables
      integer :: nszx, nszy, nszs
!
      if( nft.le.0 .or. len_trim(cft).le.0 ) then
      call wexit("plvdsk","invalid data of nft or cft")
      endif
!
      if( lmype.ne.lmspe ) then
        call wexit("plvdsk","call lmype .ne. lmspe")
      endif
!
!-----------------------------------------------------------------------
!::write
!-----------------------------------------------------------------------
      if( kk.eq.1 ) then
      open(unit=nft, file=cft, form="unformatted", action="write")
      write(n6,'(2x,"   plvdsk   kk =",i2,"  cft = ",a)')
     >  kk, trim(cft)
!
      nszx = ndx
      nszy = ndy
      nszs = ndsp
!
      write(nft) nszx, nszy, nszs
      write(nft) vna, vne, vni, vnezef, vzf,
     >           vva, vve, vti, vte, vcs, vea
      close(nft)
!
      return
      endif
!
!-----------------------------------------------------------------------
!::read
!-----------------------------------------------------------------------
      if( kk.eq.2 ) then
      open(unit=nft, file=cft, form="unformatted", action="read")
      write(n6,'(2x,"   plvdsk   kk =",i2,"  cft = ",a)')
     >  kk, trim(cft)
!
      read(nft) nszx, nszy, nszs
      write(n6,'(2x,"  ndx =",2i5,"  ndy =",2i5, "ndsp =",2i3)')
     >   ndx, nszx, ndy, nszy, ndsp, nszs
!
      read(nft) vna, vne, vni, vnezef, vzf,
     >          vva, vve, vti, vte, vcs, vea
      close(nft)
      return
      endif
!
      end
