!**********************************************************************
      subroutine pllist(nft)
!**********************************************************************
      use cplmet, only : itpvs, itsle
      use csonic, only : itim, time
      implicit none
!
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  nft
      integer, intent(in) :: nft
!
!::local variables
      integer  ia, nmax, n, it
!
      if( nft.le.0 ) return
!
      ia = 1
!
      nmax = 2
      do n = 1, nmax
      if( n.eq.1 ) it = itsle
      if( n.eq.2 ) it = itpvs
!
      write(nft,'(/2x,"*** pllist ***  itim =",i7,"  time =",1pe12.4,
     >   "  ia =",i2,"  it =",i3)') itim,time,ia,it
!x    call plprnt(nft,"plc",ia,it,1.0d0, 1)
      call plprnt(nft,"vna",ia,it,1.0d19,1)
      call plprnt(nft,"vva",ia,it,0.0d0, 1)
      call plprnt(nft,"vti",ia,it,0.0d0, 1)
      call plprnt(nft,"vte",ia,it,0.0d0, 1)
!x    call plprnt(nft,"vcs",ia,it,0.0d0, 1)
!
      enddo
!
      return
      end
