!   difference btween pldiskL.f  pldisk.f
!
!<       time = qtim       ! pldiskL.f
!> !xxx  time = qtim       ! time in topics code
!
c***********************************************************************
      subroutine pldisk(nft,kk,cft)
c***********************************************************************
      use mpi
      use csize
      use cplcom
      use cplmet
      use cgdcom
      use csonic
      use cunit
      implicit none
c
      integer nft, kk
      character cft*(*)
c
c::local variables
      character ctim*5, cmsg(2)*5, cnum*15, dsn*80, dsn2*80, comd*80
      logical lex, lopn
      integer i, nls, nle, nbf, ierr, lenx
      data cmsg(1)/"write"/, cmsg(2)/"read"/
      integer itmd
      data itmd/-1/
      integer mji, indxr
      real*8 tms
!
      if( cft(1:1).eq." " ) return
!
!x      write(n6,'(/2x,"*** pldisk ***  kk =",i2,2x,a)')
!x     >  kk, cft(1:lenx(cft))
c
c----------------------------------------------------------------------
c::in case of interrupt
c----------------------------------------------------------------------
      if( lstop.eq.1 ) then
      cft = cft(1:lenx(cft))//"E"
      write(n6,'(/2x,"*** pldisk ***  find error  lstop =",i3,
     >  "  file name =",a)') lstop,cft(1:lenx(cft))
      itmd = -1
      endif
c
c----------------------------------------------------------------------
c::write
c----------------------------------------------------------------------
      if( kk.eq.1 ) then
      if( nft.le.0 ) return
      if( itim.eq.itmd ) return
      if( lmype.ne.lmspe ) return
c
      mji = lenx(cft)
      dsn = cft(1:mji)
      qtim = time
c
c::rename for dtprf_xx
      if( dsn(1:8).eq."../dtprf" ) then
      inquire( file = dsn, exist = lex )
      if( lex .and. dsn(mji:mji).ne."E" ) then
      comd = "\mv  "//dsn(1:mji)//"  "//dsn(1:mji)//"B"
      call system( comd )
      write(n6,'(2x,"System command  ",a)') comd(1:lenx(comd))
      endif
      endif
c
      call nopen(nft,dsn,"binary write",dsn2,lex)
      write(n6,'(2x,"ndx =",i4,"  ndy =",i4,"  ndsp =",i4)')
     >   ndx, ndy, ndsp
c
      write(nft) qtim,q1a,q2a,q3,q4
      close (nft)
      itmd = itim
      endif
c
c----------------------------------------------------------------------
c::read
c----------------------------------------------------------------------
      if( kk.eq.2 ) then
      if( nft.eq.0 ) goto 900
      dsn = cft
c
c::master pe
      if( lmype.eq.lmspe ) then
c
      call nopen(nft,dsn,"binary",dsn2,lex)
      write(n6,'(2x,"ndx =",i4,"  ndy =",i4,"  ndsp =",i4)')
     >   ndx, ndy, ndsp
c
      read(nft,end=910,err=910) qtim,q1a,q2a,q3,q4
      close (nft)
c
      endif
c
!::[MPI_Bcast in pldisk]  cplcom/cmqnow/ (q1a,qtim)  10/04/21
      if( lnope.gt.1 ) then
      nls = loc( q1a(1,1,1) )
      nle = loc( qtim ) + 8
      nbf = nle - nls
      call MPI_Bcast( q1a, nbf, MPI_BYTE, lmspe, mywld, ierr ) 
      endif
      time = qtim
      endif     
c
c----------------------------------------------------------------------
c::debug write
c----------------------------------------------------------------------
      write(n6,'(2x,"*** pldisk ***   nft =",i3,"  itim =",i6,
     >  "  time =",1pe12.4,"  qtim =",1pe16.8,2x,a,2x,a)') 
     >   nft,itim,time,qtim,cmsg(kk),dsn2(1:lenx(dsn2))
      return
c
  900 continue
      write(n6,'(/2x,"*** pldisk *** initial run start.")')
      return
c
 910  continue
      call wexit( "pldisk","error occured during reading restart file")
      end
