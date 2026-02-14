!   difference btween pldiskL.f  pldisk.f
!
!<       time = qtim       ! pldiskL.f
!> !xxx  time = qtim       ! time in topics code
!
!***********************************************************************
      subroutine pldisk(nft,kk,cft)
!***********************************************************************
      use cplcom,      only : q1a, q2a, q3, q4, qtim, qtim_ini
      use csize,       only : ndsp, ndx, ndy
      use csonic,      only : itim, lstop, time
      use cunit,       only : lmspe, lmype, lnope, mywld, n6
      use mpi!,         only : mpi_bcast
      use topics_mod,  only : lgtpc
      implicit none
!
      integer,   intent(in)   :: nft, kk
      character, intent(inout) :: cft*(*)
!
!::local variables
      character cmsg(2)*5, dsn*80, dsn2*80, comd*80
      logical lex
      integer ierr, itag, ityp
      data cmsg(1)/"write"/, cmsg(2)/"read"/
      integer itmd
      data itmd/-1/
      integer    mji

! function
      integer    lenx
!
      if( cft(1:1).eq." " ) return
!
!----------------------------------------------------------------------
!::in case of interrupt
!----------------------------------------------------------------------
      if( lstop.eq.1 ) then
        cft = cft(1:lenx(cft))//"E"
        write(n6,'(/2x,"*** pldisk ***  find error  lstop =",i3,
     >    "  file name =",a)') lstop,cft(1:lenx(cft))
        itmd = -1
      endif
!
!----------------------------------------------------------------------
!::write
!----------------------------------------------------------------------
      if( kk.eq.1 ) then
        if( nft.le.0 ) return
        if( itim.eq.itmd ) return
        if( lmype.ne.lmspe ) return
!
        mji = lenx(cft)
        dsn = cft(1:mji)
        qtim = time
!
!::rename for dtprf_xx
        if( dsn(1:8).eq."../dtprf" ) then
          inquire( file = dsn, exist = lex )
          if( lex .and. dsn(mji:mji).ne."E" ) then
            comd = "\mv  "//dsn(1:mji)//"  "//dsn(1:mji)//"B"
            call system( comd )
            write(n6,'(2x,"System command  ",a)') comd(1:lenx(comd))
          endif
        endif
!
        call nopen(nft,dsn,"binary write",dsn2,lex)
        write(n6,'(2x,"ndx =",i4,"  ndy =",i4,"  ndsp =",i4)')
     >   ndx, ndy, ndsp
!
        write(nft) qtim,q1a,q2a,q3,q4
        close (nft)
        itmd = itim
      endif
!
!----------------------------------------------------------------------
!::read
!----------------------------------------------------------------------
      if( kk.eq.2 ) then
        if( nft.eq.0 ) goto 900
        dsn = cft
!
!::master pe
        if( lmype.eq.lmspe ) then
          call nopen(nft,dsn,"binary",dsn2,lex)
          write(n6,'(2x,"ndx =",i4,"  ndy =",i4,"  ndsp =",i4)')
     >     ndx, ndy, ndsp
          read(nft,end=910,err=910) qtim,q1a,q2a,q3,q4
          close (nft)
          qtim_ini = qtim ! save first value, qtim change in calculation 
        endif
!
!::[MPI_Bcast in pldisk]  cplcom/cmqnow/ (q1a,qtim)  10/04/21
        if( lnope.gt.1 ) then
          call cplcom_lsr
        endif
        if( lgtpc == 0 ) time = qtim
      endif
!
!----------------------------------------------------------------------
!::debug write
!----------------------------------------------------------------------
      write(n6,'(2x,"*** pldisk ***   nft =",i3,"  itim =",i6,
     >  "  time =",1pe12.4,"  qtim =",1pe16.8,2x,a,2x,a)')
     >   nft,itim,time,qtim,cmsg(kk),dsn2(1:lenx(dsn2))
      return
!
  900 continue
      write(n6,'(/2x,"*** pldisk *** initial run start.")')
      return
!
 910  continue
      call wexit( "pldisk","error occured during reading restart file")
      end
