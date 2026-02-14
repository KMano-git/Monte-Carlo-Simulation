!***********************************************************************
      subroutine xdt_cpfl(cenv,cvar,dnam)
!***********************************************************************
      use cdfcom, only : nxvr, xvar
      use cxdcom, only : ndxdt, ndxkn, nxs, nyp, xdnam
      implicit none
!
!::argument
! modified 1/2 lines organize local variables and include files by kamata 2021/06/16
!ik   character cenv*(*),cvar*(*),dnam*(*)
      character, intent(in)  :: cenv*(*), cvar*(*)
      character, intent(out) :: dnam*(*)

!::local variables
      character cxs*40, dsn*80
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   integer   nft, i, lenx
      integer   nft, i
      integer   mjs(110),mje(110),nk
      integer   nfl; data nfl/0/
      logical   lex
! added 2 lines organize local variables and include files by kamata 2021/06/16
! function
      integer    lenx
!-----------------------------------------------------------------------
!
      nfl = nfl + 1
!
!::clear counter variables
      if( nfl.eq.1 ) then
        nxs = 0;  nyp = 0
      endif
!
!::file open
      nft = 21
      call nopen( nft, cenv, "text", dsn, lex )
      if( lenx(dsn).le.0 ) then
        call wexit( "xdt_cpfl", "no found dsn" )
      endif
      call cdf_init( nft, dsn )
!
!::get all variables
      if( cvar(1:3).eq."all" ) then
      do i = 1, nxvr
      cxs = xvar(i)
      call xdt_copy(cxs)
      enddo
!
!::select variables
      else
      call linsep(cvar," ,",nk,mjs,mje,110)
      do i = 1, nk
      cxs = cvar(mjs(i):mje(i))
      call xdt_copy(cxs)
      enddo
      endif
!
!::dnam ( [dd_00_elastic] )
      dnam = xdnam(nxs)
!
!::file close
      close(nft)
!
!::debug write
      write(6,'(/2x,"*** xdt_cpfl ***  nfl =",i2,"  dsn = [",a,"]",
     >  "  dnam =",a)') nfl,dsn(1:lenx(dsn)),dnam(1:lenx(dnam))
      write(6,'(4x,"number of variables ",2i8)') nxs,ndxkn
      write(6,'(4x,"number of data      ",2i8)') nyp,ndxdt
!
      return
      end
