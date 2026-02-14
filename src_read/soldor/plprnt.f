!***********************************************************************
      subroutine plprnt(nft,vnam,na,nt,scl,kal)
!***********************************************************************
      use cntmnt, only : ssn, ssp, swe, swi
      use cplcom, only : vcs, vna, vte, vti, vva
      use cplmet, only : icel, itsle, jcel, jtmax, jtmin
      use csize,  only : ndx
      implicit none
!
!::argument
! modified 3/3 lines organize local variables and include files by kamata 2021/05/31
!ik   integer   nft, na, nt, kal
!ik   character vnam*(*)
!ik   real*8    var(ndx), scl
      integer,   intent(in) :: nft, na, nt, kal
      character, intent(in) :: vnam*(*)
      real(8),   intent(in) :: scl
!
!::local variables
      integer  ia, it, jtst, jten, jt, j, i, ivscl, nl
      integer  jt1, jt2, jt3, jt4
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   real*8   vmax, zscl, fmax, vscl
      real(8)  fmax, var(ndx), vmax, vscl, zscl
      character cnam*6
!
      cnam = vnam
      ia = na
      it = nt
      if( it.eq.0 ) it = itsle
      if( ia.eq.0 ) ia = 1
!
      jtst = jtmin(it)
      jten = jtmax(it)
!
!::set var
      if( vnam(1:3).eq."plc" ) then
        call plenc(it,var)
      else
!
      do 110 jt = jtst, jten
      j = jcel(jt,it)
      i = icel(jt,it)
      if( vnam(1:3).eq."vna" ) var(jt) = vna(j,i,ia)
      if( vnam(1:3).eq."vva" ) var(jt) = vva(j,i,ia)
      if( vnam(1:3).eq."vti" ) var(jt) = vti(j,i)
      if( vnam(1:3).eq."vte" ) var(jt) = vte(j,i)
      if( vnam(1:3).eq."vcs" ) var(jt) = vcs(j,i)
!
      if( vnam(1:3).eq."ssn" ) var(jt) = ssn(j,i,ia)
      if( vnam(1:3).eq."ssp" ) var(jt) = ssp(j,i,ia)
      if( vnam(1:3).eq."swi" ) var(jt) = swi(j,i)
      if( vnam(1:3).eq."swe" ) var(jt) = swe(j,i)
 110  continue
      endif
!
!::scale
      vmax = 0.0
      do 210 jt = jtst, jten
      if( dabs(var(jt)).gt.vmax ) vmax = dabs(var(jt))
 210  continue
!
      if( vmax.eq.0.0 ) then
      write(nft,602) vnam,vmax
 602  format(2x,a6,1x,1pe11.4,11(0pf10.4))
      return
      endif
!
      zscl = dabs(scl)
      if( dabs(zscl).le.1.0d-20 ) zscl = 1.0d0
      fmax = vmax/zscl
      if( fmax.lt.0.1d0 .or. fmax.gt.10.0d0 ) then
      vscl = dlog10(vmax)
      ivscl = int(vscl)
      if( vscl.lt.0.0d0 ) ivscl = ivscl - 1
      zscl = 10.0d0**ivscl
      endif
!
!::thin out
      nl  = (jten-jtst+1)/10
      jt1 = jtst
      jt2 = min0(jt1 + 19,jten)
      jt3 = max0((nl-1)*10 + 1,jt2+1)
      jt4 = jten
!
!::list
      if( kal.eq.1 ) then
      write(nft,604) cnam,zscl,(var(jt)/zscl,jt=jtst,jten)
 604  format(2x,a6,1x,1pe10.3,10(0pf10.4)
     >  :/(2x,6x,1x,10x,10(0pf10.5)))
      else
!
      write(nft,604) cnam,zscl,(var(jt)/zscl,jt=jt1,jt2)
      if( jt2.lt.jten ) then
      write(nft,606) (var(jt)/zscl,jt=jt3,jt4)
 606  format(2x,6x,1x,10x,10(0pf10.5):/(2x,6x,1x,10x,10(0pf10.5)))
      endif
      endif
!
      return
      end
