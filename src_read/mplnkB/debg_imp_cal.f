!**********************************************************************
      subroutine debg_imp_cal
!**********************************************************************
      use cimden, only : nzmx, tdnz, twne, twrd
      use cntcom, only : iplx, iply
      use cntpls, only : dene, teme
      use cplcom, only : vne, vte
      use cunit,  only : cdgr, mygrp, n6
      implicit none

      integer :: maxz, ic, i, j, iz
      real(8) :: zne, zte

      if( mygrp /= cdgr(1) .and. mygrp /= cdgr(3)
     > .and. mygrp /= cdgr(4) .and. mygrp /= cdgr(5)) goto 920

!::No use
      if( mygrp.eq.cdgr(1) ) then
      write(n6,'(2x,"call mcplas in soldor  ",i3)') mygrp
      call mcplas
      endif
!
!::debug
!xx   cmsg = "sub. lnkimp_cal  Total"
!xx   call dbg_wrad(trim(cmsg),nv,twrd,twci,tdnz)
!      goto 920

 900  continue
      maxz = min0( nzmx, 3 )
      do ic = 100, 300, 50
      zne = 0.0d0
      zte = 0.0d0
      j = iplx(ic)
      i = iply(ic)
      if( j.gt.0 .and. i.gt.0 ) then
      zne = vne(j,i)
      zte = vte(j,i)
      endif
!
      write(n6,'(2x,i6,2i5,"  Ne =",1p2e12.3,"  Te =",1p2e12.3,
     >  " twne,twrd,tdnz =",1p2e12.3,1p5e12.3)')
     >   ic, j, i, zne, dene(ic), zte, teme(ic), twne(ic),
     >   twrd(ic), (tdnz(iz,ic),iz=0,maxz)
      enddo
!
 920  continue
!xx      call trmark("lnkimp_cal","return")
!
      return
      end
