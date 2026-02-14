!***********************************************************************
      subroutine out_mon
!***********************************************************************
      use cunit, only : lmspe, lmype
      implicit none
!
      integer  lp
      data     lp/0/
      save     lp
!
      if( lmype.ne.lmspe ) return
      lp = lp + 1
      if( lp.le.3 .or. mod(lp,10).ne.0 ) return
!
      call out_monv("deni")
      call out_monv("temi")
!
      return
      end
!
!***********************************************************************
      subroutine out_monv(cnam)
!***********************************************************************
      use cntcom, only : iplx, iply, mcel
      use cntpls, only : dene, deni, teme, temi
      use cplmet, only : icaxs, jcdp1, jcdp2, jcxp1, jcxp2
      use cpmpls, only : jmd1, rohmp, vlmp
      use csize,  only : ndmc
      use csonic, only : itim, time
      implicit none
!
!::argument
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   character  cnam*(*)
      character, intent(in) :: cnam*(*)
!
!::local variables
      integer    jxtb(9)
      real*8     tmp(9)
      integer    ift, k, jx, iy, ic, icmd, ia
      character  dsn*80
      real*8     var(0:ndmc)
      real*8     dvl
      logical    lex
!
      ift = 21
      dsn = "chkM_"//trim(cnam)
      inquire(file=trim(dsn),exist=lex)
      if( lex ) then
        open(unit=ift,file=dsn,position="append")
      else
        open(unit=ift,file=dsn)
      endif
!
      ia = 1
      if( trim(cnam).eq."deni" ) var(0:ndmc) = deni(0:ndmc,ia)
      if( trim(cnam).eq."dene" ) var(0:ndmc) = dene(0:ndmc)
      if( trim(cnam).eq."teme" ) var(0:ndmc) = teme(0:ndmc)
      if( trim(cnam).eq."temi" ) var(0:ndmc) = temi(0:ndmc)
!
      jxtb(1) = jcdp1
      jxtb(2) = jcdp1 + 1
      jxtb(3) = jcxp1
      jxtb(4) = jcxp1+1
      jxtb(5) = jmd1
      jxtb(6) = jcxp2-1
      jxtb(7) = jcxp2
      jxtb(8) = jcdp2 - 1
      jxtb(9) = jcdp2
!
      write(ift,'(2x,a,2x,"time =",1pe12.4,"  itim =",i6)')
     >     trim(cnam), time, itim
      write(ift,'(5x,"icmd",2x,"iplx",2x,"iply",
     >  4x,"iy",3x,"rovh",8x,"vol",9x,"dvl",4x,9(3x,i3,6x))')
     >   (jxtb(k),k=1,9)
!
      do iy = 1, icaxs
!
      do k = 1, 9
      jx = jxtb(k)
      ic = mcel(jx,iy)
      tmp(k) = var(ic)
      if( k.eq.5 ) icmd = ic
      enddo
!
      dvl = vlmp(iy) - vlmp(iy+1)
      write(ift,'(2x,i7,3i6,1p12e12.4)')
     >  icmd, iplx(icmd), iply(icmd),
     >  iy, rohmp(iy), vlmp(iy), dvl, (tmp(k),k=1,9)
      enddo
      close(ift)
!
      return
      end
!
!***********************************************************************
      subroutine dbg_temi(cmsg)
!***********************************************************************
      use cntcom, only : mcel
      use cntpls, only : temi
      use cplmet, only : icaxs
      use cpmpls, only : jmd1
      use csize,  only : ndy
      use csonic, only : itim, time
      use cunit,  only : mype
      implicit none
!
!::argument
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   character  cmsg*(*)
      character, intent(in) :: cmsg*(*)
!
!::local variables
      real*8     tmp(ndy)
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   integer    ift, k, jx, iy, ic, icmd, ia
      integer    ift, jx, iy, ic
      character  dsn*80
! deleted 2 lines organize local variables and include files by kamata 2021/06/16
!ik   real*8     var(0:ndmc)
!ik   real*8     dvl
!
! modified 3/1 lines organize local variables and include files by kamata 2021/06/16
!ik   integer  lp
!ik   data     lp/0/
!ik   save     lp
      integer :: lp = 0
!
      lp = lp + 1
!
      ift = 21
      write(dsn,'(a,i3.3)') "dbgM_", mype
      if( lp.eq.1 ) open(unit=ift,file=dsn)
      if( lp.ne.1 ) open(unit=ift,file=dsn,position="append")
!
! deleted 2 lines organize local variables and include files by kamata 2021/06/16
!ik   ia = 1
!ik   var(0:ndmc) = temi(0:ndmc)
!
      write(ift,'(2x,a,2x,"time =",1pe12.4,"  itim =",i6)')
     >  trim(cmsg), time, itim
      do iy = 1, icaxs
      jx = jmd1
      ic = mcel(jx,iy)
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   tmp(iy) = var(ic)
      tmp(iy) = temi(ic)
      enddo
!
      write(ift,'(1x,"temi",1p10e12.4)') (tmp(iy),iy=1,icaxs)
      close(ift)
!
      return
      end
