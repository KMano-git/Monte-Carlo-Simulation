!***********************************************************************
!
!       print variables in the main plasma
!          parameter & diffusion coefficients
!
!-----------------------------------------------------------------------
!***********************************************************************
      subroutine out_sol2
!***********************************************************************
      use cplcom, only : vdda, vdet, vdxe, vdxi, vni, vte, vti
      use cplmet, only : icmpe, icmps, jcmax
      use csonic, only : itim, time
      use cunit,  only : lmspe, lmype
      implicit none
!
!::local variables
      integer  ia
      integer  jc, ic, i6, icst, icen
!
      if( lmype.ne.lmspe ) return
!
      if( itim.gt.3 .and. mod(itim,1000).ne.0 ) return
      i6 = 83000
!
      write(i6,'(/2x,"out_sol2   time =",1pe14.6,"  itim =",i6,
     >  "  icmps,icmpe =",2i5)') time, itim, icmps, icmpe

 611  format(1x,a,2x,i3,1p15e11.3)
!
      ia = 1
      icst = icmps - 2
      icen = icmps + 2
!
      do ic = icst, icen
      write(i6,611)  "Ni",ic,(vni(jc,ic),jc=1,jcmax,10)
      enddo
!
      do ic = icst, icen
      write(i6,611)  "Ti",ic,(vti(jc,ic),jc=1,jcmax,10)
      enddo
!
      do ic = icst, icen
      write(i6,611)  "Te",ic,(vte(jc,ic),jc=1,jcmax,10)
      enddo
!
      do ic = icst, icen
      write(i6,611)  "Di",ic,(vdda(jc,ic,ia),jc=1,jcmax,10)
      enddo
!
      do ic = icst, icen
      write(i6,611)  "Et",ic,(vdet(jc,ic,ia),jc=1,jcmax,10)
      enddo
!
      do ic = icst, icen
      write(i6,611)  "Xi",ic,(vdxi(jc,ic),jc=1,jcmax,10)
      enddo
!
      do ic = icst, icen
      write(i6,611)  "Xe",ic,(vdxe(jc,ic),jc=1,jcmax,10)
      enddo
!
      return

!x      ia = 1
!x      call out_solv("vni",vni)
!x      call out_solv("vti",vti)
!x      call out_solv("vte",vte)
!x      call out_solv("Dna",vdda(1,1,ia))
!x      call out_solv("Eta",vdet(1,1,ia))
!x      call out_solv("Xi", vdxi)
!x      call out_solv("Xe", vdxe)
!x
!x      return
      end
!***********************************************************************
      subroutine out_sol
!***********************************************************************
      use cplcom, only : vdda, vdet, vdxe, vdxi, vni, vte, vti
      use cplmet, only : icmpe, icmps, jcmax
      use csonic, only : itim, time
      use cunit,  only : cdgr, lmspe, lmype, mygrp
      implicit none
!
!::local variables
      integer  ia
      integer  jc, ic, i6, icst, icen
!
      if( mygrp.ne.cdgr(1) ) return
      if( lmype.ne.lmspe ) return
!
      if( itim.gt.3 .and. mod(itim,1000).ne.0 ) return
      i6 = 82000
!
      write(i6,'(/2x,"out_sol   time =",1pe14.6,"  itim =",i6,
     >  "  icmps,icmpe =",2i5)') time, itim, icmps, icmpe

 611  format(1x,a,2x,i3,1p15e11.3)
!
      ia = 1
      icst = icmps - 2
      icen = icmps + 2
!
      do ic = icst, icen
      write(i6,611)  "Ni",ic,(vni(jc,ic),jc=1,jcmax,10)
      enddo
!
      do ic = icst, icen
      write(i6,611)  "Ti",ic,(vti(jc,ic),jc=1,jcmax,10)
      enddo
!
      do ic = icst, icen
      write(i6,611)  "Te",ic,(vte(jc,ic),jc=1,jcmax,10)
      enddo
!
      do ic = icst, icen
      write(i6,611)  "Di",ic,(vdda(jc,ic,ia),jc=1,jcmax,10)
      enddo
!
      do ic = icst, icen
      write(i6,611)  "Et",ic,(vdet(jc,ic,ia),jc=1,jcmax,10)
      enddo
!
      do ic = icst, icen
      write(i6,611)  "Xi",ic,(vdxi(jc,ic),jc=1,jcmax,10)
      enddo
!
      do ic = icst, icen
      write(i6,611)  "Xe",ic,(vdxe(jc,ic),jc=1,jcmax,10)
      enddo
!
      return

!x      ia = 1
!x      call out_solv("vni",vni)
!x      call out_solv("vti",vti)
!x      call out_solv("vte",vte)
!x      call out_solv("Dna",vdda(1,1,ia))
!x      call out_solv("Eta",vdet(1,1,ia))
!x      call out_solv("Xi", vdxi)
!x      call out_solv("Xe", vdxe)
!x
!x      return
      end
!
!***********************************************************************
      subroutine out_solv(cnam,var)
!***********************************************************************
      use cplmet, only : icaxs, jcdp1, jcdp2, jcxp1, jcxp2
      use cpmpls, only : jmd1, rohmp, vlmp
      use csize,  only : ndx, ndy
      use csonic, only : itim, time
      implicit none
!
!::argument
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   character  cnam*(*)
!ik   real*8  var(ndx,ndy)
      character, intent(in) :: cnam*(*)
      real(8),   intent(in) :: var(ndx,ndy)
!
!::local variables
      integer  jxtb(9)
      real*8   tmp(9)
      integer  ift, k, jx, iy
      character  dsn*80
      logical  lex
!
      ift = 21
      dsn = "chkS_"//trim(cnam)
      inquire(file=trim(dsn),exist=lex)
      if( lex ) then
        open(unit=ift,file=dsn,position="append")
      else
        open(unit=ift,file=dsn)
      endif
!
      jxtb(1) = jcdp1
      jxtb(2) = jcdp1+1
      jxtb(3) = jcxp1
      jxtb(4) = jcxp1+1
      jxtb(5) = jmd1
      jxtb(6) = jcxp2-1
      jxtb(7) = jcxp2
      jxtb(8) = jcdp2-1
      jxtb(9) = jcdp2
!
      write(ift,'(2x,a,2x,"time =",1pe12.4,"  itim =",i6)')
     >  trim(cnam), time, itim
      write(ift,'(4x,"iy",3x,"rovh",8x,"vol",9x,"dvl",4x,9(3x,i3,
     >  6x))') (jxtb(k),k=1,9)
!
      do iy = 1, icaxs
!
      do k = 1, 9
      jx = jxtb(k)
      tmp(k) = var(jx,iy)
      enddo
!
      write(ift,'(2x,i4,1p12e12.4)')
     >  iy, rohmp(iy), vlmp(iy),vlmp(iy)-vlmp(iy+1), (tmp(k),k=1,9)
      enddo
      close(ift)
!
      return
      end
