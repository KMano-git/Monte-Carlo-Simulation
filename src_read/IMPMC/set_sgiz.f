!***********************************************************************
      subroutine set_sgiz(kk)
!***********************************************************************
!
!       ionization of neutral impurity : cdi0 (/sec)
!
!-----------------------------------------------------------------------
      use catcom, only : ip_ion
      use cimcom, only : cdi, cdi0, irmax, ndis, sgmtn
      use cntpls, only : dene, teme
      use csize,  only : ndmc
      use cunit,  only : n6
      implicit none
!
!::argument
! modified 1/1 lines organize local variables and include files by kamata 2021/06/28
!ik   integer  kk
      integer, intent(in) :: kk
!
!::local variables
      integer  ic, isno
      real*8   zne, zte, zero, acof(ndis)
!
!::local common
      real*8 fcmtn
!
      fcmtn = 1.0d0
!
!::clear
      zero = 1.0d-30
      do ic = 0, ndmc
      cdi0(ic) = zero
      enddo
!
!xxx if( emttyp(1:3).eq."Cch" ) goto 100
!
!::C ==> C+
      if( kk.eq.0 ) then
      write(n6,'(2x,"*** set_sgiz ***  C0 ==> C+")')
      do ic = 1, irmax
      cdi0(ic) = cdi(0,ic)
      enddo
      return
      endif
!
!::CD4 ==> C+
 100  continue
      write(n6,'(2x,"*** set_sgiz ***  CD4 ==> C+")')
      write(n6,'(2x,"sigv(Te) ==> sigv(Te-sgmtn)  sgmtn =",f8.3,
     >   "  fcmtn =",f8.3)') sgmtn, fcmtn
      do ic = 1, irmax
      zne = dene(ic)
      zte = teme(ic)
      zte = zte - sgmtn
      if( zte.lt.0.0d0 ) then
      cdi0(ic) = zero
      else
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   call atm_eval2(ip_ion,zne,zte,acof,isno,ndis)
      call atm_eval2( ip_ion(1), zne, zte, acof, isno, ndis, 1 )
      acof(1) = acof(1)*fcmtn
      cdi0(ic) = dmax1(zne*acof(1),zero)
      endif
      if( mod(ic,1000).eq.0 ) then
      write(n6,'(2x,"new sgv(C0=>C+) ",i7,1p3e12.3," =>",1pe12.3)')
     >  ic, dene(ic), teme(ic), cdi(0,ic), cdi0(ic)
      endif
      enddo
!
      return
      end
