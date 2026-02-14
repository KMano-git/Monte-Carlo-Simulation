!***********************************************************************
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   subroutine pyslos(nj)
      subroutine pyslos(j)
!***********************************************************************
!
!     simple loss model in scrape off layer
!
!        dq1a = - q1a/tau
!        dq2a = 0.0
!        dq3  = - gi * q3/tau
!        dq4 =  - ge * q4/tau
!
!       tau = Lpara / vpara
!       Lpara : parallel distance between stagnation point and Xp
!       vpara : 0.1 x Cs
!
!-----------------------------------------------------------------------
      use cplcom, only : ama, ff, nion, q1a, q3, q4, ss
      use cplmet, only : icspx, icwl1
      use cpmpls, only : alsfr, clsni, clste, clsti
      use cunit,  only : n6
      implicit none
!
!::argument
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  nj
      integer, intent(in) :: j
!
!::local vvariables
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer iwst, iwen, j, i, iw, ia, m1a, m2a, m3, m4
      integer iwst, iwen, i, ia, m1a, m3, m4
      real*8  zfr, tsi, twi, twe
!
      iwst = icwl1 + 1
      iwen = icspx
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   j = nj
!
      tsi = 0.0d0
      twi = 0.0d0
      twe = 0.0d0
!
! modified 3/2 lines organize local variables and include files by kamata 2021/05/31
!ik   do iw = iwst, iwen
!ik   i = iw
!ik   zfr = alsfr(iw)
      do i = iwst, iwen
      zfr = alsfr(i)
!
      do ia = 1, nion
      m1a = 2*ia - 1
! modified 3/2 lines organize local variables and include files by kamata 2021/05/31
!ik   m2a = 2*ia
!ik   ss(m1a,m1a,iw) = ss(m1a,m1a,iw)     + zfr*clsni
!ik   ff(m1a,iw) = ff(m1a,iw) - q1a(j,i,ia)*zfr*clsni
      ss(m1a,m1a,i) = ss(m1a,m1a,i)     + zfr*clsni
      ff(m1a,i) = ff(m1a,i) - q1a(j,i,ia)*zfr*clsni
      enddo
!
      m3 = 2*nion + 1
      m4 = 2*nion + 2
! modified 4/4 lines organize local variables and include files by kamata 2021/05/31
!ik   ss(m3,m3,iw) = ss(m3,m3,iw)   + zfr*clsti
!ik   ss(m4,m4,iw) = ss(m4,m4,iw)   + zfr*clste
!ik   ff(m3,iw) = ff(m3,iw) - q3(j,i)*zfr*clsti
!ik   ff(m4,iw) = ff(m4,iw) - q4(j,i)*zfr*clste
      ss(m3,m3,i) = ss(m3,m3,i)   + zfr*clsti
      ss(m4,m4,i) = ss(m4,m4,i)   + zfr*clste
      ff(m3,i) = ff(m3,i) - q3(j,i)*zfr*clsti
      ff(m4,i) = ff(m4,i) - q4(j,i)*zfr*clste
!
      tsi = tsi - q1a(j,i,1)/ama(1)*zfr*clsni
      twi = twi - q3(j,i)   *zfr*clsti
      twe = twe - q4(j,i)   *zfr*clste
!
      enddo
      write(n6,'(2x,"tsi =",1pe10.3,"  twi =",1pe10.3,"  twe =",
     >  1pe10.3)') tsi, twi, twe
!
      return
      end
