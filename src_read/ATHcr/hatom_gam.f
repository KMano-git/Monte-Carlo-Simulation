!***********************************************************************
      subroutine hatom_gam(ate,ane,saha,rz0,rz1,rz2)
!***********************************************************************
!
!       ate   [i] : electron temperature   [eV]
!       ane   [i] : electron density       [1/m3]
!       saha(40)  [o] : Z(p) Saha-Bolzmann coefficient [m3]
!       rz0 (40)  [o] : r0(p)*Z(p)         [m3]
!       rz1 (40)  [o] : r1(p)*[Z(p)/Z(1)]  [m3]
!       rz2 (40)  [o] : r2(p)*??           [m3]
!                        temporary [Z(p)/Z(1)]
!
!      n(p) = r0(p)*Z(p)*nz*ne + r1(p)*[Z(p)/Z(1)]/ne*n(1)*ne
!            <--------->         <------------------>  n0
!               rz0                   rz1
!
!       n(p)_ion = ne*ni*rz0    [1/m3]
!       n(p)_atm = ne*n0*rz1    [1/m3]
!       n(p)_mol = ne*nD2*rz2   [1/m3]
!
!         See. Plasma Spectroscopy  by T. Fujimoto  p.91
!         Consult with   Nakano  tel. 7341
!
!-----------------------------------------------------------------------
      implicit none
!
!::argument
      real*8  ate, ane
      real*8  saha(40), rz0(40), rz1(40), rz2(40)
!
!::local common
      real*8  osc(40,40), a(40,40), c(40,40), f(40,40)
      real*8  s(40), alpha(40), beta(40), h2(40)
      integer lim, lup, k
      save osc, a
!
!::local variables
      real*8  zne, zte
      real*8  sgv0, sgv1, sgv2
      integer  lp, i
      data lp/0/; save lp
!
!::MKS ==> CGS
      zne = ane*1.0d-6
      zte = ate
!
      lp = lp + 1
      if( lp.eq.1 ) then    
      call einstn(osc,a,40)
      endif
!
      call clsaha(zte,saha)

      lim=40
      lup=35
!
      call ratcof(zte,osc,saha,c,f,s,alpha,beta,k)
      call molcof(zte,h2)
      call popcof(zne,saha,c,f,h2,s,a,alpha,beta,
     >           lup,lim,rz0,rz1,rz2)
!
!x      do i = 2, 30
!x      write(6,'(2x,i3,1p3e12.4)') i,r0(i),r1(i),r2(i)
!x      enddo
!
!::CGS ==> MKS
      do i = 1, 40
      saha(i) = saha(i)*1.0d-6
      rz0(i)  = rz0(i)*1.0d-6
      rz1(i)  = rz1(i)*1.0d-6
      rz2(i)  = rz2(i)*1.0d-6
      enddo
!
      return
      end
