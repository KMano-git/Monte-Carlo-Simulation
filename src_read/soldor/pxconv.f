!**********************************************************************
      subroutine pxconv(it)
!**********************************************************************
!
!::NaN trap
!xx      if( it.eq.24 ) then
!xx      call check_NaN(m,n,jw,dd(m,n,jw),"pxconv-dd")
!xx      endif
!
!::viscocity  lvisc =1
!       ref.  Dr, K. Yamashita  Chiba univ  summer school
!
!----------------------------------------------------------------------
      use cphcns, only : cev
      use cplcom, only : ama, aza, beta, c00, c11, c12, c13, c14, c23
     >    , c32, c43, c53, cc, crod, dd, ee, faim, faip, ff, fl1a, fl2a
     >    , fl3, fl4, fl4m, fl4p, ibcyl, lmuscl, nion, nlp, q1a, q2a, q4
     >    , ss, ssnc, ssnv, vna, vea, vni, vte, vti, vva, vve
      use cplmet, only : hvol, hvsb, icel, jcel, jtmax
      use cplqcn, only : qfx_cv
      use csize,  only : ndeq, ndsp, ndx
      use csonic, only : itim, lstop
      use cunit,  only : n6
      implicit none
!
!::argument
      integer, intent(in) :: it
!
!::local variables
      integer k1, k2, nequ
      integer i, ipbc, jmax, jmax1, jw, j, jwm, jw2
      integer jp, jm, jc
      integer ia, ib, m1a, m2a, m3, m4, m1b, m2b
      real*8  zratd, zrate, zratc
      real*8  Ra1c(ndsp,ndsp), Ra2c(ndsp,ndsp), Ra3c(ndsp)
      real*8  Ra1d(ndsp,ndsp), Ra2d(ndsp,ndsp), Ra3d(ndsp)
      real*8  Ra1e(ndsp,ndsp), Ra2e(ndsp,ndsp), Ra3e(ndsp)
      real*8  Ra1s(ndsp,ndsp), Ra2s(ndsp,ndsp), Ra3s(ndsp)
      real*8  ha, ca2, ca, va, va2, zavv, zavh, zavc, zavc2
      real*8  ar11, ar12, ar13, ar21, ar22, ar23, ar31, ar32, ar33
      real*8  al11, al12, al13, al21, al22, al23, al31, al32, al33
      real*8  dsq1p, dsq2p, dsq3p, dsq4p, dsq1m, dsq2m, dsq3m, dsq4m
      real*8  h1q1p, h1q2p, h1q3p, h1q4p, h2q1m, h2q2m, h2q3m, h2q4m
      real*8  spq1m, spq2m, spq3m, smq1p, smq2p, smq3p
      real*8  bpq1m, bpq2m, bpq3m, bpq4m, bmq1p, bmq2p, bmq3p, bmq4p
      real*8  rmd1, rmd2, rmd3, rmp1, rmp2, rmp3, rmm1, rmm2, rmm3
      real*8  tap11,tap12,tap13,tap21,tap22,tap23,tap31,tap32,tap33
      real*8  tam11,tam12,tam13,tam21,tam22,tam23,tam31,tam32,tam33
      real*8  tab11,tab12,tab13,tab21,tab22,tab23,tab31,tab32,tab33
      real*8  tah44, taa44, tab44, tap44, tam44
      real*8  zf1p, zf2p, zf3p, zf4p, zf1m, zf2m, zf3m, zf4m
      real*8  zf1, zf2, zf3, zf4, zavq1, zavq2, zsqp, zsqm
      real*8  zwtp, zwtm, scm, scp, zne, zve
!
      real*8  zq1 (ndx), zq2 (ndx), zq3 (ndx), zq4 (ndx)
      real*8  zq1p(ndx), zq2p(ndx), zq3p(ndx), zq4p(ndx)
      real*8  zq1m(ndx), zq2m(ndx), zq3m(ndx), zq4m(ndx)
      real*8  zdq1(ndx), zdq2(ndx), zdq3(ndx), zdq4(ndx)
!
      real*8  aa11(ndx,ndsp),aa12(ndx,ndsp),aa13(ndx,ndsp)
      real*8  aa21(ndx,ndsp),aa22(ndx,ndsp),aa23(ndx,ndsp)
      real*8  aa31(ndx,ndsp),aa32(ndx,ndsp),aa33(ndx,ndsp)
!
      real*8  ap11(ndx,ndsp),ap12(ndx,ndsp),ap13(ndx,ndsp)
      real*8  ap21(ndx,ndsp),ap22(ndx,ndsp),ap23(ndx,ndsp)
      real*8  ap31(ndx,ndsp),ap32(ndx,ndsp),ap33(ndx,ndsp)
!
      real*8  am11(ndx,ndsp),am12(ndx,ndsp),am13(ndx,ndsp)
      real*8  am21(ndx,ndsp),am22(ndx,ndsp),am23(ndx,ndsp)
      real*8  am31(ndx,ndsp),am32(ndx,ndsp),am33(ndx,ndsp)
!
      real*8  ap44(ndx),am44(ndx)
      integer  n, m, jwb
!
!::electron velocity at j+1/2  !  %%%  2002/11/18
      real*8  xq1ap(ndx,ndsp), xq2ap(ndx,ndsp)
      real*8  xq1am(ndx,ndsp), xq2am(ndx,ndsp)
      real*8  xq1av(ndx,ndsp), xq2av(ndx,ndsp)
      real*8  wvep(ndx),wvem(ndx),wvev(ndx)
      real*8  zsmq1p,zsmq1m,zsmq1v,zsmq2p,zsmq2m,zsmq2v
!
      real*8 zsp(ndeq,ndeq),zsm(ndeq,ndeq),zdp(ndeq,ndeq),zdm(ndeq,ndeq)
      real*8 zcc(ndeq,ndeq),zee(ndeq,ndeq)
!
!::TVD scheme
      real*8  zfaip, zfaim
!
!::viscocity
      integer lvisc
      parameter (lvisc = 0)
      real(8) :: zvp, zhp, zcp2, zcp, zvm, zhm, zcm2, zcm
      real(8) :: rmd1R, rmd2R, rmd3R, rmd1L, rmd2L, rmd3L
      real(8) :: ab_rmd1, ab_rmd2, ab_rmd3, vs_rmd1, vs_rmd2, vs_rmd3
      real(8) :: rmde, rmdeR, rmdeL, ab_rmde, vs_rmde
!
      integer i6, lp
      character cdsn*80
      integer jmm
      real*8 zdvf1
      data lp/0/; save lp
!
!::minmod function
      real*8  fmmod, x, y
      fmmod(x,y)=dmax1( dmin1(x,0.0d0), dmin1(y,0.0d0), dmin1(x,y) )
!
      lp = lp + 1
!
!-----------------------------------------------------------------------
!::flux tube : it
!-----------------------------------------------------------------------
      i  = icel(1,it)
      ipbc  = ibcyl(it)
      jmax  = jtmax(it)
      jmax1 = jmax-1
!
      zfaip = faip
      zfaim = faim
      if( lmuscl.eq.0 ) then
      zfaip = 0.0d0
      zfaim = 0.0d0
      endif
!
!::debug write
      i6 = 0
!
!-----------------------------------------------------------------------
!::ion species : ia
!-----------------------------------------------------------------------
      do 100 ia = 1, nion
!
!-----------------------------------------------------------------------
!::q(jw) = Q(j) at j
!-----------------------------------------------------------------------
      do 110 jw = 1, jmax
      j = jcel(jw,it)
      zq1(jw) = q1a(j,i,ia)
      zq2(jw) = q2a(j,i,ia)
      zq3(jw) = vea(j,i,ia)
 110  continue
!
!-----------------------------------------------------------------------
!::dq(jw) = q(jw+1)-q(jw)   at j+1/2
!-----------------------------------------------------------------------
      do 120 jw = 1, jmax1
      zdq1(jw) = zq1(jw+1)-zq1(jw)
      zdq2(jw) = zq2(jw+1)-zq2(jw)
      zdq3(jw) = zq3(jw+1)-zq3(jw)
  120 continue
!
!-----------------------------------------------------------------------
!::qp(j-1/2) = Q+(j-1/2), qm(j+1/2) = Q-(j+1/2)  at cell boundary
!-----------------------------------------------------------------------
      do 130 jw = 2, jmax1
      j   = jcel(jw,it)
      jwm = jw-1
!
!::va,ca,Ha
      va  = zq2(jw)/zq1(jw)
      va2 = va*va
      ha  = c53*zq3(jw)/zq1(jw)-c13*va2
      ca2 = c23*(ha-c12*va2)
!
      if( ca2.le.0.0d0 ) then
      write(n6,'(2x,"warning pxconv  ca2 <= 0.0  ",1pe12.3,
     >  "  jw,it =",2i5,"  j,i =",2i5)') ca2, jw, it, j, i
      lstop = 1
      return
      endif
!
      ca  = sqrt(ca2)
!
!::matrix R(j)
      ar11 = c11/ca
      ar12 = c11/ca
      ar13 = c11/ca
      ar21 = va/ca
      ar22 = va/ca-c11
      ar23 = va/ca+c11
      ar31 = c12*va*va/ca
      ar32 = c12*va*va/ca+c32*ca-va
      ar33 = c12*va*va/ca+c32*ca+va
!
!::matrix L(j)
      al11 = ca-c13*va*va/ca
      al12 = c23*va/ca
      al13 = -c23/ca
      al21 = c12*va*(c11+c13*va/ca)
      al22 = -c12-c13*va/ca
      al23 = c13/ca
      al31 = c12*va*(-c11+c13*va/ca)
      al32 = c12-c13*va/ca
      al33 = c13/ca
!
!::matrix A(j)
      aa11(jw,ia) = c00
      aa12(jw,ia) = c11
      aa13(jw,ia) = c00
      aa21(jw,ia) = -c23*va*va
      aa22(jw,ia) = c43*va
      aa23(jw,ia) = c23
      aa31(jw,ia) = (-ha+c13*va*va)*va
      aa32(jw,ia) = ha-c23*va*va
      aa33(jw,ia) = c53*va
!
!::dq+(jw) = L(jw)*(Q(jw+1)-Q(jw))
      dsq1p = al11*zdq1(jw )+al12*zdq2(jw )+al13*zdq3(jw )
      dsq2p = al21*zdq1(jw )+al22*zdq2(jw )+al23*zdq3(jw )
      dsq3p = al31*zdq1(jw )+al32*zdq2(jw )+al33*zdq3(jw )
!
!::dq-(jw) = L(jw)*(Q(jw)-Q(jw-1))
      dsq1m = al11*zdq1(jwm)+al12*zdq2(jwm)+al13*zdq3(jwm)
      dsq2m = al21*zdq1(jwm)+al22*zdq2(jwm)+al23*zdq3(jwm)
      dsq3m = al31*zdq1(jwm)+al32*zdq2(jwm)+al33*zdq3(jwm)
!
!::h1_dq+(jw)
      h1q1p = fmmod( dsq1p, beta*dsq1m )
      h1q2p = fmmod( dsq2p, beta*dsq2m )
      h1q3p = fmmod( dsq3p, beta*dsq3m )
!
!::h2_dq-(jw)
      h2q1m = fmmod( dsq1m, beta*dsq1p )
      h2q2m = fmmod( dsq2m, beta*dsq2p )
      h2q3m = fmmod( dsq3m, beta*dsq3p )
!
!::sqp(j-1/2)
      spq1m = -zfaip*h2q1m-zfaim*h1q1p
      spq2m = -zfaip*h2q2m-zfaim*h1q2p
      spq3m = -zfaip*h2q3m-zfaim*h1q3p
!
!::sqm(j+1/2)
      smq1p =  zfaip*h1q1p+zfaim*h2q1m
      smq2p =  zfaip*h1q2p+zfaim*h2q2m
      smq3p =  zfaip*h1q3p+zfaim*h2q3m
!
!::bpq(j-1/2)
      bpq1m = zq1(jw) + ar11*spq1m+ar12*spq2m+ar13*spq3m
      bpq2m = zq2(jw) + ar21*spq1m+ar22*spq2m+ar23*spq3m
      bpq3m = zq3(jw) + ar31*spq1m+ar32*spq2m+ar33*spq3m
!
!::bmq(j+1/2)
      bmq1p = zq1(jw) + ar11*smq1p+ar12*smq2p+ar13*smq3p
      bmq2p = zq2(jw) + ar21*smq1p+ar22*smq2p+ar23*smq3p
      bmq3p = zq3(jw) + ar31*smq1p+ar32*smq2p+ar33*smq3p
!
!::store
      zq1p(jwm) = bpq1m
      zq2p(jwm) = bpq2m
      zq3p(jwm) = bpq3m
      zq1m(jw ) = bmq1p
      zq2m(jw ) = bmq2p
      zq3m(jw ) = bmq3p
!
  130 continue
!
!-----------------------------------------------------------------------
!::qp(jmax1+1/2),qm(1+1/2)  at cell boundary
!-----------------------------------------------------------------------
      if( ipbc.eq.1 ) then
        jw  = jmax1
        jw2 = 1
        zq1p(jw) = zq1p(jw2)
        zq2p(jw) = zq2p(jw2)
        zq3p(jw) = zq3p(jw2)
        jw  = 1
        jw2 = jmax1
        zq1m(jw) = zq1m(jw2)
        zq2m(jw) = zq2m(jw2)
        zq3m(jw) = zq3m(jw2)
      else
        jw  = jmax1
        jw2 = jmax
        zq1p(jw) = zq1(jw2)
        zq2p(jw) = zq2(jw2)
        zq3p(jw) = zq3(jw2)
        jw  = 1
        jw2 = 1
        zq1m(jw) = zq1(jw2)
        zq2m(jw) = zq2(jw2)
        zq3m(jw) = zq3(jw2)
      endif
!
!-----------------------------------------------------------------------
!::dQ=Q+(j+1/2)-Q-(j+1/2) at cell boundary
!-----------------------------------------------------------------------
      do 210 jw = 1, jmax1
      zdq1(jw) = zq1p(jw)-zq1m(jw)
      zdq2(jw) = zq2p(jw)-zq2m(jw)
      zdq3(jw) = zq3p(jw)-zq3m(jw)
  210 continue
!
!----------------------------------------------------------------------
!::jacobian matrix
!----------------------------------------------------------------------
      do 220 jw = 1, jmax1
      j = jcel(jw,it)
!
!::roe average
      if( zq1p(jw).le.0.0d0 .or. zq1m(jw).le.0.0d0 ) then
      write(n6,'(2x,"warning pxconv  zq1p,zq1m <= 0.0 ",1p2e12.3,
     >   "  jw,it =",2i5,"  j,i =",2i5)')
     >   zq1p(jw), zq1m(jw), jw, it, j, i
      write(n6,'(2x,"recalculation with lmuscl = 0")')
      lstop = 1
      return
      endif
!
      zsqp = sqrt(zq1p(jw))
      zsqm = sqrt(zq1m(jw))
      zwtp = zsqp/(zsqp+zsqm)
      zwtm = zsqm/(zsqp+zsqm)
      zavv = zwtp*(zq2p(jw)/zq1p(jw))
     >      +zwtm*(zq2m(jw)/zq1m(jw))
      zavh = zwtp*(c53*zq3p(jw)/zq1p(jw)-c13*zq2p(jw)**2/zq1p(jw)**2)
     >      +zwtm*(c53*zq3m(jw)/zq1m(jw)-c13*zq2m(jw)**2/zq1m(jw)**2)
      zavc2 = c23*(zavh-c12*zavv**2)
      if( zavc2.le.0.0d0 ) then
      write(n6,'(2x,"warning  pxconv zavc2 <= 0.0 ",1pe12.3,
     >  "  jw,it =",2i5,"  j,i =",2i5)') zavc2, jw, it, j, i
      lstop = 1
      return
      endif
      zavc=sqrt(zavc2)
!
!::ca,va
      va = zavv
      ca = zavc
      ha = zavh
!
!::matrix R(j+1/2)
      ar11 = c11/ca
      ar12 = c11/ca
      ar13 = c11/ca
      ar21 = va/ca
      ar22 = va/ca-c11
      ar23 = va/ca+c11
      ar31 = c12*va*va/ca
      ar32 = c12*va*va/ca+c32*ca-va
      ar33 = c12*va*va/ca+c32*ca+va
!
!::matrix L(j+1/2)
      al11 = ca-c13*va*va/ca
      al12 = c23*va/ca
      al13 = -c23/ca
      al21 = c12*va*(c11+c13*va/ca)
      al22 = -c12-c13*va/ca
      al23 = c13/ca
      al31 = c12*va*(-c11+c13*va/ca)
      al32 = c12-c13*va/ca
      al33 = c13/ca
!
!::eigen value
      rmd1 = va
      rmd2 = va-ca
      rmd3 = va+ca
      rmp1 = dmax1(rmd1,0.0d0)
      rmp2 = dmax1(rmd2,0.0d0)
      rmp3 = dmax1(rmd3,0.0d0)
      rmm1 = rmd1-rmp1
      rmm2 = rmd2-rmp2
      rmm3 = rmd3-rmp3
!
!::A+(j+1/2) = R*RMD+*L
      tap11 = ar11*rmp1*al11+ar12*rmp2*al21+ar13*rmp3*al31
      tap12 = ar11*rmp1*al12+ar12*rmp2*al22+ar13*rmp3*al32
      tap13 = ar11*rmp1*al13+ar12*rmp2*al23+ar13*rmp3*al33
      tap21 = ar21*rmp1*al11+ar22*rmp2*al21+ar23*rmp3*al31
      tap22 = ar21*rmp1*al12+ar22*rmp2*al22+ar23*rmp3*al32
      tap23 = ar21*rmp1*al13+ar22*rmp2*al23+ar23*rmp3*al33
      tap31 = ar31*rmp1*al11+ar32*rmp2*al21+ar33*rmp3*al31
      tap32 = ar31*rmp1*al12+ar32*rmp2*al22+ar33*rmp3*al32
      tap33 = ar31*rmp1*al13+ar32*rmp2*al23+ar33*rmp3*al33
!
!::A-(j+1/2) = R*RMD-*L
      tam11 = ar11*rmm1*al11+ar12*rmm2*al21+ar13*rmm3*al31
      tam12 = ar11*rmm1*al12+ar12*rmm2*al22+ar13*rmm3*al32
      tam13 = ar11*rmm1*al13+ar12*rmm2*al23+ar13*rmm3*al33
      tam21 = ar21*rmm1*al11+ar22*rmm2*al21+ar23*rmm3*al31
      tam22 = ar21*rmm1*al12+ar22*rmm2*al22+ar23*rmm3*al32
      tam23 = ar21*rmm1*al13+ar22*rmm2*al23+ar23*rmm3*al33
      tam31 = ar31*rmm1*al11+ar32*rmm2*al21+ar33*rmm3*al31
      tam32 = ar31*rmm1*al12+ar32*rmm2*al22+ar33*rmm3*al32
      tam33 = ar31*rmm1*al13+ar32*rmm2*al23+ar33*rmm3*al33
!
!::abs(A(j+1/2)) = A+ - A-
      tab11 = tap11 - tam11
      tab12 = tap12 - tam12
      tab13 = tap13 - tam13
      tab21 = tap21 - tam21
      tab22 = tap22 - tam22
      tab23 = tap23 - tam23
      tab31 = tap31 - tam31
      tab32 = tap32 - tam32
      tab33 = tap33 - tam33
!
!-----
!::artificial viscosity for expansion wave
      if( lvisc.eq.1 ) then
      zvp = zq2p(jw)/zq1p(jw)
      zvm = zq2m(jw)/zq1m(jw)
      zhp = c53*zq3p(jw)/zq1p(jw)-c13*zq2p(jw)**2/zq1p(jw)**2
      zhm = c53*zq3m(jw)/zq1m(jw)-c13*zq2m(jw)**2/zq1m(jw)**2
      zcp2= c23*(zhp-c12*zvp**2)
      zcm2= c23*(zhm-c12*zvm**2)
      if( zcp2.le.0.0d0 .or. zcm2.le.0.0d0 ) then
      write(n6,'(2x,"warning  pxconv zcp2,zcm2 <= 0.0 ",1p2e12.3,
     >  "  jw,it =",2i5,"  j,i =",2i5)') zcp2, zcm2, jw, it, j, i
      lstop = 1
      return
      endif
!
      zcp = sqrt(zcp2)
      zcm = sqrt(zcm2)
!
      rmd1R = zvp
      rmd2R = zvp-zcp
      rmd3R = zvp+zcp
      rmd1L = zvm
      rmd2L = zvm-zcm
      rmd3L = zvm+zcm
!
      ab_rmd1 = rmp1 - rmm1
      ab_rmd2 = rmp2 - rmm2
      ab_rmd3 = rmp3 - rmm3
      vs_rmd1 = 0.0d0
      vs_rmd2 = 0.0d0
      vs_rmd3 = 0.0d0
      if( rmd1L<0.0d0 .and. rmd1R>0.0d0 ) vs_rmd1=0.5d0*(rmd1R-rmd1L)
      if( rmd2L<0.0d0 .and. rmd2R>0.0d0 ) vs_rmd2=0.5d0*(rmd2R-rmd2L)
      if( rmd3L<0.0d0 .and. rmd3R>0.0d0 ) vs_rmd3=0.5d0*(rmd3R-rmd3L)
!
      ab_rmd1 = ab_rmd1 + vs_rmd1
      ab_rmd2 = ab_rmd2 + vs_rmd2
      ab_rmd3 = ab_rmd3 + vs_rmd3
!
      tab11 = ar11*ab_rmd1*al11+ar12*ab_rmd2*al21+ar13*ab_rmd3*al31
      tab12 = ar11*ab_rmd1*al12+ar12*ab_rmd2*al22+ar13*ab_rmd3*al32
      tab13 = ar11*ab_rmd1*al13+ar12*ab_rmd2*al23+ar13*ab_rmd3*al33
      tab21 = ar21*ab_rmd1*al11+ar22*ab_rmd2*al21+ar23*ab_rmd3*al31
      tab22 = ar21*ab_rmd1*al12+ar22*ab_rmd2*al22+ar23*ab_rmd3*al32
      tab23 = ar21*ab_rmd1*al13+ar22*ab_rmd2*al23+ar23*ab_rmd3*al33
      tab31 = ar31*ab_rmd1*al11+ar32*ab_rmd2*al21+ar33*ab_rmd3*al31
      tab32 = ar31*ab_rmd1*al12+ar32*ab_rmd2*al22+ar33*ab_rmd3*al32
      tab33 = ar31*ab_rmd1*al13+ar32*ab_rmd2*al23+ar33*ab_rmd3*al33
      endif
!-----
!
!::flux at cell boundary
!   F(j+1/2) = 1/2 * ( F(Q+)  + F(Q-) - |A|*( Q+ - Q- ) )
      zf1p = zq2p(jw)
      zf1m = zq2m(jw)
      zf1  = 0.5d0*( zf1p + zf1m
     >              -  (tab11*zdq1(jw)+tab12*zdq2(jw)+tab13*zdq3(jw)) )
      zf2p = c23*zq2p(jw)**2/zq1p(jw) + c23*zq3p(jw)
      zf2m = c23*zq2m(jw)**2/zq1m(jw) + c23*zq3m(jw)
      zf2  = 0.5d0*( zf2p + zf2m
     >              -  (tab21*zdq1(jw)+tab22*zdq2(jw)+tab23*zdq3(jw)) )
      zf3p = c53*zq2p(jw)/zq1p(jw)*zq3p(jw)-c13*zq2p(jw)**3/zq1p(jw)**2
      zf3m = c53*zq2m(jw)/zq1m(jw)*zq3m(jw)-c13*zq2m(jw)**3/zq1m(jw)**2
      zf3  = 0.5d0*( zf3p + zf3m
     >              -  (tab31*zdq1(jw)+tab32*zdq2(jw)+tab33*zdq3(jw)) )
!
!::save for electron velocity at cell boundary
      zavq1 = c14*(zsqp+zsqm)**2
      zavq2 = c14*(zsqp+zsqm)*
     >            (zsqp*zq2p(jw)/zq1p(jw)+zsqm*zq2m(jw)/zq1m(jw))
      xq1ap(jw,ia) = zq1p(jw)
      xq1am(jw,ia) = zq1m(jw)
      xq2ap(jw,ia) = zq2p(jw)
      xq2am(jw,ia) = zq2m(jw)
!-----                           ! %%% 2002/11/19
      xq1av(jw,ia) = zavq1
      xq2av(jw,ia) = zavq2
!-----
!
!::store
      fl1a(jw,ia) = zf1
      fl2a(jw,ia) = zf2
      fl3 (jw)    = fl3(jw) + zf3
!
      ap11(jw,ia) = tap11
      ap12(jw,ia) = tap12
      ap13(jw,ia) = tap13
      ap21(jw,ia) = tap21
      ap22(jw,ia) = tap22
      ap23(jw,ia) = tap23
      ap31(jw,ia) = tap31
      ap32(jw,ia) = tap32
      ap33(jw,ia) = tap33
!
      am11(jw,ia) = tam11
      am12(jw,ia) = tam12
      am13(jw,ia) = tam13
      am21(jw,ia) = tam21
      am22(jw,ia) = tam22
      am23(jw,ia) = tam23
      am31(jw,ia) = tam31
      am32(jw,ia) = tam32
      am33(jw,ia) = tam33
!
!::debug write
      if( i6.gt.0 ) then
      zdvf1 = 0.0d0
      if( jw.gt.1 ) then
      jmm = jcel(jw-1,it)
      zdvf1 = (hvsb(j,i)*fl1a(jw,ia)-hvsb(jmm,i)*fl1a(jw-1,ia))/ama(ia)
      endif
      write(i6,'(2x,3i5,f8.2,1p3e12.3,0pf8.2,1p5e12.3)') jw, j, i,
     >  dfloat(j), q2a(j,i,ia)/ama(ia), vna(j,i,ia), vva(j,i,ia),
     >  dfloat(j)+0.5d0, zf1p/ama(ia), zf1m/ama(ia), zf1/ama(ia),
     >  zdvf1, (ssnc(j,i,ia)+ssnv(j,i,ia)*vna(j,i,ia))*hvol(j,i)
      endif
!
  220 continue  ! loop(jw)
!
      if( i6.gt.0 ) then
      jw = jmax
      j  = jcel(jw,it)
      write(i6,'(2x,3i5,f8.2,1p3e12.3)') jw, j, i,
     >  dfloat(j), q2a(j,i,ia)/ama(ia), vna(j,i,ia), vva(j,i,ia)
      endif
!
  100 continue  ! loop(ia)
!
!----------------------------------------------------------------------
!::qcon
!----------------------------------------------------------------------
      m3 = 2*nion + 1
      do jw = 1, jmax1
      j = jcel(jw,it)
      do ia = 1, nion
      m1a = 2*ia - 1
      m2a = 2*ia
      qfx_cv(j,i,m1a) = hvsb(j,i)*fl1a(jw,ia)
      qfx_cv(j,i,m2a) = hvsb(j,i)*fl2a(jw,ia)
      enddo
      qfx_cv(j,i,m3)  = hvsb(j,i)*fl3(jw)
      enddo
!
!----------------------------------------------------------------------
!::ss,cc,dd,ee    ! change loop  %%% 2002/11/7
!----------------------------------------------------------------------
      m3 = 2*nion + 1
      m4 = 2*nion + 2
!
      do jw = 2, jmax1
      jwm = jw - 1
      j   = jcel(jw,it)
      jp  = jcel(jw+1,it)
      jm  = jcel(jw-1,it)
!
!::clear
      nequ = 2*nion + 2
      do k2 = 1, nequ
      do k1 = 1, nequ
      zsp(k1,k2) = 0.0d0
      zsm(k1,k2) = 0.0d0
      zdm(k1,k2) = 0.0d0
      zdp(k1,k2) = 0.0d0
      zcc(k1,k2) = 0.0d0
      zee(k1,k2) = 0.0d0
      enddo
      enddo

!
!::conversion from q3a to Q1,Q2,Q3
      jc  = j
      do ia = 1, nion
      zratd = q1a(jc,i,ia)/ama(ia)/vni(jc,i)
      zrate = q1a(jp,i,ia)/ama(ia)/vni(jp,i)
      zratc = q1a(jm,i,ia)/ama(ia)/vni(jm,i)
      if( nion.eq.1 ) then
      zratd = 1.0d0; zrate = 1.0d0; zratc = 1.0d0
      endif
      Ra3d(ia) = zratd
      Ra3e(ia) = zrate
      Ra3c(ia) = zratc
      Ra3s(ia) = Ra3d(ia)
      do ib = 1, nion
      Ra1d(ia,ib) = (crod(ia,ib)-zratd)*(c32*vti(jc,i)*cev/ama(ib)
     >               -c12*vva(jc,i,ib)**2)
      Ra1e(ia,ib) = (crod(ia,ib)-zrate)*(c32*vti(jp,i)*cev/ama(ib)
     >               -c12*vva(jp,i,ib)**2)
      Ra1c(ia,ib) = (crod(ia,ib)-zratc)*(c32*vti(jm,i)*cev/ama(ib)
     >               -c12*vva(jm,i,ib)**2)
      Ra2d(ia,ib) = (crod(ia,ib)-zratd)*vva(jc,i,ib)
      Ra2e(ia,ib) = (crod(ia,ib)-zrate)*vva(jp,i,ib)
      Ra2c(ia,ib) = (crod(ia,ib)-zratc)*vva(jm,i,ib)
      Ra1s(ia,ib) = Ra1d(ia,ib)
      Ra2s(ia,ib) = Ra2d(ia,ib)
      enddo
      enddo
!
      do ia = 1, nion
      m1a = 2*ia - 1
      m2a = 2*ia
      scm = hvsb(jm,i)
      scp = hvsb(j,i)
!-----
      zsp(m1a,m1a) = zsp(m1a,m1a) + scp*aa11(jw, ia)
      zsp(m1a,m2a) = zsp(m1a,m2a) + scp*aa12(jw, ia)
      zsp(m1a,m3 ) = zsp(m1a,m3 ) + scp*aa13(jw, ia)*Ra3s(ia)
      zsp(m2a,m1a) = zsp(m2a,m1a) + scp*aa21(jw, ia)
      zsp(m2a,m2a) = zsp(m2a,m2a) + scp*aa22(jw, ia)
      zsp(m2a,m3 ) = zsp(m2a,m3 ) + scp*aa23(jw, ia)*Ra3s(ia)
      zsp(m3, m1a) = zsp(m3, m1a) + scp*aa31(jw, ia)
      zsp(m3, m2a) = zsp(m3, m2a) + scp*aa32(jw, ia)
      zsp(m3, m3 ) = zsp(m3, m3 ) + scp*aa33(jw, ia)*Ra3s(ia)
!-----
      zsm(m1a,m1a) = zsm(m1a,m1a) + scm*aa11(jw, ia)
      zsm(m1a,m2a) = zsm(m1a,m2a) + scm*aa12(jw, ia)
      zsm(m1a,m3 ) = zsm(m1a,m3 ) + scm*aa13(jw, ia)*Ra3s(ia)
      zsm(m2a,m1a) = zsm(m2a,m1a) + scm*aa21(jw, ia)
      zsm(m2a,m2a) = zsm(m2a,m2a) + scm*aa22(jw, ia)
      zsm(m2a,m3 ) = zsm(m2a,m3 ) + scm*aa23(jw, ia)*Ra3s(ia)
      zsm(m3, m1a) = zsm(m3, m1a) + scm*aa31(jw, ia)
      zsm(m3, m2a) = zsm(m3, m2a) + scm*aa32(jw, ia)
      zsm(m3, m3 ) = zsm(m3, m3 ) + scm*aa33(jw, ia)*Ra3s(ia)
!-----
      zcc(m1a,m1a) = zcc(m1a,m1a) + scm*ap11(jwm,ia)
      zcc(m1a,m2a) = zcc(m1a,m2a) + scm*ap12(jwm,ia)
      zcc(m1a,m3 ) = zcc(m1a,m3 ) + scm*ap13(jwm,ia)*Ra3c(ia)
      zcc(m2a,m1a) = zcc(m2a,m1a) + scm*ap21(jwm,ia)
      zcc(m2a,m2a) = zcc(m2a,m2a) + scm*ap22(jwm,ia)
      zcc(m2a,m3 ) = zcc(m2a,m3 ) + scm*ap23(jwm,ia)*Ra3c(ia)
      zcc(m3, m1a) = zcc(m3, m1a) + scm*ap31(jwm,ia)
      zcc(m3, m2a) = zcc(m3, m2a) + scm*ap32(jwm,ia)
      zcc(m3, m3 ) = zcc(m3, m3 ) + scm*ap33(jwm,ia)*Ra3c(ia)
!-----
      zee(m1a,m1a) = zee(m1a,m1a) + scp*am11(jw, ia)
      zee(m1a,m2a) = zee(m1a,m2a) + scp*am12(jw, ia)
      zee(m1a,m3 ) = zee(m1a,m3 ) + scp*am13(jw, ia)*Ra3e(ia)
      zee(m2a,m1a) = zee(m2a,m1a) + scp*am21(jw, ia)
      zee(m2a,m2a) = zee(m2a,m2a) + scp*am22(jw, ia)
      zee(m2a,m3 ) = zee(m2a,m3 ) + scp*am23(jw, ia)*Ra3e(ia)
      zee(m3, m1a) = zee(m3, m1a) + scp*am31(jw, ia)
      zee(m3, m2a) = zee(m3, m2a) + scp*am32(jw, ia)
      zee(m3, m3 ) = zee(m3, m3 ) + scp*am33(jw, ia)*Ra3e(ia)
!-----
      zdp(m1a,m1a) = zdp(m1a,m1a) - scp*am11(jw,ia)
      zdp(m1a,m2a) = zdp(m1a,m2a) - scp*am12(jw,ia)
      zdp(m1a,m3 ) = zdp(m1a,m3 ) - scp*am13(jw, ia)*Ra3d(ia)
      zdp(m2a,m1a) = zdp(m2a,m1a) - scp*am21(jw,ia)
      zdp(m2a,m2a) = zdp(m2a,m2a) - scp*am22(jw,ia)
      zdp(m2a,m3 ) = zdp(m2a,m3 ) - scp*am23(jw, ia)*Ra3d(ia)
      zdp(m3, m1a) = zdp(m3, m1a) - scp*am31(jw,ia)
      zdp(m3, m2a) = zdp(m3, m2a) - scp*am32(jw,ia)
      zdp(m3, m3 ) = zdp(m3, m3 ) - scp*am33(jw, ia)*Ra3d(ia)
!-----
      zdm(m1a,m1a) = zdm(m1a,m1a) - scm*ap11(jwm,ia)
      zdm(m1a,m2a) = zdm(m1a,m2a) - scm*ap12(jwm,ia)
      zdm(m1a,m3 ) = zdm(m1a,m3 ) - scm*ap13(jwm,ia)*Ra3d(ia)
      zdm(m2a,m1a) = zdm(m2a,m1a) - scm*ap21(jwm,ia)
      zdm(m2a,m2a) = zdm(m2a,m2a) - scm*ap22(jwm,ia)
      zdm(m2a,m3 ) = zdm(m2a,m3 ) - scm*ap23(jwm,ia)*Ra3d(ia)
      zdm(m3, m1a) = zdm(m3, m1a) - scm*ap31(jwm,ia)
      zdm(m3, m2a) = zdm(m3, m2a) - scm*ap32(jwm,ia)
      zdm(m3, m3 ) = zdm(m3, m3 ) - scm*ap33(jwm,ia)*Ra3d(ia)
!-----
!
      do ib = 1, nion
      m1b = 2*ib - 1
      m2b = 2*ib
!-----
      zsp(m1a,m1b) = zsp(m1a,m1b) + scp*aa13(jw,ia)*Ra1s(ia,ib)
      zsp(m1a,m2b) = zsp(m1a,m2b) + scp*aa13(jw,ia)*Ra2s(ia,ib)
      zsp(m2a,m1b) = zsp(m2a,m1b) + scp*aa23(jw,ia)*Ra1s(ia,ib)
      zsp(m2a,m2b) = zsp(m2a,m2b) + scp*aa23(jw,ia)*Ra2s(ia,ib)
      zsp(m3, m1b) = zsp(m3, m1b) + scp*aa33(jw,ia)*Ra1s(ia,ib)
      zsp(m3, m2b) = zsp(m3, m2b) + scp*aa33(jw,ia)*Ra2s(ia,ib)
!-----
      zsm(m1a,m1b) = zsm(m1a,m1b) + scm*aa13(jw,ia)*Ra1s(ia,ib)
      zsm(m1a,m2b) = zsm(m1a,m2b) + scm*aa13(jw,ia)*Ra2s(ia,ib)
      zsm(m2a,m1b) = zsm(m2a,m1b) + scm*aa23(jw,ia)*Ra1s(ia,ib)
      zsm(m2a,m2b) = zsm(m2a,m2b) + scm*aa23(jw,ia)*Ra2s(ia,ib)
      zsm(m3, m1b) = zsm(m3, m1b) + scm*aa33(jw,ia)*Ra1s(ia,ib)
      zsm(m3, m2b) = zsm(m3, m2b) + scm*aa33(jw,ia)*Ra2s(ia,ib)
!-----
      zcc(m1a,m1b) = zcc(m1a,m1b) + scm*ap13(jwm,ia)*Ra1c(ia,ib)
      zcc(m1a,m2b) = zcc(m1a,m2b) + scm*ap13(jwm,ia)*Ra2c(ia,ib)
      zcc(m2a,m1b) = zcc(m2a,m1b) + scm*ap23(jwm,ia)*Ra1c(ia,ib)
      zcc(m2a,m2b) = zcc(m2a,m2b) + scm*ap23(jwm,ia)*Ra2c(ia,ib)
      zcc(m3, m1b) = zcc(m3, m1b) + scm*ap33(jwm,ia)*Ra1c(ia,ib)
      zcc(m3, m2b) = zcc(m3, m2b) + scm*ap33(jwm,ia)*Ra2c(ia,ib)
!-----
      zee(m1a,m1b) = zee(m1a,m1b) + scp*am13(jw,ia)*Ra1e(ia,ib)
      zee(m1a,m2b) = zee(m1a,m2b) + scp*am13(jw,ia)*Ra2e(ia,ib)
      zee(m2a,m1b) = zee(m2a,m1b) + scp*am23(jw,ia)*Ra1e(ia,ib)
      zee(m2a,m2b) = zee(m2a,m2b) + scp*am23(jw,ia)*Ra2e(ia,ib)
      zee(m3, m1b) = zee(m3, m1b) + scp*am33(jw,ia)*Ra1e(ia,ib)
      zee(m3, m2b) = zee(m3, m2b) + scp*am33(jw,ia)*Ra2e(ia,ib)
!-----
      zdp(m1a,m1b) = zdp(m1a,m1b) - scp*am13(jw, ia)*Ra1d(ia,ib)
      zdp(m1a,m2b) = zdp(m1a,m2b) - scp*am13(jw, ia)*Ra2d(ia,ib)
      zdp(m2a,m1b) = zdp(m2a,m1b) - scp*am23(jw, ia)*Ra1d(ia,ib)
      zdp(m2a,m2b) = zdp(m2a,m2b) - scp*am23(jw, ia)*Ra2d(ia,ib)
      zdp(m3, m1b) = zdp(m3, m1b) - scp*am33(jw, ia)*Ra1d(ia,ib)
      zdp(m3, m2b) = zdp(m3, m2b) - scp*am33(jw, ia)*Ra2d(ia,ib)
!-----
      zdm(m1a,m1b) = zdm(m1a,m1b) - scm*ap13(jwm,ia)*Ra1d(ia,ib)
      zdm(m1a,m2b) = zdm(m1a,m2b) - scm*ap13(jwm,ia)*Ra2d(ia,ib)
      zdm(m2a,m1b) = zdm(m2a,m1b) - scm*ap23(jwm,ia)*Ra1d(ia,ib)
      zdm(m2a,m2b) = zdm(m2a,m2b) - scm*ap23(jwm,ia)*Ra2d(ia,ib)
      zdm(m3, m1b) = zdm(m3, m1b) - scm*ap33(jwm,ia)*Ra1d(ia,ib)
      zdm(m3, m2b) = zdm(m3, m2b) - scm*ap33(jwm,ia)*Ra2d(ia,ib)
!-----
      enddo  !  loop(ib)
      enddo  !  loop(ia)
!
!::save
      do m = 1, m3
      do n = 1, m3
      ss(m,n,jw) = ss(m,n,jw) + zsp(m,n) - zsm(m,n)
      cc(m,n,jw) = cc(m,n,jw)            - zcc(m,n)
      dd(m,n,jw) = dd(m,n,jw) + zdp(m,n) - zdm(m,n)
      ee(m,n,jw) = ee(m,n,jw) + zee(m,n)
      enddo  !  loop(n)
      enddo  !  loop(m)
!
!----------------------------------------------------------------------
!::bc-condition   ! %%% 2002/11/11
!----------------------------------------------------------------------
      if( jw.eq.2 ) then
      jwb = jw-1
      do m = 1, m3
      do n = 1, m3
      dd(m,n,jwb) = dd(m,n,jwb) + zcc(m,n)
      ee(m,n,jwb) = ee(m,n,jwb) + zsm(m,n) + zdm(m,n)
      enddo
      enddo
      endif
      if( jw.eq.jmax1 ) then
      jwb = jw + 1
      do m = 1, m3
      do n = 1, m3
      cc(m,n,jwb) = cc(m,n,jwb) - zsp(m,n) - zdp(m,n)
      dd(m,n,jwb) = dd(m,n,jwb) - zee(m,n)
!-----
      enddo
      enddo
      endif
!
      enddo  !  loop(jw)
!
!----------------------------------------------------------------------
!::ff
!----------------------------------------------------------------------
      do 360 jw = 2, jmax1
      j  = jcel(jw,it)
      jm = jcel(jw-1,it)
      ff(m3,jw) = - hvsb(j,i)*fl3(jw) + hvsb(jm,i)*fl3(jw-1)
      do 370 ia = 1, nion
      m1a = 2*ia - 1
      m2a = 2*ia
      ff(m1a,jw) = - hvsb(j,i)*fl1a(jw,ia) + hvsb(jm,i)*fl1a(jw-1,ia)
      ff(m2a,jw) = - hvsb(j,i)*fl2a(jw,ia) + hvsb(jm,i)*fl2a(jw-1,ia)
!
 370  continue
 360  continue
!
!----------------------------------------------------------------------
!::electron
!----------------------------------------------------------------------
!
!::ve+, ve-, <ve> at cell boundary
      do 510 jw = 1, jmax1
      zsmq1p = 0.0d0
      zsmq1m = 0.0d0
      zsmq1v = 0.0d0
      zsmq2p = 0.0d0
      zsmq2m = 0.0d0
      zsmq2v = 0.0d0
      do ia = 1, nion
      zsmq1p = zsmq1p + aza(ia)/ama(ia)*xq1ap(jw,ia)
      zsmq1m = zsmq1m + aza(ia)/ama(ia)*xq1am(jw,ia)
      zsmq1v = zsmq1v + aza(ia)/ama(ia)*xq1av(jw,ia)
      zsmq2p = zsmq2p + aza(ia)/ama(ia)*xq2ap(jw,ia)
      zsmq2m = zsmq2m + aza(ia)/ama(ia)*xq2am(jw,ia)
      zsmq2v = zsmq2v + aza(ia)/ama(ia)*xq2av(jw,ia)
      enddo
      wvep(jw) = zsmq2p/zsmq1p
      wvem(jw) = zsmq2m/zsmq1m
      wvev(jw) = zsmq2v/zsmq1v
 510  continue
!
!::zq4(j) = qe(j) at cell center
      do 520 jw = 1, jmax
      j = jcel(jw,it)
      zq4(jw) = q4(j,i)
 520  continue
!
!::zdq4(j+1/2) = q4(j+1)-q4(j)
      do 530 jw = 1, jmax1
      zdq4(jw) = zq4(jw+1)-zq4(jw)
 530  continue
!
!::q+(j-1/2),q-(j+1/2) at cell boundary
      do 540 jw = 2, jmax1
      j   = jcel(jw,it)
      jwm = jw-1
!
!::dq+,dq-
      dsq4p = zdq4(jw)
      dsq4m = zdq4(jwm)
!
!::h1_dq+,h2_dq-
      h1q4p = fmmod( dsq4p, beta*dsq4m )
      h2q4m = fmmod( dsq4m, beta*dsq4p )
!
!::q+(j-1/2),q-(j+1/2)
      bpq4m = zq4(jw) - zfaip*h2q4m - zfaim*h1q4p
      bmq4p = zq4(jw) + zfaip*h1q4p + zfaim*h2q4m
!
!::store
      zq4p(jwm) = bpq4m
      zq4m(jw)  = bmq4p
!
 540  continue
!
!----------------------------------------------------------------------
!::q+(jmax-1+1/2),q-(1+1/2) at cell boundary
!----------------------------------------------------------------------
      if( ipbc.eq.1 ) then
        jw  = jmax1
        jw2 = 1
        zq4p(jw) = zq4p(jw2)
        jw  = 1
        jw2 = jmax1
        zq4m(jw) = zq4m(jw2)
      else
        jw  = jmax1
        jw2 = jmax
        zq4p(jw) = zq4(jw2)
        jw  = 1
        jw2 = 1
        zq4m(jw) = zq4(jw2)
      endif
!
!----------------------------------------------------------------------
!::jacobian matrix
!----------------------------------------------------------------------
      do 550 jw = 1, jmax1
      j = jcel(jw,it)
!
!::a+(j+1/2),a-(j+1/2),|a(j+1/2)|
      tah44 = c53*wvev(jw)
      taa44 = tah44
      tap44 = dmax1(taa44,0.0d0)
      tam44 = taa44 - tap44
      tab44 = tap44 - tam44
!
!::artificial viscosity
      if( lvisc.eq.1 ) then
      rmde  = c53*wvev(jw)
      rmdeR = c53*wvep(jw)
      rmdeL = c53*wvem(jw)
      ab_rmde = dabs(rmde)
      vs_rmde = 0.0d0
      if( rmdeL<0.0d0 .and. rmdeR>0.0d0 ) vs_rmde=0.5d0*(rmdeR-rmdeL)
!
      ab_rmde = ab_rmde + vs_rmde
      tab44 = ab_rmde
      endif
!
!::a(j) for ss
      zne = 0.0d0
      zve = 0.0d0
      do 555 ia = 1, nion
      zne = zne + aza(ia)/ama(ia)*q1a(j,i,ia)
      zve = zve + aza(ia)/ama(ia)*q2a(j,i,ia)
 555  continue
      zve = zve/zne
!
!::flux at cell boundary
      zf4p = c53*wvep(jw)*zq4p(jw)
      zf4m = c53*wvem(jw)*zq4m(jw)
      zf4  = 0.5d0*( zf4p + zf4m - tab44*(zq4p(jw)-zq4m(jw)) )
!
!::store
      fl4(jw)  = zf4
      ap44(jw) = tap44
      am44(jw) = tam44
!
      fl4p(jw) = zf4p
      fl4m(jw) = zf4m
!
 550  continue
!
!----------------------------------------------------------------------
!::qcon
!----------------------------------------------------------------------
      m4 = 2*nion + 2
      do jw = 1, jmax1
      j = jcel(jw,it)
      qfx_cv(j,i,m4) = hvsb(j,i)*fl4(jw)
      enddo
!
!----------------------------------------------------------------------
!::ss,cc,dd,ee
!----------------------------------------------------------------------
      m4 = 2*nion + 2
      do jw = 2, jmax1
      jwm = jw - 1
      j   = jcel(jw,it)
      jm  = jcel(jw-1,it)
      scm = hvsb(jm,i)
      scp = hvsb(j,i)
!-----
      zsp(m4,m4) = scp*c53*wvem(jw)
      zsm(m4,m4) = scm*c53*wvep(jwm)
!-----
      zcc(m4,m4) = scm*ap44(jwm)
      zee(m4,m4) = scp*am44(jw)
      zdp(m4,m4) =-scp*am44(jw)
      zdm(m4,m4) =-scm*ap44(jwm)
!-----
      ss(m4,m4,jw) = ss(m4,m4,jw) + zsp(m4,m4) - zsm(m4,m4)
      cc(m4,m4,jw) = cc(m4,m4,jw)              - zcc(m4,m4)
      dd(m4,m4,jw) = dd(m4,m4,jw) + zdp(m4,m4) - zdm(m4,m4)
      ee(m4,m4,jw) = ee(m4,m4,jw) + zee(m4,m4)
!-----
      ff(m4,jw)    = -scp*fl4(jw) + scm*fl4(jw-1)
!
!----------------------------------------------------------------------
!::bc-condition  ! %%% 2002/11/7
!----------------------------------------------------------------------
      if( jw.eq.2 ) then
      jwb = jw - 1
      dd(m4,m4,jwb) = dd(m4,m4,jwb) + zcc(m4,m4)
      ee(m4,m4,jwb) = ee(m4,m4,jwb) + zsm(m4,m4) + zdm(m4,m4)
      endif
      if( jw.eq.jmax1 ) then
      jwb = jw + 1
      cc(m4,m4,jwb) = cc(m4,m4,jwb) - zsp(m4,m4) - zdp(m4,m4)
      dd(m4,m4,jwb) = dd(m4,m4,jwb) - zee(m4,m4)
      endif
!
      enddo  ! loop(jw)
!
!----------------------------------------------------------------------
!::ff at boundary ! %%% 2002/11/11
!----------------------------------------------------------------------
      jwb = 1
      jw  = 1
      j   = jcel(jw,it)
      do ia = 1, nion
      m1a = 2*ia - 1
      m2a = 2*ia
      ff(m1a,jwb) = - hvsb(j,i)*fl1a(jw,ia)
      ff(m2a,jwb) = - hvsb(j,i)*fl2a(jw,ia)
      enddo
      ff(m3, jwb) = - hvsb(j,i)*fl3 (jw)
      ff(m4, jwb) = - hvsb(j,i)*fl4 (jw)
      jwb = jmax
      jw  = jmax1
      j   = jcel(jw,it)
      do ia = 1, nion
      m1a = 2*ia - 1
      m2a = 2*ia
      ff(m1a,jwb) = + hvsb(j,i)*fl1a(jw,ia)
      ff(m2a,jwb) = + hvsb(j,i)*fl2a(jw,ia)
      enddo
      ff(m3,jwb) =  + hvsb(j,i)*fl3 (jw)
      ff(m4,jwb) =  + hvsb(j,i)*fl4 (jw)
!
!----------------------------------------------------------------------
!::debug
!----------------------------------------------------------------------
      if( i6.gt.0 ) close(i6)
!
      return
      end
