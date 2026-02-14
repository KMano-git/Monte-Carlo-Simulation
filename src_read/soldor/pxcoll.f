!***********************************************************************
      subroutine pxcoll(nt)
!***********************************************************************
!
!     dq/dt + div(F) + Hc = S
!
!     dlt( dq/dt + divF + Hc - S ) = -( dq/dt + div(F) + Hc - S )^l
!
!      H1a = 0
!      H2a = +sum(b)[ffrab*(va-vb)]   [zf2]   qvl_cl(m2a)
!      H3  = +ne*feqp*[eV]*(Ti-Te)    [zf3]   qvl_cl(m3)
!      H4  = - H3                    [zf4]   qvl_cl(m4)
!
!          note.  V*dq/dt + H*V + ... = 0
!                 ffrab = cfrab*clg*na*nb/Ti**1.5
!                 feqp  = sum(a)[3me/ma*1/tauea]
!
!-----------------------------------------------------------------------
      use cphcns, only : cev
      use cplcom, only : alnr, ama, aza, c13, c23, cfrab, feqp, ff
     >    , mdl_eqp, nion, q1a, q2a, ss, vna, vnezef, vni, vte, vti, vva
      use cplmet, only : hvol, icel, jcel, jtmax, kreg
      use cplqcn, only : qvl_cl
      implicit none
!
!::argument
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  nt
      integer, intent(in) :: nt
!
!::local varaibles
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  it, i, jmax, jmax1, jw, j, irg
      integer  i, jmax, jmax1, jw, j, irg
      integer  ia, m1a, m2a, m3, m4, ib, m1b, m2b
      real*8   zfrab, zfeq, znei, zs31, zs32, zs33, zs34, zf2, zf3
      real*8   znezef
      real*8   fceq(10)
!
!::input data
      real(8) :: hnlmt = 5.0d19
      real(8) :: hfeqp = 30.0d0

!
! modified 3/2 lines organize local variables and include files by kamata 2021/05/31
!ik   it = nt
!ik   i  = icel(1,it)
!ik   jmax  = jtmax(it)
      i  = icel(1,nt)
      jmax  = jtmax(nt)
      jmax1 = jmax - 1
!
      m3 = 2*nion + 1
      m4 = 2*nion + 2
!
!x      write(n6,'(/2x,"***  pxcoll  ***   it =",i3)') it
!
!-----------------------------------------------------------------------
!::H2a = sum(b)[ffrab*(va-vb)]  friction force
!-----------------------------------------------------------------------
      do jw = 2, jmax1
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   j = jcel(jw,it)
      j = jcel(jw,nt)
!
      do ia = 1, nion
      m1a = 2*ia - 1
      m2a = 2*ia
!
      zf2 = 0.0d0
      do ib = 1, nion
      m1b = 2*ib - 1
      m2b = 2*ib
      zfrab = cfrab(ia,ib)*alnr(j)*vna(j,i,ia)*vna(j,i,ib)
     >         /(vti(j,i)*sqrt(vti(j,i)))*hvol(j,i)
      ss(m2a,m1a,jw) = ss(m2a,m1a,jw)- zfrab*q2a(j,i,ia)/q1a(j,i,ia)**2
      ss(m2a,m2a,jw) = ss(m2a,m2a,jw)+ zfrab/q1a(j,i,ia)
      ss(m2a,m1b,jw) = ss(m2a,m1b,jw)+ zfrab*q2a(j,i,ib)/q1a(j,i,ib)**2
      ss(m2a,m2b,jw) = ss(m2a,m2b,jw)- zfrab/q1a(j,i,ib)
      zf2  = zf2 + zfrab*(vva(j,i,ia)-vva(j,i,ib))
      enddo  !  loop(ib)
      ff(m2a,jw)     = ff(m2a,jw) - zf2
!
!::qcon  ! %%% 2002/10/25
      qvl_cl(j,i,m2a) = +zf2
!
      enddo  !  loop(ia)
      enddo  !  loop(jw)
!
!-----------------------------------------------------------------------
!:: H3 = -H4 = +ne*feqp*[ceV]*(Ti-Te)
!-----------------------------------------------------------------------
!
!::MODEL enhance equi-partition in O/I divertor region
      fceq(1:10) = 1.0d0
      if( mdl_eqp.eq.1 ) then
        irg = 1
        jw  = 2
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik     j  = jcel(jw,it)
        j  = jcel(jw,nt)
        if( vni(j,i).lt.hnlmt ) fceq(irg) = hfeqp
!
!xx     irg = 3
!xx     jw  = jmax1
!xx     j   = jcel(jw,it)
!xx     if( vni(j,i).lt.hnlmt ) fceq(irg) = hfeqp
      endif
!
      do jw = 2, jmax1
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   j  = jcel(jw,it)
      j  = jcel(jw,nt)
      znezef = vnezef(j,i)
      zfeq = feqp(j)*hvol(j,i)
!xx   znei = vne(j,i)/vni(j,i)
      znei = znezef/vni(j,i)
!
      irg = kreg(j,i)
      zfeq  = zfeq*fceq(irg)
!
      zf3 = 0.0d0
      do ia = 1, nion
      m1a = 2*ia - 1
      m2a = 2*ia
!
      zs31 =  zfeq*(  znei*(c13*vva(j,i,ia)**2-vti(j,i)*cev/ama(ia))
     >              + aza(ia)*vte(j,i)*cev/ama(ia) )
      zs32 = -zfeq*c23*znei*vva(j,i,ia)
!xx   zf3  = zf3 + vne(j,i)*zfeq*(vti(j,i)-vte(j,i))*cev
      zf3  = zf3 + znezef*zfeq*(vti(j,i)-vte(j,i))*cev
!
      ss(m3,m1a,jw) = ss(m3,m1a,jw) + zs31
      ss(m3,m2a,jw) = ss(m3,m2a,jw) + zs32
      ss(m4,m1a,jw) = ss(m4,m1a,jw) - zs31
      ss(m4,m2a,jw) = ss(m4,m2a,jw) - zs32
      enddo  ! loop(ia)
!
      zs33 =  zfeq*c23*znei
      zs34 = -zfeq*c23
      ss(m3,m3,jw) = ss(m3,m3,jw) + zs33
      ss(m3,m4,jw) = ss(m3,m4,jw) + zs34
      ss(m4,m3,jw) = ss(m4,m3,jw) - zs33
      ss(m4,m4,jw) = ss(m4,m4,jw) - zs34

      ff(m3,jw)    = ff(m3,jw)    - zf3
      ff(m4,jw)    = ff(m4,jw)    + zf3
!
!::qcon  !  %%% 2002/10/25
      qvl_cl(j,i,m3) = +zf3
      qvl_cl(j,i,m4) = -zf3
!
      enddo  ! loop(jw)
!
      return
      end
