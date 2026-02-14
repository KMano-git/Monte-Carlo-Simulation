!**********************************************************************
      subroutine pledprf
!**********************************************************************
      use cphcns, only : cev
      use cplcom, only : ama, aza, mdl_edt, nion, q1a, q2a, q3, q4, qw1a
     >    , qw2a, qw3, qw4, vna, vne, vni, vte, vti, vva
      use cplmet, only : icel, itmax, itpve, itsle, itsls, jcel, jtmax
     >    , jtmin
      use cunit,  only : n6
      implicit none
!
!::local variables
! modified 2/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  it, jt, jc, ic, ia, kbc, n6sv
!ik   integer  itst, iten
      integer  it, jt, jc, ic, ia, kbc
      real*8   zne, zni, zq3
!ik   character  dsn*20
!
      write(n6,'(2x,"*** pledprf ***  Passed mdl_edt =",i2)') mdl_edt
!
!-----------------------------------------------------------------------
!::edit profile
!-----------------------------------------------------------------------
!xx   if( mdl_edt.eq.1 ) call pledprv(itst,iten)  ! (2010/07/06)
!xx   if( mdl_edt.eq.2 ) call plednmx             ! (2010/07/11)
!xx   call pledit                !  temin, timin   (2010/01/15)
!
      if( mdl_edt.eq.0 ) then
        write(n6,'(2x,"No change of plasma parameter")')
        goto 100
      endif
!
      if( mdl_edt.eq.1 ) call pledtem  ! 2012/08/20
      if( mdl_edt.eq.2 ) call pledsol  ! 2011/05/19
      if( mdl_edt.eq.3 ) call pledmpl  ! 2011/05/20
      if( mdl_edt.eq.4 ) call pledden  ! 2011/10/17
 100  continue
!
!-----------------------------------------------------------------------
!::conservative variables & vne
!-----------------------------------------------------------------------
      do 510 it = 1, itmax
      do 520 jt = jtmin(it), jtmax(it)
      jc  = jcel(jt,it)
      ic  = icel(jt,it)
      zne = 0.0d0
      zni = 0.0d0
      zq3 = 0.0
      do 530 ia = 1,nion
      q1a(jc,ic,ia) = ama(ia)*vna(jc,ic,ia)
      q2a(jc,ic,ia) = ama(ia)*vna(jc,ic,ia)*vva(jc,ic,ia)
      zne = zne + aza(ia)*vna(jc,ic,ia)
      zni = zni + vna(jc,ic,ia)
      zq3 = zq3 + 1.5d0*vna(jc,ic,ia)*vti(jc,ic)*cev
     >          + 0.5d0*ama(ia)*vna(jc,ic,ia)*vva(jc,ic,ia)**2
 530  continue
      vne(jc,ic) = zne
      vni(jc,ic) = zni
      q3(jc,ic)  = zq3
      q4(jc,ic)  = 1.5d0*vne(jc,ic)*vte(jc,ic)*cev
 520  continue
 510  continue
!
!-----------------------------------------------------------------------
!::boundary condition at sol-wall and prv-wall
!-----------------------------------------------------------------------
      kbc = 1
      it  = itsls
      jt  = (jtmin(it)+jtmax(it))/2
      jc = jcel(jt,it)
      ic = icel(jt,it)
      do ia = 1, nion
      qw1a(ia,kbc) = q1a(jc,ic,ia)
      qw2a(ia,kbc) = q2a(jc,ic,ia)
      enddo
      qw3(kbc) = q3(jc,ic)
      qw4(kbc) = q4(jc,ic)
!
      kbc = 2
      it = itpve
      jt = (jtmin(it)+jtmax(it))/2
      jc = jcel(jt,it)
      ic = icel(jt,it)
      do ia = 1, nion
      qw1a(ia,kbc) = q1a(jc,ic,ia)
      qw2a(ia,kbc) = q2a(jc,ic,ia)
      enddo
      qw3(kbc) = q3(jc,ic)
      qw4(kbc) = q4(jc,ic)
!
!-----------------------------------------------------------------------
!::non conservative variables
!-----------------------------------------------------------------------
      call plauxv
!
      if( mdl_edt.eq.1 ) then
      write(n6,'(2x,"After edit")')
      write(n6,'(2x,"itsle =",i3)') itsle
      write(n6,'(2x,"Odiv Te =",0p15f8.3)')
     > (vte(jcel(jtmin(it),it),icel(jtmin(it),it)),it=itsls,itpve)
      write(n6,'(2x,"Odiv Ti =",0p15f8.3)')
     > (vti(jcel(jtmin(it),it),icel(jtmin(it),it)),it=itsls,itpve)
      write(n6,'(2x,"Idiv Te =",0p15f8.3)')
     > (vte(jcel(jtmax(it),it),icel(jtmax(it),it)),it=itsls,itpve)
      write(n6,'(2x,"Idiv Ti =",0p15f8.3)')
     > (vti(jcel(jtmax(it),it),icel(jtmax(it),it)),it=itsls,itpve)
      endif
!
      return
      end
