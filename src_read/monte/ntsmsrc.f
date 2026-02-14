!**********************************************************************
      subroutine ntsmsrc(kk)
!**********************************************************************
      use cntctl, only : nclnt
      use cntsms, only : i6_ntsms
      use cntwcn, only : wcsrc, wden, wisrc, wssn, wssp, wswe, wswi
      use cplcom, only : nion
      use csize,  only : ndgs, ndmc
      use csonic, only : itim
      implicit none
!
!::argument
      integer, intent(in) :: kk
!
!::local variables
      integer  i6, kout, ia
      integer  k0, k1, k2, k3, k4
!
!::local common
      real(8), dimension(0:ndmc,ndgs) :: hssn, hssp, hden
      real(8), dimension(0:ndmc)      :: hswe, hswi
!
      i6 = 0
      kout = 0
      if( nclnt.eq.5 .and. wisrc.eq.4 ) kout = 1
      if( kout.eq.1 ) i6 = i6_ntsms
!
      if( i6.gt.0 ) then
      write(i6,'(/2x,"*** ntsmsrc ***  repl =",i2,"  itim =",i8,
     >  "  nclnt =",i5,"  isrc =",i2,2x,a)')
     >  kk, itim, nclnt, wisrc, wcsrc
      endif
!
      k0 = 0*kout
      k1 = 1*kout
      k2 = 2*kout
      k3 = 3*kout
      k4 = 4*kout
!
      do ia = 1, nion
      call mksmdt(wden(0,ia),hden(0,ia),k0)
      call mksmdt(wssn(0,ia),hssn(0,ia),k1)
      call mksmdt(wssp(0,ia),hssp(0,ia),k2)
      enddo
      call mksmdt(wswi(0),hswi(0),k3)
      call mksmdt(wswe(0),hswe(0),k4)
!
      if( kk.eq.0 ) return
!
!::replace
      wden(0:ndmc,1:ndgs) = hden(0:ndmc,1:ndgs)
      wssn(0:ndmc,1:ndgs) = hssn(0:ndmc,1:ndgs)
      wssp(0:ndmc,1:ndgs) = hssp(0:ndmc,1:ndgs)
      wswe(0:ndmc) = hswe(0:ndmc)
      wswi(0:ndmc) = hswi(0:ndmc)
!
      return
      end
!
!**********************************************************************
      subroutine ntsmset
!**********************************************************************
      use cntcom, only : iplx, iply, mcel, next
      use cntctl, only : mdl_ntsm
      use cntsms, only : i6_ntsms, icsm, itsm, ncsm, ndsm, nsmax, wtsm
      use cplmet, only : icel, itpve, itpvs, itsle, itsls, jcel, jtmax
     >    , jtmin
      use cunit,  only : lmspe, lmype, n6
      implicit none
!
!::local variables
      integer ::  ii, i6
      integer ::  it, jt, j, i, ic, ic1, ic2, k
      integer ::  n, knx
!
      i6_ntsms = 80000 + 301
      if( lmype.ne.lmspe ) i6_ntsms = 0
      i6 = i6_ntsms
      i6_ntsms = 0 !KH0160802 no debug write
!
      write(n6,'(/2x,"*** ntsmset ***  mdl_ntsm =",i3,"  i6_ntsms =",
     >  i7)') mdl_ntsm, i6_ntsms
!
      if( mdl_ntsm.eq.0 ) then
        write(n6,'(5x,"Source is NOT smoothed")')
        return
      endif
!
      if( i6.gt.0 ) write(i6,'(/2x,"*** ntsmset ***")')
!
      ii = 0
      do it = itsls, itpve
        if( it.eq.itsls ) cycle
        if( it.eq.itpve ) cycle
!
        do n = 1, 2
          if( n.eq.1 ) then
            jt = jtmin(it) + 1
            j  = jcel(jt,it)
            i  = icel(jt,it)
            ic = mcel(j,i)
            knx = 2
          else
            jt = jtmax(it) - 1
            j  = jcel(jt,it)
            i  = icel(jt,it)
            ic = mcel(j,i)
            knx = 4
          endif
!
          do k = 1, 4
            ic1 = ic
            ic2 = next(ic1,knx)
            if( ic1.le.0 .or. ic2.le.0 ) goto 920
!
            ii = ii + 1
            if( ii.gt.ndsm ) goto 910
            itsm(ii) = it
            icsm(ii,1) = ic1
            icsm(ii,2) = ic2
            wtsm(ii,1) = 1.0d0
            wtsm(ii,2) = 1.0d0
            ncsm(ii) = 2
!
            ii = ii + 1
            if( ii.gt.ndsm ) goto 910
            itsm(ii) = it
            icsm(ii,1) = ic2
            icsm(ii,2) = ic1
            wtsm(ii,1) = 1.0d0
            wtsm(ii,2) = 1.0d0
            ncsm(ii) = 2
!
            ic = next(ic2,knx)
          enddo   ! loop(k)   4 cells
        enddo   ! loop(n)   outer/inner
      enddo   ! loop(it)  tube
!
      nsmax = ii
!
!::debug write
      if( i6.gt.0 ) then
      write(i6,'(/2x,"nsmax =",i4,"  ndsm =",i4 )') nsmax, ndsm
      do i = 1, nsmax
      if( itsm(i).ne.itsle .and. itsm(i).ne.itpvs ) cycle
      write(i6,'(2x,i3,2x,5(i7,2i5,0pf7.2))')
     >  i, (icsm(i,k),iplx(icsm(i,k)),iply(icsm(i,k)),
     >   wtsm(i,k),k=1,ncsm(i))
      enddo
      endif
!
      return
!
!::error
 910  continue
      call wexit("ntsmset","dimension icsm(ndsm,9) ii > ndsm")
!
 920  continue
      call wexit("ntsmset","ic1, ic2 < 0")
      end
!
!**********************************************************************
      subroutine mksmdt(sc,scsm,kout)
!**********************************************************************
!
!    real*8  wflx,wssn(0:ndmc,ndgs),wssp(0:ndmc,ndgs)
!    real*8  wswe(0:ndmc),wswi(0:ndmc),wsbr(0:ndmc,ndgs)
!
!----------------------------------------------------------------------
      use cntcom, only : iplx, iply, mcel, volm
      use cntsms, only : i6_ntsms, icsm, ncsm, nsmax, wtsm
      use cplmet, only : icel, itpve, itpvs, itsle, itsls, jcel, jtmax
     >    , jtmin
      use csize,  only : ndmc
      implicit none
!
!::argument
      real(8), dimension(0:ndmc), intent(in)  :: sc
      real(8), dimension(0:ndmc), intent(out) :: scsm
      integer,                    intent(in)  :: kout
!
!::local variables
      integer  i, j, ic, ic0, jt, it, jtst, jten, i6
      real(8) :: sums, sumv, avsc
!
!
      do ic = 0, ndmc
        scsm(ic) = sc(ic)
      enddo
!
      do i = 1, nsmax
        sums = 0.0d0
        sumv = 0.0d0
        ic0 = icsm(i,1)
        do j = 1, ncsm(i)
          ic = icsm(i,j)
          sums = sums + sc(ic)*wtsm(i,j)
          sumv = sumv + volm(ic)*wtsm(i,j)
        enddo
        avsc = sums/sumv
        scsm(ic0) = avsc*volm(ic0)
      enddo
!
      if( kout.le.0 ) return
!
!::debug write
      i6 = i6_ntsms
      if( i6.gt.0 ) then
      write(i6,'(/2x,"*** mksmdt ***  kout =",i2,2x,a)')
     >  kout, "1:wssn/2:wssp/3:wswi/4:wswe"
      do it = itsls, itpve
      if( it.ne.itsle .and. it.ne.itpvs ) cycle
      jtst = jtmin(it)
      jten = jtmax(it)
      sums = 0.0d0
      sumv = 0.0d0
      do jt = jtst, jten
      j = jcel(jt,it)
      i = icel(jt,it)
      ic = mcel(j,i)
      if( ic.eq.0 ) cycle
      sums = sums + sc(ic)
      sumv = sumv + scsm(ic)
      if( jt-jtst.le.10 .or. jten-jt.le.10 ) then
      write(i6,'(2x,2i5,i7,2i5,1p4e12.3)')
     >   jt, it, ic, iplx(ic), iply(ic), sc(ic), scsm(ic), sums, sumv
      endif
      enddo
      enddo
      endif
!
      return
      end
