!***********************************************************************
      subroutine impdt(nft,ip,comt,ptmx,ptim,pdtm)
!***********************************************************************
!
!       i6   : output file
!       comt : message
!       ptmx : time to be traced untill t = ptmx
!       ptim : current time
!       pdtm : birth time of neutral particle
!            : time step of ion   ( pdunt)
!
!::Header @trace_1
!
!    comt lpcm   ip  ptmx       ptim       tt         pdtm        ic
!  ix  iy    ln ien  ko  is ml   rr       zz      rO     wght   vz
!     Evel      vflw      Ti        lstp      tauz      dt
!    ftot       dvz       v
!
!-----------------------------------------------------------------------
      use cimcom, only : aimas, azmas, ami, amz, cfez, cfgte, cfgti
     >    , cstrg, gdez, gdte, gdti, glte, glti, ien, igi, il, iml, ir
     >    , is, lfgti, lfgtzf, lpcm, rr, tt, v, vv, vvr, vz, wght, zz
      use cimntl, only : svx, svy, svz, sx0, sy0, sz0
      use cntpls, only : temi, vflw, zefm
      use cphcns, only : cev
      use cunit,  only : n6
      implicit none
!
!::argument
! modified 3/3 lines organize local variables and include files by kamata 2021/06/28
!ik   integer   nft, ip
!ik   character comt*(*)
!ik   real*8    ptmx, ptim, pdtm
      integer,   intent(in) :: nft, ip
      character, intent(in) :: comt*(*)
      real(8),   intent(in) :: ptmx, ptim, pdtm
!
!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/06/28
!ik   integer  ic, ko, lstp, imox, imoy, ln, i6, ierr, iz
      integer  ic, ko, lstp, ln, i6, ierr, iz
      real*8   tauz, dtstp
! modified 3/3 lines organize local variables and include files by kamata 2021/06/28
!ik   real*8   wtm, xi, yi, zi, ri, vlx, vly, vlz, zvl2, zvlb
!ik   real*8   funroh, zro
!ik   real*8   hdt, hftot, hdvz
      real*8   wtm, xi, yi, zi, ri, zvl2, zvlb
      real*8   zro
      real*8   hdvz
!
!::local variables (ftot)
      real(8) :: ftot
! modified 1/1 lines organize local variables and include files by kamata 2021/06/28
!ik   real(8) :: ckp0, zvf, zti, zvi2, xx2, xxp2, zbet, zbet0
      real(8) :: ckp0, zvf, zti, zvi2, xx2, xxp2, zbet
      real(8) :: zcti1, grdti, grdte
      real(8) :: za, zfca, zcti, zcte
      integer :: ia
      integer :: nlp = 0
      save :: nlp
! added 3 lines organize local variables and include files by kamata 2021/06/28
! function
      real(8)    funroh
      integer    imox, imoy
!
!xx   return
      if( nft <= 0 ) return
      nlp = nlp + 1
      if( nlp > 500 ) return
      i6 = nft
!
!::ion
      if( is(ip).ne.0 ) then
      ic = ir(ip)
      ln = il(ip)  ! <== 2010/05/12
      iz = is(ip)
      za = iz
!---------------------------------------------------------------------
!::Fz = electric field + e-thermal + i-thermal force
      ftot = 0.0d0
      if( temi(ic).gt.0.0d0 ) then
      ckp0 = 1.0d0/0.5657d0
      ia  = 1
      zvf = vflw(ic,ia)
      zti = temi(ic)
      zvi2 = 2.0d0*cev*zti/ami
      xx2  = (vvr(ip)+(vz(ip)-zvf)**2)/zvi2
      xxp2 = (vz(ip)-zvf)**2/zvi2
      zbet = 1.5d0*ckp0*(1.0d0-2.0d0*xxp2)*exp(-xx2)
      zcti1 = (1.0d0+aimas/azmas)*za**2*zbet
! deleted 1 line organize local variables and include files by kamata 2021/06/28
!ik   zbet0 = cfgti(iz)/((1.0D0+aimas/azmas)*za**2)
!
      zfca = 1.0d0
      if(lfgtzf.eq.1)  zfca = 1.0d0/za**2*(za/zefm(ic)-1.0d0)*za
      zcti = cfgti(iz)*zfca
      zcte = cfgte(iz)*zfca
      if( lfgti.eq.1 ) zcti  = zcti1*zfca
!
!::new limit of temperature gradient
      grdti = dmin1( glti(ic), dabs(gdti(ic)) )
      grdti = dsign( grdti, gdti(ic) )
      grdte = dmin1( glte(ic), dabs(gdte(ic)) )
      grdte = dsign( grdte, gdte(ic) )
      ftot  = cfez(iz)*gdez(ic)+zcte*grdte+zcti*grdti
      endif
!---------------------------------------------------------------------
!-----
! modified 3/1 lines organize local variables and include files by kamata 2021/06/28
!ik   hdt   = pdtm
!ik   hftot = ftot
!ik   hdvz  = hdt*cev/amz*hftot
      hdvz  = pdtm*cev/amz*ftot
!-----
      call mchkin(rr(ip),zz(ip),ic,ko)
      call subdvt(ip,pdtm,tauz,dtstp,lstp)
      zro = funroh(ic,rr(ip),zz(ip))
!
      if( i6.gt.0 ) then
! modified 1/1 lines with TOPICS by kamata 2021/12/22
!ik   write(i6,'(1x,a6,1x,i5,i5,1pe10.2,1p2e12.4,1pe10.2,
      write(i6,'(1x,a6,1x,i10,i5,1pe10.2,1p2e12.4,1pe10.2,
     >  i6,2i4,i6,3i4,i3,
     >  0p2f9.4, 0p2f7.3, 1p6e10.2, 1p4e10.2)')
     >  comt, lpcm, ip, ptmx, ptim, tt(ip), pdtm,
     >  ic, imox(ic), imoy(ic), ln, ien(ip), ko, is(ip), iml(ip),
     >  rr(ip), zz(ip), zro, wght(ip), vz(ip), 0.5d0*amz*vv(ip)/cev,
     >  vflw(ic,1), temi(ic), dfloat(lstp), tauz
! modified 1/1 lines organize local variables and include files by kamata 2021/06/28
!ik  > ,hdt, hftot, hdvz, v(ip)
     > ,pdtm, ftot, hdvz, v(ip)
      else
! modified 1/1 lines with TOPICS by kamata 2021/12/22
!ik   write(cstrg,'(1x,a6,1x,i5,i5,1pe10.2,1p2e12.4,1pe10.2,
      write(cstrg,'(1x,a6,1x,i10,i5,1pe10.2,1p2e12.4,1pe10.2,
     >  i6,2i4,i6,3i4,i3,
     >  0p2f9.4, 0p2f7.3, 1p6e10.2, 1p4e10.2)')
     >  comt, lpcm, ip, ptmx, ptim, tt(ip), pdtm,
     >  ic, imox(ic), imoy(ic), ln, ien(ip), ko, is(ip), iml(ip),
     >  rr(ip), zz(ip), zro, wght(ip), vz(ip), 0.5d0*amz*vv(ip)/cev,
     >  vflw(ic,1), temi(ic), dfloat(lstp), tauz
! modified 1/1 lines organize local variables and include files by kamata 2021/06/28
!ik  > ,hdt, hftot, hdvz, v(ip)
     > ,pdtm, ftot, hdvz, v(ip)
      endif
!
!::neutral
      else
      ic = ir(ip)
      ln = il(ip)
      wtm = ptim - pdtm
      if( ln.lt.0 ) wtm = wtm*0.9999d0
      xi = sx0(ip) + wtm*svx(ip)
      yi = sy0(ip) + wtm*svy(ip)
      zi = sz0(ip) + wtm*svz(ip)
      ri = sqrt(xi**2+yi**2)
      call mchkin(ri,zi,ic,ko)
      zvl2 = svx(ip)**2+svy(ip)**2+svz(ip)**2
      call velion(xi,yi,zi,svx(ip),svy(ip),svz(ip),zvlb,ierr)
      zro = funroh(ic,ri,zi)
!
      if( ierr.ne.0 ) ien(ip) = 9
!
      if( i6.gt.0 ) then
! modified 1/1 lines with TOPICS by kamata 2021/12/22
!ik   write(i6,'(1x,a6,1x,i5,i5,1pe10.2,1p2e12.4,1pe10.2,
      write(i6,'(1x,a6,1x,i10,i5,1pe10.2,1p2e12.4,1pe10.2,
     >  i6,2i4,i6,3i4,i3,
     >  0p2f9.4, 0p2f7.3, 1p4e10.2, 2x,i6)')
     >  comt, lpcm, ip, ptmx, ptim, tt(ip), pdtm,
     >  ic, imox(ic), imoy(ic), ln, ien(ip), ko, is(ip), iml(ip),
     >  ri, zi, zro, wght(ip), zvlb, 0.5d0*amz*zvl2/cev,
     >  vflw(ic,1), temi(ic), igi(ip)
      else
! modified 1/1 lines with TOPICS by kamata 2021/12/22
!ik   write(cstrg,'(1x,a6,1x,i5,i5,1pe10.2,1p2e12.4,1pe10.2,
      write(cstrg,'(1x,a6,1x,i10,i5,1pe10.2,1p2e12.4,1pe10.2,
     >  i6,2i4,i6,3i4,i3,
     >  0p2f9.4, 0p2f7.3, 1p4e10.2, 2x,i6)')
     >  comt, lpcm, ip, ptmx, ptim, tt(ip), pdtm,
     >  ic, imox(ic), imoy(ic), ln, ien(ip), ko, is(ip), iml(ip),
     >  ri, zi, zro, wght(ip), zvlb, 0.5d0*amz*zvl2/cev,
     >  vflw(ic,1), temi(ic), igi(ip)
      endif
!
      endif
!
      call flush(n6)
      return
      end
!
!***********************************************************************
      subroutine subdvt(ip,dtunt,tauz,dtstp,lstp)
!***********************************************************************
      use cimcom, only : aimas, azmas, ir, is, slnv, slw0
      use cntpls, only : dene
      use cunit,  only : n6
      implicit none
!
!::argument
! modified 2/4 lines organize local variables and include files by kamata 2021/06/28
!ik   integer  ip, lstp
!ik   real*8   dtunt, tauz, dtstp
      integer, intent(in)  :: ip
      integer, intent(out) :: lstp
      real(8), intent(in)  :: dtunt
      real(8), intent(out) :: tauz, dtstp
!
!::local common
      real*8  facdt, cftz
! deleted 1 line replace all include files with module files by kamata 2021/08/18
!ik   common /com_subdvt/ facdt, cftz
!
!::local variables
      integer ic
      real*8  zne
!
      ic = ir(ip)
      zne = 0.0d0
      if( ic.gt.0 .and. is(ip).gt.0 ) zne = dene(ic)
!
      if( zne.le.0.0d0 ) then
        lstp = 1
        tauz = 1.0d30
      else
        facdt = 1.0d0/20.0d0
        cftz = 1.0d0/dsqrt(aimas)*azmas/(azmas+aimas)
        tauz = 1.0d0/(is(ip)**2*slw0(ic,2)*slnv(ic,2))*cftz
        dtstp = tauz*facdt
!------
        if( dtstp.eq.0.0d0 ) then
        write(n6,'(2x,"ip, is, ic, slw0, slnv, cftz =",
     >   i5,2x,i3,2x,i7,2x,1p3e12.3)') ip,is(ip),ic,slw0(ic,2),
     >   slnv(ic,2),cftz
        call wexit("impdt","zero divide")
        endif
!-------
        lstp = int(dtunt/dtstp)
        if( lstp.lt.1 ) lstp = 1
      endif
      dtstp = dtunt/dfloat(lstp)
!
      return
      end
