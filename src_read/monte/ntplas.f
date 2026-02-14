!**********************************************************************
      subroutine ntplas
!**********************************************************************
!
!      set plasma parameter
!
!  (D+,T+,C+,C2+,C3+,C4+,C5+,C6+)  nsrc = 3
!
!   isrc = 1       ngas
!     D0  ==>  D+    1   rmas(1) = 2.0*cmp = ama(1)
!              T+    2   rmas(2) = 3.0*cmp = ama(2)
!
!   isrc = 2
!     T0  ==>  D+    1   rmas(1) = 2.0*cmp = ama(1)
!              T+    2   rmas(2) = 3.0*cmp = ama(2)
!
!   isrc = 3
!     C0  ==>  C+    1   rmas(1) = 12.0*cmp = ama(3)
!
!    at present,  D0 ==> D+
!
!----------------------------------------------------------------------
      use cntcom, only : e0in, e0ps, e0wl, eips, flxin, imps, isxp
     >    , isyp, isxw, isyw, mcel, nbxp, nbyp, ncmax, ngas, noia, npep
     >    , npew, npsp, npsw, nptldp, nptlvl, nptlwl, rmas, tbsg_ei
     >    , v0in, vion, wtmin
      use cntctl, only : mdl_hrt
      use cntpls, only : dene, deni, teme, temi, vflw
      use cphcns, only : cev
      use cplcom, only : ama, nion, vte, vti
      use cplmet, only : icspx, jcdp1, jcdp2
      use csize,  only : ndgs, ndmc, ndwp
      use csonic, only : itim, lfdbg
      use cunit,  only : n6
      implicit none
!
!::local variables
      integer lon, ig, ia, nsizg, ic, icx, icy
      integer nw, iw, ixbc, iybc, j, i, iwm
      real*8  zei
      real(8) :: emti, emte, timn_emt
! function
      integer    imox, imoy
!
!::input
      timn_emt = 20.0d0
!
!::output
      lon = 0
      if( itim.eq.0 ) lon = 1
      if( lfdbg(5).ge.2 ) lon = 2
      if( lfdbg(5).eq.0 ) lon = 0
!
!-----
!::   move ntinpt                              03/4/6
!        set        ngas, rmas, e0in, wtmin
!-----
!
!::parameter set (dene,deni,teme,temi,vflw)
!::aux. parameter
      nsizg = (ndmc+1)*ndgs
      call setd( vion, nsizg, 0.0d0 )
      do ia = 1, nion
      do ic = 1, ncmax
      if( teme(ic).gt.0.0 ) then
      vion(ic,ia) = dsqrt(temi(ic)*cev/ama(ia))
      endif
      enddo
      enddo
!
!::eips
      call setd( eips, ndwp, 0.0d0 )
      call setd( e0ps, ndwp, 0.0d0 )
      do nw = 1, 4
        do iw = npsp(nw), npep(nw)-1
          j = isxp(iw) + nbxp(nw)
          i = isyp(iw) + nbyp(nw)
          zei = 0.0d0
          if( j.gt.0 .and. i.gt.0 ) then
!-----
!xx   zei = 3.0d0*vti(j,i) + 0.5d0*vte(j,i) + 3.0d0*vte(j,i)
            emti = dmin1( vti(j,i), timn_emt )
            emte = vte(j,i)
            zei  = 3.0d0*emti + 0.5d0*emte + 3.0d0*emte
!-----
          endif
          eips(iw) = zei
          if( nw.eq.1 .or. nw.eq.3) then ! for sol/prv boundary
            e0ps(iw) = 3.0d0
          else ! for diverter
            iwm = imps(iw)
            e0ps(iw) = e0wl(iwm)   ! <== Note  04/03/07
          endif
        enddo
      enddo
!
!::el-ionization
      call ntcrtb
!
!:: set H radiation trapping effects
      if(mdl_hrt > 0) call hrt_set
!
!
!::debug write
!
      if( lon.ge.1 ) then
      write(n6,'(2x,"ngas =",i2,"  noia =",10i4)')
     >     ngas,(noia(ig),ig=1,ngas)
      write(n6,'(2x,"rmas =",1p10e12.3)') (rmas(ig),ig=1,ngas)
      write(n6,'(2x,"e0in =",1pe12.3,"  v0in =",1pe12.3,"  flxin ="
     >   ,1pe12.3)') e0in,v0in,flxin
      write(n6,'(2x,"nsamp=",3i7,"  wtmin =",1pe12.3)')
     >  nptldp, nptlwl, nptlvl, wtmin
      endif
!
      if( lon.ge.2 ) then
      write(n6,'(/2x,"plasma parameter at separatrix")')
      write(n6,'(6x,"j",2x,"icx",2x,"icy",4x,"ic",3x,"dene",8x,"deni"
     >  ,8x,"teme",8x,"temi",8x,"vflw",8x,"vion",8x,"sgei")')
      icy = icspx
      do 310 j = jcdp1, jcdp2
      icx = j
      ic  = mcel(icx,icy)
      write(n6,'(2x,3i5,i6,1p7e12.3)')
     > j,icx,icy,ic,dene(ic),deni(ic,1),teme(ic),temi(ic),vflw(ic,1)
     >  ,vion(ic,1),tbsg_ei(ic)
 310  continue
      endif
!
      return
      end
