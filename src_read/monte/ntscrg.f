!**********************************************************************
      subroutine ntscrg
!**********************************************************************
!
!    soldori/cntscg
!
!      ic  : cell number
!      ia  : ion species
!      irc : reaction species (ndrmx)
!      irg : region  (10)
!
!    sn_mt(ic,ia),    sp_mt(ic,ia),    we_mt(ic),       wi_mt(ic)
!   tsn_mt(irc,irg), tsp_mt(irc,irg), twe_mt(irc,irg), twi_mt(irc,irg)
!   gsn_mt(irc),     gsp_mt(irc),     gwe_mt(irc),     gwi_mt(irc)
!   hsn_mt(irg),     hsp_mt(irg),     hwe_mt(irg),     hwi_mt(irg)
!   vsn_mt,          vsp_mt,          vwe_mt           vwi_mt
!
!    _mt : track-length estimater
!    _mc : collision estimater
!    _mf : fluid estimater   (No use)
!    _pl : in soldor (all reaction)
!
!      common /cntscr_tot/ tsn_mt, tsp_mt, twe_mt, twi_mt,
!     >                    tsn_mc, tsp_mc, twe_mc, twi_mc
!      parameter ( nwkmp_sc = ndrmx*10*8 )
!
!----------------------------------------------------------------------
      use cntcom, only : cstyp, mrgn, ncmax, ndrmx, nrcmx, volm, lbrc
      use cntmnt, only : dotn, vcsrc, vflux, vnsrc
      use cntscg, only : nwkmp_sc, tsn_mc, tsn_mt, tsp_mc, tsp_mt
     >    , twe_mc, twe_mt, twi_mc, twi_mt
      use cntwcn, only : swnrm, wssn, wssp, wswe, wswi
      use csonic, only : itim, kpcn
      use cunit,  only : lmype,lmspe, mywld, n6
      use mpi!,    only : mpi_bcast, mpi_real8, mpi_reduce, mpi_sum
      implicit none
!
!
!::local variables
      real*8, dimension(nwkmp_sc) :: wmwksc, vmwksc
!
!::track
      real*8, dimension(ndrmx) :: gsn_mt, gsp_mt, gwe_mt, gwi_mt
      real*8, dimension(10)    :: hsn_mt, hsp_mt, hwe_mt, hwi_mt
      real*8                   :: vsn_mt, vsp_mt, vwe_mt, vwi_mt
!
!::collision
      real*8, dimension(ndrmx) :: gsn_mc, gsp_mc, gwe_mc, gwi_mc
      real*8, dimension(10)    :: hsn_mc, hsp_mc, hwe_mc, hwi_mc
      real*8                   :: vsn_mc, vsp_mc, vwe_mc, vwi_mc
!
!::soldor
!xxx  real*8 wssn, wssp, wswe, wswi
      real*8, dimension(10)    :: hsn_pl, hsp_pl, hwe_pl, hwi_pl
      real*8                   :: vsn_pl, vsp_pl, vwe_pl, vwi_pl
!
      integer  :: i6, isrc, irc, irg, ic
      integer  :: ierr
      real*8   :: zfc, zsn, zsp, zwe, zwi
      integer  :: lp = 0
      save   lp
!
      lp = lp + 1
      i6 = n6
      if(lmype==lmspe)then
            i6 = 10099
            open(unit=i6,file="zscrg.txt")
      endif
!
!::[MPI_Reduce in ntscrg]  cntscg/cntscr_tot/ (tsn_mt,...)  10/12/21
!::[MPI_Bcast  in ntscrg]  cntscg/cntscr_tot/ (tsn_mt,...)  10/12/21
      vmwksc(1:nwkmp_sc) = 0.0d0
!
      call setvmwksc( 2, wmwksc )
      call MPI_Reduce( wmwksc, vmwksc, nwkmp_sc, MPI_REAL8, MPI_SUM,
     >                 lmspe, mywld, ierr )
      call MPI_Bcast( vmwksc, nwkmp_sc, MPI_REAL8, lmspe, mywld, ierr )
!
      call setvmwksc( 1, vmwksc )
!
!:debug write
      write(i6,'(/2x,"*** ntscrg ***")')
      write(i6,'(2x,"kpcn =",i3,"  itim =",i8)') kpcn, itim
      isrc = vnsrc
      write(i6,'(2x,"isrc =",i2,2x,a,2x,a,"  flxin =",1p2e12.4)')
     >     isrc,trim(vcsrc(isrc)), trim(cstyp), vflux(isrc), dotn
!
!::wssn,wssp,wswe,wswi  see plntsr
!----------------------------------------------------------------------
      hsn_pl(1:10) = 0.0d0
      hsp_pl(1:10) = 0.0d0
      hwe_pl(1:10) = 0.0d0
      hwi_pl(1:10) = 0.0d0
!
      do ic = 1, ncmax
      irg = mrgn(ic)
      zfc = dotn/swnrm/volm(ic)
      if( irg.le.0 .or. volm(ic).le.0.0d0 ) goto 910
      zsn = zfc*wssn(ic,1)
      zsp = zfc*wssp(ic,1)
      zwe = zfc*wswe(ic)
      zwi = zfc*wswi(ic)
      hsn_pl(irg) = hsn_pl(irg) + zsn*volm(ic)
      hsp_pl(irg) = hsp_pl(irg) + zsp*volm(ic)
      hwe_pl(irg) = hwe_pl(irg) + zwe*volm(ic)
      hwi_pl(irg) = hwi_pl(irg) + zwi*volm(ic)
      enddo
!
      vsn_pl = 0.0d0
      vsp_pl = 0.0d0
      vwe_pl = 0.0d0
      vwi_pl = 0.0d0
!
      do irg = 1, 7
      vsn_pl = vsn_pl + hsn_pl(irg)
      vsp_pl = vsp_pl + hsp_pl(irg)
      vwe_pl = vwe_pl + hwe_pl(irg)
      vwi_pl = vwi_pl + hwi_pl(irg)
      enddo
!
!::track
!----------------------------------------------------------------------
      gsn_mt(1:nrcmx) = 0.0d0
      gsp_mt(1:nrcmx) = 0.0d0
      gwe_mt(1:nrcmx) = 0.0d0
      gwi_mt(1:nrcmx) = 0.0d0
      hsn_mt(1:10) = 0.0d0
      hsp_mt(1:10) = 0.0d0
      hwe_mt(1:10) = 0.0d0
      hwi_mt(1:10) = 0.0d0
!
      zfc = dotn/swnrm
      do irc = 1, nrcmx
      do irg = 1, 7
      gsn_mt(irc) = gsn_mt(irc) + zfc*tsn_mt(irc,irg)
      gsp_mt(irc) = gsp_mt(irc) + zfc*tsp_mt(irc,irg)
      gwe_mt(irc) = gwe_mt(irc) + zfc*twe_mt(irc,irg)
      gwi_mt(irc) = gwi_mt(irc) + zfc*twi_mt(irc,irg)
!
      hsn_mt(irg) = hsn_mt(irg) + zfc*tsn_mt(irc,irg)
      hsp_mt(irg) = hsp_mt(irg) + zfc*tsp_mt(irc,irg)
      hwe_mt(irg) = hwe_mt(irg) + zfc*twe_mt(irc,irg)
      hwi_mt(irg) = hwi_mt(irg) + zfc*twi_mt(irc,irg)
      enddo
      enddo
!
      vsn_mt = 0.0d0
      vsp_mt = 0.0d0
      vwe_mt = 0.0d0
      vwi_mt = 0.0d0
!
      do irg = 1, 7
      vsn_mt = vsn_mt + hsn_mt(irg)
      vsp_mt = vsp_mt + hsp_mt(irg)
      vwe_mt = vwe_mt + hwe_mt(irg)
      vwi_mt = vwi_mt + hwi_mt(irg)
      enddo
!
!::collision
!----------------------------------------------------------------------
      gsn_mc(1:nrcmx) = 0.0d0
      gsp_mc(1:nrcmx) = 0.0d0
      gwe_mc(1:nrcmx) = 0.0d0
      gwi_mc(1:nrcmx) = 0.0d0
      hsn_mc(1:10) = 0.0d0
      hsp_mc(1:10) = 0.0d0
      hwe_mc(1:10) = 0.0d0
      hwi_mc(1:10) = 0.0d0
!
      zfc = dotn/swnrm
      do irc = 1, nrcmx
      do irg = 1, 7
      gsn_mc(irc) = gsn_mc(irc) + zfc*tsn_mc(irc,irg)
      gsp_mc(irc) = gsp_mc(irc) + zfc*tsp_mc(irc,irg)
      gwe_mc(irc) = gwe_mc(irc) + zfc*twe_mc(irc,irg)
      gwi_mc(irc) = gwi_mc(irc) + zfc*twi_mc(irc,irg)
!
      hsn_mc(irg) = hsn_mc(irg) + zfc*tsn_mc(irc,irg)
      hsp_mc(irg) = hsp_mc(irg) + zfc*tsp_mc(irc,irg)
      hwe_mc(irg) = hwe_mc(irg) + zfc*twe_mc(irc,irg)
      hwi_mc(irg) = hwi_mc(irg) + zfc*twi_mc(irc,irg)
      enddo
      enddo
!
      vsn_mc = 0.0d0
      vsp_mc = 0.0d0
      vwe_mc = 0.0d0
      vwi_mc = 0.0d0
!
      do irg = 1, 7
      vsn_mc = vsn_mc + hsn_mc(irg)
      vsp_mc = vsp_mc + hsp_mc(irg)
      vwe_mc = vwe_mc + hwe_mc(irg)
      vwi_mc = vwi_mc + hwi_mc(irg)
      enddo
!
!::output
!----------------------------------------------------------------------
!::Sn
      write(i6,'(/2x,"--- Sn_tr -----------------------")')
      write(i6,603)"irc","rct","tot",
     >              "odp","sol","idp","opv","ipv","edg","man"
      do irc = 1, nrcmx
      write(i6,602) irc,lbrc(irc),gsn_mt(irc),
     >                  (zfc*tsn_mt(irc,irg),irg=1,7)
      enddo
      write(i6,602) 0,"---",vsn_mt,(hsn_mt(irg),irg=1,7)
!
      write(i6,'(/2x,"--- Sn_cl -----------------------")')
      write(i6,603)"irc","rct","tot",
     >              "odp","sol","idp","opv","ipv","edg","man"
      do irc = 1, nrcmx
      write(i6,602) irc,lbrc(irc),gsn_mc(irc),
     >                  (zfc*tsn_mc(irc,irg),irg=1,7)
      enddo
      write(i6,602) 0,"---",vsn_mc,(hsn_mc(irg),irg=1,7)
!
      write(i6,'(/2x,"--- Sn_pl -----------------------")')
      write(i6,603)"irc","rct","tot",
     >              "odp","sol","idp","opv","ipv","edg","man"
      write(i6,602) 0,"---",vsn_pl,(hsn_pl(irg),irg=1,7)
!
!::We
      write(i6,'(/2x,"--- We_tr -----------------------")')
      write(i6,603)"irc","rct","tot",
     >              "odp","sol","idp","opv","ipv","edg","man"
      do irc = 1, nrcmx
      write(i6,602) irc,lbrc(irc),gwe_mt(irc),
     >                  (zfc*twe_mt(irc,irg),irg=1,7)
      enddo
      write(i6,602) 0,"---",vwe_mt,(hwe_mt(irg),irg=1,7)
!
      write(i6,'(/2x,"--- We_cl -----------------------")')
      write(i6,603)"irc","rct","tot",
     >              "odp","sol","idp","opv","ipv","edg","man"
      do irc = 1, nrcmx
      write(i6,602) irc,lbrc(irc),gwe_mc(irc),
     >                  (zfc*twe_mc(irc,irg),irg=1,7)
      enddo
      write(i6,602) 0,"---",vwe_mc,(hwe_mc(irg),irg=1,7)
!
      write(i6,'(/2x,"--- We_pl -----------------------")')
      write(i6,603)"irc","rct","tot",
     >              "odp","sol","idp","opv","ipv","edg","man"
      write(i6,602) 0,"---",vwe_pl,(hwe_pl(irg),irg=1,7)
!
!
!::Wi
      write(i6,'(/2x,"--- Wi_tr -----------------------")')
      write(i6,603)"irc","rct","tot",
     >              "odp","sol","idp","opv","ipv","edg","man"
      do irc = 1, nrcmx
      write(i6,602) irc,lbrc(irc),gwi_mt(irc),
     >                  (zfc*twi_mt(irc,irg),irg=1,7)
      enddo
      write(i6,602) 0,"---",vwi_mt,(hwi_mt(irg),irg=1,7)
!
      write(i6,'(/2x,"--- Wi_cl -----------------------")')
      write(i6,603)"irc","rct","tot",
     >              "odp","sol","idp","opv","ipv","edg","man"
      do irc = 1, nrcmx
      write(i6,602) irc,lbrc(irc),gwi_mc(irc),
     >                  (zfc*twi_mc(irc,irg),irg=1,7)
      enddo
      write(i6,602) 0,"---",vwi_mc,(hwi_mc(irg),irg=1,7)
!
      write(i6,'(/2x,"--- Wi_pl -----------------------")')
      write(i6,603)"irc","rct","tot",
     >              "odp","sol","idp","opv","ipv","edg","man"
      write(i6,602) 0,"---",vwi_pl,(hwi_pl(irg),irg=1,7)
!
!
 602  format(2x,i2,2x,a3,2x,1pe12.3,2x,1p7e12.3)
 603  format(1x,a3,2x,a3,5x,a3,2x,7(9x,a3))
!
      if(lmype==lmspe) close(i6)
!      if( lp.gt.5 ) call wexit("ntscrg","due to check")
      return
!
 910  continue
      call wexit("ntscrg","invalid irg, volm(ic)")
      end
