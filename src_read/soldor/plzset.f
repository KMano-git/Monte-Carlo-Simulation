!**********************************************************************
      subroutine plzset
!**********************************************************************
!
! cimcom:  common /cimcom_1/ aimas, azmas, ami, amz, ismax
!    --> climp
!
!       mion      : N  include impurity
!       wai(mdsp) : Ai
!       wza(mdsp) : Zi
!       wma(mdsp) : Mi
!
!        H+ D+ C+ C2+ C3+ C4+ C5+ C6+
!         nion = 2,  mion = 2 + 6 = 8 < mdsp (10)
!
!   collision frequency
!       wfeb(mdsp), wfab(mdsp,mdsp), wkab(mdsp,mdsp), wfrab(mdsp,mdsp)
!
!----------------------------------------------------------------------
      use cphcns, only : cmp
      use cplcom, only : aion, aza, nion
      use cplimp, only : aimasl, amil, amzl, azmasl, ismaxl, wmc_nty
      use cmeffz, only : dzan, mion, vdnz, vsez, vsez0, vsezg, wai, wfab
     >    , wfeb, wkab, wma, wza, xzflz, xzwtm, xzwtp
      use csize,  only : ndmis, mdsp, ndx, ndy, nzsmx
      use csonic, only : limp
      use cunit,  only : n6
      implicit none
!
!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  ia, iz, ii, ib, ic, nty
      integer  ia, iz, ii, ib, nty
!
      write(n6,'(/2x,"*** plzset (NEW) ***  limp =",i2)') limp
      if( limp.eq.0 ) then
        ismaxL = 0
      endif
!
      mion = nion
      do nty =1, wmc_nty
      write(n6,'(2x,i3," aimas =",f8.4,"  azmas =",f8.4,"  ismax =",i3,
     >  "  ami =",1pe12.4,"  amz =",1pe12.4)')
     >  nty, aimasL(nty), azmasL(nty), ismaxL(nty), amiL(nty), amzL(nty)
!
      mion = mion + ismaxL(nty)
      enddo
!
      write(n6,'(2x,"bf: nion =",i2,"  mion =",i2,"  mdsp =",i3)')
     >   nion, mion, mdsp
!
      call flush(n6)
      if( mion.gt.mdsp ) call wexit("plzset","mion.gt.mdsp")
!
      ii = 0
      do ia = 1, nion
      ii = ii + 1
      wai(ii) = aion(ia)
      wza(ii) = aza(ia)
      wma(ii) = aion(ia)*cmp
      enddo
      mion = ii
!
!::impurity (Argon)
      do nty = 1, wmc_nty
      if( ismaxL(nty).gt.0 ) then
      do iz = 1, ismaxL(nty)
      ii = ii + 1
      wai(ii) = azmasL(nty)
      wza(ii) = dfloat(iz)
      wma(ii) = azmasL(nty)*cmp
      enddo
      mion = ii
      endif
      enddo
!
!::debug write
      write(n6,'(2x,"af: nion =",i2,"  mion =",i2)') nion, mion
      write(n6,'(2x,"ia",3x,"Ai",10x,"Zi",10x,"Mi")')
      do ib = 1, mion
      write(n6,'(2x,i2,1p3e12.3)')
     >   ib, wai(ib), wza(ib), wma(ib)
      enddo
!
!::extesion of sub. plpset
      do ia = 1, mion
      wfeb(ia) = wza(ia)**2/3.4411d11
      do ib = 1, mion
      wfab(ia,ib) = wza(ia)**2*wza(ib)**2
     >  /(sqrt((wai(ia)+wai(ib))/2.0d0)*2.0853d13)
      wkab(ia,ib) = wfab(ia,ib)*sqrt(wai(ia)*wai(ib))
      enddo
      enddo
!
!::friction force for plasma ion
!     cfrab(ia,ib)
!
!::clear  /cmzeflx/
! modified 1/1 lines dynamic allocation of arrays by kamata 2022/05/29
!ik   vdnz (1:ndx,1:ndy,0:ndis2L,1:nzsmx) = 0.0d0
      vdnz (1:ndx,1:ndy,0:ndmis,1:nzsmx) = 0.0d0
      vsez (1:ndx,1:ndy) = 0.0d0
      vsez0(1:ndx,1:ndy) = 0.0d0
      vsezg(1:ndx,1:ndy) = 0.0d0
      xzwtm(1:ndx,1:ndy) = 0.0d0
      xzwtp(1:ndx,1:ndy) = 0.0d0
      xzflz(1:ndx,1:ndy) = 0.0d0
      dzan (1:ndx,1:ndy) = 0.0d0
!
      return
      end
!
!**********************************************************************
      subroutine dbg_ama(cmsg)
!**********************************************************************
      use cplcom, only : ama, nion
      use cmeffz, only : mion
      use cunit,  only : n6
      implicit none
!
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   character(*) :: cmsg
      character(*), intent(in) :: cmsg
!
      write(n6,'(2x,"### dbg_ama  ",a,2x,2i4,1p10e12.3)')
     >    trim(cmsg),nion,mion,ama(1:nion)
!
      return
      end
