!**********************************************************************
      subroutine plwrdMC
!**********************************************************************
!
!    2017/1/18 K.Hoshino
!    based on plcrad_mcN.f rev.302
!
!    energy loss due to interactions with impurity
!
!::Note
!     IMPMC/calwrd.f
!         twrd(ic)  iz + il + rc
!         twci(ic)  collision with ions
!
!----------------------------------------------------------------------
      use catcom,  only : ip_ion, ip_plt, ip_prb
      use cntcom,  only : iplx, iply, mrgn, ncmax
      use cphcns,  only : cev
      use cplcom,  only : anamp, atemp, vnezef, vte
      use cplimp,  only : eipotl, nzmxl, tdnzl, twcil, wmc_nty
      use cplwrd2, only : sradi, sradli, sradr, wmc_zwe, wmc_zwi
      use csize,   only : ndmis
      implicit none

      real*8 :: zwiz, zwil, zwrc
      integer  :: nty, ic, j, i, iz,  isno
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  :: ir, ir2
      integer  :: ir
      real*8   :: hne, hte, hwne
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   real*8,dimension(ndis2L) :: zLz_i, zLz_r, zLz_cx, zsg_i, zsg_r
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/05/08
!ik   real*8,dimension(ndis2L) :: zLz_i, zLz_r, zsg_i
      real*8,dimension(ndmis) :: zLz_i, zLz_r, zsg_i
!
!::evaluate radiation using IMPMC
!
! returned 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
      do nty = 1, wmc_nty
!ik   nty = 1
!
      do ic = 1, ncmax  ! KSFUJI
        j = iplx(ic)
        i = iply(ic)
        ir  = mrgn(ic)
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik     ir2 = mrgnp(ic)
        if( j.le.0 .or. i.le.0 ) cycle
!
!::hot core
        if( ir.eq.7 ) then
          hne = anamp(i,1)
          hte = atemp(i)
        else
          !hne = vne(j,i)
          hne = vnezef(j,i)
          hte = vte(j,i)
        endif
!
!::ele density
!xx   hwne = twne(ic) ! el density at cal. of IMPMC
        hwne = hne      ! the last el density
!
!::ionization
        zwiz = 0.0d0
        zwil = 0.0d0
! modified 2/2 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik     call atm_eval2(ip_plt,hwne,hte,zLz_i,isno,ndis2L)
!ik     call atm_eval2(ip_ion,hwne,hte,zsg_i,isno,ndis2L)
! modified 2/2 lines treat 4 or more impurities with IMPMC by kamata 2022/05/08
!ik     call atm_eval2(ip_plt(nty),hwne,hte,zLz_i,isno,ndis2L,nty)
!ik     call atm_eval2(ip_ion(nty),hwne,hte,zsg_i,isno,ndis2L,nty)
        call atm_eval2(ip_plt(nty),hwne,hte,zLz_i,isno,ndmis,nty)
        call atm_eval2(ip_ion(nty),hwne,hte,zsg_i,isno,ndmis,nty)
        do iz = 0, nzmxL(nty)-1
          sradi(iz,ic,nty) = hwne*tdnzL(iz,ic,nty)*zLz_i(iz+1)
          zwiz = zwiz + hwne*tdnzL(iz,ic,nty)*zLz_i(iz+1)
          sradli(iz,ic,nty) =
     >    hwne*tdnzL(iz,ic,nty)*zsg_i(iz+1)*eipotL(iz,nty)*cev
          zwil = zwil
     >            + hwne*tdnzL(iz,ic,nty)*zsg_i(iz+1)*eipotL(iz,nty)*cev
        enddo
!
!::recombinatin  (Note 2010/01/17)
        zwrc = 0.0d0
        if( hte.gt.0.8d0 ) then
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik       call atm_eval2(ip_prb,hwne,hte,zLz_r,isno,ndis2L)
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/05/08
!ik       call atm_eval2(ip_prb(nty),hwne,hte,zLz_r,isno,ndis2L,nty)
          call atm_eval2(ip_prb(nty),hwne,hte,zLz_r,isno,ndmis,nty)
          do iz = 1, nzmxL(nty)
            sradr(iz,ic,nty) =  hwne*tdnzL(iz,ic,nty)*zLz_r(iz)
            zwrc = zwrc + hwne*tdnzL(iz,ic,nty)*zLz_r(iz)
          enddo
        endif
!
!::energy loss due to interactions with impurity
!xxx  wmce(ic) = -twrd(ic) ! IMPMC result
      wmc_zwe(ic, nty) = -( zwiz + zwil + zwrc )  ! recal. new Te
      wmc_zwi(ic, nty) = twciL(ic,nty)  ! see imscatt
!
      enddo  ! ic
! returned 1 line treat 4 or more impurities with IMPMC by kamata 2022/04/15
      enddo  ! nty
!
      return
      end
