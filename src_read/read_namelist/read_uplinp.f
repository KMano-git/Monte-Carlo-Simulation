!**********************************************************************
      subroutine read_uplinp(nft, dtq1rg, dtq2rg, dtq3rg, dtq4rg
     >  ,rgsize, dtq1dp, dtq2dp, dtq3dp, dtq4dp, dpsize
     >  ,xplmt, xplmt_len, timeNum, interNum, vlda)
!**********************************************************************
      use cplcom, only : aion, aza, beta, bmpfn, bmpni, bmpte, bmpti
     >    , bwpfn, bwpni, bwpte, bwpti, bwsfn, bwsni, bwste, bwsti, caim
     >    , cdle, cdli, cevpr, cfth, cftp, cftr, cftw, chfcv, chfpr
     >    , cimp, cimprg, cimprg2, clmdnz, clzef, cprm, dtfn, dtlmt
     >    , dtmax, dtmin, dtrateq, dttb, edtb, edtmn, edtmx, elpmx, eltb
     >    , exni, exte, exti, exvp, factor_bal, fai, faim, faip, fcda
     >    , fcet, fcna, fcbcvl, fdps, fdsl2, flime, flimi, flimv, gcse
     >    , gcse_pv, gcse_sl, gcsi, gcsi_pv, gcsi_sl, gwpni, gwpte
     >    , gwpti, gwsni, gwste, gwsti, itfix, itsl2, ittb
     >    , jen_bal, jst_bal, lbcgd, lbcpw, lbcsw, ldluc, ldps, lfbbc
     >    , lmuscl, lordr, lttb, lttbr, lupwd, mdl_cgen, mdl_bal
     >    , mdl_bale, mdl_eqp, mdl_edt, mdl_fimp, mdl_hcv, mdl_ini
     >    , mdl_srp, mdl_srw, mdl_vis, mdl_wrd, nad0, nas0, nfth, nftp
     >    , nftr, nftw, nimin_aux, nion, nlpmn, nlpmx, nsmp_wrad, ntfn
     >    , nttb, rcydt, rcysp, relmt, rxflp, ted0, temin_aux, temn_dprv
     >    , temn_dpsl, temndf, tes0, tid0, timin_aux, timin_vis
     >    , timndf, timneta, tis0, vlet, vlxe, vlxi, wfac_lv
     >    , fdeg2
      use cplmet, only : set_hdsv
      use cplqcn, only : itqcn
      use cplvpn, only : vpn_man, vpn_sol
      use cplwrd, only : wcr_cnc, wcr_nty, wcr_typ
      use csonic, only : itend, mxdsk, mxhst, mxprf
     >    , ndsk, nhsav, nhst, npcn, nprf, tend, time
      use dbg_mod, only : dtmdbgsol, timdbgsol
      use cimcom, only : mrgn_max
      implicit none
!
!::argument
      ! nft: number for read namelist file
      integer, intent(in) :: nft, rgsize, dpsize, xplmt_len
      real(8), dimension(rgsize), intent(out) :: dtq1rg, dtq2rg
     >  , dtq3rg, dtq4rg
      real(8), dimension(dpsize), intent(out) :: dtq1dp, dtq2dp
     >  , dtq3dp, dtq4dp
      character(xplmt_len), intent(out) :: xplmt
      integer, intent(out) :: timeNum, interNum
      real(8), intent(out) :: vlda(mrgn_max)
!
!::local variable
      integer  nfmt ! dummy
!
      namelist /uplinp/
     >   caim, cprm
     >  ,nion, aion, aza
     >  ,exni, exvp, exti, exte, fcna
     >  ,nfmt
     >  ,nftr, nftw, nfth, nftp, cftr, cftw, cfth, cftp
     >  ,ndsk, nhst, nhsav, nprf, mxdsk, mxhst, mxprf, npcn
     >  ,itend, itfix, tend
     >  ,dtmax, dtmin, dtlmt
!-----
     >  ,dtq1rg, dtq2rg, dtq3rg, dtq4rg
     >  ,dtq1dp, dtq2dp, dtq3dp, dtq4dp
     >  ,itsl2,  fdsl2
     >  ,fdeg2 ! diffusion enhance factor for edge
!-----
     >  ,nttb, lttb, lttbr, dttb, eltb, edtb, dtfn
     >  ,nlpmn, nlpmx, elpmx, edtmx, edtmn, relmt
     >  ,vlda, vlet, vlxi, vlxe, fcda, fcet
     >  ,ldps, fdps, ldluc
     >  ,cdli, cdle
     >  ,lbcgd, lbcsw, lbcpw, lordr, lupwd, lmuscl
     >  ,bwsni, bwsti, bwste, bwsfn, bwpni, bwpti, bwpte, bwpfn
     >  ,bmpni, bmpti, bmpte, bmpfn
     >  ,gwsni, gwsti, gwste, gwpni, gwpti, gwpte
     >  ,clzef
     >  ,gcsi_sl, gcsi_pv, gcse_sl, gcse_pv
     >  ,cimp, cimprg, cimprg2, rxflp, rcydt, rcysp, clmdnz

     >  ,wcr_nty, wcr_typ, wcr_cnc

     >  ,itqcn, xplmt
     >  ,vpn_sol, vpn_man
     >  ,temndf, timndf, timneta, temn_dpsl, temn_dprv
     >  ,mdl_wrd, mdl_eqp, mdl_edt, mdl_srw, mdl_srp, mdl_vis
     >  ,mdl_ini, mdl_hcv, mdl_cgen, mdl_fimp
     >  ,temin_aux, timin_aux, timin_vis, nimin_aux
     >  ,nas0, tis0, tes0, nad0, tid0, ted0
     >  ,flimi, flime, flimv, fcbcvl
     >  ,lfbbc
     >  ,nsmp_wrad, wfac_lv

     >  ,mdl_bal
     >  ,jst_bal, jen_bal, factor_bal
     >  ,mdl_bale
!----- time Series of Qdpl_i,o.txt
     >  , timeNum, interNum
!----- periodic invocation of IMPMC/figdat
     >  , dtmdbgsol, timdbgsol
!----- set anomalous diffusion for poloidal direction 
     >  , set_hdsv  
!
      rewind(nft)
      read(nft,uplinp)
      end subroutine