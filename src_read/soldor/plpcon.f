!***********************************************************************
      subroutine plpcon
!***********************************************************************
!
!                  N+1  N
!      R = -V/dt*(q1a-wq1a)  - Flux + Source
!
!               qvl_dt       qfx_cv   ssnc+q1a*ssnv
!                            qfx_df
!                            qfy_df
!
!                            psm_ion  psm_ntl
!-----------------------------------------------------------------------
      use cplcom, only : ama, dtrateq, nion, nlp, q1a, ssnc, ssnv, tfvl
     >    , vna, wq1a
      use cplmet, only : hvol, icel, itmax, itmpe, itpve, jcel, jtmax
      use cplqcn, only : pcn_dtm, pcn_itm, pcn_nlp, pcn_pdt, pcn_pdt2
     >    , pcn_pfl, pcn_psi, pcn_ptb, pcn_ptb2, pcn_ptf, pcn_ptl
     >    , pcn_ptl2, pcn_rhs, pcn_tim, qfx_cv, qfx_df, qfy_df, qvl_dt
      use csize,  only : ndmfl, ndx, ndy
      use csonic, only : dtim, itim, time
      implicit none
!
!::local variables
      integer ia, j, i, m1a, it, jmax, jw
      real*8  zsum, zsum1, zsum2, zsum3, zsum4, sg
      real*8  zsum1p, zsum2p, zsum3p, zvlm
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   integer lp; data lp/0/
!
!::q(t_N-1) = q(t_N)  for N = 0
      if( itim.eq.0 ) then
      do ia = 1, nion
      do j  = 1, ndx
      do i  = 1, ndy
      wq1a(j,i,ia) = q1a(j,i,ia)
      enddo
      enddo
      enddo
      endif
!
!::clear
      pcn_ptl = 0.0d0; pcn_ptb = 0.0d0; pcn_pdt = 0.0d0
      pcn_ptl2= 0.0d0; pcn_ptb2= 0.0d0; pcn_pdt2= 0.0d0
      pcn_psi = 0.0d0; pcn_ptf = 0.0d0
      do i = 1, ndmfl
      pcn_pfl(i) = 0.0d0
      enddo
!
!::loop  ia
      do ia = 1, nion
      m1a = 2*ia - 1
!
!::integrate in plasma
      zsum1 = 0.0d0; zsum2 = 0.0d0; zsum3 = 0.0d0; zsum4 = 0.0d0
      zsum1p= 0.0d0; zsum2p= 0.0d0; zsum3p= 0.0d0
      do it = 2, itmax-1
      if( it.eq.itpve ) cycle
      jmax = jtmax(it)
      do jw = 2, jmax-1
      j = jcel(jw,it)
      i = icel(jw,it)
!xx   zvlm  = hvol(j,i)/(dtrat(j,i)*dteq(m1a))  ! KS 2011/09/26
      zvlm  = hvol(j,i)/dtrateq(j,i,m1a)
!-----
      zsum1 = zsum1 +  q1a(j,i,ia)*zvlm
      zsum2 = zsum2 + wq1a(j,i,ia)*zvlm
      zsum3 = zsum3 + (q1a(j,i,ia)-wq1a(j,i,ia))*zvlm
      zsum4 = zsum4 + ama(ia)*(ssnc(j,i,ia)+ssnv(j,i,ia)*vna(j,i,ia))
     >                           *hvol(j,i)
!-----
      zsum1p = zsum1p +  q1a(j,i,ia)*hvol(j,i)
      zsum2p = zsum2p + wq1a(j,i,ia)*hvol(j,i)
      zsum3p = zsum3p + qvl_dt(j,i,m1a)
!-----
!
      enddo  ! loop(jw)
      enddo  ! loop(it)
!-----
      pcn_ptl = pcn_ptl + zsum1/ama(ia)
      pcn_ptb = pcn_ptb + zsum2/ama(ia)
      pcn_pdt = pcn_pdt + zsum3/ama(ia)
      pcn_psi = pcn_psi + zsum4/ama(ia)
!-----
      pcn_ptl2 = pcn_ptl2 + zsum1p/ama(ia)
      pcn_ptb2 = pcn_ptb2 + zsum2p/ama(ia)
      pcn_pdt2 = pcn_pdt2 + zsum3p/ama(ia)
!
!::flux (sol)
      zsum = 0.0d0
      it = 1
      jmax = jtmax(it)
      do jw = 2, jmax-1
      j = jcel(jw,it)
      i = icel(jw,it)
      zsum = zsum + qfy_df(j,i,m1a)
      enddo
      sg = -1.0d0
      pcn_pfl(1) = pcn_pfl(1) + sg*zsum/ama(ia)
!
!::flux (idp)
      zsum = 0.0d0
      do it = 2, itpve-1
      jw = jtmax(it)-1
      j  = jcel(jw,it)
      i  = icel(jw,it)
      zsum = zsum + qfx_cv(j,i,m1a)+qfx_df(j,i,m1a)
      enddo
      sg = 1.0d0
      pcn_pfl(2) = pcn_pfl(2) + sg*zsum/ama(ia)
!
!::flux (prv)
      zsum = 0.0d0
      it = itpve-1
      jmax = jtmax(it)
      do jw = 2, jmax-1
      j = jcel(jw,it)
      i = icel(jw,it)
      zsum = zsum + qfy_df(j,i,m1a)
      enddo
      sg = 1.0d0
      pcn_pfl(3) = pcn_pfl(3) + sg*zsum/ama(ia)
!
!::flux (odp)
      zsum = 0.0d0
      do it = 2, itpve-1
      jw = 1
      j  = jcel(jw,it)
      i  = icel(jw,it)
      zsum = zsum + qfx_cv(j,i,m1a)+qfx_df(j,i,m1a)
      enddo
      sg = -1.0d0
      pcn_pfl(4) = pcn_pfl(4) + sg*zsum/ama(ia)
!
!::flux (man)
      zsum = 0.0d0
      it = itmpe-1
      jmax = jtmax(it)
      do jw = 2, jmax-1
      j = jcel(jw,it)
      i = icel(jw,it)
      zsum = zsum + qfy_df(j,i,m1a)
      enddo
      sg = 1.0d0
      pcn_pfl(5) = pcn_pfl(5) + sg*zsum/ama(ia)
!
!::flux (puf)
      pcn_pfl(6) = pcn_pfl(6) + 0.0d0
!
!::flux (vol)
      pcn_pfl(7) = pcn_pfl(7) + 0.0d0
!
!::flux(vli)
      pcn_pfl(8) = pcn_pfl(8) + tfvl(ia)
!
      enddo  !  loop(ia)
!
!::total flux
      pcn_pdt = pcn_pdt/dtim
      pcn_ptf = pcn_pfl(1)+pcn_pfl(2)+pcn_pfl(3)+pcn_pfl(4)+pcn_pfl(5)
      pcn_rhs = -pcn_pdt - pcn_ptf + pcn_psi
      pcn_tim = time
      pcn_dtm = dtim
      pcn_itm = itim
      pcn_nlp = nlp
!
!::debug write
!xx   call plpcon_out
!
      return
      end
!
!***********************************************************************
      subroutine plpcon_out
!***********************************************************************
      use cntmnt, only : i6_src, mcflx, mfmax, pfl_abs, pfl_err, pfl_ion
     >    , pfl_man, pfl_ntl, pfl_pmp, pfl_src, psm_abs, psm_err
     >    , psm_ion, psm_man, psm_ntl, psm_pmp, psm_src
      use cplqcn, only : pcn_dtm, pcn_itm, pcn_nlp, pcn_pdt, pcn_pdt2
     >    , pcn_pfl, pcn_psi, pcn_ptb, pcn_ptb2, pcn_ptf, pcn_ptl
     >    , pcn_ptl2, pcn_rhs, pcn_tim
      use cunit,  only : lmspe, lmype
      implicit none
!
!::argument
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   integer  nft
!
!::local variables
      real*8   zsm_ion
      integer  i, ift
      integer  lp; data lp/0/; save lp
!
!::debug write
      if( lmype.ne.lmspe ) return
!
!::summary (zflx)
      lp = lp + 1
      ift = 62
      if( mod(lp-1,50).eq.0 ) then
      write(ift,'(2x,"itim",2x,"time",6x,"dtim",6x,
     > "Tot_N",5x,"dlt_N",5x,"Resd",6x,"dtT_N",6x,"-Fl+S",6x,
     > "Finp",7x,"Fout",7x,"Fcore",6x,"Fpmp",7x,"Fl_mp",6x,
     > "Fl_sw",5x,"Fl_id",5x,"Fl_pw",5x,"Fl_od",5x,"Srec")')
      endif
!
      write(ift,'(i6,1p5e10.2,1p7e11.3,1p5e10.3)')
     > pcn_itm, pcn_tim, pcn_dtm,
     > pcn_ptl, pcn_ptl-pcn_ptb, pcn_rhs, pcn_pdt, -pcn_ptf+pcn_psi,
     > -pcn_pfl(5)+pfl_ntl(6), psm_man+psm_abs, psm_man, psm_pmp,
     >                                                     pcn_pfl(5),
     > pcn_pfl(1), pcn_pfl(2), pcn_pfl(3), pcn_pfl(4), pfl_ion(8)
!
      if( i6_src.le.0 ) return
      ift = i6_src
!
!::plpcon
      write(ift,'(/2x,"***  plpcon ***  ",1pe14.6,1pe11.2,i7,i3,2x,
     >  "rhs,Dt,-Fl,Si =",1p4e14.6)')
     >  pcn_tim, pcn_dtm, pcn_itm, pcn_nlp
     > ,pcn_rhs, pcn_pdt, -pcn_ptf, pcn_psi
      write(ift,'(2x,"pcn_ptl =",1p2e14.6,"  pcn_ptb =",1p2e14.6,
     > "  dlt_ptl =",1p2e14.6,"  pcn_pdt =",1p2e14.6)')
     >  pcn_ptl, pcn_ptl2, pcn_ptb, pcn_ptb2,
     >  pcn_ptl-pcn_ptb, pcn_ptl2-pcn_ptb2, pcn_pdt, pcn_pdt2
!
      write(ift,'(2x,"conserv =",1p12e14.6)')
     > pcn_ptl, pcn_ptl-pcn_ptb, pcn_rhs,
     > pcn_pdt, -pcn_ptf+pcn_psi, -pcn_ptf, pcn_psi
      write(ift,'(2x,"pcn_pfl =",1p5e14.6)')
     > pcn_pfl(1), pcn_pfl(2), pcn_pfl(3), pcn_pfl(4), pcn_pfl(5)
!
!::plsorc
      write(ift,'(2x,"i",1x,"flx",6x,"ion-plpflx",4x,"ion-ntnflx",4x,
     >  "ntl",11x,"src",11x,"man",11x,"abs",11x,"pmp",11x,"err")')
!
      zsm_ion = 0.0d0
      do i = 1, mfmax
      write(ift,'(i3,1x,a,1p10e14.6)')
     >  i,mcflx(i),pcn_pfl(i),pfl_ion(i),
     >  pfl_ntl(i),pfl_src(i),pfl_man(i),
     >  pfl_abs(i),pfl_pmp(i),pfl_err(i)
      zsm_ion = zsm_ion + pcn_pfl(i)
      enddo
      write(ift,'(i3,1x,a,1p10e14.6)') 0,"total ",
     > zsm_ion,psm_ion,psm_ntl,psm_src,psm_man,psm_abs,psm_pmp,psm_err
!
      return
      end
