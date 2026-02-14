!***********************************************************************
      subroutine pwstep(kout,dtstp,kcnv)
!***********************************************************************
!
!       function   solve q(N+1,l+1)    dtstp = t^N+1 - t^N
!
!        kout  : I  print option of information on covergence
!                    = 0 (no print),  = 1 (final result)
!                    = 2 (iteration)
!        dtstp : I  time step
!        kcnv  : O  flag of convergence
!                    = 0 (coverge),   = 1 (iterlation)
!                    = 2 (dt-step),   = 3 (negative value)
!
!        omega   relaxsation parameter (pwstep)  ! ## 2002/11/22
!
!         check linealized source   plstsrc     2009/09/16
!
!------------------------------------------------------------------------
      use cntmnt,     only : i6_src, n6_src
      use cplcom,     only : ddtq, ddtq1, ddtq2, ddtq3, ddtq4, dlpq
     >    , dlpq1, dlpq2, dlpq3, dlpq4, dq1a, dq2a, dq3, dq4, edtmx
     >    , edtq, edtq1, edtq2, edtq3, edtq4, elpmx, elpq, elpq1, elpq2
     >    , elpq3, elpq4, ibcyl, nion, nlp, nlpmn, nlpmx, q1a, q2a, q3
     >    , q4, rxflp
      use cplmet,     only : icel, itmax, jcel, jtmax, jtmin
      use csonic,     only : dtim, itim, kpcn, lfopt, lstop, time
      use cunit,      only : cdrout, n6
      use topics_mod, only : dtcal, dtduc, dtim_s, lgtpc
      implicit none
!
!::argument
      integer kout, kcnv;  real*8 dtstp

!
!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer ntr, it, ipbc, jwst, jwen, ia, jw, j, i, l, ic
      integer ntr, it, ipbc, jwst, jwen, ia, jw, j, i
      real*8  omega
      character cnum*10, cdsn*120
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer mji, indxr, lenx, iw, m
! modified 1/1 lines no used by kamata 2022/05/29
!ik   integer mji, indxr, lenx, iw
      integer mji, indxr, lenx
!
!------------------------------------------------------------------------
!::dtim
!------------------------------------------------------------------------
      rewind 81
!
! added 8 lines integrated calculation with TOPICS by kamata 2020/11/16
      dtim_s = dtstp
      if( lgtpc == 1 .and. dtim_s > dtduc ) then
        if( dtduc > 0.0_8 ) then
          dtim   = dtduc
        else
          dtim   = min( dtim_s, dtcal )
        endif
      else
        dtim  = dtim_s
! added 1 line integrated calculation with TOPICS by kamata 2020/11/16
      endif

      n6_src = 72
      nlp = 0
!
!xx   call lst_rdst("start")    ! <== 2009/07/28
!
      if( kout.eq.2 )
     > write(n6,'(2x,"*** pwstep ***   time =",1pe15.6,"  dtim ="
     >  ,1pe12.3,"  itim =",i5,"  nlpmx =",i2)') time,dtim,itim,nlpmx
!
!-----------------------------------------------------------------------
!::iterlation
!-----------------------------------------------------------------------
      do 1000 ntr = 1, nlpmx
      nlp = ntr
      omega = rxflp(nlp)
!
!::smoothing
!xx      call pdstep
!xx      call plimit
!
!::parallel diffusion for Te < temn_dprv
! deleted 4 lines no used by kamata 2022/05/29
!ik   do iw = 1, ndy
!ik   ibcdpo(iw) = 0
!ik   ibcdpi(iw) = 0
!ik   enddo
!
!-----------------------------------------------------------------------
!::particle conservation
!-----------------------------------------------------------------------
      if( i6_src.gt.100 ) then
      else
      i6_src = 0
      if( kpcn.eq.2 .and. nlp.eq.nlpmn+1 ) then
      i6_src = n6_src
      write(cnum,'(i10)') itim
      mji = indxr(cnum," ")
      cdsn = cdrout(1:lenx(cdrout))//"zpcn"//"_"//cnum(mji+1:10)
      open( unit=i6_src, file=cdsn )
      endif
      endif
!
!-----------------------------------------------------------------------
!::weit for X
!-----------------------------------------------------------------------
      call plwtfy
      call plfimp
!
!-----------------------------------------------------------------------
!::linearlized source terms  (sne => snic,snic)
!-----------------------------------------------------------------------
      call plpflx
      if( lstop.ne.0 ) goto 900
      call plsorc(1)         ! Sn, Wi, We at soldor

      if( kpcn > 0 .and. nlp > nlpmn ) then
        call pldens("Pw")    ! N0, E0 Flx
      endif

!xx   call plsorc(2)  !  abs. value of N0, E0 =>  Sn,Wi=0.0
!
!-----------------------------------------------------------------------
!::clear qcon ! %%% 2002/10/25
!-----------------------------------------------------------------------
      call plqcon_cl
!
!-----------------------------------------------------------------------
!::poloidal direction
!-----------------------------------------------------------------------
      call pxstep
      if( lstop.eq.1 ) goto 900
      call plpcon
!
!-----------------------------------------------------------------------
!::psi derection
!-----------------------------------------------------------------------
      call pystep
!
!-----------------------------------------------------------------------
!::con. variables q(N+1,l+1) including boundary
!-----------------------------------------------------------------------
      do 210 it = 1, itmax
      ipbc = ibcyl(it)
      jwst = jtmin(it)
      jwen = jtmax(it)
      if( ipbc.eq.1 ) then
      jwst = jwst + 1
      jwen = jwen - 1
      endif
      do 220 ia = 1, nion
      do 230 jw = jwst, jwen
      j = jcel(jw,it)
      i = icel(jw,it)
      q1a(j,i,ia) = q1a(j,i,ia) + dq1a(j,i,ia)*omega
      q2a(j,i,ia) = q2a(j,i,ia) + dq2a(j,i,ia)*omega
 230  continue
 220  continue
      do 240 jw = jwst, jwen
      j = jcel(jw,it)
      i = icel(jw,it)
      q3(j,i) = q3(j,i) + dq3(j,i)*omega
      q4(j,i) = q4(j,i) + dq4(j,i)*omega
 240  continue
 210  continue
!
      if( i6_src .gt. 0 ) close(i6_src)
!
!-----------------------------------------------------------------------
!::check minimum value
!-----------------------------------------------------------------------
!xx   call plkeep
!xx   call plmdfy
!
!-----------------------------------------------------------------------
!::aux. parameter except for cornor points
!-----------------------------------------------------------------------
      call plauxv
      if( lstop.eq.1 ) goto 900
!
!-----------------------------------------------------------------------
!::con. & aux variables at cornor points
!-----------------------------------------------------------------------
      if( lfopt(2).eq.1 ) then
        if( lfopt(3).eq.0 ) then
          call pycpnt0
        else
          call pycpnt
        endif
      if( lstop.eq.1 ) goto 900
      endif
!
      call plauxv
!-----
!xx      call pxbclst2
!xx      call lst_rdst("cal")    ! <== 2009/07/28

!xx      call check_pxbcsd(1)
!xx      call check_pxbcsd(2)
!-----
      if( lstop.eq.1 ) goto 900
!
!-----------------------------------------------------------------------
!::variation in iterlation
!-----------------------------------------------------------------------
      call plvrlp
!xx   call plwrlp(64)
!
!xx   dlpq  = dmax1( dlpq1,dlpq2,dlpq3,dlpq4 )
!xx   dlpq  = dmax1( dlpq1,0.5d0*dlpq2,dlpq3,dlpq4 )  !  99/11/07
      dlpq  = dmax1( dlpq1,   dlpq3,dlpq4 )  !  02/10/09
      elpq  = dmax1( elpq1,   elpq3,elpq4 )
!
      if( kout.eq.2 ) then
      write(n6,'(2x,"lp =",i3,f8.3,2x,1pe10.2,2x,1p4e10.2,2x,1p4e10.2
     >  )') nlp,omega,dlpq, dlpq1,dlpq2,dlpq3,dlpq4,
     >  elpq1,elpq2,elpq3,elpq4
      endif
!
!::minimum loop
      if( nlp.le.nlpmn ) goto 1000
!
!::convergence check
      if( dlpq.le.elpmx ) goto 1200  ! <== steady state
!xx   if( elpq.le.elpmx ) goto 1200  ! <== time evolution
!
!-----------------------------------------------------------------------
!::iterration loop
!-----------------------------------------------------------------------
 1000 continue
!
!-----------------------------------------------------------------------
!::not converge
!-----------------------------------------------------------------------
      kcnv = 1
!x      if( kout.ge.0 ) then
!x      write(n6,'(2x,"large dlpq at time",2x,1p2e12.3,2x,"lp =",i2,
!x     >  1p2e10.2)') time,dtim,nlp,dlpq,elpmx
!x      endif
      return
!
!-----------------------------------------------------------------------
!::converge
!-----------------------------------------------------------------------
 1200 continue
!x      if( kout.eq.1 ) then
!x      write(n6,'(2x,"t,dt =",1p2e12.3,"  lp",i3,1p4e10.2)')
!x     >  time,dtim,nlp,dlpq1,dlpq2,dlpq3,dlpq4
!x      endif
!
!::qcon %%% 2002/10/25
!x      if( kpcn.eq.1 ) then
!x      call plqcon_sx(jcxp1+1,jcxp2-1,icspx,icmpe-1,1)
!x      call plqcon_sx( 2,     jcxp1,  icspx+1,icwl2-1,1)
!x      call plqcon_sx(jcxp2,  jcmax-1,icspx+1,icwl2-1,1)
!x      endif
!
!
!-----------------------------------------------------------------------
!::variation during dtim
!-----------------------------------------------------------------------
      call plvrdt
!xx   edtq = dmax1( edtq1,edtq2,edtq3,edtq4 )
      edtq = dmax1( edtq1,edtq3,edtq4 )
      ddtq = dmax1( ddtq1,ddtq2,ddtq3,ddtq4 )
      if( kout.eq.2 ) then
      write(n6,'(2x,"dt =",i3,10x,1pe10.2,2x,1p4e10.2,2x,1p4e10.2)')
     >  0,ddtq, ddtq1,ddtq2,ddtq3,ddtq4, edtq1,edtq2,edtq3,edtq4
      endif
!
!::next step
      kcnv = 0
!xx   if( ddtq.le.edtmx ) return
      if( edtq.le.edtmx ) return
!
!::too large dtim
      kcnv = 2
      if( kout.ge.0 ) then
      write(n6,'(2x,"large edtq at time",2x,1pe15.6,1pe12.3,2x,"lp =",
     >  i2,1p2e10.2)') time,dtim,nlp,edtq,edtmx
      endif
      return
!
!-----------------------------------------------------------------------
!::serious error  (negative value)
!-----------------------------------------------------------------------
 900  continue
      kcnv = 3
      return
      end
