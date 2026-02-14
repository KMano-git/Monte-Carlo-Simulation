!***********************************************************************
      subroutine pldfan
!***********************************************************************
      use cplcom, only : factor_bal, fcda, fcet, fdsl2, itsl2, jen_bal
     >    , jst_bal, mdl_bal, nion, vdda, vdet, vdxe, vdxi, vlda, vlet
     >    , vlxe, vlxi, fdeg2
     >    , lerp_points_vldar, lerp_points_vletr
     >    , lerp_points_vlxir, lerp_points_vlxer
     >    , lerp_r_vldar, lerp_diff_vldar 
     >    , lerp_r_vletr, lerp_diff_vletr
     >    , lerp_r_vlxir, lerp_diff_vlxir
     >    , lerp_r_vlxer, lerp_diff_vlxer
      use cplmet, only : icel, icmpe, icmps, itmax, itmpe, itmps, itsle
     >    , jcel, jcxp1, jcxp2, jtmax, jtmin, kreg
      use cplvpn, only : vdvp
      use cpmpls, only : jmd1, jmd2
      use csize,  only : ndsp, ndx, ndy
      use cunit,  only : n6
      implicit none
!
!::local variables
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  it, jt, jts, jte, jc, ic, irg, ia, iy
!ik   integer  jchsl, icosl, icisl, i6
      integer  it, jt, jts, jte, jc, ic, irg, ia
      integer  jchsl, icosl, icisl
      real*8   fcenh, fcslw

      logical :: secondSOL_Flag = .false. ! True:Enhancement diffusion coefficient in 2nd-SOL. False:No enhanced
      logical :: balloon_Flag = .false.   ! True:With ballooning. False:No ballooning
      real*8 factor_bal_tmp ! equal to factor_bal if balloon_Flag = True
      real*8 :: edge_factor = 1.0 ! enhancement factor of diffusion coefficient in edge

      integer n_omd ! size of vldar
      real*8, allocatable :: vldar(:),vletr(:),vlxir(:),vlxer(:) ! diffsion cefficient distribution for ic direction

!
!::allocate and set value of vldar
      call set_vldar(vldar,vletr,vlxir,vlxer,n_omd)
!
!::clear
      vdda(1:ndx,1:ndy,1:ndsp) = 0.0d0
      vdet(1:ndx,1:ndy,1:ndsp) = 0.0d0
      vdxe(1:ndx,1:ndy) = 0.0d0
      vdxi(1:ndx,1:ndy) = 0.0d0
      vdvp(1:ndx,1:ndy) = 0.0d0
!
!==================== 
! 2018/01/16 Anomalous heat diffusivity poloidal setting - Homma
!::Da, Xi, Xe in main plasma
      write(n6,'(/2x,"*** pldfan ***  Ballooning model Anomalous heat
     > diffusivity poloidal setting ON - Homma")')
!
!     Ballooning flag check
      if(mdl_bal .eq. 0) then
        write(n6,'(/2x,"Balloon chi poloidal distribution OFF")')
      elseif((jst_bal.lt.0).or.(jen_bal.lt.0)) then
        write(n6,'(/2x,"jst_bal and jen_bal must be positive.    &
     & Ballooning model SWITCHED OFF!")')
      else
        write(n6,'(/2x,"Balloon chi poloidal distribution ON")')
        balloon_Flag = .true.
      endif
!==================== 
!::2nd-SOL
      if( itsl2.le.0 .or. fdsl2.le.0.0d0 ) then
        write(n6,'(/2x,"*** pldfan ***  
     >   no enhancement diff in 2nd-SOL")')
      else
        secondSOL_Flag = .true.
        jchsl = (jcxp1+jcxp2)/2
        icosl = itsl2
        icisl = itsl2
        fcenh = fdsl2
        write(n6,'(/2x,"*** pldfan (2nd-SOL) ***")')
        write(n6,'(2x,"jchsl =",i3,"  icosl =",i3,"  icisl =",i3,
     >    "  fcenh =",f8.3)') jchsl, icosl, icisl, fcenh
      endif
!====================
!::diffusion coef. in cell
      do it = 1, itmax ! radial direction
        jts = jtmin(it)
        jte = jtmax(it)
        if( it.ge.itmps .and. it.le.itmpe ) then ! edge location
          jts = jts + 1
          jte = jte - 1
        endif
!
        do jt = jts, jte ! poloidal direction 
          jc  = jcel(jt,it)
          ic  = icel(jt,it)
          irg = kreg(jc,ic)

          !::enhancement in 2nd SOL region
          fcslw = 1.0d0
          if( secondSOL_Flag ) then
            if( it.le.itsle ) then
              if( jc.le.jchsl ) then
                if( ic.le.icosl ) fcslw = fcenh
              else
                if( ic.le.icisl ) fcslw = fcenh
              endif
            endif
          endif

          !::enhancement in edge region
          edge_factor = 1.0
          if( it .ge. itmps .and. it .le. itmpe) then
            edge_factor = fdeg2
          endif

          do ia = 1, nion
            ! vldar is not used in private regine
            if((lerp_points_vldar>0).and.(irg.ne.4).and.(irg.ne.5)) then ! value set by vldar
              vdda(jc,ic,ia) = vldar(ic)*fcda(ia)
            else
              vdda(jc,ic,ia) = vlda(irg)*fcda(ia)*fcslw*edge_factor
            endif

            ! vletr is not used in private regine
            if((lerp_points_vletr>0).and.(irg.ne.4).and.(irg.ne.5)) then ! value set by vletr
              vdet(jc,ic,ia) = vletr(ic)*fcet(ia)
            else
              vdet(jc,ic,ia) = vlet(irg)*fcet(ia)*fcslw*edge_factor
            endif
          enddo

          ! balloing
          factor_bal_tmp = 1.0
          if((irg == 6) .and. balloon_Flag) then
            if((jc.ge.jst_bal).and.(jc.le.jen_bal) ) then
               factor_bal_tmp = factor_bal
            endif
          endif

          ! vlxir is not used in private regine
          if((lerp_points_vlxir>0).and.(irg.ne.4).and.(irg.ne.5)) then ! value set by vlxi
            vdxi(jc,ic) = factor_bal_tmp*vlxir(ic)
          else
            vdxi(jc,ic) = factor_bal_tmp*vlxi(irg)*fcslw*edge_factor
          endif

          ! vlxer is not used in private regine
          if((lerp_points_vlxer>0).and.(irg.ne.4).and.(irg.ne.5)) then ! value set by vlxer
            vdxe(jc,ic) = factor_bal_tmp*vlxer(ic)
          else
            vdxe(jc,ic) = factor_bal_tmp*vlxe(irg)*fcslw*edge_factor
          endif
        enddo ! jt 
      enddo !it
!
!====================
!::debug write
      write(n6,*) "homma chi poloidal setting, check output  ++"
      write(n6,*) "itmps, itmpe", itmps, itmpe
      write(n6,*) "icmps, icmpe", icmps, icmpe
      write(n6,*) "jmd1, jmd2", jmd1, jmd2
      write(n6,*) "mdl_bal, jst_bal,jen_bal,factor_bal"
      write(n6,'(3i6,1(f8.3))') mdl_bal, jst_bal,jen_bal,factor_bal

      write(n6,*)
      write(n6,*) " vdxi := Xi, vdxe := Xe"

      write(n6,'(/2x,"jc, ic, jtmin(it), jtmax(it), jcel(jt,it),
     >  icel(jt,it), kreg(jc,ic), vdxi(jc,ic), vdxe(jc,ic)")')
      write(n6,*)
      write(n6,*)

!====================
      it = itmpe -1
      jts = jtmin(it)
      jte = jtmax(it)
      if( it.ge.itmps .and. it.le.itmpe ) then
         jts = jts + 1
         jte = jte - 1
      endif
!
      do jt = jts, jte
         jc  = jcel(jt,it)
         ic  = icel(jt,it)
         irg = kreg(jc,ic)

         write(n6,'(7i6,2f8.3)') jc, ic, jtmin(it), jtmax(it),
     >        jcel(jt,it),icel(jt,it), kreg(jc,ic), vdxi(jc,ic),
     >        vdxe(jc,ic)
      enddo ! jt
                                                                                          
      write(n6,*)
      write(n6,*)  "Radial evolution at Outer Midplane"
      write(n6,'(/2x,"jc, ic, jtmin(it), jtmax(it), jcel(jt,it),                                                 
     >  icel(jt,it), kreg(jc,ic), vdxi(jc,ic), vdxe(jc,ic)")')
      do ic = icmps-2, icmpe
            write(n6,'(7i6,2f8.3)') jmd1, ic, jtmin(it), jtmax(it),
     >           jcel(jt,it),icel(jt,it), kreg(jmd1,ic), vdxi(jmd1,ic),
     >           vdxe(jmd1,ic)
      enddo
      write(n6,*)
      write(n6,*) "++++++++++++++++++++++++++++++++++++++"
!====================

      call dbg_pldfan
!
      return

      contains
!***********************************************************************
      subroutine set_vldar(vldar,vletr,vlxir,vlxer,n_omd)
!***********************************************************************
      use cplmet, only : icmax, icmin
      implicit none
! arguments
      real(8),allocatable,intent(out) :: 
     >  vldar(:),vletr(:),vlxir(:),vlxer(:)
      integer,intent(out) :: n_omd
! local
      real(8) :: rhf_omd(ndy)
      real(8) :: rst(ndy), ren(ndy)

! calculate rhf_omd for vldar
      call codrad(jmd1,rst,ren,rhf_omd,n_omd,ndy
     > , icmin(jmd1),icmax(jmd1))
      allocate(vldar(n_omd))
      allocate(vletr(n_omd))
      allocate(vlxir(n_omd))
      allocate(vlxer(n_omd))

! sort and interpolate value of vldar,vletr,vlxir,vlxer
      call interpolate_array(vldar,lerp_points_vldar
     > ,lerp_r_vldar(1:lerp_points_vldar)
     > ,lerp_diff_vldar(1:lerp_points_vldar)
     > ,rhf_omd(1:n_omd),n_omd)
      call interpolate_array(vletr,lerp_points_vletr
     > ,lerp_r_vletr(1:lerp_points_vletr)
     > ,lerp_diff_vletr(1:lerp_points_vletr)
     > ,rhf_omd(1:n_omd),n_omd)
      call interpolate_array(vlxir,lerp_points_vlxir
     > ,lerp_r_vlxir(1:lerp_points_vlxir)
     > ,lerp_diff_vlxir(1:lerp_points_vlxir)
     > ,rhf_omd(1:n_omd),n_omd)
      call interpolate_array(vlxer,lerp_points_vlxer
     > ,lerp_r_vlxer(1:lerp_points_vlxer)
     > ,lerp_diff_vlxer(1:lerp_points_vlxer)
     > ,rhf_omd(1:n_omd),n_omd)
      end subroutine set_vldar
      end subroutine pldfan
!
!***********************************************************************
      subroutine dbg_pldfan
!***********************************************************************
      use cplcom, only : vdda, vdet, vdxe, vdxi
      use cplmet, only : icmpe, icwl1, icwl2, jcxp1, jcxp2
      use cpmpls, only : jmd1
      use cunit,  only : n6
      implicit none
!
!::local variables
      integer :: ia, j, i, ist, ien
!
      ia  = 1
      j   = jcxp1 - 5
      ist = icwl1
      ien = icwl2
      write(n6,'(2x,"diff coeff (D,Et,Xi,Xe)  j =",i3," i =",2i5)')
     >   j, ist, ien
      write(n6,'(5x,"Da =",10f8.3)') (vdda(j,i,ia),i=ist,ien)
      write(n6,'(5x,"Et =",10f8.3)') (vdet(j,i,ia),i=ist,ien)
      write(n6,'(5x,"Xi =",10f8.3)') (vdxi(j,i),i=ist,ien)
      write(n6,'(5x,"Xe =",10f8.3)') (vdxe(j,i),i=ist,ien)
!
      j   = jmd1
      ist = icwl1
      ien = icmpe
      write(n6,'(2x,"diff coeff (D,Et,Xi,Xe)  j =",i3," i =",2i5)')
     >   j, ist, ien
      write(n6,'(5x,"Da =",10f8.3)') (vdda(j,i,ia),i=ist,ien)
      write(n6,'(5x,"Et =",10f8.3)') (vdet(j,i,ia),i=ist,ien)
      write(n6,'(5x,"Xi =",10f8.3)') (vdxi(j,i),i=ist,ien)
      write(n6,'(5x,"Xe =",10f8.3)') (vdxe(j,i),i=ist,ien)
!
      j   = jcxp2 + 5
      ist = icwl1
      ien = icwl2
      write(n6,'(2x,"diff coeff (D,Et,Xi,Xe)  j =",i3," i =",2i5)')
     >   j, ist, ien
      write(n6,'(5x,"Da =",10f8.3)') (vdda(j,i,ia),i=ist,ien)
      write(n6,'(5x,"Et =",10f8.3)') (vdet(j,i,ia),i=ist,ien)
      write(n6,'(5x,"Xi =",10f8.3)') (vdxi(j,i),i=ist,ien)
      write(n6,'(5x,"Xe =",10f8.3)') (vdxe(j,i),i=ist,ien)
!
      return
      end
