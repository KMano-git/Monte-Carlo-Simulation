!***********************************************************************
      subroutine pst_qdpl( cpath, ctime )
!***********************************************************************
!
!   Yacobian at boundary
!
!   boundary condition at divertor plate
!      j=1  (outer plate) / j=jmax (inner plate)
!
!  (1) na   S^sita*b^sita*Q2a = flb1a
!  (2) va   va = Csd   sound velocity
!  (3) Ti   Vj*dQ3/dt + qfx_cv + qfx_cd + qfx_vh =
!             S^sita*b^sita*{1/2*flb1ro*Vs**2+cdli*flb1ni*Ti*eV}
!  (4) Te   Vj*dQ4/dt + qfx_cv + qfx_cd =
!             S^sita*b^sita*{cdle*flb1ne*Te*eV}
!
!                        Note.  no include qfx_df term (ano-diff)
!   flux at dplate
!
!   flb1a = S^sita*b^sita * ma*na*va
!   flb2a = S*b * (ma*na*va^2+Pa + ]~[a )
!   flb3  = S*b *{ Sum_a[(1/2*ma*na*va^2+3/2*na*Ti+Pa)*va] + q//i
!                               + Sum_a[(va*]~[)]a }
!   flb4  = S*b *{ (3/2*ne*Te+Pe)*va + q//e }
!
!   flb1ro = sum_a[flb1a]
!   flb1ni = sum_a[1/ma*flb1a]
!   flb1ne = sum_a[za/ma*flb1a]
!
                                      ! %%% 2002/11/11
!
!      V*dQ/dt = -[cond+conv]1/2 - Hbc     j = 1 (o-dplate)
!      V*dQ/dt = +[cond+conv]J-1/2 - Hbc   j = J (i-dplate)
!        Hbc = sig(V)*flx*T   flx = n*V
!                                     ! 2006/08/07
!
!    flsmi(it,1) =   e-conv//
!    flsmi(it,2) = + e-cond//
!    flsmi(it,3) = + i-conv
!    flsmi(it,4) = + i-cond + i-viss
!    flsmi(it,5) = + recomb at plate
!    flsmi(it,6) = + radiation
!    flsmi(it,7) = + neutral
!
!     q_target = flsmi * h*cosw
!     sum_q_target
!
!-----------------------------------------------------------------------
      use cgdcom, only : grdx, grdy
      use cphcns, only : cev, cpi
      use cplcom, only : ama, aza, c12, cdle, cdli, fcbcvl, nion, nlp
     >    , q1a, vcs, vna, vne, vte, vti, vva, vve
      use cplmet, only : hpit, hvol, hvsb, icel, itpve, itsle, itsls
     >    , jcel, jtmax
      use cplpst, only : dwrad, farei, fareo, fbdpi, fbdpo, fldpi, fldpo
     >    , flsmi, flsmo, rhf_idp, rhf_imd, rhf_odp, rhf_omd, tflsmi
     >    , tflsmo
      use cplqcn, only : qfx_cd, qfx_cv, qfx_vh 
      use cntwfl, only : dwntl
      use csize,  only : ndy
      use csonic, only : dtim, itim, time
      use cunit,  only : n6
!     22/08/24 yamamoto, for soldor time Series
      use mod_soldorTimeSeries
      implicit none
! arguments
      character, intent(in) :: cpath*(*), ctime*(*)
! cpath : path name
! ctime : time

!
!::local variables
      integer  nt1, nt2, m, it, nt, i, jmax, jw, j, jp, jm, jsf
      integer  ia, m1a, m2a, m3, m4
      integer  i6, i7, ic, k, j1
      real*8   sigv, zcs, wdle, zsmro, zsmni, zsmne, zroi
      real*8   flb1a, flb2a, flb1ro, flb1ni, flb1ne, flb3, flb4
      real*8   z3cs, z3ti, z4te
!
      integer  mx1, mx2, mx3, mx4, my1, my2, my3, my4
      real*8   r1, r2, r3, r4, z1, z2, z3, z4, dl
      real*8   sumi, sume, sums, zsv, zsp, fcvol
      real*8   fsg, zflx, erec, hcsw
      character  dsn*200, dsn2*200

!     22/08/29 yamamoto, Time averaging before accumulating flsmo and flsmi values
      integer, parameter :: howmuchAcc = 5 !qt_ev,qt_ed,qt_iv,qt_id,qt_rec
      real*8,dimension(howmuchAcc):: flsmo_ave 
      real*8,dimension(howmuchAcc):: flsmi_ave 
!     22/08/29 yamamoto, Get the maximum and minimum values in a time series such as qt_ed
!     Once the accumulation is taken, find the maximum and minimum after the accumulation.
      real*8,dimension(timeNum,howmuchAcc):: flsmo_accumulate
      real*8,dimension(timeNum,howmuchAcc):: flsmi_accumulate
      integer i_time
!
      write(n6,'(/2x,"*** pst_qdpl ***")')
!
      m3   = 2*nion + 1
      m4   = 2*nion + 2
      nt1 = itsls + 1
      nt2 = itpve - 1
!
      fcvol = fcbcvl
      if( itim.le.50 ) fcvol = 0.02d0
!
!::clear
      do m = 1, m4
      do it = 1, ndy
      fbdpi(it,m) = 0.0d0
      fbdpo(it,m) = 0.0d0
      fldpi(it,m) = 0.0d0
      fldpo(it,m) = 0.0d0
      enddo
      enddo
      do i = 1, 10
      tflsmi(i) = 0.0d0
      tflsmo(i) = 0.0d0
      do it = 1, ndy
      flsmi(it,i) = 0.0d0
      flsmo(it,i) = 0.0d0
      enddo
      enddo
!
!::loop (it)
      do nt = nt1, nt2
      it   = nt
      jmax = jtmax(it)
!
!-----------------------------------------------------------------------
!::boundary condition j = 1
!-----------------------------------------------------------------------
      jw   = 1
      i    = icel(jw,it)
      j    = jcel(jw,it)
      jp   = jcel(jw+1,it)
      jsf  = j
      sigv = -1.0d0
      zcs  = sigv*vcs(j,i)   !  for all species
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   kbc  = 1
!
!::are
      call mcpnt(j,i,mx1,mx2,mx3,mx4,my1,my2,my3,my4)
      r2 = grdx(mx2,my2)
      r3 = grdx(mx3,my3)
      z2 = grdy(mx2,my2)
      z3 = grdy(mx3,my3)
      dl = dsqrt((r3-r2)**2+(z3-z2)**2)
      fareo(it) = 2.0d0*cpi*0.5d0*(r2+r3)*dl
!
!::cdle
      wdle = cdle
!
!::hvol of dummy cell
      hvol(j,i) = hvol(jp,i)*fcvol
!
!::flux and variables at boundary
      zsmro = 0.0d0
      zsmni = 0.0d0
      zsmne = 0.0d0
      zroi  = 0.0d0
      do ia = 1, nion
      m1a = 2*ia - 1
      m2a = 2*ia
      flb1a = qfx_cv(jsf,i,m1a)
      flb2a = qfx_cv(jsf,i,m2a) + qfx_cd(jsf,i,m2a)
      zsmro = zsmro + flb1a
      zsmni = zsmni + flb1a/ama(ia)
      zsmne = zsmne + flb1a/ama(ia)*aza(ia)
      zroi  = zroi + q1a(j,i,ia)
      fldpo(it,m1a) = -flb1a
      fldpo(it,m2a) = -flb2a
      enddo
      flb1ro = zsmro
      flb1ni = zsmni
      flb1ne = zsmne
      flb3 = qfx_cv(jsf,i,m3)+qfx_cd(jsf,i,m3)+qfx_vh(jsf,i)
      flb4 = qfx_cv(jsf,i,m4)+qfx_cd(jsf,i,m4)
      fldpo(it,m3) = -flb3
      fldpo(it,m4) = -flb4
!
!::New type
      fsg = -1.0d0
      flsmo(it,1) =               fsg*qfx_cv(jsf,i,m4)
      flsmo(it,2) = flsmo(it,1) + fsg*qfx_cd(jsf,i,m4)
      flsmo(it,3) = flsmo(it,2) + fsg*qfx_cv(jsf,i,m3)
      flsmo(it,4) = flsmo(it,3) + fsg*(qfx_cd(jsf,i,m3)+qfx_vh(jsf,i))
      ia  = 1
      m1a = 2*ia - 1
      zflx = qfx_cv(jsf,i,m1a)/ama(ia)
      erec = 13.595
      flsmo(it,5) = flsmo(it,4)+fsg*zflx*erec*cev
      flsmo(it,6) = flsmo(it,5)+dwrad(i,4)
      flsmo(it,7) = flsmo(it,6)+dwntl(i,4)
!
      tflsmo(1) = tflsmo(1) + flsmo(it,1)
      tflsmo(2) = tflsmo(2) + flsmo(it,2)
      tflsmo(3) = tflsmo(3) + flsmo(it,3)
      tflsmo(4) = tflsmo(4) + flsmo(it,4)
      tflsmo(5) = tflsmo(5) + flsmo(it,5)
      tflsmo(6) = tflsmo(6) + flsmo(it,6)
      tflsmo(7) = tflsmo(7) + flsmo(it,7)
!
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   zdcs = c23/zroi
      z3cs = sigv*c12*flb1ro
      z3ti = sigv*cdli*flb1ni
      z4te = sigv*wdle*flb1ne
!
      do ia = 1, nion
      m1a = 2*ia - 1
      m2a = 2*ia
!
!::Ne
      fbdpo(it,m1a) = +q1a(j,i,ia)-qfx_cv(jsf,i,m1a)/(hvsb(jsf,i)*zcs)
!
!::Va
      fbdpo(it,m2a) = +vva(j,i,ia)-zcs
      enddo  ! loop(ia)
!
!::Hbc_I = sig(va)*(1/2*flb1ro*Cs^2+cdli*flb1ni*Ti*eV)
      fbdpo(it,m3) = +z3cs*zcs**2 + z3ti*vti(j,i)*cev
!
!::Hbc_e = sig(va)*cdle*flb1ne*Te*eV
      fbdpo(it,m4) = +z4te*vte(j,i)*cev
!
!-----------------------------------------------------------------------
!::boundary condition j = jmax
!-----------------------------------------------------------------------
      jw   = jmax
      i    = icel(jw,it)
      j    = jcel(jw,it)
      jm   = jcel(jw-1,it)
      jsf  = jcel(jw-1,it)
      sigv = +1.0d0
      zcs  = sigv*vcs(j,i)   !  for all species
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   kbc  = 2
!
!::are
      call mcpnt(j,i,mx1,mx2,mx3,mx4,my1,my2,my3,my4)
      r4 = grdx(mx4,my4)
      r1 = grdx(mx1,my1)
      z4 = grdy(mx4,my4)
      z1 = grdy(mx1,my1)
      dl = dsqrt((r4-r1)**2+(z4-z1)**2)
      farei(it) = 2.0d0*cpi*0.5d0*(r4+r1)*dl
!
!::cdle
      wdle = cdle
!
!::hvol of dummy cell
      hvol(j,i) = hvol(jm,i)*fcvol
!
!::flux and variables at boundary
      zsmro = 0.0d0
      zsmni = 0.0d0
      zsmne = 0.0d0
      zroi  = 0.0d0
      do ia = 1, nion
      m1a = 2*ia - 1
      m2a = 2*ia
      flb1a = qfx_cv(jsf,i,m1a)
      flb2a = qfx_cv(jsf,i,m2a) + qfx_cd(jsf,i,m2a)
      zsmro = zsmro + flb1a
      zsmni = zsmni + flb1a/ama(ia)
      zsmne = zsmne + flb1a/ama(ia)*aza(ia)
      zroi  = zroi + q1a(j,i,ia)
      fldpi(it,m1a) = +flb1a
      fldpi(it,m2a) = +flb2a
      enddo
      flb1ro = zsmro
      flb1ni = zsmni
      flb1ne = zsmne
      flb3 = qfx_cv(jsf,i,m3)+qfx_cd(jsf,i,m3)+qfx_vh(jsf,i)
      flb4 = qfx_cv(jsf,i,m4)+qfx_cd(jsf,i,m4)
      fldpi(it,m3) = +flb3
      fldpi(it,m4) = +flb4
!
!::New type
      fsg = +1.0d0
      flsmi(it,1) =               fsg*qfx_cv(jsf,i,m4)
      flsmi(it,2) = flsmi(it,1) + fsg*qfx_cd(jsf,i,m4)
      flsmi(it,3) = flsmi(it,2) + fsg*qfx_cv(jsf,i,m3)
      flsmi(it,4) = flsmi(it,3) + fsg*(qfx_cd(jsf,i,m3)+qfx_vh(jsf,i))
      ia  = 1
      m1a = 2*ia - 1
      zflx = qfx_cv(jsf,i,m1a)/ama(ia)
      erec = 13.595
      flsmi(it,5) = flsmi(it,4)+fsg*zflx*erec*cev
      flsmi(it,6) = flsmi(it,5)+dwrad(i,2)
      flsmi(it,7) = flsmi(it,6)+dwntl(i,2)
!
      tflsmi(1) = tflsmi(1) + flsmi(it,1)
      tflsmi(2) = tflsmi(2) + flsmi(it,2)
      tflsmi(3) = tflsmi(3) + flsmi(it,3)
      tflsmi(4) = tflsmi(4) + flsmi(it,4)
      tflsmi(5) = tflsmi(5) + flsmi(it,5)
      tflsmi(6) = tflsmi(6) + flsmi(it,6)
      tflsmi(7) = tflsmi(7) + flsmi(it,7)
!
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   zdcs = c23/zroi
      z3cs = sigv*c12*flb1ro
      z3ti = sigv*cdli*flb1ni
      z4te = sigv*wdle*flb1ne
!
      do ia = 1, nion
      m1a = 2*ia - 1
      m2a = 2*ia
!
!::Ne
      fbdpi(it,m1a) = +q1a(j,i,ia)-qfx_cv(jsf,i,m1a)/(hvsb(jsf,i)*zcs)
!
!::Va
      fbdpi(it,m2a) = +vva(j,i,ia)-zcs
      enddo  ! loop(ia)
!
!::Hbc_I = sig(va)*(1/2*flb1ro*Cs^2+cdli*flb1ni*Ti*eV)
      fbdpi(it,m3) = +z3cs*zcs**2 + z3ti*vti(j,i)*cev
!
!::Hbc_e = sig(va)*cdle*flb1ne*Te*eV
      fbdpi(it,m4) = +z4te*vte(j,i)*cev
!
      enddo  !  loop (it)
!
!
!::outer plate
      sumi = 0.0d0
      sume = 0.0d0
      sums = 0.0d0
      do it = 2, itpve-1
      sumi = sumi + fldpo(it,m3)
      sume = sume + fldpo(it,m4)
      sums = sums + fareo(it)
      enddo
!
      i6 = 21
      dsn = "zhflx_o.txt"
      if( ctime .ne. ' ' ) then
        dsn = trim( cpath ) // "zhflx_o" 
     >   // '_' // trim( ctime ) // ".txt"
      endif
      open(unit=i6,file=trim(dsn),position="append")
!
      write(i6,'(2x,"*** plhflx_odp ***  itim =",i7,"  time =",1pe12.3,
     > "  nlp =",i3)') itim, time, nlp
      it  = itsle
      j   = jcel(1,it)
      i   = icel(1,it)
      jsf = j
      zsv = fareo(it)
      zsp = hvsb(jsf,i)
      write(i6,'(2x,"HQi/HQe/HQtot =",1p3e11.3,"  Are =",1pe11.3,
     >  "  hcsw =",1pe11.3,"  pitch =",1pe11.3)')
     >  sumi, sume, sumi+sume, sums, zsp/zsv, hpit(jsf,i)
      write(i6,'(3x,"ic",2x,"r_odp",6x,"r_idp",6x,"r_omd",6x,"r_imd",
     > 6x,"Ned",8x,"Ted",8x,"Tid",8x,"Csd",
     > 8x,"Nec",8x,"Tec",8x,"Tic",8x,"Csc",
     > 8x,"gvi_b",6x,"qve_b",6x,"qv_b",7x,"qv_f",7x,
     > 11x,1x,"Sv",9x,"hcsw",7x,"pitch",
     > 6x,"qpi_s",6x,"qpi_b",6x,"qpi_f",
     > 6x,"qpe_s",6x,"qpe_b",6x,"qpe_f")')
!
      do it = 2, itpve-1
      j   = jcel(1,it)
      i   = icel(1,it)
      jsf = j
      ic  = i
      zsv = fareo(it)    ! S^sita
      zsp = hvsb(jsf,i)  ! S^sita*h*cosw
      ia  = 1
      write(i6,'(2x,i3,1p16e11.3,11x,1x,1p9e11.3)')
     >  ic, rhf_odp(ic)*100.0d0, rhf_idp(ic)*100.0d0,
     >      rhf_omd(ic)*100.0d0, rhf_imd(ic)*100.0d0,
     >  vne(j,i), vte(j,i), vti(j,i), vcs(j,i),
     >  vne(j+1,i), vte(j+1,i), vti(j+1,i), vcs(j+1,i),
     >  fbdpo(it,m3)/zsv, fbdpo(it,m4)/zsv,
     >  (fbdpo(it,m3)+fbdpo(it,m4))/zsv,
     >  (fldpo(it,m3)+fldpo(it,m4))/zsv,
     >  zsv, zsp/zsv, hpit(jsf,i),
     >  -vna(j,i,ia)*vva(j,i,ia)*
     >      (0.5d0*ama(ia)*vva(j,i,ia)**2+cdli*vti(j,i)*cev),
     >  fbdpo(it,m3)/zsp, fldpo(it,m3)/zsp,
     >  -cdle*vne(j,i)*vve(j,i)*vte(j,i)*cev,
     >  fbdpo(it,m4)/zsp, fldpo(it,m4)/zsp
      enddo
      close(i6)
!
!::inner plate
      sumi = 0.0d0
      sume = 0.0d0
      sums = 0.0d0
      do it = 2, itpve-1
      sumi = sumi + fldpi(it,m3)
      sume = sume + fldpi(it,m4)
      sums = sums + farei(it)
      enddo
!
      dsn = "zhflx_i.txt"
      if( ctime .ne. ' ' ) then
        dsn = trim( cpath ) // "zhflx_i" 
     >   // '_' // trim( ctime ) // ".txt"
      endif
      open(unit=i6,file=trim(dsn),position="append")
!
      write(i6,'(2x,"*** plhflx_idp ***  itim =",i7,"  time =",1pe12.3,
     > "  nlp =",i3)') itim, time, nlp
      it  = itsle
      jmax = jtmax(it)
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   j   = jcel(jmax,it)
      i   = icel(jmax,it)
      jsf = jcel(jmax-1,it)
      zsv = farei(it)
      zsp = hvsb(jsf,i)
      write(i6,'(2x,"HQi/HQe/HQtot =",1p3e11.3,"  Are =",1pe11.3,
     >  "  hcsw =",1pe11.3,"  pitch =",1pe11.3)')
     >  sumi, sume, sumi+sume, sums, zsp/zsv, hpit(jsf,i)
      write(i6,'(3x,"ic",2x,"r_odp",6x,"r_idp",6x,"r_omd",6x,"r_imd",
     > 6x,"Ned",8x,"Ted",8x,"Tid",8x,"Csd",
     > 8x,"Nec",8x,"Tec",8x,"Tic",8x,"Csc",
     > 8x,"gvi_b",6x,"qve_b",6x,"qv_b",7x,"qv_f",7x,
     > 11x,1x,"Sv",9x,"hcsw",7x,"pitch",
     > 6x,"qpi_s",6x,"qpi_b",6x,"qpi_f",
     > 6x,"qpe_s",6x,"qpe_b",6x,"qpe_f")')
!
      do it = 2, itpve-1
      jmax = jtmax(it)
      j  = jcel(jmax,it)
      i  = icel(jmax,it)
      jsf = jcel(jmax-1,it)
      ic = i
      zsv = farei(it)    ! S^sita
      zsp = hvsb(jsf,i)  ! S^sita*h*cosw
      ia  = 1
      write(i6,'(2x,i3,1p16e11.3,11x,1x,1p9e11.3)')
     >  ic, rhf_odp(ic)*100.0d0, rhf_idp(ic)*100.0d0,
     >      rhf_omd(ic)*100.0d0, rhf_imd(ic)*100.0d0,
     >  vne(j,i), vte(j,i), vti(j,i), vcs(j,i),
     >  vne(j-1,i), vte(j-1,i), vti(j-1,i), vcs(j-1,i),
     >  fbdpi(it,m3)/zsv, fbdpi(it,m4)/zsv,
     >  (fbdpi(it,m3)+fbdpi(it,m4))/zsv,
     >  (fldpi(it,m3)+fldpi(it,m4))/zsv,
     >  zsv, zsp/zsv, hpit(jsf,i),
     >  vna(j,i,ia)*vva(j,i,ia)*
     >      (0.5d0*ama(ia)*vva(j,i,ia)**2+cdli*vti(j,i)*cev),
     >  fbdpi(it,m3)/zsp, fldpi(it,m3)/zsp,
     >  cdle*vne(j,i)*vve(j,i)*vte(j,i)*cev,
     >  fbdpi(it,m4)/zsp, fldpi(it,m4)/zsp
      enddo
      close(i6)
!
!::New type
      dsn = "Qdpl_o.txt"
      if( ctime .ne. ' ' ) then
        dsn = trim( cpath ) // "Qdpl_o" //
     >    '_' // trim( ctime ) // ".txt"
      endif
      open(unit=i6,file=trim(dsn),position="append")
!
      i7 = 22
      dsn2 = "Qdpl_o_stat.txt"
      if( ctime .ne. ' ' ) then
        dsn2 = trim( cpath ) // "Qdpl_o_stat" // 
     >   '_' // trim( ctime ) // ".txt"
      endif
      open(unit=i7,file=trim(dsn2),position="append")
!
      write(i6,'(//2x,"*** plhflx_odp ***  itim =",i7,"  time =",
     > 1pe14.6,"  dtim =",1pe12.3,"  nlp =",i3)') itim,time,dtim,nlp
      write(i6,'(2x,"total q_heat =",2x,"e-conV",5x,"e-conD",5x,
     > "i-conV",5x,"i-conD",5x,"recomb",5x,"radiation",2x,"neutral")')
      write(i6,'(2x,"total q_heat =",1p7e11.3)') (tflsmo(k),k=1,7)

      write(i7,'(//2x,"*** plhflx_odp ***  itim =",i7,"  time =",
     > 1pe14.6,"  dtim =",1pe12.3,"  nlp =",i3)') itim,time,dtim,nlp
      write(i7,'(2x,"total q_heat =",2x,"e-conV",5x,"e-conD",5x,
     > "i-conV",5x,"i-conD",5x,"recomb",5x,"radiation",2x,"neutral")')
      write(i7,'(2x,"total q_heat =",1p7e11.3)') (tflsmo(k),k=1,7)
 
!    220826 output time_average, maximum and minimum of qt_ev","qt_ed","qt_iv","qt_id","qt_rec","Ned","Ted","Tid","Vid", and "Vni"
      write(i6,'(1x,"ic",2x,"r_odp",5x,"r_idp",
     > 5x,"r_omd",5x,"r_imd",5x,"pitch",6x,"csw",8x,"hcsw",
     > 7x,"qt_ev",6x,"qt_ed",6x,"qt_iv",
     > 6x,"qt_id",6x,"qt_rec",
     > 5x,"qt_rad",5x,"qt_ntl",
     > 5x,"qpi_sm",5x,"qpi_bc",5x,"qpe_sm",5x,"qpe_bc",
     > 5x,"Ned",8x,"Ted",8x,"Tid",8x,"Vid",8x,"Vni",
     > 8x,"Are")')

      write(i7,'(1x,"ic",2x,"r_odp",5x,"r_idp",
     > 5x,"r_omd",5x,"r_imd",
     > 5x,"qt_ev",6x,"qt_ed",6x,"qt_iv",
     > 6x,"qt_id",6x,"qt_rec",
     > 5x,"Ned",8x,"Ted",8x,"Tid",8x,"Vid",8x,"Vni",
     > 8x,"qt_ev_max",2x,"qt_ev_min",
     > 2x,"qt_ed_max",2x,"qt_ed_min",
     > 2x,"qt_iv_max",2x,"qt_iv_min",
     > 2x,"qt_id_max",2x,"qt_id_min",
     > 2x,"qt_rec_ma",2x,"qt_rec_mi",
     > 2x,"Ned_max",4x,"Ned_min",
     > 4x,"Ted_max",4x,"Ted_min",
     > 4x,"Tid_max",4x,"Tid_min",
     > 4x,"Vid_max",4x,"Vid_min",
     > 4x,"Vni_max",4x,"Vni_min",
     > 4x, "(max and min are the maximum 
     > and minimum values of the last few time steps.)")')

      do it = 2, itpve-1
      j   = jcel(1,it)
      i   = icel(1,it)
      j1  = jcel(2,it)
      jsf = j
      ic  = i
      zsv = fareo(it)    ! S^sita
      zsp = hvsb(jsf,i)  ! S^sita*h*cosw
      hcsw = zsp/zsv
      ia  = 1

!     22/08/29 yamamoto, flsmo_ave: calculate time average before accumulating flsmo values
!     22/08/29 yamamoot, flsmo_accumulate: Once the accumulation is taken, find the maximum and minimum after the accumulation.
      do k = 1,howmuchAcc
        flsmo_ave(k) = sum(timeSeries_o
     >   (timeNum-seriesNow+1:timeNum,k,it-1))/seriesNow
        do i_time = 1,timeNum
         flsmo_accumulate(i_time,k)
     >    = sum(timeSeries_o(i_time,1:k,it-1))
        end do
      end do

      write(i6,'(i3,0p4f10.6,1p53e11.3)')
     >  ic, rhf_odp(ic), rhf_idp(ic), rhf_omd(ic), rhf_imd(ic),
     >   hpit(jsf,i), hcsw/hpit(jsf,i), hcsw,
     >  flsmo_ave(1), !qt_ev
     >  sum(flsmo_ave(1:2)), !qt_ed
     >  sum(flsmo_ave(1:3)), !qt_iv
     >  sum(flsmo_ave(1:4)), !qt_id
     >  sum(flsmo_ave(1:5)), !qt_rec
     > (flsmo(it,k)/zsp*hcsw,k=6,7),
     > (flsmo(it,4)-flsmo(it,2))/zsp,fbdpo(it,m3)/zsp,
     >  flsmo(it,2)/zsp, fbdpo(it,m4)/zsp, 
     >  (sum(timeSeries_o(timeNum-seriesNow+1:timeNum,k,it-1))
     >  /seriesNow,k=6,10),
     >  zsv

      write(i7,'(i3,0p4f10.6,1p53e11.3)')
     >  ic, rhf_odp(ic), rhf_idp(ic), rhf_omd(ic), rhf_imd(ic),
     >  flsmo_ave(1), !qt_ev
     >  sum(flsmo_ave(1:2)), !qt_ed
     >  sum(flsmo_ave(1:3)), !qt_iv
     >  sum(flsmo_ave(1:4)), !qt_id
     >  sum(flsmo_ave(1:5)), !qt_rec
     >  (sum(timeSeries_o(timeNum-seriesNow+1:timeNum,k,it-1))
     >  /seriesNow,k=6,10),
     >  (maxval(flsmo_accumulate(timeNum-seriesNow+1:timeNum,k)),
     >   minval(flsmo_accumulate(timeNum-seriesNow+1:timeNum,k)),
     >   k=1,5),
     >  (maxval(timeSeries_o(timeNum-seriesNow+1:timeNum,k,it-1)),
     >  minval(timeSeries_o(timeNum-seriesNow+1:timeNum,k,it-1)),
     >  k=6,10)
      enddo
      close(i6)
      close(i7)
!
      dsn = "Qdpl_i.txt"
      if( ctime .ne. ' ' ) then
        dsn = trim( cpath ) // "Qdpl_i" 
     >   // '_' // trim( ctime ) // ".txt"
      endif
      open(unit=i6,file=trim(dsn),position="append")
!
      dsn2 = "Qdpl_i_stat.txt"
      if( ctime .ne. ' ' ) then
        dsn2 = trim( cpath ) // "Qdpl_i_stat" 
     >   // '_' // trim( ctime ) // ".txt"
      endif
      open(unit=i7,file=trim(dsn2),position="append")
!
      write(i6,'(//2x,"*** plhflx_idp ***  itim =",i7,"  time =",
     > 1pe14.6,"  dtim =",1pe12.3,"  nlp =",i3)') itim,time,dtim,nlp
      write(i6,'(2x,"total q_heat =",2x,"e-conV",5x,"e-conD",5x,
     > "i-conV",5x,"i-conD",5x,"recomb",5x,"radiation",2x,"neutral")')
      write(i6,'(2x,"total q_heat =",1p7e11.3)') (tflsmi(k),k=1,7)

      write(i7,'(//2x,"*** plhflx_idp ***  itim =",i7,"  time =",
     > 1pe14.6,"  dtim =",1pe12.3,"  nlp =",i3)') itim,time,dtim,nlp
      write(i7,'(2x,"total q_heat =",2x,"e-conV",5x,"e-conD",5x,
     > "i-conV",5x,"i-conD",5x,"recomb",5x,"radiation",2x,"neutral")')
      write(i7,'(2x,"total q_heat =",1p7e11.3)') (tflsmi(k),k=1,7)

!    220826 output time_average, maximum and minimum of qt_ev","qt_ed","qt_iv","qt_id","qt_rec","Ned","Ted","Tid","Vid", and "Vni
      write(i6,'(1x,"ic",2x,"r_odp",5x,"r_idp",5x,"r_omd",5x,"r_imd",
     > 5x,"pitch",6x,"csw",8x,"hcsw",7x,"qt_ev",6x,"qt_ed",6x,"qt_iv",
     > 6x,"qt_id",6x,"qt_rec",5x,"qt_rad",5x,"qt_ntl",
     > 5x,"qpi_sm",5x,"qpi_bc",5x,"qpe_sm",5x,"qpe_bc",5x,
     > "Ned",8x,"Ted",8x,"Tid",8x,"Vid",8x,"Vni",
     > 8x,"Are")')

      write(i7,'(1x,"ic",2x,"r_odp",5x,"r_idp",5x,"r_omd",5x,"r_imd",
     > 5x,"qt_ev",6x,"qt_ed",6x,"qt_iv",
     > 6x,"qt_id",6x,"qt_rec",5x,
     > "Ned",8x,"Ted",8x,"Tid",8x,"Vid",8x,"Vni",
     > 8x,"qt_ev_max",2x,"qt_ev_min",
     > 2x,"qt_ed_max",2x,"qt_ed_min",
     > 2x,"qt_iv_max",2x,"qt_iv_min",
     > 2x,"qt_id_max",2x,"qt_id_min",
     > 2x,"qt_rec_ma",2x,"qt_rec_mi",
     > 2x,"Ned_max",4x,"Ned_min",
     > 4x,"Ted_max",4x,"Ted_min",
     > 4x,"Tid_max",4x,"Tid_min",
     > 4x,"Vid_max",4x,"Vid_min",
     > 4x,"Vni_max",4x,"Vni_min",
     > 4x, "(max and min are the maximum 
     > and minimum values of the last few time steps.)")')
!
      do it = 2, itpve-1
      jmax = jtmax(it)
      j   = jcel(jmax,it)
      j1  = jcel(jmax-1,it)
      i   = icel(jmax,it)
      jsf = jcel(jmax-1,it)
      ic  = i
      zsv = farei(it)    ! S^sita
      zsp = hvsb(jsf,i)  ! S^sita*h*cosw
      hcsw = zsp/zsv
      ia  = 1

!     22/08/29 yamamoto, flsmi_ave:Calculate time averages before accumulating flsmi values
!     22/08/29 yamamoto, flsmi_accumulate:Once the accumulation is taken, the maximum and minimum after the accumulation are obtained.
      do k = 1,howmuchAcc
        flsmi_ave(k) = sum(timeSeries_i
     >   (timeNum-seriesNow+1:timeNum,k,it-1))/seriesNow
        do i_time = 1,timeNum
            flsmi_accumulate(i_time,k)
     >      = sum(timeSeries_i(i_time,1:k,it-1))
        end do
      end do

      write(i6,'(i3,0p4f10.6,1p50e11.3)')
     >  ic, rhf_odp(ic), rhf_idp(ic), rhf_omd(ic), rhf_imd(ic),
     >   hpit(jsf,i), hcsw/hpit(jsf,i), hcsw, 
     >  flsmi_ave(1), !qt_ev
     >  sum(flsmi_ave(1:2)),!qt_ed
     >  sum(flsmi_ave(1:3)), !qt_iv
     >  sum(flsmi_ave(1:4)), !qt_id
     >  sum(flsmi_ave(1:5)), !qt_rec
     > (flsmi(it,k)/zsp*hcsw,k=6,7),
     > (flsmi(it,4)-flsmi(it,2))/zsp,fbdpi(it,m3)/zsp,
     >  flsmi(it,2)/zsp, fbdpi(it,m4)/zsp,
     >  (sum(timeSeries_i(timeNum-seriesNow+1:timeNum,k,it-1))
     >  /seriesNow,k=6,10),
     >  zsv

      write(i7,'(i3,0p4f10.6,1p50e11.3)')
     >  ic, rhf_odp(ic), rhf_idp(ic), rhf_omd(ic), rhf_imd(ic),
     >  flsmi_ave(1), !qt_ev
     >  sum(flsmi_ave(1:2)),!qt_ed
     >  sum(flsmi_ave(1:3)), !qt_iv
     >  sum(flsmi_ave(1:4)), !qt_id
     >  sum(flsmi_ave(1:5)), !qt_rec
     >  (sum(timeSeries_i(timeNum-seriesNow+1:timeNum,k,it-1))
     >  /seriesNow,k=6,10),
     >  (maxval(flsmi_accumulate(timeNum-seriesNow+1:timeNum,k)),
     >   minval(flsmi_accumulate(timeNum-seriesNow+1:timeNum,k)),
     >   k=1,5),
     >  (maxval(timeSeries_i(timeNum-seriesNow+1:timeNum,k,it-1)),
     >  minval(timeSeries_i(timeNum-seriesNow+1:timeNum,k,it-1)),
     >  k=6,10)
      enddo
      close(i6)
      close(i7)
!
      return
      end
