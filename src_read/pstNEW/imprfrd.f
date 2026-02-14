!***********************************************************************
      subroutine imprfrd
!***********************************************************************
      use csize
      use cimcom
      use cimden
      use cntcom
      use cntmnt
      use cntpls
      use cplcom
      use cplimp
      use cplmet
      use cplpst
      use cplwrd
      use cplwrd2
      use cpmpls
      use cgdcom
      use cunit
      implicit none

!
!::local variables
      integer ix, iy, i6, n6sav
      real*8  xaxs, yaxs
      integer jc, n, k
      real*8  rst(ndy), ren(ndy), rhf(ndy)
      data i6/6/; save i6
      logical :: is_midpl = .false.
!
!::axis
      ix = mqs1
      iy = npmx
      xaxs = grdx(ix,iy)
      yaxs = grdy(ix,iy)
!
      write(n6,'(/2x,"*** imprfrd ***  wfac =",f9.5)') wfac
      write(n6,'(2x,"axis point =",2f8.4,"  out-mid =",2f8.4,2i5,
     >  "  in-mid =",2f8.4,2i5)') xaxs,yaxs,xmd1,ymd1,jmd1,imd1
     >  ,xmd2,ymd2,jmd2,imd2
      write(n6,'(2x,"R0 =",f8.4,"  a =",f8.4)') r0mp,ramp
      write(n6,'(2x,"icwl1 =",i3,"  icspx =",i3,"  icmps =",i3
     > ,"  icmpe =",i3,"  icaxs =",i3)') icwl1,icspx,icmps,icmpe,icaxs
!
!::radial position
      jc = 1
      call codrad(jc,rst,ren,rhf_odp,n,ndy)
      jc = jcmax
      call codrad(jc,rst,ren,rhf_idp,n,ndy)
      jc = jmd1
      call codrad(jc,rst,ren,rhf_omd,n,ndy)
      jc = jmd2
      call codrad(jc,rst,ren,rhf_imd,n,ndy)
      ! jc = 18
      ! call codrad(jc,rst,ren,rhf_uodp,n,ndy)
      ! jc = 102
      ! call codrad(jc,rst,ren,rhf_uidp,n,ndy)
!
!::radial profile
      do k = 1, 6
      if( k.eq.1 ) then
        jc = 1
        is_midpl = .false.
      elseif( k.eq.2 ) then
        jc = 2
        is_midpl = .false.
      elseif( k.eq.3 ) then
        jc = jmd1
        is_midpl = .true.
      elseif( k.eq.4 ) then
        jc = jmd2
        is_midpl = .true.
      elseif( k.eq.5 ) then
        jc = jcmax-1
        is_midpl = .false.
      elseif( k.eq.6 ) then
        jc = jcmax
        is_midpl = .false.
      endif
      ! Yamoto
      ! if( k.eq.7 ) jc = 18
      ! if( k.eq.8 ) jc = 102
      call lst_rad(jc,is_midpl)
      enddo
!
      return
      end
!
!***********************************************************************
      subroutine lst_rad(kk,is_midpl)
!***********************************************************************
      use csize
      use cimcom
      use cimden
      use cntcom
      use cntmnt
      use cntpls
      use cplcom
      use cplimp
      use cplmet
      use cplpst
      use cplwrd
      use cplwrd2
      use cpmpls
      use cgdcom
      use cunit
      use cplimp_plimp
      implicit none
!
      integer, intent(in) :: kk
      logical, intent(in) :: is_midpl
!
      integer :: ix, iy, ia, iy1, iy2, i6, ic, iz, j, i, i6y
      real(8) :: zni, zne, zef, z0, ztz
      real(8), dimension(1:nzsmx) :: zav, znzi, znzt
      real(8) :: undf
      real(8) :: zef_bulk
      real(8) :: vwrt, vwrc, vwrm, vwrm0(1:nzsmx), wrdm
      real(8) :: tfrf, tthf, tfoc, tvel
      real(8) :: eff_counter
      real(8), dimension(1:nzsmx) :: zeff_ech
      real(8), dimension(1:nzsmx) :: tnfrz, tnthz, tnfoc, tnvlz
! modified 2/2 lines treat 4 or more impurities with IMPMC
!      real(8), dimension(0:ndis2L) :: hnz
!      real(8), dimension(0:ndis2L) :: hfrz, hthz, hvlz
      real(8), dimension(0:ndmis) :: hnz
      real(8), dimension(0:ndmis) :: hfrz, hthz, hvlz
      character(80) :: drimp, cdsn
      integer :: nty
      real(8) :: Nzxa(wmc_nty), Ezxa(wmc_nty), Pzxa(wmc_nty)
     >         , Nzxi(wmc_nty), Ezxi(wmc_nty), Pzxi(wmc_nty)
      real(8) :: odp_tmp, idp_tmp
!
      ix = kk
      write(n6,'(2x,"*** lst_rad ***  ix =",i3)') ix
      if( ix.le.0 .or. ix.gt.jcmax ) return
!
!::file name
      i6 = 21
      call getenv("IMP1D",drimp)
      write(cdsn,'(a,i3.3,a)') trim(drimp) // "/PRF_j", ix, ".txt"
      call delspc(cdsn)
      open(unit=i6, file=cdsn)
      write(n6,'(5x,"cdsn = ",a)') trim(cdsn)
! Yamoto additional outputs
      i6y = 504
      call getenv("IMP1D",drimp)
      write(cdsn,'(a,i3.3,a)') trim(drimp) // "/PRF_yj", ix, ".txt"
      call delspc(cdsn)
      open(unit=i6y, file=cdsn)
      write(n6,'(5x,"cdsn = ",a)') trim(cdsn)
!
!::plasma parameter
      ia = 1
      iy1 = icmin(ix)
      iy2 = icmax(ix)
      write(i6,'(2x,"*** radial profile ***  ix =",i3,"  iy =",2i5,
     >  "  wfac =",f9.5)') ix, iy1, iy2, wfac
      write(i6,'(4x,"ic",4x,"iy",2x,
     >  "r_odp",6x,"r_idp",6x,"r_omd",6x,
     >  "r_imd",6x,"r_uodp",5x,"r_uidp",5x,
     >  "dene",7x,"vnezef",5x,"zne",8x,"deni",7x,"teme",7x,
     >  "temi",7x,"Vf",9x,"vwrt",7x,"vwrm",7x,"Zeff",
     >  3((7x,a5,i1,5x,a2,i1,a4,4x,a2,i1,a4,4x,a2,i1,a1),
     >      6(7x,a2,i1,a1)))')
     >   ("wrdmz",nty, "Nz",nty,"i/Ni"
     >     ,"Nz",nty,"t/Ni", "<Z",nty,">"
     >     ,"Nz",nty,"a", "Ez",nty,"a", "Pz",nty,"a"
     >     ,"Nz",nty,"i", "Ez",nty,"i", "Pz",nty,"i"
     >     ,nty=1,wmc_nty)

      write(i6y,'(2x,"*** radial profile ***  ix =",i3,"  iy =",2i5,
     >  "  wfac =",f9.5)') ix, iy1, iy2, wfac
      write(i6y,'(4x,"ic",4x,"iy",2x,"r_odp",6x,"r_idp",6x,"r_omd",6x,
     >  "r_imd",6x,"r_uodp",5x,"r_uidp",5x,
     >  "Zeff",7x,"Zeff_bulk",2x,3(4(a5,6x)),3(4(a5,6x)))')  
     >     (("Zeff2"),nty=1,wmc_nty),
     >     (("Fric4"),nty=1,wmc_nty),
     >     (("Ther4"),nty=1,wmc_nty),
     >     (("Tot_4"),nty=1,wmc_nty),
     >     (("fri_A"),nty=1,wmc_nty),
     >     (("thf_A"),nty=1,wmc_nty),
     >     (("Tot_A"),nty=1,wmc_nty),
     >     (("vlz_A"),nty=1,wmc_nty)

!
      do iy = iy1, iy2
      ic = mcel(ix,iy)
!
!::sol/prv wall
      if( ic.eq.0 ) then
      j = ix
      i = iy
      ia = 1
      z0 = 1.0d-30
      undf = 0.0d0
!!
      if(is_midpl .and. rhf_odp(iy) < 0) then
        odp_tmp = 0.0_8
      else
        odp_tmp = rhf_odp(iy)*100.0d0
      endif
!!
      if(is_midpl .and. rhf_idp(iy) > 0) then
        idp_tmp = 0.0_8
      else
        idp_tmp = rhf_idp(iy)*100.0d0
      endif
!!
      write(i6,'(1x,i5,1x,i5,1p50e11.3)') ic, iy,
     >  odp_tmp, idp_tmp,
     >  rhf_omd(iy)*100.0d0, rhf_imd(iy)*100.0d0,
     >  rhf_uodp(iy)*100.0d0, rhf_uidp(iy)*100.0d0,
     > dmax1(z0,vne(j,i)), dmax1(z0,vnezef(j,i)), dmax1(z0,vne(j,i)), 
     > dmax1(z0,vna(j,i,ia)), dmax1(z0,vte(j,i)),dmax1(z0,vti(j,i)),
     > vva(j,i,ia)+z0,
     > dmax1(z0,undf),dmax1(z0,undf),dmax1(z0,undf),
     > (dmax1(z0,undf),dmax1(z0,undf),dmax1(z0,undf),
     >   dmax1(z0,undf),dmax1(z0,undf),dmax1(z0,undf),
     >   dmax1(z0,undf),dmax1(z0,undf),dmax1(z0,undf),dmax1(z0,undf)
     >   ,nty=1,wmc_nty)
      write(i6y,'(1x,i5,1x,i5,1p50e11.3)') ic, iy,
     >  odp_tmp, idp_tmp,
     >  rhf_omd(iy)*100.0d0, rhf_imd(iy)*100.0d0,
     >  rhf_uodp(iy)*100.0d0, rhf_uidp(iy)*100.0d0,
     > dmax1(z0,undf),dmax1(z0,undf),
     > ((dmax1(z0,undf)),nty=1,wmc_nty),
     > ((dmax1(z0,undf)),nty=1,wmc_nty),
     > ((dmax1(z0,undf)),nty=1,wmc_nty),
     > ((dmax1(z0,undf)),nty=1,wmc_nty),
     > ((dmax1(z0,undf)),nty=1,wmc_nty),
     > ((dmax1(z0,undf)),nty=1,wmc_nty),
     > ((dmax1(z0,undf)),nty=1,wmc_nty),
     > ((dmax1(z0,undf)),nty=1,wmc_nty)
      goto 150
      endif ! ic.eq.0
!
!::monte cell
      j  = iplx(ic)
      i  = iply(ic)
!
      if( j.eq.ix .and. i.eq.iy ) goto 120
      write(i6,'(2x,"wrong index  ix,iy =",2i5,"  ic =",i7," j,i ="
     >   ,2i5)') ix, iy, ic, j, i
      call wexit("imprfrd","wring index")
 120  continue
!
!::Ne, Nz, <Z>, Zef
      zne = 0.0d0
      zni = 0.0d0
      zef = 0.0d0
      zeff_ech(:) = 0.0d0
      do ia = 1, nion
      zne = zne + aza(ia)*deni(ic,ia)
      zni = zni + deni(ic,ia)
      zef = zef + aza(ia)**2*deni(ic,ia)
      enddo
      zef_bulk = zef
!
      ztz = 0.0d0
      zav = 0.0d0
!
      do nty = 1, wmc_nty
        ztz = 0.0d0 ! SY zeroset of ztz should be placed here
        zeff_ech(nty)=0.0d0 ! yamoto
        hnz(0:ismaxL(nty)) = tdnzL(0:ismaxL(nty),ic,nty) 
                        ! wfac(wfinp) already multiplyed in set_imp
        hfrz(0:ismaxL(nty)) = tfrzL(0:ismaxL(nty),ic,nty)
        hthz(0:ismaxL(nty)) = tthzL(0:ismaxL(nty),ic,nty)
        hvlz(0:ismaxL(nty)) = tvlzL(0:ismaxL(nty),ic,nty)
        tfrf = 0.0d0
        tthf = 0.0d0
        tfoc = 0.0d0
        tvel = 0.0d0
        eff_counter = 0.0d0
        do iz = 1, ismaxL(nty)
          ztz = ztz + hnz(iz)
          zne = zne + dfloat(iz)*hnz(iz)
          zav(nty) = zav(nty) + dfloat(iz)*hnz(iz)
          zef = zef + dfloat(iz)**2*hnz(iz)
          zeff_ech(nty) = zeff_ech(nty) + dfloat(iz)**2*hnz(iz)
          tfrf = tfrf + hfrz(iz)/(dfloat(iz)**2)
          tthf = tthf + hthz(iz)/(dfloat(iz)**2)
          tfoc = tfoc + hfrz(iz)/(dfloat(iz)**2)
     >      +hthz(iz)/(dfloat(iz)**2)
          tvel = tvel + hvlz(iz)
          if(hvlz(iz).ne.0.0d0) eff_counter=eff_counter+1.0d0
        enddo

        znzi(nty) = ztz             ! ion
        znzt(nty) = ztz + hnz(0)     ! total

        if( ztz.ne.0.0d0 ) zav(nty) = zav(nty)/ztz
        znzi(nty) = znzi(nty)/zni
        znzt(nty) = znzt(nty)/zni
        if(eff_counter.ne.0.0d0)then
        tnfrz(nty) = tfrf/eff_counter
        tnthz(nty) = tthf/eff_counter
        tnfoc(nty) = tfoc/eff_counter
        tnvlz(nty) = tvel/eff_counter
        else
        tnfrz(nty) = 0.0d0
        tnthz(nty) = 0.0d0
        tnfoc(nty) = 0.0d0
        tnvlz(nty) = 0.0d0
        endif
      enddo ! nty

      if( zne.ne.0.0d0 ) zef = zef/zne
      if( zne.ne.0.0d0 ) zef_bulk = zef_bulk/zne
      if( zne.ne.0.0d0 ) zeff_ech(:) = zeff_ech(:)/zne
!
!::radiation
!  v : soldor variables  t:total m:monte c:corona   vwrm = wrdm
      vwrt = -wime(ic)  ! total radiation in soldor
      vwrc = -wcre(ic)  ! Corona Model in soldor
      vwrm = -wmce(ic)  ! MC radiation in soldor
      vwrm0(1:wmc_nty)= -wmc_zwe(ic,1:wmc_nty)  ! Each MC radiation in soldor
!      wrdm = twrd(ic)   ! MC radiation in IMPMC 
!
      do nty = 1, wmc_nty
        Nzxa(nty) = tdnzL(0,ic,nty)
        Ezxa(nty) = tengzL(0,ic,nty)
        Pzxa(nty) = tprzL(0,ic,nty)
        Nzxi(nty) = 0.0d0
        Ezxi(nty) = 0.0d0
        Pzxi(nty) = 0.0d0
        do iz = 1, ndmis
          Nzxi(nty) = Nzxi(nty) + tdnzL(iz,ic,nty)
          Ezxi(nty) = Ezxi(nty) + tengzL(iz,ic,nty)
          Pzxi(nty) = Pzxi(nty) + tprzL(iz,ic,nty)
        enddo
      enddo
!
      ia = 1
      z0 = 1.0d-30
!!
      if(is_midpl .and. rhf_odp(iy) < 0) then
        odp_tmp = 0.0_8
      else
        odp_tmp = rhf_odp(iy)*100.0d0
      endif
!!
      if(is_midpl .and. rhf_idp(iy) > 0) then
        idp_tmp = 0.0_8
      else
        idp_tmp = rhf_idp(iy)*100.0d0
      endif
!!
      write(i6,'(1x,i5,1x,i5,1p50e11.3)') ic, iy,
     >  odp_tmp, idp_tmp,
     >  rhf_omd(iy)*100.0d0, rhf_imd(iy)*100.0d0,
     >  rhf_uodp(iy)*100.0d0, rhf_uidp(iy)*100.0d0,
     > dmax1(z0,dene(ic)), dmax1(z0,vnezef(j,i)), dmax1(z0,zne), 
     > dmax1(z0,deni(ic,ia)), dmax1(z0,teme(ic)),dmax1(z0,temi(ic)),
     > vflw(ic,ia)+z0,
     > dmax1(z0,vwrt), dmax1(z0,vwrm),dmax1(z0,zef),
     > (dmax1(z0,vwrm0(nty)),dmax1(z0,znzi(nty)),dmax1(z0,znzt(nty)),
     >   dmax1(z0,zav(nty)),
     >   Nzxa(nty),Ezxa(nty),Pzxa(nty),Nzxi(nty),Ezxi(nty),Pzxi(nty)
     >   ,nty=1,wmc_nty)

      write(i6y,'(1x,i5,1x,i5,1p50e11.3)') ic, iy,
     >  odp_tmp, idp_tmp,
     >  rhf_omd(iy)*100.0d0, rhf_imd(iy)*100.0d0,
     >  rhf_uodp(iy)*100.0d0, rhf_uidp(iy)*100.0d0,
     > dmax1(z0,zef), dmax1(z0,zef_bulk),
     > ((dmax1(z0,zeff_ech(nty))),nty=1,wmc_nty),
     > ((tfrzL(4,ic,nty)+z0),nty=1,wmc_nty),
     > ((tthzL(4,ic,nty)+z0),nty=1,wmc_nty),
     > ((tfrzL(4,ic,nty)+tthzL(4,ic,nty)+z0),nty=1,wmc_nty),
     > ((tnfrz(nty)+z0),nty=1,wmc_nty),
     > ((tnthz(nty)+z0),nty=1,wmc_nty),
     > ((tnfoc(nty)+z0),nty=1,wmc_nty),
     > ((tnvlz(nty)+z0),nty=1,wmc_nty)
 150  continue
      enddo !iy
      close(i6)
      close(i6y)
!
      return
      end
!
!***********************************************************************
      subroutine codrad(jc,rst,ren,rhf,n,nd)
!***********************************************************************
!
!   cell name in metric calculation
!
!                       +icaxs--------------------+
!                       !                         !
!                       !icmpe~~~~~~~icmax(jc)~~~~!
!      icwl2 /+^^^^^^^^^+                         +^^^^^^^^^+/
!            /!  (4)    !icmps  Main plasma (6)   !    (5)  !/
!            /!---------!-------------------------!---------!/
!      icspx /!  OUT    !                         !  IN     !/
!            /!   DIV   !      Scrape-Off         !   DIV   !/
!            /!  (1)    !                  (2)    !    (3)  !/
!   A  icwl1 /j+^^^^^^^j^^^^^^^^^^^^^icmin(jc)^^^^^j^^^^^^^^j j
! i !         c        c                           c        c c
! c !         d        x                           x        d m
!   !-->    1=p        p                           p        p=a
!     jc      1        1                           2        2 x
!
!                do jc = 1, jcmax                ir = kreg(jc,ic)
!                do ic = icmin(jc), icmax(jc)
!
!----------------------------------------------------------------------
      use csize
      use cimcom
      use cimden
      use cntcom
      use cntmnt
      use cntpls
      use cplcom
      use cplimp
      use cplmet
      use cplpst
      use cplwrd
      use cplwrd2
      use cpmpls
      use cgdcom
      use csonic
      implicit none
!
!::argument
      integer jc, n, nd
      real*8  rst(nd), ren(nd), rhf(nd)
!
!::local variables
      integer  ic1, ic2, ic, i6, kin
      integer  mx1, mx2, mx3, mx4, my1, my2, my3, my4
      real*8   sum, ra, rb, rc, rd, za, zb, zc, zd, rp, zp, rq, zq, dl
!
      data i6/111/; save i6
!
!xx      write(i6,'(/2x,"*** codrad ***    jc =",i4)') jc
      ic1 = icmin(jc)
      ic2 = icmax(jc)
!xx      write(i6,'(2x,"jc =",i4,"  ic1 =",i4,"  ic2 =",i4)') jc,ic1,ic2      
!
      do ic = 1, nd
      rst(ic) = 0.0d0
      ren(ic) = 0.0d0
      rhf(ic) = 0.0d0
      enddo
!
      sum = 0.0d0
      do ic = ic1, ic2
      call mcpnt(jc,ic,mx1,mx2,mx3,mx4,my1,my2,my3,my4)
!
      ra = grdx(mx1,my1)
      rb = grdx(mx2,my2)
      rc = grdx(mx3,my3)
      rd = grdx(mx4,my4)
      za = grdy(mx1,my1)
      zb = grdy(mx2,my2)
      zc = grdy(mx3,my3)
      zd = grdy(mx4,my4)
!
      rp = 0.5d0*(ra+rb)
      zp = 0.5d0*(za+zb)
      rq = 0.5d0*(rc+rd)
      zq = 0.5d0*(zc+zd)
!
      dl = dsqrt( (rq-rp)**2+(zq-zp)**2)
      rst(ic) = sum
      ren(ic) = sum + dl
      rhf(ic) = sum + 0.5d0*dl
      sum = sum + dl
      enddo
!
      kin = 0
      if( jc.ge.(jcxp1+jcxp2)/2 ) kin = 1
      do ic = ic1, ic2
      rhf(ic) = ren(icspx) - rhf(ic)
      if( kin.eq.1 ) rhf(ic) = -rhf(ic)
      enddo
      n = ic2
!
!xx      do ic = ic1, ic2
!xx      write(i6,'(2x,i4,1p3e12.4)') ic, rst(ic), ren(ic) , rhf(ic)
!xx      enddo
!
      return
      end
