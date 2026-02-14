!**********************************************************************
      subroutine imp_cutbk
!**********************************************************************
!)
!)   cutbk : reduction in number of test particles in core region
!)
!)    Note   funroh   I made a mistake. X  roh  O rho
!)
!)  chro(npmax)   ! rho of ion particle in main plasma c: cutback
!)  nhr = 20      ! division of rho [0.0,1.0]
!)  hdr = 0.05d0  ! jr = crho(ip)/hdr + 1.0d0
!)  hro(nhr)      ! hro(1) = 0.05, hro(2) = 0.10, ...  hro(nhr) = 1.0
!)  hstz(jr,iz)   ! particle number with charge state iz
!)
!)  Edit SONIC =>  MPMPD system
!)    lnope, lmype   =>  nope_grp, nrnk_grp
!)    use mpi
!)    use mod mpicomm   ! sonicV4/MPMD5/WFL/mod_mpicomm
!)   >  ,n5M => n5, n6M => n6, ngrpM => ngrp, mygrpM => mygrp
!)
!)  imp_pack exclude particles with ien(ip) >= 6
!)  modified imp_pack from Shimizu-san's version
!)  previous imp_pack only kills ions
!)  terminate trace  if(ien(ip) >= 6) goto 310  IMPTpg/imp_step.f
!)
!)  flowchart
!)       1PE  sv-data => TPE sw-data => sort,pair,roullet,new wght
!)    => 1PE  sv-data => new wght
!)
!)   chaeck nsw(sw-data), npr(pair), ntr(trim)
!)
!)  key word : algo : verify algotithm
!)           : npr : paring  : nrdc : reduction
!)
!)--------------------------------------------------------------------
      use cimcom,      only : crho, cutiz, cutjr, hdr, hroh, hstz
     >    , i6_scat, ien, ir, is, ismax, jrmax, lcutz, ndhr, ndis
     >    , ndmct
     >    , ndmp, nhr, nhsiz, npcut, npmax, pemt, rr, stimz, wght, zz
      use cunit,       only : n6
      use mod_mpicomm, only : nope_grp, nrnk_grp
      implicit none

!::local common  important !! need package due to MPI_Send/Recv
!::svXXX (1PE: GM)
      integer, parameter :: ndmsv = 500
      real(8) :: svnpn
      real(8), dimension(ndmsv) :: svrnk, svipn, svpem, svwgt, svwgc

!::swXXX (TPE: GL)
      integer, parameter :: ndmsw = 500*10
      real(8) :: swnpn
      real(8), dimension(ndmsw) :: swrnk, swipn, swpem, swwgt, swwgc
      integer, dimension(ndmsw) :: iwsrt  ! srt: sort
      real(8), dimension(ndmsw) :: wtsrt

      integer, parameter :: ndmcpe = 10000 !maximum irnk
      real(8) :: vrnkp(0:ndmcpe), vrnkw(0:ndmcpe)
      real(8) :: vcutp(0:ndmcpe), vcutw(0:ndmcpe)
      real(8), dimension(3) :: dltpw, tdlpw

!::variables forcheck
      integer :: i, ip, ic, ko, jr, iz, ik, ist, jst, j, ii
      integer :: nsv, nsw, jr0, iz0, nsw0, iw, iv, nswh
      integer :: itrm, i6
      integer :: jw1, jw2, iw1, iw2, iwL, iwD, idbg
      integer :: irnk1, irnk2, irnk
      integer :: nct, ict, totp
      real(8) :: x0, y0, zr, eps, xran, pcn
      real(8) :: totw
      integer :: ierr
      real(8) :: wgh0(ndmp)
! function
      integer    imox, imoy
      real(8)    funroh, random

!::demension check
      if( nope_grp > ndmcpe ) then
        write(n6,'(2x,"nope_grp > ndmcpe  ",2i6)') nope_grp, ndmcpe
        call wexit("imp_cutbk","nope_grp > ndmcpe")
      endif

!::input data
!xx   lcutz = 1   ! on/off cutbk model
!xx   nhr   = 50  ! histgram in jr  < ndhr = 50
!xx   jrmax = 48  ! exclude vicinity of separatrix
!xx   npcut = 300 ! target of reduction particles
!xx   nhsiz = 10  ! limitation of histgram size

      write(n6,'(/2x,"*** imp_cutbk ***  npmax =",i6,
     >  "  lcutz =",i2)') npmax, lcutz
      if( lcutz == 0 ) return

      write(n6,'(2x,"npcut =",i6,"  jrmax/nhr/ndhr =",3i5,
     >  "  nhsiz =",i5)')  npcut, jrmax, nhr, ndhr, nhsiz
      if( jrmax > nhr .or. nhr > ndhr ) then
        call wexit("imp_cutbk","jrmax > nhr .or. nhr > ndhr")
      endif

!::initialize for avoid error
      ist = 1
      jst = 1

!::check ien(ip)
      ierr = 0
      do ip = 1, npmax
        if( ien(ip) /= 0 ) ierr = ierr + 1
      enddo
      write(n6,'(2x,"number of ien(ip) /= 0  ierr =",i5)') ierr
      if( ierr > 0 ) then
        call lst_open("scat")
        i6 = i6_scat
        call lst_scat(stimz,1,i6)  ! error
        call wexit("imp_cutbk","find ien /= 0")
      endif

!::old wght before cutback
      wgh0(1:npmax) = wght(1:npmax)

!::real => integer  hstz(jr,iz), svnpn, svrnk(iv), svipn(iv)
      eps = 1.0d-3

!::def. rho
      hdr = 1.0d0/dfloat(nhr)
      do i = 1, nhr
        hroh(i) = hdr*dfloat(i)
      enddo

      crho(1:npmax) = 1.5d0
      hstz(1:ndhr,0:ndis) = 0.0d0

      write(n6,'(/2x,"particle inf.  1PE  nrnk =",i5)') nrnk_grp
      write(n6,'(2x,4x,"ip",2x,"ien",4x,"is",4x,"ic",6x,"ix",3x,"iy",
     >  2x,"ko",2x,"jro",2x,"rho",3x,"wght")')
      do ip = 1, npmax
        if( is(ip)  == 0 ) cycle
        ic = ir(ip)
        iz = is(ip)
        x0 = rr(ip)
        y0 = zz(ip)
        zr = funroh(ic,x0,y0)
        if( zr >= 1.0d0 ) cycle
        call mchkin(x0,y0,ic,ko)

        jr = int(zr/hdr + 1.0d0)
        crho(ip) = zr
        hstz(jr,iz) = hstz(jr,iz) + 1.0d0

        if( ip <= 3 .or. ip > npmax-3 ) then
        write(n6,'(2x,3i6,i7,i6,i5,i5,i5,2f11.5)') ip, ien(ip),
     >  is(ip), ic, imox(ic), imoy(ic), ko, jr, crho(ip), wght(ip)
        endif
      enddo

!::real hstz => int(hstz)
      do jr = 1, nhr
      do iz = 1, ismax
        if( hstz(jr,iz) > 0.0d0 ) hstz(jr,iz) = hstz(jr,iz) + eps
      enddo
      enddo

!::histgram
      do iz = 1, ismax
        ik = 0
        do jr = 1, nhr
          if( hstz(jr,iz) > dble(nhsiz) ) ik = ik + 1
        enddo
        if( ik > 0 ) then
          ist = iz-1
          exit
        endif
      enddo
      do jr = 1, nhr
        ik = 0
        do iz = 1, ismax
          if( hstz(jr,iz) > dble(nhsiz) ) ik = ik + 1
        enddo
        if( ik > 0 ) then
          jst = jr-1
          exit
        endif
      enddo

      write(n6,'(/2x,"histgram  hstz(jr,iz) 1PE  nrnk =",i5)') nrnk_grp
      write(n6,'(3x,"jr",5x,"hroh",5x,20(3x,"iz",i2.2,1x:))')
     >    (i,i=ist,ismax)
      do jr = jst, nhr
        write(n6,'(2x,i3,2x,f10.5,20i8)')
     >    jr, hroh(jr), (int(hstz(jr,i)),i=ist,ismax)
      enddo

!::sumup
      write(n6,'(/2x,"histgram  hstz(jr,iz) TPE  nope =",i5)') nope_grp
      write(n6,'(3x,"jr",5x,"hroh",5x,20(3x,"iz",i2.2,1x:))')
     >    (i,i=ist,ismax)
      do jr = jst, nhr
        write(n6,'(2x,i3,2x,f10.5,20i8)')
     >    jr, hroh(jr), (int(hstz(jr,i)),i=ist,ismax)
      enddo

!::select (jr,iz)
        ii = 0
        do iz = 1, ismax
        do j = jrmax, 1, -1
          jr = j
          if( hstz(jr,iz) > dfloat(npcut) ) then
            ii = ii + 1
            if( ii > ndmct ) goto 920
            cutjr(ii) = jr
            cutiz(ii) = iz
          endif
        enddo
        enddo
        nct = ii

!::cut-condition
      write(n6,'(/2x,"cut-conditon nct =",i4,"  TPE  nope =",i5)')
     >     nct, nope_grp
      if( nct == 0 ) then
        write(n6,'(2x,"not execute imp_cutbk due to nct = 0")')
        goto 200
      endif

      do ict = 1, nct
        jr = cutjr(ict)
        iz = cutiz(ict)
        write(n6,'(2x,i5,i5,i5,f11.5,i7)')
     >    ict, jr, iz, hroh(jr), int(hstz(jr,iz))
      enddo

!::rnk  sumup for cut-condition
      vcutp(0:ndmcpe) = 0.0d0
      vcutw(0:ndmcpe) = 0.0d0

!::cut condition
      do ict = 1, nct
        jr0 = cutjr(ict)
        iz0 = cutiz(ict)
        nsw0 = int(hstz(jr0,iz0))
        write(n6,'(/2x,80("-"))')
        write(n6,'(2x,"ict =",i3,"/",i3,"  jr0 =",i4,"  iz0 =",i4,
     >    "  nsw0 =",i6)') ict, nct, jr0, iz0, nsw0

!::rnk  every cut-condition
        vrnkp(0:ndmcpe) = 0.0d0
        vrnkw(0:ndmcpe) = 0.0d0

        iw = 0
        iv = 0
        do ip = 1, npmax
          zr = crho(ip)
          jr = int(zr/hdr + 1.0d0)
          if( jr == jr0 .and. is(ip) == iz0 ) then
            iv = iv + 1
            if( iv > ndmsv ) goto 940
            svrnk(iv) = nrnk_grp + eps
            svipn(iv) = ip       + eps
            svpem(iv) = pemt(ip)
            svwgt(iv) = wght(ip)
          endif
!::verify algorithm
        enddo
        svnpn = iv + eps
        nsv = int(svnpn)

!::debug write
        write(n6,'(2x,"sv-data  1PE  nrnk =",i4,"  nsv =",i5)')
     >    nrnk_grp, nsv
        write(n6,'(2x,"    iv    iv  svrnk   svipn  svpem     svwgt")')
        do iv = 1, nsv
          if( iv <= 3 .or. iv > nsv-3 ) then
            write(n6,'(2x,3i6,i8,1pe12.3,0pf10.5)') iv, iv,
     >        int(svrnk(iv)), int(svipn(iv)), svpem(iv), svwgt(iv)
          endif
        enddo
        write(n6,'(2x)')

!::sw-data
        nsv = int(svnpn)
        do j = 1, nsv
          iv = j
          iw = iw + 1
          if( iw > ndmsw ) goto 960
          swrnk(iw) = svrnk(iv)
          swipn(iw) = svipn(iv)
          swpem(iw) = svpem(iv)
          swwgt(iw) = svwgt(iv)
        enddo
        swnpn = iw + eps

!::debug write
        nsw  = int(swnpn)
        nswh = nsw/2
        write(n6,'(2x,"nsw/nsw0 =",2i6,"  nswh =",i6,"  TPE  nope =",
     >    i6)') nsw, nsw0, nswh, nope_grp
        if( nsw /= nsw0 ) then
          call wexit("imp_cutbk","nsw /= nsw0")
        endif
!
!::sort swwgt
        swwgc(1:nsw) = swwgt(1:nsw)
        wtsrt(1:nsw) = swwgt(1:nsw)

!::too many particles with wght=1  f: full
!::define wtsrt  weight for sorting
        write(n6,'(2x,"define wtsrt = wght*(1.0+0.01*random)")')
        do iw = 1, nsw
          wtsrt(iw) = wtsrt(iw)*(1.0d0+0.01d0*random(0))
        enddo
        call tsort(nsw,wtsrt,iwsrt)
        write(n6,'(2x,"sort wtsrt  nsw =",i5,"  TPE  nope =",i6)')
     >   nsw, nope_grp

!::coupling
        idbg = 0
        if( nsw < 150 ) idbg = 1
        idbg = 1
        if( idbg == 1 ) then
          write(n6,'(2x)')
          write(n6,'(2x,"coupling and Russian roulet  wght => wghtN")')
          write(n6,'(5x,a,6x,a)')
     >     "jw1   iw1   rnk1  ip1   pemt1       wght1     wghtN1",
     >     "jw2   iw2   rnk2  ip2   pemt2       wght2     wghtN2"
        endif

        jw1 = 1
        jw2 = nsw
 100    continue
        iw1 = iwsrt(jw1)
        iw2 = iwsrt(jw2)

!::russian roullete  l:life d:death
        xran = random(0)
        iwL = iw1
        iwD = iw2
        if( xran > 0.5d0 ) then
          iwL = iw2
          iwD = iw1
        endif
        pcn = swpem(iw1)*swwgt(iw1) + swpem(iw2)*swwgt(iw2)
        swwgc(iwL) = pcn/swpem(iwL)
        swwgc(iwD) = 0.0d0
      
!::balance
        irnk1 = int(swrnk(iw1))
        irnk2 = int(swrnk(iw2))
        if( swwgc(iw1) == 0.0d0 ) then
          vrnkp(irnk1) = vrnkp(irnk1) - 1
        elseif( swwgc(iw2) == 0.0d0 ) then
          vrnkp(irnk2) = vrnkp(irnk2) - 1
        endif
        vrnkw(irnk1) = vrnkw(irnk1) + swpem(iw1)*swwgc(iw1)
     >    - swpem(iw1)*swwgt(iw1)
        vrnkw(irnk2) = vrnkw(irnk2) + swpem(iw2)*swwgc(iw2)
     >    - swpem(iw2)*swwgt(iw2)

        if( idbg == 1 ) then
          if( jw1 < 5 .or. jw1 > nswh-5 ) then
            write(n6,'(2(3x,i5,i6,i6,i6,1pe12.3,0p2f10.5))')
     >        jw1, iw1, int(swrnk(iw1)), int(swipn(iw1)),
     >            swpem(iw1), swwgt(iw1), swwgc(iw1),
     >        jw2, iw2, int(swrnk(iw2)), int(swipn(iw2)),
     >            swpem(iw2), swwgt(iw2), swwgc(iw2)
          endif
        endif

        jw1 = jw1 + 1
        jw2 = jw2 - 1
        if( jw2-jw1 >= 1 ) goto 100

        if( jw2 == jw1 ) then
          iw1 = iwsrt(jw1)
          swwgc(iw1) = swwgt(iw1)
          if( idbg == 1 ) then
            write(n6,'(2(3x,i5,i6,i6,i6,1pe12.3,0p2f10.5))')
     >       jw1, iw1, int(swrnk(iw1)), int(swipn(iw1)),
     >            swpem(iw1), swwgt(iw1), swwgc(iw1)
          endif
        endif

!::sumup for cut-condition
        do i = 1, nope_grp
          irnk = i-1
          vcutp(irnk) = vcutp(irnk) + vrnkp(irnk)
          vcutw(irnk) = vcutw(irnk) + vrnkw(irnk)
        enddo

!::conservation wght
!::wght only conservation does not effective to evaluate
!::particle conservation when there are several sources (e.g. Arpf+Bk)
        write(n6,'(/2x,"Conservation of wght  TPE  nope =",i5)')
     >   nope_grp
        write(n6,'(2x," nrnk   dlpt     dlwg")')
        totp = 0
        totw = 0.0d0
        do i = 1, nope_grp
          irnk = i-1
          itrm = int(-vrnkp(irnk) + eps)
          totp = totp + itrm
          totw = totw + vrnkw(irnk)
          write(n6,'(2x,i5,i7,0pf12.5)') irnk, itrm, vrnkw(irnk)
        enddo
        write(n6,'(2x,a,i7,0pf12.5)') "tot  ", totp, totw

!xx   nsw  = number of particle with (jr0,iz0) in sw-data
!xx   nswh = nsw/2
!xx   npr  = number of paring
!xx   totp : integer

!::new wght
        nsw = int(swnpn)
        iv = 0
        do iw = 1, nsw
          if( int(swrnk(iw)) == nrnk_grp ) then
            iv = iv + 1
            ip = int(swipn(iw))
            wght(ip) = swwgc(iw)
          endif
        enddo

!::Note  recv sv-data on GL  overwrite
!::sw-data => sv-data
        iv = 0
        do iw = 1, nsw
          if( int(swrnk(iw)) == nrnk_grp ) then
            iv = iv + 1
            svrnk(iv) = swrnk(iw)
            svipn(iv) = swipn(iw)
            svpem(iv) = swpem(iw)
            svwgt(iv) = swwgt(iw)
          endif
        enddo
        nsv = iv

!::debug write
        write(n6,'(2x)')
        write(n6,'(2x,"sv-data  nrnk =",i4,"  nsv =",i5)') nrnk_grp,nsv
        write(n6,'(2x,a)')
     >   "    iv    iv  svrnk   svipn  svpem       svwgt"//
     >   "     wgh0  =>  wght"
        do iv = 1, nsv
          if( iv <= 3 .or. iv > nsv-3 ) then
            ip = int(svipn(iv))
          write(n6,'(2x,3i6,i8,1pe12.3,0p3f10.5)') iv, iv,
     >      int(svrnk(iv)), int(svipn(iv)), svpem(iv), svwgt(iv),
     >      wgh0(ip), wght(ip)
          endif
        enddo

      enddo   ! loop(ict)

!::var wght
 200  continue
      write(n6,'(2x,80("-"))')
      dltpw(1) = 0.0d0  ! ptcl
      dltpw(2) = 0.0d0  ! wght
      do ip = 1, npmax
        if( wght(ip) == 0.0d0 ) dltpw(1) = dltpw(1) + 1.0d0
        dltpw(2) = dltpw(2) + wght(ip) - wgh0(ip)
      enddo
      dltpw(3) = npmax

      tdlpw(1) = 0.0d0
      tdlpw(2) = 0.0d0
      tdlpw(3) = 0.0d0

      write(n6,'(2x,10x,6x,6x,3x,"  tot",i8,1pe12.2,3i8)')
     >   int(tdlpw(1)+eps), tdlpw(2),
     >   int(tdlpw(3)+eps), int(tdlpw(1)+eps),
     >   int(tdlpw(3)+eps) -int(tdlpw(1)+eps)

      return

!::dimension error
 920  continue
      write(n6,'(2x,"dimension error imp_cutbk  ii > ndmct",
     >  2x,2i6)') ii, ndmct
      call wexit("imp_cutbk","ii > ndmct")

 940  continue
      write(n6,'(2x,"dimension error imp_cutbk  iv > ndmsv",
     >  2x,2i6)') iv, ndmsv
      call wexit("imp_cutbk","iv > ndmsv")

 960  continue
      write(n6,'(2x,"dimension error imp_cutbk  iw > ndmsw",
     >  2x,2i6)') iw, ndmsw
      call wexit("imp_cutbk","iw > ndmsw")
      end
