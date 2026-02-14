!**********************************************************************
      subroutine listmc
!**********************************************************************
!
!      file    "MC_parm.txt"   whole region
!      region   Odiv + Oprv
!               SOL  + Main
!               Idiv + Iprv
!               Vac
!
!     ir ==> jr  because of concatinate ir(ndmp) in cplcom
!
!   cell name in Monte Carlo    gmdisk/readMe_mesh
!   include                     soldor/cplmet
!
!   region number  mrgn(ic)
!  -------------------------------
!     Odiv 1  SOL 2  Idiv 3  Oprv 4  Iprv 5  edge 6  main 7
!     Vacume (sub-divertor) 8
!
!
!                       +icaxs--------------------+
!                       !                   (7)   !
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
!
!---------------------------------------------------------------------
      use csize
      use cimcom
      use cimden
      use cntpls
      use cntsrc
      use cplcom
      use cplimp
      use cplmet
      use cplwrd
      use cunit, only : n6, wh_dir
      use cntcom ! mknd, ihwl, chwl
      use com_impvar, only : zni, zne, znz, zav
     > , zef, zrt, znzi, znzt
      use cplimp_plimp, only : twrdL
      implicit none

!::local variables
      integer :: i, ic, jr, k, ip, nft, ix, iy, j, jmax
      integer :: ia, iz, nty
      integer :: ip1, ip2, ip3, ip4, kmax
      integer :: ixmn(8), ixmx(8), iymn(8), iymx(8)
      integer :: kxmn(3), kxmx(3), kymn(3), kymx(3)
      real(8) :: xc, yc
      real(8) :: hni(ndsp), hnz(0:ndmis,1:nzsmx)
      real(8) :: p0, pg

      character(3) :: creg(8) = (/"Odv", "SOL", "Idv", "Opv", "Ipv",
     >    "Edg", "Cor", "Vac"/)
      character(9) :: ctyp(3)
      integer :: i6 = 180000 + 21
      integer :: iys, iye

      ! add wall object information in MC_parm_vac.txt
      ! iw_size must equal or greater than the max value of mseg
      integer, parameter :: iw_size = 4
      character(len=4) :: wall_obj(iw_size) ! wall information
      integer :: iw(iw_size) ! wall cell index

      character(80) :: cdsn
! function
      integer lenx
      integer imox, imoy ! gmdisk/ntwlst.f

!::wfac (normalization factor in soldor)
      write(n6,'(/2x,"***listmc ***  wfac =",1pe12.4)') wfac
      write(i6,'(/2x,"***listmc *** enh wfac=",1pe12.4)') wfac
!KH      if( wfac == 0.0d0 ) stop

!::file
      nft = 101
! cdsn => "../wxdr_XX/MC_parm.txt"
      cdsn = wh_dir(1:lenx(wh_dir)) // "/MC_parm.txt"
      open(unit=nft, file=cdsn)

!::index min and max
      ixmn(1:8) =  1000
      ixmx(1:8) = -1000
      iymn(1:8) =  1000
      iymx(1:8) = -1000

      do ic = 1, ncmax2
        jr = mrgn(ic)
        ix = imox(ic)
        iy = imoy(ic)
        if( jr == 0 ) cycle
        if( ix < ixmn(jr) ) ixmn(jr) = ix
        if( ix > ixmx(jr) ) ixmx(jr) = ix
        if( iy < iymn(jr) ) iymn(jr) = iy
        if( iy > iymx(jr) ) iymx(jr) = iy
      enddo
      do jr = 1, 8
        write(nft,'(a,3x,a,2x,4i5)') 
     >   "*", creg(jr), ixmn(jr), ixmx(jr), iymn(jr), iymx(jr)
      enddo

!::OdvOpv SOLEdgCor, IdvIpv
      ctyp(1) = creg(1)//creg(4)
      ctyp(2) = creg(2)//creg(6)//creg(7)
      ctyp(3) = creg(3)//creg(5)
      kxmn(1) = min0( ixmn(1), ixmn(4) )
      kxmx(1) = max0( ixmx(1), ixmx(4) )
      kymn(1) = min0( iymn(1), iymn(4) )
      kymx(1) = max0( iymx(1), iymx(4) )

      kxmn(2) = min0( ixmn(2), ixmn(6), ixmn(7) )
      kxmx(2) = max0( ixmx(2), ixmx(6), ixmx(7) )
      kymn(2) = min0( iymn(2), iymn(6), iymn(7) )
      kymx(2) = max0( iymx(2), iymx(6), iymx(7) )

      kxmn(3) = min0( ixmn(3), ixmn(5) )
      kxmx(3) = max0( ixmx(3), ixmx(5) )
      kymn(3) = min0( iymn(3), iymn(5) )
      kymx(3) = max0( iymx(3), iymx(5) )

      write(nft,'(a4,1x,a5,4x,a,1x,a,25(2x,a3,6x),3(2x,a5,6x))') 
     >  "* rg","ic ","ix ","iy ",
     >  "xc ","yc ",
     >  "vol", "Zav","Zef", "Nzi", "Nzt", "zne",
     >  "Ni ", "Ti ", "Te ", "N0 ", "E0 ", "P0",
     >  "Ng ", "Eg ", "Pg",
     >  "xg1", "yg1", "xg2", "yg2", "xg3", "yg3", "xg4", "yg4",
     >  "nfl0x", "nfl0y", "nfl0z"

!::exclude sub_divertor
      do ic = 1, ncmax
        jr = mrgn(ic)

!::xc, yc
        ip1 = mgrd(ic,1)
        ip2 = mgrd(ic,2)
        ip3 = mgrd(ic,3)
        ip4 = mgrd(ic,4)
        jmax = mseg(ic)
        xc = 0.0d0
        yc = 0.0d0
        do j = 1, jmax
          ip = mgrd(ic,j)
          xc = xc + xpnt(ip)
          yc = yc + ypnt(ip)
        enddo
        xc = xc/dfloat(jmax)
        yc = yc/dfloat(jmax)

!::<Z>, Zeff
!
        hni(1:nion) = deni(ic,1:nion)
        do nty = 1, wmc_nty
          hnz(0:ismaxL(nty),nty) = tdnzL(0:ismaxL(nty),ic,nty)
        enddo

        zni = 0.0d0
        zne = 0.0d0
        zef = 0.0d0
        do ia = 1, nion
          zni = zni + hni(ia)
          zne = zne + aza(ia)*hni(ia)
          zef = zef + aza(ia)**2*hni(ia)
        enddo

        znz = 0.0d0
        znzi = 0.0d0
        zav = 0.0d0
        do nty = 1, wmc_nty
          do iz = 0, ismaxL(nty)
            zne = zne + dble(iz)*hnz(iz,nty)
            zef = zef + dble(iz)**2*hnz(iz,nty)
            znz = znz + hnz(iz,nty)
            zav = zav + dble(iz)*hnz(iz,nty)
            if( iz > 0 ) znzi = znzi + hnz(iz,nty)
          enddo ! iz
        enddo ! nty
        znzt = znz
      
        if( znz > 0.0d0 ) then
          zav = zav/znz
        else
          zav = 0.0d0
        endif
        if( zne > 0.0d0 ) then
          zef = zef/zne
          zrt = znz/zne
        else
          zef = 0.0d0
          zrt = 0.0d0
        endif

!::no plasma but impurity exists at void
        if( zni == 0.0d0 ) then
          zef = 0.0d0
          zav = 0.0d0
          zne = 0.0d0
        endif

        p0 = tden0(ic,1)*teng0(ic,1)*1.60210d-19
        pg = tdeng(ic,1)*tengg(ic,1)*1.60210d-19

        write(nft,'(i4,i7,i5,i4,0p2f11.6,1p15e11.3,0p8f11.6,1p3e12.3)')
     >     jr, ic, imox(ic), imoy(ic),
     >     xc, yc, 
     >     volm(ic), zav, zef, znzi, znzt, zne,
     >     deni(ic,1), temi(ic), teme(ic),
     >     tden0(ic,1), teng0(ic,1), p0,
     >     tdeng(ic,1), tengg(ic,1), pg,
     >     xpnt(ip1), ypnt(ip1), xpnt(ip2), ypnt(ip2),
     >     xpnt(ip3), ypnt(ip3), xpnt(ip4), ypnt(ip4),
     >     nfl0x(ic,1), nfl0y(ic,1), nfl0z(ic,1)
      enddo !ic
      close(nft)

!::debug
      write(i6,'(5x,"ic",3x,"ix",3x,"iy",3x,"twrd",5x,"tdnz",
     >     10(i2,9x))') (iz,iz=0,12,2)
      ix = jcxp1 + 1
      iys = icspx - 5
      iye = icspx + 5
      do iy = iys, iye, 3
        ic = mcel(ix,iy)
!-------
        hni(1:nion) = deni(ic,1:nion)
        if(wmc_nty > 0) then
        hnz(0:ismaxl(1),1) = tdnzL(0:ismaxL(1),ic,1) !*wfac
        do nty = 2, wmc_nty
          hnz(0:ismaxl(nty),nty) = tdnzL(0:ismaxL(nty),ic,nty)
        enddo
        endif
!-------
        do nty = 1, wmc_nty
          write(i6,'(i2, i7,i5,i3,1pe11.3,1p10e11.3)')
     >     nty, ic, imox(ic), imoy(ic),
     >     twrdL(ic,nty), (hnz(iz,nty),iz=0,12,2)
        enddo
      enddo !iy

!KH191126
!vac region
! cdsn => "../wxdr_XX/MC_parm_vac.txt"
      cdsn = wh_dir(1:lenx(wh_dir)) // "/MC_parm_vac.txt"
      open(unit=nft, file=cdsn)
      write(nft,'(a6,4x,a,1x,a,17(2x,a3,6x),3(2x,a5,5x),
     >  4(2x,a4,4x,a3))')
     >  "*   ic","ix ","iy ",
     >  "xc ","yc ",
     >  "vol", 
     >  "N0 ", "E0 ", "P0 ", "Ng ", "Eg ","Pg ",
     >  "xg1", "yg1", "xg2", "yg2", "xg3", "yg3", "xg4", "yg4",
     >  "nfl0x", "nfl0y", "nfl0z",
     >  "obj1", "iw1", "obj2", "iw2","obj3", "iw3","obj4", "iw4"
!
      do ic = ncmax+1, ncmax2
!::xc, yc
        ip1 = mgrd(ic,1)
        ip2 = mgrd(ic,2)
        ip3 = mgrd(ic,3)
        ip4 = mgrd(ic,4)
        jmax = mseg(ic)
        xc = 0.0d0
        yc = 0.0d0
        do j = 1, jmax
          ip = mgrd(ic,j)
          xc = xc + xpnt(ip)
          yc = yc + ypnt(ip)
        enddo
        xc = xc/dfloat(jmax)
        yc = yc/dfloat(jmax)

        p0 = tden0(ic,1)*teng0(ic,1)*1.60210d-19
        pg = tdeng(ic,1)*tengg(ic,1)*1.60210d-19

        if(imox(ic).ne.1 .and. imoy(ic).eq.1) write(nft,*)

        call get_wall_obj(ic,wall_obj,iw,iw_size) ! set wall_obj and iw
        write(nft,'(i8,i5,i4,0p2f11.6,1p7e11.3,0p8f11.6,1p3e12.3,
     >   4(3x,a4,2x,i4))')
     >     ic, imox(ic), imoy(ic),
     >     xc, yc, 
     >     volm(ic), 
     >     tden0(ic,1), teng0(ic,1), p0,
     >     tdeng(ic,1), tengg(ic,1), pg,
     >     xpnt(ip1), ypnt(ip1), xpnt(ip2), ypnt(ip2),
     >     xpnt(ip3), ypnt(ip3), xpnt(ip4), ypnt(ip4),
     >     nfl0x(ic,1), nfl0y(ic,1), nfl0z(ic,1),
     >     (wall_obj(k), iw(k), k=1,iw_size)
      enddo
      
      close(nft)
      return

      contains
!####################################################
      subroutine get_wall_obj(ic_in,wall,iw_out,iw_size_in)
!####################################################
!     Return the wall-index(iw) and the wall-information(wall_obj), 
!       which is next to ic(fluid mesh index).
!     If ic is not next to wall, iw = 0 and wall_obj = "#   "
!####################################################
      implicit none
      ! arguments
      integer, intent(in) :: ic_in
      integer, intent(in) :: iw_size_in
      integer, intent(out) :: iw_out(iw_size_in) ! wall cell index
      character(len=4), intent(out) :: wall(iw_size)  ! wall information

      ! local variables
      integer :: k = 0, ln = 1

      ! check iw_size is large enough
      if(mseg(ic_in) .gt. iw_size_in) then
        call wexit("listmc","too small iw_size. Check the mseg value")
      endif

      ! initialization
      do k = 1, iw_size_in
        iw_out(k) = 0
        wall(k) = "#   "
      enddo

      ! if ic is next to wall, ln < 0
      do k = 1, mseg(ic_in) ! kmax = mseg(ic_in)
        ln = mknd(ic_in,k)
        if(ln .lt. 0) then
          iw_out(k) = -ln
          wall(k) = chwl(ihwl(iw_out(k)))
        endif
      enddo

      end subroutine get_wall_obj

      end subroutine listmc