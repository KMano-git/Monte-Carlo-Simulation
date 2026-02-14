!*********************************************************************
      subroutine caldne(eroh,edne,etne,iore)
!*********************************************************************
!
!          see  gmdisk/readMe_mesh
!
!      icmps                         icaxs    cell center  i-1/2
!     |     |     |     |     |     |     |
!     *--x--*--x--*--x--*--x--*--x--*--x--*   cell surface I
!     |       A      A                    |
!    romn(mpsp)      |                 romn(mpax)
!             |    croh(iy)
!             |    cdnz(iy)  iy = imps, .. icaxs   c : core
!             |
!          fls_roh(1:5)
!          fls_dnz(1:5)  interpolation from (croh,cdnz)
!
!    cdnz(iy) = sum(j)[denz(ic)*volm(ic]/ sum(j)[volm(ic)]
!           j: poloidal mesh number  ic:cell number
!        1D variables (MONTE) ic = 2D variables (SOLDOR)(j,iy)
!
!
!    donot com_msdat
!      mpsp => ./sonic/inc/size_68/cgdcom
!      romp => ./soldor/cplmet
!
!      sub. mtrdisk  include 'com_gmdisk'  Note. 2015/08_30
!            csize, cplmet, cgdcom, cntcom, cntmnt, csonic, mpif.h
!
!     broh(5), bdnz(iz,5)  e : core edge
!
!     iore = 1 ion
!     iore = 2 electron
!---------------------------------------------------------------------
      use cgdcom, only : mpsp
      use cimcom, only : sflux, sptyc, swtot
      use cntcom, only : mcel, ncmax2, volm
      use cntpls, only : dene, deni, teme, temi
      use cplmet, only : icaxs, icmpe, icmps, jcxp1, jcxp2, romn
      use csize,  only : ndmc, ndy
      use cunit,  only : n6
      implicit none

!::argument
      real(8), intent(in)  :: eroh(5)
      real(8), intent(out) :: edne(5), etne(5)
      integer, intent(in)  :: iore

!::local variables
      integer :: ic, ix, iy, mst, men, ii, n
      integer :: ixs, ixe
      real(8) :: znz, sumd, sumv, ros, roe, roh, ztemp, sumt
      real(8) :: hro, hne, hte, wro(ndy), wne(ndy), wte(ndy), croh(ndy)

      real(8) :: dne(ndmc), tne(ndmc)
      real(8) :: cdne(ndy), ctne(ndy)
! function
      real(8) :: dintp

!::index check  ixs,ixe
      ixs = jcxp1 + 1
      ixe = jcxp2 - 1
      mst = mpsp
      do iy = icmps, icaxs
        men = mst + 1
        ros = romn(mst)
        roe = romn(men)
        roh = 0.5d0*(ros+roe)
        croh(iy) = roh
        mst = men
      enddo

!::defince dnz(iz,ic)  z:charge state ic:cell number
      write(n6,'(2x,"sptyc = ",a,"  sflux =",1pe12.4,
     &      "  swtot =",1pe12.4)') trim(sptyc), sflux, swtot
      do ic = 1, ncmax2
        if(iore.eq.1) then
          znz = deni(ic,1) ! assuming single fluid
          ztemp = temi(ic)
        else
          znz = dene(ic)
          ztemp = teme(ic)
        endif
        dne(ic) = znz
        tne(ic) = ztemp
      enddo

!::averaged density cdnz(iy)
      do iy = icmps, icaxs
        sumd = 0.0d0
        sumt = 0.0d0
        sumv = 0.0d0
        do ix = ixs, ixe
          ic = mcel(ix,iy)
          sumd = sumd + dne(ic)*volm(ic)
          sumt = sumt + tne(ic)*volm(ic)
          sumv = sumv + volm(ic)
        enddo  ! loop(ix)
        cdne(iy) = sumd/sumv
        ctne(iy) = sumt/sumv
      enddo  ! loop(iy)

!:: averaged density at sep SOL
      iy = 30
      sumd = 0.0d0
      sumt = 0.0d0
      sumv = 0.0d0
      do ix = ixs, ixe
        ic = mcel(ix,iy)
        sumd = sumd + dne(ic)*volm(ic)
        sumt = sumt + tne(ic)*volm(ic)
        sumv = sumv + volm(ic)
      enddo  ! loop(ix)
      cdne(iy) = sumd/sumv
      ctne(iy) = sumt/sumv

      if(iore.eq.1) write(n6,*) 'CALDN <ni>,sep =', cdne(30)
      if(iore.eq.2) write(n6,*) 'CALDN <ne>,sep =', cdne(30)

!::interpolation (wro, wnz) => hro, hdnz
      ii = 0
      do iy = icmpe, icmps, -1   ! in the big order
        ii = ii + 1
        wro(ii) = croh(iy)
        wne(ii) = cdne(iy)
        wte(ii) = ctne(iy)
      enddo
      do n = 1, 5
        hro = eroh(n)
        hne = dintp(hro,ii,wro,wne)
        hte = dintp(hro,ii,wro,wte)
        edne(n) = hne
        etne(n) = hte
      enddo

      return
      end
