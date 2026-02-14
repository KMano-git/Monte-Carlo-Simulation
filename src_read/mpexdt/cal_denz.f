! remade calculation of impurity density for each valence to TOPICS by kamata 2022/09/05
! base : IMPMC/caldnz.f
      subroutine cal_denz( ero, ednz, nd1, nd2 )
!---------------------------------------------------------------------
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
!    donot com_msdat
!      mpsp => ./sonic/inc/size_68/cgdcom
!      romp => ./soldor/cplmet
!
!      sub. mtrdisk  include 'com_gmdisk'  Note. 2015/08_30
!            csize, cplmet, cgdcom, cntcom, cntmnt, csonic, mpif.h
!
!     broh(5), bdnz(iz,5)  e : core edge
!
!---------------------------------------------------------------------
      use cimcom, only : ndis !, sflux, swtot
      use cmeffz, only : vdnz0
      use cntcom, only : mcel, volm
      use cplimp, only : ismaxl, wmc_nty
      use cplmet, only : icmpe, icmps, jcxp1, jcxp2
      use csize,  only : ndx, ndy
      use cunit,  only : n6
      use topics_mod, only : nddnz
      implicit none

! argument
      integer, intent(in)  :: nd1, nd2
      real(8), intent(in)  :: ero(nd1)
      real(8), intent(out) :: ednz(nd1,nd2)
! nd1  : number of rho near the core edge
! nd2  : array size of valences
! ero  : mesh rho near the core edge
! ednz : impurity density for each valence near the core edge

! local variables
      integer :: ic, is, ix, ixs, ixe, iy, iz, jy, jz
!     real(8) :: flnrm, sumd, sumv
      real(8) :: sumd, sumv

! check impurity density
! no impurity data
      if( wmc_nty < 1 ) then
        nddnz(3) = 0
        ednz(1:nd1,1:nd2) = 0.0_8
        return
      endif
      
! index check  ixs,ixe
      ixs = jcxp1 + 1
      ixe = jcxp2 - 1

!     flnrm = sflux / swtot
      jz = 0
      do is = 1, wmc_nty
        do iz = 1, ismaxl(is)
          jz = jz + 1
! averaged density ednz(iy,iz)
          jy = nd1
          do iy = icmps, icmpe
            sumd = 0.0_8
            sumv = 0.0_8
            do ix = ixs, ixe
              ic = mcel(ix,iy)
              sumd = sumd + vdnz0(ix,iy,iz,is) * volm(ic)
              sumv = sumv + volm(ic)
            enddo ! ix
            ednz(jy,jz) = sumd / sumv
            jy = jy - 1
          enddo   ! iy
        enddo     ! iz
      enddo     ! is
      nddnz(3) = jz

! output impurity density for each valence
      write(n6,'(a)') 'impurity density for each valence.'
      write(n6,'(2(a,i5))') 'nddnz(1) = ', nddnz(1)
     >    , 'nddnz(3) = ', nddnz(3)
      write(n6,'(a)') 'valence  rho  den-1  den-2  den-3 ...'
      do iy = 1, nddnz(1)
        write(n6,'(i5,f8.4,80es14.5)') iy, ero(iy)
     >      , ( ednz(iy,iz), iz = 1, nddnz(3) )
      enddo

      return
      end
