!*********************************************************************
      subroutine caldnz(eroh,ednz)
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
!---------------------------------------------------------------------
      use cgdcom, only : mpsp
      use cimcom, only : denz, ismax, ndis, phdnz, sflux, sptyc, swtot
      use cntcom, only : mcel, ncmax2, volm
      use cplmet, only : icaxs, icmpe, icmps, jcxp1, jcxp2, romn
      use csize,  only : ndmc, ndy
      use cunit,  only : n6
      use mod_shexe, only : impmc_model
      implicit none

!::argument
! modified 1/2 lines organize local variables and include files by kamata 2021/06/28
!ik   real(8) :: eroh(5), ednz(0:ndis,5)
      real(8), intent(in)  :: eroh(5)
      real(8), intent(out) :: ednz(0:ndis,5)

!::local variables
! modified 2/2 lines organize local variables and include files by kamata 2021/06/28
!ik   integer :: nv, ic, iz, ix, iy, mst, men, ii, n, nw
!ik   integer :: ixs, ixe, it, jws, jwe, ixa, ixb, iya, iyb
      integer :: ic, iz, ix, iy, mst, men, ii, n
      integer :: ixs, ixe
      real(8) :: flnrm, znz, sumd, sumv, ros, roe, roh
      real(8) :: hro, hnz, wro(ndy), wnz(ndy), croh(ndy)
      integer :: i6

      real(8) :: dnz(0:ndis,ndmc)
! modified 1/2 lines organize local variables and include files by kamata 2021/06/28
!ik   real(8) :: cnz(0:ndis,ndmc), cdnz(0:ndis,ndy)
      real(8) :: cdnz(0:ndis,ndy)
! function
      real(8) :: dintp

      i6 = n6
!xx   write(i6,'(2x,"*** caldnz ***  itim =",i6)') itim

!::index check  ixs,ixe
      ixs = jcxp1 + 1
      ixe = jcxp2 - 1
! deleted 7 lines organize local variables and include files by kamata 2021/06/28
!ik   it  = itmps      ! flux tube
!ik   jws = jtmin(it)
!ik   jwe = jtmax(it)
!ik   ixa = jcel(jws+1,it)  ! dummy cell in core regin
!ik   ixb = jcel(jwe-1,it)
!ik   iya = icel(jws+1,it)
!ik   iyb = icel(jwe-1,it)

!x      write(i6,'(2x,"<pol> : ixs/ixa =",2i3,"  ixe/iyb =",2i3,
!x     >  "  iya/oyb =",2i3,"  icmps/icmpe =",2i3)')
!x     >  ixs, ixa, ixe, ixb, iya, iyb, icmps, icmpe

!::define cro(iy) and croh(iy)
!x      write(i6,'(/2x,"define croh(iy) iy=icmps,icaxs")')
!x      write(i6,'(4x,"iy",3x,"mst",3x,"men",3x,"ros",9x,"roe",
!x     >   9x,"roh")')
        mst = mpsp
        do iy = icmps, icaxs
          men = mst + 1
          ros = romn(mst)
          roe = romn(men)
          roh = 0.5d0*(ros+roe)
!x        write(i6,'(2x,3i5,3f12.4)') iy, mst, men, ros, roe, roh
          croh(iy) = roh
          mst = men
        enddo

!::defince dnz(iz,ic)  z:charge state ic:cell number
      write(i6,'(2x,"sptyc = ",a,"  sflux =",1pe12.4,
     &      "  swtot =",1pe12.4)') trim(sptyc), sflux, swtot
! added 1 line integration of ST and TD versionis of IMPMC by kamata 2022/05/14
      if( impmc_model == 0 ) then
! deleted 1 line organize local variables and include files by kamata 2021/06/28
!ik   nv = ncmax2
      flnrm = sflux/swtot
! modified 1/1 lines organize local variables and include files by kamata 2021/06/28
!ik   do ic = 1, nv
      do ic = 1, ncmax2
      do iz = 0, ismax
      znz = flnrm*denz(iz,ic)/volm(ic)
      dnz(iz,ic) = znz
      enddo
      enddo
! added 1 line integration of ST and TD versionis of IMPMC by kamata 2022/05/14
      endif

!::averaged density cdnz(iy)
      do iz = 0, ismax
      do iy = icmps, icaxs
        sumd = 0.0d0
        sumv = 0.0d0
        do ix = ixs, ixe
          ic = mcel(ix,iy)
! modified 1/2 lines integration of ST and TD versionis of IMPMC by kamata 2022/05/14
!ik       sumd = sumd + dnz(iz,ic)*volm(ic)
          sumd = sumd + ( impmc_model * phdnz(iz,ic )
     >         + ( 1 - impmc_model ) * dnz(iz,ic) ) * volm(ic)
          sumv = sumv + volm(ic)
        enddo  ! loop(ix)
        cdnz(iz,iy) = sumd/sumv
      enddo  ! loop(iy)
      enddo  ! loop(iz)

!::debug write
!x      write(i6,'(/2x,"cdnz(iz,iy)")')
!x      write(i6,'(3x,"i",3x,"iy",3x,10(i3,"+",7x:))')
!x     >  (iz,iz=0,ismax,3)
!x      ii = 0
!x      do iy = icmps, icaxs
!x        ii = ii + 1
!x        write(i6,'(2x,i3,2x,i3,2x,1p8e11.3)')
!x     >    ii, iy, (cdnz(iz,iy),iz=0,ismax,3)
!x      sumd = 0.0d0
!x      do iz = 0, ismax
!x        sumd = sumd + cdnz(iz,iy)
!x      enddo
!x      if( sumd <= 0.0d0 ) exit
!x      enddo

!::interpolation (wro, wnz) => hro, hdnz
      do iz = 0, ismax
        ii = 0
        do iy = icmpe, icmps, -1   ! in the big order
          ii = ii + 1
          wro(ii) = croh(iy)
          wnz(ii) = cdnz(iz,iy)
        enddo
! deleted 1 line organize local variables and include files by kamata 2021/06/28
!ik     nw = ii
        do n = 1, 5
          hro = eroh(n)
! modified 1/1 lines organize local variables and include files by kamata 2021/06/28
!ik       hnz = dintp(hro,nw,wro,wnz)
          hnz = dintp(hro,ii,wro,wnz)
          ednz(iz,n) = hnz
!---------
!x          write(i6,'(/2x,"iz =",i3,"  n =",i3,"  nw =",i3)') iz, n, nw
!x          write(i6,'(2x,"wro =",1p10e11.3)') wro(1:nw)
!x          write(i6,'(2x,"wnz =",1p10e11.3)') wnz(1:nw)
!x          write(i6,'(2x,"hro =",1pe11.3,"  hnz =",1pe11.3)') hro, hnz
!---------
        enddo
      enddo

      return
      end
