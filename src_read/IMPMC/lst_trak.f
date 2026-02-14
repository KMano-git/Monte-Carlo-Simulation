!***********************************************************************
      subroutine lst_trak(ip,comt,ptmx,ptim,pdtm)
!***********************************************************************
!)
!) How to use lst_trac  ==> lst_trak
!)
!) IMPHK/imshft.f:
!)   call lst_trak(ip,"imshft",ptmx,ptim,pdtm)
!)
!) IMPHK/imtrci.f:
!)   call lst_trak(ip,"ion_S",   ptmx, ptim, dtunt)
!)   call lst_trak(ip,"ion",     ptmx, ptim, dtunt)
!)   call lst_trak(ip,"ion_E",   ptmx, ptim, dtunt)
!)
!) IMPHK/imtrcn.f:
!)!  call lst_trak   comment out "return"
!)   call lst_trak(ip, "ntl_S",  ptmx, ptim, stb(ip))
!)   call lst_trak(ip, "n_hitw", ptmx, ptim, stb(ip))
!)   call lst_trak(ip, "n_psgt", ptmx, ptim, stb(ip))
!)   call lst_trak(ip, "n_wrfl", ptmx, ptim, stb(ip))
!)   call lst_trak(ip, "n_wstk", ptmx, ptim, stb(ip))
!)   call lst_trak(ip, "n_wabs", ptmx, ptim, stb(ip))
!)   call lst_trak(ip, "ntl_E",  ptmx, ptim, stb(ip))
!)
!)    comt : message
!)    ptmx : time to be traced untill t = ptmx
!)    ptim : current time
!)    pdtm : birth time of neutral particle
!)         : time step of ion   (dtunt=1.0d-6)
!)
!)      modified sub. impdt
!)
!)  ic-number
!)    ion O.K.
!)    neutal  error    ko = 1  !xx   if( ln.lt.0 ) wtm = wtm*0.9999d0
!)    neutral correct  ko = 0        if( ln > 0 )  wtm = wtm*1.001d0
!)
!)  selection
!)    kspcm: sput-type (CEpf,CEin,CEex) at imp_cal
!)    lpcm:  lt (dtout) at immont
!)    istcm: step at imstep
!)    ipcm:  ip   at imstep
!)
!-----------------------------------------------------------------------
      use cimcom, only : amz, ichcm, ien, i6_trak, il, ip_trak, ir, is
     >, istcm, kspcm, lpcm, rr, tt, wght, vv, vz, zz
      use cimntl, only : svx, svy, svz, sx0, sy0, sz0
      use cntcom, only : mrgn
      use cphcns, only : cev
      implicit none
!
!::argument
! modified 3/3 lines organize local variables and include files by kamata 2021/07/04
!ik   integer   :: ip
!ik   character :: comt*(*)
!ik   real(8)   :: ptmx, ptim, pdtm
      integer,   intent(in) :: ip
      character, intent(in) :: comt*(*)
      real(8),   intent(in) :: ptmx, ptim, pdtm
!
!::local variables
! modified 5/3 lines organize local variables and include files by kamata 2021/07/04
!ik   integer  ic, ko, lstp, imox, imoy, ln, ierr, iz, kwr, i6
!ik   real*8   tauz, dtstp
!ik   real*8   wtm, xi, yi, zi, ri, vlx, vly, vlz, zvl2, zvlb
!ik   real*8   funroh, zro
!ik   real*8   hdt, hftot, hdvz
      integer  ic, ko, ln, iz, kwr, i6
      real*8   wtm, xi, yi, zi, ri, zvl2, zvlb
      real*8   zro
!
!::reduce the number of points to plot   ds : displacement
      integer :: isbb = -1
      real(8) :: rbb = 0.0d0, zbb = 0.0d0
      real(8) :: dsr = 0.03d0, dsz = 0.03d0  ! <==
      save :: rbb, zbb, dsr, dsz
      save :: isbb
! added 3 lines organize local variables and include files by kamata 2021/07/04
! function
      integer    imox, imoy
      real(8)    funroh

      if( i6_trak <= 0 )  return
      if( ip /= ip_trak ) return

      i6 = i6_trak

!::condition
      if( comt.eq. "lnch" ) then
        write(i6,*)
      end if

      if( comt == "ionL00" ) then
        isbb = -1
        rbb  = 0.0d0
        zbb  = 0.0d0
      endif

!::ion
      if( is(ip) /= 0 ) then
      ic = ir(ip)
      ln = il(ip)  ! <== 2010/05/12
      iz = is(ip)
      call mchkin(rr(ip),zz(ip),ic,ko)
      zro = funroh(ic,rr(ip),zz(ip))

!::spin off (too many points for ion)
      kwr = 0
      if( iz /= isbb .or.
     >   dabs(rr(ip)-rbb) > dsr .or. dabs(zz(ip)-zbb) > dsz ) kwr = 1
      if( comt /= "ionLlt" ) kwr = 1

      if( kwr == 1 ) then
! modified 1/1 lines with TOPICS by kamata 2021/12/22
!ik   write(i6,'(1x,a6,i5,i8,i8,i8,i6,1pe10.2,1p2e12.4,1pe10.2,
      write(i6,'(1x,a6,i5,i10,i8,i8,i6,1pe10.2,1p2e12.4,1pe10.2,
     >  i6,2i4,i4,i6,3i4, 0p2f10.5, 0p2f8.4, 1p2e11.2)')
     >  comt, kspcm, lpcm, ichcm, istcm, ip, ptmx, ptim, tt(ip), pdtm,
     >  ic, imox(ic), imoy(ic), mrgn(ic), ln, ien(ip), ko, is(ip),
     >  rr(ip), zz(ip), zro, wght(ip), vz(ip), 0.5d0*amz*vv(ip)/cev
      rbb = rr(ip)
      zbb = zz(ip)
      isbb = iz
      endif
      call flush(i6)
!
!::neutral
      else
      ic = ir(ip)
      ln = il(ip)
      wtm = ptim - pdtm
!xx   if( ln.lt.0 ) wtm = wtm*0.9999d0
      if( ln > 0 )  wtm = wtm*1.001d0
      xi = sx0(ip) + wtm*svx(ip)
      yi = sy0(ip) + wtm*svy(ip)
      zi = sz0(ip) + wtm*svz(ip)
      ri = sqrt(xi**2+yi**2)
      call mchkin(ri,zi,ic,ko)
      zvl2 = svx(ip)**2+svy(ip)**2+svz(ip)**2

!::no need zvlb  2015/09/04
!xx   call velion(xi,yi,zi,svx(ip),svy(ip),svz(ip),zvlb,ierr)
      zvlb = 0.0d0  ! avoid to occur errors
      zro = 99.0d0
      if( ic > 0 ) zro = funroh(ic,ri,zi)
! modified 1/1 lines with TOPICS by kamata 2021/12/22
!ik   write(i6,'(1x,a6,i5,i8,i8,i8,i6,1pe10.2,1p2e12.4,1pe10.2,
      write(i6,'(1x,a6,i5,i10,i8,i8,i6,1pe10.2,1p2e12.4,1pe10.2,
     >  i6,2i4,i4,i6,3i4, 0p2f10.5, 0p2f8.4, 1p2e11.2)')
     >  comt, kspcm, lpcm, ichcm, istcm, ip, ptmx, ptim, tt(ip), pdtm,
     >  ic, imox(ic), imoy(ic), mrgn(ic), ln, ien(ip), ko, is(ip),
     >  ri, zi, zro, wght(ip), zvlb, 0.5d0*amz*zvl2/cev
      call flush(i6)
      endif

      return
      end


!***********************************************************************
      subroutine lst_trak_open(ipnm)
!***********************************************************************
      use cimcom, only : i6_trak, ip_trak
      use cunit,  only : lmype, ufl6
      implicit none

! modified 1/1 lines organize local variables and include files by kamata 2021/07/04
!ik   integer :: ipnm
      integer, intent(in) :: ipnm

      character :: dsn*(30)
      integer :: i6

      i6_trak = 0
      ip_trak = 0
      dsn = "  "
      if( ufl6 == "/dev/null" ) return

      i6_trak = 161
      ip_trak = ipnm

      i6 = i6_trak
      write(dsn,'(a,i3.3)') "lst_phst_",lmype
      open(unit=i6,file=dsn)

      write(i6,'(1x, 2x,a, 2x,a, 3x,a, 2x,a, 5x,a, 4x,a,
     >   6x,a, 8x,a, 10x,a, 6x,a,
     >   4x,a, 2x,a, 2x,a, 1x,a, 4x,a, 1x,a, 2x,a, 2x,a,
     >   8x,a, 8x,a, 6x,a, 4x,a, 9x,a, 8x,a)')
     >   "comt", "ksp", "lp", "ich", "ist", "ip",
     >       "ptmx", "ptim", "tt", "pdtm",
     >   "ic", "ix", "iy", "reg","ln", "ien", "ko", "is",
     >   "rr", "zz", "ro", "wght", "vz", "Eng"

      return
      end

!***********************************************************************
      subroutine lst_trak_head(ip,comt)
!***********************************************************************
      use cimcom, only : i6_trak, ip_trak
      implicit none

! modified 2/2 lines organize local variables and include files by kamata 2021/07/04
!ik   integer :: ip
!ik   character :: comt*(*)
      integer,   intent(in) :: ip
      character, intent(in) :: comt*(*)

      integer :: i6

      if( i6_trak <= 0 ) return
      if( ip /= ip_trak ) return

      i6 = i6_trak
      write(i6,'(1x,a6,4x,i6)') comt, ip
      return
      end
