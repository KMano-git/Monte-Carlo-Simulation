!*********************************************************************
      subroutine lst_scat(tnow,mskp,i6)
!*********************************************************************
!)
!)     lst_scat : list for scatter plot
!)
!)  How to use
!)    if( i6_scat <= 0 ) call lst_open("scat")
!)    call lst_scat(tnow,mskp)   mskp = 1 (all) = 500 (eco)
!)
!)  Note.  neutral particle information
!)       birth time : stb(ip), (sx0,sy0,sz0) (svx,svy,svz)
!)       now   time : tt(ip)  stm(ip) from stb(ip)  tt = stb + stm
!)                    x0=sx0(ip)+stm(ip)*svx(ip)
!)                    y0=sy0(ip)+stm(ip)*svy(ip)
!)
!--------------------------------------------------------------------
      use cimcom, only : amz, ien, il, ir, is, npmax, rr, tt, vv, vz
     >    , wght, wstp, zz
      use cimntl, only : stm, svx, svy, svz, sx0, sy0, sz0
      use cntcom, only : mrgn
      use cphcns, only : cev
      implicit none

! modified 2/2 lines organize local variables and include files by kamata 2021/07/04
!ik   real(8) :: tnow
!ik   integer :: mskp, i6
      real(8), intent(in) :: tnow
      integer, intent(in) :: mskp, i6

! modified 1/1 lines organize local variables and include files by kamata 2021/07/04
!ik   integer :: ip, iz, ic, imox, imoy, ko
      integer :: ip, ic, ko
      real(8) :: ztm, x1, y1, z1, r1
      real(8) :: hrr, hzz, hv2, hvl, hvz, hez
      integer :: iskp
! added 2 lines organize local variables and include files by kamata 2021/07/04
! function
      integer    imox, imoy

      if( i6 <= 0 ) return

      iskp = mskp
      if( iskp <= 0 ) iskp = npmax

!::header of lst_scat
      write(i6,'(/2x,"particle information  tnow =",1pe12.4,
     >  "  npmax =",i6,"  iskp =",i6)') tnow, npmax, iskp
      write(i6,'(4x,"ip",2x,"tt",11x,"ic",6x,"ix",3x,"iy",2x,"mrg",2x,
     >  "ien",3x,"ko",3x,"il",3x,"is",4x,"rr",10x,"zz",9x,"wgt",6x,
     >  "wst",6x,"Vz",10x,"Vl",10x,"Ez")')

      do ip = 1, npmax, iskp

!::neutral
        if( is(ip) == 0 ) then
          ztm = stm(ip)
          ztm = 0.9999d0*ztm  ! wall point => interior point
          x1  = sx0(ip) + ztm*svx(ip)
          y1  = sy0(ip) + ztm*svy(ip)
          z1  = sz0(ip) + ztm*svz(ip)
          r1  = sqrt(x1**2+y1**2)

          hrr = r1
          hzz = z1
          hv2 = svx(ip)**2+svy(ip)**2+svz(ip)**2
          hvl = sqrt(hv2)
          hvz = hvl
          hez = 0.5d0*amz*hv2/cev
!::ion
        else
          hrr = rr(ip)
          hzz = zz(ip)
          hvl = sqrt(vv(ip))
          hvz = vz(ip)
          hez = 0.5d0*amz*vv(ip)/cev
        endif

      ic = ir(ip)
      call mchkin(hrr,hzz,ic,ko)
      write(i6,'(1x,i5,1pe12.4,i7,i6,i5,5i5,0p2f12.6,0p2f9.4,
     >  1p3e12.3)') ip,tt(ip),ic,imox(ic),imoy(ic),mrgn(ic), ien(ip),
     >  ko,il(ip),is(ip),hrr,hzz,wght(ip), wstp(ip), hvz, hvl, hez
      enddo

      return
      end
