!***********************************************************************
      subroutine imgpont(cspt)
!***********************************************************************
      use cimcom, only : amz, sflux, tfbz
      use cimptl, only : ic_00, is_00, ri_00, vv_00, vz_00, zi_00
      use cntcom, only : mgrd, mseg, ncmax, xpnt, ypnt
      use cphcns, only : cev
      use cunit,  only : n6
      implicit none
!
!::argument
! modified 1/1 lines organize local variables and include files by kamata 2021/06/28
!ik   character  cspt*(*)
      character, intent(in) :: cspt*(*) ! dummy
!
!::local variables
      real*8   ri, zi, engi, cosi, xc, yc
! modified 1/3 lines organize local variables and include files by kamata 2021/06/28
!ik   integer  n, ic, ko, k, imox, imoy
      integer  n, ic, ko, k
! function
      integer    imox, imoy
!
!::input data
      ri  = 3.5300d0
      zi  = 3.5300d0
      engi = 7.0d3
      cosi = 0.85d0
!
      write(n6,'(/2x,"*** imgpont ***")')
      write(n6,'(2x,"ri,zi =",2f9.5,"  engi,cosi =",1p2e11.3)')
     >  ri, zi, engi, cosi
!
      do n = 1, ncmax
      ic = n
      call mchkin(ri,zi,ic,ko)
      if( ko.eq.0 ) goto 110
      enddo
      call wexit("imgpnt","no found cell number")
!
 110  continue
      write(n6,'(2x,"rg,zg =",2f9.5)')
     >   (xpnt(mgrd(ic,k)),ypnt(mgrd(ic,k)),k=1,mseg(ic))
!
      xc = 0.0d0
      yc = 0.0d0
      do k = 1, mseg(ic)
      xc = xc + xpnt(mgrd(ic,k))
      yc = yc + ypnt(mgrd(ic,k))
      enddo
      xc = xc/dfloat(mseg(ic))
      yc = yc/dfloat(mseg(ic))
      write(n6,'(2x,"rc,zc =",2f9.5)') xc, yc
!
      ic_00 = ic
      is_00 = 16
      ri_00 = xc
      zi_00 = yc
      vv_00 = 2.0d0*engi*cev/amz
      vz_00 = sqrt(vv_00)*cosi
!
      write(n6,'(2x,"ri,zi =",2f9.5,"  v,vz =",1p2e11.3)')
     >   ri_00, zi_00, sqrt(vv_00), vz_00
      write(n6,'(2x,"ic =",i6,"  ix =",i4,"  iy =",i3)')
     >   ic_00, imox(ic_00), imoy(ic_00)
!
!::flux
      tfbz = 4.5d21
      sflux = tfbz
      write(n6,'(2x,"tfbz =",1pe11.3)') tfbz
!
      return
      end
