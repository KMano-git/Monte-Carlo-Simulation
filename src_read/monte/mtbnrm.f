!**********************************************************************
      subroutine mtbnrm
!**********************************************************************
      use cntcom, only : fnor, fnor2, ndnrm, ndnrm2, tnor, tnor2
      use cunit,  only : n6
      implicit none
!
!::local variables
      integer :: ndim
!
      write(n6,'(2x)')
      write(n6,'(2x,"*** mtbnrm *** for CX process")')
      ndim = ndnrm
      call set_norm(ndim,fnor,tnor)
!
      write(n6,'(2x,"*** mtbnrm *** for EL process")')
      ndim = ndnrm2
      call set_norm(ndim,fnor2,tnor2)
!
      return
      end
!
!**********************************************************************
      subroutine set_norm(ndim,fnrm,tnrm)
!**********************************************************************
      use cunit, only: n6
      implicit none
!
!::argument
      integer, intent(in)  :: ndim
      real(8), intent(out) :: fnrm, tnrm(ndim)
!
!::local variables
      real(8) :: cut
      integer :: imax, i2, i, nnrm
      real(8) :: dx, x, y, z
      real(8) :: zsg, zfc
!
      cut  = 0.49977d0
!
      nnrm = ndim - 1
      fnrm = dfloat(nnrm)
!
      imax = nnrm/2
      nnrm = imax*2
      i2   = nnrm + 1
      dx   = cut/dfloat(imax)
      do  i = 1, imax
        i2 = i2 - 1
        x = dx*dfloat(i)
        y = sqrt(2.0d0*dlog(1.0d0/(0.5d0-x)))
        z = y-(2.30753d0+0.27061d0*y)/(1.0d0+y*(0.99229d0+0.04481d0*y))
        tnrm(i)  =  z
        tnrm(i2) = -z
      enddo
!
!::ave & sgm
      zsg = 0.0d0
      do i = 1, nnrm
        zsg = zsg + tnrm(i)**2
      enddo
      zsg = zsg/dfloat(nnrm)
!
!::renormalization
      zfc = dsqrt(1.0d0/zsg)
      do i = 1, nnrm
        tnrm(i) = tnrm(i)*zfc
      enddo
!
      return
      end
