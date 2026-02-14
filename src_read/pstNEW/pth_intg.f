!**********************************************************************
      subroutine pth_intg(ncol,posx,posy,ipos,pvar,vintg)
!**********************************************************************
      use csize
      use cphcns
      implicit none


!
!::argument
      integer :: ncol
      real*8,  dimension(ncol) :: posx, posy
      integer, dimension(ncol) :: ipos
      real*8,  dimension(ndmc) :: pvar
      real*8  :: vintg
!
!::local variables
      real*8   sum, xp, yp, dl
      integer  i, ic, imox, imoy
! 
      sum = 0.0d0
      do i = 1, ncol-1
      xp = posx(i)
      yp = posy(i)
      ic = ipos(i)
!
      dl = sqrt( (posx(i+1)-posx(i))**2 + (posy(i+1)-posy(i))**2 )
      sum = sum + pvar(ic)*dl
      enddo
!
      vintg = sum/(4.0d0*cpi)
!
      return
      end
