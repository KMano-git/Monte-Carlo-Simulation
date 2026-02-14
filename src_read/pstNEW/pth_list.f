!**********************************************************************
      subroutine pth_list(ncol,posx,posy,ipos,ctyp,pvar,nft)
!**********************************************************************
      use csize
      use cntpls
      use cntsrc
      use cphcns
      implicit none
!
!::argument
      integer :: ncol
      real*8,  dimension(ncol) :: posx, posy
      integer, dimension(ncol) :: ipos
      integer :: nft
!
      character*(*) :: ctyp
      real*8,  dimension(ndmc) :: pvar
!
!::local variables
      real*8   sum, xp, yp, dl
      integer  i, ic, imox, imoy
!
      if( nft.gt.0 ) then
      write(nft,'(5x,"name = [",a,"]")') trim(ctyp)
      write(nft,'(5x,"i",3x,"ix",3x,"iy",3x,"xp",7x,"yp",6x,"dl",9x,
     >  "Te",9x,"Ti",9x,"Ne",9x,"N0",9x,"var",8x,"intg")')
      endif
!
      sum = 0.0d0
      do i = 1, ncol-1
      xp = posx(i)
      yp = posy(i)
      ic = ipos(i)

      dl = sqrt( (posx(i+1)-posx(i))**2 + (posy(i+1)-posy(i))**2 )
      sum = sum + pvar(ic)*dl
!
      if( nft.gt.0 ) then
      write(nft,'(2x,i4,2i5,2f9.4,1p9e11.3)')
     >   i, imox(ic), imoy(ic), posx(i), posy(i), dl,
     >   teme(ic), temi(ic), dene(ic), tden0(ic,1),
     >   pvar(ic), sum/(4.0d0*cpi)
      endif
!
      enddo ! i
!
      return
      end
