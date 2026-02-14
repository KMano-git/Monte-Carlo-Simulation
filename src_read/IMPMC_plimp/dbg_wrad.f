!***********************************************************************
      subroutine dbg_wrad(cmsg,nv,wrd,wci,dnz)
!***********************************************************************
      use csize
      use cimcom
      use cimden
      use cntpls
      use cplcom
      use csonic
      use cunit
      implicit none

 !  twne(ndmc): Ne at IMPMC cal 

!
!       No define sptyc
!
!::argument
      character(*) :: cmsg
      integer ::  nv
      real(8) ::  wrd(ndmc), wci(ndmc), dnz(0:ndis,ndmc)
!
!::local variables
      integer :: tbic(6) = (/2675, 2677, 2679, 2681, 2683, 2685/)
      integer :: n, ia, ic, ix, iy
      integer :: imox, imoy
      real(8) :: twr, twr_rg(10), twr_it(ndy)
!
      if( lmype .ne. lmspe ) return
!
      write(n6,'(/2x,"----- dbg_wrad ----- start  itim = ",i8)') itim
      write(n6,'(2x,a,2x,a,"  sflux =",1pe12.4,"  swtot =",1pe12.4)')
     >   trim(cmsg), sptyc, sflux, swtot
      write(n6,'(4x,"ic",4x,"ix",4x,"iy",3x,"wrd",9x,"twne",8x,
     >  "dnz15",7x,"dene",8x,"deni",8x,"teme",8x,
     >  "vne",9x,"vna",9x,"vte")')
!
      ia = 1
      do n = 1, 6
        ic = tbic(n)
        ix = imox(ic)
        iy = imoy(ic)
        write(n6,'(2x,i6,i6,i5,1p9e12.4)') ic, ix, iy, 
     >   wrd(ic), twne(ic), dnz(15,ic), dene(ic), deni(ic,ia), teme(ic),
     >   vne(ix,iy), vna(ix,iy,ia), vte(ix,iy)
      enddo
!
!
!::integration
      call totrgn(wrd, twr, twr_rg, twr_it)
!
      write(n6,'(2x,"totwr =",1pe12.3,"  wr_rg =",1p6e12.3,
     >    "  wr_core =",1pe12.3)') twr, twr_rg(1:6), twr_rg(7)
      write(n6,'(2x,"----- dbg_wrad ----- end    itim = ",i8/)') itim
!
      return
      end
