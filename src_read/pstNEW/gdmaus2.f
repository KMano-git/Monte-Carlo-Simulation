!***********************************************************************
      subroutine gdmaus2( kax, key, xcrd, ycrd )
!***********************************************************************
      implicit none
!
!::local common
      integer npen
      real*4  clx,cly,cox,coy,plx,ply,
     >        cpx1,cpx2,cpy1,cpy2,
     >        cfhl,cfhl2,cfhn,cfhn2,cfhn3,cfhh,cfh2,cfhw,hsiz,
     >        axln,ayln,axlmin,aylmin
      common /cgcomm/ npen,clx,cly,cox,coy,plx,ply,
     >                cpx1,cpx2,cpy1,cpy2,
     >                cfhl,cfhl2,cfhn,cfhn2,cfhn3,cfhh,cfh2,cfhw,hsiz,
     >                axln,ayln,axlmin,aylmin

!::argument
      real*4    xcrd,  ycrd
      integer   kax
      character key*1
!
!::local variables
      real*4  plnx, plny, vlnx, vlny
      real*4  xpn, ypn, fac
      character  mode*1, cpos*80
      integer icx, icy, lenx, mji
      character  cinp*20
!
      integer ifst, ldv, lpp
      real*4  xleg, yleg
      data   ifst/0/
      save   ifst, ldv, lpp, xleg, yleg
!
      write(6,'(/2x,"*** gdmaus2 ***")')
      write(6,'(2x,"gdlib cpx1,cpx2,cpy1,cpy2 =",4f9.4)')
     >   cpx1,cpx2,cpy1,cpy2
!
!::set position of legend
!         ldv = 1 (V) / = 2 (X),  lpp = 1 (Land) / = 2 (Port)
      if( ifst.eq.0 ) then
      ifst = 1
      call gtdevinf(ldv,lpp)
      write(6,'(2x,"Term =",i2," Paper =",i2)') ldv, lpp
      call getenv("POSLG",cpos)
      read( cpos(1:lenx(cpos)),* ) xleg, yleg
      write(6,'(2x,"position of legend = ",2f8.3)') xleg, yleg
      endif
!
      if( xleg.gt.0.0 .and. yleg.gt.0.0 ) then
      key = "a"
      xcrd = xleg
      ycrd = yleg
      write(6,'(2x,"xleg,yleg =",2f9.4)') xleg, yleg
      endif
!
!::VersaTerm
 100  continue
      if( ldv.eq.1 ) then
      call gdput(" move cursor to position of legend and crick")
      plnx = 29.7
      plny = 29.7
      plnx = plnx*(11.1/12.0)
      plny = plny*(11.1/12.0)
      vlnx = 4096.0
      vlny = 4096.
      call gin( mode, icx, icy )
      write(6,'(2x,"icx,icy =",2i10)') icx,icy
!
      xcrd = plnx/vlnx*icx
      ycrd = plny/vlny*icy
      key  = mode
      else
!
!::X-term
      xcrd = cpx2 + 0.6
      ycrd = cpy2 - 3.0
      endif
!
 200  continue
      write(6,'(2x,"gdmaus  (x,y) =",2f10.3)') xcrd, ycrd
!
      return
      end

