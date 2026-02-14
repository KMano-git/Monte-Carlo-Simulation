!***********************************************************************
      subroutine gdfleg_new( xp, yp, ctbl, kclr, ntbl, iord )
!***********************************************************************
!
!                        i   i   i     i     i
!     xp   : x coordinate of left-upper corner (cm)
!     yp   : y coordinate of left-upper corner (cm)
!     ctbl : class tabel (value)
!     ntbl : number of class
!     iord : order of class ( = 1: increase, = 2: decrease )
!             if iord >= 10 then no plot the lowest class
!     kclr : colour
!
!
! ldev  : device   = 1 : VersaTerm,  = 2 : X-lib
! ltdev : same as ldev, but only trmget & trmput
! lsymfg : center symbol output flag for Macdraw
!          = 0: normal  = 1: center symbol  = 2: not output
! lpapaer : papaer = 1 : Landscape,  = 2 : Portrait
!-----------------------------------------------------------------------
      implicit none
!
!::local common
!xx   include 'clattr'
      integer ldev, ltdev, lsymfg, lpaper
      common /cmsym/ ldev, ltdev, lsymfg, lpaper
!
!::argument
      real*4   xp, yp, ctbl(ntbl)
      integer  kclr(31), ntbl, iord
!
!::local variables
      real*4   hei, he2, dy, xs, ys, xe, ye, xx(5), yy(5), xv, yv
      integer  n, jord, jplt, is, ia, i, iclr, ll, lenx, j
      character  txt*20
!
      hei = 0.2
      he2 = hei * 0.5
      dy  = 0.4
      xs  = xp + 0.2
      ys  = yp
      xe  = xs + 0.8
      n   = 5
      jord = mod( iord, 10 )
      jplt = iord / 10
      if( jord .eq. 2 ) then
        is = 0
        ia = 1
      else
        is = ntbl + 2
        ia = -1
      endif
      do i=1,ntbl+1
        is = is + ia
cik     if( is .eq. 1 ) cycle
        ye = ys - dy
!
!:: paint box
        xx(1) = xs
        yy(1) = ys
        xx(2) = xe
        yy(2) = yy(1)
        xx(3) = xx(2)
        yy(3) = ye
        xx(4) = xx(1)
        yy(4) = yy(3)
        xx(5) = xx(1)
        yy(5) = yy(1)
!
!::painting
        iclr = kclr(is)
        if( iclr.gt.0 ) then
        if( ldev .eq. 1 ) then
          call v_fillpoly( iclr, 0, xx, yy, n )
        else
          call x_fillpoly( iclr, 0, xx, yy, n )
        endif
        endif
        if( i .gt. ntbl ) exit
!
!:: plot value
        call newpen( 1 )
        xv  = xe + 0.2
        yv  = ye - he2
        write(txt,'(1pg12.5)') ctbl(is-1)
        ll = lenx( txt )
        call symbol( xv, yv, hei, txt, 0.0, ll )
        ys = ye
      enddo
!
!:: draw box
      xx(1) = xs
      yy(1) = yp
      xx(2) = xe
      yy(2) = yy(1)
      xx(3) = xx(2)
      yy(3) = yp - ( ntbl + 1 ) * dy
      xx(4) = xx(1)
      yy(4) = yy(3)
      xx(5) = xx(1)
      yy(5) = yy(1)
      call newpen( 1 )
      call plot( xx(1), yy(1), 3 )
      do j=2,n
        call plot( xx(j), yy(j), 2 )
      enddo
      ys = yp
      do i=1,ntbl
        ys =  ys - dy
        call plot( xs, ys, 3 )
        call plot( xe, ys, 2 )
      enddo
!
      return
      end
