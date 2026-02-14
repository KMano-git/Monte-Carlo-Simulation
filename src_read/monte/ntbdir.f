!**********************************************************************
      subroutine ntbdir(bdir)
!**********************************************************************
!
!       b-field  (bx,by,bz)    b is a unit vector.
!
!       b-direction must be the same as to e_sita
!         metric in soldor  b*e_sita > 0
!                                             ! %%% 2003/09/24
!----------------------------------------------------------------------
      use cgdcom, only : grdx, grdy
      use cntcom, only : bvx, bvy, mcel
      use cplmet, only : icel, itpve, itpvs, itsle, itsls, jcel, jtmax
      use cunit,  only : n6
      implicit none
!
!::argument
      real(8), intent(out) :: bdir
!
!::local variables
      integer  it, jmax1, jt, j, i, ic
      integer  mx1, mx2, mx3, mx4, my1, my2, my3, my4
      real*8   x1, y1, x2, y2, zln, esx, esy, bx, by
      real*8   ebx, eby, sgn
!
      write(n6,'(/2x,"*** ntbdir ***")')
!
!::direction of e_sita
      write(n6,'(2x,"sign of e_sita*b")')
      write(n6,'(4x,"jt",3x,"it",3x,"j",4x,"i",4x,4x,"ic",4x,"sgn",7x,
     >  "x1",8x,"y1",8x,"esx",7x,"esy",7x,"bx",8x,"by")')
!
      do it = itsls+1, itpve-1
      jmax1 = jtmax(it)-1
      do jt = 2, jmax1
      if( jt.ne.2 .and. jt.ne.jmax1 ) cycle
      j = jcel(jt,it)
      i = icel(jt,it)
      call mcpnt(j,i,mx1,mx2,mx3,mx4,my1,my2,my3,my4)
!
      x1  = grdx(mx1,my1)
      y1  = grdy(mx1,my1)
      x2  = grdx(mx2,my2)
      y2  = grdy(mx2,my2)
      zln = dsqrt((x2-x1)**2+(y2-y1)**2)
      esx = (x2-x1)/zln
      esy = (y2-y1)/zln
!
      ic  = mcel(j,i)
      bx  = bvx(ic)
      by  = bvy(ic)
!
      zln = dsqrt(bx**2+by**2)
      ebx = bx/zln
      eby = by/zln
!
      sgn = esx*ebx+esy*eby
      if( it.eq.itsle ) bdir = sgn
!
      if( it.ge.itsle-1 .and. it.le.itpvs+1 ) then
      write(n6,'(2x,4i5,i7,2x,7f10.5)')
     >  jt,it,j,i,ic,sgn,x1,y1,esx,esy,ebx,eby
      endif
!
      enddo
      enddo
!
      if( bdir.gt.0.0 ) then
      write(n6,'(4x,"b-direction is correct in ntgdsk",
     >     "  b*e_sita > 0  ",f10.5)') bdir
      else
      write(n6,'(4x,"b-direction is incorrect in ntgsk",
     >     "  b*e_sita > 0  ",f10.5)') bdir
      endif
!
      return
      end
