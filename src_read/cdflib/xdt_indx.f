!***********************************************************************
      subroutine xdt_indx(krnk,dnam,cvar,ixs)
!***********************************************************************
!
!   krnk  : (i) int   rank of variables
!                if you don't know, krnk = -1
!                you can get rank,  prnk(ixs)
!
!   dnam  : (i) char  data-name       dd_00_elastic
!   cvar  : (i) char  variable name   reaction_rate
!   ixs   : (o) int   pointer
!
!     rank  = xdrnk(ixs)
!     y-name=xwvar(1,ixs),  y-unit=xwunt(1,ixs),  y-mult=xwmlt(1,ixs)
!     x-name=xwvar(2,ixs),  x-unit=xwunt(2,ixs),
!     x-min =xwmin(2,ixs),  x-max =xwmax(2,ixs),  x_mult=xwmlt(2,ixs)
!     y-size=xwnum(1,ixs),  x-size=xwnum(2,ixs)
!     y-data=xdaty(j)  j=xjsta(ixs),...xjend(ixs)
!
!-----------------------------------------------------------------------
      use cxdcom, only : nxs, prk, xdnam, xwvar
      implicit none
!
!::argument
! modified 2/3 lines organize local variables and include files by kamata 2021/06/16
!ik   character dnam*(*), cvar*(*)
!ik   integer   krnk
      character, intent(in)  :: dnam*(*), cvar*(*)
      integer,   intent(in)  :: krnk
      integer,   intent(out) :: ixs

!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   integer   lenx, mjd, mjv, i, ixs, irnk
      integer   mjd, mjv, i, irnk
      character cmsg*80
! added 2 lines organize local variables and include files by kamata 2021/06/16
! function
      integer    lenx

!-----------------------------------------------------------------------
!
      mjd = lenx(dnam)
      mjv = lenx(cvar)
!
!::ixs
      ixs = 0
      do i = 1, nxs
      if( xdnam(i)(1:lenx(xdnam(i))).eq.dnam(1:mjd) .and.
     >    xwvar(1,i)(1:lenx(xwvar(1,i))).eq.cvar(1:mjv) ) then
      ixs = i;  exit
      endif
      enddo
!
      if( ixs.eq.0 ) then
      write(cmsg,'("no-found var. of dnam =",a,"  cvar =",a)')
     >   dnam(1:mjd),cvar(1:mjv)
      call wexit("xt_indx",cmsg)
      endif
!
!::rank
      irnk = prk(ixs)
      if( krnk.ge.0 ) then
      if( krnk.ne.irnk ) then
      write(cmsg,'("incorrect rank   cvar =",a,2i4)')
     >    cvar(1:mjv),krnk,irnk
      call wexit("xdt_indx",cmsg)
      endif
      endif
!
!::debug write
!x      write(6,'(2x,"*** xdt_indx ***  [",a,"]  ixs =",i3,"  cvar =",a
!x     >  ,"  nsiz =",i7,"  ips,ipe =",2i7)')
!x     >  cdsn(1:lenx(cdsn)),ixs,cvar(1:lenx(cvar)),jnm,jst,jen
!x      call xdt_debg(ixs)
!
      return
      end
