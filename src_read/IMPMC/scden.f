!**********************************************************************
      subroutine scden(nat,isat,tmat,nmv,icmv,tmmv,
     >       np, ptm, pdt, pis, pmv, ndsc)
!**********************************************************************
      implicit none
!
!::argument
! modified 5/6 lines organize local variables and include files by kamata 2021/06/28
!ik   integer  ndsc
!ik   integer  nat, isat(ndsc), nmv, icmv(ndsc)
!ik   real*8   tmat(ndsc), tmmv(ndsc)
!ik   integer  np, pis(ndsc), pmv(ndsc)
!ik   real*8   ptm(ndsc), pdt(ndsc)
      integer, intent(in)    :: ndsc
      integer, intent(in)    :: nat, nmv
      integer, intent(inout) :: isat(ndsc), icmv(ndsc)
      real(8), intent(inout) :: tmat(ndsc), tmmv(ndsc)
      integer, intent(out)   :: np, pis(ndsc), pmv(ndsc)
      real(8), intent(out)   :: ptm(ndsc), pdt(ndsc)
!
!::local variable
      integer  ii, i1, i2, i
! modified 5/1 lines organize local variables and include files by kamata 2021/06/28
!ik   integer  j, idbg, imox, imoy
!ik   real*8   tm, tmb, tm1, tm2
!ik   integer  lp
!ik   data lp/0/
!ik   save lp
      real*8   tm, tm1, tm2
!
!::KS_100319
      tmat(nat+1) = tmat(nat)
      isat(nat+1) = isat(nat)
      tmmv(nmv+1) = tmmv(nmv)
      icmv(nmv+1) = icmv(nmv)
!::KS_100319
!
      tm = 0.0d0
! deleted 1 line organize local variables and include files by kamata 2021/06/28
!ik   tmb = 0.0d0
      ii = 0
      i1 = 1
      i2 = 1
      tm1 = tmat(i1)
      tm2 = tmmv(i2)
 100  continue
      if( tm1.lt.tm2 ) then
        tm = tm1
        i1 = i1 + 1
      elseif( tm1.eq.tm2 ) then
        tm = tm1
        i1 = i1 + 1
        i2 = i2 + 1
      else
        tm = tm2
        i2 = i2 + 1
      endif
!
!::Note  i1 > nat,  i2 > nmv
!
      ii  = ii + 1
      tm1 = tmat(i1)
      tm2 = tmmv(i2)
      pis(ii) = isat(i1)
      pmv(ii) = icmv(i2)
      if( tm1.gt.tm ) pis(ii) = isat(i1-1)
      if( tm2.gt.tm ) pmv(ii) = icmv(i2-1)
      ptm(ii) = tm
! deleted 1 line organize local variables and include files by kamata 2021/06/28
!ik   tmb = tm
      if(i1.gt.nat .or. i2.gt.nmv ) goto 150
      goto 100
 150  continue
      np = ii
!
      do i = 1, np-1
      pdt(i) = ptm(i+1)-ptm(i)
      enddo
      pdt(np) = 0.0d0
      pis(np) = 0
      pmv(np) = 0
!
!::debug write
!x      lp = lp + 1
!x      write(n6,'(/2x,"lp =",i5)') lp
!x      write(n6,'(2x,"SCDEN nat, nmv, np =",3i5)')   nat, nmv, np
!x      write(n6,'(2x,"SCDEN tmat =",1p10e12.3)')    (tmat(j),j=1,nat)
!x      write(n6,'(2x,"SCDEN tmmv =",1p10e12.3)')    (tmmv(j),j=1,nmv)
!x      write(n6,'(2x,"SCDEN ptm  =",1p10e12.3)')    (ptm(j),j=1,np)
!x      write(n6,'(2x,"SCDEN pdt  =",1p10e12.3)')    (pdt(j),j=1,np)
!x      write(n6,'(2x,"SCDEN pis  =",10(3x,i5,4x))') (pis(j),j=1,np)
!x      write(n6,'(2x,"SCDEN pix  =",10(3x,i5,4x))') (imox(pmv(j)),j=1,np)
!x      write(n6,'(2x,"SCDEN piy  =",10(3x,i5,4x))') (imoy(pmv(j)),j=1,np)
!x      if(lp.gt.100) call wexit("scden","too many loop")
!
      return
      end
