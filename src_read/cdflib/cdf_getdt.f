!***********************************************************************
      subroutine cdf_getdt(ckey,ibas,iinc,num,y,ndt)
!***********************************************************************
      use cdfcom, only : nfcd
      implicit none
!
!::arguments
! modified 4/5 lines organize local variables and include files by kamata 2021/06/16
!ik   integer ndt
!ik   integer ibas,iinc,num
!ik   real*8  y(ndt)
!ik   character ckey*(*)
      integer,   intent(in)  :: ndt
      integer,   intent(in)  :: ibas, iinc
      integer,   intent(out) :: num
      real(8),   intent(out) :: y(ndt)
      character, intent(in)  :: ckey*(*)

!::local varaibales
! modified 3/3 lines organize local variables and include files by kamata 2021/06/16
!ik   integer  ii, jj, kk, ml, ns, n, lenx, me, i, nk, nlst
!ik   integer  mjs(110), mje(110), nrnk
!ik   character  clin*150, cmsg*120, ctyp
      integer  ii, jj, kk, ml, ns, n, me, i, nk
      integer  mjs(110), mje(110)
      character  clin*150, cmsg*120
      integer  i1, i2, i3, i4
! added 2 lines organize local variables and include files by kamata 2021/06/16
! function
      integer    lenx
!
      ml = lenx(ckey)
      ii = 0
      jj = 0
      kk = 0
!
 100  continue
      read(nfcd,'(a)',end=910) clin
      call linsep(clin,'=",;',nk,mjs,mje,110)
      if( nk.le.0 ) goto 100
      if( clin(mjs(1):mje(1)).eq.ckey(1:ml) ) kk = 1
      if( kk.gt.0 ) then
      ns = 1; if(kk.eq.1) ns = 2
!--
      do n = ns, nk
      ii = ii + 1
      if( ii.gt.ibas ) then
        jj = jj + 1
        if( jj.gt.ndt ) goto 920
        read(clin(mjs(n):mje(n)),*) y(jj)
        if( ii.eq.ibas+iinc ) goto 110
      endif
      enddo
!--
      kk = kk + 1
      me = lenx(clin)
      if( clin(me:me).eq.";" ) goto 110
      endif
      goto 100
!
 110  continue
      num = jj
!
      return
      write(6,'(2x,"*** cdf_getdt ***  num =",i7,"  ibas =",i7,
     >   "  iinc =",i7)') num,ibas,iinc
      i1 = 1; i2 = min0(i1+5,num)
      i3 = max0( 1, num-5 ); i4 = num
      write(6,'(4x,2i7,2x,1p6e12.3)') i1,i2,(y(i),i=i1,i2)
      write(6,'(4x,2i7,2x,1p6e12.3)') i3,i4,(y(i),i=i3,i4)
!
      return
!
!::error
 910  continue
      call wexit("cdf_getky","no found "//ckey(1:lenx(ckey)))
 920  continue
      write(cmsg,'("too many data  jj > ndt  ",2i7)') jj,ndt
      call wexit("cdf_getdt",cmsg)
      end
