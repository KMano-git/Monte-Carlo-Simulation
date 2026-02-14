!**********************************************************************
      subroutine plenc(it,alnx)
!**********************************************************************
!
!       |.*.|..x..|..x..|..x..|..x..|..x..|
!   J     1 |  2     3    j-1    j    j+1      jhm = j-1/2
!           |              <--|-->             dxm = hdxm(jhm)
!           |              dxm dxp             dxp = hdxp(jhm)
!           |-------------------->  alnx(j)
!
!----------------------------------------------------------------------
      use cplmet, only : hdxm, hdxp, icel, itmps, jcel, jtmax, jtmin
      use csize,  only : ndx
      implicit none
!
!::argument
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  it
!ik   real*8   alnx(ndx)
      integer, intent(in)  :: it
      real(8), intent(out) :: alnx(ndx)
!
!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  j, jts, jte, jsm, ic, jt, jc
      integer  j, jts, jte, jsm, ic, jt
      real*8  zln, zdx
!
      do 110 j = 1, ndx
      alnx(j) = 0.0
 110  continue
!
      zln = 0.0d0
      jts = jtmin(it) + 1
      jte = jtmax(it)
      if( it.ge.itmps ) then
        jte = jte - 1
        jsm = jcel(jts-1,it)
        ic  = icel(jts,it)
        zln = -hdxm(jsm,ic)
      endif
      do 120 jt = jts, jte
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   jc  = jcel(jt,it)
      ic  = icel(jt,it)
      jsm = jcel(jt-1,it)
      zdx = hdxm(jsm,ic) + hdxp(jsm,ic)
      zln = zln + zdx
      alnx(jt) = zln
 120  continue
!
      return
      end
!
!**********************************************************************
      subroutine plenb(it,alnx)
!**********************************************************************
!
!       |.*.|..x..|..x..|..x..|..x..|..x..|
!   J     1 |  2     3    j-1    j    j+1         jhm = j-1/2
!           |                 |-><--|             dxp = hdxm(jhm)
!           |                 dxp dxm             dxm = hdxp(jhp)
!           |----------------------->  alnx(j)
!
!----------------------------------------------------------------------
      use cplmet, only : hdxm, hdxp, icel, jcel, jtmax, jtmin
      use csize,  only : ndx
      implicit none
!
!::argument
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  it
!ik   real*8   alnx(ndx)
      integer, intent(in)  :: it
      real(8), intent(out) :: alnx(ndx)
!
!::local variables
      integer j, jts, jte, jt, jc, ic, jsm, jsp
      real*8  zln
!
      do 110 j = 1, ndx
      alnx(j) = 0.0
 110  continue
!
      zln = 0.0d0
      jts = jtmin(it) + 1
      jte = jtmax(it) - 1
      do 120 jt = jts, jte
      jc  = jcel(jt,it)
      ic  = icel(jt,it)
      jsm = jcel(jt-1,it)
      jsp = jc
      zln = zln + hdxp(jsm,ic) + hdxm(jsp,ic)
      alnx(jt) = zln
 120  continue
!
      return
      end
!
!**********************************************************************
      subroutine slenc(it,alnx)
!**********************************************************************
!
!       |.*.|..x..|..x..|..x..|..x..|..x..|
!   J     1 |  2     3    j-1    j    j+1      jhm = j-1/2
!           |              <--|-->             dxm = hdxm(jhm)
!           |              dxm dxp             dxp = hdxp(jhm)
!           |-------------------->  alnx(j)
!
!----------------------------------------------------------------------
      use cplmet, only : hdxm, hdxp, hpit, icel, itmps, jcel, jtmax
     >    , jtmin
      use csize,  only : ndx
      implicit none
!
!::argument
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  it
!ik   real*8   alnx(ndx)
      integer, intent(in)  :: it
      real(8), intent(out) :: alnx(ndx)
!
!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  j, jts, jte, jsm, ic, jt, jc
      integer  j, jts, jte, jsm, ic, jt
      real*8   zln, zdx
!
      do 110 j = 1, ndx
      alnx(j) = 0.0
 110  continue
!
      zln = 0.0d0
      jts = jtmin(it) + 1
      jte = jtmax(it)
      if( it.ge.itmps ) then
        jte = jte - 1
        jsm = jcel(jts-1,it)
        ic  = icel(jts,it)
        zln = -hdxm(jsm,ic)/hpit(jsm,ic)
      endif
      do 120 jt = jts, jte
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   jc  = jcel(jt,it)
      ic  = icel(jt,it)
      jsm = jcel(jt-1,it)
      zdx = (hdxm(jsm,ic) + hdxp(jsm,ic))/hpit(jsm,ic)
      zln = zln + zdx
      alnx(jt) = zln
 120  continue
!
      return
      end
!
!**********************************************************************
      subroutine slenb(it,alnx)
!**********************************************************************
!
!      |.*.|..x..|..x..|..x..|..x..|..x..|
!   J     1 |  2     3    j-1    j    j+1         jhm = j-1/2
!           |                 |-><--|             dxp = hdxm(jhm)
!           |                 dxp dxm             dxm = hdxp(jhp)
!           |----------------------->  alnx(j)
!
!----------------------------------------------------------------------
      use cplmet, only : hdxm, hdxp, hpit, icel, jcel, jtmax, jtmin
      use csize,  only : ndx
      implicit none
!
!::argument
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  it
!ik   real*8   alnx(ndx)
      integer, intent(in)  :: it
      real(8), intent(out) :: alnx(ndx)
!
!::local variables
      integer j, jts, jte, jt, jc, ic, jsm, jsp
      real(8) hpitm, hpitp
      real*8  zln
!
      do j = 1, ndx
      alnx(j) = 0.0
      enddo
!
      zln = 0.0d0
      jts = jtmin(it) + 1
      jte = jtmax(it) - 1
      do jt = jts, jte
      jc  = jcel(jt,it)
      ic  = icel(jt,it)
      jsm = jcel(jt-1,it)
      jsp = jc
      !FIXME:: 0div reported SYamoto
!      if(hpit(jsm,ic)==0.0d0 .or. hpit(jsp,ic)==0.0d0)then
!      zln = zln + hdxp(jsm,ic)/hpit(jsm,ic) + hdxm(jsp,ic)/hpit(jsp,ic)
      ! 0-div occur at corner of domain
      ! if hpit=0, then set hpit=1 tentatively
      hpitm = hpit(jsm,ic)
      hpitp = hpit(jsp,ic)
      if(hpitm == 0.0d0) hpitm = 1.0d0
      if(hpitp == 0.0d0) hpitp = 1.0d0
      zln = zln + hdxp(jsm,ic)/hpitm + hdxm(jsp,ic)/hpitp
      alnx(jt) = zln
      enddo
!
      return
      end
