!***********************************************************************
      subroutine cdf_getvr(ivar,nvdt,cvdt,mvdt,ndvd)
!***********************************************************************
      use cdfcom, only : ctab, mtab, mxvr, nxvr, tbrnk, xvar
      implicit none
!
!::arguments
! modified 3/3 lines organize local variables and include files by kamata 2021/06/16
!ik   integer   ndvd
!ik   character cvdt(ndvd)*(*)
!ik   integer   ivar,nvdt,mvdt(ndvd)
      integer,   intent(in)  :: ivar, ndvd
      character, intent(out) :: cvdt(ndvd)*(*)
      integer,   intent(out) :: nvdt, mvdt(ndvd)

!::local varaibales
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   integer    lenx, i0, ii, i, ist, ien
      integer    i0, ii, i, ist, ien
      character  cmsg*80
!
      if( ivar.le.0 .or. ivar.gt.nxvr ) goto 910
!
      i0 = tbrnk*(ivar-1)
      ist = i0 + 1
      ien = i0 + tbrnk
      if( tbrnk.le.0 ) ien = ist
!
      ii = 0
      do i = ist, ien
      ii = ii + 1
      if( ii.gt.ndvd ) goto 920
      cvdt(ii) = ctab(i)
      mvdt(ii) = mtab(i)
      enddo
      nvdt = ii
!
      return
      write(6,'(2x,"*** cdf_getvr ***  ",i2,"  [",a,"]")')
     >   ivar,xvar(ivar)(1:mxvr(ivar))
      write(6,'(4(4x,a20))') (cvdt(i)(1:mvdt(i)),i=1,nvdt)
!
      return
!
!::error
 910  continue
      write(cmsg,'("wrong ivar  ivar<=nxvr ",2i4)') ivar,nxvr
      call wexit("cdf_getvr",cmsg)
 920  continue
! modified 1/1 lines bug by kamata 2022/05/26
!ik   write(cmsg,'("dimension error  ii.gt.ndvd")') ii,ndvd
      write(cmsg,'("dimension error  ii.gt.ndvd",2i4)') ii,ndvd
      call wexit("cdf_getvr",cmsg)
      end
