!***********************************************************************
      subroutine plupwd
!***********************************************************************
!
!       lupwd = 0   no up-wind
!       lupwd = 1   up-wind in private region
!       lupwd = 2   up-wind in whole region
!
!-----------------------------------------------------------------------
      use cplcom, only : iupwd, lupwd
      use cplmet, only : icel, icmpe, icspx, icwl2, itmax, jcel, jcxp1
     >    , jcxp2, jtmax
      use csize,  only : ndx, ndy
      use cunit,  only : n6
      implicit none
!
!::local viariables
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  nsiz, it, jmax, jw, j, i, iupw, ipmx
      integer  nsiz, it, jw, j, i, iupw, ipmx
!
      write(n6,'(/2x,"*** plupwd ***   lupwd =",i3)') lupwd
!
      nsiz = ndx*ndy
      call seti( iupwd, nsiz, 9 )
!
      do it = 1, itmax
! modified 2/1 lines organize local variables and include files by kamata 2021/05/31
!ik   jmax = jtmax(it)
!ik   do jw = 1, jmax
      do jw = 1, jtmax(it)
      j = jcel(jw,it)
      i = icel(jw,it)
      iupw = 0
!-----
      if( lupwd.eq.1 ) then
      if( j.le.jcxp1 .or. j.ge.jcxp2 ) then
      if( i.ge.icspx ) iupw = 1
      endif
      endif
      if( lupwd.eq.2 ) iupw = 1
!-----
      iupwd(j,i) = iupw
      enddo
      enddo
!
!::debug write
      ipmx = max0( icwl2, icmpe ) + 3
      do i = 1, ipmx
      write(n6,'(2x,i3,2x,170i1)') i,(iupwd(j,i),j=1,ndx)
      enddo
!
      return
      end
