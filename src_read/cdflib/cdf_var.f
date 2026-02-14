!**********************************************************************
      subroutine cdf_var(nvar,cvar,ndvr)
!**********************************************************************
      use cdfcom, only : ctab, mtab, nrnk_max, ntab, xnam
      implicit none
!
!::argument
! modified 2/3 lines organize local variables and include files by kamata 2021/06/16
!ik   integer nvar,ndvr
!ik   character cvar(ndvr)*(*)
      integer,   intent(in)  :: ndvr
      integer,   intent(out) :: nvar
      character, intent(out) :: cvar(ndvr)*(*)

!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   integer   i6, ndat, ii, i, lenx
      integer   i6, ndat, ii, i
      character ckey*40
! added 2 lines organize local variables and include files by kamata 2021/06/16
! function
      integer    lenx
!
      i6 = 6
!
      ckey = "xs_data_size"
      call cdf_getky(ckey)
      if( ntab.ne.1 ) goto 910
      read(ctab(1)(1:mtab(1)),*) ndat
!
      ckey = "xs_num_dep_var"
      call cdf_getky(ckey)
      if( ntab.ne.1 ) goto 920
      read(ctab(1)(1:mtab(1)),*) nvar
!
      ckey = "xs_name"
      call cdf_getky(ckey)
      if( ntab.ne.1 ) goto 930
      xnam = ctab(1)(1:mtab(1))
!
      ckey = "xs_var"
      call cdf_getky(ckey)
      if( ntab.le.0 ) goto 940
      ii = 0
      do i = 1, ntab, nrnk_max
      ii = ii + 1
      cvar(ii) = ctab(i)(1:mtab(i))
      if( ii.eq.nvar ) exit
      enddo
!
      write(i6,'(/2x,"*** cdf_var ***")')
      write(i6,'(4x,"xnam = [",a,"]")') xnam(1:lenx(xnam))
      write(i6,'(4x,"nvar =",i3,"  ndat =",i7)') nvar,ndat
      write(i6,'(5(4x,a))') (cvar(i)(1:lenx(cvar(i))),i=1,nvar)
!
      return
!
!::error
 910  continue
 920  continue
 930  continue
 940  continue
      call wexit("cdf_var","no found "//ckey(1:lenx(ckey)))
!
      end
