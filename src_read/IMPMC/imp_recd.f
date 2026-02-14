!**********************************************************************
      subroutine imp_recd
!**********************************************************************
      use cimctl, only : cdirz, icalz, kimp, nimp, tmimp
      use csonic, only : itim, time
      use cunit,  only : n6
      implicit none
!
!::local variables
      real*8   tmb
!
!::define kimp
      if( kimp /= 1 ) return
!
!::calculation
      icalZ = icalZ + 1
      tmb   = tmimp
      tmimp = time
      cdirZ = "dIMP_tmp"

!
      write(n6,'(/2x,30("*"),2x,a,2x,"itim =",i6,"  time =",1pe14.6,
     >   "  timeb =",1pe14.6,"  icalZ =",i4/)')
     >   "IMPMC-START", itim, time, tmb, icalZ

!
      if( kimp > 0 ) then
        write(n6,'(2x,"  itim/nimp/kimp =",i8,2i5)') itim, nimp, kimp
      endif
      return
      end
