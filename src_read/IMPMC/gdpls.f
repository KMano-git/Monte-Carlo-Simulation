!**********************************************************************
      subroutine gdpls(j,i,ane,ani,ate,ati,acs)
!**********************************************************************
      use cplcom, only : vcs, vne, vni, vte, vti
      implicit none
!
!::argument
      integer, intent(in)  :: j, i
      real(8), intent(out) :: ane, ani, ate, ati, acs
!
      if( j.le.0 .or. i.le.0 ) then
        ane = 0.0d0
        ani = 0.0d0
        ate = 0.0d0
        ati = 0.0d0
        acs = 0.0d0
      else
        ane = vne(j,i)
        ani = vni(j,i)
        ate = vte(j,i)
        ati = vti(j,i)
        acs = vcs(j,i)
      endif
!
      return
      end
