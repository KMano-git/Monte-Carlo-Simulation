!**********************************************************************
      subroutine dtyplnk(dtyp,sgrp,rgrp,kcon)
!**********************************************************************

      implicit none

      character(*), intent(in) :: dtyp
      integer, intent(in)  :: sgrp, rgrp
      integer, intent(out) :: kcon

      kcon = 0
      return
      end subroutine dtyplnk
