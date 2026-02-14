! added replace all include files with module files by kamata 2021/08/18
      module cdfcom
      implicit none

!::file number & dimension size
      integer :: nfcd = 0, nrnk_max = 0, nvar_max = 0, ndat_max = 0

!::list of xs_var
      integer, parameter :: ndxvr = 20

      character :: xnam*40 = ' ', xvar(ndxvr)*40 = ' '
      integer   :: nxvr = 0, mxvr(ndxvr) = 0

!::list of variables
      integer, parameter :: ndhvr = 20

      character :: hvar(ndhvr)*20 = ' ', hvty(ndhvr)*6 = ' '
      integer   :: hvrk(ndhvr) = 0, nhvr = 0, mhvr(ndhvr) = 0

!::table
      integer, parameter :: ndtb = 100

      integer   :: ntab = 0, mtab(ndtb) = 0, tbrnk = 0
      character :: ctab(ndtb)*40 = ' ', tbkey*40 = ' ', tbtyp*1 = ' '

      end module cdfcom
