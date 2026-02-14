! added get the code number from the run-time parameters by kamata 2021/08/18
! remaked treat 4 or more impurities with IMPMC by kamata 2022/04/21
! get run-time argument parameter
      subroutine get_par( carg, marg, narg )
      implicit none
! arguments
      integer,   intent(in)  :: marg
      integer,   intent(out) :: narg
      character, intent(out) :: carg(marg)*10
! carg : run-time argument parameters
! marg : dimension size of cpara
! narg : number of parameters ( < 0 : too many parameters )

! local variables
      integer    i, nd
! i  : counter
! nd : number of daat to set

! functions
      integer    iargc

! get number of command line arguments
      narg = iargc( )

! check number of parameters
      nd = narg
      if( narg == 0 ) then
        carg(1:marg) = ' '
        return
      elseif( narg > marg ) then
        narg = -1
        nd = marg
      endif

! repeated number of arguments
      do i = 1, nd
! get command line argument
        call getarg( i, carg(i) )
      enddo 

      return
      end
