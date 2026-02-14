!*********************************************************************
      subroutine imp_lnch(temt,tmend)
!*********************************************************************
      use cimcom, only : igemt, mypcn, myptl, mywgt, ndmp, npemt, npmax
     >    , ntags, pemt, tt, wght
      use cimden, only : nsrc_spt
      use cimntl, only : stb
      use cunit,  only : n6
      implicit none

!::argument
      real(8), intent(in) :: temt, tmend

!::local variables
      real(8) :: ptim, ptmx
      integer :: i, ip

      ip = npmax
      do i = 1, npemt
        ip = ip + 1
        if( ip > ndmp ) goto 910

        call imemt_wall(ip)
        tt(ip)  = temt   ! time of particle at the present time
        stb(ip) = temt   ! birth time of neutral  <IMPORTANT>

        myptl(igemt) = myptl(igemt) + 1.0d0
        mywgt(igemt) = mywgt(igemt) + wght(ip)
        mypcn(igemt) = mypcn(igemt) + wght(ip)*pemt(ip)

        ptmx = tmend
        ptim = tt(ip)
        ntags(ip) = nsrc_spt
        call lst_trak(ip,"lnch",ptmx,ptim,stb(ip))  ! TRAC
      enddo

      npmax = ip

      return

 910  continue
      write(n6,'(2x,"stop at sub. imp_lnch 11/24 ",2i6)') ip, ndmp
      call wexit("imp_lnch  11/24","ip > ndmp")
      stop
      end
