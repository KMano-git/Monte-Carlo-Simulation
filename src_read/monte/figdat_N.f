!**********************************************************************
      subroutine figdat( cpath, ctime )
!**********************************************************************
!
!                  IMPMC                soldor (call plsorc)
!                                       /cntsrc/
!
!        a005128/PRJ/sonicV3/src/soldor     2013/06/05
!
!  group in charge of the code   2-GRP type
!    [PLS]: 1 / NTL]: 2 /[IMP]: 2 /  <== cdgr(i) i=1,2,3(ncode)
!
!----------------------------------------------------------------------
      use cunit, only : lmspe, lmype, mype, n6
      implicit none
! arguments
      character, intent(in) :: cpath*(*), ctime*(*)
! cpath : path name
! ctime : time

!::local variables
      integer :: i60
      integer :: kio, nft, n6sv
      character :: cdsn*200, ched*20, ctyp*10
!
      if( lmype.ne.lmspe ) return
      write(n6,'(/2x,"*** figdat ***",2x,a)') "DTNTL  w"
!
      i60 = 180000 + mype
!
!---------------------------------------------------------------------
!::DTNTL : neutral data
!---------------------------------------------------------------------
!   vtime, vflux, vmwork, vitim, vitnt, vnsrc, vsmty, vsmno
!   visrc, vkflx, vksty,  vnsmp, vcsrc
!  --> NTL_2
!
      ctyp="DTNTL"
      kio = 1
!
      n6sv = n6
      n6   = i60
      nft = 21 ! KH150906
!
      if( ctime == ' ' ) then
        cdsn = "./" // trim(ctyp)
      else
        cdsn = trim( cpath ) // trim( ctyp ) // '_' // trim( ctime )
      endif
      call ntdisk(nft,kio,cdsn)
!
      n6 = n6sv
!
      return
      end
