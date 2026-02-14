!***********************************************************************
      subroutine imedrfl(ip,dtstp)
!***********************************************************************
      use cimcom, only : frf, ien, ir, is, ivd, ivs, rf2_ro, rf1_c
     >    , rf1_r, rf1_z, rr, zz
      use cntcom, only : mrgn
      implicit none
!
!::argument
      integer, intent(in) :: ip
      real(8), intent(in) :: dtstp ! dummy
!
!::local variables
      integer  ic, ix
      real*8   roh
! function
      real(8)    funroh, random
!
!::ion impuirty in the main region
      if( ien(ip).ge.6 ) return
      if( is(ip).eq.0  ) return
!
      ic = ir(ip)
      if( mrgn(ic).ne.6 .and. mrgn(ic).ne.7 ) return
!
!::rho
      roh = funroh(ic,rr(ip),zz(ip))
!
!::flag of penetration
      if( ivs(ip).eq.0 .and. roh.le.0.998d0 ) ivs(ip) = 1
      if( ivd(ip).eq.0 .and. roh.le.0.900d0 ) ivd(ip) = 1
      if( roh.ge.rf2_ro ) return
!
!::new position
      ix = int(frf*random(0) + 1.0d0)
      rr(ip) = rf1_r(ix)
      zz(ip) = rf1_z(ix)
      ir(ip) = rf1_c(ix)
      ic = ir(ip)
      roh = funroh(ic,rr(ip),zz(ip))
      end
