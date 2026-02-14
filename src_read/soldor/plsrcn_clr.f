!***********************************************************************
      subroutine plsrcn_clr
!***********************************************************************
      use cntmnt, only : sn0, ssn, ssp, swe, swi
      use cpmpls, only : an0mp, asnmp
      use csize,  only : ndsp, ndx, ndy
      implicit none
!
!::local variables
      integer  nsza, nsiz
      real*8   zero
!
!::clear
      nsza = ndx*ndy*ndsp
      nsiz = ndx*ndy
      zero = 0.0d0
      call setd( sn0, nsza, zero )
      call setd( ssn, nsza, zero )
      call setd( ssp, nsza, zero )
      call setd( swi, nsiz, zero )
      call setd( swe, nsiz, zero )
!
      nsiz = ndy
      call setd( an0mp, nsiz, 0.0d0 )
      call setd( asnmp, nsiz, 0.0d0 )
!
!::plasma parameter
      call plmdprf(1,0)
!
      return
      end
