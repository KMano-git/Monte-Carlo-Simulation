!**********************************************************************
      subroutine plcler
!**********************************************************************
      use cmeffz, only : dzan, vdnz, vsez, vsez0, vsezg, vsezz, xzflz
     >    , xzwtm, xzwtp
      use cplcom, only : dq1a, dq2a, dq3, dq4, q1a, q2a, q3, q4, vcs
     >    , vea, vna, vne, vnezef, vni, vte, vti, vva, vve, vzf, wq1a
     >    , wq2a, wq3, wq4
      use cplimp, only : ismaxl, wmc_nty
      use cplmet, only : gdsv, gwtm, gwtp, hdsp, hdsv, hdxm, hdxp, hgdx
     >    , hgdy, hvol, hwtm, hwtp
      use csize,  only : ndmis, ndsp, ndx, ndy, nzsmx
      use cunit,  only : n6
      implicit none
!
!::local variables
      integer  nsiz, nsz4, nsza
      real*8   zero
!
      write(n6,'(/2x,"*** plcler ***")')

!:::zero
      zero = 0.0d0
!xx   zero = 1.0d40
      write(n6,'(5x,"zero =",1pe11.3)') zero
!
!::metric
      write(n6,'(5x,"/cmetrc/: metric quantities")')
      nsiz = ndx*ndy
      nsz4 = ndx*ndy*4
      call setd( hvol, nsiz, zero )
      call setd( hdsp, nsz4, zero )
      call setd( hdsv, nsz4, zero )
      call setd( gdsv, nsz4, zero )
      call setd( hwtp, nsiz, zero )
      call setd( hwtm, nsiz, zero )
      call setd( gwtp, nsiz, zero )
      call setd( gwtm, nsiz, zero )
!
      write(n6,'(5x,"/cmgpnt/: cell center")')
      nsiz = ndx*ndy
      call setd( hgdx, nsiz, zero )
      call setd( hgdy, nsiz, zero )
      call setd( hdxp, nsiz, zero )
      call setd( hdxm, nsiz, zero )
!
!::conservative variables
      write(n6,'(5x,"/cmqnow/,/cmqprv/,/cmqdlt/: q(N+1),q(N),dq(N+1)")')
      nsza = ndx*ndy*ndsp
      nsiz = ndx*ndy
      call setd(  q1a, nsza, zero )
      call setd(  q2a, nsza, zero )
      call setd(  q3,  nsiz, zero )
      call setd(  q4,  nsiz, zero )
      call setd( wq1a, nsza, zero )
      call setd( wq2a, nsza, zero )
      call setd( wq3,  nsiz, zero )
      call setd( wq4,  nsiz, zero )
      call setd( dq1a, nsza, zero )
      call setd( dq2a, nsza, zero )
      call setd( dq3,  nsiz, zero )
      call setd( dq4,  nsiz, zero )
!
!::aux. parameter
      write(n6,'(5x,"/cmauxv/: aux.parameter")')
      nsza = ndx*ndy*ndsp
      nsiz = ndx*ndy
      call setd( vna, nsza, zero )
      call setd( vne, nsiz, zero )
      call setd( vni, nsiz, zero )
      call setd( vnezef, nsiz, zero )  ! 2012/11/24
      call setd( vzf, nsiz, zero )
      call setd( vva, nsza, zero )
      call setd( vve, nsiz, zero )
      call setd( vti, nsiz, zero )
      call setd( vte, nsiz, zero )
      call setd( vcs, nsiz, zero )
      call setd( vea, nsza, zero )
!
!::impurity effect
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/05/08
!ik   vdnz (1:ndx,1:ndy,0:ndis2L, 1:nzsmx) = 0.0d0
      vdnz (1:ndx,1:ndy,0:ndmis, 1:nzsmx) = 0.0d0
      dzan (1:ndx,1:ndy) = 0.0d0
      vsez0(1:ndx,1:ndy) = 0.0d0
      vsez (1:ndx,1:ndy) = 0.0d0
      vsezz(1:ndx,1:ndy) = 0.0d0
      vsezg(1:ndx,1:ndy) = 0.0d0
      xzwtm(1:ndx,1:ndy) = 0.0d0
      xzwtp(1:ndx,1:ndy) = 0.0d0
      xzflz(1:ndx,1:ndy) = 0.0d0
!
      wmc_nty = 0
      ismaxL(1:nzsmx) = -1
!
      return
      end
