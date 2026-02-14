!***********************************************************************
      subroutine impvar(ic)
!***********************************************************************
      use csize
      use cimcom
      use cimden
      use cntpls
      use cplcom
      use cplimp
      use com_impvar, only : zni, zne, znz, zav
     > , zef, zrt, znzi, znzt
      implicit none

      integer  ic

!::local variables
! modified 1/1 lines treat 4 or more impurities with IMPMC
!      real(8) :: hni(ndsp), hnz(0:ndis2L,1:nzsmx)
      real(8) :: hni(ndsp), hnz(0:ndmis,1:nzsmx)
      integer :: ia, iz, nty
!
      hni(1:nion) = deni(ic,1:nion)
      do nty = 1, wmc_nty
      hnz(0:ismaxL(nty),nty) = tdnzL(0:ismaxL(nty),ic,nty)
      enddo    
 
      zni = 0.0d0
      zne = 0.0d0
      zef = 0.0d0
      do ia = 1, nion
      zni = zni + hni(ia)
      zne = zne + aza(ia)*hni(ia)
      zef = zef + aza(ia)**2*hni(ia)
      enddo

      znz = 0.0d0
      znzi = 0.0d0
      zav = 0.0d0
      do nty = 1, wmc_nty
      do iz = 0, ismaxL(nty)
        zne = zne + dble(iz)*hnz(iz,nty)
        zef = zef + dble(iz)**2*hnz(iz,nty)
        znz = znz + hnz(iz,nty)
        zav = zav + dble(iz)*hnz(iz,nty)
        if( iz > 0 ) znzi = znzi + hnz(iz,nty)
      enddo
      enddo
      znzt = znz

      if( znz > 0.0d0 ) then
        zav = zav/znz
      else
        zav = 0.0d0
      endif
      if( zne > 0.0d0 ) then
        zef = zef/zne
        zrt = znz/zne 
      else
        zef = 0.0d0
        zrt = 0.0d0
      endif

!::no plasma but impurity exists at void
      if( zni == 0.0d0 ) then
        zef = 0.0d0
        zav = 0.0d0
        zne = 0.0d0
       endif

      return
      end
