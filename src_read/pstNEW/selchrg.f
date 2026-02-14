!***********************************************************************
      subroutine selchrg(cnam,chrg,nmc,var,kcn,nty)
!***********************************************************************
!
!     impmc     nmc = ncmax or ncmax2
!               Nz   :  tdnz(0:ndis2,ndmc)   ==> var
!        
!-----------------------------------------------------------------------
      use csize
      use cimcom
      use cimden
      use cplimp
      use cplimp_plimp, only : tvlzL,tionZL,trecZL
      implicit none
!
!::argument
      character(10)::cnam
      character(*) :: chrg
      integer :: nmc, nty
      real(8), dimension(ndmc) :: var
      integer :: kcn
!
!::local variables
      character(80) :: cond
      integer, parameter :: ndim = 100
      integer, dimension(ndim) :: kval
      integer :: nval
!
      integer :: mji, ic, k, iz
      real(8) :: zsm
!
      if( len_trim(chrg).le.0 ) then
      call gdget(">> enetr is (0/a:all/i:ion/2-5) ==> ",cond,mji)
      chrg = cond
      else
      cond = chrg
      endif
!
      if( cond(1:1).eq."a" ) then
        write(cond,'(i1,"-",i2)') 0, ismaxL(nty)
      elseif( cond(1:1).eq."i" ) then
        write(cond,'(i1,"-",i2)') 1, ismaxL(nty)
      endif
!
      call interp(cond,nval,kval,ndim)
      if( nval.le.0 ) then
        kcn = 1
        return
      endif

      chrg = trim(cond)
      if(trim(cnam(1:2)).eq."Nz")then
      do ic = 1, nmc
      zsm = 0.0d0
      do k = 1, nval
      iz = kval(k)
      zsm = zsm + tdnzL(iz,ic,nty)
      enddo
      var(ic) = zsm
      enddo

      else if(trim(cnam(1:2)).eq."Vz")then
      do ic = 1, nmc
      zsm = 0.0d0
      do k = 1, nval
      iz = kval(k)
      zsm = zsm + tvlzL(iz,ic,nty)
      enddo
      var(ic) = zsm
      enddo

      else if(trim(cnam(1:3)).eq."ion")then
      do ic = 1, nmc
      zsm = 0.0d0
      do k = 1, nval
      iz = kval(k)
      zsm = zsm + tionZL(iz,ic,nty)
      enddo
      var(ic) = zsm
      enddo

      else if(trim(cnam(1:3)).eq."rec")then
      do ic = 1, nmc
      zsm = 0.0d0
      do k = 1, nval
      iz = kval(k)
      zsm = zsm + trecZL(iz,ic,nty)
      enddo
      var(ic) = zsm
      enddo
      end if

      kcn = 0
!
      return
      end
