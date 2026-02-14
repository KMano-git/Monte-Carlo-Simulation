!***********************************************************************
      subroutine selchrg_avg(cnam,chrg,nmc,var,kcn,nty)
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
      use cplimp_plimp, only : tvlzL, tfrzL,tthzL
      implicit none

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
      real(8) :: cnd
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

      if(trim(cnam(1:2)).eq."Vz") then
      do ic = 1, nmc
      zsm = 0.0d0
      cnd = 0.0d0
      do k = 1, nval
      iz = kval(k)
      zsm = zsm + tvlzL(iz,ic,nty)
      if(tvlzL(iz,ic,nty).ne.0.0d0) cnd=cnd+1.0d0
      enddo
      if(cnd == 0.0d0 .and. zsm == 0.0d0) then
         var(ic) = 0.0d0
      else
         var(ic) = zsm/cnd
      endif
      enddo

      elseif(trim(cnam(1:3)).eq."Ffr") then
      do ic = 1, nmc
      zsm = 0.0d0
      cnd = 0.0d0
      do k = 1, nval
      iz = kval(k)
      zsm = zsm + tfrzL(iz,ic,nty)
      if(tfrzL(iz,ic,nty).ne.0.0d0) cnd=cnd+1.0d0
      enddo
      if(cnd == 0.0d0 .and. zsm == 0.0d0) then
         var(ic) = 0.0d0
      else
         var(ic) = zsm/cnd
      endif
      enddo

      elseif(trim(cnam(1:4)).eq."Ffrt") then
      do ic = 1, nmc
      zsm = 0.0d0
      cnd = 0.0d0
      do k = 1, nval
      iz = kval(k)
      zsm = zsm + tfrzL(iz,ic,nty)/(iz**2)
      if(tfrzL(iz,ic,nty).ne.0.0d0) cnd=cnd+1.0d0
      enddo
      if(cnd == 0.0d0 .and. zsm == 0.0d0) then
         var(ic) = 0.0d0
      else
         var(ic) = zsm/cnd
      endif
      enddo

      elseif(trim(cnam(1:3)).eq."Fth") then
      do ic = 1, nmc
      zsm = 0.0d0
      cnd = 0.0d0
      do k = 1, nval
      iz = kval(k)
      zsm = zsm + tthzL(iz,ic,nty)
      if(tthzL(iz,ic,nty).ne.0.0d0) cnd=cnd+1.0d0
      enddo
      if(cnd == 0.0d0 .and. zsm == 0.0d0) then
         var(ic) = 0.0d0
      else
         var(ic) = zsm/cnd
      endif
      enddo

      elseif(trim(cnam(1:4)).eq."Ftht") then
      do ic = 1, nmc
      zsm = 0.0d0
      cnd = 0.0d0
      do k = 1, nval
      iz = kval(k)
      zsm = zsm + tthzL(iz,ic,nty)/(iz**2)
      if(tthzL(iz,ic,nty).ne.0.0d0) cnd=cnd+1.0d0
      enddo
      if(cnd == 0.0d0 .and. zsm == 0.0d0) then
         var(ic) = 0.0d0
      else
         var(ic) = zsm/cnd
      endif
      enddo

      elseif(trim(cnam(1:3)).eq."FBl") then
      do ic = 1, nmc
      zsm = 0.0d0
      cnd = 0.0d0
      do k = 1, nval
      iz = kval(k)
      zsm = zsm + tfrzL(iz,ic,nty) + tthzL(iz,ic,nty)
      if(tfrzL(iz,ic,nty).ne.0.0d0) cnd=cnd+1.0d0
      enddo
      if(cnd == 0.0d0 .and. zsm == 0.0d0) then
         var(ic) = 0.0d0
      else
         var(ic) = zsm/cnd
      endif
      var(ic) = dsign(1.0d0,var(ic))
      enddo

      elseif(trim(cnam(1:4)).eq."Fnet") then
      do ic = 1, nmc
      zsm = 0.0d0
      cnd = 0.0d0
      do k = 1, nval
      iz = kval(k)
      zsm = zsm + tfrzL(iz,ic,nty)/(iz**2) + tthzL(iz,ic,nty)/(iz**2)
      if(tfrzL(iz,ic,nty).ne.0.0d0) cnd=cnd+1.0d0
      enddo
      if(cnd == 0.0d0 .and. zsm == 0.0d0) then
         var(ic) = 0.0d0
      else
         var(ic) = zsm/cnd
      endif
      enddo

      endif
      kcn = 0
!
      return
      end
