!**********************************************************************
      subroutine imp_pls
!**********************************************************************
      use cimcom, only : cdif, cdifrg, dfcf
     > ,lerp_points_cdif, lerp_r_cdif, lerp_diff_cdif
      use cimden, only : twne
      use cntcom, only : mrgn, ncmax2, iply
      use cntmnt, only : mcflx, mfmax, pfl_ion, pfl_ntl, visrc
      use cntpls, only : dene
      use csize,  only : ndmc
      use cunit,  only : n6
      use mod_externalgrid, only : use_exdata
      use cimctl, only : nprg
      implicit none
!
!::local variables
      integer iflx, isrc
      integer :: ic, irg
      real(8), allocatable :: dfcf_r(:)
!
      write(n6,'(30("-"),"  imp_pls(sta)  twne = dene")')

!::plasma parameter in monte code   2015/04/27
      call mcplas
!
!::radiation in soldor  Ne*Nz*Lz(Te)  Ne:twne Nz:tdnz
      twne(1:ndmc) = dene(1:ndmc)
!
!::set diffusion coeffientdistribution
      call set_dfcf_dist(dfcf_r, lerp_points_cdif(nprg)
     > , lerp_diff_cdif(1:100,nprg), lerp_r_cdif(1:100,nprg))
!
!::ion flux (pfl_ion,pfl_ntl)   see sub. plpflx
      call set_nflx(1)
      write(n6,'(/2x,"*** flux ***  total neutral density")')
      do iflx = 1, mfmax
      isrc = visrc(iflx)
      write(n6,'(i3,2x,i3,2x,a,2x,1p2e14.6)')
     >  iflx, isrc, mcflx(iflx), pfl_ion(iflx), pfl_ntl(iflx)
      enddo
!
!::diffusion coefficient after set plasma
      write(n6,'(2x,"dfcf = 0.0 where vac and void in prv")')
      do ic = 1, ncmax2
        irg = mrgn(ic)
        !::set diffusion coefficient
        if( (irg.eq.4 .or. irg.eq.5) .and. dene(ic).eq.0.0d0 ) then ! private region
          dfcf(:,ic) = 0.0d0
        elseif(irg.ge.8 .and. .not.use_exdata) then ! vacume region
          dfcf(:,ic) = 0.0d0
        elseif((lerp_points_cdif(nprg).ne.0)  !::radial distribution (except for private region)
     >      .and. (irg.ne.4) .and. (irg.ne.5) .and. (irg.lt.8)
     >      .and. iply(ic).ne.0) then
          dfcf(:,ic) = dfcf_r(iply(ic))
        elseif(cdifrg(irg).lt.0.0d0) then ! using isotropic value
          dfcf(:,ic) = cdif
        else
          dfcf(:,ic) = cdifrg(irg) ! change value for each region
        endif
      enddo
!
!
!::various table
      call set_scatt
      call set_ratez
!
!::set variables which depend on plasma & neutral
      call set_forc
      call set_gcoef
      write(n6,'(30("-"),"  imp_pls(end)"/)')
      call flush(n6)
!
      return

!**********************************************************************
      contains
!**********************************************************************
      subroutine set_dfcf_dist(dfcf_r,points,value_r,r_dist)
!**********************************************************************
      use cpmpls, only : jmd1
      use cplmet, only : icwl1, icaxs
      use csize,  only : ndy
      implicit none
! arguments
      integer,intent(in) :: points
      real(8),intent(in) :: value_r(100), r_dist(100)
      real(8),allocatable,intent(out) :: dfcf_r(:)
! local
      real(8) :: rhf_omd(ndy)
      real(8) :: rst(ndy), ren(ndy)
      integer n_omd

! calculate rhf_omd for vldar
      call codrad(jmd1,rst,ren,rhf_omd,n_omd,ndy,icwl1,icaxs)
      allocate(dfcf_r(n_omd))
      call interpolate_array(dfcf_r,points,r_dist(1:points)
     > ,value_r(1:points),rhf_omd(1:n_omd),n_omd)

      end subroutine set_dfcf_dist

!**********************************************************************
!:: debug write for check dicf distribution
      subroutine debug_dfcf(jc)
!**********************************************************************
      use cimcom, only : dfcf
      use cpmpls, only : jmd1
      use cplmet, only : icmax, icmin, jcmax, icwl1, icaxs, icwl2
     > ,jcxp1,jcxp2
      use cntcom, only : mcel,xpnt,ypnt, mgrd
      use cunit,  only : lmspe, lmype
      use csize,  only : ndy
      use cimctl, only : nprg
      implicit none

      integer jc, ic1, ic2, i6, ic_MC

      real(8) :: rhf_omd(ndy)
      real(8) :: rst(ndy), ren(ndy)
      integer n_omd
      character(80) :: cdsn

      if( lmype.ne.lmspe ) return

! calculate rhf_omd
      ic1 = icwl1      
      if(jcxp1 < jc .and. jc < jcxp2) then
        ic2 = icaxs
      else
        ic2 = icwl2
      endif
      i6 = 21

      call codrad(jc,rst,ren,rhf_omd,n_omd,ndy,ic1,ic2)

      write(cdsn,'(a,i3.3,a,i2.2,a)') "test",jc,"_",nprg,".txt"
      open(unit=i6, file=cdsn)
      write(i6,'(2x,"*** radial profile ***  jc =",i3)') jc
      write(i6,'(5x,"i",3x,"r_omd",6x,"Da",9x,"ic_MC")')

      do ic = ic1, ic2
        ic_MC = mcel(jc,ic)
        if(ic_MC.le.0) cycle
        write(i6,'(1x,i5,1p2e11.3,i5,1p8e11.3)') ic
     >    ,rhf_omd(ic)*100.0d0
     >    ,dfcf(0,ic_MC)
     >    ,ic_MC
     >    ,xpnt(mgrd(ic_MC,1)),xpnt(mgrd(ic_MC,2))
     >    ,xpnt(mgrd(ic_MC,3)),xpnt(mgrd(ic_MC,4))
     >    ,ypnt(mgrd(ic_MC,1)),ypnt(mgrd(ic_MC,2))
     >    ,ypnt(mgrd(ic_MC,3)),ypnt(mgrd(ic_MC,4))
      enddo
      close(i6)
      end subroutine debug_dfcf

      end subroutine imp_pls
