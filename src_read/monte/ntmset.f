!**********************************************************************
      subroutine ntmset
!**********************************************************************
!
!      neutral transport code
!
!         mfmax : number of ion/neutral flux
!         mcflx : name of flux
!         mnknd : 1(surface)/2(puff)/3(volume)
!         mnsmp : number of sample particle
!         mnvwk : index of vwork
!
!         mnvwk  flux  ==> source
!         vkflx  flux  <== source
!
!----------------------------------------------------------------------
      use cntcom, only : lemwl, nptldp, nptlvl, nptlwl
      use cntmnt, only : mcflx, mfmax, mnknd, mnsmp, mnvwk, pfl_ion
     >    , pfl_ntl, vkflx, vnsrc
      use csize,  only : ndbr, ndmsr, ndwp
      use cunit,  only : n6
      implicit none
!
!::local variables
      integer  ii, i, iflx
      real*8   tfion, tfntl
!
      mfmax = 8
      mcflx(1) = "sol"; mnknd(1) = 1;  mnsmp(1) = nptlwl
      mcflx(2) = "idp"; mnknd(2) = 1;  mnsmp(2) = nptldp
      mcflx(3) = "prv"; mnknd(3) = 1;  mnsmp(3) = nptlwl
      mcflx(4) = "odp"; mnknd(4) = 1;  mnsmp(4) = nptldp
      mcflx(5) = "man"; mnknd(5) = 1;  mnsmp(5) = 0
      mcflx(6) = "puf"; mnknd(6) = 2;  mnsmp(6) = nptlwl
      mcflx(7) = "vol"; mnknd(7) = 3;  mnsmp(7) = nptlvl
      mcflx(8) = "vli"; mnknd(8) = 3;  mnsmp(8) = 0
!
      if(lemwl.eq.3)then
        mnsmp(1)=0 ! sol
        mnsmp(3)=0 ! prv
      endif
!
!::flux
      ii = 0
      do iflx = 1, mfmax
!-----
      call test_ntnflx(iflx,tfion,tfntl)
!-----
      if( mnsmp(iflx).le.0 ) tfntl = 0.0d0
      pfl_ion(iflx) = tfion
      pfl_ntl(iflx) = tfntl
      mnvwk(iflx)  = 0
      if( tfntl.le.0.0d0 ) cycle
      ii = ii + 1
      mnvwk(iflx)  = ii
      vkflx(ii) = iflx
      enddo
      vnsrc = ii
!
!::recomb  H+ ==> H0
      mnvwk(8) = mnvwk(7)
!
!::debug write
      write(n6,'(/2x,"*** ntmset ***  dummy flux")')
      write(n6,'(4x,"ndwp =",i4,"  ndbr =",i5)') ndwp,ndbr
      write(n6,'(4x,"i",2x,"flx",3x,"knd",4x,"nsmp",2x,"iv",5x,
     >   "fion",11x,"fntl",12x,"dif")')
      do i = 1, mfmax
      write(n6,'(2x,i3,2x,a,i3,i8,i4,1p2e15.6,1pe11.2)')
     >  i, mcflx(i), mnknd(i), mnsmp(i), mnvwk(i),
     >  pfl_ion(i), pfl_ntl(i), pfl_ion(i)-pfl_ntl(i)
      enddo
      write(n6,'(2x,"vnsrc =",i2,"  ndmsr =",i2)') vnsrc,ndmsr
!
      if( vnsrc.gt.ndmsr .or. vnsrc.le.0 ) then
      call wexit("ntmset","dimension error  vnsrc > ndmsr")
      endif
      if( ndwp.gt.ndbr ) then
      call wexit("ntmset","dimension error  ndwp > ndbr")
      endif
!
      return
      end
!
!**********************************************************************
      subroutine test_ntnflx(kflx,tfion,tfntl)
!**********************************************************************
!
!       dummy flux
!
!----------------------------------------------------------------------
      use cntcom, only : npf, pfmag
      implicit none
!
!::argument
      integer, intent(in)  :: kflx
      real(8), intent(out) :: tfion, tfntl
!
!   kflx       : 1:(sol)/2:(idp)/3:(prv)/4:(odp)/5:(man)
!
      if( kflx.eq.1 ) then
        tfion = 1.0d22
        tfntl = tfion
      elseif( kflx.eq.2 ) then
        tfion = 2.0d24
        tfntl = tfion
      elseif( kflx.eq.3 ) then
        tfion = 3.0d22
        tfntl = tfion
      elseif( kflx.eq.4 ) then
        tfion = 4.0d24
        tfntl = tfion
      elseif( kflx.eq.5 ) then
        tfion = 5.0d22
        tfntl = 0.0d0
      elseif( kflx.eq.6 ) then
        tfion = 0.0d0
        if( npf.gt.0 ) tfntl = pfmag(1)
      elseif( kflx.eq.7 ) then
        tfion = -7.0d22
        tfntl = -tfion
      elseif( kflx.eq.8 ) then
        tfion = -7.0d22
        tfntl = 0.0d0
      else
        tfion = 0.0d0
        tfntl = 0.0d0
      endif
!
      return
      end
