!**********************************************************************
      subroutine pledmpl
!**********************************************************************
      use cplcom, only : fcna, nion, vna, vni
      use cplmet, only : icaxs, icmpe, icwl1, jcxp1, jcxp2
      use cpmpls, only : bniedg, jmd1, lnedg
      use csize,  only : ndsp
      use cunit,  only : n6
      implicit none
!
!::local variables
      integer :: jc, ic, ia
      real*8  :: zdni
      real*8  :: zfac(ndsp)
      real*8  :: bnip
!
      write(n6,'(/2x,"*** pledmpl ***")')
!
      jc = jmd1
      ic = icmpe
      bnip = bniedg
!
      write(n6,'(2x,"lnedg =",i2,"  jc,ic =",2i5,"  bnip =",1pe11.3,
     >   "  vni(jc,ic) =",1pe11.3)') lnedg, jc, ic, bnip, vni(jc,ic)
      write(n6,'(2x,"nion =",i2)') nion
      write(n6,'(2x,"fcna =",1p10e11.3)') (fcna(ia),ia=1,nion)
!
      if( lnedg .eq. 0 ) then
      write(n6,'(2x,"no effect of sub. pledmpl due to lnedg = 0")')
      return
      endif
!
      zdni = bnip - vni(jc,ic)
      write(n6,'(2x,"zdni =",1pe11.3)') zdni
! KH20121119
!      if( zdni.lt.0.0d0 ) call wexit("pledmpl","zdni < 0")
!
      do ic = icwl1, icaxs
      do jc = jcxp1+1, jcxp2-1
      do ia = 1, nion
      zfac(ia) = vna(jc,ic,ia)/vni(jc,ic)
      enddo
      vni(jc,ic) = vni(jc,ic) + zdni
      do ia = 1, nion
      vna(jc,ic,ia) = vni(jc,ic)*zfac(ia)
      enddo
      enddo
      enddo
!
      call plmdprf(1,1)
      call plmdprf(2,1)
!
      return
      end
