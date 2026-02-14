!***********************************************************************
      subroutine plmmtr
!***********************************************************************
!
!     metric in the hot core
!
!     hvol(i)     = sum_(j)[hvol(j,i)]
!     gdsv(i)_E/W = 0.0d0
!     gdsv(i)_N/S = sum_(j)[gdsv(j,i)_N/S]
!     gwtmp(i)    = 1/2*dlr(i)
!
!------------------------------------------------------------------------
      use cplmet, only : gdsv, gwtm, gwtp, hvol, icmpe, icspx, icwl1
     >    , jcxp1, jcxp2, kce, kcn, kcs, kcw
      use cpmpls, only : arhmp, dvmp, jmd1, sfmp, wdmp, wfmp
      use cunit,  only : n6
      implicit none
!
!::local variables
      integer  jcst, jcen, ic, jc, jcm, icst, icen
      real*8   zvol, zgde, zgdw, zgdn, zgds, zgd0, zgd1, dr0, dr1
      real*8   drp, drm, zwp, zwm
!
      icst = icwl1
      icen = icmpe
      jcst = jcxp1 + 1
      jcen = jcxp2 - 1
      jcm = jmd1
!
      write(n6,'(/2x,"*** plmmtr **")')
      write(n6,'(2x,"icspx =",i3,"  icmpe =",i3,
     >  "  icst =",i3,"  icen =",i3,"  jcst =",i3,"  jcen =",i3)')
     >  icspx,icmpe,icst,icen,jcst,jcen
!
!::hvol,gdsv
      write(n6,'(2x,1x,"ic",4x,"dvl",9x,"dvl",9x,"dr0",9x,"dr1",9x
     >  ,"gdsv_E",6x,"gdsv_W",6x,"gdsv_N",6x,"gdsv_S",6x,"Sps/dr0"
     >  ,5x,"Sps/dr1")')
      do ic = icst, icen
      zvol = 0.0d0
      zgde = 0.0d0
      zgdw = 0.0d0
      zgdn = 0.0d0
      zgds = 0.0d0
      do jc = jcst, jcen
      zvol = zvol + hvol(jc,ic)
      zgde = zgde + gdsv(jc,ic,kce)
      zgdw = zgdw + gdsv(jc,ic,kcw)
      zgdn = zgdn + gdsv(jc,ic,kcn)
      zgds = zgds + gdsv(jc,ic,kcs)
      enddo
      dr0  = -(arhmp(ic+1)-arhmp(ic))
      dr1  = wfmp(ic)
      zgd0 = +sfmp(ic)/dr0
      zgd1 = +sfmp(ic)/dr1
      write(n6,'(2x,i3,1p10e12.3)') ic
     >  ,zvol, dvmp(ic), dr0, dr1
     >  ,zgde, zgdw, zgdn, zgds, zgd0, zgd1
!
      hvol(jcm,ic) = zvol
      gdsv(jcm,ic,kce) = zgde
      gdsv(jcm,ic,kcw) = zgdw
      gdsv(jcm,ic,kcn) = zgdn
      gdsv(jcm,ic,kcs) = zgds
      enddo
!
!::gwtm,gwtp
      write(n6,'(2x,1x,"ic",4x,"gwtm",8x,"gwtp",8x,"drm/dr",6x
     >   ,"drp/dr")')
      do ic = icst, icen
      drp = 0.5d0*wdmp(ic+1)
      drm = 0.5d0*wdmp(ic)
      zwp = drp/(drp+drm)
      zwm = drm/(drp+drm)
      write(n6,'(2x,i3,1p4e12.3)')
     >   ic,gwtm(jcm,ic),gwtp(jcm,ic),zwm,zwp
      gwtm(jcm,ic) = zwm
      gwtp(jcm,ic) = zwp
      enddo
!
      return
      end
