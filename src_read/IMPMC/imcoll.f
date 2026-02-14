!**********************************************************************
      subroutine slwtime(j,i,iz,tauZe,tauZi,kout)
!**********************************************************************
!
!       C-log ramda = 16.0d0  in IMPMC
!       C-log ramda = 12.0d0  in soldor
!
!----------------------------------------------------------------------
      use cimcom, only : aimas, azmas, slnv, slw0
      use cntcom, only : mcel
      use cntpls, only : dene, teme, temi
      use cphcns, only : cev, cme, cmp, cpi
      use cplcom, only : vne, vni, vte, vti
      use cunit,  only : n6
      implicit none
!
!::argument
! modified 2/2 lines organize local variables and include files by kamata 2021/06/28
!ik   integer  j, i, iz, kout
!ik   real*8   tauZe, tauZi
      integer, intent(in)  :: j, i, iz, kout
      real(8), intent(out) :: tauZe, tauZi
!
!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/06/28
!ik   integer ic, imox, imoy
      integer ic
      real*8  hev, hme, hmp, hee, hrmd, hmi, hmz, hcn0, hcne, hcni
      real*8  zi, zne, zni, zte, zti, cftz, tauz
!
      integer lp
      data    lp/0/
      save
! added 2 lines organize local variables and include files by kamata 2021/06/28
! function
      integer    imox, imoy
!
!::physical const MKS ==> CGS
      if( kout.eq.-1 ) then
      lp = 0
      hev  = cev*1.0d7
      hme  = cme*1.0d3
      hmp  = cmp*1.0d3
      hee  = 4.80298d-10
!
      hrmd = 16.0d0
      hmi  = aimas*hmp
      hmz  = azmas*hmp
      hcn0 = 3.0d0/4.0d0/sqrt(2.0d0*cpi)
      hcne = hcn0*hev**1.5*hmz/sqrt(hme)/hee**4/hrmd
      hcni = hcn0*hev**1.5*hmz*hmz/(hmi+hmz)/sqrt(hmi)/hee**4/hrmd
!
!::debg write
      write(n6,'(/2x,"*** slwtime ***  (CGS_unit)")')
      write(n6,'(2x,"heV  =",1pe12.4,"  hme  =",1pe12.4,"  hmp  =",
     > 1pe12.4,"  hee  =",1pe12.4)') hev, hme, hmp, hee
      write(n6,'(2x,"hrmd =",1pe12.4,"  hmi  =",1pe12.4,"  hmz  =",
     >  1pe12.4,"  hcne =",1pe12.4,"  hcni =",1pe12.4)')
     >  hrmd, hmi, hmz, hcne, hcni
      write(n6,'(/2x,4x,"lp",4x,"j",3x,"i",3x,"ic",6x,"ix",3x,"iy",3x,
     >  "Ne",10x,"Ni",10x,"Te",10x,"Ti",10x,"dene",8x,"teme",8x,
     >  "temi",8x,"tauZe",7x,"tauZi",7x,"tauz")')
      return
      endif
!
!::Ne,Ni MKS ==> CGS
      zne = vne(j,i)*1.0d-6
      zni = vni(j,i)*1.0d-6
      zte = vte(j,i)
      zti = vti(j,i)
!
!::including charge state
      zi = iz
      tauZe = hcne*zte*sqrt(zte)/zne/zi**2
      tauZi = hcni*zti*sqrt(zti)/zni/zi**2
!
      if( kout.ne.0 ) then
      tauz = 0.0d0
      ic = mcel(j,i)
      cftz = 1.0d0/dsqrt(aimas)*azmas/(azmas+aimas)
      if(ic.gt.0 )
     >   tauz = 1.0d0/(zi**2*slw0(ic,2)*slnv(ic,2))*cftz
!
      write(n6,'(2x,i6,i5,i4,i7,i6,i5,1p10e12.3)')
     >  lp, j, i, ic, imox(ic), imoy(ic),
     >  zne, zni, zte, zti, dene(ic), teme(ic), temi(ic),
     >  tauZe, tauZi, tauz
      endif
!
      return
      end
!
!**********************************************************************
      subroutine mfpath(j,i,rmfpe,rmfpi)
!**********************************************************************
!
!  rmfpe = vthe*taue
!     = sqrt(Te*cev/ame)*3.4411e11/Lrmd/zi**2*Te**1.5/Ni
!
!  rmfpi = vthi*taui
!     = sqrt(Ti*cev/ami)*2.0853e13*sqrt(ami/amp)/Lrmd/zi**4*Ti**1.5/Ni
!
!----------------------------------------------------------------------
      use cimcom, only : aimas
      use cphcns, only : cev, cme, cmp
      use cplcom, only : vni, vte, vti
      implicit none
!
!::argument
! modified 2/2 lines organize local variables and include files by kamata 2021/06/28
!ik   integer  j, i
!ik   real*8   rmfpe, rmfpi
      integer, intent(in)  :: j, i
      real(8), intent(out) :: rmfpe, rmfpi
!
!::local variables
      real*8  hrmd, hzi, hzi2, hzi4, hcne, hcni
!
      hrmd = 12.0d0
      hzi  = 1.0d0   ! D+ ion
      hzi2 = hzi**2
      hzi4 = hzi**4
!
      hcne = 3.4411d11/hrmd/hzi2
      hcni = 2.0853d13/hrmd/hzi4*sqrt(aimas)
!
      hcne = hcne*sqrt(cev/cme)
      hcni = hcni*sqrt(cev/(aimas*cmp))
!
      rmfpe = hcne*vte(j,i)**2/vni(j,i)
      rmfpi = hcni*vti(j,i)**2/vni(j,i)
!
      return
      end
