!**********************************************************************
      subroutine pycpnt
!**********************************************************************
!
!       new value at cornor points
!
!          e^psi*grad(na)|w = (1/lna)*na|w
!          e^psi*grad(Ti)|w = (1/lTi)*Ti|w
!          e^psi*grad(Te)|w = (1/lTe)*Te|w
!          va|w = 0.0
!
!       input parameter : decay length with sign
!             lna, lTi, lTe =  0.01 [m]  at sol side wall
!                           = -0.01 [m]  at prv side wall
!
!         1 cm  ==> 5 cm   2009/07/24
!----------------------------------------------------------------------
      use cphcns, only : cev
      use cplcom, only : ama, aza, c12, gcse, gcsi, nion, q1a, q2a, q3
     >    , q4, vcs, vea, vna, vne, vni, vte, vti, vva, vve, vzf
      use cplmet, only : gare, gdsv, icmax, jcmax, kce, kcn, kcs, kcw
      use csize,  only : ndsp
      implicit none
!
!::local varaiables
      integer  jc, ic, js, is, j, i, k, jp, ip, ia, jerr
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer ldbg, jctb(4), ictb(4)
      integer jctb(4), ictb(4)
      real*8  gnaw1, gtiw1, gtew1, gnaw2, gtiw2, gtew2, zrlx
      real*8  zln, zdm, zne, zq3, zni, zzf, zve, zki, zro, zpe, zpi
      real*8  zge, zgi
! modified 2/1 lines organize local variables and include files by kamata 2021/05/31
!ik   real*8  bna(4,ndsp),bva(4,ndsp),bti(4),bte(4)
!ik   data ldbg/0/
      real*8  bna(4,ndsp),bti(4),bte(4)
!
!x      write(n6,'(/2x,"***  pycpnt  ***  time =",1pe12.3)') time
!
!::input (absolute value)  1 cm ==> 5 cm
      gnaw1 = 0.05
      gtiw1 = 0.05
      gtew1 = 0.05
      gnaw2 = 0.05
      gtiw2 = 0.05
      gtew2 = 0.05
!
!::corner cell
      jctb(1) = 1
      ictb(1) = 1
      jctb(2) = 1
      ictb(2) = icmax(1)
      jctb(3) = jcmax
      ictb(3) = 1
      jctb(4) = jcmax
      ictb(4) = icmax(jcmax)
!
!----------------------------------------------------------------------
!::cornor cell (j=1,i=1)
!----------------------------------------------------------------------
! modified 8/4 lines organize local variables and include files by kamata 2021/05/31
!ik   jc = 1
!ik   ic = 1
!ik   js = jc
!ik   is = ic
!ik   j  = jc
!ik   i  = ic
!ik   jp = j + 1
!ik   ip = i + 1
      js = 1
      is = 1
      jp = js + 1
      ip = is + 1
!
      zln = gnaw1
      zdm = c12*(gdsv(js,is,kcw)+gdsv(js,is,kcs))+gare(js,is)/zln
      do 110 ia = 1, nion
      bna(1,ia)= c12/zdm*
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik  >    ( (gdsv(js,is,kce)-gdsv(js,is,kcs))*vna(jp,i ,ia)
     >    ( (gdsv(js,is,kce)-gdsv(js,is,kcs))*vna(jp,is,ia)
     >     +(gdsv(js,is,kce)+gdsv(js,is,kcn))*vna(jp,ip,ia)
! modified 2/1 lines organize local variables and include files by kamata 2021/05/31
!ik  >     -(gdsv(js,is,kcw)-gdsv(js,is,kcn))*vna(j, ip,ia) )
!ik   bva(1,ia) = 0.0d0
     >     -(gdsv(js,is,kcw)-gdsv(js,is,kcn))*vna(js,ip,ia) )
 110  continue
      zln = gtiw1
      zdm = c12*(gdsv(js,is,kcw)+gdsv(js,is,kcs))+gare(js,is)/zln
      bti(1) = c12/zdm*
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik  >    ( (gdsv(js,is,kce)-gdsv(js,is,kcs))*vti(jp,i )
     >    ( (gdsv(js,is,kce)-gdsv(js,is,kcs))*vti(jp,is)
     >     +(gdsv(js,is,kce)+gdsv(js,is,kcn))*vti(jp,ip)
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik  >     -(gdsv(js,is,kcw)-gdsv(js,is,kcn))*vti(j, ip) )
     >     -(gdsv(js,is,kcw)-gdsv(js,is,kcn))*vti(js,ip) )
      zln = gtew1
      zdm = c12*(gdsv(js,is,kcw)+gdsv(js,is,kcs))+gare(js,is)/zln
      bte(1) = c12/zdm*
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik  >    ( (gdsv(js,is,kce)-gdsv(js,is,kcs))*vte(jp,i )
     >    ( (gdsv(js,is,kce)-gdsv(js,is,kcs))*vte(jp,is)
     >     +(gdsv(js,is,kce)+gdsv(js,is,kcn))*vte(jp,ip)
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik  >     -(gdsv(js,is,kcw)-gdsv(js,is,kcn))*vte(j, ip) )
     >     -(gdsv(js,is,kcw)-gdsv(js,is,kcn))*vte(js,ip) )
!
!----------------------------------------------------------------------
!::cornor cell (j=1,i=imax)
!----------------------------------------------------------------------
! modified 8/4 lines organize local variables and include files by kamata 2021/05/31
!ik   jc = 1
!ik   ic = icmax(1)
!ik   js = jc
!ik   is = ic - 1
!ik   j  = jc
!ik   i  = ic - 1
!ik   jp = j + 1
!ik   ip = i + 1
      js = 1
      ip = icmax(1)
      jp = js + 1
      is = ip - 1
      zln = -gnaw2
      zdm = c12*(gdsv(js,is,kcw)-gdsv(js,is,kcn))+gare(js,is)/zln
      do 120 ia = 1, nion
      bna(2,ia)= c12/zdm*
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik  >    (-(gdsv(js,is,kcw)+gdsv(js,is,kcs))*vna(j, i ,ia)
!ik  >     +(gdsv(js,is,kce)-gdsv(js,is,kcs))*vna(jp,i, ia)
     >    (-(gdsv(js,is,kcw)+gdsv(js,is,kcs))*vna(js,is,ia)
     >     +(gdsv(js,is,kce)-gdsv(js,is,kcs))*vna(jp,is,ia)
     >     +(gdsv(js,is,kce)+gdsv(js,is,kcn))*vna(jp,ip,ia) )
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   bva(2,ia) = 0.0d0
 120  continue
      zln = -gtiw2
      zdm = c12*(gdsv(js,is,kcw)-gdsv(js,is,kcn))+gare(js,is)/zln
      bti(2) = c12/zdm*
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik  >    (-(gdsv(js,is,kcw)+gdsv(js,is,kcs))*vti(j, i )
!ik  >     +(gdsv(js,is,kce)-gdsv(js,is,kcs))*vti(jp,i )
     >    (-(gdsv(js,is,kcw)+gdsv(js,is,kcs))*vti(js,is)
     >     +(gdsv(js,is,kce)-gdsv(js,is,kcs))*vti(jp,is)
     >     +(gdsv(js,is,kce)+gdsv(js,is,kcn))*vti(jp,ip) )
      zln = -gtew2
      zdm = c12*(gdsv(js,is,kcw)-gdsv(js,is,kcn))+gare(js,is)/zln
      bte(2) = c12/zdm*
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik  >    (-(gdsv(js,is,kcw)+gdsv(js,is,kcs))*vte(j, i )
!ik  >     +(gdsv(js,is,kce)-gdsv(js,is,kcs))*vte(jp,i )
     >    (-(gdsv(js,is,kcw)+gdsv(js,is,kcs))*vte(js,is)
     >     +(gdsv(js,is,kce)-gdsv(js,is,kcs))*vte(jp,is)
     >     +(gdsv(js,is,kce)+gdsv(js,is,kcn))*vte(jp,ip) )
!
!----------------------------------------------------------------------
!::cornor cell (j=jmax,i=1)
!----------------------------------------------------------------------
! modified 8/4 lines organize local variables and include files by kamata 2021/05/31
!ik   jc = jcmax
!ik   ic = 1
!ik   js = jc
!ik   is = ic
!ik   j  = jc - 1
!ik   i  = ic
!ik   jp = j + 1
!ik   ip = i + 1
      js = jcmax
      is = 1
      j  = js - 1
      ip = is + 1
      zln = gnaw1
      zdm = c12*(gdsv(js,is,kce)-gdsv(js,is,kcs))-gare(js,is)/zln
      do 130 ia = 1, nion
      bna(3,ia)= c12/zdm*
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik  >    ( (gdsv(js,is,kcw)+gdsv(js,is,kcs))*vna(j, i ,ia)
!ik  >     -(gdsv(js,is,kce)+gdsv(js,is,kcn))*vna(jp,ip,ia)
     >    ( (gdsv(js,is,kcw)+gdsv(js,is,kcs))*vna(j, is,ia)
     >     -(gdsv(js,is,kce)+gdsv(js,is,kcn))*vna(js,ip,ia)
     >     +(gdsv(js,is,kcw)-gdsv(js,is,kcn))*vna(j, ip,ia) )
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   bva(3,ia) = 0.0d0
 130  continue
      zln = gtiw1
      zdm = c12*(gdsv(js,is,kce)-gdsv(js,is,kcs))-gare(js,is)/zln
      bti(3) = c12/zdm*
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik  >    ( (gdsv(js,is,kcw)+gdsv(js,is,kcs))*vti(j, i )
!ik  >     -(gdsv(js,is,kce)+gdsv(js,is,kcn))*vti(jp,ip)
     >    ( (gdsv(js,is,kcw)+gdsv(js,is,kcs))*vti(j, is )
     >     -(gdsv(js,is,kce)+gdsv(js,is,kcn))*vti(js,ip)
     >     +(gdsv(js,is,kcw)-gdsv(js,is,kcn))*vti(j, ip) )
      zln = gtew1
      zdm = c12*(gdsv(js,is,kce)-gdsv(js,is,kcs))-gare(js,is)/zln
      bte(3) = c12/zdm*
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik  >    ( (gdsv(js,is,kcw)+gdsv(js,is,kcs))*vte(j, i )
!ik  >     -(gdsv(js,is,kce)+gdsv(js,is,kcn))*vte(jp,ip)
     >    ( (gdsv(js,is,kcw)+gdsv(js,is,kcs))*vte(j, is )
     >     -(gdsv(js,is,kce)+gdsv(js,is,kcn))*vte(js,ip)
     >     +(gdsv(js,is,kcw)-gdsv(js,is,kcn))*vte(j, ip) )
!
!----------------------------------------------------------------------
!::cornor cell (j=jmax,i=imax)
!----------------------------------------------------------------------
! modified 8/4 lines organize local variables and include files by kamata 2021/05/31
!ik   jc = jcmax
!ik   ic = icmax(jcmax)
!ik   js = jc
!ik   is = ic - 1
!ik   j  = jc - 1
!ik   i  = ic - 1
!ik   jp = j + 1
!ik   ip = i + 1
      js = jcmax
      ip = icmax(jcmax)
      j  = js - 1
      is = ip - 1
      zln = -gnaw2
      zdm = c12*(gdsv(js,is,kce)+gdsv(js,is,kcn))-gare(js,is)/zln
      do 140 ia = 1, nion
      bna(4,ia)= c12/zdm*
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik  >    ( (gdsv(js,is,kcw)+gdsv(js,is,kcs))*vna(j, i ,ia)
!ik  >     -(gdsv(js,is,kce)-gdsv(js,is,kcs))*vna(jp,i, ia)
     >    ( (gdsv(js,is,kcw)+gdsv(js,is,kcs))*vna(j, is,ia)
     >     -(gdsv(js,is,kce)-gdsv(js,is,kcs))*vna(js,is,ia)
     >     +(gdsv(js,is,kcw)-gdsv(js,is,kcn))*vna(j, ip,ia) )
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   bva(4,ia) = 0.0d0
 140  continue
      zln = -gtiw2
      zdm = c12*(gdsv(js,is,kce)+gdsv(js,is,kcn))-gare(js,is)/zln
      bti(4) = c12/zdm*
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik  >    ( (gdsv(js,is,kcw)+gdsv(js,is,kcs))*vti(j, i )
!ik  >     -(gdsv(js,is,kce)-gdsv(js,is,kcs))*vti(jp,i )
     >    ( (gdsv(js,is,kcw)+gdsv(js,is,kcs))*vti(j, is )
     >     -(gdsv(js,is,kce)-gdsv(js,is,kcs))*vti(js,is )
     >     +(gdsv(js,is,kcw)-gdsv(js,is,kcn))*vti(j, ip) )
      zln = -gtew2
      zdm = c12*(gdsv(js,is,kce)+gdsv(js,is,kcn))-gare(js,is)/zln
      bte(4) = c12/zdm*
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik  >    ( (gdsv(js,is,kcw)+gdsv(js,is,kcs))*vte(j, i )
!ik  >     -(gdsv(js,is,kce)-gdsv(js,is,kcs))*vte(jp,i )
     >    ( (gdsv(js,is,kcw)+gdsv(js,is,kcs))*vte(j, is )
     >     -(gdsv(js,is,kce)-gdsv(js,is,kcs))*vte(js,is )
     >     +(gdsv(js,is,kcw)-gdsv(js,is,kcn))*vte(j, ip) )
!
!----------------------------------------------------------------------
!::relaxation
!----------------------------------------------------------------------
      do 150 k = 1, 4
      jc = jctb(k)
      ic = ictb(k)
      zrlx = 1.0d0
      do 160 ia = 1, nion
      vna(jc,ic,ia) = (1.0d0-zrlx)*vna(jc,ic,ia) + zrlx*bna(k,ia)
      vva(jc,ic,ia) = 0.0d0
 160  continue
      vti(jc,ic) = (1.0d0-zrlx)*vti(jc,ic) + zrlx*bti(k)
      vte(jc,ic) = (1.0d0-zrlx)*vte(jc,ic) + zrlx*bte(k)
 150  continue
!
!----------------------------------------------------------------------
!::convervative variables
!----------------------------------------------------------------------
      do 210 k = 1, 4
      jc  = jctb(k)
      ic  = ictb(k)
      zne = 0.0d0
      zq3 = 0.0
      do 220 ia = 1,nion
      q1a(jc,ic,ia) = ama(ia)*vna(jc,ic,ia)
      q2a(jc,ic,ia) = ama(ia)*vna(jc,ic,ia)*vva(jc,ic,ia)
      zne = zne + aza(ia)*vna(jc,ic,ia)
      zq3 = zq3 + 1.5d0*vna(jc,ic,ia)*vti(jc,ic)*cev
     >          + 0.5d0*ama(ia)*vna(jc,ic,ia)*vva(jc,ic,ia)**2
 220  continue
      vne(jc,ic) = zne
      q3(jc,ic)  = zq3
      q4(jc,ic)  = 1.5d0*vne(jc,ic)*vte(jc,ic)*cev
 210  continue
!
!-----------------------------------------------------------------------
!::plasma parameter at cell center
!-----------------------------------------------------------------------
      jerr = 0
      do 310 k = 1, 4
      j = jctb(k)
      i = ictb(k)
      zni = 0.0
      zne = 0.0
      zzf = 0.0
      zve = 0.0
      zki = 0.0
      zro = 0.0
!
      do 315 ia = 1,nion
      vna(j,i,ia) = q1a(j,i,ia)/ama(ia)
      vva(j,i,ia) = q2a(j,i,ia)/q1a(j,i,ia)
      zni = zni + vna(j,i,ia)
      zne = zne + aza(ia)*vna(j,i,ia)
      zzf = zzf + aza(ia)**2*vna(j,i,ia)
      zve = zve + aza(ia)*q2a(j,i,ia)/ama(ia)
      zki = zki + 0.5d0*q2a(j,i,ia)*vva(j,i,ia)
      zro = zro + q1a(j,i,ia)
 315  continue
!
      zpe = q4(j,i)/1.5d0
      zpi = (q3(j,i)-zki)/1.5d0
!
!::check of sound velocity
      if( zpi.le.0.0d0 .or. zpe.le.0.0d0 .or. zro.le.0.0d0 ) then
        jerr = jerr + 1
      else
        zge = gcse(i)
        zgi = gcsi(i)
        vcs(j,i) = sqrt((zge*zpe+zgi*zpi)/zro)
      endif
!
      vti(j,i) = zpi/(zni*cev)
      vve(j,i) = zve/zne
      vzf(j,i) = zzf/zne
      vte(j,i) = zpe/(zne*cev)
      vne(j,i) = zne
      vni(j,i) = zni
!
      do 325 ia = 1, nion
      vea(j,i,ia) = 1.5d0*vna(j,i,ia)*vti(j,i)*cev
     >             +0.5d0*ama(ia)*vna(j,i,ia)*vva(j,i,ia)**2
 325  continue
!
!::check of Te & Ti
      if( vte(j,i).le.0.0 .or. vti(j,i).le.0.0 ) then
        jerr = jerr + 1
      endif
!
 310  continue
!
!::nagative value
      if( jerr.gt.0 ) then
        call plauxv
      endif
!
      return
      end
!
!**********************************************************************
      subroutine pycpnt0
!**********************************************************************
!
!       new value at cornor points
!
!          na(jc,ic) = na(jc,id)
!          Ti(jc,ic) = Ti(jc,id)   (jc,ic)  corner point
!          Te(jc,ic) = Te(jc,id)   (jc,id)  divertor plate
!          vp(jc,ic) = 0.0
!
!----------------------------------------------------------------------
      use cphcns, only : cev
      use cplcom, only : ama, aza, gcse, gcsi, nion, q1a, q2a, q3, q4
     >    , vcs, vea, vna, vne, vni, vte, vti, vva, vve, vzf
      use cplmet, only : icmax, jcmax
      implicit none
!
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  jctb(4),ictb(4),jdtb(4),idtb(4), ldbg
      integer  jctb(4),ictb(4),jdtb(4),idtb(4)
      integer  k, jc, ic, jd, id, ia, jerr, j, i
      real*8   zne, zq3, zni, zzf, zve, zki, zro, zpe, zpi, zge, zgi
      real*8   facp
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   data ldbg/0/
!
!x      write(n6,'(/2x,"***  pycpnt0  ***  time =",1pe12.3)') time
!
!::corner cell
      jctb(1) = 1
      ictb(1) = 1
      jctb(2) = 1
      ictb(2) = icmax(1)
      jctb(3) = jcmax
      ictb(3) = 1
      jctb(4) = jcmax
      ictb(4) = icmax(jcmax)
!
!::d-plate
      idtb(1) = ictb(1) + 1
      idtb(2) = ictb(2) - 1
      idtb(3) = ictb(3) + 1
      idtb(4) = ictb(4) - 1
      jdtb(1) = jctb(1)
      jdtb(2) = jctb(2)
      jdtb(3) = jctb(3)
      jdtb(4) = jctb(4)
!
!::facp
      facp = 0.95d0
!
!----------------------------------------------------------------------
!::plasma parameter & convervative variables
!----------------------------------------------------------------------
      do 210 k = 1, 4
      jc  = jctb(k)
      ic  = ictb(k)
      jd  = jdtb(k)
      id  = idtb(k)
      zne = 0.0d0
      zq3 = 0.0
      vti(jc,ic) = vti(jd,id) * facp
      vte(jc,ic) = vte(jd,id) * facp
      do 220 ia = 1,nion
      vna(jc,ic,ia) = vna(jd,id,ia) * facp
      vva(jc,ic,ia) = 0.0d0
      q1a(jc,ic,ia) = ama(ia)*vna(jc,ic,ia)
      q2a(jc,ic,ia) = ama(ia)*vna(jc,ic,ia)*vva(jc,ic,ia)
      zne = zne + aza(ia)*vna(jc,ic,ia)
      zq3 = zq3 + 1.5d0*vna(jc,ic,ia)*vti(jc,ic)*cev
     >          + 0.5d0*ama(ia)*vna(jc,ic,ia)*vva(jc,ic,ia)**2
 220  continue
      vne(jc,ic) = zne
      q3(jc,ic)  = zq3
      q4(jc,ic)  = 1.5d0*vne(jc,ic)*vte(jc,ic)*cev
 210  continue
!
!-----------------------------------------------------------------------
!::plasma parameter at cell center
!-----------------------------------------------------------------------
      jerr = 0
      do 310 k = 1, 4
      j = jctb(k)
      i = ictb(k)
      zni = 0.0
      zne = 0.0
      zzf = 0.0
      zve = 0.0
      zki = 0.0
      zro = 0.0
!
      do 315 ia = 1,nion
      vna(j,i,ia) = q1a(j,i,ia)/ama(ia)
      vva(j,i,ia) = q2a(j,i,ia)/q1a(j,i,ia)
      zni = zni + vna(j,i,ia)
      zne = zne + aza(ia)*vna(j,i,ia)
      zzf = zzf + aza(ia)**2*vna(j,i,ia)
      zve = zve + aza(ia)*q2a(j,i,ia)/ama(ia)
      zki = zki + 0.5d0*q2a(j,i,ia)*vva(j,i,ia)
      zro = zro + q1a(j,i,ia)
 315  continue
!
      zpe = q4(j,i)/1.5d0
      zpi = (q3(j,i)-zki)/1.5d0
!
!::check of sound velocity
      if( zpi.le.0.0d0 .or. zpe.le.0.0d0 .or. zro.le.0.0d0 ) then
        jerr = jerr + 1
      else
        zge = gcse(i)
        zgi = gcsi(i)
        vcs(j,i) = sqrt((zge*zpe+zgi*zpi)/zro)
      endif
!
      vti(j,i) = zpi/(zni*cev)
      vve(j,i) = zve/zne
      vzf(j,i) = zzf/zne
      vte(j,i) = zpe/(zne*cev)
      vne(j,i) = zne
      vni(j,i) = zni
!
      do 325 ia = 1, nion
      vea(j,i,ia) = 1.5d0*vna(j,i,ia)*vti(j,i)*cev
     >             +0.5d0*ama(ia)*vna(j,i,ia)*vva(j,i,ia)**2
 325  continue
!
!::check of Te & Ti
      if( vte(j,i).le.0.0 .or. vti(j,i).le.0.0 ) then
        jerr = jerr + 1
      endif
!
 310  continue
!
!::nagative value
      if( jerr.gt.0 ) then
        call plauxv
      endif
!
      return
      end
