!***********************************************************************
      subroutine plprof(kt)
!***********************************************************************
      use cntcom, only : mcel
      use cntmnt, only : sn0
      use cntpls, only : dene, teme, temi, vflw
      use cntsrc, only : tden0, tdeng, teng0, tengg, tssn, tssp, tswe
     >    , tswi !  local common (plntsr_plt.f)
      use cphcns, only : cev
      use cplcom, only : ama, q2a, q3, q4, ssnc, ssnv, sspc, sspv, swec
     >    , swev, swic, swiv, vcs, vna, vne, vni, vte, vti, vva, wrad
      use cplmet, only : icel, jcel, jtmax, jtmin
      use csize,  only : ndx
      use cunit,  only : lmspe, lmype
      implicit none
!
!KH111108   edit format  add slen
!
!::argument
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  kt
      integer, intent(in) :: kt
!
!::local variables
      character  dsn*80
      integer    ift, ia, it, jtst, jten, jt, j, i, ig, ic
      real*8     alnx(ndx), zlmax, zmch, zpcn
      real*8     slnx(ndx), zlmax2
!
      if( lmype.ne.lmspe ) return
!
      ift = 21
      ia = 1
      ig = 1
      it = kt
!
!xx   write(n6,'(2x,"*** plprof ***   it =",i3)') it
!
!::plasma
      write(dsn,'("mflpls_",i2.2)') it
      open(unit=ift,file=dsn)
!
      jtst = jtmin(it)
      jten = jtmax(it)
      call plenc(it,alnx)
      call slenc(it,slnx)
      zlmax = alnx(jten)
      zlmax2 = slnx(jten)
!
      write(ift,'(1x,3x,"jt",4x,"j",4x,"pl",7x,"pl2",7x,
     >  "sl",6x,"sl2",4x,
     >  "Ne",9x,"Ti",9x,"Te",9x,"Va",9x,"Mach",7x,"pcon",7x
     >  "N0",9x,"Sn",9x,"Sp",9x,"Wi",9x,"We",9x,"Wrad")')
      do jt = jtst, jten
      j = jcel(jt,it)
      i = icel(jt,it)
      zmch = vva(j,i,ia)/vcs(j,i)
      zpcn = vne(j,i)*vte(j,i)+vni(j,i)*vti(j,i)
     >         +ama(ia)*vni(j,i)*vva(j,i,ia)**2/cev
      zpcn = zpcn*cev
      ic = mcel(j,i)
      write(ift,'(1x,2i5,0p2f9.5,0p2f9.4,1p15e11.3)')
     >  jt, j, alnx(jt), zlmax-alnx(jt),
     >  slnx(jt), zlmax2-slnx(jt),
     >  vne(j,i), vti(j,i), vte(j,i), vva(j,i,ia), zmch, zpcn,
     >  sn0(j,i,ia),
     >  ssnc(j,i,ia)+ssnv(j,i,ia)*vna(j,i,ia),
     >  (sspc(j,i,ia)+sspv(j,i,ia)*q2a(j,i,ia)),
     >  (swic(j,i)+swiv(j,i)*q3(j,i)),
     >  (swec(j,i)+swev(j,i)*q4(j,i)),
     >  wrad(j,i)
      enddo
!
      close(ift)
!      return
!
!::neutral
      write(dsn,'("mflntl_",i2.2)') it
      open(unit=ift,file=dsn)
!
      write(ift,'(1x,3x,"jt",4x,"j",4x,"pl",8x,"pl2",8x,
     >  "sl",8x,"sl2",8x,
     >   "Ne",9x,"Ti",9x,"Te",9x,"Vf",9x,"N0",9x,"ND2",8x,"ND2+",7x,
     >   "T0",9x,"Tg",9x,
     >   "Sn",9x,"Sp",9x,"Wi",9x,"We",9x,"Wrad")')
      do jt = jtst, jten
      j = jcel(jt,it)
      i = icel(jt,it)
      ic = mcel(j,i)
      write(ift,'(1x,2i5,0p4f11.5,1p15e11.3)')
     >  jt, j, alnx(jt), zlmax-alnx(jt),
     >  slnx(jt), zlmax2-slnx(jt),
     >  dene(ic), temi(ic), teme(ic), vflw(ic,ia),
     >  tden0(ic,ig), tdeng(ic,1), tdeng(ic,2),
     >  teng0(ic,ig), tengg(ic,1),
     >  tssn(ic,ia), tssp(ic,ia), tswi(ic), tswe(ic),
     >  wrad(j,i)
      enddo
!
      close(ift)
!
      return
      end
