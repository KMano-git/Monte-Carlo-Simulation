!***********************************************************************
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   subroutine plqcon_pcn(nt)
      subroutine plqcon_pcn(it)
!***********************************************************************
!
!    pressure constant
!
!-----------------------------------------------------------------------
      use cphcns, only : cev
      use cplcom, only : ama, aza, nlp, vna, vne, vte, vti, vva
      use cplmet, only : hdxm, hdxp, hpit, hvol, hvsb, hwtm, hwtp, icel
     >    , jcel, jtmax
      use cplqcn, only : qfx_cd, qfx_cv, qfx_df, qfx_vh, qfy_cd, qfy_df
     >    , qfy_vh, qvl_cl, qvl_dt, qvl_pe, qvl_pi, qvl_sc
      use csize,  only : nqfx => ndx
      use csonic, only : itim, time
      use cunit,  only : cdrout
      implicit none
!
!::argument
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  nt
      integer, intent(in) :: it
!
!::local variables
! modified 3/2 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  it, i, im, jmax, jmax1, ia, mj, k, mm
!ik   integer  ift, jw, j, jm, lenx, jp
!ik   integer  m1a, m2a, m3, m4
      integer  i, im, jmax, jmax1, ia, mj, k, mm
      integer  ift, jw, j, jm, jp
      character  ctim*5, ctub*2, ckey*3, fnm*80
      real*8   zvx_cv, zvx_df, zvy_df, zvx_cd, zvy_cd, zvl_pe
      real*8   zvl_sc, zvl_pi, zvx_vh, zvy_vh, zvl_dt, zvl_cl, zsum
!
      real*8   zfpa, zlb, zpe, zpep, zpem, zf22, zf22p
      real*8   zpehp, zpehm, zhvsb, zptot
      real*8   pfx_cv(nqfx), wvl_pe(nqfx), wfx_cd(nqfx)
! added 2 lines organize local variables and include files by kamata 2021/05/31
! function
      integer    lenx
!
!::index
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   it = nt
      i  = icel(1,it)
      im = i - 1
      jmax  = jtmax(it)
      jmax1 = jmax-1
      ia = 1
! modified 4/1 lines organize local variables and include files by kamata 2021/05/31
!ik   m1a = 2*ia - 1
!ik   m2a = 2*ia
!ik   m3  = 2*nion + 1
!ik   m4  = 2*nion + 2
      mm  = 2*ia
!
!::ctim
      mj = lenx(cdrout)
      if( itim.le.999  ) then
        write(ctim,'(i3.3)') itim
      elseif( itim.le.9999 ) then
        write(ctim,'(i4.4)') itim
      else
        write(ctim,'(i5.5)') itim
      endif
      write(ctub,'(i2.2)') it
      write(6,'(2x,"ctim =",a)') ctim
      write(6,'(2x,"ctub =",a)') ctub
      write(6,'(2x,"cdrout =",a)') cdrout(1:mj)
!
!::q1a,q2a,q3,q4
      k = 2
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   mm = m2a
!
!::file name
      ift = 21
      write(ckey,'(a2,i1)') "Pc", k
      write(6,'(2x,"ckey =",a)') ckey
      fnm = "it"//ctub//"_"//ckey//"_"//ctim(1:lenx(ctim))
      write(6,'(2x,"fnm =",a)') fnm(1:lenx(fnm))
      open( unit=ift, file=cdrout(1:mj)//fnm )
!
      write(ift,'(/2x,"*** plqcon_pcn ***  ",a,"  it =",i3,
     > "  itim =",i6,"  nlp =",i2,"  time =",1pe12.3,"  mm =",i2)')
     >  ckey,it,itim,nlp,time,mm
!
      write(ift,'(2x,"jw",2x,"j",3x,"i",4x,
     >   "dvl_pe",6x,"dvl_pe",6x,"dvl_pe2",6x,"hvsb_V",6x,
     >   "j-1/2", 7x,"j+1/2", 7x, "zfpa", 8x,
     >   "pe_j-1/2", 4x, "pe_j+1/2",4x,"pconB",7x,"pcon")')
!
!::flux at j+1/2
      do jw = 2, jmax1
      j  = jcel(jw,it)
      jm = jcel(jw-1,it)
      jp = jcel(jw+1,it)
!
      zfpa = aza(ia)*vna(j,i,ia)/vne(j,i)
      zlb  = hdxp(jm,i)/hpit(jm,i) + hdxm(j,i)/hpit(j,i)
      zpe  = vne(j,i)*vte(j,i)*cev
      zpep = vne(jp,i)*vte(jp,i)*cev
      zpem = vne(jm,i)*vte(jm,i)*cev
      zpehp = hwtp(j,i)*zpe+hwtm(j,i)*zpep
      zpehm = hwtp(jm,i)*zpem+hwtm(jm,i)*zpe
      zf22  = +zfpa*hvol(j,i)/zlb*(zpehp-zpehm)
      zf22p = +zfpa*hvsb(j,i)*zpehp - zfpa*hvsb(jm,i)*zpehm
      wvl_pe(j) = zf22p
!
      pfx_cv(j) = qfx_cv(j,i,mm) + zfpa*hvsb(j,i)*zpehp
     >                           + qfx_cd(j,i,mm)
      if( jw.eq.2 ) then
      pfx_cv(jm) = qfx_cv(jm,i,mm) + zfpa*hvsb(jm,i)*zpehm
     >                             + qfx_cd(jm,i,mm)
      endif
!
      write(ift,'(3i4,1p13e12.4)') jw,j,i,
     >  qvl_pe(j,i,mm), zf22, zf22p,
     >  hvol(j,i)/zlb, hvsb(jm,i), hvsb(j,i),
     >  zfpa, zpehm, zpehp, pfx_cv(jm), pfx_cv(j)
     > ,hvol(j,i)/zlb*(zpehp-zpehm)
     > ,hvsb(j,i)*zpehp - zfpa*hvsb(jm,i)*zpehm
      enddo
!
!::flux |~| = -eta*grad(v//) at j
      do jw = 2, jmax1
      j  = jcel(jw,it)
      jm = jcel(jw-1,it)
      wfx_cd(j) = 0.5d0*(qfx_cd(j,i,mm)+qfx_cd(jm,i,mm))
      enddo
      jw = 1; j = jcel(jw,it)
      wfx_cd(j) = qfx_cd(j,i,mm)
      jw = jmax; j = jcel(jw,it); jm = jcel(jw-1,it)
      wfx_cd(j) = qfx_cd(jm,i,mm)
!
!::debug write
      write(ift,'(2x)')
      write(ift,'(/2x,"*** plqcon_pcn (ptot) ***  ",a,"  it =",i3,
     > "  itim =",i6,"  nlp =",i2,"  time =",1pe12.3,"  mm =",i2)')
     >  ckey,it,itim,nlp,time,mm
      write(ift,'(2x,"jw",2x,"j",3x,"i",4x,
     >  "pe",9x,"pi",9x,"mnv2",7x,"visc",7x,
     >  "j",10x,"hvs_j",6x,"S*ptot_j",2x,"ptot_j",5x,
     >  "jh",9x,"hvs_jh",6x,"S*ptot_jh",2x,"ptot_jh")')
      do jw = 2, jmax
      j  = jcel(jw,it)
      jm = jcel(jw-1,it)
      jp = jcel(jw+1,it)
!
      if( jw.eq.1 ) then
      zhvsb = hvsb(j,i)
      elseif( jw.eq.jmax ) then
      zhvsb = hvsb(jm,i)
      else
      zlb   = hdxp(jm,i)/hpit(jm,i) + hdxm(j,i)/hpit(j,i)
      zhvsb = hvol(j,i)/zlb
      endif
!
      zptot = vne(j,i)*vte(j,i)*cev + vna(j,i,ia)*vti(j,i)*cev
     >      + vna(j,i,ia)*ama(ia)*vva(j,i,ia)**2
     >      + wfx_cd(j)/zhvsb
!
      if( jw.ne.jmax ) then
      write(ift,'(3i4,1p12e11.3)') jw, j, i,
     >  vne(j,i)*vte(j,i)*cev, vna(j,i,ia)*vti(j,i)*cev,
     >    vna(j,i,ia)*ama(ia)*vva(j,i,ia)**2, wfx_cd(j)/zhvsb,
     >  dfloat(j), zhvsb,  zhvsb*zptot, zptot,
     >  dfloat(j)+0.5d0, hvsb(j,i), pfx_cv(j), pfx_cv(j)/hvsb(j,i)
      else
      write(ift,'(3i4,1p12e11.3)') jw, j, i,
     >  vne(j,i)*vte(j,i)*cev, vna(j,i,ia)*vti(j,i)*cev,
     >    vna(j,i,ia)*ama(ia)*vva(j,i,ia)**2, wfx_cd(j)/zhvsb,
     >  dfloat(j), zhvsb,  zhvsb*zptot, zptot
      endif
!
      enddo
!
      write(ift,'(2x)')
      write(ift,'(/2x,"*** plqcon_pcn (inc pe) ***  ",a,"  it =",i3,
     > "  itim =",i6,"  nlp =",i2,"  time =",1pe12.3,"  mm =",i2)')
     >  ckey,it,itim,nlp,time,mm
      write(ift,'(2x,"jw",2x,"j",3x,"i",3x,
     >  "vx_cvB", 4x, "vx_cv", 5x, "dvx_cv",4x,"dvx_cv2",3x,"err",6x,
     >  "-dvx_df", 3x, "-dvy_df", 3x, "-dvx_cd", 3x, "-dvy_cd", 3x,
     >  "-dvl_pe", 3x, "-dvl_cl", 3x, "-dvl_pi", 3x, "-dvx_vh", 3x
     >  "-dvy_vh", 3x, "+dvl_sc", 3x, "-dvl_dt")')
!
      do jw = 2, jmax1
      j  = jcel(jw,it)
      jm = jcel(jw-1,it)
      jp = jcel(jw+1,it)
!
!xx   zvx_cv = qfx_cv(j,i,mm)-qfx_cv(jm,i,mm)
      zvx_cv = pfx_cv(j)-pfx_cv(jm)
      zvx_df = qfx_df(j,i,mm)-qfx_df(jm,i,mm)
      zvy_df = qfy_df(j,i,mm)-qfy_df(j,im,mm)
      zvx_cd = qfx_cd(j,i,mm)-qfx_cd(jm,i,mm)
      zvx_cd = 0.0d0
      zvy_cd = qfy_cd(j,i,mm)-qfy_cd(j,im,mm)
!xx   zvl_pe =  qvl_pe(j,i,mm)
      zvl_pe =  qvl_pe(j,i,mm) - wvl_pe(j)
      zvl_cl =  qvl_cl(j,i,mm)
      zvl_sc =  qvl_sc(j,i,mm)
      zvl_pi =  0.0d0
      zvx_vh =  0.0d0
      zvy_vh =  0.0d0
      if( k.eq.2 ) then
        zvl_pi =  qvl_pi(j,i,ia)
      elseif( k.eq.3 ) then
        zvx_vh = qfx_vh(j,i)-qfx_vh(jm,i)
        zvy_vh = qfy_vh(j,i)-qfy_vh(j,im)
      endif
      zvl_dt =  qvl_dt(j,i,mm)
      zsum   =  - zvx_df - zvy_df - zvx_cd - zvy_cd
     >          - zvl_pe - zvl_cl
     >          - zvl_pi - zvx_vh - zvy_vh
     >          + zvl_sc - zvl_dt
!
      write(ift,'(3i4,1p17e10.2)') jw,j,i,
!xx     >  qfx_cv(jm,i,mm), qfx_cv(j,i,mm), zvx_cv, zsum, zvx_cv-zsum,
     >  pfx_cv(jm), pfx_cv(j), zvx_cv, zsum, zvx_cv-zsum,
     >  -zvx_df, -zvy_df, -zvx_cd, -zvy_cd,
     >  -zvl_pe, -zvl_cl, -zvl_pi, -zvx_vh, -zvy_vh,
     >  +zvl_sc, -zvl_dt
!
      enddo  ! loop (jw)
!
      close(ift)
!
      return
      end
