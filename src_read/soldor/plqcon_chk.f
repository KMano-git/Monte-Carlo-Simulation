!***********************************************************************
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   subroutine plqcon_chk(nt)
      subroutine plqcon_chk(it)
!***********************************************************************
!
!      dQ/dt + div(F) + div(G) + Hp + Hc = S(j)
!
!      dQ/dt = -div(F) -div(G) -Hp -Hc +S
!            = -F(j+1/2)+F(j-1/2) -G(i+1/2)+G(i-1/2) -Hp -Hc + S
!
!-----------------------------------------------------------------------
      use cplcom, only : nion, nlp
      use cplmet, only : icel, jcel, jtmax
      use cplqcn, only : qfx_cd, qfx_cv, qfx_df, qfx_vh, qfy_cd, qfy_df
     >    , qfy_vh, qvl_al, qvl_cl, qvl_dt, qvl_pe, qvl_pi, qvl_sc
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
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  it, i, im, jmax, jmax1, jw, j, jm, mj, lenx
      integer  i, im, jmax, jmax1, jw, j, jm, mj, lenx
      integer  k, ift1, ift2
      integer  ia, m1a, m2a, m3, m4, mm
      real*8   zvx_cv, zvx_df, zvy_df, zvx_cd, zvy_cd
      real*8   zvl_pe, zvl_cl, zvl_sc, zvl_pi, zvx_vh, zvy_vh
      real*8   zvl_al, zvl_dt, zsum
      real*8   fvx_cv, fvx_df, fvy_df, fvx_cd, fvy_cd
      real*8   fvl_pe, fvl_cl, fvl_sc, fvl_pi, fvx_vh, fvy_vh
      real*8   fvl_al, fvl_dt, fsum
      character  ctim*5, ctub*2, ckey1*3, ckey2*3, fnm1*80, fnm2*80
!
!::index
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   it = nt
      i  = icel(1,it)
      im = i - 1
      jmax  = jtmax(it)
      jmax1 = jmax-1
      ia = 1
      m1a = 2*ia - 1
      m2a = 2*ia
      m3  = 2*nion + 1
      m4  = 2*nion + 2
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
!
!::q1a,q2a,q3,q4
      do k = 1, 4
      if( k.eq.1 ) mm = m1a
      if( k.eq.2 ) mm = m2a
      if( k.eq.3 ) mm = m3
      if( k.eq.4 ) mm = m4
!
!::file name
      ift1 = 21
      ift2 = 22
      write(ckey1,'(a2,i1)') "dq", k
      write(ckey2,'(a2,i1)') "Fq", k
      fnm1 = "it"//ctub//"_"//ckey1//"_"//ctim(1:lenx(ctim))
      fnm2 = "it"//ctub//"_"//ckey2//"_"//ctim(1:lenx(ctim))
      open( unit=ift1, file=cdrout(1:mj)//fnm1 )
      open( unit=ift2, file=cdrout(1:mj)//fnm2 )
!
      write(ift1,'(/2x,"*** plqcon_chk(dq/dt) ***  ",a,"  it =",i3,
     > "  itim =",i6,"  nlp =",i2,"  time =",1pe12.3,"  mm =",i2)')
     >  ckey1,it,itim,nlp,time,mm
      write(ift1,'(3x,"jw",3x,"j",4x,"i",4x,"-vx_cv",5x,"-vx_df",5x,
     >  "-vy_df",5x,"-vx_cd",5x,"-vy_cd",5x,"-vl_pe",5x,"-vl_cl",5x,
     >  "vl_sc",6x,"-vl_pi",5x,"-vx_vh",5x,"-vy_vh",5x,"sum",8x,
     >  "vl_al",6x,"vl_dt",6x,"res")')
!
      write(ift2,'(/2x,"*** plqcon_chk(flx) ***  ",a,"  it =",i3,
     > "  itim =",i6,"  nlp =",i2,"  time =",1pe12.3,"  mm =",i2)')
     >  ckey2,it,itim,nlp,time,mm
      write(ift2,'(3x,"jw",3x,"j",4x,"i",4x,"vx_cv",6x,"vx_df",6x,
     >  "vy_df",6x,"vx_cd",6x,"vy_cd",6x,"vl_pe",6x,"vl_cl",6x,
     >  "vl_sc",6x,"vl_pi",6x,"vx_vh",6x,"vy_vh",6x,"sum",8x,
     >  "vl_al",6x,"vl_dt",6x,"res")')
!
      do jw = 1, jmax1
      j  = jcel(jw,it)
      jm = j
      if( jw.ne.1 ) jm = jcel(jw-1,it)
!
      zvx_cv = qfx_cv(j,i,mm)-qfx_cv(jm,i,mm)
      zvx_df = qfx_df(j,i,mm)-qfx_df(jm,i,mm)
      zvy_df = qfy_df(j,i,mm)-qfy_df(j,im,mm)
      zvx_cd = qfx_cd(j,i,mm)-qfx_cd(jm,i,mm)  ! zero
      zvy_cd = qfy_cd(j,i,mm)-qfy_cd(j,im,mm)
      zvl_pe =  qvl_pe(j,i,mm)
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
      zsum   =- zvx_cv - zvx_df - zvy_df - zvx_cd - zvy_cd
     >        - zvl_pe - zvl_cl + zvl_sc
     >        - zvl_pi - zvx_vh - zvy_vh
      zvl_al =  qvl_al(j,i,mm)
      zvl_dt =  qvl_dt(j,i,mm)
!
      if( jw.eq.1 ) then
      fvx_cv = +qfx_cv(j,i,mm)
      fvx_df = +qfx_df(j,i,mm)
      fvy_df = +qfy_df(j,i,mm)
      fvx_cd = +qfx_cd(j,i,mm)
      fvy_cd = +qfy_cd(j,i,mm)
      fvl_pe =  0.0d0
      fvl_cl =  0.0d0
      fvl_sc =  0.0d0
      fvl_pi =  0.0d0
      fvx_vh =  0.0d0
      fvy_vh =  0.0d0
      if( k.eq.3 ) then
        fvx_vh = +qfx_vh(j,i)
        fvy_vh = +qfy_vh(j,i)
      endif
      fsum   =  0.0d0
      fvl_al =  0.0d0
      fvl_dt =  0.0d0
      else
      fvx_cv = fvx_cv + zvx_cv
      fvx_df = fvx_df + zvx_df
      fvy_df = fvy_df + zvy_df + qfy_df(j,im,mm)
      fvx_cd = fvx_cd + zvx_cd
      fvy_cd = fvy_cd + zvy_cd + qfy_cd(j,im,mm)
      fvl_pe = fvl_pe + zvl_pe
      fvl_cl = fvl_cl + zvl_cl
      fvl_sc = fvl_sc + zvl_sc
      fvl_pi = fvl_pi + zvl_pi
      fvx_vh = fvx_vh + zvx_vh
      if( k.ne.3 ) then
      fvy_vh = fvy_vh + zvy_vh
      else
      fvy_vh = fvy_vh + zvy_vh + qfy_vh(j,im)
      endif
      fsum   = fsum
     >        + zvx_cv + zvx_df + zvy_df + zvx_cd + zvy_cd
     >        + zvl_pe + zvl_cl + zvl_sc
     >        + zvl_pi + zvx_vh + zvy_vh
      fvl_al =  fvl_al + qvl_al(j,i,mm)
      fvl_dt =  fvl_dt + qvl_dt(j,i,mm)
      endif
!
      write(ift1,'(3i5,1p15e11.2)') jw,j,i,
     >  -zvx_cv, -zvx_df, -zvy_df, -zvx_cd, -zvy_cd,
     >  -zvl_pe, -zvl_cl, +zvl_sc, -zvl_pi, -zvx_vh, -zvy_vh,
     >   zsum,    zvl_al, +zvl_dt, +zvl_dt-zvl_al
!
      write(ift2,'(3i5,1p15e11.2)') jw,j,i,
     >  fvx_cv, fvx_df, fvy_df, fvx_cd, fvy_cd,
     >  fvl_pe, fvl_cl, fvl_sc, fvl_pi, fvx_vh, fvy_vh,
     >  fsum,   fvl_al, fvl_dt, -fvl_dt+fvl_al
!
      enddo  ! loop (jw)
!
      close(ift1)
      close(ift2)
!
      enddo  ! loop (k)
!
      return
      end
