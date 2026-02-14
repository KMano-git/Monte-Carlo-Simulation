!***********************************************************************
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   subroutine plqcon_flx(nt)
      subroutine plqcon_flx(it)
!***********************************************************************
!
!     dq/dt + div(F) + div(G) + Hp = S
!
!      F^sita(j+1/2) - F^sita(j-1/2)
!            = - Vj*div(F) - Vj*div(G) + Vj*Source - Vj*dq/dt
!
!-----------------------------------------------------------------------
      use cplcom, only : nion, nlp, vna, vte, vti, vva
      use cplmet, only : icel, jcel, jtmax
      use cplqcn, only : qfx_cd, qfx_cv, qfx_df,qfx_vh,  qfy_cd, qfy_df
     >    , qfy_vh, qvl_cl, qvl_dt, qvl_pe, qvl_pi, qvl_sc
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
!ik   integer  it, i, im, jmax, jmax1, ia, mj, k, mm
      integer  i, im, jmax, jmax1, ia, mj, k, mm
      integer  ift, jw, j, jm, lenx
      integer  m1a, m2a, m3, m4
      character  ctim*5, ctub*2, ckey*3, fnm*80
      real*8   zvx_cv, zvx_df, zvy_df, zvx_cd, zvy_cd, zvl_pe
      real*8   zvl_sc, zvl_pi, zvx_vh, zvy_vh, zvl_dt, zvl_cl, zsum
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
      write(ctub,'(i2.2)') it   ! KS 2011/09/22
!xx      write(6,'(2x,"ctim =",a)') ctim
!xx      write(6,'(2x,"ctub =",a)') ctub
!xx      write(6,'(2x,"cdrout =",a)') cdrout(1:mj)
!
!::q1a,q2a,q3,q4
      do k = 1, 4
      if( k.eq.1 ) mm = m1a
      if( k.eq.2 ) mm = m2a
      if( k.eq.3 ) mm = m3
      if( k.eq.4 ) mm = m4
!
!::file name
      ift = 21
      write(ckey,'(a2,i1)') "Fv", k
!xx      write(6,'(2x,"ckey =",a)') ckey
      fnm = "it"//ctub//"_"//ckey//"_"//ctim(1:lenx(ctim))
!xx      write(6,'(2x,"fnm =",a)') fnm(1:lenx(fnm))
      open( unit=ift, file=cdrout(1:mj)//fnm )
!
      write(ift,'(/2x,"*** plqcon_flx ***  ",a,"  it =",i3,
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
      zvl_dt =  qvl_dt(j,i,mm)
      zsum   =  - zvx_df - zvy_df - zvx_cd - zvy_cd
     >          - zvl_pe - zvl_cl
     >          - zvl_pi - zvx_vh - zvy_vh
     >          + zvl_sc - zvl_dt
!
      write(ift,'(3i4,1p16e10.2)') jw,j,i,
     >  qfx_cv(jm,i,mm), qfx_cv(j,i,mm), zvx_cv, zsum, zvx_cv-zsum,
     >  -zvx_df, -zvy_df, -zvx_cd, -zvy_cd,
     >  -zvl_pe, -zvl_cl, -zvl_pi, -zvx_vh, -zvy_vh,
     >  +zvl_sc, -zvl_dt
!
      enddo  ! loop (jw)
!
      write(ift,'(//)')
      write(ift,'(2x,"jw",2x,"j",3x,"i",3x,
     >   "Na",10x,"Va",10x,"Ti",10x,"Te",10x,"fx_cv",7x,"fx_cd")')
!
      do jw = 1, jmax
      j  = jcel(jw,it)
      write(ift,'(3i4,1p10e12.3)') jw, j, i
     >  ,vna(j,i,ia), vva(j,i,ia), vti(j,i), vte(j,i)
     >  ,qfx_cv(j,i,mm), qfx_cd(j,i,mm)
      enddo
!
      close(ift)
!
      enddo  ! loop (k)
!
      return
      end
