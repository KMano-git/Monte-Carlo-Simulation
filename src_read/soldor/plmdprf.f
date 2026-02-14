!**********************************************************************
      subroutine plmdprf(kset,kout)
!**********************************************************************
!)
!)   plasma profile at mid-plane
!)       kset =  1 : vna,vti,vte ==> mdp_ni, mdp_ti, mdp_te
!)            =  2 : mdp_ni, ... ==> anamp,  atimp, atemp
!)            =  0 : No effect for debug check
!)
!)::Old
!)   *--- prg_init
!)   *--- prg_parm
!)          *--- pls_ini
!)                 *--- plrprf
!)                        *--- plmdprf(2) => PLS2
!)   *--- prg_prof
!)   *--- prg_exec
!)          *--- plmdprf(2) => PLS2
!)   *--- prg_term
!)
!)::New
!)   *--- prg_exec
!)          *--- pls_cal
!)          *--- plmdprf(2)
!)
!)
!) plmd 45 8.74E-01 8.65E-01 4.6020E+19 1.6570E+03 1.6120E+03 parm 0
!) plmd 45 8.74E-01 8.65E-01 2.4760E+19 3.6759E+03 3.3033E+03 parm 0 X
!) plmd 45 8.74E-01 8.65E-01 2.4760E+19 3.6759E+03 3.3033E+03 exec 0 X
!) plmd 45 8.74E-01 8.65E-01 2.4760E+19 3.6759E+03 1.8795E+03 exec 1
!) plmd 45 8.74E-01 8.65E-01 2.4760E+19 3.6759E+03 1.8785E+03 exec 2
!) plmd 45 8.74E-01 8.65E-01 2.4760E+19 3.6759E+03 1.8779E+03 exec 3
!)
!)--------------------------------------------------------------------
      use cplcom, only : anamp, animp, atemp, atimp, fcna, nion, vna
     >    , vni, vte, vti
      use cplmet, only : icaxs, icmpe, icmps, icwl1
      use cpmpls, only : imd1, jmd1, lnedg,  mdp_ni, mdp_rh, mdp_ro
     >    , mdp_te, mdp_ti, prfni, prfte
     >    , prfti, rohmp, romp, xmd1, ymd1
      use csize,  only : ndsp
      use csonic, only : itim, time
      use cunit,  only : n6
      implicit none
!
!::argument
      integer, intent(in) :: kset, kout
!
!::local variables
      integer ia, jc, ic, ic1, iout
      real*8  zro, zro1, fnib, ftib, fteb
      real*8  zdni, zdti, zdte, zfac(ndsp)
      real*8  zni0, zti0, zte0, zni1, zti1, zte1
      logical :: is_Error = .false.
! function
      real(8)    fprfmp

      if( jmd1.eq.0 .or. imd1.eq.0 ) then
        call wexit("plmdprf","undefined index of outer mid-plane")
      endif

      iout = kout
      if( itim < 5 ) iout = 1
!
!----------------------------------------------------------------------
!::vni,vti,vte ==> mdp_ni, mdp_ti, mdp_te
!----------------------------------------------------------------------
      if( kset.eq.1 ) then
        do ic = icwl1, icaxs
          mdp_rh(ic) = rohmp(ic)
          mdp_ro(ic) = romp(ic)
        enddo
        jc = jmd1
        do ic = icwl1, icmpe
          mdp_ni(ic) = vni(jc,ic)
          mdp_ti(ic) = vti(jc,ic)
          mdp_te(ic) = vte(jc,ic)
        enddo
        ic = icmpe
        zni0 = vni(jc,ic)
        zti0 = vti(jc,ic)
        zte0 = vte(jc,ic)
!
!::hot core region
        ic1  = icmpe
        zro1 = rohmp(ic1)
        fnib = (1.0d0-zro1**2)**prfni(4)
        ftib = (1.0d0-zro1**2)**prfti(4)
        fteb = (1.0d0-zro1**2)**prfte(4)
!-----
        prfni(2) = (zni0-prfni(1)*fnib)/(1.0d0-fnib)
        prfti(2) = (zti0-prfti(1)*ftib)/(1.0d0-ftib)
        prfte(2) = (zte0-prfte(1)*fteb)/(1.0d0-fteb)
!-----
!
        zni1 = fprfmp(prfni,zro1)
        zti1 = fprfmp(prfti,zro1)
        zte1 = fprfmp(prfte,zro1)
!
        do ic = ic1, icaxs
          zro = rohmp(ic)
          mdp_ni(ic) = fprfmp(prfni,zro) - zni1 + zni0
          mdp_ti(ic) = fprfmp(prfti,zro) - zti1 + zti0
          mdp_te(ic) = fprfmp(prfte,zro) - zte1 + zte0
        enddo
!
        if( iout.eq.1 ) then
          write(n6,'(/2x,"*** plmdprf ***  kset =",i2,i7,1pe12.4)')
     >       kset, itim, time
          write(n6,'(2x,"out-mid =",2f8.4,2i5)') xmd1,ymd1,jmd1,imd1
          write(n6,'(2x,"icmpe =",i3,"  vni,vti,vte =",1p3e11.3,
     >       "  func =",1p3e11.3,"  mdp =",1p3e11.3)') icmpe
     >      ,zni0, zti0, zte0, zni1, zti1, zte1
     >      ,mdp_ni(icmpe), mdp_ti(icmpe), mdp_te(icmpe)
          write(n6,'(2x,3x,"ic",3x,"rh",10x,"ro",10x,"Ni",10x,"Ti",10x,
     >      "Te")')
          do ic = icwl1, icaxs
            write(n6,'(2x,i5,1p5e12.3)')
     >        ic, mdp_rh(ic), mdp_ro(ic), mdp_ni(ic)
     >        , mdp_ti(ic), mdp_te(ic)
          enddo
        endif  ! lout = 1
      endif  ! kset = 1
!
!----------------------------------------------------------------------
!::mdp_ni, mdp_ti, mdp_te ==> anamp, atimp, atemp
!----------------------------------------------------------------------
      if( kset.eq.2 ) then
        jc = jmd1
        ic = icmpe
!
        if( vni(jc,ic).le.0.0d0 ) then
          zdni = 0.0d0
          zdte = 0.0d0
          zdti = 0.0d0
          do ia = 1, nion
            zfac(ia) = fcna(ia)
          enddo
        elseif( lnedg.eq.1 )then
          zdni = 0.0d0
          zdte = 0.0d0
          zdti = 0.0d0
          do ia = 1, nion
            zfac(ia) = vna(jc,ic,ia)/vni(jc,ic)
         enddo
        else
          zdni = vni(jc,ic) - mdp_ni(ic)
          zdti = vti(jc,ic) - mdp_ti(ic)
          zdte = vte(jc,ic) - mdp_te(ic)
          do ia = 1, nion
            zfac(ia) = vna(jc,ic,ia)/vni(jc,ic)
          enddo
        endif
!
        if( iout.eq.1 ) then
          write(n6,'(/2x,"*** plmdprf ***  kset =",i2,i7,2x,1pe12.4)')
     >     kset, itim, time
          write(n6,'(2x,"core edge  jc,ic =",2i5,"  Ni =",1pe11.3,
     >     "  Ti =",1pe11.3,"  Te =",1pe11.3)')
     >     jc,ic,vni(jc,ic),vti(jc,ic),vte(jc,ic)
          write(n6,'(2x,"smooth  dNi =",1pe11.3,"  dTi =",1pe11.3,
     >     "  dTe =",1pe11.3,"  fac =",0p10f8.4)')
     >      zdni, zdti, zdte, (zfac(ia),ia=1,nion)
        endif
!
        do ic = icmps, icaxs
          animp(ic) = mdp_ni(ic) + zdni
          atimp(ic) = mdp_ti(ic) + zdti
          atemp(ic) = mdp_te(ic) + zdte
          if((animp(ic)<0 .or. atimp(ic)<0 .or. atemp(ic)<0)
     >      .and. ic > icmpe) then
            is_Error = .true.
          endif
          do ia = 1, nion
            anamp(ic,ia) = animp(ic)*zfac(ia)
          enddo
        enddo

        if(is_Error) then
          write(n6,'(/2x,"*** Error in plmdprf ***")')
          write(n6,'(2x,"zdni=",1pe11.3,2x,"zdti=",1pe11.3
     >     ,2x,"zdte=",1pe11.3,2x,"icmpe=",i5)')zdni,zdti,zdte,icmpe
          write(n6,'(2x,"ic",4x,"animp(ic)",3x,"mdp_ni(ic)"
     >     ,4x,"atimp(ic)",3x,"mdp_ti(ic)",4x
     >     ,"atemp(ic)",3x,"mdp_te(ic)")')
          write(n6,'(2x,i2,2x,1pe11.3,2x,1pe11.3,
     >     2x,1pe11.3,2x,1pe11.3
     >     ,2x,1pe11.3,2x,1pe11.3)')(ic,animp(ic)
     >     ,mdp_ni(ic),atimp(ic),mdp_ti(ic),atemp(ic),mdp_te(ic)
     >     ,ic=icmps,icaxs)
          call wexit("plmdprf","animp(ic).or.atimp(ic).or.atemp(ic)<0")
        endif
!
        if( iout.eq.1 ) then
          write(n6,'(2x,"*** plmdprf ***  itim =",i6)') itim

          do ic = icmpe, icmpe+5
            write(n6,'(2x,a,2x,i5,1p5e15.6)') "plmdprf",
     >      ic, rohmp(ic), romp(ic), anamp(ic,1), atimp(ic), atemp(ic)
          enddo
        endif  ! iout = 1
      endif  ! kset = 2
!
!----------------------------------------------------------------------
!::No effect when kset = 0  for debug check   2016/05/19
!----------------------------------------------------------------------
      if( kset.eq.0 .and. iout.eq.1) then
        write(n6,'(2x,"*** plmdprf ***  itim =",i6)') itim
        do ic = icmpe, icmpe+5
          write(n6,'(2x,a,2x,i5,1p5e15.6)') "plmdprf",
     >      ic, rohmp(ic), romp(ic), anamp(ic,1), atimp(ic), atemp(ic)
        enddo
      endif  ! iout = 1
!
      return
      end
!
!**********************************************************************
      real(8) function fprfmp( fit, x )
!**********************************************************************
      real(8), intent(in) :: fit(10), x
!
      fprfmp = (fit(1)-fit(2))*(1.0d0-x**2)**fit(4) + fit(2)
!
      return
      end
