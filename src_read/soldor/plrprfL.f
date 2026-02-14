!**********************************************************************
      subroutine plrprf(lscal)
!**********************************************************************
!
!    meaning of sub-index     pointer of ic
!      c : center               icaxs : plasma axis
!      p : edge condition       icmpe : plasma edge
!      q : = ip-1               icmps : plasma surface
!      b : plasma surface       icspx : flux tube close to separatrix
!                               icwl1 : sol wall
!
!    radial profile in sol
!     fac(ic) = (dmid(ic)-dmid(icspx))/(dmid(icwl1)-dmid(icspx))
!      na(ic) = nasx - (nasx-naws)*fac(ic)
!      Ti(ic) = Tisx - (Tisx-Tiws)*fac(ic)
!      Te(ic) = Tesx - (Tesx-Tews)*fac(ic)
!
!       where Tisx = Ti(icspx), Tiws = Ti(icwl1)
!
!       fprna ==> fcna    2005/7/10
!
!-----------------------------------------------------------------------
      use cplcom, only : anamp, anapv, anasl, animp, anipv, anisl, atemp
     >    , atepv, atesl, atimp, atipv, atisl, bwpni, bwpte, bwpti
     >    , bwsni, bwste, bwsti, cfti, dmid, fcna, nfti, nion, trad_max
     >    , trad_min, vna, vni, vte, vti
      use cplmet, only : icaxs, icmpe, icmps, icspx, icwl1, icwl2, itmps
     >    , itpvs
      use cpmpls, only : arhmp, bniedg, clsav, clsfc, clsni, clste
     >    , clsti, clsxp, dn0c, dn0s, fflna, flxni, flxqe, flxqi, fprna
     >    , jmd1, lnedg, mdp_ni, mdp_rh, mdp_ro, mdp_te, mdp_ti, prfna
     >    , prfni, prfte, prfti, pufs, rmdc, rmds, rohmp, taup, wlose
     >    , wlosi, ltedg
     >    , nty_ctl, rad_const_step, rad_const_fac
      use csize,  only : ndsp, ndy
      use csonic, only : lmstd, lntl
      use cunit,  only : lmspe, lmype, lnope, mywld, n5, n6, wh_dir
      use mpi!,    only : mpi_bcast, mpi_byte
      implicit none
!
!::argument
      integer, intent(in) :: lscal
!
!::local variables
      real*8  bnip, btip, btep, bnix, btix, btex
      real*8  zni1, zti1, zte1, zni2, zti2, zte2, zfac
      integer jc, ic, ia, ic1, ic2, it, i6, icxx
      integer icmpsx, icmpex, icaxsx, mj
      integer ktmax, ntrmax; real*8  dtstep
      integer nls, nle, nbf, ierr
      character  clin*80, cdsn*80, cmesh*20
      logical lex
      integer    itag, ityp
! function
      integer    lenx
!
      namelist /umpinp/ flxni, flxqi, flxqe, fflna,
     >  trad_min, trad_max,
     >  prfni, prfti, prfte, fprna,
     >  bnip, btip, btep, bnix, btix, btex,
     >  lnedg, ltedg, nfti, cfti,
     >  ktmax, ntrmax, dtstep, rmdc, dn0c, taup, wlosi, wlose,
     >  rmds, dn0s, pufs,
     >  clsav, clsxp, clsfc, clsni, clsti, clste
!::radiation power constraint
     >  , nty_ctl, rad_const_step, rad_const_fac
!
      write(n6,'(/2x,"***  plrprf  ***")')
!
!::KSFUJI  zero clear in plrprf
      rmdc(1:2) = 0.0d0
      dn0c(1:2) = 0.0d0
      rmds(1:2) = 0.0d0
      dn0s(1:2) = 0.0d0
      call deflt_umpinp
!
!-----------------------------------------------------------------------
!::input data
!-----------------------------------------------------------------------
      read(n5,umpinp)
!
!::lscal (=0:read file/=1:calculation)
!
      if( nion.eq.1 ) fflna(1) = 1.0d0
      do ia = 1, nion
      prfna(1,ia) = prfni(1)*fcna(ia)
      prfna(2,ia) = prfni(2)*fcna(ia)
      prfna(3,ia) = prfni(3)
      prfna(4,ia) = prfni(4)
      enddo
!
      if( trad_min.le.0.0d0 ) trad_min = 0.3d0*(flxqi+flxqe)
      if( trad_max.le.0.0d0 ) trad_max = 0.8d0*(flxqi+flxqe)
!
      bniedg = bnip
!
!-----------------------------------------------------------------------
!::debug write
!-----------------------------------------------------------------------
!
!::flux from the main plasma
      write(n6,'(5x,"lnedg =",i2)') lnedg
      write(n6,'(5x,"ltedg =",i2)') ltedg
      write(n6,'(5x,"fflna =",10f7.3)') (fflna(ia),ia=1,nion)
      write(n6,'(5x,"flxni =",1pe11.3,"  flxqi =",1pe11.3,"  flxqe =",
     >  1pe11.3)') flxni,flxqi,flxqe
      write(n6,'(5x,"trad_min =",1pe12.3,"  trad_max =",1pe12.3)')
     >    trad_min, trad_max
!
!::boudary value
      write(n6,'(5x,"bnip  =",1pe11.3,"  btip  =",1pe11.3,"  btep  =",
     >  1pe11.3)') bnip,btip,btep
!
!
!::initial profile
      write(n6,'(5x,"prfni =",1p6e11.3)') (prfni(ia),ia=1,6)
      write(n6,'(5x,"prfti =",1p6e11.3)') (prfti(ia),ia=1,6)
      write(n6,'(5x,"prfte =",1p6e11.3)') (prfte(ia),ia=1,6)
!
!::main plasma parameter
      write(n6,'(5x,"cfti =",a)') cfti(1:30)
!
      write(n6,'(5x,"ktmax =",i4,"  ntrmax -",i2
     >  ,"  dtstep =",1pe12.3)') ktmax,ntrmax,dtstep
      write(n6,'(5x,"rmdc =",0p2f9.4,"  dn0c =",1p2e12.3,"  rmds =",
     >  0p2f9.4,"  dn0s =",1p2e12.3)') rmdc,dn0c,rmds,dn0s
      write(n6,'(5x,"taup  =",1pe12.3,"  wlosi =",1pe12.3,"  wlose =",
     >  1pe12.3)') taup,wlosi,wlose
      write(n6,'(5x,"clsav =",f7.3,"  clsxp =",f7.3,"  clsfc =",f7.3,
     >  "  clsni/clsti/clste =",3f7.3)') clsav,clsxp,clsfc,
     >   clsni,clsti,clste
!
!-----------------------------------------------------------------------
!::preparation for analytical neutral transport
!           and recombination (plrecm)
!-----------------------------------------------------------------------
      if( lscal.eq.1 .or. lntl.eq.0 ) then
!::input data
        call ntinpt(n5)
      endif
!
!-----------------------------------------------------------------------
!::read file data of main & SOL profile ( mdp_ni,Ti,Te )
!-----------------------------------------------------------------------
      if( lscal.eq.0 ) then
        write(n6,'(/2x,"*** plrprf ***  read file")')
        write(n6,'(2x,"icwl1 =",i4,"  icaxs =",i4)') icwl1, icaxs
        write(n6,'(2x,"cfti =",a)') trim(cfti)
!
!::master pe
        if( lmype.eq.lmspe ) then
          i6 = 21
          call nopen(i6,cfti,"text",cdsn,lex)
          open(unit=i6, file=cdsn)
          read(i6,'(a)') clin
          write(n6,'(a)') clin(1:lenx(clin))
          mj = index(clin,"ic=")
          cmesh = clin(mj+3:)
          read(i6,'(a)') clin
!-----
          if( index(clin,"Rh").eq.0 ) then
            do ic = icwl1, icaxs
              read(i6,'(2x,i5,1p5e12.3)')
     >         icxx, mdp_rh(ic), mdp_ro(ic),
     >         mdp_ni(ic), mdp_ti(ic), mdp_te(ic)
            enddo
          else
            do ic = icwl1, icaxs
              read(i6,'(2x,i5,1p3e14.5,1p3e12.3)')
     >         icxx, mdp_rh(ic), mdp_ro(ic), arhmp(ic),
     >         mdp_ni(ic), mdp_ti(ic), mdp_te(ic)
            enddo
          endif
!-----
          close( i6 )
          read(cmesh,*) icmpsx, icmpex, icaxsx
          write(n6,'(2x,"icmps,icmpe,icaxs =",3i4,"  file =",3i4)')
     >      icmps, icmpe, icaxs, icmpsx, icmpex, icaxsx
          if(icmpsx.ne.icmps.or.icmpex.ne.icmpe.or.icaxsx.ne.icaxs)then
            call wexit("plrprf","mesh number of inpmpr is invalid")
          endif
        endif
!
!::[MPI_Bcast in cpmpls]  cpmpls/com_plmdprf/ (mdp_ro,emrk)  10/04/21
        if( lnope.gt.1 ) then
          call cpmpls_lsr
        endif
!
!::debug write
        write(n6,'(/2x,"initial profile  mdp_rh, mdp_ni")')
        do ic = icwl1, icaxs
          write(n6,'(2x,i5,1p3e14.5,1p3e12.3)')
     >     ic, mdp_rh(ic), mdp_ro(ic), arhmp(ic),
     >        mdp_ni(ic), mdp_ti(ic), mdp_te(ic)
        enddo
!
!-----------------------------------------------------------------------
!::calculation to get steady profile  ( mdp_ni,Ti,Te )
!-----------------------------------------------------------------------
      else
        write(n6,'(2x,"*** plrprf ***  calculation")')
        write(n6,'(5x,"bnix  =",1pe11.3,"  btix  =",1pe11.3,"  btex =",
     >    1pe11.3)') bnix,btix,btex
        zni1 = bnix;  zti1 = btix;  zte1 = btex;  ic1 = icspx
        zni2 = bwsni; zti2 = bwsti; zte2 = bwste; ic2 = icwl1
        jc = jmd1
        do ic = ic1, ic2, -1
          zfac = (dmid(ic)-dmid(ic1))/(dmid(ic2)-dmid(ic1))
          vni(jc,ic) = zni1 + (zni2-zni1)*zfac
          vti(jc,ic) = zti1 + (zti2-zti1)*zfac
          vte(jc,ic) = zte1 + (zte2-zte1)*zfac
          do ia = 1, nion
            vna(jc,ic,ia) = vni(jc,ic)*fcna(ia)
          enddo
        enddo
!
!--in core edge
        zni1 = bnip; zti1 = btip; zte1 = btep; ic1 = icmpe
        zni2 = bnix; zti2 = btix; zte2 = btex; ic2 = icspx
        jc = jmd1
        do ic = ic2+1, ic1
          zfac = (rohmp(ic)-rohmp(ic1))/(rohmp(ic2)-rohmp(ic1))
          if( ic.gt.ic1 ) zfac = 0.0d0
          vni(jc,ic) = zni1 + (zni2-zni1)*zfac
          vti(jc,ic) = zti1 + (zti2-zti1)*zfac
          vte(jc,ic) = zte1 + (zte2-zte1)*zfac
          do ia = 1, nion
            vna(jc,ic,ia) = vni(jc,ic)*fcna(ia)
          enddo
        enddo
!
!--vna,vti,vte ==> mdp_ni, mdp_ti, mdp_te
        call plmdprf(1,1)
      endif
!
!
!-----------------------------------------------------------------------
!::mdp_ni,mdp_ti,mdp_te ==> animp, anisl, anipv
!-----------------------------------------------------------------------
!::KSFUJI  set value for undefined variables
      anamp(1:ndy,1:ndsp) = 0.0d0
      animp(1:ndy) = 0.0d0
      atimp(1:ndy) = 0.0d0
      atemp(1:ndy) = 0.0d0
      anasl(1:ndy,1:ndsp) = 0.0d0
      anisl(1:ndy) = 0.0d0
      atisl(1:ndy) = 0.0d0
      atesl(1:ndy) = 0.0d0
      anapv(1:ndy,1:ndsp) = 0.0d0
      anipv(1:ndy) = 0.0d0
      atipv(1:ndy) = 0.0d0
      atepv(1:ndy) = 0.0d0
!
!::main
      call plmdprf(2,1)
!
!::SOL
      do ic = icwl1, icspx
      anisl(ic) = mdp_ni(ic)
      atisl(ic) = mdp_ti(ic)
      atesl(ic) = mdp_te(ic)
      do ia = 1, nion
      anasl(ic,ia) = anisl(ic)*fcna(ia)
      enddo
      enddo
!
!::prv
      do ic = icspx+1, icwl2
      zfac = dfloat(ic-icwl2)/dfloat(icspx-icwl2)
      anipv(ic) = bwpni + (mdp_ni(icspx)-bwpni)*zfac
      atipv(ic) = bwpti + (mdp_ti(icspx)-bwpti)*zfac
      atepv(ic) = bwpte + (mdp_te(icspx)-bwpte)*zfac
      do ia = 1, nion
      anapv(ic,ia) = anipv(ic)*fcna(ia)
      enddo
      enddo
!
!-----------------------------------------------------------------------
!::calculation to get steady profile
!-----------------------------------------------------------------------
      if( lscal.eq.1 ) then
!
!--initial profile   anamp, anasl, anapv ==> vna,vti,vte
        call plinit
!
!--steady profile
        lmstd = 1
        call plmstdy(ktmax,ntrmax,dtstep)
        lmstd = 0
!
!--in hot core   vna,vti,vte ==> mdp_ni, mdp_ti, mdp_te
        call plmdprf(1,1)
!
!--output to file
        if( lmype.eq.lmspe ) then
          i6 = 21
          call nopen(i6,cfti,"text write",cdsn,lex)
          ia = 1
          write(i6,'(2x,"*** plmpdrf ***  ",a,2x,"ic=",3i5)')
     >       wh_dir(1:lenx(wh_dir)), icmps, icmpe, icaxs
          write(i6,'(2x,3x,"ic",3x,"rh",12x,"ro",12x,"Rh",12x,"Ni",10x,
     >       "Ti",10x,"Te")')
          do ic = icwl1, icaxs
            write(i6,'(2x,i5,1p3e14.5,1p3e12.3)')
     >       ic, mdp_rh(ic), mdp_ro(ic), arhmp(ic),
     >       mdp_ni(ic), mdp_ti(ic), mdp_te(ic)
          enddo
          close( i6 )
        endif
      endif
!
!----------------------------------------------------------------------
!::debug write
!----------------------------------------------------------------------
      write(n6,'(/2x,"lscal =",i3)') lscal
      write(n6,'(/2x," anamp, anasl, anapv")')
      do ic = icwl1, icaxs
      write(n6,'(2x,i5,3(1p3e11.3,2x))')
     >  ic, anamp(ic,1), atimp(ic), atemp(ic)
     >     ,anasl(ic,1), atisl(ic), atesl(ic)
     >     ,anapv(ic,1), atipv(ic), atepv(ic)
      enddo
!
      write(n6,'(2x,"plasma parameter in main plasma")')
      write(n6,'(4x,"it",2x,"ic",4x,"roh",8x,"ti",9x,"te",9x,"na")')
      ia = 1
      do ic = icmps, icaxs
      it = itmps + ic - icmps
      write(n6,'(2x,2i4,1p4e11.3)')
     >  it, ic, rohmp(ic), atimp(ic),atemp(ic),anamp(ic,ia)
      enddo
!
      write(n6,'(2x,"plasma parameter in sol plasma")')
      write(n6,'(4x,"it",2x,"ic",4x,"dmid",7x,"ti",9x,"te",9x,"na")')
      ia = 1
      do ic = icwl1, icspx
      it = ic
      write(n6,'(2x,2i4,1p4e11.3)')
     >  it, ic, dmid(ic), atisl(ic),atesl(ic),anasl(ic,ia)
      enddo
!
      write(n6,'(2x,"plasma parameter in prv plasma")')
      write(n6,'(4x,"it",2x,"ic",4x,"dmid",7x,"ti",9x,"te",9x,"na")')
      ia = 1
      do ic = icspx+1, icwl2
      it = itpvs + ic - (icspx+1)
      write(n6,'(2x,2i4,1p4e11.3)')
     >  it,ic, 0.0d0, atipv(ic),atepv(ic),anapv(ic,ia)
      enddo
!
!----------------------------------------------------------------------
!::ending
!----------------------------------------------------------------------
      write(n6,'(2x,"---  normal end  plrprf  ---")')
      end
!
!**********************************************************************
      subroutine deflt_umpinp
!**********************************************************************
      use cpmpls, only : clsav, clsfc, clsni, clste, clsti, clsxp
      implicit none
!
!xx   data  clsav,clsxp,clsfc/0.1d0, 2.0d0, 0.6d0/
!xx   data  clsni,clsti,clste/1.0d0, 1.0d0, 1.0d0/
!
      clsav = 0.1d0
      clsxp = 2.0d0
      clsfc = 0.6d0
      clsni = 1.0d0
      clsti = 1.0d0
      clste = 1.0d0
!
      end
