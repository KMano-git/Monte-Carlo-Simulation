!**********************************************************************
      subroutine mtrdsk(nft,cact)
!**********************************************************************
      use cgdcom,      only : grdx, grdy, hbr, hbt, hbz, mpax, mpman
     >    , mpprv, mpsol, mpsp, mpw1, mpw2, mqd1, mqd2, mqh1, mqh2, mqs1
     >    , mqs2, mqx1, mqx2, npmx, nqmx, psman, psprv, pssol
      use cplmet,      only : gare, gdsv, gwtm, gwtp, hare, hdsp, hdsv
     >    , hdxm, hdxp, hgdx, hgdy, hpit, hvol, hvsb, hwtm, hwtp, icaxs
     >    , icel, icmax, icmin, icmpe, icmps, icspx, icwl1, icwl2, itmax
     >    , itmpe, itmps, itpve, itpvs, itsle, itsls, jcdp1, jcdp2, jcel
     >    , jcmax, jcxp1, jcxp2, jnxm, jnxp, jtmax, jtmin, kce, kcn, kcs
     >    , kcw, kgdx, kgdy, kreg, nompl, noprv, nosol, romn, vlmn
      use csize,       only : ndp => ndy, ndq => ndx, ndx, ndy
      use cunit,       only : lmspe, lmype, lnope, mywld, n6
      use mod_dtypdef, only : nddt, typ_dnam, typ_itag
      use mpi!,         only : mpi_bcast
      implicit none
!
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  nft
!ik   character cact*(*)
      integer,   intent(in) :: nft
      character, intent(in) :: cact*(*)
!
!::local variables
      character  cver*80
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  nndq, nndp, nndx, nndy, ldbg, lenx
!ik   integer  it, jts, jte, jt, jc, ic, i
      integer  nndq, nndp, nndx, nndy, ldbg
      integer  jts, jte, jt, jc, ic, i
! modified 1/1 lines replace all include files with module files by kamata 2021/08/18
!ik   integer  nls, nle, nbf, ierr
      integer  ierr, itag, ityp
! added 2 lines organize local variables and include files by kamata 2021/05/31
! function
      integer    lenx

! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   integer  jmax, jw
!
!::version
      cver = "/home/shimizu/Gmesh/soc/gmtrc/mtrdsk   04/4/28"
!
      ldbg = 1
!
      write(n6,'(/2x,"*** mtrdsk ***   ",a)') cact(1:lenx(cact))
      write(n6,'(2x,"version : ",a)') cver(1:lenx(cver))

      if( lmype.ne.lmspe ) goto 1000
!
!----------------------------------------------------------------------
!::write
!----------------------------------------------------------------------
      if( cact(1:1).eq."w" ) then
      write(nft) cver
!
      write(nft) ndq, ndp
      write(nft)
     >   mqd1, mqx1, mqs1, mqs2, mqh1, mqh2, mqx2, mqd2, nqmx
     >  ,mpw1, mpsp, mpw2, mpax, npmx
     >  ,mpsol, mpprv, mpman
     >  ,grdx, grdy, hbr, hbz, hbt
     >  ,pssol, psprv, psman
!
      write(nft) ndx, ndy
      write(nft)
     >   nosol, noprv, nompl
     >  ,jcdp1, jcdp2, jcxp1, jcxp2, jcmax
     >  ,icwl1, icspx, icwl2, icmps, icmpe, icaxs
     >  ,itsls, itsle, itpvs, itpve, itmps, itmpe, itmax
     >  ,icmin, icmax, jtmin, jtmax
     >  ,kgdx,  kgdy,  kreg
     >  ,jcel,  icel,  jnxp, jnxm
     >  ,kce, kcw, kcn, kcs
     >  ,vlmn, romn
     >  ,hvol, hgdx, hgdy
     >  ,hdxp, hdxm
     >  ,hdsp, hdsv, hvsb
     >  ,hare, hpit
     >  ,hwtp, hwtm
     >  ,gdsv, gare
     >  ,gwtp, gwtm
      endif
!
!----------------------------------------------------------------------
!::read
!----------------------------------------------------------------------
      if( cact(1:1).eq."r" ) then
      read(nft) cver
      write(n6,'(2x,"version : ",a)') cver(1:lenx(cver))
!-----
      read(nft) nndq, nndp
!
      write(n6,'(2x,"dimension size  ndq =",2i5,"  ndp =",2i5)')
     >  ndq, nndq, ndp, nndp
      if( ndq.ne.nndq .or. ndp.ne.nndp ) then
      call wexit("mtrdsk","ndq.ne.nndq .or. ndq.ne.nndp")
      endif
!
      read(nft)
     >   mqd1, mqx1, mqs1, mqs2, mqh1, mqh2, mqx2, mqd2, nqmx
     >  ,mpw1, mpsp, mpw2, mpax, npmx
     >  ,mpsol, mpprv, mpman
     >  ,grdx, grdy, hbr, hbz, hbt
     >  ,pssol, psprv, psman
!
      read(nft) nndx, nndy
!
      write(n6,'(2x,"dimension size  ndx =",2i5,"  ndy =",2i5)')
     >  ndx, nndx, ndy, nndy
      if( ndx.ne.nndx .or. ndy.ne.nndy ) then
      call wexit("mtrdsk","ndx.ne.nndx .or. ndy.ne.nndy")
      endif
!
      read(nft)
     >   nosol, noprv, nompl
     >  ,jcdp1, jcdp2, jcxp1, jcxp2, jcmax
     >  ,icwl1, icspx, icwl2, icmps, icmpe, icaxs
     >  ,itsls, itsle, itpvs, itpve, itmps, itmpe, itmax
     >  ,icmin, icmax, jtmin, jtmax
     >  ,kgdx,  kgdy,  kreg
     >  ,jcel,  icel,  jnxp, jnxm
     >  ,kce, kcw, kcn, kcs
     >  ,vlmn, romn
     >  ,hvol, hgdx, hgdy
     >  ,hdxp, hdxm
     >  ,hdsp, hdsv, hvsb
     >  ,hare, hpit
     >  ,hwtp, hwtm
     >  ,gdsv, gare
     >  ,gwtp, gwtm
      endif
!
!----------------------------------------------------------------------
!::send ntgdsk
!----------------------------------------------------------------------
 1000 continue
      if( lnope.gt.1 ) then
      write(n6,'(5x,"passed MPI_Bcast")')
!
!::[MPI_Bcast in mtrdsk]   in cgdcom (grdx,emrk)  10/04/21
! modified 4/3 lines replace all include files with module files by kamata 2021/08/18
!ik   nls = loc( grdx(1,1) )
!ik   nle = loc( cgmsh_emrk ) + 1
!ik   nbf = nle - nls
!ik   call MPI_Bcast( grdx, nbf, MPI_BYTE, lmspe, mywld, ierr )
      call tbfind( nddt, typ_dnam, ityp, 'CGMSH' )
      itag = typ_itag(ityp)
! modified 1/2 lines dynamic allocation of arrays by kamata 2022/05/29
!ik   call MPI_Bcast( grdx, 1, itag, lmspe, mywld, ierr )
      call MPI_Bcast( mqd1, 1, itag, lmspe, mywld, ierr )
      call cgdcom_sr( 1, lmspe )
!
!::[MPI_Bcast in mtrdsk]   cplmet (vlmn,emrk)    10/04/21
! modified 4/3 lines replace all include files with module files by kamata 2021/08/18
!ik   nls = loc( vlmn(1) )
!ik   nle = loc( cmetrc_emrk ) + 1
!ik   nbf = nle - nls
!ik   call MPI_Bcast( vlmn, nbf, MPI_BYTE, lmspe, mywld, ierr )
      call tbfind( nddt, typ_dnam, ityp, 'CMETRC' )
      itag = typ_itag(ityp)
! modified 1/2 lines dynamic allocation of arrays by kamata 2022/05/29
!ik   call MPI_Bcast( vlmn, 1, itag, lmspe, mywld, ierr )
      call MPI_Bcast( nosol, 1, itag, lmspe, mywld, ierr )
      call cplmet_sr( 1, lmspe )
!
      endif
!
!----------------------------------------------------------------------
!::debug write
!----------------------------------------------------------------------
      write(n6,'(2x,"psi-mesh")')
      write(n6,'(5x,"mpw1 =",i3,"  mpsp =",i3,"  mpw2 =",i3,"  mpax =",
     >  i3,"  npmx =",i3,"  ndp =",i3)')
     >  mpw1, mpsp, mpw2, mpax, npmx, ndp
      write(n6,'(2x,"xi-mesh")')
      write(n6,'(5x,"mqd1 =",i3,"  mqx1 =",i3,"  mqs1 =",i3,"  mqs2 =",
     >  i3,"  mqx2 =",i3,"  mqd2 =",i3,"  nqmx =",i3,"  ndq =",i3)')
     >  mqd1, mqx1, mqs1, mqs2, mqx2, mqd2, nqmx, ndq
      write(n6,'(5x,"mqh1 =",i3,"  mqh2 =",i3)') mqh1, mqh2
      write(n6,'(2x,"sol  ",i3)') mpsol
      write(n6,'(5x,"pssol =",1p8e12.3)') (pssol(i),i=1,mpsol)
      write(n6,'(2x,"prv  ",i3)') mpprv
      write(n6,'(5x,"psprv =",1p8e12.3)') (psprv(i),i=1,mpprv)
      write(n6,'(2x,"main ",i3)') mpman
      write(n6,'(5x,"psman =",1p8e12.3)') (psman(i),i=1,mpman)
!
      write(n6,'(2x,"nosol =",i3,"  noprv =",i3,"  nompl =",i3)')
     >   nosol, noprv, nompl
      write(n6,'(2x,"itsls =",i3,"  itsle =",i3,"  itpvs =",i3
     >  ,"  itpve =",i3,"  itmps =",i3,"  itmpe =",i3)')
     >   itsls, itsle, itpvs, itpve, itmps, itmpe
      write(n6,'(2x,"icwl1 =",i3,"  icspx =",i3,"  icwl2 =",i3
     >  ,"  icmps =",i3,"  icmpe =",i3,"  icaxs =",i3)')
     >   icwl1, icspx, icwl2, icmps, icmpe, icaxs
      write(n6,'(2x,"jcdp1 =",i3,"  jcxp1 =",i3,"  jcxp2 =",i3
     >  ,"  jcdp2 =",i3)') jcdp1, jcxp1, jcxp2, jcdp2
!
      if( ldbg.eq.1 ) then
! modified 4/3 lines organize local variables and include files by kamata 2021/05/31
!ik   it  = itsle
!ik   jts = jtmin(it)
!ik   jte = jtmax(it)
!ik   write(n6,'(2x,"flux tube it =",i3)') it
      jts = jtmin(itsle)
      jte = jtmax(itsle)
      write(n6,'(2x,"flux tube it =",i3)') itsle
      write(n6,'(4x,"jt",2x,"jc",3x,"hgdx",4x,"hgdy",4x,"hvol",7x
     >  ,"hwtm",4x,"hwtp")')
      do jt = jts, jte
      if( jt.ge.jts+3 .and. jt.le.jte-3 ) cycle
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   jc = jcel(jt,it)
!ik   ic = icel(jt,it)
      jc = jcel(jt,itsle)
      ic = icel(jt,itsle)
      write(n6,'(2x,2i4,2f8.3,1pe11.2,0p2f8.3)')
     >  jt,jc,hgdx(jc,ic),hgdy(jc,ic),hvol(jc,ic)
     > ,hwtm(jc,ic),hwtp(jc,ic)
      enddo
      endif
!
!::check
      write(n6,'(2x)')
      write(n6,'(2x,"vlmn(ndx)  mpsp,mpax,ndx =",3i5)') mpsp,mpax,ndx
      write(n6,'(2x,"romn(ndy)  mpsp,mpax,ndy =",3i5)') mpsp,mpax,ndy
      write(n6,'(2x,"vlmn =",1p10e11.3)') (vlmn(i),i=mpsp,mpax)
      write(n6,'(2x,"romn =",1p10e11.3)') (romn(i),i=mpsp,mpax)
      write(n6,'(2x)')
!
      return
      end
