!**********************************************************************
!
!       neut2d : 2d neutral transport code
!
!       multi-source
!
!     ntl_ini   read input data
!     ntl_cal   calculation
! XX  ntl_src   calculation source terms  ! cal. source in plsorc
!     ntl_out   output
!     ntl_end   epilog
!
!**********************************************************************
      subroutine ntl_ini
!**********************************************************************
      use cntcom,      only : cftnr, ichti, ichto, iplx, iply, lntmd
     >    , mcel, nftnr
      use cntctl,      only : nclnt, itmnt, itmntb, nhfl, nhunt, tmntl
      use cntmnt,      only : iseed, mfmax, mnsmp, nsmpmx, vcsrc, vflux
     >    , visrc, vitim, vitnt, vkflx, vksty, vmwork, vnsmp, vnsrc
     >    , vsmno, vsmty, vtime, cntmnt_emrk
      use cntxwk,      only : ndhfl, ndhun, nhno, nhtm, xmwork
      use cplmet,      only : icel, itsle, jcel, jtmax
      use csize,       only : ndmfl, ndmsr, ndptl, nwkmp_nt
      use csonic,      only : lrand
      use cunit,       only : mygrp, mype, n5, n6, stcpy
      use mod_shexe,   only : use_namelist, INP_NTL, MESH_NTL
      use mpi!,         only : mpi_bcast, mpi_real8
      implicit none
!
!::local variables
      integer nft, i, nszv, nszx
      integer jt, ic, j, m, l
      integer, save :: ntset = 0
      logical lex
      character cdsn*80
!
      call trmark("ntl_ini","start")
!
!::avoid error at atm_copy (nxs>ndxs=20)
!   when  call ntl_ini 2 times
      ntset = ntset + 1
      write(n6,'(/2x,"=== passed ntl_ini  mype =",i5,"  mygrp =",i2
     >  "  ntset =",i2)') mype, mygrp, ntset
      if( ntset > 1 ) then
      write(n6,'(2x,"return because of ntset > 1")')
      return
      endif
!
!-----------------------------------------------------------------------
!
!::input data
      call trmark("ntl_ini","PRE input")
      nft = 21
      if(use_namelist) then
!..   use variable in inppls
         call ncopy(nft,INP_NTL,"text",cdsn,lex)
      else
!..   use environment variable
         call ncopy(nft,"INP_NTL","text",cdsn,lex)
      end if
      open(unit=n5,file=cdsn)
      call ntinpt(n5)
!
!::set table data
      call ntcnst
!
!::grid data
      nft = 21
      if(use_namelist) then
!..   use variable in inppls
         call nopen(nft,MESH_NTL,"binary",cdsn,lex)
      else
!..   use environment variable
         call nopen(nft,"MESH_NTL","binary",cdsn,lex)
      end if
      call ntgdsk(nft,"read")
      close (nft)
      call ntwlty  ! <== redefine wall type
      call ntmrgn  ! <== tube(itmpe) out of system
!
!::recycling data
      call ntwinp(n5)
      close (n5, status=stcpy)
!
!::recycling and temeperatute
      call ntwset
      call ntwlst
!
!::clear loop
      nclnt = 0
      itmnt = 0
      itmntb = 0
      tmntl = -1.0d0
      vitnt = 0
!
!::plasma surface
      call ntwpls
!
!::gass puff
      call ntpuff
!
!::Analytical model
      if( lntmd.eq.0 ) return
!
!::set degree of wall element
      call ntwdeg
!
!::EDF on wall
      call ntwedf_init
!
!::clear scoreing variables
      nszv = nwkmp_nt*ndmsr
      nszx = nwkmp_nt*ndmsr*ndhfl
      call setvmwork( 3, 0 )
      call setd( vmwork,  nszv, 0.0d0 )
      call setd( xmwork,  nszx, 0.0d0 )
!
!::clean 2012/01/19
      nhno(1:ndhfl) = 0
      nhtm(1:ndhun,1:ndhfl) = 0
!
!::clear in cntmnt
      vtime = 0.0d0
      vflux(1:ndmfl) = 0.0d0
      vitim = 0
      vitnt = 0
      vnsrc = 0
      vsmty = 0
      vsmno = 0
      visrc(1:ndmfl) = 0
      vkflx(1:ndmfl) = 0
      vksty(1:ndmfl) = 0
      vnsmp(1:ndmfl) = 0
      vcsrc(1:ndmfl) = "      "
      cntmnt_emrk = "#"

!::random seed for lrand = 1
      if( lrand.eq.1 ) then
      nsmpmx = 0
      do m = 1, mfmax
      nsmpmx = max0( nsmpmx, mnsmp(m) )
      enddo
      write(n6,'(/2x,"Note. random seed for lrand =1   nsmpmx =",i7,
     >  "  ndptl =",i7)') nsmpmx,ndptl
      if( nsmpmx.gt.ndptl )
     >   call wexit("ntl_ini","if lrand = 1, nsamp < ndptl")
      call rnseed(nsmpmx,iseed)
      endif
!
      if( nftnr.gt.0 ) then
      call ntdisk(nftnr,2,cftnr)
      endif
!
!::index check
      call ntindx
!
!::strike point (not dummy cell)
      do l = 1, 2
      if( l.eq.1 ) jt = 2
      if( l.eq.2 ) jt = jtmax(itsle) - 1
      j  = jcel(jt,itsle)
      i  = icel(jt,itsle)
      ic = mcel(j,i)
      if( l.eq.1 ) ichto = ic
      if( l.eq.2 ) ichti = ic
      write(n6,'(2x,"strike point  ",i2,"  ic =",i6,"  ix,iy =",2i5,
     >  "  ne,te,ti =",1p3e11.3)') l,ic,iplx(ic),iply(ic)
      enddo
!
!::input data for smothing results
      if( nhunt.gt.1 ) then
      write(n6,'(/2x,"smothing results for neut2d")')
      write(n6,'(2x,"nhunt =",i3,"  ndhun =",i3,"  nhfl =",i3,
     >  "  ndhfl =",i3)') nhunt, ndhun, nhfl, ndhfl
      if( nhfl.gt.ndhfl )  call wexit("ntl_ini","nhfl.gt.ndhfl")
      if( nhunt.gt.ndhun ) call wexit("ntl_ini","nhunt.gt.ndhun")
      call seti( nhno, ndhfl, 0 )
      endif
!
!::spatially smoothing
      call ntsmset
!
!::clear particle variables
      call ntclnpv
!
 100  continue
      call trmark("ntl_ini","return")

      return
      end
