!**********************************************************************
      subroutine ntmont
!**********************************************************************
!
!        neutral transport code
!
!::2011/10/31
!        nhunt = 3,  nhfl = 5      avarage (11-15)
!        nhunt = 1,  nhfl = 20     average 20
!        nhunt = 0,  nhfl = 1      average  1
!
!----------------------------------------------------------------------
      use BGProf, only : ntbgprof, ntwprof, swnrmForEachSource ! monte/ntbgprf.f90
      use cntcom, only : cstyp, flxin, iptl, lnnel, lscrg, ltrc
     > , xpos, ypos, weit, icpo
     > , is_atom_or_mole,ntmont_flg
      use cntctl, only : mdl_ntsm, nclnt, nhfl, nhunt
      use cntmnt, only : dotn, iseed, mcflx, mfmax, mnknd, mnsmp, mnvwk
     >    , nsmpmx, pfl_ion, pfl_ntl, vcsrc, vflux, visrc, vitim, vitnt
     >    , vkflx, vmwork, vnsmp, vnsrc, vsmno, vsmty, vtime
      use cntpfctl, only : fpfctl, lpfctl
      use cntwcn, only : swnrm, wabs, wcsrc, werr, wion, wisrc, wsum
     >    , wtot
      use cntxwk, only : ihfl, jhfl, ndhfl, ndhun, nhno, nhtm, xmwork
      use cplcom, only : tfps
      use cputim, only : cntcal, cntsum
      use csize,  only : ndmsr, nwkmp_nt
      use csonic, only : itim, lrand, time
      use cunit,  only : lmceq, lmspe, lmype, lnope, mywld, n6
      use mpi!,    only : mpi_bcast, mpi_real8, mpi_reduce, mpi_sum
!    >    , mpi_wtime
      use ntpfctl, only : pfctl
      use mod_externalgrid, only : use_exdata
      implicit none
!
!::local vaiables
      real*8  wmwork(nwkmp_nt)
      integer ido1, ido2, ido3, mpsum
      integer i, m, ii, jseed, ia
!
!::mpi variables
      integer  ierr
!
!::local variables
      character  csrc*6
      integer ih, jh, iflx, isrc, nknd, nsmp
      real*8  tfion, tfntl, zsum
      real(8) :: cptm1, cptm2, cptm3, cptm4
!::for ntfolw 
      real(8) :: zincx(0:1), zint1(0:1), zint2(0:1), ztim1(0:1)
     > , ztim2(0:1)

! function
      real(8)    rnlast

      call ran_init
!
!::number of calculation
      nclnt = nclnt + 1
!
!::dimension check
      if( nhfl.gt.ndhfl .or. nhunt.gt.ndhun ) then
        write(n6,'(/2x,"dimension err at sub. ntmont ")')
        write(n6,'(2x,"nhfl,ndhfl =",2i4,"  nhunt,ndhun =",2i5)')
     >   nhfl, ndhfl, nhunt, ndhun
        call wexit( 'ntmont', 'dimension err' )
      endif
!
!::smothing
      vsmty = 1  ! when vsmty = 2, set dtim = 2.0d-7
      vsmno = 1  ! number of smothing (max = nhfl*nhunt = 3*5 )
!
      if( nhunt.eq.0 ) then
        nhfl = 1
        ihfl = 1
        jhfl = 1
      else
        ihfl = mod((nclnt-1)/nhunt,nhfl)+1
        jhfl = mod(nclnt-1,nhunt)+1
      endif
!
      if( jhfl.eq.1 ) then
        nhno(ihfl) = 0
        do i = 1, nwkmp_nt
          do m = 1, ndmsr
            xmwork(i,m,ihfl) = 0.0d0
          enddo
        enddo
      endif
      nhno(ihfl) = nhno(ihfl) + 1
      nhtm(nhno(ihfl),ihfl) = itim
!
      vsmty = 1
      if( jhfl.eq.1 )  vsmty = 2
      if( nhunt.le.1 ) vsmty = 1
!
      ii = 0
      do ih = 1, nhfl
        ii = ii + nhno(ih)
      enddo
      vsmno = ii
!
!-----
      write(n6,'(/2x,"*** itim =",i6,"  vsmty =",i2,"  vsmno =",i5)')
     >    itim, vsmty, vsmno
      do ih = 1, nhfl
        write(n6,'(2x,"xmwork  itim =",i5,"  ihfl =",i4,"  nhno =",i4,
     >    :"  nhtm =",20i7)') itim,ih,nhno(ih),
     >   (nhtm(jh,ih),jh=1,nhno(ih))
      enddo
!-----
!
!
!::ion & neutral flux
      write(n6,'(/2x,"=== ntmont ===  itim =",i5,"  vsmty =",i2,
     >     "  vsmno =",i3,"  nhfl =",i2,"  ran =",0pf12.9)')
     >   itim, vsmty, vsmno, nhfl, rnlast(0)
      write(n6,'(3x,"isr",2x,"ifl",2x,"name",5x,"dotn",8x,"wtot",10x,
     > "w(ion+abs)",4x,"werr",10x,"wion",10x,"wabs")')
!
      call ntpump(1)
      call ntbgprof   ! 2015.10.06, by toku
!
      isrc = 0
      do iflx = 1, mfmax
        if(iflx==6) then
          fpfctl = 1.0d0
          if(lpfctl>0) call pfctl ! set fpfctl
        endif

        call ntnflx(iflx,tfion,tfntl)
!
        if( mnsmp(iflx).le.0 ) tfntl = 0.0d0
        pfl_ion(iflx) = tfion
        pfl_ntl(iflx) = tfntl
        mnvwk(iflx) = 0
        visrc(iflx) = 0
        if( tfntl.le.0.0d0 ) cycle
        isrc = isrc + 1
        if( isrc.gt.ndmsr ) then
          call wexit("ntmont","isrc.gt.ndmsr")
        endif
        mnvwk(iflx) = isrc
        visrc(iflx) = isrc
        vkflx(isrc) = iflx
!
!::loop ksrc
        cptm1 = MPI_WTIME()
        csrc = mcflx(iflx)
        nknd = mnknd(iflx)
        nsmp = mnsmp(iflx)
        flxin = pfl_ntl(iflx)
        cstyp = csrc
!
        wisrc = isrc
        wcsrc = csrc
!
!::save source type
        vnsrc = isrc
        vcsrc(isrc) = csrc
        vnsmp(isrc) = nsmp
        vflux(isrc) = flxin
!
        if( lrand.eq.1 .and. nsmp.gt.nsmpmx) then
          call wexit("ntmont","nsmp.gt.nsmpmx")
        endif
!
!::MPI
        ido1 = 1; ido2 = nsmp; ido3 = 1; mpsum = 0
        if( lmceq(2).eq.0 .and. lnope.gt.1 ) then
          ido1 = lmype+1; ido2 = nsmp; ido3 = lnope;  mpsum = 1
        endif
!
!::clear scoreing variables
        call ntcler   ! set dotn
!
!::surface/volume source
        do i = ido1, ido2, ido3
          is_atom_or_mole = 0 ! atom
          ntmont_flg = 0
          iptl = i
          if(ltrc.gt.0) then
            call ntrc680(0,0.0d0,0.0d0,0.0d0,0,"ntrc_"//csrc)
          endif
          if( lrand.eq.1 ) then
            jseed = iseed(i)
            call srandom(jseed)
          endif
          if( nknd.ne.3  ) then
            call ntstaw
          else
            call ntstav
          endif
!::
          ! initialize value
          zincx = 0.0d0
          zint1 = 0.0d0
          zint2 = 0.0d0
          ztim1 = 0.0d0
          ztim2 = 0.0d0
          do while(.true.)
          !------------------------
          ! folow sample particle
          !------------------------
            call ntfolw(zincx,zint1,zint2,ztim1,ztim2)
          !------------------------
            if(is_atom_or_mole==0) then
              exit
            else ! molecular
              ! from ntstaw d-plate or puff
              if(ntmont_flg == 1) then
                call ntrc680( 1,xpos,ypos,weit,icpo,"staw=end" )
                is_atom_or_mole = 0
                ntmont_flg = 0
              ! return from ntrefl in hiwall_atom
              elseif(ntmont_flg == 2) then
                ntmont_flg = 3
              ! done ntrefl mole calculation and back to ntfolw as atom
              elseif(ntmont_flg == 3) then
                is_atom_or_mole = 0
              endif
            endif
!::
          enddo ! ntfolw loop
        enddo
!
!::track / collision estimator
        call ntscpy
        call ntwcon
!
!::MPIsumup
!::[MPI_Reduce in ntmont]  cntwcn/cntvar/ (wflx,werr)  10/04/21
!::[MPI_Bcast  in ntmont]  cntwcn/cntvar/ (wflx,werr)  10/04/21
        cptm3 = MPI_WTIME()
        if( mpsum.eq.1 ) then
          call setvmwork( 2, isrc )
          do i= 1,  nwkmp_nt
            wmwork(i) = vmwork(i,isrc)
          enddo
          call MPI_Reduce( wmwork, vmwork(1,isrc), nwkmp_nt
     >                , MPI_REAL8, MPI_SUM, lmspe, mywld, ierr )
          call MPI_Bcast( vmwork(1,isrc), nwkmp_nt
     >                , MPI_REAL8, lmspe,mywld,ierr )
        else
          call setvmwork( 2, isrc )
        endif
        cptm4 = MPI_WTIME()
        cntsum(isrc) = cptm4-cptm3
        cntsum(0) = cntsum(0) + cntsum(isrc)
!
!::check
        call setvmwork( 1, isrc )
!
        if(.not.use_exdata) then
          if( mdl_ntsm.eq.1 ) call ntsmsrc(1)   ! smoothing
        endif
        call ntwcon
!
!::replace
        call setvmwork( 2, isrc )
!
!::source for each reaction
        if(lscrg.eq.1) then
          call ntscrg   ! check neutral source 2010/12/28
        endif
!
!::smothing
        do i = 1, nwkmp_nt
          xmwork(i,isrc,ihfl) = xmwork(i,isrc,ihfl) + vmwork(i,isrc)
        enddo
        do i = 1, nwkmp_nt
          zsum = 0.0d0
          do ih = 1, nhfl
            zsum = zsum + xmwork(i,isrc,ih)
          enddo
          vmwork(i,isrc) = zsum
        enddo
!
!::scoreing variables
        call setvmwork( 1, isrc )
        call ntwcon
        swnrmForEachSource(isrc) = swnrm
!
!::debug write
        write(n6,'(4x,i2,3x,i2,3x,a,1pe12.4,1p5e14.6,
     >    "  ran =",0pf12.9)')
     >    isrc,iflx,csrc,dotn,wtot,wsum,werr,wion,wabs, rnlast(0)
!
!::pump
        call ntpump(2)
!
!::loop
        cptm2 = MPI_WTIME()
        cntcal(isrc) = cptm2-cptm1
        cntcal(0) = cntcal(0) + cntcal(isrc)
      enddo  ! loop(iflx)
!
      call ntwprof   ! 2015.10.06, by toku
!
!::volume source ( H+ ==> H0 )
      if( mfmax.ge.8 ) then
        mnvwk(8) = mnvwk(7)
      endif
!
!::debug write
      ia = 1
      write(n6,'(2x)')
      write(n6,'(3x,"isr",2x,"ifl",2x,"name",3x,"knd",4x,"nsmp",3x,
     >  "fion",11x,"fntl",11x,"diff",8x,"plpflx")')
      do m = 1, mfmax
        write(n6,'(4x,i2,3x,i2,3x,a,i3,i8,1p2e15.6,1pe12.3,1pe15.6)')
     >  mnvwk(m), m, mcflx(m), mnknd(m), mnsmp(m), pfl_ion(m),
     >  pfl_ntl(m), pfl_ion(m)-pfl_ntl(m), tfps(m,ia)
      enddo
      if( vnsrc.gt.ndmsr .or. vnsrc.le.0 ) then
        call wexit("ntmont","dimension error  vnsrc > ndmsr")
      endif
!
      call ntpump(3)
!
!::time
      vtime = time
      vitim = itim
      vitnt = vitnt + 1
!
!::debug write
      write(n6,'(2x,"=== plsorc in ntmont ===  itim =",i6,"  time =",
     >  1pe14.6,"  loop =",i5)') vitim,vtime,vitnt
      write(n6,'(2x)')
!
      if( ltrc>0 ) then
        if(lnnel==1)then
          if(itim>0) call wexit("ntmont","ltrc = 1(lnnel==1)")
        else
          call wexit("ntmont","ltrc = 1")
        endif
      endif
!
      return
      end