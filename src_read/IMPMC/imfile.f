!***********************************************************************
      subroutine imfile(nft)
!***********************************************************************
!
!
!
!
!-----------------------------------------------------------------------
      use cimcom, only : i6_epos, i6_ipos, i6_repo, i6_trac, i6_wrfi
     >, i6_wrfn, ip1tr, ip2tr, isw2tr, iswtr, lp1tr, lp2tr, mpetr
      use cimctl, only : nsvtm, svdsn, svitm
      use csonic, only : itend
      use cunit,  only : lmype, mype, n6
      implicit none
!
!::argument
! modified 1/1 lines organize local variables and include files by kamata 2021/06/28
!ik   integer  nft
      integer, intent(in) :: nft
!
!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/06/28
!ik   character  cline*80, csvtm*80, cno*8
      character  csvtm*80, cno*8
      integer  mjs(30), mje(30), kmx, ii, i, no, jnl
      logical  lbck
      character(30) :: dstrc, dspos
      character(5)  :: clmpe
!
      namelist /uimfil/
     >  csvtm,
     >  mpetr, ip1tr, ip2tr, lp1tr, lp2tr, iswtr, isw2tr,
     >  i6_repo, i6_trac, i6_ipos, i6_wrfn, i6_wrfi
!
!::default value
      i6_repo = 0          ! sub. impcal  dummy
      i6_wrfn = 0          ! sub. imwrfn  dummy
      i6_wrfi = 0          ! sub. imwrfi  dummy
      i6_trac = 0
      i6_ipos = 0
      i6_epos = 0
      dstrc = "-----"
      dspos = "-----"
!
      write(n6,'(/2x,"*** imfile ***")')
!
!::input
      read(nft,uimfil)
!
!::save directory
      write(n6,'(2x,"csvtm = [",a,"]")') trim(csvtm)
      call linsep( csvtm, " ,", kmx, mjs, mje, 30 )
      ii = 0
      do i = 1, kmx
      read( csvtm(mjs(i):mje(i)), * ) no
      if( no.gt.itend ) cycle
      ii = ii + 1
      svitm(ii) = no
      enddo
!
      if( ii.eq.0 ) then
      ii = ii + 1
      svitm(ii) = itend
      endif
      if( svitm(ii).ne.itend ) then
      ii = ii + 1
      svitm(ii) = itend
      endif
      nsvtm = ii
!
      lbck = .TRUE.
      do i = 1, nsvtm
      write(cno,'(i8)') svitm(i)
      jnl = index(cno,' ',lbck)
      svdsn(i) = "dIMP_"//cno(jnl+1:8)
      enddo
      write(n6,'(2x)')
      write(n6,'(2x,"svitm =",10i12)') (svitm(i),i=1,nsvtm)
      write(n6,'(2x,"svdsn =",10a12)') (svdsn(i),i=1,nsvtm)
!
!::trace
      write(clmpe,'(i5.3)') lmype
      if( i6_trac.gt.0 ) then
      i6_trac = 0
      iswtr = 0
      dstrc = "-----"
      if( lmype.eq.mpetr ) then
        i6_trac = 300000 + 61000 + mype
        iswtr = 1
        dstrc = "@trac_"//trim(clmpe)
        open(unit=i6_trac,file=dstrc)
!
      write(i6_trac,
! modified 1/1 lines with TOPICS by kamata 2021/12/22
!ik  >  '(3x,"comt",2x,"lpcm",3x,"ip",2x,"ptmx",6x,"ptim",8x,
     >  '(3x,"comt",7x,"lpcm",3x,"ip",2x,"ptmx",6x,"ptim",8x,
     >  "tt",10x,"pdtm",8x,"ic",2x,"ix",2x,"iy",4x,"ln",1x,"ien",2x,
     >  "ko",2x,"is",1x,"ml",3x,"rr",7x,"zz",6x,"rO",5x,"wght",3x,
     >  "vz",8x,"Evel",6x,"vflw",6x,"Ti",9x,"lstp",5x,"tauz",5x,
     >  "dt",7x,"Ftot",5x,"dVz",6x,"V")')

      endif
      endif
!
!::ionization point
      if( i6_ipos.gt.0 ) then
      i6_ipos = 0
      dspos = "-----"
      if( mype.eq.mpetr ) then
      i6_ipos = 300000 + 62000 + mype
      dspos = "@ipos_"//clmpe
!
      open(unit=i6_ipos,file=dspos)
      write(i6_ipos,
! modified 1/1 lines with TOPICS by kamata 2021/12/22
!ik  >  '(3x,"comt",1x,"lpcm",2x,"ip",2x,"ptmx",7x,"ptim",7x,
     >  '(3x,"comt",6x,"lpcm",2x,"ip",2x,"ptmx",7x,"ptim",7x,
     >  "tt",9x,"pdtm",9x,"ic",2x,"ix",2x,"iy",4x,"ln",1x,"ien",2x,
     >  "ko",2x,"is",1x,"ml",3x,"rr",7x,"zz",6x,"rO",5x,"wght",3x,
     >  "vz",8x,"Evel",6x,"vflw",6x,"Ti",8x,"lstp",6x,"tauz")')
      endif
      endif

!::debug
      write(n6,'(2x,"trace ",a,2x,"i6_trac =",i6,"  mpetr =",i3,
     >  "  iptr =",2i6,"  lptr =",2i4,"  iswtr,iswtr2 =",2i3)')
     >  trim(dstrc), i6_trac, mpetr, ip1tr, ip2tr, lp1tr, lp2tr,
     >  iswtr, isw2tr
      write(n6,'(2x,"mype/lmype =",2i4," trc : ",a,"  pos : ",a)')
     >    mype, lmype, trim(dstrc), trim(dspos)
!
!::debug for imwrfn & imwrfi
!
      return
      end
