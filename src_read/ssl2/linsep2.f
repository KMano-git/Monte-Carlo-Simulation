!%%mkspt +
!***********************************************************************
      subroutine test_linsep
!***********************************************************************
      implicit none
! modified 4/2 lines organize local variables and include files by kamata 2021/05/31
!ik   integer mjs(20),mje(20),ne,nk
!ik   integer i, ii, j, num, itab(50)
!ik   character csep*5, ctab(50)*20
!ik   data  csep/" =,();"/
      integer mjs(20),mje(20),nk
      integer i
!
      character ctxt(8)*120, clin*120
      ctxt(1) = "xs_data_base = 0, 101, 2702, 5303, 7904, 10505,"
      ctxt(2) = "   23315, 23316, 23317, 23318, 23319 ; "
      ctxt(3) = "xs_var = "
      ctxt(4) = '"cross_section   "'
      ctxt(5) = '"specific_energy  ",'
      ctxt(6) = '"unknown     " ; '
      ctxt(7) = 'xs_unit = "cm^2", "    ", "ev/amu" ;'
      ctxt(7) = '  "       ",  '
      ctxt(7) = 'cross_section, I_1_0, sigv_max'
      ctxt(8) = '* formatI  (2i5)'
!
!xx   call cdf_geti("xs_data_base",num,itab,50)
!
      clin = ctxt(8)
!xx   call linsep(clin,csep,nk,mjs,mje,20)
      call linsep(clin," ,",nk,mjs,mje,20)
      write(6,'(2x,"nk =",i3)') nk
      write(6,'(10(2x,"[",a,"]"))') (clin(mjs(i):mje(i)),i=1,nk)

!
      stop
      end
!
!***********************************************************************
      subroutine linsep(ckey,csep,nk,mjs,mje,ndm)
!***********************************************************************
      implicit none
!
!::arguments
! modified 2/3 lines organize local variables and include files by kamata 2021/05/31
!ik   character ckey*(*), csep*(*)
!ik   integer   ndm, nk, mjs(ndm), mje(ndm)
      character, intent(in)  :: ckey*(*), csep*(*)
      integer,   intent(in)  :: ndm
      integer,   intent(out) :: nk, mjs(ndm), mje(ndm)
!
!::local varaibales
      integer nmax, ii, js, je, i, is
      character ctab*1
!
      ctab = char(9)
      nmax = len_trim(ckey)
!
!xx   write(6,*) "ckey = [",trim(ckey),"]","  csep = [",trim(csep),"]"
!
      ii = 0
!
!::start
      is = 0
      do i = 1, nmax
      is = i
      if( ckey(i:i).eq.ctab ) cycle
      if( ckey(i:i).ne." " )  exit
      enddo
      if( is.eq.0 ) goto 100
!
!::separateor
      if( index(csep,ckey(i:i)).gt.0 ) is = is+1
      js = is
      do i = is, nmax
      if( index(csep,ckey(i:i)).gt.0 ) then
        ii = ii + 1
        if( ii.gt.ndm ) goto 910
        je = i - 1
        je = len_trim(ckey(1:je))
!
!--Note.    xs_units = "      ",  ==> " "
      if( ckey(je:je).eq.'"' .and. ckey(i:i).eq.'"' ) je = js+ 1
!--
        mjs(ii) = js
        mje(ii) = je
         if( js.gt.je ) then
            ii = ii - 1
         endif
        js = i+1
      endif
      enddo
!
      if( js.le.nmax ) then
      ii = ii + 1
      if( ii.gt.ndm ) goto 910
      mjs(ii) = js
      mje(ii) = nmax
      endif
!
 100  continue
      nk = ii
!
!::debug write
!xx      write(6,'(2x,"[",a,"]")') trim(ckey)
!xx      write(6,'(2x,10(i2,":","/"a,"/",2x))')
!xx     >  (i,ckey(mjs(i):mje(i)),i=1,nk)
!
      return
!
 910  continue
      write(6,'(/2x,"*** linsep ***  too many word  ",2i5)') ii,ndm
      call wexit("linsep2","too many word")
      end
