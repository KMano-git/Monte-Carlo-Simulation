!***********************************************************************
!
!     comment out   10 lines
!
!     setenv BTYPE  big  ! little
!     setenv BTYPE <rtn>
!
!x      call test_nopen
!x      stop
!x      end
!
!x      subroutine wexit
!x      return
!x      end
!
!***********************************************************************
      subroutine test_nopen
!***********************************************************************
      use cunit, only : lmspe, lmype, n6
      implicit none
!
      integer    nft
      character  cdsn*80, cinp*80
      logical    lex
!
!::local variables
      real*8   time, x(10), y(10)
      integer  itim
      integer  n
!
      integer  i
!
      lmype = 0
      lmspe = 0
      n6 = 6
      nft = 21
!
      write(n6,'(2x,"setenv BTYPE  big (/little/<rtn>)")')
      write(n6,'(">> enter (w/r) ==> ",$)')
      read(5,*) cinp
!
!---------------------------------------------------------------------
      if( cinp(1:1).eq."w" ) then
!---------------------------------------------------------------------
!
!::data
      time = 25.5e-3
      itim = 12000
      n = 10
      do i = 1, n
      x(i) = 0.1*dfloat(i-1)
      y(i) = x(i)**2
      enddo
!
!::data
      write(n6,'(2x,50("-"))')
      write(n6,'(2x,"time =",1pe12.4,"  itim =",i6,"  n =",i4)')
     >   time, itim, n
      write(n6,'(2x,"x =",1p10e12.3)') x(1:10)
      write(n6,'(2x,"y =",1p10e12.3)') y(1:10)
      write(n6,'(2x,50("-"))')
!
!:write
      call nopen(nft, "MTRC", "write,binary", cdsn, lex)
      write(nft) time, itim, n, x, y
      close(nft)
!
!---------------------------------------------------------------------
      elseif( cinp(1:1).eq."r" ) then
!---------------------------------------------------------------------
!
!::read
      call nopen(nft, "MTRC", "read,binary", cdsn, lex)
      read(nft) time, itim, n, x, y
      close(nft)
!
      write(n6,'(2x,"time =",1pe12.4,"  itim =",i6,"  n =",i4)')
     >   time, itim, n
      write(n6,'(2x,"x =",1p10e12.3)') x(1:10)
      write(n6,'(2x,"y =",1p10e12.3)') y(1:10)
!
!---------------------------------------------------------------------
!::error
!---------------------------------------------------------------------
      else
      write(n6,'(2x,"invalid input")')
      endif
!
      return
      end
!
!***********************************************************************
      subroutine nopen(nft,ckey,cfrm,cdsn,lex)
!***********************************************************************
!
!   alpha station version
!-----------------------------------------------------------------------
!     cfrm = 'text'   : ascii character
!            'binary' : binary
!
!         big endian ( ibm, sun, vpp, origin )
!         little endian binary( alpha station & pc & intel ) standard
!
!   SGI version  (Error in open statement)
!-----------------------------------------------------------------------
!         open( unit=nft, ... convert="..." )
!
!-----------------------------------------------------------------------
!   mplbry
!      common variables   cunit :  lmype,lmspe,n6
!      procedure          wexit
!-----------------------------------------------------------------------
!   shell
!     setenv BTYPE  big   !  little
!     call   getenv("BTYPE",cbty)
!-----------------------------------------------------------------------
!
!   nopen( unit,  ckey,    cfrm,   cdsn,  lex )
!          21    "METRIC", "read,binary,append", "exe/run05/mtrc", "T"
!                "exe/run05/mtxt", "read,text",  "exe/run05/mtxt", "T"
!
!-----------------------------------------------------------------------
      use cunit,     only : lmspe, lmype, mype, n6
      use mod_shexe, only : use_namelist, BTYPE
      implicit none
!
!::argument
! modified 3/4 lines organize local variables and include files by kamata 2021/05/31
!ik   integer    nft
!ik   character  ckey*(*), cdsn*(*), cfrm*(*)
!ik   logical    lex
      integer,   intent(in)  :: nft
      character, intent(in)  :: ckey*(*), cfrm*(*)
      character, intent(out) :: cdsn*(*)
      logical,   intent(out) :: lex
!
!::local variables
      character*(20) ::  cfmt, cact, cpos
      character*(20) ::  cbty, ccnv, ctrm
!
      integer    ist
      data       ist/0/
      save       ist, ccnv
!
!----------------------------------------------------------------------
!::endian type   all for binary file   initial set
!----------------------------------------------------------------------
      if( ist.eq.0 ) then
      ist = 1
      if(use_namelist) then
         cbty = BTYPE
      else
         call getenv("BTYPE",cbty)
      end if

      ccnv = " "
      if( index(cbty,"big").gt.0 ) ccnv ="big_endian"
      if( index(cbty,"little").gt.0 ) ccnv ="little_endian"
      write(n6,'(/2x,"*** nopen ***  type of binary data BTYPE ",a,
     >  "  ccnv ",a)') "["//trim(cbty)//"]", "["//trim(ccnv)//"]"
      endif
!
!----------------------------------------------------------------------
!::file name
!----------------------------------------------------------------------
      if(use_namelist) then
         cdsn = ckey
      else
         call getenv( ckey, cdsn )
!.. if ckey is not environment variable, get file name from ckey
         if(len_trim(cdsn).le.0) cdsn = ckey
      end if
!
!
!----------------------------------------------------------------------
!::MPI
!----------------------------------------------------------------------
!::slave PE
      if( lmype.ne.lmspe ) then
        return
      endif
!
!----------------------------------------------------------------------
!::condition
!----------------------------------------------------------------------
!::action
      cact = "read"
      if( index(cfrm,"read") .gt.0  ) cact = "read"
      if( index(cfrm,"write").gt.0 )  cact = "write"
      if( index(cfrm,"read") .gt.0 .and.
     >    index(cfrm,"write").gt.0 )  cact = "readwrite"
!
!::formatt
      cfmt = "formatted"
      if( index(cfrm,"text").gt.0 )   cfmt = "formatted"
      if( index(cfrm,"binary").gt.0 ) cfmt = "unformatted"
!
!::position
      cpos = "asis"
      if( index(cfrm,"rewind").gt.0 ) cpos = "rewind"
      if( index(cfrm,"append").gt.0 ) cpos = "append"
      if( index(cfrm,"asis")  .gt.0 ) cpos = "asis"
!
!----------------------------------------------------------------------
!::check of file
!----------------------------------------------------------------------
!::read
      if( cact(1:1).eq."r" ) then
        if( len_trim(cdsn) .gt. 0 ) then
          inquire( file=cdsn, exist=lex )
        else
          cdsn = '_'
          lex = .false.
        endif
!
!::not find file
        if( .not.lex ) then
          write(n6,'(2x,"=== file allocation error (nopen) ===")')
          write(n6,'(2x,"nft =",i3,"  ckey =",a,"  cfrm =",a,
     >      "  cdsn =",a)') nft,trim(ckey),trim(cfrm),trim(cdsn)
          call wexit("nopen","not find file : "//trim(cdsn))
        endif
      endif
!
      lex = .true.
!
!----------------------------------------------------------------------
!::read/write
!----------------------------------------------------------------------
      if( cact(1:1).eq."r" ) ctrm = "R-cdsn"
      if( cact(1:1).eq."w" ) ctrm = "W-cdsn"
!
!::open of text file
      if( cfmt(1:1).eq."f" ) then
      write(n6,602) "=== open file (text) ",
     >  mype, nft, trim(ckey), trim(ctrm), trim(cdsn), trim(cpos)
      open( unit=nft,file=cdsn,form=cfmt,action=cact,position=cpos )
      else
!
!::open of binary file
      write(n6,602) "=== open file(binary)",
     >  mype, nft, trim(ckey), trim(ctrm), trim(cdsn), trim(cpos)
      if( len_trim(ccnv).eq.0 ) then
        open( unit=nft,file=cdsn,form=cfmt,action=cact,position=cpos )
        write(n6,'(a)') '  (NO convert endian)'
      else
        open( unit=nft,file=cdsn,form=cfmt,action=cact,position=cpos
     >       ,convert=ccnv )
        write(n6,'(a)') "  (W convert "//trim(ccnv)//")"
      endif
!
      endif
!
      return
!
!
!::format
 602  format(2x,a,"  mype =",i5,"  nft =",i3,"  ckey =",a,2x,
     >  a,2x,a,2x,a)
      end
