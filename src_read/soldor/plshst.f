!***********************************************************************
      subroutine plshst
!***********************************************************************
      use cplcom, only : cfth, nfth
      use cplhst, only : cvnm, hstt, hstv, iphs, jphs, nhpn, nhtm, nhvl
      use cunit,  only : n6
      implicit none
!
!::local variables
! modified 2/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  nfw, nft, irec, ll, ih, k, jtm
!ik   character  cft*80
      integer  nfw, irec, ll, ih, k, jtm
      logical  lex
!
      if( nfth.le.0 ) return
!
      write(n6,'(/2x,"*** plshst ***")')
      nfw = 21
      open(unit=nfw,file="cmphst")
      write(nfw,'(2x,"itim",2x,"time",9x,"dtim",11x,
     >   "nedO",10x,"vadO",10x,"TidO",10x,"TedO",10x,
     >   "nedI",10x,"vadI",10x,"TidI",10x,"TedI",10x,
     >   "Finp",10x,"DotN",10x,"Fabs",10x,"Fpmp")')
!
!::initialization of file
! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik   nft = nfth
!ik   cft = cfth
      irec = 0
      jtm = 0
!
! modified 4/4 lines organize local variables and include files by kamata 2021/05/31
!ik   inquire(file=cft, exist=lex)
!ik   open(unit=nft,file=cft,form='unformatted')
!ik   write(n6,'(5x,"hist-file exist lex =",l2,"  nft =",i3,
!ik  >  "  file =",a)') lex,nft,cft(1:20)
      inquire(file=cfth, exist=lex)
      open(unit=nfth,file=cfth,form='unformatted')
      write(n6,'(5x,"hist-file exist lex =",l2,"  nfth =",i3,
     >  "  file =",a)') lex,nfth,trim( cfth )
      if( .not.lex ) then
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   write(n6,'(2x,"no found file  ",a)') cft(1:20)
      write(n6,'(2x,"no found file  ",a)') trim( cfth )
      return
      endif
!
      do ll = 1,10000
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   read(nft,end=215,err=212)
      read(nfth,end=215,err=212)
     >   cvnm,iphs,jphs,nhpn,hstt,hstv,nhtm,nhvl
      irec = irec + 1
      if( irec.lt.2 )
     >  write(n6,'(5x,"irec =",i3,"  nhtm =",i6)') irec,nhtm
!
!::output
      do ih = 1, nhtm
      jtm = jtm + 1
      write(nfw,'(i6,2f13.7,15f14.4)')
     >  jtm, hstv(ih,1)*1.0d3, hstv(ih,2)*1.0d3,
     >  hstv(ih, 9)/1.0d19,hstv(ih,10)/1.0d4, hstv(ih,11),hstv(ih,12),
     >  hstv(ih,13)/1.0d19,hstv(ih,14)/1.0d4, hstv(ih,15),hstv(ih,16),
     >  (hstv(ih,k)/1.0d22,k=57,60)
      enddo
!
      enddo
 215  continue
      write(n6,'(5x,"irec =",i3,"  nhtm =",i6)') irec,nhtm
!
      return
!
 212  continue
      call wexit("plshst","error occured during read file")
      end
