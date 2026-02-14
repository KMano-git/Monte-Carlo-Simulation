!***********************************************************************
      subroutine rflini(ctyp)
!***********************************************************************
      use cntrfl,      only : dlge, mpe, npe, rfelg, rfeng, rfer, rfrn
      use mod_shexe,   only : use_namelist, DEGAS_REFL
      use cunit,       only : lmspe, lmype, n6
      implicit none
!
!::argument
      character, intent(in) :: ctyp*(*)
!
!::local variables
      character clin*120,dsn*80
      integer   nft, je, kon, js, i, j
      real*8    zdel, emin, emax
      logical lex
!
!
      if( lmype.eq.lmspe ) then
      nft = 21
      if(use_namelist) then
!..   use variable in inppls
         call nopen(nft,DEGAS_REFL,"text",dsn,lex)
      else
!..   use environment variable
         call nopen(nft,"DEGAS_REFL","text",dsn,lex)
      end if
      write(n6,'( 2x,"nft  =",i3,"  dsn =",a)') nft,dsn
      if( .not.lex ) then
      call wexit("rflini","no refl-data")
      endif
      endif
!
      write(n6,'(/2x,"*** rflini ***")')
      write(n6,'( 2x,"ctyp =",a)') ctyp
!
!----------------------------------------------------------------------
!::energy
!----------------------------------------------------------------------
      mpe = npe
      zdel=0.1*log(10.0)
      do 10 je=1,mpe
      rfelg(je)=real(je-1)*zdel
      rfeng(je)=exp(rfelg(je))
   10 continue
      emin=rfeng(1)
      emax=rfeng(mpe)
      dlge=real(mpe-1)/rfelg(mpe)
!
      write(n6,'( 2x,"mpe  =",i5,2x,"emin =",1pe11.3,"  emax ="
     >  ,1pe11.3)')  mpe,emin,emax
!
!----------------------------------------------------------------------
!::read reflection data
!----------------------------------------------------------------------
      if( lmype.eq.lmspe ) then
!
      kon=0
  100 continue
      read(nft,end=190,fmt='(a)') clin
      if(index(clin,ctyp).gt.0) then
      if(index(clin,'RN').gt.0) then
      kon = kon + 1
      js=0
      do 110 i=1,8
      js=js+1
      read(nft,64) (rfrn(j),j=js,js+5)
      js=js+5
  110 continue
      endif
      if(index(clin,'ER').gt.0) then
      kon = kon + 1
      js=0
      do 120 i=1,8
      js=js+1
      read(nft,64) (rfer(j),j=js,js+5)
      js=js+5
  120 continue
      endif
      endif
      go to 100
  190 continue
!
      if( kon.ne.2 ) then
      call wexit("refle","no found data of "//ctyp)
      endif
!
      write(n6,'(2x,"read refle data of ",a)') ctyp
!
!
   64 format(6(1x,1pe12.5))
!
      close (nft)
      endif

!::debug write
      write(n6,'(2x,"eng  =",1p10e11.3:/(2x,6x,1p10e11.3))')
     >  (rfeng(i),i=1,mpe)
      write(n6,'(2x,"rn   =",1p10e11.3:/(2x,6x,1p10e11.3))')
     >  (rfrn(i),i=1,mpe)
      write(n6,'(2x,"er   =",1p10e11.3:/(2x,6x,1p10e11.3))')
     >  (rfer(i),i=1,mpe)
!
      return
      end
!
!***********************************************************************
      subroutine refcof(ein,erfl,rn,re)
!***********************************************************************
      use cntrfl, only : dlge, mpe, rfelg, rfer, rfrn
      implicit none
      real(8), intent(in)  :: ein
      real(8), intent(out) :: erfl, rn, re
!
!::local variables
      real*8 te,zle,fce
      integer nte
!
!---ein < 1.0
      if( ein.le.1.0d0 ) then
      erfl = ein
      rn = 0.0d0
      re = 0.0d0
      return
      endif
!
!--intepolate
      te=ein
      zle=log(te)
      nte=int(zle*dlge)+1
      nte=max(1,nte)
      nte=min(nte,mpe-1)
      fce=(zle-rfelg(nte))/(rfelg(nte+1)-rfelg(nte))
!
      rn  =rfrn(nte)+fce*(rfrn(nte+1)-rfrn(nte))
      erfl=rfer(nte)+fce*(rfer(nte+1)-rfer(nte))
!-----
      if( erfl.gt.ein ) erfl = ein
!-----
      re  =erfl/ein*rn
      return
      end
