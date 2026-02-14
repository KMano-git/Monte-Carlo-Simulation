!x      call test_gpsta
!x      stop
!x      end
!
!***********************************************************************
      subroutine test_gpsta
!***********************************************************************
      use csonic
      implicit none
      integer :: nft
!
      nft = 21
      open(unit=nft,file="inpfig")
      call gpsta(nft)
      close(nft)
!
      return
      end     
!
!***********************************************************************
      subroutine gpsta(nft)
!***********************************************************************
!
! &upfgsz
!   csiz(1) = "Main = 1.6   4.2   -2.4   2.4",
!   csiz(2) = "Div  = 1.8   3.2   -2.4  -1.4",
!   csiz(3) = "UpM  = 1.6   3.8    0.9   2.6",
!   cleg(1) = "5.0, 1.0",
!   clge(2) = "2.8, -1.2",
!   cleg(3) = "  ",
! &end
!             cgnm    cgsz        ngsz = 3
!
! &upcntv
!  cnvl( 1)= "Ni   = 0.5 1.0 1.5 2.0 3.0 4.0 5.0 10.0 20.0 40.0",
!  cnvl( 2)= "Mach = 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9",
!  cnvl( 3)= "Vp   = -4.0 -3.0 -2.0 -1.0 -0.5 0.0 0.5 1.0 2.0 3.0 4.0",
!  cnvl( 5)= "Ti = 1.0 2.0 5.0 10.0 20.0 50.0 100.0 200.0 300.0 400.0",
!  cnvl( 6)= "N0  = 99.0",
!  cnvl(16)= "Wrad = 0.1 0.2 0.5 1.0 2.0 3.0 4.0 5.0 6.0 7.0",
! &end
!
!     ctnm(1) = "Ni"
!     ctcv(1) = "0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 30.0, 50.0"
!     ctnm(2) = "log(Te)"
!     ctcv(2) = "0.0,3.0,(0.2)"
!
!-----------------------------------------------------------------------
      use cunit
      use com_gmsiz
      use com_gpsta
      implicit none

!
      integer  nft
!
!::local variables
      integer   lenx, mji, ii, i, no, j, lst
      integer   nk, mjs(20), mje(20)
      character cinp*80, clin*80, cnam*4
      character(80) :: clgp
      save no
!
!::input data
      character(80), dimension(10) :: csiz, cleg
      character(80), dimension(30) :: cnvl
      namelist /upfgsz/ csiz, cleg
!
      write(6,'(/2x,"*** gpsta ***")')
!
!----------------------------------------------------------------------
!::gsiz
!----------------------------------------------------------------------
      if( nft.gt.0 ) then
      csiz(1:10) = "  "
      cleg(1:10) = "  "
      cnvl(1:30) = "  "
      rewind nft
      read(nft,upfgsz)
!
      ii = 0
      do i = 1, 10
      clin = csiz(i)
      if( len_trim(clin).le.0 ) cycle
      call linsep(clin,":",nk,mjs,mje,20)
      if( nk.ne.2 ) cycle
      ii = ii + 1
      cgnm(ii) = trim(clin(mjs(1):mje(1)))
      cgsz(ii) = trim(clin(mjs(2):mje(2)))
!xx      write(6,'(2x,"ii =",i2,"  cgnm = [",a,"]  cgsz = [",a,"]")')
!xx     >   ii, cgnm(ii), trim(cgsz(ii))
      enddo
      ngsz = ii
!
      write(6,'(2x,"ngsz =",i2)') ngsz
      if( ngsz.le.0 ) then
        write(6,'(2x,"invalid input data upfgsz")')
        write(6,'(2x,"correct example  [",a,"]")')
     >       'csiz(2) = "Div  : 1.8,  3.2,  -2.7, -1.7"'
        stop
      endif

      do i = 1, ngsz
      write(6,'(2x,i2,2x,a,2x,a)') i, cgnm(i), trim(cgsz(i))
      enddo
!
!----------------------------------------------------------------------
!::ctnm, ctcv
!----------------------------------------------------------------------
      ii = 0
      lst = 0
      rewind nft
 110  continue
      read(nft,'(a)',end=190) clin
      call linsep(clin, "=", nk, mjs, mje, 20)
      if( nk.eq.0 ) goto 110
      if( index(clin(mjs(1):mje(1)),"&upcntv").gt.0 ) then
        lst = 1
        goto 110
      endif
      if( lst.eq.0 ) goto 110
      if( index(clin(mjs(1):mje(1)),"&end").gt.0 ) goto 120
      ii = ii + 1
      if( ii.gt.ndvr ) goto 910
      ctnm(ii) = clin(mjs(1):mje(1))
      ctcv(ii) = clin(mjs(2):mje(2))
      goto 110
!
 120  continue
 190  continue
      rewind nft
      nvar = ii
!
      write(6,'(2x)')
      do i = 1, nvar
      write(6,'(2x,i3,2x,"[",a,"]",2x,"[",a,"]")')
     >   i, ctnm(i), trim(ctcv(i))
      enddo
!
      endif      
!------------------------------- nft > 0
!
      do i = 1, ngsz
      write(clin,'(i2,2x,a,2x,a)') i, cgnm(i), trim(cgsz(i))
      call gdput(trim(clin))
      enddo
!
      call gdget("enter plot size (1/2/-1.8,2.4,-2.1,-0.0/q) ==> ",
     >  cinp,mji)
      no = 0
      if( mji.eq.0 ) goto 210
      if( mji.eq.1 ) then
      if( cinp(1:1).eq."q" ) goto 210
      read(cinp,*) no
      if( no.ge.1 .and. no.le.10 ) then
      cinp = cgsz(no)
      else
      no = 1
      cinp = cgsz(1)
      endif
      write(6,'(2x,"gmsta  no =",i2,"  [",a,"]")') no,cinp(1:lenx(cinp))
      endif
      read(cinp,*) rmi,rmx,zmi,zmx
      write(6,'(2x,"gpsta  rmi,rmx,zmi,zmx =",4f9.4)') rmi,rmx,zmi,zmx
!
!::legend postion
      clgp = "  "
      if( no.gt.0 ) clgp = cleg(no)
      if( len_trim(clgp).gt.0 ) then
      read(clgp,*) pslgx, pslgy
      else
      pslgx = rmx - 0.2d0*(rmx-rmi)
      pslgy = zmx - 0.2d0*(zmx-zmi)
      endif
      write(6,'(2x,"gpsta  clgp = [",a,"]","  pslgx,pslgy =",2f9.4)')
     >   clgp(1:len_trim(clgp)+1), pslgx, pslgy
!
      call flush(6)
      mfig = 1
cx      write(6,'(/2x,"*** gpsta ***   lfig =",i2,"  plot region =",
cx     >  4f10.4)') mfig,rmi,rmx,zmi,zmx
c
      return
c
 210  continue
      mfig = 0
!KH      call gdpag(0)
      call gdmod(0)
      stop
      return
!
 910  continue
      write(n6,'(/2x,"Stop at sub. gpsta ii.gt.ndvr ",2i4)') ii, ndvr
      stop
      end
