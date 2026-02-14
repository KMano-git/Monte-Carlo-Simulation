!**********************************************************************
      subroutine figdat( cpath, ctime )
!**********************************************************************
!
!                  IMPMC                soldor (call plsorc)
!                                       /cntsrc/
!
!        a005128/PRJ/sonicV3/src/soldor     2013/06/05
!
!  group in charge of the code   2-GRP type
!    [PLS]: 1 / NTL]: 2 /[IMP]: 2 /  <== cdgr(i) i=1,2,3(ncode)
!
!----------------------------------------------------------------------
      use cntcom, only : ncmax, iplx, iply, mrgn
      use cplcom, only : qtim, vne, vnezef, vni, vte, vti
      use csonic, only : time
      use cunit,  only : lmspe, lmype, mype, n6
      implicit none
! arguments
      character, intent(in) :: cpath*(*), ctime*(*)
! cpath : path name
! ctime : time

!::local variables
      integer :: i60
      integer :: nskp = 500
      integer :: kio, nft = 21, ic, ir
      integer :: icx, icy
      character :: cdsn*200, ctyp*10
!
      if( lmype.ne.lmspe ) return
!
!::soldor
      i60 = 180000 + mype
      kio = 1
!
      write(n6,'(/2x,"*** figdat ***",2x,a)') "DTWRD  w"
!---------------------------------------------------------------------
!::DTPLS : conservative variables in soldor (q1a,q2a,q3,q4)
!::DTPLV : non-conservative variables       (vna,vva,vti,vte)
!---------------------------------------------------------------------
      ctyp = "DTPLV"
      if( ctime == ' ' ) then
        cdsn = "./" // trim(ctyp)
      else
        cdsn = trim( cpath ) // trim( ctyp ) // '_' // trim( ctime )
      endif
      call plvdsk(nft,kio,cdsn)
!
      write(i60,'(/2x,"*** plvdsk *** ",a,"  vne,vna,vti,vte")')
     >   trim(cdsn)
      write(i60,'(2x,"time,qtim =",1p2e14.6)') time, qtim
!
      write(i60,'(6x,"ic",2x,"icx",1x,"icy",1x,"ir",2x,"vne",9x,"vne2",
     >  8x,"vni",9x,"vte",9x,"vti")')
      do ic = nskp, ncmax, nskp
        icx = iplx(ic)
        icy = iply(ic)
        ir  = mrgn(ic)
        if( icx.le.0 .or. icy.le.0 ) cycle
        write(i60,'(2x,i6,i5,i4,i3,1p5e12.4)')
     >   ic, icx, icy, ir, vne(icx,icy), vnezef(icx,icy), vni(icx,icy),
     >   vte(icx,icy), vti(icx,icy)
      enddo
!
!
      write(n6,'(/2x,"*** figdat ***",2x,a)') "DTPLV  w"
!---------------------------------------------------------------------
!::DTWRD : radiation in soldor
!---------------------------------------------------------------------
!     wtime, wfac, wcre, wmce, wime  (wmci, wimi)
!
      ctyp = "DTWRD"
      if( ctime == ' ' ) then
        cdsn = "./" // trim(ctyp)
      else
        cdsn =  trim( cpath ) // trim( ctyp ) // '_' // trim( ctime )
      endif
      open(unit=nft,file=cdsn,form="unformatted")
      call plwrdt(nft,kio,cdsn)   ! close
!
      return
      end
