!***********************************************************************
      subroutine ntdisk(nft,kk,cft)
!***********************************************************************
!
!::No confusion wwork
!     csize   nwksz = ....        ==> nwkmp_nt
!     cntwcn  wwork(nwksz)        ==> comment  real*8 wmwork(nwkmp_nt)
!     cntmnt  vwork(nwksz,ndmsr)  ==> vmwork(nwkmp_nt,ndmsr)
!     cntxwk  nxksz=nwksz         ==> nxksz = nwkmp_nt
!
!-----------------------------------------------------------------------
      use cntmnt,      only : dotn, vcsrc, vflux, visrc, vitim, vitnt
     >    , vkflx, vksty, vmwork, vnsmp, vnsrc, vsmno, vsmty, vtime
     >    , cntmnt_emrk
      use cntwcn,      only : wabs, werr, wion, wsum, wtot
      use cunit,       only : lmspe, lmype, lnope, mywld, n6
      use mod_dtypdef, only : nddt, typ_dnam, typ_itag
      use mpi!,         only : mpi_bcast
      implicit none
!
      logical lex
! modified 2/2 lines organize local variables and include files by kamata 2021/06/16
!ik   integer  nft, kk
!ik   character cft*(*)
      integer,   intent(in) :: nft, kk
      character, intent(in) :: cft*(*)
!
!::local variables
! deleted 3 lines replace all include files with module files by kamata 2021/08/18
!ik   real*8 wmwork(nwkmp_nt)
!ik   equivalence (wflx,wmwork(1))
!ik   integer  nwmsz
!
! modified 3/1 lines organize local variables and include files by kamata 2021/06/16
!ik   character ctim*5, cnum*15, dsn*80, cact*5, csrc*6, dsn2*80
!ik   integer mji, indxr, itmn;  real*8 tms
!ik   real*8  timnt, ritnt
      character csrc*6, dsn2*80
! modified 1/1 lines replace all include files with module files by kamata 2021/08/18
!ik   integer  nls, nle, nbf, ierr
      integer  ierr, itag, ityp
      integer  isrc, i
!
      if( nft.le.0 ) return
! deleted 2 lines organize local variables and include files by kamata 2021/06/16
!ik   timnt = time
!ik   ritnt = itim
! deleted 1 line replace all include files with module files by kamata 2021/08/18
!ik   nwmsz = nwkmp_nt
!
!-----------------------------------------------------------------------
!::write
!-----------------------------------------------------------------------
      if( kk.eq.1 ) then
      if( lmype.ne.lmspe ) return
! deleted 1 line organize local variables and include files by kamata 2021/06/16
!ik   cact = "write"
!
!::ctim
! deleted 4 lines organize local variables and include files by kamata 2021/06/16
!ik   tms = time*1.0e3
!ik   write(cnum,'(f15.3)') tms
!ik   mji  = indxr(cnum," ")
!ik   ctim = cnum(mji+1:mji+5)
!
!xx      mji = index(cft," ")
!xx      dsn = cft(1:mji-1)//"_"//ctim
! deleted 1 line organize local variables and include files by kamata 2021/06/16
!ik   dsn = cft
!xx   open(unit=nft,file=dsn,form="unformatted",action="write")
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   call nopen(nft,dsn,"binary write",dsn2,lex)
      call nopen(nft,cft,"binary write",dsn2,lex)
!
      write(nft)
     >     vtime,vflux,vmwork
     >    ,vitim,vitnt,vnsrc,vsmty,vsmno
     >    ,visrc,vkflx,vksty,vnsmp,vcsrc
     >    ,cntmnt_emrk
      close (nft)
      endif
!
!-----------------------------------------------------------------------
!::read
!-----------------------------------------------------------------------
      if( kk.eq.2 ) then
! deleted 2 lines organize local variables and include files by kamata 2021/06/16
!ik   cact = "read"
!ik   dsn = cft
!
!::master pe
      if( lmype.eq.lmspe ) then
!xx   open(unit=nft,file=dsn,form="unformatted",action="read")
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   call nopen(nft,dsn,"binary read",dsn2,lex)
      call nopen(nft,cft,"binary read",dsn2,lex)
      read(nft,end=910,err=910)
     >     vtime,vflux,vmwork
     >    ,vitim,vitnt,vnsrc,vsmty,vsmno
     >    ,visrc,vkflx,vksty,vnsmp,vcsrc
     >    ,cntmnt_emrk
      close (nft)
      endif
!
!::send/receive
!::[MPI_Bcast in ntdisk]  cntmnt /cntmnt/ (vtime,emrk)  10/04/21
      if( lnope.gt.1 ) then
! modified 4/3 lines replace all include files with module files by kamata 2021/08/18
!ik   nls = loc( vtime )
!ik   nle = loc( cntmnt_emrk ) + 1
!ik   nbf = nle - nls
!ik   call MPI_Bcast( vtime, nbf, MPI_BYTE, lmspe, mywld, ierr)
      call tbfind( nddt, typ_dnam, ityp, 'CNTMNT' )
      itag = typ_itag(ityp)
      call MPI_Bcast( vtime, 1, itag, lmspe, mywld, ierr )
! added 1 line dynamic allocation of arrays by kamata 2022/05/29
      call cntmnt_sr( 1, lmspe )
      endif
!
      endif
!
!----------------------------------------------------------------------
!::debug write
!----------------------------------------------------------------------
! deleted 1 line organize local variables and include files by kamata 2021/06/16
!ik   itmn = ritnt
      do isrc = 1, vnsrc
      csrc = vcsrc(isrc)
      dotn = vflux(isrc)
!
! modified 3/1 lines replace all include files with module files by kamata 2021/08/18
!ik   do i = 1, nwmsz
!ik   wmwork(i) = vmwork(i,isrc)
!ik   enddo
      call setvmwork( 1, isrc )
      call ntwcon
!
      write(n6,'(2x,"=== ntdisk ===",i2,2x,a,"  dotn =",1pe14.6,
     > "  wtot =",1p2e14.6,"  werr =",1pe14.6,"  wion/wabs =",
     > 1p2e14.6)') isrc,csrc,dotn,wtot,wsum,werr,wion,wabs
      enddo
!
      return
!
  910 continue
      call wexit("ntdisk","error occured during read file")
      end
