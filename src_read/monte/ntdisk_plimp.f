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
      use cunit,       only : lmspe, lmype, n6
      use cntwedf,     only : nene_mx, nsmx, nwmx
      use csize,       only : ndgs, ndmc, ndwp, ndmsr
      implicit none
!
      integer,   intent(in) :: nft, kk
      character, intent(in) :: cft*(*)
!
!::local variables
      logical lex, is_err
      character csrc*6, dsn2*80
      integer  isrc, i, nwkmp_nt_old, j
      real*8, allocatable :: vmwork_tmp(:,:)
!
      if( nft.le.0 ) return
!-----------------------------------------------------------------------
!::write
!-----------------------------------------------------------------------
      if( kk.eq.1 ) then
      if( lmype.ne.lmspe ) return
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
!
!::master pe
      if( lmype.eq.lmspe ) then
        call nopen(nft,cft,"binary read",dsn2,lex)
        is_err = .true.
        read(nft,end=10,err=10)
     >     vtime,vflux,vmwork
     >    ,vitim,vitnt,vnsrc,vsmty,vsmno
     >    ,visrc,vkflx,vksty,vnsmp,vcsrc
     >    ,cntmnt_emrk
        close (nft)
        is_err = .false.
!
! for old version DTNTL(before Wed Jul 16 2025)
 10     continue
        if(is_err) then
          ! size of vmwork in old version
          nwkmp_nt_old = 1 + (ndmc+1)*(ndgs*2+2+ndgs) 
     >   + (ndmc+1)*(ndgs*3) + (ndmc+1)*(ndgs*3)
     >   + (ndmc+1)*(ndgs+2*2) + (ndmc+1)*(2*3) 
     >   + ndwp*2+10+1+30*3*2  + ndwp*4 
     >   + 11+6 + nene_mx*nwmx*nsmx  + ndwp*3
!
          allocate(vmwork_tmp(nwkmp_nt_old,ndmsr))
          ! read binray
          rewind(nft)
          read(nft,end=910,err=910)
     >      vtime,vflux,vmwork_tmp
     >     ,vitim,vitnt,vnsrc,vsmty,vsmno
     >     ,visrc,vkflx,vksty,vnsmp,vcsrc
     >     ,cntmnt_emrk
          close(nft)
          ! set value for vmwork
          do i = 1, nwkmp_nt_old
            do j = 1, ndmsr
              vmwork(i,j) = vmwork_tmp(i,j)
            enddo
          enddo
          deallocate(vmwork_tmp)
        endif
      endif
!
      endif
!
!----------------------------------------------------------------------
!::debug write
!----------------------------------------------------------------------
      do isrc = 1, vnsrc
      csrc = vcsrc(isrc)
      dotn = vflux(isrc)
!
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
