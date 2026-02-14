!***********************************************************************
      subroutine ncopy(nft,ckey,cfrm,cdsn,lex)
!***********************************************************************
      use csize,     only : nzsmx
      use cunit,     only : cmype, lmspe, lmype, lnope, mjpe, mygrp
     >    , mype, mywld, n6
      use mod_shexe, only : use_namelist
      use mpi!,       only : mpi_bcast, mpi_character, mpi_integer
      implicit none

!::argument
      integer, intent(in)    :: nft
      logical, intent(in)    :: lex ! dummy
      character, intent(in)  :: cfrm*(*), ckey*(*)
      character, intent(out) :: cdsn*80

!::local variables
      integer :: istat, ndmx = 0
      character  clin*120, cmsg*80, fmw*80
      character, allocatable :: ctbl(:)*120
      integer i, ierr, ii, mj, ncha, nmax
! function
      integer  lenx

      write(n6,'(2x,"=== passed ncopy  ",a,2x,"mype =",i5)') ckey,mype

!::check
      if( cfrm(1:1).ne."t" ) then
      call wexit("ncopy","binary file can not be copied.")
      endif

      if(use_namelist) then
         cdsn = ckey
      else
         call getenv( ckey, cdsn )
      end if
      mj = lenx(cdsn)
      if( mj.le.0 ) goto 920

! get number of records in csdn
      if( ndmx == 0 .and. lmype == lmspe ) then
        write(fmw,'(i10)') mygrp
        fmw = '@rec_' // adjustl( fmw )
        call getnrec( cdsn(1:mj), fmw, ndmx )
        write(n6,*) 'ncopy ndmx = ', ndmx
      endif
! allocate
      if( lnope > 1 )
     >  call MPI_Bcast( ndmx, 1, MPI_INTEGER, lmspe, mywld, ierr )

      allocate( ctbl(ndmx), stat = istat )
      if( istat /= 0 ) then
        write(cmsg,'(a,3i8)')
     >    'ctbl allocate error in ncopy, istat,mygrp,lmype = '
     >  , istat, mygrp, lmype
        call wexit( 'ncopy', trim( cmsg ) )
      endif

!::copy file
      write(n6,602) "=== copy file (text)  ",
     >  mype,nft,ckey(1:lenx(ckey)),cdsn(1:mj)

!::master pe
      if( lmype.eq.lmspe ) then
      open(unit=nft,file=cdsn(1:mj),form="formatted")
      ii = 0
 110  continue
      read(nft,'(a)',end=120) clin
      ii = ii + 1
      if( ii.gt.ndmx ) goto 910
      ctbl(ii) = clin
      goto 110
 120  continue
      close (nft)
      nmax = ii
!-----
      write(n6,'(2x,"nmax = ",i5,"  nft =",i5)') nmax, nft
      endif

!::send data
      if( lnope.gt.1 ) then
      ncha = 120*ndmx
      call MPI_Bcast( ctbl, ncha, MPI_CHARACTER, lmspe, mywld, ierr)
      endif

!::write
      cdsn = cdsn(1:mj)//"_tmp_"//cmype(mjpe:8)
      mj   = lenx(cdsn)
      write(n6,602) "=== write temp file   ",mype,nft,ckey,cdsn(1:mj)
      open(unit=nft,file=cdsn)
      do i = 1, ndmx
      write(nft,'(a)') ctbl(i)(1:lenx(ctbl(i)))
      write(n6,'(2x,a)') ctbl(i)(1:lenx(ctbl(i)))
      enddo
      close (nft)
      return

!::error
 910  continue
      call wexit("ncopy","too many lines")
 920  continue
      call wexit("ncopy","no found setenv "//ckey)

 602  format(2x,a,"  mype =",i5,"  nft =",i3,"  ckey =",a,"  R-cdsn ="
     >  ,a)

      end
