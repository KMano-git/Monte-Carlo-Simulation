!      module mpi
!
!      integer    MPI_INTEGER, MPI_BYTE
!      integer    MPI_CHARACTER
!      real(8)    MPI_REAL8, MPI_SUM!, MPI_WTIME
!
!      integer    MPI_STATUS_SIZE
!      parameter (MPI_STATUS_SIZE=1)
!
!      contains
c***********************************************************************
      subroutine MPI_Bcast!(*, no, dtype, root, comm, ierr)
c***********************************************************************
      implicit none
!      real*8 
      character buf
      integer no, dtype, root, comm, ierr
c
      return
      end
!
      subroutine MPI_RECV
      return
      end
!
      subroutine MPI_REDUCE
      return
      end
!
      subroutine MPI_SEND
      return
      end
!
c**********************************************************************
      subroutine mpiend
c**********************************************************************
      return
      end
!
      subroutine mpifile
      return
      end
!
c**********************************************************************
      subroutine mpista
c**********************************************************************
      use cunit
c
      mype  = 0
      mspe  = 0
      nope  = 1
      lmype = 0
      lmspe = 0
      lnope = 1
!
      cdgr(1) = 1
      cdgr(2) = 1
      cdgr(3) = 1
c
      ipg = 0
      write(6,'(2x,"*** mpista ***   mype =",i5,"  mspe =",i3,
     >  "  nope =",i3)') mype,mspe,nope
c
      return
      end
!
c***********************************************************************
      subroutine ncopy(nft,ckey,cfrm,cdsn,lex)
c***********************************************************************
      use mpi
      use cunit, only : cmype, lmspe, lmype, lnope, mjpe, mygrp
     >    , mype, n6
      implicit none
c
      integer nft,lenx,l,mj
      logical lex
      character ckey*(*),cdsn*80,cfrm*(*)
      character cenv*80,cact*5
c
c::local variables
      integer :: istat, ndmx = 0
      character  clin*120, cmsg*80, fmw*80
      character, allocatable :: ctbl(:)*120

      integer ii,nmax,ncha,i
c
      write(n6,'(2x,"=== passed ncopy  ",a,2x,"mype =",i5)') ckey,mype
c
c::check
      if( cfrm(1:1).ne."t" ) then
      call wexit("ncopy","binary file can not be copied.")
      endif
c
      call getenv( ckey, cdsn )
      mj = lenx(cdsn)
      if( mj.le.0 ) goto 920

! added lines dynamic allocation of arrays
! get number of records in csdn
      if( ndmx == 0 .and. lmype == lmspe ) then
        write(fmw,'(i10)') mygrp
        fmw = '@rec_' // adjustl( fmw )
        call getnrec( cdsn(1:mj), fmw, ndmx )
        write(n6,*) 'ncopy ndmx = ', ndmx
      endif

! allocate
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
c
c::master pe
      if( lmype.eq.lmspe ) then
      open(unit=nft,file=cdsn(1:mj),form="formatted")
      ii = 0
 110    continue
      read(nft,'(a)',end=120) clin
      ii = ii + 1
      if( ii.gt.ndmx ) goto 910
      ctbl(ii) = clin
      goto 110
 120  continue
      write(6,'(2x,"##1111 mpidmy master_pe nft = ",i6)') nft
      close (nft)
      nmax = ii
!-----
      write(n6,'(2x,"nmax = ",i5,"  nft =",i5)') nmax, nft
      endif
c
c::send data
      if( lnope.gt.1 ) then
      ncha = 120*ndmx
      endif
c
c::write
      cdsn = cdsn(1:mj)//"_tmp_"//cmype(mjpe:8)
      mj   = lenx(cdsn)
      write(n6,602) "=== write temp file   ",mype,nft,ckey,cdsn(1:mj)
      open(unit=nft,file=cdsn)
      do i = 1, ndmx
      write(nft,'(a)') ctbl(i)(1:lenx(ctbl(i)))
      write(n6,'(2x,a)') ctbl(i)(1:lenx(ctbl(i)))
      enddo
      write(6,'(2x,"##1111 mpidmy send_data nft = ",i6)') nft
      close (nft)
      return
c
c::error
 910    continue
      call wexit("ncopy","too many lines")
 920    continue
      call wexit("ncopy","no found setenv "//ckey)
c
 602    format(2x,a,"  mype =",i5,"  nft =",i3,"  ckey =",a,"  R-cdsn ="
     >  ,a)
c
      end
!
c***********************************************************************
      subroutine wexit(csub,cmsg)
c***********************************************************************
      implicit none
c
      character csub*(*),cmsg*(*)
c
      character dsn*6, wpg*6
      integer n7, mype, lenx
c
      n7 = 7
      dsn = "erstop"
      wpg = "PROG"
      mype = 0
c
      open( unit=n7, file=dsn )
      write(n7,'(/2x,"*** wexit ***    PRG =",a,"  mype =",i4)')
     >   wpg,mype
      write(n7,'(5x,"stop at sub. ",a)') csub(1:lenx(csub))
      write(n7,'(5x,"because of   ",a)') cmsg(1:lenx(cmsg))
c
      write(6,'(/2x,"*** wexit ***    PRG =",a,"  mype =",i4)')
     >   wpg,mype
      write(6,'(5x,"stop at sub. ",a)') csub(1:lenx(csub))
      write(6,'(5x,"because of   ",a)') cmsg(1:lenx(cmsg))
c
      stop  ! Abnormal end (wexit)
      end
!
!
      subroutine mpi_file_open
      return
      end
!
      subroutine mpi_file_seek
      return
      end
!
      subroutine mpi_file_write
      return
      end
!
      subroutine mpi_file_sync
      return
      end
!
      subroutine mpi_file_close
      return
      end
!
      subroutine mpi_barrier
      return
      end

!      end module mpi
