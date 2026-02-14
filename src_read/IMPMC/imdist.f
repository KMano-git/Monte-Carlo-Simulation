!***********************************************************************
      subroutine imdist(lptm,dtunt)
!***********************************************************************
!
!       lptm : loop number lt in imp_cal
!
!  comt lpcm  ip  ptmx       ptim,      tt         pdtm        ic  ix  iy
!   ln  en  ko  is ml   rr       zz      ro     wght   vz        Evel
!     vflw      Ti        lstp      tauz
!
!      lmype = 0            = 1                 = 2
!         50 100 150 200     50 100 150 200      50 100 150 200
!   |-|  ---|---|---|---|  |---|---|---|---|   |---|---|---|---|
!        A       A    isiz   <----ksiz------>    1   2   3   4
!        |       |            nctp = nchd + kcsz*lmype
!       nctp    offset        offset = nctp + isiz*(k-1)
!
!                              icsz = 256*50 line  200 PE
!                              kcsz = icsz*kmax
!
!-----------------------------------------------------------------------
      use cimcom, only : cstrg, is, lprf, ndclng, npmax, sptyc, tt
      use cimctl, only : cdirz
      use cimntl, only : stb
      use cunit,  only : lmspe, lmype, mype, mywld, n6
      use mpi!,    only : mpi_barrier, mpi_character,  mpi_file_close
!    >    , mpi_file_open, mpi_file_seek, mpi_file_sync, mpi_file_write
!    >    , mpi_mode_create, mpi_info_null, mpi_mode_rdwr, mpi_seek_set
!    >    , mpi_status_size
      implicit none
!
!::argument
! modified 2/2 lines organize local variables and include files by kamata 2021/06/28
!ik   integer lptm
!ik   real*8  dtunt
      integer, intent(in) :: lptm
      real(8), intent(in) :: dtunt
!
!::local variables
! modified 3/3 lines organize local variables and include files by kamata 2021/06/28
!ik   integer  ip, nft, l, lpe, nsiz, ninf, i, ii, fh
!ik   integer  ic, ko, lstp, imox, imoy, k, lenx, i6_save, md
!ik   real*8   tauz, dtstp, ptim, pdtm
      integer  ip, nft, i, fh
      integer  k, md
      real*8   ptim, pdtm
      character  cdsn*80, comt*6
      integer    kmax, nchd, icsz, kcsz, nctp
      integer*8  offset   ! <===
! modified 1/1 lines organize local variables and include files by kamata 2021/06/28
!ik   integer  nls, nle, nbf, spe, rpe, ierr
      integer  ierr
      integer  status(MPI_STATUS_SIZE)
!
!::text file
      integer     ndunp
      parameter  (ndunp=5)
      integer     nunp
      character  cln0*(ndclng), clin(ndunp)*(ndclng)
      character  cLF*1
!
!::select time
      if( lprf.eq.0 ) return
!
!::file name
      nft = 50000 + 1
      write(cdsn,'(a,i6.6,a)') trim(cdirZ)//"/"//trim(sptyc)//"_dist_",
     >   lptm,".txt"
      write(n6,'(2x,"imdist check cdsn =",a)') trim(cdsn)
      cLF = char(10)
!
!::header
      if( lmype.eq.lmspe ) then
      open(unit=nft,file=trim(cdsn))
      write(cln0,'(3x,"comt",1x,"lpcm",2x,"ip",2x,"ptmx",7x,"ptim",7x,
     >  "tt",9x,"pdtm",9x,"ic",2x,"ix",2x,"iy",4x,"ln",1x,"ien",2x,
     >  "ko",2x,"is",1x,"ml",3x,"rr",7x,"zz",6x,"rO",5x,"wght",3x,
     >  "vz",8x,"Evel",6x,"vflw",6x,"Ti",8x,"lstp",6x,"tauz",6x,
     >  "dt",8x,"Ftot",6x,"dVz",7x,"V")')
      cln0(ndclng:ndclng) = cLF
      write(nft,'(a)') cln0
      close(nft)
      endif
!
!::loop
      nunp = ndunp
      if( nunp.gt.npmax ) nunp = npmax
      kmax = npmax/nunp
      if( npmax.ne.nunp*kmax ) then
        write(n6,'(2x,"imdist   npmax =",i6,"  nunp/ndunp =",2i5,
     >  "  kmax =",i3)') npmax, nunp, ndunp, kmax
        call wexit("imdist","npmax.ne.nunp*kmax")
      endif
      write(comt,'(i6.6)') mype   !  outimp_mype
!
!::header
      nchd = ndclng         ! header
      icsz = ndclng*nunp    ! 256*50 line
      kcsz = icsz*kmax      ! 256*50*4
      nctp = nchd + kcsz*lmype  ! top
!
! modified 1/1 lines organize local variables and include files by kamata 2021/06/28
!ik   ii = 0
      ip = 0
      do k = 1, kmax
      offset = nctp + icsz*(k-1)
!
      call MPI_FILE_OPEN(mywld, cdsn,
     >     MPI_MODE_CREATE+MPI_MODE_RDWR, MPI_INFO_NULL, fh, ierr)
      call MPI_FILE_SEEK(fh, offset, MPI_SEEK_SET, ierr)
!
      do i = 1, nunp
! modified 2/1 lines organize local variables and include files by kamata 2021/06/28
!ik   ii = ii + 1
!ik   ip = ii
      ip = ip + 1
      cstrg = " "
      if( ip.le.npmax ) then
      ptim = tt(ip)
      pdtm = dtunt
      if( is(ip).eq.0 ) pdtm = stb(ip)
      call impdt(0,ip,comt,ptim,ptim,pdtm)
      endif
      cstrg(ndclng:ndclng) = cLF
      clin(i) = cstrg
      enddo
!
      call MPI_FILE_WRITE(fh, clin, icsz, MPI_CHARACTER, status, ierr)
!
      call MPI_FILE_SYNC(fh,ierr)
      call MPI_FILE_CLOSE(fh,ierr)
      call MPI_BARRIER(mywld,ierr)
      enddo   ! loop(k)
!
      return
!
      md = npmax/5
      if( md.eq.0 ) md = 1
      do ip = 1, npmax
      if( mod(ip,md).ne.0 ) cycle
      ptim = tt(ip)
      if( is(ip).ne.0 ) then
      call impdt(n6,ip,"imdist",ptim,ptim,dtunt)    ! ion
      else
      call impdt(n6,ip,"imdist",ptim,ptim,stb(ip))  ! neutral
      endif
      enddo
!
 190  continue
      return
      end
