!***********************************************************************
      subroutine atm_set_x(drnm,catm,kis)
!***********************************************************************
!
!         drnm : full path of dir which contains adas-data
!
!-----------------------------------------------------------------------
      use catcom, only : ched, dratm, nfatm, rtyp, tbdsn
      use cunit,  only : n6
      implicit none
!
!::argument
! modified 1/2 lines organize local variables and include files by kamata 2021/06/28
!ik   character(*) :: drnm, catm
      character(*), intent(in)  :: drnm
      character(*), intent(out) :: catm
      integer, intent(in) :: kis
!
!::local variables
      integer, parameter :: ndm = 20
! modified 2/1 lines organize local variables and include files by kamata 2021/06/28
!ik   character(80), dimension(ndm) :: clst
!ik   integer :: nls, i, id, m1, m2
      integer :: i, m1, m2
!
      m1 = index( drnm, "/", .true. )
      m2 = index( drnm, "_", .true. )
      catm = drnm(m1+1:m2-1)
!
      write(n6,'(/2x,"*** atm_set ***")')
      write(n6,'(2x,"drnm = ",a,"  catm = ",a)')
     >    trim(drnm), trim(catm)
!
      if( len_trim(drnm).le.0 ) then
      write(n6,'(2x,"stop at sub. atm_set due to NULL file")')
! added 1 line correction of stop statement only by kamata 2021/12/28
      call wexit( 'atm_set', 'due to NULL file' )
      stop
      endif
!
!::tbdsn
! modified 1/1 lines multiple IMPMC execution bug by kamata 2022/04/15
!ik   call getls(drnm,"f",nfatm,tbdsn,10)
      call getls(drnm,"f",nfatm,tbdsn,10,' ')
!
!::table
      rtyp = (/"ION","REC","CXR", "PLT", "RAD", "PRB", "PRC"/)
      ched = (/"scd","acd","ccd", "plt", "plt", "prb", "prc"/)
!
      write(n6,'(2x,"rtyp = ",10(a,2x))') (rtyp(i),i=1,7)
      write(n6,'(2x,"ched = ",10(a,2x))') (ched(i),i=1,7)
!
      write(n6,'(2x)')
      write(n6,'(2x,"adas-file = ",10(a,2x))')
     >   (trim(tbdsn(i)),i=1,nfatm)
!
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   dratm = drnm
      dratm(kis) = drnm
!
!::how to get dsn name
!      call atm_fnam("ATOM_ION",dsn)
!
      call flush(n6)
      return
      end
!
!***********************************************************************
! modified 1/1 lines multiple IMPMC execution bug by kamata 2022/04/15
!ik   subroutine getls(cdir,cknd,nls,clst,ndm)
!       subroutine getls(cdir,cknd,nls,clst,ndm,cid)
! !***********************************************************************
! !miyamoto>>>
!       use, intrinsic :: iso_c_binding
! !<<<
!       use cunit, only : n6
!       implicit none

! ! modified 3/4 lines organize local variables and include files by kamata 2021/06/28
! !ik   integer :: ndm, nls
! !ik   character(*) :: cdir, cknd
! !ik   character(*), dimension(ndm) :: clst
!       integer,      intent(in)  :: ndm
!       integer,      intent(out) :: nls
! ! modified 1/1 lines multiple IMPMC execution bug by kamata 2022/04/15
! !ik   character(*), intent(in)  :: cdir, cknd
!       character(*), intent(in)  :: cdir, cid, cknd
!       character(*), dimension(ndm), intent(out) :: clst
! !
! !::local variables
! ! modified 3/2 lines organize local variables and include files by kamata 2021/06/28
! !ik   character(80) :: cmsg, clin, cdsn, ccd
! !ik   character(10) :: iam_char
! !ik   integer :: ii, nft, i, iam, ierr
!       character(80) :: cmsg, clin, cdsn
!       integer :: ii, nft, i, ierr
! !miyamoto>>>
!       interface
!          integer(c_int) function system(command) bind(c)
!            import
!            character(kind=c_char),intent(in) :: command(*)
!          end function system
!       end interface
! !<<<
! !
! !::dir  ~/sonicV2/adas2/
! !      call MPI_Comm_rank(MPI_COMM_WORLD, iam,ierr)
! !      write(*,'("GETLS: called by rank ",i0," of MPMD program.")') iam
! !      write(nprg_char,*) nprg
! !      write(iam_char,'(i0)') iam
! !      cdsn = "@ls_"//trim(nprg_char)
! !      cdsn = "@ls_"//trim(iam_char)
! ! modified 1/1 lines multiple IMPMC execution bug by kamata 2022/04/15
! !ik   cdsn = "@ls"
!       cdsn = "@ls" // cid
! !      call system("rm -f @ls")
! !miyamoto     call system("rm -f "//trim(cdsn))
! !miyamoto>>>
!       ierr=system("rm -f "//trim(cdsn)//c_null_char)
!       if(ierr/=0) then
!          write(*,'(a,i0)')
!      &        "'rm -f "//trim(cdsn)//"' returned an error: ",ierr
!       end if
! !<<<
! ! deleted 1 line organize local variables and include files by kamata 2021/06/28
! !ik   ccd  = "cd "//trim(cdir)
!       if( cknd(1:1).eq."d" ) then
!       cmsg = "ls "//trim(cdir)//" -1F | grep \/$  > "//trim(cdsn)
!       elseif( cknd(1:1).eq."f" ) then
!       cmsg = "ls "//trim(cdir)//" -1F | grep -v \/$ > "//trim(cdsn)
!       else
!       cmsg = "ls "//trim(cdir)//" -1F > "//trim(cdsn)
!       endif
! !miyamoto     call system(trim(cmsg))
! !miyamoto>>>
!       ierr=system(trim(cmsg)//c_null_char)
!       if(ierr/=0) then
!          write(*,'(a,i0)')
!      &        "'"//trim(cmsg)//"' returned an error: ",ierr
!       end if
! !<<<
! !
!       ii = 0
!       nft = 21
! !
!       open(nft,file=trim(cdsn),form='formatted')
!  100  continue
!       read(nft,'(a)',end=190) clin
!       if( len_trim(clin).eq.0 ) goto 100
!       ii = ii + 1
!       if( ii.gt.ndm ) goto 910
!       clst(ii) = trim(clin)
!       goto 100
! !
!  190  continue
!       close(nft)
! !
! !::wrong name cdir
!       if( ii.le.0 ) then
!       write(n6,'(/2x,"Stop at sub. getls due to  No found dir ",a)')
!      >    trim(cdir)
! ! added 1 line correction of stop statement only by kamata 2021/12/28
!       call wexit( 'getls', 'due to No found dir' )
!       stop
!       endif
! !
!       nls = ii
!       write(n6,'(/2x,"*** getls ***  cdir = ",a)') trim(cdir)
!       if( cknd(1:1).eq."d" ) then
!       write(n6,'(2x,"dir  = ",10(a,2x))') (trim(clst(i)),i=1,nls)
!       elseif( cknd(1:1).eq."f" ) then
!       write(n6,'(2x,"file  = ",10(a,2x))') (trim(clst(i)),i=1,nls)
!       else
!       write(n6,'(2x,"dr/fl = ",10(a,2x))') (trim(clst(i)),i=1,nls)
!       endif
! !xx   write(n6,'(2x)')
!       return
! !
!       if( nls.le.0 ) then
! !xx     call wexit("atm_fnam","no found data")
!         write(n6,'(2x,"error at sub. getls  cdir = ",a)')
!      >    trim(cdir)
! ! added 1 line correction of stop statement only by kamata 2021/12/28
!         call wexit( 'getls', 'due to No found data' )
!         stop
!       endif
! !
!  910  continue
!       write(n6,'(2x,"error at sub.getls  ii > ndm ",2i5)')
!      >   ii, ndm
! ! added 1 line correction of stop statement only by kamata 2021/12/28
!       call wexit( 'getls', 'due to array size obver' )
!       stop
!       end
!
!***********************************************************************
      subroutine atm_fnam_x(cxsp,dsn,kis)
!***********************************************************************
      use catcom, only : ched, dratm, nfatm, rtyp, tbdsn
      use cunit,  only : n6
      implicit none
!
! modified 1/2 lines organize local variables and include files by kamata 2021/06/28
!ik   character(*) :: cxsp, dsn
      character(*), intent(in)  :: cxsp
      character(*), intent(out) :: dsn
      integer, intent(in) :: kis
!
!
!::local variables
      integer       :: mj, iknd, i
      character(3)  :: cknd
      character(20) :: cmem
!
!::find file
      mj = index(cxsp,"_",.true.)
      cknd = cxsp(mj+1:)
!
      iknd = 0
      do i = 1, 7
      iknd = i
      if( cknd.eq.rtyp(i) ) goto 210
      enddo
      goto 910
 210  continue
!
      do i = 1, nfatm
      cmem = tbdsn(i)
      if( cmem(1:3).eq.ched(iknd) ) goto 220
      enddo
      goto 910
 220  continue
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   dsn = trim(dratm)//"/"//trim(cmem)
      dsn = trim(dratm(kis))//"/"//trim(cmem)
!
      write(n6,'(/2x,"cxsp = ",a,"  dsn = ",a)') trim(cxsp), trim(dsn)
      return
!
 910  continue
      dsn = ""
      write(n6,'(/2x,"cxsp = ",a,"  dsn = [",a,"]  No found file")')
     >   trim(cxsp), trim(dsn)
      return
      end
