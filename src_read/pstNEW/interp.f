!***********************************************************************
      subroutine test_interp
!***********************************************************************
      implicit none
!
      character cond*80
      integer    ndim
      parameter (ndim=20)
      integer  nval, kval(ndim)
!
      cond = "2 - 5,8, 10-15, 21"
      call interp(cond,nval,kval,ndim)
!
      stop
      end
!
!
!***********************************************************************
      subroutine interp(cond,nval,kval,ndim)
!***********************************************************************
!
!     "5-8"
!     "5, 15, 16, 17, 18"
!     "0"
!
!-----------------------------------------------------------------------
      implicit none
!
!::argument
      character  cond*(*)
      integer    ndim
      integer    nval, kval(ndim)
!
!::local variables
      character  clin*80, cwrd*20
      integer    i, ii, j, n, jsta, jend, jinc, jnum
      integer    nk, mjs(30), mje(30)
      integer    nk2, mjs2(30), mje2(30)
!
      clin = cond
      call linsep(clin,",",nk,mjs,mje,30)
!
      ii = 0
      do i = 1, nk
      cwrd = clin(mjs(i):mje(i))
      if( index(cwrd,"-").gt.0 ) then
        call linsep(cwrd,"-",nk2,mjs2,mje2,10)
        if( nk2.ne.2 ) goto 910
!
      read(cwrd(mjs2(1):mje2(1)),*) jsta
      read(cwrd(mjs2(2):mje2(2)),*) jend
!
      jinc = 1
      if( jend.lt.jsta ) jinc = -1
!
      do n = jsta, jend, jinc
        ii = ii + 1
        kval(ii) = n
      enddo
!
      else      
        read(cwrd,*) jnum
        ii = ii + 1
        kval(ii) = jnum
      endif
      enddo
!
      nval = ii
!
!::debug write
!x      write(6,'(2x,"cond = [",a,"]","  nval =",i3)')
!x     >    trim(cond), nval
!x      write(6,'(2x,"kval =",15i5)') (kval(i),i=1,nval)
!      
      return
!
 910  continue
      write(6,'(2x,"error at sub. interp  [",a,"]")') trim(cond)
      write(6,'(2x,"word can not be separated [",a,"]")') cwrd
      stop
      end
