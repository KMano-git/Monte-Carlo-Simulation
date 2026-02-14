!***********************************************************************
      subroutine gdfps_new(ncl,kclr,ndcl)
!***********************************************************************
!
!     ncl    (i) : number of paint colour
!     kclr   (o) : colour type
!     ndcl   (i) : dimension size
!
!-----------------------------------------------------------------------
      implicit none
!
!::argument
      integer ndcl
      integer ncl, kclr(ndcl)
!
!::local variables
      integer maxclr, i, mji, lenx
      character clin*80, cinp*80, cdef*80
!
      cdef = "1,5,10,14,19,23,28,31"
!
      write(clin,'("0(white),1(blue),31(red)  no =",i3)') ncl+1
      call gdput(clin(1:lenx(clin)))
      call gdget(">> enter clr-no. (0,1,6,12,24,28,31/d) ==> ",
     >  cinp,mji)
!
      if( mji.le.0 ) then
        call gdfps(ncl)
        do i = 1, ncl+1
        kclr(i) = i
        enddo
        return
      endif
!
      if( cinp(1:1).eq."d" ) then
        cinp = cdef
        mji = lenx(cinp)
        write(6,'(2x,"cinp =",a)') cinp(1:lenx(cinp))
        write(6,'(2x,"ncl+1 =",i3)') ncl+1
      endif
!
      maxclr = 30
      call gdfps( maxclr )
      read(cinp(1:mji),*) (kclr(i),i=1,ncl+1)
      write(6,'(2x,"kclr =",20i4)') (kclr(i),i=1,ncl+1)
!
      return
      end
