!**********************************************************************
      subroutine plwrdcrG(swrd)
!**********************************************************************
!
!     pl(plasma) wrd(radiation) cr(corona) G(general)
!     swrd : total radiation
!
!----------------------------------------------------------------------
      use cplcom, only : nlp, nlpmn
      use cplwrd, only : wcr_cnc, wcr_nty, wcr_typ, wcr_wrd
      use csize,  only : ndmc
      use csonic, only : itim, lpost
      use cunit,  only : n6
      implicit none

!::argument
      real(8) :: swrd(ndmc)

!::local variables
      real(8) :: wrd(ndmc)
      real(8) :: cimz(10)
      real(8) :: wrsm, wrgn(10)
      integer :: i, iz, ldbg, ir
      character(4) :: cnam

      ldbg = 0
      if( mod(itim,100).eq.0  .and. nlp.eq.nlpmn ) ldbg = 1
      if( lpost.eq.1 ) ldbg = 1

      swrd(1:ndmc) = 0.0d0

      do i = 1, wcr_nty
        wrd(1:ndmc) = 0.0d0
        iz = i
        cimz(1:10) = wcr_cnc(1:10,iz)
        cnam = wcr_typ(iz)

            if( trim(cnam) == "C" ) then
                    call wrdcr_C(cimz,wrd)
        elseif( trim(cnam) == "W" ) then
                    call wrdcr_W(cimz,wrd)
        elseif( trim(cnam) == "W2" ) then
                    call wrdcr_W2(cimz,wrd)
        elseif( trim(cnam) == "Be" ) then
                    call wrdcr_Be(cimz,wrd)
        elseif( trim(cnam) == "Ar" ) then
                    call wrdcr_Ar(cimz,wrd)
        elseif( trim(cnam) == "Ne" ) then
                    call wrdcr_Ne(cimz,wrd)
        elseif( trim(cnam) == "Kr" ) then
                    call wrdcr_Kr(cimz,wrd)
        elseif( trim(cnam) == "Xe" ) then
                    call wrdcr_Xe(cimz,wrd)
        else
          call wexit("plwrdcrG","no found imp-name "//trim(cnam))
        endif
        wcr_wrd(1:ndmc,iz) = wrd(1:ndmc)
        swrd(1:ndmc) = swrd(1:ndmc) + wrd(1:ndmc)

!::debug check
        if( ldbg == 1 ) then
          call wrintg(wrd,wrsm,wrgn)
          write(n6,'(2x,"wrdcrG_",a2,2x,1pe14.6,1x,1p7e11.3,
     >   "  cimz(1) =",1pe11.3)')
     >      cnam(1:2),wrsm,(wrgn(ir),ir=1,7),cimz(1)
        endif

      enddo

      return
      end
