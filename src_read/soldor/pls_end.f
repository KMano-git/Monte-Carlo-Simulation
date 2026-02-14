!***********************************************************************
      subroutine pls_end
!***********************************************************************
      use cplcom, only : cftw, nftw
      use csonic, only : lntl, lqick
      use cunit,  only : cdrout, lmspe, lmype, n6
      use dbg_mod,     only : cpaths

      implicit none
!
      integer nf
      integer krad, kntl
! function
      integer    lenx
!
!
      krad = 1
      kntl = 1
      if( lqick.eq.1 ) then
        krad = 0
        kntl = 0
      endif
!
      write(n6,'(/2x,"*** epilog ***   lqick =",i2,"  krad =",i2,
     >  "  kntl =",i2)')   lqick, krad, kntl
!
!::flux at core edge
      call plqcon_edge
      call plmflx
!
!::mapping data
      call pst_rprf
      call pst_bprf
!
!::heat flux onto the plate
      if( lntl.eq.3 ) then
      call ntwdeg
      call pst_wrad(krad)
      call pst_wntl(kntl)
      call pst_qdpl( cpaths, ' ' )
      call pst_hlod
      endif
!
      call pldisk(nftw,1,cftw)
      call plhist(3)
!
      if( lmype.eq.lmspe ) then
        nf = 21
        open( unit=nf, file=cdrout(1:lenx(cdrout))//"outend" )
        call pllist(nf)
        close (nf)
        call plsprf    ! <== compare prf data
      endif

      call write_plasma_source_terms
      end subroutine pls_end

!***********************************************************************
      subroutine write_plasma_source_terms
!***********************************************************************
      use cntmnt, only: ssn, ssp, swi, swe
      use cntcom, only: iplx, iply
      use cntcom, only: ncmax
      implicit none
!
!::local variables
      integer i6, ic, j, i
!
      i6 = 21
      open(unit=i6,file="pla_source.txt")
      write(i6,'(" *   ic",4x,"j",5x,"i",4x,
     > "ssn",9x,"ssp",9x,"swi",9x,"swe")')
      do ic = 1, ncmax
        j = iplx(ic)
        i = iply(ic)
        if(j>0 .and. i>0) then
          write(i6,'(i6,1x,i5,1x,i5,4(1x,1pe11.3))') 
     >     ic,j,i,ssn(j,i,1),ssp(j,i,1),swi(j,i),swe(j,i)
        else
          write(i6,'(i6,1x,i5,1x,i5,4(1x,1pe11.3))') 
     >     ic,j,i,0.0d0,0.0d0,0.0d0,0.0d0
        endif
      enddo
      close(i6)
      end subroutine write_plasma_source_terms
