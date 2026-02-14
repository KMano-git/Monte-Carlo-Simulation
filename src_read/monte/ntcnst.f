!**********************************************************************aion
      subroutine ntcnst
!**********************************************************************
!
!   physical constants to be used in neut2d code
!
!----------------------------------------------------------------------
      use cntcom, only : fnfi, fnor, fnth, iransd, tcfi, tcth, tnor
     >    , tsfi, tsth, wmate
      use cntctl, only : mdl_hrt
      use cplcom, only : aion
      use csonic, only : mstop
      use cunit,  only : n6
      implicit none
!
      integer  i, ith, ifi, ivx
! modified 2/1 lines organize local variables and include files by kamata 2021/06/16
!ik   real*8   x, random, ctx, stx, cfx, sfx, vlx
!ik   real*8   rnlast
      real*8   x, ctx, stx, cfx, sfx, vlx
      integer  ldbg
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   character clin*120, crfl*20
      character crfl*20
      data ldbg/1/
! added 2 lines organize local variables and include files by kamata 2021/06/16
! function
      real(8)    random, rnlast
!
!----------------------------------------------------------------------
!::pi,cev,cmp
!----------------------------------------------------------------------
      call phcnst
      call phycnst
!
!----------------------------------------------------------------------
!::random
!----------------------------------------------------------------------
      do 110 i = 1, iransd
      x = random(0)
 110  continue
      write(n6,'(4x,"initial random number iransd =",i10,2f20.12)')
     >   iransd, x, rnlast(0)
!
!----------------------------------------------------------------------
!::table
!----------------------------------------------------------------------
      call mtbnrm
      call mtbtri
      write(n6,'(4x,"call mtbnrm : normal random number")')
      write(n6,'(4x,"call mtbtri : tri function")')
!
!----------------------------------------------------------------------
!::reflection
!----------------------------------------------------------------------
      crfl = " "
      if( dabs(aion(1)-2.0).le.0.001 ) crfl = " D -> "
      if( dabs(aion(1)-1.0).le.0.001 ) crfl = " H -> "
      if( wmate(1:1).eq."C"  ) crfl(7:) = "C"
      if( wmate(1:2).eq."FE" ) crfl(7:) = "FE"
      if( wmate(1:1).eq."W"  ) crfl(7:) = "W"
      write(n6,'(4x,"call rflini : reflection data =",a)') crfl
      if( index(crfl,">").eq.0 .or. crfl(7:7).eq." " ) then
        mstop = 1
        return
!xx     call wexit("ntcnst","invarid plasma or wall material")
      endif
      call rflini( crfl )
!xxx  call test_refl
!
!----------------------------------------------------------------------
!::ATOMIC/MOLECULAR data
!----------------------------------------------------------------------
      call ntcros
!
!::atomic hudrogen cross section
      write(n6,'(4x)')
      write(n6,'(4x,"call atomhy : ionization & recombination",
     >  " data of Hydrogen  degas2  00/08/29")')
      call atomhy
!
!::molecular hydrogen cross section
      write(n6,'(4x,"call molehy : cross section ",
     >  "data of molecular hydrogen  degas2  00/09/20")')
      call molehy
!
!::elastic collision
      write(n6,'(4x,"call elinit : elastic collision  02/01/28")')
      call elinit
!
!:: H radiation trapping by Sawada-sensei
      if(mdl_hrt .eq. 1) then
        write(n6, '(4x,"call hrt_read : H radiation trapping 150814")')
        call hrt_read
      endif
!
!----------------------------------------------------------------------
!::debug
!----------------------------------------------------------------------
      if( ldbg.eq.1 ) then
      write(n6,'(/2x,"*** ntcnst ***  mtbnrm & mtbtri")')
      write(n6,'(2x,2x,"i",3x,"vlx",9x,"cth",9x,"sth",9x,"cfi",9x,
     >  "sfi")')
      do 120 i = 1, 10
      ith = int(fnth*random(0) + 1.0d0)
      ctx = tcth(ith)
      stx = tsth(ith)
      ifi = int(fnfi*random(0) + 1.0d0)
      cfx = tcfi(ifi)
      sfx = tsfi(ifi)
      ivx = int(fnor*random(0) + 1.0d0)
      vlx = tnor(ivx)
      write(n6,'(2x,i3,1p5e12.3)')
     >  i,vlx,ctx,stx,cfx,sfx
 120  continue
      endif
!
!----------------------------------------------------------------------
!::cross section data ( mfp : mean free path )
!----------------------------------------------------------------------
!xx   call ntrmfp
!xx   call test_mtbnrm
!
      return
      end
