!***********************************************************************
      subroutine pst_hlod
!***********************************************************************
      use cplmet, only : itmpe, itpve, itsle, itsls, jcdp2, jcxp1, jcxp2
      use cpmpls, only : jmd1, jmd2
      use cunit,  only : n6
      implicit none
!
!::local variables
      integer  ix, it
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   integer  i6, jmax, jw
!
      write(n6,'(/2x,"*** pst_hlod *** ")')
      write(n6,'(2x,"type",6x,"ix",2x,"kdr",4x,"Fi",10x,"Qt",
     >   10x,"Qi",10x,"Qe")')
!
!::o-mid
      ix = jmd1
      call plxflx("sol",ix,+1)
!::i-mid
      ix = jmd2
      call plxflx("sol",ix,+1)
!::o-div
      ix = 1
      call plxflx("odv",ix,+1)
!::i-div  (error)
      ix = jcdp2
      call plxflx("idv",ix,-1)  ! <==
!::o-Xp (SOL)
      ix = jcxp1+1
      call plxflx("sol",ix,-1)
!::o-Xp (Div)
      ix = jcxp1-1
      call plxflx("odv",ix,-1)
!::i-Xp (SOL)
      ix = jcxp2-1
      call plxflx("sol",ix,+1)
!::i-Xp (Div)
      ix = jcxp2+1
      call plxflx("idv",ix,+1)
!
      write(n6,'(2x)')
      write(n6,'(2x,"type",6x,"it",2x,"kdr",4x,"Fi",10x,"Qt",
     >   10x,"Qi",10x,"Qe")')
!
!::edge
      it = itmpe
      call plyflx("edg",it,-1)
!::spx
      it = itsle
      call plyflx("sol",it,+1)
!::sol
      it = itsls
      call plyflx("sol",it,+1)
!::o-div
      it = itsls
      call plyflx("odv",it,+1)
!::i-div
      it = itsls
      call plyflx("idv",it,+1)
!::prv
      it = itpve
      call plyflx("prv",it,-1)
!
!PSI2014
      it = itsle
      call plyflx("odv",it,  +1)
      call plyflx("odv",it-1,+1)
      call plyflx("idv",it,  +1)
      call plyflx("idv",it-1,+1)
!
!**************************************************
!Heat flux tube plot data(2018/12/17 Y.HOMMA)
      call heat_flux_along_flux_tube(24)
!**************************************************


      write(n6,'(2x)')
      return
      end
