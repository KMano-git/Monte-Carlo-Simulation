!**********************************************************************
      subroutine ntcler
!**********************************************************************
!
!     clear the scoreing variables
!
!     delete fluid scoring  sn_mf, sp_mf, we_mf, wi_mf   2010/12/21
!
!----------------------------------------------------------------------
      use cntcom,  only : flxin
      use cntmnt,  only : dotn
      use cntscg,  only : ndrmx, sn_mc, sn_mt, sp_mc, sp_mt, tsn_mc
     >    , tsn_mt, tsp_mc, tsp_mt, twe_mc, twe_mt, twi_mc, twi_mt
     >    , we_mc, we_mt, wi_mc, wi_mt
      use cntwedf, only : wedf, wedf_line
      use cntwcn,  only : wden, weema, weemm, wees, wehta, wehtm, wemt
     >    , wend, weng, wfhta, wfhtm, wgden, wgeng, whta, whtm, wnfl0x
     >    , wnfl0y, wnfl0z, wnflgx, wnflgy, wnflgz, wsbr, wtion, wtot
     >    , wvlp, wwal
      use csize,   only : ndgs, ndmc, ndsp, ndwp
      implicit none
!
!::local variables
      integer  nsiz, nsizg
!
!::neutal flux into the plasma when calculete neutral profile
      dotn  = flxin
!
!::source at recombination
      nsiz = (ndmc+1)*ndgs
      call setd( wsbr, nsiz, 0.0d0 )
!
!::neutral density and temperature
      nsiz = (ndmc+1)*ndgs
      call setd( wden,   nsiz, 0.0d0 )
      call setd( weng,   nsiz, 0.0d0 )
      call setd( wvlp,   nsiz, 0.0d0 )
      call setd( wnfl0x, nsiz, 0.0d0 ) ! neutral flow 160623 toku
      call setd( wnfl0y, nsiz, 0.0d0 )
      call setd( wnfl0z, nsiz, 0.0d0 )
      nsiz = (ndmc+1)*2
      call setd( wgden,  nsiz, 0.0d0 )
      call setd( wgeng,  nsiz, 0.0d0 )
      call setd( wnflgx, nsiz, 0.0d0 )
      call setd( wnflgy, nsiz, 0.0d0 )
      call setd( wnflgz, nsiz, 0.0d0 )
!
!::weit
      wtot = 0.0d0
      nsiz = 10
      call setd( wend, nsiz, 0.0d0 )
      nsiz = (ndmc+1)*ndgs
      call setd( wtion, nsiz, 0.0d0 )
      nsiz = ndwp
      call setd( wemt, nsiz, 0.0d0 )  ! No use ==> wsbr
      call setd( wwal, nsiz, 0.0d0 )
      nsiz = 30*3
      call setd( whta, nsiz, 0.0d0 )
      call setd( whtm, nsiz, 0.0d0 )
!
!::source terms
      nsizg = (ndmc+1)*ndsp
      nsiz = ndmc + 1
!
!xx   call setd( sn_mf, nsizg, 0.0d0 )
!xx   call setd( sp_mf, nsizg, 0.0d0 )
!xx   call setd( we_mf, nsiz,  0.0d0 )
!xx   call setd( wi_mf, nsiz,  0.0d0 )
!
      call setd( sn_mt, nsizg, 0.0d0 )
      call setd( sp_mt, nsizg, 0.0d0 )
      call setd( we_mt, nsiz,  0.0d0 )
      call setd( wi_mt, nsiz,  0.0d0 )
!
      call setd( sn_mc, nsizg, 0.0d0 )
      call setd( sp_mc, nsizg, 0.0d0 )
      call setd( we_mc, nsiz,  0.0d0 )
      call setd( wi_mc, nsiz,  0.0d0 )
!
!::total source for each reaction
      tsn_mt(1:ndrmx,1:10) = 0.0d0
      tsp_mt(1:ndrmx,1:10) = 0.0d0
      twe_mt(1:ndrmx,1:10) = 0.0d0
      twi_mt(1:ndrmx,1:10) = 0.0d0
      tsn_mc(1:ndrmx,1:10) = 0.0d0
      tsp_mc(1:ndrmx,1:10) = 0.0d0
      twe_mc(1:ndrmx,1:10) = 0.0d0
      twi_mc(1:ndrmx,1:10) = 0.0d0
!
!::neutral flux onto the wall
      nsiz = ndwp
      call setd( wfhta, nsiz, 0.0d0 )
      call setd( wfhtm, nsiz, 0.0d0 )
      call setd( wehta, nsiz, 0.0d0 )
      call setd( wehtm, nsiz, 0.0d0 )
      call setd( weema, nsiz, 0.0d0 )
      call setd( weemm, nsiz, 0.0d0 )
!
      wees = 0.0d0
      wedf = 0.0d0
      wedf_line = 0.0d0
!
      return
      end
