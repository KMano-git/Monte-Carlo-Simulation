!**********************************************************************
      subroutine plwrdt(nft,kio,cdsn)
!**********************************************************************
!
!xx      common /cmwrad/ wtime, witim, wnlp, wfac,
!xx     >   wcre, wmci, wmce, wimi, wime
!xx
!        soldor  wcre, wmce  < 0.0
!        IMPMC   twrd, wfac  > 0.0
!
!----------------------------------------------------------------------
      use cplwrd,  only : wcre, wfac, wime, wimi, witim, wmce, wmci
     >    , wnlp, wtime
      use cpmpls,  only : nty_ctl
      use cplwrd2, only : wmc_zwe, wmc_zwi
      use csize,   only : nzsmx, ndmc
      use cunit,   only : n6
      implicit none
!
!::argument
      integer, intent(in)      :: nft, kio
      character(*), intent(in) :: cdsn     ! dummy
!
!::local variables
      real*8 :: wmc_zwi_tmp(ndmc,nzsmx), wmc_zwe_tmp(ndmc,nzsmx)
      integer i
!
      if( kio.eq.0 ) return
!
      if( kio.eq.1 ) then
        do i = 1, nzsmx
          if(i == nty_ctl) then
            wmc_zwi_tmp(:,i) = wmc_zwi(:,i)*wfac
            wmc_zwe_tmp(:,i) = wmc_zwe(:,i)*wfac
          else
            wmc_zwi_tmp(:,i) = wmc_zwi(:,i)
            wmc_zwe_tmp(:,i) = wmc_zwe(:,i)
          endif
        enddo
        write(n6,'("wtime, witim, wnlp, wfac = ",(f13.8,f13.4,2f13.8))')
     >            wtime, witim, wnlp, wfac
        write(nft)  wtime, witim, wnlp, wfac
        write(nft)  wcre, wmci, wmce, wimi, wime
        write(nft)  wmc_zwi_tmp(:,:)
        write(nft)  wmc_zwe_tmp(:,:)
      endif
!
      if( kio.eq.2 ) then
        read(nft)  wtime, witim, wnlp, wfac
        write(n6,*) 'wtime, witim, wnlp, wfac',
     >     wtime, witim, wnlp, wfac
        read(nft)  wcre, wmci, wmce, wimi, wime
        read(nft)  wmc_zwi(:,1),wmc_zwi(:,2:nzsmx)
        read(nft)  wmc_zwe(:,1),wmc_zwe(:,2:nzsmx)
      endif
      close(nft)
!
      return
      end
