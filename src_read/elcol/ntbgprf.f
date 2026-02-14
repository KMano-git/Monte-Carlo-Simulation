      MODULE BGProf
      use csize,  only : ndgs, ndmc, ndmsr
      use cntmnt, only : pfl_ntl, vkflx, vmwork
      use cunit,  only : n6
      IMPLICIT NONE

      integer,parameter :: nMolType=2
      real(8), private, allocatable :: wden0BG(:,:,:), wdengBG(:,:,:)
     >    , weng0BG(:,:,:), wenggBG(:,:,:), wnfl0BGx(:,:,:)
     >    , wnfl0BGy(:,:,:), wnfl0BGz(:,:,:), wnflgBGx(:,:,:)
     >    , wnflgBGy(:,:,:), wnflgBGz(:,:,:),  wvlp0BG(:,:,:)
     >    , wvlpgBG(:,:,:)

! BackGround Profile from previous calculation
      real(8), public, allocatable :: den0BG(:,:), dengBG(:,:)
     >    , eng0BG(:,:), enggBG(:,:), nfl0BGx(:,:), nfl0BGy(:,:)
     >    , nfl0BGz(:,:), nflgBGx(:,:), nflgBGy(:,:), nflgBGz(:,:)
     >    , vlp0BG(:,:), vlpgBG(:,:)

      real(8),public :: swnrmForEachSource(ndmsr)

      PUBLIC :: ntwprof,ntbgprof

      CONTAINS

! --------
! conserve weight profile
      Subroutine ntwprof
        IMPLICIT NONE

! BackGroung profile prepared for n-n collision & photon transport, 2015.07.27, by toku
!ik     INTEGER,parameter :: nwkmp_wden_sta=1 + (ndmc+1)*(ndgs*2+2+ndgs) ! starting point of wden. see monte/inc/cntwcn
        integer    nwkmp_wden_sta ! starting point of wden. see monte/inc/cntwcn
        INTEGER :: ivmwork, iSource, iSpecies, iCell, ivmwork0

        nwkmp_wden_sta = 1 + (ndmc+1) * (ndgs*2+2+ndgs)
        if( .not. allocated( wden0BG ) )
     >    allocate( wden0BG(0:ndmc,ndgs,ndmsr)
     >      , wdengBG(0:ndmc,nMolType,ndmsr), weng0BG(0:ndmc,ndgs,ndmsr)
     >      , wenggBG(0:ndmc,nMolType,ndmsr)
     >      , wnfl0BGx(0:ndmc,ndgs,ndmsr), wnfl0BGy(0:ndmc,ndgs,ndmsr)
     >      , wnfl0BGz(0:ndmc,ndgs,ndmsr)
     >      , wnflgBGx(0:ndmc,nMolType,ndmsr)
     >      , wnflgBGy(0:ndmc,nMolType,ndmsr)
     >      , wnflgBGz(0:ndmc,nMolType,ndmsr)
     >      , wvlp0BG(0:ndmc,ndgs,ndmsr)
     >      , wvlpgBG(0:ndmc,nMolType,ndmsr) )

!::Pickup wden0, wdeng etc. from vmwork

       Lp_Source: do iSource = 1, ndmsr

       ivmwork0 = nwkmp_wden_sta
       ivmwork  = nwkmp_wden_sta

! wden
       do iSpecies = 1, ndgs
          do iCell = 0, ndmc
             ivmwork = ivmwork + 1
             wden0BG(iCell,iSpecies,iSource) = vmwork(ivmwork,iSource)
          enddo
       enddo

! weng
       do iSpecies = 1, ndgs
          do iCell = 0, ndmc
             ivmwork = ivmwork + 1
             weng0BG(iCell,iSpecies,iSource) = vmwork(ivmwork,iSource)
          enddo
       enddo

! wvlp
       do iSpecies = 1, ndgs
          do iCell = 0, ndmc
             ivmwork = ivmwork + 1
             wvlp0BG(iCell,iSpecies,iSource) = vmwork(ivmwork,iSource)
          enddo
       enddo

! nflow 160623
       do iSpecies = 1, ndgs
          do iCell = 0, ndmc
             ivmwork = ivmwork + 1
             wnfl0BGx(iCell,iSpecies,iSource) = vmwork(ivmwork,iSource)
          enddo
       enddo
       do iSpecies = 1, ndgs
          do iCell = 0, ndmc
             ivmwork = ivmwork + 1
             wnfl0BGy(iCell,iSpecies,iSource) = vmwork(ivmwork,iSource)
          enddo
       enddo
       do iSpecies = 1, ndgs
          do iCell = 0, ndmc
             ivmwork = ivmwork + 1
             wnfl0BGz(iCell,iSpecies,iSource) = vmwork(ivmwork,iSource)
          enddo
       enddo

! wtion (SKIP)
       ivmwork = ivmwork+ndgs*(ndmc+1)

! wdeng
       do iSpecies = 1, nMolType
          do iCell = 0, ndmc
             ivmwork = ivmwork + 1
             wdengBG(iCell,iSpecies,iSource) = vmwork(ivmwork,iSource)
          enddo
       enddo

! wengg
       do iSpecies = 1, nMolType
          do iCell = 0, ndmc
             ivmwork = ivmwork + 1
             wenggBG(iCell,iSpecies,iSource) = vmwork(ivmwork,iSource)
          enddo
       enddo

! nflow 160623
       do iSpecies = 1, nMolType
          do iCell = 0, ndmc
             ivmwork = ivmwork + 1
             wnflgBGx(iCell,iSpecies,iSource) = vmwork(ivmwork,iSource)
          enddo
       enddo
       do iSpecies = 1, nMolType
          do iCell = 0, ndmc
             ivmwork = ivmwork + 1
             wnflgBGy(iCell,iSpecies,iSource) = vmwork(ivmwork,iSource)
          enddo
       enddo
       do iSpecies = 1, nMolType
          do iCell = 0, ndmc
             ivmwork = ivmwork + 1
             wnflgBGz(iCell,iSpecies,iSource) = vmwork(ivmwork,iSource)
          enddo
       enddo
! wvlpg (Currently Unavailable)
       wvlpgBG=0.d0

      enddo Lp_Source
      Return
      END Subroutine ntwprof

! ============================================================

      Subroutine ntbgprof
      use cntcom,      only : rmas, volm
      use com_phycnst, only : cev
      IMPLICIT NONE

      REAL(8) :: tfsrc=0.d0, zfce, swnrm1=0.d0, zmas
      real(8) :: vol, denBG
      INTEGER :: iSource, iSpecies, iCell

!::Calc. profiles

!::reset
      den0BG =0.d0;eng0BG =0.d0;vlp0BG =0.d0
      dengBG =0.d0;enggBG =0.d0;vlpgBG =0.d0
      nfl0BGx=0.d0;nfl0BGy=0.d0;nfl0BGz=0.d0
      nflgBGx=0.d0;nflgBGy=0.d0;nflgBGz=0.d0

! den0, eng0, vlp0
      do iSource = 1, ndmsr
       swnrm1 = swnrmForEachSource(iSource)
       if(swnrm1==0.0d0)then
         write(n6,'(/2x,"*** ntbgprf  skip iSource ",i3)') iSource
         cycle
       endif
       tfsrc  = pfl_ntl(vkflx(iSource))
       do iSpecies = 1, ndgs
          zfce=2.d0/3.d0*0.5d0*rmas(iSpecies)/cev
          do iCell=1,ndmc
             vol = volm(iCell)
             denBG = wden0BG(iCell,iSpecies,iSource)

             if( vol>0.0d0 .and. denBG > 0.0d0) then

               den0BG(iCell,iSpecies) = den0BG(iCell,iSpecies)          
     &               + tfsrc/swnrm1/vol * denBG

               eng0BG(iCell,iSpecies) = eng0BG(iCell,iSpecies)          
     &               + zfce * weng0BG(iCell,iSpecies,iSource) / denBG

               vlp0BG(iCell,iSpecies) = vlp0BG(iCell,iSpecies)          
     &               + wvlp0BG(iCell,iSpecies,iSource) / denBG

               nfl0BGx(iCell,iSpecies) = nfl0BGx(iCell,iSpecies)        
     &               + wnfl0BGx(iCell,iSpecies,iSource) / denBG

               nfl0BGy(iCell,iSpecies) = nfl0BGy(iCell,iSpecies)        
     &               + wnfl0BGy(iCell,iSpecies,iSource) / denBG

               nfl0BGz(iCell,iSpecies) = nfl0BGz(iCell,iSpecies)        
     &               + wnfl0BGz(iCell,iSpecies,iSource) / denBG
             endif

          enddo
       enddo
      enddo

! deng, engg, vlpg
      do iSource = 1, ndmsr
         swnrm1 = swnrmForEachSource(iSource)
         if(swnrm1==0.0d0)then
            write(n6,'(2x,"*** ntbgprf  skip iSource ",i3)') iSource
            cycle
         endif
         tfsrc  = pfl_ntl(vkflx(iSource))
         do iSpecies = 1, ndgs
            zmas=2.d0*rmas(1)   !!!!! This must be changed when Tritium is taken into account.
            zfce=2.d0/3.d0*0.5d0*zmas/cev
            do iCell=1,ndmc
               vol = volm(iCell)
               denBG = wdengBG(iCell,iSpecies,iSource)

               if( vol>0.0d0 .and. denBG > 0.0d0) then

                  dengBG(iCell,iSpecies) = dengBG(iCell,iSpecies)       
     &                 + tfsrc/swnrm1/vol * denBG

                  enggBG(iCell,iSpecies) = enggBG(iCell,iSpecies)       
     &                 + zfce * wenggBG(iCell,iSpecies,iSource) / denBG

                  vlpgBG(iCell,iSpecies) = vlpgBG(iCell,iSpecies)       
     &                 + wvlpgBG(iCell,iSpecies,iSource) / denBG

                  nflgBGx(iCell,iSpecies) = nflgBGx(iCell,iSpecies)     
     &                 + wnflgBGx(iCell,iSpecies,iSource) / denBG

                  nflgBGy(iCell,iSpecies) = nflgBGy(iCell,iSpecies)     
     &                 + wnflgBGy(iCell,iSpecies,iSource) / denBG

                  nflgBGz(iCell,iSpecies) = nflgBGz(iCell,iSpecies)     
     &                 + wnflgBGz(iCell,iSpecies,iSource) / denBG
               endif

            enddo
         enddo
      enddo

      Return
      END Subroutine ntbgprof

      subroutine allocate_ntbgprof
      if( .not. allocated( den0BG ) )
     >  allocate( den0BG(0:ndmc,ndgs), dengBG(0:ndmc,nMolType)
     >    , eng0BG(0:ndmc,ndgs), enggBG(0:ndmc,nMolType)
     >    , nfl0BGx(0:ndmc,ndgs), nfl0BGy(0:ndmc,ndgs)
     >    , nfl0BGz(0:ndmc,ndgs), nflgBGx(0:ndmc,nMolType)
     >    , nflgBGy(0:ndmc,nMolType), nflgBGz(0:ndmc,nMolType)
     >    , vlp0BG(0:ndmc,ndgs), vlpgBG(0:ndmc,nMolType) )
      den0BG =0.d0;eng0BG =0.d0;vlp0BG =0.d0
      dengBG =0.d0;enggBG =0.d0;vlpgBG =0.d0
      nfl0BGx=0.d0;nfl0BGy=0.d0;nfl0BGz=0.d0
      nflgBGx=0.d0;nflgBGy=0.d0;nflgBGz=0.d0
      end subroutine allocate_ntbgprof

! --------
      END Module BGProf
