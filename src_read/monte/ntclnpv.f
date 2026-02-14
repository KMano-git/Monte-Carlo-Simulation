!**********************************************************************
      subroutine ntclnpv
!**********************************************************************
      use cntcom, only : evel, evelb, icpo, ievt, igas, igasb, ihit
     >    , ikon, ilon, iptl, istp, ixpo, iypo, izpo, tmas, tmasb, vel
     >    , vel2, vel2b, velb, velp, velpb, velx, velxb, vely, velyb
     >    , velz, velzb, weit, weitb, xpos, ypos, zpos
      use cunit,  only : n6
      implicit none
!
!----------------------------------------------------------------------
!::particle variables
!----------------------------------------------------------------------
!x      real*8    xpos,ypos,zpos
!x     >         ,velx,vely,velz,velp,vel,vel2,evel,tmas,weit
!x      integer   icpo,ixpo,iypo,izpo,ilon,igas,istp,ievt,ihit,iptl
!x      real*8    velxb,velyb,velzb,velpb,velb,vel2b,evelb,tmasb,weitb
!x      integer   igasb
!
!x      common /cntptl/ xpos,ypos,zpos
!x     >          ,velx,vely,velz,velp,vel,vel2,evel,tmas,weit
!x     >          ,icpo,ixpo,iypo,izpo,ilon,igas,istp,ievt,ihit,iptl
!x      common /cntptb/
!x     >     velxb,velyb,velzb,velpb,velb,vel2b,evelb,tmasb,weitb,igasb
!----------------------------------------------------------------------
!
      write(n6,'(/2x,"*** ntclnpv ***  clear particle variables")')
!
      xpos=0.0d0; ypos=0.0d0; zpos=0.0d0
      velx=0.0d0; vely=0.0d0; velz=0.0d0
      velp=0.0d0; vel =0.0d0; vel2=0.0d0; evel=0.0d0
      tmas=0.0d0; weit=0.0d0
      icpo=0;     ixpo=0;     iypo=0;     izpo=0
      ilon=0;     ikon=0;     igas=0;     istp=0;     ievt=0
      ihit=0;     iptl=0
      velxb=0.0d0; velyb=0.0d0; velzb=0.0d0
      velpb=0.0d0; velb =0.0d0; vel2b=0.0d0; evelb=0.0d0
      tmasb=0.0d0; weitb=0.0d0; igasb=0
!
      return
      end
