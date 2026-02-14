      subroutine hrt_set
      use BGProf, only : den0bg, eng0bg ! monte/ntbgprf.f90
      use chrt,   only : sn0min_hrt
      use cntcom, only : ncmax, tbel_ei, tbsg_ei, iplx, iply
      use cntpls, only : dene, teme
      use cunit,  only : n6, lmype, lmspe
      use cntctl, only : hrt_dbg
      implicit none

      integer  ic
      real*8   zne, zte
      real*8   sghrt, bg_den0, bg_deng, bg_eng0, sghrt_div
      real*8   sgei
! function
      real(8)    hrt_eval
      logical  ldbg 
      integer :: ixsv
      ldbg = .false.

      write(n6,*)
      write(n6,'(1x,"HRT: set table")')
!      write(n6,*)
!      write(n6,'(1x,3x,"ic",3x,"sgei",12x,"sghrt",11x,"sghrt/sgei")')

      if( lmype.eq.lmspe ) then
        if(hrt_dbg>0) ldbg=.true.
      endif
      
      if(ldbg) open(999,file="hrt.txt")

      do ic = 1, ncmax
        zne = dene(ic)
        zte = teme(ic)
        bg_eng0 = eng0BG(ic,1)
        bg_den0 = den0BG(ic,1)
        !if(zte .gt. stemax) cycle
        !if(bg_eng0 .gt. st0max) cycle
        if(bg_den0 .lt. sn0min_hrt .and. .not. ldbg) cycle
!        if(bg_den0 .gt. sn0max) cycle
        sghrt = hrt_eval("alp",zte,zne,bg_eng0, bg_den0, bg_deng)
        sghrt = dmax1(zne*sghrt,1.0d-99)
        if(tbsg_ei(ic).ne.0.0d0) then
          sghrt_div = sghrt/tbsg_ei(ic)
        else
          sghrt_div = 0.0d0
        endif
        if(mod(ic,50)==0)then
           write(n6,'(1x,"HRT",3i5,2e16.5,f16.5)')
     >       ic,iplx(ic),iply(ic),tbsg_ei(ic),sghrt,sghrt_div
        endif
        if(ldbg)then
        if(iplx(ic)/=0)then 
          if(ixsv/=iplx(ic)) write(999,*)
          ixsv=iplx(ic)
          write(999,'(3i5, 2(1pe16.5), 0pf16.7, 4(1pe16.5))')
     >            ic,iplx(ic),iply(ic)
     >            ,tbsg_ei(ic),sghrt,sghrt_div
     >            ,zne, zte, bg_den0, bg_eng0
        endif ! iplx/=0
        endif ! ldbg
        sgei = tbsg_ei(ic)
        tbsg_ei(ic) = sgei + sghrt
        tbel_ei(ic) = (tbel_ei(ic)*sgei + 3.398d0*sghrt )/(sgei+sghrt)
        !tbel_ei(ic) = (tbel_ei(ic)*sgei + 13.595*sghrt )/(sgei+sghrt)
      enddo

      if(ldbg)close(999)

      return
      end subroutine
