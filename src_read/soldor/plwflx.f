!**********************************************************************
      subroutine plwflx
!**********************************************************************
!
!    plasma flux onto the wall
!
!----------------------------------------------------------------------
      use cntcom, only : imwl, npew, npsw
      use cntwfl, only : xare, xfhti, xtfi
      use cplcom, only : flps
      use csize,  only : ndmsr, ndwp
      use cunit,  only : cprg
      use mod_externalgrid, only : use_exdata
      implicit none
!
!::local variables
      integer nsiz, ia, nw, iw, iws, iwe, iwp
      real*8  zfi, zsm
!
!::check flps/tfps(flux on boundary) to be positive
      if(cprg /= 'plimp') then
        call check_flps
      endif
!
!::clear
      nsiz = ndwp*(ndmsr+1)
      call setd( xfhti, nsiz, 0.0d0 )
      nsiz = ndmsr+1
      call setd( xtfi, nsiz, 0.0d0 )
!
      ia = 1
!
      do nw = 2, 4, 2 ! only for inner diverter(nw=2) or outer diverter(nw=4)
        zsm = 0.0d0
        iws = npsw(nw)
        iwe = npew(nw)
        do iw = iws, iwe
          if(.not.use_exdata .and. iw==iwe) cycle
          iwp = imwl(iw) !  wall => plasma surface
          if( iwp.gt.0 ) then
            zfi = flps(iwp,ia)/xare(iw)
            xfhti(iw,nw) = zfi
            xfhti(iw,0) = xfhti(iw,0) + zfi
            zsm = zsm + xfhti(iw,nw)*xare(iw)
          endif
        enddo
        xtfi(nw) = zsm
        xtfi(0) = xtfi(0) + zsm
      enddo
!
      return
      end subroutine plwflx
!
!**********************************************************************
      subroutine check_flps
!**********************************************************************
!
!    check flux on the boundary have postive value.
!    calculate in soldor/plpflx.f
!
!----------------------------------------------------------------------
      use cplcom, only : flps, tfps, nion ! flps:flux each cell / tfps:total flux of boundary
      use cntcom, only : npsp, npep
      use cunit,  only : n6
      use csonic, only : time
      implicit none
!
!::local variables
      integer ia,nw,iw
      character cend(5)*3
      data cend/"sol","idp","prv","odp","man"/

!:: ia(fluid species)
!:: nw(boundary): 1(sol),2(idp),3(prv),4(odp),5(man)
!:: iw(each cell)

      do ia = 1, nion
        !:: flux must positive for nw = 1(sol),2(idp),3(prv),4(odp)

        !:: nw = 1 (SOL)
        nw = 1
        ! for SOL boundary, flux=0 at initial condition.
        if(time.gt.0.0d0) then
          do iw = npsp(nw), npep(nw)-1
            if(flps(iw,ia).le.0.0d0) then
              ! you may call wexit
              write(n6,'("Waring:flux on boundary <= 0 :iw"
     >         ,i6,",",a,",",1pe12.3)') iw,cend(nw),flps(iw,ia)
            endif
          enddo
          !:: check total flux
          if(tfps(nw,ia).le.0.0d0) then
            write(n6,'("Error:total flux on boundary <=0 :",a,":"
     >        ,1pe12.3)') cend(nw),tfps(nw,ia)
            call wexit("check_flps","total flux on boundary <= 0")
          endif
        else
          write(n6,'("total flux on SOL boundary =",1pe12.3,
     >      ", time=",1pe12.3)') tfps(nw,ia),time
        endif

        do nw = 2, 4
          do iw = npsp(nw), npep(nw)-1
            if(flps(iw,ia).le.0.0d0) then
              ! you may call wexit
              write(n6,'("Waring:flux on boundary <= 0 :iw"
     >         ,i6,",",a,",",1pe12.3)') iw,cend(nw),flps(iw,ia)
            endif
          enddo
          !:: check total flux
          if(tfps(nw,ia).le.0.0d0) then
            write(n6,'("Error:total flux on boundary <=0 :",a,":"
     >        ,1pe12.3)') cend(nw),tfps(nw,ia)
            call wexit("check_flps","total flux on boundary <= 0")
          endif
        enddo ! nw
        !::
        !:: flux must negative for nw = 5 (core to SOL)
        nw = 5
        do iw = npsp(nw), npep(nw)-1
          if(flps(iw,ia).ge.0.0d0) then
            ! you may call wexit
            write(n6,'("Waring:flux on core boundary => 0 :iw"
     >         ,i6,",",a,",",1pe12.3)') iw,cend(nw),flps(iw,ia)
          endif
        enddo !iw
        !:: check total flux
        if(tfps(nw,ia).ge.0.0d0) then
          write(n6,'("Waring:total flux on core boundary =>0 :",a,":"
     >        ,1pe12.3)') cend(nw),tfps(nw,ia)
        endif
      !::
      enddo ! ia
      end subroutine check_flps
