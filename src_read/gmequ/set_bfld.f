!***********************************************************************
      subroutine set_bfld
!***********************************************************************
      use com_eqdat, only : nr, nz, ubx, uby
      use cunit,     only : n6
      use mod_shexe, only : use_namelist, reqil, refit
      implicit none
!
!::local variables
      integer  nft, nmax
      character  dsn*80,cdsn*80
      logical lex,leq
      real*8  bdir
! The lack of the initial settings for logical parameters may lead to
! the mismatch of the Bfield parameters
      leq = .false.
!
!::equil data
      nft = 21
!::EQDISK
      if(use_namelist) then
!..   use variable in inppls
         cdsn = REQIL
      else
!..   use environment variable
         call getenv( "REQIL", cdsn )
      end if
      if( len_trim(cdsn) .gt. 0 ) then
        leq = .true.
        if(use_namelist) then
!..   use variable in inppls
           call nopen(nft,REQIL,"text read",dsn,lex)
        else
!..   use environment variable
           call nopen(nft,"REQIL","text read",dsn,lex)
        end if
        call eqdisk(nft,"read")
      endif
!::EFIT
      if( .not.leq ) then
         if(use_namelist) then
!..   use variable in inppls
            cdsn = REFIT
         else
!..   use environment variable
            call getenv( "REFIT", cdsn )
         end if
        if( len_trim(cdsn) .gt. 0 ) then
          leq = .true.
          if(use_namelist) then
!..   use variable in inppls
             call nopen(nft,REFIT,"text read",dsn,lex)
          else
!..   use environment variable
             call nopen(nft,"REFIT","text read",dsn,lex)
          end if
          call eqdsk_efit(nft,"read")
        endif
      endif
!::no equil data
      if(.not.leq) call wexit("set_bfld","no definition of equil data")
!
      close(nft)
!
!::unit b-direction
      call bfield
!
!::direction
      call bfdir(bdir)
      write(n6,'(2x,"bdir =",f10.5)') bdir
!
!::correction
      if( bdir.lt.0.0d0 ) then
      write(n6,'(2x,"change b-direction of data (ubx,uby,ubz)")')
!
      nmax = nr*nz   ! <== 2012/01/19
      ubx(1:nmax) = -ubx(1:nmax)
      uby(1:nmax) = -uby(1:nmax)
      endif
!
!::plot
!xx      call plt_bfld
!
      return
      end
