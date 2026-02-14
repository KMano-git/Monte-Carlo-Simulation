!**********************************************************************
      subroutine pledden
!**********************************************************************
      use cplcom, only : fcna, ibcyl, mdl_edt, nion, vna
      use cplmet, only : icel, icmpe, icspx, icwl1, icwl2, itmax, itsle
     >    , jcel, jtmax, jtmin, kreg
      use cpmpls, only : jmd1
      use cunit,  only : n6
      implicit none
!
!::local variables
      integer :: it, jwst, jwen, ipbc, jw, j, i, ia, irg
      integer :: ist, ien
      real(8) :: dnae, dnas, dna
!
!::input data
      real(8) :: dnied = 1.0d19
      real(8) :: dnisl = 0.5d19
!
      write(n6,'(/2x,"*** pledden ***  mdl_edt(3) = ",i2)') mdl_edt
      write(n6,'(2x,"dnied =",1pe12.3,"  dnisl =",1pe12.3)')
     >    dnied, dnisl
!
!::debug write
      ia  = 1
!
      it  = itsle
      j   = jtmin(it)
      ist = icwl1
      ien = icwl2
      write(n6,'(5x,"vna(j,i,ia)  j = ",i4,"  icwl1,icwl2 =",2i4)')
     >   j, icwl1, icwl2
      write(n6,'(5x,"vna =",1p10e11.3)') (vna(j,i,ia),i=ist,ien)
!
      it  = itsle
      j   = jtmax(it)
      ist = icwl1
      ien = icwl2
      write(n6,'(5x,"vna(j,i,ia)  j = ",i4,"  icwl1,icwl2 =",2i4)')
     >   j, icwl1, icwl2
      write(n6,'(5x,"vna =",1p10e11.3)') (vna(j,i,ia),i=ist,ien)
!
      j   = jmd1
      ist = icwl1
      ien = icmpe
      write(n6,'(5x,"vna(j,i,ia)  j = ",i4,"  icwl1,icmpe =",2i4)')
     >   j, icwl1, icmpe
      write(n6,'(5x,"vna =",1p10e11.3)') (vna(j,i,ia),i=ist,ien)
!
!::new density
      do it = 1, itmax
      jwst = jtmin(it)
      jwen = jtmax(it)
      ipbc = ibcyl(it)
      if( ipbc.eq.1 ) then
      jwst = jwst+1
      jwen = jwen-1
      endif
!
      do jw = jwst, jwen
      j = jcel(jw,it)
      i = icel(jw,it)
      do ia = 1, nion
      dnae = dnied*fcna(ia)  ! extra density at icmpe
      dnas = dnisl*fcna(ia)  ! extra density at icspx
      dna  = dnas + (dnae-dnas)/dfloat(icmpe-icspx)*(dfloat(i)-icspx)
      irg  = kreg(j,i)
      if( irg.ne.6 ) dna = dnas
      vna(j,i,ia) = vna(j,i,ia) + dna
      enddo
      enddo
!
      enddo
!
!::debug write
      write(n6,'(5x,"new density")')
      ia  = 1
!
      it  = itsle
      j   = jtmin(it)
      ist = icwl1
      ien = icwl2
      write(n6,'(5x,"vna(j,i,ia)  j = ",i4,"  icwl1,icwl2 =",2i4)')
     >   j, icwl1, icwl2
      write(n6,'(5x,"vna =",1p10e11.3)') (vna(j,i,ia),i=ist,ien)
!
      it  = itsle
      j   = jtmax(it)
      ist = icwl1
      ien = icwl2
      write(n6,'(5x,"vna(j,i,ia)  j = ",i4,"  icwl1,icwl2 =",2i4)')
     >   j, icwl1, icwl2
      write(n6,'(5x,"vna =",1p10e11.3)') (vna(j,i,ia),i=ist,ien)
!
      j   = jmd1
      ist = icwl1
      ien = icmpe
      write(n6,'(5x,"vna(j,i,ia)  j = ",i4,"  icwl1,icmpe =",2i4)')
     >   j, icwl1, icmpe
      write(n6,'(5x,"vna =",1p10e11.3)') (vna(j,i,ia),i=ist,ien)
!
      return
      end
