!***********************************************************************
      subroutine plout_reg
!***********************************************************************
      use cntmnt, only : totsn, totsp, totwe, totwi
      use cplcom, only : ama, nion, q1a, q2a, q3, q4
      use cplmet, only : hvol, icel, itmax, itpve, jcel, jtmax, kreg
      use cplqcn, only : qvl_dt
      use csize,  only : ndsp
      use csonic, only : dtim, itim, time
      use cunit,  only : cdrout, lmspe, lmype
      implicit none
!
!::local variables
      real*8  tdtna(10,ndsp),tdtpa(10,ndsp),tdtei(10),tdtee(10)
      real*8  tsmna(10,ndsp),tsmpa(10,ndsp),tsmei(10),tsmee(10)
      integer ir, mj, nsza, nsiz, it, lenx, jmax, jt, j, i, ia
      integer m1a, m2a, m3, m4
!
      integer noft(10), ift, ism
      character creg(10)*3
      data  ift/0/
      data  (noft(i),i=1,7)/ 121,  122,  123,  124,  125,  126,  127/
      data  (creg(i),i=1,7)/"odp","sol","idp","opv","ipv","edg","man"/
      parameter (ism=8)
      data  noft(ism)/128/, creg(ism)/"all"/
      save  noft, ift, creg
!
!xx   return
!
      if( lmype.ne.lmspe ) return
!
!::file allocation
      if( ift.eq.0 ) then
      do ir = 1, ism
      ift = noft(ir)
      mj = lenx(cdrout)
      open( unit=ift, file=cdrout(1:mj)//"zreg_"
     >           //creg(ir)(1:lenx(creg(ir))) )
      write(ift,'(3x,"itim",4x,"time",8x,"dtim",8x,"tmNa",8x,"tmPa",8x,
     >  "tmEi",8x,"tmEe",8x,"dtNa",8x,"dtPa",8x,"dtEi",8x,"dtEe",8x,
     >  "Sn",10x,"Sp",10x,"Wi",10x,"We")')
      enddo
      endif
!
!::clear
      nsza = 10*ndsp
      nsiz = 10
      call setd( tdtna, nsza, 0.0d0 )
      call setd( tdtpa, nsza, 0.0d0 )
      call setd( tdtei, nsiz, 0.0d0 )
      call setd( tdtee, nsiz, 0.0d0 )
      call setd( tsmna, nsza, 0.0d0 )
      call setd( tsmpa, nsza, 0.0d0 )
      call setd( tsmei, nsiz, 0.0d0 )
      call setd( tsmee, nsiz, 0.0d0 )
!
!::sumup
      m3 = 2*nion + 1
      m4 = 2*nion + 2
      do it = 2, itmax-1
      if( it.eq.itpve ) cycle
      jmax = jtmax(it)
      do jt = 2, jmax-1
      j  = jcel(jt,it)
      i  = icel(jt,it)
      ir = kreg(j,i)
      tdtei(ir) = tdtei(ir) + qvl_dt(j,i,m3)
      tdtee(ir) = tdtee(ir) + qvl_dt(j,i,m4)
      tsmei(ir) = tsmei(ir) + q3(j,i)*hvol(j,i)
      tsmee(ir) = tsmee(ir) + q4(j,i)*hvol(j,i)
      do ia = 1, nion
      m1a = 2*ia - 1
      m2a = 2*ia 
      tdtna(ir,ia) = tdtna(ir,ia) + qvl_dt(j,i,m1a)
      tdtpa(ir,ia) = tdtpa(ir,ia) + qvl_dt(j,i,m2a)
      tsmna(ir,ia) = tsmna(ir,ia) + q1a(j,i,ia)*hvol(j,i)
      tsmpa(ir,ia) = tsmpa(ir,ia) + q2a(j,i,ia)*hvol(j,i)
      enddo
      enddo
      enddo
!
!::particle number
      do ia = 1, nion
      do ir = 1, 6
      tdtna(ir,ia) = tdtna(ir,ia)/ama(ia)
      tsmna(ir,ia) = tsmna(ir,ia)/ama(ia)
      enddo
      enddo
!
!::hot core  (avoid Zero divide)
      ir = 7
      do ia = 1, nion
      tdtna(ir,ia) = 1.0
      tdtpa(ir,ia) = 1.0
      enddo
      tdtei(ir) = 1.0
      tdtee(ir) = 1.0
!
!::total
      do ir = 1, 6
      tdtei(ism) = tdtei(ism) + tdtei(ir)
      tdtee(ism) = tdtee(ism) + tdtee(ir)
      tsmei(ism) = tsmei(ism) + tsmei(ir)
      tsmee(ism) = tsmee(ism) + tsmee(ir)
      do ia = 1, nion
      tdtna(ism,ia) = tdtna(ism,ia) + tdtna(ir,ia)
      tdtpa(ism,ia) = tdtpa(ism,ia) + tdtpa(ir,ia)
      tsmna(ism,ia) = tsmna(ism,ia) + tsmna(ir,ia)
      tsmpa(ism,ia) = tsmpa(ism,ia) + tsmpa(ir,ia)
      enddo
      enddo 
!
!::total source in system (except for core)
      do ir = 1, 6
      totwi(ism) = totwi(ism) + totwi(ir)
      totwe(ism) = totwe(ism) + totwe(ir)
      do ia = 1, nion
      totsn(ism,ia) = totsn(ism,ia) + totsn(ir,ia)
      totsp(ism,ia) = totsp(ism,ia) + totsp(ir,ia)
      enddo
      enddo
!
!::output
      ia = 1
      do ir = 1, ism
      ift = noft(ir)
      write(ift,'(i7,1p14e12.3)')
     >  itim,time,dtim,
     >  tsmna(ir,ia)/tdtna(ir,ia), tsmpa(ir,ia)/tdtpa(ir,ia),
     >  tsmei(ir)/tdtei(ir), tsmee(ir)/tdtee(ir),
     >  tdtna(ir,ia),tdtpa(ir,ia),tdtei(ir),tdtee(ir),
     >  totsn(ir,ia),totsp(ir,ia),totwi(ir),totwe(ir)
      enddo
!
      return
      end
