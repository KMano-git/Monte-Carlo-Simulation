!***********************************************************************
      subroutine plauxg
!***********************************************************************
!
!     plasma parameter at grid point
!
!    (j-1/2,i+1/2)   (j+1/2,i+1/2)       (j-1,i)         (j,i)
!        *----------------*                *----------------*  
!        |                |                |                |  
!        |      (j,i)     |                |      (j,i)     |
!        |                |                |                |  
!        *----------------*                *----------------*  
!    (j-1/2,i-1/2)   (j+1/2,i-1/2)       (j-1,i-1)       (j,i-1)
!
!
!      F^sita(j+1/2,i)                 F^psi(j,i+1/2)
!                                                   N
!                       Qg(j,i)          Qg(jm,i)   /|   Qg(j,i) 
!        *----------------N                N--------||------E  
!        |                |                |        ||      |  
!        |    E (j,i)   ====>  W           |     S(j,i)     |
!        |                |                |                |  
!        *----------------S                *----------------*  
!                       Qg(j,i-1)                                 
!                                        Note  jm#j-1
!
!     Qg(j,i) = Q(j+1/2,i+1/2)
!             = 1/4*( Q(j,i)+Q(jp,i)+Q(j,ip)+Q(jp,ip) )
!
!    Note.  No need special treatment at the X-point
!
!      different interpolation at X-point
!            Qg(jcxp1,ispcx) # Qg(jcxp2-1,ispx)
!
!      No effect on cal. result even if  Qg(Xp) = 1.0d40
!            because of orthogonal mesh at X-point
!
!     j  = 1   2   3 ...      J-1   J
!     jg =   1   2  ...         J-1
!               im main    same point  j=1 and J-1  
!
!CHK_ia  plauxg [do ia] = 1, nion  Done
!----------------------------------------------------------------------
      use cplcom, only : c12, c14, nion, vna, vnag, vte, vteg, vti, vtig
     >    , vva, vvag
      use cplmet, only : icel, icspx, itmax, itmps, itpve, jcdp1, jcdp2
     >    , jcel, jcxp1, jcxp2, jtmax
      use csonic, only : itim
      implicit none
!
      integer  it, jmax, jwst, jwen, jw, j, jp, i, ip, ia, j1, j2
      integer  ic, jc1, jc2
      real*8   zna, zva
!
!x      write(n6,'(2x,"*** plauxg ***")')
!
      do it = 1, itmax
      jmax  = jtmax(it)
      jwst  = 1
      jwen  = jmax-1
      if( it.ge.itmps ) jwst = 2
!
      do jw = jwst, jwen
      j  = jcel(jw,it)
      jp = jcel(jw+1,it)
      i  = icel(jw,it)
      ip = i + 1
!
      do ia = 1, nion
      zna = c14*(vna(j,i,ia)+vna(j,ip,ia)+vna(jp,i,ia)+vna(jp,ip,ia))
      zva = c14*(vva(j,i,ia)+vva(j,ip,ia)+vva(jp,i,ia)+vva(jp,ip,ia))
      vnag(j,i,ia) = zna
      vvag(j,i,ia) = zva
      enddo
!
      vtig(j,i) = c14*(vti(j,i)+vti(j,ip)+vti(jp,i)+vti(jp,ip))
      vteg(j,i) = c14*(vte(j,i)+vte(j,ip)+vte(jp,i)+vte(jp,ip))
      enddo
!
      enddo
!
!::d-plate (j=3/2, j=jmax-1/2)
      j1 = jcdp1
      j2 = jcdp2
! modified 2/1 lines organize local variables and include files by kamata 2021/05/31
!ik   do it = 1, itpve-1
!ik   i  = it
      do i = 1, itpve-1
      ip = i + 1
      do ia = 1, nion
      vnag(j1,  i,ia) = c12*(vna(j1,i,ia)+vna(j1,ip,ia))
      vvag(j1,  i,ia) = c12*(vva(j1,i,ia)+vva(j1,ip,ia))
      vnag(j2-1,i,ia) = c12*(vna(j2,i,ia)+vna(j2,ip,ia))
      vvag(j2-1,i,ia) = c12*(vva(j2,i,ia)+vva(j2,ip,ia))
      enddo
      vtig(j1,  i) = c12*(vti(j1,i)+vti(j1,ip))
      vteg(j1,  i) = c12*(vte(j1,i)+vte(j1,ip))
      vtig(j2-1,i) = c12*(vti(j2,i)+vti(j2,ip))
      vteg(j2-1,i) = c12*(vte(j2,i)+vte(j2,ip))
      enddo
!
!::X-point  Qg(Xp) = 0.5*(Qg(Xp_o)+Qg(Xp_i))
!x      ic  = icspx
!x      jc1 = jcxp1
!x      jc2 = jcxp2-1
!x      do ia = 1, nion
!x      zna = c12*(vnag(jc1,ic,ia)+vnag(jc2,ic,ia))
!x      zva = c12*(vvag(jc1,ic,ia)+vvag(jc2,ic,ia))
!x      vnag(jc1,ic,ia) = zna
!x      vnag(jc2,ic,ia) = zna
!x      vvag(jc1,ic,ia) = zva
!x      vvag(jc2,ic,ia) = zva
!x      enddo
!x      zti = c12*(vtig(jc1,ic)+vtig(jc2,ic))
!x      zte = c12*(vteg(jc1,ic)+vteg(jc2,ic))
!x      vtig(jc1,ic) = zti
!x      vtig(jc2,ic) = zti
!x      vteg(jc1,ic) = zte
!x      vteg(jc2,ic) = zte
!
!::no depend X-point value
      if( itim.eq.0 ) then
!x    write(n6,'(2x,"set value to 1.0d40 at X-point")')
      endif
!
      ic  = icspx
      jc1 = jcxp1
      jc2 = jcxp2-1
      do ia = 1, nion
      vnag(jc1,ic,ia) = 1.0d40
      vnag(jc2,ic,ia) = 1.0d40
      vvag(jc1,ic,ia) = 1.0d40
      vvag(jc2,ic,ia) = 1.0d40
      enddo
      vtig(jc1,ic) = 1.0d40
      vtig(jc2,ic) = 1.0d40
      vteg(jc1,ic) = 1.0d40
      vteg(jc2,ic) = 1.0d40
!
      return
      end
