!**********************************************************************
      subroutine set_forc
!**********************************************************************
!
!     set thermal force for impmc calculaion
!
!----------------------------------------------------------------------
      use cimcom, only : cfez, cfgte, cfgti, flmxe_at_impmc
     >    , flmxi_at_impmc, gdez, gdte, gdti, glte, glti
     >    , heat_flux_para_at_ic_center, heat_flux_para_e_ic_center
     >    , ndis
      use cntcom, only : mcel
      use cplmet, only : icel, itmax, itmps, itpve, itsls, jcel, jtmax
      use csize,  only : ndmc, ndx
      use cunit,  only : n6
      implicit none
!
!::local variables
! deleted 1 line organize local variables and include files by kamata 2021/06/28
!ik   character  clin*120
      integer  ic, iz, it, jt
      integer  jte, j, i
      real*8   frez(ndx), frte(ndx), frti(ndx), flge(ndx), flgi(ndx)
      real*8   q_para_temp_along_j_center(ndx)
      real*8   q_para_temp_elec_j_center(ndx)
      real*8   gdzero, taue, taui
      real*8   flmxi_temp, flmxe_temp
!
      write(n6,'(2x,"*** set_forc ***")')
!
!----------------------------------------------------------------------
!::zero clear
!----------------------------------------------------------------------
      gdzero = 0.0d0
      gdez(0:ndmc) = 0.0d0
      gdte(0:ndmc) = 0.0d0
      gdti(0:ndmc) = 0.0d0
      glte(0:ndmc) = gdzero
      glti(0:ndmc) = gdzero
!
      cfez (0:ndis) = 0.0d0
      cfgte(0:ndis) = 0.0d0
      cfgti(0:ndis) = 0.0d0
!
      q_para_temp_along_j_center(1:ndx) = 0.d0
      q_para_temp_elec_j_center(1:ndx) = 0.d0

      heat_flux_para_at_ic_center(0:ndmc) = 0.d0
      heat_flux_para_e_ic_center(0:ndmc) = 0.d0
!
!::header
      call slwtime(0,0,0,taue,taui,-1)
!
      do it = 1, itmax
      if( it.ne.itsls .and. it.ne.itmps ) cycle
      write(n6,'(2x,"it =",i3)') it
      do jt = 1, jtmax(it)
      if( jt.le.3 .or. jt.ge.jtmax(it)-2 .or. mod(jt,20).eq.0 ) then
      j = jcel(jt,it)
      i = icel(jt,it)
      iz = 1
      call slwtime(j,i,iz,taue,taui,1)
      endif
      enddo
      enddo
!
!----------------------------------------------------------------------
!::Ez,gradTi,grdTe force
!----------------------------------------------------------------------
      do it = 1, itmax
      if( it.eq.itsls .or. it.eq.itpve ) cycle
      call gdfrc(it,frez,frte,frti,flge,flgi)

! 2018/05/08 receive_q_para from SOLDOR (Y.Homma)
      call receive_q_para(it,q_para_temp_along_j_center,
     >     q_para_temp_elec_j_center)
!
      jte = jtmax(it)-1
      do jt = 2, jte
      i = icel(jt,it)
      j = jcel(jt,it)
      ic  = mcel(j,i)
      if( ic.le.0 ) then
        write(n6,'("jt,it =",2i5,"  j,i =",2i5,"  ic =",i7)')
     >     jt, it, j, i, ic
        call wexit("gdfrc","ic.le.0")
      endif
      gdez(ic) = frez(jt)
      gdte(ic) = frte(jt)
      gdti(ic) = frti(jt)
      glte(ic) = flge(jt)
      glti(ic) = flgi(jt)

      heat_flux_para_at_ic_center(ic)
     &     = q_para_temp_along_j_center(jt)
      heat_flux_para_e_ic_center(ic)
     &     = q_para_temp_elec_j_center(jt)

      call receive_flmxi(j, i, flmxi_temp, flmxe_temp )
      flmxi_at_impmc(ic) = flmxi_temp
      flmxe_at_impmc(ic) = flmxe_temp
      enddo
      enddo
!
      return
      end
!
!**********************************************************************
      subroutine dbg_forc(cmsg)
!**********************************************************************
      use cimcom, only : gdez, gdte, gdti
      use cntcom, only : ncmax
      use cunit,  only : n6
      implicit none
!
!::argument
! modified 1/1 lines organize local variables and include files by kamata 2021/06/28
!ik   character  cmsg*(*)
      character, intent(in) :: cmsg*(*)
!
!::local variables
! modified 1/3 lines organize local variables and include files by kamata 2021/06/28
!ik   integer  ic, imox, imoy
      integer  ic
! function
      integer    imox, imoy
!
      write(n6,'(2x,"*** dbg_forc ***  ",a)') trim(cmsg)
      do ic = 1, ncmax, 200
      write(n6,'(2x,"DBG_FORC",2x,i7,i6,i5,2x,1p3e12.3)')
     >  ic, imox(ic), imoy(ic), gdez(ic), gdte(ic), gdti(ic)
      enddo
!
      return
      end
