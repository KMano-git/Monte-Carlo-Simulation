!
!**********************************************************************
      subroutine receive_q_para(it,q_para_temp_along_j_center,
     >     q_para_temp_elec_j_center)
!**********************************************************************
!        2018/05/08 Yuki HOMMA:
!                              Receive conductive heat flux density q//
!             from SOLDOR, to pass it to the array
!             heat_flux_para_by_kappa_para(ic over all IMPMC mesh) in set_forc.f.
!
!             2018/10/3  q// ion & elec are both received in this routine,
!             to avoid "include cplcom(SOLDOR variables)" at "set_forc.f" routine
!             where this "receive_q_para" routine is called.
!
!        it   :  tube number
!        q_para_temp_along_j :  heat flux density [ joule/(m2*s) ]
!                             given from SOLDOR.
!                           q// containing only //diffusive heat transfer
!                           at cell boundary  j--jp, i.e. j+1/2.
!                           Contribution from q_perp at mesh boundary j+1/2
!                           due to mesh slant is excluded (cf. soldor/pxviscN.f).
!
!----------------------------------------------------------------------
      use cplcom, only : heat_flux_para_by_kappa_para
     >    , heat_flux_para_elec
      use cplmet, only : hdxm, hdxp, icel, itmax, jcel, jtmax
      use csize,  only : ndx
      implicit none
!
!::argument
! modified 1/1 lines organize local variables and include files by kamata 2021/06/28
!ik   integer   it
      integer, intent(in)  :: it
      real(8), intent(out) :: q_para_temp_along_j_center(ndx)
      real(8), intent(out) :: q_para_temp_elec_j_center(ndx)
!
!::local common
! deleted 2 lines replace all include files with module files by kamata 2021/08/18
!ik   real*8  dlnx(ndx), zpes(ndx), ztes(ndx), ztis(ndx)
!ik   common /com_gdfrc/ dlnx, zpes, ztes, ztis
!
!::local
!      real*8    gpe, gte, gti, zne, fjt, rmfpe, rmfpi
! modified 1/1 lines organize local variables and include files by kamata 2021/06/28
!ik   integer  jte, jt, j, i, jp, jm, js
      integer  jte, jt, j, i, jm
      real*8  q_para_temp_along_j_east_side(ndx)
      real*8  q_para_temp_elec_j_east_side(ndx)
      real*8 weight_east, weight_west
!
!::length in cell along the magnetic field line
!      call gdlen(it,dlnx)
!
      if( it.le.0 .or. it.gt.itmax ) return
!
!::Zero set
      q_para_temp_along_j_center(1:ndx) = 0.d0
      q_para_temp_elec_j_center(1:ndx) = 0.d0

      q_para_temp_along_j_east_side(1:ndx) = 0.0d0
      q_para_temp_elec_j_east_side(1:ndx) = 0.0d0
!
!
!     Loop set
      jte = jtmax(it)
!
!::receieve q// from SOLDOR
      do jt = 1, jte
         j  = jcel(jt,it)
         i  = icel(jt,it)
!
!     substitute from SOLDOR variablaes
         q_para_temp_along_j_east_side(jt) =
     &        heat_flux_para_by_kappa_para(j,i)
         q_para_temp_elec_j_east_side(jt) =
     &        heat_flux_para_elec(j,i)
!      heat flux at cell boundary j+1/2, i.e. East side
!      Unit [joule/(m2 * s)]
      enddo                     ! do jt = 1, jte

!::heat flux at cell center   interpolation
      ! jt = 1 extreme
      jt = 1
      q_para_temp_along_j_center(jt)
     &     = q_para_temp_along_j_east_side(jt)
      q_para_temp_elec_j_center(jt)
     &     = q_para_temp_elec_j_east_side(jt)

      ! jt = jte extreme
      jt = jte
      q_para_temp_along_j_center(jt)
     &     = q_para_temp_along_j_east_side(jt-1)
      q_para_temp_elec_j_center(jt)
     &     = q_para_temp_elec_j_east_side(jt-1)

      ! jt = 2 <--> (jte-1)
      do jt = 2, jte-1

         j  = jcel(jt,it)
         i  = icel(jt,it)
         jm = jcel(jt-1, it)

         weight_west = hdxm(j,i) / ( hdxp(jm,i)+hdxm(j,i) )
         weight_east = hdxp(jm,i) / ( hdxp(jm,i)+hdxm(j,i) )
!
         q_para_temp_along_j_center(jt)
     &        = (weight_west*q_para_temp_along_j_east_side(jt-1) )
     &        + (weight_east*q_para_temp_along_j_east_side(jt) )

         q_para_temp_elec_j_center(jt)
     &        = (weight_west*q_para_temp_elec_j_east_side(jt-1) )
     &        + (weight_east*q_para_temp_elec_j_east_side(jt) )
      enddo                     ! do jt = 2, jte-1
!
!     if( jt.eq.1 .or. jt.eq.jte ) then
!      else
!      enddo
!
      return
      end subroutine
!!!
!!!
!!!
!
!**********************************************************************
      subroutine receive_flmxi(j, i, flmxi_received, flmxe_received)
!**********************************************************************
!        2018/06/23 Yuki HOMMA:
!                              receive flmxi&flmxe from SOLDOR,
!       to pass it to the array
!       flmxi&flmxe_at_impmc(ic over all IMPMC mesh) in set_forc.f.
!       Since this flmxi is used for test with switch "mdl_fthi_wo_limiter",
!       interpolation to cell center is not considered.
!
!       Once the thermal force development phase is over,
!       flmxi & flmxe may not be necessary anymore, and can be removed.
!----------------------------------------------------------------------
      use cplcom, only : flmxe, flmxi
      implicit none
!
!::argument
! modified 1/1 lines organize local variables and include files by kamata 2021/06/28
!ik   integer   j, i
      integer, intent(in)  :: j, i
      real(8), intent(out) :: flmxi_received, flmxe_received

!::Zero set
      flmxi_received = 0.0d0
      flmxe_received = 0.0d0
      flmxi_received = flmxi(j,i)
      flmxe_received = flmxe(j,i)
      return
      end subroutine
