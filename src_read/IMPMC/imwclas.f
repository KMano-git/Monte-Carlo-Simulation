!***********************************************************************
      subroutine imwclas(rgpa,rgpc,cnpt,cnwt)
!***********************************************************************
!
!   wght(ip) => rgpa, rgpc, cnpt, cnwt
!
!   rgpa(0:10), rgpc(0:10)  rg:region, pa:allpartcile  pc:cal-partcile
!   cnpt(0:10), cnwt(0:10)  cn:condition, pt:partcile, wt:weight
!
!   cnpt(i)  idemt=1, idcal=2, idexh = 3, idstk = 4
!            iderr=5, idetc=6
!   idexh : exhausted (contains stick)
!   idstk : stick and wght < wrfmin (No trace)
!   idetc : et cetra
!
!
!::refer common variables
!      integer :: npmax, ien(1), ihwl(1), il(1), ir(1), mrgn(1)
!      character(1) :: chwl(1)
!      real(8) :: wght(1)
!      integer :: idcal, idexh, idstk, iderr, idetc
!
!    difference of krfl = 5 and 6
!      krfl = 5  lwstk(ih) = 1                     ien(ip) = 6
!      krfl = 6  lwstk(ih) = 0 and wght < wrfmin   ien(ip) = 7
!
!    ien(ip) = -1(run), 0(tmend), 1(puf),  2(rec),  3(ion),
!              4(wln),  5(wli),   6(stk),  7(abs),  9(err)
!     No trace in case of  ien = 6(stk), 7(abs) 9(err)
!
!     wght(ip) > 0 in those cases  (flux to wall)
!
!     elseif( icnd.eq.6 .or. icnd.eq.7 ) then
!       hcnpt(idstk) = hcnpt(idstk) + 1
!       hcnwt(idstk) = hcnwt(idstk) + zwt
!       hcnwt(idexh) = hcnwt(idexh) + zwt  ! Not wgt0-zwt   No trace
!
!   Note  hcnpt(idexh) = 0  hcnpt(idstk) > 0  last condition
!         hcnwt(idexh) > 0  hcnwt(idstk)
!
!-----------------------------------------------------------------------
      use cimcom, only : idcal, idemt, iderr, idetc, idexh, idstk, ien
     >    , ir, npmax, wght, wpemt, wtemt
      use cntcom, only : mrgn
      implicit none
!
!
!::local variables
      integer :: ip, icnd, ic, irg
      real(8) :: wgt0, zwt
!
!::argument
      real(8), intent(out) :: cnpt(0:10), cnwt(0:10)
     >                      , rgpa(0:10), rgpc(0:10)
!
!
      rgpa(0:10) = 0.0d0   ! region for all particle
      rgpc(0:10) = 0.0d0   ! region for cal-particle  ien = 0
      cnpt(0:10) = 0.0d0   ! condition of all particles
      cnwt(0:10) = 0.0d0   ! condition of all particles
!
!::wght when emitted
      wgt0 = 1.0d0    ! see imemit
!
      do ip = 1, npmax
        icnd = ien(ip)
        zwt  = wght(ip)
        cnwt(idexh) = cnwt(idexh) + wgt0-zwt    ! exhausted wght
        if( icnd.eq.0 ) then
           cnpt(idcal) = cnpt(idcal) + 1
           cnwt(idcal) = cnwt(idcal) + zwt
        elseif( icnd.eq.6 .or. icnd.eq.7 ) then
           cnpt(idstk) = cnpt(idstk) + 1
           cnwt(idstk) = cnwt(idstk) + zwt
           cnwt(idexh) = cnwt(idexh) + zwt      ! No trace
        elseif( icnd.eq.9 ) then
           cnpt(iderr) = cnpt(iderr) + 1
           cnwt(iderr) = cnwt(iderr) + zwt
        else
           cnpt(idetc) = cnpt(idetc) + 1
           cnwt(idetc) = cnwt(idetc) + zwt
        endif
        ic  = ir(ip)
        irg = mrgn(ic)
        rgpa(irg) = rgpa(irg) + zwt
        if( icnd.eq.0 ) rgpc(irg) = rgpc(irg) + zwt
      enddo
!
      cnpt(idemt) = wpemt
      cnwt(idemt) = wtemt
!
      return
      end
