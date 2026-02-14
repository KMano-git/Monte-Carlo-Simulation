!***********************************************************************
      subroutine pst_rprf
!***********************************************************************
      use cgdcom, only : grdx, grdy, mqs1, npmx
      use cplmet, only : icaxs, icmpe, icmps, icspx, icwl1, jcmax
     >   , icmax, icmin 
      use cplpst, only : rhf_idp, rhf_imd, rhf_odp, rhf_omd
      use cpmpls, only : imd1, imd2, jmd1, jmd2, r0mp, ramp, xmd1, xmd2
     >    , ymd1, ymd2
      use csize,  only : ndy
      use cunit,  only : n6
      implicit none
!
!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer ix, iy, i6, n6sav
      integer ix, iy
      real*8  xaxs, yaxs
      integer jc, n, k
! modified 2/1 lines organize local variables and include files by kamata 2021/05/31
!ik   real*8  rst(ndy), ren(ndy), rhf(ndy)
!ik   data i6/6/; save i6
      real*8  rst(ndy), ren(ndy)
      logical :: is_midpl, is_diverter
!
!::axis
      ix = mqs1
      iy = npmx
      xaxs = grdx(ix,iy)
      yaxs = grdy(ix,iy)
!
      write(n6,'(2x,"axis point =",2f8.4,"  out-mid =",2f8.4,2i5,
     >  "  in-mid =",2f8.4,2i5)') xaxs,yaxs,xmd1,ymd1,jmd1,imd1
     >  ,xmd2,ymd2,jmd2,imd2
      write(n6,'(2x,"R0 =",f8.4,"  a =",f8.4)') r0mp,ramp
      write(n6,'(2x,"icwl1 =",i3,"  icspx =",i3,"  icmps =",i3
     > ,"  icmpe =",i3,"  icaxs =",i3)') icwl1,icspx,icmps,icmpe,icaxs
!
!::radial position
      jc = 1
      call codrad(jc,rst,ren,rhf_odp,n,ndy,icmin(jc),icmax(jc))
      jc = jcmax
      call codrad(jc,rst,ren,rhf_idp,n,ndy,icmin(jc),icmax(jc))
      jc = jmd1
      call codrad(jc,rst,ren,rhf_omd,n,ndy,icmin(jc),icmax(jc))
      jc = jmd2
      call codrad(jc,rst,ren,rhf_imd,n,ndy,icmin(jc),icmax(jc))
!
!::radial profile
      do k = 1, 6
        is_midpl = .false.
        is_diverter = .false.
        if( k.eq.1 ) then
          jc = 1
          is_diverter = .true.
        elseif( k.eq.2 ) then
          jc = 2
          is_diverter = .true.
        elseif( k.eq.3 ) then
          jc = jmd1
          is_midpl = .true.
        elseif( k.eq.4 ) then
          jc = jmd2
          is_midpl = .true.
        elseif( k.eq.5 ) then
          jc = jcmax-1
          is_diverter = .true.
        elseif( k.eq.6 ) then
          jc = jcmax
          is_diverter = .true.
        endif
        call lst_rad(jc,is_midpl,is_diverter)
      enddo
!
      return
      end
!
!***********************************************************************
      subroutine lst_rad(jc,is_midpl,is_diverter)
!***********************************************************************
      use cntmnt, only : sn0
      use cplcom, only : vdda, vdxe, vdxi, vne, vni, vte, vti, vva
     > ,tN0, tT0, tNg, tTg
      use cplmet, only : icmax, icmin, jcmax
      use cplpst, only : rhf_idp, rhf_imd, rhf_odp, rhf_omd
      use cntcom, only : mcel
      implicit none
!
      integer, intent(in) :: jc
      logical, intent(in) :: is_midpl, is_diverter
!
      integer ia, ic1, ic2, i6, ic, ic_MC
      character  cdsn*12
      real(8) :: p0, pg, odp_tmp, idp_tmp, omd_tmp, imd_tmp
!
      if( jc.le.0 .or. jc.gt.jcmax ) return
!
      i6 = 21
!
      write(cdsn,'(a,i3.3,a)') "PRF_j", jc, ".txt"
      open(unit=i6, file=cdsn)
!
!::plasma parameter
      ia = 1
      ic1 = icmin(jc)
      ic2 = icmax(jc)
      write(i6,'(2x,"*** radial profile ***  jc =",i3)') jc
      write(i6,'(3x,"ic",4x,"r_odp",6x,"r_idp",6x,"r_omd",6x,"r_imd",
     >   9x,"Ne",9x,"Ni",9x,"Te",9x,"Ti",9x,"Va",9x,"N0",
     >   9x,"E0",9x,"P0",9x,"Ng",9x,"Eg",9x,"Pg",
     >   9x,"Da",9x,"xi",9x,"xe")')
!
      do ic = ic1, ic2
        ic_MC = mcel(jc,ic)
        p0 = tN0(ic_MC,ia)*tT0(ic_MC,ia)*1.60210d-19
        pg = tNg(ic_MC,ia)*tTg(ic_MC,ia)*1.60210d-19
!!
        if(is_midpl .and. rhf_odp(ic) < 0) then
          odp_tmp = 0.0_8
        else
          odp_tmp = rhf_odp(ic)*100.0d0
        endif
!!
        if(is_midpl .and. rhf_idp(ic) > 0) then
          idp_tmp = 0.0_8
        else
          idp_tmp = rhf_idp(ic)*100.0d0
        endif
!!
        if(is_diverter .and. rhf_omd(ic) < 0) then
          omd_tmp = 0.0_8
        else
          omd_tmp = rhf_omd(ic)*100.0d0
        endif
!!
        if(is_diverter .and. rhf_imd(ic) > 0) then
          imd_tmp = 0.0_8
        else
          imd_tmp = rhf_imd(ic)*100.0d0
        endif
!!
        write(i6,'(1x,i5,1p18e11.3,i8)') ic,
     >    odp_tmp, idp_tmp,
     >    omd_tmp, imd_tmp,
     >    vne(jc,ic), vni(jc,ic), vte(jc,ic), vti(jc,ic),
     >    vva(jc,ic,ia), sn0(jc,ic,ia),
     >    tT0(ic_MC,1), p0, tNg(ic_MC,1), tTg(ic_MC,1), pg,
     >    vdda(jc,ic,ia), vdxi(jc,ic), vdxe(jc,ic)
      enddo
      close(i6)
!
      return
      end
