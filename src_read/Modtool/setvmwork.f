! added replace all include files with module files by kamata 2021/08/18
! data copy of vmwork and cntwcn ( cntwedf ) ( monte/inc/cntwcn /cntvar/ )
      subroutine setvmwork( kk, isrc )
      use cntmnt,  only : vmwork
      use cntwcn,  only : wabs, wden, weema, weemm, wees, wehta, wehtm
     >    , wemt, wend, weng, werr, wfhta, wfhtm, wflx, wgden, wgeng
     >    , whta, whtm, wion, wnfl0x, wnfl0y, wnfl0z, wnflgx, wnflgy
     >    , wnflgz, wnrm, wpmp, wreg, wsbr, wssn, wssp, wsum, wswe, wswi
     >    , wtion, wtot, wvlp, wwal
      use cntwedf, only : nene_mx, nsmx, nwmx, wedf, nene_mx_line
     > , wedf_line
      use csize,   only : ndgs, ndmc, ndwp 
      implicit none
! arguments
      integer, intent(in) :: isrc, kk
! isrc : second index number of vmwork
! kk   : flag to set data, = 1 : vmwork to 'cntwcn', = 2 : 'cntwcn' to vmwork
!                          = 3 : initialization of 'cntwcn'

! local variables
      integer    i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, ia, ib, ic, id
     >         , ie, ig, ih, ii, ij, ik, il, im, in, io, ip, iq, ir, is
     >         , it, iu, iv, iw, ix, iy, iz, j0, j1, j2, j3, j4, j5

! set start posiotion of each variable
      i0 = 1                          ! wflx
      i1 = i0 + 1                     ! wssn
      i2 = i1 + ( ndmc + 1 ) * ndgs   ! wssp
      i3 = i2 + ( ndmc + 1 ) * ndgs   ! wswe
      i4 = i3 + ( ndmc + 1 )          ! wswi
      i5 = i4 + ( ndmc + 1 )          ! wsbr
      i6 = i5 + ( ndmc + 1 ) * ndgs   ! wden
      i7 = i6 + ( ndmc + 1 ) * ndgs   ! weng
      i8 = i7 + ( ndmc + 1 ) * ndgs   ! wvlp

      i9 = i8 + ( ndmc + 1 ) * ndgs   ! wnfl0x
      ia = i9 + ( ndmc + 1 ) * ndgs   ! wnfl0y
      ib = ia + ( ndmc + 1 ) * ndgs   ! wnfl0z

      ic = ib + ( ndmc + 1 ) * ndgs   ! wtion
      id = ic + ( ndmc + 1 ) * ndgs   ! wgden
      ie = id + ( ndmc + 1 ) * 2      ! wgeng

      ig = ie + ( ndmc + 1 ) * 2      ! wnflgx
      ih = ig + ( ndmc + 1 ) * 2      ! wnflgy
      ii = ih + ( ndmc + 1 ) * 2      ! wnflgz

      ij = ii + ( ndmc + 1 ) * 2      ! wemt
      ik = ij + ndwp                  ! wwal
      il = ik + ndwp                  ! wend
      im = il + 10                    ! wtot
      in = im + 1                     ! whta
      io = in + 30 * 3                ! whtm
      ip = io + 30 * 3                ! wfhta
      iq = ip + ndwp                  ! wfhtm
      ir = iq + ndwp                  ! wehta
      is = ir + ndwp                  ! wehtm
      it = is + ndwp                  ! wreg
      iu = it + 11                    ! wsum
      iv = iu + 1                     ! wnrm
      iw = iv + 1                     ! wion
      ix = iw + 1                     ! wabs
      iy = ix + 1                     ! wpmp
      iz = iy + 1                     ! werr
      j0 = iz + 1                     ! wedf
      j1 = j0 + nene_mx * nwmx * nsmx ! weema
      j2 = j1 + ndwp                  ! weemm
      j3 = j2 + ndwp                  ! wees
      j4 = j3 + ndwp                  ! wedf_line
      j5 = j4 + nene_mx_line * nwmx * nsmx ! (last)+1
!
      select case ( kk )
      case( 1 )
! vmwork to 'cntwcn'
        wflx = vmwork(i0,isrc)
        wtot = vmwork(im,isrc)
        wsum = vmwork(iu,isrc)
        wnrm = vmwork(iv,isrc)
        wion = vmwork(iw,isrc)
        wabs = vmwork(ix,isrc)
        wpmp = vmwork(iy,isrc)
        werr = vmwork(iz,isrc)

        wswe(0:ndmc)  = vmwork(i3:i4-1,isrc)
        wswi(0:ndmc)  = vmwork(i4:i5-1,isrc)

        wemt(1:ndwp)  = vmwork(ij:ik-1,isrc)
        wwal(1:ndwp)  = vmwork(ik:il-1,isrc)
        wfhta(1:ndwp) = vmwork(ip:iq-1,isrc)
        wfhtm(1:ndwp) = vmwork(iq:ir-1,isrc)
        wehta(1:ndwp) = vmwork(ir:is-1,isrc)
        wehtm(1:ndwp) = vmwork(is:it-1,isrc)
        weema(1:ndwp) = vmwork(j1:j2-1,isrc)
        weemm(1:ndwp) = vmwork(j2:j3-1,isrc)
        wees(1:ndwp)  = vmwork(j3:j4-1,isrc)

        wend(1:10)    = vmwork(il:im-1,isrc)

        wreg(0:10)    = vmwork(it:iu-1,isrc)

        call datcp( vmwork(i1,isrc), wssn,   i2-i1 )
        call datcp( vmwork(i2,isrc), wssp,   i3-i2 )
        call datcp( vmwork(i5,isrc), wsbr,   i6-i5 )
        call datcp( vmwork(i6,isrc), wden,   i7-i6 )
        call datcp( vmwork(i7,isrc), weng,   i8-i7 )
        call datcp( vmwork(i8,isrc), wvlp,   i9-i8 )
        call datcp( vmwork(i9,isrc), wnfl0x, ia-i9 )
        call datcp( vmwork(ia,isrc), wnfl0y, ib-ia )
        call datcp( vmwork(ib,isrc), wnfl0z, ic-ib )
        call datcp( vmwork(ic,isrc), wtion,  id-ic )

        call datcp( vmwork(id,isrc), wgden,  ie-id )
        call datcp( vmwork(ie,isrc), wgeng,  ig-ie )
        call datcp( vmwork(ig,isrc), wnflgx, ih-ig )
        call datcp( vmwork(ih,isrc), wnflgy, ii-ih )
        call datcp( vmwork(ii,isrc), wnflgz, ij-ii )

        call datcp( vmwork(in,isrc), whta,   io-in )
        call datcp( vmwork(io,isrc), whtm,   ip-io )

        call datcp( vmwork(j0,isrc), wedf,   j1-j0 )
        call datcp( vmwork(j4,isrc), wedf_line, j5-j4 )
      case( 2 )
! 'cntwcn' to vmwork
        vmwork(i0,isrc) = wflx
        vmwork(im,isrc) = wtot
        vmwork(iu,isrc) = wsum
        vmwork(iv,isrc) = wnrm
        vmwork(iw,isrc) = wion
        vmwork(ix,isrc) = wabs
        vmwork(iy,isrc) = wpmp
        vmwork(iz,isrc) = werr

        vmwork(i3:i4-1,isrc) = wswe(0:ndmc)
        vmwork(i4:i5-1,isrc) = wswi(0:ndmc)

        vmwork(ij:ik-1,isrc) = wemt(1:ndwp)
        vmwork(ik:il-1,isrc) = wwal(1:ndwp)
        vmwork(ip:iq-1,isrc) = wfhta(1:ndwp)
        vmwork(iq:ir-1,isrc) = wfhtm(1:ndwp)
        vmwork(ir:is-1,isrc) = wehta(1:ndwp)
        vmwork(is:it-1,isrc) = wehtm(1:ndwp)
        vmwork(j1:j2-1,isrc) = weema(1:ndwp)
        vmwork(j2:j3-1,isrc) = weemm(1:ndwp)
        vmwork(j3:j4-1,isrc) = wees(1:ndwp)

        vmwork(il:im-1,isrc) = wend(1:10)

        vmwork(it:iu-1,isrc) = wreg(0:10)

        call datcp( wssn,   vmwork(i1,isrc), i2-i1 )
        call datcp( wssp,   vmwork(i2,isrc), i3-i2 )
        call datcp( wsbr,   vmwork(i5,isrc), i6-i5 )
        call datcp( wden,   vmwork(i6,isrc), i7-i6 )
        call datcp( weng,   vmwork(i7,isrc), i8-i7 )
        call datcp( wvlp,   vmwork(i8,isrc), i9-i8 )
        call datcp( wnfl0x, vmwork(i9,isrc), ia-i9 )
        call datcp( wnfl0y, vmwork(ia,isrc), ib-ia )
        call datcp( wnfl0z, vmwork(ib,isrc), ic-ib )
        call datcp( wtion,  vmwork(ic,isrc), id-ic )

        call datcp( wgden,  vmwork(id,isrc), ie-id )
        call datcp( wgeng,  vmwork(ie,isrc), ig-ie )
        call datcp( wnflgx, vmwork(ig,isrc), ih-ig )
        call datcp( wnflgy, vmwork(ih,isrc), ii-ih )
        call datcp( wnflgz, vmwork(ii,isrc), ij-ii )

        call datcp( whta,   vmwork(in,isrc), io-in )
        call datcp( whtm,   vmwork(io,isrc), ip-io )

        call datcp( wedf,   vmwork(j0,isrc), j1-j0 )
        call datcp( wedf_line, vmwork(j4,isrc), j5-j4 )
      case( 3 )
! initialization of 'cntwcn'
        wflx = 0.0d0
        wtot = 0.0d0
        wsum = 0.0d0
        wnrm = 0.0d0
        wion = 0.0d0
        wabs = 0.0d0
        wpmp = 0.0d0
        werr = 0.0d0

        wswe(0:ndmc)  = 0.0d0
        wswi(0:ndmc)  = 0.0d0

        wemt(1:ndwp)  = 0.0d0
        wwal(1:ndwp)  = 0.0d0
        wfhta(1:ndwp) = 0.0d0
        wfhtm(1:ndwp) = 0.0d0
        wehta(1:ndwp) = 0.0d0
        wehtm(1:ndwp) = 0.0d0
        weema(1:ndwp) = 0.0d0
        weemm(1:ndwp) = 0.0d0
        wees(1:ndwp)  = 0.0d0

        wend(1:10)    = 0.0d0

        wreg(0:10)    = 0.0d0

        wssn(0:ndmc, 1:ndgs)  = 0.0d0
        wssp(0:ndmc, 1:ndgs)  = 0.0d0
        wsbr(0:ndmc, 1:ndgs)  = 0.0d0
        wden(0:ndmc, 1:ndgs)  = 0.0d0
        weng(0:ndmc, 1:ndgs)  = 0.0d0
        wvlp(0:ndmc, 1:ndgs)  = 0.0d0
        wnfl0x(0:ndmc,1:ndgs) = 0.0d0
        wnfl0y(0:ndmc,1:ndgs) = 0.0d0
        wnfl0z(0:ndmc,1:ndgs) = 0.0d0
        wtion(0:ndmc,1:ndgs)  = 0.0d0

        wgden(0:ndmc,1:2)  = 0.0d0
        wgeng(0:ndmc,1:2)  = 0.0d0
        wnflgx(0:ndmc,1:2) = 0.0d0
        wnflgy(0:ndmc,1:2) = 0.0d0
        wnflgz(0:ndmc,1:2) = 0.0d0

        whta(1:30,1:3) = 0.0d0
        whtm(1:30,1:3) = 0.0d0

        wedf(1:nene_mx,1:nwmx,1:nsmx) = 0.0d0
        wedf_line(1:nene_mx_line,1:nwmx,1:nsmx) = 0.0d0
      end select

      return
      end
