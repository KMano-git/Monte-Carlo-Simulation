!**********************************************************************
      subroutine mcplas
!**********************************************************************
!
!      set plasma parameter for monte-carlo calculaion
!
!          neutral(ig)      plasma ion(ia)
!
!----------------------------------------------------------------------
      use cntcom, only : ncmax, iplx, iply, mrgn, void_ne, void_ni
     >    , void_te, void_ti, volm
      use cntpls, only : dene, deni, teme, temi, vflw, zefm
      use cplcom, only : anamp, atemp, atimp, aza, nion, vnezef, vna
     >    , vte, vti, vva, vzf
      use csize,  only : ndgs, ndmc
      use csonic, only : kpcn
      implicit none
!
!::local variables
      integer  nsiz, nsizg, ic, icx, icy, ia
      real*8   zne, zvl, zni
      integer  ir
!
!::zero clear
      nsiz  = ndmc + 1
      nsizg = (ndmc+1)*ndgs
      call setd( dene, nsiz,  void_ne )
      call setd( deni, nsizg, void_ni )
      call setd( teme, nsiz,  void_te )
      call setd( temi, nsiz,  void_ti )
      call setd( vflw, nsizg, 0.0d0 )
!
!::zeff in thermal force  avoid zero division  see imorbit
      call setd( zefm, nsiz,  1.0d0 )
!
!::plasma parameter
      do ic = 1, ncmax
        icx = iplx(ic)
        icy = iply(ic)
        ir  = mrgn(ic)
        if( icx.le.0 .or. icy.le.0 ) cycle
!
!::hot core
        if( ir.eq.7 ) then
          teme(ic) = atemp(icy)
          temi(ic) = atimp(icy)
          zne = 0.0d0
          do ia = 1, nion
            deni(ic,ia) = anamp(icy,ia)
            vflw(ic,ia) = 0.0d0
            zne = zne + aza(ia)*anamp(icy,ia)
          enddo
          dene(ic) = zne
        else
!::soldor
          dene(ic) = vnezef(icx,icy)
          teme(ic) = vte(icx,icy)
          temi(ic) = vti(icx,icy)
          zefm(ic) = vzf(icx,icy)
          do ia = 1, nion
            deni(ic,ia) = vna(icx,icy,ia)
            vflw(ic,ia) = vva(icx,icy,ia)
          enddo
        endif
      enddo
!
!::vol & <Ni>
      zvl = 0.0d0
      zni = 0.0d0
      ia  = 1
      do ic = 1, ncmax
        if( mrgn(ic).ne.6 .and. mrgn(ic).ne.7 ) cycle
        zvl = zvl + volm(ic)
        zni = zni + deni(ic,ia)*volm(ic)
      enddo
      zni = zni/zvl
!
      if( kpcn.eq.0 ) return
!
!----------------------------------------------------------------------
!::debug write
!----------------------------------------------------------------------
      call trmark("mcplas","return")
      return
      end
