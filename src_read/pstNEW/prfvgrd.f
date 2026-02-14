!***********************************************************************
      subroutine prfvgrd(nt,np,pls,pne,pni,pte,pti,gne,gni,gte,gti,ndm)
!***********************************************************************
!
!    prfvgrd : profile of value and gradient of plasma parameter
!
!      Lp  : poloidal length            (j)
!
!      Ls  : position along B-field     (j)
!      Ti  : temperature at cell center (j)
!      gTi : ion temp. gradient at cell center (j) 
!           gTi(j) = (Ti(j+1/2)-Ti(j-1/2))/{Ls(j+1/2)-Ls(j-1/2)}
!                    (Tih(j)-Tih(j-1))/{Lsh(j)-Lsh(j-1)}
!
!      lsh : position of cell boundary (j+1/2)
!      Tih : density at cell boundary  (j+1/2)
!
!      Ls(j), Ti(j), gTi(j)
!           gTi(j=1) Not defined at dummy cell in SOL
!           gTi(j=1) = gTi(jte-1)  gTi(jte) = gTi(2) in Main plasma
!
!
!
!        dp                                          dp
!        |-o-|--o--|--o--|--o--|--o--|--o--|--o--|-o-|
!          1    2     3           j          J-1   J
!            X     X     X     X     X     X     X
!                            j-1/2  j+1/2
!                    dLs(j):   <---->
!            |---------------------->  Lsh(j+1/2) = Lsh(j)
!            |------------------->     Ls(j) = Ls(j)
!
!         dLs(j) : Basic distance
!         Lsh(3/2) = 0.0
!         Lsh(j+1/2) = Lsh(j-1/2) + dLs(j)  j=2,3...J-1
!
!         Ls(1) = Lsh(3/2) - dLs(1) 
!         Ls(j) = Lsh(j-1/2) + dLs(j)       j=2,3...J
!
!          j+1/2 => j,  j-1/2 => j-1    
!
!-----------------------------------------------------------------------
      use csize
      use cplcom
      use cplmet
      use csonic
      use cunit
      implicit none
!

  ! jtmin, jtmax
  ! vne, vte, vti
  ! n6

!
!::argument
      integer :: nt, np, ndm
      real(8) :: pls(ndm), pne(ndm), pni(ndm), pte(ndm), pti(ndm)
      real(8) :: gne(ndm), gni(ndm), gte(ndm), gti(ndm)
!
!::local variables
      integer :: it,  nft
      integer :: jwst, jwen, jw, j, i
      real(8) :: alnx(ndm), dls(ndm), plsh2(ndm)
      real(8) :: pstate(ndm), pstati(ndm), pdyni(ndm)
      real(8) :: pva(ndm), pvai(ndm)

      real(8) :: plsh(ndm), pneh(ndm), pnih(ndm), pteh(ndm), ptih(ndm)
      character(3)  :: cno
      character(80) :: dsn, drimp
      real(8) :: dsm, dsp, wsm, wsp

      integer :: lp = 0
      save :: lp, drimp

!::dsn name
      lp = lp + 1
      if( lp.eq.1 ) then
        call getenv("IMP1D",drimp)
      endif
!
      it = nt
      write(cno,'(i3)') it
      dsn = trim(drimp) // "/it" // trim(cno) // "PLS.txt"
      call delspc(dsn)
      write(n6,'(2x,"*** prfvgrd ***  dsn = ",a)') trim(dsn)
!
!::position and distance between cell boundaries
      call plenc(it,alnx)  ! Lp(j)
      call slenb(it,plsh2) ! Lsh(j+1/2)
      call gdlen(it,dls)   ! dLs(j)
!
      jwst = jtmin(it)
      jwen = jtmax(it)

!xx      if( it.gt.itmps ) then
!xx        jwst = jwst + 1
!xx        jwen = jwen - 1
!xx      endif
!
!::Lsh(j+1/2),Ls(j)
      plsh(jwst) = 0.0d0
      do jw = jwst+1, jwen-1
        plsh(jw) = plsh(jw-1) + dls(jw)
      enddo
      plsh(jwen) = 0.0d0

      pls(jwst) = plsh(jwst) - 0.5d0*dls(jwst)
      do jw = jwst+1, jwen
        pls(jw) = plsh(jw-1) + 0.5d0*dls(jw)
      enddo
!
!::plasma parameter at cell center (j)
      do jw = jwst, jwen
        j = jcel(jw,it)
        i = icel(jw,it)
        pne(jw) = vne(j,i)
        pni(jw) = vni(j,i)
        pte(jw) = vte(j,i)
        pti(jw) = vti(j,i)
        pstate(jw) = 1.602d-19 * pne(jw) * pte(jw)
        pstati(jw) = 1.602d-19 * pni(jw) * pti(jw)
        pdyni(jw) = 3.35d-27 * pni(jw) * vva(j,i,1) * vva(j,i,1)
        pva(jw) = vva(j,i,1)
        pvai(jw) = -1.d0 * vva(j,i,1)
      enddo
!
!::plasma parameter at cell boundary (j+1/2)
      do jw = jwst, jwen-1
        dsm = 0.5d0*dls(jw)
        dsp = 0.5d0*dls(jw+1)
        wsm = dsm/(dsm+dsp)
        wsp = dsp/(dsm+dsp) 
!
        pneh(jw) = wsp*pne(jw)+wsm*pne(jw+1)
        pnih(jw) = wsp*pni(jw)+wsm*pni(jw+1)
        pteh(jw) = wsp*pte(jw)+wsm*pte(jw+1)
        ptih(jw) = wsp*pti(jw)+wsm*pti(jw+1)
      enddo
      jw = jwen
      pneh(jw) = 0.0d0
      pnih(jw) = 0.0d0
      pteh(jw) = 0.0d0
      ptih(jw) = 0.0d0
!
!::gradient
      do jw = jwst+1, jwen-1
        gne(jw) = (pneh(jw)-pneh(jw-1))/(plsh(jw)-plsh(jw-1))
        gni(jw) = (pnih(jw)-pnih(jw-1))/(plsh(jw)-plsh(jw-1))
        gte(jw) = (pteh(jw)-pteh(jw-1))/(plsh(jw)-plsh(jw-1))
        gti(jw) = (ptih(jw)-ptih(jw-1))/(plsh(jw)-plsh(jw-1))
      enddo
!
!::Not define gradient in dummy cell
      if( it.lt.itmps ) then
        gne(jwst) = 0.0d0
        gni(jwst) = 0.0d0
        gte(jwst) = 0.0d0
        gti(jwst) = 0.0d0
        gne(jwen) = 0.0d0
        gni(jwen) = 0.0d0
        gte(jwen) = 0.0d0
        gti(jwen) = 0.0d0
!
!::cyclic condition in the main plasma
      else
        gne(jwst) = gne(jwen-1)
        gni(jwst) = gni(jwen-1)
        gte(jwst) = gte(jwen-1)
        gti(jwst) = gti(jwen-1)
        gne(jwen) = gne(jwst+1)
        gni(jwen) = gni(jwst+1)
        gte(jwen) = gte(jwst+1)
        gti(jwen) = gti(jwst+1)
      endif
!
!::output
      nft = 101
      open(unit=nft,file=dsn)
      write(nft,'(2x,a,2x,i4)') trim(dsn), it

      write(nft,'(4x,"it",3x,"jw",3x,"Lp",10x,"Ls",10x,"dLs",9x,
     >  "pne",9x,"pni",9x,"pte",9x,"pti",9x,
     >  "pva",9x,"pvai",8x,"pse",9x,"psi",9x,"pdi",9x,
     >  "Gne",9x,"Gni",9x,"Gti",9x,"Gte",9x,
     >  "jwh",4x,"Lsh",9x,"Lsh2",8x,"pneh",8x,"pnih",8x,"ptih")')

      np = jwen
      do jw = jwst, jwen
      write(nft,'(2x,2i5,1p16e12.4,0pf7.2,1p2e12.5,1p3e12.4)')
     >  it, jw, alnx(jw), pls(jw), dls(jw),
     >  pne(jw), pni(jw), pte(jw), pti(jw), 
     >  pva(jw), pvai(jw),
     >  pstate(jw), pstati(jw), pdyni(jw),
     >  gne(jw), gni(jw), gti(jw), gte(jw),
     >  dfloat(jw)+0.5d0, plsh(jw), plsh2(jw),
     >  pneh(jw), pnih(jw), ptih(jw)
      enddo

      close(nft)

      return
      end
