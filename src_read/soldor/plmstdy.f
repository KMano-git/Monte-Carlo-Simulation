!***********************************************************************
      subroutine plmstdy(ktmax,ntrmax,dtstep)
!***********************************************************************
!
!      j-direction  only jc = jmd1
!
!-----------------------------------------------------------------------
      use cplcom, only : dq1a, dq2a, dq3, dq4, dteq, lbcgd, nion, nlp
     >    , q1a, q2a, q3, q4, vna, vva, vte, vti, wq1a, wq2a, wq3, wq4
      use cplmet, only : icel, icmax, icmpe, icspx, icwl1, itmax, itmpe
     >    , itpve, itpvs, itsls, jcel, jcxp1, jcxp2, jtmax
     >    , jtmin
      use cplqcn, only : qsx_cd, qsx_df, qsx_vh
      use cplvpn, only : qsx_vp
      use cpmpls, only : jmd1, lnedg
      use csonic, only : dtim, itim, kpcn, time
      use cunit,  only : n6
      implicit none
!
!::argument
      integer, intent(in) :: ktmax, ntrmax
      real(8), intent(in) :: dtstep
!
!::local variables
      integer  ktm, jc, j, imax, i, ntr
      integer  it
      integer  jm
      integer  ia, m1a, m2a, m3, m4
!
      write(n6,'(/2x,"*** plmstdy ***")')
!
!-----------------------------------------------------------------------
!::input data
!-----------------------------------------------------------------------
      lbcgd = 0   !  fixed value
      dtim  = dtstep
!
      write(n6,'(2x,"lnedg =",i2,"  lbcgd =",i2,"  dtim =",1pe12.4)')
     >   lnedg, lbcgd, dtim
!
!-----------------------------------------------------------------------
!::metric data  sum up along j-direction
!-----------------------------------------------------------------------
      call plmmtr
!
!::source terms in main plasma
      call plsrcn_clr
      call plsrcn_man
!
!-----------------------------------------------------------------------
!::dQ/dt = 0 in SOL & no viscous heating
!-----------------------------------------------------------------------
!
!::j-direction
      jm = jmd1
      jcxp1 = jm - 1
      jcxp2 = jm + 1
      do it = 1, itmax
        jtmin(it) = 1
        jtmax(it) = 3
        if( it.ge.itpvs .and. it.le.itpve ) cycle
        jcel(1,it) = jm
        jcel(2,it) = jm
        jcel(3,it) = jm
      enddo
!
!::set dteq ( due to no visocus heating in main plasma )
      do ia = 1, nion
        m2a = 2*ia
        dteq(m2a) = 1.0d-6
      enddo
!
!::set va = 0 in separatrix
      do ia = 1, nion
        i = icspx
        do jc = jcxp1+1, jcxp2-1
          j = jc
          q2a(j,i,ia) = 0.0d0
        enddo
      enddo
      call plauxv
!
!-----------------------------------------------------------------------
!::steady state in main plasma
!-----------------------------------------------------------------------
      itim = 0
      time = 0.0d0
      write(n6,'(/2x,"plmstdy itim =",i3,"  time =",1pe12.3)') itim,time
      write(n6,'(6x,"j",4x,"i",4x,"Na",9x,"Va",9x,"Ti",9x,"Te")')
      do i = 1, icmpe-1
        write(n6,'(2x,2i5,1p13e11.2)')
     >   jm,i,vna(jm,i,1),vva(jm,i,1),vti(jm,i),vte(jm,i)
      enddo
!
      do ktm = 1, ktmax
!
!::source terms in main plasma
        call plsrcn_clr
        call plsrcn_man
!
!::store q(N)
        do ia = 1, nion
          do j  = jcxp1+1, jcxp2-1
            do i  = icwl1, icmpe
              wq1a(j,i,ia) = q1a(j,i,ia)
              wq2a(j,i,ia) = q2a(j,i,ia)
            enddo
          enddo
        enddo
        do j = jcxp1+1, jcxp2-1
          do i = icwl1, icmpe
            wq3(j,i) = q3(j,i)
            wq4(j,i) = q4(j,i)
          enddo
        enddo
!
!::iterlation
        do ntr = 1, ntrmax
          nlp = ntr
!
!::solve q(N+1,l+1)
          call plwtfy    !   weit for X
          call pmpstp
!
          do j = jcxp1+1, jcxp2-1
            imax  = icmax(j)
            do i = 1, imax
              do ia = 1, nion
                q1a(j,i,ia) = q1a(j,i,ia) + dq1a(j,i,ia)
                q2a(j,i,ia) = q2a(j,i,ia) + dq2a(j,i,ia)
              enddo     ! loop(ia)
              q3(j,i) = q3(j,i) + dq3(j,i)
              q4(j,i) = q4(j,i) + dq4(j,i)
            enddo     ! loop(i)
          enddo     ! loop(j)
!
          call plauxv
        enddo
!
        itim = ktm
        time = time + dtim
!
        kpcn = 0
        if( itim.le.5 ) kpcn = 1
        if( mod(itim,200).eq.0 ) kpcn = 1
!
!::flux
        call plqcon_cl
        do it = itsls+1, itmpe-1   ! <== KS 2011/09/22
          if( it.ge.itpvs .and. it.le.itpve ) cycle
          call pxrdyv(it)
          call pxrdvp(it)
        enddo
!
!::debug write
        if( kpcn.eq.1 ) then
          write(n6,'(/2x,"plmstdy itim =",i3,"  time =",1pe12.3)')
     >      itim,time
          write(n6,'(6x,"j",4x,"i",4x,"Na",8x,"Va",8x,"Ti",
     >    8x,"Te",8x,"Fi",8x,"Qi",8x,"Qe",8x,"Qi_cnv",4x,
     >    "Qi_cnd",4x,"Qi_vis",4x,"Qe_cnv",4x,"Qe_cnd",4x,"Qe_vis")')
          ia  = 1
          m1a = 2*ia - 1
          m2a = 2*ia
          m3  = 2*nion + 1
          m4  = 2*nion + 2
          do i = icspx, icmpe-1
            call plqcon_sx(jcxp1+1,jcxp2-1,i,i,0)
            write(n6,'(2x,2i5,1p16e10.2)')
     >       jm,i,vna(jm,i,1),vva(jm,i,1),vti(jm,i),vte(jm,i),
     >       -qsx_df(m1a),-(qsx_df(m3)+qsx_cd(m3)+qsx_vh),
     >       -(qsx_df(m4)+qsx_cd(m4)),
     >       -qsx_df(m3), -qsx_cd(m3), -qsx_vh,
     >       -qsx_df(m4), -qsx_cd(m4), 0.0d0,
     >       -qsx_vp(m1a), -qsx_vp(m3), -qsx_vp(m4)
          enddo
          i = icmpe
          write(n6,'(2x,2i5,1p13e11.2)')
     >     jm,i,vna(jm,i,1),vva(jm,i,1),vti(jm,i),vte(jm,i)
        endif
!
      enddo     ! loop(k)
      end
!
!**********************************************************************
      subroutine pmpstp
!**********************************************************************
!
!   dq/dt*vl + viscos + pressure = source
!
!                                               copy pystep
!----------------------------------------------------------------------
      use cntmnt, only : ssn, ssp, swe, swi
      use cplcom, only : ama, c00, c11, cc, dd, dq1a, dq2a, dq3, dq4
     >    , dtrateq, ee, ff, lbcgd, nion, q1a, q2a, q3, q4, ss, wq1a
     >    , wq2a, wq3, wq4
      use cplmet, only : hvol, icmax, jcxp1, jcxp2
      use csize,  only : ndeq, ndxy
      use csonic, only : dtim, lfopt
      implicit none
!
!::local variables
      integer  j, imax, imax1, nequ, n, m, i, ia
      integer  m1a, m2a, m3, m4, nd1, nd2, neq, nmx
      real*8   dtvl, dtyy
!
      dtyy = dtim
      if( lfopt(2).eq.0 ) dtyy = 0.0d0
!
!----------------------------------------------------------------------
!::plasma parameter at the main plasma
!----------------------------------------------------------------------
      call plbman
!
!----------------------------------------------------------------------
!::colum : vertical direction
!----------------------------------------------------------------------
      do j = jcxp1+1, jcxp2-1
        imax  = icmax(j)
        imax1 = imax-1
!
!----------------------------------------------------------------------
!::yacobian matrix
!----------------------------------------------------------------------
!
!::number of equations
        nequ = 2*nion + 2
        m3   = 2*nion + 1
        m4   = 2*nion + 2
!
!::clear
        call pycler(j)
!
!::viscosity term
        call pyvisc(j)
        call pyvpin(j)
!
!::loss in scrape-off layer
        call pyslos(j)
!
!::source terms
        do i = 2, imax1
          do ia = 1, nion
            m1a = 2*ia - 1
            m2a = 2*ia
            ff(m1a,i) = ff(m1a,i) + ssn(j,i,ia)*ama(ia)*hvol(j,i)
            ff(m2a,i) = ff(m2a,i) + ssp(j,i,ia)*hvol(j,i)
          enddo
          ff(m3, i) = ff(m3, i) + swi(j,i)*hvol(j,i)
          ff(m4, i) = ff(m4, i) + swe(j,i)*hvol(j,i)
        enddo
!
!-----------------------------------------------------------------------
!::factor dtim/hvol(j,i) & term I
!-----------------------------------------------------------------------
        do n = 1, nequ
          do m = 1, nequ
            do i = 2, imax1
              dtvl = dtyy/hvol(j,i)*dtrateq(j,i,n)
              dd(n,m,i) = dd(n,m,i) + ss(n,m,i)
              cc(n,m,i) = cc(n,m,i)*dtvl
              dd(n,m,i) = dd(n,m,i)*dtvl
              ee(n,m,i) = ee(n,m,i)*dtvl
            enddo
          enddo
        enddo
!
        do n  = 1, nequ
          do i = 2, imax1
            dtvl = dtyy/hvol(j,i)*dtrateq(j,i,n)
            dd(n,n,i) = dd(n,n,i) + c11
            ff(n,i)   = ff(n,i)*dtvl
          enddo
        enddo
!
!::dt_V*R = -(q(N+1,l)-q(N)) + ...
        do ia = 1, nion
          m1a = 2*ia - 1
          m2a = 2*ia
          do i = 2, imax1
            ff(m1a,i) = ff(m1a,i)-(q1a(j,i,ia)-wq1a(j,i,ia))
            ff(m2a,i) = ff(m2a,i)-(q2a(j,i,ia)-wq2a(j,i,ia))
          enddo
        enddo
        do i = 2, imax1
          ff(m3,i) = ff(m3,i)-(q3(j,i)-wq3(j,i))
          ff(m4,i) = ff(m4,i)-(q4(j,i)-wq4(j,i))
        enddo
!
!-----------------------------------------------------------------------
!::dq^bar in RHS
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!::boundary condition : sol region
!-----------------------------------------------------------------------
!      cc = 0.0, dd = 1.0, ee =0.0, ff = 0.0
!      dtrat(j,i) = 0.0 in sol & prv region
!
!-----------------------------------------------------------------------
!::boundary condition iw = imax
!-----------------------------------------------------------------------
        do n = 1, nequ
          do m = 1, nequ
            cc(n,m,imax) = c00
            dd(n,m,imax) = c00
            ee(n,m,imax) = c00
          enddo
          dd(n,n,imax) = c11
          ff(n,imax)   = c00
        enddo
!
!-----------------------------------------------------------------------
!::boundary condition iw = 1 & imax
!-----------------------------------------------------------------------
        if( lbcgd.eq.1 ) then
          call pybcon(j)
        else
          call pybcon2(j)
        endif
!
!-----------------------------------------------------------------------
!::solve
!-----------------------------------------------------------------------
        nd1  = ndeq
        nd2  = ndxy
        neq  = nequ
        nmx  = imax
        call tdmslv( nd1, nd2, neq, nmx, cc, dd, ee, ff )
!
!-----------------------------------------------------------------------
!::store the solution
!-----------------------------------------------------------------------
        do i = 1, imax
          do ia = 1, nion
            m1a = 2*ia - 1
            m2a = 2*ia
            dq1a(j,i,ia) = ff(m1a,i)
            dq2a(j,i,ia) = ff(m2a,i)
          enddo
          dq3(j,i) = ff(m3,i)
          dq4(j,i) = ff(m4,i)
        enddo
!
      enddo
      end
