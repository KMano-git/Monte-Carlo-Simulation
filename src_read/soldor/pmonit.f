!***********************************************************************
      subroutine pmonit
!***********************************************************************
      use cntmnt, only : pfl_ion, psm_abs, psm_man, psm_pmp, tnpuf
      use cplcom, only : trad, vna, vne, vte, vti
      use cplmet, only : icspx, jcdp1, jcdp2
      use cplqcn, only : pcn_pdt
      use cplwrd, only : wfac
      use cpmpls, only : jmd1, jmd2
      use csonic, only : dtim, itim, time
      use cunit,  only : lmspe, lmype, wh_job
      implicit none
!
!::local variables
      character cjob*7
      integer   lp, nft, i, jx(10), iy
      real*8    unt, unp, unf, und, uns
      data lp/0/; save lp
!
      if( lmype.ne.lmspe ) return
!
      lp = lp + 1
!
      nft = 61   !  see sub. opnfrt
!
      jx(1) = jcdp1
      jx(2) = jcdp2
      jx(3) = jcdp1 + 1
      jx(4) = jcdp2 - 1
      jx(5) = jmd1
      jx(6) = jmd2
      iy    = icspx
!
      cjob = wh_job
!
      unt = 1.0d-3   ! unit of time
      uns = 1.0d-6   ! unit of dtim
      unp = 1.0d6    ! unit of power
      unf = 1.0d22   ! unit of flux
      und = 1.0d19   ! unit of density
!
      if( lp.eq.1 ) then
      write(nft,'(1x,64x,"i =",i3,6(4x,i3,4x,i3,4x,i3))')
     >   iy,(jx(i),jx(i),jx(i),i=1,6)
      write(nft,'(a,3x,a,5x,a,7x,a,5x,25(a,7x),a)')
     >    'job',  'itim', 'time', 'dtim', 'trad', 'wfac', 'DotN'
     >  , 'Finp', 'Fout', 'Fpmp', 'nedO', 'TedO', 'TidO', 'nedI'
     >  , 'TedI', 'TidI', 'necO', 'TecO', 'TicO', 'necI', 'TecI'
     >  , 'TicI', 'nesO', 'nisO', 'TesO', 'TisO', 'nesI', 'nisl'
     >  , 'TesI', 'TisI'
      endif
!
      write(nft,'(a3,1x,i6,1x,f11.3,1x,f7.3,26(1x,f10.3))')
     >  cjob(1:3), itim, time/unt, dtim/uns, trad/unp
     > ,wfac
     > ,pcn_pdt/unf, (-pfl_ion(5)+tnpuf)/unf
     > ,(psm_abs+psm_man)/unf, psm_pmp/unf
     > ,(vne(jx(i),iy)/und,vte(jx(i),iy),vti(jx(i),iy),i=1,4)
     > ,(vne(jx(i),iy)/und,vna(jx(i),iy,1)/und
     > ,vte(jx(i),iy),vti(jx(i),iy),i=5,6)
!
!::KSFUJI
      if( mod(lp,100).eq.0 ) call flush(nft)
!
      return
      end
