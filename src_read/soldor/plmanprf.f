!**********************************************************************
      subroutine plmanprf(kset)
!**********************************************************************
!
!     plasma profile in main plasma  at mid-plane
!
!   kset =  1 : anamp, atimp(tokrcv) ==> man_ni, man_ti : set
!        =  2 : man_ni, man_ti ==> anamp, atimp         : correction
!        =  3 : profile at mid plane
!
!      ift = 21
!      dsn = "chk_manprf"
!      inquire(file=trim(dsn),exist=lex)
!      if( lex ) then
!        open(unit=ift,file=dsn,position="append")
!      else
!        open(unit=ift,file=dsn)
!      endif
!
!----------------------------------------------------------------------
      use clocal, only : man_na, man_ne, man_ni, man_te, man_ti
      use cplcom, only : anamp, anemp, animp, atemp, atimp, vna, vne
     >    , vni, vte, vti
      use cplmet, only : icaxs, icmpe, icmps, icwl1
      use cpmpls, only : jmd1, rohmp, romp
      use csize,  only : ndsp, ndy
      use csonic, only : itim, time
      implicit none
!
!::argument
      integer    kset
!
!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer :: ia, i, icy, jc, ic, ii
      integer :: ia, icy, jc, ic, ii
      real(8), dimension(ndsp) :: zdna
      real(8) :: zdni, zdti, zdne, zdte
      integer :: ift
! deleted 4 lines organize local variables and include files by kamata 2021/05/31
!ik   character(20) :: dsn
!ik   logical :: lex
!ik   integer :: lp = 0
!ik   save :: lp
!
!::important variables
! deleted 3 lines dynamic allocation of arrays by kamata 2022/05/29
!ik   real(8), dimension(ndy,ndsp) :: man_na
!ik   real(8), dimension(ndy) :: man_ni, man_ti, man_ne, man_te
!ik   save :: man_na, man_ni, man_ti, man_ne, man_te
!
      ift = 7000 + 20
!
!----------------------------------------------------------------------
!::anamp,atimp(tokrcv) ==> man_na, man_ti  : set
!----------------------------------------------------------------------
      if( kset.eq.1 ) then
      ia = 1
! modified 2/1 lines organize local variables and include files by kamata 2021/05/31
!ik   do i = icaxs, icmps, -1
!ik   icy = i
      do icy = icaxs, icmps, -1
      man_na(icy,ia) = anamp(icy,ia)
      man_ni(icy) = animp(icy)
      man_ti(icy) = atimp(icy)
      man_ne(icy) = anemp(icy)
      man_te(icy) = atemp(icy)
      enddo
!
      icy = icmpe
      write(ift,'(2x,"*** plmanprf ***  set man_ti ",i6,1pe14.6,
     >   "  animp,atimp,atemp =",1p3e12.3)')
     >    itim, time, animp(icy), atimp(icy), atemp(icy)
      return
      endif
!
!----------------------------------------------------------------------
!::man_ti ==> atimp         : correction
!----------------------------------------------------------------------
      if( kset.eq.2 ) then
      ia = 1
      jc = jmd1
      ic = icmpe
      zdna(ia) = vna(jc,ic,ia) - man_na(ic,ia)
      zdni     = vni(jc,ic)    - man_ni(ic)
      zdti     = vti(jc,ic)    - man_ti(ic)
      zdne     = vne(jc,ic)    - man_ne(ic)
      zdte     = vte(jc,ic)    - man_te(ic)
!
      do ic = icmps, icaxs
      anamp(ic,ia) = man_na(ic,ia) + zdna(ia)
      animp(ic)    = man_ni(ic)    + zdni
      atimp(ic)    = man_ti(ic)    + zdti
      anemp(ic)    = man_ne(ic)    + zdne
      atemp(ic)    = man_te(ic)    + zdte
      enddo
!
      icy = icmpe
      write(ift,'(2x,"*** plmanprf ***  set atimp ",i6,1pe14.6,
     >   "  animp,atimp,atemp =",1p3e12.3)')
     >    itim, time, animp(ic), atimp(ic), atemp(ic)
      return
      endif
!
!----------------------------------------------------------------------
!::profile at mid-plane
!----------------------------------------------------------------------
      if( kset.eq.3 ) then
!
      ia = 1
      jc = jmd1
      ic = icmpe
      write(ift,'(2x,"*** plmanprf ***  out atimp ",i6,1pe14.6,
     >   "  animp,atimp,atemp =",1p3e12.3)')
     >    itim, time, animp(ic), atimp(ic), atemp(ic)
      write(ift,'(2x,"jc, ic =",2i5,"  vni,vti,vte =",1p3e12.3)')
     >  jc, ic, vni(jc,ia),vti(jc,ic),vte(jc,ic)
      write(ift,'(5x,"ic",3x,"i",4x,"roh",9x,"ro",10x,"vni",9x,
     >  "vti",9x,"vte",9x,"animp",7x,"atimp",7x,"atemp")')
!
      ii = 0
      do ic = icaxs, icwl1, -1
      ii = ii + 1
      write(ift,'(2x,2i5,1p8e12.3)')
     >  ic, ii, rohmp(ic), romp(ic), vni(jc,ic), vti(jc,ic), vte(jc,ic),
     >    animp(ic), atimp(ic), atemp(ic)
      enddo
!
      return
      endif
!
      end
