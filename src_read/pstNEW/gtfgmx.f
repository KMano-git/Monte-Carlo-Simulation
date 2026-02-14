c***********************************************************************
      subroutine gtfgmx
c***********************************************************************
      use cunit
      implicit none

! move from cplfig2 (to be local valuable) 23/6/19
! for parallel profile plot
      real*4, save :: fpnemx_sol = 0.0_4, fpnemx_prv = 0.0_4
     >, fpnemx_man = 0.0_4
     >, fpnemn_sol = 0.0_4, fpnemn_prv = 0.0_4, fpnemn_man = 0.0_4
     >, fptemx_sol = 0.0_4, fptemx_prv = 0.0_4, fptemx_man = 0.0_4
     >, fptemn_sol = 0.0_4, fptemn_prv = 0.0_4, fptemn_man = 0.0_4
     >, fpvpmx_sol = 0.0_4, fpvpmx_prv = 0.0_4, fpvpmx_man = 0.0_4
     >, fpvpmn_sol = 0.0_4, fpvpmn_prv = 0.0_4, fpvpmn_man = 0.0_4
! for vertical profile plot
      real*4, save :: fvnemx_sol = 0.0_4, fvnemx_prv = 0.0_4
     >, fvnemx_man = 0.0_4
     >, fvnemn_sol = 0.0_4, fvnemn_prv = 0.0_4, fvnemn_man = 0.0_4
     >, fvtemx_sol = 0.0_4, fvtemx_prv = 0.0_4, fvtemx_man = 0.0_4
     >, fvtemn_sol = 0.0_4, fvtemn_prv = 0.0_4, fvtemn_man = 0.0_4
     >, fvvpmx_sol = 0.0_4, fvvpmx_prv = 0.0_4, fvvpmx_man = 0.0_4
     >, fvvpmn_sol = 0.0_4, fvvpmn_prv = 0.0_4, fvvpmn_man = 0.0_4
! for history plot
      real*4, save :: fhnemx_sol = 0.0_4, fhnemx_prv = 0.0_4
     >, fhnemx_man = 0.0_4
     >, fhnemn_sol = 0.0_4, fhnemn_prv = 0.0_4, fhnemn_man = 0.0_4
     >, fhtemx_sol = 0.0_4, fhtemx_prv = 0.0_4, fhtemx_man = 0.0_4
     >, fhtemn_sol = 0.0_4, fhtemn_prv = 0.0_4, fhtemn_man = 0.0_4
     >, fhvpmx_sol = 0.0_4, fhvpmx_prv = 0.0_4, fhvpmx_man = 0.0_4
     >, fhvpmn_sol = 0.0_4, fhvpmn_prv = 0.0_4, fhvpmn_man = 0.0_4
     >, fhflmn = 0.0_4, fhflmx = 0.0_4
     >, fhitmn = 0.0_4, fhitmx = 0.0_4

      character cdsn*80
      logical   lex
      integer   nft
c
      namelist /upfgmx/
     >       fpnemx_sol, fpnemx_prv, fpnemx_man
     >     , fpnemn_sol, fpnemn_prv, fpnemn_man
     >     , fptemx_sol, fptemx_prv, fptemx_man
     >     , fptemn_sol, fptemn_prv, fptemn_man
     >     , fpvpmx_sol, fpvpmx_prv, fpvpmx_man
     >     , fpvpmn_sol, fpvpmn_prv, fpvpmn_man
     >     , fvnemx_sol, fvnemx_prv, fvnemx_man
     >     , fvnemn_sol, fvnemn_prv, fvnemn_man
     >     , fvtemx_sol, fvtemx_prv, fvtemx_man
     >     , fvtemn_sol, fvtemn_prv, fvtemn_man
     >     , fvvpmx_sol, fvvpmx_prv, fvvpmx_man
     >     , fvvpmn_sol, fvvpmn_prv, fvvpmn_man
     >     , fhnemx_sol, fhnemx_prv, fhnemx_man
     >     , fhnemn_sol, fhnemn_prv, fhnemn_man
     >     , fhtemx_sol, fhtemx_prv, fhtemx_man
     >     , fhtemn_sol, fhtemn_prv, fhtemn_man
     >     , fhvpmx_sol, fhvpmx_prv, fhvpmx_man
     >     , fhvpmn_sol, fhvpmn_prv, fhvpmn_man
     >     , fhflmn, fhflmx
     >     , fhitmn, fhitmx
c
c:: default
      fpnemn_sol =  0.0;  fpnemn_prv =  0.0;  fpnemn_man =  0.0
      fptemn_sol =  0.0;  fptemn_prv =  0.0;  fptemn_man =  0.0
      fpvpmn_sol =  0.0;  fpvpmn_prv =  0.0;  fpvpmn_man =  0.0
      fvnemn_sol =  0.0;  fvnemn_prv =  0.0;  fvnemn_man =  0.0
      fvtemn_sol =  0.0;  fvtemn_prv =  0.0;  fvtemn_man =  0.0
      fvvpmn_sol =  0.0;  fvvpmn_prv =  0.0;  fvvpmn_man =  0.0
      fhnemn_sol =  0.0;  fhnemn_prv =  0.0;  fhnemn_man =  0.0
      fhtemn_sol =  0.0;  fhtemn_prv =  0.0;  fhtemn_man =  0.0
      fhvpmn_sol =  0.0;  fhvpmn_prv =  0.0;  fhvpmn_man =  0.0
c
      fpnemx_sol = 99.0;  fpnemx_prv = 99.0;  fpnemx_man = 99.0
      fptemx_sol = 99.0;  fptemx_prv = 99.0;  fptemx_man = 99.0
      fpvpmx_sol = 99.0;  fpvpmx_prv = 99.0;  fpvpmx_man = 99.0
      fvnemx_sol = 99.0;  fvnemx_prv = 99.0;  fvnemx_man = 99.0
      fvtemx_sol = 99.0;  fvtemx_prv = 99.0;  fvtemx_man = 99.0
      fvvpmx_sol = 99.0;  fvvpmx_prv = 99.0;  fvvpmx_man = 99.0
      fhnemx_sol = 99.0;  fhnemx_prv = 99.0;  fhnemx_man = 99.0
      fhtemx_sol = 99.0;  fhtemx_prv = 99.0;  fhtemx_man = 99.0
      fhvpmx_sol = 99.0;  fhvpmx_prv = 99.0;  fhvpmx_man = 99.0
c
      fhflmn = -0.5d22;   fhflmx = 2.0d22
      fhitmn = 0.0;  fhitmx = 99.0
c
      nft = 21
      call nopen(nft,"inpfig","text",cdsn,lex)
c
c::input /upfgmx/
      write(6,'(2x,"##1111 gtfgmx upfgmx nft = ",i6)') nft
      read(nft,upfgmx,end=900,err=900)
      close(nft)
c
c::debug
      write(n6,'(/2x,"*** gtfgmx ***")')
      write(n6,'(4x,6x,4x,"sol",9x,"prv",9x,"man")')
 602  format(4x,a,1p3e12.3)
      write(n6,602) "para-profle"
      write(n6,602) "fpnemx",  fpnemx_sol, fpnemx_prv, fpnemx_man
      write(n6,602) "fpnemn",  fpnemn_sol, fpnemn_prv, fpnemn_man
      write(n6,602) "fptemx",  fptemx_sol, fptemx_prv, fptemx_man
      write(n6,602) "fptemn",  fptemn_sol, fptemn_prv, fptemn_man
      write(n6,602) "fpvpmx",  fpvpmx_sol, fpvpmx_prv, fpvpmx_man
      write(n6,602) "fpvpmx",  fpvpmn_sol, fpvpmn_prv, fpvpmn_man
      write(n6,602) "vart-profle"
      write(n6,602) "fvnemx",  fvnemx_sol, fvnemx_prv, fvnemx_man
      write(n6,602) "fvnemn",  fvnemn_sol, fvnemn_prv, fvnemn_man
      write(n6,602) "fvtemx",  fvtemx_sol, fvtemx_prv, fvtemx_man
      write(n6,602) "fvtemn",  fvtemn_sol, fvtemn_prv, fvtemn_man
      write(n6,602) "fvvpmx",  fvvpmx_sol, fvvpmx_prv, fvvpmx_man
      write(n6,602) "fvvpmx",  fvvpmn_sol, fvvpmn_prv, fvvpmn_man
      write(n6,602) "history"
      write(n6,602) "fhflmx",  fhflmn, fhflmx
cx      write(n6,602) "fhnemx",  fhnemx_sol, fhnemx_prv, fhnemx_man
cx      write(n6,602) "fhnemn",  fhnemn_sol, fhnemn_prv, fhnemn_man
cx      write(n6,602) "fhtemx",  fhtemx_sol, fhtemx_prv, fhtemx_man
cx      write(n6,602) "fhtemn",  fhtemn_sol, fhtemn_prv, fhtemn_man
cx      write(n6,602) "fhvpmx",  fhvpmx_sol, fhvpmx_prv, fhvpmx_man
cx      write(n6,602) "fhvpmx",  fhvpmn_sol, fhvpmn_prv, fhvpmn_man
      return
c
  900 continue
      call wexit("gtfgmx","during reading file")
c
      return
      end
