!***********************************************************************
subroutine test_git_mesg   ! 2012/07/01)
!***********************************************************************
integer :: n6
call git_mesg(n6)
stop
end
!
!***********************************************************************
subroutine git_mesg(n6)
!***********************************************************************
implicit none
integer :: n6
!
 61   format(1x,a)
!
!::Message from git
write(n6,61)  "   "
write(n6,61)  "=== Load module ==="
write(n6,61)  "prj : /work/NSSDIV/mano/sonic/sonicV4"
write(n6,61)  "lod : ../LOD/plimp"
write(n6,61)  "dat : Fri Oct 17 17:29:52 JST 2025"
write(n6,61)  "   "
write(n6,61)  "=== git branch  ==="
write(n6,61)  "* develop"
write(n6,61)  "  master"
write(n6,61)  "   "
write(n6,61)  "=== git log  ==="
write(n6,61)  "commit 4893ff3713ef13cb9106a848da486b9ed7c71292"
write(n6,61)  "Author: Tatsuto Yamamoto <yamamoto.tatsuto@qst.go.jp>"
write(n6,61)  "Date:   Mon Sep 29 15:46:58 2025 +0900"
write(n6,61)  "    bugfix: set icalZ = 0 for set bkflw distribution at first."
write(n6,61)  "   "
write(n6,61)  "=== git status  ==="
write(n6,61)  "   "
write(n6,61)  "   "
   
return
end
