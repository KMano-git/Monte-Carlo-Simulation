c
      subroutine pack (n,st)
c                      o io
      character*(*) st
      j=0
      do 100 i=1,len(st)
      if (st(i:i).ne.' ') then
      j=j+1
      if (i.ne.j) st(j:j)=st(i:i)
      endif
  100 continue
      st(j+1:)=' '
      n=j
c
      return
      end
