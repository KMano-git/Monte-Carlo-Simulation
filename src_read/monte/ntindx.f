!**********************************************************************
      subroutine ntindx
!**********************************************************************
!
!      index check
!
!----------------------------------------------------------------------
      use cntcom, only : iaxs, impl, iplx, iply, ispx, ivs1, ivs2, iwl1
     >    , iwl2, mcel, migx, migy
      use cplmet, only : icel, icmpe, icmps, icspx, icwl1, icwl2, itmpe
     >    , itmps, itpve, itpvs, itsle, itsls, jcel
      use cunit,  only : n6
      implicit none
!
!::local variables
      integer  it, jt, j, i, ic
!
      write(n6,'(/2x,"*** ntindx ***")')
!
      write(n6,602) "icwl1 =",icwl1,  "itsls =",itsls,icel(1,itsls)
      write(n6,602) "icspx =",icspx,  "itsle =",itsle,icel(1,itsle)
      write(n6,602) "icwl2 =",icwl2,  "itpve =",itpve,icel(1,itpve)
      write(n6,602) "icmps =",icmps,  "itmps =",itmps,icel(1,itmps)
      write(n6,602) "icmpe =",icmpe,  "itmpe =",itmpe,icel(1,itmpe)
 602  format(2x,"monte: ",a,i5,2x,"soldor :",a,2i5)
!-----
      if( icmpe.eq.icel(1,itmpe) ) then
      write(n6,'(2x,"O.K. icmpe.eq.icel(1,itmpe)")')
      else
      write(n6,'(2x,"???  icmpe.ne.icel(1,itmpe)")')
      call wexit("ntindx","error")
      endif
!-----
!
      write(n6,'(2x)')
      write(n6,604) "itsls =",itsls,"  itsle =",itsle
      write(n6,604) "itpvs =",itpvs,"  itpve =",itpve
      write(n6,604) "itmps =",itmps,"  itmpe =",itmpe
      write(n6,604) "ivs1  =",ivs1, "  iwl1  =",iwl1, "  ispx  =",ispx
      write(n6,604) "iwl2  =",iwl2, "  ivs2  =",ivs2
      write(n6,604) "impl  =",impl, "  iaxs  =",iaxs
 604  format(2x,a,i5,a,i5,a,i5)
      write(n6,'(/2x,"jt",2x,"it",4x,"j",3x,"i",4x,"ic",4x,
     >  "mgx",2x,"mgy",2x,"ipx",2x,"ipy")')
      do it = itsls, itmpe
      jt = 2
      j  = jcel(jt,it)
      i  = icel(jt,it)
      ic = mcel(j,i)
      if( it.ne.itmpe ) cycle
      write(n6,'(2i4,2x,2i4,2x,i5,2i5,2x,2i5)')
     >  jt,it, j,i, ic, migx(ic), migy(ic), iplx(ic), iply(ic)
!-----
      if( it.eq.itmpe ) then
      if( migy(ic).eq.impl-1 ) then
        write(n6,'(2x,"O.K.  when it=itmpe, migy(ic) = impl-1")')
      else
        write(n6,'(2x,"???   when it=itmpe, migy(ic) .ne. impl-1")')
        call wexit("ntindx","error")
      endif
      endif
!-----
      enddo
!
      return
      end
