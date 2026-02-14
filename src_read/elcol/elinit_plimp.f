!***********************************************************************
      subroutine elinit
!***********************************************************************
      use celcom,      only : agmi, el_agmi, el_angl, el_atm, el_cros
     >    , el_dnam, el_i_10, el_i_11, el_i_12, el_mol, el_rate, el_sgmx
     >    , ndel, sgmx
      use cunit,       only : lmspe, lmype, mype, mywld, n6
      use cxdcom,      only : nxs, xdnam, xwmlt
      use mod_shexe,   only : use_namelist, DEGAS_ELATM, DEGAS_ELMOL
      use nnIintg,     only : set_nnIintgCoef
      implicit none

!
!::local variables
      character  cenv*80, dnam*80
      integer    i, ixs
      real*8     zval, fmks
      character  csub*6; data csub/"elinit"/
! function
      real(8)    xdt_eval0
!
!-----------------------------------------------------------------------
!
      write(n6,'(/2x,"*** ",a," ***  start")') csub
      el_atm = 1
      el_mol = 2
!
!-----------------------------------------------------------------------
!::copy data of cdf-file
!-----------------------------------------------------------------------
      if( lmype.eq.lmspe ) then
!
         if(use_namelist) then
!..   use variable in inppls
            cenv = DEGAS_ELATM
         else
!..   use environment variable
            call getenv( "DEGAS_ELATM", cenv )
         end if
      call xdt_cpfl( cenv, "all", el_dnam(el_atm) )
!
      if(use_namelist) then
!..   use variable in inppls
         cenv = DEGAS_ELMOL
      else
!..   use environment variable
         call getenv( "DEGAS_ELMOL", cenv )
      end if
      call xdt_cpfl( cenv, "all", el_dnam(el_mol) )
!
      endif

!-----------------------------------------------------------------------
!::set interp-variables
!-----------------------------------------------------------------------
      call xdt_pset
!
!-----------------------------------------------------------------------
!::set index of variables
!-----------------------------------------------------------------------
      do i = 1, 2
      dnam = el_dnam(i)
      call xdt_indx( 1, dnam, "cross_section",     el_cros(i) )
      call xdt_indx( 2, dnam, "reaction_rate",     el_rate(i) )
      call xdt_indx( 2, dnam, "I_1_0",             el_I_10(i) )
      call xdt_indx( 2, dnam, "I_1_1*up",          el_I_11(i) )
      call xdt_indx( 2, dnam, "I_1_2*up^2",        el_I_12(i) )
      call xdt_indx( 2, dnam, "scattering_angle",  el_angl(i) )
      call xdt_indx( 0, dnam, "sigv_max",          el_sgmx(i) )
      call xdt_indx( 0, dnam, "angle_min",         el_agmi(i) )
      enddo
!
!-----------------------------------------------------------------------
!::debug write
!-----------------------------------------------------------------------
      write(n6,'(/2x,"*** ",a,"***  mype =",i4,"   index check")')
     >   csub,mype
      do i = 1, 2
      call el_indx( "cross_section",     el_cros(i) )
      call el_indx( "reaction_rate",     el_rate(i) )
      call el_indx( "I_1_0",             el_I_10(i) )
      call el_indx( "I_1_1*up",          el_I_11(i) )
      call el_indx( "I_1_2*up^2",        el_I_12(i) )
      call el_indx( "scattering_angle",  el_angl(i) )
      call el_indx( "sigv_max",          el_sgmx(i) )
      call el_indx( "angle_min",         el_agmi(i) )
      enddo

!-----------------------------------------------------------------------
!::set polynomial fit coefficients for n-n elastic collision
!-----------------------------------------------------------------------
      call set_nnIintgCoef

!-----------------------------------------------------------------------
!::sigv_max & angle_min
!-----------------------------------------------------------------------
      write(n6,'(/2x,"sigv_max & angle_min in MKS-unit")')
      do i = 1, 2
      ixs = el_sgmx(i)
      zval = xdt_eval0(ixs)
      fmks = xwmlt(1,ixs)
      sgmx(i) = zval*fmks
      write(n6,'(4x,"sgmx(",i1,")=",1pe12.3,2x,i3,2x,a16,2x,1p2e12.3)')
     >  i,sgmx(i),ixs,xdnam(ixs),zval,fmks
      ixs = el_agmi(i)
      zval = xdt_eval0(ixs)
      fmks = xwmlt(1,ixs)
      agmi(i) = zval*fmks
      write(n6,'(4x,"agmi(",i1,")=",1pe12.3,2x,i3,2x,a16,2x,1p2e12.3)')
     >  i,agmi(i),ixs,xdnam(ixs),zval,fmks
      enddo
!
      write(n6,'(/2x,"*** ",a," ***  end")') csub
!
      return
      end
!
!***********************************************************************
      subroutine el_indx(cvar,n)
!***********************************************************************
      use cunit,  only : n6
      use cxdcom, only : xdnam, xdrnk, xwmlt, xwunt, xwvar
      implicit none
!
!::argument
      character, intent(in) :: cvar*(*)
      integer,   intent(in) :: n
!
! function
      integer    lenx
!
      write(n6,'(2x,i3,2x,a18,2x,a18,2x,i2,1pe12.3,2x,a8)')
     >   n,xdnam(n),xwvar(1,n),xdrnk(n),xwmlt(1,n),xwunt(1,n)
!
      if( xwvar(1,n)(1:lenx(xwvar(1,n))).ne.cvar ) then
        call wexit("el_indx","not matched cvar")
      endif
!
      return
      end
