! added dynamic allocation of arrays by kamata 2022/05/29
! allocate variables in cplqcn
      subroutine cplqcn_alc( kk )
      use cplqcn, only : mrgnp, qfx_cd, qfx_cv, qfx_df, qfx_vh, qfy_cd
     >    , qfy_df, qfy_vh, qvl_al, qvl_cl, qvl_dt, qvl_pe, qvl_pi
     >    , qvl_sc
      use csize,  only : nqfm => ndeq, ndmc, nqfa => ndsp, nqfx => ndx
     >    , nqfy => ndy
      implicit none
! arguments
      integer, intent(in) :: kk
! kk : flag = 1: allocate,    = 2: deallocate

! local variables
      integer    istat
      character  cmsg*80

      select case( kk )
      case( 1 )
! allocate
        if( .not. allocated( mrgnp ) ) then
          allocate( mrgnp(ndmc), qfx_cd(nqfx,nqfy,nqfm)
     >      , qfx_cv(nqfx,nqfy,nqfm), qfx_df(nqfx,nqfy,nqfm)
     >      , qfx_vh(nqfx,nqfy), qfy_cd(nqfx,nqfy,nqfm)
     >      , qfy_df(nqfx,nqfy,nqfm), qfy_vh(nqfx,nqfy)
     >      , qvl_al(nqfx,nqfy,nqfm), qvl_cl(nqfx,nqfy,nqfm)
     >      , qvl_dt(nqfx,nqfy,nqfm), qvl_pe(nqfx,nqfy,nqfm)
     >      , qvl_pi(nqfx,nqfy,nqfa), qvl_sc(nqfx,nqfy,nqfm)
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'mrgnp allocate error in cplqcn_alc, istat = ', istat
            call wexit( 'cplqcn', trim( cmsg ) )
          endif

! initial set
          mrgnp(1:ndmc) = 0
          qfx_cd(1:nqfx,1:nqfy,1:nqfm) = 0.0_8 
          qfx_cv(1:nqfx,1:nqfy,1:nqfm) = 0.0_8 
          qfx_df(1:nqfx,1:nqfy,1:nqfm) = 0.0_8 
          qfx_vh(1:nqfx,1:nqfy) = 0.0_8 
          qfy_cd(1:nqfx,1:nqfy,1:nqfm) = 0.0_8 
          qfy_df(1:nqfx,1:nqfy,1:nqfm) = 0.0_8 
          qfy_vh(1:nqfx,1:nqfy) = 0.0_8
          qvl_al(1:nqfx,1:nqfy,1:nqfm) = 0.0_8 
          qvl_cl(1:nqfx,1:nqfy,1:nqfm) = 0.0_8 
          qvl_dt(1:nqfx,1:nqfy,1:nqfm) = 0.0_8 
          qvl_pe(1:nqfx,1:nqfy,1:nqfm) = 0.0_8 
          qvl_pi(1:nqfx,1:nqfy,1:nqfa) = 0.0_8 
          qvl_sc(1:nqfx,1:nqfy,1:nqfm) = 0.0_8 
        endif
      case( 2 )
! deallocate
        if( allocated( mrgnp ) ) then
          deallocate( mrgnp, qfx_cd, qfx_cv, qfx_df, qfx_vh, qfy_cd
     >      , qfy_df, qfy_vh, qvl_al, qvl_cl, qvl_dt, qvl_pe, qvl_pi
     >      , qvl_sc
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'mrgnp deallocate error in cplqcn_alc, istat = ', istat
            call wexit( 'cplqcn', trim( cmsg ) )
          endif
        endif
      end select

      return
      end
