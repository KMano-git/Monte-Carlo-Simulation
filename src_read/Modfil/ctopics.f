! added addition and reconsider of passed data by kamata 2022/02/23
      module ctopics
      implicit none
      integer :: ktop_cr  = 0, ktop_di  = 1, ktop_ese = 1, ktop_esi = 1
     >         , ktop_fni = 1, ktop_fqe = 1, ktop_fqi = 1, ktop_hp  = 1
     >         , ktop_ne  = 1, ktop_ni  = 1, ktop_pp  = 1, ktop_ps  = 1
     >         , ktop_sn  = 0, ktop_te  = 1, ktop_ti  = 1, ktop_we  = 0
     >         , ktop_wi  = 0, ktop_xe  = 1, ktop_xi  = 1
! ktop_cr  : radiation loss due to impurities
! ktop_di  : particle diffusion coefficient
! ktop_ese : energy   source by external, jules, nuclear fusion and arbitrary ( electron )
! ktop_esi : energy   source by external, jules, nuclear fusion and arbitrary ( ion )
! ktop_fni : particle flux ( ion )
! ktop_fqe : heat     flux ( electoron )
! ktop_fqi : heat     flux ( ion )
! ktop_hp  : heat     pinch
! ktop_ne  : electron density
! ktop_ni  : ion      density
! ktop_pp  : particle pinch
! ktop_ps  : particle source by NB, nuclear fusion and arbitrary
! ktop_sn  : particle source by neutral
! ktop_te  : electron trmperature
! ktop_ti  : ion      trmperature
! ktop_we  : energy source by neutral ( electron )
! ktop_wi  : energy source by neutral ( ion )
! ktop_xe  : thermal diffusivity ( electron )
! ktop_xi  : thermal diffusivity ( ion )
! the above variable is set to 1 when used

      end module ctopics
