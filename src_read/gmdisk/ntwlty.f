!**********************************************************************
      subroutine ntwlty
!**********************************************************************
!
!   define name of <tywl> to sum up incident flux to wall
!
!           not distinguish divertor from wall
!
!                                                      <--
!          W3o        +-------------------------+  W3i  (3)
!           gat2      !                         !      gat1
!         ##    ##### !~~~~~~~~~~~~~~~~~~~~~~~~~! #####    ##
!       | #---------+ !                         ! +---------#
!     W | # Void    ! !     Main Plasma         ! !         #
!     4 V #---------! !-------------------------! !---------#  W2
!         #  OUT    ! !                         ! !  IN     #
!         #   DIV   ! !      Scrape-Off         ! !   DIV   #
!         #         ! !                         ! !         # A
!    (4)  #---------- --------------------------- ----------# | (2)
!         #         ! !      Void               ! !         # |
!         ########### ##################     #### ###########
!            W1o     --->   W1     (1)   gat3        W1i
!
!      when no vacume region,   gato = albo,  gati = albi
!
!
!                g1(9)    g2(10)     V3   (7)
!         ######    #####    ################################
!         #           |          |       #                  #
!         #         a1|        a2|       #p1                # g V2
!     V4  #           |          |       #                  # t
!         #         (13)       (14)      #p2                # v
!    (8)  #           |          |       #(11,12)           # 3 (6)
!         ####################################################
!                                 V1   (5)
!
!----------------------------------------------------------------------
      use cntcom, only : cwal, icwl, igtw, inwl, iplx, iply, nowl2
     >    , npew, npsw, tywl
      use csize,  only : ndwp
      use cunit,  only : n6
      use mod_externalgrid, only:use_exdata
      implicit none
!
!::local variables
      integer  i, nw, iws, iwe, iw, ic
      character cty*4, ctyb*4
!
      if(use_exdata) then
         call ntwlty_ex
         return
      endif
!
!::redefine tywl (temporary)
      write(n6,'(/2x,"*** ntwlty ***    set tywl")')
!
!::clear
      do i = 1, ndwp
      tywl(i) = "----"
      enddo
!
!----------------------------------------------------------------------
!::label for plasma wall and vacume wall
!----------------------------------------------------------------------
      write(n6,'(2x,1x,"nw",4x,"iw",2x,"inwl",2x,"igtw",2x,"tywl")')
!
      do nw = 1, nowl2
          iws = npsw(nw)
          iwe = npew(nw)
          ctyb = "xxxx"
          do iw = iws, iwe-1
            ic = icwl(iw)
            cty = cwal(nw)
            !--W2
            if( inwl(iw).eq.2 
     >          .and. iplx(ic).gt.0 .and. iply(ic).gt.0) then
                cty = "W2d"
            endif
            !--W4
            if( inwl(iw).eq.4 
     >          .and. iplx(ic).gt.0 .and. iply(ic).gt.0) then
                cty = "W4d"
            endif
!
!--gate in W and V
            if( cty(1:1).eq."W" .and. igtw(iw).gt.0) then
                 write( cty, '(a2,a,i1)') cty(1:2),"g",igtw(iw)
            endif
!
!::error
            if( cty.eq."????" ) then
              write(n6,'(2x,"find  unidentified plasma-wall   nw,iw ="
     >           ,2i6)') nw,iw
              call wexit("ntwlty","find Not identified plasma-wall")
            endif
!
!::define
            tywl(iw) = cty
            if( cty.ne.ctyb .or. iw.eq.iwe-1 ) then
                write(n6,'(2x,i3,3i6,2x,a)')
     >            nw, iw, inwl(iw), igtw(iw), tywl(iw)
            endif
            ctyb = cty
          enddo    !  loop (iw)
      enddo    !  loop (nw)
!
      return
      end

!**********************************************************************
      subroutine ntwlty_ex
!**********************************************************************
      use cntcom, only : cwal, icwl, igtw, inwl, nowl2
     >    , npew, npsw, tywl, ncmax
      use csize,  only : ndwp
      use cunit,  only : n6
      implicit none
!
!::local variables
      integer  i, nw, iws, iwe, iw, ic
      character cty*4, ctyb*4
!
!::redefine tywl (temporary)
      write(n6,'(/2x,"*** ntwlty_ex ***    set tywl")')
!
!::clear
      do i = 1, ndwp
      tywl(i) = "----"
      enddo
!
!----------------------------------------------------------------------
!::label for plasma wall and vacume wall
!----------------------------------------------------------------------
      write(n6,'(2x,1x,"nw",4x,"iw",2x,"inwl",2x,"igtw",2x,"tywl")')
!
      do nw = 1, nowl2
          iws = npsw(nw)
          iwe = npew(nw)
          ctyb = "xxxx"
          do iw = iws, iwe
            ic = icwl(iw)
            cty = cwal(nw)
            !--W2
            if( inwl(iw).eq.2 .and. ic.le.ncmax) then
              cty = "W2d"
            endif
            !--W4
            if( inwl(iw).eq.4 .and. ic.le.ncmax) then
              cty = "W4d"
            endif
!
!--gate in W and V
            if( cty(1:1).eq."W" .and. igtw(iw).gt.0) then
              write( cty, '(a2,a,i1)') cty(1:2),"g",igtw(iw)
            endif
!
!::error
            if( cty.eq."????" ) then
              write(n6,'(2x,"find  unidentified plasma-wall   nw,iw ="
     >           ,2i6)') nw,iw
              call wexit("ntwlty","find Not identified plasma-wall")
            endif
!
!::define
            tywl(iw) = cty
            if( cty.ne.ctyb .or. iw.eq.iwe ) then
                write(n6,'(2x,i3,3i6,2x,a)')
     >            nw, iw, inwl(iw), igtw(iw), tywl(iw)
            endif
            ctyb = cty
          enddo    !  loop (iw)
      enddo    !  loop (nw)
!
      return
      end subroutine ntwlty_ex