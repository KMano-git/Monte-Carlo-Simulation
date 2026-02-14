!****************************************************************************
      subroutine interpolate_array(val_out,points,r_in,val_in,r_norm
     > ,n_norm)
!****************************************************************************
      implicit none
! arguments
      integer,intent(in)    :: points ! effective size of r_in and val_in
      integer,intent(in)    :: n_norm  ! effective size of r_norm
      real(8),intent(inout) :: r_in(points)
      real(8),intent(out)   :: val_out(n_norm), val_in(points)
      real(8),intent(in)    :: r_norm(n_norm)
! local
      integer ic, ic_point

      if(points <= 0) then
        return
      elseif(points == 1) then
        val_out(1:n_norm) = val_in(1)
      else
        call sort_array(r_in,val_in,points)
        do ic = 1, n_norm
          if(r_norm(ic) < r_in(1)) then
            val_out(ic) = val_in(1)
          elseif(r_norm(ic) > r_in(points)) then
            val_out(ic) = val_in(points)
          else
            do ic_point = 1, points-1
              ! interpolation of val_out(ic)
              if((r_in(ic_point)  <= r_norm(ic)) .and. 
     >           (r_in(ic_point+1) > r_norm(ic)) ) then
                call ievalue(val_out(ic),ic_point,val_in,r_norm,r_in
     >                       ,points,ic,n_norm)
              endif
            enddo
          endif
        enddo
      endif
!***********************************************************************   
      contains
!***********************************************************************      
      subroutine sort_array(r_in,val_in,points)
!***********************************************************************
      implicit none
! arguments
      integer,intent(in)    :: points ! effective size of r_in and val_in
      real(8),intent(inout) :: r_in(points)
      real(8),intent(out)   :: val_in(points)
! local
      integer i, j
      real(8) temp_omd, temp_diff

      ! bubble sort
      do i = 1, points-1
        do j = 1, points-i
          if (r_in(j) > r_in(j+1)) then
            ! exchange r_in
            temp_omd = r_in(j)
            r_in(j) = r_in(j+1)
            r_in(j+1) = temp_omd
            ! exchange temp
            temp_diff = val_in(j)
            val_in(j)   = val_in(j+1)
            val_in(j+1) = temp_diff
          endif
        enddo
      enddo
      end subroutine sort_array

!****************************************************************************
      subroutine ievalue(inter_out,ic_point,val_in,r_norm,r_in,points,ic
     > ,n_norm)
!****************************************************************************
! liner interpolation or extrapolation
      implicit none
! arguments
      integer,intent(in) :: ic_point, points, ic, n_norm
      real(8),intent(in) :: val_in(points)
     > ,r_norm(n_norm), r_in(points)
      real(8),intent(out) :: inter_out
!
      inter_out = val_in(ic_point+1)+
     >             (r_norm(ic)-r_in(ic_point+1))
     >            /(r_in(ic_point)-r_in(ic_point+1))
     >            *(val_in(ic_point)-val_in(ic_point+1))
      end subroutine ievalue

      end subroutine interpolate_array
