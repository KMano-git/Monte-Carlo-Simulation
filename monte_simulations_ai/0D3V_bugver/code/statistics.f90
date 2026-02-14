!===============================================================================
! Module: statistics
! 統計量の計算と出力
!===============================================================================
module statistics
   use constants, only: dp, M_D, J_TO_EV
   use particle_data, only: particle_0d
   implicit none

contains

   !-----------------------------------------------------------------------------
   ! エネルギー統計量の計算
   !-----------------------------------------------------------------------------
   subroutine compute_energy_stats(particles, n_particles, E_mean, E_std, T_eff)
      type(particle_0d), intent(in) :: particles(:)
      integer, intent(in) :: n_particles
      real(dp), intent(out) :: E_mean, E_std, T_eff

      real(dp) :: E_i, E_sum, E2_sum
      integer :: i, n_alive

      E_sum = 0.0d0
      E2_sum = 0.0d0
      n_alive = 0

      do i = 1, n_particles
         if (particles(i)%alive) then
            E_i = 0.5d0 * M_D * (particles(i)%vx**2 + particles(i)%vy**2 + &
               particles(i)%vz**2) * J_TO_EV
            E_sum = E_sum + E_i
            E2_sum = E2_sum + E_i**2
            n_alive = n_alive + 1
         end if
      end do

      if (n_alive > 0) then
         E_mean = E_sum / n_alive
         E_std = sqrt(max(E2_sum / n_alive - E_mean**2, 0.0d0))
         ! 有効温度: 3D理想気体 <E> = (3/2) k T → T = (2/3) <E>
         T_eff = (2.0d0 / 3.0d0) * E_mean
      else
         E_mean = 0.0d0
         E_std = 0.0d0
         T_eff = 0.0d0
      end if

   end subroutine compute_energy_stats

   !-----------------------------------------------------------------------------
   ! エネルギーヒストグラムの計算
   !-----------------------------------------------------------------------------
   subroutine compute_energy_histogram(particles, n_particles, n_bins, &
      E_min, E_max, hist, E_centers)
      type(particle_0d), intent(in) :: particles(:)
      integer, intent(in) :: n_particles, n_bins
      real(dp), intent(in) :: E_min, E_max
      integer, intent(out) :: hist(n_bins)
      real(dp), intent(out) :: E_centers(n_bins)

      real(dp) :: E_i, dE
      integer :: i, ibin

      ! 初期化
      hist(:) = 0
      dE = (E_max - E_min) / n_bins

      ! ビン中心のエネルギー
      do i = 1, n_bins
         E_centers(i) = E_min + (i - 0.5d0) * dE
      end do

      ! ヒストグラム計算
      do i = 1, n_particles
         if (particles(i)%alive) then
            E_i = 0.5d0 * M_D * (particles(i)%vx**2 + particles(i)%vy**2 + &
               particles(i)%vz**2) * J_TO_EV

            ! ビン番号を計算
            ibin = int((E_i - E_min) / dE) + 1

            ! 範囲内ならカウント
            if (ibin >= 1 .and. ibin <= n_bins) then
               hist(ibin) = hist(ibin) + 1
            end if
         end if
      end do

   end subroutine compute_energy_histogram

   !-----------------------------------------------------------------------------
   ! 緩和データの書き込み
   !-----------------------------------------------------------------------------
   subroutine write_relaxation_data(unit, time, E_mean, T_eff)
      integer, intent(in) :: unit
      real(dp), intent(in) :: time, E_mean, T_eff

      write(unit, '(3(ES15.6E3,","))') time, E_mean, T_eff

   end subroutine write_relaxation_data

   !-----------------------------------------------------------------------------
   ! ヒストグラムデータの書き込み
   !-----------------------------------------------------------------------------
   subroutine write_histogram_data(filename, hist, E_centers, n_bins)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: n_bins
      integer, intent(in) :: hist(n_bins)
      real(dp), intent(in) :: E_centers(n_bins)

      integer :: unit, i

      open(newunit=unit, file=trim(filename), status='replace')
      write(unit, '(A)') '# E_center[eV], count'

      do i = 1, n_bins
         write(unit, '(ES15.6E3,",",I0)') E_centers(i), hist(i)
      end do

      close(unit)

   end subroutine write_histogram_data

end module statistics
