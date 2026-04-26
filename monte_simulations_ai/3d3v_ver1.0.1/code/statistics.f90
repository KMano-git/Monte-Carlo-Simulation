!===============================================================================
! Module: statistics
! Running statistics for source terms.
!===============================================================================
module statistics
   use constants, only: dp
   use data_types, only: running_stats_t, source_stats_t, estimator_stats_t, &
      score_stats_data, source_terms_t, estimator_score_t, score_data, N_REACTIONS
   implicit none

   private
   public :: reset_running_stats, update_running_stats
   public :: running_variance, running_stddev, running_stderr
   public :: reset_score_stats, update_score_stats

contains

   subroutine reset_running_stats(stats)
      type(running_stats_t), intent(out) :: stats

      stats%n = 0
      stats%mean = 0.0d0
      stats%m2 = 0.0d0
   end subroutine reset_running_stats

   subroutine update_running_stats(stats, value)
      type(running_stats_t), intent(inout) :: stats
      real(dp), intent(in) :: value

      real(dp) :: delta, delta2

      stats%n = stats%n + 1
      delta = value - stats%mean
      stats%mean = stats%mean + delta / real(stats%n, dp)
      delta2 = value - stats%mean
      stats%m2 = stats%m2 + delta * delta2
   end subroutine update_running_stats

   real(dp) function running_variance(stats)
      type(running_stats_t), intent(in) :: stats

      if (stats%n > 1) then
         running_variance = max(0.0d0, stats%m2 / real(stats%n - 1, dp))
      else
         running_variance = 0.0d0
      end if
   end function running_variance

   real(dp) function running_stddev(stats)
      type(running_stats_t), intent(in) :: stats

      running_stddev = sqrt(running_variance(stats))
   end function running_stddev

   real(dp) function running_stderr(stats)
      type(running_stats_t), intent(in) :: stats

      if (stats%n > 0) then
         running_stderr = running_stddev(stats) / sqrt(real(stats%n, dp))
      else
         running_stderr = 0.0d0
      end if
   end function running_stderr

   subroutine reset_source_stats(stats)
      type(source_stats_t), intent(out) :: stats
      integer :: k

      call reset_running_stats(stats%sn_plus)
      call reset_running_stats(stats%sn_minus)
      call reset_running_stats(stats%sn_net)
      do k = 1, 3
         call reset_running_stats(stats%sp(k))
      end do
      call reset_running_stats(stats%we)
      call reset_running_stats(stats%wi)
   end subroutine reset_source_stats

   subroutine reset_estimator_stats(stats)
      type(estimator_stats_t), intent(out) :: stats
      integer :: ir

      do ir = 1, N_REACTIONS
         call reset_source_stats(stats%reaction(ir))
      end do
   end subroutine reset_estimator_stats

   subroutine reset_score_stats(stats)
      type(score_stats_data), intent(out) :: stats

      call reset_estimator_stats(stats%cl)
      call reset_estimator_stats(stats%cl_avg)
      call reset_estimator_stats(stats%tr)
      call reset_estimator_stats(stats%tr_pretab)
   end subroutine reset_score_stats

   subroutine update_source_stats(stats, src, norm)
      type(source_stats_t), intent(inout) :: stats
      type(source_terms_t), intent(in) :: src
      real(dp), intent(in) :: norm

      call update_running_stats(stats%sn_plus, src%sn_plus * norm)
      call update_running_stats(stats%sn_minus, src%sn_minus * norm)
      call update_running_stats(stats%sn_net, src%sn_net * norm)
      call update_running_stats(stats%sp(1), src%sp(1) * norm)
      call update_running_stats(stats%sp(2), src%sp(2) * norm)
      call update_running_stats(stats%sp(3), src%sp(3) * norm)
      call update_running_stats(stats%we, src%we * norm)
      call update_running_stats(stats%wi, src%wi * norm)
   end subroutine update_source_stats

   subroutine update_estimator_stats(stats, score, norm)
      type(estimator_stats_t), intent(inout) :: stats
      type(estimator_score_t), intent(in) :: score
      real(dp), intent(in) :: norm
      integer :: ir

      do ir = 1, N_REACTIONS
         call update_source_stats(stats%reaction(ir), score%reaction(ir), norm)
      end do
   end subroutine update_estimator_stats

   subroutine update_score_stats(stats, score, norm, include_pretab)
      type(score_stats_data), intent(inout) :: stats
      type(score_data), intent(in) :: score
      real(dp), intent(in) :: norm
      logical, intent(in) :: include_pretab

      call update_estimator_stats(stats%cl, score%cl, norm)
      call update_estimator_stats(stats%cl_avg, score%cl_avg, norm)
      call update_estimator_stats(stats%tr, score%tr, norm)
      if (include_pretab) call update_estimator_stats(stats%tr_pretab, score%tr_pretab, norm)
   end subroutine update_score_stats

end module statistics
