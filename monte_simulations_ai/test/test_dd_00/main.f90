program main
   use constants, only: dp, PI
   use cdf_reader

   implicit none
   integer :: ierr
   integer :: i, e
   real(dp) :: energy, sigma
   real(dp) :: p1, p2, chi1, chi2, dchi, dP_dchi
   real(dp) :: chi_mid, dsigma_dOmega
   real(dp) :: energies(4) = [0.1_dp, 1.0_dp, 10.0_dp, 100.0_dp]
   integer :: file_unit

   ! Initialize CDF
   call load_elastic_cdf('../../3d3v_event-driven/run/dd_00_elastic.cdf', ierr)
   if (ierr /= 0) then
      write(*,*) "Failed to load CDF file. Error: ", ierr
      stop
   end if

   ! Open CSV file
   file_unit = 10
   open(unit=file_unit, file='dcs_output.csv', status='replace')
   write(file_unit, '(A)') "Energy_eV,Angle_rad,Angle_deg,DCS_m2_sr"

   ! Loop over energies
   do e = 1, 4
      energy = energies(e)
      sigma = get_sigma_elastic(energy)

      write(*,*) "Processing energy: ", energy, " eV, Sigma: ", sigma, " m^2"

      ! Use exact probability grid spacing from cdf_reader (N_PROB=251 -> 250 intervals)
      do i = 1, 250
         p1 = (i - 1.0_dp) / 250.0_dp
         p2 = (i) / 250.0_dp

         p1 = max(1.0e-6_dp, p1)
         p2 = min(1.0_dp - 1.0e-6_dp, p2)

         chi1 = sample_scattering_angle(energy, p1)
         chi2 = sample_scattering_angle(energy, p2)

         dchi = chi2 - chi1
         if (i == 125) then
            write(*,*) "Debug i=125, p1=", p1, " p2=", p2, " chi1=", chi1, " chi2=", chi2, " dchi=", dchi
         end if
         if (abs(dchi) > 1.0e-12_dp) then
            dchi = abs(dchi)
            chi_mid = (chi1 + chi2) / 2.0_dp
            dP_dchi = (p2 - p1) / dchi

            if (sin(chi_mid) > 1.0e-12_dp) then
               dsigma_dOmega = (sigma / (2.0_dp * PI * sin(chi_mid))) * dP_dchi
               write(file_unit, '(F10.4, ",", F10.6, ",", F10.6, ",", E15.6)') &
                  energy, chi_mid, chi_mid * 180.0_dp / PI, dsigma_dOmega
            end if
         end if
      end do
   end do

   close(file_unit)
   write(*,*) "Output written to dcs_output.csv"
end program main
