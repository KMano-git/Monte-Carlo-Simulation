Krstic D+ + D cross-section working files

This directory stores reproducible intermediate data derived from the
Krstic fit coefficients summarized in `../DpD_fit_memo_v2.md`.

Recent DCS-based update

- The Krstic workflow is now split into two layers:
  - Existing `Krstic/` scripts keep the earlier elastic-total / composite workflow.
  - New root-level helpers in `../` rebuild pure-elastic tables directly from
    `g_pure(theta, E) = g_el_total(theta, E) - g_se(theta, E)`.
- New helper scripts:
  - `../krstic_dcs.py`
  - `../build_krstic_angle_cdf.py`
  - `../calc_I_kernel.py`
- New defaults write their outputs into this `Krstic/` directory and keep the
  filenames explicitly Krstic-specific, so they do not overwrite the generic
  `el_data/` outputs.
- `../calc_I_kernel.py` no longer uses surrogate `R_theta(E)` from the angle
  CDF. It now rebuilds `sigma_mt(E_cm)` from direct Krstic DCS moments and
  uses `E_cm = 0.5 * E_lab` for D + D.

Scope

- Reliable, structured data are provided for the integral cross-sections:
  - elastic
  - momentum-transfer
  - viscosity
  - spin-exchange
- `elastic` and `momentum-transfer` are sampled from the analytic fit
  coefficients.
- `viscosity` and `spin-exchange` are sampled from log-log interpolation of
  the page-104 tabulated values, because the current OCR snapshot of the memo
  does not preserve those analytic coefficients cleanly enough to regenerate a
  trustworthy fit.
- Output tables are sampled on the same logarithmic energy grids that the
  current `dd_00_elastic.cdf` workflow uses:
  - 101-point `cross_section` grid
  - 51-point `scattering_angle` energy grid
- Differential-fit coefficients are now curated from the manual-entry
  section in `../DpD_fit_memo_v2.md`.

Generated files

- `krstic_dd_integral_coeffs.json`
- `krstic_dd_integral_reference_points.csv`
- `krstic_dd_integral_cdf_grid_101.csv`
- `krstic_dd_integral_angle_grid_51.csv`
- `krstic_dd_dcs_coeffs.json`
- `krstic_dd_dcs_integral_checks.csv`
- `krstic_dd_scattering_angle_tabulated_31.csv`
- `krstic_dd_composite_compat.cdf`
- `krstic_dd_scattering_angle_compat.csv`
- `krstic_dd_transport_compat.csv`
- `krstic_dd_angle_validation.json`

New pure-elastic DCS outputs

- `krstic_pure_dcs_coeffs.json`
- `krstic_pure_dcs_validation.csv`
- `krstic_pure_dcs_validation.json`
- `krstic_pure_scattering_angle_compat.csv`
- `krstic_pure_transport_from_dcs.csv`
- `krstic_dd_pure_dcs_compat.cdf`
- `krstic_dd_pure_el_angle_fixed.cdf`

Generator scripts

- `generate_krstic_integral_data.py`
- `krstic_dcs.py`
- `build_krstic_composite_cdf.py`

New helper scripts outside this directory

- `../krstic_dcs.py`
- `../build_krstic_angle_cdf.py`
- `../calc_I_kernel.py`

Regeneration

```bash
python3 generate_krstic_integral_data.py
python3 build_krstic_composite_cdf.py
python3 ../krstic_dcs.py
python3 ../build_krstic_angle_cdf.py
python3 ../calc_I_kernel.py
```

Notes

- The fit is tabulated in the memo for 0.1-100 eV (center-of-mass).
- The generated CSV files are sampled on the exact logarithmic energy grids
  encoded in the repository `dd_00_elastic.cdf` metadata, and each row carries
  an `in_tabulated_fit_range` flag so that extrapolated values can be
  identified explicitly.
- `build_krstic_composite_cdf.py` now uses the manually curated Krstic
  differential-fit coefficients directly. It:
  - parses the `Elastic` and `Spin Exchange` coefficient tables from
    `../DpD_fit_memo_v2.md`
  - writes a machine-readable coefficient dump to `krstic_dd_dcs_coeffs.json`
  - integrates the direct Krstic `Elastic` DCS at the 31 tabulated energies
  - interpolates the inverse CDF in log-energy onto the runtime 51-point
    `scattering_angle` energy grid
  - recomputes `reaction_rate`, `I_1_0`, `I_1_1*up`, `I_1_2*up^2`, and
    `sigv_max` so that the file remains internally consistent
- The runtime angle grid extends below the tabulated DCS range. For
  `E < 0.1 eV`, the angle generator clamps to the 0.1 eV direct-DCS slice.
- `krstic_dd_dcs_integral_checks.csv` compares the direct DCS integrals against
  the Krstic integral elastic / momentum-transfer fits at the 31 tabulated
  energies.
- `krstic_dd_angle_validation.json` records the remaining mismatch between:
  - direct DCS and integral-fit transport ratios at the tabulated energies
  - the runtime 51-point angle table and the log-energy-interpolated direct DCS
- The new pure-elastic DCS workflow stores both:
  - a full rebuilt CDF file, `krstic_dd_pure_dcs_compat.cdf`
  - an I-kernel-only regenerated file, `krstic_dd_pure_el_angle_fixed.cdf`
- Validation for the new workflow is written to
  `krstic_pure_dcs_validation.json/csv`.
- Current limitation:
  - many tabulated energies still produce locally negative
    `g_pure = g_el_total - g_se` from the present OCR/manual coefficients.
  - the current default is `warn-clip`, so these tables are useful for
    implementation plumbing and comparison work, but the coefficients still
    need a tighter physical/data review before being treated as final
    production tables.
