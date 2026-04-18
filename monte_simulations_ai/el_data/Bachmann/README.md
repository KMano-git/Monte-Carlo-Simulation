Bachmann working files

This directory stores two kinds of Bachmann-related assets.

1. Operational D + D+ data extracted from the current
   `../dd_00_elastic.cdf`, which the repository documentation describes as
   Bachmann-based.
2. Raw memo excerpts for the Bachmann coefficient tables, copied out of
   `../DpD_fit_memo.md` into a local file so that later cleanup can happen in
   isolation.

Generated files

- `bachmann_dd_cross_section_from_cdf.csv`
- `bachmann_dd_transport_tables_from_cdf.csv`
- `bachmann_dd_scattering_angle_from_cdf.csv`
- `bachmann_dd_scalars_from_cdf.json`
- `bachmann_raw_coefficients_excerpt.md`
- `bachmann_dd_from_split_tables_compat.cdf`

Generator scripts

- `extract_bachmann_cdf_data.py`
- `build_bachmann_compat_cdf.py`

Regeneration

```bash
python3 extract_bachmann_cdf_data.py
python3 build_bachmann_compat_cdf.py
```

Notes

- The extraction script intentionally does not implement the original
  Bachmann H.1 fit evaluator yet. The raw OCR coefficients in the memo still
  need a clean primary-source transcription before a trustworthy evaluator can
  be written.
- Until that cleanup is done, the safest reproducible Bachmann dataset in this
  repository is the one already embedded in `../dd_00_elastic.cdf`.
- `build_bachmann_compat_cdf.py` reconstructs a `.cdf`-compatible file from
  the extracted split CSV/JSON assets so that the repository has a verified
  round-trip path between the operational data tables and the compact CDF/CDL
  layout that the simulation code reads.
