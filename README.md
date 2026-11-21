# Child stunting in Nigeria, 2018 to 2023-24: R code

Analysis code for "Population shift, survey design and early-life shocks:
explaining the rise in child stunting in Nigeria, 2018 to 2023-24", using the
2018 and 2023-24 Nigeria Demographic and Health Surveys (NDHS).

**No data are included.** The DHS terms do not allow the survey files, or
anything derived from them at the child or household level, to be shared. To
run the code, register at dhsprogram.com, request the Nigeria 2018 and 2023-24
surveys (KR, PR, HR and GE files, and the 2021 MIS for the optional checks), and
download the free external data listed in section 3. The `.gitignore` keeps
`data/`, `derived/` and `output/` out of the repository.

## 1. Folder layout

Set `project_dir` at the top of `00_settings.R` to your project folder and put
the scripts in it. You do not need to unzip or rename anything:

1. Create `data/downloads` in the project folder (or run `00_organize_files.R`
   once and it creates it for you).
2. Copy everything you downloaded into it, still zipped: the DHS zips for 2018
   and 2023-24, the UCDP zip, the World Bank food price CSV and the ERA5-Land
   files.
3. Run `00_organize_files.R`. It unzips everything (including zips inside zips)
   and files each dataset as below, then prints a checklist of what the
   analysis will find. Rerun it whenever you add downloads.

```
project/
  00_organize_files.R, 00_settings.R ... 09_survey_checks.R, run_all.R, README.md
  data/
    downloads/       what you downloaded, untouched
    dhs/
      NGKR7BDT/  NGPR7BDT/  NGHR7BDT/  NGGE7AFL/     (2018)
      NGKR8BDT/  NGPR8BDT/  NGHR8BDT/  NGGE8AFL/     (2023-24)
    external/
      ucdp/      GEDEvent_v26_1.csv
      prices/    NGA_RTFP_mkt_*.csv
      era5/      the .nc files, under any names
      floods/    cluster_month_floods.csv            (optional)
  derived/   (created by the scripts)
  output/    (created by the scripts)
```

The external files are found by name pattern, so they keep the names they were
downloaded with. For ERA5-Land, the three variables can arrive in one file or
several; the code finds t2m, d2m and tp inside them.

Keep `data/dhs` out of any shared or synced location that people outside the
registered team can open: the DHS terms do not allow redistribution of the
microdata.

## 2. Packages

```r
install.packages(c("haven", "survey", "terra", "foreign"))
install.packages("SPEI")   # recommended: gives SPEI-3; without it the code uses an SPI-type index
install.packages("ncdf4")  # recommended: reads the dates stored in Climate Data Store NetCDF files
```

The fixed-effects models use a base-R estimator (in `00_functions.R`), so
`fixest` is not required. Its coefficients and cluster-robust (CR1) standard
errors were checked against `lm()` with dummy variables and match exactly. If
you prefer `fixest` for the final tables, model B is
`feols(haz ~ x | cl_key + round_zone_age + round_zone_bmoy, weights = ~wt, cluster = ~psu)`.

## 3. External data (all free)

| Domain | Source | What to download |
|---|---|---|
| Conflict | UCDP Georeferenced Event Dataset (GED) 26.1, ucdp.uu.se/downloads | Global CSV (ged261-csv.zip), covering 1989 to 2025. License CC BY 4.0. Set `ucdp_last_month` in `00_settings.R` to the last month it covers. |
| Food prices | World Bank Real Time Food Prices (RTFP), microdata.worldbank.org, catalog 4503 | Nigeria market-level CSV. Needs columns `mkt_name`, `lat`, `lon`, `DATES`, `c_food_price_index`. |
| Climate | ERA5-Land monthly averaged reanalysis, Copernicus Climate Data Store (free account) | Variables `2m_temperature`, `2m_dewpoint_temperature`, `total_precipitation`; years 1991 to 2025; all months; area North 14.5, West 2.5, South 4, East 15; NetCDF. One request for all three is fine (the Store may send a zip of two .nc files). |
| Floods (optional) | Copernicus Global Flood Monitoring (GFM) | A CSV with columns `round` ("2018" or "2023-24"), `v001` (cluster number), `month` ("YYYY-MM") and `flooded` (1 if the cluster was flooded that month). Skip if not built; the code runs without it. |

## 4. Run order

Open `run_all.R` and source it, or run the scripts one at a time:

| Script | What it does | Main outputs |
|---|---|---|
| `01_prepare_dhs.R` | Replicates the official stunting figures and re-estimates them with corrected weights, harmonizes the two rounds, adds GPS and household buffers, sets the analysis weights | `01_official_replication.csv`, `01b_zone_mix_all_vs_measured.csv`, `02_data_quality.csv` |
| `02_decomposition.R` | Pooled asset index, descriptive change, shift-share, sampling-frame check and Oaxaca-Blinder; key results again with the published weights | `03a` to `07b` tables, including `04b_frame_constant.csv` |
| `03_exposure_conflict.R` | Conflict deaths within 10, 25 and 50 km of each cluster, by month | `derived/conflict_panels.rds` |
| `04_exposure_prices.R` | Year-on-year food inflation at the nearest market | `derived/price_panels.rds` |
| `05_exposure_climate.R` | Humid heat (wet-bulb), SPEI-3 drought, very wet months, hot-dry months | `derived/climate_panels.rds` |
| `06_child_exposures.R` | Each child's exposure in pregnancy, year 1, year 2, and the placebo windows | `derived/children_exposures.rds` |
| `07_models.R` | Models A to H (fixed effects: cluster, round by zone by age, round by zone by birth month) | `08_models.csv`, `fig_window_effects.png` |
| `08_attribution.R` | How much of the 2018 to 2023-24 change each shock domain explains | `09_attribution.csv`, `fig_attribution.png` |
| `09_survey_checks.R` | Where children under 5 live and anemia in 2018, 2021 (MIS, if present) and 2023-24; households picked for measuring versus the rest | `10` to `12` tables |
| `10_figures.R` | The six manuscript figures (shock timeline, zones, age, survey design, Oaxaca-Blinder, shocks by window), 300 dpi and 174 mm wide | `output/figures/fig1_timeline.png` to `fig6_shocks.png` |

Step 01 starts by printing which DHS files it found and runs on every round that
is complete. Any step that cannot run yet is marked "Skipped" with the reason,
and the run carries on:

- With 2023-24 only: 01 runs in full, 02 saves the one-round descriptive tables,
  and the decomposition (02) and attribution (08) wait for the 2018 files.
- Without the conflict, food price and climate files: 03 to 05 note what is
  missing, and 06 to 08 wait for at least one of them.

## 5. Checks before the first full run

1. **DHS files found**: step 01 (and `00_organize_files.R`) lists the file it
   found for each recode and round. The files are located by type and phase,
   so the 2018 version letters do not need to be typed.
2. **2018 zone and state variables**: `zone_var` uses `v024` for 2018 and
   `szone` for 2023-24; `state_var` uses `sstate` and `sstate1`. After step 01,
   the zone table must show the same six zone names in both rounds, and the
   "state names present in only one round" message should list none (or only
   spelling differences to fix, such as "fct-abuja").
3. **2018 response codes** (`codes` in `00_settings.R`): confirm clean cooking
   fuels (v161) for DHS-7, plus the water, toilet and facility-birth codes.
   Set `show_code_labels <- TRUE` (the default) to print the labels.
4. **Official replication**: `01_official_replication.csv` should give about
   37% stunting for 2018 and 39.5% for 2023-24 (2023-24 already replicates).
5. **Exposure coverage**: step 06 prints the share of children with a value for
   each exposure. A domain enters the models only if at least 80% of children in
   each round have it.

## 6. Weighting

Children were measured only in the one-in-three households picked for the men's
interview. In both rounds, the picked households in the South hold far more young
children (and women) than the other households, while in the North the two look
alike (see `12_measured_households_check.csv`; mothers' birth histories, collected
in every household, show the same gap). With the published weights the measured
children therefore over-represent the South, where stunting is lower.

`weighting` in `00_settings.R` sets the analysis weights:

- `"state_residence"` (default, main analysis): measured children are
  post-stratified to the state and urban/rural mix of all children under 5 in
  the same survey.
- `"zone"`: the same, by zone only.
- `"published"`: DHS weights unchanged.

`01_official_replication.csv` reports national stunting under all three, and
`07b_published_weights.csv` repeats the key decomposition results with the
published weights. Post-stratification treats the all-children totals as known,
so intervals for composition terms are somewhat narrow.

**Sampling frames.** The 2018 survey sampled from the 2006 census enumeration
areas; the 2023-24 survey (like the 2021 MIS) used the new frame prepared for the
planned census, and the North's share of households is higher under it.
`04b_frame_constant.csv` repeats the change in stunting and HAZ with both rounds
put on one round's household distribution across states (each state's child
weights times the ratio of the two rounds' household shares), so that only
children per household and stunting within states change.

## 7. Notes

- The placebo after the interview (model D) needs exposure data that reach 12
  months past the last 2023-24 interview, and is skipped until the files cover
  that period. With the cluster fixed effects it carries little information,
  because exposure after the interview varies mostly between clusters.
- Rural residence is constant within a cluster, so in model F only its
  interactions with the shocks are estimated.
- Send a PDF of any resulting publication to references@dhsprogram.com, as the
  DHS terms require.

## 8. License and citation

The code is released under the MIT License (see `LICENSE`). If you use it,
please cite the article. The DHS data must be obtained from the DHS Program and
used under its terms.
