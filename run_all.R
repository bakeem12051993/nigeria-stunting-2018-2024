# =============================================================================
#  run_all.R: runs the whole analysis in order.
#  First time: copy your downloads into data/downloads and run
#  00_organize_files.R, which files them where these scripts look.
#  Then open this project folder in RStudio (or setwd() to it) and run this file.
#  Each step can also be run on its own once the steps before it have run.
#  Steps that cannot run yet (missing 2018 files or exposure data) are skipped
#  with a message, and the run carries on.
# =============================================================================

steps <- c("01_prepare_dhs.R",        # DHS files, official replication, quality, GPS
           "02_decomposition.R",      # descriptive change, shift-share, Oaxaca-Blinder
           "03_exposure_conflict.R",  # UCDP conflict around each cluster
           "04_exposure_prices.R",    # food inflation at the nearest market
           "05_exposure_climate.R",   # ERA5-Land climate anomalies, drought, floods
           "06_child_exposures.R",    # each child's exposure in each window
           "07_models.R",             # fixed-effects models
           "08_attribution.R",        # share of the change each shock explains
           "09_survey_checks.R",      # zone shares, anemia, measured households
           "10_figures.R")            # manuscript figures (output/figures)

for (s in steps) {
  cat("\n==============================", s, "==============================\n")
  tryCatch(source(s, echo = FALSE),
           skip_step = function(e) message("Skipped: ", conditionMessage(e)))
}
cat("\nFinished. Results are in the output folder; any step marked Skipped is waiting for data.\n")
