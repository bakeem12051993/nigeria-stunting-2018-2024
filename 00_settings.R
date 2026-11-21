# =============================================================================
#  Population shift, survey design and early-life shocks: explaining the
#  rise in child stunting in Nigeria, 2018 to 2023-24
#
#  00_settings.R: folders, file names, DHS codes and analysis choices.
#  Every other script reads this file first. Edit settings here only.
#
#  Author : Akeem A. Bakare
#  Version: 2.0, November 2026
# =============================================================================


# ---- Project folder ----------------------------------------------------------
# Set this to the folder that holds the scripts and data/ (forward slashes,
# also on Windows). If it does not exist, the current working folder is used.
project_dir <- "path/to/project"
if (dir.exists(project_dir)) {
  setwd(project_dir)
} else message("project_dir not found; using the current working folder: ", getwd())

dhs_dir      <- "data/dhs"        # unzipped DHS folders (NGKR8BDT/, NGGE8AFL/, ...)
external_dir <- "data/external"   # UCDP, food prices, ERA5-Land, floods
derived_dir  <- "derived"         # intermediate files written by the scripts
out_dir      <- "output"          # tables and figures
for (d in c(derived_dir, out_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)


# ---- DHS files ------------------------------------------------------------------
# Found inside dhs_dir by recode type and DHS phase (7 = 2018, 8 = 2023-24), so
# the version letter (e.g. NGKR7BFL, NGGE8AFL) does not need to be typed.
# Standard DHS files have a letter version (7B, 8A); malaria and other interim
# surveys have a number (NGKR81FL, the 2021 MIS) and are ignored.
# To force a file, add a line after this block, e.g.
#   dhs_files[["2018"]]$kr <- "NGKR7BDT/NGKR7BFL.dta"
dhs_phase <- c("2018" = "7", "2024" = "8")
dhs_find <- function(type, phase) {
  ext <- if (type == "GE") "dbf" else "dta"
  f <- list.files(dhs_dir, pattern = sprintf("^NG%s%s[A-Z]FL\\.%s$", type, phase, ext),
                  recursive = TRUE, ignore.case = TRUE)
  if (length(f)) sort(f, decreasing = TRUE)[1] else sprintf("NG%s%s?FL.%s", type, phase, ext)
}
dhs_files <- lapply(dhs_phase, function(p)
  list(kr = dhs_find("KR", p), pr = dhs_find("PR", p), hr = dhs_find("HR", p), ge = dhs_find("GE", p)))
round_label <- c("2018" = "2018", "2024" = "2023-24")

# Zone and state variables (matched across rounds by label, not code).
# In 2023-24, v024 numbers the zones differently from 2018, so szone is used.
# CHECK 2018: confirm v024 and sstate.
zone_var  <- c("2018" = "v024",   "2024" = "szone")
state_var <- c("2018" = "sstate", "2024" = "sstate1")
# The same in the household member (PR) and household (HR) recodes
pr_zone_var  <- c("2018" = "hv024",   "2024" = "shzone")
pr_state_var <- c("2018" = "shstate", "2024" = "shstate")

# ---- Weighting -----------------------------------------------------------------------
# Children were measured only in the one-in-three households picked for the men's
# interview. In both rounds those households hold more young children in the South
# than the other households do, so with the published weights the measured children
# over-represent the South (lower stunting).
#   "state_residence"  weight measured children to the state and urban/rural mix of
#                      all children under 5 in the same survey (main analysis)
#   "zone"             the same, by zone only
#   "published"        DHS weights unchanged (kept as a sensitivity analysis)
weighting <- "state_residence"

# DHS response codes. 2023-24 confirmed against the data. CHECK 2018.
codes <- list(
  water_unimproved = c(32, 42, 43),
  water_na         = c(96, 97, 98, 99),
  toilet_improved  = c(11, 12, 13, 15, 21, 22, 41),
  toilet_open      = c(31),
  toilet_na        = c(96, 97, 98, 99),
  finished         = 30:39,
  clean_fuel       = list("2018" = c(1, 2, 3, 4),            # CHECK for DHS-7
                          "2024" = c(1, 2, 3, 4, 5, 6)),     # DHS-8, confirmed
  facility_birth   = 20:49,
  dob_complete     = c(0, 1)
)
show_code_labels <- TRUE   # print value labels so codes can be checked


# ---- External data files --------------------------------------------------------
# See README.md for where to download each one. 00_organize_files.R puts each
# download in its folder; files are then found by name pattern, so they can
# keep the names they were downloaded with.
first_file <- function(dir, pattern, shown) {
  f <- list.files(dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  if (length(f)) sort(f, decreasing = TRUE)[1] else file.path(dir, shown)
}

# UCDP Georeferenced Event Dataset, CSV version (ged261-csv.zip holds GEDEvent_v26_1.csv)
ucdp_file       <- first_file(file.path(external_dir, "ucdp"), "^(GEDEvent|ged).*\\.csv$", "GEDEvent*.csv")
ucdp_last_month <- "2025-12"     # last month the file covers (GED 26.1 ends Dec 2025)

# World Bank Real-Time Food Prices, Nigeria market-level CSV (NGA_RTFP_mkt_*.csv)
prices_file  <- first_file(file.path(external_dir, "prices"), "RTFP.*mkt.*\\.csv$", "NGA_RTFP_mkt_*.csv")
prices_index <- "c_food_price_index"   # close estimate of the market food price index
prices_max_km <- 100   # clusters farther than this from any market take their zone's
                       # average series from the same file (NULL: always nearest market)

# ERA5-Land monthly averaged reanalysis: the .nc files go in this folder, under
# any names. Each variable (t2m, d2m, tp) is found inside them by name, whether
# the three come in one file or in several.
era5_dir         <- file.path(external_dir, "era5")
era5_first_month <- "1991-01"    # used only if the files carry no time information

# Optional: cluster-month flood flags built from Copernicus Global Flood
# Monitoring (columns: round, v001, month "YYYY-MM", flooded 0/1)
flood_file <- file.path(external_dir, "floods/cluster_month_floods.csv")


# ---- Exposure choices --------------------------------------------------------------
conflict_radius_km   <- 25          # main radius; sensitivity radii below
conflict_radius_sens <- c(10, 50)
ucdp_max_precision   <- 3           # where_prec 1-3: exact, within 25 km, or ADM2
baseline_years       <- 1991:2020   # climate normal for anomalies
heat_threshold_sd    <- 1           # hot month: temperature anomaly above 1 SD
heat_detrend         <- TRUE        # heat and temperature anomalies against each
                                    # cluster's own trend (FALSE: fixed 1991-2020 normal)
drought_threshold    <- -1          # drought month: SPEI-3 below -1
wet_threshold_sd     <- 1.5         # extremely wet month: rainfall anomaly above 1.5 SD

# Exposure windows, in months relative to the birth month (0 = birth month).
# Conception is taken as 9 months before birth.
windows <- list(
  preg = c(-9, -1),      # pregnancy
  y1   = c(0, 11),       # first year of life
  y2   = c(12, 23),      # second year of life
  pre  = c(-21, -10)     # placebo: the 12 months before conception
)
post_window <- c(1, 12)  # placebo: the 12 months after the interview


# ---- Inference ---------------------------------------------------------------------
n_boot <- 500            # bootstrap replicates for the decomposition
set.seed(20261119)
options(survey.lonely.psu = "adjust", stringsAsFactors = FALSE)
