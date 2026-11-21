# =============================================================================
#  06_child_exposures.R: each child's exposure in each window
#  Windows (months relative to birth): pregnancy (-9 to -1), first year (0-11),
#  second year (12-23); placebo windows before conception (-21 to -10) and
#  after the interview (+1 to +12).
#  Only months the child had lived by the interview count. A window not yet
#  reached counts as no exposure; the age fixed effects in 07 absorb this.
#
#  Exposure measures (suffix _preg, _y1, _y2, _pre, _post, and _1000 for
#  pregnancy through age 2 combined):
#    conf      log(1 + conflict deaths within 25 km)          [conf10, conf50: 10/50 km]
#    food      mean food inflation (percentage points / 10) at the nearest market
#    heat      mean humid-heat anomaly (SD)                     [temp: air temperature]
#    dry       mean drought index (SD; positive = drier)
#    wet       number of extremely wet months
#    hotdry    number of months both hot and in drought
#    flood     number of flooded months (only if flood data were supplied)
#  Writes derived/children_exposures.rds
# =============================================================================

source("00_settings.R"); source("00_functions.R")
ch <- readRDS(file.path(derived_dir, "children.rds"))
conflict <- readRDS(file.path(derived_dir, "conflict_panels.rds"))
prices   <- readRDS(file.path(derived_dir, "price_panels.rds"))
climate  <- readRDS(file.path(derived_dir, "climate_panels.rds"))
if (is.null(conflict) && is.null(prices) && is.null(climate))
  skip_step("06: no exposure data yet; add the conflict, food price or climate files and run 03 to 05 first")

# Window start and end (CMC) for each child; lived months only
win <- list()
for (w in names(windows)) {
  a <- ch$birth_cmc + windows[[w]][1]
  b <- ch$birth_cmc + windows[[w]][2]
  if (w != "pre") b <- pmin(b, ch$interview_cmc)
  win[[w]] <- list(a = a, b = b, len = diff(windows[[w]]) + 1)
}
win$post <- list(a = ch$interview_cmc + post_window[1], b = ch$interview_cmc + post_window[2],
                 len = diff(post_window) + 1)
win$`1000` <- list(a = ch$birth_cmc + windows$preg[1],
                   b = pmin(ch$birth_cmc + windows$y2[2], ch$interview_cmc),
                   len = windows$y2[2] - windows$preg[1] + 1)

# Add one exposure for every window: transform(sum over window, window length)
add_exposure <- function(M, name, transform) {
  if (is.null(M)) return(invisible())
  P <- prep_cum(M); row <- match(ch$cl_key, rownames(M))
  for (w in names(win)) {
    s <- window_sum(P, row, win[[w]]$a, win[[w]]$b)
    ch[[paste0(name, "_", w)]] <<- transform(s, win[[w]]$len)
  }
}
log_count <- function(s, len) log1p(s)
mean_per  <- function(s, len) s / len
count     <- function(s, len) s

if (!is.null(conflict)) {
  cp <- conflict$panels
  add_exposure(cp[[paste0("deaths_r", conflict_radius_km)]], "conf", log_count)
  for (r in conflict_radius_sens) add_exposure(cp[[paste0("deaths_r", r)]], paste0("conf", r), log_count)
}
if (!is.null(prices)) add_exposure(prices$panels$food_infl / 10, "food", mean_per)
if (!is.null(climate)) {
  cl <- climate$panels
  add_exposure(cl$heat,     "heat",   mean_per)
  add_exposure(cl$temp,     "temp",   mean_per)
  add_exposure(cl$drought,  "dry",    mean_per)
  add_exposure(cl$wet_m,    "wet",    count)
  add_exposure(cl$hotdry_m, "hotdry", count)
  if (!is.null(cl$flood_m)) add_exposure(cl$flood_m, "flood", count)
}

expo <- grep("_(preg|y1|y2|pre|post|1000)$", names(ch), value = TRUE)
cat("\nExposure coverage (share of children with a value), by round:\n")
cov <- sapply(expo, function(v) tapply(is.finite(ch[[v]]), ch$round, mean))
print(round(t(cov), 2))

saveRDS(ch, file.path(derived_dir, "children_exposures.rds"))
step_done(sprintf("06 done: %d exposure variables for %d children", length(expo), nrow(ch)))
