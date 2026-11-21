# =============================================================================
#  05_exposure_climate.R: monthly climate anomalies at each DHS cluster
#  Source: ERA5-Land monthly averaged reanalysis (2m temperature, 2m dewpoint,
#  total precipitation), 1991 onward. Optional: cluster-month flood flags.
#  Builds, for every cluster and month:
#    heat     humid-heat (wet-bulb) anomaly, in SD from the 1991-2020 normal
#    temp     air temperature anomaly (SD), for a sensitivity analysis
#    precip   rainfall anomaly (SD)
#    drought  SPEI-3 (SPEI package if installed; otherwise a standardized
#             3-month rainfall index, which ignores evaporation)
#    flags    hot, drought, extremely wet and hot-and-dry months
#  Needs the terra package. Writes derived/climate_panels.rds
# =============================================================================

source("00_settings.R"); source("00_functions.R")
children <- readRDS(file.path(derived_dir, "children.rds"))
clusters <- readRDS(file.path(derived_dir, "clusters.rds"))
cmc_range <- c(min(children$birth_cmc) + windows$pre[1],
               max(children$interview_cmc) + post_window[2])

# Read one ERA5-Land variable and return a cluster x month matrix
read_era5 <- function(src) {
  f <- src$file
  r <- if (src$several) terra::rast(f, subds = src$var) else terra::rast(f)
  tm <- tryCatch(terra::time(r), error = function(e) NULL)
  if (is.null(tm) || all(is.na(tm))) tm <- nc_valid_time(f)
  if (!is.null(tm) && !all(is.na(tm)) && length(tm) == terra::nlyr(r)) {
    cmc <- cmc_from_date(as.Date(tm))
  } else {
    message("  No dates readable in ", basename(f), " (install ncdf4); assuming consecutive months from ",
            era5_first_month)
    cmc <- cmc_from_text(era5_first_month) + seq_len(terra::nlyr(r)) - 1
  }
  pts <- terra::vect(data.frame(lon = clusters$lon, lat = clusters$lat),
                     geom = c("lon", "lat"), crs = "EPSG:4326")
  v <- as.matrix(terra::extract(r, pts, ID = FALSE))
  # ERA5-Land covers land only: move coastal clusters to the nearest land cell
  miss <- which(is.na(v[, 1]))
  if (length(miss)) {
    cells <- terra::as.data.frame(r[[1]], xy = TRUE, cells = TRUE, na.rm = TRUE)
    for (i in miss) {
      k <- which.min(haversine_km(clusters$lat[i], clusters$lon[i], cells$y, cells$x))
      v[i, ] <- as.numeric(unlist(r[cells$cell[k]]))
    }
    message("  ", length(miss), " clusters moved to the nearest land cell (", src$var, ")")
  }
  dimnames(v) <- list(clusters$cl_key, as.character(cmc))
  v
}

# Standardized anomaly against each cluster's calendar-month normal
anomaly <- function(M) {
  cmc <- as.integer(colnames(M)); moy <- cmc_month(cmc)
  base <- cmc_year(cmc) %in% baseline_years
  Z <- M
  for (m in 1:12) {
    j <- moy == m; jb <- j & base
    mu <- rowMeans(M[, jb, drop = FALSE], na.rm = TRUE)
    s  <- apply(M[, jb, drop = FALSE], 1, sd, na.rm = TRUE)
    Z[, j] <- (M[, j, drop = FALSE] - mu) / s
  }
  Z
}

# Anomaly against each cluster's own linear trend, by calendar month: a "shock"
# is an unusual month for that place and time, not the slow warming trend.
# Without this, heat in the year before conception "predicted" growth (a failed
# placebo), because later-born children face a warmer climate in every window.
anomaly_detrended <- function(M) {
  cmc <- as.integer(colnames(M)); moy <- cmc_month(cmc); yr <- cmc_year(cmc)
  Z <- M
  for (m in 1:12) {
    j <- which(moy == m)
    x <- yr[j] - mean(yr[j])
    Y <- M[, j, drop = FALSE]
    mu <- rowMeans(Y)
    b  <- as.vector((Y - mu) %*% x) / sum(x^2)          # slope for each cluster
    R  <- Y - mu - outer(b, x)
    s  <- sqrt(rowSums(R^2) / (length(x) - 2))
    Z[, j] <- R / s
  }
  Z
}

# Relative humidity (Magnus) and wet-bulb temperature (Stull 2011), both in deg C
rel_humidity <- function(t, td) 100 * exp(17.625 * td / (243.04 + td)) / exp(17.625 * t / (243.04 + t))
wet_bulb <- function(t, rh) {
  t * atan(0.151977 * sqrt(rh + 8.313659)) + atan(t + rh) - atan(rh - 1.676331) +
    0.00391838 * rh^1.5 * atan(0.023101 * rh) - 4.686035
}

# Drought index: SPEI-3 if the SPEI package is installed, else a standardized
# 3-month rainfall index
# Standardized 3-month sum of a monthly series, by calendar month
standardized_3m <- function(M) {
  M3 <- M
  for (j in seq_len(ncol(M))) M3[, j] <- if (j < 3) NA else rowSums(M[, (j - 2):j, drop = FALSE])
  anomaly(M3)
}

# SPEI can return no value, or an infinite one, for extreme months that fall
# outside its fitted distribution. Those cells take the standardized 3-month
# water balance instead, and all values are capped at -3 and 3.
repair_spei <- function(S, balance) {
  Zb  <- standardized_3m(balance)
  bad <- !is.finite(S) & is.finite(Zb)
  S[bad] <- Zb[bad]
  S[!is.finite(S)] <- NA
  S <- pmin(pmax(S, -3), 3)
  attr(S, "n_replaced") <- sum(bad)
  S
}

drought_index <- function(P, Tm) {
  cmc <- as.integer(colnames(P)); start <- c(cmc_year(cmc[1]), cmc_month(cmc[1]))
  ref <- list(c(min(baseline_years), 1), c(max(baseline_years), 12))
  if (requireNamespace("SPEI", quietly = TRUE)) {
    quiet <- function(expr) {                 # SPEI prints progress messages
      out <- NULL
      invisible(capture.output(out <- suppressMessages(suppressWarnings(expr))))
      out
    }
    res <- lapply(seq_len(nrow(P)), function(i) {
      pet <- as.numeric(quiet(SPEI::thornthwaite(ts(Tm[i, ], frequency = 12, start = start),
                                                 lat = clusters$lat[i])))
      bal <- P[i, ] - pet
      fit <- quiet(SPEI::spei(ts(bal, frequency = 12, start = start), scale = 3,
                              ref.start = ref[[1]], ref.end = ref[[2]]))
      list(spei = as.numeric(fit$fitted), bal = bal)
    })
    S <- do.call(rbind, lapply(res, `[[`, "spei"))
    B <- do.call(rbind, lapply(res, `[[`, "bal"))
    dimnames(S) <- dimnames(B) <- dimnames(P)
    out <- repair_spei(S, B)
    attr(out, "method") <- sprintf(paste0("SPEI-3 (Thornthwaite evapotranspiration); %d cluster-months ",
                                          "outside SPEI's fitted range took the standardized water balance"),
                                   attr(out, "n_replaced"))
  } else {
    out <- pmin(pmax(standardized_3m(P), -3), 3)
    attr(out, "method") <- "Standardized 3-month rainfall (SPI-type); install SPEI for SPEI-3"
  }
  dimnames(out) <- dimnames(P)
  out
}

build_climate <- function() {
  era5_files <- find_era5(era5_dir)
  miss <- setdiff(c("t2m", "d2m", "tp"), names(era5_files))
  if (length(miss)) {
    message("ERA5-Land variables not found in ", era5_dir, ": ", paste(miss, collapse = ", "),
            "; climate exposures skipped")
    return(NULL)
  }
  for (v in names(era5_files)) cat("  ", v, "from", basename(era5_files[[v]]$file), "\n")
  t2m <- read_era5(era5_files$t2m) - 273.15                 # K to deg C
  d2m <- read_era5(era5_files$d2m) - 273.15
  tp  <- read_era5(era5_files$tp)                           # m per day (monthly mean)
  cmc <- as.integer(colnames(tp))
  days <- as.integer(format(as.Date(paste(cmc_year(cmc) + (cmc_month(cmc) == 12),
                                          cmc_month(cmc) %% 12 + 1, 1, sep = "-")) - 1, "%d"))
  precip <- sweep(tp * 1000, 2, days, `*`)                  # mm per month
  tw <- wet_bulb(t2m, pmin(pmax(rel_humidity(t2m, d2m), 1), 100))

  drought <- drought_index(precip, t2m)
  cat("Drought index:", attr(drought, "method"), "\n")
  temp_anomaly <- if (isTRUE(heat_detrend)) anomaly_detrended else anomaly
  panels <- list(heat = temp_anomaly(tw), temp = temp_anomaly(t2m), precip = anomaly(precip),
                 drought = -drought)                         # positive = drier
  panels$hot_m     <- (panels$heat > heat_threshold_sd) * 1
  panels$drought_m <- (drought < drought_threshold) * 1
  panels$wet_m     <- (panels$precip > wet_threshold_sd) * 1
  panels$hotdry_m  <- panels$hot_m * panels$drought_m

  if (!is.null(flood_file) && file.exists(flood_file)) {
    fl <- read.csv(flood_file)
    key <- paste(ifelse(fl$round == "2018", "2018", "2024"), fl$v001, sep = "_")
    panels$flood_m <- make_panel(key, cmc_from_text(fl$month), fl$flooded,
                                 clusters$cl_key, range(cmc), fill = 0)
  }
  # Keep only the months the analysis needs
  keep <- as.character(seq(cmc_range[1], cmc_range[2]))
  panels <- lapply(panels, function(M) {
    out <- matrix(NA_real_, nrow(M), length(keep), dimnames = list(rownames(M), keep))
    both <- intersect(keep, colnames(M)); out[, both] <- M[, both]; out
  })
  list(panels = panels, cmc_range = cmc_range, drought_method = attr(drought, "method"))
}

climate <- build_climate()
saveRDS(climate, file.path(derived_dir, "climate_panels.rds"))
step_done("05 done: climate exposures")
