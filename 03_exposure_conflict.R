# =============================================================================
#  03_exposure_conflict.R: monthly conflict exposure around each DHS cluster
#  Source: UCDP Georeferenced Event Dataset (GED), CSV version.
#  For every cluster and month: deaths and events from organized violence
#  within 25 km (and 10 and 50 km for sensitivity analyses).
#  Writes derived/conflict_panels.rds
# =============================================================================

source("00_settings.R"); source("00_functions.R")
children <- readRDS(file.path(derived_dir, "children.rds"))
clusters <- readRDS(file.path(derived_dir, "clusters.rds"))
cmc_range <- c(min(children$birth_cmc) + windows$pre[1],
               max(children$interview_cmc) + post_window[2])

build_conflict <- function() {
  if (is.null(ucdp_file) || !file.exists(ucdp_file)) {
    message("UCDP file not found (", ucdp_file, "); conflict exposures skipped")
    return(NULL)
  }
  # Read only the columns needed (the file is large, mostly source text)
  hdr  <- names(read.csv(ucdp_file, nrows = 1))
  keep <- c("country", "where_prec", "latitude", "longitude", "date_start", "best", "best_est",
            "type_of_violence")
  ged <- read.csv(ucdp_file, colClasses = ifelse(tolower(hdr) %in% keep, NA, "NULL"))
  names(ged) <- tolower(names(ged))
  deaths_col <- intersect(c("best", "best_est"), names(ged))[1]   # name differs by version
  ng <- ged[ged$country == "Nigeria" & ged$where_prec <= ucdp_max_precision &
              !is.na(ged$latitude) & !is.na(ged$longitude), ]
  ng$cmc    <- cmc_from_text(ng$date_start)
  ng$deaths <- ng[[deaths_col]]
  ng <- ng[ng$cmc >= cmc_range[1] & ng$cmc <= cmc_range[2], ]
  cat(sprintf("UCDP: %d events in Nigeria, %s to %s (precision <= %d)\n", nrow(ng),
              cmc_label(min(ng$cmc)), cmc_label(max(ng$cmc)), ucdp_max_precision))

  radii <- sort(unique(c(conflict_radius_km, conflict_radius_sens)))
  near <- lapply(seq_len(nrow(clusters)), function(i) {
    dkm <- haversine_km(clusters$lat[i], clusters$lon[i], ng$latitude, ng$longitude)
    j <- which(dkm <= max(radii))
    if (!length(j)) return(NULL)
    data.frame(cl_key = clusters$cl_key[i], cmc = ng$cmc[j], deaths = ng$deaths[j], dist = dkm[j])
  })
  near <- do.call(rbind, near)

  last_cmc <- cmc_from_text(ucdp_last_month)
  panels <- list()
  for (r in radii) {
    sub <- if (is.null(near)) near else near[near$dist <= r, ]
    for (what in c("deaths", "events")) {
      val <- if (is.null(sub)) numeric(0) else if (what == "deaths") sub$deaths else rep(1, nrow(sub))
      M <- make_panel(sub$cl_key, sub$cmc, val, clusters$cl_key, cmc_range, fill = 0)
      M[, as.integer(colnames(M)) > last_cmc] <- NA   # months the file does not cover
      panels[[paste0(what, "_r", r)]] <- M
    }
  }
  main <- panels[[paste0("deaths_r", conflict_radius_km)]]
  yrs <- as.integer(colnames(main)) >= cmc_from_ym(2015, 1)
  cat(sprintf("Clusters with any conflict death within %d km since 2015: %.1f%%\n",
              conflict_radius_km, 100 * mean(rowSums(main[, yrs, drop = FALSE], na.rm = TRUE) > 0)))
  list(panels = panels, cmc_range = cmc_range, radii = radii)
}

conflict <- build_conflict()
saveRDS(conflict, file.path(derived_dir, "conflict_panels.rds"))
step_done("03 done: conflict exposures")
