# =============================================================================
#  04_exposure_prices.R: monthly food inflation at the nearest market
#  Source: World Bank Real-Time Food Prices (RTFP), Nigeria market-level CSV.
#  For every market: year-on-year change in the log food price index.
#  Each cluster takes the series of its nearest market. Clusters farther than
#  prices_max_km from any market (most of North Central, which has no monitored
#  market) take their zone's average series from the same file instead.
#  Writes derived/price_panels.rds
# =============================================================================

source("00_settings.R"); source("00_functions.R")
children <- readRDS(file.path(derived_dir, "children.rds"))
clusters <- readRDS(file.path(derived_dir, "clusters.rds"))
cmc_range <- c(min(children$birth_cmc) + windows$pre[1],
               max(children$interview_cmc) + post_window[2])

build_prices <- function() {
  if (is.null(prices_file) || !file.exists(prices_file)) {
    message("Food price file not found (", prices_file, "); food price exposures skipped")
    return(NULL)
  }
  rt <- read.csv(prices_file)
  if ("ISO3" %in% names(rt)) rt <- rt[rt$ISO3 == "NGA", ]
  # The file also carries summary series (national, zone and overall market
  # averages); keep real markets only
  zone_rows <- if ("adm1_name" %in% names(rt)) rt[rt$adm1_name %in% "Geopolitical Zone", ] else rt[0, ]
  if ("adm1_name" %in% names(rt))
    rt <- rt[!rt$adm1_name %in% c("National Average", "Geopolitical Zone", "Market Average"), ]
  rt <- rt[!is.na(rt$lat) & !is.na(rt$lon), ]
  need <- c("mkt_name", "lat", "lon", prices_index)
  if (length(setdiff(need, names(rt)))) stop("CHECK price file columns; missing: ",
                                             paste(setdiff(need, names(rt)), collapse = ", "))
  rt$cmc <- if ("DATES" %in% names(rt)) cmc_from_text(rt$DATES) else cmc_from_ym(rt$year, rt$month)
  rt$market <- paste(rt$mkt_name, round(rt$lat, 4), round(rt$lon, 4))
  rt <- rt[!is.na(rt[[prices_index]]) & rt[[prices_index]] > 0, ]

  mk <- unique(rt[, c("market", "mkt_name", "lat", "lon")])
  full <- c(cmc_range[1] - 12, cmc_range[2])                 # 12 extra months for the lag
  lp <- make_panel(rt$market, rt$cmc, log(rt[[prices_index]]), mk$market, full, fill = NA)
  infl <- 100 * (lp[, -(1:12), drop = FALSE] - lp[, 1:(ncol(lp) - 12), drop = FALSE])
  colnames(infl) <- colnames(lp)[-(1:12)]                    # % change on 12 months earlier

  # Nearest market to each cluster
  nearest <- t(sapply(seq_len(nrow(clusters)), function(i) {
    dkm <- haversine_km(clusters$lat[i], clusters$lon[i], mk$lat, mk$lon)
    c(which.min(dkm), min(dkm))
  }))
  M <- infl[nearest[, 1], , drop = FALSE]
  rownames(M) <- clusters$cl_key
  near_tab <- data.frame(cl_key = clusters$cl_key, market = mk$mkt_name[nearest[, 1]],
                         dist_km = nearest[, 2], source = "nearest market")

  # Zone averages for clusters far from every market
  zone_of_cl <- unique(children[, c("cl_key", "zone")])
  near_tab$zone <- zone_of_cl$zone[match(near_tab$cl_key, zone_of_cl$cl_key)]
  if (!is.null(prices_max_km) && nrow(zone_rows)) {
    zone_rows$cmc <- if ("DATES" %in% names(zone_rows)) cmc_from_text(zone_rows$DATES) else
      cmc_from_ym(zone_rows$year, zone_rows$month)
    zone_rows <- zone_rows[!is.na(zone_rows[[prices_index]]) & zone_rows[[prices_index]] > 0, ]
    zl <- make_panel(zone_rows$mkt_name, zone_rows$cmc, log(zone_rows[[prices_index]]),
                     unique(zone_rows$mkt_name), full, fill = NA)
    zi <- 100 * (zl[, -(1:12), drop = FALSE] - zl[, 1:(ncol(zl) - 12), drop = FALSE])
    colnames(zi) <- colnames(zl)[-(1:12)]
    far <- which(near_tab$dist_km > prices_max_km & near_tab$zone %in% rownames(zi))
    if (length(far)) {
      M[far, ] <- zi[near_tab$zone[far], , drop = FALSE]
      near_tab$source[far] <- "zone average"
    }
    cat("Clusters using their zone's average series (more than", prices_max_km, "km from a market):\n")
    print(table(near_tab$zone, near_tab$source))
  }
  cat(sprintf("Food prices: %d markets; distance from cluster to nearest market (km): median %.0f, 90th percentile %.0f, max %.0f\n",
              nrow(mk), median(near_tab$dist_km), quantile(near_tab$dist_km, 0.9), max(near_tab$dist_km)))
  list(panels = list(food_infl = M), nearest = near_tab, cmc_range = cmc_range)
}

prices <- build_prices()
saveRDS(prices, file.path(derived_dir, "price_panels.rds"))
step_done("04 done: food price exposures")
