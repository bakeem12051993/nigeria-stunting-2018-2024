# =============================================================================
#  00_functions.R: helper functions used by the other scripts.
#  Base R only (plus haven and survey), so nothing here needs extra packages.
# =============================================================================

suppressPackageStartupMessages({
  library(haven)
  library(survey)
})


# ---- Reading DHS files -------------------------------------------------------

# Read only the needed columns of a DHS Stata file
read_dhs <- function(path, vars) {
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
  df <- haven::read_dta(path, col_select = tidyselect::any_of(vars))
  miss <- setdiff(vars, names(df))
  if (length(miss)) message("  Not in ", basename(path), ": ", paste(miss, collapse = ", "))
  df
}

# A column as plain numbers (NA if the column is absent)
getv <- function(df, v) {
  if (!v %in% names(df)) return(rep(NA_real_, nrow(df)))
  as.numeric(haven::zap_labels(df[[v]]))
}

# A column's value labels as text (NA if the column is absent)
lab_of <- function(df, v) {
  if (!v %in% names(df)) return(rep(NA_character_, nrow(df)))
  x <- df[[v]]
  if (inherits(x, "haven_labelled")) as.character(haven::as_factor(x, levels = "labels"))
  else as.character(x)
}

# Normalize place names so they match across rounds ("FCT-Abuja" = "fct abuja")
key_of <- function(x) gsub("[^a-z]", "", tolower(x))

zone_display <- c(northcentral = "North Central", northeast = "North East",
                  northwest = "North West", southeast = "South East",
                  southsouth = "South South", southwest = "South West")

# State names to keys that match across surveys ("FCT Abuja" -> "fctabuja")
state_key <- function(x) {
  s <- key_of(x)
  s[s %in% c("fct", "abuja", "fctabuja")] <- "fctabuja"
  s[s %in% c("abroad", "")] <- NA
  s
}
# Geopolitical zone of each state (for surveys that record state only, e.g. the MIS)
zone_of_state <- c(
  benue = "northcentral", kogi = "northcentral", kwara = "northcentral", nasarawa = "northcentral",
  niger = "northcentral", plateau = "northcentral", fctabuja = "northcentral",
  adamawa = "northeast", bauchi = "northeast", borno = "northeast", gombe = "northeast",
  taraba = "northeast", yobe = "northeast",
  jigawa = "northwest", kaduna = "northwest", kano = "northwest", katsina = "northwest",
  kebbi = "northwest", sokoto = "northwest", zamfara = "northwest",
  abia = "southeast", anambra = "southeast", ebonyi = "southeast", enugu = "southeast", imo = "southeast",
  akwaibom = "southsouth", bayelsa = "southsouth", crossriver = "southsouth", delta = "southsouth",
  edo = "southsouth", rivers = "southsouth",
  ekiti = "southwest", lagos = "southwest", ogun = "southwest", ondo = "southwest",
  osun = "southwest", oyo = "southwest")
# Zone display name from a zone or a state label
zone_from <- function(x) {
  k <- key_of(x)
  k <- ifelse(k %in% names(zone_display), k, zone_of_state[state_key(x)])
  unname(zone_display[k])
}

# Weight factor that makes measured children match all children in each cell
# (post-stratification). Cells with children but none measured keep factor 1.
match_factor <- function(wt, measured, all, cell) {
  pop <- tapply(wt[all], cell[all], sum)
  got <- tapply(wt[measured], cell[measured], sum)
  lost <- setdiff(names(pop), names(got))
  if (length(lost)) message("  ", length(lost), " weighting cell(s) have no measured children: ",
                            paste(lost, collapse = ", "))
  f <- as.numeric(pop[cell] / got[cell])
  f[is.na(f)] <- 1
  f
}
weighting_label <- c(published       = "Published DHS weights",
                     zone            = "Matched to all children by zone",
                     state_residence = "Matched to all children by state and urban/rural")

# Print value labels so codes can be checked against the recode manual
show_codes <- function(df, vars, title) {
  cat("\n================ Value labels:", title, "================\n")
  for (v in intersect(vars, names(df))) {
    cat("\n--", v, ":", attr(df[[v]], "label"), "\n")
    lab <- attr(df[[v]], "labels")
    if (is.null(lab)) cat("   (no value labels)\n") else print(lab)
  }
}

# 1/0 indicator with NA for missing codes
yes_no <- function(x, yes, na_codes = c(7, 8, 9, 97, 98, 99)) {
  out <- ifelse(x %in% yes, 1, 0)
  out[is.na(x) | x %in% na_codes] <- NA
  out
}


# ---- Dates: DHS century month codes (CMC) --------------------------------------
# CMC = (year - 1900) * 12 + month. All monthly data are indexed this way.

cmc_from_ym   <- function(year, month) (year - 1900) * 12 + month
cmc_from_text <- function(x) {            # "YYYY-MM" or "YYYY-MM-DD"
  cmc_from_ym(as.integer(substr(x, 1, 4)), as.integer(substr(x, 6, 7)))
}
cmc_from_date <- function(d) cmc_from_ym(as.integer(format(d, "%Y")), as.integer(format(d, "%m")))
cmc_year      <- function(cmc) 1900 + (cmc - 1) %/% 12
cmc_month     <- function(cmc) (cmc - 1) %% 12 + 1
cmc_label     <- function(cmc) sprintf("%04d-%02d", cmc_year(cmc), cmc_month(cmc))


# ---- Geography -------------------------------------------------------------------

# Great-circle distance in km between one point and many points
haversine_km <- function(lat1, lon1, lat2, lon2) {
  rad <- pi / 180
  dlat <- (lat2 - lat1) * rad; dlon <- (lon2 - lon1) * rad
  a <- sin(dlat / 2)^2 + cos(lat1 * rad) * cos(lat2 * rad) * sin(dlon / 2)^2
  6371 * 2 * atan2(sqrt(a), sqrt(1 - a))
}


# ---- Cluster-month panels and exposure windows -------------------------------------
# A panel is a matrix: one row per cluster ("round_v001"), one column per month (CMC).

# Turn a long table (key, cmc, value) into a cluster x month matrix.
# Values falling in the same cell are summed; empty cells get `fill`.
make_panel <- function(key, cmc, value, keys, cmc_range, fill = 0) {
  months <- seq(cmc_range[1], cmc_range[2])
  M <- matrix(fill, nrow = length(keys), ncol = length(months),
              dimnames = list(keys, as.character(months)))
  ok <- key %in% keys & !is.na(cmc) & cmc >= cmc_range[1] & cmc <= cmc_range[2] & !is.na(value)
  if (!any(ok)) return(M)
  lin <- match(key[ok], keys) + (cmc[ok] - cmc_range[1]) * nrow(M)   # column-major cell index
  s <- tapply(value[ok], lin, sum)
  M[as.numeric(names(s))] <- s
  M
}

# Cumulative sums of a panel, prepared once and reused for every window
prep_cum <- function(M) {
  M0 <- M; M0[is.na(M0)] <- 0
  list(S  = cbind(0, t(apply(M0, 1, cumsum))),
       Sn = cbind(0, t(apply(is.na(M), 1, cumsum))),
       first = as.integer(colnames(M)[1]), last = as.integer(colnames(M)[ncol(M)]))
}

# Sum of a panel over each child's months [a, b] (CMC), for the child's cluster.
# A window not yet reached (b < a) sums to 0. Months outside the panel, or NA
# months, make the result NA.
window_sum <- function(P, row, a, b) {
  S <- P$S; Sn <- P$Sn; first <- P$first; last <- P$last
  out <- rep(NA_real_, length(row))
  empty <- !is.na(a) & !is.na(b) & b < a                  # window not reached yet
  out[empty] <- 0
  ok <- !is.na(row) & !is.na(a) & !is.na(b) & b >= a & a >= first & b <= last
  ia <- a[ok] - first + 1; ib <- b[ok] - first + 1
  s  <- S[cbind(row[ok], ib + 1)]  - S[cbind(row[ok], ia)]
  sn <- Sn[cbind(row[ok], ib + 1)] - Sn[cbind(row[ok], ia)]
  s[sn > 0] <- NA
  out[ok] <- s
  out
}


# ---- Fixed-effects regression (base R) ---------------------------------------------
# Weighted least squares with several sets of fixed effects absorbed by
# alternating projections, and cluster-robust (CR1) standard errors.
# Gives the same coefficients as lm() with dummy variables for each effect.

demean_fe <- function(Z, fe, w, tol = 1e-9, maxit = 5000) {
  g <- lapply(fe, function(f) as.integer(factor(f)))
  for (it in seq_len(maxit)) {
    change <- 0
    for (k in seq_along(g)) {
      m <- rowsum(Z * w, g[[k]], reorder = TRUE) / as.vector(rowsum(w, g[[k]], reorder = TRUE))
      adj <- m[g[[k]], , drop = FALSE]
      Z <- Z - adj
      change <- max(change, max(abs(adj)))
    }
    if (change < tol) return(structure(Z, iterations = it))
  }
  warning("Fixed-effect demeaning did not fully converge")
  structure(Z, iterations = maxit)
}

fe_reg <- function(data, y, x, fe, w = "wt", cluster = "psu", drop_singletons = TRUE) {
  vars <- unique(c(y, x, fe, w, cluster))
  d <- data[complete.cases(data[, vars]), vars, drop = FALSE]
  num <- vapply(d, is.numeric, logical(1))
  if (any(num)) d <- d[Reduce(`&`, lapply(d[num], is.finite)), , drop = FALSE]   # no Inf
  d <- d[d[[w]] > 0, , drop = FALSE]
  if (drop_singletons) {                       # groups of one carry no information
    repeat {
      keep <- rep(TRUE, nrow(d))
      for (f in fe) { n <- table(d[[f]]); keep <- keep & n[as.character(d[[f]])] > 1 }
      if (all(keep)) break
      d <- d[keep, , drop = FALSE]
    }
  }
  Z <- demean_fe(as.matrix(d[, c(y, x)]), d[fe], d[[w]])
  yd <- Z[, 1]; Xd <- Z[, -1, drop = FALSE]; wt <- d[[w]]
  # Drop regressors with no variation left after the fixed effects
  keepx <- apply(Xd, 2, function(v) sum(abs(v)) > 1e-8)
  if (any(!keepx)) warning("Dropped (no variation within fixed effects): ",
                           paste(colnames(Xd)[!keepx], collapse = ", "))
  Xd <- Xd[, keepx, drop = FALSE]
  XtWX <- crossprod(Xd, Xd * wt)
  q <- qr(XtWX)
  if (q$rank < ncol(Xd)) {
    bad <- colnames(Xd)[q$pivot[(q$rank + 1):ncol(Xd)]]
    warning("Dropped (collinear): ", paste(bad, collapse = ", "))
    Xd <- Xd[, setdiff(colnames(Xd), bad), drop = FALSE]
    XtWX <- crossprod(Xd, Xd * wt)
  }
  beta <- solve(XtWX, crossprod(Xd, yd * wt))[, 1]
  e <- yd - as.vector(Xd %*% beta)
  bread <- solve(XtWX)
  cl <- d[[cluster]]
  meat <- crossprod(rowsum(Xd * (wt * e), cl))
  G <- length(unique(cl)); N <- nrow(d); K <- ncol(Xd)
  V <- bread %*% meat %*% bread * (G / (G - 1)) * ((N - 1) / (N - K))
  dimnames(V) <- list(names(beta), names(beta))
  list(coef = beta, vcov = V, se = sqrt(diag(V)), n = N, clusters = G,
       y = y, fe = fe, iterations = attr(Z, "iterations"))
}

tidy_fe <- function(fit, model) {
  est <- fit$coef; se <- fit$se
  data.frame(model = model, term = names(est), estimate = est, se = se,
             lower = est - 1.96 * se, upper = est + 1.96 * se,
             p_value = 2 * pnorm(-abs(est / se)), n = fit$n, clusters = fit$clusters,
             row.names = NULL)
}

# Linear combination of coefficients (estimate and delta-method SE)
lincom <- function(fit, weights) {
  b <- fit$coef[names(weights)]; V <- fit$vcov[names(weights), names(weights), drop = FALSE]
  est <- sum(weights * b); se <- sqrt(as.numeric(t(weights) %*% V %*% weights))
  c(estimate = est, se = se, lower = est - 1.96 * se, upper = est + 1.96 * se)
}


# ---- Survey summaries ------------------------------------------------------------------

est_by <- function(var, by, design) {
  x  <- svyby(as.formula(paste0("~", var)),
              as.formula(paste0("~", paste(by, collapse = " + "))),
              design, svymean, na.rm = TRUE)
  ci <- confint(x)
  out <- as.data.frame(x)[, by, drop = FALSE]
  out$outcome  <- var
  out$estimate <- as.numeric(coef(x))
  out$lower    <- ci[, 1]
  out$upper    <- ci[, 2]
  rownames(out) <- NULL
  out
}

tidy_rep <- function(x, label) {
  est <- as.numeric(coef(x)); se <- as.numeric(SE(x))
  data.frame(analysis = label, term = names(coef(x)), estimate = est, se = se,
             lower = est - 1.96 * se, upper = est + 1.96 * se, row.names = NULL)
}

wmean <- function(x, w) { ok <- !is.na(x) & !is.na(w); sum(x[ok] * w[ok]) / sum(w[ok]) }

step_done <- function(msg) cat(sprintf("\n[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

# ---- Steps that cannot run yet --------------------------------------------------
# skip_step() stops the current script with a clear reason. run_all.R reports it
# as "Skipped" and moves on to the next step.
skip_step <- function(msg) stop(structure(class = c("skip_step", "error", "condition"),
                                          list(message = msg, call = NULL)))

# ---- Which DHS rounds are on disk? ------------------------------------------------
dhs_path <- function(yr, f) file.path(dhs_dir, dhs_files[[yr]][[f]])

available_rounds <- function(report = FALSE) {
  found <- sapply(names(dhs_files), function(yr)
    sapply(names(dhs_files[[yr]]), function(f) file.exists(dhs_path(yr, f))))
  colnames(found) <- round_label[colnames(found)]
  if (report) {
    cat("Working folder:", getwd(), "\n")
    cat("DHS files in", dhs_dir, ":\n")
    for (yr in names(dhs_files)) for (f in names(dhs_files[[yr]]))
      cat(sprintf("  %-8s %-3s %-8s %s\n", round_label[[yr]], f,
                  if (file.exists(dhs_path(yr, f))) "found" else "MISSING", dhs_files[[yr]][[f]]))
  }
  ok <- names(dhs_files)[colSums(!found) == 0]
  partial <- names(dhs_files)[colSums(found) > 0 & colSums(!found) > 0]
  if (length(partial))
    message("Some files are missing for ", paste(round_label[partial], collapse = ", "),
            ", so that round is left out. Check its folder and file names in dhs_files (00_settings.R).")
  if (!length(ok))
    stop("No complete DHS round found in ", normalizePath(dhs_dir, mustWork = FALSE),
         ". Check project_dir and dhs_dir in 00_settings.R.", call. = FALSE)
  if (report && length(ok) < length(dhs_files))
    message("Running with ", paste(round_label[ok], collapse = " and "), " only. ",
            "Steps that compare rounds (decomposition, attribution) wait until every round is in place.")
  ok
}

# ---- ERA5-Land files: which variable is in which .nc file? ------------------------
nc_vars <- function(f) {
  v <- tryCatch(names(terra::sds(f)), error = function(e) NULL)          # several variables
  if (is.null(v)) v <- tryCatch(terra::varnames(terra::rast(f)), error = function(e) character(0))
  unique(v)
}
# Dates of the layers in a NetCDF file. Climate Data Store files keep them in
# "valid_time", which terra may not read; ncdf4 can.
nc_valid_time <- function(f) {
  if (!requireNamespace("ncdf4", quietly = TRUE)) return(NULL)
  nc <- ncdf4::nc_open(f); on.exit(ncdf4::nc_close(nc))
  for (v in c("valid_time", "time")) {
    if (!v %in% c(names(nc$dim), names(nc$var))) next
    x <- ncdf4::ncvar_get(nc, v)
    u <- ncdf4::ncatt_get(nc, v, "units")$value
    if (!is.character(u)) next
    step <- c(seconds = 1, hours = 3600, days = 86400)[sub(" since.*", "", u)]
    if (is.na(step)) next
    return(as.POSIXct(x * step, origin = sub(".* since ", "", u), tz = "UTC"))
  }
  NULL
}

find_era5 <- function(dir, need = c("t2m", "d2m", "tp")) {
  out <- list()
  for (f in list.files(dir, pattern = "\\.nc$", full.names = TRUE, ignore.case = TRUE)) {
    vs <- nc_vars(f)
    for (v in intersect(need, vs)) {
      if (!is.null(out[[v]])) message("  ", v, " is in more than one file; using ", basename(out[[v]]$file))
      else out[[v]] <- list(file = f, var = v, several = length(vs) > 1)
    }
  }
  out
}
