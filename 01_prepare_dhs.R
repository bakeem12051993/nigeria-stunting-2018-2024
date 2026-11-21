# =============================================================================
#  01_prepare_dhs.R
#  1. Replicate the official stunting estimates (household member recode)
#  2. Build a harmonized children's file for each round (children's recode,
#     with livestock and farmland from the household recode)
#  3. Attach cluster GPS locations
#  4. Data quality checks
#  Writes: derived/children.rds, derived/clusters.rds,
#          output/01_official_replication.csv (published and corrected weights),
#          output/01b_zone_mix_all_vs_measured.csv, output/02_data_quality.csv
# =============================================================================

source("00_settings.R"); source("00_functions.R")
rounds <- available_rounds(report = TRUE)   # rounds with all four files on disk


# ---- 1. Official stunting estimates, and corrected for where children live --------------
# DHS method: de facto children (hv103 = 1) with a valid hc70, published weights.
# Targets: 37% (2018) and 40% (2023-24).
# Children were measured only in households picked for the men's interview, and the
# measured children over-represent the South (see 09_survey_checks.R). The corrected
# rows weight measured children to the zone, or state and urban/rural, mix of all
# de facto children under 5 in the same survey.

pr_children <- function(yr) {
  pr <- read_dhs(dhs_path(yr, "pr"), c("hv005", "hv021", "hv022", "hv025", "hv103", "hv105",
                                       "hc70", pr_zone_var[[yr]], pr_state_var[[yr]]))
  d <- data.frame(round = round_label[[yr]], wt = getv(pr, "hv005") / 1e6,
                  psu = paste(yr, getv(pr, "hv021")), strata = paste(yr, getv(pr, "hv022")),
                  residence = getv(pr, "hv025"), defacto = getv(pr, "hv103"),
                  age_y = getv(pr, "hv105"), hc70 = getv(pr, "hc70"),
                  zone = zone_from(lab_of(pr, pr_zone_var[[yr]])),
                  state = state_key(lab_of(pr, pr_state_var[[yr]])))
  d$measured <- d$defacto %in% 1 & !is.na(d$hc70) & d$hc70 < 9996
  d$all_u5   <- d$defacto %in% 1 & d$age_y %in% 0:4
  d <- d[d$measured | d$all_u5, ]
  d$stunted  <- ifelse(d$measured, as.numeric(d$hc70 < -200), NA)
  d
}
pr <- do.call(rbind, lapply(rounds, pr_children))

stunting_by_weighting <- function(k) {
  pr$cell <- switch(k, published = "all", zone = paste(pr$round, pr$zone),
                    state_residence = paste(pr$round, pr$state, pr$residence))
  m <- pr[pr$measured, ]
  des <- svydesign(ids = ~psu, strata = ~strata, weights = ~wt, data = m, nest = TRUE)
  if (k != "published") {                       # post-stratify to all children under 5
    pop <- tapply(pr$wt[pr$all_u5], pr$cell[pr$all_u5], sum)
    des <- postStratify(des, ~cell, data.frame(cell = names(pop), Freq = as.numeric(pop)),
                        partial = TRUE)
  }
  r <- svyby(~stunted, ~round, des, svymean); ci <- confint(r)
  out <- data.frame(weighting = weighting_label[[k]], round = r$round,
                    children = as.numeric(table(m$round)[r$round]),
                    stunted_pct = 100 * r$stunted, lower = 100 * ci[, 1], upper = 100 * ci[, 2],
                    row.names = NULL)
  if (length(unique(m$round)) == 2) {
    g <- svyglm(stunted ~ round, design = des)
    b <- coef(g)[[2]]; se <- SE(g)[[2]]
    out$change_pp    <- 100 * b
    out$change_lower <- 100 * (b - 1.96 * se)
    out$change_upper <- 100 * (b + 1.96 * se)
  }
  out
}
official <- do.call(rbind, lapply(names(weighting_label), stunting_by_weighting))
print(official, digits = 3)
write.csv(official, file.path(out_dir, "01_official_replication.csv"), row.names = FALSE)

# Where children live: all de facto children under 5 versus measured children
share <- function(keep, label) {
  t <- tapply(pr$wt[keep], list(pr$zone[keep], pr$round[keep]), sum)
  t <- 100 * sweep(t, 2, colSums(t, na.rm = TRUE), "/")
  colnames(t) <- paste(label, colnames(t)); t
}
zs_all <- share(pr$all_u5, "all children"); zs_meas <- share(pr$measured, "measured")
zone_mix <- data.frame(zone = rownames(zs_all), zs_all, zs_meas[rownames(zs_all), , drop = FALSE],
                       check.names = FALSE, row.names = NULL)
cat("\nShare of children by zone (%), published weights: all under-5s vs measured\n")
print(zone_mix, digits = 3)
write.csv(zone_mix, file.path(out_dir, "01b_zone_mix_all_vs_measured.csv"), row.names = FALSE)


# ---- 2. Harmonized children's file ---------------------------------------------------

kr_vars <- c("caseid", "v001", "v002", "v005", "v008", "v012", "v021", "v022",
             "v024", "v025", "v106", "v113", "v116", "v119", "v120", "v121",
             "v122", "v123", "v124", "v125", "v127", "v128", "v129", "v160",
             "v161", "v438", "v445", "b0", "b3", "b4", "b5", "b9", "b10", "b11",
             "bord", "hw13", "hw15", "hw70", "m15", "m18", "h11", "h22",
             "szone", "sstate", "sstate1", "v042")
code_check_vars <- c("v024", "szone", "v025", "v106", "v113", "v116", "v127",
                     "v161", "b10", "hw13", "hw15", "m15", "m18", "h11")

prep_kr <- function(yr) {
  kr <- read_dhs(dhs_path(yr, "kr"), kr_vars)
  if (show_code_labels) show_codes(kr, code_check_vars, paste("KR", round_label[[yr]]))
  g <- function(v) getv(kr, v)

  d <- data.frame(
    year    = as.integer(yr),
    round   = round_label[[yr]],
    v001    = g("v001"), v002 = g("v002"),
    cl_key  = paste(yr, g("v001"), sep = "_"),       # cluster id, unique across rounds
    mother  = paste(yr, as.character(kr$caseid)),    # mother id, for sibling models
    psu     = paste(yr, g("v021")),
    strata  = paste(yr, g("v022")),                  # v022 = v023 in recent DHS
    wt      = g("v005") / 1e6,
    alive   = g("b5"), with_mother = g("b9"),
    birth_cmc = g("b3"), interview_cmc = g("v008"),
    hw70 = g("hw70"), hw13 = g("hw13"), hw15 = g("hw15"), b10 = g("b10"),
    v042 = g("v042")                                 # household selected for measuring
  )
  d$age_m <- d$interview_cmc - d$birth_cmc
  d$birth_moy <- cmc_month(d$birth_cmc)              # calendar month of birth

  # Outcome (9996-9998 implausible or flagged, 9999 missing)
  d$haz     <- ifelse(!is.na(d$hw70) & d$hw70 < 9996, d$hw70 / 100, NA)
  d$stunted <- as.numeric(d$haz < -2)

  # Child
  d$age_band <- cut(d$age_m, c(-1, 5, 11, 23, 35, 47, 59),
                    labels = c("0-5", "6-11", "12-23", "24-35", "36-47", "48-59"))
  d$female <- yes_no(g("b4"), 2)
  bord <- g("bord"); b11 <- g("b11")
  d$birth_order <- cut(bord, c(0, 1, 3, 5, Inf), labels = c("1", "2-3", "4-5", "6+"))
  d$short_interval <- ifelse(bord == 1, 0, ifelse(is.na(b11), NA, as.numeric(b11 < 24)))
  d$multiple <- as.numeric(g("b0") > 0)
  m18 <- g("m18")
  d$small_at_birth <- ifelse(m18 %in% c(4, 5), 1, ifelse(m18 %in% 1:3, 0, NA))

  # Mother
  d$mother_age  <- cut(g("v012"), c(14, 19, 34, 49), labels = c("15-19", "20-34", "35-49"))
  d$mother_educ <- factor(g("v106"), levels = 0:3,
                          labels = c("none", "primary", "secondary", "higher"))
  v438 <- g("v438"); v445 <- g("v445")
  d$mother_short <- ifelse(v438 < 9990, as.numeric(v438 < 1450), NA)
  d$mother_thin  <- ifelse(v445 < 9990, as.numeric(v445 < 1850), NA)

  # WASH
  v113 <- g("v113")
  d$water_unimproved <- ifelse(v113 %in% codes$water_unimproved, 1, 0)
  d$water_unimproved[is.na(v113) | v113 %in% codes$water_na] <- NA
  v116 <- g("v116"); v160 <- g("v160")
  san <- ifelse(v116 %in% codes$toilet_open, "open defecation",
         ifelse(v116 %in% codes$toilet_improved,
                ifelse(v160 %in% 1, "limited (shared)", "basic"), "unimproved"))
  san[is.na(v116) | v116 %in% codes$toilet_na] <- NA
  d$sanitation <- factor(san, levels = c("basic", "limited (shared)", "unimproved", "open defecation"))

  # Assets (pooled index built in 02)
  d$electricity <- yes_no(g("v119"), 1); d$radio      <- yes_no(g("v120"), 1)
  d$tv          <- yes_no(g("v121"), 1); d$fridge     <- yes_no(g("v122"), 1)
  d$bicycle     <- yes_no(g("v123"), 1); d$motorcycle <- yes_no(g("v124"), 1)
  d$car         <- yes_no(g("v125"), 1)
  fin <- function(x) { o <- as.numeric(x %in% codes$finished); o[is.na(x) | x %in% c(96, 97, 99)] <- NA; o }
  d$finished_floor <- fin(g("v127")); d$finished_wall <- fin(g("v128"))
  d$finished_roof  <- fin(g("v129"))
  v161 <- g("v161")
  d$clean_fuel <- ifelse(v161 %in% codes$clean_fuel[[yr]], 1, 0)
  d$clean_fuel[is.na(v161) | v161 %in% c(95, 96, 97, 99)] <- NA

  # Health care and recent illness
  m15 <- g("m15")
  d$facility_birth <- ifelse(m15 %in% codes$facility_birth, 1,
                             ifelse(m15 %in% c(10:19, 90:96), 0, NA))
  d$diarrhea <- ifelse(g("h11") %in% 1:2, 1, ifelse(g("h11") %in% 0, 0, NA))
  d$fever    <- yes_no(g("h22"), 1)

  # Place
  d$residence <- factor(g("v025"), levels = 1:2, labels = c("urban", "rural"))
  zkey <- key_of(lab_of(kr, zone_var[[yr]]))
  d$zone_key <- zkey
  d$zone  <- ifelse(zkey %in% names(zone_display), zone_display[zkey], zkey)
  d$state <- state_key(lab_of(kr, state_var[[yr]]))

  # Household recode: livestock and farmland
  hr <- read_dhs(dhs_path(yr, "hr"), c("hv001", "hv002", "hv244", "hv246"))
  idx <- match(paste(d$v001, d$v002), paste(getv(hr, "hv001"), getv(hr, "hv002")))
  d$livestock <- yes_no(getv(hr, "hv246"), 1)[idx]
  d$farmland  <- yes_no(getv(hr, "hv244"), 1)[idx]

  # GPS (displaced cluster locations). 0,0 means no coordinates.
  ge <- foreign::read.dbf(dhs_path(yr, "ge"), as.is = TRUE)
  ge <- ge[!(ge$LATNUM == 0 & ge$LONGNUM == 0), ]
  gi <- match(d$v001, ge$DHSCLUST)
  d$lat <- ge$LATNUM[gi]; d$lon <- ge$LONGNUM[gi]
  d
}

kr_all <- do.call(rbind, lapply(rounds, prep_kr))
kr_all$round <- factor(kr_all$round, levels = unname(round_label[rounds]))
kr_all$round_age   <- paste(kr_all$round, kr_all$age_m)      # fixed effects used later
kr_all$round_bmoy  <- paste(kr_all$round, kr_all$birth_moy)

cat("\nChildren by zone and round (labels must match across rounds):\n")
print(table(kr_all$zone, kr_all$round, useNA = "ifany"))
st <- table(kr_all$state, kr_all$round)
if (ncol(st) == 2) {
  unmatched <- rownames(st)[st[, 1] == 0 | st[, 2] == 0]
  if (length(unmatched)) message("CHECK state names present in only one round: ",
                                 paste(unmatched, collapse = ", "))
}


# ---- 3. Data quality checks ---------------------------------------------------------------

# Only households selected for measuring count as eligible. In 2018, children in
# other households carry hw13 = 7; in 2023-24 their hw13 is blank. v042 marks
# selected households in both rounds.
quality_check <- function(d) {
  selected <- if (any(!is.na(d$v042))) d$v042 %in% 1 else !is.na(d$hw13)
  e <- d[d$alive %in% 1 & d$with_mother %in% 0 & d$age_m >= 0 & d$age_m < 60 & selected, ]
  data.frame(
    round                  = as.character(e$round[1]),
    selected_for_measuring = nrow(e),
    pct_measured           = 100 * mean(e$hw13 %in% 0),
    pct_valid_haz          = 100 * mean(!is.na(e$haz)),
    pct_flagged_haz        = 100 * mean(e$hw70 %in% 9996:9998),
    sd_haz                 = sd(e$haz, na.rm = TRUE),
    pct_dob_imputed        = 100 * mean(!(e$b10 %in% codes$dob_complete)),
    pct_under2_standing    = 100 * mean(e$hw15[e$age_m < 24] %in% 2, na.rm = TRUE),
    pct_no_gps             = 100 * mean(is.na(e$lat))
  )
}
quality <- do.call(rbind, lapply(split(kr_all, kr_all$round), quality_check))
print(quality, digits = 3)
write.csv(quality, file.path(out_dir, "02_data_quality.csv"), row.names = FALSE)


# ---- 4. Analysis sample and cluster list ---------------------------------------------------

# Analysis weights. wt_dhs keeps the published weight (sensitivity analysis); wt is
# the weight chosen in 00_settings.R, matching the measured children to all living
# children under 5 who live with their mother in the same survey.
eligible <- kr_all$alive %in% 1 & kr_all$with_mother %in% 0 & kr_all$age_m >= 0 & kr_all$age_m < 60
measured <- eligible & !is.na(kr_all$haz)
kr_all$wt_dhs <- kr_all$wt
cell <- switch(weighting, published = NULL,
               zone = paste(kr_all$round, kr_all$zone),
               state_residence = paste(kr_all$round, kr_all$state, kr_all$residence))
if (is.null(cell) && weighting != "published") stop("weighting must be published, zone or state_residence")
ps_totals <- NULL                                  # population totals for post-stratification
kr_all$ps_cell <- "all"
if (!is.null(cell)) {
  kr_all$ps_cell <- cell
  kr_all$wt <- kr_all$wt_dhs * match_factor(kr_all$wt_dhs, measured, eligible, cell)
  pt <- tapply(kr_all$wt_dhs[eligible], cell[eligible], sum)
  ps_totals <- data.frame(ps_cell = names(pt), Freq = as.numeric(pt))
}
saveRDS(ps_totals, file.path(derived_dir, "ps_totals.rds"))
cat("\nAnalysis weights:", weighting_label[[weighting]], "\n")

children <- kr_all[measured, ]

clusters <- unique(children[!is.na(children$lat),
                            c("cl_key", "year", "round", "v001", "lat", "lon")])
clusters$first_cmc <- tapply(children$interview_cmc, children$cl_key, min)[clusters$cl_key]
clusters$last_cmc  <- tapply(children$interview_cmc, children$cl_key, max)[clusters$cl_key]

saveRDS(children, file.path(derived_dir, "children.rds"))
saveRDS(clusters, file.path(derived_dir, "clusters.rds"))
step_done(sprintf("01 done: %d children in %d clusters with GPS", nrow(children), nrow(clusters)))
