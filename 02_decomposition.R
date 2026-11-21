# =============================================================================
#  02_decomposition.R: the descriptive base of the paper
#  1. Pooled asset index (comparable across rounds)
#  2. Change in stunting and HAZ, nationally and by zone, residence and age
#  3. Shift-share: composition versus within-zone change
#  3b. Sampling-frame check: the change with the household distribution across
#      states held fixed (the 2018 survey used the 2006 census frame, the
#      2023-24 survey the new frame prepared for the planned census)
#  4. Oaxaca-Blinder decomposition of the change in mean HAZ and stunting
#  5. The key results again with the published DHS weights (sensitivity)
#  Reads derived/children.rds and the household recodes. Writes output/03 to 07b tables,
#  and adds the asset index to derived/children.rds.
# =============================================================================

source("00_settings.R"); source("00_functions.R")
ana <- readRDS(file.path(derived_dir, "children.rds"))
cat("Weights:", weighting_label[[weighting]], "\n")


# ---- 1. Pooled asset index -------------------------------------------------------
# DHS wealth quintiles are relative within each survey, so they cannot explain
# change between surveys. One principal component across both rounds can.
asset_vars <- c("electricity", "radio", "tv", "fridge", "bicycle", "motorcycle",
                "car", "finished_floor", "finished_wall", "finished_roof", "clean_fuel")
ok  <- complete.cases(ana[, asset_vars])
pca <- prcomp(ana[ok, asset_vars], center = TRUE, scale. = TRUE)
sgn   <- sign(cor(pca$x[, 1], ana$electricity[ok]))   # higher score = more assets
score <- sgn * pca$x[, 1]
ana$asset_score <- NA_real_
ana$asset_score[ok] <- score
ana$asset_q <- cut(ana$asset_score, quantile(ana$asset_score, 0:5 / 5, na.rm = TRUE),
                   include.lowest = TRUE, labels = paste0("Q", 1:5))
cat("\nAsset index loadings (first component):\n"); print(round(sgn * pca$rotation[, 1], 3))
saveRDS(ana, file.path(derived_dir, "children.rds"))

ana$north    <- grepl("^north", ana$zone_key)
ana$south    <- grepl("^south", ana$zone_key)
ana$urban    <- ana$residence %in% "urban"
ana$rural    <- ana$residence %in% "rural"
ana$age24up  <- ana$age_m >= 24
ana$under36  <- ana$age_m < 36
ana$dob_full <- ana$b10 %in% codes$dob_complete

ana$zone_f   <- factor(ana$zone)

# Survey design. With corrected weighting, the published-weight design is
# post-stratified to all children under 5 (totals saved by step 01), so the
# standard errors reflect the correction.
ps_totals <- readRDS(file.path(derived_dir, "ps_totals.rds"))
calibrate_to_all <- function(d) if (is.null(ps_totals)) d else
  postStratify(d, ~ps_cell, ps_totals, partial = TRUE)
des_dhs <- svydesign(ids = ~psu, strata = ~strata, weights = ~wt_dhs, data = ana, nest = TRUE)
des     <- calibrate_to_all(des_dhs)


# ---- 2. Descriptive results, by round ----------------------------------------------------
desc_nat  <- rbind(est_by("stunted", "round", des), est_by("haz", "round", des))
desc_zone <- est_by("stunted", c("round", "zone"), des)
desc_res  <- est_by("stunted", c("round", "residence"), des)
desc_age  <- est_by("stunted", c("round", "age_band"), des)
# Where children live: weighted share of measured children in each zone, by round.
# A shift between rounds feeds the composition part of the shift-share below.
zs <- svyby(~zone_f, ~round, des, svymean, na.rm = TRUE)
zci <- confint(zs)
zone_shares <- data.frame(round = sub(":.*", "", names(coef(zs))),
                          zone  = sub(".*:zone_f", "", names(coef(zs))),
                          share = as.numeric(coef(zs)), lower = zci[, 1], upper = zci[, 2],
                          row.names = NULL)
print(desc_nat, digits = 3)
cat("\nShare of children by zone (%):\n")
print(round(100 * xtabs(share ~ zone + round, zone_shares), 1))
write.csv(zone_shares, file.path(out_dir, "03f_zone_shares.csv"), row.names = FALSE)
write.csv(desc_nat,  file.path(out_dir, "03a_change_national.csv"),       row.names = FALSE)
write.csv(desc_zone, file.path(out_dir, "03b_stunting_by_zone.csv"),      row.names = FALSE)
write.csv(desc_res,  file.path(out_dir, "03c_stunting_by_residence.csv"), row.names = FALSE)
write.csv(desc_age,  file.path(out_dir, "03d_stunting_by_age.csv"),       row.names = FALSE)

# Everything below compares the two rounds
if (length(unique(ana$round)) < 2)
  skip_step("02: descriptive tables saved for one round; the change, shift-share and Oaxaca-Blinder parts run once the 2018 files are added")
rdes_dhs <- as.svrepdesign(des_dhs, type = "subbootstrap", replicates = n_boot)
rdes     <- calibrate_to_all(rdes_dhs)

change_in <- function(design, label) {
  fit <- svyglm(stunted ~ round, design = design)
  ci  <- confint(fit)["round2023-24", ]
  data.frame(group = label, change_pp = 100 * coef(fit)[["round2023-24"]],
             lower = 100 * ci[1], upper = 100 * ci[2], row.names = NULL)
}
change <- rbind(change_in(des, "Nigeria"),
                do.call(rbind, lapply(sort(unique(ana$zone)), function(z)
                  change_in(subset(des, zone == z), z))))
print(change, digits = 3)


# ---- 3. Shift-share --------------------------------------------------------------------
shift_share <- function(w, data, y = "stunted", g = "zone") {
  ok <- !is.na(data[[y]]) & !is.na(data[[g]])
  w <- w[ok]; yv <- data[[y]][ok]; gv <- factor(data[[g]][ok]); rv <- data$round[ok]
  share <- function(k) { s <- tapply(w[rv == k], gv[rv == k], sum); s[is.na(s)] <- 0; s / sum(s) }
  prev  <- function(k) tapply(w[rv == k] * yv[rv == k], gv[rv == k], sum) /
                       tapply(w[rv == k], gv[rv == k], sum)
  s1 <- share("2018"); s2 <- share("2023-24"); p1 <- prev("2018"); p2 <- prev("2023-24")
  both <- !is.na(p1) & !is.na(p2)
  c(total       = sum(s2 * p2, na.rm = TRUE) - sum(s1 * p1, na.rm = TRUE),
    composition = sum(((s2 - s1) * (p1 + p2) / 2)[both]),
    within      = sum(((s1 + s2) / 2 * (p2 - p1))[both]))
}
ss <- rbind(
  tidy_rep(withReplicates(rdes, shift_share, y = "stunted", g = "zone"),  "Stunting, by zone"),
  tidy_rep(withReplicates(rdes, shift_share, y = "stunted", g = "state"), "Stunting, by state"),
  tidy_rep(withReplicates(rdes, shift_share, y = "haz",     g = "zone"),  "HAZ, by zone"))
print(ss, digits = 3)


# ---- 3b. Sampling-frame check ----------------------------------------------------------
# Each survey's weighted household recode gives the share of households in each
# state under that survey's frame. Children under 5 per household = the child
# share divided by the household share. Holding the household distribution fixed
# means multiplying one round's child weights in each state by the ratio of the
# two rounds' household shares; children per household and stunting within each
# state stay as surveyed. The household shares are treated as known, like the
# post-stratification totals.
hh_share <- lapply(c("2018", "2024"), function(yr) {
  h <- read_dhs(dhs_path(yr, "hr"), c("hv005", pr_state_var[[yr]]))
  s <- tapply(getv(h, "hv005") / 1e6, state_key(lab_of(h, pr_state_var[[yr]])), sum)
  s / sum(s)
})
names(hh_share) <- c("2018", "2023-24")
ana$skey <- state_key(ana$state)
miss <- setdiff(unique(ana$skey), intersect(names(hh_share[["2018"]]), names(hh_share[["2023-24"]])))
if (length(miss)) message("Frame check: no household share for ", paste(miss, collapse = ", "),
                          "; their weights are left unchanged")
frame_factor <- function(from, to) {
  f <- hh_share[[to]][ana$skey] / hh_share[[from]][ana$skey]
  f[is.na(f)] <- 1
  as.numeric(f)
}
# hold = the round whose household distribution both rounds take
to_2018 <- ifelse(ana$round == "2023-24", frame_factor("2023-24", "2018"), 1)
to_2024 <- ifelse(ana$round == "2018",    frame_factor("2018", "2023-24"), 1)
frame_change <- function(w, data, y, factor) {
  ok <- !is.na(data[[y]]); w <- w * factor
  m  <- function(k) { i <- ok & data$round == k; sum(w[i] * data[[y]][i]) / sum(w[i]) }
  c(rise = m("2023-24") - m("2018"))
}
frame_tab <- do.call(rbind, lapply(c("stunted", "haz"), function(y) rbind(
  tidy_rep(withReplicates(rdes, frame_change, y = y, factor = rep(1, nrow(ana))),
           paste0(y, ": as surveyed")),
  tidy_rep(withReplicates(rdes, frame_change, y = y, factor = to_2018),
           paste0(y, ": both rounds on the 2018 household distribution")),
  tidy_rep(withReplicates(rdes, frame_change, y = y, factor = to_2024),
           paste0(y, ": both rounds on the 2023-24 household distribution")))))
cat("\nNorth's share of households (%): 2018", round(100 * sum(hh_share[["2018"]][names(zone_of_state)[grepl("^north", zone_of_state)]], na.rm = TRUE), 1),
    "| 2023-24", round(100 * sum(hh_share[["2023-24"]][names(zone_of_state)[grepl("^north", zone_of_state)]], na.rm = TRUE), 1), "\n")
cat("Change with the household distribution held fixed:\n"); print(frame_tab, digits = 3)


# ---- 4. Oaxaca-Blinder -------------------------------------------------------------------
domain_of <- c(
  age_band = "Child", female = "Child", birth_order = "Child",
  short_interval = "Child", multiple = "Child", small_at_birth = "Child",
  mother_age = "Mother", mother_educ = "Mother", mother_short = "Mother",
  mother_thin = "Mother", asset_q = "Household assets",
  water_unimproved = "WASH", sanitation = "WASH", facility_birth = "Health care",
  diarrhea = "Recent illness", fever = "Recent illness",
  livestock = "Livestock and farmland", farmland = "Livestock and farmland",
  residence = "Place", zone = "Place")

# Terms observed for at least 90% of children in each round. Place of delivery
# and size at birth were asked only under about 36 months in 2023-24, so they
# enter only the under-36-months sensitivity analysis.
usable <- sapply(names(domain_of), function(v)
  all(tapply(!is.na(ana[[v]]), ana$round, mean) >= 0.9))
if (any(!usable)) message("Left out of the decomposition (mostly missing in a round): ",
                          paste(names(domain_of)[!usable], collapse = ", "))
model_terms <- names(domain_of)[usable]

ob_decompose <- function(w, data, y, terms, reference = "pooled", keep = NULL) {
  sel <- if (is.null(keep)) rep(TRUE, nrow(data)) else data[[keep]] %in% TRUE
  cc  <- sel & complete.cases(data[, c(y, terms)])
  d <- data[cc, , drop = FALSE]; w <- w[cc]
  X <- model.matrix(as.formula(paste("~", paste(terms, collapse = " + "))), d)
  r <- d$round == "2023-24"; yv <- d[[y]]
  wm <- function(idx) colSums(X[idx, , drop = FALSE] * w[idx]) / sum(w[idx])
  xbar18 <- wm(!r); xbar24 <- wm(r)
  ybar18 <- sum(yv[!r] * w[!r]) / sum(w[!r]); ybar24 <- sum(yv[r] * w[r]) / sum(w[r])
  beta <- switch(reference,
    pooled    = lm.wfit(cbind(X, round_2023_24 = as.numeric(r)), yv, w)$coefficients[colnames(X)],
    "2018"    = lm.wfit(X[!r, , drop = FALSE], yv[!r], w[!r])$coefficients,
    "2023-24" = lm.wfit(X[r, , drop = FALSE],  yv[r],  w[r])$coefficients)
  beta[is.na(beta)] <- 0
  contrib  <- (xbar24 - xbar18) * beta
  term_col <- c("(Intercept)", terms)[attr(X, "assign") + 1]
  dom_col  <- factor(domain_of[term_col], levels = unique(domain_of[terms]))
  by_dom   <- tapply(contrib, dom_col, sum); by_dom[is.na(by_dom)] <- 0
  total <- ybar24 - ybar18; explained <- sum(contrib)
  c(n_children = sum(cc), total = total, explained = explained,
    unexplained = total - explained,
    setNames(by_dom, paste0("explained: ", names(by_dom))))
}
run_ob <- function(y, label, terms = model_terms, design = rdes, ...)
  tidy_rep(withReplicates(design, ob_decompose, y = y, terms = terms, ...), label)

ob_main <- rbind(run_ob("haz", "HAZ, all children"),
                 run_ob("stunted", "Stunting (linear probability), all children"))
ob_sub  <- rbind(run_ob("haz", "HAZ, North", keep = "north"),
                 run_ob("haz", "HAZ, South", keep = "south"),
                 run_ob("haz", "HAZ, urban", keep = "urban"),
                 run_ob("haz", "HAZ, rural", keep = "rural"))
ob_sens <- rbind(
  run_ob("haz", "HAZ, 2018 coefficients as reference",    reference = "2018"),
  run_ob("haz", "HAZ, 2023-24 coefficients as reference", reference = "2023-24"),
  run_ob("haz", "HAZ, children 24 to 59 months",          keep = "age24up"),
  run_ob("haz", "HAZ, complete date of birth only",       keep = "dob_full"),
  run_ob("haz", "HAZ, under 36 months, adding facility birth and size at birth",
         terms = union(model_terms, c("facility_birth", "small_at_birth")), keep = "under36"))
print(ob_main, digits = 3)


# ---- 5. Sensitivity: published DHS weights ----------------------------------------------
# The main analysis uses the weights chosen in 00_settings.R; this repeats the key
# results with the published DHS weights.
if (weighting != "published") {
  des_pub  <- des_dhs
  rdes_pub <- rdes_dhs
  ch_pub   <- change_in(des_pub, "Nigeria")
  pub <- rbind(
    data.frame(analysis = "Change in stunting (percentage points)", term = "change",
               estimate = ch_pub$change_pp, se = NA, lower = ch_pub$lower, upper = ch_pub$upper),
    tidy_rep(withReplicates(rdes_pub, shift_share, y = "stunted", g = "zone"), "Stunting, by zone"),
    run_ob("haz",     "HAZ, all children",                           design = rdes_pub),
    run_ob("stunted", "Stunting (linear probability), all children", design = rdes_pub))
  pub$weights <- weighting_label[["published"]]
  write.csv(pub, file.path(out_dir, "07b_published_weights.csv"), row.names = FALSE)
}


# ---- Save ---------------------------------------------------------------------------------
write.csv(change,    file.path(out_dir, "03e_change_by_zone.csv"),       row.names = FALSE)
write.csv(ss,        file.path(out_dir, "04_shift_share.csv"),           row.names = FALSE)
write.csv(frame_tab, file.path(out_dir, "04b_frame_constant.csv"),       row.names = FALSE)
write.csv(ob_main,   file.path(out_dir, "05_oaxaca_main.csv"),           row.names = FALSE)
write.csv(ob_sub,    file.path(out_dir, "06_oaxaca_subgroups.csv"),      row.names = FALSE)
write.csv(ob_sens,   file.path(out_dir, "07_oaxaca_sensitivity.csv"),    row.names = FALSE)
step_done("02 done: descriptive change and decomposition saved")
