# =============================================================================
#  07_models.R: effects of shocks in each window on child height-for-age
#
#  Design: survey-weighted regressions with fixed effects for
#    - DHS cluster (compares children in the same community),
#    - round x zone x age in months (removes everything common to children of
#      the same age in the same survey and zone, including the zone's own age
#      pattern of growth and the zone-wide path of prices and conflict), and
#    - round x zone x calendar month of birth (birth seasonality by zone).
#  So effects are identified only from local differences in the timing and size
#  of shocks. Standard errors are clustered by DHS cluster.
#  With a national age pattern instead (round x age), conflict before
#  conception "predicted" height-for-age (a failed placebo), because growth
#  falls faster with age in the North, where conflict also rose over time.
#  That design is kept as a sensitivity analysis (model H).
#
#  Models
#    A  each shock domain alone           E  compound shocks (interactions)
#    B  all domains together (main)        F  household buffers (livestock,
#    C  siblings (mother fixed effects)       farmland, rural residence)
#    D  placebo windows                    G  climate effects with and without
#    H  sensitivity: 10 and 50 km             food prices (pathway)
#       conflict radius; air temperature
#  Writes output/08_models.csv, output/fig_window_effects.png,
#  derived/main_fits.rds
# =============================================================================

source("00_settings.R"); source("00_functions.R")
f_exp <- file.path(derived_dir, "children_exposures.rds")
if (!file.exists(f_exp)) skip_step("07: no exposure file yet; run 03 to 06 first")
ch <- readRDS(f_exp)

ch$round_zone_age      <- paste(ch$round, ch$zone, ch$age_m)
ch$round_zone_bmoy     <- paste(ch$round, ch$zone, ch$birth_moy)
ch$round_state_ageband <- paste(ch$round, ch$state, cut(ch$age_m, seq(-1, 59, 6)))
fe_main   <- c("cl_key", "round_zone_age", "round_zone_bmoy")
fe_mother <- c("mother", "round_zone_age", "round_zone_bmoy")
main_win  <- c("preg", "y1", "y2")

# Domains with data (each needs at least 80% coverage in both rounds)
has <- function(v) v %in% names(ch) && all(tapply(is.finite(ch[[v]]), ch$round, mean) >= 0.8)
dom <- list(conflict = "conf", food = "food", climate = c("heat", "dry", "wet"))
dom <- lapply(dom, function(p) p[sapply(p, function(x) has(paste0(x, "_y1")))])
dom <- dom[lengths(dom) > 0]
if (!length(dom)) skip_step("07: no exposure reaches 80% coverage in every round; check the coverage table from 06")
cat("Domains with data:", paste(names(dom), collapse = ", "), "\n")
terms_for <- function(prefixes, wins = main_win) as.vector(outer(prefixes, wins, paste, sep = "_"))
x_main <- terms_for(unlist(dom))

results <- list()
fit_and_store <- function(label, y, x, fe = fe_main, data = ch) {
  fit <- tryCatch(fe_reg(data, y, x, fe), error = function(e) { message(label, ": ", conditionMessage(e)); NULL })
  if (!is.null(fit)) results[[label]] <<- tidy_fe(fit, label)
  fit
}

# A. Each domain alone
for (d in names(dom)) fit_and_store(paste0("A. ", d, " only (HAZ)"), "haz", terms_for(dom[[d]]))

# B. All domains together: the main model
fit_haz  <- fit_and_store("B. All shocks (HAZ)", "haz", x_main)
fit_stnt <- fit_and_store("B. All shocks (stunting, linear probability)", "stunted", x_main)

# C. Siblings: mother fixed effects
fit_and_store("C. All shocks, mother fixed effects (HAZ)", "haz", x_main, fe = fe_mother)

# D. Placebo: shocks before conception and after the interview should not matter.
#    Before conception is the informative test. Children in the same cluster were
#    interviewed in the same weeks, so their after-interview exposure barely
#    differs and the cluster fixed effect absorbs most of it (wide intervals).
x_pre  <- terms_for(unlist(dom), "pre")
x_post <- terms_for(unlist(dom), "post")
fit_and_store("D. Placebo: before conception (HAZ)", "haz", c(x_main, x_pre))
if (all(sapply(x_post, has))) {
  fit_and_store("D. Placebo: after the interview (HAZ)", "haz", c(x_main, x_post))
} else message("Placebo after the interview skipped: exposure data do not reach 12 months past the interviews")

# E. Compound shocks: first-1,000-days totals and their interactions
k1000 <- paste0(unlist(dom), "_1000")
pairs <- list(c("food_1000", "conf_1000"), c("food_1000", "dry_1000"),
              c("heat_1000", "dry_1000"), c("conf_1000", "dry_1000"))
pairs <- Filter(function(p) all(p %in% k1000), pairs)
for (p in pairs) ch[[paste(p, collapse = "_x_")]] <- ch[[p[1]]] * ch[[p[2]]]
x_comp <- c(k1000, sapply(pairs, paste, collapse = "_x_"))
if ("hotdry_1000" %in% names(ch) && has("hotdry_1000")) x_comp <- c(x_comp, "hotdry_1000")
fit_and_store("E. Compound shocks, first 1,000 days (HAZ)", "haz", x_comp)

# F. Household buffers: does a shock hurt less with livestock, farmland, or in towns?
ch$rural_hh <- as.numeric(ch$residence == "rural")
# Rural residence is fixed within a cluster, so the cluster fixed effect absorbs
# its main effect; only its interactions with the shocks are estimated.
for (m in c("livestock", "farmland", "rural_hh")) {
  inter <- paste0(k1000, "_x_", m)
  for (k in k1000) ch[[paste0(k, "_x_", m)]] <- ch[[k]] * ch[[m]]
  main_m <- if (m == "rural_hh") character(0) else m
  fit_and_store(paste0("F. Buffer: ", m, " (HAZ)"), "haz", c(k1000, main_m, inter))
}

# G. Pathway: climate effects with and without local food prices
if (all(c("climate", "food") %in% names(dom))) {
  clim1000 <- paste0(dom$climate, "_1000")
  others <- setdiff(k1000, c(clim1000, "food_1000"))
  fit_and_store("G. Climate without food prices (HAZ)", "haz", c(clim1000, others))
  fit_and_store("G. Climate with food prices (HAZ)",    "haz", c(clim1000, others, "food_1000"))
}

# H. Sensitivity: conflict radius and air temperature instead of humid heat
for (r in conflict_radius_sens) {
  alt <- sub("^conf_", paste0("conf", r, "_"), x_main)
  if (all(sapply(alt, has))) fit_and_store(paste0("H. Conflict within ", r, " km (HAZ)"), "haz", alt)
}
fit_and_store("H. State-specific age pattern, 6-month bands (HAZ)", "haz", x_main,
              fe = c("cl_key", "round_state_ageband", "round_zone_bmoy"))
fit_and_store("H. National age pattern, original design (HAZ)", "haz", x_main,
              fe = c("cl_key", "round_age", "round_bmoy"))
if ("heat" %in% unlist(dom) && has("temp_y1"))
  fit_and_store("H. Air temperature instead of humid heat (HAZ)", "haz", sub("^heat_", "temp_", x_main))

models <- do.call(rbind, results)
write.csv(models, file.path(out_dir, "08_models.csv"), row.names = FALSE)
print(models[grepl("^B\\. All shocks \\(HAZ", models$model), c("term", "estimate", "lower", "upper", "p_value")], digits = 3)
saveRDS(list(haz = fit_haz, stunted = fit_stnt, x_main = x_main, dom = dom),
        file.path(derived_dir, "main_fits.rds"))


# ---- Figure: effect of each shock by window (main model) ---------------------------
b <- results[["B. All shocks (HAZ)"]]
if (!is.null(b)) {
  b$prefix <- sub("_(preg|y1|y2)$", "", b$term)
  b$window <- factor(sub("^.*_", "", b$term), levels = main_win)
  labels <- c(conf = "Conflict (log deaths + 1)", food = "Food inflation (per 10 pp)",
              heat = "Humid heat (SD)", dry = "Drought (SD)", wet = "Very wet months")
  b <- b[order(match(b$prefix, names(labels)), b$window), ]
  png(file.path(out_dir, "fig_window_effects.png"), width = 2000, height = 1300, res = 220)
  par(mar = c(4.5, 13, 2.5, 1))
  y <- rev(seq_len(nrow(b)))
  plot(b$estimate, y, xlim = range(c(b$lower, b$upper, 0)), yaxt = "n", pch = 19,
       col = c("#7f7f7f", "#4f81bd", "#1f3d7a")[as.integer(b$window)],
       xlab = "Change in HAZ (SD units) per unit of exposure", ylab = "",
       main = "Shock effects by window: pregnancy, year 1, year 2")
  arrows(b$lower, y, b$upper, y, angle = 90, code = 3, length = 0.03,
         col = c("#7f7f7f", "#4f81bd", "#1f3d7a")[as.integer(b$window)])
  abline(v = 0, col = "grey40")
  axis(2, at = y, labels = paste(labels[b$prefix], c("pregnancy", "year 1", "year 2")[as.integer(b$window)], sep = ", "),
       las = 1, cex.axis = 0.75)
  dev.off()
}
step_done("07 done: models saved")
