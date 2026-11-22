# =============================================================================
#  13_food_markets_bounds.R: two added analyses, specified before running
#  (plan Doc, "Part 3", 22 November 2026). Both are reported whatever they show.
#
#  A. Food inflation near markets. The median cluster is about 89 km from its
#     nearest monitored market, so many children's inflation values are
#     approximations, and error of that kind pulls estimates toward zero.
#     Primary: the main model (all shocks, HAZ, same fixed effects and weights)
#     refit on children in clusters within 50 km of a market, with the placebo
#     before conception refit on the same children.
#     Secondary: 25 km and 100 km (a threshold with fewer than 100 clusters is
#     reported but not interpreted); the full sample with each food inflation
#     window interacted with being within 50 km; stunting (linear probability
#     model) at 50 km.
#  B. What the data rule out. For each shock and window in the main model: the
#     minimum detectable effect (80% power, two-sided 5%: 2.80 x SE), per unit
#     and per SD of exposure, and the most harmful effect inside the 95% CI.
#     For each domain's contribution to the 2018 to 2023-24 change (step 08):
#     the most harmful contribution inside the 95% CI and the minimum
#     detectable contribution, for HAZ and for stunting.
#
#  Reads derived/children_exposures.rds, price_panels.rds, main_fits.rds and
#  output/09_attribution.csv (run 07 and 08 first).
#  Writes output/13a_food_near_markets.csv, 13b_detectable_effects.csv,
#  13c_contribution_bounds.csv
# =============================================================================

source("00_settings.R"); source("00_functions.R")
f_exp <- file.path(derived_dir, "children_exposures.rds")
f_pr  <- file.path(derived_dir, "price_panels.rds")
f_fit <- file.path(derived_dir, "main_fits.rds")
f_att <- file.path(out_dir, "09_attribution.csv")
for (f in c(f_exp, f_pr, f_fit, f_att))
  if (!file.exists(f)) skip_step(paste0("13: ", basename(f), " not found; run steps 04 to 08 first"))

ch     <- readRDS(f_exp)
prices <- readRDS(f_pr)
fits   <- readRDS(f_fit)
if (is.null(prices)) skip_step("13: no food price exposure; nothing to test")

# Same design as the main model in step 07
ch$round_zone_age  <- paste(ch$round, ch$zone, ch$age_m)
ch$round_zone_bmoy <- paste(ch$round, ch$zone, ch$birth_moy)
fe_main <- c("cl_key", "round_zone_age", "round_zone_bmoy")
x_main  <- fits$x_main
x_pre   <- paste0(unlist(fits$dom), "_pre")
food_terms <- grep("^food_", x_main, value = TRUE)
if (!length(food_terms)) skip_step("13: food inflation is not in the main model")

ch$mkt_km <- prices$nearest$dist_km[match(ch$cl_key, prices$nearest$cl_key)]
ch$mkt_source <- prices$nearest$source[match(ch$cl_key, prices$nearest$cl_key)]
z95 <- qnorm(0.975); z_mde <- qnorm(0.975) + qnorm(0.80)   # 2.80


# ---- A. Food inflation near markets ------------------------------------------------------
cat("\nClusters by distance to the nearest monitored market:\n")
cl_tab <- unique(ch[, c("cl_key", "round", "zone", "mkt_km")])
cl_tab$band <- cut(cl_tab$mkt_km, c(-Inf, 25, 50, 100, Inf),
                   labels = c("within 25 km", "25 to 50 km", "50 to 100 km", "over 100 km"))
print(table(cl_tab$zone, cl_tab$band))
print(table(cl_tab$round, cl_tab$band))

food_rows <- function(fit, label, km, terms) {
  if (is.null(fit)) return(NULL)
  t <- tidy_fe(fit, label)
  t <- t[t$term %in% terms, ]
  if (!nrow(t)) return(NULL)
  t$max_km <- km
  t
}
fit_safe <- function(d, y, x) tryCatch(fe_reg(d, y, x, fe_main),
                                       error = function(e) { message("  ", conditionMessage(e)); NULL })

res_a <- list()
for (km in c(Inf, 100, 50, 25)) {
  d <- if (is.finite(km)) ch[!is.na(ch$mkt_km) & ch$mkt_km <= km, ] else ch
  tag <- if (is.finite(km)) paste0("within ", km, " km") else "all clusters"
  cat("\n", tag, ": ", length(unique(d$cl_key)), " clusters, ", nrow(d), " children\n", sep = "")
  res_a[[length(res_a) + 1]] <- food_rows(fit_safe(d, "haz", x_main), paste("HAZ,", tag), km, food_terms)
  res_a[[length(res_a) + 1]] <- food_rows(fit_safe(d, "haz", c(x_main, x_pre)),
                                          paste("HAZ placebo before conception,", tag), km, "food_pre")
  if (km %in% c(Inf, 50))
    res_a[[length(res_a) + 1]] <- food_rows(fit_safe(d, "stunted", x_main),
                                            paste("Stunting (linear probability),", tag), km, food_terms)
}

# Interaction: does the food inflation association differ within 50 km?
ch$near50 <- as.numeric(ch$mkt_km <= 50)
x_int <- paste0(food_terms, "_x_near50")
for (k in seq_along(food_terms)) ch[[x_int[k]]] <- ch[[food_terms[k]]] * ch$near50
res_a[[length(res_a) + 1]] <- food_rows(fit_safe(ch, "haz", c(x_main, x_int)),
                                        "HAZ, all clusters, interaction with within 50 km", 50,
                                        c(food_terms, x_int))

tab_a <- do.call(rbind, res_a)
tab_a$interpret <- ifelse(tab_a$clusters >= 100, "yes", "no (fewer than 100 clusters)")
write.csv(tab_a, file.path(out_dir, "13a_food_near_markets.csv"), row.names = FALSE)
cat("\nA. Food inflation (per 10 points) by distance to a market:\n")
print(tab_a[, c("model", "term", "estimate", "lower", "upper", "p_value", "n", "clusters")],
      digits = 3, row.names = FALSE)


# ---- B1. Detectable effects and bounds for each shock and window -------------------------
wsd <- function(x, w) { ok <- is.finite(x) & is.finite(w); m <- sum(w[ok] * x[ok]) / sum(w[ok])
                        sqrt(sum(w[ok] * (x[ok] - m)^2) / sum(w[ok])) }
bounds_for <- function(fit, y, scale, harm_low) {
  if (is.null(fit)) return(NULL)
  x <- names(fit$coef)
  d <- ch[complete.cases(ch[, c(y, x, "wt")]), ]
  sdx <- sapply(x, function(v) wsd(d[[v]], d$wt))
  est <- fit$coef * scale; se <- fit$se * scale
  data.frame(outcome = y, term = x, sd_exposure = sdx,
             estimate = est, lower = est - z95 * se, upper = est + z95 * se,
             mde_per_unit = z_mde * se, mde_per_sd = z_mde * se * sdx,
             most_harmful_per_unit = if (harm_low) est - z95 * se else est + z95 * se,
             most_harmful_per_sd = (if (harm_low) est - z95 * se else est + z95 * se) * sdx,
             row.names = NULL)
}
tab_b <- rbind(bounds_for(fits$haz, "haz", 1, harm_low = TRUE),
               bounds_for(fits$stunted, "stunted", 100, harm_low = FALSE))  # stunting in points
write.csv(tab_b, file.path(out_dir, "13b_detectable_effects.csv"), row.names = FALSE)
cat("\nB1. Main model: minimum detectable effect (80% power) and most harmful effect in the 95% CI\n",
    "    HAZ in SD; stunting in percentage points; per unit and per SD of exposure\n", sep = "")
print(tab_b[, c("outcome", "term", "estimate", "mde_per_unit", "mde_per_sd", "most_harmful_per_sd")],
      digits = 2, row.names = FALSE)


# ---- B2. Bounds on each domain's contribution to the change ------------------------------
att <- read.csv(f_att)
att$scale <- ifelse(grepl("^Stunting", att$outcome), 100, 1)    # stunting in percentage points
tab_c <- data.frame(outcome = att$outcome, domain = att$domain,
                    contribution = att$contribution * att$scale,
                    lower = att$lower * att$scale, upper = att$upper * att$scale,
                    observed_change = att$observed_change * att$scale,
                    mde_contribution = z_mde * att$se * att$scale)
tab_c$most_harmful <- ifelse(att$scale == 100, tab_c$upper, tab_c$lower)
write.csv(tab_c, file.path(out_dir, "13c_contribution_bounds.csv"), row.names = FALSE)
cat("\nB2. Contribution of each domain to the change, 2018 to 2023-24\n",
    "    HAZ in SD; stunting in percentage points\n", sep = "")
print(tab_c, digits = 3, row.names = FALSE)

step_done("13 done: food inflation near markets and detectable effects saved")
