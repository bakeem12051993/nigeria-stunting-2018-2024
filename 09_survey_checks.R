# =============================================================================
#  09_survey_checks.R: checks that bear on how the change is read
#  1. Where children under 5 live, by zone: 2018 DHS, 2021 MIS (if present) and
#     2023-24 DHS, from the household member recodes (all de facto children)
#  2. Anemia in children 6 to 59 months in the same surveys
#  3. Households picked for measuring (one in three, those picked for the men's
#     interview) versus the other households: young children, women and births
#     per household, by zone. Under random selection the two should look alike.
#  The 2021 Malaria Indicator Survey did not measure height, so it adds no
#  stunting estimate; it only helps here. Its files are found automatically in
#  data/dhs (NGPR81FL and similar) if present.
#  Writes output/10_zone_shares_three_surveys.csv, output/11_anemia_three_surveys.csv,
#         output/12_measured_households_check.csv
# =============================================================================

source("00_settings.R"); source("00_functions.R")
rounds <- available_rounds()

# Surveys to compare: the DHS rounds on disk, plus a malaria survey if present
mis_pr <- list.files(dhs_dir, pattern = "^NGPR[0-9][0-9]FL\\.dta$", recursive = TRUE, ignore.case = TRUE)
srv <- lapply(rounds, function(yr) list(label = paste(round_label[[yr]], "DHS"),
                                        pr = dhs_path(yr, "pr"), zone = pr_zone_var[[yr]]))
for (f in mis_pr) {
  yr <- getv(read_dhs(file.path(dhs_dir, f), "hv007")[1, ], "hv007")
  srv[[length(srv) + 1]] <- list(label = paste(yr, "MIS"), pr = file.path(dhs_dir, f), zone = "hv024")
}
srv <- srv[order(sapply(srv, function(x) substr(x$label, 1, 4)))]
cat("Surveys compared:", paste(sapply(srv, `[[`, "label"), collapse = ", "), "\n")


# ---- 1 and 2. Zone shares and anemia ----------------------------------------------------
shares <- list(); anemia <- list()
for (s in srv) {
  p <- read_dhs(s$pr, c("hv005", "hv021", "hv022", "hv103", "hv105", "hc1", "hc57", s$zone))
  d <- data.frame(wt = getv(p, "hv005") / 1e6, psu = getv(p, "hv021"), strata = getv(p, "hv022"),
                  defacto = getv(p, "hv103"), age_y = getv(p, "hv105"), age_m = getv(p, "hc1"),
                  hc57 = getv(p, "hc57"), zone = zone_from(lab_of(p, s$zone)))
  if (any(is.na(d$zone))) message("  ", s$label, ": ", sum(is.na(d$zone)), " people with no zone")
  u5 <- d[d$defacto %in% 1 & d$age_y %in% 0:4 & !is.na(d$zone), ]
  u5$zone <- factor(u5$zone, levels = unname(zone_display))
  des <- svydesign(ids = ~psu, strata = ~strata, weights = ~wt, data = u5, nest = TRUE)
  m <- svymean(~zone, des); ci <- confint(m)
  shares[[s$label]] <- data.frame(survey = s$label, zone = sub("^zone", "", names(coef(m))),
                                  share_pct = 100 * coef(m), lower = 100 * ci[, 1],
                                  upper = 100 * ci[, 2], row.names = NULL)
  # Anemia (hc57: 1 severe, 2 moderate, 3 mild, 4 not anemic), children 6-59 months
  a <- u5[u5$age_m %in% 6:59 & u5$hc57 %in% 1:4, ]
  if (nrow(a)) {
    a$any_anemia <- as.numeric(a$hc57 %in% 1:3)
    a$moderate_or_severe <- as.numeric(a$hc57 %in% 1:2)
    dA <- svydesign(ids = ~psu, strata = ~strata, weights = ~wt, data = a, nest = TRUE)
    nat <- svymean(~any_anemia + moderate_or_severe, dA); cn <- confint(nat)
    bz  <- svyby(~any_anemia, ~zone, dA, svymean); cz <- confint(bz)
    anemia[[s$label]] <- rbind(
      data.frame(survey = s$label, group = "Nigeria", measure = names(coef(nat)), children = nrow(a),
                 pct = 100 * coef(nat), lower = 100 * cn[, 1], upper = 100 * cn[, 2], row.names = NULL),
      data.frame(survey = s$label, group = as.character(bz$zone), measure = "any_anemia",
                 children = as.numeric(table(a$zone)[as.character(bz$zone)]),
                 pct = 100 * bz$any_anemia, lower = 100 * cz[, 1], upper = 100 * cz[, 2], row.names = NULL))
  }
}
shares <- do.call(rbind, shares); anemia <- do.call(rbind, anemia)
cat("\nShare of children under 5 living in each zone (%):\n")
print(round(xtabs(share_pct ~ zone + survey, shares), 1))
cat("\nAnemia, children 6-59 months (%):\n")
print(anemia[anemia$group == "Nigeria", ], digits = 3, row.names = FALSE)
write.csv(shares, file.path(out_dir, "10_zone_shares_three_surveys.csv"), row.names = FALSE)
write.csv(anemia, file.path(out_dir, "11_anemia_three_surveys.csv"), row.names = FALSE)


# ---- 3. Households picked for measuring versus the rest ------------------------------------
hh_check <- function(yr) {
  h <- read_dhs(dhs_path(yr, "hr"), c("hv001", "hv002", "hv010", "hv014", "hv042", pr_zone_var[[yr]]))
  k <- read_dhs(dhs_path(yr, "kr"), c("v001", "v002", "v008", "b3", "b5"))
  h <- data.frame(id = paste(getv(h, "hv001"), getv(h, "hv002")), picked = getv(h, "hv042") %in% 1,
                  women = getv(h, "hv010"), under5_listed = getv(h, "hv014"),
                  zone = zone_from(lab_of(h, pr_zone_var[[yr]])))
  u5 <- getv(k, "b5") %in% 1 & (getv(k, "v008") - getv(k, "b3")) < 60
  births <- tapply(u5, paste(getv(k, "v001"), getv(k, "v002")), sum)
  h$births_u5 <- ifelse(is.na(births[h$id]), 0, births[h$id])
  out <- do.call(rbind, lapply(split(h, h$zone), function(z) data.frame(
    round = round_label[[yr]], zone = z$zone[1], households = nrow(z),
    pct_picked = 100 * mean(z$picked),
    under5_listed_picked = mean(z$under5_listed[z$picked]), under5_listed_other = mean(z$under5_listed[!z$picked]),
    births_u5_picked = mean(z$births_u5[z$picked]),         births_u5_other = mean(z$births_u5[!z$picked]),
    women_picked = mean(z$women[z$picked]),                 women_other = mean(z$women[!z$picked]))))
  rownames(out) <- NULL
  out
}
hh <- do.call(rbind, lapply(rounds, hh_check))
cat("\nHouseholds picked for measuring versus the rest (means per household, unweighted).\n",
    "Under random selection, picked and other should be close in every zone:\n", sep = "")
print(hh, digits = 2, row.names = FALSE)
write.csv(hh, file.path(out_dir, "12_measured_households_check.csv"), row.names = FALSE)
step_done("09 done: survey checks saved")
