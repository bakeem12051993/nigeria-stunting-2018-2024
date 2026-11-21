# =============================================================================
#  08_attribution.R: how much of the change between rounds each shock explains
#  Counterfactual: if 2023-24 children had faced the average 2018 exposures,
#  how different would mean HAZ (and stunting) be? For each domain:
#    contribution = sum over windows of beta x (mean exposure 2023-24 - mean 2018)
#  with delta-method 95% CIs from the main model (07). The exposure means are
#  treated as fixed, so the CIs are somewhat too narrow.
#  Writes output/09_attribution.csv and output/fig_attribution.png
# =============================================================================

source("00_settings.R"); source("00_functions.R")
f_fit <- file.path(derived_dir, "main_fits.rds")
if (!file.exists(f_fit)) skip_step("08: no model results yet; run 07 first")
ch   <- readRDS(file.path(derived_dir, "children_exposures.rds"))
fits <- readRDS(f_fit)
if (length(unique(ch$round)) < 2)
  skip_step("08: attribution compares the two rounds; it runs once the 2018 files are added")

attribute <- function(fit, y, label) {
  if (is.null(fit)) return(NULL)
  x <- names(fit$coef)
  d <- ch[complete.cases(ch[, c(y, x, "wt", "cl_key")]), ]
  r <- d$round == "2023-24"
  dx <- sapply(x, function(v) wmean(d[[v]][r], d$wt[r]) - wmean(d[[v]][!r], d$wt[!r]))
  observed <- wmean(d[[y]][r], d$wt[r]) - wmean(d[[y]][!r], d$wt[!r])
  prefix <- sub("_(preg|y1|y2)$", "", x)
  domain <- names(fits$dom)[sapply(prefix, function(p) which(sapply(fits$dom, function(v) p %in% v)))]
  rows <- lapply(c(unique(domain), "All shocks"), function(g) {
    use <- if (g == "All shocks") x else x[domain == g]
    est <- lincom(fit, setNames(dx[use], use))
    data.frame(outcome = label, domain = g, contribution = est[["estimate"]], se = est[["se"]],
               lower = est[["lower"]], upper = est[["upper"]])
  })
  out <- do.call(rbind, rows)
  out$observed_change <- observed
  out$pct_of_observed <- 100 * out$contribution / observed
  out
}

attr_tab <- rbind(attribute(fits$haz, "haz", "HAZ (SD units)"),
                  attribute(fits$stunted, "stunted", "Stunting (proportion)"))
if (is.null(attr_tab)) skip_step("08: the main models did not run in 07; see the messages there")
write.csv(attr_tab, file.path(out_dir, "09_attribution.csv"), row.names = FALSE)
print(attr_tab, digits = 3)

# Figure: contributions to the change in mean HAZ
a <- attr_tab[attr_tab$outcome == "HAZ (SD units)", ]
if (nrow(a)) {
  png(file.path(out_dir, "fig_attribution.png"), width = 1800, height = 1000, res = 220)
  par(mar = c(4.5, 9, 3, 1))
  lim <- range(pretty(c(a$lower, a$upper, 0, a$observed_change[1])))
  dom_lab <- c(conflict = "Conflict", food = "Food prices", climate = "Climate", "All shocks" = "All shocks")
  bp <- barplot(a$contribution, names.arg = dom_lab[a$domain], horiz = TRUE, las = 1, xlim = lim,
                col = ifelse(a$domain == "All shocks", "#1f3d7a", "#4f81bd"), border = NA,
                xlab = "Contribution to the change in mean HAZ, 2018 to 2023-24 (SD units)",
                main = "How much of the change do the shocks explain?")
  arrows(a$lower, bp, a$upper, bp, angle = 90, code = 3, length = 0.04)
  abline(v = 0, col = "grey40")
  abline(v = a$observed_change[1], lty = 2, col = "#c0504d")
  mtext("Dashed line: observed change", side = 3, line = 0.2, cex = 0.7, col = "#c0504d")
  dev.off()
}
step_done("08 done: attribution saved")
