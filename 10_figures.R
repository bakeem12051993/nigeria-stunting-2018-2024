# =============================================================================
#  10_figures.R: the manuscript figures
#  Figure 1  Shock timeline: conflict, food inflation and humid heat by month,
#            with survey fieldwork and the birth cohort behind the rise marked
#  Figure 2  Stunting by zone, 2018 and 2023-24
#  Figure 3  Stunting by age, 2018 and 2023-24
#  Figure 4  Survey design: (A) children under 5 per household, households picked
#            for measuring versus others; (B) the North's share of households and
#            of children under 5 in three surveys
#  Figure 5  Oaxaca-Blinder: (A) contributions by domain; (B) by subgroup
#  Figure 6  Shock associations with HAZ by window, with the placebo window
#  Reads the output tables from steps 01 to 09, the exposure panels in derived/
#  and the household recodes (for Figure 4B). Base R graphics only.
#  Writes output/figures/fig1_timeline.png ... fig6_shocks.png (300 dpi, 174 mm wide)
# =============================================================================

source("00_settings.R"); source("00_functions.R")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
need <- function(f) {
  p <- file.path(out_dir, f)
  if (!file.exists(p)) skip_step(paste0("10: ", f, " not found; run steps 01 to 09 first"))
  read.csv(p)
}

# ---- Style -------------------------------------------------------------------------
col_a    <- "#2a78d6"   # first series (2018; picked households; households)
col_b    <- "#eb6834"   # second series (2023-24; other households; children)
col_c    <- "#1baf7a"   # third series
win_cols <- c(pre = "#8c8c8c", preg = "#86b6ef", y1 = "#2a78d6", y2 = "#104281")
ink      <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e1e0d9"; axis_col <- "#c3c2b7"
band_fw  <- "#e8e7e1"; band_coh <- "#dbe9fa"
width_in <- 174 / 25.4

open_fig <- function(name, height) {
  png(file.path(fig_dir, name), width = width_in, height = height, units = "in",
      res = 300, bg = "white", type = "cairo")
  par(family = "sans", col.axis = ink2, col.lab = ink, fg = axis_col, cex = 0.8,
      mgp = c(2, 0.5, 0), tcl = -0.25, las = 1)
}
panel_title <- function(txt, adj = 0) title(main = txt, adj = adj, cex.main = 0.95,
                                              font.main = 2, col.main = ink, line = 0.6)
hgrid <- function(at) abline(h = at, col = grid_col, lwd = 0.6)
vgrid <- function(at) abline(v = at, col = grid_col, lwd = 0.6)
pct   <- function(x, d = 1) formatC(100 * x, format = "f", digits = d)
num1   <- function(x) ifelse(abs(x) < 0.05 & x != 0, sprintf("%.2f", x), sprintf("%.1f", x))
fmt_ci <- function(e, l, u) paste0(ifelse(e >= 0, "+", ""), num1(e), " (", num1(l), " to ", num1(u), ")")
set_cex <- function() par(cex = 0.8)   # layout() and mfrow reset cex


# ---- Figure 1: shock timeline ---------------------------------------------------------
conflict <- readRDS(file.path(derived_dir, "conflict_panels.rds"))
prices   <- readRDS(file.path(derived_dir, "price_panels.rds"))
climate  <- readRDS(file.path(derived_dir, "climate_panels.rds"))
if (!is.null(conflict) && !is.null(prices) && !is.null(climate)) {
  roll <- function(x, k) as.numeric(stats::filter(x, rep(1 / k, k), sides = 1))
  cm  <- as.integer(colnames(conflict$panels[[paste0("deaths_r", conflict_radius_km)]]))
  yr  <- cmc_year(cm) + (cmc_month(cm) - 0.5) / 12
  shown <- cm >= cmc_from_ym(2013, 1) & cm <= cmc_from_ym(2024, 12)
  d25 <- conflict$panels[[paste0("deaths_r", conflict_radius_km)]]
  conf_share <- 100 * colMeans(d25 > 0, na.rm = TRUE)
  food <- colMeans(prices$panels$food_infl[, as.character(cm)], na.rm = TRUE)
  heat <- colMeans(climate$panels$heat[, as.character(cm)], na.rm = TRUE)
  series <- list(
    list(y = roll(conf_share, 3), raw = NULL, ylab = "% of clusters",
         title = sprintf("A. Clusters with a conflict death within %d km (3-month average)", conflict_radius_km)),
    list(y = food, raw = NULL, ylab = "% (log change)",
         title = "B. Food price inflation at the nearest market (year on year)"),
    list(y = roll(heat, 12), raw = heat, ylab = "SD",
         title = "C. Humid-heat anomaly (bars: monthly; line: 12-month average)"))
  fw   <- list(c(cmc_from_ym(2018, 8), cmc_from_ym(2018, 12)), c(cmc_from_ym(2023, 12), cmc_from_ym(2024, 5)))
  coh  <- c(cmc_from_ym(2019, 1), cmc_from_ym(2021, 5))
  dec  <- function(c) cmc_year(c) + (cmc_month(c) - 1) / 12
  open_fig("fig1_timeline.png", 6.2)
  par(mfrow = c(3, 1), mar = c(2.2, 4, 2.2, 1), oma = c(1.2, 0, 2.8, 0)); set_cex()
  for (k in seq_along(series)) {
    s <- series[[k]]
    yv <- c(s$y[shown], if (!is.null(s$raw)) s$raw[shown])
    yl <- range(c(0, yv), na.rm = TRUE); yl <- yl + c(-0.04, 0.08) * diff(yl)
    plot(NA, xlim = c(2013, 2025), ylim = yl, xlab = "", ylab = s$ylab, axes = FALSE, xaxs = "i")
    for (b in fw) rect(dec(b[1]), yl[1], dec(b[2]) + 1 / 12, yl[2], col = band_fw, border = NA)
    rect(dec(coh[1]), yl[1], dec(coh[2]) + 1 / 12, yl[2], col = band_coh, border = NA)
    hgrid(pretty(yl)); abline(h = 0, col = axis_col, lwd = 0.8)
    if (!is.null(s$raw)) segments(yr[shown], 0, yr[shown], s$raw[shown], col = "#b7d3f6", lwd = 1.2, lend = 1)
    lines(yr[shown], s$y[shown], col = col_a, lwd = 1.8)
    axis(1, at = 2013:2025, labels = if (k == 3) 2013:2025 else FALSE, col = axis_col)
    axis(2, at = pretty(yl), col = axis_col)
    panel_title(s$title)
    if (k == 1)                    # band labels go in the outer margin above panel A
      lab_x <- grconvertX(c(dec(fw[[1]][2]) + 1 / 12, mean(dec(coh)) + 1 / 24,
                            mean(dec(fw[[2]])) + 1 / 24), "user", "ndc")
  }
  mtext(c("2018 survey", "Births of children aged\n36-59 months in 2023-24", "2023-24 survey"),
        side = 3, outer = TRUE, at = lab_x, adj = c(1, 0.5, 0.5), line = 0.3, cex = 0.62,
        col = ink2)
  dev.off()
} else message("Figure 1 skipped: exposure panels not found (run steps 03 to 05)")


# ---- Figure 2: stunting by zone ------------------------------------------------------
bz <- need("03b_stunting_by_zone.csv"); ch <- need("03e_change_by_zone.csv"); nat <- need("03a_change_national.csv")
nat <- nat[nat$outcome == "stunted", ]
zones <- unique(bz$zone)
z24 <- bz[bz$round == "2023-24", ]; zones <- z24$zone[order(z24$estimate)]
rows <- c("Nigeria", zones)                     # bottom to top
get_row <- function(r, rnd) {
  if (r == "Nigeria") nat[nat$round == rnd, c("estimate", "lower", "upper")]
  else bz[bz$zone == r & bz$round == rnd, c("estimate", "lower", "upper")]
}
ypos <- c(1, seq_along(zones) + 1.6)
open_fig("fig2_zones.png", 3.6)
par(mar = c(3.2, 8, 2.6, 11))
plot(NA, xlim = c(0.10, 0.62), ylim = c(0.4, max(ypos) + 0.6), axes = FALSE,
     xlab = "Children under 5 stunted (%)", ylab = "")
vgrid(seq(0.1, 0.6, 0.1)); abline(h = 1.8, col = axis_col, lwd = 0.6, lty = 3)
for (i in seq_along(rows)) {
  a <- get_row(rows[i], "2018"); b <- get_row(rows[i], "2023-24"); y <- ypos[i]
  segments(a$estimate, y, b$estimate, y, col = axis_col, lwd = 2.5)
  segments(a$lower, y + 0.16, a$upper, y + 0.16, col = col_a, lwd = 1)
  segments(b$lower, y - 0.16, b$upper, y - 0.16, col = col_b, lwd = 1)
  points(a$estimate, y, pch = 21, bg = col_a, col = "white", cex = 1.5, lwd = 1)
  points(b$estimate, y, pch = 21, bg = col_b, col = "white", cex = 1.5, lwd = 1)
  cg <- if (rows[i] == "Nigeria") ch[ch$group == "Nigeria", ] else ch[ch$group == rows[i], ]
  sig <- cg$lower > 0 | cg$upper < 0
  mtext(fmt_ci(cg$change_pp, cg$lower, cg$upper), side = 4, at = y, las = 1, line = 0.6,
        cex = 0.75, col = ink, font = if (sig) 2 else 1)
}
axis(1, at = seq(0.1, 0.6, 0.1), labels = paste0(seq(10, 60, 10)), col = axis_col)
axis(2, at = ypos, labels = rows, tick = FALSE, col.axis = ink, font.axis = 1)
mtext("Change, points (95% CI)", side = 4, at = max(ypos) + 0.75, las = 1, line = 0.6, cex = 0.75, col = ink2)
legend("top", inset = c(0, -0.13), xpd = NA, horiz = TRUE, bty = "n", pch = 21,
       pt.bg = c(col_a, col_b), col = "white", pt.cex = 1.4, legend = c("2018", "2023-24"),
       text.col = ink, cex = 0.85)
dev.off()


# ---- Figure 3: stunting by age ----------------------------------------------------------
ba <- need("03d_stunting_by_age.csv")
bands <- c("0-5", "6-11", "12-23", "24-35", "36-47", "48-59")
open_fig("fig3_age.png", 3.4)
par(mar = c(3.4, 4, 2.4, 1))
plot(NA, xlim = c(0.6, 6.4), ylim = c(0.1, 0.58), axes = FALSE,
     xlab = "Age (months)", ylab = "Children stunted (%)")
rect(4.5, 0.1, 6.5, 0.58, col = band_coh, border = NA)
hgrid(seq(0.1, 0.5, 0.1))
for (k in 1:2) {
  rnd <- c("2018", "2023-24")[k]; cl <- c(col_a, col_b)[k]; off <- c(-0.07, 0.07)[k]
  s <- ba[ba$round == rnd, ]; s <- s[match(bands, s$age_band), ]
  segments(1:6 + off, s$lower, 1:6 + off, s$upper, col = cl, lwd = 1)
  lines(1:6 + off, s$estimate, col = cl, lwd = 1.8)
  points(1:6 + off, s$estimate, pch = 21, bg = cl, col = "white", cex = 1.3)
}
axis(1, at = 1:6, labels = bands, col = axis_col)
axis(2, at = seq(0.1, 0.5, 0.1), labels = seq(10, 50, 10), col = axis_col)
text(5.5, 0.115, "In 2023-24, born from\nlate 2018 to mid-2021", cex = 0.72, col = ink2, adj = c(0.5, 0))
legend("topleft", bty = "n", pch = 21, pt.bg = c(col_a, col_b), col = c(col_a, col_b), pt.cex = 1.3,
       lwd = 1.8, lty = 1, legend = c("2018", "2023-24"), text.col = ink, cex = 0.85)
dev.off()


# ---- Figure 4: survey design ---------------------------------------------------------------
mh  <- need("12_measured_households_check.csv")
zs3 <- need("10_zone_shares_three_surveys.csv")
north <- c("North Central", "North East", "North West")
# Households: weighted share in the North, from each survey's household recode
hh_files <- list("2018 DHS" = list(f = dhs_path("2018", "hr"), z = pr_zone_var[["2018"]]),
                 "2023-24 DHS" = list(f = dhs_path("2024", "hr"), z = pr_zone_var[["2024"]]))
mis_hr <- list.files(dhs_dir, pattern = "^NGHR[0-9][0-9]FL\\.dta$", recursive = TRUE, ignore.case = TRUE)
if (length(mis_hr)) hh_files[["2021 MIS"]] <- list(f = file.path(dhs_dir, mis_hr[1]), z = "hv024")
hh_north <- sapply(hh_files, function(s) {
  h <- read_dhs(s$f, c("hv005", s$z))
  z <- zone_from(lab_of(h, s$z)); w <- getv(h, "hv005")
  100 * sum(w[z %in% north]) / sum(w[!is.na(z)])
})
u5_north <- tapply(zs3$share_pct * (zs3$zone %in% north), zs3$survey, sum)
surveys <- intersect(c("2018 DHS", "2021 MIS", "2023-24 DHS"), names(u5_north))

open_fig("fig4_survey_design.png", 3.5)
layout(matrix(1:3, 1), widths = c(2.44, 1.51, 2.9)); set_cex(); par(oma = c(1.4, 0, 0, 0))
zone_order <- rev(c("North West", "North East", "North Central", "South East", "South South", "South West"))
for (k in 1:2) {
  rnd <- c("2018", "2023-24")[k]
  s <- mh[mh$round == rnd, ]; s <- s[match(zone_order, s$zone), ]
  par(mar = c(2.2, if (k == 1) 6.6 else 0.8, 2.4, 0.8))
  plot(NA, xlim = c(0.2, 1.5), ylim = c(0.5, 6.5), axes = FALSE, xlab = "", ylab = "")
  vgrid(seq(0.4, 1.4, 0.2))
  segments(s$under5_listed_other, 1:6, s$under5_listed_picked, 1:6, col = axis_col, lwd = 2.5)
  points(s$under5_listed_picked, 1:6, pch = 21, bg = col_a, col = "white", cex = 1.4)
  points(s$under5_listed_other, 1:6, pch = 21, bg = col_b, col = "white", cex = 1.4)
  axis(1, at = seq(0.4, 1.4, 0.2), labels = c("0.4", "", "0.8", "", "1.2", ""), col = axis_col)
  if (k == 1) axis(2, at = 1:6, labels = zone_order, tick = FALSE, col.axis = ink)
  panel_title(paste0(c("A. ", "")[k], rnd))
  if (k == 2) legend("topleft", bty = "n", pch = 21, pt.bg = c(col_a, col_b), col = "white",
                     pt.cex = 1.2, legend = c("Picked for\nmeasuring", "Other"),
                     text.col = ink, cex = 0.72, y.intersp = 1.3)
}
mtext("Children under 5 per household", side = 1, outer = TRUE, at = 0.36, line = 0.3, cex = 0.7, col = ink)
par(mar = c(2.2, 4, 2.4, 5.2))
xs <- seq_along(surveys)
yl <- c(45, 80)
plot(NA, xlim = c(0.6, length(xs) + 0.4), ylim = yl, axes = FALSE, xlab = "", ylab = "Share in the North (%)")
hgrid(seq(45, 80, 5))
if ("2021 MIS" %in% surveys) {
  abline(v = 1.5, col = ink2, lty = 2, lwd = 0.8)
  text(1.5, yl[1] + 1, "New sampling\nframe", pos = 4, cex = 0.68, col = ink2, offset = 0.2)
}
hn <- hh_north[surveys]; un <- u5_north[surveys]
lines(xs, hn, col = col_a, lwd = 1.8); points(xs, hn, pch = 21, bg = col_a, col = "white", cex = 1.3)
lines(xs, un, col = col_b, lwd = 1.8); points(xs, un, pch = 21, bg = col_b, col = "white", cex = 1.3)
text(xs, hn, sprintf("%.1f", hn), pos = 1, cex = 0.7, col = ink)
text(xs, un, sprintf("%.1f", un), pos = 3, cex = 0.7, col = ink)
mtext(c("Households", "Children\nunder 5"), side = 4, at = c(tail(hn, 1), tail(un, 1)),
      las = 1, line = 0.4, cex = 0.62, col = ink)
axis(1, at = xs, labels = sub(" ", "\n", surveys), col = axis_col, padj = 0.4, cex.axis = 0.8, mgp = c(2, 0.3, 0))
axis(2, at = seq(45, 80, 5), col = axis_col)
panel_title("B. The North's share")
dev.off()


# ---- Figure 5: Oaxaca-Blinder ---------------------------------------------------------------
ob  <- need("05_oaxaca_main.csv"); obs <- need("06_oaxaca_subgroups.csv")
m <- ob[ob$analysis == "HAZ, all children", ]
dom <- m[grepl("^explained: ", m$term), ]; dom$label <- sub("^explained: ", "", dom$term)
summ <- m[match(c("explained", "unexplained", "total"), m$term), ]
summ$label <- c("All characteristics", "Unexplained", "Total change")
A <- rbind(dom[, c("label", "estimate", "lower", "upper")], summ[, c("label", "estimate", "lower", "upper")])
yA <- c(rev(seq_len(nrow(dom))) + 3.6, 3, 2, 1)
groups <- c("HAZ, all children", "HAZ, North", "HAZ, South", "HAZ, urban", "HAZ, rural")
glab   <- c("All children", "North", "South", "Urban", "Rural")
B <- rbind(m, obs)
open_fig("fig5_oaxaca.png", 3.8)
layout(matrix(1:2, 1), widths = c(1.35, 1)); set_cex()
par(mar = c(3.4, 9.5, 2.4, 1))
xl <- c(-0.16, 0.05)
plot(NA, xlim = xl, ylim = c(0.5, max(yA) + 0.5), axes = FALSE,
     xlab = "Contribution to change in mean HAZ (SD)", ylab = "")
vgrid(seq(-0.15, 0.05, 0.05)); abline(v = 0, col = ink2, lwd = 0.8)
abline(h = 3.8, col = axis_col, lty = 3, lwd = 0.6)
fillc <- ifelse(seq_len(nrow(A)) > nrow(dom), "#104281", col_a)
rect(pmin(0, A$estimate), yA - 0.3, pmax(0, A$estimate), yA + 0.3, col = fillc, border = NA)
segments(A$lower, yA, A$upper, yA, col = ink, lwd = 0.9)
axis(1, at = seq(-0.15, 0.05, 0.05), labels = sprintf("%.2f", seq(-0.15, 0.05, 0.05)), col = axis_col)
axis(2, at = yA, labels = A$label, tick = FALSE, col.axis = ink)
panel_title("A. All children, by domain")
par(mar = c(3.4, 5.2, 2.4, 1))
plot(NA, xlim = c(-0.25, 0.17), ylim = c(0.5, length(groups) + 0.5), axes = FALSE,
     xlab = "Change in mean HAZ (SD)", ylab = "")
vgrid(seq(-0.2, 0.15, 0.1)); abline(v = 0, col = ink2, lwd = 0.8)
terms3 <- c("total", "explained", "unexplained"); cols3 <- c(ink, col_a, col_b); off3 <- c(0.22, 0, -0.22)
pch3 <- c(23, 21, 21)
for (g in seq_along(groups)) {
  y <- length(groups) + 1 - g
  for (t in 1:3) {
    r <- B[B$analysis == groups[g] & B$term == terms3[t], ]
    segments(r$lower, y + off3[t], r$upper, y + off3[t], col = cols3[t], lwd = 1)
    points(r$estimate, y + off3[t], pch = pch3[t], bg = cols3[t], col = "white", cex = 1.2)
  }
}
axis(1, at = seq(-0.2, 0.1, 0.1), labels = sprintf("%.1f", seq(-0.2, 0.1, 0.1)), col = axis_col)
axis(2, at = length(groups):1, labels = glab, tick = FALSE, col.axis = ink)
legend("topright", bty = "n", pch = pch3, pt.bg = cols3, col = "white", pt.cex = 1.1,
       legend = c("Total", "Explained", "Unexplained"), text.col = ink, cex = 0.75)
panel_title("B. By subgroup")
dev.off()


# ---- Figure 6: shock associations by window ----------------------------------------------
mods <- need("08_models.csv")
main <- mods[mods$model == "B. All shocks (HAZ)", ]
plac <- mods[mods$model == "D. Placebo: before conception (HAZ)" & grepl("_pre$", mods$term), ]
doms <- list(conf = c("A. Conflict", "per log(1 + deaths within 25 km)"),
             food = c("B. Food inflation", "per 10 points"),
             heat = c("C. Humid heat", "per SD of the anomaly"),
             dry  = c("D. Drought (SPEI-3)", "per SD"),
             wet  = c("E. Very wet months", "per month"))
wins <- c(pre = "Before conception\n(placebo)", preg = "Pregnancy", y1 = "Year 1", y2 = "Year 2")
open_fig("fig6_shocks.png", 4.6)
layout(matrix(1:6, 2, byrow = TRUE), widths = c(2.92, 1.96, 1.96)); set_cex()
par(oma = c(1.6, 0, 0, 0))
for (k in seq_along(doms)) {
  p <- names(doms)[k]
  r <- rbind(plac[plac$term == paste0(p, "_pre"), ],
             main[main$term %in% paste0(p, c("_preg", "_y1", "_y2")), ])
  if (!nrow(r)) { plot.new(); next }
  r$win <- sub(".*_", "", r$term); r <- r[match(names(wins), r$win), ]
  y <- 4:1
  par(mar = c(3.4, if (k %in% c(1, 4)) 7.2 else 1.2, 2.2, 0.8))
  xl <- range(c(0, r$lower, r$upper), na.rm = TRUE); xl <- xl + c(-0.08, 0.08) * diff(xl)
  plot(NA, xlim = xl, ylim = c(0.5, 4.5), axes = FALSE, xlab = "", ylab = "")
  vgrid(pretty(xl, 4)); abline(v = 0, col = ink2, lwd = 0.8)
  abline(h = 3.5, col = axis_col, lty = 3, lwd = 0.6)
  segments(r$lower, y, r$upper, y, col = win_cols[r$win], lwd = 1.2)
  points(r$estimate, y, pch = ifelse(r$win == "pre", 21, 19), cex = 1.25,
         col = win_cols[r$win], bg = "white", lwd = 1.2)
  axis(1, at = pretty(xl, 4), col = axis_col)
  if (k %in% c(1, 4)) axis(2, at = y, labels = wins, tick = FALSE, col.axis = ink, cex.axis = 0.95)
  panel_title(doms[[p]][1]); mtext(doms[[p]][2], side = 1, line = 2, cex = 0.62, col = ink2)
}
par(mar = c(0, 0, 0, 0)); plot.new()
legend("center", bty = "n", pch = c(21, 19, 19, 19), col = win_cols, pt.bg = "white", pt.cex = 1.2,
       legend = c("Before conception (placebo)", "Pregnancy", "Year 1", "Year 2"),
       text.col = ink, cex = 0.8, lwd = 1.2, lty = 1, seg.len = 1.6,
       title = "Window", title.col = ink2)
mtext("Change in HAZ (SD) per unit of exposure, with 95% CI; main model with all shocks together",
      side = 1, outer = TRUE, line = 0.4, cex = 0.65, col = ink2)
dev.off()

step_done(paste("10 done: figures saved in", fig_dir))
