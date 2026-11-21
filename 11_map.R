# =============================================================================
#  11_map.R: map of stunting by state, 2018 and 2023-24, and the change
#  (A) 2018, (B) 2023-24, (C) change in percentage points. Zone boundaries are
#  drawn over the states. Estimates use the analysis weights (measured children
#  matched to all children by state and urban/rural), as in step 02.
#  State boundaries: GADM version 4.1, downloaded once with the geodata package
#  into data/external/boundaries (needs internet the first time).
#  Writes output/figures/fig_map_stunting.png (300 dpi, 174 mm wide) and
#  output/11_state_stunting.csv
# =============================================================================

source("00_settings.R"); source("00_functions.R")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
for (p in c("terra", "survey")) if (!requireNamespace(p, quietly = TRUE))
  skip_step(paste0("11: package ", p, " is not installed"))


# ---- State boundaries ------------------------------------------------------------------
bnd_dir  <- file.path(external_dir, "boundaries")
bnd_file <- file.path(bnd_dir, "nga_states_gadm41.gpkg")
if (!file.exists(bnd_file)) {
  if (!requireNamespace("geodata", quietly = TRUE))
    skip_step("11: install the geodata package once (install.packages(\"geodata\")) to download state boundaries")
  dir.create(bnd_dir, showWarnings = FALSE, recursive = TRUE)
  g <- geodata::gadm(country = "NGA", level = 1, path = bnd_dir, version = "4.1")
  if (is.null(g)) skip_step("11: the GADM download failed; check the internet connection and rerun")
  terra::writeVector(g, bnd_file, overwrite = TRUE)
}
states <- terra::vect(bnd_file)
states$state <- state_key(states$NAME_1)
states$state[states$state %in% c("federalcapitalterritory", "abujafederalcapitalterritory")] <- "fctabuja"
states$state[states$state %in% "nassarawa"] <- "nasarawa"
states <- states[states$state %in% names(zone_of_state), ]
miss <- setdiff(names(zone_of_state), states$state)
if (length(miss)) stop("CHECK boundary names; no polygon for: ", paste(miss, collapse = ", "),
                       ". Names in the file: ", paste(sort(unique(terra::vect(bnd_file)$NAME_1)), collapse = ", "))
states$zone <- zone_of_state[states$state]
# For display: simplify the coastline, then build zone outlines with a small
# buffer (closes the slivers between neighboring states in the boundary file),
# filled holes, and islands under 100 km2 dropped (so the delta coast does not
# draw as dark blobs). The buffer is 500 m, invisible at this scale.
states <- terra::simplifyGeom(states, tolerance = 0.005)
zones  <- terra::fillHoles(terra::aggregate(terra::buffer(states, 500), by = "zone"))
zparts <- terra::disagg(zones)
zones  <- terra::aggregate(zparts[terra::expanse(zparts, unit = "km") > 100, ], by = "zone")


# ---- Stunting by state -----------------------------------------------------------------
ana <- readRDS(file.path(derived_dir, "children.rds"))
ps_totals <- readRDS(file.path(derived_dir, "ps_totals.rds"))
des <- svydesign(ids = ~psu, strata = ~strata, weights = ~wt_dhs, data = ana, nest = TRUE)
if (!is.null(ps_totals)) des <- postStratify(des, ~ps_cell, ps_totals, partial = TRUE)
est <- svyby(~stunted, ~state + round, des, svymean, na.rm = TRUE)
st  <- data.frame(state = est$state, round = est$round, stunted = est$stunted, se = SE(est), row.names = NULL)
n   <- aggregate(stunted ~ state + round, ana, length); names(n)[3] <- "n_children"
st  <- merge(st, n)
wide <- reshape(st, idvar = "state", timevar = "round", direction = "wide")
names(wide) <- sub("\\.2023-24$", "_2024", sub("\\.2018$", "_2018", names(wide)))
wide$change_pp    <- 100 * (wide$stunted_2024 - wide$stunted_2018)
wide$change_lower <- wide$change_pp - 1.96 * 100 * sqrt(wide$se_2018^2 + wide$se_2024^2)
wide$change_upper <- wide$change_pp + 1.96 * 100 * sqrt(wide$se_2018^2 + wide$se_2024^2)
wide$zone <- zone_display[zone_of_state[wide$state]]
wide <- wide[order(wide$zone, wide$state), ]
write.csv(wide, file.path(out_dir, "11_state_stunting.csv"), row.names = FALSE)
cat("\nStunting by state (%), analysis weights:\n")
print(data.frame(state = wide$state, `2018` = round(100 * wide$stunted_2018, 1),
                 `2023-24` = round(100 * wide$stunted_2024, 1), change = round(wide$change_pp, 1),
                 n_2018 = wide$n_children_2018, n_2024 = wide$n_children_2024, check.names = FALSE),
      row.names = FALSE)


# ---- Map ---------------------------------------------------------------------------------
ink <- "#0b0b0b"; ink2 <- "#52514e"
seq_breaks <- c(0, 0.2, 0.3, 0.4, 0.5, 0.6, 1.0001)
seq_labels <- c("<20", "20-30", "30-40", "40-50", "50-60", "60+")
seq_cols   <- c("#eaf2fc", "#c3daf6", "#8fbcef", "#4f94e2", "#1f66c2", "#0d3570")
div_breaks <- c(-Inf, -10, -5, -2, 2, 5, 10, Inf)
div_labels <- c("\u2264 -10", "-10 to -5", "-5 to -2", "-2 to 2", "2 to 5", "5 to 10", "\u2265 10")
div_cols   <- c("#2a78d6", "#86b6ef", "#cfe1f8", "#e8e7e1", "#f8cdb9", "#f1946a", "#c9501f")
fill_of <- function(x, breaks, cols) cols[findInterval(x, breaks, rightmost.closed = TRUE, left.open = FALSE)]

m <- match(states$state, wide$state)
fills <- list(fill_of(wide$stunted_2018[m], seq_breaks, seq_cols),
              fill_of(wide$stunted_2024[m], seq_breaks, seq_cols),
              fill_of(wide$change_pp[m],    div_breaks, div_cols))
titles <- c("A. Stunting, 2018 (%)", "B. Stunting, 2023-24 (%)", "C. Change (percentage points)")

ext <- terra::ext(states)
png(file.path(fig_dir, "fig_map_stunting.png"), width = 174 / 25.4, height = 3.15, units = "in",
    res = 300, bg = "white", type = "cairo")
layout(matrix(c(1, 2, 3, 4, 4, 5), nrow = 2, byrow = TRUE), heights = c(1, 0.36))
par(family = "sans", cex = 0.8)
for (k in 1:3) {
  par(mar = c(0.2, 0.2, 1.6, 0.2))
  terra::plot(states, col = fills[[k]], border = "white", lwd = 0.4, axes = FALSE,
       legend = FALSE, mar = NA, xlim = ext[1:2], ylim = ext[3:4])
  terra::lines(zones, col = ink2, lwd = 1.1)
  title(main = titles[k], adj = 0, cex.main = 0.95, font.main = 2, col.main = ink, line = 0.4)
}
legend_block <- function(labels, cols, title, ncol, cex = 0.82) {
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("top", legend = labels, fill = cols, border = ink2, ncol = ncol, bty = "n",
         title = title, title.adj = 0.5, cex = cex, x.intersp = 0.3, y.intersp = 0.9,
         text.col = ink2, title.col = ink, xpd = NA)
}
legend_block(seq_labels, seq_cols, "Stunting (%)", ncol = 6)
legend_block(div_labels, div_cols, "Change (points)", ncol = 4, cex = 0.74)
dev.off()
cat("\nMap written to", file.path(fig_dir, "fig_map_stunting.png"), "\n")
step_done("11 done: state map")
