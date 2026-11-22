# =============================================================================
#  12_submission_figures.R: figure files for journal submission
#  Copies the figures from steps 10 and 11 to output/submission_figures under
#  their manuscript numbers, checks that each is at least 600 dpi at 174 mm
#  wide, and zips them into output/submission_figures.zip.
#  Also writes output/session_info.txt (R and package versions, for Section 2.8).
#
#  Manuscript order:
#    Figure 1  fig1_timeline.png        shock timeline
#    Figure 2  fig2_zones.png           stunting by zone
#    Figure 3  fig_map_stunting.png     stunting by state (step 11)
#    Figure 4  fig3_age.png             stunting by age
#    Figure 5  fig4_survey_design.png   survey design
#    Figure 6  fig5_oaxaca.png          Oaxaca-Blinder
#    Figure 7  fig6_shocks.png          shocks by window
# =============================================================================

source("00_settings.R"); source("00_functions.R")
fig_dir <- file.path(out_dir, "figures")
sub_dir <- file.path(out_dir, "submission_figures")

fig_map <- c(fig1_timeline.png      = "Figure1.png",
             fig2_zones.png         = "Figure2.png",
             fig_map_stunting.png   = "Figure3.png",
             fig3_age.png           = "Figure4.png",
             fig4_survey_design.png = "Figure5.png",
             fig5_oaxaca.png        = "Figure6.png",
             fig6_shocks.png        = "Figure7.png")

missing <- names(fig_map)[!file.exists(file.path(fig_dir, names(fig_map)))]
if (length(missing))
  skip_step(paste("12: run steps 10 and 11 first; not found:", paste(missing, collapse = ", ")))

# ---- Copy under manuscript numbers -------------------------------------------------
dir.create(sub_dir, showWarnings = FALSE, recursive = TRUE)
ok <- file.copy(file.path(fig_dir, names(fig_map)), file.path(sub_dir, fig_map), overwrite = TRUE)
if (!all(ok)) stop("12: could not copy ", paste(names(fig_map)[!ok], collapse = ", "))

# ---- Check resolution: pixel width and height from the PNG header -------------------
png_size <- function(f) {
  con <- file(f, "rb"); on.exit(close(con))
  b <- readBin(con, "raw", 24)
  if (!identical(b[2:4], charToRaw("PNG"))) return(c(NA, NA))
  to_int <- function(x) sum(as.integer(x) * 256^(3:0))
  c(width = to_int(b[17:20]), height = to_int(b[21:24]))
}
width_in <- 174 / 25.4
chk <- do.call(rbind, lapply(names(fig_map), function(src) {
  wh <- png_size(file.path(sub_dir, fig_map[[src]]))
  data.frame(figure = fig_map[[src]], source = src,
             width_px = wh[1], height_px = wh[2],
             dpi_at_174mm = round(wh[1] / width_in),
             height_mm = round(wh[2] / (wh[1] / width_in) * 25.4))
}))
rownames(chk) <- NULL
print(chk, row.names = FALSE)
low <- chk$figure[is.na(chk$dpi_at_174mm) | chk$dpi_at_174mm < 590]
if (length(low))
  warning("12: below 600 dpi: ", paste(low, collapse = ", "),
          ". Rerun steps 10 and 11 (both now save at 600 dpi), then this step.")

# ---- Zip -----------------------------------------------------------------------------
zip_file <- file.path(normalizePath(out_dir), "submission_figures.zip")
if (file.exists(zip_file)) file.remove(zip_file)
files_abs <- normalizePath(file.path(sub_dir, fig_map))
zipped <- tryCatch({
  if (requireNamespace("zip", quietly = TRUE)) {
    zip::zip(zip_file, files = files_abs, mode = "cherry-pick")
  } else {
    if (utils::zip(zip_file, files = files_abs, flags = "-j9Xq") != 0) stop("zip program failed")
  }
  file.exists(zip_file)
}, error = function(e) { message("12: could not zip (", conditionMessage(e), ")"); FALSE })
if (!zipped)
  message("12: install the zip package (install.packages(\"zip\")) and rerun, or zip the ",
          sub_dir, " folder by hand.")

# ---- R and package versions ------------------------------------------------------------
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
cat("\n", R.version.string, "\n", sep = "")

step_done(paste0("12 done: figures in ", sub_dir,
                 if (zipped) paste0(" and ", zip_file) else "",
                 "; versions in ", file.path(out_dir, "session_info.txt")))
