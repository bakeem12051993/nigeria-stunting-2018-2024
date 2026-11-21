# =============================================================================
#  00_organize_files.R: put every downloaded file where the pipeline reads it
#
#  1. Copy everything you downloaded, unchanged and still zipped, into
#     data/downloads inside the project folder. Any mix works:
#       - DHS downloads (2018 and 2023-24), as one bundle zip or one zip per file
#       - UCDP GED (ged261-csv.zip or the CSV inside it)
#       - World Bank food prices (NGA_RTFP_mkt_*.csv, zipped or not)
#       - ERA5-Land monthly data (.nc files, or the .zip the Climate Data Store gives)
#  2. Run this script. It unzips everything, including zips inside zips, and files:
#       DHS recodes and GPS  -> data/dhs/<dataset>/      e.g. data/dhs/NGKR7BDT/NGKR7BFL.dta
#       UCDP GED csv          -> data/external/ucdp/
#       World Bank prices csv -> data/external/prices/
#       ERA5-Land .nc files   -> data/external/era5/
#  3. It ends with a checklist of what the analysis will find.
#
#  Safe to rerun. Files in data/downloads are never changed or deleted, so you
#  can delete that folder yourself once the checklist looks right.
# =============================================================================

source("00_settings.R"); source("00_functions.R")
downloads_dir <- "data/downloads"

if (!dir.exists(downloads_dir)) {
  dir.create(downloads_dir, recursive = TRUE)
  stop("Created ", normalizePath(downloads_dir), ". Copy your downloads into it and run this script again.",
       call. = FALSE)
}
inbox <- list.files(downloads_dir, full.names = TRUE, recursive = TRUE)
if (!length(inbox)) stop("No files in ", normalizePath(downloads_dir), call. = FALSE)


# ---- 1. Unzip everything into a scratch folder -------------------------------------
stage <- file.path(tempdir(), "organize")
unlink(stage, recursive = TRUE); dir.create(stage)

stem_of <- function(f) tools::file_path_sans_ext(basename(f))
is_zip  <- function(f) grepl("\\.zip$", f, ignore.case = TRUE)

expand <- function(zipfile, into) {              # unzip, then any zips inside, recursively
  dir.create(into, recursive = TRUE, showWarnings = FALSE)
  got <- tryCatch(utils::unzip(zipfile, exdir = into), error = function(e) character(0),
                  warning = function(w) character(0))
  if (!length(got)) { message("  Could not unzip ", basename(zipfile), " (damaged or incomplete download?)"); return(invisible()) }
  for (z in got[is_zip(got)]) {
    expand(z, file.path(dirname(z), stem_of(z)))
    file.remove(z)
  }
}

cat("Unzipping", length(inbox), "downloaded file(s)...\n")
for (f in inbox) {
  if (is_zip(f)) expand(f, file.path(stage, stem_of(f)))
  else { dir.create(file.path(stage, "_loose"), showWarnings = FALSE)
         file.copy(f, file.path(stage, "_loose", basename(f)), overwrite = TRUE) }
}
staged <- list.files(stage, full.names = TRUE, recursive = TRUE)
used   <- rep(FALSE, length(staged))

put <- function(from, to_dir, to_name = basename(from)) {
  dir.create(to_dir, recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(from, file.path(to_dir, to_name), overwrite = TRUE)
  if (!ok) message("  Could not copy ", basename(from), " to ", to_dir)
  ok
}
lower_ext <- function(f) {                       # NGKR7BFL.DTA -> NGKR7BFL.dta
  e <- tools::file_ext(f)
  if (nzchar(e)) paste0(stem_of(f), ".", tolower(e)) else basename(f)
}


# ---- 2. DHS datasets ----------------------------------------------------------------
# A DHS file is named like NGKR7BFL.dta: country, recode type, phase, version, FL.
# Survey recodes live in a ...DT folder (NGKR7BDT); GPS and covariates keep FL
# (NGGE7AFL, NGGC8AFL).
dhs_re <- "^[A-Z]{2}[A-Z]{2}[0-9][A-Z0-9]FL$"
stems  <- toupper(stem_of(staged))
is_dhs <- grepl(dhs_re, stems)
dhs_folder <- function(stem) if (substr(stem, 3, 4) %in% c("GE", "GC")) stem else sub("FL$", "DT", stem)

n_dhs <- 0
for (i in which(is_dhs)) {
  folder <- dhs_folder(stems[i])
  if (put(staged[i], file.path(dhs_dir, folder), lower_ext(staged[i]))) n_dhs <- n_dhs + 1
  used[i] <- TRUE
}
# Documentation that sits beside a single dataset (README, PDF) goes with it
for (d in unique(dirname(staged[is_dhs]))) {
  here <- dirname(staged) == d
  ds <- unique(stems[here & is_dhs])
  if (length(ds) == 1)
    for (i in which(here & !used)) { put(staged[i], file.path(dhs_dir, dhs_folder(ds))); used[i] <- TRUE }
}
dhs_sets <- sort(unique(sapply(stems[is_dhs], dhs_folder)))
cat("DHS:", n_dhs, "files in", length(dhs_sets), "datasets:", paste(dhs_sets, collapse = ", "), "\n")


# ---- 3. Conflict, food prices, climate -------------------------------------------------
b <- basename(staged)

i <- which(!used & grepl("^(GEDEvent|ged).*\\.csv$", b, ignore.case = TRUE))
for (k in i) put(staged[k], file.path(external_dir, "ucdp"))
used[i] <- TRUE
cat("UCDP:", if (length(i)) paste(b[i], collapse = ", ") else "none", "\n")

i <- which(!used & grepl("RTFP.*\\.csv$", b, ignore.case = TRUE))
for (k in i) put(staged[k], file.path(external_dir, "prices"))
used[i] <- TRUE
cat("Food prices:", if (length(i)) paste(b[i], collapse = ", ") else "none", "\n")

# Climate Data Store zips often hold files with the same names (e.g.
# data_stream-moda_stepType-avgua.nc), so files from a zip get the zip's name
# in front to keep them apart.
i <- which(!used & grepl("\\.nc$", b, ignore.case = TRUE))
for (k in i) {
  top <- strsplit(substring(staged[k], nchar(stage) + 2), "[/\\\\]")[[1]][1]
  nm  <- if (top == "_loose") b[k] else paste0(top, "_", b[k])
  put(staged[k], era5_dir, nm)
}
used[i] <- TRUE
cat("ERA5-Land:", if (length(i)) paste(b[i], collapse = ", ") else "none", "\n")

i <- which(!used & grepl("floods?.*\\.csv$", b, ignore.case = TRUE))
for (k in i) put(staged[k], dirname(flood_file), basename(flood_file))
used[i] <- TRUE

if (any(!used)) cat("Not used by the analysis (left in data/downloads):",
                    paste(unique(b[!used]), collapse = ", "), "\n")
unlink(stage, recursive = TRUE)


# ---- 4. Checklist: what the analysis will find -----------------------------------------
source("00_settings.R")                          # picks up the new file names
cat("\n==================== Checklist ====================\n")
rounds <- tryCatch(available_rounds(report = TRUE), error = function(e) { message(conditionMessage(e)); character(0) })
cat("\nDHS rounds ready:", if (length(rounds)) paste(round_label[rounds], collapse = " and ") else "none", "\n")

chk <- function(label, ok, detail) cat(sprintf("%-22s %-8s %s\n", label, if (ok) "ready" else "MISSING", detail))
chk("Conflict (UCDP)", file.exists(ucdp_file), ucdp_file)
chk("Food prices (RTFP)", file.exists(prices_file), prices_file)
e5 <- find_era5(era5_dir)
for (v in c("t2m", "d2m", "tp"))
  chk(paste("ERA5-Land", v), !is.null(e5[[v]]), if (is.null(e5[[v]])) era5_dir else e5[[v]]$file)
chk("Floods (optional)", file.exists(flood_file), flood_file)
cat("\nNext: run run_all.R. Steps whose data are MISSING are skipped with a message.\n")
