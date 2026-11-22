# =============================================================================
#  14_multiple_testing.R: false discovery rate across the shock estimates
#  Specified before running (plan Doc, "Part 4", 22 November 2026).
#  Benjamini-Hochberg q-values in two families:
#    Primary    the 15 main-model HAZ coefficients (5 shock measures x 3 windows)
#    Secondary  every shock coefficient estimated in steps 07 and 13 (all models,
#               terms and outcomes, placebo and interaction terms included).
#               Household characteristics entered as main effects (livestock,
#               farmland) are not shock coefficients and are left out; a model
#               refit that repeats an estimate already in 07 is counted once.
#  Reads output/08_models.csv and output/13a_food_near_markets.csv.
#  Writes output/14_fdr.csv
# =============================================================================

source("00_settings.R"); source("00_functions.R")
f07 <- file.path(out_dir, "08_models.csv")
f13 <- file.path(out_dir, "13a_food_near_markets.csv")
for (f in c(f07, f13)) if (!file.exists(f)) skip_step(paste0("14: ", basename(f), " not found; run steps 07 and 13 first"))

m07 <- read.csv(f07); m07$step <- "07"
m13 <- read.csv(f13); m13$step <- "13"
cols <- c("step", "model", "term", "estimate", "se", "lower", "upper", "p_value", "n", "clusters")
all <- rbind(m07[, cols], m13[, cols])

shock <- "^(conf|conf10|conf50|food|heat|temp|dry|wet|hotdry)_"
all <- all[grepl(shock, all$term), ]
# The same estimate from a refit in step 13 is the same test: count it once
key <- paste(signif(all$estimate, 10), signif(all$se, 10))
all <- all[!duplicated(key), ]

primary <- all$step == "07" & all$model == "B. All shocks (HAZ)" &
           grepl("^(conf|food|heat|dry|wet)_(preg|y1|y2)$", all$term)
if (sum(primary) != 15) warning("14: expected 15 main-model coefficients, found ", sum(primary))

all$q_primary   <- NA
all$q_primary[primary] <- p.adjust(all$p_value[primary], method = "BH")
all$q_secondary <- p.adjust(all$p_value, method = "BH")
all$in_primary  <- primary
write.csv(all, file.path(out_dir, "14_fdr.csv"), row.names = FALSE)

cat("\nPrimary family (main model, HAZ):", sum(primary), "coefficients;",
    sum(all$p_value[primary] < 0.05), "with p < 0.05;",
    sum(all$q_primary[primary] < 0.05, na.rm = TRUE), "with q < 0.05\n")
print(all[primary, c("term", "estimate", "p_value", "q_primary")], digits = 3, row.names = FALSE)

cat("\nSecondary family (every shock coefficient in steps 07 and 13):", nrow(all), "coefficients;",
    sum(all$p_value < 0.05), "with p < 0.05;", sum(all$q_secondary < 0.05), "with q < 0.05\n")
sig <- all[all$p_value < 0.05, c("model", "term", "estimate", "p_value", "q_secondary")]
print(sig[order(sig$p_value), ], digits = 3, row.names = FALSE)

step_done("14 done: false discovery rate saved")
