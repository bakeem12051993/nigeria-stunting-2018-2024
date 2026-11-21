# =============================================================================
#  00_design_tables.R: sample design figures from the two final reports, used by
#  the check on clusters that could not be surveyed (02_decomposition.R, 3c).
#
#  alloc: clusters allocated to each state, urban and rural (Appendix A,
#         Table A.3, of NPC and ICF 2019 [FR359] and FMoHSW, NPC and ICF 2025
#         [FR395]). Clusters actually surveyed are counted from the data.
#  frame_hh_2024: residential households in each state, urban and rural, in the
#         2023 NPHC enumeration-area frame used by the 2023-24 survey (FR395,
#         Table A.1). The 2018 report gives the 2006 census population rather
#         than households, so the weight check is run for 2023-24 only.
#  Reported reasons: insecurity kept teams out of 11 clusters in 2018 and 20 in
#  2023-24. In 2018, 11 Borno LGAs holding about 39% of Borno households were
#  also dropped and their clusters replaced from the other 16 LGAs (FR359, A.3).
# =============================================================================

design_states <- c(
  "fctabuja", "benue", "kogi", "kwara", "nasarawa", "niger", "plateau",
  "adamawa", "bauchi", "borno", "gombe", "taraba", "yobe",
  "jigawa", "kaduna", "kano", "katsina", "kebbi", "sokoto", "zamfara",
  "abia", "anambra", "ebonyi", "enugu", "imo",
  "akwaibom", "bayelsa", "crossriver", "delta", "edo", "rivers",
  "ekiti", "lagos", "ogun", "ondo", "osun", "oyo")

alloc <- rbind(
  data.frame(round = "2018", state = design_states,
    urban = c(22, 5, 12, 24, 9, 10, 10,   10, 6, 13, 9, 6, 8,   5, 20, 22, 9, 6, 8, 7,
              8, 32, 30, 25, 18,   3, 9, 6, 18, 20, 19,   26, 50, 18, 17, 27, 30),
    rural = c(13, 33, 24, 11, 26, 28, 25,   25, 33, 25, 26, 29, 27,   34, 23, 31, 33, 29, 29, 29,
              28, 7, 6, 11, 21,   34, 26, 29, 20, 15, 22,   9, 3, 19, 19, 9, 12)),
  data.frame(round = "2023-24", state = design_states,
    urban = c(23, 12, 13, 16, 17, 15, 9,   7, 7, 29, 7, 10, 11,   5, 16, 19, 11, 8, 5, 8,
              19, 23, 13, 25, 13,   20, 19, 21, 32, 30, 34,   25, 49, 39, 24, 32, 35),
    rural = c(13, 29, 23, 21, 19, 25, 30,   29, 34, 9, 24, 25, 22,   31, 27, 26, 28, 28, 32, 27,
              16, 15, 21, 12, 26,   19, 10, 14, 8, 6, 8,   6, 2, 7, 13, 6, 8)))
# The 2018 state rows sum to 577 urban and 823 rural clusters (1,400 in all); the
# report's South West subtotal (171 urban, 68 rural) does not match its own rows.
stopifnot(sum(alloc$urban[alloc$round == "2018"] + alloc$rural[alloc$round == "2018"]) == 1400,
          sum(alloc$urban[alloc$round == "2023-24"]) == 701,
          sum(alloc$rural[alloc$round == "2023-24"]) == 699)

frame_hh_2024 <- data.frame(state = design_states,
  urban = c(684948, 528895, 392007, 514629, 551004, 631458, 347334,
            243394, 341540, 1071927, 143710, 303267, 274373,
            165352, 878427, 1169947, 484847, 272906, 169138, 257730,
            588309, 881956, 331472, 784841, 514353,
            745989, 305362, 572179, 1314029, 887174, 1572081,
            528627, 4140649, 2402970, 760710, 1156336, 1821459),
  rural = c(390287, 1382367, 740795, 701091, 616127, 1086781, 1253371,
            1019666, 1737304, 356534, 507825, 783176, 538375,
            1063317, 1566184, 1661023, 1225827, 970459, 1183068, 863561,
            519781, 583902, 567036, 395845, 1067329,
            753019, 173824, 380952, 316117, 196122, 385483,
            139144, 140000, 445414, 411622, 230430, 440541))
stopifnot(sum(frame_hh_2024$urban) == 28735329, sum(frame_hh_2024$rural) == 26793699)
