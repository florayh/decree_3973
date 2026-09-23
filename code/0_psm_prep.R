#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Project: Bolivia Fires
# Purpose: PSM matching for DID fire analysis. Computes pre-treatment grid
#          characteristics, runs Propensity Score Matching to construct a
#          comparable control group, and produces the pre-trends figure.
#
# Inputs:
#   - data/grid_month_panel_2km.rds
#
# Outputs:
#   - data/grid_psm_weights_2km.rds  (PSM match indicators + weights)
#   - output/figures/did_pretrends_fire.pdf
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(data.table)
library(MatchIt)
library(ggplot2)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1. Load data -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_root <- "data"

message("=== Loading grid panel ===")
dt <- readRDS(file.path(data_root, "grid_month_panel_2km.rds"))
setDT(dt)
message("Panel: ", format(nrow(dt), big.mark = ","), " rows x ", ncol(dt), " cols")
message("Grids: ", uniqueN(dt$ID), " | Months: ", uniqueN(dt$year_month))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2. Compute pre-treatment grid characteristics -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== Computing pre-treatment characteristics ===")

# Pre-treatment: all months before decree (event_time < 0 = before July 2019)
pre <- dt[event_time < 0]

# Fire-related and climate means over pre-treatment period
fire_chars <- pre[, .(
  mean_fires_pre     = mean(n_fires, na.rm = TRUE),
  pct_months_fire    = mean(n_fires > 0, na.rm = TRUE), # mean of a logical value: proportion of months that have one fire
  mean_precip_pre    = mean(chirps_precip_mm, na.rm = TRUE),
  mean_temp_pre      = mean(temperature_2m, na.rm = TRUE)
), by = ID]

# Static covariates (time-invariant -- first observation per grid)
static <- dt[, .(
  treated              = treated[1],
  elevation            = elevation[1],
  slope                = slope[1],
  pct_forest2014       = pct_forest2014[1],
  population           = population[1],
  pct_private_small    = pct_private_small[1],
  pct_private_large    = pct_private_large[1],
  pct_campesino        = pct_campesino[1],
  pct_indigenous       = pct_indigenous[1],
  pct_military         = pct_military[1],
  pct_within_npa       = pct_within_npa[1],
  pct_within_forest_reserve = pct_within_forest_reserve[1],
  pct_tierra_fiscal    = pct_tierra_fiscal[1],
  ADM1_ES              = ADM1_ES[1],
  ADM3_PCODE           = ADM3_PCODE[1],
  centroid_lon         = centroid_lon[1],
  centroid_lat         = centroid_lat[1]
), by = ID]

grid_chars <- merge(fire_chars, static, by = "ID")

# Elevation filter already applied in bundled panel (controls <= 1800m).

# Dominant tenure type (highest share among 6 analytical categories; exclude military)
# Protected combines NPA and forest reserve coverage before classification.
# Must exceed 40% of grid to be considered dominant; grids below threshold are dropped
stopifnot(!anyNA(grid_chars[, .(pct_within_npa, pct_within_forest_reserve)]))
grid_chars[, pct_protected := pmin(100,
  pct_within_npa + pct_within_forest_reserve)]
tenure_cols <- c("pct_private_small", "pct_private_large", "pct_campesino",
                 "pct_indigenous", "pct_protected", "pct_tierra_fiscal")
grid_chars[, dominant_tenure := tenure_cols[max.col(.SD, ties.method = "first")],
           .SDcols = tenure_cols]
grid_chars[, dominant_tenure := gsub("pct_", "", dominant_tenure)]
grid_chars[, max_tenure_share := do.call(pmax, .SD), .SDcols = tenure_cols]
stopifnot(!any(grid_chars$dominant_tenure %in% c("npa", "forest_reserve")))

n_before_tenure <- nrow(grid_chars)
grid_chars <- grid_chars[max_tenure_share >= 40]
message("Tenure share filter (>= 40%): dropped ", n_before_tenure - nrow(grid_chars),
        " grids (", nrow(grid_chars), " remain)")

grid_chars[, dominant_tenure := factor(dominant_tenure)]

message("Treatment distribution:")
print(table(Treated = grid_chars$treated))
message("\nDominant tenure x treatment:")
print(table(grid_chars$dominant_tenure, grid_chars$treated, dnn = c("Tenure", "Treated")))


# Drop rows with NA in matching variables
match_vars <- c("elevation", "slope", "pct_forest2014", "population",
                 "mean_precip_pre", "mean_temp_pre", "dominant_tenure")
n_na <- sum(!complete.cases(grid_chars[, ..match_vars]))
if (n_na > 0) message("Dropping ", n_na, " grids with NA in matching variables")
grid_chars_complete <- grid_chars[complete.cases(grid_chars[, ..match_vars])]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3b. PSM matching -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== PSM matching ===")

# Nearest-neighbor 1:1 matching with caliper = 0.2 SD of logit propensity score
m_psm <- matchit(
  treated ~ elevation + slope + pct_forest2014 + population +
    mean_precip_pre + mean_temp_pre + dominant_tenure,
  data = grid_chars_complete,
  method = "nearest",
  distance = "glm",
  caliper = 0.2,
  ratio = 1
)

message("\nPSM summary:")
print(summary(m_psm))

# Extract matched data with PSM weights
psm_data <- match.data(m_psm)
psm_matched_ids <- psm_data$ID

message("\nPSM matched sample: ", nrow(psm_data), " grids")
message("  PSM treated: ", sum(psm_data$treated == 1))
message("  PSM control: ", sum(psm_data$treated == 0))
message("  Dropped: ", nrow(grid_chars_complete) - nrow(psm_data))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5. Pre-trend visualization -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== Pre-trend visualization ===")

# Fire season subset (Jul-Oct, matching DR-DID sample); elevation filter pre-applied
dt_lowland <- dt
dt_lowland_fire <- dt_lowland[month %in% 7:10 & year %in% 2015:2022]
trends <- dt_lowland_fire[, .(mean_fires = mean(n_fires, na.rm = TRUE)),
                           by = .(year_month, year, month, treated)]
trends[, treated_label := fifelse(treated == 1,
                                   "Treated (Santa Cruz + Beni)",
                                   "Control")]

# Sequential index for equal spacing (no gaps between Oct and next Jul)
setorder(trends, treated, year_month)
month_seq <- sort(unique(trends$year_month))
trends[, seq_idx := match(year_month, month_seq)]

# X-axis labels: show year at each July
jul_idx <- trends[month == 7 & treated == 1, .(seq_idx, year)]

# Decree line: between Oct 2018 (last pre-treatment) and Jul 2019 (first treated)
decree_idx <- mean(c(trends[year_month == as.Date("2018-10-01") & treated == 1, seq_idx],
                     trends[year_month == as.Date("2019-07-01") & treated == 1, seq_idx]))

p_pretrends <- ggplot(trends, aes(x = seq_idx, y = mean_fires,
                                   color = treated_label)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.2) +
  geom_vline(xintercept = decree_idx, linetype = "dashed",
             color = "red", alpha = 0.7) +
  annotate("text", x = decree_idx + 0.3, y = Inf, label = "Decree 3973",
           hjust = 0, vjust = 1.8, size = 2.8, color = "red") +
  scale_x_continuous(breaks = jul_idx$seq_idx, labels = jul_idx$year) +
  scale_color_manual(values = c("Treated (Santa Cruz + Beni)" = "#d62728",
                                "Control" = "#1f77b4")) +
  labs(x = NULL, y = "Mean monthly fire count per grid cell",
       color = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

dir.create(file.path("output", "figures"), recursive = TRUE, showWarnings = FALSE)
ggsave(file.path("output", "figures", "did_pretrends_fire.pdf"),
       p_pretrends, width = 8, height = 4, device = cairo_pdf)
message("Pre-trends figure exported: output/figures/did_pretrends_fire.pdf")

# PSM match status (used by Section 7b save)
grid_chars[, psm_matched := ID %in% psm_matched_ids]
grid_chars[, psm_match_status := fifelse(
  !psm_matched, "Dropped",
  fifelse(treated == 1, "Matched treated", "Matched control")
)]


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 7b. Save PSM weights -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== Saving PSM weights ===")

psm_weights <- data.table(
  ID         = psm_data$ID,
  psm_weight = psm_data$weights,
  psm_subclass = as.character(psm_data$subclass)
)

psm_out <- merge(
  grid_chars[, .(ID, treated, mean_fires_pre, pct_months_fire,
                 mean_precip_pre,
                 dominant_tenure, psm_matched, psm_match_status)],
  psm_weights,
  by = "ID", all.x = TRUE
)
psm_out[is.na(psm_weight), psm_weight := 0]

saveRDS(psm_out, file.path("data", "grid_psm_weights_2km.rds"))
message("PSM weights saved: ", nrow(psm_out), " grids")
message("  Matched (weight > 0): ", sum(psm_out$psm_weight > 0))
message("  Unmatched (weight = 0): ", sum(psm_out$psm_weight == 0))

message("\n=== Script 0 complete ===")
