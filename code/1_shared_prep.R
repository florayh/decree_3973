#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Project: Bolivia Fires
# Purpose: Shared data prep for DID fire/land-use analysis scripts.
#          Loads grid-month panel, merges PSM weights, creates derived variables,
#          builds fire-season and annual subsets, and defines common outcome
#          vectors and helper functions. Sourced by scripts 2, 4, and 5.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(data.table)
library(fixest)
library(ggplot2)
library(tidyverse)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1. Paths -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_root <- "data"
tab_dir   <- file.path("output", "tables")
fig_dir   <- file.path("output", "figures")

dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2. Load and prepare panel -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("=== Loading data ===")
dt <- readRDS(file.path(data_root, "grid_month_panel_2km.rds"))
setDT(dt)


# Merge PSM weights
psm <- readRDS(file.path(data_root, "grid_psm_weights_2km.rds"))
setDT(psm)
dt <- merge(dt, psm[, .(ID, psm_weight, psm_matched)],
            by = "ID", all.x = TRUE)
rm(psm)

# Fill NAs for unmatched
dt[is.na(psm_weight), psm_weight := 0]
dt[is.na(psm_matched), psm_matched := FALSE]

# Derived variables
dt[, after_rescission := as.integer(post_decree == 1 & decree_period == 0)]
dt[, period := fifelse(post_decree == 0, "pre",
                fifelse(decree_period == 1, "during", "after"))]

message("Panel: ", format(nrow(dt), big.mark = ","), " rows")
message("PSM matched grids: ", uniqueN(dt$ID[dt$psm_matched == TRUE]))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3. Fire-season subset (Jul-Oct) -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt_fire <- dt[month %in% 7:10]

message("Fire-season panel: ", format(nrow(dt_fire), big.mark = ","),
        " rows (", uniqueN(dt_fire$ID), " grids x ", uniqueN(dt_fire$year_month), " months)")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4. Annual panel -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

build_annual_panel <- function(dt_input) {
  dt_yr <- dt_input[, .(
    # Land use (annual -- take first month)
    defor_area           = defor_area[1],
    primary_loss_km2     = primary_loss_km2[1],
    secondary_loss_km2   = secondary_loss_km2[1],
    secondary_growth_km2 = secondary_growth_km2[1],
    all_ag_sum           = all_ag_sum[1],
    soy_area_km          = soy_area_km[1],
    pasture              = pasture[1],
    forest               = forest[1],
    burned_area_km2       = burned_area_km2[1],
    # Climate (annual mean)
    temperature_2m   = mean(temperature_2m, na.rm = TRUE),
    chirps_precip_mm = mean(chirps_precip_mm, na.rm = TRUE),
    wind_speed_10m   = mean(wind_speed_10m, na.rm = TRUE),
    # Identifiers (static)
    treated          = treated[1],
    ADM1_ES          = ADM1_ES[1],
    ADM3_PCODE       = ADM3_PCODE[1],
    psm_matched      = psm_matched[1],
    psm_weight       = psm_weight[1],
    cluster_20km     = cluster_20km[1],
    elevation        = elevation[1],
    slope            = slope[1],
    pct_forest2014   = pct_forest2014[1],
    population       = population[1]
  ), by = .(ID, year)]

  # Treatment timing
  dt_yr[, during_year     := as.integer(year %in% c(2019, 2020))]
  dt_yr[, after_year      := as.integer(year >= 2021)]
  dt_yr[, post_year       := as.integer(year >= 2019)]
  dt_yr[, event_time_year := year - 2019]
  dt_yr[, period := fifelse(year < 2019, "pre",
                     fifelse(during_year == 1, "during", "after"))]

  return(dt_yr)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5. Constants -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Outcome vectors
fire_outcomes <- c("n_fires",
                   "n_fires_grassland", "n_fires_other_natural",
                   "n_fires_agricultural")
fire_labels   <- c("Total fires",
                   "Natural grass fires", "Forest/other veg. fires",
                   "Ag. land fires")

lu_outcomes <- c("primary_loss_km2", "all_ag_sum", "soy_area_km",
                 "pasture", "forest", "burned_area_km2")
lu_labels   <- c("Primary forest loss",
                 "Agriculture (total)", "Soy area",
                 "Pasture", "Forest", "Burned area")

# Climate control formula fragment
climate_controls <- "temperature_2m + chirps_precip_mm + wind_speed_10m"

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 6. Helper functions -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Extract event study coefficients from fixest model
extract_es <- function(model, ref_value = -1) {
  ct <- coeftable(model)
  idx <- grep("event_time", rownames(ct))
  if (length(idx) == 0) return(NULL)

  es <- data.table(
    term     = rownames(ct)[idx],
    estimate = ct[idx, "Estimate"],
    se       = ct[idx, "Std. Error"]
  )
  es[, ci_low  := estimate - 1.96 * se]
  es[, ci_high := estimate + 1.96 * se]
  es[, event_time := as.numeric(gsub(".*::(-?[0-9]+):.*", "\\1", term))]

  ref_row <- data.table(term = "ref", estimate = 0, se = 0,
                        ci_low = 0, ci_high = 0, event_time = ref_value)
  es <- rbind(es, ref_row)
  setorder(es, event_time)
  return(es)
}

# Event study ggplot
plot_es <- function(es_data, title = "", ylab = "Estimate",
                    xlab = "Periods relative to decree",
                    decree_time = 0, rescission_time = 16) {
  ggplot(es_data, aes(x = event_time, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "solid", color = "grey50") +
    geom_vline(xintercept = decree_time - 0.5, linetype = "dashed",
               color = "red", alpha = 0.7) +
    geom_vline(xintercept = rescission_time - 0.5, linetype = "dotted",
               color = "darkred", alpha = 0.5) +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15,
                fill = "#1f77b4") +
    geom_point(size = 1, color = "#1f77b4") +
    geom_line(linewidth = 0.4, color = "#1f77b4") +
    labs(title = title, x = xlab, y = ylab) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12))
}

# Compute dominant land tenure per grid (plurality, >= min_share threshold).
# Combines NPA and forest reserve into "protected" before determining plurality.
# Returns data.table(ID, dominant_tenure) for grids passing the threshold.
compute_dominant_tenure <- function(dt, min_share = 40) {
  tenure_source_cols <- c("pct_private_small", "pct_private_large", "pct_campesino",
                          "pct_indigenous", "pct_within_npa",
                          "pct_within_forest_reserve", "pct_tierra_fiscal")
  tenure_cols <- c("pct_private_small", "pct_private_large", "pct_campesino",
                   "pct_indigenous", "pct_protected", "pct_tierra_fiscal")

  grid_tenure <- dt[, .SD[1], by = ID, .SDcols = tenure_source_cols]

  stopifnot(!anyNA(grid_tenure[, .(pct_within_npa, pct_within_forest_reserve)]))
  grid_tenure[, pct_protected := pmin(100, pct_within_npa + pct_within_forest_reserve)]

  grid_tenure[, dominant_tenure := tenure_cols[max.col(.SD, ties.method = "first")],
              .SDcols = tenure_cols]
  grid_tenure[, dominant_tenure := gsub("pct_", "", dominant_tenure)]
  stopifnot(!any(grid_tenure$dominant_tenure %in% c("npa", "forest_reserve")))

  grid_tenure[, max_share := do.call(pmax, .SD), .SDcols = tenure_cols]

  n_before <- nrow(grid_tenure)
  grid_tenure <- grid_tenure[max_share >= min_share]
  message(sprintf("compute_dominant_tenure: %d / %d grids pass >= %d%% share threshold",
                  nrow(grid_tenure), n_before, min_share))

  grid_tenure[, .(ID, dominant_tenure)]
}

# DR-DID event study extractor (from aggte(type="dynamic"), bootstrap critical value)
extract_drdid_es <- function(agg_obj) {
  if (is.null(agg_obj)) return(NULL)
  cval <- agg_obj$crit.val.egt
  if (is.null(cval) || !is.finite(cval)) cval <- 1.96
  data.table(
    event_time = agg_obj$egt,
    estimate   = agg_obj$att.egt,
    se         = agg_obj$se.egt,
    ci_low     = agg_obj$att.egt - cval * agg_obj$se.egt,
    ci_high    = agg_obj$att.egt + cval * agg_obj$se.egt
  )
}

# Significance stars (p-value -> LaTeX star notation)
add_stars <- function(p_val) {
  if (is.na(p_val)) return("")
  if (p_val < 0.001) return("***")
  if (p_val < 0.01)  return("**")
  if (p_val < 0.05)  return("*")
  if (p_val < 0.1)   return("$^{+}$")
  return("")
}

message("=== Shared prep complete ===")
