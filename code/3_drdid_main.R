#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Project: Bolivia Fires
# Purpose: Doubly Robust DID (Sant'Anna & Zhao 2020) as robustness check.
#          Uses the `did` package (Callaway & Sant'Anna 2021) with est_method = "dr" with symmetric
#          propensity score trimming (Crump et al. 2009, threshold [0.01, 0.99])
#          to enforce common support. After elevation <= 1800m filter, control
#          grids with PS < 0.01 and treated grids with PS > 0.99 are dropped. Fire outcomes use fire-season-only panel (Jul-Oct, 2015-2022).
#          Land use outcomes use annual panel (2015-2022). Produces event study
#          figures and DR-DID summary tables.
#
# Reference: Sant'Anna, P.H.C. and Zhao, J. (2020). "Doubly Robust
#            Difference-in-Differences Estimators." Journal of Econometrics.
#            Callaway, B. and Sant'Anna, P.H.C. (2021). "Difference-in-
#            Differences with Multiple Time Periods." Journal of Econometrics.
#
# Inputs:
#   - data/grid_month_panel_2km.rds
#   - data/grid_psm_weights_2km.rds (loaded by 1_shared_prep.R)
#
# Outputs:
#   Figures (output/figures/):
#     - did_drdid_eventstudy_n_fires.pdf
#     - did_drdid_eventstudy_{primary_loss_km2,all_ag_sum,soy_area_km}.pdf
#     - did_drdid_trimmed_sample.pdf
#   Tables (output/tables/):
#     - did_balance_ps_trimmed.tex
#     - did_drdid_fire.tex
#     - did_drdid_landuse.tex
#
# RUNTIME NOTE: DR-DID estimation with bootstrap SEs is slow.
#   Expect several hours for the full run (fire + land use outcomes).
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

source("code/1_shared_prep.R")
library(did)
library(sf)
library(kableExtra)

set.seed(20260414)

message("Full panel: ", format(nrow(dt), big.mark = ","), " rows")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1-PS. Propensity score trimming (Crump et al. 2009) -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== Propensity score trimming ===")

# Elevation filter already applied in bundled panel (controls <= 1800m).

# Collapse pre-treatment climate covariates to one row per grid
grid_covs <- dt[event_time < 0, .(
  mean_temp_pre   = mean(temperature_2m, na.rm = TRUE),
  mean_precip_pre = mean(chirps_precip_mm, na.rm = TRUE)
), by = ID]

# Collapse static covariates to one row per grid
grid_static <- dt[, .(
  treated        = treated[1],
  elevation      = elevation[1],
  slope          = slope[1],
  pct_forest2014 = pct_forest2014[1],
  log_pop        = log1p(population[1])
), by = ID]

grid_ps <- merge(grid_static, grid_covs, by = "ID", all.x = TRUE)

# Estimate propensity score via logistic regression
ps_model <- glm(treated ~ elevation + slope + pct_forest2014 + log_pop +
                   mean_precip_pre + mean_temp_pre,
                 data = grid_ps, family = binomial(link = "logit"))

grid_ps[, ps := predict(ps_model, newdata = grid_ps, type = "response")]

message("PS distribution by treatment status:")
message("  Treated:  ", paste(sprintf("%.3f", quantile(grid_ps[treated == 1]$ps,
        probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = TRUE)),
        collapse = " | "))
message("  Control:  ", paste(sprintf("%.3f", quantile(grid_ps[treated == 0]$ps,
        probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = TRUE)),
        collapse = " | "))
message("  (Quantiles: 0% | 5% | 25% | 50% | 75% | 95% | 100%)")

# Symmetric trimming (Crump et al. 2009): drop units outside [0.01, 0.99]
ps_threshold <- 0.01

# Drop control grids with PS < 0.01 or NA
n_control_before <- sum(grid_ps$treated == 0)
n_na_ps <- sum(is.na(grid_ps[treated == 0]$ps))
if (n_na_ps > 0) message("  ", n_na_ps, " control grids have NA PS (missing covariates) -- dropped")
drop_ctrl <- grid_ps[treated == 0 & (is.na(ps) | ps < ps_threshold), ID]

message(sprintf("PS trimming controls (PS < %.2f): dropped %d / %d control grids (%.1f%%)",
                ps_threshold, length(drop_ctrl), n_control_before,
                100 * length(drop_ctrl) / n_control_before))

# Drop treated grids with PS > 1 - threshold (= 0.99)
n_treated_before <- sum(grid_ps$treated == 1)
drop_treat <- grid_ps[treated == 1 & (is.na(ps) | ps > (1 - ps_threshold)), ID]

message(sprintf("PS trimming treated (PS > %.2f): dropped %d / %d treated grids (%.1f%%)",
                1 - ps_threshold, length(drop_treat), n_treated_before,
                100 * length(drop_treat) / n_treated_before))

drop_ids <- c(drop_ctrl, drop_treat)
n_remaining <- uniqueN(grid_ps$ID) - length(drop_ids)
message("Remaining after symmetric trimming: ", n_remaining, " grids (",
        n_treated_before - length(drop_treat), " treated + ",
        n_control_before - length(drop_ctrl), " control)")

# Collapse to grid-level pre-treatment means for balance table
grid_bal_static <- dt[, .(
  treated          = treated[1],
  elevation        = elevation[1],
  slope            = slope[1],
  pct_forest2014   = pct_forest2014[1],
  log_pop          = log1p(population[1]),
  pct_private_small       = pct_private_small[1],
  pct_private_large       = pct_private_large[1],
  pct_campesino           = pct_campesino[1],
  pct_indigenous          = pct_indigenous[1],
  pct_within_npa          = pct_within_npa[1],
  pct_within_forest_reserve = pct_within_forest_reserve[1],
  pct_tierra_fiscal       = pct_tierra_fiscal[1]
), by = ID]

grid_bal_timevar <- dt[year < 2019, .(
  mean_precip  = mean(chirps_precip_mm, na.rm = TRUE),
  mean_wind    = mean(wind_speed_10m, na.rm = TRUE),
  mean_fires   = mean(n_fires, na.rm = TRUE),
  mean_primary_loss  = mean(primary_loss_km2, na.rm = TRUE),
  mean_all_ag        = mean(all_ag_sum, na.rm = TRUE),
  mean_soy           = mean(soy_area_km, na.rm = TRUE),
  mean_pasture       = mean(pasture, na.rm = TRUE),
  mean_forest        = mean(forest, na.rm = TRUE),
  mean_burned_area   = mean(burned_area_km2, na.rm = TRUE),
  mean_defor         = mean(defor_area, na.rm = TRUE)
), by = ID]

grid_bal <- merge(grid_bal_static, grid_bal_timevar, by = "ID", all.x = TRUE)
grid_bal[, trimmed := ID %in% drop_ids]
rm(grid_bal_static, grid_bal_timevar)

# Filter panel to surviving grids
dt <- dt[!(ID %in% drop_ids)]

# Show trimmed PS distribution
trimmed_ps <- grid_ps[!(ID %in% drop_ids)]
message("Trimmed PS distribution:")
message("  Treated:  ", paste(sprintf("%.3f", quantile(trimmed_ps[treated == 1]$ps,
        probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = TRUE)),
        collapse = " | "))
message("  Control:  ", paste(sprintf("%.3f", quantile(trimmed_ps[treated == 0]$ps,
        probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = TRUE)),
        collapse = " | "))

# Store trim status for map before cleaning up
grid_ps[, trim_status := fifelse(
  ID %in% drop_ids, "Trimmed",
  fifelse(treated == 1, "Retained treated", "Retained control")
)]

# Clean up PS objects (keep grid_ps for map)
rm(grid_covs, grid_static, trimmed_ps, drop_ctrl, drop_treat)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1-MAP. Map of PS-trimmed sample -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== Map of PS-trimmed sample ===")

grid_sf <- st_read(file.path("data", "shapefiles/grids/grid_bolivia_2km.shp"),
                   quiet = TRUE)
dept_sf <- st_read(file.path("data", "shapefiles/adm1",
                              "bol_admbnda_adm1_gov_2020514.shp"),
                   quiet = TRUE)

hex_map <- merge(grid_sf, grid_ps[, .(ID, trim_status)],
                 by = "ID", all.x = FALSE)
hex_map$trim_status <- factor(hex_map$trim_status,
                               levels = c("Retained treated", "Retained control", "Trimmed"))

p_trim_map <- ggplot() +
  geom_sf(data = hex_map, aes(fill = trim_status), color = NA, linewidth = 0) +
  geom_sf(data = dept_sf, fill = NA, color = "grey40", linewidth = 0.3) +
  scale_fill_manual(
    values = c("Retained treated" = "#d62728",
               "Retained control" = "#1f77b4",
               "Trimmed"          = "grey85"),
    name = NULL
  ) +
  labs(title = "DR-DID Sample: Symmetric PS Trimming [0.01, 0.99]") +
  theme_minimal(base_size = 11) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom")

ggsave(file.path(fig_dir, "did_drdid_trimmed_sample.pdf"),
       p_trim_map, width = 6, height = 8, device = cairo_pdf)
message("Map exported: ", file.path(fig_dir, "did_drdid_trimmed_sample.pdf"))

rm(grid_sf, dept_sf, hex_map, p_trim_map, grid_ps, ps_model, drop_ids)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1-BAL. Post-trimming balance table -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== Post-trimming balance table ===")

# SMD function
compute_smd <- function(data, var, treat_var = "treated") {
  t_idx <- data[[treat_var]] == 1
  c_idx <- data[[treat_var]] == 0
  x_t <- data[[var]][t_idx]
  x_c <- data[[var]][c_idx]
  m_t <- mean(x_t, na.rm = TRUE)
  m_c <- mean(x_c, na.rm = TRUE)
  sd_pool <- sqrt((var(x_t, na.rm = TRUE) + var(x_c, na.rm = TRUE)) / 2)
  smd <- ifelse(sd_pool > 0, (m_t - m_c) / sd_pool, 0)
  return(c(mean_t = m_t, mean_c = m_c, smd = smd))
}

# Balance variables and labels
bal_vars <- c(
  "elevation", "slope", "pct_forest2014", "log_pop",
  "mean_precip", "mean_wind",
  "pct_private_small", "pct_private_large", "pct_campesino",
  "pct_indigenous", "pct_within_npa", "pct_within_forest_reserve",
  "pct_tierra_fiscal",
  "mean_fires",
  "mean_primary_loss", "mean_all_ag", "mean_soy",
  "mean_pasture", "mean_forest", "mean_burned_area", "mean_defor"
)

bal_labels <- c(
  "Elevation (m)", "Slope (degrees)", "Forest cover 2014 (\\%)",
  "log(1 + Population)",
  "Precipitation (mm)", "Wind speed (m/s)",
  "Private small (\\%)", "Private large (\\%)",
  "Campesino (\\%)", "Indigenous (\\%)",
  "Protected area (\\%)", "Forest reserve (\\%)",
  "Tierra fiscal (\\%)",
  "Fire count",
  "Primary forest loss (km\\textsuperscript{2})",
  "Agriculture total (km\\textsuperscript{2})",
  "Soy area (km\\textsuperscript{2})",
  "Pasture (km\\textsuperscript{2})",
  "Forest (km\\textsuperscript{2})",
  "Burned area (km\\textsuperscript{2})",
  "Deforestation (km\\textsuperscript{2})"
)

# Before trimming: all grids
bal_before <- t(sapply(bal_vars, function(v) compute_smd(grid_bal, v)))

# After trimming: exclude trimmed grids
bal_after <- t(sapply(bal_vars, function(v)
  compute_smd(grid_bal[trimmed == FALSE], v)))

bal_table <- data.frame(
  Variable      = bal_labels,
  Treat_Before  = round(bal_before[, "mean_t"], 3),
  Ctrl_Before   = round(bal_before[, "mean_c"], 3),
  SMD_Before    = round(bal_before[, "smd"], 3),
  Treat_After   = round(bal_after[, "mean_t"], 3),
  Ctrl_After    = round(bal_after[, "mean_c"], 3),
  SMD_After     = round(bal_after[, "smd"], 3),
  row.names     = NULL
)

message("\nBalance table (before vs. after PS trimming):")
print(bal_table)

# Export LaTeX
kbl_out <- kbl(bal_table, format = "latex", booktabs = TRUE,
               col.names = c("", "Treated", "Control", "SMD",
                              "Treated", "Control", "SMD"),
               escape = FALSE, align = c("l", rep("c", 6)),
               caption = "Covariate balance before and after propensity score trimming",
               label = "did_balance_ps_trimmed") |>
  kable_styling(latex_options = c("scale_down", "hold_position"),
                font_size = 9) |>
  add_header_above(c(" " = 1,
                      "Before PS trimming" = 3,
                      "After PS trimming" = 3)) |>
  footnote(general = "Pre-treatment means for treated and control grids. SMD = standardized mean difference (pooled SD). Sample starts with elevation $\\\\leq$ 1800m filter on controls, then symmetric PS trimming [0.01, 0.99] (Crump et al.\\\\ 2009). Time-varying covariates averaged over the pre-treatment period (2015--2018). Climate and fire variables are monthly; land cover variables are annual (constant within year).",
           general_title = "",
           escape = FALSE,
           threeparttable = TRUE)

writeLines(kbl_out, file.path(tab_dir, "did_balance_ps_trimmed.tex"))
message("Balance table exported: ", file.path(tab_dir, "did_balance_ps_trimmed.tex"))

rm(grid_bal, bal_before, bal_after, bal_table, kbl_out)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1a. Fire panel: fire season (Jul-Oct), 2015-2022 -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== Building fire season panel ===")

dt_fire <- dt[month %in% 7:10 & year %in% 2015:2022]

fire_month_indices <- sort(unique(dt_fire$month_index))
dt_fire[, fire_period := match(month_index, fire_month_indices)]

treatment_fire_period <- match(54L, fire_month_indices)
message("Treatment fire_period: ", treatment_fire_period,
        " (month_index 54 = Jul 2019)")

dt_fire[, treatment_period := fifelse(treated == 1, as.numeric(treatment_fire_period), 0)]

dt_fire[, log_pop := log1p(population)]
dt_fire[, month_f := as.factor(month)]

# Check panel balance
n_grids_fire <- uniqueN(dt_fire$ID)
n_periods_fire <- uniqueN(dt_fire$month_index)
message("Fire panel: ", n_grids_fire, " grids x ", n_periods_fire,
        " fire-season months = ", format(nrow(dt_fire), big.mark = ","), " obs")
message("Expected balanced: ", format(n_grids_fire * n_periods_fire, big.mark = ","))

# Verify balance
obs_per_grid <- dt_fire[, .N, by = ID]
if (all(obs_per_grid$N == n_periods_fire)) {
  message("Panel is balanced")
} else {
  message("WARNING: Panel is UNBALANCED. ",
          sum(obs_per_grid$N != n_periods_fire), " grids have incomplete obs")
  complete_grids <- obs_per_grid[N == n_periods_fire, ID]
  dt_fire <- dt_fire[ID %in% complete_grids]
  message("After balancing: ", uniqueN(dt_fire$ID), " grids x ",
          n_periods_fire, " periods")
}

message("Treated grids: ", uniqueN(dt_fire$ID[dt_fire$treated == 1]))
message("Control grids: ", uniqueN(dt_fire$ID[dt_fire$treated == 0]))
message("Month_index values: ", paste(sort(unique(dt_fire$month_index)), collapse = ", "))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1b. Land use panel: annual 2015-2022 -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== Building annual land use panel ===")

dt_year <- build_annual_panel(dt[year %in% 2015:2022])

dt_year[, treatment_year := fifelse(treated == 1, 2019, 0)]

dt_year[, log_pop := log1p(population)]

n_grids_yr <- uniqueN(dt_year$ID)
n_years <- uniqueN(dt_year$year)
message("Land use panel: ", n_grids_yr, " grids x ", n_years,
        " years = ", format(nrow(dt_year), big.mark = ","), " obs")
message("Treated: ", uniqueN(dt_year$ID[dt_year$treated == 1]),
        " | Control: ", uniqueN(dt_year$ID[dt_year$treated == 0]))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2. DR-DID: Fire outcomes -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== DR-DID: Fire outcomes ===")

xformla_fire <- ~ elevation + slope + pct_forest2014 + log_pop
stopifnot(identical(all.vars(xformla_fire),
                    c("elevation", "slope", "pct_forest2014", "log_pop")))
message("DR-DID covariates: ", paste(all.vars(xformla_fire), collapse = ", "))

df_fire <- as.data.frame(dt_fire)

fire_attgt <- list()
fire_agg_dyn <- list()
fire_agg_grp <- list()

for (i in seq_along(fire_outcomes)) {
  y <- fire_outcomes[i]
  message("  Estimating att_gt for: ", y)

  fire_attgt[[y]] <- tryCatch(
    att_gt(
      yname         = y,
      tname         = "fire_period",
      idname        = "ID",
      gname         = "treatment_period",
      xformla       = xformla_fire,
      control_group = "nevertreated",
      anticipation  = 0,
      est_method    = "dr",
      clustervars   = "ADM3_PCODE",
      data          = df_fire,
      print_details = FALSE
    ),
    error = function(e) {
      message("    ERROR: ", e$message)
      NULL
    }
  )

  if (!is.null(fire_attgt[[y]])) {
    message("    att_gt succeeded: ", length(fire_attgt[[y]]$att), " ATT(g,t) estimates")

    fire_agg_dyn[[y]] <- tryCatch(
      aggte(fire_attgt[[y]], type = "dynamic"),
      error = function(e) { message("    Dynamic agg error: ", e$message); NULL }
    )

    fire_agg_grp[[y]] <- tryCatch(
      aggte(fire_attgt[[y]], type = "group"),
      error = function(e) { message("    Group agg error: ", e$message); NULL }
    )

    if (!is.null(fire_agg_grp[[y]])) {
      message("    Overall ATT: ", round(fire_agg_grp[[y]]$overall.att, 4),
              " (SE: ", round(fire_agg_grp[[y]]$overall.se, 4), ")")
    }
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3. DR-DID: Land use outcomes -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== DR-DID: Land use outcomes ===")

xformla_lu <- ~ elevation + slope + pct_forest2014 + log_pop
stopifnot(identical(all.vars(xformla_lu),
                    c("elevation", "slope", "pct_forest2014", "log_pop")))

df_year <- as.data.frame(dt_year)

lu_attgt <- list()
lu_agg_dyn <- list()
lu_agg_grp <- list()

for (i in seq_along(lu_outcomes)) {
  y <- lu_outcomes[i]
  message("  Estimating att_gt for: ", y)

  lu_attgt[[y]] <- tryCatch(
    att_gt(
      yname         = y,
      tname         = "year",
      idname        = "ID",
      gname         = "treatment_year",
      xformla       = xformla_lu,
      control_group = "nevertreated",
      anticipation  = 0,
      est_method    = "dr",
      clustervars   = "ADM3_PCODE",
      data          = df_year,
      print_details = FALSE
    ),
    error = function(e) {
      message("    ERROR: ", e$message)
      NULL
    }
  )

  if (!is.null(lu_attgt[[y]])) {
    message("    att_gt succeeded: ", length(lu_attgt[[y]]$att), " ATT(g,t) estimates")

    lu_agg_dyn[[y]] <- tryCatch(
      aggte(lu_attgt[[y]], type = "dynamic"),
      error = function(e) { message("    Dynamic agg error: ", e$message); NULL }
    )

    lu_agg_grp[[y]] <- tryCatch(
      aggte(lu_attgt[[y]], type = "group"),
      error = function(e) { message("    Group agg error: ", e$message); NULL }
    )

    if (!is.null(lu_agg_grp[[y]])) {
      message("    Overall ATT: ", round(lu_agg_grp[[y]]$overall.att, 4),
              " (SE: ", round(lu_agg_grp[[y]]$overall.se, 4), ")")
    }
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4. Event study figures -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== Plotting DR-DID event studies ===")

plot_drdid_es <- function(es_data, title = "", ylab = "Estimate",
                          xlab = "Periods relative to treatment",
                          vline_pos = -0.5) {
  if (is.null(es_data) || nrow(es_data) == 0) return(NULL)

  ggplot(es_data, aes(x = event_time, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "solid", color = "grey50") +
    geom_vline(xintercept = vline_pos, linetype = "dashed",
               color = "red", alpha = 0.7) +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15,
                fill = "#1f77b4") +
    geom_point(size = 1.5, color = "#1f77b4") +
    geom_line(linewidth = 0.4, color = "#1f77b4") +
    labs(title = title, x = xlab, y = ylab) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12))
}

# Fire event studies
for (i in seq_along(fire_outcomes)) {
  y <- fire_outcomes[i]
  es_data <- extract_drdid_es(fire_agg_dyn[[y]])
  if (is.null(es_data)) {
    message("  Skipping ", y, " -- no event study data")
    next
  }

  p <- plot_drdid_es(
    es_data,
    title = paste0("DR-DID Event Study (PS-Trimmed): ", fire_labels[i]),
    ylab = paste0("ATT on ", tolower(fire_labels[i])),
    xlab = "Fire-season months relative to Jul 2019"
  )

  out_file <- file.path(fig_dir, paste0("did_drdid_eventstudy_", y, ".pdf"))
  ggsave(out_file, p, width = 9, height = 5, device = cairo_pdf)
  message("  Exported: ", out_file)
}

# Land use event studies
for (i in seq_along(lu_outcomes)) {
  y <- lu_outcomes[i]
  es_data <- extract_drdid_es(lu_agg_dyn[[y]])
  if (is.null(es_data)) {
    message("  Skipping ", y, " -- no event study data")
    next
  }

  p <- plot_drdid_es(
    es_data,
    title = paste0("DR-DID Event Study (PS-Trimmed): ", lu_labels[i]),
    ylab = paste0("ATT on ", tolower(lu_labels[i])),
    xlab = "Years relative to 2019"
  )

  out_file <- file.path(fig_dir, paste0("did_drdid_eventstudy_", y, ".pdf"))
  ggsave(out_file, p, width = 7, height = 5, device = cairo_pdf)
  message("  Exported: ", out_file)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5. DR-DID summary tables -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== Building DR-DID summary tables ===")

# --- Fire table ---
fire_tex <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\small",
  "\\caption{Doubly Robust DID: Fire Outcomes (PS-Trimmed, Fire Season)}",
  "\\label{tab:did_drdid_fire}",
  "\\begin{tabular}{lccc}",
  "\\toprule",
  " & DR-DID & Treated mean (pre) & $N$ \\\\",
  "\\midrule"
)

for (i in seq_along(fire_outcomes)) {
  y <- fire_outcomes[i]

  if (!is.null(fire_agg_grp[[y]])) {
    dr_att <- fire_agg_grp[[y]]$overall.att
    dr_se  <- fire_agg_grp[[y]]$overall.se
    dr_p <- 2 * pnorm(-abs(dr_att / dr_se))
    dr_cell <- paste0(sprintf("%.3f", dr_att), add_stars(dr_p))
    dr_se_cell <- paste0("(", sprintf("%.3f", dr_se), ")")
  } else {
    dr_cell <- "---"; dr_se_cell <- ""
  }

  treat_mean <- sprintf("%.3f",
    mean(dt_fire[treated == 1 & event_time < 0][[y]], na.rm = TRUE))
  n_obs <- nrow(dt_fire)

  fire_tex <- c(fire_tex,
    paste0(fire_labels[i], " & ", dr_cell,
           " & ", treat_mean,
           " & ", format(n_obs, big.mark = ","), " \\\\"),
    paste0(" & ", dr_se_cell, " & & \\\\")
  )
  if (i < length(fire_outcomes)) fire_tex <- c(fire_tex, "[0.5em]")
}

fire_tex <- c(fire_tex,
  "\\midrule",
  "Covariates in OR/PS & Yes & & \\\\",
  "Sample & PS-trimmed & & \\\\",
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}[flushleft]",
  "\\small",
  "\\item \\textit{Notes:} DR-DID = Sant'Anna \\& Zhao (2020) doubly robust estimator, implemented via the \\texttt{did} package (Callaway \\& Sant'Anna, 2021); ATT aggregated across post-treatment periods. The trimming propensity score includes elevation, slope, 2014 forest share, log population, and pre-treatment mean precipitation and temperature. The DR-DID propensity-score and outcome-regression models include elevation, slope, 2014 forest share, and log population. Sample: control grids with elevation $>$ 1800m dropped, then symmetric PS trimming (Crump et al.\\ 2009): controls with PS $<$ 0.01 and treated with PS $>$ 0.99 dropped. Fire season months (Jul--Oct) only, 2015--2022. SE clustered at municipality level (bootstrap). $^{+}p<0.1$, $^{*}p<0.05$, $^{**}p<0.01$, $^{***}p<0.001$.",
  "\\end{tablenotes}",
  "\\end{table}"
)

writeLines(fire_tex, file.path(tab_dir, "did_drdid_fire.tex"))
message("Fire DR-DID table exported: ", file.path(tab_dir, "did_drdid_fire.tex"))

# --- Land use table ---
lu_tex <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\small",
  "\\caption{Doubly Robust DID: Land Use Outcomes (PS-Trimmed, Annual)}",
  "\\label{tab:did_drdid_landuse}",
  "\\begin{tabular}{lccc}",
  "\\toprule",
  " & DR-DID & Treated mean (pre) & $N$ \\\\",
  "\\midrule"
)

for (i in seq_along(lu_outcomes)) {
  y <- lu_outcomes[i]

  if (!is.null(lu_agg_grp[[y]])) {
    dr_att <- lu_agg_grp[[y]]$overall.att
    dr_se  <- lu_agg_grp[[y]]$overall.se
    dr_p <- 2 * pnorm(-abs(dr_att / dr_se))
    dr_cell <- paste0(sprintf("%.3f", dr_att), add_stars(dr_p))
    dr_se_cell <- paste0("(", sprintf("%.3f", dr_se), ")")
  } else {
    dr_cell <- "---"; dr_se_cell <- ""
  }

  treat_mean <- sprintf("%.3f",
    mean(dt_year[treated == 1 & year < 2019][[y]], na.rm = TRUE))
  n_obs <- nrow(dt_year)

  lu_tex <- c(lu_tex,
    paste0(lu_labels[i], " & ", dr_cell,
           " & ", treat_mean,
           " & ", format(n_obs, big.mark = ","), " \\\\"),
    paste0(" & ", dr_se_cell, " & & \\\\")
  )
  if (i < length(lu_outcomes)) lu_tex <- c(lu_tex, "[0.5em]")
}

lu_tex <- c(lu_tex,
  "\\midrule",
  "Covariates in OR/PS & Yes & & \\\\",
  "Sample & PS-trimmed & & \\\\",
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}[flushleft]",
  "\\small",
  "\\item \\textit{Notes:} DR-DID = Sant'Anna \\& Zhao (2020) doubly robust estimator, implemented via the \\texttt{did} package (Callaway \\& Sant'Anna, 2021); ATT aggregated across post-treatment periods. The trimming propensity score includes elevation, slope, 2014 forest share, log population, and pre-treatment mean precipitation and temperature. The DR-DID propensity-score and outcome-regression models include elevation, slope, 2014 forest share, and log population. Sample: control grids with elevation $>$ 1800m dropped, then symmetric PS trimming (Crump et al.\\ 2009): controls with PS $<$ 0.01 and treated with PS $>$ 0.99 dropped. Annual panel 2015--2022. SE clustered at municipality level (bootstrap). $^{+}p<0.1$, $^{*}p<0.05$, $^{**}p<0.01$, $^{***}p<0.001$.",
  "\\end{tablenotes}",
  "\\end{table}"
)

writeLines(lu_tex, file.path(tab_dir, "did_drdid_landuse.tex"))
message("Land use DR-DID table exported: ", file.path(tab_dir, "did_drdid_landuse.tex"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 6. Summary diagnostics -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== Summary diagnostics ===")

message("\n--- Pre-treatment ATT(g,t) check ---")
for (y in fire_outcomes) {
  if (is.null(fire_agg_dyn[[y]])) next
  es <- extract_drdid_es(fire_agg_dyn[[y]])
  pre <- es[event_time < 0]
  if (nrow(pre) > 0) {
    n_sig <- sum(abs(pre$estimate) > 1.96 * pre$se, na.rm = TRUE)
    message(sprintf("  %s: %d/%d pre-treatment coefficients significant at 5%%",
                    y, n_sig, nrow(pre)))
    if (n_sig > 0) {
      sig_rows <- pre[abs(estimate) > 1.96 * se]
      for (j in seq_len(nrow(sig_rows))) {
        message(sprintf("    t=%d: ATT=%.4f (SE=%.4f)",
                        sig_rows$event_time[j], sig_rows$estimate[j], sig_rows$se[j]))
      }
    }
  }
}

for (y in lu_outcomes) {
  if (is.null(lu_agg_dyn[[y]])) next
  es <- extract_drdid_es(lu_agg_dyn[[y]])
  pre <- es[event_time < 0]
  if (nrow(pre) > 0) {
    n_sig <- sum(abs(pre$estimate) > 1.96 * pre$se, na.rm = TRUE)
    message(sprintf("  %s: %d/%d pre-treatment coefficients significant at 5%%",
                    y, n_sig, nrow(pre)))
  }
}

message("\n=== Script 3 complete ===")
