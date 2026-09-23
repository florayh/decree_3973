#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Project: Bolivia Fires
# Purpose: Heterogeneous DR-DID by land tenure type. Recomputes plurality
#          dominant tenure (40% share threshold), applies symmetric PS trimming
#          [0.01, 0.99] (Crump et al. 2009), and re-runs DR-DID (C&S 2021)
#          on each of 6 tenure subsets. Produces combined comparison tables.
#
# Design decisions:
#   - Elevation <= 1800m for controls (same as 3_drdid_main.R)
#   - Symmetric PS trimming [0.01, 0.99] (Crump et al. 2009)
#   - Dominant tenure requires >= 40% share (consistent with 0_psm_prep.R)
#   - 6 categories: private_large, private_small, campesino, tierra_fiscal,
#     indigenous, protected (NPA + forest reserve)
#
# Inputs:
#   - data/grid_month_panel_2km.rds (loaded by 1_shared_prep.R)
#   - data/grid_psm_weights_2km.rds (loaded by 1_shared_prep.R)
#
# Outputs:
#   Tables (output/tables/):
#     - did_drdid_fire_byland.tex
#     - did_drdid_landuse_byland.tex
#   Figures (output/figures/):
#     - did_drdid_byland_sample_map.pdf
#     - did_drdid_byland_lu_coefplot.pdf
#     - did_drdid_eventstudy_byland_{n_fires,primary_loss_km2,all_ag_sum,soy_area_km}.pdf
#
# RUNTIME NOTE: Expect several hours due to bootstrap SEs across 6 tenure subsets.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

source("code/1_shared_prep.R")
library(did)
library(sf)

set.seed(20260506)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 0. Settings -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Tenure categories (order = table column order)
tenure_cats <- c("private_large", "private_small", "campesino", "tierra_fiscal",
                  "indigenous", "protected")
tenure_labels <- c("Lg. private", "Sm. private", "Communal", "Untitled public",
                    "Indigenous", "Protected")
names(tenure_labels) <- tenure_cats

# Minimum grid threshold per tenure subset
MIN_GRIDS <- 50

# Outcomes for which to produce by-tenure event study figures
es_outcomes <- c("n_fires", "primary_loss_km2", "all_ag_sum", "soy_area_km")
es_labels   <- c("Total fires", "Primary forest loss", "Agriculture (total)", "Soy area")
names(es_labels) <- es_outcomes
es_fire <- intersect(es_outcomes, fire_outcomes)  # "n_fires"
es_lu   <- intersect(es_outcomes, lu_outcomes)     # the land use ones

# Symmetric PS trimming threshold [0.01, 0.99]
ps_threshold <- 0.01

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1. Elevation filter + pre-treatment climate means -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== Pre-treatment climate means ===")
# Elevation filter already applied in bundled panel (controls <= 1800m).

# Pre-treatment climate means (for PS trimming within tenure subsets)
grid_covs <- dt[event_time < 0, .(
  mean_temp_pre   = mean(temperature_2m, na.rm = TRUE),
  mean_precip_pre = mean(chirps_precip_mm, na.rm = TRUE)
), by = ID]
dt <- merge(dt, grid_covs, by = "ID", all.x = TRUE)
rm(grid_covs)

# Load shapefiles for by-tenure map below
grid_sf <- st_read(file.path("data", "shapefiles/grids/grid_bolivia_2km.shp"),
                   quiet = TRUE)
dept_sf <- st_read(file.path("data", "shapefiles/adm1",
                              "bol_admbnda_adm1_gov_2020514.shp"),
                   quiet = TRUE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2. Compute dominant tenure (40% share threshold) -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== Computing dominant tenure (40% threshold) ===")

dt <- merge(dt, compute_dominant_tenure(dt), by = "ID")

message("Panel after tenure assignment (6 categories): ", uniqueN(dt$ID), " grids, ",
        format(nrow(dt), big.mark = ","), " obs")

# Print final counts per tenure
message("\nFinal grid counts per tenure subset:")
for (t in tenure_cats) {
  n_tr <- uniqueN(dt[dominant_tenure == t & treated == 1, ID])
  n_ct <- uniqueN(dt[dominant_tenure == t & treated == 0, ID])
  message("  ", tenure_labels[t], ": ", n_tr, " treated, ", n_ct, " control")
}

#--- Map of treated vs control by tenure ---
message("\n=== Map of trimmed sample by land tenure ===")

grid_map_dt <- dt[, .(treated = treated[1], dominant_tenure = dominant_tenure[1]), by = ID]
grid_map_dt[, status := fifelse(treated == 1, "Treated", "Control")]
grid_map_dt[, tenure_label := tenure_labels[dominant_tenure]]
grid_map_dt[, tenure_label := factor(tenure_label, levels = tenure_labels)]

hex_tenure <- merge(grid_sf, grid_map_dt[, .(ID, status, tenure_label)],
                    by = "ID", all.x = FALSE)

p_tenure_map <- ggplot() +
  geom_sf(data = dept_sf, fill = "grey95", color = "grey70", linewidth = 0.2) +
  geom_sf(data = hex_tenure, aes(fill = status), color = NA, linewidth = 0) +
  geom_sf(data = dept_sf, fill = NA, color = "grey40", linewidth = 0.3) +
  facet_wrap(~ tenure_label, nrow = 2) +
  scale_fill_manual(
    values = c("Treated" = "#d62728", "Control" = "#1f77b4"),
    name = NULL
  ) +
  labs(title = "DR-DID Sample by Land Tenure (PS-trimmed [0.01, 0.99])") +
  theme_minimal(base_size = 10) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 10))

ggsave(file.path(fig_dir, "did_drdid_byland_sample_map.pdf"),
       p_tenure_map, width = 12, height = 12, device = cairo_pdf)
message("Map exported: ", file.path(fig_dir, "did_drdid_byland_sample_map.pdf"))

rm(grid_sf, dept_sf, hex_tenure, grid_map_dt, p_tenure_map)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3. Helpers -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4. Loop over tenure categories -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xformla <- ~ elevation + slope + pct_forest2014 + log_pop
stopifnot(identical(all.vars(xformla),
                    c("elevation", "slope", "pct_forest2014", "log_pop")))
message("DR-DID covariates: ", paste(all.vars(xformla), collapse = ", "))

# Storage for combined tables
all_fire_atts <- list()
all_lu_atts   <- list()

# Storage for dynamic event study aggregations (by-tenure, selected outcomes only)
all_es_dyn <- list()  # keyed by paste(tenure, outcome)

# Storage for analysis-sample counts (after all filters/balancing)
fire_n_treated <- setNames(rep(NA_integer_, length(tenure_cats)), tenure_cats)
fire_n_control <- setNames(rep(NA_integer_, length(tenure_cats)), tenure_cats)
lu_n_treated   <- setNames(rep(NA_integer_, length(tenure_cats)), tenure_cats)
lu_n_control   <- setNames(rep(NA_integer_, length(tenure_cats)), tenure_cats)

for (tenure in tenure_cats) {
  message("\n", strrep("=", 70))
  message("=== Tenure: ", tenure_labels[tenure], " (", tenure, ") ===")
  message(strrep("=", 70))

  dt_sub <- dt[dominant_tenure == tenure]
  n_treated <- uniqueN(dt_sub$ID[dt_sub$treated == 1])
  n_control <- uniqueN(dt_sub$ID[dt_sub$treated == 0])
  message("Grids: ", n_treated, " treated, ", n_control, " control")

  # Within-subsample PS trimming [0.01, 0.99]
  grid_sub <- dt_sub[, .(
    treated = treated[1], elevation = elevation[1], slope = slope[1],
    pct_forest2014 = pct_forest2014[1], log_pop = log1p(population[1]),
    mean_precip_pre = mean_precip_pre[1], mean_temp_pre = mean_temp_pre[1]
  ), by = ID]
  ps_sub <- glm(treated ~ elevation + slope + pct_forest2014 + log_pop +
                  mean_precip_pre + mean_temp_pre,
                data = grid_sub, family = binomial(link = "logit"))
  grid_sub[, ps := predict(ps_sub, newdata = grid_sub, type = "response")]
  drop_sub <- grid_sub[
    (treated == 0 & (is.na(ps) | ps < ps_threshold)) |
    (treated == 1 & (is.na(ps) | ps > (1 - ps_threshold))), ID]
  if (length(drop_sub) > 0) {
    dt_sub <- dt_sub[!(ID %in% drop_sub)]
    n_treated <- uniqueN(dt_sub$ID[dt_sub$treated == 1])
    n_control <- uniqueN(dt_sub$ID[dt_sub$treated == 0])
    message("  Within-subsample PS trimming: dropped ", length(drop_sub),
            " grids -> ", n_treated, " treated, ", n_control, " control")
  }
  rm(grid_sub, ps_sub, drop_sub)

  # Check minimum sample
  if (n_treated < MIN_GRIDS || n_control < MIN_GRIDS) {
    message("  WARNING: Skipping -- insufficient grids (need >= ", MIN_GRIDS,
            " treated AND control)")
    all_fire_atts[[tenure]] <- data.table(
      tenure = tenure, outcome = fire_outcomes, att = NA_real_, se = NA_real_
    )
    all_lu_atts[[tenure]] <- data.table(
      tenure = tenure, outcome = lu_outcomes, att = NA_real_, se = NA_real_
    )
    next
  }

  #--- 4a. Fire season panel (Jul-Oct, 2015-2022) ---
  message("\n  Building fire season panel...")
  dt_fire_sub <- dt_sub[month %in% 7:10 & year %in% 2015:2022]

  # Consecutive fire_period indexing
  fire_month_indices <- sort(unique(dt_fire_sub$month_index))
  dt_fire_sub[, fire_period := match(month_index, fire_month_indices)]

  # Treatment timing
  treatment_fire_period <- match(54L, fire_month_indices)
  dt_fire_sub[, treatment_period := fifelse(treated == 1,
                                             as.numeric(treatment_fire_period), 0)]

  dt_fire_sub[, log_pop := log1p(population)]

  # Balance check
  n_periods_fire <- uniqueN(dt_fire_sub$fire_period)
  obs_per_grid <- dt_fire_sub[, .N, by = ID]
  if (!all(obs_per_grid$N == n_periods_fire)) {
    complete_grids <- obs_per_grid[N == n_periods_fire, ID]
    dt_fire_sub <- dt_fire_sub[ID %in% complete_grids]
    message("  Balanced fire panel: ", uniqueN(dt_fire_sub$ID), " grids x ",
            n_periods_fire, " periods")
  } else {
    message("  Fire panel: ", uniqueN(dt_fire_sub$ID), " grids x ",
            n_periods_fire, " periods (balanced)")
  }

  df_fire <- as.data.frame(dt_fire_sub)

  # Record analysis-sample counts (after balancing)
  fire_n_treated[tenure] <- uniqueN(dt_fire_sub$ID[dt_fire_sub$treated == 1])
  fire_n_control[tenure] <- uniqueN(dt_fire_sub$ID[dt_fire_sub$treated == 0])

  #--- 4b. Annual land use panel (2015-2022) ---
  message("  Building annual land use panel...")
  dt_year_sub <- build_annual_panel(dt_sub[year %in% 2015:2022])

  # Treatment timing
  dt_year_sub[, treatment_year := fifelse(treated == 1, 2019, 0)]

  dt_year_sub[, log_pop := log1p(population)]

  message("  Land use panel: ", uniqueN(dt_year_sub$ID), " grids x ",
          uniqueN(dt_year_sub$year), " years")

  df_year <- as.data.frame(dt_year_sub)

  # Record analysis-sample counts
  lu_n_treated[tenure] <- uniqueN(dt_year_sub$ID[dt_year_sub$treated == 1])
  lu_n_control[tenure] <- uniqueN(dt_year_sub$ID[dt_year_sub$treated == 0])

  #--- 4c. DR-DID: Fire outcomes ---
  message("\n  DR-DID: Fire outcomes")
  tenure_fire_atts <- data.table(tenure = tenure, outcome = fire_outcomes,
                                  att = NA_real_, se = NA_real_)

  for (i in seq_along(fire_outcomes)) {
    y <- fire_outcomes[i]
    message("    ", y, "...")

    attgt_obj <- tryCatch(
      att_gt(
        yname         = y,
        tname         = "fire_period",
        idname        = "ID",
        gname         = "treatment_period",
        xformla       = xformla,
        control_group = "nevertreated",
        anticipation  = 0,
        est_method    = "dr",
        clustervars   = "cluster_20km",
        data          = df_fire,
        print_details = FALSE
      ),
      error = function(e) { message("      ERROR: ", e$message); NULL }
    )

    if (!is.null(attgt_obj)) {
      agg_grp <- tryCatch(aggte(attgt_obj, type = "group"),
                           error = function(e) NULL)

      if (!is.null(agg_grp)) {
        tenure_fire_atts[outcome == y, `:=`(att = agg_grp$overall.att,
                                             se = agg_grp$overall.se)]
        message("      ATT: ", round(agg_grp$overall.att, 4),
                " (SE: ", round(agg_grp$overall.se, 4), ")")
      }

      # Store dynamic aggregation for selected event study outcomes
      if (y %in% es_fire) {
        agg_dyn <- tryCatch(aggte(attgt_obj, type = "dynamic"),
                             error = function(e) NULL)
        if (!is.null(agg_dyn)) all_es_dyn[[paste(tenure, y, sep = "|")]] <- agg_dyn
      }
    }
  }
  all_fire_atts[[tenure]] <- tenure_fire_atts

  #--- 4d. DR-DID: Land use outcomes ---
  message("\n  DR-DID: Land use outcomes")
  tenure_lu_atts <- data.table(tenure = tenure, outcome = lu_outcomes,
                                att = NA_real_, se = NA_real_)

  for (i in seq_along(lu_outcomes)) {
    y <- lu_outcomes[i]
    message("    ", y, "...")

    attgt_obj <- tryCatch(
      att_gt(
        yname         = y,
        tname         = "year",
        idname        = "ID",
        gname         = "treatment_year",
        xformla       = xformla,
        control_group = "nevertreated",
        anticipation  = 0,
        est_method    = "dr",
        clustervars   = "cluster_20km",
        data          = df_year,
        print_details = FALSE
      ),
      error = function(e) { message("      ERROR: ", e$message); NULL }
    )

    if (!is.null(attgt_obj)) {
      agg_grp <- tryCatch(aggte(attgt_obj, type = "group"),
                           error = function(e) NULL)

      if (!is.null(agg_grp)) {
        tenure_lu_atts[outcome == y, `:=`(att = agg_grp$overall.att,
                                           se = agg_grp$overall.se)]
        message("      ATT: ", round(agg_grp$overall.att, 4),
                " (SE: ", round(agg_grp$overall.se, 4), ")")
      }

      # Store dynamic aggregation for selected event study outcomes
      if (y %in% es_lu) {
        agg_dyn <- tryCatch(aggte(attgt_obj, type = "dynamic"),
                             error = function(e) NULL)
        if (!is.null(agg_dyn)) all_es_dyn[[paste(tenure, y, sep = "|")]] <- agg_dyn
      }
    }
  }
  all_lu_atts[[tenure]] <- tenure_lu_atts
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5. Combined comparison tables -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== Building combined comparison tables ===")

fire_atts <- rbindlist(all_fire_atts)
lu_atts   <- rbindlist(all_lu_atts)

build_byland_tex <- function(atts_dt, outcome_vec, label_vec, caption, label_tag,
                             filename, panel_type, treated_means = NULL,
                             n_treated_vec = NULL, n_control_vec = NULL) {
  n_tenures <- length(tenure_cats)

  # Header row
  col_header <- paste(paste0("& ", tenure_labels[tenure_cats]), collapse = " ")
  tex <- c(
    "\\begin{table*}[t!]",
    "\\centering",
    "\\small",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{tab:", label_tag, "}"),
    paste0("\\begin{tabular}{l", paste(rep("c", n_tenures), collapse = ""), "}"),
    "\\toprule",
    paste0(" ", col_header, " \\\\"),
    "\\midrule"
  )

  for (i in seq_along(outcome_vec)) {
    y <- outcome_vec[i]
    cells <- character(n_tenures)
    se_cells <- character(n_tenures)
    cm_cells <- character(n_tenures)

    for (j in seq_along(tenure_cats)) {
      t <- tenure_cats[j]
      row <- atts_dt[tenure == t & outcome == y]
      if (nrow(row) == 0 || is.na(row$att)) {
        cells[j] <- "---"
        se_cells[j] <- ""
      } else {
        p_val <- 2 * pnorm(-abs(row$att / row$se))
        cells[j] <- paste0(sprintf("%.3f", row$att), add_stars(p_val))
        se_cells[j] <- paste0("(", sprintf("%.3f", row$se), ")")
      }
      # Treated pre-treatment mean
      if (!is.null(treated_means) && y %in% rownames(treated_means) &&
          t %in% colnames(treated_means)) {
        cm_cells[j] <- sprintf("%.3f", treated_means[y, t])
      } else {
        cm_cells[j] <- ""
      }
    }

    tex <- c(tex,
      paste0(label_vec[i], " & ", paste(cells, collapse = " & "), " \\\\"),
      paste0(" & ", paste(se_cells, collapse = " & "), " \\\\"),
      paste0("Treated mean (pre) & ", paste(cm_cells, collapse = " & "), " \\\\")
    )
    if (i < length(outcome_vec)) tex <- c(tex, "[0.5em]")
  }

  notes_text <- paste0(
    "\\item \\textit{Notes:} DR-DID ATT (Sant'Anna \\& Zhao 2020; \\texttt{did} package, Callaway \\& Sant'Anna 2021) ",
    "estimated separately for each tenure category. Dominant tenure = ",
    "plurality (highest share $\\geq$ 40\\%). ",
    "The trimming propensity score includes elevation, slope, 2014 forest share, ",
    "log population, and pre-treatment mean precipitation and temperature. ",
    "DR-DID nuisance models include elevation, slope, 2014 forest share, and log population. ",
    "Symmetric PS trimming (Crump et al.\\ 2009): within each tenure subset, ",
    "controls with PS $<$ 0.01 and treated with PS $>$ 0.99 dropped; ",
    "elevation $\\leq$ 1800m for controls. ",
    panel_type, ". ",
    "SE clustered at 20km super-grid level (bootstrap). ",
    "$^{+}p<0.1$, $^{*}p<0.05$, $^{**}p<0.01$, $^{***}p<0.001$."
  )

  # Sample size rows
  if (!is.null(n_treated_vec) && !is.null(n_control_vec)) {
    tr_cells <- sapply(tenure_cats, function(t) {
      if (is.na(n_treated_vec[t])) "---"
      else format(n_treated_vec[t], big.mark = ",")
    })
    ct_cells <- sapply(tenure_cats, function(t) {
      if (is.na(n_control_vec[t])) "---"
      else format(n_control_vec[t], big.mark = ",")
    })
    tex <- c(tex,
      "\\midrule",
      paste0("Treated & ", paste(tr_cells, collapse = " & "), " \\\\"),
      paste0("Control & ", paste(ct_cells, collapse = " & "), " \\\\")
    )
  }

  tex <- c(tex,
    "\\bottomrule",
    "\\end{tabular}",
    "\\begin{tablenotes}[flushleft]",
    "\\small",
    notes_text,
    "\\end{tablenotes}",
    "\\end{table*}"
  )

  stopifnot(
    identical(sum(tex == "\\begin{table*}[t!]"), 1L),
    identical(sum(tex == "\\end{table*}"), 1L)
  )
  writeLines(tex, file.path(tab_dir, filename))
  message("Exported: ", file.path(tab_dir, filename))
}

# Compute treated pre-treatment means per tenure x outcome
fire_treat_means <- matrix(NA_real_, nrow = length(fire_outcomes),
                           ncol = length(tenure_cats),
                           dimnames = list(fire_outcomes, tenure_cats))
lu_treat_means <- matrix(NA_real_, nrow = length(lu_outcomes),
                         ncol = length(tenure_cats),
                         dimnames = list(lu_outcomes, tenure_cats))

for (tenure in tenure_cats) {
  dt_sub <- dt[dominant_tenure == tenure]
  # Fire: treated grids, pre-treatment, fire season
  dt_fire_pre <- dt_sub[treated == 1 & month %in% 7:10 & year %in% 2015:2018]
  for (y in fire_outcomes) {
    if (nrow(dt_fire_pre) > 0) {
      fire_treat_means[y, tenure] <- mean(dt_fire_pre[[y]], na.rm = TRUE)
    }
  }
  # Land use: treated grids, pre-treatment (annual, year < 2019)
  dt_lu_pre <- build_annual_panel(dt_sub[treated == 1 & year %in% 2015:2018])
  for (y in lu_outcomes) {
    if (nrow(dt_lu_pre) > 0) {
      lu_treat_means[y, tenure] <- mean(dt_lu_pre[[y]], na.rm = TRUE)
    }
  }
}

build_byland_tex(
  fire_atts, fire_outcomes, fire_labels,
  caption = "DR-DID by Land Tenure: Fire Outcomes",
  label_tag = "did_drdid_fire_byland",
  filename = "did_drdid_fire_byland.tex",
  panel_type = "Fire season (Jul--Oct), 2015--2022",
  treated_means = fire_treat_means,
  n_treated_vec = fire_n_treated,
  n_control_vec = fire_n_control
)

build_byland_tex(
  lu_atts, lu_outcomes, lu_labels,
  caption = "DR-DID by Land Tenure: Land Use Outcomes",
  label_tag = "did_drdid_landuse_byland",
  filename = "did_drdid_landuse_byland.tex",
  panel_type = "Annual panel, 2015--2022",
  treated_means = lu_treat_means,
  n_treated_vec = lu_n_treated,
  n_control_vec = lu_n_control
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5b. Land use coefficient plot by tenure -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== Land use coefficient plot by tenure ===")

lu_plot <- lu_atts[!is.na(att)]
lu_plot[, `:=`(
  ci_low  = att - 1.96 * se,
  ci_high = att + 1.96 * se
)]

lu_label_map <- setNames(lu_labels, lu_outcomes)
lu_plot[, outcome_label := factor(lu_label_map[outcome], levels = lu_labels)]
lu_plot[, tenure_label := factor(tenure_labels[tenure], levels = rev(tenure_labels))]

p_lu_coef <- ggplot(lu_plot, aes(x = att, y = tenure_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(xmin = ci_low, xmax = ci_high), size = 0.4) +
  facet_wrap(~ outcome_label, scales = "free_x", ncol = 3) +
  labs(x = "ATT estimate (95% CI)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(fig_dir, "did_drdid_byland_lu_coefplot.pdf"),
       p_lu_coef, width = 10, height = 6, device = cairo_pdf)
message("Exported: ", file.path(fig_dir, "did_drdid_byland_lu_coefplot.pdf"))

rm(lu_plot, p_lu_coef)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5c. By-tenure event study figures (selected outcomes) -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== By-tenure event study figures ===")

for (y in es_outcomes) {
  # Collect event study data across tenure types
  es_all <- rbindlist(lapply(tenure_cats, function(t) {
    key <- paste(t, y, sep = "|")
    es <- extract_drdid_es(all_es_dyn[[key]])
    if (is.null(es)) return(NULL)
    es[, tenure := t]
    es
  }))

  if (nrow(es_all) == 0) {
    message("  Skipping ", y, " -- no event study data for any tenure")
    next
  }

  es_all[, tenure_label := factor(tenure_labels[tenure], levels = tenure_labels)]

  # X-axis label depends on fire vs land use outcome
  if (y %in% es_fire) {
    x_lab <- "Fire-season months relative to Jul 2019"
  } else {
    x_lab <- "Years relative to 2019"
  }

  p <- ggplot(es_all, aes(x = event_time, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "solid", color = "grey50") +
    geom_vline(xintercept = -0.5, linetype = "dashed", color = "red", alpha = 0.7) +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15, fill = "#1f77b4") +
    geom_point(size = 1, color = "#1f77b4") +
    geom_line(linewidth = 0.3, color = "#1f77b4") +
    facet_wrap(~ tenure_label, nrow = 3, scales = "free_y") +
    labs(x = x_lab, y = paste0("ATT on ", tolower(es_labels[y]))) +
    theme_minimal(base_size = 10) +
    theme(strip.text = element_text(face = "bold", size = 9))

  out_file <- file.path(fig_dir, paste0("did_drdid_eventstudy_byland_", y, ".pdf"))
  ggsave(out_file, p, width = 8, height = 9, device = cairo_pdf)
  message("  Exported: ", basename(out_file))
}

rm(all_es_dyn)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 6. Pre-treatment diagnostics -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== Pre-treatment ATT(g,t) diagnostics ===")
message("(Checking for systematic pre-trend violations per tenure)")

message("\n--- Quick pre-treatment check (n_fires by tenure) ---")

for (tenure in tenure_cats) {
  key <- paste(tenure, "n_fires", sep = "|")
  agg_dyn <- all_es_dyn[[key]]

  if (is.null(agg_dyn)) {
    message("  ", tenure_labels[tenure], ": no stored event study (skipped or failed)")
    next
  }

  es <- extract_drdid_es(agg_dyn)
  pre <- es[event_time < 0]
  n_sig <- sum(abs(pre$estimate) > 1.96 * pre$se, na.rm = TRUE)
  message(sprintf("  %s: %d/%d pre-treatment coefficients significant at 5%%",
                  tenure_labels[tenure], n_sig, nrow(pre)))
  if (n_sig > 0) {
    sig_rows <- pre[abs(estimate) > 1.96 * se]
    for (j in seq_len(nrow(sig_rows))) {
      message(sprintf("    t=%d: ATT=%.4f (SE=%.4f)",
                      sig_rows$event_time[j], sig_rows$estimate[j], sig_rows$se[j]))
    }
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 7. Summary -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== Summary ===")
message("Tenure categories processed:")
for (t in tenure_cats) {
  fire_row <- fire_atts[tenure == t]
  n_estimated <- sum(!is.na(fire_row$att))
  lu_row <- lu_atts[tenure == t]
  n_estimated_lu <- sum(!is.na(lu_row$att))
  message("  ", tenure_labels[t], ": ", n_estimated, "/",
          length(fire_outcomes), " fire + ",
          n_estimated_lu, "/", length(lu_outcomes), " land use outcomes estimated")
}

message("\n=== Script 4 (drdid_byland) complete ===")
