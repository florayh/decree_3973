#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Project: Bolivia Fires
# Purpose: DR-DID robustness tables restricted to the decree period only
#          (2015-2020). The main DR-DID results (scripts 3 and 4) use panels
#          through 2022, so the aggregated ATT pools over both the decree period
#          and the post-rescission period. This script isolates just the decree
#          period (Jul 2019 - Sep 2020; SD 3973 revoked 16 Sep 2020).
#          Full 2020 fire season (Jul-Oct) included because the decree was active
#          for ~85% of Jul-Sep 2020, and October fires largely continue clearing
#          started earlier. For annual land use, MapBiomas captures the full year
#          and decisions were made while the decree was active.
#
# Inputs:
#   - data/grid_month_panel_2km.rds (via 1_shared_prep.R)
#
# Outputs:
#   Tables (output/tables/):
#     - did_drdid_fire_decree.tex
#     - did_drdid_landuse_decree.tex
#     - did_drdid_fire_byland_decree.tex
#     - did_drdid_landuse_byland_decree.tex
#
# RUNTIME NOTE: Expect several hours due to bootstrap SEs.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

source("code/1_shared_prep.R")
library(did)

set.seed(20260625)

message("Full panel: ", format(nrow(dt), big.mark = ","), " rows")

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

# Symmetric PS trimming threshold (used within tenure subsets below)
ps_threshold <- 0.01

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2. Fire panel: fire season (Jul-Oct), 2015-2020 -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== Building fire season panel (2015-2020) ===")

dt_fire <- dt[month %in% 7:10 & year %in% 2015:2020]

# Consecutive fire_period indexing for att_gt()
fire_month_indices <- sort(unique(dt_fire$month_index))
dt_fire[, fire_period := match(month_index, fire_month_indices)]

# Treatment timing: Jul 2019 = month_index 54
treatment_fire_period <- match(54L, fire_month_indices)
message("Treatment fire_period: ", treatment_fire_period,
        " (month_index 54 = Jul 2019)")

dt_fire[, treatment_period := fifelse(treated == 1,
                                       as.numeric(treatment_fire_period), 0)]

dt_fire[, log_pop := log1p(population)]
dt_fire[, month_f := as.factor(month)]

# Balance check
n_grids_fire <- uniqueN(dt_fire$ID)
n_periods_fire <- uniqueN(dt_fire$fire_period)
message("Fire panel: ", n_grids_fire, " grids x ", n_periods_fire,
        " fire-season months = ", format(nrow(dt_fire), big.mark = ","), " obs")

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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3. Land use panel: annual 2015-2020 -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== Building annual land use panel (2015-2020) ===")

dt_year <- build_annual_panel(dt[year %in% 2015:2020])

# Treatment timing: gname = 2019 for treated, 0 for never-treated
dt_year[, treatment_year := fifelse(treated == 1, 2019, 0)]

dt_year[, log_pop := log1p(population)]

n_grids_yr <- uniqueN(dt_year$ID)
n_years <- uniqueN(dt_year$year)
message("Land use panel: ", n_grids_yr, " grids x ", n_years,
        " years = ", format(nrow(dt_year), big.mark = ","), " obs")
message("Treated: ", uniqueN(dt_year$ID[dt_year$treated == 1]),
        " | Control: ", uniqueN(dt_year$ID[dt_year$treated == 0]))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4. DR-DID: Fire outcomes -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== DR-DID: Fire outcomes (decree period only) ===")

xformla <- ~ elevation + slope + pct_forest2014 + log_pop
stopifnot(identical(all.vars(xformla),
                    c("elevation", "slope", "pct_forest2014", "log_pop")))
message("DR-DID covariates: ", paste(all.vars(xformla), collapse = ", "))

df_fire <- as.data.frame(dt_fire)

fire_agg_grp <- list()

for (i in seq_along(fire_outcomes)) {
  y <- fire_outcomes[i]
  message("  Estimating att_gt for: ", y)

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
      clustervars   = "ADM3_PCODE",
      data          = df_fire,
      print_details = FALSE
    ),
    error = function(e) {
      message("    ERROR: ", e$message)
      NULL
    }
  )

  if (!is.null(attgt_obj)) {
    message("    att_gt succeeded: ", length(attgt_obj$att), " ATT(g,t) estimates")

    fire_agg_grp[[y]] <- tryCatch(
      aggte(attgt_obj, type = "group"),
      error = function(e) { message("    Group agg error: ", e$message); NULL }
    )

    if (!is.null(fire_agg_grp[[y]])) {
      message("    Overall ATT: ", round(fire_agg_grp[[y]]$overall.att, 4),
              " (SE: ", round(fire_agg_grp[[y]]$overall.se, 4), ")")
    }
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5. DR-DID: Land use outcomes -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== DR-DID: Land use outcomes (decree period only) ===")

df_year <- as.data.frame(dt_year)

lu_agg_grp <- list()

for (i in seq_along(lu_outcomes)) {
  y <- lu_outcomes[i]
  message("  Estimating att_gt for: ", y)

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
      clustervars   = "ADM3_PCODE",
      data          = df_year,
      print_details = FALSE
    ),
    error = function(e) {
      message("    ERROR: ", e$message)
      NULL
    }
  )

  if (!is.null(attgt_obj)) {
    message("    att_gt succeeded: ", length(attgt_obj$att), " ATT(g,t) estimates")

    lu_agg_grp[[y]] <- tryCatch(
      aggte(attgt_obj, type = "group"),
      error = function(e) { message("    Group agg error: ", e$message); NULL }
    )

    if (!is.null(lu_agg_grp[[y]])) {
      message("    Overall ATT: ", round(lu_agg_grp[[y]]$overall.att, 4),
              " (SE: ", round(lu_agg_grp[[y]]$overall.se, 4), ")")
    }
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 6. Main DR-DID summary tables -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== Building DR-DID summary tables (decree period only) ===")

decree_note_fire <- paste0(
  "\\item \\textit{Notes:} DR-DID = Sant'Anna \\& Zhao (2020) doubly robust ",
  "estimator, implemented via the \\texttt{did} package (Callaway \\& Sant'Anna, ",
  "2021); ATT aggregated across post-treatment periods. ",
  "SD~3973 decree period only: panel truncated at 2020. ",
  "Fire season months (Jul--Oct), 2015--2020. ",
  "Full 2020 fire season included (decree active for $\\sim$85\\% of Jul--Sep 2020; ",
  "revoked 16 Sep 2020 by SD~4389). ",
  "DR-DID nuisance models include elevation, slope, 2014 forest share, and log population. ",
  "Sample: control grids with elevation $>$ 1800m dropped. ",
  "SE clustered at municipality level (bootstrap). ",
  "$^{+}p<0.1$, $^{*}p<0.05$, $^{**}p<0.01$, $^{***}p<0.001$."
)

decree_note_lu <- paste0(
  "\\item \\textit{Notes:} DR-DID = Sant'Anna \\& Zhao (2020) doubly robust ",
  "estimator, implemented via the \\texttt{did} package (Callaway \\& Sant'Anna, ",
  "2021); ATT aggregated across post-treatment periods. ",
  "SD~3973 decree period only: panel truncated at 2020. ",
  "Annual panel 2015--2020. ",
  "Full 2020 included (decree active for $\\sim$85\\% of Jul--Sep 2020; ",
  "revoked 16 Sep 2020 by SD~4389; MapBiomas captures full year). ",
  "DR-DID nuisance models include elevation, slope, 2014 forest share, and log population. ",
  "Sample: control grids with elevation $>$ 1800m dropped. ",
  "SE clustered at municipality level (bootstrap). ",
  "$^{+}p<0.1$, $^{*}p<0.05$, $^{**}p<0.01$, $^{***}p<0.001$."
)

# --- Fire table ---
fire_tex <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\small",
  "\\caption{Doubly Robust DID: Fire Outcomes (Decree Period Only)}",
  "\\label{tab:did_drdid_fire_decree}",
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
  decree_note_fire,
  "\\end{tablenotes}",
  "\\end{table}"
)

writeLines(fire_tex, file.path(tab_dir, "did_drdid_fire_decree.tex"))
message("Fire DR-DID table exported: ", file.path(tab_dir, "did_drdid_fire_decree.tex"))

# --- Land use table ---
lu_tex <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\small",
  "\\caption{Doubly Robust DID: Land Use Outcomes (Decree Period Only)}",
  "\\label{tab:did_drdid_landuse_decree}",
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
  decree_note_lu,
  "\\end{tablenotes}",
  "\\end{table}"
)

writeLines(lu_tex, file.path(tab_dir, "did_drdid_landuse_decree.tex"))
message("Land use DR-DID table exported: ", file.path(tab_dir, "did_drdid_landuse_decree.tex"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 7. Dominant tenure assignment (40% threshold) -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== Computing dominant tenure (40% threshold) ===")

tenure_cats <- c("private_large", "private_small", "campesino", "tierra_fiscal",
                  "indigenous", "protected")
tenure_labels_ten <- c("Lg. private", "Sm. private", "Communal", "Untitled public",
                        "Indigenous", "Protected")
names(tenure_labels_ten) <- tenure_cats

MIN_GRIDS <- 50

# Compute dominant tenure on the elevation-filtered panel.
dt_ten <- merge(dt, compute_dominant_tenure(dt), by = "ID")

message("Panel after tenure assignment: ", uniqueN(dt_ten$ID), " grids, ",
        format(nrow(dt_ten), big.mark = ","), " obs")

for (t in tenure_cats) {
  n_tr <- uniqueN(dt_ten[dominant_tenure == t & treated == 1, ID])
  n_ct <- uniqueN(dt_ten[dominant_tenure == t & treated == 0, ID])
  message("  ", tenure_labels_ten[t], ": ", n_tr, " treated, ", n_ct, " control")
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 8. DR-DID by tenure: loop over categories -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== DR-DID by tenure (decree period only) ===")

# Storage
all_fire_atts <- list()
all_lu_atts   <- list()

fire_n_treated <- setNames(rep(NA_integer_, length(tenure_cats)), tenure_cats)
fire_n_control <- setNames(rep(NA_integer_, length(tenure_cats)), tenure_cats)
lu_n_treated   <- setNames(rep(NA_integer_, length(tenure_cats)), tenure_cats)
lu_n_control   <- setNames(rep(NA_integer_, length(tenure_cats)), tenure_cats)

for (tenure in tenure_cats) {
  message("\n", strrep("=", 70))
  message("=== Tenure: ", tenure_labels_ten[tenure], " (", tenure, ") ===")
  message(strrep("=", 70))

  dt_sub <- dt_ten[dominant_tenure == tenure]
  n_treated <- uniqueN(dt_sub$ID[dt_sub$treated == 1])
  n_control <- uniqueN(dt_sub$ID[dt_sub$treated == 0])
  message("Grids: ", n_treated, " treated, ", n_control, " control")

  # Within-subsample PS trimming [0.01, 0.99] (same as script 4)
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

  #--- Fire season panel (Jul-Oct, 2015-2020) ---
  message("\n  Building fire season panel (2015-2020)...")
  dt_fire_sub <- dt_sub[month %in% 7:10 & year %in% 2015:2020]

  fire_month_indices_sub <- sort(unique(dt_fire_sub$month_index))
  dt_fire_sub[, fire_period := match(month_index, fire_month_indices_sub)]

  treatment_fire_period_sub <- match(54L, fire_month_indices_sub)
  dt_fire_sub[, treatment_period := fifelse(treated == 1,
                                             as.numeric(treatment_fire_period_sub), 0)]

  dt_fire_sub[, log_pop := log1p(population)]

  # Balance check
  n_periods_sub <- uniqueN(dt_fire_sub$fire_period)
  obs_per_grid_sub <- dt_fire_sub[, .N, by = ID]
  if (!all(obs_per_grid_sub$N == n_periods_sub)) {
    complete_grids <- obs_per_grid_sub[N == n_periods_sub, ID]
    dt_fire_sub <- dt_fire_sub[ID %in% complete_grids]
    message("  Balanced fire panel: ", uniqueN(dt_fire_sub$ID), " grids x ",
            n_periods_sub, " periods")
  } else {
    message("  Fire panel: ", uniqueN(dt_fire_sub$ID), " grids x ",
            n_periods_sub, " periods (balanced)")
  }

  df_fire_sub <- as.data.frame(dt_fire_sub)

  fire_n_treated[tenure] <- uniqueN(dt_fire_sub$ID[dt_fire_sub$treated == 1])
  fire_n_control[tenure] <- uniqueN(dt_fire_sub$ID[dt_fire_sub$treated == 0])

  #--- Annual land use panel (2015-2020) ---
  message("  Building annual land use panel (2015-2020)...")
  dt_year_sub <- build_annual_panel(dt_sub[year %in% 2015:2020])

  dt_year_sub[, treatment_year := fifelse(treated == 1, 2019, 0)]

  dt_year_sub[, log_pop := log1p(population)]

  message("  Land use panel: ", uniqueN(dt_year_sub$ID), " grids x ",
          uniqueN(dt_year_sub$year), " years")

  df_year_sub <- as.data.frame(dt_year_sub)

  lu_n_treated[tenure] <- uniqueN(dt_year_sub$ID[dt_year_sub$treated == 1])
  lu_n_control[tenure] <- uniqueN(dt_year_sub$ID[dt_year_sub$treated == 0])

  #--- DR-DID: Fire outcomes ---
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
        data          = df_fire_sub,
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
    }
  }
  all_fire_atts[[tenure]] <- tenure_fire_atts

  #--- DR-DID: Land use outcomes ---
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
        data          = df_year_sub,
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
    }
  }
  all_lu_atts[[tenure]] <- tenure_lu_atts
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 9. By-tenure comparison tables -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== Building by-tenure comparison tables (decree period only) ===")

fire_atts <- rbindlist(all_fire_atts)
lu_atts   <- rbindlist(all_lu_atts)

build_byland_tex <- function(atts_dt, outcome_vec, label_vec, caption, label_tag,
                             filename, panel_type, treated_means = NULL,
                             n_treated_vec = NULL, n_control_vec = NULL) {
  n_tenures <- length(tenure_cats)

  col_header <- paste(paste0("& ", tenure_labels_ten[tenure_cats]), collapse = " ")
  tex <- c(
    "\\begin{table}[htbp]",
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
    "\\item \\textit{Notes:} DR-DID ATT (Sant'Anna \\& Zhao 2020; \\texttt{did} package, ",
    "Callaway \\& Sant'Anna 2021) estimated separately for each tenure category. ",
    "Dominant tenure = plurality (highest share $\\geq$ 40\\%). ",
    "SD~3973 decree period only: panel truncated at 2020 (revoked 16 Sep 2020 by SD~4389). ",
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
    "\\end{table}"
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
  dt_sub <- dt_ten[dominant_tenure == tenure]
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
  caption = "DR-DID by Land Tenure: Fire Outcomes (Decree Period Only)",
  label_tag = "did_drdid_fire_byland_decree",
  filename = "did_drdid_fire_byland_decree.tex",
  panel_type = "Fire season (Jul--Oct), 2015--2020",
  treated_means = fire_treat_means,
  n_treated_vec = fire_n_treated,
  n_control_vec = fire_n_control
)

build_byland_tex(
  lu_atts, lu_outcomes, lu_labels,
  caption = "DR-DID by Land Tenure: Land Use Outcomes (Decree Period Only)",
  label_tag = "did_drdid_landuse_byland_decree",
  filename = "did_drdid_landuse_byland_decree.tex",
  panel_type = "Annual panel, 2015--2020",
  treated_means = lu_treat_means,
  n_treated_vec = lu_n_treated,
  n_control_vec = lu_n_control
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 10. Summary diagnostics -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== Summary ===")

message("\n--- Main DR-DID (decree period only) ---")
for (i in seq_along(fire_outcomes)) {
  y <- fire_outcomes[i]
  if (!is.null(fire_agg_grp[[y]])) {
    message(sprintf("  %s: ATT=%.4f (SE=%.4f)",
                    fire_labels[i], fire_agg_grp[[y]]$overall.att,
                    fire_agg_grp[[y]]$overall.se))
  } else {
    message("  ", fire_labels[i], ": FAILED")
  }
}

for (i in seq_along(lu_outcomes)) {
  y <- lu_outcomes[i]
  if (!is.null(lu_agg_grp[[y]])) {
    message(sprintf("  %s: ATT=%.4f (SE=%.4f)",
                    lu_labels[i], lu_agg_grp[[y]]$overall.att,
                    lu_agg_grp[[y]]$overall.se))
  } else {
    message("  ", lu_labels[i], ": FAILED")
  }
}

message("\n--- By-tenure DR-DID (decree period only) ---")
for (t in tenure_cats) {
  fire_row <- fire_atts[tenure == t]
  n_est_fire <- sum(!is.na(fire_row$att))
  lu_row <- lu_atts[tenure == t]
  n_est_lu <- sum(!is.na(lu_row$att))
  message("  ", tenure_labels_ten[t], ": ", n_est_fire, "/",
          length(fire_outcomes), " fire + ",
          n_est_lu, "/", length(lu_outcomes), " land use outcomes estimated")
}

message("\n=== Script 6 (drdid_decree_only) complete ===")
