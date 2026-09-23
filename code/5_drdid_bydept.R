#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Project: Bolivia Fires
# Purpose: Heterogeneous DR-DID by department x land tenure. Applies elevation
#          filter (controls <= 1800m), within-subsample PS trimming [0.01, 0.99]
#          per tenure subset, then restricts treated grids to Santa Cruz or Beni
#          separately. Tests differential treatment intensity (SC private was
#          already permitted pre-decree; Beni was entirely new).
#
# Design decisions:
#   - Elevation <= 1800m for controls, within-subsample PS trimming [0.01, 0.99]
#     per tenure subset, dominant tenure >= 40% share, MIN_GRIDS = 50
#   - Protected tenure combines NPA and forest reserve shares before assignment
#   - PS estimated within each tenure subset (NOT on the full panel)
#   - Controls: all non-treated departments pooled (same control set for both)
#
# Inputs:
#   - data/grid_month_panel_2km.rds (via 1_shared_prep.R)
#
# Outputs:
#   Tables (output/tables/):
#     - did_drdid_fire_bydept_sc.tex
#     - did_drdid_fire_bydept_beni.tex
#     - did_drdid_landuse_bydept_sc.tex
#     - did_drdid_landuse_bydept_beni.tex
#
# RUNTIME NOTE: Expect several hours due to bootstrap SEs across 2 departments
#   x 6 tenure subsets x fire + land use outcomes.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

source("code/1_shared_prep.R")
library(did)

set.seed(20260525)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 0. Settings -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tenure_cats <- c("private_large", "private_small", "campesino", "tierra_fiscal",
                  "indigenous", "protected")
tenure_labels <- c("Lg. private", "Sm. private", "Communal", "Untitled public",
                    "Indigenous", "Protected")
names(tenure_labels) <- tenure_cats

MIN_GRIDS <- 50
ps_threshold <- 0.01

dept_loop <- c("Santa Cruz", "Beni")
dept_short <- c("Santa Cruz" = "sc", "Beni" = "beni")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1. Elevation filter + pre-treatment climate means -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Elevation filter already applied in bundled panel (controls <= 1800m).

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2. Dominant tenure assignment -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== Computing dominant tenure (40% threshold) ===")

dt <- merge(dt, compute_dominant_tenure(dt), by = "ID")

message("Panel after tenure assignment: ", uniqueN(dt$ID), " grids, ",
        format(nrow(dt), big.mark = ","), " obs")

# Print final counts per tenure
message("\nFinal grid counts per tenure subset:")
for (t in tenure_cats) {
  n_tr <- uniqueN(dt[dominant_tenure == t & treated == 1, ID])
  n_ct <- uniqueN(dt[dominant_tenure == t & treated == 0, ID])
  message("  ", tenure_labels[t], ": ", n_tr, " treated, ", n_ct, " control")
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3. Helpers -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xformla <- ~ elevation + slope + pct_forest2014 + log_pop
stopifnot(identical(all.vars(xformla),
                    c("elevation", "slope", "pct_forest2014", "log_pop")))
message("DR-DID covariates: ", paste(all.vars(xformla), collapse = ", "))

build_bydept_tex <- function(atts_dt, outcome_vec, label_vec, caption, label_tag,
                             filename, panel_type, treated_means = NULL,
                             n_treated_vec = NULL, n_control_vec = NULL,
                             dept_name = NULL) {
  n_tenures <- length(tenure_cats)
  col_header <- paste(paste0("& ", tenure_labels[tenure_cats]), collapse = " ")

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

  dept_note <- if (!is.null(dept_name)) {
    paste0("Treated grids restricted to ", dept_name, " department; ",
           "controls pool all non-treated departments. ")
  } else ""

  notes_text <- paste0(
    "\\item \\textit{Notes:} DR-DID ATT (Sant'Anna \\& Zhao 2020; \\texttt{did} package, Callaway \\& Sant'Anna 2021) ",
    "estimated separately for each tenure category. ",
    dept_note,
    "Dominant tenure = plurality (highest share $\\geq$ 40\\%); protected share ",
    "is the capped sum of NPA and forest reserve coverage. ",
    "Propensity score trimming and DR-DID nuisance models both include elevation, slope, 2014 forest share, and log population. ",
    "Symmetric PS trimming (Crump et al.\\ 2009): within each tenure subset, ",
    "controls with PS $<$ 0.01 and treated with PS $>$ 0.99 dropped; ",
    "elevation $\\leq$ 1800m for controls. ",
    panel_type, ". ",
    "SE clustered at 20km super-grid level (bootstrap). ",
    "$^{+}p<0.1$, $^{*}p<0.05$, $^{**}p<0.01$, $^{***}p<0.001$."
  )

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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4. Department loop -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for (dept in dept_loop) {

  dept_tag <- dept_short[dept]

  message("\n", strrep("#", 70))
  message("### DEPARTMENT: ", dept, " ###")
  message(strrep("#", 70))

  #--- 4a. Subset: treated from this dept + all controls ---
  dt_dept <- dt[treated == 0 | (treated == 1 & ADM1_ES == dept)]
  message("Grids: ",
          uniqueN(dt_dept$ID[dt_dept$treated == 1]), " treated (",
          dept, "), ",
          uniqueN(dt_dept$ID[dt_dept$treated == 0]), " control")

  # Print counts per tenure
  message("\nGrid counts per tenure:")
  for (t in tenure_cats) {
    n_tr <- uniqueN(dt_dept$ID[dt_dept$dominant_tenure == t & dt_dept$treated == 1])
    n_ct <- uniqueN(dt_dept$ID[dt_dept$dominant_tenure == t & dt_dept$treated == 0])
    message("  ", tenure_labels[t], ": ", n_tr, " treated, ", n_ct, " control")
  }

  #--- 4b. Tenure loop: DR-DID ---
  all_fire_atts <- list()
  all_lu_atts   <- list()
  fire_n_treated <- setNames(rep(NA_integer_, length(tenure_cats)), tenure_cats)
  fire_n_control <- setNames(rep(NA_integer_, length(tenure_cats)), tenure_cats)
  lu_n_treated   <- setNames(rep(NA_integer_, length(tenure_cats)), tenure_cats)
  lu_n_control   <- setNames(rep(NA_integer_, length(tenure_cats)), tenure_cats)

  for (tenure in tenure_cats) {
    message("\n", strrep("=", 60))
    message("  ", dept, " | ", tenure_labels[tenure], " (", tenure, ")")
    message(strrep("=", 60))

    dt_sub <- dt_dept[dominant_tenure == tenure]
    n_treated <- uniqueN(dt_sub$ID[dt_sub$treated == 1])
    n_control <- uniqueN(dt_sub$ID[dt_sub$treated == 0])
    message("  Grids: ", n_treated, " treated, ", n_control, " control")

    # Within-subsample PS trimming [0.01, 0.99] (same as script 4)
    grid_sub <- dt_sub[, .(
      treated = treated[1], elevation = elevation[1], slope = slope[1],
      pct_forest2014 = pct_forest2014[1], log_pop = log1p(population[1])
    ), by = ID]
    ps_sub <- glm(treated ~ elevation + slope + pct_forest2014 + log_pop,
                  data = grid_sub, family = binomial(link = "logit"))
    grid_sub[, ps := predict(ps_sub, newdata = grid_sub, type = "response")]
    drop_sub <- grid_sub[
      (treated == 0 & (is.na(ps) | ps < ps_threshold)) |
      (treated == 1 & (is.na(ps) | ps > (1 - ps_threshold))), ID]
    if (length(drop_sub) > 0) {
      dt_sub <- dt_sub[!(ID %in% drop_sub)]
      n_treated <- uniqueN(dt_sub$ID[dt_sub$treated == 1])
      n_control <- uniqueN(dt_sub$ID[dt_sub$treated == 0])
      message("    Within-subsample PS trimming: dropped ", length(drop_sub),
              " grids -> ", n_treated, " treated, ", n_control, " control")
    }
    rm(grid_sub, ps_sub, drop_sub)

    if (n_treated < MIN_GRIDS || n_control < MIN_GRIDS) {
      message("  WARNING: Skipping -- insufficient grids (need >= ", MIN_GRIDS, ")")
      all_fire_atts[[tenure]] <- data.table(
        tenure = tenure, outcome = fire_outcomes, att = NA_real_, se = NA_real_
      )
      all_lu_atts[[tenure]] <- data.table(
        tenure = tenure, outcome = lu_outcomes, att = NA_real_, se = NA_real_
      )
      next
    }

    #--- Fire season panel (Jul-Oct, 2015-2022) ---
    dt_fire_sub <- dt_sub[month %in% 7:10 & year %in% 2015:2022]

    fire_month_indices <- sort(unique(dt_fire_sub$month_index))
    dt_fire_sub[, fire_period := match(month_index, fire_month_indices)]

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
    fire_n_treated[tenure] <- uniqueN(dt_fire_sub$ID[dt_fire_sub$treated == 1])
    fire_n_control[tenure] <- uniqueN(dt_fire_sub$ID[dt_fire_sub$treated == 0])

    #--- Annual land use panel (2015-2022) ---
    dt_year_sub <- build_annual_panel(dt_sub[year %in% 2015:2022])
    dt_year_sub[, treatment_year := fifelse(treated == 1, 2019, 0)]

    dt_year_sub[, log_pop := log1p(population)]

    message("  Land use panel: ", uniqueN(dt_year_sub$ID), " grids x ",
            uniqueN(dt_year_sub$year), " years")

    df_year <- as.data.frame(dt_year_sub)
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
      }
    }
    all_lu_atts[[tenure]] <- tenure_lu_atts
  }

  #--- 4c. Build tables for this department ---
  message("\n=== Building tables for ", dept, " ===")

  fire_atts <- rbindlist(all_fire_atts)
  lu_atts   <- rbindlist(all_lu_atts)

  # Treated pre-treatment means
  fire_treat_means <- matrix(NA_real_, nrow = length(fire_outcomes),
                             ncol = length(tenure_cats),
                             dimnames = list(fire_outcomes, tenure_cats))
  lu_treat_means <- matrix(NA_real_, nrow = length(lu_outcomes),
                           ncol = length(tenure_cats),
                           dimnames = list(lu_outcomes, tenure_cats))

  for (tenure in tenure_cats) {
    dt_sub <- dt_dept[dominant_tenure == tenure]
    dt_fire_pre <- dt_sub[treated == 1 & month %in% 7:10 & year %in% 2015:2018]
    for (y in fire_outcomes) {
      if (nrow(dt_fire_pre) > 0) {
        fire_treat_means[y, tenure] <- mean(dt_fire_pre[[y]], na.rm = TRUE)
      }
    }
    dt_lu_pre <- build_annual_panel(dt_sub[treated == 1 & year %in% 2015:2018])
    for (y in lu_outcomes) {
      if (nrow(dt_lu_pre) > 0) {
        lu_treat_means[y, tenure] <- mean(dt_lu_pre[[y]], na.rm = TRUE)
      }
    }
  }

  dept_label <- dept

  build_bydept_tex(
    fire_atts, fire_outcomes, fire_labels,
    caption = paste0("DR-DID by Land Tenure: Fire Outcomes (", dept_label,
                     " Treated Only)"),
    label_tag = paste0("did_drdid_fire_bydept_", dept_tag),
    filename = paste0("did_drdid_fire_bydept_", dept_tag, ".tex"),
    panel_type = "Fire season (Jul--Oct), 2015--2022",
    treated_means = fire_treat_means,
    n_treated_vec = fire_n_treated,
    n_control_vec = fire_n_control,
    dept_name = dept_label
  )

  build_bydept_tex(
    lu_atts, lu_outcomes, lu_labels,
    caption = paste0("DR-DID by Land Tenure: Land Use Outcomes (", dept_label,
                     " Treated Only)"),
    label_tag = paste0("did_drdid_landuse_bydept_", dept_tag),
    filename = paste0("did_drdid_landuse_bydept_", dept_tag, ".tex"),
    panel_type = "Annual panel, 2015--2022",
    treated_means = lu_treat_means,
    n_treated_vec = lu_n_treated,
    n_control_vec = lu_n_control,
    dept_name = dept_label
  )

  #--- Summary for this department ---
  message("\n=== ", dept, " summary ===")
  for (t in tenure_cats) {
    fire_row <- fire_atts[tenure == t]
    n_est_fire <- sum(!is.na(fire_row$att))
    lu_row <- lu_atts[tenure == t]
    n_est_lu <- sum(!is.na(lu_row$att))
    message("  ", tenure_labels[t], ": ", n_est_fire, "/",
            length(fire_outcomes), " fire + ",
            n_est_lu, "/", length(lu_outcomes), " land use outcomes estimated")
  }
}

message("\n=== Script 5 (drdid_bydept) complete ===")
