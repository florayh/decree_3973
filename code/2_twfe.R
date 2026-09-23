#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Project: Bolivia Fires
# Purpose: TWFE and event study estimation of Decree 3973 effects on fires
#          and land use. Runs full and PSM-matched specifications with and
#          without climate controls. Fire outcomes restricted to fire season
#          months (Jul-Oct). Includes Poisson robustness and heterogeneity by
#          land tenure.
#
# Inputs:
#   - data/grid_month_panel_2km.rds
#   - data/grid_psm_weights_2km.rds
#
# Outputs:
#   Tables (output/tables/):
#     - did_fire_n_fires.tex
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

source("code/1_shared_prep.R")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3. TWFE -- Fire outcomes (monthly) -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== TWFE: Fire count (n_fires) ===")

fml1 <- n_fires ~ treated:post_decree | ID + year_month
fml2 <- as.formula(paste0("n_fires ~ treated:post_decree + ",
  climate_controls, " | ID + year_month"))

mods <- list(
  "(1)" = feols(fml1, data = dt_fire, cluster = ~ADM3_PCODE),
  "(2)" = feols(fml2, data = dt_fire, cluster = ~ADM3_PCODE),
  "(3)" = feols(fml1, data = dt_fire[psm_matched == TRUE],
                weights = ~psm_weight, cluster = ~ADM3_PCODE),
  "(4)" = feols(fml2, data = dt_fire[psm_matched == TRUE],
                weights = ~psm_weight, cluster = ~ADM3_PCODE)
)

treat_mean_full <- sprintf("%.3f",
  mean(dt_fire[treated == 1 & post_decree == 0, n_fires], na.rm = TRUE))
treat_mean_psm <- sprintf("%.3f",
  mean(dt_fire[treated == 1 & post_decree == 0 & psm_matched == TRUE, n_fires],
       na.rm = TRUE))

etable(mods,
       tex = TRUE,
       file = file.path(tab_dir, "did_fire_n_fires.tex"),
       replace = TRUE,
       digits = 3,
       title = "TWFE: Total fires",
       label = "tab:did_fire_n_fires",
       dict = c(
         "treated:post_decree" = "Treated $\\times$ Post",
         "temperature_2m"      = "Temperature",
         "chirps_precip_mm"    = "Precipitation",
         "wind_speed_10m"      = "Wind speed"
       ),
       keep = c("Treated"),
       fitstat = ~ n + r2,
       style.tex = style.tex("aer",
                             depvar.title = "",
                             fixef.title = "",
                             fixef.suffix = " FE",
                             yesNo = c("Yes", "No")),
       signif.code = c("***" = 0.001, "**" = 0.01, "*" = 0.05, "+" = 0.1),
       postprocess.tex = function(x) sub("\\\\centering", "\\\\centering\n\\\\small", x),
       extralines = list(
         "Treated mean (pre)" = c(treat_mean_full, treat_mean_full,
                                  treat_mean_psm, treat_mean_psm),
         "Climate controls"   = c("No", "Yes", "No", "Yes"),
         "Sample"             = c("Full", "Full", "PSM", "PSM")
       ),
       notes = "Fire season months (Jul--Oct) only. SE clustered at municipality level. PSM specifications use nearest-neighbor 1:1 matching with caliper = 0.2 SD.")

message("Fire table exported: ", file.path(tab_dir, "did_fire_n_fires.tex"))

message("\n=== Script 2 complete ===")
