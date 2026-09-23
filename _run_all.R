#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Project: Bolivia Fires
# Purpose: Master replication script. Runs all locally-executable analysis
#          scripts in order and produces all tables and figures.
# INSTRUCTIONS:
#   1. Set your working directory to this folder (replication_package/)
#      before running:
#        setwd("/path/to/replication_package")
#   2. Run this script (or individual code/ scripts) from that directory.
#   3. Outputs land in output/figures/ and output/tables/.
#
# REQUIREMENTS:
#   R packages: data.table, fixest, MatchIt, did, sf, ggplot2, patchwork,
#               tidyverse, kableExtra, scales
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Confirm working directory
cat("Working directory:", getwd(), "\n")
if (!file.exists("data/grid_month_panel_2km.rds")) {
  stop("ERROR: Working directory must be set to the replication_package/ folder. ",
       "data/grid_month_panel_2km.rds not found.")
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step 0: PSM matching + pre-trends figure -----
#         Outputs: data/grid_psm_weights_2km.rds
#                  output/figures/did_pretrends_fire.pdf
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cat("\n=== Running 0_psm_prep.R ===\n")
source("code/0_psm_prep.R")
cat("Script 0 complete.\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step 2: TWFE -----
#         Outputs: output/tables/did_fire_n_fires.tex
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cat("\n=== Running 2_twfe.R ===\n")
source("code/2_twfe.R")
cat("Script 2 complete.\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step 3: DR-DID main -----
#         Outputs: output/figures/did_drdid_trimmed_sample.pdf
#                  output/figures/did_drdid_eventstudy_*.pdf
#                  output/tables/did_balance_ps_trimmed.tex
#                  output/tables/did_drdid_fire.tex
#                  output/tables/did_drdid_landuse.tex
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cat("\n=== Running 3_drdid_main.R ===\n")
source("code/3_drdid_main.R")
cat("Script 3 complete.\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step 4: DR-DID by land tenure -----
#         Outputs: output/figures/did_drdid_byland_sample_map.pdf
#                  output/figures/did_drdid_byland_lu_coefplot.pdf
#                  output/figures/did_drdid_eventstudy_byland_*.pdf
#                  output/tables/did_drdid_fire_byland.tex
#                  output/tables/did_drdid_landuse_byland.tex
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cat("\n=== Running 4_drdid_byland.R ===\n")
source("code/4_drdid_byland.R")
cat("Script 4 complete.\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step 5: DR-DID by department x land tenure -----
#         Outputs: output/tables/did_drdid_fire_bydept_{sc,beni}.tex
#                  output/tables/did_drdid_landuse_bydept_{sc,beni}.tex
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cat("\n=== Running 5_drdid_bydept.R ===\n")
source("code/5_drdid_bydept.R")
cat("Script 5 complete.\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step 6: DR-DID decree-period only -----
#         Outputs: output/tables/did_drdid_fire_decree.tex
#                  output/tables/did_drdid_landuse_decree.tex
#                  output/tables/did_drdid_fire_byland_decree.tex
#                  output/tables/did_drdid_landuse_byland_decree.tex
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cat("\n=== Running 6_drdid_decree_only.R ===\n")
source("code/6_drdid_decree_only.R")
cat("Script 6 complete.\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step 7: Maps -----
#         Outputs: output/figures/did_fire_deforestation_map.pdf
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cat("\n=== Running 7_maps.R ===\n")
source("code/7_maps.R")
cat("Script 7 complete.\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summary -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cat("\n", strrep("=", 70), "\n")
cat("All scripts complete. Outputs in:\n")
cat("  output/figures/  --", length(list.files("output/figures", pattern = "\\.pdf$")),
    "figures\n")
cat("  output/tables/   --", length(list.files("output/tables", pattern = "\\.tex$")),
    "tables\n")
cat(strrep("=", 70), "\n")

cat("\nNOTE: SI Figure S3 (fire ignition by tenure) is not included in this package.\n")
cat("  It requires ~20GB of fire detection data processed on a 256GB cloud instance.\n")
