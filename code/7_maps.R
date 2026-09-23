#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Project: Bolivia Fires
# Purpose: Two-panel maps for DID descriptive section.
#          Map 1 -- Fire density: left = mean fires per fire season (2015-2018,
#            Jul-Oct), right = 2019 fire season.
#          Map 2 -- Primary forest loss: left = mean annual loss (2015-2018),
#            right = 2019.
#          Shows lowland grids (elevation <= 1800m), no PS trimming.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

source("code/1_shared_prep.R")
library(sf)
library(patchwork)

log1p_trans <- scales::trans_new("log1p", transform = log1p, inverse = expm1,
                                 domain = c(0, Inf))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1. Load panel (elevation filter pre-applied) -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Elevation filter already applied in bundled panel (controls <= 1800m).

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2. Compute fire density by grid -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== Computing fire density ===")

# Fire season subset (Jul-Oct)
dt_fs <- dt[month %in% 7:10]

# Pre-period: mean annual fire count across 2015-2018 fire seasons
pre <- dt_fs[year %in% 2015:2018, .(total_fires = sum(n_fires, na.rm = TRUE)), by = .(ID, year)]
pre_mean <- pre[, .(fires_pre = mean(total_fires, na.rm = TRUE)), by = ID]

# 2019 fire season
post <- dt_fs[year == 2019, .(fires_2019 = sum(n_fires, na.rm = TRUE)), by = ID]

# Merge
fire_grid <- merge(pre_mean, post, by = "ID", all = TRUE)
fire_grid[is.na(fires_pre), fires_pre := 0]
fire_grid[is.na(fires_2019), fires_2019 := 0]

# Add treatment indicator
fire_grid <- merge(fire_grid,
                   unique(dt[, .(ID, treated)]),
                   by = "ID")

message("Grids: ", nrow(fire_grid))
message("Pre-period mean fires -- treated median: ",
        round(median(fire_grid[treated == 1]$fires_pre), 2),
        ", control median: ",
        round(median(fire_grid[treated == 0]$fires_pre), 2))
message("2019 fires -- treated median: ",
        round(median(fire_grid[treated == 1]$fires_2019), 2),
        ", control median: ",
        round(median(fire_grid[treated == 0]$fires_2019), 2))

rm(dt_fs, pre, pre_mean, post)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3. Load shapefiles -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== Loading shapefiles ===")

grid_sf <- st_read(file.path("data", "shapefiles/grids/grid_bolivia_2km.shp"),
                   quiet = TRUE)
dept_sf <- st_read(file.path("data", "shapefiles/adm1",
                              "bol_admbnda_adm1_gov_2020514.shp"),
                   quiet = TRUE)

# Mark treated departments for thicker boundary
treated_dept_names <- c("Santa Cruz", "Beni")
dept_sf$is_treated <- dept_sf$ADM1_ES %in% treated_dept_names

# Common theme
map_theme <- theme_minimal(base_size = 11) +
  theme(axis.text       = element_blank(),
        axis.ticks      = element_blank(),
        axis.title      = element_blank(),
        panel.grid      = element_blank(),
        legend.position = "bottom",
        legend.key.width  = unit(1.5, "cm"),
        legend.key.height = unit(0.3, "cm"),
        plot.title      = element_text(face = "bold", size = 11))

# Shared department boundary layers
dept_layers <- list(
  geom_sf(data = dept_sf[!dept_sf$is_treated, ],
          fill = NA, color = "grey40", linewidth = 0.3),
  geom_sf(data = dept_sf[dept_sf$is_treated, ],
          fill = NA, color = "white", linewidth = 0.8)
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4. Fire density map -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== Building fire density map ===")

fire_map_sf <- merge(grid_sf, fire_grid, by = "ID", all.x = FALSE)

max_fire <- max(c(fire_grid$fires_pre, fire_grid$fires_2019), na.rm = TRUE)

p_fire_pre <- ggplot() +
  geom_sf(data = fire_map_sf, aes(fill = fires_pre),
          color = NA, linewidth = 0) +
  dept_layers +
  scale_fill_viridis_c(option = "inferno",
                       limits = c(0, max_fire),
                       trans = log1p_trans,
                       name = "Fire count",
                       breaks = c(0, 1, 5, 20, 100, 500)) +
  labs(title = "A. Mean fires per fire season (2015\u20132018)") +
  map_theme

p_fire_post <- ggplot() +
  geom_sf(data = fire_map_sf, aes(fill = fires_2019),
          color = NA, linewidth = 0) +
  dept_layers +
  scale_fill_viridis_c(option = "inferno",
                       limits = c(0, max_fire),
                       trans = log1p_trans,
                       name = "Fire count",
                       breaks = c(0, 1, 5, 20, 100, 500)) +
  labs(title = "B. Total fires during 2019 fire season") +
  map_theme

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5. Primary forest loss map -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message("\n=== Building primary forest loss map ===")

# primary_loss_km2 is annual (constant within year), so take one value per grid-year
loss_yr <- dt[, .(primary_loss_km2 = primary_loss_km2[1]), by = .(ID, year)]

# Pre-period: mean annual loss across 2015-2018
loss_pre <- loss_yr[year %in% 2015:2018,
                    .(loss_pre = mean(primary_loss_km2, na.rm = TRUE)), by = ID]

# 2019
loss_post <- loss_yr[year == 2019,
                     .(loss_2019 = primary_loss_km2), by = ID]

loss_grid <- merge(loss_pre, loss_post, by = "ID", all = TRUE)
loss_grid[is.na(loss_pre), loss_pre := 0]
loss_grid[is.na(loss_2019), loss_2019 := 0]

# Add treatment indicator
loss_grid <- merge(loss_grid,
                   unique(dt[, .(ID, treated)]),
                   by = "ID")

message("Grids: ", nrow(loss_grid))
message("Pre-period mean loss (km2) -- treated median: ",
        round(median(loss_grid[treated == 1]$loss_pre), 4),
        ", control median: ",
        round(median(loss_grid[treated == 0]$loss_pre), 4))
message("2019 loss (km2) -- treated median: ",
        round(median(loss_grid[treated == 1]$loss_2019), 4),
        ", control median: ",
        round(median(loss_grid[treated == 0]$loss_2019), 4))

loss_map_sf <- merge(grid_sf, loss_grid, by = "ID", all.x = FALSE)

max_loss <- max(c(loss_grid$loss_pre, loss_grid$loss_2019), na.rm = TRUE)

p_loss_pre <- ggplot() +
  geom_sf(data = loss_map_sf, aes(fill = loss_pre),
          color = NA, linewidth = 0) +
  dept_layers +
  scale_fill_viridis_c(option = "mako",
                       limits = c(0, max_loss),
                       trans = log1p_trans,
                       name = expression("km"^2),
                       breaks = c(0, 0.1, 0.5, 2, 4)) +
  labs(title = "C. Mean annual primary forest loss (2015\u20132018)") +
  map_theme

p_loss_post <- ggplot() +
  geom_sf(data = loss_map_sf, aes(fill = loss_2019),
          color = NA, linewidth = 0) +
  dept_layers +
  scale_fill_viridis_c(option = "mako",
                       limits = c(0, max_loss),
                       trans = log1p_trans,
                       name = expression("km"^2),
                       breaks = c(0, 0.1, 0.5, 2, 4)) +
  labs(title = "D. Primary forest loss in 2019") +
  map_theme

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 6. Combine and export -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Top row: fire (shared inferno legend), bottom row: forest loss (shared mako legend)
p_fire_row <- (p_fire_pre + p_fire_post) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

p_loss_row <- (p_loss_pre + p_loss_post) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

p_combined <- p_fire_row / p_loss_row

out_file <- file.path(fig_dir, "did_fire_deforestation_map.pdf")
ggsave(out_file, p_combined, width = 12, height = 12, device = cairo_pdf)
message("Exported: ", out_file)

message("\n=== Script 7 (maps) complete ===")
