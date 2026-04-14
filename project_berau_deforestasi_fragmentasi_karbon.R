# =====================================================================
# MINI PROJECT BERAU
# Analisis Deforestasi, Fragmentasi Hutan, dan Karbon
# =====================================================================

rm(list = ls())

# Install packages
install.packages("landscapemetrics")

library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(tidyr)
library(landscapemetrics)

# =====================================================================
# 01. IMPORT DAN CLIP DATA
# =====================================================================
# Working directory & output folder
setwd("E:/geosoftware - 2026/R - Bundle Spasial/Kelas 2 - spasial tematik R/project")

data_dir   <- "data"
output_dir <- "output"

dir.create(output_dir, showWarnings = FALSE)
dir.create(file.path(output_dir, "aoi"),    showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "raster"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "maps"),   showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "graphs"), showWarnings = FALSE, recursive = TRUE)

shp_path       <- file.path(data_dir, "IDN_ADM2_fixedInternalTopology.shp")
treecover_path <- file.path(data_dir, "Hansen_GFC-2024-v1.12_treecover2000_10N_110E.tif")
lossyear_path  <- file.path(data_dir, "Hansen_GFC-2024-v1.12_lossyear_10N_110E.tif")
datamask_path  <- file.path(data_dir, "Hansen_GFC-2024-v1.12_datamask_10N_110E.tif")
agb_path       <- file.path(data_dir, "N10E110_ESACCI-BIOMASS-L4-AGB-MERGED-100m-2022-fv6.0.tif")

# AOI Kabupaten Berau
adm2 <- st_read(shp_path)
print(names(adm2))

adm2_no_geom <- st_drop_geometry(adm2)
char_cols <- names(adm2_no_geom)[sapply(adm2_no_geom, is.character)]

match_matrix <- sapply(char_cols, function(col) {
  grepl("berau", adm2_no_geom[[col]], ignore.case = TRUE)
})

row_match <- apply(match_matrix, 1, any)
berau <- adm2[row_match, ][1, ]

st_write(
  berau,
  file.path(output_dir, "aoi", "berau_aoi.shp"),
  delete_dsn = TRUE,
  quiet = TRUE
)

# Preview AOI
par(mfrow = c(1, 1), mar = c(2, 2, 2, 2))
plot(st_geometry(berau), col = "lightgreen", main = "AOI Kabupaten Berau")

# Simpan AOI
png(file.path(output_dir, "maps", "01_aoi_berau.png"), width = 1200, height = 900, res = 150)
par(mfrow = c(1, 1), mar = c(2, 2, 2, 2))
plot(st_geometry(berau), col = "lightgreen", main = "AOI Kabupaten Berau")
dev.off()

# Baca & clip raster 
treecover <- rast(treecover_path)
lossyear  <- rast(lossyear_path)
datamask  <- rast(datamask_path)
agb       <- rast(agb_path)

berau_vect_gfc <- vect(st_transform(berau, crs(treecover)))
berau_vect_agb <- vect(st_transform(berau, crs(agb)))

treecover_berau <- crop(treecover, berau_vect_gfc) |> mask(berau_vect_gfc)
lossyear_berau  <- crop(lossyear, berau_vect_gfc)  |> mask(berau_vect_gfc)
datamask_berau  <- crop(datamask, berau_vect_gfc)  |> mask(berau_vect_gfc)
agb_berau       <- crop(agb, berau_vect_agb)       |> mask(berau_vect_agb)

writeRaster(treecover_berau, file.path(output_dir, "raster", "treecover2000_berau.tif"), overwrite = TRUE)
writeRaster(lossyear_berau,  file.path(output_dir, "raster", "lossyear_berau.tif"), overwrite = TRUE)
writeRaster(datamask_berau,  file.path(output_dir, "raster", "datamask_berau.tif"), overwrite = TRUE)
writeRaster(agb_berau,       file.path(output_dir, "raster", "agb_berau_2022.tif"), overwrite = TRUE)

cat("\nRaster hasil clip berhasil disimpan.\n")

# Preview raster awal
par(mfrow = c(2, 2), mar = c(3, 3, 3, 5))
plot(treecover_berau, main = "Tree Cover 2000 - Kabupaten Berau")
plot(lossyear_berau,  main = "Loss Year - Kabupaten Berau")
plot(datamask_berau,  main = "Data Mask - Kabupaten Berau")
plot(agb_berau,       main = "AGB 2022 - Kabupaten Berau")

# Threshold hutan
forest_30 <- ifel(datamask_berau == 1 & treecover_berau >= 30, 1, NA)
forest_75 <- ifel(datamask_berau == 1 & treecover_berau >= 75, 1, NA)
forest_90 <- ifel(datamask_berau == 1 & treecover_berau >= 90, 1, NA)

par(mfrow = c(2, 2), mar = c(3, 3, 3, 3))
plot(treecover_berau, main = "Persentase Tutupan Pohon 2000")
plot(forest_30, col = "darkgreen", legend = FALSE, main = "Forest Cover >= 30%")
plot(forest_75, col = "darkgreen", legend = FALSE, main = "Forest Cover >= 75%")
plot(forest_90, col = "darkgreen", legend = FALSE, main = "Forest Cover >= 90%")

# =====================================================================
# 01. ANALISIS DEFORESTASI TAHUNAN
# =====================================================================
forest_threshold <- 30
cell_area_ha <- cellSize(treecover_berau, unit = "ha")

forest_2000 <- ifel(datamask_berau == 1 & treecover_berau >= forest_threshold, 1, NA)
writeRaster(
  forest_2000, file.path(output_dir, "raster", "forest_mask_2000_berau.tif"), 
  overwrite = TRUE)

forest_2002 <- ifel(forest_2000 == 1 & (is.na(lossyear_berau) | lossyear_berau == 0 | lossyear_berau > 2), 1, NA)
forest_2012 <- ifel(forest_2000 == 1 & (is.na(lossyear_berau) | lossyear_berau == 0 | lossyear_berau > 12), 1, NA)
forest_2022 <- ifel(forest_2000 == 1 & (is.na(lossyear_berau) | lossyear_berau == 0 | lossyear_berau > 22), 1, NA)

writeRaster(forest_2002, file.path(output_dir, "raster", "forest_mask_2002_berau.tif"), overwrite = TRUE)
writeRaster(forest_2012, file.path(output_dir, "raster", "forest_mask_2012_berau.tif"), overwrite = TRUE)
writeRaster(forest_2022, file.path(output_dir, "raster", "forest_mask_2022_berau.tif"), overwrite = TRUE)

# Preview hutan lintas waktu
par(mfrow = c(1, 3), mar = c(3, 3, 3, 3))
plot(forest_2002, col = "darkgreen", legend = FALSE, main = "Forest 2002")
plot(forest_2012, col = "darkgreen", legend = FALSE, main = "Forest 2012")
plot(forest_2022, col = "darkgreen", legend = FALSE, main = "Forest 2022")

# Tabel deforestasi tahunan 2001-2024
# Luas hutan awal tahun 2000 (ha)
forest_2000_area_ha <- as.numeric(
  global(ifel(forest_2000 == 1, cell_area_ha, NA), "sum", na.rm = TRUE)[1, 1]
)

lossyear_forest <- ifel(forest_2000 == 1 & lossyear_berau > 0, lossyear_berau, NA)
loss_area_ha <- ifel(!is.na(lossyear_forest), cell_area_ha, NA)

# Hitung luas deforestasi per tahun dengan zonal
loss_zonal <- as.data.frame(
  zonal(loss_area_ha, lossyear_forest, fun = "sum", na.rm = TRUE)
)

if (nrow(loss_zonal) > 0) {
  names(loss_zonal)[1:2] <- c("loss_code", "annual_loss_ha")
}

defo_tbl <- data.frame(
  loss_code = 1:24,
  year = 2001:2024
)

defo_tbl <- left_join(defo_tbl, loss_zonal, by = "loss_code")
defo_tbl$annual_loss_ha[is.na(defo_tbl$annual_loss_ha)] <- 0
defo_tbl$cumulative_loss_ha <- cumsum(defo_tbl$annual_loss_ha)
defo_tbl$remaining_forest_ha <- forest_2000_area_ha - defo_tbl$cumulative_loss_ha

dir.create(file.path(output_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
write.csv(
  defo_tbl,
  file.path(output_dir, "tables", "deforestation_annual_berau.csv"),
  row.names = FALSE
)

defo_year_map <- ifel(!is.na(lossyear_forest), lossyear_forest + 2000, NA)
writeRaster(
  defo_year_map,
  file.path(output_dir, "raster", "deforestation_year_berau.tif"),
  overwrite = TRUE
)

# Peta tahun kehilangan hutan
defo_year_map <- ifel(forest_2000 == 1 & lossyear_berau > 0, lossyear_berau + 2000, NA)
writeRaster(defo_year_map, file.path(output_dir, "raster", "deforestation_year_berau.tif"), overwrite = TRUE)

par(mfrow = c(1, 1), mar = c(4, 4, 3, 5))
plot(defo_year_map,
     main = "Peta Deforestasi Kabupaten Berau (Tahun Kehilangan Hutan)",
     col = rainbow(24))
plot(berau_vect_gfc, add = TRUE, border = "black", lwd = 0.7)

png(file.path(output_dir, "maps", "02_deforestation_year_map_berau.png"),
    width = 1400, height = 1000, res = 150)
par(mfrow = c(1, 1), mar = c(4, 4, 3, 5))
plot(defo_year_map,
     main = "Peta Deforestasi Kabupaten Berau (Tahun Kehilangan Hutan)",
     col = rainbow(24))
plot(berau_vect_gfc, add = TRUE, border = "black", lwd = 0.7)
dev.off()

# Grafik tren deforestasi
p_loss <- ggplot(defo_tbl, aes(x = year, y = annual_loss_ha)) +
  geom_col(fill = "firebrick2") +
  theme_minimal() +
  labs(
    title = "Deforestasi Tahunan Kabupaten Berau",
    x = "Tahun",
    y = "Kehilangan hutan (ha)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_cover <- ggplot(defo_tbl, aes(x = year, y = remaining_forest_ha)) +
  geom_line(color = "darkgreen", linewidth = 1.2) +
  geom_point(color = "darkgreen", size = 1.8) +
  theme_minimal() +
  labs(
    title = "Luas Hutan Tersisa Kabupaten Berau",
    x = "Tahun",
    y = "Luas hutan (ha)"
  )

print(p_loss)
ggsave(
  filename = file.path(output_dir, "graphs", "03_deforestation_annual_loss_berau.png"),
  plot = p_loss, width = 10, height = 6, dpi = 150
)

print(p_cover)
ggsave(
  filename = file.path(output_dir, "graphs", "04_remaining_forest_trend_berau.png"),
  plot = p_cover, width = 10, height = 6, dpi = 150
)

cat("\nTahap deforestasi selesai.\n")

# =====================================================================
# 02. ANALISIS FRAGMENTASI
# =====================================================================
crs_utm <- "EPSG:32650"

# Klasifikasi hutan vs non-hutan
forest_cat_2002 <- ifel(datamask_berau == 1, ifel(!is.na(forest_2002), 1, 2), NA)
forest_cat_2012 <- ifel(datamask_berau == 1, ifel(!is.na(forest_2012), 1, 2), NA)
forest_cat_2022 <- ifel(datamask_berau == 1, ifel(!is.na(forest_2022), 1, 2), NA)

# Proyeksi ke UTM
forest_cat_2002_utm <- as.int(project(forest_cat_2002, crs_utm, method = "near"))
forest_cat_2012_utm <- as.int(project(forest_cat_2012, crs_utm, method = "near"))
forest_cat_2022_utm <- as.int(project(forest_cat_2022, crs_utm, method = "near"))

# Lima metrik fragmentasi
metric_forest <- c(
  "lsm_c_np",
  "lsm_c_pd",
  "lsm_c_area_mn",
  "lsm_c_lpi",
  "lsm_c_ed"
)

# Hitung metrik tiap tahun
frag_2002_metrics <- calculate_lsm(
  forest_cat_2002_utm,
  what = metric_forest,
  directions = 4,
  classes = 1
) %>%
  mutate(year = 2002)

frag_2012_metrics <- calculate_lsm(
  forest_cat_2012_utm,
  what = metric_forest,
  directions = 4,
  classes = 1
) %>%
  mutate(year = 2012)

frag_2022_metrics <- calculate_lsm(
  forest_cat_2022_utm,
  what = metric_forest,
  directions = 4,
  classes = 1
) %>%
  mutate(year = 2022)

fragmentation_metrics_long <- bind_rows(
  frag_2002_metrics,
  frag_2012_metrics,
  frag_2022_metrics
)

fragmentation_metrics_summary <- fragmentation_metrics_long %>%
  select(year, metric, value) %>%
  group_by(year, metric) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = metric, values_from = value)

write.csv(
  fragmentation_metrics_summary,
  file.path(output_dir, "tables", "fragmentation_metrics_summary_berau.csv"),
  row.names = FALSE
)
print(fragmentation_metrics_summary)

# Peta patch untuk visualisasi
forest_only_2002_utm <- ifel(forest_cat_2002_utm == 1, 1, NA)
forest_only_2022_utm <- ifel(forest_cat_2022_utm == 1, 1, NA)

forest_only_2002_map <- aggregate(forest_only_2002_utm, fact = 2, fun = "modal", na.rm = TRUE)
forest_only_2022_map <- aggregate(forest_only_2022_utm, fact = 2, fun = "modal", na.rm = TRUE)

patch_2002 <- patches(forest_only_2002_map, directions = 4)
patch_2022 <- patches(forest_only_2022_map, directions = 4)

writeRaster(patch_2002, file.path(output_dir, "raster", "patches_2002_berau.tif"), overwrite = TRUE)
writeRaster(patch_2022, file.path(output_dir, "raster", "patches_2022_berau.tif"), overwrite = TRUE)

par(mfrow = c(1, 2), mar = c(3, 3, 3, 3))
plot(patch_2002, legend = FALSE, col = terrain.colors(60), main = "Patch Hutan 2002")
plot(patch_2022, legend = FALSE, col = terrain.colors(60), main = "Patch Hutan 2022")

cat("\nTahap fragmentasi selesai.\n")

# =====================================================================
# 03. Analisis Karbon
# =====================================================================
agb_utm <- project(agb_berau, crs_utm, method = "bilinear")
forest_2022_utm <- project(forest_2022, crs_utm, method = "near")

# Samakan grid AGB dengan forest 2022
agb_match <- project(agb_utm, forest_2022_utm, method = "bilinear")

# Mask hanya hutan 2022
agb_forest_2022 <- mask(agb_match, forest_2022_utm)

# Konversi biomassa ke karbon
carbon_2022 <- agb_forest_2022 * 0.47
names(carbon_2022) <- "Carbon_MgC_ha"

# Luas piksel untuk total karbon
pixel_area_carbon_ha <- cellSize(carbon_2022, unit = "ha")
carbon_total_per_cell <- carbon_2022 * pixel_area_carbon_ha

total_carbon_ton <- r_sum(carbon_total_per_cell)
mean_carbon_density_ton_ha <- global(carbon_2022, "mean", na.rm = TRUE)[1, 1]

carbon_summary_tbl <- data.frame(
  year = 2022,
  forest_area_ha = r_sum(ifel(forest_2022 == 1, cell_area_ha, NA)),
  mean_carbon_density_ton_ha = as.numeric(mean_carbon_density_ton_ha),
  total_carbon_ton = total_carbon_ton
)

write.csv(
  carbon_summary_tbl,
  file.path(output_dir, "tables", "carbon_summary_berau.csv"),
  row.names = FALSE
)

writeRaster(agb_match,             file.path(output_dir, "raster", "agb_2022_match_berau.tif"), overwrite = TRUE)
writeRaster(carbon_2022,           file.path(output_dir, "raster", "carbon_density_2022_berau.tif"), overwrite = TRUE)
writeRaster(carbon_total_per_cell, file.path(output_dir, "raster", "carbon_total_per_cell_2022_berau.tif"), overwrite = TRUE)

cat("\nRingkasan karbon:\n")
print(carbon_summary_tbl)

# Area karbon tinggi > 35 MgC/ha
carbon_high <- ifel(carbon_2022 > 35, 1, NA)
carbon_high_patch <- patches(carbon_high, directions = 4)

carbon_high_area_table <- as.data.frame(
  zonal(pixel_area_carbon_ha, carbon_high_patch, fun = "sum", na.rm = TRUE)
)

if (nrow(carbon_high_area_table) > 0) {
  names(carbon_high_area_table)[1:2] <- c("patch_id", "patch_area_ha")
  carbon_patch_ge500 <- carbon_high_area_table$patch_id[carbon_high_area_table$patch_area_ha >= 500]
  carbon_high_large <- ifel(carbon_high_patch %in% carbon_patch_ge500, 1, NA)
} else {
  carbon_high_large <- carbon_high_patch
}

writeRaster(carbon_high,       file.path(output_dir, "raster", "carbon_high_gt35_berau.tif"), overwrite = TRUE)
writeRaster(carbon_high_large, file.path(output_dir, "raster", "carbon_high_largepatch_berau.tif"), overwrite = TRUE)

par(mfrow = c(2, 2), mar = c(3, 3, 3, 5))
plot(agb_match,         main = "Biomassa 2022 (AGB)")
plot(forest_2022_utm,   col = "darkgreen", legend = FALSE, main = "Hutan 2022")
plot(carbon_2022,       main = "Carbon Stock (MgC/ha)")
plot(carbon_high_large, main = "Patch Karbon Tinggi >= 500 ha")
save_current_plot(file.path(output_dir, "maps", "11_carbon_maps_berau.png"), width = 1800, height = 1200)

par(mfrow = c(1, 1), mar = c(4, 4, 3, 5))
plot(carbon_2022, main = "Peta Stok Karbon Hutan Tersisa Berau", col = hcl.colors(20, "YlGn"))
save_current_plot(file.path(output_dir, "maps", "12_carbon_stock_2022_berau.png"), width = 1400, height = 1000)

cat("\nTahap karbon selesai.\n")

# =====================================================================
# 10. RINGKASAN AKHIR MINI PROJECT
# =====================================================================
forest_area_2002_ha <- r_sum(ifel(forest_2002 == 1, cell_area_ha, NA))
forest_area_2012_ha <- r_sum(ifel(forest_2012 == 1, cell_area_ha, NA))
forest_area_2022_ha <- r_sum(ifel(forest_2022 == 1, cell_area_ha, NA))
total_deforestation_ha <- sum(defo_tbl$annual_loss_ha, na.rm = TRUE)
high_carbon_area_ha <- r_sum(ifel(carbon_high_large == 1, pixel_area_carbon_ha, NA))

final_summary_tbl <- data.frame(
  wilayah_studi = "Berau, Kalimantan Timur",
  forest_threshold_pct = forest_threshold,
  forest_area_2002_ha = forest_area_2002_ha,
  forest_area_2012_ha = forest_area_2012_ha,
  forest_area_2022_ha = forest_area_2022_ha,
  total_deforestation_2001_2024_ha = total_deforestation_ha,
  total_carbon_ton_2022 = total_carbon_ton,
  mean_carbon_density_ton_ha_2022 = as.numeric(mean_carbon_density_ton_ha),
  high_carbon_area_ge500ha_ha = high_carbon_area_ha
)

write.csv(
  final_summary_tbl,
  file.path(output_dir, "tables", "final_summary_berau.csv"),
  row.names = FALSE
)

cat("\n==================== SELESAI ====================\n")
cat("Mini project selesai diproses.\n")
cat("Cek folder output untuk tabel, grafik, peta, dan raster hasil.\n\n")

cat("Output utama:\n")
cat("- output/tables/deforestation_annual_berau.csv\n")
cat("- output/tables/fragmentation_metrics_summary_berau.csv\n")
cat("- output/tables/carbon_summary_berau.csv\n")
cat("- output/tables/final_summary_berau.csv\n")
cat("- output/maps/05_deforestation_year_map_berau.png\n")
cat("- output/maps/10_patch_map_2002_2022_berau.png\n")
cat("- output/maps/12_carbon_stock_2022_berau.png\n")
cat("================================================\n")