library(sp)
library(raster)
library(tidyverse)
library(ggspatial)
library(rgdal)

# fence ####
fence_sh <- readOGR("data/shapes/merged_fences_old and 2021 AA5FEB2021.shp") 
fence_sh %>% plot()

empty_raster <- ld_cov_crop
values(empty_raster) <- 1
plot(empty_raster)

# buffer
fence_raster <- buffer(fence_sh, width = .1) %>% 
  rasterize(., empty_raster, field = 0, background = 1)

plot(fence_raster)
fence_sh %>% plot(add = T)
scalebar(100, type = "bar", divs = 4)

writeRaster(fence_raster, filename = "data/fence_raster.tif")

# roads ####
roads_sh <- readOGR("data/shapes/roads_main_GCS.shp")
roads_sh %>% plot()

# buffer
road_raster <- buffer(roads_sh, width = .1) %>% 
  rasterize(., empty_raster, field = 0, background = 1)
plot(road_raster)
roads_sh %>% plot(add = T, col = "red" )

writeRaster(road_raster, filename = "data/road_raster.tif")

# rivers ####
rivers_sh <- readOGR("data/shapes/primary_trunk_rds_5countries.shp")
rivers_sh %>% plot(col = "steelblue")
rivers_sh_crop <- crop(rivers_sh, extent(ld_cov_crop))
rivers_sh_crop %>% plot(col = "steelblue", lwd = 2)
plot(rworldmap::countriesLow, add = T, lty = 2)

# BioClim ####
bc_data <- getData("worldclim", var = "bio", res = 0.5, lon = 22, lat = -18.5)
bc_data_crop <- crop(bc_data, extent(ld_cov_crop))

writeRaster(bc_data_crop, 
            filename = "data/bioclim_rasters.tif", format = "GTiff",
            options = c("INTERLEAVE=BAND","COMPRESS=LZW"), overwrite = T)

# MODIS rasters ####
library(MODIStsp)
MODIStsp_get_prodlayers("M*D13A2")
spat_file <- system.file("data/shapes/merged_fences_old and 2021 AA5FEB2021.shp", 
                         package = "MODIStsp")
MODIStsp(
  gui = FALSE,
  out_folder = "data/",
  selprod = "Vegetation_Indexes_16Days_1Km (M*D13A2)",
  bandsel = c("EVI", "NDVI"),
  quality_bandsel = "QA_usef",
  indexes_bandsel = "SR",
  spatmeth = "file",
  spafile = "data/shapes/merged_fences_old and 2021 AA5FEB2021.shp",
  user = "mstp_test" ,
  password = "MSTP_test_01",
  start_date = "2020.06.01",
  end_date = "2020.06.15",
  verbose = FALSE,
  parallel = FALSE
)

ndvi_15 <- raster("data/merged_fences_old and 2021 AA5FEB2021/VI_16Days_1Km_v6/NDVI/MOD13A2_NDVI_2020_161.tif")
ndvi_15_rp <- projectRaster(ndvi_15, crs = crs(bc_data_crop))
plot(ndvi_15_rp)
