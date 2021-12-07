library(sp)
library(sf)
library(raster)
library(tidyverse)
library(ggspatial)
library(rgdal)

# land cover ####
land_cover <- raster("data/land-cover-data_2021-10-15_1503/kaza_landCover_2016_project.tif")
e <- extent(fence_raster) %>% as("SpatialPolygons")
proj4string(e) <- crs(fence_raster)
e_rp <- spTransform(e, crs(land_cover))
ld_cov_crop <- crop(land_cover, e_rp)
ld_cov_crop_rp <- projectRaster(ld_cov_crop, crs = crs(fence_raster))
writeRaster(ld_cov_crop_rp, filename = "data/kaza_croped.tif")
ld_cov_crop_rp <- raster("data/kaza_croped.tif")

plot(ld_cov_crop_rp) 
ant_devel <- c(9:13, 21, 24, 29) # anthropogenic development
water <- c(1, 2, 3, 28) # water
wood <- c(7, 8) # woodland/forest
bare <- c(19, 20, 25, 26) # bare ground
open <- c(14:16) # open woodland/thicket

m_change <- cbind(is = c(water, wood, ant_devel, open, bare), 
                  becomes = c( rep(min(water), length(water)), rep(min(wood), length(wood)), 
                               rep(min(ant_devel), length(ant_devel)), rep(min(open), length(open)), 
                               rep(min(bare), length(bare)) ))
ld_cov_crop_rp_reclf <- reclassify(ld_cov_crop_rp, m_change)
ld_cov_crop_rp_reclf[is.na(ld_cov_crop_rp_reclf[])] <- 0
writeRaster(ld_cov_crop_rp_reclf, filename = "data/kaza_croped_reclaf.tif"
            , overwrite = T
            )

plot(ld_cov_crop_rp_reclf)

# fence ####
fence_sh <- readOGR("data/shapes/merged_fences_old and 2021 AA5FEB2021.shp") 
fence_sh %>% plot()

empty_raster <- ld_cov_crop_rp_reclf
values(empty_raster) <- 1

# buffer
fence_raster <- buffer(fence_sh, width = .005) %>% # 0.001 degrees ~ 111m
  rasterize(., empty_raster, field = 0, background = 1)

plot(fence_raster)
fence_sh %>% plot(add = T, col = "red", lty = 2)
scalebar(50, type = "bar", divs = 4, below = "km")

writeRaster(fence_raster, filename = "data/fence_raster.tif"
            , overwrite = T
            )

# roads ####
roads_sh <- readOGR("data/shapes/roads_main_GCS.shp")
roads_sh %>% plot()

# buffer
road_raster <- buffer(roads_sh, width = .005) %>% 
  rasterize(., empty_raster, field = 0, background = 1)
plot(road_raster)
roads_sh %>% plot(add = T, col = "lightgreen", lty = 3)
scalebar(50, type = "bar", divs = 4)

writeRaster(road_raster, filename = "data/road_raster.tif"
            #, overwrite = T
            )
# primary trunk ####
trunk_sh <- readOGR("data/shapes/primary_trunk_rds_5countries.shp")
trunk_sh_crop <- crop(trunk_sh, extent(ld_cov_crop_rp_reclf))
trunk_sh_crop %>% plot(col = "darkviolet")
plot(rworldmap::countriesLow, add = T, lty = 2)

road_trunk <- buffer(trunk_sh_crop, width = .005) %>% 
  rasterize(., empty_raster, field = 0, background = 1)
plot(road_trunk)
trunk_sh_crop %>% plot(add = T, col = "coral", lty = 3)
scalebar(50, type = "bar", divs = 4)

writeRaster(road_trunk, filename = "data/trunk_raster.tif"
            , overwrite = T
            )

# rivers ####
rivers_sh <- readOGR("data/shapes/rivers_main_KAZA_dig.shp")
rivers_sh %>% plot(col = "steelblue")
rivers_sh_crop <- crop(rivers_sh, extent(ld_cov_crop_rp_reclf))
rivers_sh_crop %>% plot(col = "steelblue", lwd = 2)
plot(rworldmap::countriesLow, add = T, lty = 2)

# protected areas ####
pas_sh <- readOGR("data/shapes/PAs_main_KAZA.shp")
pas_sh %>% plot(col = "forestgreen", add = T)
e <- extent(fence_raster) %>% as("SpatialPolygons")
r <- raster(ncol = fence_raster@ncols, nrow = fence_raster@nrows)
extent(r) <- e
pas_sh@data$Text_20 <- 1
pas_ras <- rasterize(pas_sh, r, field = "Text_20", background = 0)
plot(pas_ras)
writeRaster(pas_ras, filename = "data/pas_raster.tif"
            #, overwrite = T
)
pas_ras <- raster("data/pas_raster.tif")

# all elements ####
plot(ld_cov_crop_rp_reclf)
pas_sh %>% plot(add = T, col = "plum")
fence_sh %>% plot(add = T, col = "red", lty = 2)
#roads_sh %>% plot(add = T, col = "lightgreen", lty = 3)
trunk_sh %>% plot(col = "darkviolet", add = T)
rivers_sh_crop %>% plot(col = "steelblue", lwd = 2, add = T)
plot(rworldmap::countriesLow, add = T, lty = 2)
scalebar(250, type = "bar", divs = 4)
