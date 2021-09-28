library(tidyverse)

raw_data <- read.csv("data/Namibia_elephant_2010-Mar2021.csv")
raw_data %>% head()

mod_data <- raw_data %>%
  mutate(date.time = as.POSIXct(date.time, format = "%Y-%m-%d %H:%M", tz = "UTC"))

# to check data and numbers ####
# number of registers by sex through years
mod_data %>% group_by(year, sex) %>% count() %>% 
  ggplot(aes(x = year, y = log(`n`), col = sex)) + geom_line() + 
  theme_classic() +
  scale_color_brewer(palette = "Set1")

# number of registers by id and sex through years
mod_data %>% group_by(year, id.new, sex) %>% count() %>% 
  ggplot(aes(x = year, y = `n`, col = id.new)) + geom_line() + 
  theme_classic() +
  #scale_color_brewer(palette = "Set1") + 
  theme(legend.position = "none") +
    facet_wrap(~sex, nrow = 2)

# maps ####
library(ggmap)
library(ggspatial)
library(raster)
library(rgdal)

# fence shape and raster
fence_sh <- readOGR("data/shapes/merged_fences_old and 2021 AA5FEB2021.shp") 
fence_sh %>% plot

fence_raster <- raster("data/fence_raster.tif")

# roads shape and raster
roads_sh <- readOGR("data/shapes/roads_main_GCS.shp")
roads_sh %>% plot()

road_raster <- raster("data/road_raster.tif")

# land cover
land_cover <- raster("data/land_cover/kaza_landCover_2005.img")

#ld_cov_crop <- crop(land_cover, area_borders)
#writeRaster(ld_cov_crop, filename = "data/kaza_croped.tiff")
ld_cov_crop <- raster("data/kaza_croped.tif")
  
plot(ld_cov_crop)
plot(fence_sh, col = "magenta", add = T, lty = 'dotted', lwd = 3)
plot(roads_sh, col = "yellow", add = T , lty = 'dotdash', lwd = 3)

# whole area map 
height <- max(mod_data$lat) - min(mod_data$lat)
width  <- max(mod_data$long) - min(mod_data$long)

area_borders <- c(left    = min(mod_data$long) - 0.025 * width, 
                  right   = max(mod_data$long) + 0.025 * width,
                  bottom  = min(mod_data$lat)  - 0.025 * height, 
                  top     = max(mod_data$lat)  + 0.025 * height)

area_map <- get_stamenmap(area_borders, zoom = 12, maptype = "terrain")
save(area_map, file = "data/stamenmap_z12_terrain.Rdata")
load("data/stamenmap_z12_terrain.Rdata")

# plot all points, color by category
ggmap(area_map, extent = 'panel') + 
  mod_data %>% 
  geom_point(data = ., 
             mapping = aes(x = long, y = lat, color = as.factor(year)), size = .25) +
  geom_path(data = fortify(fence_sh), aes(x = long, y = lat, group = group), colour = "red") +  
  labs(title = "", x = "Longitude", y = "Latitude") + 
  theme(axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(color = "black", size = 10, face = "bold")) +
  theme(legend.position = "top") + scale_colour_viridis_d(option = "inferno")

# plot by individual's id
i=3
subbs_dt <- filter(mod_data, id == unique(mod_data$id)[i])

ggmap(area_map, extent = 'panel') + 
  geom_point(data = subbs_dt, 
             mapping = aes(x = long, y = lat, color = id), size = .25) + 
  #geom_line(data = subbs_dt, aes(x = long, y = lat, color = as.factor(month)), linetype = 3) +
  geom_path(data = fortify(fence_sh), 
            aes(x = long, y = lat, group = group), 
            colour = "red", linetype = "dotdash") +  
  labs(title = "", subtitle = paste0("id: ", unique(raw_data$id)[i]), 
       x = "Longitude", y = "Latitude") + 
  xlim(c(min(subbs_dt$long) - 0.01, max(subbs_dt$long) + 0.01)) + 
  ylim(c(min(subbs_dt$lat) - 0.01, max(subbs_dt$lat) + 0.01)) +
  theme(axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(color = "black", size = 10, face = "bold")) +
  theme(legend.position = "right") + 
  scale_colour_viridis_d(option = "inferno")

## plot by category with levels
#creates a folder and saves the graphs

cath <- "year"
var <- mod_data %>% pull(cath) %>% unique()

dir.create(paste0("figures/", cath))
for(i in 1:length(var)){
  subbs_dt <- filter(mod_data, year == var[i])
  ggmap(area_map, extent = 'panel') + 
    geom_point(data = subbs_dt, 
               mapping = aes(x = long, y = lat, color = id), size = .25) + 
    geom_path(data = fortify(fence_sh), 
              aes(x = long, y = lat, group = group), 
              colour = "red", linetype = "dotdash") +  
    labs(title = paste0(cath,": ", var[i]), 
         x = "Longitude", y = "Latitude", cex = 14) + 
    xlim(c(min(subbs_dt$long) - 0.01, max(subbs_dt$long) + 0.01)) + 
    ylim(c(min(subbs_dt$lat) - 0.01, max(subbs_dt$lat) + 0.01)) +
    theme(axis.title.x = element_text(color = "black", size = 10, face = "bold"),
          axis.title.y = element_text(color = "black", size = 10, face = "bold")) +
    theme(legend.position = "right", legend.text = element_text(size = 12),
          axis.text = element_text(size = 12)) + 
    scale_colour_viridis_d(option = "inferno") +  
    annotation_north_arrow(style = north_arrow_orienteering)
  ggsave(filename = paste0("figures/", cath, "/", paste0(cath,"_", var[i]), ".png"), 
         device = "png",
         width = 50, height = 30, units = "cm", dpi = 300)
}


  