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

# add Khaudum IDs
khaudum_ids <- 
c("MET.kh.1", "MET.kh.2", "MET.kh.4", "MET.kh.3", "MET.kh.5", 
  "MET.kh.6", "MET.kh.7", "MET.kh.8", "MET.kh.9", "MET.kh.10", 
  "01108738SKY8D", "01108752SKYC5", "01108753SKY49", "01108735SKY01", 
  "01108739SKY11", "01108741SKY19", "01108744SKYA5", "01108740SKY95", "01108743SKY21", 
  "01108734SKY7D", "01145133SKY80", "01108734SKY7D", "01145131SKY78", "01145141SKYA0", 
  "01226982SKY22", "01177691SKYB4", "01178599SKYF5", "01183261SKY28", "01226977SKY8E", 
  "01401510SKY73", "01426026SKY62", "01401410SKYE187")

mod_data <- mod_data %>% #sample_n(100) %>% 
  mutate(collared = ifelse(test = id.new %in% khaudum_ids, yes = "khadum", no = "no_khadum"))

# add remote sensing data ####

# regularize data ####
# add spatial vars ####
library(ggmap)
library(ggspatial)
library(raster)
library(rgdal)
fence_raster <- raster("data/fence_raster.tif") # fence buffer
road_trunk <- raster("data/trunk_raster.tif") # trunk buffer
ld_cov_crop_rp_reclf <- raster("data/kaza_croped_reclaf.tif") # land cover reclasified
pas_raster <- raster("data/pas_raster.tif") # protected areas raster
