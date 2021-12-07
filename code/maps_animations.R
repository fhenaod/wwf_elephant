library(tidyverse)

raw_data <- read.csv("data/Namibia_elephant_2010-Mar2021.csv")
raw_data %>% head()

mod_data <- raw_data %>%
  mutate(date.time = as.POSIXct(date.time, format = "%Y-%m-%d %H:%M", tz = "UTC"))

# mapswith gps points ####
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





# animation with path ####
library(gganimate)
library(gifski)

subbs_dt <- filter(mod_data, id == unique(mod_data$id)[i]) # by individual id
p <- ggmap(area_map, extent = 'panel') +
  geom_point(data = subbs_dt, aes(x = long, y = lat, fill = as.factor(month)), 
             alpha = .7, shape = 21) +
  geom_path(data = subbs_dt, aes(x = long, y = lat, color = as.factor(month)),
            alpha = .7, size = .25) +
  labs(title = "", subtitle = paste0("id: ", unique(raw_data$id)[i]), 
       x = "Longitude", y = "Latitude") + 
  xlim(c(min(subbs_dt$long) - 0.01, max(subbs_dt$long) + 0.01)) + 
  ylim(c(min(subbs_dt$lat) - 0.01, max(subbs_dt$lat) + 0.01)) +
  #scale_fill_viridis_c(option = "inferno")  +
  #scale_color_viridis_c(option = "inferno") +
  scale_size_continuous(range = c(0.1,10))  +
  theme_dark() +
  theme(panel.grid = element_blank(),
        legend.position = "none")

anim <- p + transition_reveal(along = date.time) + 
  ease_aes("linear") + ggtitle("Date: {frame_along}")

animate(anim, nframes = subbs_dt$date %>% unique() %>% length(), fps = 4,
        end_pause = 10, rewind = F)  
dev.off()

# animation with tail ####
#library(OpenStreetMap)
library(moveVis)
library(move)

var2iter <- mod_data %>% pull(year) %>% unique() # by variable
  
for(i in 1:length(var2iter)){
  subbs_dt <- filter(mod_data, year == var2iter[i])
  subbs_dt %>% dim()
  #subbs_dt <- subbs_dt %>% filter(id.new %in% c("MET.kh.1", "MET.kh.10"))
  
  var2plt <- "id.new"
  mov_dat <- subbs_dt %>% 
    df2move(proj = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", 
            x = "long", y = "lat", time = "date.time", track_id = var2plt)
  mov_dat %>% head()
  ninds <- subbs_dt %>% pull(var2plt) %>% unique() %>% length()
  
  mm <- align_move(mov_dat, res = 6, unit = "hours")
  mm %>% dim()
  
  pal <- colorRampPalette(c("red", "blue", "orange", "darkgreen", 
                            "brown", "yellow", "purple", "pink")) 
  fram <- frames_spatial(mm, 
                         path_colours = pal(ninds), path_size = .7, 
                         tail_size = .7, 
                         map_service = "osm", map_res = 1, 
                         path_legend_title = "Individual" 
                         #, map_type = "terrain_bg" 
                         #, alpha = 0.2, trace_show = T
  ) %>%
    add_labels(x = "Longitude", y = "Latitude") %>% # add some customizations, such as axis labels
    add_northarrow() %>%
    add_scalebar() %>%
    add_timestamps(mm, type = "label") %>%
    add_progress()
  
  fram[[length(fram)-1]]
  
  animate_frames(fram, out_file = paste0("move_",var2iter[i],".gif"), fps = 5, 
                 overwrite = T, end_pause = 3)  
}
