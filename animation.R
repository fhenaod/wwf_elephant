library(tidyverse)

# animation with path ####
library(gganimate)
library(gifski)

subbs_dt <- filter(mod_data, id == unique(mod_data$id)[i])
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

var2iter <- mod_data %>% pull(year) %>% unique()
  
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
