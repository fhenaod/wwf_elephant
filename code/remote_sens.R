library(tidyverse)

# BioClim ####
bc_data <- getData("worldclim", var = "bio", res = 0.5, 
                   lon = 22, lat = -18.5)
bc_data_crop <- crop(bc_data, extent(ld_cov_crop))

writeRaster(bc_data_crop, 
            filename = "data/bioclim_rasters.tif", format = "GTiff",
            options = c("INTERLEAVE=BAND","COMPRESS=LZW"), overwrite = T)

# MODIS rasters ####
library(MODIStsp)
MODIStsp_get_prodlayers("M*D13A2")
spat_file <- system.file("data/shapes/merged_fences_old and 2021 AA5FEB2021.shp", 
                         package = "MODIStsp")

raw_data <- read.csv("data/Namibia_elephant_2010-Mar2021.csv")
hmm_data <- raw_data %>%
  mutate(date.time = as.POSIXct(date.time, format = "%Y-%m-%d %H:%M", tz = "UTC"))

modis_dat_tsp <- 
  MODIStsp(
    gui = FALSE,
    out_folder = "data/",
    selprod = "Vegetation_Indexes_16Days_1Km (M*D13A2)",
    bandsel = c("EVI"),
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

#MODIStsp(
#  gui = FALSE,
#  out_folder = "data/",
#  selprod = "Vegetation_Indexes_16Days_1Km (M*D13A2)",
#  bandsel = c("EVI"),
#  quality_bandsel = "QA_usef",
#  indexes_bandsel = "SR",
#  spatmeth = "file",
#  spafile = "data/shapes/merged_fences_old and 2021 AA5FEB2021.shp",
#  user = "mstp_test" ,
#  password = "MSTP_test_01",
#  start_date = gsub("-", "\\.", min(hmm_data$date)),
#  end_date = gsub("-", "\\.", max(hmm_data$date)),
#  verbose = FALSE,
#  parallel = FALSE
#)

ndvi_15 <- raster("data/merged_fences_old and 2021 AA5FEB2021/VI_16Days_1Km_v6/NDVI/MOD13A2_NDVI_2020_161.tif")
ndvi_15_rp <- projectRaster(ndvi_15, crs = crs(bc_data_crop))
plot(ndvi_15_rp)

# MODISTools
library(MODISTools)
products <- mt_products() # product names
bands <- mt_bands(product = "MOD13Q1") # NDVI product bands

# product spatial-temporal resolution
dates <- mt_dates(product = "MOD13Q1", 
                  lat = extent(fence_sh)[3], lon = extent(fence_sh)[1])
dates$calendar_date[c(1,length(dates$calendar_date))]

namibia_evi  <- mt_subset(product = "MOD13Q1",
                          lat = mean(raw_data$lat),
                          lon =  mean(raw_data$long),
                          band = "250m_16_days_EVI",
                          start = min(hmm_data$date),
                          end = max(hmm_data$date),
                          km_lr = 20,
                          km_ab = 20,
                          site_name = "namibia",
                          internal = TRUE,
                          progress = T)

saveRDS(namibia_evi, "data/namibia_250m_16_evi.rds")
namibia_evi <- readRDS("data/namibia_250m_16_evi.rds")

# viz plots evi per day
namibia_evi %>%
  group_by(calendar_date) %>% 
  summarize(doy = as.numeric(format(as.Date(calendar_date)[1],"%j")), 
            evi_mean = median(value * as.double(scale))) %>% 
  ggplot(aes(x = doy, y = evi_mean)) +
  geom_point() +
  geom_smooth(span = 0.3, method = "loess") +
  labs(x = "day of year (DOY)",
       y = "EVI") +
  theme_minimal() 

# viz plots evi per year
namibia_evi %>% group_by(calendar_date) %>% 
  mutate(doy = as.numeric(format(as.Date(calendar_date)[1],"%Y"))) %>% 
  group_by(doy) %>% #summarize(evi_mean = mean(value*as.double(scale))) %>% 
  ggplot(aes(x = doy, y = value)) +
  geom_point() +
  geom_smooth(span = 0.3, method = "loess") +
  labs(x = "Year",
       y = "Mean EVI") +
  theme_minimal() 

# convert to raster
df <- namibia_evi
dates <- unique(df$calendar_date)
df$scale[df$scale == "Not Available"] <- 1
rr <- do.call("stack", lapply(dates, function(date) {
  m <- matrix(as.numeric(df$value[df$calendar_date == date]) * 
                as.numeric(df$scale[df$calendar_date == date]), df$nrows[1], 
              df$ncols[1], byrow = TRUE)
  return(raster::raster(m))
}))
bb <- MODISTools::mt_bbox(xllcorner = df$xllcorner[1], yllcorner = df$yllcorner[1], 
                          cellsize = df$cellsize[1], 
                          nrows = df$nrows[1], ncols = df$ncols[1], 
                          transform = T)
r <- rr
bb <- st_sf(bb) 
raster::extent(r) <- raster::extent(bb)
raster::projection(r) <- raster::projection(bb)
names(r) <- as.character(dates)
r <- raster::projectRaster(r, crs = "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
plot(r$X2010.10.16)
fence_sh %>% plot(add = T, col = "red", lty = 2)
roads_sh %>% plot(add = T, col = "lightgreen", lty = 3)

# chirps ####
library(chirps)

mod_data %>% head()
chirps_dat_sf <- get_chirps(sample_n(mod_data, 3) %>% select(lon = long, lat), 
                         c(min(mod_data$date), max(mod_data$date))
                         , as.sf = T
)

mo_tb <- sample_n(mod_data, 3) 
chirps_dat <- get_chirps(mo_tb %>% select(lon = long, lat),
                         c(min(mo_tb$date), max(mo_tb$date)), as.sf = F )
chirps_dat %>% data.frame()
p_ind <- precip_indices(chirps_dat, timeseries = T, intervals = 15)
p_ind %>% data.frame() %>% spread(index, value)

mm <- c("R95p","R99p","Rx5day","Rtotal")
days <- c("MLDS","MLWS","R10mm","R20mm")

p_ind %>% 
  filter(index %in% mm) %>% 
  group_by(index) %>% 
  mutate(ab = mean(value)) %>% 
  ggplot() +
  geom_line(aes(x = date, y = value, group = id)) + 
  geom_hline(aes(yintercept = ab), colour = "red", lwd = 0.7) +
  geom_smooth(aes(x = date, y = value), method = "loess") +
  facet_wrap(. ~ index) +
  labs(x = "Year", y = "Index (mm)") +
  theme_classic()

p_ind %>% 
  filter(index %in% days) %>% 
  group_by(index) %>% 
  mutate(ab = mean(value)) %>% 
  ggplot() +
  geom_line(aes(x = date, y = value, group = id)) + 
  geom_hline(aes(yintercept = ab), colour = "red", lwd = 0.7) + 
  geom_smooth(aes(x = date, y = value), method = "loess") +
  facet_wrap(. ~ index) +
  labs(x = "Year", y = "Index (days)") +
  theme_classic()

#
library(heavyRain)
url <- "ftp://ftp.chc.ucsb.edu/pub/org/chc/products/CHIRPS-2.0"
destfile <- "New Chirps Data.R"

gzs <- getCHIRPS(region = "africa", format = "tifs", tres = "monthly"
                 , begin = as.Date("1982-01-01"), end = as.Date("1983-12-31")
                 , sres = .25
                 , dsn = file.path(getwd(), "data")
                 , download.file(url, destfile)
)

data("tapajos", package = "chirps")
tapajos %>% head()
tapajos %>% plot()

# sample three points within the Tapajos area
set.seed(1234)
tp_point <- st_sample(tapajos, 1)

# coerce as sf points
tp_point <- st_as_sf(tp_point)

dt <- get_chirps(tp_point, dates = c("2013-01-01","2013-01-31"))
dt
p_ind <- precip_indices(dt, timeseries = TRUE, intervals = 15)
p_ind %>% data.frame() %>% spread(index, value)

p_ind %>% 
  filter(index %in% mm) %>% 
  group_by(index) %>% 
  mutate(ab = mean(value)) %>% 
  ggplot() +
  geom_line(aes(x = date, y = value, group = id)) + 
  geom_hline(aes(yintercept = ab), colour = "red", lwd = 0.7) +
  geom_smooth(aes(x = date, y = value), method = "loess") +
  facet_wrap(. ~ index) +
  labs(x = "Year", y = "Index (mm)") +
  theme_classic()


lonlat <- data.frame(lon = c(-55.0281,-54.9857),
                     lat = c(-2.8094, -2.8756))

dates <- c("2017-12-15", "2017-12-31")

dt <- get_chirps(lonlat, dates)

dt %>% data.frame()
p_ind <- precip_indices(dt, timeseries = TRUE, intervals = 15)
p_ind %>% data.frame() %>% spread(index, value)



