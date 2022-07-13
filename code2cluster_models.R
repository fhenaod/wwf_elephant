library(momentuHMM)
library(sp)
library(raster)
library(tidyverse)

# load data #####
raw_data <- read.csv("data/Namibia_elephant_2010-Mar2021.csv")
hmm_data <- raw_data %>%
  mutate(date.time = as.POSIXct(date.time, format = "%Y-%m-%d %H:%M", tz = "UTC"))

# add Khaudum IDs
khaudum_ids <- 
  c("MET.kh.1", "MET.kh.2", "MET.kh.4", "MET.kh.3", "MET.kh.5", 
    "MET.kh.6", "MET.kh.7", "MET.kh.8", "MET.kh.9", "MET.kh.10", 
    "01108738SKY8D", "01108752SKYC5", "01108753SKY49", "01108735SKY01", 
    "01108739SKY11", "01108741SKY19", "01108744SKYA5", "01108740SKY95", "01108743SKY21", 
    "01108734SKY7D", "01145133SKY80", "01108734SKY7D", "01145131SKY78", "01145141SKYA0", 
    "01226982SKY22", "01177691SKYB4", "01178599SKYF5", "01183261SKY28", "01226977SKY8E", 
    "01401510SKY73", "01426026SKY62", "01401410SKYE187")

hmm_data <- hmm_data %>% 
  mutate(collared = ifelse(test = id.new %in% khaudum_ids, 
                           yes = "khadum", no = "no_khadum"))

library(raster)
library(rgdal)
fence_raster <- raster("data/fence_raster.tif") # fence buffer
road_trunk <- raster("data/trunk_raster.tif") # trunk buffer
ld_cov_crop_rp_reclf <- raster("data/kaza_croped_reclaf.tif") # land cover reclasified
pas_raster <- raster("data/pas_raster.tif") # protected areas raster

#crs(fence_raster)         <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
#crs(road_trunk)           <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
#crs(ld_cov_crop_rp_reclf) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
#crs(pas_raster)           <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"

load_prep_data <- T
if(load_prep_data==F){
# regularize data ####
reg_data = T # regularize data by indiv

data_reg_ls <- list()
if(reg_data == T){
for(i in 1:length(unique(hmm_data$id)) ){
  hmm_data_f <- hmm_data %>% 
    filter(id.new == unique(hmm_data$id)[i])

  ti <- seq(hmm_data_f$date.time[1], hmm_data_f$date.time[length(hmm_data_f$date.time)],
            by = 2*60*60)
  
  iLoc <- as.data.frame(cbind(lon = approx(hmm_data_f$date.time, hmm_data_f$lon, xout = ti)$y, 
                              lat = approx(hmm_data_f$date.time, hmm_data_f$lat, xout = ti)$y))
  
  data_reg_ls[[i]] <- cbind(ID = unique(hmm_data$id)[i], date = ti, iLoc, 
                    sex = unique(hmm_data_f$sex), collared = unique(hmm_data_f$collared))
}
  data2prep <- map_df(data_reg_ls, rbind) #  for serial version
} else {
  data2prep <- hmm_data %>% 
    select(ID = id, lat, lon = long, sex, date, collared)
}

data2prep %>% head()
data2prep %>% distinct(ID)
data2prep %>% dim()

# prepare data ####
paral_prep = F
prep_ss_ls <- c()
if(paral_prep == T){
  # parallel version not working jet !!!!
  library(parallel)
  x <- data_reg_ls
  prep_ss_ls <- 
  mclapply(x, function(x) prepData(x, type = "LL", 
                                coordNames = c("lon", "lat"),
                                covNames = c("sex", "collared"), 
                                spatialCovs = list(cover = ld_cov_crop_rp_reclf, 
                                                   fence = fence_raster, 
                                                   road = road_trunk, 
                                                   pas = pas_raster) 
                                                 #, prec = bc_data_crop$bioclim_rasters.12
  ), mc.cores = 2)
  
  #prep_ss_ls <- 
  #mclapply(function(x) prepData(x, type = "LL", coordNames = c("lon", "lat"),
  #                                               covNames = c("sex", "collared"),
  #                                               spatialCovs = list(cover = ld_cov_crop_rp_reclf, 
  #                                                                  fence = fence_raster, 
  #                                                                  road = road_trunk,
  #                                                                  pas = pas_raster) 
  #                                               #, prec = bc_data_crop$bioclim_rasters.12
  #), X = x, mc.cores = 3)
  
  data_prep_cv2 <- map_df(prep_ss_ls, rbind) # convert prep_data list in one data.frame
} else {
  data_prep_cv2 <- prepData(data2prep, type = "LL", 
                            coordNames = c("lon", "lat"),
                            covNames = c("sex", "collared"),
                            spatialCovs = list(cover = ld_cov_crop_rp_reclf, 
                                               fence = fence_raster, 
                                               road = road_trunk,
                                               pas = pas_raster) 
                            #, prec = bc_data_crop$bioclim_rasters.12
  )
}

if(sum(data_prep_cv2$step == 0, na.rm = T)!=0) print("step_len has zeros")

dir.create("output/")

  data_prep_cv2$hour  <- as.integer(strftime(data_prep_cv2$date, format = "%H", tz = "UTC"))
  data_prep_cv2$month <- as.integer(strftime(data_prep_cv2$date, format = "%m", tz = "UTC"))
  data_prep_cv2$year  <- as.integer(strftime(data_prep_cv2$date, format = "%Y", tz = "UTC"))
  data_prep_cv2$cover <- data_prep_cv2$cover %>% as.integer()
  data_prep_cv2$sex <- data_prep_cv2$sex %>% as.factor()
  data_prep_cv2$collared <- data_prep_cv2$collared %>% as.factor()
  
saveRDS(data_prep_cv2, "output/data_prep2hmm.rds")
write.csv(data_prep_cv2, "output/data_prep2hmm.csv")
 } else {
   dir.create("output/")
  data_prep_cv2 <- readRDS("output/data_prep2hmm.rds")
  data_prep_cv <- data_prep_cv2 #%>% filter(collared == "no_khadum")
}

  print("COVARS:")
  table(data_prep_cv$fence) ; table(data_prep_cv$road); table(data_prep_cv$pas)
  table(data_prep_cv$cover) ; 
  
# detect multiple modes in continuous data ####
  modes_ts <- c(2, 3, 4, 6)
  mean_stp <- list() # step
  for(i in 1:length(modes_ts)){
    mean_stp[[i]] <- multimode::locmodes(data_prep_cv$step, mod0 = modes_ts[i])
  }
  
  mean_ang <- list() # angle
  for(i in 1:length(modes_ts)){
    mean_ang[[i]] <- multimode::locmodes(data_prep_cv$angle, mod0 = modes_ts[i])
  }
  
  par(mfrow = c(2,2))
  png(filename = "output/mean_step_plot.png")
  sapply(mean_stp, function(x) plot(x, addLegend = F))
  dev.off()
  png(filename = "output/mean_angle_plot.png")
  sapply(mean_ang, function(x) plot(x, addLegend = F))
  dev.off()
  par(mfrow = c(1,1))
  
  png(filename = "output/acf_step_.png")
  acf(data_prep_cv$step[!is.na(data_prep_cv$step)], lag.max = 500)
  dev.off()
  png(filename = "output/acf_angle_.png")
  acf(data_prep_cv$step[!is.na(data_prep_cv$angle)], lag.max = 500)
  dev.off()
  
  #step_locs <- sapply(mean_stp, "[", 1, simplify = T)
  #ang_locs <- sapply(mean_ang, "[", 1, simplify = T)

# model set ups ####
  step_dist <- "gamma" # step distribution
  angle_dist <- "vm" # turning angle distribution
  retry_fits <- 0
  outfile <- "output/"

  int_pars <- data_prep_cv %>% dplyr::summarize(mean_step = mean(step, na.rm = T),
                                                median_step = median(step, na.rm = T),
                                                sd_step = sd(step, na.rm = T),
                                                mean_ang = mean(angle, na.rm = T),
                                                sd_ang = sd(angle, na.rm = T))
  
  mu0 <- c(mean_stp[[1]]$locations[seq(from = 1, to = length(mean_ang[[1]]$locations), 2)])
  sigma0 <- c(int_pars$sd_step, int_pars$sd_step*2)
  kappa0 <- c(int_pars$sd_ang, int_pars$sd_ang)

# run simple models ####
run_simp_mods <- F
  if(run_simp_mods==T){
  # Model 1
  #if (model == 1){
    model1 <- try(fitHMM(data_prep_cv, nbState = 2, 
                         dist = list(step = step_dist, angle = angle_dist),
                         Par0 = list(step = c(mu0, sigma0, rep(0.01, 2)), angle = kappa0),
                         formula = ~1, estAngleMean = NULL, retryFits = retry_fits))
    saveRDS(model1, paste0(outfile, "model1.rds"))
    capture.output(print(paste0("Analysis- model 1 done.")), 
                   file = paste0(outfile, "results.log"), append = TRUE)
  #}
  
  mu03s <- c(mean_stp[[2]]$locations[seq(from = 1, to = length(mean_ang[[2]]$locations), 2)]) 
  
  sigma03s <- c(model1$mle$step["sd",1],
                model1$mle$step["sd",2],
                model1$mle$step["sd",2]*2) 
  
  kappa03s <- c(model1$mle$angle["concentration",1],
                model1$mle$angle["concentration",2],
                model1$mle$angle["concentration",2]*10)
  
  # Model 2
  #if (model == 2){
    model2 <- try(fitHMM(data_prep_cv, nbState = 3, 
                         dist = list(step = step_dist, angle = angle_dist),
                         Par0 = list(step = c(mu03s, sigma03s, rep(0.01, 3)), angle = kappa03s), 
                         formula = ~1, estAngleMean = NULL, retryFits = retry_fits))
    saveRDS(model2, paste0(outfile, "model2.rds"))
    capture.output(print(paste0("Analysis- model  - model 2 done.")), 
                   file = paste0(outfile, "results.log"), append = TRUE)
  #}
  
  mu04s <- c(mean_stp[[3]]$locations[seq(from = 1, to = length(mean_ang[[3]]$locations), 2)]) 
  
  sigma04s <- c(model1$mle$step["sd",1],
                model1$mle$step["sd",2],
                model1$mle$step["sd",2]*2,
                model1$mle$step["sd",2]*.5) 
  
  kappa04s <- c(model1$mle$angle["concentration",1],
                model1$mle$angle["concentration",2],
                model1$mle$angle["concentration",2]*10,
                model1$mle$angle["concentration",2]*.5)
  
  # Model 3
  #if (model == 3){
    model3 <- try(fitHMM(data_prep_cv, nbState = 4, 
                         dist = list(step = step_dist, angle = angle_dist),
                         Par0 = list(step = c(mu04s, sigma04s, rep(0.01, 4)), angle = kappa04s),
                         formula = ~1, estAngleMean = NULL, retryFits = retry_fits))
    saveRDS(model3, paste0(outfile, "model3.rds"))
    capture.output(print(paste0("Analysis- model  - model 3 done.")), 
                   file = paste0(outfile, "results.log"), append = TRUE)
  #}
  
  mu06s <- c(mean_stp[[4]]$locations[seq(from = 1, to = length(mean_ang[[4]]$locations), 2)]) 
  
  sigma06s <- c(model1$mle$step["sd",1],
                model1$mle$step["sd",2],
                model1$mle$step["sd",2]*2,
                model1$mle$step["sd",2]*.5,
                model1$mle$step["sd",1]*2,
                model1$mle$step["sd",1]*.5) 
  
  kappa06s <- c(model1$mle$angle["concentration",1],
                model1$mle$angle["concentration",2],
                model1$mle$angle["concentration",2]*10,
                model1$mle$angle["concentration",2]*.5,
                model1$mle$angle["concentration",1]*10,
                model1$mle$angle["concentration",1]*.5)
  
  # Model 4
  #if (model == 4){
    model4 <- try(fitHMM(data_prep_cv, nbState = 6, 
                         dist = list(step = step_dist, angle = angle_dist),
                         Par0 = list(step = c(mu06s, sigma06s, rep(0.01, 6)), angle = kappa06s),
                         formula = ~1, estAngleMean = NULL, retryFits = retry_fits))
    saveRDS(model4, paste0(outfile, "model4.rds"))
    capture.output(print(paste0("Analysis- model  - model 4 done.")), 
                   file = paste0(outfile, "results.log"), append = TRUE)
  #}  
  } else {
   
    model1 <- readRDS(paste0(outfile, "model1.rds"))
    model2 <- readRDS(paste0(outfile, "model2.rds"))
    model3 <- readRDS(paste0(outfile, "model3.rds"))
   # model4 <- readRDS(paste0(outfile, "model4.rds"))
    
 }
# function to further models ####      
fit_HMM_mods <- function(model, data_prep_cv, outfile, ...){
  # spatial models
  # Model 5
  if (model == 5){
    par0b <- getPar0(model = model1, formula = ~cover)
    model5 <- try(fitHMM(data_prep_cv, nbState = 2, 
                         dist = list(step = step_dist, angle = angle_dist),
                         Par0 = par0b$Par, beta0 = par0b$beta, 
                         formula = ~cover))
    saveRDS(model5, paste0(outfile, "model5.rds"))
    capture.output(print(paste0("Analysis- model  - model 5 done.")), 
                   file = paste0(outfile, "results.log"), append = TRUE)
  }
  
  # Model 6
  if (model == 6){
    par0b <- getPar0(model = model2, formula = ~cover)
    model6 <- try(fitHMM(data_prep_cv, nbState = 3, 
                         dist = list(step = step_dist, angle = angle_dist),
                         Par0 = par0b$Par, beta0 = par0b$beta, 
                         formula = ~cover))
    saveRDS(model6, paste0(outfile, "model6.rds"))
    capture.output(print(paste0("Analysis- model  - model 6 done.")), 
                   file = paste0(outfile, "results.log"), append = TRUE)
  }
  
  # Model 7
  if (model == 7){
    par0b <- getPar0(model = model3, formula = ~cover)
    model7 <- try(fitHMM(data_prep_cv, nbState = 4, 
                         dist = list(step = step_dist, angle = angle_dist),
                         Par0 = par0b$Par, beta0 = par0b$beta, 
                         formula = ~cover))
    saveRDS(model7, paste0(outfile, "model7.rds"))
    capture.output(print(paste0("Analysis- model  - model 7 done.")), 
                   file = paste0(outfile, "results.log"), append = TRUE)
  }
  
  # Model 8
  if (model == 8){
    par0b <- getPar0(model = model1, formula = ~fence)
    model8 <- try(fitHMM(data_prep_cv, nbState = 2, 
                         dist = list(step = step_dist, angle = angle_dist),
                         Par0 = par0b$Par, beta0 = par0b$beta, 
                         formula = ~fence))
    saveRDS(model8, paste0(outfile, "model8.rds"))
    capture.output(print(paste0("Analysis- model  - model 8 done.")), 
                   file = paste0(outfile, "results.log"), append = TRUE)
  }
  
  # Model 9
  if (model == 9){
    par0b <- getPar0(model = model2, formula = ~fence)
    model9 <- try(fitHMM(data_prep_cv, nbState = 3, 
                         dist = list(step = step_dist, angle = angle_dist),
                         Par0 = par0b$Par, beta0 = par0b$beta, 
                         formula = ~fence))
    saveRDS(model9, paste0(outfile, "model9.rds"))
    capture.output(print(paste0("Analysis- model  - model 9 done.")), 
                   file = paste0(outfile, "results.log"), append = TRUE)
  }
  
  # Model 10
  if (model == 10){
    par0b <- getPar0(model = model3, formula = ~fence)
    model10 <- try(fitHMM(data_prep_cv, nbState = 4, 
                          dist = list(step = step_dist, angle = angle_dist),
                          Par0 = par0b$Par, beta0 = par0b$beta, 
                          formula = ~fence))
    saveRDS(model10, paste0(outfile, "model10.rds"))
    capture.output(print(paste0("Analysis- model  - model 10 done.")), 
                   file = paste0(outfile, "results.log"), append = TRUE)
  }
  
  # Model 11
  if (model == 11){
    par0b <- getPar0(model = model1, formula = ~road)
    model11 <- try(fitHMM(data_prep_cv, nbState = 2, 
                          dist = list(step = step_dist, angle = angle_dist),
                          Par0 = par0b$Par, beta0 = par0b$beta, 
                          formula = ~road))
    saveRDS(model11, paste0(outfile, "model11.rds"))
    capture.output(print(paste0("Analysis- model  - model 11 done.")), 
                   file = paste0(outfile, "results.log"), append = TRUE)
  }
  
  # Model 12
  if (model == 12){
    par0b <- getPar0(model = model2, formula = ~road)
    model12 <- try(fitHMM(data_prep_cv, nbState = 3, 
                          dist = list(step = step_dist, angle = angle_dist),
                          Par0 = par0b$Par, beta0 = par0b$beta, 
                          formula = ~road))
    saveRDS(model12, paste0(outfile, "model12.rds"))
    capture.output(print(paste0("Analysis- model  - model 12 done.")), 
                   file = paste0(outfile, "results.log"), append = TRUE)
  }
  
  # Model 13
  if (model == 13){
    par0b <- getPar0(model = model3, formula = ~road)
    model13 <- try(fitHMM(data_prep_cv, nbState = 4, 
                          dist = list(step = step_dist, angle = angle_dist),
                          Par0 = par0b$Par, beta0 = par0b$beta, 
                          formula = ~road))
    saveRDS(model13, paste0(outfile, "model13.rds"))
    capture.output(print(paste0("Analysis- model  - model 13 done.")), 
                   file = paste0(outfile, "results.log"), append = TRUE)
  }
  
  # Model 14
  if (model == 14){
    par0b <- getPar0(model = model1, formula = ~cover+fence)
    model14 <- try(fitHMM(data_prep_cv, nbState = 2, 
                          dist = list(step = step_dist, angle = angle_dist),
                          Par0 = par0b$Par, beta0 = par0b$beta, 
                          formula = ~cover+fence))
    saveRDS(model14, paste0(outfile, "model14.rds"))
    capture.output(print(paste0("Analysis- model  - model 14 done.")), 
                   file = paste0(outfile, "results.log"), append = TRUE)
  }
  
  # Model 15
  if (model == 15){
    par0b <- getPar0(model = model2, formula = ~cover+fence)
    model15 <- try(fitHMM(data_prep_cv, nbState = 3, 
                          dist = list(step = step_dist, angle = angle_dist),
                          Par0 = par0b$Par, beta0 = par0b$beta, 
                          formula = ~cover+fence))
    saveRDS(model15, paste0(outfile, "model15.rds"))
    capture.output(print(paste0("Analysis- model  - model 15 done.")), 
                   file = paste0(outfile, "results.log"), append = TRUE)
  }
  
  # Model 16
  if (model == 16){
    par0b <- getPar0(model = model3, formula = ~cover+fence)
    model16 <- try(fitHMM(data_prep_cv, nbState = 4, 
                          dist = list(step = step_dist, angle = angle_dist),
                          Par0 = par0b$Par, beta0 = par0b$beta, 
                          formula = ~cover+fence))
    saveRDS(model16, paste0(outfile, "model16.rds"))
    capture.output(print(paste0("Analysis- model  - model 16 done.")), 
                   file = paste0(outfile, "results.log"), append = TRUE)
  }
  
  # Model 17
  if (model == 17){
    par0b <- getPar0(model = model1, formula = ~cover+fence+road)
    model17 <- try(fitHMM(data_prep_cv, nbState = 2, 
                          dist = list(step = step_dist, angle = angle_dist),
                          Par0 = par0b$Par, beta0 = par0b$beta, 
                          formula = ~cover+fence+road))
    saveRDS(model17, paste0(outfile, "model17.rds"))
    capture.output(print(paste0("Analysis- model  - model 17 done.")), 
                   file = paste0(outfile, "results.log"), append = TRUE)
  }
  
  # Model 18
  if (model == 18){
    par0b <- getPar0(model = model2, formula = ~cover+fence+road)
    model18 <- try(fitHMM(data_prep_cv, nbState = 3, 
                          dist = list(step = step_dist, angle = angle_dist),
                          Par0 = par0b$Par, beta0 = par0b$beta, 
                          formula = ~cover+fence+road))
    saveRDS(model18, paste0(outfile, "model18.rds"))
    capture.output(print(paste0("Analysis- model  - model 18 done.")), 
                   file = paste0(outfile, "results.log"), append = TRUE)
  }
  
  # Model 19
  if (model == 19){
    par0b <- getPar0(model = model3, formula = ~cover+fence+road)
    model19 <- try(fitHMM(data_prep_cv, nbState = 4, 
                          dist = list(step = step_dist, angle = angle_dist),
                          Par0 = par0b$Par, beta0 = par0b$beta, 
                          formula = ~cover+fence+road))
    saveRDS(model19, paste0(outfile, "model19.rds"))
    capture.output(print(paste0("Analysis- model  - model 19 done.")), 
                   file = paste0(outfile, "results.log"), append = TRUE)
  }
  
  # Model 20
  if (model == 20){
    par0b <- getPar0(model = model1, formula = ~fence+road)
    model20 <- try(fitHMM(data_prep_cv, nbState = 2, 
                          dist = list(step = step_dist, angle = angle_dist),
                          Par0 = par0b$Par, beta0 = par0b$beta, 
                          formula = ~fence+road))
    saveRDS(model20, paste0(outfile, "model20.rds"))
    capture.output(print(paste0("Analysis- model  - model 20 done.")), 
                   file = paste0(outfile, "results.log"), append = TRUE)
  }
  
  # Model 21
  if (model == 21){
    par0b <- getPar0(model = model2, formula = ~fence+road)
    model21 <- try(fitHMM(data_prep_cv, nbState = 3, 
                          dist = list(step = step_dist, angle = angle_dist),
                          Par0 = par0b$Par, beta0 = par0b$beta, 
                          formula = ~fence+road))
    saveRDS(model21, paste0(outfile, "model21.rds"))
    capture.output(print(paste0("Analysis- model  - model 21 done.")), 
                   file = paste0(outfile, "results.log"), append = TRUE)
  }
  
  # temporal models
  # Model 22
  if (model == 22){
    par0b <- getPar0(model = model1, formula = ~cosinor(hour, period = 24))
    model22 <- try(fitHMM(data_prep_cv, nbState = 2, 
                          dist = list(step = step_dist, angle = angle_dist),
                          Par0 = par0b$Par, beta0 = par0b$beta, 
                          formula = ~cosinor(hour, period = 24)))
    saveRDS(model22, paste0(outfile, "model22.rds"))
    capture.output(print(paste0("Analysis- model  - model 22 done.")), 
                   file = paste0(outfile, "results.log"), append = TRUE)
  }
  
  # Model 23
  if (model == 23){
    par0b <- getPar0(model = model1, formula = ~cosinor(month, period = 6))
    model23 <- try(fitHMM(data_prep_cv, nbState = 2, 
                          dist = list(step = step_dist, angle = angle_dist),
                          Par0 = par0b$Par, beta0 = par0b$beta, 
                          formula = ~cosinor(month, period = 6)))
    saveRDS(model23, paste0(outfile, "model23.rds"))
    capture.output(print(paste0("Analysis- model  - model 23 done.")), 
                   file = paste0(outfile, "results.log"), append = TRUE)
  }
  
  # Model 24
  if (model == 24){
    par0b <- getPar0(model = model2, formula = ~cosinor(month, period = 6))
    model24 <- try(fitHMM(data_prep_cv, nbState = 3, 
                          dist = list(step = step_dist, angle = angle_dist),
                          Par0 = par0b$Par, beta0 = par0b$beta, 
                          formula = ~cosinor(month, period = 6)))
    saveRDS(model24, paste0(outfile, "model24.rds"))
    capture.output(print(paste0("Analysis- model  - model 24 done.")), 
                   file = paste0(outfile, "results.log"), append = TRUE)
  }
  # demographic models
  # Model 25
  if (model == 25){
    par0b <- getPar0(model = model1)
    model25 <- try(fitHMM(data_prep_cv, nbState = 2, 
                          dist = list(step = step_dist, angle = angle_dist),
                          Par0 = par0b$Par, 
                          formula = ~ sex + 0))
    saveRDS(model25, paste0(outfile, "model25.rds"))
    capture.output(print(paste0("Analysis- model  - model 25 done.")), 
                   file = paste0(outfile, "results.log"), append = TRUE)
  }
  # Model 26
    if (model == 26){
      par0b <- getPar0(model = model2)
      model26 <- try(fitHMM(data_prep_cv, nbState = 3, 
                            dist = list(step = step_dist, angle = angle_dist),
                            Par0 = par0b$Par, 
                            formula = ~ sex + 0))
      saveRDS(model26, paste0(outfile, "model26.rds"))
      capture.output(print(paste0("Analysis- model  - model 26 done.")), 
                     file = paste0(outfile, "results.log"), append = TRUE)
    }
  # Model 27
    if (model == 27){
      par0b <- getPar0(model = model3)
      model27 <- try(fitHMM(data_prep_cv, nbState = 4, 
                            dist = list(step = step_dist, angle = angle_dist),
                            Par0 = par0b$Par, 
                            formula = ~ sex + 0))
      saveRDS(model27, paste0(outfile, "model27.rds"))
      capture.output(print(paste0("Analysis- model  - model 27 done.")), 
                     file = paste0(outfile, "results.log"), append = TRUE)
    }
    # Model 28
    if (model == 28){
      par0b <- getPar0(model = model1)
      model28 <- try(fitHMM(data_prep_cv, nbState = 2, 
                            dist = list(step = step_dist, angle = angle_dist),
                            Par0 = par0b$Par, 
                            formula = ~ collared + 0))
      saveRDS(model28, paste0(outfile, "model28.rds"))
      capture.output(print(paste0("Analysis- model  - model 28 done.")), 
                     file = paste0(outfile, "results.log"), append = TRUE)
    }  
    # Model 29
    if (model == 29){
      par0b <- getPar0(model = model2)
      model29 <- try(fitHMM(data_prep_cv, nbState = 3, 
                            dist = list(step = step_dist, angle = angle_dist),
                            Par0 = par0b$Par, 
                            formula = ~ collared + 0))
      saveRDS(model29, paste0(outfile, "model29.rds"))
      capture.output(print(paste0("Analysis- model  - model 29 done.")), 
                     file = paste0(outfile, "results.log"), append = TRUE)
    }  
    # Model 30
    if (model == 30){
      par0b <- getPar0(model = model3)
      model30 <- try(fitHMM(data_prep_cv, nbState = 4, 
                            dist = list(step = step_dist, angle = angle_dist),
                            Par0 = par0b$Par, 
                            formula = ~ collared + 0))
      saveRDS(model30, paste0(outfile, "model30.rds"))
      capture.output(print(paste0("Analysis- model  - model 30 done.")), 
                     file = paste0(outfile, "results.log"), append = TRUE)
    }
  # Model 31
  if (model == 31){
    par0b <- getPar0(model = model1)
    model31 <- try(fitHMM(data_prep_cv, nbState = 2, 
                          dist = list(step = step_dist, angle = angle_dist),
                          Par0 = par0b$Par, 
                          formulaPi = ~ sex + 0, mixtures = 2))
    saveRDS(model31, paste0(outfile, "model31.rds"))
    capture.output(print(paste0("Analysis- model  - model 31 done.")), 
                   file = paste0(outfile, "results.log"), append = TRUE)
  } 
  # Model 32
  if (model == 32){
    par0b <- getPar0(model = model1)
    model32 <- try(fitHMM(data_prep_cv, nbState = 2, 
                          dist = list(step = step_dist, angle = angle_dist),
                          Par0 = par0b$Par, 
                          formulaPi = ~ collared + 0, mixtures = 2))
    saveRDS(model32, paste0(outfile, "model32.rds"))
    capture.output(print(paste0("Analysis- model  - model 32 done.")), 
                   file = paste0(outfile, "results.log"), append = TRUE)
  } 
  print("Analysis completed")
}    

# running models ####
run_type <- c("parallel", "serial")
n_models <- 32
if(run_type[1] == "serial"){
  # serial or individual model fits
  for(i in 5:n_models){
    fit_HMM_mods(i, data_prep_cv, "output/"
                 , retryFits = retry_fits
                 )
    print(paste0("serial_run: ", i))
  }
} else {
  # parallel model fits
  library(parallel)
  print("parallel run")
  mclapply(X = c(5:n_models), 
           function(X) fit_HMM_mods(X, data_prep_cv, "output/", retryFits = retry_fits), mc.cores = n_models)
}

# load results an compute model fits ####
path <- paste0("output/")
res_fs <- dir(path = path, pattern = "model") %>% jamba::mixedSort()
model <- sapply(strsplit(res_fs, split = "\\."), "[", 1) #%>%strsplit(split = "\\.") %>% sapply("[", 1)
 
models <- list()
n_covars <- c(); n_iteras <- c()
form_mod <- c(); n_states <- c()
 for(f in 1:length(res_fs)){
   models[[f]] <- readRDS(paste0(path, res_fs[f]))
   if(class(models[[f]][1]) == "try-error" | class(models[[f]])[1] == "try-error"){
     n_covars[f] <- -999
     n_iteras[f] <- -999
     form_mod[f] <- -999
     n_states[f] <- -999
   } else {
     n_covars[f] <- dim(models[[f]]$rawCovs)[2]
     n_iteras[f] <- models[[f]]$mod$iterations
     form_mod[f] <- as.character(models[[f]]$conditions$formula)[2]
     n_states[f] <- models[[f]]$stateNames %>% length()
    }
 }

if(sum(n_covars==-999, na.rm = T)>0){
  m2f <- which(n_covars==-999)
  data.frame(model = model[-m2f], n_covars = n_covars[-m2f], n_iter = n_iteras[-m2f], 
             formula = form_mod[-m2f], n_states = n_states[-m2f], 
             AIC = sapply(models[-m2f], AIC), 
             geiger::aicw(sapply(models[-m2f], AIC))) %>% select(-fit) %>% 
    mutate(across(where(is.numeric), round, 4)) %>% arrange(desc(w)) %>% 
    write.csv("output/fit_mods_res.csv")
  } else {
  data.frame(model = model, n_covars = n_covars, n_iter = n_iteras, 
             formula = form_mod, n_states = n_states, 
             AIC = sapply(models, 
                          geiger::aicw(sapply(models, AIC)))) %>% select(-fit) %>% 
  mutate(across(where(is.numeric), round, 4)) %>% arrange(desc(w)) %>%  
  write.csv("output/fit_mods_res.csv")
}
