library(momentuHMM)
library(sp)
library(raster)
library(tidyverse)

raw_data <- read.csv("data/Namibia_elephant_2010-Mar2021.csv")
hmm_data <- raw_data %>%
  mutate(date.time = as.POSIXct(date.time, format = "%Y-%m-%d %H:%M", tz = "UTC"))

head(hmm_data)

fence_raster <- raster("data/fence_raster.tif")
road_raster <- raster("data/road_raster.tif")
ld_cov_crop <- raster("data/kaza_croped.tif")
bc_data_crop <- brick("data/bioclim_rasters.tif") 

for(i in 1:length(unique(mod_data$id))){
  hmm_data_f <- hmm_data %>% filter(id.new == unique(mod_data$id)[i])
  
  ti <- seq(hmm_data_f$date.time[1], hmm_data_f$date.time[length(hmm_data_f$date.time)],
            by = 2*60*60)
  
  iLoc <- as.data.frame(cbind(lon = approx(hmm_data_f$date.time, hmm_data_f$lon, xout = ti)$y, 
                              lat = approx(hmm_data_f$date.time, hmm_data_f$lat, xout = ti)$y))
  data_reg <- cbind(date = ti, iLoc)
  
  data_prep <- prepData(data_reg, type = "LL", coordNames = c("lon", "lat"))
  
  if(sum(data_prep$step == 0, na.rm = T)!=0) stop("ERROR !!!")
  
  data_prep_cv2 <- prepData(data_reg, type = "LL", coordNames = c("lon", "lat"), 
                           spatialCovs = list(cover = ld_cov_crop, 
                                              fence = fence_raster, 
                                              road = road_raster), 
                                              prec = bc_data_crop$bioclim_rasters.12)
  
  data_prep_cv2 %>% select(-c(ID:y)) %>% cor(method = "spearman") %>% print()
  
  data_prep_cv <- data_prep_cv2 %>% sample_n(150)
  
  print("COVARS:")
  table(data_prep_cv$cover) ; table(data_prep_cv$fence) ; table(data_prep_cv$road)

  int_pars <- data_prep %>% dplyr::summarize(mean_step = mean(step, na.rm = T),
                                             med_step = median(step, na.rm = T),
                                             sd_step = sd(step, na.rm = T),
                                             mean_ang = mean(angle, na.rm = T),
                                             sd_ang = sd(angle, na.rm = T))
  
  #getParDM to get initial pars
  
  mu0 <- c(int_pars$mean_step, int_pars$med_step)
  sigma0 <- c(int_pars$sd_step, int_pars$sd_step*2)
  kappa0 <- c(int_pars$sd_ang, int_pars$sd_ang*10)
  
  hmm_2st <- fitHMM(data_prep_cv, nbState = 2, 
                        dist = list(step = "gamma", angle = "vm"),
                        Par0 = list(step = c(mu0, sigma0), angle = kappa0),
                        formula = ~1, estAngleMean = NULL)
  
  mu03s <- c(hmm_2st$mle$step["mean",1], 
             hmm_2st$mle$step["mean",2], 
             hmm_2st$mle$step["mean",2]*2) 
  
  sigma03s <- c(hmm_2st$mle$step["sd",1],
                hmm_2st$mle$step["sd",2],
                hmm_2st$mle$step["sd",2]*2) 
  
  kappa03s <- c(hmm_2st$mle$angle["concentration",1],
                hmm_2st$mle$angle["concentration",2],
                hmm_2st$mle$angle["concentration",2]*10)
  
  hmm_3st <- fitHMM(data_prep_cv, nbState = 3, 
                         dist = list(step = "gamma", angle = "vm"),
                         Par0 = list(step = c(mu03s, sigma03s), angle = kappa03s), 
                         formula = ~1, estAngleMean = NULL)
  
  mu04s <- c(hmm_2st$mle$step["mean",1], 
             hmm_2st$mle$step["mean",2], 
             hmm_2st$mle$step["mean",2]*2,
             hmm_2st$mle$step["mean",2]*.5) 
  
  sigma04s <- c(hmm_2st$mle$step["sd",1],
                hmm_2st$mle$step["sd",2],
                hmm_2st$mle$step["sd",2]*2,
                hmm_2st$mle$step["sd",2]*.5) 
  
  kappa04s <- c(hmm_2st$mle$angle["concentration",1],
                hmm_2st$mle$angle["concentration",2],
                hmm_2st$mle$angle["concentration",2]*10,
                hmm_2st$mle$angle["concentration",2]*.5)
  
  hmm_4st <- fitHMM(data_prep_cv, nbState = 4, 
                         dist = list(step = "gamma", angle = "vm"),
                         Par0 = list(step = c(mu04s, sigma04s), angle = kappa04s),
                         formula = ~1, estAngleMean = NULL)
  
  mu06s <- c(hmm_2st$mle$step["mean",1], 
             hmm_2st$mle$step["mean",2], 
             hmm_2st$mle$step["mean",2]*2,
             hmm_2st$mle$step["mean",2]*.5,
             hmm_2st$mle$step["mean",1]*2,
             hmm_2st$mle$step["mean",1]*.5) 
  
  sigma06s <- c(hmm_2st$mle$step["sd",1],
                hmm_2st$mle$step["sd",2],
                hmm_2st$mle$step["sd",2]*2,
                hmm_2st$mle$step["sd",2]*.5,
                hmm_2st$mle$step["sd",1]*2,
                hmm_2st$mle$step["sd",1]*.5) 
  
  kappa06s <- c(hmm_2st$mle$angle["concentration",1],
                hmm_2st$mle$angle["concentration",2],
                hmm_2st$mle$angle["concentration",2]*10,
                hmm_2st$mle$angle["concentration",2]*.5,
                hmm_2st$mle$angle["concentration",1]*10,
                hmm_2st$mle$angle["concentration",1]*.5)
  
  hmm_6st <- fitHMM(data_prep_cv, nbState = 6, 
                         dist = list(step = "gamma", angle = "vm"),
                         Par0 = list(step = c(mu06s, sigma06s), angle = kappa06s),
                         formula = ~1, estAngleMean = NULL)
  
  par0b <- getPar0(model = hmm_2st, formula = ~cover)
  cov_2st_cover <- fitHMM(data_prep_cv, nbState = 2, 
                        dist = list(step = "gamma", angle = "vm"),
                        Par0 = par0b$Par, beta0 = par0b$beta, 
                        formula = ~cover)
  
  par0b <- getPar0(model = hmm_3st, formula = ~cover)
  cov_3st_cover <- fitHMM(data_prep_cv, nbState = 3, 
                          dist = list(step = "gamma", angle = "vm"),
                          Par0 = par0b$Par, beta0 = par0b$beta, 
                          formula = ~cover)
  
  par0b <- getPar0(model = hmm_4st, formula = ~cover)
  cov_4st_cover <- fitHMM(data_prep_cv, nbState = 4, 
                          dist = list(step = "gamma", angle = "vm"),
                          Par0 = par0b$Par, beta0 = par0b$beta, 
                          formula = ~cover)
  
  par0b <- getPar0(model = hmm_2st, formula = ~fence)
  cov_2st_fence <- fitHMM(data_prep_cv, nbState = 2, 
                          dist = list(step = "gamma", angle = "vm"),
                          Par0 = par0b$Par, beta0 = par0b$beta, 
                          formula = ~fence)
  
  par0b <- getPar0(model = hmm_3st, formula = ~fence)
  cov_3st_fence <- fitHMM(data_prep_cv, nbState = 3, 
                          dist = list(step = "gamma", angle = "vm"),
                          Par0 = par0b$Par, beta0 = par0b$beta, 
                          formula = ~fence)
  
  par0b <- getPar0(model = hmm_4st, formula = ~fence)
  cov_4st_fence <- fitHMM(data_prep_cv, nbState = 4, 
                          dist = list(step = "gamma", angle = "vm"),
                          Par0 = par0b$Par, beta0 = par0b$beta, 
                          formula = ~fence)
  
  par0b <- getPar0(model = hmm_2st, formula = ~road)
  cov_2st_road <- fitHMM(data_prep_cv, nbState = 2, 
                         dist = list(step = "gamma", angle = "vm"),
                         Par0 = par0b$Par, beta0 = par0b$beta, 
                         formula = ~road)
  
  par0b <- getPar0(model = hmm_3st, formula = ~road)
  cov_3st_road <- fitHMM(data_prep_cv, nbState = 3, 
                         dist = list(step = "gamma", angle = "vm"),
                         Par0 = par0b$Par, beta0 = par0b$beta, 
                         formula = ~road)
  
  par0b <- getPar0(model = hmm_4st, formula = ~road)
  cov_4st_road <- fitHMM(data_prep_cv, nbState = 4, 
                         dist = list(step = "gamma", angle = "vm"),
                         Par0 = par0b$Par, beta0 = par0b$beta, 
                         formula = ~road)
  
  par0b <- getPar0(model = hmm_2st, formula = ~cover+fence)
  cov_2st_cov_fen <- fitHMM(data_prep_cv, nbState = 2, 
                            dist = list(step = "gamma", angle = "vm"),
                            Par0 = par0b$Par, beta0 = par0b$beta, 
                            formula = ~cover+fence)
  
  par0b <- getPar0(model = hmm_3st, formula = ~cover+fence)
  cov_3st_cov_fen <- fitHMM(data_prep_cv, nbState = 3, 
                            dist = list(step = "gamma", angle = "vm"),
                            Par0 = par0b$Par, beta0 = par0b$beta, 
                            formula = ~cover+fence)
  
  par0b <- getPar0(model = hmm_4st, formula = ~cover+fence)
  cov_4st_cov_fen <- fitHMM(data_prep_cv, nbState = 4, 
                            dist = list(step = "gamma", angle = "vm"),
                            Par0 = par0b$Par, beta0 = par0b$beta, 
                            formula = ~cover+fence)
  
  par0b <- getPar0(model = hmm_2st, formula = ~cover+fence+road)
  cov_2st_cov_fen_rd <- fitHMM(data_prep_cv, nbState = 2, 
                               dist = list(step = "gamma", angle = "vm"),
                               Par0 = par0b$Par, beta0 = par0b$beta, 
                               formula = ~cover+fence+road)
  
  par0b <- getPar0(model = hmm_2st, formula = ~fence+road)
  cov_2st_fen_rd <- fitHMM(data_prep_cv, nbState = 2, 
                               dist = list(step = "gamma", angle = "vm"),
                               Par0 = par0b$Par, beta0 = par0b$beta, 
                               formula = ~fence+road)
  
  par0b <- getPar0(model = hmm_3st, formula = ~cover+fence+road)
  cov_3st_cov_fen_rd <- fitHMM(data_prep_cv, nbState = 3, 
                               dist = list(step = "gamma", angle = "vm"),
                               Par0 = par0b$Par, beta0 = par0b$beta, 
                               formula = ~cover+fence+road)
  
  par0b <- getPar0(model = hmm_3st, formula = ~fence+road)
  cov_3st_fen_rd <- fitHMM(data_prep_cv, nbState = 3, 
                               dist = list(step = "gamma", angle = "vm"),
                               Par0 = par0b$Par, beta0 = par0b$beta, 
                               formula = ~fence+road)
  
  par0b <- getPar0(model = hmm_4st, formula = ~cover+fence+road)
  cov_4st_cov_fen_rd <- fitHMM(data_prep_cv, nbState = 4, 
                               dist = list(step = "gamma", angle = "vm"),
                               Par0 = par0b$Par, beta0 = par0b$beta, 
                               formula = ~cover+fence+road)
  
  fit_mods <- list(hmm_2st, hmm_3st, #hmm_4st, 
                   hmm_6st,
                   cov_2st_cover, cov_3st_cover, #cov_4st_cover,
                   cov_2st_fence, cov_3st_fence, #cov_4st_fence,
                   cov_2st_road, cov_3st_road, #cov_4st_road,
                   cov_2st_cov_fen, cov_3st_cov_fen, #cov_4st_cov_fen,
                   cov_2st_cov_fen_rd, cov_3st_cov_fen_rd , #cov_4st_cov_fen_rd
                   cov_2st_fen_rd, cov_3st_fen_rd
                   )
  
  #saveRDS(fit_mods, paste0("models_", unique(mod_data$id)[i], ".rds"))
}

sapply(fit_mods, AIC)
AICweights(hmm_2st, hmm_3st, #hmm_4st, 
           hmm_6st,
           cov_2st_cover, cov_3st_cover, #cov_4st_cover,
           cov_2st_fence, cov_3st_fence, #cov_4st_fence,
           cov_2st_road, cov_3st_road, #cov_4st_road,
           cov_2st_cov_fen, cov_3st_cov_fen, #cov_4st_cov_fen,
           cov_2st_cov_fen_rd, cov_3st_cov_fen_rd,#, cov_4st_cov_fen_rd
           cov_2st_fen_rd, cov_3st_fen_rd) %>% mutate(across(where(is.numeric), round, 3))
