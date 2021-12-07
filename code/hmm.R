library(momentuHMM)
library(sp)
library(raster)
library(tidyverse)

raw_data <- read.csv("data/Namibia_elephant_2010-Mar2021.csv")
raw_data %>% head()

hmm_data <- raw_data %>%
  mutate(date.time = as.POSIXct(date.time, format = "%Y-%m-%d %H:%M", tz = "UTC")) %>% 
  filter(id.new == unique(mod_data$id)[i])

t_int <- as.numeric(diff(hmm_data$date.time), units = "hours")
plot(t_int, ylab = "Time interval (hr)", main = unique(hmm_data$id.new),
     pch = 19, cex = 0.5, las = 1
     , ylim = c(0,10)
     )
hist((t_int), xlab = "Time interval (hr)", main = unique(hmm_data$id.new))

# regularize data, creates 2hr sequence, from firts to last data point
ti <- seq(hmm_data$date.time[1], hmm_data$date.time[length(hmm_data$date.time)],
          by = 2*60*60)

# Interpolate location 
iLoc <- as.data.frame(cbind(lon = approx(hmm_data$date.time, hmm_data$lon, xout = ti)$y, 
                            lat = approx(hmm_data$date.time, hmm_data$lat, xout = ti)$y))
data_reg <- cbind(date = ti, iLoc)

plot(hmm_data$lon, hmm_data$lat, pch = 19, cex = 0.5, 
     xlab = "Lon", ylab = "Lat", las = 1, main = unique(hmm_data$id.new))
points(data_reg$lon, data_reg$lat, pch = 19, cex = 0.5, col = rgb(1,0, 0, 0.5))

# step and angles estimation
data_prep <- prepData(data_reg, type = "LL", coordNames = c("lon", "lat"))
data_prep %>% head()
plot(data_prep)

sum(data_prep$step == 0, na.rm = T)
int_pars <- data_prep %>% dplyr::summarize(mean_step = mean(step, na.rm = T), 
                                           med_step = median(step, na.rm = T), 
                                           sd_step = sd(step, na.rm = T), 
                                           mean_ang = mean(angle, na.rm = T), 
                                           sd_ang = sd(angle, na.rm = T))




# starting values
mu0 <- c(int_pars$mean_step, int_pars$med_step)
sigma0 <- c(int_pars$sd_step, int_pars$sd_step*2)
kappa0 <- c(int_pars$sd_ang, int_pars$sd_ang*10)

# fit two-state hhm
hmm_2st <- fitHMM(data_prep, nbState = 2, 
                  dist = list(step = "gamma", angle = "vm"),
                  Par0 = list(step = c(mu0, sigma0), angle = kappa0),
                  formula = ~1, estAngleMean = NULL)
plot(hmm_2st)
data_states <- viterbi(hmm_2st)
table(data_states)

sts_probs <- stateProbs(hmm_2st)
plotStates(hmm_2st)  

hmm_2st_rf <- fitHMM(data_prep, nbState = 2, 
                  dist = list(step = "gamma", angle = "vm"),
                  Par0 = list(step = c(mu0, sigma0), angle = kappa0),
                  formula = ~1, estAngleMean = NULL, retryFits = 1)
c(original = hmm_2st$mod$minimum, retryFits = hmm_2st_rf$mod$minimum)

plotPR(hmm_2st) # check pseudoresiduals 

# # fit three-state hhm
mu03s <- c(hmm_2st_rf$mle$step["mean",1], 
           hmm_2st_rf$mle$step["mean",2], 
           hmm_2st_rf$mle$step["mean",2]*2) 

sigma03s <- c(hmm_2st_rf$mle$step["sd",1],
              hmm_2st_rf$mle$step["sd",2],
              hmm_2st_rf$mle$step["sd",2]*2) 

kappa03s <- c(hmm_2st_rf$mle$angle["concentration",1],
              hmm_2st_rf$mle$angle["concentration",2],
              hmm_2st_rf$mle$angle["concentration",2]*10)

hmm_3st <- fitHMM(data_prep, nbState = 3, 
                  dist = list(step = "gamma", angle = "vm"),
                  Par0 = list(step = c(mu03s, sigma03s), angle = kappa03s),
                  formula = ~1, estAngleMean = NULL)

hmm_3st_rf <- fitHMM(data_prep, nbState = 3, 
                  dist = list(step = "gamma", angle = "vm"),
                  Par0 = list(step = c(mu03s, sigma03s), angle = kappa03s),
                  formula = ~1, estAngleMean = NULL, retryFits = 1)

data.frame(model = c("two_st", "two_st_rf", "three_st", "hmm_3st_rf"),
           logLik = c(hmm_2st$mod$minimum, 
                      hmm_2st_rf$mod$minimum, 
                      hmm_3st$mod$minimum,
                      hmm_3st_rf$mod$minimum),
           AIC = AIC(hmm_2st, hmm_2st_rf, hmm_3st, hmm_3st_rf)$AIC,
           AICw = AICweights(hmm_2st, hmm_2st_rf, hmm_3st, hmm_3st_rf)$weight) %>% 
  arrange(desc(AICw))

data.frame(
  mean_step = hmm_2st$mle$step["mean",],
  sd_step = hmm_2st$mle$step["sd",],
  turn_con = hmm_2st$mle$angle["concentration",]
) %>% round(3)

# Simulate data
data_sim <- simData(model = hmm_2st, nbAnimals = 1, 
                    obsPerAnimal = dim(data_prep)[1])
model_sim_fit <- fitHMM(data_sim, nbState = length(hmm_2st$stateNames), 
                       dist = list(step = "gamma", angle = "vm"), 
                       Par0 = list(step = c(hmm_2st$mle$step["mean",],
                                            hmm_2st$mle$step["sd",]), 
                                   angle = hmm_2st$mle$angle["concentration",]), 
                       formula = ~1, estAngleMean = NULL)
plot(model_sim_fit)

# Spatial covariables #####
plot(ld_cov_crop, 
     xlim = c(min(data_reg$lon), max(data_reg$lon)), 
     ylim = c(min(data_reg$lat), max(data_reg$lat)))
points(hmm_data$lon, hmm_data$lat, pch = 20, cex = 0.5, las = 1)
#points(data_reg$lon, data_reg$lat, pch = 19, cex = 0.5, col = rgb(1,0, 0, 0.5))
plot(fence_sh, col = "magenta", add = T, lty = 'dotted', lwd = 2)
plot(roads_sh, col = "yellow", add = T , lty = 'dotdash', lwd = 2)

# prep data with covars
data_prep_cv <- prepData(data_reg, type = "LL", coordNames = c("lon", "lat"), 
                         spatialCovs = list(cover = ld_cov_crop, 
                                            fence = fence_raster, 
                                            road = road_raster))
data_prep_cv %>% head()
table(data_prep_cv$cover) ; table(data_prep_cv$fence) ; table(data_prep_cv$road)

hmm_2st_cov <- fitHMM(data_prep_cv, nbState = 2, 
                  dist = list(step = "gamma", angle = "vm"),
                  Par0 = list(step = c(mu0, sigma0), angle = kappa0),
                  formula = ~1, estAngleMean = NULL)

# fit the model with the covariates and new start pars
par0b <- getPar0(model = hmm_2st_cov, formula = ~cover)
cov_2st_mod <- fitHMM(data_prep_cv, nbState = 2, 
                      dist = list(step = "gamma", angle = "vm"),
                      Par0 = par0b$Par, beta0 = par0b$beta, 
                      formula = ~cover)

AICweights(hmm_2st_cov, cov_2st_mod)

plot(cov_2st_mod)
sim_cov_dat <- simData(model = cov_2st_mod, nbAnimals = 1, 
                    obsPerAnimal = dim(data_prep_cv)[1])

p1 <- data_prep_cv %>% ggplot(aes(x = x, y = y, color = as.factor(cover))) +
  geom_point() + theme_minimal() + 
  theme(legend.position = "none", 
        axis.text.x = element_blank(), axis.text.y = element_blank())
p2 <- sim_cov_dat %>% ggplot(aes(x = x, y = y, color = as.factor(cover))) + 
  geom_point() + theme_minimal() + 
  theme(legend.position = "none", 
        axis.text.x = element_blank(), axis.text.y = element_blank())

ggpubr::ggarrange(p2, p1, 
                  labels = c("A", "B"),
                  ncol = 2, nrow = 1)
ggsave(filename = "figures/panel_emp_sim.png", plot = last_plot(), 
       width = 500, height = 400, units = "mm")


