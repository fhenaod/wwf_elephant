library(tidyverse)

raw_data <- read.csv("data/Namibia_elephant_2010-Mar2021.csv")
raw_data %>% head()

mod_data <- raw_data %>%
  mutate(date.time = as.POSIXct(date.time, format = "%Y-%m-%d %H:%M", tz = "UTC"),
         date = as.Date(date.time),
         year = as.numeric(strftime(as.POSIXlt(date.time), format = "%Y")),
         month = as.numeric(strftime(as.POSIXlt(date.time), format = "%m")),
         week = as.numeric(strftime(as.POSIXlt(date.time), format = "%W")),
         day = as.numeric(strftime(as.POSIXlt(date.time), format = "%j")))

mod_data %>% head()

hist(sort(table(mod_data$id.new)), xlab = "N locations", 
     main = paste0("Distribution of individual sample size \n N = ", 
                   length(unique(mod_data$id.new)), " individuals"))

mod_data <- mod_data %>% mutate(numID = group_indices(., id.new))

# plot sampling times per individual
mod_data %>% ggplot(aes(x = date.time, y = numID, color = as.factor(id.new))) + 
  geom_point() + theme_classic() +
  theme(legend.position = "none") +
  scale_color_viridis_d()

# check for potential duplicated records
which(mod_data$date.time[1:(nrow(mod_data)-1)] == mod_data$date.time[2:(nrow(mod_data))])

# first register id
foo <- which(mod_data$numID[1:(nrow(mod_data)-1)] != mod_data$numID[2:nrow(mod_data)])  
mod_data <- mod_data %>% mutate(first_reg = rep(0, nrow(mod_data)))
mod_data$first_reg[foo+1] <- 1
mod_data$first_reg[1] <- 1

length(unique(mod_data$numID)) == sum(mod_data$first_reg)
mod_data[sort(c(foo-1,foo,foo+1)),c('numID','date.time','first_reg')]

# step length #### 
# euclidian distance
euc_dist <- sqrt((mod_data$lat[1:(nrow(mod_data)-1)] - mod_data$lat[2:(nrow(mod_data))])^2 + 
                 (mod_data$long[1:(nrow(mod_data)-1)] - mod_data$long[2:(nrow(mod_data))])^2)

summary(euc_dist)
mod_data$euc_dis <- c(NA, euc_dist) # add NA to first register in df
mod_data$euc_dis <- 
  ifelse(mod_data$first_reg == 1, NA, mod_data$euc_dis) # add NA to first register per indiv

# step lengths - if non projected coordinates
library(fossil)
cor_dist <- deg.dist(mod_data$long[1:(nrow(mod_data)-1)], mod_data$lat[1:(nrow(mod_data)-1)], 
                     mod_data$long[2:nrow(mod_data)], mod_data$lat[2:nrow(mod_data)])

summary(cor_dist)
mod_data$cor_dist <- c(NA, cor_dist) # add NA to first register in df
mod_data$cor_dist <- 
  ifelse(mod_data$first_reg == 1, NA, mod_data$cor_dist) # add NA to first register per indiv

mod_data <- mod_data %>% mutate(cor_dist = cor_dist*1000) # convert to m from km

hist(mod_data$euc_dis - mod_data$cor_dist)

# cummulative step length hist by sex
mod_data %>% group_by(id.new, sex) %>% 
  summarize(cum_d = sum(cor_dist, na.rm = T)) %>% 
  ggplot(aes(x = cum_d)) + 
  geom_histogram() + theme_classic() +
  labs(x = "Cummulative step distance (m)", y = "Count") +
  facet_wrap(~sex, nrow = 2, scales = "free")

mod_data %>% group_by(id.new, sex, year) %>% 
  summarize(cum_d = sum(cor_dist, na.rm = T)) %>% 
  mutate(cum_d = cum_d/1000) %>% 
  arrange(desc(cum_d)) %>% head()

# mean step length hist by sex
mod_data %>% group_by(id.new, sex) %>% 
  summarize(m = mean(cor_dist, na.rm = T)) %>% 
  ggplot(aes(x = m)) + geom_histogram() + theme_classic() +
  labs(x = "Mean step distance (m)", y = "Count") +
    facet_wrap(~sex, nrow = 2)

# mean step length boxplot by year and sex
mod_data %>% group_by(id.new, sex, year, month) %>% 
  summarize(m = mean(cor_dist, na.rm = T)) %>% 
  ggplot(aes(x = as.factor(year), y = m)) + geom_violin() + 
  theme_classic() +
  labs(x = "Year", y = "Mean step distance (m)") +
  facet_wrap(~sex, nrow = 2)

# time between steps ####
tm_bs <- difftime(mod_data$date.time[2:(nrow(mod_data))], 
                  mod_data$date.time[1:(nrow(mod_data)-1)], 
                  units = "mins")
mod_data$dt <- c(NA, tm_bs)
mod_data$dt <- ifelse(mod_data$first_reg == 1, NA, mod_data$dt)
summary(mod_data$dt)
hist(mod_data$dt[mod_data$dt < 300], 
     main = "", xlab = "Time between steps (minutes)")

# step length and time-lag between steps ####
library(mgcv)
plot(gam(cor_dist ~ s(dt), data = mod_data))
plot(gam(cor_dist ~ s(dt), data = mod_data, subset = dt < 300))

# Turn angle ####
angs <- earth.bear(mod_data$long[1:(nrow(mod_data)-1)], 
                        mod_data$lat[1:(nrow(mod_data)-1)], 
                   mod_data$long[2:nrow(mod_data)], 
                            mod_data$lat[2:nrow(mod_data)])
mod_data$angs <- c(NA, angs)
mod_data$angs <- ifelse(mod_data$first_reg == 1, NA, mod_data$angs)

mod_data %>% ggplot(aes(x = angs)) +
  geom_histogram() + theme_classic() +
  labs(x = "Angles", y = "Count")  +
  facet_wrap(~sex, nrow = 2, scales = "free")

# estimate trajectories
library(adehabitatLT)
traj <- as.ltraj(xy = mod_data[,c("long","lat")], 
                 date = mod_data$date.time, id = mod_data$numID)
summary(traj[[1]])
plot(traj[[1]])

traj_df <- ld(traj)
library(CircStats)

# angles (radians) histogram, just 4 to check
par(mfrow = c(2,2))
indvs <- levels(as.factor(traj_df$id))
for(i in 1:4){
  hist((traj_df[traj_df$id==indvs[i] & traj_df$dt < 6000,'rel.angle']), 
       main = indvs[i], xlab = "Turning angle (radians)")
}

# angles (degrees) histogram, just 4 to check
par(mfrow = c(2,2))
for(i in 1:4){
  hist(na.omit(traj_df[traj_df$id==indvs[i] & traj_df$dt < 6000,'rel.angle'])*180/pi, 
       main = indvs[i], xlab = "Turning angle (degrees)")
}

# Rose diagram
par(mfrow = c(2,2))
for(i in 1:4){
  rose.diag(na.omit(traj_df[traj_df$id==indvs[i] & traj_df$dt < 6000,'rel.angle'])*180/pi, 
            bins = 30, prop = 2, main = indvs[i])
}

saveRDS(mod_data, "mod_data.rds")
saveRDS(traj_df, "traj_df.rds")

mod_data <- readRDS("mod_data.rds")

# change point analysis ######
library(bcpa)

i = 1

X <- mod_data %>% filter(numID == i) %>% pull(long)
Y <- mod_data %>% filter(numID == i) %>% pull(lat)
Time <- mod_data %>% filter(numID == i) %>% pull(date.time)

# step lengths, anges and movement stats
mytrack <- MakeTrack(X,Y,Time)
plot(mytrack)
mytrack.VT <- GetVT(mytrack)
head(mytrack.VT)
summary(mytrack.VT)

par(mfrow = c(1,2))
hist(mytrack.VT$V, breaks = 20, col = "grey")
hist(mytrack.VT$Theta, breaks = 20, col = "grey")
par(mfrow = c(1,1))

# window sweep  to find break points
mytrack.ws <- WindowSweep(mytrack.VT, "V*cos(Theta)", # Vcos (persistent velocity)
                          windowsize = 50, progress = T, K = 2)

plot(mytrack.ws, type = "smooth") 
plot(mytrack.ws, type = "flat")

PathPlot(mytrack, mytrack.ws, type = "flat", clusterwidth = 3, main = "Flat BCPA")
PathPlot(mytrack, mytrack.ws, type = "smooth", main = "Smooth BCPA")
DiagPlot(mytrack.ws) # Diagnostics

PhasePlot(mytrack.ws, type = "smooth", clusterwidth = 3)
