#------------loading the library------------------

rm(list=ls(all=TRUE))
library(tidyverse)
library(ggplot2)
library(plyr)
library(dplyr)
library(spatialEco)# for finding the peaks and valleys

### Make sure you run the code for each dotted line sections 
##don't run the full code at once


#---------------------------------------------------- load the waveform dataset
acclima <- read_csv("finite_305N_30_percent_sat_calib_verif_ready_3.csv")

###First, calculating the first derivative maxima, in the first raising limb, dv1
max1_derv <- data.frame(Amplitude=integer(),
                        max_1_derv=double())

for(i in 1:11){
  dv_max1 <- max(acclima[,i+16][acclima[,1]>2300 & acclima[,1]<2800])
  dv_max1_time <- min(acclima[,1][acclima[,i+16]==dv_max1&acclima[,1]>2300 & acclima[,1]<2800])
  max1_derv <- rbind(max1_derv,c(i,dv_max1_time,dv_max1))
}
#------------------------------------------------------------------------------
colnames(max1_derv)<- c('Amplitude_number','dv_max1_time','dv_max1_value')

####Now calculating first derivative maxima in the second raising limb, dv2
max2_derv <- data.frame(Amplitude=integer(),
                        max_2_derv=double())

for(i in 1:11){
  if(i<10){
    dv_max2 <- max(acclima[,i+12][acclima[,1]>2800 & acclima[,1]<3600])
    dv_max2_time <- min(acclima[,1][acclima[,i+12]==dv_max2 & acclima[,1]>2800& acclima[,1]<3600])
  }
  if(i>9){
    dv_max2 <- max(acclima[,i+12][acclima[,1]>2800 & acclima[,1]<3600])
    dv_max2_time <- max(acclima[,1][acclima[,i+12]==dv_max2 & acclima[,1]>2800& acclima[,1]<3600])
    
  }
  max2_derv <- rbind(max2_derv,c(i,dv_max2_time,dv_max2))
}
#-----------------------------------------------------------------------------------
colnames(max2_derv)<- c('Amplitude_number','dv_max2_time','dv_max2_value')

max1_derv <- cbind(max1_derv,max2_derv$dv_max2_time)
#----------------------------------------------------------------------------------
colnames(max1_derv)<- c('Amplitude_number','dv_max1_time','dv_max1_value','dv_max2_time')

#### Calculating m, i.e. the point of reference for local minima in low permittivity mediums
max1_derv$lm_last_swath_point <- max1_derv$dv_max1_time + (0.7*(max1_derv$dv_max2_time-
                                                                  max1_derv$dv_max1_time))

max1_derv$rounded_lm_last_swath <- round_any(max1_derv$lm_last_swath_point,5,f=ceiling)
max1_derv$lm_swath_last_index <- (max1_derv$rounded_lm_last_swath/5)+1

# Number of swath points to be included for tangent from 'm'
# Default Schwartz algorithm 7
# Modeified and Model based algorithm, use Table 10 as reference
for (i in 1:11) {
  if(i<3||i==11)
    max1_derv[i,8] <- max1_derv[i,7]-7 #for depth 5cm and 0, 0.5cm
  if(i==8)
    max1_derv[i,8] <- max1_derv[i,7]-7 #for depth 3.5cm
  if(i==9)
    max1_derv[i,8] <- max1_derv[i,7]-7 #for depth 4 cm
  if(i==10||((i<5)&(i>2)))
    max1_derv[i,8] <- max1_derv[i,7]-7 #for depth 4.5cm and 1, 1.5cm 
  if(i==5)
    max1_derv[i,8] <- max1_derv[i,7]-7 #for depth 2cm
  if(i>5&i<8)
    max1_derv[i,8] <- max1_derv[i,7]-7 #for depth 3 and 2.5cm 
}
#-------------------------------------------------------------------------
names(max1_derv)[names(max1_derv) == 'V8'] <- 'lm_swath_first_index'


### Generating the tangent line equations for each depth of insertion
acclima_needed <- list()
#for drawing the linear regression lines
for (i in 1:11) {
  acclima_needed[[i]] <- data.frame(acclima[max1_derv[i,8]:max1_derv[i,7],])
}

#model <- lapply(acclima_needed,function(x) {lm(data.frame(acclima_needed[[x]][x+1])~data.frame(acclima_needed[[x]][1]))})

model <- list()

for (i in 1:11) {
  model[[i]] <- lm(unlist(acclima_needed[[i]][i+1])~unlist(acclima_needed[[i]][1]))
}

#acclima_needed <- acclima[547:554,]

#creating regression equations for each waveform in the dataset

#model <- lapply(acclima_needed[2:12], function(x)lm(x~acclima_needed$X_Time_0))

#minima_filter <- lapply(df_minima_list, function(x) { x <- filter(x, unique_minima>400 & unique_minima<1700)})

lm_data <- data.frame(Amplitude=integer(),
                      intercept=double())
for (j in 1:11) {
  row_len <- length(model[[j]]$coefficients)
  # unique_minima <- rbind(unique_minima,c(j,min(minima_filter[[j]][[row_len,1]],minima_filter[[j]][[row_len-1,1]])))
  lm_data <- rbind(lm_data,c(j-1,model[[j]][[1]]))
}
###-------------------------------------------------------------------------------
colnames(lm_data)<- c('Amplitude_number','intercept','slope')

# Calculating maximum second derivative in second rasing limb
local_second_maxima_time <- data.frame(Amplitude=integer(),
                                       local_minima_time=double())
for (i in 1:11) {
  # we are selecting the second derivative after local minima
  local_second_derivative <- acclima[,i+23][acclima[,1]>2700&acclima[,1]<3600]
  # selecting the time as well after local minima
  local_second_time <- acclima[,1][acclima[,1]>2700&acclima[,1]<3600]
  #creating data frmae with the above inputs to ease calculation
  local_second_maxima <- data.frame(local_second_time)
  local_second_maxima <- cbind(local_second_maxima,local_second_derivative)
  # calculating local maxima of the second derivative
  max_second <- max(local_second_derivative)
  #now selecting the local second maxima time
  second_maxima_time <- local_second_maxima[,1][local_second_maxima[,2]==max_second&local_second_maxima[,1]>2750&local_second_maxima[,1]<3600]
  local_second_maxima_time <- rbind(local_second_maxima_time,c(i-1,max_second,max(second_maxima_time)))
}

###---------
colnames(local_second_maxima_time)<- c('Amplitude','local_second_maxima','local_second_maxima_time')

##This is an alternative try to draw tangent from dv''2, you can run it but nothing much to infer
for (i in 1:11) {
  local_second_maxima_time[i,4] <- acclima[,i+1][acclima[,1]==local_second_maxima_time[i,3]]
}

for (i in 1:11){
  local_second_maxima_time[i,5] <- acclima[,i+12][acclima[,1]==local_second_maxima_time[i,3]]
}
#---------
colnames(local_second_maxima_time) <- c('Amplitude_number','local_second_maxima','local_second_maxima_time',
                                        'max_sec_derv_Amplitude','max_sec_derv_slope')

local_second_maxima_time$max_sec_derv_slope <- local_second_maxima_time$max_sec_derv_slope/200

local_second_maxima_time$max_sec_derv_intercept <- 
  local_second_maxima_time$max_sec_derv_Amplitude-
  (local_second_maxima_time$max_sec_derv_slope*local_second_maxima_time$local_second_maxima_time)

# calculating local maximum first derivative to draw tangent lines
local_first_maxima_time <- max2_derv
####--------------------
colnames(local_first_maxima_time)<- c('Amplitude','local_first_maxima_time','local_first_maxima')

#now selecting the respective amplitudes for local maximum derivative to draw tangent lines

for (i in 1:11) {
  local_first_maxima_time[i,4] <- acclima[,i+1][acclima[,1]==local_first_maxima_time[i,2]]
}

#-----------------------
## this tangent is just using the dv2 point, no swath points included
## this is avialable to show the difference when we include swath points
colnames(local_first_maxima_time)<- c('Amplitude','local_first_maxima_time','local_first_maxima','local_first_maxima_amplitude')

local_first_maxima_time$first_derv_slope <- local_first_maxima_time$local_first_maxima/200

local_first_maxima_time$first_derv_intercept <- 
  local_first_maxima_time$local_first_maxima_amplitude - (local_first_maxima_time$first_derv_slope * local_first_maxima_time$local_first_maxima_time)

intersect_time <- data.frame(Amplitude=integer(),
                             local_minima_time=double())
for (i in 1:11) {
  A <- rbind(c(1,-local_first_maxima_time[i,5]),
             c(1,-lm_data[i,3]))
  B <- c(local_first_maxima_time[i,6],lm_data[i,2])
  intersect_time <- rbind(intersect_time,c(i-1,solve(A,B))) 
}
#------------------------------------
colnames(intersect_time) <- c('Amplitude_number','Amplitude_intersect','intersect_time')
### Calculating the theoretical time for diffrent depths of insertion
intersect_time$air_start_time <- seq(2415,2415,0)

intersect_time$air_length <- seq(0.05,0,-0.005)

intersect_time$air_end_time <- ((20000*intersect_time$air_length)/3)+intersect_time$air_start_time

intersect_time$rounded_air_end_time <- round_any(intersect_time$air_end_time,5,f=ceiling)

intersect_time$rounded_probe_end_time <- round_any(intersect_time$intersect_time,5,f=ceiling)

intersect_time$soil_perm <- ((intersect_time$rounded_probe_end_time-intersect_time$rounded_air_end_time)*3/
                               (2*(0.05-intersect_time$air_length)*10000))^2

intersect_time$soil_length <- 0.05 - intersect_time$air_length
#### Calculating the index number, to apply the Schwarz algortihm for including swath points
local_first_maxima_time$index_maxima <- (local_first_maxima_time$local_first_maxima_time/5)+1#change it to 15 for arduino dat
#--------------------------------------------------------------------------------
### For default SChwartz Algortihm, the value is 5
### FOr modified Schwartz algortihm, use Table 9, according the mediuma and depth
### FOr model algortihm, use Table 13, Column Low
for(i in 1:11){
  if(i==5){
    local_first_maxima_time[i,8] <- local_first_maxima_time[i,7]+5 #Depth 2 cm
    local_first_maxima_time[i,9] <- local_first_maxima_time[i,7]-5
  }
  if(i==6){
    local_first_maxima_time[i,8] <- local_first_maxima_time[i,7]+5 # Depth 2.5 cm
    local_first_maxima_time[i,9] <- local_first_maxima_time[i,7]-5
  }
  if(i==7){
    local_first_maxima_time[i,8] <- local_first_maxima_time[i,7]+5 #Depth 3cm
    local_first_maxima_time[i,9] <- local_first_maxima_time[i,7]-5
  }
  if(i==11){
    local_first_maxima_time[i,8] <- local_first_maxima_time[i,7]+5 # Depth 5 cm
    local_first_maxima_time[i,9] <- local_first_maxima_time[i,7]-5
  }
  if(i==8){
    local_first_maxima_time[i,8] <- local_first_maxima_time[i,7]+5 # Depth 3.5 cm
    local_first_maxima_time[i,9] <- local_first_maxima_time[i,7]-5
  }
  if(i==9){
    local_first_maxima_time[i,8] <- local_first_maxima_time[i,7]+5  #Depth 4 cm
    local_first_maxima_time[i,9] <- local_first_maxima_time[i,7]-5
  }
  if(i==10){
    local_first_maxima_time[i,8] <- local_first_maxima_time[i,7]+5 # Depth 4.5 cm
    local_first_maxima_time[i,9] <- local_first_maxima_time[i,7]-5
  }
  if(i==1){
    local_first_maxima_time[i,8] <- local_first_maxima_time[i,7]+5 #Depth 0 cm
    local_first_maxima_time[i,9] <- local_first_maxima_time[i,7]-5
  }
  if(i==2||i==3||i==4){
    local_first_maxima_time[i,8] <- local_first_maxima_time[i,7]+5 #Depth 05,1,1.5 cm
    local_first_maxima_time[i,9] <- local_first_maxima_time[i,7]-5
  }
  
}
names(local_first_maxima_time)[names(local_first_maxima_time) == 'V8'] <- 'index_tangent_last'
names(local_first_maxima_time)[names(local_first_maxima_time) == 'V9'] <- 'index_tangent'



acclima_needed_maxima_tangent <- list()
#for drawing the linear regression lines
for (i in 1:11) {
  acclima_needed_maxima_tangent[[i]] <- data.frame(acclima[local_first_maxima_time[i,9]:local_first_maxima_time[i,8],])
}

#model <- lapply(acclima_needed,function(x) {lm(data.frame(acclima_needed[[x]][x+1])~data.frame(acclima_needed[[x]][1]))})

model_maxima <- list()

for (i in 1:11) {
  model_maxima[[i]] <- lm(unlist(acclima_needed_maxima_tangent[[i]][i+1])~unlist(acclima_needed_maxima_tangent[[i]][1]))
}

#acclima_needed <- acclima[547:554,]

#creating regression equations for each waveform in the dataset

#model <- lapply(acclima_needed[2:12], function(x)lm(x~acclima_needed$X_Time_0))

#minima_filter <- lapply(df_minima_list, function(x) { x <- filter(x, unique_minima>400 & unique_minima<1700)})

lm_data_maxima <- data.frame(Amplitude=integer(),
                             intercept=double())
for (j in 1:11) {
  #row_len <- length(model_maxima[[j]]$coefficients)
  # unique_minima <- rbind(unique_minima,c(j,min(minima_filter[[j]][[row_len,1]],minima_filter[[j]][[row_len-1,1]])))
  lm_data_maxima <- rbind(lm_data_maxima,c(j-1,model_maxima[[j]][[1]]))
}
#-----------------------------------------------------------------------
colnames(lm_data_maxima)<- c('Amplitude_number','intercept','slope')
### Calculating the correct t3, based on intersection with tangent from 'm'
intersect_time_maxima <- data.frame(Amplitude=integer(),
                                    local_minima_time=double())
for (i in 1:11) {
  A <- rbind(c(1,-lm_data[i,3]),
             c(1,-lm_data_maxima[i,3]))
  B <- c(lm_data[i,2],lm_data_maxima[i,2])
  intersect_time_maxima <- rbind(intersect_time_maxima,c(i,solve(A,B))) 
}
#------------------------------------
colnames(intersect_time_maxima) <- c('Amplitude_number','Amplitude_intersect','intersect_time')
### Calculating theoretical probe end time, to compare t3 from the above method
intersect_time$probe_end_swath_inters <- intersect_time_maxima$intersect_time

intersect_time$swath_perm <- ((intersect_time$probe_end_swath_inters-intersect_time$air_end_time)*3/
                                (2*(intersect_time$soil_length)*10000))^2

intersect_time$theoretical_probe_end_time <- (((7)^0.5)*2*(0.05-intersect_time$air_length)*10000/3)+intersect_time$rounded_air_end_time
intersect_time$swath_last_point <- max1_derv$lm_last_swath_point
intersect_time$offset_swath_theo <- intersect_time$theoretical_probe_end_time - intersect_time$probe_end_swath_inters
intersect_time$local_Second_derv_maxima <- local_second_maxima_time$local_second_maxima_time
#----------------------------------------------------------------------------------------------
write_csv(intersect_time,"model_variable_swath_proper_algo_finite_2.5_perm_TDR_305N_air_dry_3_analysis.csv")
