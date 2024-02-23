rm(list=ls(all=TRUE))
library(tidyverse)
library(ggplot2)
library(plyr) #to round numbers
library(dplyr)
library(spatialEco)# for finding the peaks and valleys

#Follow the steps as written below
#Make sure you run the code till where the lines are drawn in the code
#The codes are compartmentalized based on the lines, if run completely at a 
#single time you will not get accurate results
####-----------------------------------------------------------
#loading acclima amplitude, first, and second derivatives
#Make sure you load the correct file to avoid confusion
#FOr example, don't load water and use methanol permittivity values
acclima <- read_csv("finite_305N_saturation_calib_verif_ready_3.csv")

#since we are interested in only in first reflection for higher permittivity 

# we are selecting the first 1001 rows for analysis

# This is based on water, which has a permittivity of 78 

acclima <- acclima[1:1101, ]

# calculating the maxima and minima using spatial eco functions
#creating minima list for the given waveforms

# first selecting the amplitude columns from acclima for calculating minima

acclima_minima <- acclima %>% select(contains('Ampli'))
# ignoring the air amplitude data to avoid confusion, because air is lower permittivity
acclima_minima <- acclima_minima[,2:11]#11 for lab data

#creating minima list for the given waveforms to draw tangent from minima
# using the column names to run and find the local minima in each column, ignoring air
df_minima_list <- lapply(1:10,function(j) if(j%%1==0) data.frame(unique_minima = unique(local.min.max(acclima_minima[[paste0("Y_Amplitude_", j)]])$minima)))

# now we have many unique local minima values for each column

# to avoid confusion and select the correct minima

# first we are going to remove values less than 300
minima_filter <- lapply(df_minima_list, function(x) { x <- filter(x, unique_minima>400 & unique_minima<1700)})#for acclima it is 400 1700


#minima_filter <- lapply(7:10,function(x) { x <- filter(x, minima_filter[[x]][1]>600 & minima_filter[[x]][1]<800)})#for acclima it is 400 1700

#now, we have to select the required local minima for tangent line approach
unique_minima <- data.frame(Amplitude=integer(),
                            Minima=double())
for (j in 2:10) {
  row_len <- length(minima_filter[[j]]$unique_minima)
  unique_minima <- rbind(unique_minima,c(j,minima_filter[[j]][[row_len,1]]))
}
#----------------------------------------------------------------------
colnames(unique_minima)<- c('Amplitude','minima')

#now we have identified the unique minima, the next step is to identify the correct time

# in many cases unique minima can repeat, so choosing the right one is a challenge

# for that we check if the first derivative is zero and second derivative is positive

# in many cases the above condition can be valid multiple times, in those situations

# we chose the one with the maximum positive second derivative out of the avaialble


medium_minima <- data.frame(Amplitude=integer(),
                            local_minima_time=double())
for (i in 2:10) {
  #minima_time_filter <- data.frame(min_time= double(),
  #min_first=double(),
  #min_second=double())
  #extracting first derivative values
  minima_first_derivative <- acclima[,i+13][acclima[,i+2]==unique_minima[i-1,2]]
  # extracting second derivative values
  minima_second_derivative <- acclima[,i+24][acclima[,i+2]==unique_minima[i-1,2]]
  # extracting the time where local minima was found
  minima_time_select <- acclima[,1][acclima[,i+2]==unique_minima[i-1,2]]
  # combining the above three to create a data frame to apply the condition
  minima_time_filter <- data.frame(minima_time_select)
  minima_time_filter <- cbind(minima_time_filter,minima_first_derivative,minima_second_derivative)
  #selecting the minima time,verifying mathematically
  minima_time <- minima_time_filter[,1][minima_time_filter[,2]==0 & minima_time_filter[,3]>0]
  medium_minima <- rbind(medium_minima,c(i,unique_minima[i-1,2],max(minima_time)))
  if(length(minima_time)==0){
    minima_time <-minima_time_filter[,1][minima_time_filter[,3]>0]
    medium_minima <- rbind(medium_minima,c(i,unique_minima[i-1,2],max(minima_time)))
  }
}
#At this point you will get warning with Inf values, please ignore it and run the next section
#----------------------------------------------------------------------------
colnames(medium_minima)<- c('Amplitude','local_minima','local_minima_time')

#removing the -inf rows

medium_minima <- medium_minima[!is.infinite(rowSums(medium_minima)),]

#Now we have the local minima for the tangent ready.

#now we have to find local maxima of the second derivative

local_second_maxima_time <- data.frame(Amplitude=integer(),
                                       local_minima_time=double())
for (i in 2:10) {
  # we are selecting the second derivative after local minima
  local_second_derivative <- acclima[,i+24][acclima[,1]>medium_minima[i-1,3]]
  # selecting the time as well after local minima
  local_second_time <- acclima[,1][acclima[,1]>medium_minima[i-1,3]]
  #creating data frmae with the above inputs to ease calculation
  local_second_maxima <- data.frame(local_second_time)
  local_second_maxima <- cbind(local_second_maxima,local_second_derivative)
  # calculating local maxima of the second derivative
  max_second <- max(local_second_derivative)
  #now selecting the local second maxima time
  second_maxima_time <- local_second_maxima[,1][local_second_maxima[,2]==max_second]
  local_second_maxima_time <- rbind(local_second_maxima_time,c(i,max_second,max(second_maxima_time)))
}
#----------------------------------------------------------------------------------
colnames(local_second_maxima_time)<- c('Amplitude','local_second_maxima','local_second_maxima_time')

# calculating local maximum first derivative to draw tangent lines
local_first_maxima_time <- data.frame(Amplitude=integer(),
                                      local_minima_time=double())
for (i in 2:10) {
  # we are selecting the first derivative after local minima
  local_first_derivative <- acclima[,i+13][acclima[,1]>medium_minima[i-1,3]]
  # selecting the time as well after local minima
  local_first_time <- acclima[,1][acclima[,1]>medium_minima[i-1,3]]
  #creating data frmae with the above inputs to ease calculation
  local_first_maxima <- data.frame(local_first_time)
  local_first_maxima <- cbind(local_first_maxima,local_first_derivative)
  # calculating local maxima of the second derivative
  max_first <- max(local_first_derivative)
  #now selecting the local second maxima time
  first_maxima_time <- local_first_maxima[,1][local_first_maxima[,2]==max_first]
  local_first_maxima_time <- rbind(local_first_maxima_time,c(i,max_first,max(first_maxima_time)))
}
#--------------------------------------------------------------------------------------
colnames(local_first_maxima_time)<- c('Amplitude','local_first_maxima','local_first_maxima_time')

#now selecting the respective amplitudes for local maximum derivative to draw tangent lines

for (i in 2:10) {
  local_first_maxima_time[i-1,4] <- acclima[,i+2][acclima[,1]==local_first_maxima_time[i-1,3]]
}

# selecting and adding the required columns into one data frame 
local_first_maxima_time[,5] <- local_second_maxima_time[,3]

local_first_maxima_time[,c(6,7)] <- medium_minima[,c(2,3)] 
#------------------------------------------------------------------------------------------
#### Tangent from dv2′(local first derivative maxima in second raising limb)

colnames(local_first_maxima_time)<- c('number','local_first_derivative_maximum',
                                      'local_first_derivative_maximum_time',
                                      'amplitude_first_derivative_maximum',
                                      'local_second_serivative_maximum_time',
                                      'local_minima_amplitude',
                                      'local_minima_time')

#alternate_Approach_using_swath_of_points_from_first_Derivative_maximum for drawing tangents
local_first_maxima_time$index_maxima <- (local_first_maxima_time$local_first_derivative_maximum_time/5)+1 #change it to 15 for arduino dat

#Note depth 0 and 0.5 cm are omitted for all high perm mediums to avoid confusion
#Right now 5 is used for all depths of insertion
#for modified Schwartz Algorithm follow Table 9 results from paper
#for model based estimation follow Table 13, High column data

for (i in 1:9) {
  if(i==9){
    local_first_maxima_time[i,9] <- local_first_maxima_time[i,8]+5#Depth 5 cm
    local_first_maxima_time[i,10] <- local_first_maxima_time[i,8]-5
  }
  if(i==8){
    local_first_maxima_time[i,9] <- local_first_maxima_time[i,8]+5 #Depth 4.5 cm
    local_first_maxima_time[i,10] <- local_first_maxima_time[i,8]-5
  }
  if(i==7){
    local_first_maxima_time[i,9] <- local_first_maxima_time[i,8]+5 #Depth 4 cm
    local_first_maxima_time[i,10] <- local_first_maxima_time[i,8]-5
  }
  if(i==6){
    local_first_maxima_time[i,9] <- local_first_maxima_time[i,8]+5 #Depth 3.5 cm
    local_first_maxima_time[i,10] <- local_first_maxima_time[i,8]-5
  }
  if(i==5){
    local_first_maxima_time[i,9] <- local_first_maxima_time[i,8]+5 #Depth 3 cm
    local_first_maxima_time[i,10] <- local_first_maxima_time[i,8]-5
  }
  if(i==4){
    local_first_maxima_time[i,9] <- local_first_maxima_time[i,8]+5 #Depth 2.5 cm
    local_first_maxima_time[i,10] <- local_first_maxima_time[i,8]-5
  }
  
  if(i<4){
    local_first_maxima_time[i,9] <- local_first_maxima_time[i,8]+5#Depth 2 cm to 1cm
    local_first_maxima_time[i,10] <- local_first_maxima_time[i,8]-5
  }
  
}
names(local_first_maxima_time)[names(local_first_maxima_time) == 'V9'] <- 'index_tangent_last'
names(local_first_maxima_time)[names(local_first_maxima_time) == 'V10'] <- 'index_tangent'


##Storing the swath points for each depth if insertion
acclima_needed <- list()
#for drawing the linear regression lines
for (i in 1:9) {
  acclima_needed[[i]] <- data.frame(acclima[local_first_maxima_time[i,10]:local_first_maxima_time[i,9],])
}


model <- list()

for (i in 1:9) {
  model[[i]] <- lm(unlist(acclima_needed[[i]][i+3])~unlist(acclima_needed[[i]][1]))
}


lm_data <- data.frame(Amplitude=integer(),
                      intercept=double())
for (j in 1:9) {
  row_len <- length(model[[j]]$coefficients)
  # unique_minima <- rbind(unique_minima,c(j,min(minima_filter[[j]][[row_len,1]],minima_filter[[j]][[row_len-1,1]])))
  lm_data <- rbind(lm_data,c(j-1,model[[j]][[1]]))
}
###------------------------------------------------------------------------
colnames(lm_data)<- c('Amplitude_number','intercept','slope')

# getting the intercept for the tangent
local_first_maxima_time$tangent_intercept_Swath <- lm_data$intercept

# getting the slope for the tangent
local_first_maxima_time$tangent_slope_Swath <- lm_data$slope

##Identifying the point of intersection for t3, using intersection of two lines concept
local_first_maxima_time$intersection_time_swath <- 
  ((local_first_maxima_time$local_minima_amplitude - 
      local_first_maxima_time$tangent_intercept_Swath)/(local_first_maxima_time$tangent_slope_Swath))


#calculating tangent intercept for local first derivative maximum
local_first_maxima_time$tangent_intercept <- 
  local_first_maxima_time$amplitude_first_derivative_maximum - (((local_first_maxima_time$local_first_derivative_maximum)/200)*local_first_maxima_time$local_first_derivative_maximum_time)

# calculating tangent slope at the first derivative maximum
local_first_maxima_time$tangent_slope <- (local_first_maxima_time$local_first_derivative_maximum)/200

# now we have to calculate the intersection time of tangent and minima
local_first_maxima_time$intersection_time <- 
  ((local_first_maxima_time$local_minima_amplitude - 
      local_first_maxima_time$tangent_intercept)/(local_first_maxima_time$local_first_derivative_maximum*0.005))

# now finalizing the travel end time for Schwartz algorithm
# we check if second order maxima time is between local minima time and intersection time
local_first_maxima_time$probe_end_time <- ifelse(local_first_maxima_time$local_second_serivative_maximum_time
                                                 < local_first_maxima_time$intersection_time , local_first_maxima_time$local_second_serivative_maximum_time,
                                                 ifelse(local_first_maxima_time$intersection_time < local_first_maxima_time$local_second_serivative_maximum_time ,local_first_maxima_time$intersection_time, 0))

# rounding the numbers to match with the resolution of waveform
local_first_maxima_time$probe_end_time <- round_any(local_first_maxima_time$probe_end_time, 5, f = ceiling)

#depth of insertion in medium in meters
local_first_maxima_time$travel_length_medium <- seq(0.01,0.05,0.005)

#Calculating travel time in medium using Equation 5 in paper
#Use correct permittivity value for the medium from the Table 6 in paper
local_first_maxima_time$travel_time_media <- ((35.4^0.5)*
                                                (2*local_first_maxima_time$travel_length_medium)*(10000))/3

# Initializing t1 as constant for all depths of insertion
local_first_maxima_time$air_start_time <- seq(2415,2415,0)

#Calculating t2 using Equation 3 in paper
local_first_maxima_time$air_end_time <- ((20000*(0.05-local_first_maxima_time$travel_length_medium))/3)+ local_first_maxima_time$air_start_time

# Calculating Theoretical probe end time using Equation 5 in paper
local_first_maxima_time$theoretical_probe_end_time <- local_first_maxima_time$air_end_time + local_first_maxima_time$travel_time_media

# ROunding the theoretical probe end time to match the resolution of the sensor
local_first_maxima_time$rounded_theo_probe_end_time <- round_any(local_first_maxima_time$theoretical_probe_end_time,5)

# Calculating the difference between Theoretical and calculated t3
local_first_maxima_time$swath_offset_time <- 
  local_first_maxima_time$theoretical_probe_end_time - 
  local_first_maxima_time$intersection_time_swath

##Schwartz used 5 point swath for tangent from minima as well
##However, in our study we did not use it but the program has that part to depict it
local_first_maxima_time$index_minima <- (local_first_maxima_time$local_minima_time/5)+1

local_first_maxima_time$index_minima_tangent_first <- local_first_maxima_time$index_minima-5

local_first_maxima_time$index_minima_tangent_last <- local_first_maxima_time$index_minima+5

acclima_needed_minima <- list()
#for drawing the linear regression lines
for (i in 1:9) {
  acclima_needed_minima[[i]] <- data.frame(acclima[local_first_maxima_time[i,26]:local_first_maxima_time[i,27],])
}

#model <- lapply(acclima_needed,function(x) {lm(data.frame(acclima_needed[[x]][x+1])~data.frame(acclima_needed[[x]][1]))})

model_minima <- list()

for (i in 1:9) {
  model_minima[[i]] <- lm(unlist(acclima_needed_minima[[i]][i+3])~unlist(acclima_needed_minima[[i]][1]))
}


lm_data_minima <- data.frame(Amplitude=integer(),
                             intercept=double())
for (j in 1:9) {
  #row_len <- length(model[[j]]$coefficients)
  # unique_minima <- rbind(unique_minima,c(j,min(minima_filter[[j]][[row_len,1]],minima_filter[[j]][[row_len-1,1]])))
  lm_data_minima <- rbind(lm_data_minima,c(j-1,model_minima[[j]][[1]]))
}
###------------------------------------------------------------------------------
colnames(lm_data_minima)<- c('Amplitude_number','intercept','slope')

intersect_time <- data.frame(Amplitude=integer(),
                             local_minima_time=double())
for (i in 1:9) {
  A <- rbind(c(1,-local_first_maxima_time[i,12]),
             c(1,-lm_data_minima[i,3]))
  B <- c(local_first_maxima_time[i,11],lm_data_minima[i,2])
  intersect_time <- rbind(intersect_time,c(i+1,solve(A,B))) 
}
#-------------------------------------------------------------------------------

colnames(intersect_time) <- c('Amplitude_number','Amplitude_intersect','intersect_time')
local_first_maxima_time$minima_swath_intersection <- intersect_time$intersect_time

local_first_maxima_time$rounded_minima_swath_inters <- round_any(local_first_maxima_time$minima_swath_intersection,5,f= ceiling)

for (i in 1:9) {
  hema <- summary(model[[i]])
  local_first_maxima_time[i,30] <- hema$r.squared
}
names(local_first_maxima_time)[names(local_first_maxima_time) == 'V30'] <- 'r_Square_maxima'
##Exporting the final results for future analysis
write_csv(local_first_maxima_time,"proper_algo_finite_swath_35.4_finite_saturation_3.csv")
##################----------------------------------------------------------------
#Please follow the procedure carefully to visualize the waveform and results.
#To visualize the waveform analysis please use the below code
#For example, In the below code we visualize for 5 cm depth of insertion, Y_Amplitude_10
#For 4.5 cm, Y_Amplitude_9, all the row values will start at 8. For 4 cm, Y_Amplitude_8, row values at 7
#keep this mind to avoid unnecessary confusion and make sure you change it and verify it
#So all the rows in array will be 9, don't change the column [row,column]
# for plotting the swath points included in tangent
#To select for 4.5 cm [[8]][1] for times, [[8]][11] for selecting the amplitude
#To select for 4 cm [[7]][1] for times, [[7]][10] for selecting the amplitude

ggplot(data=acclima, aes(x=X_Time_0, y=Y_Amplitude_10, group=1)) +
  geom_line(size=2,color = "#FFDB6D")+
  geom_hline(aes(yintercept = local_first_maxima_time[9,6], color= "local minima tangent"),size = 1.5 )+#local minma amplitude line
  geom_abline(aes(intercept = lm_data[9,2], slope = lm_data[9,3],color="first derivative maxima tangent"),size = 1.5)+ #tangent at first maxima derivative
  geom_vline(aes(xintercept = local_first_maxima_time[9,13], color= "Probe end time (t3,Schwartz)"),size = 1.5 ,linetype="dashed")+
  geom_vline(aes(xintercept = local_first_maxima_time[9,21], color= "air end time (t2)"),size = 1.5 ,linetype="dashed")+
  geom_vline(aes(xintercept = local_first_maxima_time[9,20], color= "air start time (t1)"),size = 1.5 ,linetype="dashed")+
  theme_bw()+
  geom_vline(aes(xintercept = local_first_maxima_time[9,22], color= "theoretical probe end time"),size = 1.5 ,linetype="dashed")+
  theme(panel.grid = element_blank())+
  #Here we plot the number of swath points included
  geom_point(data = acclima_needed[[9]],
             aes(x = unlist(acclima_needed[[9]][1]),#for [[9][1] selects the time for 5 cm depth of insertion
                 y = unlist(acclima_needed[[9]][12]),# for [[9]][12] selects the amplitude to be plotted for 5 cm depth of insertion
                 color = "Swath points used"),  # Map the label to the color aesthetic
             shape = 16)+
  labs(
    x = "Time in picoseconds",
    y = "Signal amplitude"
  )+
  scale_color_manual(
    values = c(
      "local minima tangent" = "blue", 
      "first derivative maxima tangent" = "red",
      "Probe end time (t3,Schwartz)" = "black",
      "air end time (t2)" = "purple",
      "air start time (t1)" = "grey",
      "theoretical probe end time" = "green",
      "Swath points used"="#56B4E9"
    )
  )+
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(face = "bold", size = 15)) +
  theme(axis.title.x = element_text(size = 18)) +
  theme(axis.title.y = element_text(size = 18)) +
  theme(plot.title = element_text(hjust = 0.5, size = 20))+
  #Need to adjust this accordingly for aesthetics and not to miss information in plot
  geom_segment(aes(x=4550,
                   y=1250,
                   xend=5300,
                   yend=1450),arrow = arrow(length = unit(0.5,"cm")))+
  geom_text(aes(x = 3600, y = 1200, label = "Error = 35 ps"), 
            hjust = -0.2, vjust = 0.5, size = 6)+
  labs(color = "5 cm in saturated soil Trial 3")+
  theme(legend.position = c(0.215, 0.7), 
        legend.text = element_text(size = 11),  # Adjust legend text size
        legend.title = element_text(size = 14))  # Adjust legend title size
