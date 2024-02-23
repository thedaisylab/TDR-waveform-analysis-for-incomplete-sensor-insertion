# Now we are generalizing the script for any given waveform dataset

#------------loading the library------------------

rm(list=ls(all=TRUE))
library(tidyverse)
library(ggplot2)
library(dplyr)
library(spatialEco)# for finding the peaks and valleys

#---------------------------------------------------- load the waveform dataset
#-------here we transpose the datset, calculate the first derivative, and append it to transposed dataset---------

# loading original dataset, before transpose
acclima_original <- read_csv("TDR_methanol_305N_verification_for_R.csv")

# transposing the original dataset
# first remember the names
n <- acclima_original$`data type`
# this will convert the data for each depth of insertion
acclima_transpose <- as.data.frame(t(acclima_original[,-1]))
colnames(acclima_transpose) <- n

colnames(acclima_transpose)[colnames(acclima_transpose)      # Rename two variable names
                            %in% c("X_Time","Y_Amplitude")] <- c("X_Time_0","Y_Amplitude_0")

##### Now we are calculating the first derivative of the signal amplitude
### Note the values are mutiplied by 200 for better visualization
m <- length(n)

for (i in 1:m) {
  for (j in  1:nrow(acclima_transpose)) {
    if(i%%2==0){
      if(j==1){
        acclima_transpose[j,m+(i/2)] <- ((acclima_transpose[2,i] - acclima_transpose[1,i])/5)*200#no average for the first row
      }
      if(j>1 & j!=4096){
        acclima_transpose[j,m+(i/2)] <- ((acclima_transpose[j+1,i]- acclima_transpose[j-1,i])/10)*200 #taking the average of two rows
      }
      if(j==4096){
        acclima_transpose[j,m+(i/2)] <- ((acclima_transpose[4096,i]-acclima_transpose[4095,i])/5)*200 #no average for the last row
        names(acclima_transpose)[ncol(acclima_transpose)] <- paste0("derivative_", i/2-1)
      }  
    }
  }
}

#now changing the NA to zero in first rows of derivative

#acclima_transpose[is.na(acclima_transpose)] = 0

p <- ncol(acclima_transpose)# to calculate the total number of columns involved in analysis

#----now we have to calculate the second derivative
## by taking the average of two rows in common
##----now calulating the second derivative
##Note the values are multiplied by 200 for better visualization
for (i in 1:m) {
  for (j in  1:nrow(acclima_transpose)) {
    if(i%%2==0){
      if(j==1){
        acclima_transpose[j,p+(i/2)] <- ((acclima_transpose[2,i] - 2*acclima_transpose[1,i]+ acclima_transpose[1,i])/25)*200#no average for the first row
      }
      if(j>1 & j!=4096){
        acclima_transpose[j,p+(i/2)] <- ((acclima_transpose[j+1,i]- 2*acclima_transpose[j,i]+ acclima_transpose[j-1,i])/25)*200 #taking the average of two rows
      }
      if(j==4096){
        acclima_transpose[j,p+(i/2)] <- ((acclima_transpose[4096,i]-acclima_transpose[4095,i])/25)*200 #no average for the last row
        names(acclima_transpose)[ncol(acclima_transpose)] <- paste0("seccond_derivative_", i/2-1)
      }  
    }
  }
}




##---- now we are selecting the only columns of interest for further analysis

# first selecting the time for all the analysis
acclima_interest <- acclima_transpose %>% select(contains('X_Time_0'))
# next selecting the columns with amplitude
acclima_interest <- cbind(acclima_interest,acclima_transpose %>% select(contains('Ampli')))
# next selecting the columns with derivative
acclima_interest <- cbind(acclima_interest,acclima_transpose %>% select(contains('derivative')))
# creating the csv file for future analysis for medium and length analysis
write_csv(acclima_interest,"finite_305N_methanol_calib_verif_ready.csv")
