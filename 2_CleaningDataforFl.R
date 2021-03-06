#Determining the cells for which there are fluorescence measurements. 
#This is based on the time that the traps are fully occupied and the last observation time
#Only cells that are born before time 43 (for young experiments), and cells in the old experiments are considered. Young cells have to also divide at least once prior to doxycycline removal
#The resulting list of cells is saved in "./CombinedData/infoall2.csv"

setwd('/Users/thomasyoung/Dropbox/MovieProcessing/Whi5Localization_maxips_git')
source('./functions/timeseries_func.Rd')
source('./functions/func.Rd')
source('./functions/Preprocessing_func.Rd')

library(dplyr)
library(ggplot2)
library(reshape2)
library(gridExtra)

info = read.csv('./CombinedData/infoall.csv')


#Also make sure young cells are born at or before time point 43.
condition = (info$birth <= 43 & info$doxtime<50) | (info$doxtime>100);
info = filter(info,condition)

#Make sure young cells divide at least once before dox removal
bt = getbudtimes(info)
ageatdoxremoval = ageattimes(bt,info$doxtime+24)
condition = ageatdoxremoval >= 1;
info = filter(info,condition);


write.csv(info,"./CombinedData/infoall2.csv",row.names=FALSE);
