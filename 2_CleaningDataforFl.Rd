setwd('/Users/thomasyoung/Dropbox/MovieProcessing/Whi5Localization_maxips')
#Determining the cells for which there are fluorescence measurements. 
#This is based on the time that the traps are fully occupied and the last observation time

source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/timeseries_func.Rd')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/func.Rd')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/Preprocessing_func.Rd')
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
