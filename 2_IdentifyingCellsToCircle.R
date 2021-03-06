#Based on the combined data from multiple experiments, make a list of the cells that are of the appropriate ages for study. These cells will be circled to measure Whi5-YFP fluorescence. 
#We only include cells from the old experiments >= age 15 at the time of doxycycline exposure
#Tables listing the cells to circle are saved in './TrapsToMeasure_oldover15/'. The smallest and largest fluorescent snapshot indices for which to circle the cells (when the cells are alive and uncensored) is also saved in the table.

setwd('/Users/thomasyoung/Dropbox/MovieProcessing/Whi5Localization_maxips_git')
source('./functions/timeseries_func.Rd')
source('./functions/func.Rd')
source('./functions/Preprocessing_func.Rd')

library(dplyr)
library(ggplot2)
library(reshape2)
library(gridExtra)


#Making sure the old cells that we pick are over age 15
info = read.csv('./CombinedData/infoall2.csv')

lastobsflindex = mapply(lastflindexbefore,info$lastobservationtime,info$flfreq,info$lastibeforepart2+1)
firstflindexafterdox = 15;

condition = lastobsflindex >= firstflindexafterdox;
info = filter(info,condition)

bt = getbudtimes(info)
ageatdox = ageattimes(bt,info$doxtime)
condition = (info$doxtime<100) | (info$doxtime>100 & ageatdox>=15)
info = filter(info,condition)
lastobsflindex = mapply(lastflindexbefore,info$lastobservationtime,info$flfreq,info$lastibeforepart2+1)
firstobsflindex = mapply(firstflindexafter,info$birth,info$flfreq,info$lastibeforepart2+1)
firstobsflindex[firstobsflindex<1]=1

trapstoinspect = select(info,xy,trap)
trapstoinspect = cbind(trapstoinspect,firstobsflindex,lastobsflindex);


trapstoinspect = split(trapstoinspect,paste0(info$date,info$strain))
trapstoinspect = trapstoinspect[lapply(trapstoinspect,nrow)!=0]


lapply(1:length(trapstoinspect),function(i) write.csv(trapstoinspect[[i]],file=paste0('./TrapsToMeasure_oldover15/','aliveafterdoxadded_',names(trapstoinspect[i]),'oldover15.csv'),row.names=FALSE))
