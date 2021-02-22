#Based on the combined data from multiple experiments, make a list of the cells that are of the appropriate ages for study
#We only include cells from the old experiments >= age 15 at the time of doxycycline exposure

setwd('/Users/thomasyoung/Dropbox/MovieProcessing/Whi5Localization_maxips')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/timeseries_func.Rd')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/func.Rd')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/Preprocessing_func.Rd')
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
