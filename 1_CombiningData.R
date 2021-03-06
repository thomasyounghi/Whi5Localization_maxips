#Make 1 big data frame each from the budding time data in the 'Data' folder
#Identifies and labels replicates (experiments sharing the same strain and age of dox induction)
#Add labels to be used for publication (redundant with the official strain names, doxtime, etc)
#Cell data saved in "./CombinedData/infoall.csv"


setwd('/Users/thomasyoung/Dropbox/MovieProcessing/Whi5Localization_maxips_git')
source('./functions/timeseries_func.Rd')
source('./functions/func.Rd')
source('./functions/Preprocessing_func.Rd')

library(dplyr)
library(stringr)

#Open up all the files in data, to get one big file of info, YFP, RFP
infofilestoopen = list.files("./Data/",pattern="[[:alnum:]]*revised2.csv")
infofilestoopen = paste('./Data/',infofilestoopen,sep="")


info = concatenatefiles(infofilestoopen)

#Check for gaps in the budding time matrix
bt = getbudtimes(info)
btgaps = apply(bt,1,checkgaps)
btdecreasing = apply(bt,1,checkdecreasing)
btrepeats = apply(bt,1,checkrepeats)
table(btgaps)
table(btdecreasing)
table(btrepeats)
write.csv(info[btgaps,],'./ProblematicBTs/NAbetweenbts.csv',row.names=FALSE)
write.csv(info[btdecreasing,],'./ProblematicBTs/decreasingbts.csv',row.names=FALSE)
write.csv(info[btrepeats,],'./ProblematicBTs/repeatedbts.csv',row.names=FALSE)


#In situations where there are repeated buddingtime, remove the repeated buddingtimes and shift the later budding times over
condition = btrepeats & !btdecreasing & !btgaps
bttofix = bt[condition,]
fixedbt = t(apply(bttofix,1,removerepeatsinincreasing))
birthcol = which(names(info)=='birth')
info[condition,(birthcol+1):ncol(info)]=fixedbt

#Eliminate all rows for which budding times are decreasing or there are gaps
condition = !(btdecreasing | btgaps)
table(condition)
info = filter(info,condition)




#Plot the birth distributions (boxplot). Also calculate the doxtime, which is important for all experiments.
doxtime = info$lastibeforepart2 + info$firstpart2ipostdoxstart;
info = cbind(doxtime,info);
neworder = order(info$strain,doxtime,info$date,info$lane)
info = info[neworder,]

#Creating the experiment field: 'strain timeofdoxycylinetreatment date lane'
experiment = paste(info$strain,info$doxtime, info$date,info$lane)
experiment = factor(experiment)

#Providing an id to each cell
id = 1:nrow(info)

#identifying replicates. Experiments sharing the same strain and doxycycline treatment time, but derived from different colonies. Order of the labeling specified by (date lane)
replabel = identifyreplicates(info);


#Labels to use for plotting
expage = factor(info$doxtime)
levels(expage) = c('young','old')
expstrain = factor(info$strain,levels=c('yTY159b','yTY160a'))
levels(expstrain) = c('SSAcontrol','SSAcontrol+RFPdegron')


#adding the additional features to the info data frame
#labels that make things easier to understand.  
bt = getbudtimes(info)
ageatdox = ageattimes(bt,info$doxtime)

info = data.frame(id,expage,ageatdox,expstrain,replabel,experiment,info)
date = as.character(info$date)
date = str_match(date,"(\\w+)_")
date = date[,2]
info$date=date;


#Checking the time between birth and the first bud. 
#If it is less than 20, then shift budding time so the time of the first bud is the new birth time. 
#In all cases it is 30 minutes which is reasonable
table(info$X1 -info$birth)
timetofb = (info$X1-info$birth)<=2
timetofb1 = timetofb & (!is.na(timetofb))
bt = getbudtimes(info)
bt1=bt;
bt1[which(timetofb1),1:(ncol(bt)-1)]= bt[which(timetofb1),2:ncol(bt)]
bt1[which(timetofb1),ncol(bt1)] = NA
info[,((which(colnames(info)=='birth'))+1):ncol(info)]=bt1

write.csv(info,"./CombinedData/infoall.csv",row.names=FALSE);

