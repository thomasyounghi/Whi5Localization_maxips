#Calculating Whi5-YFP localization scores for imaged cells as the difference between the 
#maximum 5x5 pixel value minus the mean 5x5 pixel value
#YFP measurements are annotated with non fluorescent data for the same cell (found in ./CombinedData/infoall2.csv)
#Adjusted fluorescent indices for cells collected on 3/12/19 by adding 1. This is because this movie was started 10 minutes later than usual
#Annotated fluorescent values are saved in './CombinedData/flall.csv'


setwd('/Users/thomasyoung/Dropbox/MovieProcessing/Whi5Localization_maxips_git')
source('./functions/timeseries_func.Rd')
source('./functions/func.Rd')
source('./functions/Preprocessing_func.Rd')

library(dplyr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(scales)
library(cowplot)
library(stringr)

#figure settings:
theme_set(theme_cowplot())
ylogscale = scale_y_continuous(trans=log10_trans(),breaks=trans_breaks("log10",function(x) 10^x),labels=trans_format("log10",math_format(10^.x)),limits=c(0.2,2000))
xscale = scale_x_continuous(limits = c(0,60),breaks = seq(0,60,12),labels=as.character(seq(0,10,2)))
fonts = theme(axis.text=element_text(size=4), axis.title=element_text(size=4),strip.text.x = element_text(size = 1),plot.title=element_text(size=6))
ylims = ylim(-10,1000)
xlims = xlim(0,12.5)
doxwindow = data.frame(15,39);
colnames(doxwindow) = c('start','end')
doxrect = geom_rect(data=doxwindow,inherit.aes=FALSE,aes(xmin=start,xmax=end,ymin=0.0003,ymax=Inf),color = "transparent",fill="green3",alpha=0.3)
doxlabel = annotate("text",label ="dox treatment",x = 4.25,y=1500,colour='red',size=5)
ylabel = ylab('YFP (a.u.)')
xlabel = xlab('time(h)')
removelegend = theme(legend.position="none")
border = theme(panel.border=element_rect(colour="black",linetype="solid",fill=NA,size=1))
linewidth = 0.2
pointsize = 0.2
hlinewidth = 0.35

#Getting the non-fluorescence cell data to join to the fluorescence data
info = read.csv('./CombinedData/infoall2.csv')

#Reading in the fluorescence data
flfilestoopen = list.files("./FLData/",pattern="[[:alnum:]]*manualroimeasurements.csv")
flfilestoopen = paste('./FlData/',flfilestoopen,sep="")
fldata = concatenatefiles(flfilestoopen)

#generating a date variable from the filename
fldate = fldata$filenames
date = str_match(fldate,"./FlData/(\\w+?)_(\\w+?)_(\\w+?)_yTY")
date = paste(date[,2],'_',date[,3],'_',date[,4],sep='')
fldata = data.frame(date,fldata)

#Assigning each row of the fldata to a strain based on the matches to data in info
infodatatoadd = select(info,c('date','xy','trap','strain','replabel'))
fldata = inner_join(fldata,infodatatoadd)

#Adding a cell id variable to the fldata
cellid = paste(fldata$date,fldata$xy,fldata$trap)
fldata = data.frame(cellid,fldata)

#Subtracting yfpmean5by5 intensity from yfpmax5by5
fldata = mutate(fldata,maxminusmean = yfpmax5by5 - yfpmean5by5)

#Correct the indices from 3/12/19.  Add 1 to each of the time indices since the movie was started 10 minutes later than usual
timecorrection  = rep(0,nrow(fldata))
timecorrection[fldata$date=='3_12_19'] = 1
fldata$time = fldata$time + timecorrection
write.csv(fldata,'./CombinedData/flall.csv',row.names=FALSE)


