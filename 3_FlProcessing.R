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
library(scales)
library(cowplot)
library(stringr)

#figure settings:
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



#Making sure the old cells that we pick are over age 15
info = read.csv('./CombinedData/infoall2.csv')

#
bt = getbudtimes(info)
cellid = paste(info$date,info$xy,info$trap)
bts = data.frame(cellid,bt)
btsmelted = melt(bts,id.vars=c('cellid'))
btsmelted = filter(btsmelted,value>=154)
btsmelted = mutate(btsmelted,bt = value-155)

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

#Correct the indices from 3/12/19.  Add 1 to each of the time indices
timecorrection  = rep(0,nrow(fldata))
timecorrection[fldata$date=='3_12_19'] = 1
fldata$time = fldata$time + timecorrection
write.csv(fldata,'./CombinedData/flall.csv',row.names=FALSE)

#Plot the yfpmax5x5 and yfpmean on the same plots (20 cells for each of the lanes)
strainstoplot = c('yTY159b','yTY160a')
for(i in 1:length(strainstoplot)){
	currfl = fldata[fldata$strain==strainstoplot[i],]
	print(nrow(currfl))
	currflsplit = split(currfl,currfl$cellid)
	currflsplit = currflsplit[lapply(currflsplit,nrow)!=0]
	listofplots = list();
	print(length(currflsplit))
	for(j in 1:30){
		print(j)
		currcell = currflsplit[[j]]
		currcellmelted = melt(currcell,id.vars=c('time'),measure.vars=c('yfpmean','maxminusmean','yfpmax5by5'))
		currbts = btsmelted[as.character(btsmelted$cellid)==as.character(currcell[1,]$cellid),]
		title = currcell[1,]$cellid
		p1 <- ggplot(currcellmelted,aes(time,value))+geom_line(aes(colour=variable),size=0.3) + ggtitle(title) + ylab('YFP (a.u.)') + xlab('time(h)') + doxrect + border +fonts + removelegend + xscale  + geom_vline(data=currbts,aes(xintercept=bt),linetype='dashed',size=0.2) + scale_y_continuous(limits=c(0.0032,0.0035))
		listofplots[[j]] = p1;
	}
	p1 = plot_grid(plotlist = listofplots,nrow=8,ncol=4)
	plotfilename = paste('./Figures/FlTrajectories/','yfp_',strainstoplot[i],'.pdf',sep="")
	print(plotfilename)
	ggsave(file = plotfilename,p1,dpi=600,width=6,height=12)
	
}



#See how mean fluorescence varies between the two experiments at different times
p1 = ggplot(fldata,aes(time,yfpmean))+ geom_smooth(aes(colour=strain))
timestrainlevels = paste(1:51,c('yTY159b','yTY160a'))
timestrain = factor(paste(fldata$time,fldata$strain),levels=timestrainlevels);
fldata = data.frame(timestrain,fldata)
p1 = ggplot(fldata[!is.na(timestrain),],aes(time,yfpmean))+ geom_boxplot(notch=FALSE,aes(group=timestrain,colour=strain),position=position_dodge(width=0.1))


