#Plotting single cell Whi5-YFP localization scores for strains yTY159b and yTY160a
#These are SSA (control) strains with Whi5-mCitrine, and either stable RFP (yTY159b) or RFPdegron (yTY160a)in the SSA (control) cassette
#Whi5 localization score is defined to be max5x5 - mean5x5 YFP pixel value
#Single cell trajectories of this score with overlaid budding times are plotted in './Figures/FlTrajectories/'
#Single cell trajectories of the max5x5 and mean5x5 YFP pixels are also saved


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
doxrect = geom_rect(data=doxwindow,inherit.aes=FALSE,aes(xmin=start,xmax=end,ymin=-Inf,ymax=Inf),color = "transparent",fill="green3",alpha=0.3)
doxlabel = annotate("text",label ="dox treatment",x = 4.25,y=1500,colour='red',size=5)
ylabel = ylab('YFP (a.u.)')
xlabel = xlab('time(h)')
removelegend = theme(legend.position="none")
border = theme(panel.border=element_rect(colour="black",linetype="solid",fill=NA,size=1))
linewidth = 0.2
pointsize = 0.2
hlinewidth = 0.35



#Making sure the old cells that we pick are over age 15
outfolder = './Figures/FlTrajectories/'
info = read.csv('./CombinedData/infoall2.csv')

#
bt = getbudtimes(info)
cellid = paste(info$date,info$xy,info$trap)
bts = data.frame(cellid,bt)
btsmelted = melt(bts,id.vars=c('cellid'))
btsmelted = filter(btsmelted,value>=154)
btsmelted = mutate(btsmelted,bt = value-155)

#Reading in the fluorescence data
fldata = read.csv('./CombinedData/flall.csv')
cutoff = 0.0000425

#Plot the yfpmax5x5 and yfpmean on the same plots (20 cells for each of the lanes)
strainstoplot = c('yTY159b','yTY160a')
doxrect = geom_rect(data=doxwindow,inherit.aes=FALSE,aes(xmin=start,xmax=end,ymin=0.00315,ymax=Inf),color = "transparent",fill="green3",alpha=0.3)
for(i in 1:length(strainstoplot)){
	currfl = fldata[fldata$strain==strainstoplot[i],]
	print(nrow(currfl))
	currflsplit = split(currfl,currfl$cellid)
	currflsplit = currflsplit[lapply(currflsplit,nrow)!=0]
	listofplots = list();
	print(length(currflsplit))
	for(j in 1:length(currflsplit)){
		print(j)
		currcell = currflsplit[[j]]
		currcellmelted = melt(currcell,id.vars=c('time'),measure.vars=c('yfpmean','yfpmax5by5'))
		currbts = btsmelted[as.character(btsmelted$cellid)==as.character(currcell[1,]$cellid),]
		title = currcell[1,]$cellid
		p1 <- ggplot(currcellmelted,aes(time,value))+geom_line(aes(colour=variable),size=0.3) + ggtitle(title) + ylab('YFP (a.u.)') + xlab('time(h)') + doxrect + border +fonts + removelegend + xscale  + geom_vline(data=currbts,aes(xintercept=bt),linetype='dashed',size=0.2) #+ #scale_y_continuous(limits=c(0.0032,0.0038))
		listofplots[[j]] = p1;
	}
	p1 = plot_grid(plotlist = listofplots,nrow=20,ncol=4)
	plotfilename = paste('./Figures/FlTrajectories/','yfp_max5x5_mean5x5',strainstoplot[i],'.pdf',sep="")
	print(plotfilename)
	ggsave(file = plotfilename,p1,dpi=600,width=6,height=30)
	
}




#Plots of only max5x5minusmean5x5
doxrect = geom_rect(data=doxwindow,inherit.aes=FALSE,aes(xmin=start,xmax=end,ymin=0,ymax=Inf),color = "transparent",fill="green3",alpha=0.3)
for(i in 1:length(strainstoplot)){
	currfl = fldata[fldata$strain==strainstoplot[i],]
	print(nrow(currfl))
	currflsplit = split(currfl,currfl$cellid)
	currflsplit = currflsplit[lapply(currflsplit,nrow)!=0]
	listofplots = list();
	print(length(currflsplit))
	for(j in 1:length(currflsplit)){
		print(j)
		currcell = currflsplit[[j]]
		currcellmelted = melt(currcell,id.vars=c('time'),measure.vars=c('maxminusmean'))
		currbts = btsmelted[as.character(btsmelted$cellid)==as.character(currcell[1,]$cellid),]
		title = currcell[1,]$cellid
		p1 <- ggplot(currcellmelted,aes(time,value))+geom_line(aes(colour=variable),size=0.3) + ggtitle(title) + ylab('YFP (a.u.)') + xlab('time(h)') + doxrect + border +fonts + removelegend + xscale  + geom_vline(data=currbts,aes(xintercept=bt),linetype='dashed',size=0.2) + geom_hline(yintercept = cutoff,colour="blue",size=0.2)
		listofplots[[j]] = p1;
	}
	p1 = plot_grid(plotlist = listofplots,nrow=20,ncol=4)
	plotfilename = paste('./Figures/FlTrajectories/','yfp_maxminusmean5x5',strainstoplot[i],'.pdf',sep="")
	print(plotfilename)
	ggsave(file = plotfilename,p1,dpi=600,width=6,height=30)
	
}


#Plotting a single example trajectory. 3/25/19, xy 24, cell 4 is a good example
fonts = theme(axis.text=element_text(size=18), axis.title=element_text(size=20),strip.text.x = element_text(size = 18))
currcell = filter(fldata,date=='3_25_19' & xy ==24 & trap==4)
currcellmelted = melt(currcell,id.vars=c('time'),measure.vars=c('maxminusmean'))
currbts = btsmelted[as.character(btsmelted$cellid)=='3_25_19 24 4',]
p1 <- ggplot(currcellmelted,aes(time,value))+geom_line(aes(colour=variable),size=1) + geom_point(aes(colour=variable))+ggtitle(title) + ylab('YFP (a.u.)') + xlab('time(h)') + doxrect + border +fonts + removelegend + xscale  + geom_vline(data=currbts,aes(xintercept=bt),linetype='dashed',size=0.8) + geom_hline(yintercept = cutoff,colour="blue",size=0.8) + scale_y_continuous(expand=c(0,0))
plotfilename = paste(outfolder,'ExampleCell_3_25_19_xy24_trap4.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=5,height=4)
removeggplottextleaveticks(p1,plotfilename,4,600,5,4)



#See how mean fluorescence varies between the two experiments at different times
p1 = ggplot(fldata,aes(time,yfpmean))+ geom_smooth(aes(colour=strain))
timestrainlevels = paste(1:51,c('yTY159b','yTY160a'))
timestrain = factor(paste(fldata$time,fldata$strain),levels=timestrainlevels);
fldata = data.frame(timestrain,fldata)
p1 = ggplot(fldata[!is.na(timestrain),],aes(time,yfpmean))+ geom_boxplot(notch=FALSE,aes(group=timestrain,colour=strain),position=position_dodge(width=0.1))


