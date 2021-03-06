#Assessing thresholds for the Whi5-YFP localization score used to call Whi5-YFP nuclear localized.
#Whi5-YFP localization score = max5x5 - mean5x5 YFP pixel intensity
#Different cutoffs are applied to each cell budding interval. Within each budding interval, the number of contiguous blocks for which the cell is above or below the cutoff are assessed for agreement with the prior literature showing Whi5 is nuclear localized once per cell cycle
#For each cutoff, the fraction of budding intervals with 0 blocks below, exactly 1 block above, and > 1 block above the cutoff are calculated.  A desireable threshold will yield a low fraction with > 1 block above, low fraction with 0 blocks below, and relatively high fraction with exactly 1 block above.
#The cutoff value of 0.0000425 
#Results are saved in './Figures/ThresholdsForYFP+/'


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
doxrect = geom_rect(data=doxwindow,inherit.aes=FALSE,aes(xmin=start,xmax=end,ymin=0.005,ymax=Inf),color = "transparent",fill="green3",alpha=0.3)
doxlabel = annotate("text",label ="dox treatment",x = 4.25,y=1500,colour='red',size=5)
ylabel = ylab('YFP (a.u.)')
xlabel = xlab('time(h)')
removelegend = theme(legend.position="none")
border = theme(panel.border=element_rect(colour="black",linetype="solid",fill=NA,size=1))
linewidth = 0.2
pointsize = 0.2
hlinewidth = 0.35
themes = theme(axis.text=element_text(size=14), axis.title=element_text(size=16),strip.text.x = element_text(size = 14),panel.border = element_rect(colour="black",linetype="solid",fill=NA,size=1))
themes = theme(axis.text=element_text(size=17), axis.title=element_text(size=16),strip.text.x = element_text(size = 14),panel.border = element_rect(colour="black",linetype="solid",fill=NA,size=1))

outfolder = './Figures/ThresholdsForYFP+/'
fl = read.csv('./CombinedData/flall.csv')
info = read.csv('./CombinedData/infoall2.csv')

flid = fl$cellid
fl = data.frame(flid,fl)
flwide <- dcast(fl,flid~time,value.var = 'maxminusmean')

infoid = paste(info$date,info$xy,info$trap)
flid = as.character(flwide$flid)
matches = match(infoid,flid)
matchesinfoindex = which(!is.na(matches))
matchesflindex = matches[!is.na(matches)]
bt = getbudtimes(info)
bt = cbind(bt,info$lastobservationtime)
bt = bt-154;


flwide = flwide[2:ncol(flwide)]
flwidelist = mattolol(flwide)
btlist = mattolol(bt)
btlist = btlist[matchesinfoindex]
flwidelist = flwidelist[matchesflindex]


#cutoffs to test.  We apply each cutoff, then acess the fraction
#of 'bad' results (lots of fluctuations)
#also bad is everthing below
cutoffs = mean(unlist(flwidelist),na.rm=TRUE) + sd(unlist(flwidelist),na.rm=TRUE)*seq(-1,5,0.05)
numblocksabove = list()
numblocksbelow = list()
for(i in 1:length(cutoffs)){
	currco = cutoffs[i]
	coapplied = lapply(flwidelist,greaterthan,currco)
	biwindows = mapply(extractwindows1, coapplied, btlist)
	biwindows = unlist(biwindows,recursive=FALSE)
	biwindows = lapply(biwindows,f<-function(x){return(x[!is.na(x)])})
	rles = lapply(biwindows,rle)
	rlelengths = lapply(rles, getrlelengths)
	rlevalues = lapply(rles, getrlevalues)
	numblocks = lapply(rlelengths,length)
	numblocksabove[[i]] = lapply(rlevalues,f<-function(x){return(sum(x))})
	numblocksbelow[[i]] = lapply(rlevalues,f<-function(x){return(sum(!x))})
}

#For each choice of a cutoff, we plot the fraction of budding intervals with numblocks #above > 2 (this is bad), numblocksbelow = 0 (this is also bad)
fractoomanyabove = lapply(numblocksabove,f<-function(x){frac = sum(unlist(x)>=2,na.rm=TRUE)/length(unlist(x));return(frac)})
fracnobelow = lapply(numblocksbelow,f<-function(x){frac = sum(unlist(x)==0,na.rm=TRUE)/length(unlist(x));return(frac)})
fraconeabove = lapply(numblocksabove,f<-function(x){frac = sum(unlist(x)==1,na.rm=TRUE)/length(unlist(x));return(frac)})
plot(unlist(fracnobelow),unlist(fractoomanyabove))


toplot= data.frame(unlist(fracnobelow),unlist(fractoomanyabove),unlist(fraconeabove),cutoffs)
colnames(toplot) = c('fracnobelow','fractoomanyabove','fraconeabove','cutoffs')
p1 = ggplot(toplot,aes(fracnobelow,fractoomanyabove,colour=cutoffs)) + geom_point(size=1) + scale_colour_gradientn(colours=terrain.colors(10)) + ylab('fraction with >= 2 blocks above') + xlab('fraction with 0 blocks below') + scale_y_continuous(breaks=seq(0,1,0.05))
ggsave('./Figures/ThresholdsForYFP+/fractionofbi_greaterthan2blocksoverco_vs_0blocksbelowcutoff.pdf',p1,dpi=100,width=5,height=4)


#Plotting fraction to many above vs cutoff
p1 = ggplot(toplot,aes(cutoffs,fractoomanyabove)) + geom_point(size=1) + scale_colour_gradientn(colours=terrain.colors(10)) + ylab('fraction with >= 2 blocks above') + xlab('cutoffs') + scale_y_continuous(breaks=seq(0,1,0.05))
ggsave('./Figures/ThresholdsForYFP+/fractionofbi_greaterthan2blocksoverco_vs_cutoff.pdf',p1,dpi=100,width=5,height=4)


p1 = ggplot(toplot,aes(cutoffs,fraconeabove)) + geom_point(size=1) + scale_colour_gradientn(colours=terrain.colors(10)) + ylab('fraction with exactly 1 block above') + xlab('cutoffs') + scale_y_continuous(breaks=seq(0,1,0.05))
ggsave('./Figures/ThresholdsForYFP+/fractionofbi_exactly1blocksoverco_vs_cutoff.pdf',p1,dpi=100,width=5,height=4)


#Plotting all fractions vs the cutoffs along with the final chosen cutoff
cutoff = 0.0000425
toplotmelted = melt(toplot,id.vars = c('cutoffs'))
toplotmelted$variable = factor(toplotmelted$variable,levels = c('fracnobelow','fraconeabove','fractoomanyabove'))
themes1 = theme(axis.text=element_text(size=12), axis.title=element_text(size=16),strip.text.x = element_text(size = 14),panel.border = element_rect(colour="black",linetype="solid",fill=NA,size=1))
p1 = ggplot(toplotmelted,aes(cutoffs,value,colour=variable)) + geom_line(size=0.5) + geom_point()  + ylab('Fraction') + xlab('Cutoffs [x10^-5]') + scale_y_continuous(breaks=seq(0,1,0.1))  +scale_x_continuous(breaks=seq(0.00001,0.00009,0.00001),labels = seq(1,9,1),limits = c(0.00001,0.00009))  + themes1 + border + geom_vline(xintercept = cutoff)
plotfilename = plotfilename = paste(outfolder,'AllFractionsVsCutoff.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=5,height=4.5)
removeggplottext(p1,plotfilename,4,600,5,4.5)
p1 = p1 + theme(axis.title=element_blank(),text=element_blank())
ggsave(file = paste(outfolder,'AllFractionsVsCutoff_onlynumbers.pdf',sep=''),p1,dpi=600,width=5,height=4.5)