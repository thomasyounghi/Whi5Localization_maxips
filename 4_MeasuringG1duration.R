#Measuring G1 durations, comparing between strains yTY159b and yTY160a
#We apply the threshold to each cell budding interval and assess the number of contiguous blocks for which the cell is above the
#threshold and below it.

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
doxrect = geom_rect(data=doxwindow,inherit.aes=FALSE,aes(xmin=start,xmax=end,ymin=0.005,ymax=Inf),color = "transparent",fill="green3",alpha=0.3)
doxlabel = annotate("text",label ="dox treatment",x = 4.25,y=1500,colour='red',size=5)
ylabel = ylab('YFP (a.u.)')
xlabel = xlab('time(h)')
removelegend = theme(legend.position="none")
border = theme(panel.border=element_rect(colour="black",linetype="solid",fill=NA,size=1))
linewidth = 0.2
pointsize = 0.2
hlinewidth = 0.35
#figure settings:
themes = theme(axis.text=element_text(size=14), axis.title=element_text(size=16),strip.text.x = element_text(size = 14),panel.border = element_rect(colour="black",linetype="solid",fill=NA,size=1))
themes = theme(axis.text=element_text(size=17), axis.title=element_text(size=16),strip.text.x = element_text(size = 14),panel.border = element_rect(colour="black",linetype="solid",fill=NA,size=1))
ylims = ylim(-0,25)
ebwidth = 0.25;
bwidth = 0.4;
fillc = "gray66"
#cartcoordx= coord_cartesian(ylim = c(0,1),clip='off')
xscale = scale_x_discrete(expand = c(0.6,0))
yscale = scale_y_continuous(breaks = seq(0,1,0.25),labels=seq(0,100,25),limits = c(0,1.1))
border = theme(panel.border=element_rect(colour="black",linetype="solid",fill=NA,size=1))


outfolder = './Figures/g1durationcomparison/'
fl = read.csv('./CombinedData/flall.csv')
info = read.csv('./CombinedData/infoall2.csv')
info = filter(info,ageatdox>=15)

flid = fl$cellid
fl = data.frame(flid,fl)
fl = mutate(fl,yfp = maxminusmean)
flwide <- dcast(fl,flid~time,value.var = 'yfp')

infoid = paste(info$date,info$xy,info$trap)
flid = as.character(flwide$flid)
matches = match(infoid,flid)
matchesinfoindex = which(!is.na(matches))
matchesflindex = matches[!is.na(matches)]
bt = getbudtimes(info)
ageatfl = ageattimes(bt,rep(155,nrow(bt)))
info = data.frame(ageatfl,info)
bt = bt-154;


#We only measure full budding intervals 
flwide = flwide[2:ncol(flwide)]
flwidelist = mattolol(flwide)
btlist = mattolol(bt)
btlist = btlist[matchesinfoindex]
flwidelist = flwidelist[matchesflindex]
infosmall = info[matchesflindex,]
cellrepcount = table(infosmall$strain,infosmall$date)
write.csv(cellrepcount,'./Figures/NumberOfCellsPerRep_G1comparison.csv')


#We apply this threshold and for each budding interval record the pattern of blocks
cutoff = 0.0000425
coapplied = lapply(flwidelist,greaterthan,cutoff)
names(coapplied) = flid[matchesflindex]
biwindows = mapply(extractwindows1, coapplied, btlist)
biwindows = unlist(biwindows,recursive=FALSE)
biwindows = lapply(biwindows,f<-function(x){return(x[!is.na(x)])})
bilengths = lapply(biwindows,length)
rles = lapply(biwindows,rle)
rlelengths = lapply(rles, getrlelengths)
rlevalues = lapply(rles, getrlevalues)
numblocks = lapply(rlelengths,length)
numblocksabove = lapply(rlevalues,f<-function(x){return(sum(x))})
numblocksbelow = lapply(rlevalues,f<-function(x){return(sum(!x))})
abovelengths = mapply(gettrurlelengths,rlevalues,rlelengths)
maxabovelength = lapply(abovelengths,max,na.rm=TRUE)

#Problem: 3/12/19 26 6 has only 2 measured budding intervals. when the is a 3rd
#corresponds to entry 377 and 378 out of 616 entries.  The 3rd budding interval
#is the last, incomplete one (since we didn't measure fl for all time points in the budding interval)
rlelengths[377:379]
biwindows[377:379]


#Checking whether are data makes sense.  1. There should be no Whi5 localization 
#when we observe a small bud so we check whether the first block of the rle
#is true or false. 2. We also assess whether the number of blocks follows a reasonable distribution.  3.  In the cells with more than 1 'above' block, is there a clear block that corresponds to the real G1?

firstblockvalue = lapply(rlevalues,f<-function(x){return(x[1])})
table(unlist(firstblockvalue))
table(unlist(numblocks))
table(unlist(bilengths),unlist(numblocks))
#Graph to summarize the number of problematic results (like >=2 blocks above threshold)


#extracting the id information for each budding interval so we can better
#understand how budding intervals are distributed among cells and how they are 
#affected by strain.
#We add strain information by matching to the info table
binames = names(numblocks);
binamesparts = str_match(binames,"(\\w+) (\\w+) (\\w+).bi:(\\d+)")
binamesparts = data.frame(binamesparts)
colnames(binamesparts) = c('biname','date','xy','trap','bi')
binamesparts$xy = as.numeric(as.character(binamesparts$xy))
binamesparts$trap = as.numeric(as.character(binamesparts$trap))
idsforbi = binamesparts
cellid = paste(idsforbi$date,idsforbi$xy,idsforbi$trap)
idsforbi = data.frame(cellid,idsforbi)


idsforinfo = data.frame(info$strain,info$date,info$replabel,info$xy,info$trap,info$ageatdox,info$ageatfl)
colnames(idsforinfo) = c('strain','date','replabel','xy','trap','ageatdox','ageatfl')
allids = left_join(idsforbi,idsforinfo)
nrow(allids)
nrow(idsforbi)






#Making a data frame with the relevant cell and budding interval id information with the extracted G1 data
bilength = unlist(bilengths)
maxabovelength = unlist(maxabovelength)
numblocks = unlist(numblocks)
firstblockabove = unlist(firstblockvalue)
g1data = data.frame(allids,bilength,maxabovelength,numblocks,firstblockabove)
#change bilengths and maxabovelength to units of minutes
g1data$bilength = 10*g1data$bilength
g1data$maxabovelength = 10*g1data$maxabovelength
g1data = filter(g1data,!is.na(xy))
g1data = filter(g1data,bilength!=0)
g1data$replabel = factor(g1data$replabel)



#For each cell we take the mean g1 across all observed buds, including the final partial buds.  Then we average across all cell specific means within each experimental replicate
g1cellsummary <- g1data %>% group_by(cellid,strain,replabel) %>% summarize(meang1 = mean(maxabovelength), meanbi = mean(bilength))
meang1bystrain <- g1cellsummary %>% group_by(strain) %>% summarize(meang1overallcells = mean(meang1))
meang1bystrainnool <- filter(g1cellsummary,meang1<100) %>% group_by(strain) %>% summarize(meang1overallcells = mean(meang1))


repsummary <- g1cellsummary %>% group_by(strain,replabel) %>% summarize(repmeang1 = mean(meang1),repmeanbi = mean(meanbi))
repsummary1 <- repsummary %>% group_by(strain) %>% summarize(meang1 =mean(repmeang1), semg1 = sd(repmeang1)/sqrt(n()),meanbi = mean(repmeanbi),sembi=sd(repmeanbi)/sqrt(n()))
p1 = ggplot(repsummary1,aes(strain,y=meang1)) + geom_bar(stat="identity",width = bwidth,fill=fillc)+ geom_point() + geom_errorbar(aes(ymin=meang1-semg1,ymax=meang1+semg1),width=ebwidth)+ylim(0,30)  +xlab('') + ylab('Mean G1 (min)') + ggtitle('')+themes 
plotfilename = plotfilename = paste(outfolder,'MeanG1_overreplicates.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=3.5,height=4.5)
removeggplottext(p1,plotfilename,4,600,3.5,4.5)

bwidth1 = 0.5
ebwidth1 = 0.25
p1 = ggplot(g1cellsummary,aes(strain,y=meang1)) + geom_boxplot(width=bwidth1,outlier.shape=NA)+ geom_jitter(width=ebwidth1,height = 0,aes(colour=replabel),alpha=0.5) +ylim(0,320)  +xlab('') + ylab('Mean G1 (min)') + ggtitle('')+themes 
plotfilename = plotfilename = paste(outfolder,'MeanG1_boxplots.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=4,height=4.5)
removeggplottext(p1,plotfilename,4,600,4,4.5)

#Now splitting the G1 plot into two different plots with different scale. Will manually overlay one on top of the other

# From the bottom up to the top whisker.
themes1 = theme(axis.text=element_text(size=18), axis.title=element_text(size=15),strip.text.x = element_text(size = 14))
p1 = ggplot(g1cellsummary,aes(strain,y=meang1)) + geom_boxplot(width=bwidth1,outlier.shape=NA)+ geom_jitter(width=ebwidth1,height = 0,aes(colour=replabel),alpha=0.5) +xlab('') + ylab('Mean G1 (min)') + ggtitle('')+themes1 + theme(legend.position='none') + coord_cartesian(ylim = c(0,30)) + scale_y_continuous(breaks = seq(0,30,10),labels = c('0','10','20','30'))
plotfilename = plotfilename = paste(outfolder,'MeanG1_boxplots_ymax50.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=4,height=2.8)
removeggplottext(p1,plotfilename,4,600,4,2.8)


#We filter out the g1duration at 30 since it is already present in the lower part of the graph
panel.border = element_rect(colour="black",linetype="solid",fill=NA,size=1)
themes2 = theme(axis.text=element_text(size=15), axis.title=element_text(size=15),strip.text.x = element_text(size = 14),axis.line.x=element_blank(), axis.ticks.x=element_blank(),axis.text.x=element_blank())
p2 = ggplot(filter(g1cellsummary,meang1>30),aes(strain,y=meang1)) + geom_jitter(width=ebwidth1,height = 0,aes(colour=replabel),alpha=0.5)  +xlab('') + ylab('Mean G1 (min)') + ggtitle('')+themes2  + coord_cartesian(ylim=c(30,340)) + scale_y_continuous(breaks=c(30,100,200,300),expand=c(0,0),labels=c('',100,200,300))
plotfilename = plotfilename = paste(outfolder,'MeanG1_boxplots_yrange40to320.pdf',sep="")
ggsave(file = plotfilename,p2,dpi=600,width=4,height=2.4)
removeggplottext(p2,plotfilename,4,600,4,2.4)




p1 = ggplot(repsummary1,aes(strain,y=meanbi)) + geom_bar(stat="identity",width=bwidth,fill=fillc) + geom_point() + geom_errorbar(aes(ymin=meanbi-sembi,ymax=meanbi+sembi),width=ebwidth)+scale_y_continuous(expand = c(0,0),limits=c(0,140),breaks=seq(0,140,20))  +xlab('') + ylab('Mean Budding Interval (min)') + ggtitle('')+themes  
plotfilename = plotfilename = paste(outfolder,'MeanBuddingInterval_overreplicates.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=3.5,height=4.5)
removeggplottext(p1,plotfilename,4,600,3.5,4.5)


#Now we look at the relation ship between g1 and total budding interval between the two strains. 

#First for smaller values of bilength
g1datanonzero = filter(g1data,maxabovelength>0);
p1 = ggplot(g1datanonzero,aes(bilength,maxabovelength,colour=strain,shape=replabel),alpha=0.01) + geom_jitter() + xlim(0,160) + ylim(0,80) + ylab('G1 duration (min)') + xlab('Budding Interval (min)') + themes+ guides(colour=guide_legend(override.aes = list(shape=15,size=5)),shape = guide_legend(override.aes=list(shape = c(1,2),size=4)))
plotfilename = plotfilename = paste(outfolder,'G1vsBuddingInterval_onlynonzeroG1s_smallerrange.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=5,height=4)
removeggplottextleaveticks(p1,plotfilename,4,600,5,4)


#Now the larger plot
p1 = ggplot(g1datanonzero,aes(bilength,maxabovelength,colour=strain,shape=replabel),alpha=0.01) + geom_jitter() + ylab('G1 duration (min)') + xlab('Budding Interval (min)') + themes+ guides(colour=guide_legend(override.aes = list(shape=15,size=5)),shape = guide_legend(override.aes=list(shape = c(1,2),size=4)))
plotfilename = plotfilename = paste(outfolder,'G1vsBuddingInterval_onlynonzeroG1s_largerrange.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=5,height=4)
removeggplottextleaveticks(p1,plotfilename,4,600,5,4)


#Comparing the cell specific g1s and budding intervals betwen the replicates
replicate = paste(g1cellsummary$strain,g1cellsummary$replabel)
g1cellsummary = data.frame(replicate,g1cellsummary)
p1 = ggplot(g1cellsummary,aes(replicate,meang1)) + geom_boxplot() + ylab('Mean G1 (min)') + xlab('Replicate') + themes 
plotfilename = paste(outfolder,'MeanG1_singlecell_byreplicate.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=5,height=4)
removeggplottext(p1,plotfilename,4,600,5,4)
p1 = ggplot(g1cellsummary,aes(replicate,meanbi)) + geom_boxplot() + ylab('Mean Budding Interval (min)') + xlab('Replicate') + themes 
plotfilename = plotfilename = paste(outfolder,'MeanBI_singlecell_byreplicate.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=5,height=4)
removeggplottext(p1,plotfilename,4,600,5,4)



#Now back checking that our data makes sense. Looking at the distribution of number of blocks, and the order of low and high blocks
numblocks1 = g1data$numblocks
numblocks1[numblocks1>=5] = '>=5'
numblocks1 =factor(numblocks1,levels = c(1,2,3,4,'>=5'))
numblocks1 = unlist(numblocks1)
g1data = data.frame(g1data,numblocks1)
g1summary = g1data %>% group_by(strain,numblocks1) %>% summarize(n = n())
g1summarybystrain = g1summary %>% group_by(strain) %>% mutate(fraction=n/sum(n))

p1 = ggplot(g1summarybystrain,aes(x=numblocks1,fill=strain,y=fraction)) + geom_bar(stat="identity",position = position_dodge(width=0.5),width = 0.5) + ggtitle('numblocks for threshold of 0.00045') + ylab('Fraction') + xlab('Number of blocks')
plotfilename = paste(outfolder,'Numberofblocksinbuddinginterval_distribution.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=5,height=4)
removeggplottextleaveticks(p1,plotfilename,4,600,5,4)

#Saving for Ping to have a look
write.csv(g1data,file='./Figures/g1data_long.csv',row.names=FALSE)

#plotting age at dox
p1 = ggplot(info,aes(strain,ageatfl)) + geom_boxplot()
plotfilename = './Figures/ageatflbystrain.pdf'
ggsave(file = plotfilename,p1,dpi=600,width=5,height=4)


strainrep = paste(g1data$strain,g1data$date)
generation = as.numeric(g1data$bi)+g1data$ageatfl
g1data = data.frame(g1data,generation,strainrep)
write.csv(g1data,file='./Figures/g1data_long.csv',row.names=FALSE)

#Restricting ourselves to budding intervals between generatin 20 and 25
g1datasmall = filter(g1data,generation <= 25 & generation >=20)
g1cellsummary <- g1datasmall %>% group_by(cellid,strain,replabel) %>% summarize(meang1 = mean(maxabovelength), meanbi = mean(bilength))
repsummary <- g1cellsummary %>% group_by(strain,replabel) %>% summarize(repmeang1 = mean(meang1),repmeanbi = mean(meanbi))
repsummary1 <- repsummary %>% group_by(strain) %>% summarize(meang1 =mean(repmeang1), semg1 = sd(repmeang1)/sqrt(n()),meanbi = mean(repmeanbi),sembi=sd(repmeanbi)/sqrt(n()))
p1 = ggplot(repsummary1,aes(strain,y=meang1)) + geom_bar(stat="identity",width = bwidth,fill=fillc)+ geom_point() + geom_errorbar(aes(ymin=meang1-semg1,ymax=meang1+semg1),width=ebwidth)+ylim(0,30)  +xlab('') + ylab('Mean G1 (min)') + ggtitle('')+themes 
plotfilename = plotfilename = paste(outfolder,'MeanG1_overreplicates_gen20-25.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=3.5,height=4.5)
removeggplottext(p1,plotfilename,4,600,3.5,4.5)

p1 = ggplot(repsummary1,aes(strain,y=meanbi)) + geom_bar(stat="identity",width=bwidth,fill=fillc) + geom_point() + geom_errorbar(aes(ymin=meanbi-sembi,ymax=meanbi+sembi),width=ebwidth)+ylim(0,150)  +xlab('') + ylab('Mean Budding Interval (min)') + ggtitle('')+themes  
plotfilename = plotfilename = paste(outfolder,'MeanBuddingInterval_overreplicates_gen20-25.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=3.5,height=4.5)
removeggplottext(p1,plotfilename,4,600,3.5,4.5)




