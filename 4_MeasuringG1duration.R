#Measuring G1 durations, comparing between strains yTY159b and yTY160a. Both strains have Whi5 tagged to mCitrine, and SSAcontrol (non-cuttable) cassettes. The difference is that the yTY159b cassette as RFP while the yTY160a cassete has RFPdegron.
#The threshold of 0.0000425 is applied to Whi5-YFP localization scores for each cell budding interval (see 4_AssesingWhi5Thresholds.R)
#The number of localization scores above the threshold in a given budding interval is the estimated G1 duration. This only applies to cells for which there is there is 
#Plot the mean G1 durations single cells each of the two strains. Replicate averages of the single-cell averages, and single-cell averages are plotted. Results are saved in './Figures/g1durationcomparison/'.  

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

#maxminusmean is the Whi5-YFP localization score
flid = fl$cellid
fl = data.frame(flid,fl)
fl = mutate(fl,yfp = maxminusmean)
flwide <- dcast(fl,flid~time,value.var = 'yfp')

#Matching the single cell localization scores to other annotations
infoid = paste(info$date,info$xy,info$trap)
flid = as.character(flwide$flid)
matches = match(infoid,flid)
matchesinfoindex = which(!is.na(matches))
matchesflindex = matches[!is.na(matches)]

#Getting budding times. These will be used to determine the budding intervals for which to measure G1 durations
bt = getbudtimes(info)
ageatfl = ageattimes(bt,rep(155,nrow(bt)))
info = data.frame(ageatfl,info)
bt = bt-154;

#Matching the single cell localization scores for each cell to budding time measurements for the same cell. This generates this list of budding times for each cell (btlist), and the list of localization scores (flwidelist). The same index in both lists corresponds to the same cell. 
flwide = flwide[2:ncol(flwide)]
flwidelist = mattolol(flwide)
btlist = mattolol(bt)
btlist = btlist[matchesinfoindex]
flwidelist = flwidelist[matchesflindex]
infosmall = info[matchesflindex,]

#Saving cells counts per replicate
cellrepcount = table(infosmall$strain,infosmall$date)
write.csv(cellrepcount,'./Figures/NumberOfCellsPerRep_G1comparison.csv')


#Apply the localization threshold. For each budding interval record the pattern of measurement blocks above or below the threshold.
#The pattern of blocks is described by an 'rle' or 'run length encoding' that stores the lengths of continguous measurements above or below the cutoff in order.

#applying the cutoff
cutoff = 0.0000425
coapplied = lapply(flwidelist,greaterthan,cutoff)
names(coapplied) = flid[matchesflindex]

#Creating a list of the 'above cutoff' measurements for each budding interval. Each entry in the list corresponds to the budding interval of a given cell. The entry is itself a list of booleans corresponding to whether each localization score in the budding interval is above the cutoff
biwindows = mapply(extractwindows1, coapplied, btlist)
biwindows = unlist(biwindows,recursive=FALSE)
biwindows = lapply(biwindows,f<-function(x){return(x[!is.na(x)])})

#Getting the length, rle, number of blocks above and below for each entry in the list.
#maxabovelength is the longest period of time the Whi5 localization score is above the cutoff. maxabovelength is used as the metric for G1 duration.
bilengths = lapply(biwindows,length)
rles = lapply(biwindows,rle)
rlelengths = lapply(rles, getrlelengths)
rlevalues = lapply(rles, getrlevalues)
numblocks = lapply(rlelengths,length)
numblocksabove = lapply(rlevalues,f<-function(x){return(sum(x))})
numblocksbelow = lapply(rlevalues,f<-function(x){return(sum(!x))})
abovelengths = mapply(gettrurlelengths,rlevalues,rlelengths)
maxabovelength = lapply(abovelengths,max,na.rm=TRUE)

#Example: 3/12/19 26 6 has only 2 measured budding intervals. when the is a 3rd
#corresponds to entry 377 and 378 out of 616 entries.  The 3rd budding interval
#is the last, incomplete one (since we didn't measure fl for all time points in the budding interval). Compare with the plot in './Figures/FlTrajectories/yfp_maxminusmean5x5yTY160a
rlelengths[387:389]
biwindows[387:389]


#Checking how sensible the G1 duration measurements are.
#	1. The first block value corresponds to the appearance of a bud in the phase g image of the cell. By this time Whi5 should no longer be nuclear localized since such cells are in S phase. Most of these block values should be False (localization score < threshold). This is true.
#	2. The number of blocks should ideally be 2 or 3. This corresponds to Whi5-YFP being (below,above,below) or (below, above) the threshold in each budding interval. Actually, the data shows that most budding intervals have no blocks above the threshold (430/627 = 69%). There are ((14+143)/627 = 0.25). 
firstblockvalue = lapply(rlevalues,f<-function(x){return(x[1])})
table(unlist(firstblockvalue))
table(unlist(numblocks))



#extracting the id information for each budding interval to match to the strain and replicate information.  Once matched, see how budding intervals are distributed among cells by strain
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


#Making a dataframe with the relevant cell and budding interval id information and the extracted G1 data. maxabovelength (the longest run of Whi5 localization scores > cutoff) is the estimated G1 duration
bilength = unlist(bilengths)
maxabovelength = unlist(maxabovelength)
numblocks = unlist(numblocks)
firstblockabove = unlist(firstblockvalue)
g1data = data.frame(allids,bilength,maxabovelength,numblocks,firstblockabove)

#change bilengths and maxabovelength to units of minutes
g1data$bilength = 10*g1data$bilength
g1data$maxabovelength = 10*g1data$maxabovelength

#Removing cells for which there is the xy value is na, or the length of the budding interval is 0
g1data = filter(g1data,!is.na(xy))
g1data = filter(g1data,bilength!=0)
g1data$replabel = factor(g1data$replabel)



#For each cell take the mean g1 across all observed buds, including the final partial buds.  Then average across all cell specific means within each experimental replicate
g1cellsummary <- g1data %>% group_by(cellid,strain,replabel) %>% summarize(meang1 = mean(maxabovelength), meanbi = mean(bilength))
meang1bystrain <- g1cellsummary %>% group_by(strain) %>% summarize(meang1overallcells = mean(meang1))

#Table of g1 durations for cells where g1 duration was less than 100
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
# From the bottom up to the top whisker. This is will be the lower part of the boxplot figure
themes1 = theme(axis.text=element_text(size=18), axis.title=element_text(size=15),strip.text.x = element_text(size = 14))
p1 = ggplot(g1cellsummary,aes(strain,y=meang1)) + geom_boxplot(width=bwidth1,outlier.shape=NA)+ geom_jitter(width=ebwidth1,height = 0,aes(colour=replabel),alpha=0.5) +xlab('') + ylab('Mean G1 (min)') + ggtitle('')+themes1 + theme(legend.position='none') + coord_cartesian(ylim = c(0,30)) + scale_y_continuous(breaks = seq(0,30,10),labels = c('0','10','20','30'))
plotfilename = plotfilename = paste(outfolder,'MeanG1_boxplots_ymax50.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=4,height=2.8)
removeggplottext(p1,plotfilename,4,600,4,2.8)

#Making the upper part of the G1 duration boxplot:
#Filter out the g1duration at 30 since it is already present in the lower part of the graph
panel.border = element_rect(colour="black",linetype="solid",fill=NA,size=1)
themes2 = theme(axis.text=element_text(size=15), axis.title=element_text(size=15),strip.text.x = element_text(size = 14),axis.line.x=element_blank(), axis.ticks.x=element_blank(),axis.text.x=element_blank())
p2 = ggplot(filter(g1cellsummary,meang1>30),aes(strain,y=meang1)) + geom_jitter(width=ebwidth1,height = 0,aes(colour=replabel),alpha=0.5)  +xlab('') + ylab('Mean G1 (min)') + ggtitle('')+themes2  + coord_cartesian(ylim=c(30,340)) + scale_y_continuous(breaks=c(30,100,200,300),expand=c(0,0),labels=c('',100,200,300))
plotfilename = plotfilename = paste(outfolder,'MeanG1_boxplots_yrange40to320.pdf',sep="")
ggsave(file = plotfilename,p2,dpi=600,width=4,height=2.4)
removeggplottext(p2,plotfilename,4,600,4,2.4)



#Plotting mean budding interval for each strain over the replicates
p1 = ggplot(repsummary1,aes(strain,y=meanbi)) + geom_bar(stat="identity",width=bwidth,fill=fillc) + geom_point() + geom_errorbar(aes(ymin=meanbi-sembi,ymax=meanbi+sembi),width=ebwidth)+scale_y_continuous(expand = c(0,0),limits=c(0,140),breaks=seq(0,140,20))  +xlab('') + ylab('Mean Budding Interval (min)') + ggtitle('')+themes  
plotfilename = plotfilename = paste(outfolder,'MeanBuddingInterval_overreplicates.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=3.5,height=4.5)
removeggplottext(p1,plotfilename,4,600,3.5,4.5)


#Now look at the relationship between g1 and total budding interval between the two strains. 
#Plot g1 duration vs budding interval length for each budding interval observed
#Only considering budding intervals less than 160 and g1 durations less than 80
g1datanonzero = filter(g1data,maxabovelength>0);
p1 = ggplot(g1datanonzero,aes(bilength,maxabovelength,colour=strain,shape=replabel),alpha=0.01) + geom_jitter() + xlim(0,160) + ylim(0,80) + ylab('G1 duration (min)') + xlab('Budding Interval (min)') + themes+ guides(colour=guide_legend(override.aes = list(shape=15,size=5)),shape = guide_legend(override.aes=list(shape = c(1,2),size=4)))
plotfilename = plotfilename = paste(outfolder,'G1vsBuddingInterval_onlynonzeroG1s_smallerrange.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=5,height=4)
removeggplottextleaveticks(p1,plotfilename,4,600,5,4)


#Plotting g1 duration vs budding interval length for each budding interval observed.
#Plotting the full range of data, though it may be harder to see the points
p1 = ggplot(g1datanonzero,aes(bilength,maxabovelength,colour=strain,shape=replabel),alpha=0.01) + geom_jitter() + ylab('G1 duration (min)') + xlab('Budding Interval (min)') + themes+ guides(colour=guide_legend(override.aes = list(shape=15,size=5)),shape = guide_legend(override.aes=list(shape = c(1,2),size=4)))
plotfilename = plotfilename = paste(outfolder,'G1vsBuddingInterval_onlynonzeroG1s_largerrange.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=5,height=4)
removeggplottextleaveticks(p1,plotfilename,4,600,5,4)


#Comparing the cell specific g1s and budding intervals between the replicates
#One boxplot for each replicate.
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



#Plotting the fraction of budding intervals that have 1,2,3,4,>=5 blocks relative to the Whi5-YFP localization cutoffnumblocks1 = g1data$numblocks
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


#plotting age at dox (boxplot)
p1 = ggplot(info,aes(strain,ageatfl)) + geom_boxplot()
plotfilename = './Figures/ageatflbystrain.pdf'
ggsave(file = plotfilename,p1,dpi=600,width=5,height=4)

#Saving a table of relevant data for each budding interval observed
strainrep = paste(g1data$strain,g1data$date)
generation = as.numeric(g1data$bi)+g1data$ageatfl
g1data = data.frame(g1data,generation,strainrep)
write.csv(g1data,file='./Figures/g1data_long.csv',row.names=FALSE)

#Looking at g1 duration, only for budding intervals between generation 20 and 25
g1datasmall = filter(g1data,generation <= 25 & generation >=20)
g1cellsummary <- g1datasmall %>% group_by(cellid,strain,replabel) %>% summarize(meang1 = mean(maxabovelength), meanbi = mean(bilength))
repsummary <- g1cellsummary %>% group_by(strain,replabel) %>% summarize(repmeang1 = mean(meang1),repmeanbi = mean(meanbi))
repsummary1 <- repsummary %>% group_by(strain) %>% summarize(meang1 =mean(repmeang1), semg1 = sd(repmeang1)/sqrt(n()),meanbi = mean(repmeanbi),sembi=sd(repmeanbi)/sqrt(n()))
p1 = ggplot(repsummary1,aes(strain,y=meang1)) + geom_bar(stat="identity",width = bwidth,fill=fillc)+ geom_point() + geom_errorbar(aes(ymin=meang1-semg1,ymax=meang1+semg1),width=ebwidth)+ylim(0,30)  +xlab('') + ylab('Mean G1 (min)') + ggtitle('')+themes 
plotfilename = plotfilename = paste(outfolder,'MeanG1_overreplicates_gen20-25.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=3.5,height=4.5)
removeggplottext(p1,plotfilename,4,600,3.5,4.5)

#Looking at budding interval length only for budding intervals between generation 20 and 25
p1 = ggplot(repsummary1,aes(strain,y=meanbi)) + geom_bar(stat="identity",width=bwidth,fill=fillc) + geom_point() + geom_errorbar(aes(ymin=meanbi-sembi,ymax=meanbi+sembi),width=ebwidth)+ylim(0,150)  +xlab('') + ylab('Mean Budding Interval (min)') + ggtitle('')+themes  
plotfilename = plotfilename = paste(outfolder,'MeanBuddingInterval_overreplicates_gen20-25.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=3.5,height=4.5)
removeggplottext(p1,plotfilename,4,600,3.5,4.5)




