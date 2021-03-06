setwd('/Users/thomasyoung/Dropbox/MovieProcessing/Whi5Localization_maxips/')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/Preprocessing_func.Rd')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/timeseries_func.Rd')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/func.Rd')
library(dplyr)
library(cowplot)
library(reshape2)
library(grid)
library(gridExtra)
library(grid)

#original ebwidth is 0.25
#figure settings:
themes = theme(axis.text=element_text(size=14), axis.title=element_text(size=16),strip.text.x = element_text(size = 14))
themes = theme(axis.text=element_text(size=17), axis.title=element_text(size=20),strip.text.x = element_text(size = 20))
themes1 = theme(axis.text=element_text(size=18), axis.title=element_text(size=16),strip.text.x = element_text(size = 16,margin = margin(.3,0,0.3,0,"cm")),strip.background=element_rect(fill=NA))
ylims = ylim(-0,25)
jitterw = 0.12
ebwidth = 0.5;
ylabel = ylab("Age (generations)");
xscale = scale_x_discrete(expand = c(0.6,0))
border = theme(panel.border=element_rect(colour="black",linetype="solid",fill=NA,size=1))



#folder to save in
outfolder = './Figures/ageatdox/'


info = read.csv('./CombinedData/infoall2.csv')
info = filter(info,ageatdox>=15)

#Plotting mean+/-sd of the age at dox by experiment  
stats <- info %>% group_by(strain,doxtime,replabel) %>% summarize(meanageatdox = mean(ageatdox),total = n())
tograph <- stats %>% group_by(strain,doxtime) %>% summarize(meandoxage = mean(meanageatdox),sem = sd(meanageatdox)/sqrt(n()))
write.csv(tograph,paste(outfolder,'ageatdox_meanandsemofreplicateaverages.csv',sep=""))









#Jitter plots of cell ages at dox with replicate based mean +/-sem overlaid
#Now plotting the raw age at dox data pooled across replicates sharing the same strain and age
stats <- info %>% group_by(strain,doxtime,replabel,expstrain) %>% summarize(meanageatdox = mean(ageatdox),total = n())
tograph <- stats %>% group_by(strain,doxtime,expstrain) %>% summarize(ageatdox= mean(meanageatdox),sem = sd(meanageatdox)/sqrt(n()))
outtable = paste(outfolder,'mean+-sem_ageatdox.csv',sep="")
write.csv(tograph,outtable,row.names=FALSE)


#SSA (old cells)
xlabel = xlab('')
p1 = ggplot(info,aes(strain,ageatdox)) + geom_jitter(width=jitterw,height=0,colour='red',size=1)+ xlab('strain') + ylims + ggtitle('')+themes+geom_point(data = tograph,aes(strain,ageatdox)) + geom_errorbar(data = tograph,aes(ymin = ageatdox-sem,ymax = ageatdox+sem),width=ebwidth) + ylabel + xscale + border
plotfilename = paste(outfolder,'rawageatdox_meansem_oldSSA_mpaornot.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=3.5,height=4.5)
removeggplottext(p1,plotfilename,4,600,3.5,4.5)

p1 <- ggplot(info,aes(birth,ageatdox,colour=strain)) + geom_point()
p1 <- ggplot(info,aes(paste(strain,date),birth)) + geom_boxplot()



