#Load libraries
library(tidyverse)
library(ggforce)
library(here)
conflict_prefer("here", "here")
library(cowplot)

#Load DDA Posaconazole evolved data - contains data from the replicates evolved in 72h transfers in POS for "FLC" "CTR" "NYT" "MCZ" "5FC" "VCZ" "POS"; all strains in POS & FLC, other drugs on SC5314 evolved replicates

DDA <- read_csv(here("Data_In", "DDA", "POSevolved", "POS", "211218_POSevolved_allDDA.csv"))
DDA$line <-ifelse(DDA$person =="M", DDA$replicate, DDA$replicate+12)

#dark2 colours
Dark2  <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")

################################
#POS DDA
################################
DDA_POS <- subset(DDA, drug=="POS")
DDA_POS_anc <- subset(DDA_POS, time =="0")

DDA_POStemp <-  DDA_POS %>%
  #  filter(line %in% Evol) %>%
  group_by(bioreplicate, strain, line,time) %>%
  summarise(nbiorep = n())

DDA_POSm <-  DDA_POS %>%
  group_by(strain, line, time) %>%
  summarise(avgRAD20 = mean(RAD20, na.rm=TRUE), avgFoG20 = mean(FoG20, na.rm=TRUE))

DDA_POSm$ST <- paste(DDA_POSm$strain, DDA_POSm$time, sep="_")

#POS Plot with all ancestral and evolved replicates side by side
RAD20 <- ggplot(data=DDA_POSm, mapping=aes(x = strain, y=avgRAD20, col=as.factor(time)))+
  geom_sina(alpha=0.6,size=2.5, position = position_dodge(width=0.4))+
#  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.5) +
  xlab("Strain")+
  ylab("Susceptibility (RAD20)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title = element_text(size=12), plot.margin = unit(c(0.2, 0.2, 0, 0.2), "cm"))+
  ylim(23,0)+
  scale_color_manual(values = c("grey","darkmagenta"),name="",labels=c("Ancestral","Evolved"))
  #annotate("text",x=c(2,3,7,8),y=c(9.5,8.5,11,9.5),label = "*")
#annotate("segment",x=c(0.75,3.75,4.75,5.75),xend = c(1.25,4.25,5.25,6.25),y = c(10,9,11.5,10), yend =c(10,9,11.5,10))

FoG20 <- ggplot(data=DDA_POSm, mapping=aes(x = strain, y=avgFoG20, col=as.factor(time)))+
  geom_sina(alpha=0.6,size=2.5, position = position_dodge(width=0.4))+
  #  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.5) +
  xlab("Strain")+
  ylab("Tolerance (FoG20)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=12), axis.title = element_text(size=12), plot.margin = unit(c(0, 0.2, 0.2, 0.2), "cm"))+
  scale_color_manual(values = c("grey","darkmagenta"),name="",labels=c("Ancestral","Evolved"))+
  ylim(0,1)
# annotate("text",x=c(1,3,4,5,6,7),y=c(0.67,0.62,0.7,0.64,0.71,0.51),label = "*")
#annotate("segment",x=c(0.75,2.75,3.75,4.75,5.75,6.75),xend = c(1.25,3.25,4.25,5.25,6.25,7.25),y = c(0.66,0.61,0.69,0.63,0.7,0.5), yend =c(0.66,0.61,0.69,0.63,0.7,0.5))

DDA20 <- plot_grid(RAD20,FoG20,nrow = 2,align = "hv")
DDA20
ggsave(here("Figures", "PDF", "POS_72_POS-DDA_20.pdf"), DDA20, width = 5, height=5)

DDA_POSgm <-  DDA_POS %>%
  group_by(strain, time) %>%
  summarise(avgRAD20 = median(RAD20, na.rm=TRUE), avgFoG20 = median(FoG20, na.rm=TRUE))

DDA_POSgm_anc <- subset(DDA_POSgm, time == "0")
DDA_POSm_anc <- subset(DDA_POSm, time == "0")
DDA_POSm_evol <- subset(DDA_POSm, time == "5")

radDiff <- c()
FoGDiff <- c()
strain <- c()
RADresults <- data.frame()
FoGresults <- data.frame()
j <- 0
for (i in unique(DDA_POSm_evol$strain)){
  print(i)
 evol <- subset(DDA_POSm_evol, strain==i)
 anc <- subset(DDA_POSm_anc, strain==i)
 anc_mean <- subset(DDA_POSgm_anc, strain==i)

 radDiff <- append(radDiff, evol$avgRAD20 - anc_mean$avgRAD20)
 FoGDiff <- append(FoGDiff, evol$avgFoG20 - anc_mean$avgFoG20)
 strain <- append(strain, evol$strain)

 if(length(evol$avgRAD20) > 2){
 t <- t.test(anc$avgRAD20, evol$avgRAD20)
 results_RAD <- c(i, round(t$estimate[2]-t$estimate[1],2), round(t$statistic,2), round(t$parameter,2), round(t$p.value,4))
 RADresults <- rbind(RADresults, results_RAD)

 t2 <- t.test(anc$avgFoG20, evol$avgFoG20)
 results_FoG <- c(i, round(t2$estimate[2]-t2$estimate[1],2), round(t2$statistic,2), round(t2$parameter,2), round(t2$p.value,4))
 FoGresults <- rbind(FoGresults, results_FoG)
 }
 else{
   RADresults <- rbind(RADresults, c(i, NA, NA, NA, NA))
   FoGresults <- rbind(FoGresults, c(i, NA, NA, NA, NA))
 }
}
colnames(RADresults) <- paste(c("Strain","Evol - Anc","t-statistic","Degrees of freedom","p-value"))
colnames(FoGresults) <- paste(c("Strain","Evol - Anc","t-statistic","Degrees of freedom","p-value"))

write.csv(RADresults, file = here("Data_out","DDA", "POS_72_POS_RADstatistics.csv"))
write.csv(FoGresults, file = here("Data_out", "DDA", "POS_72_POS_FoGstatistics.csv"))

diffs <- data.frame(strain, radDiff, FoGDiff)

DDA_POSgm <-  DDA_POS %>%
  group_by(strain, time) %>%
  summarise(avgRAD20 = median(RAD20, na.rm=TRUE), avgFoG20 = median(FoG20, na.rm=TRUE))

DDA_POSgm_anc <- subset(DDA_POSgm, time == "0")
DDA_POSgm_evol <- subset(DDA_POSgm, time == "5")

DDA_POS_diff <- data.frame(strain=character(), RAD20diff = numeric(), FoG20diff = numeric())
for (i in unique(DDA_POSm_evol$strain)){
  evolTemp <- subset(DDA_POSm_evol, strain==i)
  ancTemp <- subset(DDA_POSgm_anc, strain==i)
  tempdf <- data.frame(strain =  evolTemp$strain, line = evolTemp$line, RAD20diff = evolTemp$avgRAD20-ancTemp$avgRAD20, FoG20diff = evolTemp$avgFoG20 - ancTemp$avgFoG20)
  DDA_POS_diff <- rbind(DDA_POS_diff, tempdf)
}

plot(as.numeric(as.factor(DDA_POS_diff$strain)), DDA_POS_diff$RAD20diff, ylim=c(-10, 10))
abline(h=0)

plot(as.numeric(as.factor(DDA_POS_diff$strain)), DDA_POS_diff$FoG20diff, ylim=c(-1, 1))
abline(h=0)

######################################
#FLC
######################################
DDA_FLC <- subset(DDA, drug=="FLC")
DDA_FLC_anc <- subset(DDA_FLC, time =="0")

DDA_FLCtemp <-  DDA_FLC %>%
  #  filter(line %in% Evol) %>%
  group_by(bioreplicate, strain, line,time) %>%
  summarise(nbiorep = n())

DDA_FLCm <-  DDA_FLC %>%
  group_by(strain, line, time) %>%
  summarise(avgRAD20 = mean(RAD20, na.rm=TRUE), avgFoG20 = mean(FoG20, na.rm=TRUE))

DDA_FLCm$ST <- paste(DDA_FLCm$strain, DDA_FLCm$time, sep="_")

RAD20_FLC <- ggplot(data=DDA_FLCm, mapping=aes(x = strain, y=avgRAD20, col=as.factor(time)))+
  geom_sina(alpha=0.6,size=2.5, position = position_dodge(width=0.4))+
  #  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.5) +
  xlab("Strain")+
  ylab("Susceptibility (RAD20)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title = element_text(size=12), plot.margin = unit(c(0.2, 0.2, 0, 0.2), "cm"))+
  ylim(23,0)+
  scale_color_manual(values = c("grey","darkmagenta"),name="",labels=c("Ancestral","Evolved"))
#annotate("text",x=c(2,3,7,8),y=c(9.5,8.5,11,9.5),label = "*")
#annotate("segment",x=c(0.75,3.75,4.75,5.75),xend = c(1.25,4.25,5.25,6.25),y = c(10,9,11.5,10), yend =c(10,9,11.5,10))

FoG20_FLC <- ggplot(data=DDA_FLCm, mapping=aes(x = strain, y=avgFoG20, col=as.factor(time)))+
  geom_sina(alpha=0.6,size=2.5, position = position_dodge(width=0.4))+
  #  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.5) +
  xlab("Strain")+
  ylab("Tolerance (FoG20)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=12), axis.title = element_text(size=12), plot.margin = unit(c(0, 0.2, 0.2, 0.2), "cm"))+
  scale_color_manual(values = c("grey","darkmagenta"),name="",labels=c("Ancestral","Evolved"))+
  ylim(0,1)
# annotate("text",x=c(1,3,4,5,6,7),y=c(0.67,0.62,0.7,0.64,0.71,0.51),label = "*")
#annotate("segment",x=c(0.75,2.75,3.75,4.75,5.75,6.75),xend = c(1.25,3.25,4.25,5.25,6.25,7.25),y = c(0.66,0.61,0.69,0.63,0.7,0.5), yend =c(0.66,0.61,0.69,0.63,0.7,0.5))

DDA20_FLC <- plot_grid(RAD20_FLC,FoG20_FLC,nrow = 2,align = "hv")
DDA20_FLC
ggsave(here("Figures", "PDF", "FLC_72_FLC-DDA_20.pdf"), DDA20, width = 5, height=5)

DDA_FLCgm <-  DDA_FLC %>%
  group_by(strain, time) %>%
  summarise(avgRAD20 = mean(RAD20, na.rm=TRUE), avgFoG20 = mean(FoG20, na.rm=TRUE))

DDA_FLCgm_anc <- subset(DDA_FLCgm, time == "0")
DDA_FLCm_anc <- subset(DDA_FLCm, time == "0")
DDA_FLCm_evol <- subset(DDA_FLCm, time == "5")

radDiff_FLC <- c()
FoGDiff_FLC <- c()
strain_FLC <- c()
RADresults_FLC <- data.frame()
FoGresults_FLC <- data.frame()
j <- 0
for (i in unique(DDA_FLCm_evol$strain)){
  evol <- subset(DDA_FLCm_evol, strain==i)
  anc <- subset(DDA_FLCm_anc, strain==i)
  anc_mean <- subset(DDA_FLCgm_anc, strain==i)

  radDiff_FLC <- append(radDiff_FLC, evol$avgRAD20 - anc_mean$avgRAD20)
  FoGDiff_FLC <- append(FoGDiff_FLC, evol$avgFoG20 - anc_mean$avgFoG20)
  strain_FLC <- append(strain_FLC, evol$strain)

  if(length(evol$avgFoG20) > 2){
    if(sum(!is.na(evol$avgFoG20)) > 2){
      t <- t.test(anc$avgRAD20, evol$avgRAD20)
      results_RAD_FLC <- c(i, round(t$estimate[2]-t$estimate[1],2), round(t$statistic,2), round(t$parameter,2), round(t$p.value,4))

      t2 <- t.test(anc$avgFoG20, evol$avgFoG20)
      results_FoG_FLC <- c(i, round(t2$estimate[2]-t2$estimate[1],2), round(t2$statistic,2), round(t2$parameter,2), round(t2$p.value,4))
    }
    else{
      results_RAD_FLC <- c(i,round(t$estimate[2]-t$estimate[1],2), NA, NA, NA)
      results_FoG_FLC <- c(i,round(t2$estimate[2]-t2$estimate[1],2), NA, NA, NA)
    }
  }
  else{
    results_RAD_FLC <- c(i,round(t$estimate[2]-t$estimate[1],2), NA, NA, NA)
    results_FoG_FLC <- c(i,round(t2$estimate[2]-t2$estimate[1],2), NA, NA, NA)
  }
    print(i)
    RADresults_FLC <- rbind(RADresults_FLC, results_RAD_FLC)
    FoGresults_FLC <- rbind(FoGresults_FLC, results_FoG_FLC)
}
colnames(RADresults_FLC) <- paste(c("Strain","Evol - Anc","t-statistic","Degrees of freedom","p-value"))
colnames(FoGresults_FLC) <- paste(c("Strain","Evol - Anc","t-statistic","Degrees of freedom","p-value"))

write.csv(RADresults_FLC, file = here("Data_out", "DDA", "FLC_72_FLC_RADstatistics.csv"))
write.csv(FoGresults_FLC, file = here("Data_out", "DDA", "FLC_72_FLC_FoGstatistics.csv"))


DDA_FLCgm <-  DDA_FLC %>%
  group_by(strain, time) %>%
  summarise(avgRAD20 = median(RAD20, na.rm=TRUE), avgFoG20 = median(FoG20, na.rm=TRUE))

DDA_FLCgm_anc <- subset(DDA_FLCgm, time == "0")

DDA_FLC_diff <- data.frame(strain=character(), RAD20diff = numeric(), FoG20diff = numeric())
for (i in unique(DDA_FLCm_evol$strain)){
  evolTemp <- subset(DDA_FLCm_evol, strain==i)
  ancTemp <- subset(DDA_FLCgm_anc, strain==i)
  tempdf <- data.frame(strain =  evolTemp$strain, line = evolTemp$line, RAD20diff = evolTemp$avgRAD20-ancTemp$avgRAD20, FoG20diff = evolTemp$avgFoG20 - ancTemp$avgFoG20)
  DDA_FLC_diff <- rbind(DDA_FLC_diff, tempdf)
}

##############################
#POS & FLC Differences combined
##############################
#run 211215-YPD-DDAanalysis.R first to get YPD data

order <- c("A08", "A02",  "A03", "A17", "A04",  "A18", "A10", "A12")
realnames <- c("P75016", "P87",  "GC75",  "SC5314", "P78048", "FH1", "P76055", "T101")

DDA_POS_diff$strain <- factor(DDA_POS_diff$strain, levels=order)
DDA_FLC_diff$strain <- factor(DDA_FLC_diff$strain, levels=order)

pdf("Figures/Figure3-DDA-YPD-POS-FLC.pdf", width=7.5, height=5)
par(mfrow=c(2, 1),mar=c(1,1,1,1), oma=c(3, 3.5, 1, 1))
plot(jitter(as.numeric(as.factor(DDA_POS_diff$strain)), factor=0.2), DDA_POS_diff$RAD20diff, ylim=c(10, -10), xlim=c(0.8,8.2), col="darkblue", ylab = "", xaxt="n", yaxt="n")
points(jitter(as.numeric(as.factor(DDA_FLC_diff$strain))+0.2, factor=0.2), DDA_FLC_diff$RAD20diff, col=Dark2[6])
points(jitter(as.numeric(as.factor(diffs72_YPD$strain))-0.2, factor=0.2), diffs72_YPD$radDiff72, col="forestgreen")

abline(h=0, lty=2)
axis(1, at=1:8, labels=FALSE)
axis(2, las=2)
txt <- expression(paste(Delta," Susceptibility (", RAD[20], ")"))
mtext(txt, side=2, line=3)

plot(jitter(as.numeric(as.factor(DDA_POS_diff$strain)), 0.2), DDA_POS_diff$FoG20diff, ylim=c(-1, 1), col="darkblue", xlim=c(0.8,8.2),ylab = "Change in Tolerance (FoG20)", xaxt="n", yaxt="n")
points(jitter(as.numeric(as.factor(DDA_FLC_diff$strain))+0.2, factor= 0.2), DDA_FLC_diff$FoG20diff, col=Dark2[6])
points(jitter(as.numeric(as.factor(diffs72_YPD$strain))-0.2, factor= 0.2), diffs72_YPD$FoGDiff72, col="forestgreen")
abline(h=0, lty=2)
axis(1, at=1:8, labels=realnames)
axis(2, las=2)
txt2 <- expression(paste(Delta," Tolerance (", FoG[20], ")"))
mtext(txt2, side=2, line=3)
legend("bottomright", pch=21, col=c("forestgreen", Dark2[8], Dark2[6]), legend=c("YPD-evolved, POS", "POS-evoled, POS", "POS-evolved, FLC"), cex=0.8)
dev.off()


###################
#SC5314
##################
DDA_POS <- subset(DDA, drug=="POS")
DDA_POSm_anc <- subset(DDA_POS, time =="0")
DDA_POSm <-  DDA_POS %>%
  group_by(bioreplicate, strain, line,time) %>%
  summarise(nbiorep = n(), avgRAD20 = mean(RAD20, na.rm=TRUE),avgFoG20 = mean(FoG20, na.rm=TRUE), avgRAD50 = mean(RAD50, na.rm=TRUE),avgFoG50 = mean(FoG50, na.rm=TRUE), deltaRAD50 = avgRAD50 - median(DDA_POSm_anc$RAD50, na.rm=TRUE), deltaFoG50 = avgFoG50 - median(DDA_POSm_anc$FoG50, na.rm=TRUE), deltaRAD20 = avgRAD20 - median(DDA_POSm_anc$RAD20, na.rm=TRUE), deltaFoG20 = avgFoG20 - median(DDA_POSm_anc$FoG20, na.rm=TRUE))
DDA_POSm_evol <- subset(DDA_POSm, time == "5")

DDA_FLC <- subset(DDA, drug=="FLC")
DDA_FLCm_anc <- subset(DDA_FLC, time =="0")
DDA_FLCm <-  DDA_FLC %>%
  group_by(bioreplicate, strain, line,time) %>%
  summarise(nbiorep = n(), avgRAD20 = mean(RAD20, na.rm=TRUE),avgFoG20 = mean(FoG20, na.rm=TRUE), avgRAD50 = mean(RAD50, na.rm=TRUE),avgFoG50 = mean(FoG50, na.rm=TRUE), deltaRAD50 = avgRAD50 - median(DDA_FLCm_anc$RAD50, na.rm=TRUE), deltaFoG50 = avgFoG50 - median(DDA_FLCm_anc$FoG50, na.rm=TRUE), deltaRAD20 = avgRAD20 - median(DDA_FLCm_anc$RAD20, na.rm=TRUE), deltaFoG20 = avgFoG20 - median(DDA_FLCm_anc$FoG20, na.rm=TRUE))
DDA_FLCm_evol <- subset(DDA_FLCm, time == "5")

DDA_CTR <- subset(DDA, drug=="CTR")
DDA_CTRm_anc <- subset(DDA_CTR, time =="0")
DDA_CTRm <-  DDA_CTR %>%
  group_by(bioreplicate, strain, line,time) %>%
  summarise(avgRAD20 = mean(RAD20, na.rm=TRUE),avgFoG20 = mean(FoG20, na.rm=TRUE), avgRAD50 = mean(RAD50, na.rm=TRUE),avgFoG50 = mean(FoG50, na.rm=TRUE), deltaRAD50 = avgRAD50 - median(DDA_CTRm_anc$RAD50, na.rm=TRUE), deltaFoG50 = avgFoG50 - median(DDA_CTRm_anc$FoG50, na.rm=TRUE), deltaRAD20 = avgRAD20 - median(DDA_CTRm_anc$RAD20, na.rm=TRUE), deltaFoG20 = avgFoG20 - median(DDA_CTRm_anc$FoG20, na.rm=TRUE))
DDA_CTRm_evol <- subset(DDA_CTRm, time == "5")

DDA_NYT <- subset(DDA, drug=="NYT")
DDA_NYTm_anc <- subset(DDA_NYT, time =="0")
DDA_NYTm <-  DDA_NYT %>%
  group_by(strain, line,time) %>%
  summarise(avgRAD20 = mean(RAD20, na.rm=TRUE),avgFoG20 = mean(FoG20, na.rm=TRUE), avgRAD50 = mean(RAD50, na.rm=TRUE),avgFoG50 = mean(FoG50, na.rm=TRUE), deltaRAD50 = avgRAD50 - median(DDA_NYTm_anc$RAD50, na.rm=TRUE), deltaFoG50 = avgFoG50 - median(DDA_NYTm_anc$FoG50, na.rm=TRUE), deltaRAD20 = avgRAD20 - median(DDA_NYTm_anc$RAD20, na.rm=TRUE), deltaFoG20 = avgFoG20 - median(DDA_NYTm_anc$FoG20, na.rm=TRUE))
DDA_NYTm_evol <- subset(DDA_NYTm, time == "5")

DDA_MCZ <- subset(DDA, drug=="MCZ")
DDA_MCZm_anc <- subset(DDA_MCZ, time =="0")
DDA_MCZm <-  DDA_MCZ %>%
  group_by(strain, line,time) %>%
  summarise(avgRAD20 = mean(RAD20, na.rm=TRUE),avgFoG20 = mean(FoG20, na.rm=TRUE), avgRAD50 = mean(RAD50, na.rm=TRUE),avgFoG50 = mean(FoG50, na.rm=TRUE), deltaRAD50 = avgRAD50 - median(DDA_MCZm_anc$RAD50, na.rm=TRUE), deltaFoG50 = avgFoG50 - median(DDA_MCZm_anc$FoG50, na.rm=TRUE), deltaRAD20 = avgRAD20 - median(DDA_MCZm_anc$RAD20, na.rm=TRUE), deltaFoG20 = avgFoG20 - median(DDA_MCZm_anc$FoG20, na.rm=TRUE))
DDA_MCZm_evol <- subset(DDA_MCZm, time == "5")

DDA_5FC <- subset(DDA, drug=="5FC")
DDA_5FCm_anc <- subset(DDA_5FC, time =="0")
DDA_5FCm <-  DDA_5FC %>%
  group_by(strain, line,time) %>%
  summarise(avgRAD20 = mean(RAD20, na.rm=TRUE),avgFoG20 = mean(FoG20, na.rm=TRUE), avgRAD50 = mean(RAD50, na.rm=TRUE),avgFoG50 = mean(FoG50, na.rm=TRUE), deltaRAD50 = avgRAD50 - median(DDA_5FCm_anc$RAD50, na.rm=TRUE), deltaFoG50 = avgFoG50 - median(DDA_5FCm_anc$FoG50, na.rm=TRUE), deltaRAD20 = avgRAD20 - median(DDA_5FCm_anc$RAD20, na.rm=TRUE), deltaFoG20 = avgFoG20 - median(DDA_5FCm_anc$FoG20, na.rm=TRUE))
DDA_5FCm_evol <- subset(DDA_5FCm, time == "5")

DDA_VCZ <- subset(DDA, drug=="VCZ")
DDA_VCZm_anc <- subset(DDA_VCZ, time =="0")
DDA_VCZm <-  DDA_VCZ %>%
  group_by(strain, line,time) %>%
  summarise(avgRAD20 = mean(RAD20, na.rm=TRUE),avgFoG20 = mean(FoG20, na.rm=TRUE), avgRAD50 = mean(RAD50, na.rm=TRUE),avgFoG50 = mean(FoG50, na.rm=TRUE), deltaRAD50 = avgRAD50 - median(DDA_VCZm_anc$RAD50, na.rm=TRUE), deltaFoG50 = avgFoG50 - median(DDA_VCZm_anc$FoG50, na.rm=TRUE), deltaRAD20 = avgRAD20 - median(DDA_VCZm_anc$RAD20, na.rm=TRUE), deltaFoG20 = avgFoG20 - median(DDA_VCZm_anc$FoG20, na.rm=TRUE))
DDA_VCZm_evol <- subset(DDA_VCZm, time == "5")

Evol <- c(3:11, 14, 22, 23)
DDA_POSm_evol_A17 <- DDA_POSm_evol %>%
  filter(line %in% Evol, strain == "A17")
DDA_POSm_evol_A17$drug <- "POS"
DDA_FLCm_evol_A17 <- DDA_FLCm_evol %>%
  filter(line %in% Evol, strain == "A17")
DDA_FLCm_evol_A17$drug <- "FLC"
DDA_CTRm_evol_A17 <- DDA_CTRm_evol %>%
  filter(line %in% Evol, strain == "A17")
DDA_CTRm_evol_A17$drug <- "CTR"
DDA_NYTm_evol_A17 <- DDA_NYTm_evol %>%
  filter(line %in% Evol, strain == "A17")
DDA_NYTm_evol_A17$drug <- "NYT"
DDA_MCZm_evol_A17 <- DDA_MCZm_evol %>%
  filter(line %in% Evol, strain == "A17")
DDA_MCZm_evol_A17$drug <- "MCZ"
DDA_5FCm_evol_A17 <- DDA_5FCm_evol %>%
  filter(line %in% Evol, strain == "A17")
DDA_5FCm_evol_A17$drug <- "5FCm"
DDA_VCZm_evol_A17 <- DDA_VCZm_evol %>%
  filter(line %in% Evol, strain == "A17")
DDA_VCZm_evol_A17$drug <- "VCZm"

evolA17 <- rbind(DDA_POSm_evol_A17, DDA_FLCm_evol_A17, DDA_CTRm_evol_A17, DDA_NYTm_evol_A17, DDA_MCZm_evol_A17, DDA_5FCm_evol_A17, DDA_VCZm_evol_A17)

evolA17m <- evolA17 %>%
  group_by(drug, line) %>%
  summarize(avgRAD20 = mean(avgRAD20), avgFoG20 = mean(avgFoG20))

write.csv(evolA17m, here("Data_Out", "DDA", "A17evol_DDA_long.csv"), row.names=FALSE)

evolA17m_wide <- pivot_wider(evolA17m, names_from = drug, values_from = c(avgRAD20, avgFoG20))

write.csv(evolA17m_wide, here("Data_Out", "DDA", "A17evol_DDA_wide.csv"), row.names=FALSE)

evolA17m_sd <- evolA17m %>%
  group_by(drug) %>%
  summarize(sd_RAD = sd(avgRAD20, na.rm=TRUE), sd_FoG = sd(avgFoG20, na.rm=TRUE))


# numerical order
plot(DDA_POSm_evol_A17$deltaRAD50, 1:12, xlim=c(6, -14), pch=19, col=Dark2[8], yaxt="n", xlab = "Change in susceptibility", ylab = "Evolved replicate", type="b")
axis(2, las=2, at=1:13, labels=13:1)
abline(v=0, lty=2)
points(DDA_CTRm_evol_A17$deltaRAD50, 1:12, col=Dark2[1], pch=21, type="b")
points(DDA_NYTm_evol_A17$deltaRAD50, 1:12, col=Dark2[2], pch=21, type="b")
points(DDA_MCZm_evol_A17$deltaRAD50, 1:12, col=Dark2[3], pch=21, type="b")
points(DDA_5FCm_evol_A17$deltaRAD50, 1:12, col=Dark2[4], pch=21, type="b")
points(DDA_VCZm_evol_A17$deltaRAD50, 1:12, col=Dark2[5], pch=21, type="b")
points(DDA_FLCm_evol_A17$deltaRAD50, 1:12, col=Dark2[6], pch=21, type="b")
abline(v=0, lty=2)
#legend(-1, 12, legend=c("POS", "FLC", "CTR", "MCZ", "VCZ", "NYT", "5FC"), col = c(Dark2[8], Dark2[6], Dark2[1], Dark2[3], Dark2[5], Dark2[2], Dark2[4]), pch=c(19, rep(21, 6)))

# karyotype order
Korder<- c(2, 12, 10, 5, 11, 4, 7, 8, 1, 3, 9, 6)
Korder_sub <-c(NA, NA, NA, NA, NA, NA, NA, NA, 1, 3, NA, NA)

plot(DDA_POSm_evol_A17$deltaRAD50[Korder], 1:12, xlim=c(6, -14), pch=19, col=Dark2[8], yaxt="n", xlab = "Change in susceptibility", ylab = "Evolved replicate", type="p")
points(DDA_POSm_evol_A17$deltaRAD50[Korder], 1:12, type="l")
axis(2, las=2, at=1:12, labels=DDA_POSm_evol_A17$line[Korder])
abline(v=0, lty=2)
points(DDA_CTRm_evol_A17$deltaRAD50[Korder], 1:12, col=Dark2[1], pch=21, type="b")
points(DDA_CTRm_evol_A17$deltaRAD50[c(6, 3)], c(12,10), type="l", col=Dark2[1])
points(DDA_NYTm_evol_A17$deltaRAD50[Korder], 1:12, col=Dark2[2], pch=21, type="b")
points(DDA_NYTm_evol_A17$deltaRAD50[c(6, 3)], c(12,10), type="l", col=Dark2[2])
points(DDA_MCZm_evol_A17$deltaRAD50[Korder], 1:12, col=Dark2[3], pch=21, type="b")
points(DDA_MCZm_evol_A17$deltaRAD50[c(6, 3)], c(12,10), type="l", col=Dark2[3])
points(DDA_5FCm_evol_A17$deltaRAD50[Korder], 1:12, col=Dark2[4], pch=21, type="b")
points(DDA_VCZm_evol_A17$deltaRAD50[Korder], 1:12, col=Dark2[5], pch=21, type="b")
points(DDA_FLCm_evol_A17$deltaRAD50[Korder], 1:12, col=Dark2[6], pch=21, type="b")
abline(v=0, lty=2)

plot(DDA_POSm_evol_A17$deltaFoG50[Korder], 1:12, xlim=c(-1, 1), pch=19, col=Dark2[8], yaxt="n", xlab = "Change in tolerance", ylab = "Evolved replicate", type="b")
#axis(2, las=2, at=1:12, labels=DDA_POSm_evol_A17$line[Korder])
axis(2, las=2, at=1:12, labels=paste0("S", 12:1))
abline(v=0, lty=2)
points(DDA_CTRm_evol_A17$deltaFoG50[Korder], 1:12, col=Dark2[1], pch=21, type="b")
points(DDA_CTRm_evol_A17$deltaFoG50[c(6, 3)], c(12,10), type="l", col=Dark2[1])
points(DDA_NYTm_evol_A17$deltaFoG50[Korder], 1:12, col=Dark2[2], pch=21, type="b")
points(DDA_NYTm_evol_A17$deltaFoG50[c(6, 3)], c(12,10), type="l", col=Dark2[2])
points(DDA_MCZm_evol_A17$deltaFoG50[Korder], 1:12, col=Dark2[3], pch=21, type="b")
points(DDA_MCZm_evol_A17$deltaFoG50[c(6, 3)], c(12,10), type="l", col=Dark2[3])
points(DDA_5FCm_evol_A17$deltaFoG50[Korder], 1:12, col=Dark2[4], pch=21, type="b")
points(DDA_VCZm_evol_A17$deltaFoG50[Korder], 1:12, col=Dark2[5], pch=21, type="b")
points(DDA_VCZm_evol_A17$deltaFoG50[c(7, 11)], c(7, 5), type="l", col=Dark2[5])
points(DDA_FLCm_evol_A17$deltaFoG50[Korder], 1:12, col=Dark2[6], pch=21, type="b")
points(DDA_FLCm_evol_A17$deltaFoG50[c(6, 3)], c(12,10), type="l", col=Dark2[6])
legend(-1, 12, legend=c("POS", "FLC", "CTR", "MCZ", "VCZ", "NYT", "5FC"), col = c(Dark2[8], Dark2[6], Dark2[1], Dark2[3], Dark2[5], Dark2[2], Dark2[4]), pch=c(19, rep(21, 6)))


# RAD/FoG20
pdf("2110Figures/FigureX-A17_CrossDrug.pdf", width=7.5, height=4)
par(mfrow=c(1, 2), mar= c(3, 5, 1, 1), oma=c(3, 3, 1, 1))
plot(DDA_POSm_evol_A17$deltaRAD20[Korder], 1:12, xlim=c(6, -12), pch=19, col=Dark2[8], yaxt="n", xlab = "Change in susceptibility", ylab = "Evolved replicate", type="p")
points(DDA_POSm_evol_A17$deltaRAD20[Korder], 1:12, type="l")
#axis(2, las=2, at=1:12, labels=DDA_POSm_evol_A17$line[Korder])
axis(2, las=2, at=1:12, labels=paste0("S", 12:1))
abline(v=0, lty=2)
points(DDA_CTRm_evol_A17$deltaRAD20[Korder], 1:12, col=Dark2[1], pch=21, type="b")
points(DDA_CTRm_evol_A17$deltaRAD20[c(6, 3)], c(12,10), type="l", col=Dark2[1])
points(DDA_NYTm_evol_A17$deltaRAD20[Korder], 1:12, col=Dark2[2], pch=21, type="b")
points(DDA_NYTm_evol_A17$deltaRAD20[c(6, 3)], c(12,10), type="l", col=Dark2[2])
points(DDA_MCZm_evol_A17$deltaRAD20[Korder], 1:12, col=Dark2[3], pch=21, type="b")
points(DDA_MCZm_evol_A17$deltaRAD20[c(6, 3)], c(12,10), type="l", col=Dark2[3])
points(DDA_5FCm_evol_A17$deltaRAD20[Korder], 1:12, col=Dark2[4], pch=21, type="b")
points(DDA_VCZm_evol_A17$deltaRAD20[Korder], 1:12, col=Dark2[5], pch=21, type="b")
points(DDA_FLCm_evol_A17$deltaRAD20[Korder], 1:12, col=Dark2[6], pch=21, type="b")
abline(v=0, lty=2)
mtext(side=1, "Change in susceptibility", line=2)

par(mar= c(3, 0, 1, 6))
plot(DDA_POSm_evol_A17$deltaFoG20[Korder], 1:12, xlim=c(-0.25, 0.6), pch=19, col=Dark2[8], yaxt="n", xlab = "Change in tolerance", ylab = "Evolved replicate", type="b")
#axis(2, las=2, at=1:12, labels=DDA_POSm_evol_A17$line[Korder])
axis(2, las=2, at=1:12, labels=FALSE)
abline(v=0, lty=2)
points(DDA_CTRm_evol_A17$deltaFoG20[Korder], 1:12, col=Dark2[1], pch=21, type="b")
points(DDA_CTRm_evol_A17$deltaFoG20[c(6, 3)], c(12,10), type="l", col=Dark2[1])
points(DDA_NYTm_evol_A17$deltaFoG20[Korder], 1:12, col=Dark2[2], pch=21, type="b")
points(DDA_NYTm_evol_A17$deltaFoG20[c(6, 3)], c(12,10), type="l", col=Dark2[2])
points(DDA_MCZm_evol_A17$deltaFoG20[Korder], 1:12, col=Dark2[3], pch=21, type="b")
points(DDA_MCZm_evol_A17$deltaFoG20[c(6, 3)], c(12,10), type="l", col=Dark2[3])
points(DDA_5FCm_evol_A17$deltaFoG20[Korder], 1:12, col=Dark2[4], pch=21, type="b")
points(DDA_VCZm_evol_A17$deltaFoG20[Korder], 1:12, col=Dark2[5], pch=21, type="b")
points(DDA_VCZm_evol_A17$deltaFoG20[c(7, 11)], c(7, 5), type="l", col=Dark2[5])
points(DDA_FLCm_evol_A17$deltaFoG20[Korder], 1:12, col=Dark2[6], pch=21, type="b")
legend(0.7, 12.5, legend=c("POS", "FLC", "CTR", "MCZ", "VCZ", "NYT", "5FC"), col = c(Dark2[8], Dark2[6], Dark2[1], Dark2[3], Dark2[5], Dark2[2], Dark2[4]), pch=c(19, rep(21, 6)), xpd=TRUE)
mtext(side=1, "Change in tolerance", line=2)
dev.off()
