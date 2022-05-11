#Load libraries
library(tidyverse)
library(ggforce)
library(here)
library(cowplot)
library(broom)

#Load DDA_YPD Posaconazole evolved data - contains data from the replicates evolved in 72h transfers in POS for "FLC" "CTR" "NYT" "MCZ" "5FC" "VCZ" "POS"; all strains in POS & FLC, other drugs on SC5314 evolved replicates

#A02 - P87
#A03 - GC75
#A04 - P78048
#A08 - P75016
#A10 - P76055
#A12 - T101
#A17 - SC5314
#A18 - FH1

#dark2 colours
Dark2  <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")

DDA_YPD <- read_csv(here("Data_In", "DDA", "YPDevolved", "YPD_24hT_72hT-DDA_POS.csv"))

names(DDA_YPD)[2] <- "strain"
DDA_YPD$line <- unlist(lapply(DDA_YPD$name, function(s) strsplit(s,"_")[[1]][3]))
DDA_YPD$line[DDA_YPD$line=="5a"] <- 5

################################
#POS DDA_YPD
################################
DDA_YPDm <-  DDA_YPD %>%
  group_by(strain, line, transfer) %>%
  summarise(avgRAD20 = mean(RAD20, na.rm=TRUE), avgFoG20 = mean(FoG20, na.rm=TRUE))

DDA_YPDm$ST <- paste(DDA_YPDm$strain, DDA_YPDm$transfer, sep="_")

DDA_YPD_anc <- subset(DDA_YPD, transfer =="0")
DDA_YPD_t24 <- subset(DDA_YPD, transfer =="24")
DDA_YPD_t72 <- subset(DDA_YPD, transfer =="72")

#POS Plot with all ancestral and evolved replicates side by side
RAD20 <- ggplot(data=DDA_YPDm, mapping=aes(x = strain, y=avgRAD20, col=as.factor(transfer)))+
  geom_sina(alpha=0.6,size=2.5, position = position_dodge(width=0.4))+
  #  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.5) +
  xlab("Strain")+
  ylab("Susceptibility (RAD20)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title = element_text(size=12), plot.margin = unit(c(0.2, 0.2, 0, 0.2), "cm"))+
  ylim(23,0)+
  scale_color_manual(values = c("grey","darkmagenta", "orange"),name="",labels=c("Ancestral","t24", "t72"))
#annotate("text",x=c(2,3,7,8),y=c(9.5,8.5,11,9.5),label = "*")
#annotate("segment",x=c(0.75,3.75,4.75,5.75),xend = c(1.25,4.25,5.25,6.25),y = c(10,9,11.5,10), yend =c(10,9,11.5,10))

FoG20 <- ggplot(data=DDA_YPDm, mapping=aes(x = strain, y=avgFoG20, col=as.factor(transfer)))+
  geom_sina(alpha=0.6,size=2.5, position = position_dodge(width=0.4))+
  #  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.5) +
  xlab("Strain")+
  ylab("Tolerance (FoG20)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=12), axis.title = element_text(size=12), plot.margin = unit(c(0, 0.2, 0.2, 0.2), "cm"))+
  scale_color_manual(values = c("grey","darkmagenta", "orange"),name="",labels=c("Ancestral","t24", "t72"))+
  ylim(0,1)
# annotate("text",x=c(1,3,4,5,6,7),y=c(0.67,0.62,0.7,0.64,0.71,0.51),label = "*")
#annotate("segment",x=c(0.75,2.75,3.75,4.75,5.75,6.75),xend = c(1.25,3.25,4.25,5.25,6.25,7.25),y = c(0.66,0.61,0.69,0.63,0.7,0.5), yend =c(0.66,0.61,0.69,0.63,0.7,0.5))

DDA_YPD20 <- plot_grid(RAD20,FoG20,nrow = 2,align = "hv")
DDA_YPD20
ggsave(here("Figures", "PDF", "YPD_24_72_POS-DDA_YPD_20.pdf"), DDA_YPD20, width = 5, height=5)

strain24 <- c()
strain72 <- c()
radDiff24 <- c()
radDiff72 <- c()
FoGDiff24 <- c()
FoGDiff72 <- c()
for (i in unique(DDA_YPDm$strain)){
  sub <- subset(DDA_YPDm, strain==i)
  tRAD <- aov(sub$avgRAD20~as.factor(sub$transfer))
  print(i)
  print(summary(tRAD))
  tFoG <- aov(sub$avgFoG20~as.factor(sub$transfer))
  print(summary(tFoG))

  anc <- subset(sub, transfer=="0")
  evol24 <- subset(sub, transfer=="24")
  evol72 <- subset(sub, transfer=="72")

  strain24 <- append(strain24, evol24$strain)
  strain72 <- append(strain72, evol72$strain)
  radDiff24 <- append(radDiff24, evol24$avgRAD20 - mean(anc$avgRAD20))
  radDiff72 <- append(radDiff72, evol72$avgRAD20 - mean(anc$avgRAD20))
  FoGDiff24 <- append(FoGDiff24, evol24$avgFoG20 - mean(anc$avgFoG20))
  FoGDiff72 <- append(FoGDiff72, evol72$avgFoG20 - mean(anc$avgFoG20))
}

diffs24_YPD <- data.frame(strain24, radDiff24, FoGDiff24)
diffs72_YPD <- data.frame(strain72, radDiff72, FoGDiff72)

# [1] "A02"
# Df Sum Sq Mean Sq F value     Pr(>F)
# as.factor(sub$transfer)  2  10.89   5.444   21.14 0.00000123 ***
#   Residuals               33   8.50   0.258
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# Df  Sum Sq   Mean Sq F value Pr(>F)
# as.factor(sub$transfer)  2 0.00211 0.0010549   1.606  0.216
# Residuals               33 0.02167 0.0006567
# [1] "A03"
# Df Sum Sq Mean Sq F value Pr(>F)
# as.factor(sub$transfer)  2  0.931  0.4653   1.871   0.17
# Residuals               33  8.208  0.2487
# Df  Sum Sq  Mean Sq F value   Pr(>F)
# as.factor(sub$transfer)  2 0.03279 0.016397   10.73 0.000257 ***
#   Residuals               33 0.05043 0.001528
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# [1] "A04"
# Df Sum Sq Mean Sq F value Pr(>F)
# as.factor(sub$transfer)  2  1.792  0.8958   0.995  0.381
# Residuals               33 29.708  0.9003
# Df  Sum Sq   Mean Sq F value Pr(>F)
# as.factor(sub$transfer)  2 0.00104 0.0005215   0.461  0.634
# Residuals               33 0.03730 0.0011302
# [1] "A08"
# Df Sum Sq Mean Sq F value Pr(>F)
# as.factor(sub$transfer)  2  0.014  0.0069    0.02   0.98
# Residuals               33 11.563  0.3504
# Df  Sum Sq   Mean Sq F value Pr(>F)
# as.factor(sub$transfer)  2 0.00096 0.0004778   0.368  0.695
# Residuals               33 0.04290 0.0013001
# [1] "A10"
# Df Sum Sq Mean Sq F value    Pr(>F)
# as.factor(sub$transfer)  2  18.28   9.141   12.47 0.0000988 ***
#   Residuals               32  23.45   0.733
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# Df  Sum Sq  Mean Sq F value Pr(>F)
# as.factor(sub$transfer)  2 0.01658 0.008291   2.342  0.112
# Residuals               32 0.11327 0.003540
# [1] "A12"
# Df Sum Sq Mean Sq F value   Pr(>F)
# as.factor(sub$transfer)  2  19.97   9.983   10.82 0.000244 ***
#   Residuals               33  30.46   0.923
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# Df  Sum Sq  Mean Sq F value Pr(>F)
# as.factor(sub$transfer)  2 0.01711 0.008556    3.56 0.0398 *
#   Residuals               33 0.07931 0.002403
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# [1] "A17"
# Df Sum Sq Mean Sq F value Pr(>F)
# as.factor(sub$transfer)  2  0.094  0.0469   0.048  0.953
# Residuals               31 30.296  0.9773
# Df  Sum Sq  Mean Sq F value Pr(>F)
# as.factor(sub$transfer)  2 0.00325 0.001626   0.405   0.67
# Residuals               31 0.12440 0.004013
# [1] "A18"
# Df Sum Sq Mean Sq F value Pr(>F)
# as.factor(sub$transfer)  2  1.339  0.6694   1.744  0.193
# Residuals               29 11.128  0.3837
# Df  Sum Sq  Mean Sq F value   Pr(>F)
# as.factor(sub$transfer)  2 0.03347 0.016736   10.97 0.000283 ***
#   Residuals               29 0.04423 0.001525

A02 <- subset(DDA_YPDm, strain=="A02") #rad
A03 <- subset(DDA_YPDm, strain=="A03") #FoG
A10 <- subset(DDA_YPDm, strain=="A10") #rad
A12 <- subset(DDA_YPDm, strain=="A12") #rad, Fog
A18 <- subset(DDA_YPDm, strain=="A18") #Fog

tA02 <- aov(A02$avgRAD20~as.factor(A02$transfer))
TukeyHSD(tA02)
# diff          lwr      upr     p adj
# 24-0  0.8333333  0.324922191 1.341744 0.0009002
# 72-0  1.3333333  0.824922191 1.841744 0.0000008
# 72-24 0.5000000 -0.008411142 1.008411 0.0546772
# t24 and t72 lower than ancestral

tA03 <- aov(A03$avgFoG20~as.factor(A03$transfer))
TukeyHSD(tA03)
# $`as.factor(A03$transfer)`
# diff         lwr         upr     p adj
# 24-0  -0.07083333 -0.10999368 -0.03167299 0.0002751
# 72-0  -0.01708333 -0.05624368  0.02207701 0.5387705
# 72-24  0.05375000  0.01458966  0.09291034 0.0053619
# t24 lower than ancestral or t72

tA10 <- aov(A10$avgRAD20~as.factor(A10$transfer))
TukeyHSD(tA10)
# $`as.factor(A10$transfer)`
# diff       lwr        upr     p adj
# 24-0  -1.3800505 -2.258157 -0.5019437 0.0014569
# 72-0  -1.6856061 -2.563713 -0.8074993 0.0001307
# 72-24 -0.3055556 -1.164361  0.5532498 0.6601338
# t24 and t72 higher than ancestral

tA12_rad <- aov(A12$avgRAD20~as.factor(A12$transfer))
TukeyHSD(tA12_rad)
# $`as.factor(A12$transfer)`
# diff        lwr        upr     p adj
# 24-0  -1.6736111 -2.6360087 -0.7112135 0.0004496
# 72-0  -0.2083333 -1.1707309  0.7540643 0.8566428
# 72-24  1.4652778  0.5028802  2.4276754 0.0019922

tA12_fog <- aov(A12$avgFoG20~as.factor(A12$transfer))
TukeyHSD(tA12_fog)
# t24 higher than ancestral or t72

tA18 <- aov(A18$avgFoG20~as.factor(A18$transfer))
TukeyHSD(tA18)
# diff         lwr         upr     p adj
# 24-0   0.07763889  0.03510775  0.12017003 0.0002849
# 72-0   0.02585859 -0.01749317  0.06921034 0.3183077
# 72-24 -0.05178030 -0.09204146 -0.01151915 0.0095776
# t24 higher than ancestral or t72


order <- c("A08", "A02",  "A03", "A17", "A04",  "A18", "A10", "A12")
realnames <- c("P75016", "P87",  "GC75",  "SC5314", "P78048", "FH1", "P76055", "T101")

diffs24_YPD$strain <- factor(diffs24_YPD$strain, levels=order)
diffs72_YPD$strain <- factor(diffs72_YPD$strain, levels=order)

pdf("Figures/PDF/FigureSx-DDA_YPD-YPD-24h-72h.pdf", width=5.5, height=6)
par(mfrow=c(2, 1),mar=c(1,1,1,1), oma=c(3, 3.5, 1, 1))
plot(jitter(as.numeric(as.factor(diffs24_YPD$strain)), factor=0.2), diffs24_YPD$radDiff24, ylim=c(10, -10), xlim=c(0.8,8.2), col="darkmagenta", ylab = "", xaxt="n", yaxt="n")
points(jitter(as.numeric(as.factor(diffs72_YPD$strain))+0.2, factor=0.2), diffs72_YPD$radDiff72, col="forestgreen")
abline(h=0, lty=2)
axis(1, at=1:8, labels=FALSE)
axis(2, las=2)
txt <- expression(paste(Delta," Susceptibility (", RAD[20], ")"))
mtext(txt, side=2, line=3)

plot(jitter(as.numeric(as.factor(diffs24_YPD$strain)), 0.2), diffs24_YPD$FoGDiff24, ylim=c(-1, 1), col="darkmagenta", xlim=c(0.8,8.2),ylab = "Change in Tolerance (FoG20)", xaxt="n", yaxt="n")
points(jitter(as.numeric(as.factor(diffs72_YPD$strain))+0.2, factor= 0.2), diffs72_YPD$FoGDiff72, col="forestgreen")
abline(h=0, lty=2)
axis(1, at=1:8, labels=realnames)
axis(2, las=2)
txt2 <- expression(paste(Delta," Tolerance (", FoG[20], ")"))
mtext(txt2, side=2, line=3)
legend("bottomright", pch=21, col=c("darkmagenta", "forestgreen"), legend=c("24 h transfers", "72 h transfers"), cex=0.8)
dev.off()







diffs <- data.frame(strain, radDiff, FoGDiff)

DDA_YPD_gm <-  DDA_YPD %>%
  group_by(strain, transfer) %>%
  summarise(avgRAD20 = median(RAD20, na.rm=TRUE), avgFoG20 = median(FoG20, na.rm=TRUE))

DDA_YPD_gm_anc <- subset(DDA_YPD_gm, transfer == "0")
DDA_YPD_gm_evol24 <- subset(DDA_YPD_gm, transfer == "24")
DDA_YPD_gm_evol72 <- subset(DDA_YPD_gm, transfer == "72")

DDA_YPD_POS_diff <- data.frame(strain=character(), RAD20diff = numeric(), FoG20diff = numeric())
for (i in unique(DDA_YPD_POSm_evol$strain)){
  evolTemp <- subset(DDA_YPD_POSm_evol, strain==i)
  ancTemp <- subset(DDA_YPD_POSgm_anc, strain==i)
  tempdf <- data.frame(strain =  evolTemp$strain, line = evolTemp$line, RAD20diff = evolTemp$avgRAD20-ancTemp$avgRAD20, FoG20diff = evolTemp$avgFoG20 - ancTemp$avgFoG20)
  DDA_YPD_POS_diff <- rbind(DDA_YPD_POS_diff, tempdf)
}

plot(as.numeric(as.factor(DDA_YPD_POS_diff$strain)), DDA_YPD_POS_diff$RAD20diff, ylim=c(-10, 10))
abline(h=0)

plot(as.numeric(as.factor(DDA_YPD_POS_diff$strain)), DDA_YPD_POS_diff$FoG20diff, ylim=c(-1, 1))
abline(h=0)

######################################
#FLC
######################################
DDA_YPD_FLC <- subset(DDA_YPD, drug=="FLC")
DDA_YPD_FLC_anc <- subset(DDA_YPD_FLC, time =="0")

DDA_YPD_FLCtemp <-  DDA_YPD_FLC %>%
  #  filter(line %in% Evol) %>%
  group_by(bioreplicate, strain, line,time) %>%
  summarise(nbiorep = n())

DDA_YPD_FLCm <-  DDA_YPD_FLC %>%
  group_by(strain, line, time) %>%
  summarise(avgRAD20 = mean(RAD20, na.rm=TRUE), avgFoG20 = mean(FoG20, na.rm=TRUE))

DDA_YPD_FLCm$ST <- paste(DDA_YPD_FLCm$strain, DDA_YPD_FLCm$time, sep="_")

RAD20_FLC <- ggplot(data=DDA_YPD_FLCm, mapping=aes(x = strain, y=avgRAD20, col=as.factor(time)))+
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

FoG20_FLC <- ggplot(data=DDA_YPD_FLCm, mapping=aes(x = strain, y=avgFoG20, col=as.factor(time)))+
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

DDA_YPD20_FLC <- plot_grid(RAD20_FLC,FoG20_FLC,nrow = 2,align = "hv")
DDA_YPD20_FLC
ggsave(here("Figures", "PDF", "FLC_72_FLC-DDA_YPD_20.pdf"), DDA_YPD20, width = 5, height=5)

DDA_YPD_FLCgm <-  DDA_YPD_FLC %>%
  group_by(strain, time) %>%
  summarise(avgRAD20 = mean(RAD20, na.rm=TRUE), avgFoG20 = mean(FoG20, na.rm=TRUE))

DDA_YPD_FLCgm_anc <- subset(DDA_YPD_FLCgm, time == "0")
DDA_YPD_FLCm_anc <- subset(DDA_YPD_FLCm, time == "0")
DDA_YPD_FLCm_evol <- subset(DDA_YPD_FLCm, time == "5")

radDiff <- c()
FoGDiff <- c()
strain <- c()
RADresults_FLC <- data.frame(nrow = 8,ncol=4)
FoGresults_FLC <- data.frame(nrow = 8,ncol=4)
j <- 0
for (i in unique(DDA_YPD_FLCm_evol$strain)){
  evol <- subset(DDA_YPD_FLCm_evol, strain==i)
  anc <- subset(DDA_YPD_FLCm_anc, strain==i)
  anc_mean <- subset(DDA_YPD_FLCgm_anc, strain==i)

  radDiff <- append(radDiff, evol$avgRAD20 - anc_mean$avgRAD20)
  FoGDiff <- append(FoGDiff, evol$avgFoG20 - anc_mean$avgFoG20)
  strain <- append(strain, evol$strain)

  j <- j+1
  if(length(evol$avgRAD20) > 1){
    t <- t.test(anc$avgRAD20, evol$avgRAD20)
    RADresults_FLC[j,1] <- i
    RADresults_FLC[j,2] <- round(t$estimate[2]-t$estimate[1],2)
    RADresults_FLC[j,3] <- round(t$statistic,2)
    RADresults_FLC[j,4] <- round(t$parameter,2)
    if(t$p.value >= 0.0001) RADresults_FLC[j,5] <- round(t$p.value,4)
    else RADresults_FLC[j,5] <- "< 0.0001"

    t2 <- t.test(anc$avgFoG20, evol$avgFoG20)
    FoGresults_FLC[j,1] <- i
    FoGresults_FLC[j,2] <- round(t2$estimate[2]-t2$estimate[1],2)
    FoGresults_FLC[j,3] <- round(t2$statistic,2)
    FoGresults_FLC[j,4] <- round(t2$parameter,2)
    if(t2$p.value >= 0.0001) FoGresults_FLC[j,5] <- round(t2$p.value,4)
    else FoGresults_FLC[j,5] <- "< 0.0001"
  }
  else{
    RADresults_FLC[j,1] <- i
    RADresults_FLC[j,2] <- NA
    RADresults_FLC[j,3] <- NA
    RADresults_FLC[j,4] <- NA
    FoGresults_FLC[j,1] <- i
    FoGresults_FLC[j,2] <- NA
    FoGresults_FLC[j,3] <- NA
    FoGresults_FLC[j,4] <- NA
  }
}
colnames(RADresults_FLC) <- paste(c("Strain","Evol - Anc","t-statistic","Degrees of freedom","p-value"))
colnames(FoGresults_FLC) <- paste(c("Strain","Evol - Anc","t-statistic","Degrees of freedom","p-value"))

write.csv(RADresults_FLC, file = here("Data_out", "DDA_YPD", "FLC_72_FLC_RADstatistics.csv"))
write.csv(FoGresults_FLC, file = here("Data_out", "DDA_YPD", "FLC_72_FLC_FoGstatistics.csv"))


DDA_YPD_FLCgm <-  DDA_YPD_FLC %>%
  group_by(strain, time) %>%
  summarise(avgRAD20 = median(RAD20, na.rm=TRUE), avgFoG20 = median(FoG20, na.rm=TRUE))

DDA_YPD_FLCgm_anc <- subset(DDA_YPD_FLCgm, time == "0")

DDA_YPD_FLC_diff <- data.frame(strain=character(), RAD20diff = numeric(), FoG20diff = numeric())
for (i in unique(DDA_YPD_FLCm_evol$strain)){
  evolTemp <- subset(DDA_YPD_FLCm_evol, strain==i)
  ancTemp <- subset(DDA_YPD_FLCgm_anc, strain==i)
  tempdf <- data.frame(strain =  evolTemp$strain, line = evolTemp$line, RAD20diff = evolTemp$avgRAD20-ancTemp$avgRAD20, FoG20diff = evolTemp$avgFoG20 - ancTemp$avgFoG20)
  DDA_YPD_FLC_diff <- rbind(DDA_YPD_FLC_diff, tempdf)
}

##############################
#POS & FLC Differences combined
##############################
order <- c("A08", "A02",  "A03", "A17", "A04",  "A18", "A10", "A12")
realnames <- c("P75016", "P87",  "GC75",  "SC5314", "P78048", "FH1", "P76055", "T101")

DDA_YPD_POS_diff$strain <- factor(DDA_YPD_POS_diff$strain, levels=order)
DDA_YPD_FLC_diff$strain <- factor(DDA_YPD_FLC_diff$strain, levels=order)

pdf("Figures/PDF/Figure3-DDA_YPD-POS-FLC.pdf", width=7.5, height=6)
par(mfrow=c(2, 1),mar=c(1,1,1,1), oma=c(3, 3.5, 1, 1))
plot(jitter(as.numeric(as.factor(DDA_YPD_POS_diff$strain)), factor=0.2), DDA_YPD_POS_diff$RAD20diff, ylim=c(10, -10), xlim=c(0.8,8.2), col="darkblue", ylab = "", xaxt="n", yaxt="n")
points(jitter(as.numeric(as.factor(DDA_YPD_FLC_diff$strain))+0.2, factor=0.2), DDA_YPD_FLC_diff$RAD20diff, col=Dark2[6])
abline(h=0, lty=2)
axis(1, at=1:8, labels=FALSE)
axis(2, las=2)
txt <- expression(paste(Delta," Susceptibility (", RAD[20], ")"))
mtext(txt, side=2, line=3)

plot(jitter(as.numeric(as.factor(DDA_YPD_POS_diff$strain)), 0.2), DDA_YPD_POS_diff$FoG20diff, ylim=c(-1, 1), col="darkblue", xlim=c(0.8,8.2),ylab = "Change in Tolerance (FoG20)", xaxt="n", yaxt="n")
points(jitter(as.numeric(as.factor(DDA_YPD_FLC_diff$strain))+0.2, factor= 0.2), DDA_YPD_FLC_diff$FoG20diff, col=Dark2[6])
abline(h=0, lty=2)
axis(1, at=1:8, labels=realnames)
axis(2, las=2)
txt2 <- expression(paste(Delta," Tolerance (", FoG[20], ")"))
mtext(txt2, side=2, line=3)
legend("bottomright", pch=21, col=c(Dark2[8], Dark2[6]), legend=c("POS", "FLC"), cex=0.8)
dev.off()


###################
#SC5314
##################


DDA_YPD_FLC <- subset(DDA_YPD, drug=="FLC")
DDA_YPD_FLCm_anc <- subset(DDA_YPD_FLC, time =="0")
DDA_YPD_FLCm <-  DDA_YPD_FLC %>%
  group_by(bioreplicate, strain, line,time) %>%
  summarise(nbiorep = n(), avgRAD20 = mean(RAD20, na.rm=TRUE),avgFoG20 = mean(FoG20, na.rm=TRUE), avgRAD50 = mean(RAD50, na.rm=TRUE),avgFoG50 = mean(FoG50, na.rm=TRUE), deltaRAD50 = avgRAD50 - median(DDA_YPD_FLCm_anc$RAD50, na.rm=TRUE), deltaFoG50 = avgFoG50 - median(DDA_YPD_FLCm_anc$FoG50, na.rm=TRUE), deltaRAD20 = avgRAD20 - median(DDA_YPD_FLCm_anc$RAD20, na.rm=TRUE), deltaFoG20 = avgFoG20 - median(DDA_YPD_FLCm_anc$FoG20, na.rm=TRUE))
DDA_YPD_FLCm_evol <- subset(DDA_YPD_FLCm, time == "5")

DDA_YPD_CTR <- subset(DDA_YPD, drug=="CTR")
DDA_YPD_CTRm_anc <- subset(DDA_YPD_CTR, time =="0")
DDA_YPD_CTRm <-  DDA_YPD_CTR %>%
  group_by(bioreplicate, strain, line,time) %>%
  summarise(avgRAD20 = mean(RAD20, na.rm=TRUE),avgFoG20 = mean(FoG20, na.rm=TRUE), avgRAD50 = mean(RAD50, na.rm=TRUE),avgFoG50 = mean(FoG50, na.rm=TRUE), deltaRAD50 = avgRAD50 - median(DDA_YPD_CTRm_anc$RAD50, na.rm=TRUE), deltaFoG50 = avgFoG50 - median(DDA_YPD_CTRm_anc$FoG50, na.rm=TRUE), deltaRAD20 = avgRAD20 - median(DDA_YPD_CTRm_anc$RAD20, na.rm=TRUE), deltaFoG20 = avgFoG20 - median(DDA_YPD_CTRm_anc$FoG20, na.rm=TRUE))
DDA_YPD_CTRm_evol <- subset(DDA_YPD_CTRm, time == "5")

DDA_YPD_NYT <- subset(DDA_YPD, drug=="NYT")
DDA_YPD_NYTm_anc <- subset(DDA_YPD_NYT, time =="0")
DDA_YPD_NYTm <-  DDA_YPD_NYT %>%
  group_by(strain, line,time) %>%
  summarise(avgRAD20 = mean(RAD20, na.rm=TRUE),avgFoG20 = mean(FoG20, na.rm=TRUE), avgRAD50 = mean(RAD50, na.rm=TRUE),avgFoG50 = mean(FoG50, na.rm=TRUE), deltaRAD50 = avgRAD50 - median(DDA_YPD_NYTm_anc$RAD50, na.rm=TRUE), deltaFoG50 = avgFoG50 - median(DDA_YPD_NYTm_anc$FoG50, na.rm=TRUE), deltaRAD20 = avgRAD20 - median(DDA_YPD_NYTm_anc$RAD20, na.rm=TRUE), deltaFoG20 = avgFoG20 - median(DDA_YPD_NYTm_anc$FoG20, na.rm=TRUE))
DDA_YPD_NYTm_evol <- subset(DDA_YPD_NYTm, time == "5")

DDA_YPD_MCZ <- subset(DDA_YPD, drug=="MCZ")
DDA_YPD_MCZm_anc <- subset(DDA_YPD_MCZ, time =="0")
DDA_YPD_MCZm <-  DDA_YPD_MCZ %>%
  group_by(strain, line,time) %>%
  summarise(avgRAD20 = mean(RAD20, na.rm=TRUE),avgFoG20 = mean(FoG20, na.rm=TRUE), avgRAD50 = mean(RAD50, na.rm=TRUE),avgFoG50 = mean(FoG50, na.rm=TRUE), deltaRAD50 = avgRAD50 - median(DDA_YPD_MCZm_anc$RAD50, na.rm=TRUE), deltaFoG50 = avgFoG50 - median(DDA_YPD_MCZm_anc$FoG50, na.rm=TRUE), deltaRAD20 = avgRAD20 - median(DDA_YPD_MCZm_anc$RAD20, na.rm=TRUE), deltaFoG20 = avgFoG20 - median(DDA_YPD_MCZm_anc$FoG20, na.rm=TRUE))
DDA_YPD_MCZm_evol <- subset(DDA_YPD_MCZm, time == "5")

DDA_YPD_5FC <- subset(DDA_YPD, drug=="5FC")
DDA_YPD_5FCm_anc <- subset(DDA_YPD_5FC, time =="0")
DDA_YPD_5FCm <-  DDA_YPD_5FC %>%
  group_by(strain, line,time) %>%
  summarise(avgRAD20 = mean(RAD20, na.rm=TRUE),avgFoG20 = mean(FoG20, na.rm=TRUE), avgRAD50 = mean(RAD50, na.rm=TRUE),avgFoG50 = mean(FoG50, na.rm=TRUE), deltaRAD50 = avgRAD50 - median(DDA_YPD_5FCm_anc$RAD50, na.rm=TRUE), deltaFoG50 = avgFoG50 - median(DDA_YPD_5FCm_anc$FoG50, na.rm=TRUE), deltaRAD20 = avgRAD20 - median(DDA_YPD_5FCm_anc$RAD20, na.rm=TRUE), deltaFoG20 = avgFoG20 - median(DDA_YPD_5FCm_anc$FoG20, na.rm=TRUE))
DDA_YPD_5FCm_evol <- subset(DDA_YPD_5FCm, time == "5")

DDA_YPD_VCZ <- subset(DDA_YPD, drug=="VCZ")
DDA_YPD_VCZm_anc <- subset(DDA_YPD_VCZ, time =="0")
DDA_YPD_VCZm <-  DDA_YPD_VCZ %>%
  group_by(strain, line,time) %>%
  summarise(avgRAD20 = mean(RAD20, na.rm=TRUE),avgFoG20 = mean(FoG20, na.rm=TRUE), avgRAD50 = mean(RAD50, na.rm=TRUE),avgFoG50 = mean(FoG50, na.rm=TRUE), deltaRAD50 = avgRAD50 - median(DDA_YPD_VCZm_anc$RAD50, na.rm=TRUE), deltaFoG50 = avgFoG50 - median(DDA_YPD_VCZm_anc$FoG50, na.rm=TRUE), deltaRAD20 = avgRAD20 - median(DDA_YPD_VCZm_anc$RAD20, na.rm=TRUE), deltaFoG20 = avgFoG20 - median(DDA_YPD_VCZm_anc$FoG20, na.rm=TRUE))
DDA_YPD_VCZm_evol <- subset(DDA_YPD_VCZm, time == "5")

Evol <- c(3:11, 14, 22, 23)
DDA_YPD_POSm_evol_A17 <- DDA_YPD_POSm_evol %>%
  filter(line %in% Evol, strain == "A17")
DDA_YPD_FLCm_evol_A17 <- DDA_YPD_FLCm_evol %>%
  filter(line %in% Evol, strain == "A17")
DDA_YPD_CTRm_evol_A17 <- DDA_YPD_CTRm_evol %>%
  filter(line %in% Evol, strain == "A17")
DDA_YPD_NYTm_evol_A17 <- DDA_YPD_NYTm_evol %>%
  filter(line %in% Evol, strain == "A17")
DDA_YPD_MCZm_evol_A17 <- DDA_YPD_MCZm_evol %>%
  filter(line %in% Evol, strain == "A17")
DDA_YPD_5FCm_evol_A17 <- DDA_YPD_5FCm_evol %>%
  filter(line %in% Evol, strain == "A17")
DDA_YPD_VCZm_evol_A17 <- DDA_YPD_VCZm_evol %>%
  filter(line %in% Evol, strain == "A17")

# numerical order
plot(DDA_YPD_POSm_evol_A17$deltaRAD50, 1:12, xlim=c(6, -14), pch=19, col=Dark2[8], yaxt="n", xlab = "Change in susceptibility", ylab = "Evolved replicate", type="b")
axis(2, las=2, at=1:13, labels=13:1)
abline(v=0, lty=2)
points(DDA_YPD_CTRm_evol_A17$deltaRAD50, 1:12, col=Dark2[1], pch=21, type="b")
points(DDA_YPD_NYTm_evol_A17$deltaRAD50, 1:12, col=Dark2[2], pch=21, type="b")
points(DDA_YPD_MCZm_evol_A17$deltaRAD50, 1:12, col=Dark2[3], pch=21, type="b")
points(DDA_YPD_5FCm_evol_A17$deltaRAD50, 1:12, col=Dark2[4], pch=21, type="b")
points(DDA_YPD_VCZm_evol_A17$deltaRAD50, 1:12, col=Dark2[5], pch=21, type="b")
points(DDA_YPD_FLCm_evol_A17$deltaRAD50, 1:12, col=Dark2[6], pch=21, type="b")
abline(v=0, lty=2)
#legend(-1, 12, legend=c("POS", "FLC", "CTR", "MCZ", "VCZ", "NYT", "5FC"), col = c(Dark2[8], Dark2[6], Dark2[1], Dark2[3], Dark2[5], Dark2[2], Dark2[4]), pch=c(19, rep(21, 6)))

# karyotype order
Korder<- c(2, 12, 10, 5, 11, 4, 7, 8, 1, 3, 9, 6)
Korder_sub <-c(NA, NA, NA, NA, NA, NA, NA, NA, 1, 3, NA, NA)

plot(DDA_YPD_POSm_evol_A17$deltaRAD50[Korder], 1:12, xlim=c(6, -14), pch=19, col=Dark2[8], yaxt="n", xlab = "Change in susceptibility", ylab = "Evolved replicate", type="p")
points(DDA_YPD_POSm_evol_A17$deltaRAD50[Korder], 1:12, type="l")
axis(2, las=2, at=1:12, labels=DDA_YPD_POSm_evol_A17$line[Korder])
abline(v=0, lty=2)
points(DDA_YPD_CTRm_evol_A17$deltaRAD50[Korder], 1:12, col=Dark2[1], pch=21, type="b")
points(DDA_YPD_CTRm_evol_A17$deltaRAD50[c(6, 3)], c(12,10), type="l", col=Dark2[1])
points(DDA_YPD_NYTm_evol_A17$deltaRAD50[Korder], 1:12, col=Dark2[2], pch=21, type="b")
points(DDA_YPD_NYTm_evol_A17$deltaRAD50[c(6, 3)], c(12,10), type="l", col=Dark2[2])
points(DDA_YPD_MCZm_evol_A17$deltaRAD50[Korder], 1:12, col=Dark2[3], pch=21, type="b")
points(DDA_YPD_MCZm_evol_A17$deltaRAD50[c(6, 3)], c(12,10), type="l", col=Dark2[3])
points(DDA_YPD_5FCm_evol_A17$deltaRAD50[Korder], 1:12, col=Dark2[4], pch=21, type="b")
points(DDA_YPD_VCZm_evol_A17$deltaRAD50[Korder], 1:12, col=Dark2[5], pch=21, type="b")
points(DDA_YPD_FLCm_evol_A17$deltaRAD50[Korder], 1:12, col=Dark2[6], pch=21, type="b")
abline(v=0, lty=2)

plot(DDA_YPD_POSm_evol_A17$deltaFoG50[Korder], 1:12, xlim=c(-1, 1), pch=19, col=Dark2[8], yaxt="n", xlab = "Change in tolerance", ylab = "Evolved replicate", type="b")
#axis(2, las=2, at=1:12, labels=DDA_YPD_POSm_evol_A17$line[Korder])
axis(2, las=2, at=1:12, labels=paste0("S", 12:1))
abline(v=0, lty=2)
points(DDA_YPD_CTRm_evol_A17$deltaFoG50[Korder], 1:12, col=Dark2[1], pch=21, type="b")
points(DDA_YPD_CTRm_evol_A17$deltaFoG50[c(6, 3)], c(12,10), type="l", col=Dark2[1])
points(DDA_YPD_NYTm_evol_A17$deltaFoG50[Korder], 1:12, col=Dark2[2], pch=21, type="b")
points(DDA_YPD_NYTm_evol_A17$deltaFoG50[c(6, 3)], c(12,10), type="l", col=Dark2[2])
points(DDA_YPD_MCZm_evol_A17$deltaFoG50[Korder], 1:12, col=Dark2[3], pch=21, type="b")
points(DDA_YPD_MCZm_evol_A17$deltaFoG50[c(6, 3)], c(12,10), type="l", col=Dark2[3])
points(DDA_YPD_5FCm_evol_A17$deltaFoG50[Korder], 1:12, col=Dark2[4], pch=21, type="b")
points(DDA_YPD_VCZm_evol_A17$deltaFoG50[Korder], 1:12, col=Dark2[5], pch=21, type="b")
points(DDA_YPD_VCZm_evol_A17$deltaFoG50[c(7, 11)], c(7, 5), type="l", col=Dark2[5])
points(DDA_YPD_FLCm_evol_A17$deltaFoG50[Korder], 1:12, col=Dark2[6], pch=21, type="b")
points(DDA_YPD_FLCm_evol_A17$deltaFoG50[c(6, 3)], c(12,10), type="l", col=Dark2[6])
legend(-1, 12, legend=c("POS", "FLC", "CTR", "MCZ", "VCZ", "NYT", "5FC"), col = c(Dark2[8], Dark2[6], Dark2[1], Dark2[3], Dark2[5], Dark2[2], Dark2[4]), pch=c(19, rep(21, 6)))


# RAD/FoG20
pdf("2110Figures/FigureX-A17_CrossDrug.pdf", width=7.5, height=4)
par(mfrow=c(1, 2), mar= c(3, 5, 1, 1), oma=c(3, 3, 1, 1))
plot(DDA_YPD_POSm_evol_A17$deltaRAD20[Korder], 1:12, xlim=c(6, -12), pch=19, col=Dark2[8], yaxt="n", xlab = "Change in susceptibility", ylab = "Evolved replicate", type="p")
points(DDA_YPD_POSm_evol_A17$deltaRAD20[Korder], 1:12, type="l")
#axis(2, las=2, at=1:12, labels=DDA_YPD_POSm_evol_A17$line[Korder])
axis(2, las=2, at=1:12, labels=paste0("S", 12:1))
abline(v=0, lty=2)
points(DDA_YPD_CTRm_evol_A17$deltaRAD20[Korder], 1:12, col=Dark2[1], pch=21, type="b")
points(DDA_YPD_CTRm_evol_A17$deltaRAD20[c(6, 3)], c(12,10), type="l", col=Dark2[1])
points(DDA_YPD_NYTm_evol_A17$deltaRAD20[Korder], 1:12, col=Dark2[2], pch=21, type="b")
points(DDA_YPD_NYTm_evol_A17$deltaRAD20[c(6, 3)], c(12,10), type="l", col=Dark2[2])
points(DDA_YPD_MCZm_evol_A17$deltaRAD20[Korder], 1:12, col=Dark2[3], pch=21, type="b")
points(DDA_YPD_MCZm_evol_A17$deltaRAD20[c(6, 3)], c(12,10), type="l", col=Dark2[3])
points(DDA_YPD_5FCm_evol_A17$deltaRAD20[Korder], 1:12, col=Dark2[4], pch=21, type="b")
points(DDA_YPD_VCZm_evol_A17$deltaRAD20[Korder], 1:12, col=Dark2[5], pch=21, type="b")
points(DDA_YPD_FLCm_evol_A17$deltaRAD20[Korder], 1:12, col=Dark2[6], pch=21, type="b")
abline(v=0, lty=2)
mtext(side=1, "Change in susceptibility", line=2)

par(mar= c(3, 0, 1, 6))
plot(DDA_YPD_POSm_evol_A17$deltaFoG20[Korder], 1:12, xlim=c(-0.25, 0.6), pch=19, col=Dark2[8], yaxt="n", xlab = "Change in tolerance", ylab = "Evolved replicate", type="b")
#axis(2, las=2, at=1:12, labels=DDA_YPD_POSm_evol_A17$line[Korder])
axis(2, las=2, at=1:12, labels=FALSE)
abline(v=0, lty=2)
points(DDA_YPD_CTRm_evol_A17$deltaFoG20[Korder], 1:12, col=Dark2[1], pch=21, type="b")
points(DDA_YPD_CTRm_evol_A17$deltaFoG20[c(6, 3)], c(12,10), type="l", col=Dark2[1])
points(DDA_YPD_NYTm_evol_A17$deltaFoG20[Korder], 1:12, col=Dark2[2], pch=21, type="b")
points(DDA_YPD_NYTm_evol_A17$deltaFoG20[c(6, 3)], c(12,10), type="l", col=Dark2[2])
points(DDA_YPD_MCZm_evol_A17$deltaFoG20[Korder], 1:12, col=Dark2[3], pch=21, type="b")
points(DDA_YPD_MCZm_evol_A17$deltaFoG20[c(6, 3)], c(12,10), type="l", col=Dark2[3])
points(DDA_YPD_5FCm_evol_A17$deltaFoG20[Korder], 1:12, col=Dark2[4], pch=21, type="b")
points(DDA_YPD_VCZm_evol_A17$deltaFoG20[Korder], 1:12, col=Dark2[5], pch=21, type="b")
points(DDA_YPD_VCZm_evol_A17$deltaFoG20[c(7, 11)], c(7, 5), type="l", col=Dark2[5])
points(DDA_YPD_FLCm_evol_A17$deltaFoG20[Korder], 1:12, col=Dark2[6], pch=21, type="b")
legend(0.7, 12.5, legend=c("POS", "FLC", "CTR", "MCZ", "VCZ", "NYT", "5FC"), col = c(Dark2[8], Dark2[6], Dark2[1], Dark2[3], Dark2[5], Dark2[2], Dark2[4]), pch=c(19, rep(21, 6)), xpd=TRUE)
mtext(side=1, "Change in tolerance", line=2)
dev.off()
