#Load libraries
library(tidyverse)
library(ggforce)
library(here)
library(Hmisc)

se <- function(x, ...) sqrt(var(x)/(length(x) - 1))

#strain converstion:
#A02 - P87
#A03 - GC75
#A04 - P78048
#A08 - P75016
#A10 - P76055
#A12 - T101 
#A17 - SC5314
#A18 - FH1

order <- c("A08", "A02",  "A03", "A17", "A04",  "A18", "A10", "A12")
realnames <- c("P75016", "P87",  "GC75",  "SC5314", "P78048", "FH1", "P76055", "T101")
LA <- read_csv(here("Data_In", "Fitness", "211119_POS-LA_FullDataSet.csv"))
LA$line <-ifelse(LA$person =="M", LA$replicate, LA$replicate+12)

#dark2 colours
Dark2  <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")

LA_POS <- subset(LA, drug==0.5)

LA_POStemp <-  LA_POS %>%
  #  filter(line %in% Evol) %>%
  group_by(bioRep, strain, line,time) %>%
  summarise(nbiorep = n())

LA_POSm <-  LA_POS %>%
  group_by(strain, line, time) %>%
  summarise(OD24 = mean(OD_24h, na.rm=TRUE), OD48 = mean(OD_48h, na.rm=TRUE), OD72mix = mean(OD_72h_Mixed), OD72 = mean(OD_72h))

cor.test(LA_POSm$OD72mix, LA_POSm$OD72)
#t = 13.663, df = 194, p-value < 2.2e-16

LA_POSm$col <- ifelse(LA_POSm$time == 0, "darkred", "darkblue")
LA_POSm$SLT <- paste(LA_POSm$strain, LA_POSm$line, LA_POSm$time, sep="_")

t24_results <-data.frame()
t72_results <-data.frame()
for(i in unique(LA_POSm$strain)){
  print(i)
  sub <- subset(LA_POSm, strain==i)
  if(length(subset(sub, time==5)$OD24) > 2){
    test24<-t.test(subset(sub, time==0)$OD24, subset(sub, time==5)$OD24)
    test72<-t.test(subset(sub, time==0)$OD72, subset(sub, time==5)$OD72)
    t24_results <-rbind(t24_results, c(i, round(test24$estimate[2]-test24$estimate[1], 3), round(test24$statistic, 3), round(test24$parameter, 2), round(test24$p.value, 4)))
    t72_results <-rbind(t72_results, c(i, round(test72$estimate[2]-test72$estimate[1], 3), round(test72$statistic, 3), round(test72$parameter, 2), round(test72$p.value, 4)))
    #print(test24)
    #print(test72)
  }
  else{
    t24_results <-rbind(t24_results, c(i, NA, NA, NA, NA))
    t72_results <-rbind(t72_results, c(i, NA, NA, NA, NA))
  }
}

names(t24_results) <- c("strain", "effectsize", "stat", "df", "p")
names(t72_results) <- c("strain", "effectsize", "stat", "df", "p")

write.csv(t24_results, "Data_out/LA-POS-24hOD_ttest.csv", row.names=FALSE)
write.csv(t72_results, "Data_out/LA-POS-72hOD_ttest.csv", row.names=FALSE)

pdf("Figures/PDF/Figure2-LA-POS.pdf", width=7.5, height=4)
k <- 0
par(mfrow=c(2, 4), mar=c(0.5,1,0.5,1), oma=c(3.5, 3.5, 1, 1))
for(j in order){
  k <- k+1
  sub<- subset(LA_POSm, strain==j)
  plot(jitter(rep(1, nrow(sub))), sub$OD24, ylim=c(0, 1.6), xlim=c(0.8, 2.2), col=sub$col, xaxt="n", yaxt="n")
  #points(jitter(rep(2, nrow(sub))), sub$OD48, col=sub$col)
  points(jitter(rep(2, nrow(sub))), sub$OD72, col=sub$col)
  for(i in sub$SLT){
    temp <- subset(sub, SLT==i)
    points(c(1,2), c(temp$OD24,  temp$OD72), type="l", col=temp$col)
  }
  if(k > 4) axis(1, at=c(1, 2), labels=c("24h", "72h"))
  else axis(1, at=c(1, 2), labels=FALSE)
 
  if(k %% 4 ==1) axis(2, las=2)
  else axis(2, labels=FALSE)
  
  if(k==8) legend("topright", pch=21, col=c("red", "darkblue"), legend=c("Anc", "Evol"))
  if(j %nin% c("A02", "A03", "A08")) text(1, max(sub$OD24)+0.1, "*", cex=2)
  if(j %nin% c("A02", "A03", "A08", "A12")) text(2, max(sub$OD72)+0.1, "*", cex=2)
  text(0.8, 1.5, realnames[k], cex=1.2, pos=4)
}
mtext("Optical density (OD)", side=2, line=2, outer=TRUE)
mtext("Time of reading (hr)", side=1, line=2, outer=TRUE)
dev.off()

write.csv(LA_POSgm, here("Data_out", "Fitness", "AncestralFitness_OD.csv"), row.names=FALSE)

#############################
#YPD
#############################

LA_YPD <- subset(LA, drug==0)
LA_YPD_anc <- subset(LA_YPD, time =="0")

LA_YPDtemp <-  LA_YPD %>%
  #  filter(line %in% Evol) %>%
  group_by(bioRep, strain, line,time) %>%
  summarise(nbiorep = n())

LA_YPDm <-  LA_YPD %>%
  group_by(strain, line, time) %>%
  summarise(OD24 = mean(OD_24h, na.rm=TRUE), OD48 = mean(OD_48h, na.rm=TRUE), OD72mix = mean(OD_72h_Mixed), OD72 = mean(OD_72h))

cor.test(LA_YPDm$OD72mix, LA_YPDm$OD72)
#t = 18.755, df = 194, p-value < 2.2e-16

LA_YPDm$col <- ifelse(LA_YPDm$time == 0, "red", "deepskyblue")
LA_YPDm$SLT <- paste(LA_YPDm$strain, LA_YPDm$line, LA_YPDm$time, sep="_")

t24_results_YPD <-data.frame()
t72_results_YPD <-data.frame()
for(i in unique(LA_YPDm$strain)){
  print(i)
  sub <- subset(LA_YPDm, strain==i)
  if(length(subset(sub, time==5)$OD24) > 2){
    test24<-t.test(subset(sub, time==0)$OD24, subset(sub, time==5)$OD24)
    test72<-t.test(subset(sub, time==0)$OD72, subset(sub, time==5)$OD72)
    t24_results_YPD <-rbind(t24_results_YPD, c(i, round(test24$estimate[2]-test24$estimate[1], 3), round(test24$statistic, 3), round(test24$parameter, 2), round(test24$p.value, 4)))
    t72_results_YPD <-rbind(t72_results_YPD, c(i, round(test72$estimate[2]-test72$estimate[1], 3), round(test72$statistic, 3), round(test72$parameter, 2), round(test72$p.value, 4)))
  }
  else{
    t24_results <-rbind(t24_results, c(i, NA, NA, NA, NA))
    t72_results <-rbind(t72_results, c(i, NA, NA, NA, NA))
  }
}

names(t24_results_YPD) <- c("strain", "effectsize", "stat", "df", "p")
names(t72_results_YPD) <- c("strain", "effectsize", "stat", "df", "p")

write.csv(t24_results_YPD, "Data_out/LA-YPD-24hOD_ttest.csv", row.names=FALSE)
write.csv(t72_results_YPD, "Data_out/LA-YPD-72hOD_ttest.csv", row.names=FALSE)


pdf("Figures/PDF/Figure2-LA-YPD.pdf", width=7.5, height=4)
k <- 0
par(mfrow=c(2, 4), mar=c(0.5,1,0.5,1), oma=c(3.5, 3.5, 1, 1))
for(j in order){
  k <- k+1
  sub<- subset(LA_YPDm, strain==j)
  plot(jitter(rep(1, nrow(sub))), sub$OD24, ylim=c(0, 2), xlim=c(0.8, 2.2), col=sub$col, xaxt="n", yaxt="n")
  #points(jitter(rep(2, nrow(sub))), sub$OD48, col=sub$col)
  points(jitter(rep(2, nrow(sub))), sub$OD72, col=sub$col)
  for(i in sub$SLT){
    temp <- subset(sub, SLT==i)
    #points(c(1,2,3), c(temp$OD24, temp$OD48, temp$OD72), type="l", col=temp$col)
    points(c(1,2), c(temp$OD24, temp$OD72), type="l", col=temp$col)
  }
  #if(k > 4) axis(1, at=c(1, 2, 3), labels=c("24h", "48h", "72h"))
  if(k > 4) axis(1, at=c(1,  2), labels=c("24h", "72h"))
  #else axis(1, at=c(1, 2, 3), labels=FALSE)
  else axis(1, at=c(1, 2), labels=FALSE)
  
  if(k %% 4 ==1) axis(2, las=2)
  else axis(2, labels=FALSE)
  
  if(k==8) legend("bottomright", pch=21, col=c("red", "darkblue"), legend=c("Anc", "Evol"))
  text(1, 1.9, realnames[k], cex=1.2)
  
  if(j %in% c("A04", "A12", "A17", "A18")) text(1, max(sub$OD24)+0.1, "*", cex=2)
  if(j %in% c("A17", "A18")) text(2, max(sub$OD72)+0.1, "*", cex=2)
}
mtext("Optical density (OD)", side=2, line=2, outer=TRUE)
mtext("Time of reading (hr)", side=1, line=2, outer=TRUE)
dev.off()


#############
#MIC
#############
LA_POS_anc <- subset(LA_POS, time =="0")

LA_POSgm_anc <-  LA_POS_anc %>%
  group_by(strain, time) %>%
  summarise(OD24 = mean(OD_24h, na.rm=TRUE), OD48 = mean(OD_48h, na.rm=TRUE), OD72mix = mean(OD_72h_Mixed), OD72 = mean(OD_72h), OD24_se = se(OD_24h, na.rm=TRUE), OD48_se = se(OD_48h, na.rm=TRUE), OD72mix_se = se(OD_72h_Mixed), OD72_se = se(OD_72h))

LA_YPD_anc <- subset(LA_YPD, time =="0")

LA_YPDgm_anc <-  LA_YPD_anc %>%
  group_by(strain, time) %>%
  summarise(OD24 = mean(OD_24h, na.rm=TRUE), OD48 = mean(OD_48h, na.rm=TRUE), OD72mix = mean(OD_72h_Mixed), OD72 = mean(OD_72h), OD24_se = se(OD_24h, na.rm=TRUE), OD48_se = se(OD_48h, na.rm=TRUE), OD72mix_se = se(OD_72h_Mixed), OD72_se = se(OD_72h))

LA_POSgm_anc$OD24/LA_YPDgm_anc$OD24

# 0.09456347 0.13834895 0.19549214 0.15237744 0.47171114 0.14758398 0.08462935 0.07867496

LA_POS_evol <- subset(LA_POS, time =="5")

LA_POSgm_evol <-  LA_POS_evol %>%
  group_by(strain, replicate, time) %>%
  summarise(OD24 = mean(OD_24h, na.rm=TRUE), OD48 = mean(OD_48h, na.rm=TRUE), OD72mix = mean(OD_72h_Mixed), OD72 = mean(OD_72h), OD24_se = se(OD_24h, na.rm=TRUE), OD48_se = se(OD_48h, na.rm=TRUE), OD72mix_se = se(OD_72h_Mixed), OD72_se = se(OD_72h))

LA_YPD_evol <- subset(LA_YPD, time =="5")

LA_YPDgm_evol <-  LA_YPD_evol %>%
  group_by(strain, replicate, time) %>%
  summarise(OD24 = mean(OD_24h, na.rm=TRUE), OD48 = mean(OD_48h, na.rm=TRUE), OD72mix = mean(OD_72h_Mixed), OD72 = mean(OD_72h), OD24_se = se(OD_24h, na.rm=TRUE), OD48_se = se(OD_48h, na.rm=TRUE), OD72mix_se = se(OD_72h_Mixed), OD72_se = se(OD_72h))

MIC5 <- LA_POSgm_evol$OD24/LA_YPDgm_evol$OD24

LA_POSgm_evol$MIC5 <- MIC5
t<- subset(LA_POSgm_evol, MIC5 > 0.5)
# Only 12 evolved replicates had MIC50 above the evolutionary level of drug, A4 x3, A8 x1, A10 x7
#################
#STAT
#################

for(i in unique(LA_POSm$strain)){
    print(i)
    sub <- subset(LA_POSm, strain==i)
  sub_long <- pivot_longer(sub, cols=c("OD24", "OD48", "OD72", "OD72mix"), names_to= "OD")
  sub_long_sub <- subset(sub_long, OD %in% c("OD24", "OD72"))
  test<-aov(sub_long_sub$value~sub_long_sub$time+sub_long_sub$OD)
  print(summary(test))
}

# [1] "A02"
# Df  Sum Sq Mean Sq F value       Pr(>F)    
# sub_long_sub$time  1 0.00610 0.00610    2.85       0.0977 .  
# sub_long_sub$OD    1 0.09153 0.09153   42.78 0.0000000342 ***
#   Residuals         49 0.10484 0.00214                         
# 
# [1] "A03"
# Df Sum Sq Mean Sq F value   Pr(>F)    
# sub_long_sub$time  1  0.036   0.036   0.984    0.325    
# sub_long_sub$OD    1  3.981   3.981 109.100 2.23e-15 ***
#   Residuals         63  2.299   0.036                     
# 
# [1] "A04"
# Df Sum Sq Mean Sq F value       Pr(>F)    
# sub_long_sub$time  1  1.627  1.6275   49.67 0.0000000245 ***
#   sub_long_sub$OD    1  1.726  1.7261   52.67 0.0000000129 ***
#   Residuals         37  1.212  0.0328                         
# 
# 
# [1] "A08"
# Df Sum Sq Mean Sq F value        Pr(>F)    
# sub_long_sub$time  1 0.5484  0.5484   35.23 0.00000474137 ***
#   sub_long_sub$OD    1 1.2980  1.2980   83.37 0.00000000413 ***
#   Residuals         23 0.3581  0.0156                          
# 
# [1] "A10"
# Df Sum Sq Mean Sq F value             Pr(>F)    
# sub_long_sub$time  1 0.2448  0.2448    14.9           0.000376 ***
#   sub_long_sub$OD    1 1.9988  1.9988   121.6 0.0000000000000407 ***
#   Residuals         43 0.7066  0.0164                               
# 
# [1] "A12"
# Df Sum Sq Mean Sq F value           Pr(>F)    
# sub_long_sub$time  1 0.1262  0.1262   7.351          0.00959 ** 
#   sub_long_sub$OD    1 1.4838  1.4838  86.458 0.00000000000749 ***
#   Residuals         43 0.7379  0.0172                             
# 
# [1] "A17"
# Df Sum Sq Mean Sq F value     Pr(>F)    
# sub_long_sub$time  1 0.4998  0.4998   20.43 0.00002416 ***
#   sub_long_sub$OD    1 0.6681  0.6681   27.30 0.00000167 ***
#   Residuals         71 1.7372  0.0245                       
# 
# [1] "A18"
# Df Sum Sq Mean Sq F value           Pr(>F)    
# sub_long_sub$time  1 0.1248  0.1248   84.85 0.00000000000983 ***
#   sub_long_sub$OD    1 0.3515  0.3515  239.07          < 2e-16 ***
#   Residuals         43 0.0632  0.0015                             



############
# Combined
############
pdf("Figures/PDF/Figure2-LA-POS-YPD.pdf", width=7.5, height=4.5)
k <- 0
par(mfrow=c(2, 4), mar=c(0.5,1,0.5,1), oma=c(3.5, 3.5, 1, 1))
for(j in order){
  k <- k+1
  sub<- subset(LA_POSm, strain==j)
  plot(jitter(rep(1, nrow(sub))), sub$OD24, ylim=c(0, 2), xlim=c(0.8, 2.2), col=sub$col, xaxt="n", yaxt="n")
  #points(jitter(rep(2, nrow(sub))), sub$OD48, col=sub$col)
  points(jitter(rep(2, nrow(sub))), sub$OD72, col=sub$col)
  for(i in sub$SLT){
    temp <- subset(sub, SLT==i)
    points(c(1,2), c(temp$OD24,  temp$OD72), type="l", col=temp$col)
  }
  if(k > 4) axis(1, at=c(1, 2), labels=c("24h", "72h"))
  else axis(1, at=c(1, 2), labels=FALSE)
  
  if(k %% 4 ==1) axis(2, las=2)
  else axis(2, labels=FALSE)
  
  if(k==8) legend("topright", pch=21, col=c("red4", "blue3", "black", "black"), legend=c("Anc", "Evol", "YPD", "POS"), lty=c(0, 0, 2, 1), cex=0.8)
  if(j %nin% c("A02", "A03", "A08")) text(0.85, mean(sub$OD24), "*", cex=2.5)
  if(j %nin% c("A02", "A03", "A08", "A12")) text(2.15, mean(sub$OD72), "*", cex=2.5)
  text(1, 1.9, realnames[k], cex=1.2, pos=4)
  
  sub_YPD<- subset(LA_YPDm, strain==j)
  points(jitter(rep(1, nrow(sub_YPD))), sub_YPD$OD24, col=sub_YPD$col)
  points(jitter(rep(2, nrow(sub_YPD))), sub_YPD$OD72, col=sub_YPD$col)
  for(i in sub_YPD$SLT){
    temp_YPD <- subset(sub_YPD, SLT==i)
    points(c(1,2), c(temp_YPD$OD24, temp_YPD$OD72), type="l", col=temp_YPD$col, lty=2)
  }
  if(j %in% c("A04", "A12", "A17", "A18")) text(0.85, mean(sub_YPD$OD24), "*", cex=2.5)
  if(j %in% c("A17", "A18")) text(2.15, mean(sub_YPD$OD72), "*", cex=2.5)
}
mtext("Optical density (OD)", side=2, line=2, outer=TRUE)
mtext("Time of reading (hr)", side=1, line=2, outer=TRUE)
dev.off()
