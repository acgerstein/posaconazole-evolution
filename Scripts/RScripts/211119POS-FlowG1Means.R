#Load libraries
library(tidyverse)
library(ggforce)
library(here)
library(Hmisc)

#strain converstion:
#A02 - P87
#A03 - GC75
#A04 - P78048
#A08 - P75016
#A10 - P76055
#A12 - T101 
#A17 - SC5314
#A18 - FH1

order <- c("A02", "A08", "A03", "A04", "A17", "A18", "A10", "A12")
realnames <- c("P87", "P75016", "GC75", "P78048", "SC5314", "FH1", "P76055", "T101")

flow <- read_csv(here("Data_In", "Genomics", "FlowCytometry", "211101_POS_G1Means.csv"))
names(flow)[1] <- "datasource"
flow$strain <- unlist(lapply(flow$datasource, function(x) strsplit(x, "_")[[1]][3]))
flow$replicate <- as.numeric(unlist(lapply(flow$datasource, function(x) strsplit(x, "_")[[1]][4])))
flow$time <- unlist(lapply(flow$datasource, function(x) strsplit(x, "_")[[1]][2]))
flow$person <- unlist(lapply(flow$datasource, function(x) strsplit(x, "_")[[1]][1]))

flow$line <-ifelse(flow$person =="M", flow$replicate, flow$replicate+12)

flow$strain[flow$strain=="A2"] <- "A02"
flow$strain[flow$strain=="A3"] <- "A03"
flow$strain[flow$strain=="A4"] <- "A04"
flow$strain[flow$strain=="A8"] <- "A08"

flow$strain <- as.factor(factor(flow$strain))
#flow$strainFac <- as.factor(factor(flow$strain, levels=order))

flow_t0 <- subset(flow, time=="t0")
flow_t5 <- subset(flow, time=="t5")

#plot(as.numeric(flow_t0$strain), flow_t0$G1_mean, ylim=c(0, 250000))
#points(as.numeric(flow_t5$strain)+0.2, flow_t5$G1_mean)

pdf("Figures/PDF/Figure4b-FlowPoints.pdf", width=2.5, height=5)
par(mfrow=c(4, 2), mar=c(2, 4, 1,2), mgp=c(3,0.5,0))
for(i in order){
  print(i)
  sub_0 <- subset(flow_t0, strain == i)
  sub_5 <- subset(flow_t5, strain == i)
  df <- rbind(sub_0, sub_5)
  plot(as.numeric(as.factor(df$time)), df$G1_mean, ylim=c(100000, 250000), xlim=c(0.8, 2.2), xaxt="n", yaxt="n", xlab="", ylab="")
  axis(1, at=1:2, labels=c("anc", "evol"))
  axis(2, las=2, at=c(100000, 150000, 200000, 250000), labels=c("10", "15","20", "25"))
  mtext("FL1 intensity", side=2, line=2, cex=0.8)
}
dev.off()

#stat
t_results <-data.frame()
for(i in unique(flow$strain)){
  print(i)
  sub <- subset(flow, strain==i)
  if(length(subset(sub, time=="t5")$G1_mean) > 3){
    test<-t.test(subset(sub, time=="t0")$G1_mean, subset(sub, time=="t5")$G1_mean)
    t_results <-rbind(t_results, c(i, round(test$estimate[2]-test$estimate[1], 3), round(test$statistic, 3), round(test$parameter, 2), round(test$p.value, 4)))
  }
}

names(t_results) <- c("strain", "effectsize", "stat", "df", "p")

write.csv(t_results, "Data_out/Flow-G1mean_ttest.csv", row.names=FALSE)
