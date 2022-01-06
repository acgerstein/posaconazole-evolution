library(tidyverse)
library(ggforce)
library(here)
library(cowplot)
library(Hmisc)

G2_Anc <- read_delim(here("Data_In", "Genomics", "CoverageAnalysis", "G2_Anc.gff.txt"), col_names = c("Chr", "Ymap", "CNV", "Start", "End", "Copy", "X7", "X9", "X9"))
G2_Anc$Chr <- unlist(lapply(G2_Anc$Chr, function(s) strsplit(s,"_")[[1]][1]))
G2_Anc$Strain <- "G2_Anc"

highG2_Anc_2.5 <- subset(G2_Anc, Copy > 2.4)
highG2_Anc_2.6 <- subset(G2_Anc, Copy > 2.5)
nrow(highG2_Anc_2.5)/nrow(G2_Anc) #0.063
nrow(highG2_Anc_2.6)/nrow(G2_Anc) #0.042

G2_Anc.list <- split(G2_Anc, G2_Anc$Chr)
G2_Anc.avg <- unlist(lapply(G2_Anc.list, function(x) median(x$Copy)))

G2_Anc.sub <- data.frame()
par(mfrow=c(2, 4), mar=c(3, 1, 1, 1))
for(i in 1:8){
  g <- G2_Anc.list[[i]]$Copy
  h <- hist(g, breaks = 10, density = 10,
            col = "lightgray", xlab = "Accuracy", main = "Overall") 
  xfit <- seq(min(g), max(g), length = 400) 
  yfit <- dnorm(xfit, mean = mean(g), sd = sd(g)) 
  yfit <- yfit * diff(h$mids[1:2]) * length(g) 
  lines(xfit, yfit, col = "black", lwd = 2)
  sub <- subset(G2_Anc.list[[i]], Copy > 2.5)
  G2_Anc.sub <- rbind(G2_Anc.sub, sub)
}
G2_Anc.sub$Strain <- "G2_Anc"
write.csv(G2_Anc.sub, here("Data_out", "Genomics", "Coverage", "G2Q_highCNV.csv"))

G4_Anc <- read_delim(here("Data_In", "Genomics", "CoverageAnalysis", "G4_Anc.gff.txt"), col_names = c("Chr", "Ymap", "CNV", "Start", "End", "Copy", "X7", "X9", "X9"))
G4_Anc$Chr <- unlist(lapply(G4_Anc$Chr, function(s) strsplit(s,"_")[[1]][1])) 
G4_Anc$Strain <- "G4_Anc"
G4_Anc.list <- split(G4_Anc, G4_Anc$Chr)
G4_Anc.avg <- unlist(lapply(G4_Anc.list, function(x) median(x$Copy)))

highG4_Anc_2.5 <- subset(G4_Anc, Copy > 2.4)
highG4_Anc_2.6 <- subset(G4_Anc, Copy > 2.5)
nrow(highG4_Anc_2.5)/nrow(G4_Anc) #0.080
nrow(highG4_Anc_2.6)/nrow(G4_Anc) #0.049

G4_Anc.sub <- data.frame()
par(mfrow=c(2, 4), mar=c(3, 1, 1, 1))
for(i in 1:8){
  g <- G4_Anc.list[[i]]$Copy
  h <- hist(g, breaks = 10, density = 10,
            col = "lightgray", xlab = "Accuracy", main = "Overall") 
  xfit <- seq(min(g), max(g), length = 400) 
  yfit <- dnorm(xfit, mean = mean(g), sd = sd(g)) 
  yfit <- yfit * diff(h$mids[1:2]) * length(g) 
  lines(xfit, yfit, col = "black", lwd = 2)
  sub <- subset(G4_Anc.list[[i]], Copy > 2.5)
  G4_Anc.sub <- rbind(G4_Anc.sub, sub)
}
G4_Anc.sub$Strain <- "G4_Anc"
write.csv(G4_Anc.sub, here("Data_out", "Genomics", "Coverage", "G2Q_highCNV.csv"))


G2_Q <- read_delim(here("Data_In", "Genomics", "CoverageAnalysis", "G2_Q.gff.txt"), col_names = c("Chr", "Ymap", "CNV", "Start", "End", "Copy", "X7", "X9", "X9"))
G2_Q$Chr <- unlist(lapply(G2_Q$Chr, function(s) strsplit(s,"_")[[1]][1])) 
G2_Q$Strain <- "G2_Q"
G2_Q.list <- split(G2_Q, G2_Q$Chr)
S10.avg <- unlist(lapply(G2_Q.list, function(x) median(x$Copy)))

G2_Q.sub <- data.frame()
par(mfrow=c(2, 4), mar=c(3, 1, 1, 1))
for(i in 1:8){
  g <- G2_Q.list[[i]]$Copy
  h <- hist(g, breaks = 10, density = 10,
            col = "lightgray", xlab = "Accuracy", main = "Overall") 
  xfit <- seq(min(g), max(g), length = 400) 
  yfit <- dnorm(xfit, mean = mean(g), sd = sd(g)) 
  yfit <- yfit * diff(h$mids[1:2]) * length(g) 
  lines(xfit, yfit, col = "black", lwd = 2)
  sub <- subset(G2_Q.list[[i]], Copy >= round(mean(g))+1)
  G2_Q.sub <- rbind(G2_Q.sub, sub)
}
G2_Q.sub$Strain <- "G2_Q"
write.csv(G2_Q.sub, here("Data_out", "Genomics", "Coverage", "G2Q_highCNV.csv"))


G3 <- read_delim(here("Data_In", "Genomics", "CoverageAnalysis", "G3_Ev.gff.txt"), col_names = c("Chr", "Ymap", "CNV", "Start", "End", "Copy", "X7", "X9", "X9"))
G3$Chr <- unlist(lapply(G3$Chr, function(s) strsplit(s,"_")[[1]][1])) 
G3$Strain <- "G3"
G3.list <- split(G3, G3$Chr)
S4.avg <- unlist(lapply(G3.list, function(x) median(x$Copy)))

G3.sub <- data.frame()
par(mfrow=c(2, 4), mar=c(3, 1, 1, 1))
for(i in 1:8){
  g <- G3.list[[i]]$Copy
  h <- hist(g, breaks = 10, density = 10,
            col = "lightgray", xlab = "Accuracy", main = "Overall") 
  xfit <- seq(min(g), max(g), length = 400) 
  yfit <- dnorm(xfit, mean = mean(g), sd = sd(g)) 
  yfit <- yfit * diff(h$mids[1:2]) * length(g) 
  lines(xfit, yfit, col = "black", lwd = 2)
  sub <- subset(G3.list[[i]], Copy >= round(mean(g))+1)
  G3.sub <- rbind(G3.sub, sub)
}
G3.sub$Strain <- "G3"
write.csv(G3.sub, here("Data_out", "Genomics", "Coverage", "G3_highCNV.csv"))

G4 <- read_delim(here("Data_In", "Genomics", "CoverageAnalysis", "G4_Ev.gff.txt"), col_names = c("Chr", "Ymap", "CNV", "Start", "End", "Copy", "X7", "X9", "X9"))
G4$Strain <- "G4"
G4$Chr <- unlist(lapply(G4$Chr, function(s) strsplit(s,"_")[[1]][1])) 
G4.list <- split(G4, G4$Chr)
S12.avg <- unlist(lapply(G4.list, function(x) median(x$Copy)))

G4.sub <- data.frame()
par(mfrow=c(2, 4), mar=c(3, 1, 1, 1))
for(i in 1:8){
  g <- G4.list[[i]]$Copy
  h <- hist(g, breaks = 10, density = 10,
            col = "lightgray", xlab = "Accuracy", main = "Overall") 
  xfit <- seq(min(g), max(g), length = 400) 
  yfit <- dnorm(xfit, mean = mean(g), sd = sd(g)) 
  yfit <- yfit * diff(h$mids[1:2]) * length(g) 
  lines(xfit, yfit, col = "black", lwd = 2)
  sub <- subset(G4.list[[i]], Copy >= round(mean(g))+1)
  G4.sub <- rbind(G4.sub, sub)
}
G4.sub$Strain <- "G4"
write.csv(G4.sub, here("Data_out", "Genomics", "Coverage", "G4_highCNV.csv"))


G5 <- read_delim(here("Data_In", "Genomics", "CoverageAnalysis", "G5_Ev.gff.txt"), col_names = c("Chr", "Ymap", "CNV", "Start", "End", "Copy", "X7", "X9", "X9"))
G5$Chr <- unlist(lapply(G5$Chr, function(s) strsplit(s,"_")[[1]][1])) 
G5$Strain <- "G5"
G5.list <- split(G5, G5$Chr)
S3.avg <- unlist(lapply(G5.list, function(x) median(x$Copy)))

G5.sub <- data.frame()
par(mfrow=c(2, 4), mar=c(3, 1, 1, 1))
for(i in 1:8){
  g <- G5.list[[i]]$Copy
  h <- hist(g, breaks = 10, density = 10,
            col = "lightgray", xlab = "Accuracy", main = "Overall") 
  xfit <- seq(min(g), max(g), length = 400) 
  yfit <- dnorm(xfit, mean = mean(g), sd = sd(g)) 
  yfit <- yfit * diff(h$mids[1:2]) * length(g) 
  lines(xfit, yfit, col = "black", lwd = 2)
  sub <- subset(G5.list[[i]], Copy >= round(mean(g))+1)
  G5.sub <- rbind(G5.sub, sub)
}
G5.sub$Strain <- "G5"
write.csv(G5.sub, here("Data_out", "Genomics", "Coverage", "G5_highCNV.csv"))

G6 <- read_delim(here("Data_In", "Genomics", "CoverageAnalysis", "G6_Ev.gff.txt"), col_names = c("Chr", "Ymap", "CNV", "Start", "End", "Copy", "X7", "X9", "X9"))
G6$Chr <- unlist(lapply(G6$Chr, function(s) strsplit(s,"_")[[1]][1])) 
G6$Strain <- "G6"
G6.list <- split(G6, G6$Chr)
S7.avg <- unlist(lapply(G6.list, function(x) median(x$Copy)))

G6.sub <- data.frame()
par(mfrow=c(2, 4), mar=c(3, 1, 1, 1))
for(i in 1:8){
  g <- G6.list[[i]]$Copy
  h <- hist(g, breaks = 10, density = 10,
            col = "lightgray", xlab = "Accuracy", main = "Overall") 
  xfit <- seq(min(g), max(g), length = 400) 
  yfit <- dnorm(xfit, mean = mean(g), sd = sd(g)) 
  yfit <- yfit * diff(h$mids[1:2]) * length(g) 
  lines(xfit, yfit, col = "black", lwd = 2)
  sub <- subset(G6.list[[i]], Copy >= round(mean(g))+1)
  G6.sub <- rbind(G6.sub, sub)
}
G6.sub$Strain <- "G6"
write.csv(G6.sub, here("Data_out", "Genomics", "Coverage", "G6_highCNV.csv"))


G7 <- read_delim(here("Data_In", "Genomics", "CoverageAnalysis", "G7_Ev.gff.txt"), col_names = c("Chr", "Ymap", "CNV", "Start", "End", "Copy", "X7", "X9", "X9"))
G7$Chr <- unlist(lapply(G7$Chr, function(s) strsplit(s,"_")[[1]][1])) 
G7$Strain <- "G7"
G7.list <- split(G7, G7$Chr)
S9.avg <- unlist(lapply(G7.list, function(x) median(x$Copy)))

G7.sub <- data.frame()
par(mfrow=c(2, 4), mar=c(3, 1, 1, 1))
for(i in 1:8){
  g <- G7.list[[i]]$Copy
  h <- hist(g, breaks = 10, density = 10,
            col = "lightgray", xlab = "Accuracy", main = "Overall") 
  xfit <- seq(min(g), max(g), length = 400) 
  yfit <- dnorm(xfit, mean = mean(g), sd = sd(g)) 
  yfit <- yfit * diff(h$mids[1:2]) * length(g) 
  lines(xfit, yfit, col = "black", lwd = 2)
  sub <- subset(G7.list[[i]], Copy >= round(mean(g))+1)
  G7.sub <- rbind(G7.sub, sub)
}
G7.sub$Strain <- "G7"
write.csv(G7.sub, here("Data_out", "Genomics", "Coverage", "G7_highCNV.csv"))


G8 <- read_delim(here("Data_In", "Genomics", "CoverageAnalysis", "G8_Ev.gff.txt"), col_names = c("Chr", "Ymap", "CNV", "Start", "End", "Copy", "X7", "X9", "X9"))
G8$Chr <- unlist(lapply(G8$Chr, function(s) strsplit(s,"_")[[1]][1])) 
G8$Strain <- "G8"
G8.list <- split(G8, G8$Chr)
S1.avg <- unlist(lapply(G8.list, function(x) median(x$Copy)))

G8.sub <- data.frame()
par(mfrow=c(2, 4), mar=c(3, 1, 1, 1))
for(i in 1:8){
  g <- G8.list[[i]]$Copy
  h <- hist(g, breaks = 10, density = 10,
            col = "lightgray", xlab = "Accuracy", main = "Overall") 
  xfit <- seq(min(g), max(g), length = 400) 
  yfit <- dnorm(xfit, mean = mean(g), sd = sd(g)) 
  yfit <- yfit * diff(h$mids[1:2]) * length(g) 
  lines(xfit, yfit, col = "black", lwd = 2)
  sub <- subset(G8.list[[i]], Copy >= round(mean(g))+1)
  G8.sub <- rbind(G8.sub, sub)
}
G8.sub$Strain <- "G8"
write.csv(G8.sub, here("Data_out", "Genomics", "Coverage", "G8_highCNV.csv"))


G9 <- read_delim(here("Data_In", "Genomics", "CoverageAnalysis", "G9_Ev.gff.txt"), col_names = c("Chr", "Ymap", "CNV", "Start", "End", "Copy", "X7", "X9", "X9"))
G9$Chr <- unlist(lapply(G9$Chr, function(s) strsplit(s,"_")[[1]][1])) 
G9$Strain <- "G9"
G9.list <- split(G9, G9$Chr)
S6.avg <- unlist(lapply(G9.list, function(x) median(x$Copy)))

G9.sub <- data.frame()
par(mfrow=c(2, 4), mar=c(3, 1, 1, 1))
for(i in 1:8){
  g <- G9.list[[i]]$Copy
  h <- hist(g, breaks = 10, density = 10,
            col = "lightgray", xlab = "Accuracy", main = "Overall") 
  xfit <- seq(min(g), max(g), length = 400) 
  yfit <- dnorm(xfit, mean = mean(g), sd = sd(g)) 
  yfit <- yfit * diff(h$mids[1:2]) * length(g) 
  lines(xfit, yfit, col = "black", lwd = 2)
  sub <- subset(G9.list[[i]], Copy >= round(mean(g))+1)
  G9.sub <- rbind(G9.sub, sub)
}
G9.sub$Strain <- "G9"
write.csv(G9.sub, here("Data_out", "Genomics", "Coverage", "G9_highCNV.csv"))


G10 <- read_delim(here("Data_In", "Genomics", "CoverageAnalysis", "G10_Ev.gff.txt"), col_names = c("Chr", "Ymap", "CNV", "Start", "End", "Copy", "X7", "X9", "X9"))
G10$Chr <- unlist(lapply(G10$Chr, function(s) strsplit(s,"_")[[1]][1])) 
G10$Strain <- "G10"
G10.list <- split(G10, G10$Chr)
S5.avg <- unlist(lapply(G10.list, function(x) median(x$Copy)))

G10.sub <- data.frame()
par(mfrow=c(2, 4), mar=c(3, 1, 1, 1))
for(i in 1:8){
  g <- G10.list[[i]]$Copy
  h <- hist(g, breaks = 10, density = 10,
            col = "lightgray", xlab = "Accuracy", main = "Overall") 
  xfit <- seq(min(g), max(g), length = 400) 
  yfit <- dnorm(xfit, mean = mean(g), sd = sd(g)) 
  yfit <- yfit * diff(h$mids[1:2]) * length(g) 
  lines(xfit, yfit, col = "black", lwd = 2)
  sub <- subset(G10.list[[i]], Copy >= round(mean(g))+1)
  G10.sub <- rbind(G10.sub, sub)
}
G10.sub$Strain <- "G10"
write.csv(G10.sub, here("Data_out", "Genomics", "Coverage", "G10_highCNV.csv"))


G10_Q <- read_delim(here("Data_In", "Genomics", "CoverageAnalysis", "G10_Q.gff.txt"), col_names = c("Chr", "Ymap", "CNV", "Start", "End", "Copy", "X7", "X9", "X9"))
G10_Q$Chr <- unlist(lapply(G10_Q$Chr, function(s) strsplit(s,"_")[[1]][1]))
G10_Q$Strain <- "G10_Q"
G10_Q.list <- split(G10_Q, G10_Q$Chr)
S8.avg <- unlist(lapply(G10_Q.list, function(x) median(x$Copy)))

G10_Q.sub <- data.frame()
par(mfrow=c(2, 4), mar=c(3, 1, 1, 1))
for(i in 1:8){
  g <- G10_Q.list[[i]]$Copy
  h <- hist(g, breaks = 10, density = 10,
            col = "lightgray", xlab = "Accuracy", main = "Overall") 
  xfit <- seq(min(g), max(g), length = 400) 
  yfit <- dnorm(xfit, mean = mean(g), sd = sd(g)) 
  yfit <- yfit * diff(h$mids[1:2]) * length(g) 
  lines(xfit, yfit, col = "black", lwd = 2)
  sub <- subset(G10_Q.list[[i]], Copy >= round(mean(g))+1)
  G10_Q.sub <- rbind(G10_Q.sub, sub)
}
G10_Q.sub$Strain <- "G10_Q"
write.csv(G10_Q.sub, here("Data_out", "Genomics", "Coverage", "G10_Q_highCNV.csv"))


G11 <- read_delim(here("Data_In", "Genomics", "CoverageAnalysis", "G11_Ev.gff.txt"), col_names = c("Chr", "Ymap", "CNV", "Start", "End", "Copy", "X7", "X9", "X9"))
G11$Chr <- unlist(lapply(G11$Chr, function(s) strsplit(s,"_")[[1]][1])) 
G11$Strain <- "G11"
G11.list <- split(G11, G11$Chr)
S2.avg <- unlist(lapply(G11.list, function(x) median(x$Copy)))

G11.sub <- data.frame()
par(mfrow=c(2, 4), mar=c(3, 1, 1, 1))
for(i in 1:8){
  g <- G11.list[[i]]$Copy
  h <- hist(g, breaks = 10, density = 10,
            col = "lightgray", xlab = "Accuracy", main = "Overall") 
  xfit <- seq(min(g), max(g), length = 400) 
  yfit <- dnorm(xfit, mean = mean(g), sd = sd(g)) 
  yfit <- yfit * diff(h$mids[1:2]) * length(g) 
  lines(xfit, yfit, col = "black", lwd = 2)
  sub <- subset(G11.list[[i]], Copy >= round(mean(g))+1)
  G11.sub <- rbind(G11.sub, sub)
}
G11.sub$Strain <- "G11"
write.csv(G11.sub, here("Data_out", "Genomics", "Coverage", "G11_highCNV.csv"))


G11_Q <- read_delim(here("Data_In", "Genomics", "CoverageAnalysis", "G11_Q.gff.txt"), col_names = c("Chr", "Ymap", "CNV", "Start", "End", "Copy", "X7", "X9", "X9"))
G11_Q$Chr <- unlist(lapply(G11_Q$Chr, function(s) strsplit(s,"_")[[1]][1]))
G11_Q$Strain <- "G11_Q"
G11_Q.list <- split(G11_Q, G11_Q$Chr)
S11.avg <- unlist(lapply(G11_Q.list, function(x) median(x$Copy)))

G11_Q.sub <- data.frame()
par(mfrow=c(2, 4), mar=c(3, 1, 1, 1))
for(i in 1:8){
  g <- G11_Q.list[[i]]$Copy
  h <- hist(g, breaks = 10, density = 10,
            col = "lightgray", xlab = "Accuracy", main = "Overall") 
  xfit <- seq(min(g), max(g), length = 400) 
  yfit <- dnorm(xfit, mean = mean(g), sd = sd(g)) 
  yfit <- yfit * diff(h$mids[1:2]) * length(g) 
  lines(xfit, yfit, col = "black", lwd = 2)
  sub <- subset(G11_Q.list[[i]], Copy >= round(mean(g))+1)
  G11_Q.sub <- rbind(G11_Q.sub, sub)
}
G11_Q.sub$Strain <- "G11_Q"
write.csv(G11_Q.sub, here("Data_out", "Genomics", "Coverage", "G11_q_highCNV.csv"))


copyNums <- rbind(S1.avg, S2.avg, S3.avg, S4.avg, S5.avg, S6.avg, S7.avg, S8.avg, S9.avg, S10.avg, S11.avg, S12.avg)

write.csv(copyNums, here("Data_out", "Genomics", "Coverage", "ChrCopyNumber.csv"))

all <- rbind(G2_Anc, G4_Anc, G2_Q, G3, G4, G5, G6, G7, G8, G9, G10, G10_Q, G11, G11_Q)
all$ChrStart <- paste(all$Chr, all$Start, sep="_")

CNVs <- rbind(G2_Anc.sub, G4_Anc.sub, G2_Q.sub, G3.sub, G4.sub, G5.sub, G6.sub, G7.sub, G8.sub, G9.sub, G10.sub, G10_Q.sub, G11.sub, G11_Q.sub)
CNVs$ChrStart <- paste(CNVs$Chr, CNVs$Start, sep="_")

CNVsAnc <- subset(CNVs, Strain %in% c("G2_Anc", "G4_Anc"))
length(unique(CNVsAnc$ChrStart)) #169 positions
CNVsEvol <- subset(CNVs, ChrStart %nin% unique(CNVsAnc$ChrStart)) 
length(unique(CNVsEvol$ChrStart)) #26 positions

write.csv(CNVsEvol, here("Data_out", "Genomics", "Coverage", "EvolvedCNVs.csv"))

CNVStrainTable <- table(CNVsEvol$ChrStart)
t <- hist(CNVStrainTable, breaks=3) 

pos3 <- subset(CNVStrainTable, CNVStrainTable==3)  #elevated in all
subset(all, ChrStart == names(pos3)[1]) #episemon-6a LTR

pos2 <- subset(CNVStrainTable, CNVStrainTable==2)  #
subset(all, ChrStart == names(pos2)[1]) # MRS Repeat sequence
subset(all, ChrStart == names(pos2)[2]) # telomere

plotFunction <- function(chr, X.list, pos, topy=4){
  AncChr <- subset(CNVsAnc, Chr==paste0("Ca21chr",chr))
  plot(G2_Anc.list[[chr]]["Copy"], type="l", col="grey", ylim=c(1, topy), xlim=c((which(G10.list[[chr]]["Start"]==pos)-20), (which(X.list[[chr]]["Start"]==pos)+20)))
  axis(2, las=2)
  points(G4_Anc.list[[chr]]["Copy"], type="l", col="grey")
  points(X.list[[chr]]["Copy"], type="l", col="purple")
  points(which(X.list[[chr]]["Start"]==pos), topy-0.25, cex=1)
#  points(which(X.list[[chr]]["Start"] %in% unique(Anc_Chr["Start"])), rep(3.5, length(unique(Anc_Chr["Start"]))), cex=0.5, col="grey", pch=19)
#  abline(h=median(X.list[[chr]]["Copy"]), type="l", col="blue")
#  abline(h=median(G2_Anc.list[[chr]]["Copy"]), type="l", col="black", lty=2)
#  abline(h=median(G4_Anc.list[[chr]]["Copy"]), type="l", col="black", lty=2)
}

plotFunction(3, G10.list, 1079536)


pos1 <- subset(CNVStrainTable, CNVStrainTable==1)  

#positions 1-16 are G10 (S5) chr3
pos1_G10_chr3 <- unlist(lapply(names(pos1)[1:16], function(s) strsplit(s,"_")[[1]][2]))
Anc_Chr3 <- subset(CNVsAnc, Chr=="Ca21chr3")

plot(G2_Anc.list[[3]]$Copy, type="l", col="grey", ylim=c(1, 4))
points(G4_Anc.list[[3]]$Copy, type="l", col="grey")
points(G10.list[[3]]$Copy, type="l", col="purple")
#abline(h=median(G10.list[[3]]$Copy, type="l"), col="blue")
#abline(h=median(G2_Anc.list[[3]]$Copy, type="l"), col="black", lty=2)
#abline(h=median(G4_Anc.list[[3]]$Copy, type="l"), col="black", lty=2)
points(which(G10.list[[3]]$Start %in% pos1_G10_chr3), rep(3.5, 16), cex=1)
points(which(G10.list[[3]]$Start %in% unique(Anc_Chr3$Start)), rep(3.5, length(unique(Anc_Chr3$Start))), cex=0.5, col="grey", pch=19)

subset(all, ChrStart == names(pos1)[17]) #11, 5, 3 already triploid, elevated in 11 (x4.1), 5 (x3.9, but not 3) #omega-6a LTR
plotFunction(6, G11.list, 778906, topy=6)
points(G5.list[[6]]$Copy, type="l", col="navy")
points(G3.list[[6]]$Copy, type="l", col="orange")
points(G2_Q.list[[6]]$Copy, type="l", col="grey")

#pos1 18-23 all G3 x 4
subset(all, ChrStart == names(pos1)[18]) # G3 x4
subset(all, ChrStart == names(pos1)[19]) # 
subset(all, ChrStart == names(pos1)[20]) # 
subset(all, ChrStart == names(pos1)[21]) # 
subset(all, ChrStart == names(pos1)[22]) # 
subset(all, ChrStart == names(pos1)[23]) # 

pos1_G3_chrR <- unlist(lapply(names(pos1)[18:23], function(s) strsplit(s,"_")[[1]][2]))
Anc_ChrR <- subset(CNVsAnc, Chr=="Ca21chrR")

plot(G2_Anc.list[[8]]$Copy, type="l", col="grey", ylim=c(1, 5))
points(G4_Anc.list[[8]]$Copy, type="l", col="grey")
points(G3.list[[8]]$Copy, type="l", col="purple")
#abline(h=median(G10.list[[3]]$Copy, type="l"), col="blue")
#abline(h=median(G2_Anc.list[[3]]$Copy, type="l"), col="black", lty=2)
#abline(h=median(G4_Anc.list[[3]]$Copy, type="l"), col="black", lty=2)
points(which(G3.list[[8]]$Start %in% pos1_G10_chr3), rep(4.5, 16), cex=1)
#points(which(G10.list[[3]]$Start %in% unique(Anc_Chr3$Start)), rep(3.5, length(unique(Anc_Chr3$Start))), cex=0.5, col="grey", pch=19)
