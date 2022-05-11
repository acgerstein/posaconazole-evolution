#Load libraries
library(tidyverse)
library(ggforce)
library(here)
library(cowplot)
library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot)


cors <- read_csv("Data_Out/DDA/A17evol_DDA.csv")
names(cors)[3:16] <- c("RAD_5FC", "RAD_CTR", "RAD_FLC", "RAD_MCZ", "RAD_NYT", "RAD_POS", "RAD_VCZ", "FoG_5FM", "FoG_CTR", "FoG_FLC", "FoG_MCZ", "FoG_NYT", "FoG_POS", "FoG_VCZ")
cors$chr3Q <- ifelse(cors$chr3 <2.21, "2N", "3N")
cors$chr6Q <- ifelse(cors$chr6 <2.21, "2N", "3N")
cors$chrRQ <- ifelse(cors$chrR <2.21, "2N", "3N")

# singular value decomposition of the centered and scaled to have unit variance data
cors.pca <- prcomp(cors[,3:16], scale. = T, center=TRUE)

summary(cors.pca)

ggbiplot(cors.pca, varname.size = 3, varname.adjust = 1, groups = cors$chrRQ, choices = c(1, 2)) +
  xlim(-3, 3.5) +
  ylim(-2, 2) +
  theme_bw()

# add chr3, chr6 aneup in post


# singular value decomposition of the centered and scaled to have unit variance data
azole.pca <- prcomp(cors[,c(4:6, 8, 9, 11:13, 15, 16)], scale. = T, center=TRUE)

summary(azole.pca)

ggbiplot(ir.pca, obs.scale = 1, var.scale = 1, alpha=0) + geom_point(aes(shape=factor(ir.organ), colour=factor(ir.species)), size=4)

ggbiplot(azole.pca, varname.size = 4, varname.adjust = 1, groups = cors$chrRQ, choices = c(1, 2)) +
ggbiplot(azole.pca, varname.size = 4, varname.adjust = 1, alpha = 0) +
  geom_point(aes(colour=as.factor(cors$chr6Q)), size = 1) +
  #geom_point(aes(colour=as.factor(cors$chrRQ)), size = 1) +
ggbiplot(azole.pca, varname.size = 4, varname.adjust = 1) +
  xlim(-3, 3.5) +
  ylim(-2, 2) +
  theme_bw() +
  labs(col = "ChrR") +
  scale_color_manual(values = c("black", "goldenrod"))+
  xlab ("PC1 (63.9% explained var.)") +
  ylab ("PC2 (13.1% explained var.)")  +
  theme(text = element_text(size = 14))

ggbiplot(azole.pca, varname.size = 3, varname.adjust = 1, groups = cors$chrRQ, choices = c(3, 4)) +
  xlim(-3, 3.5) +
  ylim(-2, 2) +
  theme_bw() +
  xlab ("PC3 (9.7% explained var.") +
  ylab ("PC4 (6.8% explained var.")


# add chr3, chr6 aneup in post

#x = the value of the rotated data (the centred (and scaled if requested) data multiplied by the rotation matrix) is returned. Hence, cov(x) is the diagonal matrix diag(sdev^2)
PCAvar <-data.frame(cors[,3:16], prin_comp$x)
write.csv(PCAvar, "/Users/acgerstein/Documents/Postdoc/Research/CryptoClinical/UgCl/tables/171026PCAloadings.csv", row.names=FALSE)
#PCAvar.rand <-data.frame(df[,1:5], prin_comp.rand$x)

#compute standard deviation of each principal component
std_dev <- prin_comp$sdev

#compute variance (eigenvalues)
pr_var <- std_dev^2

prop_varex <- pr_var/sum(pr_var)

plot(pr_var, xlab = "Principal Component",
     ylab = "Eigenvalues",
     type = "b", pch=19, yaxt="n", xaxt="n", cex=1.2)
axis(2, las=2, cex.axis=0.8)
axis(1, at=1:21, cex.axis=0.8)
abline(h=1)

plot(1:12, PCAvar[,1], ylim=c(0, 20))
points(PCAvar[,2], col="red", pch=21)
par(new=TRUE)
text(1:12, 0, c("", "", "X", "X", "X", "","","","","","",""))

#######################
# Growth ability
#######################
LA <- read_csv(here("Data_In", "Fitness", "211119_POS-LA_FullDataSet.csv"))
LA$line <-ifelse(LA$person =="M", LA$replicate, LA$replicate+12)

LA_POS <- subset(LA, drug==0.5)

LA_POStemp <-  LA_POS %>%
  #  filter(line %in% Evol) %>%
  group_by(bioRep, strain, line,time) %>%
  summarise(nbiorep = n())

LA_POSm <-  LA_POS %>%
  group_by(strain, line, time) %>%
  summarise(OD24 = mean(OD_24h, na.rm=TRUE), OD48 = mean(OD_48h, na.rm=TRUE), OD72mix = mean(OD_72h_Mixed), OD72 = mean(OD_72h))

LA17 <- subset(LA_POSm, strain=="A17")

LA17_wide <- pivot_wider(LA17, names_from = time, values_from = c(OD24, OD48, OD72mix, OD72))
LA17_wide <- LA17_wide[!is.na(LA17_wide[,"OD24_5"]),]
LA17_wide_sub <- subset(LA17_wide, line != "21")

cors <- cors[order(cors$line),]

cor.test(c(LA17_wide_sub$OD24_5-LA17_wide_sub$OD24_0), cors$FoG_POS)
cor.test(c(LA17_wide_sub$OD72_5-LA17_wide_sub$OD72_0), cors$FoG_POS)
cor.test(LA17_wide_sub$OD24_5, cors$FoG_POS, method= "spearman")
cor.test(LA17_wide_sub$OD72_5, cors$FoG_POS, method= "spearman")
