#libraries
library(tidyverse)
library(ggforce)
library(here)
#library(conflicted)
#conflict_prefer("here", "here")

var_all <- read_tsv(here("Data_In", "Genomics", "Variants", "220420_Posaconazole.OnePerLine_ANN.tsv"))

table(var_all$HET)
# 0      1      2      3      4      5      6      7      8      9     10     11     12     13
# 7973   1505    687    682    495    556    488    462    605    621    748    948   1417   5133
# 14
# 166523

var <- subset(var_all, HET!="14" & `HOM-VAR` != "14") #14688

table(var$Annotation)

var <- subset(var, Annotation!="intergenic_region") #8762

# Gene	GeneID
# orf19.3188	TAC1
# orf19.3269	GSL2
# orf19.6000	CDR1
# orf19.5958	CDR2
# orf19.922	ERG11
# orf19.7359	CRZ1
# orf19.7372	MRR1
# orf19.391	UPC2
# orf19.5604	MDR1
# orf19.58	RRP6
# orf19.101	RIM9
# orf19.6577	FLU1 -
# orf19.328	NPR 2.00
# orf19.2842	GZF3in 1

var.TAC1 <- subset(var, Gene_ID == "orf19.3188")
# no variants differentiate evolved from ancestral (2 variants are in all backgrounds including ANC, presumable an LOH)

var.GSL2 <- subset(var, Gene_ID == "orf19.3269")
# nothing

var.CDR1 <- subset(var, Gene_ID == "orf19.6000")
# nothing that looks well-supported as anything other than wonky alignment

var.CDR2 <- subset(var, Gene_ID == "orf19.5958")
write_tsv(var.CDR2, here("Data_Out", "Genomics", "Variants", "POS_var_CDR2.tsv"))
# one heterozygous synonymous SNP in S11 (Leu1023Leu)

var.ERG11 <- subset(var, Gene_ID == "orf19.922")
# nothing

var.CRZ1 <- subset(var, Gene_ID == "orf19.7359")
# nothing well supported

var.MRR1 <- subset(var, Gene_ID == "orf19.7372")
#nothing

var.UPC2 <- subset(var, Gene_ID == "orf19.391")
# nothing

var.MDR1 <- subset(var, Gene_ID == "orf19.5604")
# nothing

var.RRP6 <- subset(var, Gene_ID == "orf19.58")
#nothing

var.RIM9 <- subset(var, Gene_ID == "orf19.101")
# nothing

var.FLU1 <- subset(var, Gene_ID == "orf19.6577")
# nothing

var.NPR <- subset(var, Gene_ID == "orf19.328")
# nothing in the gene, weird As in a downstream region

var.GZF3 <- subset(var, Gene_ID == "orf19.2842")
# LOH in positions 646493, 646559, 646560, 646608, 647030, 649877, 649963 in S4 which is a mix of upstream and downstream but there are many other het positions so not full LOH or anything like that

var_missense <- subset(var, Annotation=="missense_variant") #701

hist(table(var_missense$Gene_ID))
tab_Genes <- table(var_missense$Gene_ID)

tab_Genes10 <- subset(tab_Genes, tab_Genes > 10) # 9 genes

orf19.104 <- subset(var_all, Gene_ID == "orf19.104") #256
orf19.104_missense <- subset(var_missense, Gene_ID == "orf19.104") #27
# unknown function

#IFD6: Aldo-keto reductase; similar to aryl alcohol dehydrogenases; protein increase correlates with MDR1 overexpression (not CDR1 or CDR2) in fluconazole-resistant clinical isolates; farnesol regulated; possibly essential; Spider biofilm induced (6, 9, 10, 11, 12, 13, 14)
orf19.1048 <- subset(var_all, Gene_ID == "orf19.1048") #138
orf19.1048_missense <- subset(var_missense, Gene_ID == "orf19.1048") #11
write_tsv(orf19.1048_missense, here("Data_Out", "Genomics", "Variants", "POS_var_IFD6-orf19.1048.tsv"))

#FGR28; Protein lacking an ortholog in S. cerevisiae; transposon mutation affects filamentous growth; possibly an essential gene, disruptants not obtained by UAU1 method
orf19.1596 <- subset(var_all, Gene_ID == "orf19.1596") #177
orf19.1596_missense <- subset(var_missense, Gene_ID == "orf19.1596") #17
write_tsv(orf19.1048, here("Data_Out", "Genomics", "Variants", "POS_ALLvar_FGR28-orf19.1596.tsv"))
write_tsv(orf19.1596_missense, here("Data_Out", "Genomics", "Variants", "POS_var_FGR28-orf19.1596.tsv"))

#uncharacterized - looks like weird alignments
orf19.2449 <- subset(var_all, Gene_ID == "orf19.2449") #
orf19.2449_missense <- subset(var_missense, Gene_ID == "orf19.2449") #

# IFD1 - Protein with a NADP-dependent oxidoreductase domain; transcript induced by ketoconazole; rat catheter and Spider biofilm induced
orf19.4476 <- subset(var_all, Gene_ID == "orf19.4476") #153
orf19.4476_missense <- subset(var_missense, Gene_ID == "orf19.4476") #18
write_tsv(orf19.4476, here("Data_Out", "Genomics", "Variants", "POS_ALLvar_IFD1-orf19.4476.tsv"))
write_tsv(orf19.4476_missense, here("Data_Out", "Genomics", "Variants", "POS_var_IFD1-orf19.4476.tsv"))

orf19.5775 <- subset(var_all, Gene_ID == "orf19.5775") #73
orf19.5775_missense <- subset(var_missense, Gene_ID == "orf19.5775") #11
#Predicted ORF overlapping the Major Repeat Sequence on chromosome 6

