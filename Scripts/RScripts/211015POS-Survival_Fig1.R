# load libraries
library(tidyverse)
library(here)
library(forcats)


### FIGURE 1: Survival ###
#Loading in the data 

ex_rates <- read.csv(here("Data_In", "ExtinctionRates", "211119_ExtinctionData.csv"))
ex_rates$Strain_Name <- as.factor(ex_rates$Strain_Name) 

anc_fitness <- read.csv(here("Data_out", "Fitness", "AncestralFitness_OD.csv"))

ex_rates <- cbind(ex_rates, anc_fitness[,4:11])

survival_ancFit_24 <- ex_rates %>% 
  mutate(Strain_Name = fct_reorder(Strain_Name, Percent_Survived)) %>% 
  ggplot(aes(x=Strain_Name,y = Percent_Survived)) +
  geom_bar(stat = "identity", width=.6) + 
  geom_point(aes(y = avgOD_24h*100), size = 3, colour = "red") +
  geom_errorbar(aes(ymin=avgOD_24h*100-seOD_24h*100, ymax=avgOD_24h*100+seOD_24h*100), width=.2,
                position=position_dodge(.9)) +
  scale_y_continuous("Survival (% replicates)", 
                     sec.axis = sec_axis(~./100, name = "Ancestral Fitness (24 h OD)")) +
  xlab("Strain") +
  ylab("Percent survived") +
  theme_bw() 
#theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave("Figures/PDF/FigureSx-Survival_AncestralFitness_24.pdf", survival_ancFit_24, width = 6, height=3.5)

# Correlation bewteen survival and fitness
cor.test(ex_rates$avgOD_24h, ex_rates$Percent_Survived)
  #t = 0.90764, df = 6, p-value = 0.3991

survival_ancFit_72 <- ex_rates %>% 
  mutate(Strain_Name = fct_reorder(Strain_Name, Percent_Survived)) %>% 
  ggplot(aes(x=Strain_Name,y = Percent_Survived)) +
  geom_bar(stat = "identity", width=.6) + 
  geom_point(aes(y = avgOD_72h*90), size = 3, colour = "red") +
  geom_errorbar(aes(ymin=avgOD_72h*90-seOD_72h*90, ymax=avgOD_72h*90+seOD_72h*90), width=.2,
                position=position_dodge(.9)) +
  scale_y_continuous("Survival (% replicates)", 
                     sec.axis = sec_axis(~./90, name = "Ancestral Fitness (72 h OD)")) +
  xlab("Strain") +
  ylab("Percent survived") +
  theme_bw() 
#theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave("Figures/PDF/Figure1-Survival_AncestralFitness_72.pdf", survival_ancFit_72, width = 6, height=3.5)

# Correlation bewteen survival and fitness
cor.test(ex_rates$avgOD_72h, ex_rates$Percent_Survived)
#t = 0.60661, df = 6, p-value = 0.5663
