#This is the test file about RIR project
#author:sangyimeng
#data:2023-03-25

# Analysis prepare using arg_ranker
system("/PATH/arg_ranker -i /metagenome/read/or/contig/path/ -o /output/path/ -t Thread_num -kkdbtype 16s")
system("/bin/bash /output/path/arg_ranking/script_output/arg_ranker.sh")
# This step will used arg_ranker v2.0 from https://github.com/caozhichongchong/arg_ranker
# arg_ranker will output arg_ranking/Sample_ranking_results.txt that we need.
# Sample_ranking_results.txt is sample_data/samlpe_data1.csv

# 1.Loading R packages----
rm(list=ls())
library(patchwork)
library(ggplot2)
library(tidyverse)
library(reshape2)

# 2.Data collation----
# The results of arg_ranker and the results of ARGs risk index were used to obtain data for model1, model2
# Model1:import and collate arg_ranker for annotation and calculate risk indices and relative abundance results for 36 samples and 72 files
# Model1 data from doi:10.1038/s41467-021-25096-3
samlple_risk_rank <- read.csv("sample_data/samlpe_data1.csv",sep="\t") %>% 
  arrange(Sample)
group=read.csv("sample_data/group_information.csv",sep = ",",header = F)
group1 <- c("spd","spd","spd","spd","spd","spd","spd","spd","spd",
            "spw","spw","spw","spw","spw","spw","spw","spw","spw",
            "vfd","vfd","vfd","vfd","vfd","vfd","vfd","vfd","vfd",
            "vfw","vfw","vfw","vfw","vfw","vfw","vfw","vfw","vfw")
group2 <- c("spd1-3","spd1-3","spd1-3","spd4-6","spd4-6","spd4-6","spd7-9","spd7-9","spd7-9",
            "spw1-3","spw1-3","spw1-3","spw4-6","spw4-6","spw4-6","spw7-9","spw7-9","spw7-9",
            "vfd1-3","vfd1-3","vfd1-3","vfd4-6","vfd4-6","vfd4-6","vfd7-9","vfd7-9","vfd7-9",
            "vfw1-3","vfw1-3","vfw1-3","vfw4-6","vfw4-6","vfw4-6","vfw7-9","vfw7-9","vfw7-9")
model1 <- merge(x =samlple_risk_rank, y = group ,by.x = "Sample",by.y = 'V1') %>%
  select(sample=V2,Rank_I_per:Total_abu,Rank_I_risk:Unassessed_risk) %>% 
  group_by(sample) %>%
  summarise(Rank_I_per=sum(Rank_I_per),
            Rank_II_per=sum(Rank_II_per),
            Rank_III_per=sum(Rank_III_per),
            Rank_IV_per=sum(Rank_IV_per),
            Unassessed_per=sum(Unassessed_per),
            Total_abu=sum(Total_abu),
            Rank_I_risk=sum(Rank_I_risk),
            Rank_II_risk=sum(Rank_II_risk),
            Rank_III_risk=sum(Rank_III_risk),
            Rank_IV_risk=sum(Rank_IV_risk),
            Unassessed_risk=sum(Unassessed_risk)) %>%
  mutate(RIsample1=(1-122/4050)*Rank_I_risk+(1-145/4050)*Rank_II_risk+(1-763/4050)*Rank_III_risk+(1-2579/4050)*Rank_IV_risk+(1-3996/4050)*Unassessed_risk) %>%
  mutate(risk_1_2=Rank_I_risk + Rank_II_risk,group=group1) %>% arrange(sample) #Double-ended merging
# model2
# Model2 data from doi:10.1038/s41467-022-29283-8
risk_index <- read.csv("sample_data/ARG_riskmodel2.csv")
ARG <- read.csv("sample_data/ARG_anno.csv",sep = "\t") %>%
  select(ARO.Name,vfw1:spd9) %>%
  group_by(ARO.Name)%>%
  summarise(vfw1=sum(vfw1),vfw2=sum(vfw2),vfw3=sum(vfw3),vfw4=sum(vfw4),vfw5=sum(vfw5),vfw6=sum(vfw6),vfw7=sum(vfw7),vfw8=sum(vfw8),vfw9=sum(vfw9),
            vfd1=sum(vfd1),vfd2=sum(vfd2),vfd3=sum(vfd3),vfd4=sum(vfd4),vfd5=sum(vfd5),vfd6=sum(vfd6),vfd7=sum(vfd7),vfd8=sum(vfd8),vfd9=sum(vfd9),
            spw1=sum(spw1),spw2=sum(spw2),spw3=sum(spw3),spw4=sum(spw4),spw5=sum(spw5),spw6=sum(spw6),spw7=sum(spw7),spw8=sum(spw8),spw9=sum(spw9),
            spd1=sum(spd1),spd2=sum(spd2),spd3=sum(spd3),spd4=sum(spd4),spd5=sum(spd5),spd6=sum(spd6),spd7=sum(spd7),spd8=sum(spd8),spd9=sum(spd9))
model2 <- merge(x = ARG,y = risk_index,by.x = "ARO.Name",by.y = "ARO.Term") %>%
  mutate(vfw1_risk=vfw1*risk,vfw2_risk=vfw2*risk,vfw3_risk=vfw3*risk,vfw4_risk=vfw4*risk,vfw5_risk=vfw5*risk,vfw6_risk=vfw6*risk,vfw7_risk=vfw7*risk,vfw8_risk=vfw8*risk,vfw9_risk=vfw9*risk,
         vfd1_risk=vfd1*risk,vfd2_risk=vfd2*risk,vfd3_risk=vfd3*risk,vfd4_risk=vfd4*risk,vfd5_risk=vfd5*risk,vfd6_risk=vfd6*risk,vfd7_risk=vfd7*risk,vfd8_risk=vfd8*risk,vfd9_risk=vfd9*risk,
         spw1_risk=spw1*risk,spw2_risk=spw2*risk,spw3_risk=spw3*risk,spw4_risk=spw4*risk,spw5_risk=spw5*risk,spw6_risk=spw6*risk,spw7_risk=spw7*risk,spw8_risk=spw8*risk,spw9_risk=spw9*risk,
         spd1_risk=spd1*risk,spd2_risk=spd2*risk,spd3_risk=spd3*risk,spd4_risk=spd4*risk,spd5_risk=spd5*risk,spd6_risk=spd6*risk,spd7_risk=spd7*risk,spd8_risk=spd8*risk,spd9_risk=spd9*risk) %>%
  select(ARO.Name,vfw1_risk:spd9_risk)
rownames(model2) = model2$ARO.Name
model2=model2[,-1]
model2=rbind(model2,Total=colSums(model2))
# Consolidated results
result1 = data.frame(sample1=model1$sample,RIsample1=model1$RIsample1) %>% arrange(sample1)
model2[nrow(model2),] %>% 
  t() %>% 
  data.frame() %>% 
  write.csv(file = "~/tem")
result2=read.csv("~/tem") %>% arrange(X)
system("rm ~/tem")
result=cbind(result1,result2) %>% 
  select(sample=sample1,RIsample1,RIsample2=Total) %>%
  mutate(RIR=RIsample1 * RIsample2,group=group1) 
# Normalization
standrir <- result %>% 
  select(sample,RIsample1,RIsample2,RIR) %>% 
  cbind(group=group1)
standrir$totla_abun <- (standrir$totla_abun-min(standrir$totla_abun))/(max(standrir$totla_abun)-min(standrir$totla_abun))
standrir$RIsample1 <- (standrir$RIsample1-min(standrir$RIsample1))/(max(standrir$RIsample1)-min(standrir$RIsample1))  
standrir$RIsample2 <- (standrir$RIsample2-min(standrir$RIsample2))/(max(standrir$RIsample2)-min(standrir$RIsample2))
standrir$RIR <- (standrir$RIR-min(standrir$RIR))/(max(standrir$RIR)-min(standrir$RIR))
standrir <- cbind(standrir,RIR_round2 = round(standrir$RIR,3)) %>% 
  arrange(sample) %>% 
  cbind(totla_abun=model1$Total_abu,rank_1_2=model1$Rank_I_risk+model1$Rank_II_risk,group2)
# Only the result variable is retained
rm(list=ls()[c(1:11)])

# 3.ANOVA----
spw <- standrir %>% select(RIsample1,RIsample2,group) %>% melt(id.var="group") %>% filter(group=="spw")
spd <- standrir %>% select(RIsample1,RIsample2,group) %>% melt(id.var="group") %>% filter(group=="spd")
vfw <- standrir %>% select(RIsample1,RIsample2,group) %>% melt(id.var="group") %>% filter(group=="vfw")
vfd <- standrir %>% select(RIsample1,RIsample2,group) %>% melt(id.var="group") %>% filter(group=="vfd")

# 4.Visualization----
# Data preparation
data_pic <- standrir %>% 
  group_by(group2) %>% 
  summarise(risample1=mean(RIsample1),
            risample2=mean(RIsample2),
            rir=mean(RIR),
            abun=mean(totla_abun),
            risample1sd=sd(RIsample1),
            risample2sd=sd(RIsample2),
            rirsd=sd(RIR),
            abunsd=sd(totla_abun),
            group=group[1])
tem <- data_pic %>% select(group2,group,risample1,risample2,rir,abun) %>% melt(id.var=c("group2","group"),value.name="riskindex")
ten <- data_pic %>% select(group2,group,risample1sd,risample2sd,rirsd,abunsd) %>% melt(id.var=c("group2","group"),value.name="risksd")
dat_pic <- cbind(tem,risksd=ten$risksd)
rm(tem,ten)
dat_pic$variable <-factor(dat_pic$variable,ordered=TRUE,levels=c("abun","risample1","risample2","rir"))
head(dat_pic)
# ggplot
pic <- ggplot(data = dat_pic, aes(x = group2,y = riskindex, fill=variable))+
  geom_col(position = position_dodge2(),color = "black")+
  geom_errorbar(aes(ymin=riskindex, ymax=riskindex+risksd), width=.2, position=position_dodge(.9))+
  labs(x = "", y = "Risk Index") +
  theme_classic()+
  scale_fill_manual(values=c("#725A7A","#C56C86","#355C7D","#613030"))+
  theme(legend.position = 'top',legend.title=element_blank())
pic

# 5.Save result
ggsave(filename = "output/Abundance_vs_risk_index.pdf",plot = pic,height = 5,width = 10)
write.csv(x = standrir,file = "outoput/RIR.csv",quote = F,row.names = FALSE)
