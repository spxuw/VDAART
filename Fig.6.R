library(ggplot2)
library(ggpubr)
library(forcats)
library(dplyr)
library(gridExtra)
library(reshape2)
library(ggrepel)

rm(list=ls())
setwd("/Users/xu-wenwang/Dropbox/Projects/Disease_prediction/code")

imputation = c("MissRF","Median","TOBMI")


for (i in 2:2){
  g1 = list()
  Performance = read.csv(file =paste('../results/',imputation[i],'_Python_Metrics_external.csv',sep = ""), header = T, sep=",")
  Performance[is.na(Performance)] = 0
  Performance = Performance[Performance$Omics.used%in%c("(2, 3)"),c(1,4,5,6,7)]
  rownames(Performance) = Performance$Method
  Performance = Performance[c(1,2),-1]
  
  Performance_covariate = read.csv(file =paste('../results/',imputation[i],'_Python_Metrics_external_covariate_23.csv',sep = ""), header = T, sep=",")
  Performance_covariate[is.na(Performance_covariate)] = 0
  rownames(Performance_covariate) = Performance_covariate$Method
  Performance_covariate = Performance_covariate[c(1,2),c(4,5,6,7)]
  
  Performance = rbind(melt(as.matrix(Performance)),melt(as.matrix(Performance_covariate)))
  Performance['types'] = c(rep("w/o covariates",8),rep("w/ covariates",8))
  
  g1 = ggplot(Performance, aes(x=Var2,y=value,fill=types)) +
    geom_bar(stat="identity", position=position_dodge())+facet_wrap(~Var1)+
    scale_fill_manual(values = c("#E1BE6A","#40B0A6"))+
    theme_bw()+xlab("")+ylab("Performance")+
    theme(
      line = element_line(size = 0.5), # Thickness of all lines, mm
      rect = element_rect(size = 0.5), # Thickness of all rectangles, mm
      text = element_text(size = 8), # Size of all text, points
      axis.text.x = element_text(size = 8,color = 'black'), # Angle x-axis text
      axis.text.y = element_text(size = 8,color = 'black'), # Angle x-axis text
      axis.title.x = element_text(size = 10), # Size x-axis title
      axis.title.y = element_text(size = 10), # Size y-axis title
      panel.grid.minor = element_blank(), # No minor grid lines
      panel.grid.major.x = element_blank(), # No major x-axis grid lines
      panel.grid.major.y = element_blank(),
      legend.title = element_blank(),
      #legend.position = 'none',
      legend.box.background = element_rect(size = 0.2), # Box for legend
      legend.key.size = unit(4, unit = 'mm'),
      legend.text = element_text(size = 8),
    )
}

########################## 24 ########################
for (i in 2:2){
  g2 = list()
  Performance = read.csv(file =paste('../results/',imputation[i],'_Python_Metrics_external.csv',sep = ""), header = T, sep=",")
  Performance[is.na(Performance)] = 0
  Performance = Performance[Performance$Omics.used%in%c("(2, 4)"),c(1,4,5,6,7)]
  rownames(Performance) = Performance$Method
  Performance = Performance[c(1,2),-1]
  
  Performance_covariate = read.csv(file =paste('../results/',imputation[i],'_Python_Metrics_external_covariate_24.csv',sep = ""), header = T, sep=",")
  Performance_covariate[is.na(Performance_covariate)] = 0
  rownames(Performance_covariate) = Performance_covariate$Method
  Performance_covariate = Performance_covariate[c(1,2),c(4,5,6,7)]
  
  Performance = rbind(melt(as.matrix(Performance)),melt(as.matrix(Performance_covariate)))
  Performance['types'] = c(rep("w/o covariates",8),rep("w/ covariates",8))
  
  g2 = ggplot(Performance, aes(x=Var2,y=value,fill=types)) +
    geom_bar(stat="identity", position=position_dodge())+facet_wrap(~Var1)+
    scale_fill_manual(values = c("#E1BE6A","#40B0A6"))+
    theme_bw()+xlab("")+ylab("Performance")+
    theme(
      line = element_line(size = 0.5), # Thickness of all lines, mm
      rect = element_rect(size = 0.5), # Thickness of all rectangles, mm
      text = element_text(size = 8), # Size of all text, points
      axis.text.x = element_text(size = 8,color = 'black'), # Angle x-axis text
      axis.text.y = element_text(size = 8,color = 'black'), # Angle x-axis text
      axis.title.x = element_text(size = 10), # Size x-axis title
      axis.title.y = element_text(size = 10), # Size y-axis title
      panel.grid.minor = element_blank(), # No minor grid lines
      panel.grid.major.x = element_blank(), # No major x-axis grid lines
      panel.grid.major.y = element_blank(),
      legend.title = element_blank(),
      #legend.position = 'none',
      legend.box.background = element_rect(size = 0.2), # Box for legend
      legend.key.size = unit(4, unit = 'mm'),
      legend.text = element_text(size = 8),
    )
}
p1 = grid.arrange(g1, g2, nrow=2)

# feature importance
feature_name = c("masthma","fasthma","mhayfever","meczema","feczema","fhayfever","mmaritalstatus",
                  "meducation","methnicity","mrace","site","assignment","trmt","basevitdng",
                  "mage","gestdays","gestage","treatment","trtmt","sitename","educ","mstat","raceeth","basevitdnmol")

Importance_23 = read.csv(file =paste('../results/Feature_importance_23.csv',sep = ""), header = F, sep=",")
Importance_24 = read.csv(file =paste('../results/Feature_importance_24.csv',sep = ""), header = F, sep=",")
dat_importance = data.frame(i2=as.numeric(Importance_23[1,601:624]),i3=as.numeric(Importance_24[1,601:624]))
dat_importance['feature_name'] = feature_name 
dat_importance = dat_importance[-c(17,18,19,20,24),]
dat_importance['feature_num'] = 1:19


pal4 <- c("#bf0000", "#f29d3d", "#a3d9b1", "#bfd9ff", "#f780ff", "#ff0000", "#ffe1bf", "#00d957", "#0020f2", "#e60099", 
          "#730f00", "#7f6600", "#336655", "#293aa6", "#a6538a", "#8c4f46", "#e5c339", "#00ffcc", "#333a66", "#40202d", 
          "#f29979", "#d8e600", "#00b3a7", "#8091ff")

g3 = ggplot(data=dat_importance,aes(x=i2,y=i3,label=feature_num))+geom_point(aes(color=factor(feature_num)),size=3)+
  stat_cor(method = "pearson", label.x = -0.25, label.y = 0.9)+
  geom_label_repel(aes(label = 1:19,fill = factor(feature_num)), color = 'white',size = 3.5) +
  scale_color_manual(values = pal4)+scale_fill_manual(values = pal4)+
  geom_smooth(method=lm, se=TRUE)+
  theme_bw()+xlab("Coefficient (miRNA-mRNA)")+ylab("Coefficient (miRNA-Microbiome)")+
  theme(
    line = element_line(size = 0.5), # Thickness of all lines, mm
    rect = element_rect(size = 0.5), # Thickness of all rectangles, mm
    text = element_text(size = 8), # Size of all text, points
    axis.text.x = element_text(size = 8,color = 'black'), # Angle x-axis text
    axis.text.y = element_text(size = 8,color = 'black'), # Angle x-axis text
    axis.title.x = element_text(size = 10), # Size x-axis title
    axis.title.y = element_text(size = 10), # Size y-axis title
    panel.grid.minor = element_blank(), # No minor grid lines
    panel.grid.major.x = element_blank(), # No major x-axis grid lines
    panel.grid.major.y = element_blank(),
    legend.title = element_blank(),
    #legend.position = 'none',
    legend.box.background = element_rect(size = 0.2), # Box for legend
    legend.key.size = unit(4, unit = 'mm'),
    legend.text = element_text(size = 8),
  )

p2 = grid.arrange(p1, g3, nrow=1)
ggsave(p2,file=paste("../figures/Fig6.pdf",sep = ""),width=11, height=4.5)

