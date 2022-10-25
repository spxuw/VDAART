library(ggplot2)
library(ggpubr)
library(forcats)
library(dplyr)  
library(reshape2)
library(gplots)

rm(list=ls())
setwd("/Users/xu-wenwang/Dropbox/Projects/Disease_prediction/code")

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

Metric=c("Accuracy","F1","AUROC","AUPRC")
imputation = c("MissRF","Median","TOBMI")

p2 = list()
for (i in 1:3){
  pdf(file=paste('../figures/Fig4_',imputation[i],'.pdf',sep = ''),height = 5 ,width = 12)
  par(mfrow = c(4,1))
  g2 = list()
  index = 1
  for (j in 1:4){
    Performance = read.csv(file =paste('../results/',imputation[i],'_Python_Metrics.csv',sep = ""), header = T, sep=",")
    Performance$Method[Performance$Method=="MONGNET"] = "MOGONET"
    Performance$Method[Performance$Method=="LR_VAE"] = "LR-VAE"
    Performance$Method[Performance$Method=="LRCV_VAE"] = "LRCV-VAE"
    
    Performance = Performance[,c(1,2,3,j+3)]
    Performance[is.na(Performance)] = 0
    Performance_summary = data_summary(Performance, varname=Metric[j], groupnames=c("Method", "Omics.used"))
    
    Performance_summary$Omics.used <- factor(Performance_summary$Omics.used, levels = unique(Performance$Omics.used))
    colnames(Performance_summary)[3] = 'performance'
    
    Performance_matrix=dcast(Performance_summary, Method ~ Omics.used,value.var = "performance")
    rownames(Performance_matrix) = Performance_matrix$Method
    Performance_matrix=Performance_matrix[,-1]

    heatmap.2(as.matrix((Performance_matrix)),col=colorRampPalette(c("navy", "white", "red"))(50), density.info="none", trace="none")
    
    g2[[index]] = ggplot(Performance_summary,aes(x=Method,y=performance))+
      geom_boxplot(lwd=1,color="#2f357c")+xlab('')+ylab(Metric[j])+
      theme_bw()+
      theme(
        line = element_line(size = 0.5), # Thickness of all lines, mm
        rect = element_rect(size = 0.5), # Thickness of all rectangles, mm
        text = element_text(size = 8), # Size of all text, points
        axis.text.x = element_text(size = 8,color = 'black',angle = 90,hjust=1,vjust=0.2), # Angle x-axis text
        axis.text.y = element_text(size = 8,color = 'black'), # Angle x-axis text
        axis.title.x = element_text(size = 10), # Size x-axis title
        axis.title.y = element_text(size = 10), # Size y-axis title
        panel.grid.minor = element_blank(), # No minor grid lines
        panel.grid.major.x = element_blank(), # No major x-axis grid lines
        panel.grid.major.y = element_blank(),
        legend.title = element_blank(),
        legend.position = 'none',
        legend.box.background = element_rect(size = 0.2), # Box for legend
        legend.key.size = unit(4, unit = 'mm'),
        legend.text = element_text(size = 8),
      )
    index = index+1
  }
  dev.off() 
  p2[[i]] = do.call(ggarrange, c(g2, list(labels = c("a", "b", "c","d"),ncol = 4, nrow = 1)))
}
p3 = do.call(ggarrange, c(p2, list(ncol = 1, nrow = 3)))
ggsave(p3,file=paste("../figures/FigS1.pdf",sep = ""),width=14, height=12)

############ ranking of omics
p1 = list()
for (i in 1:3){
  g1 = list()
  index = 1
  for (j in 1:4){
    Performance = read.csv(file =paste('../results/',imputation[i],'_Python_Metrics.csv',sep = ""), header = T, sep=",")
    Performance$Method[Performance$Method=="MONGNET"] = "MOGONET"
    Performance$Method[Performance$Method=="LR_VAE"] = "LR-VAE"
    Performance$Method[Performance$Method=="LRCV_VAE"] = "LRCV-VAE"
    Performance = Performance[,c(1,2,3,j+3)]
    Performance[is.na(Performance)] = 0

    Performance$Omics.used <- factor(Performance$Omics.used, levels = unique(Performance$Omics.used))
    colnames(Performance)[4] = 'performance'
    
    Performance_summary <- Performance %>%                                        
      group_by(Omics.used) %>%                        
      summarise_at(vars(performance),list(name = median)) 
    Performance_summary = Performance_summary[order(Performance_summary$name,decreasing = T),]
    
    g1[[index]] =  ggdotchart(Performance_summary[1:10,], x = "Omics.used", y = "name",
                              color = "#4876ab",                                
                              sorting = "descending",
                              dot.size = 3,
                              add = "segments",
                              add.params = list(color = "lightgray", size = 0.5))+
      geom_hline(yintercept = 0, linetype = 1, color = "lightgray")+
      xlab('')+ylab(Metric[j])+
      theme_bw()+
      theme(
        line = element_line(size = 0.5), # Thickness of all lines, mm
        rect = element_rect(size = 0.5), # Thickness of all rectangles, mm
        text = element_text(size = 8), # Size of all text, points
        axis.text.x = element_text(size = 6,color = 'black',angle=90,hjust=0.95,vjust=0.2), # Angle x-axis text
        axis.text.y = element_text(size = 8,color = 'black'), # Angle x-axis text
        axis.title.x = element_text(size = 10), # Size x-axis title
        axis.title.y = element_text(size = 10), # Size y-axis title
        panel.grid.minor = element_blank(), # No minor grid lines
        panel.grid.major.x = element_blank(), # No major x-axis grid lines
        panel.grid.major.y = element_blank(),
        legend.title = element_blank(),
        legend.position = 'none',
        legend.box.background = element_rect(size = 0.2), # Box for legend
        legend.key.size = unit(4, unit = 'mm'),
        legend.text = element_text(size = 8)
      )
    index = index+1
  }
  p1[[i]] = do.call(ggarrange, c(g1, list(labels = c("a", "b", "c","d"),ncol = 4, nrow = 1)))
}
p2 = do.call(ggarrange, c(p1, list(ncol = 1, nrow = 3)))
ggsave(p2,file=paste("../figures/FigS4.pdf",sep = ""),width=14, height=10)

######################### different imputations #####################
Performance_all= c()
for (i in 1:3){
    Performance = read.csv(file =paste('../results/',imputation[i],'_Python_Metrics.csv',sep = ""), header = T, sep=",")
    Performance[is.na(Performance)] = 0
    dat_Performance = data.frame(Method=rep(Performance$Method,4),
                                 Metric=c(rep('Accuracy',nrow(Performance)),rep('F1',nrow(Performance)),rep('AUROC',nrow(Performance)),rep('AUPRC',nrow(Performance))),
                                 performance=c(Performance$Accuracy,Performance$F1,Performance$AUROC,Performance$AUPRC),
                                 imputation=rep(imputation[i],4*nrow(Performance)))
    Performance_all = rbind(Performance_all,dat_Performance)
}
Performance_all$imputation[Performance_all$imputation=="MissRF"] = "missForest"
g1 =  ggplot(Performance_all,aes(x=imputation,y=performance,color=imputation))+
  geom_boxplot(lwd=0.5)+scale_color_manual(values = c("#708573", "#b43749", "#4876ab", "#ffce21"))+
  theme_bw()+stat_compare_means(comparisons = list(c("missForest", "Median"), c("missForest", "TOBMI"),  c("Median", "TOBMI")),
                                    label = "p.signif",paired=TRUE)+facet_wrap(~Metric,nrow = 1)+
  theme(
    line = element_line(size = 0.5), # Thickness of all lines, mm
    rect = element_rect(size = 0.5), # Thickness of all rectangles, mm
    text = element_text(size = 8), # Size of all text, points
    axis.text.x = element_text(size = 8,color = 'black',angle = 90,hjust=1,vjust=0.2), # Angle x-axis text
    axis.text.y = element_text(size = 8,color = 'black'), # Angle x-axis text
    axis.title.x = element_text(size = 10), # Size x-axis title
    axis.title.y = element_text(size = 10), # Size y-axis title
    panel.grid.minor = element_blank(), # No minor grid lines
    panel.grid.major.x = element_blank(), # No major x-axis grid lines
    panel.grid.major.y = element_blank(),
    legend.title = element_blank(),
    legend.position = 'none',
    legend.box.background = element_rect(size = 0.2), # Box for legend
    legend.key.size = unit(4, unit = 'mm'),
    legend.text = element_text(size = 8),
  )

ggsave(g1,file=paste("../figures/FigS5.pdf",sep = ""),width=6, height=3)
  