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
for (i in 2:2){
  pdf(file=paste('../figures/Fig5_',imputation[i],'.pdf',sep = ''),height = 5 ,width = 12)
  par(mfrow = c(4,1))
  g2 = list()
  index = 1
  for (j in 1:4){
    Performance = read.csv(file =paste('../results/',imputation[i],'_Python_Metrics_external.csv',sep = ""), header = T, sep=",")
    Performance$Method[Performance$Method=="MONGNET"] = "MOGONET"
    Performance$Method[Performance$Method=="LR_VAE"] = "LR-VAE"
    Performance$Method[Performance$Method=="LRCV_VAE"] = "LRCV-VAE"
    Performance = Performance[,c(1,2,3,j+3)]
    Performance[is.na(Performance)] = 0
    colnames(Performance)[4] = 'performance'
    Performance$Omics.used <- factor(Performance$Omics.used, levels = unique(Performance$Omics.used))
    
    Performance_matrix=dcast(Performance, Method ~ Omics.used,value.var = "performance")
    rownames(Performance_matrix) = Performance_matrix$Method
    Performance_matrix=Performance_matrix[,-1]
    
    heatmap.2(as.matrix((Performance_matrix)),col=colorRampPalette(c("navy", "white", "red"))(50), density.info="none", trace="none")
    
    
    g2[[index]] = ggplot(Performance,aes(x=Method,y=performance))+
      geom_boxplot(lwd=1,color="#2f357c")+
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
ggsave(p3,file=paste("../figures/FigS6.pdf",sep = ""),width=14, height=12)

############ ranking of omics
for (i in 2:2){
  g1 = list()
  index = 1
  for (j in 1:4){
    Performance = read.csv(file =paste('../results/',imputation[i],'_Python_Metrics_external.csv',sep = ""), header = T, sep=",")
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
      xlab('')+ylab('Performance')+
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
  p1 = do.call(ggarrange, c(g1, list(labels = c("a", "b", "c","d"),ncol = 4, nrow = 1)))
  ggsave(p1,file=paste("../figures/FigS7.pdf",sep = ""),width=10, height=2.5)
}
  