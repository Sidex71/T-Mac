singletData <- readRDS("/mnt/8TB/users/shameed/shameed/Zheng/singletData.rds")
ZhengData <- readRDS("/mnt/8TB/users/shameed/shameed/Zheng/ZhengData.rds")

###########################################Cell proportion plots##########################################
library(ggplot2)
library(tidyverse)
sigData <- subset(singletData, subset = refined_celltype != 'Proliferative cells' & refined_celltype != 'HSC')
sigData$Groups <- sigData$Groups %>%
  str_replace('PBMC', 'blood') %>% 
  str_replace('Metastatic Tumor', 'metastasis') %>%
  str_replace('Lymph Node', 'lymphnode') %>% 
  str_replace('Ascites', 'ascites') %>%
  str_replace('Primary Tumor', 'primary tumor')

saveRDS(sigData, '/mnt/8TB/users/shameed/shameed/Zheng/sigData.rds')

p1<-DimPlot(sigData, reduction = 'umap', group.by = 'refined_celltype', label=T) + 
  ggtitle('')


png("/mnt/8TB/users/shameed/shameed/Zheng/figures/umap.png", width = 8.5, height = 5.5, units = 'in', res = 600)
p1
dev.off()

#############################
tubo_bar<-sigData@meta.data %>% group_by(`Response to chemotherapy`, refined_celltype) %>%
  tally() %>% mutate(cell_percentage= n*100/sum(n)) %>% 
  filter(`Response to chemotherapy`=='Sensitive' | `Response to chemotherapy`=='Resistant') 

p6<-ggplot(tubo_bar, aes(x=`Response to chemotherapy`, y=cell_percentage, fill= refined_celltype)) + 
  geom_bar(stat = 'identity', width = 0.2) + 
  labs(title = 'Cell distribution by response',
       y= 'Cell proportion (%)', x=NULL)+
  theme_bw() +
  theme(axis.text.x = element_text(size=18, hjust = 0.5, face = 'bold'),
        axis.text.y = element_text(size = 18, hjust = 0.5, face = 'bold'),
        legend.text = element_text(size = 17, face = 'bold'),
        legend.key.size = unit(1.2, "cm"),
        axis.title.y = element_text(size = 20, hjust = 0.5, face = 'bold'),
        plot.title = element_text(size = 25, hjust = 0.5, face = 'bold')) +
  scale_fill_discrete(name=NULL)

saveRDS(tubo_bar, '/mnt/8TB/users/shameed/shameed/Zheng/tubo_bar_cell_response.rds')

png("/mnt/8TB/users/shameed/shameed/Zheng/figures/cell_response.png", width = 10, height = 8.5, units = 'in', res = 600)
p6
dev.off()  
#################################
tubo_bar<-sigData@meta.data %>% group_by(`Response to chemotherapy`, Groups) %>%
  tally() %>% mutate(cell_percentage= n*100/sum(n)) %>% 
  filter(`Response to chemotherapy`=='Sensitive' | `Response to chemotherapy`=='Resistant') 
p6<-ggplot(tubo_bar, aes(x=`Response to chemotherapy`, y=cell_percentage, fill= Groups)) + 
  geom_bar(stat = 'identity', width = 0.2) + 
  labs(title = 'Site distribution by response',
       y= 'Cell proportion (%)', x=NULL)+
  theme_bw() +
  theme(axis.text.x = element_text(size=18, hjust = 0.5, face = 'bold'),
        axis.text.y = element_text(size = 18, hjust = 0.5, face = 'bold'),
        legend.text = element_text(size = 17, face = 'bold'),
        legend.key.size = unit(1.2, "cm"),
        axis.title.y = element_text(size = 20, hjust = 0.5, face = 'bold'),
        plot.title = element_text(size = 25, hjust = 0.5, face = 'bold')) +
  scale_fill_discrete(name=NULL)
#saveRDS(tubo_bar, 'tubo_bar.rds')

png("/mnt/8TB/users/shameed/shameed/Zheng/figures/site_response.png", width = 10, height = 8.5, units = 'in', res = 600)
p6
dev.off()  

####################################
tubo_bar<-sigData@meta.data %>% group_by(Groups, refined_celltype) %>%
  tally() %>% mutate(cell_percentage= n*100/sum(n)) #%>% 
#filter(`Response to chemotherapy`=='Sensitive' | `Response to chemotherapy`=='Resistant') 
p6<-ggplot(tubo_bar, aes(x=Groups, y=cell_percentage, fill= refined_celltype)) + 
  geom_bar(stat = 'identity', width = 0.2) + 
  labs(title = 'Cell distribution by site',
       y= 'Cell proportion (%)', x=NULL)+
  theme_bw() +
  theme(axis.text.x = element_text(size=18, hjust = 0.5, face = 'bold'),
        axis.text.y = element_text(size = 18, hjust = 0.5, face = 'bold'),
        legend.text = element_text(size = 17, face = 'bold'),
        legend.key.size = unit(1.2, "cm"),
        axis.title.y = element_text(size = 20, hjust = 0.5, face = 'bold'),
        plot.title = element_text(size = 25, hjust = 0.5, face = 'bold')) +
  scale_fill_discrete(name=NULL)
#saveRDS(tubo_bar, 'tubo_bar.rds')

png("/mnt/8TB/users/shameed/shameed/Zheng/figures/cell_site.png", width = 15, height = 8.5, units = 'in', res = 600)
p6
dev.off()  

##########################
table(sigData$stage)
sigData$stage <- str_replace(sigData$stage, 'Ⅵ', 'IV')
table(sigData$stage)
tubo_bar<-sigData@meta.data %>% group_by(stage, `Response to chemotherapy`) %>%
  tally() %>% mutate(cell_percentage= n*100/sum(n)) %>% 
  filter(`Response to chemotherapy`=='Sensitive' | `Response to chemotherapy`=='Resistant') 
p6<-ggplot(tubo_bar, aes(x=stage , y=cell_percentage, fill= `Response to chemotherapy`)) + 
  geom_bar(stat = 'identity', width = 0.2) + 
  labs(title = 'Response distribution by stage',
       y= 'Cell proportion (%)', x=NULL)+
  theme_bw() +
  theme(axis.text.x = element_text(size=18, hjust = 0.5, face = 'bold'),
        axis.text.y = element_text(size = 18, hjust = 0.5, face = 'bold'),
        legend.text = element_text(size = 17, face = 'bold'),
        legend.key.size = unit(1.2, "cm"),
        axis.title.y = element_text(size = 20, hjust = 0.5, face = 'bold'),
        plot.title = element_text(size = 25, hjust = 0.5, face = 'bold')) +
  scale_fill_discrete(name=NULL)
#saveRDS(tubo_bar, 'tubo_bar.rds')

png("/mnt/8TB/users/shameed/shameed/Zheng/figures/cell_response_stage.png", width = 15, height = 8.5, units = 'in', res = 600)
p6
dev.off()  

#######################################################multiplet proportions#########################################################
T_MacData <- subset(ZhengData, ident = c('CD4+ T', 'CD4+ T_Macrophage',
                                         'CD8+ T', 'CD8+ T_Macrophage', 
                                         'Macrophage'
) )

TMacDoub <- subset(ZhengData, ident = c( 'CD4+ T_Macrophage',
                                         'CD8+ T_Macrophage') )

############################
tubo_bar<-TMacDoub@meta.data %>% group_by(`Response to chemotherapy`, celltype_ulm) %>%
  tally() %>% mutate(cell_percentage= n*100/sum(n)) %>% 
  filter(`Response to chemotherapy`=='Sensitive' | `Response to chemotherapy`=='Resistant') 

p6<-ggplot(tubo_bar, aes(x=`Response to chemotherapy`, y=cell_percentage, fill= celltype_ulm)) + 
  geom_bar(stat = 'identity', width = 0.2) + 
  labs(title = 'T-Mac doublet distribution by response',
       y= 'Cell proportion (%)', x=NULL)+
  theme_bw() +
  theme(axis.text.x = element_text(size=18, hjust = 0.5, face = 'bold'),
        axis.text.y = element_text(size = 18, hjust = 0.5, face = 'bold'),
        legend.text = element_text(size = 17, face = 'bold'),
        legend.key.size = unit(1.2, "cm"),
        axis.title.y = element_text(size = 20, hjust = 0.5, face = 'bold'),
        plot.title = element_text(size = 25, hjust = 0.5, face = 'bold')) +
  scale_fill_discrete(name=NULL)
#saveRDS(tubo_bar, 'tubo_bar.rds')

png("/mnt/8TB/users/shameed/shameed/Zheng/figures/TMac_prop.png", width = 10, height = 8.5, units = 'in', res = 600)
p6
dev.off()  

##################
tubo_bar<-TMacDoub@meta.data %>% group_by(celltype_ulm, sample_site) %>%
  tally() %>% mutate(cell_percentage= n*100/sum(n)) #%>% 
# filter(`Response to chemotherapy`=='Sensitive' | `Response to chemotherapy`=='Resistant') 
p6<-ggplot(tubo_bar, aes(x=celltype_ulm, y=cell_percentage, fill= sample_site)) + 
  geom_bar(stat = 'identity', width = 0.2) + 
  labs(title = 'T-Mac doublet distribution by site',
       y= 'Cell proportion (%)', x=NULL)+
  theme_bw() +
  theme(axis.text.x = element_text(size=18, hjust = 0.5, face = 'bold'),
        axis.text.y = element_text(size = 18, hjust = 0.5, face = 'bold'),
        legend.text = element_text(size = 17, face = 'bold'),
        legend.key.size = unit(1.2, "cm"),
        axis.title.y = element_text(size = 20, hjust = 0.5, face = 'bold'),
        plot.title = element_text(size = 25, hjust = 0.5, face = 'bold')) +
  scale_fill_discrete(name=NULL)
#saveRDS(tubo_bar, 'tubo_bar.rds')

png("/mnt/8TB/users/shameed/shameed/Zheng/figures/TMac_site.png", width = 10, height = 8.5, units = 'in', res = 600)
p6
dev.off() 