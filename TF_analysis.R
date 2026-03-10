###########################################TF Activities########################
library('decoupleR')
library(tidyverse)
library(cowplot)
#remotes::install_github("cvarrichio/Matrix.utils")
library(Matrix.utils)
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(data.table)
ZhengData <- readRDS("/mnt/8TB/users/shameed/shameed/Zheng/ZhengData.rds")

T_MacData <- subset(ZhengData, ident = c('CD4+ T', 'CD4+ T_Macrophage',
                                         'CD8+ T', 'CD8+ T_Macrophage', 
                                         'Macrophage'
) )
T_MacData$response <- T_MacData$`Response to chemotherapy`
T_MacData<- subset(T_MacData, subset =response!= "/")
T_MacData$response_cell <- paste(T_MacData$celltype_ulm, T_MacData$response, sep = '_')
table(T_MacData$celltype_ulm)

##########################
Idents(T_MacData) <- T_MacData$celltype_ulm
cluster_names <- levels(Idents(T_MacData)) 
groups<- T_MacData@meta.data[, c("celltype_ulm", "response")]
head(T_MacData@assays$RNA@counts)
head(groups)
# Aggregate across cluster-sample groups
# transposing row/columns to have cell_ids as row names matching those of groups
aggr_counts <- aggregate.Matrix(t(T_MacData@assays$RNA@data), 
                                groupings = groups, fun = "sum") 
agg_data<- as.data.frame(t(aggr_counts)) #%>% head()

################run TF activity
net <- get_collectri(organism='human', split_complexes=FALSE)

acts <- run_ulm(mat=agg_data, net=net, .source='source', .target='target',
                .mor='mor', minsize = 5)
head(acts)

saveRDS(acts, '/mnt/8TB/users/shameed/shameed/Zheng/acts.rds')

###viewing the top 30

my_TF <- acts %>% filter(p_value <0.05) %>% group_by(condition) %>% 
  slice_max(abs(score), n=25) %>% pull(source)

# Transform to wide matrix
sample_acts_mat <- acts %>% filter(source %in% my_TF) %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()

#saveRDS(sample_acts_mat, 'tf_25_cluster.rds')

# Scale per sample
sample_acts_mat <- scale(sample_acts_mat)

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))
library(pheatmap)
# Plot
p1 <- pheatmap(sample_acts_mat, border_color = NA, color=my_color, 
               breaks = my_breaks, main = 'Transcriptional factor activity', treeheight_col = 0)

png("/mnt/8TB/users/shameed/shameed/Zheng/figures/TF_all.png", width = 10, height = 4.5, units = 'in', res = 600)
p1
dev.off() 

#############DEG TF##########
set.seed(062425)
##CD8_Mac
Idents(T_MacData)<- T_MacData$response_cell
deg <- FindMarkers(T_MacData, ident.1 = "CD8+ T_Macrophage_Sensitive",ident.2 = "CD8+ T_Macrophage_Resistant",min.pct = 0.25)
#####CD4_Mac
deg <- FindMarkers(T_MacData, ident.1 = "CD4+ T_Macrophage_Sensitive",ident.2 = "CD4+ T_Macrophage_Resistant",min.pct = 0.25)

########
deg<- deg[deg$p_val_adj <= 0.05,]
# Run ulm
contrast_acts <- run_ulm(mat=deg[, 'avg_log2FC', drop=FALSE], net=net, .source='source', .target='target',
                         .mor='mor', minsize = 5)

saveRDS(contrast_acts, '/mnt/8TB/users/shameed/shameed/Zheng/contrast_acts.CD8.rds')
saveRDS(contrast_acts, '/mnt/8TB/users/shameed/shameed/Zheng/contrast_acts.CD4.rds')

###viewing the top 30

# Filter top TFs in both signs
TF_up <- contrast_acts %>% slice_max(., n=15, order_by= score)
TF_down <- contrast_acts%>% slice_min(., n=15, order_by= score)
TF_temp <- rbind(TF_up, TF_down)

# Plot
p1<-ggplot(TF_temp, aes(x = reorder(source, score), y = score)) + 
  geom_bar(aes(fill = score), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Transcription Factors") + ggtitle('CD8+ T_Macrophage: Sensitive Vs Resistant')

png("/mnt/8TB/users/shameed/shameed/Zheng/figures/TF_DEG_CD8.png", width = 12, height = 4.5, units = 'in', res = 600)
p1
dev.off() 

p2<-ggplot(TF_temp, aes(x = reorder(source, score), y = score)) + 
  geom_bar(aes(fill = score), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Transcription Factors") + ggtitle('CD4+ T_Macrophage: Sensitive Vs Resistant')

png("/mnt/8TB/users/shameed/shameed/Zheng/figures/TF_DEG_CD4.png", width = 12, height = 4.5, units = 'in', res = 600)
p2
dev.off() 

png("/mnt/8TB/users/shameed/shameed/Zheng/figures/TF_DEG.png", width = 12, height = 7.5, units = 'in', res = 600)
p1/p2
dev.off() 

#########Doublet Vs Singlet##################
#######CD8 sensitive
deg <- FindMarkers(T_MacData, ident.1 = "CD8+ T_Macrophage_Sensitive",ident.2 = c("CD8+ T_Sensitive", 'Macrophage_Sensitive'),min.pct = 0.25)
#######CD8 Resistant
deg <- FindMarkers(T_MacData, ident.1 = "CD8+ T_Macrophage_Resistant",ident.2 = c("CD8+ T_Resistant", 'Macrophage_Resistant'),min.pct = 0.25)

########
deg<- deg[deg$p_val_adj <= 0.05,]
# Run ulm
contrast_acts <- run_ulm(mat=deg[, 'avg_log2FC', drop=FALSE], net=net, .source='source', .target='target',
                         .mor='mor', minsize = 5)

saveRDS(contrast_acts, '/mnt/8TB/users/shameed/shameed/Zheng/contrast_acts.CD8_doubSing_Sensitive.rds')
saveRDS(contrast_acts, '/mnt/8TB/users/shameed/shameed/Zheng/contrast_acts.CD8_doubSing_Resistant.rds')

###viewing the top 30

# Filter top TFs in both signs
TF_up <- contrast_acts %>% slice_max(., n=15, order_by= score)
TF_down <- contrast_acts%>% slice_min(., n=15, order_by= score)
TF_temp <- rbind(TF_up, TF_down)

# Plot
p1<-ggplot(TF_temp, aes(x = reorder(source, score), y = score)) + 
  geom_bar(aes(fill = score), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Transcription Factors") + ggtitle('CD8+ Doublet Vs Singlet (Sensitive)')

png("/mnt/8TB/users/shameed/shameed/Zheng/figures/TF_DEG_CD8_DoubSing_Sensit.png", width = 12, height = 4.5, units = 'in', res = 600)
p1
dev.off() 

p2<-ggplot(TF_temp, aes(x = reorder(source, score), y = score)) + 
  geom_bar(aes(fill = score), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Transcription Factors") + ggtitle('CD8+ Doublet Vs Singlet (Resistant)')

png("/mnt/8TB/users/shameed/shameed/Zheng/figures/TF_DEG_CD8_DoubSing_Resistant.png", width = 12, height = 4.5, units = 'in', res = 600)
p2
dev.off() 

png("/mnt/8TB/users/shameed/shameed/Zheng/figures/TF_DEG_CD8_DoubSing.png", width = 12, height = 7.5, units = 'in', res = 600)
p1/p2
dev.off() 

######################CD4
#######CD4 sensitive
deg <- FindMarkers(T_MacData, ident.1 = "CD4+ T_Macrophage_Sensitive",ident.2 = c("CD4+ T_Sensitive", 'Macrophage_Sensitive'),min.pct = 0.25)
#######CD4 Resistant
deg <- FindMarkers(T_MacData, ident.1 = "CD4+ T_Macrophage_Resistant",ident.2 = c("CD4+ T_Resistant", 'Macrophage_Resistant'),min.pct = 0.25)

########
deg<- deg[deg$p_val_adj <= 0.05,]
# Run ulm
contrast_acts <- run_ulm(mat=deg[, 'avg_log2FC', drop=FALSE], net=net, .source='source', .target='target',
                         .mor='mor', minsize = 5)

saveRDS(contrast_acts, '/mnt/8TB/users/shameed/shameed/Zheng/contrast_acts.CD4_doubSing_Sensitive.rds')
saveRDS(contrast_acts, '/mnt/8TB/users/shameed/shameed/Zheng/contrast_acts.CD4_doubSing_Resistant.rds')

###viewing the top 30

# Filter top TFs in both signs
TF_up <- contrast_acts %>% slice_max(., n=15, order_by= score)
TF_down <- contrast_acts%>% slice_min(., n=15, order_by= score)
TF_temp <- rbind(TF_up, TF_down)

# Plot
p3<-ggplot(TF_temp, aes(x = reorder(source, score), y = score)) + 
  geom_bar(aes(fill = score), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Transcription Factors") + ggtitle('CD4+ Doublet Vs Singlet (Sensitive)')

png("/mnt/8TB/users/shameed/shameed/Zheng/figures/TF_DEG_CD4_DoubSing_Sensit.png", width = 12, height = 4.5, units = 'in', res = 600)
p3
dev.off() 

p4<-ggplot(TF_temp, aes(x = reorder(source, score), y = score)) + 
  geom_bar(aes(fill = score), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Transcription Factors") + ggtitle('CD4+ Doublet Vs Singlet (Resistant)')

png("/mnt/8TB/users/shameed/shameed/Zheng/figures/TF_DEG_CD4_DoubSing_Resistant.png", width = 12, height = 4.5, units = 'in', res = 600)
p4
dev.off() 

png("/mnt/8TB/users/shameed/shameed/Zheng/figures/TF_DEG_CD4_DoubSing.png", width = 12, height = 7.5, units = 'in', res = 600)
p3/p4
dev.off() 

