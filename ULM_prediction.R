###############################ulm modeling####################################################
singletData <- readRDS("/mnt/8TB/users/shameed/shameed/Zheng/singletData.rds")
ZhengData <- readRDS("/mnt/8TB/users/shameed/shameed/Zheng/ZhengData.rds")
table(singletData$refined_celltype)
VlnPlot(singletData, features = c('CD14', 'CD86', 'CD4', 'CD8A', 'HLA-DRA'), group.by = 'maintypes_2', raster = F)
singletData$refined_celltype <- singletData$maintypes_2 %>%
  str_replace('Fibroblast', 'Stromal') %>% 
  str_replace('Other stromal cells', 'Stromal') %>%
  str_replace('Epithelial cells', 'Cancer cells') %>% 
  str_replace('Mesothelial cells', 'Stromal') 

#singletData@meta.data <- singletData@meta.data %>%
# mutate(refined_celltype = as.character(refined_celltype),
#       refined_celltype = if_else(refined_celltype %in% c('CD4+ T', 'CD8+ T'),
#                                 'T cells', refined_celltype))
table(singletData$refined_celltype)
class(singletData$refined_celltype)
sigData <- subset(singletData, subset = refined_celltype != 'Proliferative cells' & refined_celltype != 'HSC')
#Idents(singletData) <- singletData$celltype_main
#sigData <- subset(singletData, idents = c('B-cells', 'Epithelial cells',
#                                         'CD4+ T-cells', 'CD8+ T-cells',
#                                        'DC', 'Endothelial cells',
#                                       'Fibroblasts', 'Macrophages',
#                                      'Monocytes', 'NK cells'))

#DimPlot(sigData, reduction = 'umap', group.by = 'celltype_main', label=T)

table(sigData$refined_celltype)

saveRDS(singletData, '/mnt/8TB/users/shameed/shameed/Zheng/singletData.rds')

library(ULMnet)
set.seed(180225)
zheng_sig <- GetSignature(sigData, ident_col = sigData$refined_celltype)
saveRDS(zheng_sig, 'zheng_sig.rds')
saveRDS(zheng_sig, '/mnt/8TB/users/shameed/shameed/Zheng/zheng_sig.rds')

set.seed(180225)
my_scores_1 <- GetCellScores(seurat_obj = ZhengData[, 1:180000], signatures = zheng_sig, slot = 'data')
saveRDS(my_scores_1, 'my_scores_1.rds')
my_scores_2 <- GetCellScores(seurat_obj = ZhengData[, 180001:300000], signatures = zheng_sig, slot = 'data')
saveRDS(my_scores_2, 'my_scores_2.rds')
my_scores_3 <- GetCellScores(seurat_obj = ZhengData[, 300001: 400000], signatures = zheng_sig, slot = 'data')
saveRDS(my_scores_3, 'my_scores_3.rds')
my_scores_4 <- GetCellScores(seurat_obj = ZhengData[, 400001: ncol(ZhengData)], signatures = zheng_sig, slot = 'data')
saveRDS(my_scores_4, 'my_scores_4.rds')
#######
#my_scores_1 <- GetCellScores(seurat_obj = ZhengData[, 1:180000], signatures = zheng_sig, slot = 'data')
#saveRDS(my_scores_1, 'my_scores_11.rds')
#my_scores_2 <- GetCellScores(seurat_obj = ZhengData[, 180001:300000], signatures = zheng_sig, slot = 'data')
#saveRDS(my_scores_2, 'my_scores_22.rds')
#my_scores_3 <- GetCellScores(seurat_obj = ZhengData[, 300001: 400000], signatures = zheng_sig, slot = 'data')
#saveRDS(my_scores_3, 'my_scores_33.rds')
#my_scores_4 <- GetCellScores(seurat_obj = ZhengData[, 400001: ncol(ZhengData)], signatures = zheng_sig, slot = 'data')
#saveRDS(my_scores_4, 'my_scores_44.rds')
#my_scores <- rbind(my_scores_11, my_scores_22, my_scores_33, my_scores_44)

my_scores <- rbind(my_scores_1, my_scores_2, my_scores_3, my_scores_4)

my_ass <- GetCellAssignments(score_data = my_scores, cut_off = 1)
head(my_ass)
table(my_ass$count_ulm)
table(my_ass$celltype_ulm)
#my_scores_zheng <- my_scores
#my_ass_zheng <- my_ass
colnames(ZhengData@meta.data)
#ZhengData@meta.data<- ZhengData@meta.data %>% dplyr:: select(-c("count_ulm" , "celltype_ulm", "avg_pvalue", "avg_score"))
ZhengData <- AddMetaObject(ZhengData, cell_class_df = my_ass)
colnames(ZhengData@meta.data)
my_mult <- GetMultiplet(ZhengData)
my_mult_filt <- FilterMultiplet(ZhengData)
my_network_df <- GetNodeDF(mat = my_mult_filt$multSummaryFilt)
PlotNetwork(my_network_df)
p1<-PlotNetwork(my_network_df)
p1
#dir.create('plots')
png("/mnt/8TB/users/shameed/shameed/Zheng/figures/ulm_net.png", width = 15, height = 10.5, units = 'in', res = 600)
p1
dev.off()

pdf("/mnt/8TB/users/shameed/shameed/Zheng/figures/ulm_net.pdf",
    width = 15, height = 10.5)

p1
dev.off()

saveRDS(ZhengData, 'ZhengData.rds')

####Marker expression########################
table(ZhengData$celltype_ulm)
Idents(ZhengData) <- ZhengData$celltype_ulm

T_MacData <- subset(ZhengData, ident = c('CD4+ T', 'CD4+ T_Macrophage',
                                         'CD8+ T', 'CD8+ T_Macrophage', 
                                         'Macrophage'
) )

VlnPlot(T_MacData, features = c('CD4', 'CD8A', 'CD8B', 'CD3E', 'CD3D', 'CD3G'), raster = F)
VlnPlot(T_MacData, features = c('PARP14'), raster = F)

TMacDoub <- subset(ZhengData, ident = c( 'CD4+ T_Macrophage',
                                         'CD8+ T_Macrophage') )

p1<- VlnPlot(TMacDoub, features = c('CD4', 'CD8A', 'CD8B', 'CD3E', 'CD3D', 'CD3G', 'TRAC'), raster = F)

png("/mnt/8TB/users/shameed/shameed/Zheng/figures/T_mark.png", width = 15, height = 15.5, units = 'in', res = 600)
p1
dev.off()

p2<-VlnPlot(TMacDoub, features = c('CD14', 'CD68', 'MARCO', 'CD163', 'HLA-DRA', 'CD74', 'FCGR3A'), raster = F)
png("/mnt/8TB/users/shameed/shameed/Zheng/figures/Mac_mark.png", width = 15, height = 15.5, units = 'in', res = 600)
p2
dev.off()

p3 <- DoHeatmap(TMacDoub, features = c('CD3D', 'CD3E', 'CD3G', 'CD4', 'CD8A','CD8B', 'TRAC', 
                                       'CD14', 'CD68', 'MARCO', 'CD163', 'HLA-DRA', 'CD74', 'FCGR3A'), 
                slot = 'data',label = F)

png("/mnt/8TB/users/shameed/shameed/Zheng/figures/TMac_heat.png", width = 10, height = 5.5, units = 'in', res = 600)
p3
dev.off()

saveRDS(T_MacData, '/mnt/8TB/users/shameed/shameed/Zheng/T_MacData.rds')
saveRDS(TMacDoub, '/mnt/8TB/users/shameed/shameed/Zheng/TMacDoub.rds')

ZhengData@meta.data$`Response to chemotherapy` %>% table()
table(T_MacData$celltype_ulm, T_MacData$`Response to chemotherapy`)

