#############################Pathway analysis doublets#####################################

ZhengData <- readRDS("/mnt/8TB/users/shameed/shameed/Zheng/ZhengData.rds")

T_MacData <- subset(ZhengData, ident = c('CD4+ T', 'CD4+ T_Macrophage',
                                         'CD8+ T', 'CD8+ T_Macrophage', 
                                         'Macrophage'
) )

set.seed(020725)
##CD8_Mac
Idents(T_MacData)<- T_MacData$response_cell
deg <- FindMarkers(T_MacData, ident.1 = "CD8+ T_Macrophage_Sensitive",ident.2 = "CD8+ T_Macrophage_Resistant",min.pct = 0.25)
#####CD4_Mac
deg <- FindMarkers(T_MacData, ident.1 = "CD4+ T_Macrophage_Sensitive",ident.2 = "CD4+ T_Macrophage_Resistant",min.pct = 0.25)

########
#deg$stat <- deg$avg_log2FC * (-log10(deg$p_val_adj))

saveRDS(deg, '/mnt/8TB/users/shameed/shameed/Zheng/deg.CD8_Mac.rds')
saveRDS(deg, '/mnt/8TB/users/shameed/shameed/Zheng/deg.CD4_Mac.rds')

############
library(msigdbr)
library(fgsea)
C2_REACTOME <-msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "REACTOME")

C.symbol<-C2_REACTOME%>% dplyr:: select(c(gs_name, gene_symbol)) %>% group_by(gs_name) %>%
  summarise(all.genes= list(gene_symbol)) %>% deframe()

gsea_gene<- deg %>% #### repeat by replacing with DEGs of the other groups
  rownames_to_column(var = 'genes') %>% 
  dplyr::select(c( genes, avg_log2FC))

gsea_gene<- gsea_gene[order(-gsea_gene$avg_log2FC),] #arrange(desc(gsea_gene$log2FoldChange))
gene_list<- gsea_gene$avg_log2FC
#gene_list<- jitter(gsea_gene$avg_log2FC, factor = 0.01) ##if there is error due to many matches in avg_log2FC
names(gene_list)<- gsea_gene$genes
gsea_path<- fgseaSimple(pathways = C.symbol, stats = gene_list, nperm = 1000)
gsea_path <- gsea_path%>% filter(pval <0.05)

CD8_Mac_reac <- gsea_path
CD8_Mac_reac$pathway <- gsub('REACTOME_', '', CD8_Mac_reac$pathway)
CD8_Mac_reac$pathway <- gsub('REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_', '', CD8_Mac_reac$pathway)
CD8_Mac_reac$pathway <- gsub('AND_EIFS_AND_SUBSEQUENT_BINDING_TO_43S_', '', CD8_Mac_reac$pathway)

#CD8_Mac_reac$group <- 'CD8_Mac'
p1 <- ggplot(CD8_Mac_reac[1:60], aes(x= reorder(pathway, NES), y= NES, fill= pval)) + 
  geom_col() + coord_flip() + labs(y='NORMALISED ENRICHMENT SCORES',
                                   x= 'REACTOME PATHWAYS') + ggtitle('CD8_Macrophage Doublet: Sensitive vs Resistant')
############
CD4_Mac_reac <- gsea_path
CD4_Mac_reac$pathway <- gsub('REACTOME_', '', CD4_Mac_reac$pathway)
CD4_Mac_reac$pathway <- gsub('AND_EIFS_AND_SUBSEQUENT_BINDING_TO_43S_', '', CD4_Mac_reac$pathway)

png("/mnt/8TB/users/shameed/shameed/Zheng/figures/pathways_CD8_Mac.png", width = 15, height = 10.5, units = 'in', res = 600)
p1
dev.off()

p2 <- ggplot(CD4_Mac_reac[1:60], aes(x= reorder(pathway, NES), y= NES, fill= pval)) + 
  geom_col() + coord_flip() + labs(y='NORMALISED ENRICHMENT SCORES',
                                   x= 'REACTOME PATHWAYS') + ggtitle('CD4_Macrophage Doublet: Sensitive vs Resistant')

png("/mnt/8TB/users/shameed/shameed/Zheng/figures/pathways_CD4_Mac.png", width = 15, height = 10.5, units = 'in', res = 600)
p2
dev.off()

#################doublet vs singlet#############
######CD8 sensitive
deg <- FindMarkers(T_MacData, ident.1 = "CD8+ T_Macrophage_Sensitive",ident.2 = c("CD8+ T_Sensitive", 'Macrophage_Sensitive'),min.pct = 0.25)
#######CD8 Resistant
deg <- FindMarkers(T_MacData, ident.1 = "CD8+ T_Macrophage_Resistant",ident.2 = c("CD8+ T_Resistant", 'Macrophage_Resistant'),min.pct = 0.25)
######CD4 sensitive
deg <- FindMarkers(T_MacData, ident.1 = "CD4+ T_Macrophage_Sensitive",ident.2 = c("CD4+ T_Sensitive", 'Macrophage_Sensitive'),min.pct = 0.25)
#######CD4 Resistant
deg <- FindMarkers(T_MacData, ident.1 = "CD4+ T_Macrophage_Resistant",ident.2 = c("CD4+ T_Resistant", 'Macrophage_Resistant'),min.pct = 0.25)

########
#deg$p_val_adj <- pmax(deg$p_val_adj, 1e-300)
#deg$stat <- deg$avg_log2FC * (-log10(deg$p_val_adj))

gsea_gene<- deg %>% #### repeat by replacing with DEGs of the other groups
  rownames_to_column(var = 'genes') %>% 
  dplyr::select(c( genes, avg_log2FC))

gsea_gene<- gsea_gene[order(-gsea_gene$avg_log2FC),] #arrange(desc(gsea_gene$log2FoldChange))

gene_list<- gsea_gene$avg_log2FC
#gene_list<- jitter(gsea_gene$avg_log2FC, factor = 0.01) ##if there is error due to many matches in avg_log2FC
names(gene_list)<- gsea_gene$genes

gsea_path<- fgseaSimple(pathways = C.symbol, stats = gene_list, nperm = 1000)
gsea_path <- gsea_path%>% filter(pval <0.05)

CD8_DoubSing_path_sensi <- gsea_path
CD8_DoubSing_path_sensi$group <- 'CD8 Sensitive'
CD8_DoubSing_path_resis <- gsea_path
CD8_DoubSing_path_resis$group <- 'CD8 Resistant'
CD4_DoubSing_path_sensi <- gsea_path
CD4_DoubSing_path_sensi$group <- 'CD4 Sensitive'
CD4_DoubSing_path_resis <- gsea_path
CD4_DoubSing_path_resis$group <- 'CD4 Resistant'
ggplot(CD4_DoubSing_path_resis[1:50], aes(x= reorder(pathway, NES), y= NES, fill= pval)) + 
  geom_col() + coord_flip() + labs(y='NORMALISED ENRICHMENT SCORES',
                                   x= 'REACTOME PATHWAYS') #+ ggtitle('CD4_Macrophage Doublet: Sensitive vs Resistant')

Pathways_Reac <- rbind(CD8_DoubSing_path_sensi,CD8_DoubSing_path_resis,
                       CD4_DoubSing_path_sensi, CD4_DoubSing_path_resis)
Pathways_filt<- Pathways_Reac %>% group_by(group) %>% slice_max(n=30, order_by = NES) %>% pull(pathway)
Pathways_Reac <- Pathways_Reac[Pathways_Reac$pathway %in% unique(Pathways_filt),]

Pathways_Reac$pathway <- gsub('REACTOME_', '', Pathways_Reac$pathway)
Pathways_Reac$pathway <- gsub('_RECEPTORS_SHR_IN_THE_PRESENCE_OF_LIGAND', '', Pathways_Reac$pathway)
Pathways_Reac$pathway <- gsub("_BY_CHEMIOSMOTIC_COUPLING_AND_HEAT_PRODUCTION_BY_UNCOUPLING_PROTEINS", '', Pathways_Reac$pathway)
Pathways_Reac$pathway <- gsub('AND_EIFS_AND_SUBSEQUENT_BINDING_TO_43S_', '', Pathways_Reac$pathway)
Pathways_Reac$pathway <- gsub('REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_', '', Pathways_Reac$pathway)


Pathways_Reac$group <- factor(Pathways_Reac$group,
                              levels = c('CD8 Sensitive', 'CD4 Sensitive',
                                         'CD8 Resistant', 'CD4 Resistant'))

# Plot
p3<-ggplot(Pathways_Reac, aes(x = group, y = pathway)) + 
  geom_point(aes(fill = NES, size = -log10(padj)), shape = 21) + 
  labs(x = "", y = "REACTOME PATHWAYS", fill = "NES", size = "-log10(padj)")  + 
  scale_fill_gradient2(low="#3B4992FF", mid="white", high="#EE0000FF", oob = scales::squish) +
  theme(
    legend.key = element_blank(), 
    panel.background = element_blank(), 
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), 
    legend.position = "right", 
    axis.text = element_text(face = "bold", size = 9), 
    legend.text = element_text(face = "bold", size = 10), 
    legend.title = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 90),
    title = element_text(face = "bold", size = 12)
  ) + 
  ggtitle('Doublets vs Singlets')


png("/mnt/8TB/users/shameed/shameed/Zheng/figures/pathways_DoubSing.png", width = 13, height = 12.5, units = 'in', res = 600)
p3
dev.off()
