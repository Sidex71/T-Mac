####################################################projecTILs######################################################################################
singletData <- readRDS("/mnt/8TB/users/shameed/shameed/Zheng/singletData.rds")
ZhengData <- readRDS("/mnt/8TB/users/shameed/shameed/Zheng/ZhengData.rds")
T_MacData <- subset(ZhengData, ident = c('CD4+ T', 'CD4+ T_Macrophage',
                                         'CD8+ T', 'CD8+ T_Macrophage', 
                                         'Macrophage'
) )

TMacDoub <- subset(ZhengData, ident = c( 'CD4+ T_Macrophage',
                                         'CD8+ T_Macrophage') )

library(remotes)

remotes::install_github("carmonalab/STACAS")
remotes::install_github("carmonalab/ProjecTILs")
library(ProjecTILs)
ref <- load.reference.map()
ref$Tcell.exhaustion
#data(query_example_seurat)
#ref <- CD8T_human_ref_v1
query.projected <- Run.ProjecTILs(query = TMacDoub, ref=ref)
saveRDS(query.projected, 'query.projected.rds')
p1<-plot.projection(ref, query.projected, linesize = 0.5, pointsize = 0.5)
png("/mnt/8TB/users/shameed/shameed/Zheng/figures/umap_ref.png", width = 10, height = 5.5, units = 'in', res = 600)
p1
dev.off()
p2<-plot.statepred.composition(ref, query.projected, metric = "Percent")
png("/mnt/8TB/users/shameed/shameed/Zheng/figures/ref_bar.png", width = 10, height = 5.5, units = 'in', res = 600)
p2
dev.off()

genes4radar = c("Cd4", "Cd8a", "Tcf7", "Ccr7", "Sell", 'Gzmm', "Gzmb", "Gzmk",'Ifng',
                "Havcr2", "Tox", "Mki67", 'Il2ra', 'Tigit', 'Ctla4', 'Lag3',"Foxp3","Pdcd1")

#genes4radar = c("FOXP3", "CD4", "CD8A", "TCF7", "CCR7", "SELL", 'GZMM', "GZMB", "GZMK", "PDCD1",
#               "HAVCR2", "TOX", "MKI67", 'IL2RA', 'IFNG')

p3<-plot.states.radar(ref, query = query.projected, genes4radar = genes4radar, min.cells = 15)

png("/mnt/8TB/users/shameed/shameed/Zheng/figures/Radar_TMAC.png", width = 15, height = 10.5, units = 'in', res = 600)
p3
dev.off()

##################################compare conditions#######################################################################
CD4<- subset(query.projected, subset =celltype_ulm== "CD4+ T_Macrophage")
CD8 <- subset(query.projected, subset =celltype_ulm== "CD8+ T_Macrophage")

p4<-plot.states.radar(ref, query = list(CD4_Macrophage =CD4 , CD8_Macrophage = CD8),genes4radar = genes4radar)
png("/mnt/8TB/users/shameed/shameed/Zheng/figures/Radar_TMAC_sep.png", width = 15, height = 10.5, units = 'in', res = 600)
p4
dev.off()

###############
query.projected$response <- query.projected$`Response to chemotherapy`
Resist<- subset(query.projected, subset =response== "Resistant")
Sensit<- subset(query.projected, subset =response== "Sensitive")


p5<-plot.states.radar(ref, query = list(Sensitive =Sensit , Resistant = Resist),genes4radar = genes4radar)
png("/mnt/8TB/users/shameed/shameed/Zheng/figures/Radar_TMAC_res.png", width = 15, height = 10.5, units = 'in', res = 600)
p5
dev.off()

query.projected$response_cell <- paste(query.projected$celltype_ulm, query.projected$response, sep = '_')
Resist_CD4<- subset(query.projected, subset =response_cell=="CD4+ T_Macrophage_Resistant")
Resist_CD8<- subset(query.projected, subset =response_cell=="CD8+ T_Macrophage_Resistant")
Sensit_CD4<- subset(query.projected, subset =response_cell=="CD4+ T_Macrophage_Sensitive")
Sensit_CD8<- subset(query.projected, subset =response_cell=="CD8+ T_Macrophage_Sensitive")

p6<-plot.states.radar(ref, query = list(CD4_Macrophage_Sensitive =Sensit_CD4 , CD4_Macrophage_Resistant = Resist_CD4,
                                        CD8_Macrophage_Sensitive =Sensit_CD8 , CD8_Macrophage_Resistant = Resist_CD8),genes4radar = genes4radar)

png("/mnt/8TB/users/shameed/shameed/Zheng/figures/Radar_TMAC_res_cell.png", width = 15, height = 10.5, units = 'in', res = 600)
p6
dev.off()

#######################discriminant genes##########################
discriminantGenes <- find.discriminant.genes(ref = ref, query =Sensit_CD8 , CD8_Resistant = Resist_CD8,
                                             state = "CD8_EffectorMemory")
saveRDS(discriminantGenes, 'discriminantGenes_CD8Eff.rds')
head(discriminantGenes, n = 10)
library(EnhancedVolcano)
p7<-EnhancedVolcano(discriminantGenes, lab = rownames(discriminantGenes), x = "avg_log2FC",
                    y = "p_val", pCutoff = 1e-09, FCcutoff = 0.5, labSize = 5, legendPosition = "none",
                    drawConnectors = F, title = "C8_Mac Sentive vs. Resistant (CD8 effector memory)")

png("/mnt/8TB/users/shameed/shameed/Zheng/figures/volcano_CD8EffMemo.png", width = 18, height = 15.5, units = 'in', res = 600)
p7
dev.off()
###############################################################
discriminantGenes <- find.discriminant.genes(ref = ref, query =Sensit_CD8 , CD8_Resistant = Resist_CD8,
                                             state = "CD8_NaiveLike")
saveRDS(discriminantGenes, 'discriminantGenes_CD8naiveLike.rds')
head(discriminantGenes, n = 10)
library(EnhancedVolcano)
p8<-EnhancedVolcano(discriminantGenes, lab = rownames(discriminantGenes), x = "avg_log2FC",
                    y = "p_val", pCutoff = 1e-09, FCcutoff = 0.5, labSize = 5, legendPosition = "none",
                    drawConnectors = F, title = "C8_Mac Sentive vs. Resistant (CD8 Naive-like)")

png("/mnt/8TB/users/shameed/shameed/Zheng/figures/volcano_CD8naiveLike.png", width = 18, height = 15.5, units = 'in', res = 600)
p8
dev.off()

#######################################
discriminantGenes_CD4Naive <- find.discriminant.genes(ref = ref, query =Sensit_CD4 , CD4_Resistant = Resist_CD4,
                                                      state = "CD8_NaiveLike")
saveRDS(discriminantGenes, '/mnt/8TB/users/shameed/shameed/Zheng/discriminantGenes_CD4naiveLike.rds')
head(discriminantGenes, n = 10)
library(EnhancedVolcano)
p8<-EnhancedVolcano(discriminantGenes_CD4Naive, lab = rownames(discriminantGenes_CD4Naive), x = "avg_log2FC",
                    y = "p_val", pCutoff = 1e-09, FCcutoff = 0.5, labSize = 5, legendPosition = "none",
                    drawConnectors = F, title = "CD4_Mac Sentive vs. Resistant (CD4 Naive-like)")

png("/mnt/8TB/users/shameed/shameed/Zheng/figures/volcano_CD4naiveLike.png", width = 18, height = 15.5, units = 'in', res = 600)
p8
dev.off()
######################################signature enrichment#####################
TMacDoub <- TMacDoub |>
  FindVariableFeatures() |>
  ScaleData() |>
  RunPCA(npcs = 20) |>
  RunUMAP(dims = 1:10)
saveRDS(TMacDoub, 'TMacDoub.rds')
ElbowPlot(TMacDoub)
DimPlot(TMacDoub)

#remotes::install_github("carmonalab/SignatuR")
library(SignatuR)
programs <-SignatuR::GetSignature(SignatuR$Mm$Programs)
gene.sets <- lapply(programs, str_to_upper)
gene.sets$Tcell.cytotoxicity <- c(gene.sets$Tcell.cytotoxicity, 'GZMM', 'GZMK', 'TNF','GNLY', 'IFNG', 'FASLG')
gene.sets$immune.checkpoint<- gene.sets$Tcell.exhaustion 
gene.sets$M1.macrophage <- c('TNF', 'IL6', 'IL12A', 'IL1B', 'CD80', 'CD86', 'NOS2')
gene.sets$M2.macrophage <- c('IL10', 'TGFBI', 'MRC1','CD163','ARG1', 'PPARG')
gene.sets$Tcell.exhaustion <- NULL
gene.sets$Tcell.exhaustion <- c('ENTPD1','LAYN','ITGAE','BATF', 'TOX')
gene.sets$Tcell.proliferation <- c('PCNA', 'MCM2', 'IL2RA', 'CD71', 'MKI67')
library(UCell)
TMacDoub <- AddModuleScore_UCell(TMacDoub, features = gene.sets, assay = "RNA", name = NULL)
#ref <- AddModuleScore_UCell(ref, features = sidex, assay = "RNA", name = NULL)
#plot.states.radar(ref, query = query.projected, meta4radar = names(gene.sets),min.cells = 15)
TMacDoub$M2.macrophage %>% head()
TMacDoub$Tcell.exhaustion

TMacDoub$IFN <- NULL
TMacDoub$HeatShock <- NULL
TMacDoub$Tcell.proliferation <- NULL
TMacDoub$cellCycle.G2M <- NULL
#VlnPlot(TMacDoub, features = c('cellCycle.G1S', 'cellCycle.G2M', 'M1.macrophage'), pt.size=0)
VlnPlot(TMacDoub, features = names(gene.sets),stack = T, flip = F) + NoLegend() #+ coord_flip()

TMacDoub$response <- TMacDoub$`Response to chemotherapy`
TMacDoub <- subset(TMacDoub, subset =response!= "/")
TMacDoub$response_cell <- paste(TMacDoub$celltype_ulm, TMacDoub$response, sep = '_')
p1<-VlnPlot(TMacDoub, features = names(gene.sets), group.by = 'response', 
            stack = T, flip = T) + NoLegend() +
  theme(
    axis.text.x = element_text(size = 25, face = "bold"),
    axis.text.y = element_text(size = 25, face = 'bold'),
    axis.title.x = element_text(size = 25, face = "bold"),
    axis.title.y = element_text(size = 25, face = "bold"),
    plot.title = element_text(size = 25, face = "bold"),
    strip.text = element_text(size = 25, face = "bold")
  )


p2<-VlnPlot(TMacDoub, features = names(gene.sets), group.by = 'response_cell', 
            stack = T, flip = T) + NoLegend() +
  theme(
    axis.text.x = element_text(size = 25, face = "bold"),
    axis.text.y = element_text(size = 25, face = 'bold'),
    axis.title.x = element_text(size = 25, face = "bold"),
    axis.title.y = element_text(size = 25, face = "bold"),
    plot.title = element_text(size = 25, face = "bold"),
    strip.text = element_text(size = 25, face = "bold")
  )

p3<-VlnPlot(TMacDoub, features = names(gene.sets),
            stack = T, flip = T) + NoLegend() +
  theme(
    axis.text.x = element_text(size = 25, face = "bold"),
    axis.text.y = element_text(size = 25, face = 'bold'),
    axis.title.x = element_text(size = 25, face = "bold"),
    axis.title.y = element_text(size = 25, face = "bold"),
    plot.title = element_text(size = 25, face = "bold"),
    strip.text = element_text(size = 25, face = "bold")
  )


png("/mnt/8TB/users/shameed/shameed/Zheng/figures/signature_response.png", width = 15, height = 20.5, units = 'in', res = 600)
p1
dev.off()
png("/mnt/8TB/users/shameed/shameed/Zheng/figures/signature_response_cell.png", width = 22, height = 20.5, units = 'in', res = 600)
p2
dev.off()
png("/mnt/8TB/users/shameed/shameed/Zheng/figures/signature_cell.png", width = 15, height = 20.5, units = 'in', res = 600)
p3
dev.off()

#########################pathway analysis- discriminant genes#################################################
library(enrichR)
#BiocManager::install('pathview')
library(pathview)
library(org.Mm.eg.db)
#BiocManager::install('org.Hs.eg.db')
library(org.Hs.eg.db)
library(EBImage)
#BiocManager::install('EBImage')
library(EBImage)
library('clusterProfiler')
#install.packages("msigdbr")
library(msigdbr)
#BiocManager::install('fgsea')
library(fgsea)

H<- msigdbr::msigdbr(species = 'Homo sapiens', category = 'H')
#C2_KEGG <-msigdbr:: msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")

C2_REACTOME <-msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "REACTOME")
#C5 <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
#C7 <- msigdbr(species = "Homo sapiens", category = "C7") #immunologic signature gene sets
#C8 <- msigdbr(species = "Homo sapiens", category = "C8") #cell type signature gene sets
H.symbol<- C2_REACTOME %>% dplyr:: select(c(gs_name, gene_symbol)) %>% group_by(gs_name) %>%
  summarise(all.genes= list(gene_symbol)) %>% deframe()

###############################convert mouse to human genes#
library(homologene)

mouse_genes <- rownames(discriminantGenes_CD8Eff)
conversion <- homologene(mouse_genes, inTax = 10090, outTax = 9606)
print(conversion)
conversion <- conversion[, 1:2]
colnames(conversion) <- c('Gene', 'human_gene')
#######################################
gsea_gene<- discriminantGenes_CD8Eff %>% rownames_to_column('Gene') %>% 
  filter(p_val <=0.05) %>% dplyr::select(c(Gene, avg_log2FC))
gsea_gene<- gsea_gene[order(-gsea_gene$avg_log2FC),] #arrange(desc(gsea_gene$log2FoldChange))
#gsea_gene <- right_join(gsea_gene, conversion, 'Gene')
gene_list<- gsea_gene$avg_log2FC
#gene_list<- jitter(gsea_gene$avg_log2FC, factor = 0.01) ##if there is error due to many matches in avg_log2FC
names(gene_list)<- gsea_gene$Gene
names(gene_list) <- str_to_upper(names(gene_list))
gsea_path<- fgseaSimple(pathways = H.symbol, stats = gene_list, nperm = 1000)
gsea_path$pathway<- str_replace(gsea_path$pathway, 'REACTOME_', '')
gsea_path <- gsea_path%>% filter(pval <0.05)

p1<-ggplot(gsea_path, aes(x= reorder(pathway, NES), y= NES, fill= pval)) + 
  geom_col() + coord_flip() + labs(y='normalised enrichment scores',
                                   x= 'Reactome pathways') + ggtitle('Sensitive vs Resistant (CD8 Effector Memory)')
saveRDS(gsea_path, 'gsea_pathways_CD8EffMem.rds')
png("/mnt/8TB/users/shameed/shameed/Zheng/figures/pathways_CD8EffMem.png", width = 15, height = 10.5, units = 'in', res = 600)
p1
dev.off()

###########################
gsea_gene<- discriminantGenes_CD8naiveLike %>% rownames_to_column('Gene') %>% 
  filter(p_val <=0.05) %>% dplyr::select(c(Gene, avg_log2FC))
gsea_gene<- gsea_gene[order(-gsea_gene$avg_log2FC),] #arrange(desc(gsea_gene$log2FoldChange))
#gsea_gene <- right_join(gsea_gene, conversion, 'Gene')
gene_list<- gsea_gene$avg_log2FC
#gene_list<- jitter(gsea_gene$avg_log2FC, factor = 0.01) ##if there is error due to many matches in avg_log2FC
names(gene_list)<- gsea_gene$Gene
names(gene_list) <- str_to_upper(names(gene_list))
gsea_path<- fgseaSimple(pathways = H.symbol, stats = gene_list, nperm = 1000)
gsea_path$pathway<- str_replace(gsea_path$pathway, 'REACTOME_', '')
gsea_path <- gsea_path%>% filter(pval <0.05)

p1<-ggplot(gsea_path, aes(x= reorder(pathway, NES), y= NES, fill= pval)) + 
  geom_col() + coord_flip() + labs(y='normalised enrichment scores',
                                   x= 'Reactome pathways') + ggtitle('Sensitive vs Resistant (CD8 Naive-like)')
saveRDS(gsea_path, 'gsea_pathways_CD8Naive.rds')
png("/mnt/8TB/users/shameed/shameed/Zheng/figures/pathways_CD8Naive.png", width = 15, height = 10.5, units = 'in', res = 600)
p1
dev.off()

#######################
gsea_gene<- discriminantGenes_CD4Naive %>% rownames_to_column('Gene') %>% 
  filter(p_val <=0.05) %>% dplyr::select(c(Gene, avg_log2FC))
gsea_gene<- gsea_gene[order(-gsea_gene$avg_log2FC),] #arrange(desc(gsea_gene$log2FoldChange))
#gsea_gene <- right_join(gsea_gene, conversion, 'Gene')
gene_list<- gsea_gene$avg_log2FC
#gene_list<- jitter(gsea_gene$avg_log2FC, factor = 0.01) ##if there is error due to many matches in avg_log2FC
names(gene_list)<- gsea_gene$Gene
names(gene_list) <- str_to_upper(names(gene_list))
gsea_path<- fgseaSimple(pathways = H.symbol, stats = gene_list, nperm = 1000)
gsea_path$pathway<- str_replace(gsea_path$pathway, 'REACTOME_', '')
gsea_path <- gsea_path%>% filter(pval <0.05)

p1<-ggplot(gsea_path, aes(x= reorder(pathway, NES), y= NES, fill= pval)) + 
  geom_col() + coord_flip() + labs(y='normalised enrichment scores',
                                   x= 'Reactome pathways') + ggtitle('Sensitive vs Resistant (CD8 Naive-like)')

saveRDS(gsea_path, '/mnt/8TB/users/shameed/shameed/Zheng/gsea_pathways_CD4NaiveLike.rds')

png("/mnt/8TB/users/shameed/shameed/Zheng/figures/pathways_CD4NaiveLike.png", width = 15, height = 10.5, units = 'in', res = 600)
p1
dev.off()