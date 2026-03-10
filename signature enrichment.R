######################################signature enrichment#####################
ZhengData <- readRDS("/mnt/8TB/users/shameed/shameed/Zheng/ZhengData.rds")

T_MacData <- subset(ZhengData, ident = c('CD4+ T', 'CD4+ T_Macrophage',
                                         'CD8+ T', 'CD8+ T_Macrophage', 
                                         'Macrophage'
) )

TMacDoub <- subset(ZhengData, ident = c( 'CD4+ T_Macrophage',
                                         'CD8+ T_Macrophage') )

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


gene_table <- as.data.frame(
  lapply(gene.sets, function(x) {
    length(x) <- max(lengths(gene.sets))
    x
  })
)

################################signature: doublet and singlet###############################
library(UCell)
T_MacData <- AddModuleScore_UCell(T_MacData, features = gene.sets, assay = "RNA", name = NULL)

T_MacData$IFN <- NULL
T_MacData$HeatShock <- NULL
T_MacData$Tcell.proliferation <- NULL
T_MacData$cellCycle.G2M <- NULL

VlnPlot(T_MacData, features = names(gene.sets),stack = T, flip = T) + NoLegend() #+ coord_flip()

T_MacData$response <- T_MacData$`Response to chemotherapy`

T_MacData$response_cell <- paste(T_MacData$celltype_ulm, T_MacData$response, sep = '_')

p1<-VlnPlot(T_MacData, features = names(gene.sets), group.by = 'response', stack = T, flip = T) + NoLegend()
p2<-VlnPlot(T_MacData, features = names(gene.sets), group.by = 'response_cell', pt.size = 0, stack = T, flip = T) + NoLegend()
p3<-VlnPlot(T_MacData, features = names(gene.sets), pt.size = 0, stack = T, flip = T) + NoLegend()

p1<-VlnPlot(T_MacData, features = names(gene.sets), group.by = 'response', 
            stack = T, flip = T) + NoLegend() +
  theme(
    axis.text.x = element_text(size = 25, face = "bold"),
    axis.text.y = element_text(size = 25, face = 'bold'),
    axis.title.x = element_text(size = 25, face = "bold"),
    axis.title.y = element_text(size = 25, face = "bold"),
    plot.title = element_text(size = 25, face = "bold"),
    strip.text = element_text(size = 25, face = "bold")
  )


p2<-VlnPlot(T_MacData, features = names(gene.sets), group.by = 'response_cell', 
            stack = T, flip = T) + NoLegend() +
  theme(
    axis.text.x = element_text(size = 25, face = "bold"),
    axis.text.y = element_text(size = 25, face = 'bold'),
    axis.title.x = element_text(size = 25, face = "bold"),
    axis.title.y = element_text(size = 25, face = "bold"),
    plot.title = element_text(size = 25, face = "bold"),
    strip.text = element_text(size = 25, face = "bold")
  )

p3<-VlnPlot(T_MacData, features = names(gene.sets),
            stack = T, flip = T) + NoLegend() +
  theme(
    axis.text.x = element_text(size = 25, face = "bold"),
    axis.text.y = element_text(size = 25, face = 'bold'),
    axis.title.x = element_text(size = 25, face = "bold"),
    axis.title.y = element_text(size = 25, face = "bold"),
    plot.title = element_text(size = 25, face = "bold"),
    strip.text = element_text(size = 25, face = "bold")
  )

png("/mnt/8TB/users/shameed/shameed/Zheng/figures/signature_response_all.png", width = 15, height = 20.5, units = 'in', res = 600)
p1
dev.off()
png("/mnt/8TB/users/shameed/shameed/Zheng/figures/signature_response_cell_all.png", width = 25, height = 20.5, units = 'in', res = 600)
p2
dev.off()
png("/mnt/8TB/users/shameed/shameed/Zheng/figures/signature_cell_all.png", width = 15, height = 20.5, units = 'in', res = 600)
p3
dev.off()
