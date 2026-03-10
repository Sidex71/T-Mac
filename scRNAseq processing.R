###from processed count matrx#############
raw_object <- readRDS("/mnt/8TB/users/shameed/shameed/Zheng/raw_object.rds")
raw_object@images <- list()
raw_object@assays$RNA <- CreateAssayObject(raw_object@assays$RNA@counts)
raw_object@assays$RNA
raw_object = UpdateSeuratObject(object = raw_object)
raw_object <- NormalizeData(raw_object)
raw_object <- FindVariableFeatures(raw_object)
raw_object <- ScaleData(raw_object)
#raw_object@assays
raw_object <- RunPCA(raw_object)
ElbowPlot(raw_object, ndims = 50)
raw_object<-FindNeighbors(raw_object,dims = 1:40)
raw_object<-FindClusters(raw_object,resolution= 0.3)
raw_object<-RunUMAP(raw_object, dims = 1:40)
DimPlot(raw_object, reduction = 'umap', group.by = 'Groups', label=T)
raw_obj$Groups %>% table()
VlnPlot(ZhengData, features = c('ADIPOQ'),group.by ='celltype_main', raster = F)
saveRDS(raw_object, 'singletData.rds')
VlnPlot(singletData,features = c("nFeature_RNA",
                                 "nCount_RNA",
                                 "percent.mt"),ncol = 3)
DimPlot(singletData, reduction = 'umap', group.by = 'maintypes_2', label=T)

singletData@reductions
########harmony correction
singletData <- RunHarmony(singletData, 'Patients')
singletData@reductions
DimPlot(singletData, reduction = 'harmony', group.by = 'Patients', label=T)
singletData<-FindNeighbors(singletData,dims = 1:40,reduction = 'harmony')
singletData<-FindClusters(singletData,resolution= 0.3, reduction = 'harmony')
singletData<-RunUMAP(singletData, dims = 1:40, reduction = 'harmony')
DimPlot(singletData, reduction = 'umap', group.by = 'Patients', label=T)


saveRDS(singletData, 'singletData.rds')

########################from raw count ##########################################
setwd("/mnt/8TB/users/shameed/shameed/Zheng/raw_feature_bc_matrix-P1-P14")
list.files(getwd())

files <-list.files(getwd()) 
file_2 <- str_replace(files, 'raw_feature_bc_matrix-', '')
obj_list<- as.list(rep(NA, length(files)))
for (i in 1:length(files)){
  options(Seurat.object.assay.version = "v3")
  C<- ReadMtx(mtx = paste0(files[i], "/matrix.mtx.gz"), cells =paste0(files[i], '/barcodes.tsv.gz'),features =paste0(files[i], "/features.tsv.gz"))
  C_obj<- CreateSeuratObject(C, min.cells = 3, min.features = 200)
  C_obj@meta.data$sample <- rep(file_2[i], nrow(C_obj@meta.data))
  obj_list[[i]] <- C_obj
}

options(Seurat.object.assay.version = "v3")
ZhengData<- merge(x = obj_list[[1]], y = obj_list[-1])
#saveRDS(obj_list, 'obj_list.rds')
meta.data <- separate(ZhengData@meta.data, col = 'sample',
                      into = c('Patient', 'sample_site'), remove = F, sep = '-')
txtDat <- read.table('/mnt/8TB/users/shameed/shameed/Zheng/sample_information.txt')
colnames(txtDat) <- c('Patient', 'Tumor_type')
meta.data <- left_join(meta.data, txtDat, by ='Patient')
#info <- read_csv('Zheng_meta.csv')
rownames(meta.data) <- rownames(ZhengData@meta.data)
ZhengData@meta.data <- meta.data
ZhengData[["percent.mt"]] <-PercentageFeatureSet(ZhengData, pattern = "^MT-")
VlnPlot(ZhengData,features = c("nFeature_RNA",
                               "nCount_RNA",
                               "percent.mt"),ncol = 3)
ZhengData<-subset(ZhengData,subset = percent.mt <10 & nFeature_RNA > 200)
ZhengData <- NormalizeData(ZhengData)
ZhengData <- FindVariableFeatures(ZhengData)
ZhengData <- ScaleData(ZhengData)
ZhengData <- RunPCA(ZhengData)
ElbowPlot(ZhengData, ndims = 50)
options(future.globals.maxSize = 1024^3)
ZhengData<-FindNeighbors(ZhengData,dims = 1:40)
ZhengData<-FindClusters(ZhengData,resolution= 0.3)
ZhengData<-RunUMAP(ZhengData, dims = 1:40)
DimPlot(ZhengData, reduction = 'umap', group.by = 'Tumor_type', label=T)
#VlnPlot(ZhengData, features = c('CD3D', 'PECAM1', 'EPCAM'),group.by ='maintypes_2', raster = F)
#saveRDS(raw_object, 'ZhengData.rds')
saveRDS(ZhengData, 'ZhengData.rds')

#####harmony batch correction##################
ZhengData <- RunHarmony(ZhengData, 'Tumor_type')
ZhengData@reductions
DimPlot(ZhengData, reduction = 'harmony', group.by = 'Tumor_type', label=T)
ZhengData<-FindNeighbors(ZhengData,dims = 1:40,reduction = 'harmony')
ZhengData<-FindClusters(ZhengData,resolution= 0.3, reduction = 'harmony')
ZhengData<-RunUMAP(ZhengData, dims = 1:40, reduction = 'harmony')
DimPlot(ZhengData, reduction = 'umap', group.by = 'Tumor_type', label=T)
DimPlot(ZhengData, reduction = 'umap', group.by = 'refined_celltype', label=T)
saveRDS(ZhengData, 'ZhengData.rds')
#########################################################################

singletData <- readRDS("/mnt/8TB/users/shameed/shameed/Zheng/singletData.rds")
ZhengData <- readRDS("/mnt/8TB/users/shameed/shameed/Zheng/ZhengData.rds")
ZhengData_annotated <- readRDS("/mnt/8TB/users/shameed/shameed/Zheng/ZhengData_annotated.rds")
library(readxl)
Zheng_meta <- read_excel("/mnt/8TB/users/shameed/shameed/Zheng/Zheng_meta.xlsx")
singletData@meta.data %>% head()
ZhengData@meta.data %>% head()
ZhengData_annotated@meta.data %>% head()
ZhengData$Tumor_type %>% unique()
table(singletData$Patients)
colnames(Zheng_meta)[1] <- 'Tumor_type'
meta_zheng <- ZhengData@meta.data %>% left_join(Zheng_meta, 'Tumor_type')
head(meta_zheng)
rownames(meta_zheng) <- rownames(ZhengData@meta.data)
ZhengData@meta.data <- meta_zheng
saveRDS(ZhengData, 'ZhengData.rds')

colnames(Zheng_meta)[1] <- 'Patients'
meta_sing <- singletData@meta.data %>% left_join(Zheng_meta, 'Patients')
head(meta_sing)
rownames(meta_sing) <- rownames(singletData@meta.data)
singletData@meta.data <- meta_sing
saveRDS(singletData, 'singletData.rds')
