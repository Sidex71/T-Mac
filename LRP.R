#############################LR analysis##############################################################
library(readxl)
library(tidyverse)
library(Seurat)
LR_Rowman <- read_excel("/mnt/8TB/users/shameed/shameed/LR_Rowman.xlsx", 
                        sheet = "All.Pairs")
table(TMacDoub$res)
TMacDoub$response  <- TMacDoub$`Response to chemotherapy`
Resist<- subset(TMacDoub, subset =response== "Resistant")
Sensit<- subset(TMacDoub, subset =response== "Sensitive")

###############################
predf <- LR_Rowman %>% filter(Ligand.ApprovedSymbol %in% rownames(Sensit)&
                                Receptor.ApprovedSymbol %in% rownames(Sensit))
LRP <- unique(predf$Pair.Name)
length(LRP)
LRdf <- sapply(LRP, function(x){
  LRSplit <- str_split(x,'_')
  mylig <- LRSplit[[1]][1]
  myrec <- LRSplit[[1]][2]
  ligcell <- Sensit@assays$RNA@data[mylig,]
  ligmean <- mean(ligcell)
  reccell <- Sensit@assays$RNA@data[myrec,]
  recmean <- mean(reccell)
  ligexp <-as.numeric(ligcell > ligmean)
  recexp <- as.numeric(reccell > recmean)
  LRexp <- ligexp + recexp
})

LRdf <- LRdf ==2
LRP_enriched <- as.data.frame(colSums(LRdf))
colnames(LRP_enriched) <- 'Freq'
LRP_enriched$percent <- LRP_enriched$Freq * 100/length(Sensit$orig.ident)
LRP_filt <- LRP_enriched[LRP_enriched$Freq >= 5,]  #280 LRPs

saveRDS(LRP_enriched, '/mnt/8TB/users/shameed/shameed/Zheng/LRP_enriched_Sensit.rds')

LRP_filt <-LRP_filt[order(LRP_filt$percent, decreasing = T),]
LRP_filt <- LRP_filt %>% rownames_to_column('LRP')

saveRDS(LRP_filt, '/mnt/8TB/users/shameed/shameed/Zheng/LRP_filt_Sensit.rds')

intersect(LRP_filt_Resist$LRP, LRP_filt_Sensit$LRP)
p1<-ggplot(LRP_filt[1:50,], aes(x = reorder(LRP, percent), y = percent, fill = percent)) + 
  geom_bar(stat = 'identity') + 
  labs(title = ' ',
       y= 'cell proportion (%)', x=NULL)+
  theme_bw() + coord_flip() #+ theme(legend.position = "none")
png("/mnt/8TB/users/shameed/shameed/Zheng/figures/LRP_Sensitive.png", width = 15, height = 10.5, units = 'in', res = 600)
p1
dev.off()

########################Resistant
predf <- LR_Rowman %>% filter(Ligand.ApprovedSymbol %in% rownames(Resist)&
                                Receptor.ApprovedSymbol %in% rownames(Resist))
LRP <- unique(predf$Pair.Name)
length(LRP)
LRdf <- sapply(LRP, function(x){
  LRSplit <- str_split(x,'_')
  mylig <- LRSplit[[1]][1]
  myrec <- LRSplit[[1]][2]
  ligcell <- Resist@assays$RNA@data[mylig,]
  ligmean <- mean(ligcell)
  reccell <- Resist@assays$RNA@data[myrec,]
  recmean <- mean(reccell)
  ligexp <-as.numeric(ligcell > ligmean)
  recexp <- as.numeric(reccell > recmean)
  LRexp <- ligexp + recexp
})

LRdf <- LRdf ==2
LRP_enriched <- as.data.frame(colSums(LRdf))
colnames(LRP_enriched) <- 'Freq'
LRP_enriched$percent <- LRP_enriched$Freq * 100/length(Resist$orig.ident)
LRP_filt <- LRP_enriched[LRP_enriched$Freq >= 5,]  #398 LRPs

saveRDS(LRP_enriched, '/mnt/8TB/users/shameed/shameed/Zheng/LRP_enriched_Resist.rds')

LRP_filt <-LRP_filt[order(LRP_filt$percent, decreasing = T),]
LRP_filt <- LRP_filt %>% rownames_to_column('LRP')

saveRDS(LRP_filt, '/mnt/8TB/users/shameed/shameed/Zheng/LRP_filt_Resist.rds')

p2<-ggplot(LRP_filt[1:50,], aes(x = reorder(LRP, percent), y = percent, fill = percent)) + 
  geom_bar(stat = 'identity') + 
  labs(title = ' ',
       y= 'cell proportion (%)', x=NULL)+
  theme_bw() + coord_flip() #+ theme(legend.position = "none")
png("/mnt/8TB/users/shameed/shameed/Zheng/figures/LRP_Resist.png", width = 15, height = 10.5, units = 'in', res = 600)
p2
dev.off()

##########################spatial LRP##########################################
library(spacexr)
LRP_filt_Resist <- readRDS("/mnt/8TB/users/shameed/shameed/Zheng/LRP_filt_Resist.rds")
LRP_filt_Sensit <- readRDS("/mnt/8TB/users/shameed/shameed/Zheng/LRP_filt_Sensit.rds")

spatialObj_list <- readRDS("/mnt/8TB/users/shameed/shameed/Doublet predictions/spatial/ovarian cancer/spatialObj_list.rds")
weight_list_2 <- readRDS("/mnt/8TB/users/shameed/shameed/Doublet predictions/spatial/ovarian cancer/weight_list_2.rds")
paired_lists <- Map(list, spatialObj_list, weight_list_2)
names(paired_lists) <- paste0('SP', 1:8)

LRE <- unique(c(LRP_filt_Resist$LRP, LRP_filt_Sensit$LRP))

Spa_LRdf_list <- lapply(paired_lists, function(x){
  temp_spa <- x[[1]]
  temp_spa <- UpdateSeuratObject(temp_spa)
  tempw <- x[[2]]
  tempw <- as.data.frame(tempw)
  tempw <- normalize_weights(tempw)
  T_Mac_spots <- tempw[, c('T cells', 'Macrophages')]
  T_Mac_spots <- T_Mac_spots %>% dplyr:: filter( `T cells` > 0.1 & Macrophages >0.1)
  
  if (nrow(T_Mac_spots !=0)) {
    spot_BC <- temp_spa[, rownames(T_Mac_spots)]
    predf <- LR_Rowman %>% filter(Ligand.ApprovedSymbol %in% rownames(temp_spa)&
                                    Receptor.ApprovedSymbol %in% rownames(temp_spa))
    my_LRE <- intersect(rownames(LRE), predf$Pair.Name)
    DefaultAssay(spot_BC) <- 'SCT'
    predf <- LR_Rowman %>% filter(Ligand.ApprovedSymbol %in% rownames(temp_spa)&
                                    Receptor.ApprovedSymbol %in% rownames(temp_spa))
    
    my_LRE <- intersect(LRE, predf$Pair.Name)
    my_enrich <-sapply(my_LRE, function(x){
      LRSplit <- str_split(x,'_')
      mylig <- LRSplit[[1]][1]
      myrec <- LRSplit[[1]][2]
      ligcell <- spot_BC@assays$SCT@data[mylig,]
      ligmean <- mean(ligcell)
      reccell <- spot_BC@assays$SCT@data[myrec,]
      recmean <- mean(reccell)
      ligexp <-as.numeric(ligcell > ligmean)
      recexp <- as.numeric(reccell > recmean)
      LRexp <- ligexp + recexp
    })
    rownames(my_enrich) <- colnames(spot_BC)
    #my_enrich$ <- unique(spot_BC$origin)
    return(my_enrich)  
  }
  
})


library(purrr)
Spa_LRdf_list <- compact(Spa_LRdf_list)

saveRDS(Spa_LRdf_list, '/mnt/8TB/users/shameed/shameed/Doublet predictions/spatial/ovarian cancer/Spa_LRdf_list.rds')


Spa_df_list <- lapply(Spa_LRdf_list, function(x){
  x <- x ==2
  LRP_enriched <- as.data.frame(colSums(x))
  colnames(LRP_enriched) <- 'Freq'
  LRP_enriched$percent <- LRP_enriched$Freq * 100/length(rownames(x))
  LRP_enriched <- LRP_enriched %>% rownames_to_column('LRP')
  return(LRP_enriched)
}) 


spa_response <- data.frame(origin = c(paste0('SP', 1:8)), 
                           response = c('poor', 'good', 'good', 'partial', 'good', 'partial', 'poor', 'poor'))

poor_list <- c(Spa_df_list['SP1'], Spa_df_list['SP8'])
good_list <- c(Spa_df_list['SP2'], Spa_df_list['SP3'], Spa_df_list['SP5'])
partial_list <- c(Spa_df_list['SP6'])
bad_list <- c(poor_list, partial_list)
good_LRdf <- do.call(rbind, good_list)
bad_LRdf <- do.call(rbind, bad_list)



good_spots <- sum(sapply(Spa_LRdf_list[c('SP2', 'SP3', 'SP5')], 
                         function(x) length(rownames(x)))) #495
good_LRdf <- good_LRdf %>% group_by(LRP) %>% 
  mutate(total_freq = sum(Freq),
         final_percent = sum(Freq)*100/good_spots)

good_LRdf_final <- good_LRdf %>% select(c(LRP, total_freq, final_percent))
good_LRdf_final <- unique(good_LRdf_final)

good_LRdf_final <-good_LRdf_final[order(good_LRdf_final$final_percent, decreasing = T),]
good_LRdf_final <- good_LRdf_final[good_LRdf_final$LRP %in% LRP_filt_Sensit$LRP,]

saveRDS(good_LRdf_final, '/mnt/8TB/users/shameed/shameed/Doublet predictions/spatial/ovarian cancer/good_LRdf_final.rds')

good_df <- left_join(good_LRdf_final, LRP_filt_Sensit, 'LRP')
colnames(good_df) <- c('LRP', 'Spot_Freq', 'Spot_%', 'Doub_Freq', 'Doub_%')

ct <- cor.test(good_df$`Spot_%`, good_df$`Doub_%`, method = "pearson")

r_val <- round(ct$estimate, 3)
p_val <- signif(ct$p.value, 3)

p1<- ggplot(good_df, aes(x = `Spot_%`, y = `Doub_%`)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("r = ", r_val, "\n p = ", p_val),
           hjust = 1.1, vjust = 1.5, size = 5) +
  labs(x = 'Spots (%)',y = 'Doublets (%)',
       title = 'Good Responders: Ligand–Receptor Correlations') +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 22),  
        axis.title = element_text(face = "bold", size = 19)                 
  )

png("/mnt/8TB/users/shameed/shameed/Zheng/figures/corr_good.png", width = 8, height = 6.5, units = 'in', res = 600)
p1
dev.off()

###########################################bad
bad_spots <- sum(sapply(Spa_LRdf_list[c('SP1', 'SP8', 'SP6')], 
                        function(x) length(rownames(x)))) #79
bad_LRdf <- bad_LRdf %>% group_by(LRP) %>% 
  mutate(total_freq = sum(Freq),
         final_percent = sum(Freq)*100/bad_spots)

bad_LRdf_final <- bad_LRdf %>% select(c(LRP, total_freq, final_percent))
bad_LRdf_final <- unique(bad_LRdf_final)

bad_LRdf_final <-bad_LRdf_final[order(bad_LRdf_final$final_percent, decreasing = T),]
bad_LRdf_final <- bad_LRdf_final[bad_LRdf_final$LRP %in% LRP_filt_Resist$LRP,]

saveRDS(bad_LRdf_final, '/mnt/8TB/users/shameed/shameed/Doublet predictions/spatial/ovarian cancer/bad_LRdf_final.rds')

bad_df <- left_join(bad_LRdf_final, LRP_filt_Resist, 'LRP')
colnames(bad_df) <- c('LRP', 'Spot_Freq', 'Spot_%', 'Doub_Freq', 'Doub_%')

ct <- cor.test(bad_df$`Spot_%`, bad_df$`Doub_%`, method = "pearson")

r_val <- round(ct$estimate, 3)
p_val <- signif(ct$p.value, 3)


p2 <- ggplot(bad_df, aes(x = `Spot_%`, y = `Doub_%`)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("r = ", r_val, "\n p = ", p_val),
           hjust = 1.1, vjust = 1.5, size = 5) +
  labs(x = 'Spots (%)',y = 'Doublets (%)',
       title = 'Non Responders: Ligand–Receptor Correlations') +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 22),  
        axis.title = element_text(face = "bold", size = 19)                 
  )

png("/mnt/8TB/users/shameed/shameed/Zheng/figures/corr_bad.png", width = 8, height = 6.5, units = 'in', res = 600)
p2
dev.off()
