#############spatial analysis################################
spatialObj <- readRDS("/mnt/8TB/users/shameed/shameed/Doublet predictions/spatial/ovarian cancer/spatialObj.rds")
spatialObj$response %>% table()
table(spatialObj$origin, spatialObj$response)
spot_labels <- readRDS("/mnt/8TB/users/shameed/shameed/Doublet predictions/spatial/ovarian cancer/spot_labels_2.rds")
head(spot_labels)
library(spacexr)
spot_labels <- normalize_weights(spot_labels)
head(spot_labels)
spot_labels <- as.data.frame(spot_labels)
T_Mac_spots <- spot_labels[, c('T cells', 'Macrophages')]
T_Mac_spots <- T_Mac_spots %>% dplyr:: filter( `T cells` > 0.1 & Macrophages >0.1)

########
spatialObj_list <- readRDS("/mnt/8TB/users/shameed/shameed/Doublet predictions/spatial/ovarian cancer/spatialObj_list.rds")
weight_list_2 <- readRDS("/mnt/8TB/users/shameed/shameed/Doublet predictions/spatial/ovarian cancer/weight_list_2.rds")
paired_lists <- Map(list, spatialObj_list, weight_list_2)

#lapply(spatialObj_list, function(x){unique(x$origin)})
T_Mac_list <- lapply(paired_lists, function(x){
  temp_spa <- x[[1]]
  tempw <- x[[2]]
  tempw <- as.data.frame(tempw)
  tempw <- normalize_weights(tempw)
  tempw <- as.data.frame(tempw)
  T_Mac_spots <- tempw[, c('T cells', 'Macrophages')]
  T_Mac_spots <- T_Mac_spots %>% dplyr:: filter( `T cells` > 0.1 & Macrophages >0.1)
  #head(T_Mac_spots)
  if (nrow(T_Mac_spots !=0)){
    T_Mac_spots$sum <- T_Mac_spots$`T cells` + T_Mac_spots$Macrophages
    T_Mac_spots$origin <- unique(temp_spa$origin)
    return(T_Mac_spots)
    #spot_TMac <- temp_spa[, rownames(T_Mac_spots)]
    #emp_exp <- t(spot_TMac@assays$SCT@data) 
    #T_Mac_spots <- cbind(T_Mac_spots, temp_exp)
  }
})

T_Mac_spots <- do.call(rbind, T_Mac_list)
saveRDS(T_Mac_spots, "/mnt/8TB/users/shameed/shameed/Doublet predictions/spatial/ovarian cancer/T_Mac_spots.rds")
saveRDS(T_Mac_list, "/mnt/8TB/users/shameed/shameed/Doublet predictions/spatial/ovarian cancer/T_Mac_list.rds")

T_Mac_spots$origin <- str_split(T_Mac_spots$origin, '_', simplify = T)[,2] 
spa_response <- data.frame(origin = c(paste0('SP', 1:8)), 
                           response = c('poor', 'good', 'good', 'partial', 'good', 'partial', 'poor', 'poor'))
T_Mac_spots <- T_Mac_spots %>% left_join(spa_response, 'origin')

spot_det <- table(T_Mac_spots$response)
spot_all <- table(spatialObj$response)
spot_norm <- (spot_det / spot_all) * 100
spot_norm <- as.data.frame(spot_norm)
colnames(spot_norm) <- c('response', 'proportion')

p1<- ggplot(spot_norm, aes(x=response, y=proportion, fill= response)) + 
  geom_bar(stat = 'identity', width = 0.2) + 
  labs(title = 'Proportion of T-Mac spots',
       y= 'Proportion (%)', x= 'Response')+
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 18, hjust = 0.5, face = 'bold'),
        legend.text = element_text(size = 17, face = 'bold'),
        legend.key.size = unit(1.2, "cm"),
        axis.title.y = element_text(size = 20, hjust = 0.5, face = 'bold'),
        axis.title.x = element_text(size = 20, hjust = 0.5, face = 'bold'),
        plot.title = element_text(size = 25, hjust = 0.5, face = 'bold')) +
  scale_fill_discrete(name=NULL)
#saveRDS(tubo_bar, 'tubo_bar.rds')

png("/mnt/8TB/users/shameed/shameed/Zheng/figures/TMac_spot.png", width = 10, height = 8.5, units = 'in', res = 600)
p1
dev.off()  

##########spatial gene expression#################
paired_lists <- Map(list, spatialObj_list, weight_list_2)

#lapply(spatialObj_list, function(x){unique(x$origin)})
my_genes <- c('CD4', 'CD8A', 'CD8B')
gene_list <- lapply(paired_lists, function(x){
  temp_spa <- x[[1]]
  
  tempw <- x[[2]]
  tempw <- as.data.frame(tempw)
  tempw <- normalize_weights(tempw)
  tempw <- as.data.frame(tempw)
  T_Mac_spots <- tempw[, c('T cells', 'Macrophages')]
  T_Mac_spots <- T_Mac_spots %>% dplyr:: filter( `T cells` > 0.1 & Macrophages >0.1)
  #head(T_Mac_spots)
  if (nrow(T_Mac_spots !=0)){
    T_Mac_spots$sum <- T_Mac_spots$`T cells` + T_Mac_spots$Macrophages
    T_Mac_spots$origin <- unique(temp_spa$origin)
    #return(T_Mac_spots)
    spot_TMac <- temp_spa[, rownames(T_Mac_spots)]
    temp_exp <- t(spot_TMac@assays$SCT@data[my_genes,]) 
    T_Mac_spots <- cbind(T_Mac_spots, temp_exp)
  }
})

spot_exp <- do.call(rbind, gene_list)
spot_exp$origin <- str_split(spot_exp$origin, '_', simplify = T)[,2] 
spa_response <- data.frame(origin = c(paste0('SP', 1:8)), 
                           response = c('poor', 'good', 'good', 'partial', 'good', 'partial', 'poor', 'poor'))
spot_exp <- spot_exp %>% left_join(spa_response, 'origin')
spot_exp$CD8 <- spot_exp$CD8A + spot_exp$CD8B
spot_exp <- spot_exp %>% select(c(CD8, CD4, response))
spot_exp <- spot_exp %>% pivot_longer(cols = c(CD4, CD8), 
                                      names_to = 'marker', values_to = 'expression')

my_comp <- list(c("CD4", "CD8"),
                c("good", "partial"),
                c("partial", "poor"),
                c("good", "poor"))
p1<-ggplot(spot_exp, aes(x = response, y = expression, fill = marker)) +
  geom_boxplot() +
  theme_minimal() +geom_signif(comparisons = my_comp, test = "wilcox.test", textsize = 4, 
                               step_increase = 0.1, map_signif_level = T) +
  theme_bw() + theme(axis.text.x = element_text(size = 18, face = 'bold'), 
                     axis.text.y = element_text(size = 18, face = 'bold'),
                     axis.title.y = element_text(size = 20, hjust = 0.5, face = 'bold'),
                     axis.title.x = element_text(size = 20, hjust = 0.5, face = 'bold')) +
  ggtitle('T-Mac spot marker expression') + theme(title = element_text(face = 'bold', hjust =0.5 ),
                                                  plot.title = element_text(face = 'bold', hjust = 0.5, size = 25),
                                                  legend.text = element_text(size = 17, face= 'bold'),
                                                  legend.title = element_text(size = 20, hjust = 0.5, face = 'bold'))

png("/mnt/8TB/users/shameed/shameed/Zheng/figures/TMac_spot_marker.png", width = 10, height = 8.5, units = 'in', res = 600)
p1
dev.off()  