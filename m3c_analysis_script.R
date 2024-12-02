library(Cardinal)
# library(rayshader)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(M3C)
library(spatstat.geom)
library(RColorBrewer)
library(keras)
library(reshape2)

res_ppm <- 20
mass_range = c(50, 650)

setCardinalBPPARAM(BPPARAM=SerialParam())

normal_peak_2 <- 157.0771
normal_peak <- 208.1144
#### data loading ####
collect_cell_info <- function(min_size=1,max_size=20,imzmlPath,path_spotlist,std_dev=3,bg_ratio=100,include_bg=F,cell_peak = 281.2485){
  
  msData_cardinal <- Cardinal::readMSIData(imzmlPath, attach.only=TRUE, resolution = res_ppm, units = "ppm", mass.range = mass_range)
  
  prefix <- c()
  
  if(grepl("EOC",imzmlPath)){
    prefix <- "EOC"
  }
  if(grepl("SIMA9",imzmlPath)){
    prefix <- "SIMA9"
  }
  if(grepl("IPSC",imzmlPath)){
    prefix <- "hiPSC2"
  }
  if(grepl("35er",imzmlPath)){
    prefix <- "hiPSC1"
  }
  if(grepl("68er",imzmlPath)){
    prefix <- "hiPSC2"
  }
  if(grepl("DMA",imzmlPath)){
    prefix <- "SIMA9"
  }
  
  tmp_tab <- read.table(path_spotlist)
  
  simple_labels_spotlist <- as.character(tmp_tab$V4)
  for(simp in seq(simple_labels_spotlist)){
    simple_labels_spotlist[simp] <- paste0("0",simple_labels_spotlist[simp])
  }
  
  simple_labels <- as.factor(simple_labels_spotlist)
  
  
  simple_labels_2 <- as.character(simple_labels)
  
  for(simp in seq(simple_labels_2)){
    
    str_tmp <- paste0(strsplit(simple_labels_2[simp],split="0")[[1]][-1],collapse="0")
    
    simple_labels_2[simp] <- paste0(strsplit(str_tmp,split="_")[[1]][1:3],collapse="_")
    
  }
  
  simple_labels <- as.factor(simple_labels_2)
  
  ms_mat <- as.matrix(t(Cardinal::spectra(msData_cardinal)))
  
  norm_peak_idx <- which.min(abs(normal_peak-msData_cardinal@featureData@mz))#
  norm_peak_idx_2 <- which.min(abs(normal_peak_2-msData_cardinal@featureData@mz))#
  
  ms_mat_normed <- (ms_mat)
  
  gc()
  ms_mat_normed <- ms_mat/ms_mat[,norm_peak_idx]
  
  
  ms_mat_normed[is.nan(ms_mat_normed)] <- 0
  ms_mat_normed[is.infinite(ms_mat_normed)] <- 0
  
  ms_mat_normed <- ms_mat_normed[,-c(norm_peak_idx,norm_peak_idx_2)]
  gc()
  
  ms_mat_normed <- scale(ms_mat_normed)
  
  ms_mat_normed[is.nan(ms_mat_normed)] <- 0
  ms_mat_normed[is.infinite(ms_mat_normed)] <- 0
  
  gc()
  
  normalize_peak <- normal_peak#208.1140
  tol <- 0.01
  
  cell_peak_idx <- which.min(abs(cell_peak-msData_cardinal@featureData@mz[-c(norm_peak_idx,norm_peak_idx_2)]))
  
  thresh <- mean(ms_mat_normed[,cell_peak_idx])+std_dev*sd(ms_mat_normed[,cell_peak_idx])
  
  msData_norm_2 <- subsetFeatures(msData_cardinal, mz > cell_peak-tol, mz < cell_peak+tol)
  norm_img_2 <- image(msData_norm_2)
  norm_mat_2 <- norm_img_2[["facets"]][[1]][[1]][["values"]]
  
  cell_idx <- which(ms_mat_normed[,cell_peak_idx]>thresh)
  
  if(include_bg==T){
    bg_index <- sample(which(ms_mat_normed[,cell_peak_idx]<thresh),size=round(length(cell_idx)/bg_ratio))
    cell_idx <- sort(c(cell_idx,bg_index))
  }
  
  print(length(cell_idx))
  msData_norm <- subsetFeatures(msData_cardinal[cell_idx], mz > cell_peak-tol, mz < cell_peak+tol)
  norm_img <- image(msData_norm)
  norm_mat <- norm_img[["facets"]][[1]][[1]][["values"]]
  
  x_coord <- coord(msData_cardinal)$x-min(coord(msData_cardinal)$x)+1
  y_coord <- coord(msData_cardinal)$y-min(coord(msData_cardinal)$y)+1
  
  mapper_matrix <-matrix(NA,nrow=nrow(norm_mat_2),ncol=ncol(norm_mat_2))
  
  t <- 1
  
  for(ind in c(1:length(x_coord))){
    
    mapper_matrix[x_coord[ind],y_coord[ind]] <- t
    t<-t+1
    
  }
  
  mapper_matrix[which(!mapper_matrix %in% cell_idx)] <- NA
  
  mapper_matrix_bin <- mapper_matrix
  mapper_matrix_bin[!is.na(mapper_matrix)] <- 1
  
  cell_im <- im((mapper_matrix_bin))
  
  cell_con <- connected(cell_im)
  
  cell_candidates <- levels(cell_con[["v"]])
  
  cell_list_idxs <- list()
  
  for(cand in cell_candidates){
    
    idxs <- which(cell_con[["v"]]==cand)
    
    if(length(idxs)<min_size){
      next
    }
    if(length(idxs)>max_size){
      next
    }
    
    cell_list_idxs[[length(cell_list_idxs)+1]] <- idxs
    
  }
  
  cell_means <- c()
  poi_means <- c()
  cell_lab <- c()
  
  cells_all_mz_mat_sum <- matrix(0,nrow=1,ncol=dim(ms_mat_normed)[2])
  cells_all_mz_mat_mean <- matrix(0,nrow=1,ncol=dim(ms_mat_normed)[2])
  cells_all_mz_mat_sky <- matrix(0,nrow=1,ncol=dim(ms_mat_normed)[2])
  
  for(cell_idxs in cell_list_idxs){
    
    
    idx_lab <- mapper_matrix[cell_idxs]
    
    idx_lab <- idx_lab[!is.na(idx_lab)]
    
    cell_lab <- c(cell_lab,as.character(simple_labels[idx_lab[1]]))
    
    if(length(idx_lab)==1){
      
      
      
      cells_all_mz_mat_sum <- rbind(cells_all_mz_mat_sum,(ms_mat_normed[idx_lab,]))
      cells_all_mz_mat_mean <- rbind(cells_all_mz_mat_mean,(ms_mat_normed[idx_lab,]))
      cells_all_mz_mat_sky <- rbind(cells_all_mz_mat_sky,(ms_mat_normed[idx_lab,]))
      
      next
    }
    
    
    cell_mat_tmp_sum <- colSums(ms_mat_normed[idx_lab,])#
    cell_mat_tmp_mean <- colMeans(ms_mat_normed[idx_lab,])#
    cell_mat_tmp_sky <- apply(ms_mat_normed[idx_lab,],2,max)#
    
    cells_all_mz_mat_sum <- rbind(cells_all_mz_mat_sum,cell_mat_tmp_sum)
    cells_all_mz_mat_mean <- rbind(cells_all_mz_mat_mean,cell_mat_tmp_mean)
    cells_all_mz_mat_sky <- rbind(cells_all_mz_mat_sky,cell_mat_tmp_sky)
    
  }
  
  cells_all_mz_mat_sum <- cells_all_mz_mat_sum[-1,]
  cells_all_mz_mat_mean <- cells_all_mz_mat_mean[-1,]
  cells_all_mz_mat_sky <- cells_all_mz_mat_sky[-1,]
  
  cells_all_mz_mat <- list(cells_all_mz_mat_sum,cells_all_mz_mat_mean,cells_all_mz_mat_sky)
  
  print(dim(cells_all_mz_mat_mean))
  
  return(list(cells_all_mz_mat,cell_lab,mapper_matrix,msData_cardinal,cell_list_idxs,prefix))
  
}
min_max <- function(x) {
  return((x - min(x,na.rm=T)) / (max(x,na.rm=T) - min(x,na.rm=T)))
}

table_raw <- read.csv("Z://10-James_Cairns//HMDBE_KEGG_10.csv")

cell_info_list <- list()
cell_info_list[[1]] <- collect_cell_info(imzmlPath="U:/James/HMDBE_KEGG_10/230601_EOC_20_1_FINAL.imzML",path_spotlist = "U:/James/spotlists/230601_EOC_20_1.25_spotlist.txt")
cell_info_list[[2]] <- collect_cell_info(imzmlPath="U:/James/HMDBE_KEGG_10/230615_SIMA9_lyo_20_1_FINAL.imzML",path_spotlist = "U:/James/spotlists/230615_SIMA9_lyo_20_1_spotlist.txt")
cell_info_list[[3]] <- collect_cell_info(imzmlPath="U:/James/HMDBE_KEGG_10/230616_SIMA9_lyo_20_1_FINAL.imzML",path_spotlist = "U:/James/spotlists/230616_SIMA9_lyo_20_1_spotlist.txt")
cell_info_list[[4]] <- collect_cell_info(imzmlPath="U:/James/HMDBE_KEGG_10/230618_SIMA9_lyo_20_1_FINAL.imzML",path_spotlist = "U:/James/spotlists/230618_SIMA9_lyo_20_1_spotlist.txt")
cell_info_list[[5]] <- collect_cell_info(imzmlPath="U:/James/HMDBE_KEGG_10/231027_68er_20_1_FINAL.imzML",path_spotlist = "U:/James/spotlists/231027_68er_20_1_spotlist.txt")
cell_info_list[[6]] <- collect_cell_info(imzmlPath="U:/James/HMDBE_KEGG_10/231027_35er_20_1_FINAL.imzML",path_spotlist = "U:/James/spotlists/231027_35er_20_1_spotlist.txt")


o <- 2

indices <- which(sapply(cell_info_list, function(x) length(x) > 0))

cell_lab <- unlist(sapply(cell_info_list[indices], function(x) x[[2]]))
prefixes <- unlist(sapply(cell_info_list[indices], function(x) x[[6]]))
sublist_lengths <- sapply(cell_info_list[indices], function(x) length(x[[2]]))
summarize_to_mat <- function(x) x[[1]][[o]]
cells_all_mz_mat <- do.call(rbind, lapply(cell_info_list[indices], summarize_to_mat))

repeated_prefixes <- c()

for(i in 1:length(indices)) {
  repeated_prefixes <- c(repeated_prefixes, rep(prefixes[i], times=sublist_lengths[i]))
}

combined_labels <- paste(repeated_prefixes, cell_lab, sep="_")


#### m3C ####

m3c_clustering <- function(mz_input=129.02){
  
  best_clus_list <- list()
  m3c_list <- list()
  plot_list <- list()
  img_lists <- list()
  
  for(indic_m3 in indices){
    gc()
    print(indic_m3)
    cell_info_x <- cell_info_list[[indic_m3]]
    cell_list_idxs <- cell_info_list[[indic_m3]][[5]]
    msData_cardinal <- cbind(cell_info_list[[indic_m3]][[4]])
    
    norm_peak_idx <- which.min(abs(normal_peak-msData_cardinal@featureData@mz))#
    norm_peak_idx_2 <- which.min(abs(normal_peak_2-msData_cardinal@featureData@mz))
    
    mz <- msData_cardinal@featureData@mz[-c(norm_peak_idx,norm_peak_idx_2)]
    # break
    mapper_matrix <- cell_info_list[[indic_m3]][[3]]
    
    
    tmp_df <- as.data.frame(t((cell_info_list[[indic_m3]][[1]][[o]])))
    
    
    idx_rm <- -c(which.min(abs(mz-mz_input)))
    
    if(length(idx_rm)!=0){
      tmp_df <- tmp_df[-idx_rm,]
      
      dimnames(tmp_df)[[1]] <- msData_cardinal@featureData@mz[-c(norm_peak_idx,norm_peak_idx_2)][-idx_rm]
    }else{
      
      
      dimnames(tmp_df)[[1]] <- msData_cardinal@featureData@mz[-c(norm_peak_idx,norm_peak_idx_2)]
      
    }
    
    dimnames(tmp_df)[[2]] <- cell_info_list[[indic_m3]][[2]] 
    
    tmp_df <- data.frame(tmp_df)
    
    med <- cell_info_list[[indic_m3]][[2]]
    
    for (idx in c(1:length(med))){
      med[idx] <- paste(strsplit(med[idx],split="_")[[1]][-3],collapse="_")
      if(is.na(med[idx])){
        print(idx)
      }
      
    }
    
    
    med <- as.factor(med)
    desx_tmp <- data.frame(ID=dimnames(tmp_df)[[2]],Cell_Label=med)
    filtered_data <- featurefilter((rbind(tmp_df,tmp_df)))#tmp_df, 100,method="var")
    
    
    mz_tmp <- as.numeric(row.names(filtered_data$filtered_data))
    
    remove(tmp_df)
    gc()
    set.seed(42)
    
    m3c_clusters_hc <- M3C::M3C(filtered_data$filtered_data, cores=1, iters = 5, maxK=8,seed = 42, des = desx_tmp,repsref = 200,repsreal = 200,clusteralg = "hc")
    
    #first we check where a pvalue is under 0.05 rejecting the null hypothesis, meaning that the data is not gaussian distributed --> k=1
    idx_relevant <- which(m3c_clusters_hc[["plots"]][[3]][["data"]][["NORM_P"]]<0.05)
    
    if(length(idx_relevant)==0){
      clus_num <- 1
      best_clus_list[[indic_m3]] <- clus_num
    }else{
      #then we find the cluster with the highest cluster stability
      clus_num <- idx_relevant[which.max(m3c_clusters_hc[["plots"]][[4]][["data"]][["RCSI"]][idx_relevant])]+1
      
      best_clus_list[[indic_m3]] <- clus_num
    }
    print(clus_num)
    
    m3c_list[[indic_m3]] <- m3c_clusters_hc
    
    if(clus_num==1){
      next
    }
    
    cell_peak <- 281.2485#281.2485#303.2329#281.2482#303.2329 #281.2482#281.2482#283.2642#281.2482
    goal_mz <-  129.0193#124.0074#231.0986#124.0074
    taurine_mz <- 124.007#208.1140
    adenine_mz <- 306.08
    tol <- 0.01
    # 
    cell_peak_idx <- which.min(abs(cell_peak-msData_cardinal@featureData@mz[-c(norm_peak_idx,norm_peak_idx_2)]))
    goal_peak_idx <- which.min(abs(goal_mz-msData_cardinal@featureData@mz[-c(norm_peak_idx,norm_peak_idx_2)]))
    taur_peak_idx <- which.min(abs(taurine_mz-msData_cardinal@featureData@mz[-c(norm_peak_idx,norm_peak_idx_2)]))
    ad_peak_idx <- which.min(abs(adenine_mz-msData_cardinal@featureData@mz[-c(norm_peak_idx,norm_peak_idx_2)]))
    
    
    df <- data.frame(
      Clusters = as.factor(m3c_clusters_hc$realdataresults[[clus_num]]$assignments),
      Treatment = as.factor(desx_tmp$Cell_Label),
      `Itaconic Acid` = (cell_info_list[[indic_m3]][[1]][[o]])[,goal_peak_idx],
      `Fatty Acid 18:1` = (cell_info_list[[indic_m3]][[1]][[o]])[,cell_peak_idx],
      `Taurine` = (cell_info_list[[indic_m3]][[1]][[o]])[,taur_peak_idx],
      `Glutathione` = (cell_info_list[[indic_m3]][[1]][[o]])[,ad_peak_idx]
    )
    
    if(all(levels(desx_tmp$Cell_Label) == c("LPS_1", "LPS_2", "LPS_3", "LPS_4", "VEH_1", "VEH_2", "VEH_3", "VEH_4"))){
      levels(df$Treatment) <- c("VEH","VEH","VEH","VEH","LPS 500","LPS 500","LPS 500","LPS 500")
    }else{
      df$Treatment <- factor(df$Treatment, levels=c("VEH_1","VEH_2","LPS_0.1", "LPS_1", "LPS_2.5", "LPS_10", "LPS_100", "LPS_500"))
      levels(df$Treatment) <- c("VEH","VEH","LPS 0.1","LPS 1","LPS 2.5","LPS 10","LPS 100","LPS 500")
    }
    
    conc_vec <- as.character(df$Treatment)
    
    conc_vec[which(conc_vec=="VEH")] <- 0.001
    conc_vec[which(conc_vec=="LPS 0.1")] <- 0.1
    conc_vec[which(conc_vec=="LPS 1")] <- 1
    conc_vec[which(conc_vec=="LPS 2.5")] <- 2.5
    conc_vec[which(conc_vec=="LPS 10")] <- 10
    conc_vec[which(conc_vec=="LPS 100")] <- 100
    conc_vec[which(conc_vec=="LPS 500")] <- 500
    
    df$Concentrations <- as.numeric(conc_vec)
    
    p1 <- ggplot(data = df, aes(x = log10(Concentrations), y = (`Itaconic.Acid`), group = Concentrations)) +
      theme(text = element_text(size=30), legend.position = "none") +
      geom_violin(width = 1.3) + # Adjusted width here
      geom_boxplot(width = 0.1, outlier.shape = NA) +
      scale_x_continuous(breaks = log10(c(0, 0.1, 1, 2.5, 10, 100, 500)), labels = c("0", "0.1", "1", "2.5", "10", "100", "500")) +
      stat_summary(fun = median, geom = "point", color = "gray50", size = 3) +
      stat_summary(fun = median, geom = "line", aes(group = 1), color = "gray50", size = 1)+
      # scale_y_continuous(limits = c(-2.5, 1.5)) +
      labs(x = "Concentration [ng/mL]")+
      labs(y = "Normalized Intensity")
    
    p2 <- ggplot(data = df, aes(x = log10(Concentrations), y = `Taurine`, group = Concentrations)) +
      theme(text = element_text(size=30), legend.position = "none") +
      geom_violin(width = 1.3) + # Adjusted width here
      geom_boxplot(width = 0.1, outlier.shape = NA) +
      scale_x_continuous(breaks = log10(c(0, 0.1, 1, 2.5, 10, 100, 500)), labels = c("0", "0.1", "1", "2.5", "10", "100", "500")) +
      stat_summary(fun = median, geom = "point", color = "gray50", size = 3) +
      stat_summary(fun = median, geom = "line", aes(group = 1), color = "gray50", size = 1)+
      # scale_y_continuous(limits = c(-2.5, 1.5)) +
      labs(x = "Concentration [ng/mL]")+
      labs(y = "Normalized Intensity")
    
    p3 <- ggplot(data = df, aes(x = log10(Concentrations), y = `Fatty.Acid.18.1`, group = Concentrations)) +
      theme(text = element_text(size=30), legend.position = "none") +
      geom_violin(width = 1.3) + # Adjusted width here
      geom_boxplot(width = 0.1, outlier.shape = NA) +
      scale_x_continuous(breaks = log10(c(0, 0.1, 1, 2.5, 10, 100, 500)), labels = c("0", "0.1", "1", "2.5", "10", "100", "500")) +
      stat_summary(fun = median, geom = "point", color = "gray50", size = 3) +
      stat_summary(fun = median, geom = "line", aes(group = 1), color = "gray50", size = 1)+
      # scale_y_continuous(limits = c(-2.5, 1.5)) +
      labs(x = "Concentration [ng/mL]")+
      labs(y = "Normalized Intensity")
    
    p4 <- ggplot(data = df, aes(x = log10(Concentrations), y = `Glutathione`, group = Concentrations)) +
      theme(text = element_text(size=30), legend.position = "none") +
      geom_violin(width = 1.3) + # Adjusted width here
      geom_boxplot(width = 0.1, outlier.shape = NA) +
      scale_x_continuous(breaks = log10(c(0, 0.1, 1, 2.5, 10, 100, 500)), labels = c("0", "0.1", "1", "2.5", "10", "100", "500")) +
      stat_summary(fun = median, geom = "point", color = "gray50", size = 3) +
      stat_summary(fun = median, geom = "line", aes(group = 1), color = "gray50", size = 1)+
      # scale_y_continuous(limits = c(-2.5, 1.5)) +
      labs(x = "Concentration [ng/mL]")+
      labs(y = "Normalized Intensity")
    
    gridExtra::grid.arrange(p1,p2,p3,p4,ncol=2)
    
    
    plot_list[[indic_m3]] <- list(p1,p2,p3)#,pca,pca_c,tsne,tsne_c)
    
    rm(p1,p2,p3,p4)#,pca,pca_c,tsne,tsne_c)
    # next()
    
    
    # m1 <- matrix(NA,nrow=nrow(mapper_matrix),ncol=ncol(mapper_matrix))
    mclus_2 <- matrix(NA,nrow=nrow(mapper_matrix),ncol=ncol(mapper_matrix))
    mclus_3 <- matrix(NA,nrow=nrow(mapper_matrix),ncol=ncol(mapper_matrix))
    mclus_4 <- matrix(NA,nrow=nrow(mapper_matrix),ncol=ncol(mapper_matrix))
    
    gc()
    # Create assignments as vectors
    assigns_2 <- as.vector(m3c_clusters_hc$realdataresults[[2]]$assignments)
    assigns_3 <- as.vector(m3c_clusters_hc$realdataresults[[3]]$assignments)
    assigns_4 <- as.vector(m3c_clusters_hc$realdataresults[[4]]$assignments)
    
    for(i in seq_along(cell_list_idxs)){
      
      x_range <-  dim(mapper_matrix)[1]
      
      if(i<=length(mapper_matrix)){
        x_inds <- c(1:x_range)
      } else {
        x_inds <- c((x_range+1):dim(mapper_matrix)[1])
      }
      
      idx_lab <- mapper_matrix[x_inds,][cell_list_idxs[[i]]]
      
      if(i%%500==0){
        print(i)
      }
      
      for(idx_tmp in idx_lab){
        idx_matrix <- which((mapper_matrix[x_inds,])==idx_tmp)
        
        mclus_2[x_inds,][idx_matrix] <- assigns_2[i]
        mclus_3[x_inds,][idx_matrix] <- assigns_3[i]
        mclus_4[x_inds,][idx_matrix] <- assigns_4[i]
      }
    }
    
    
    img_lists[[indic_m3]] <- list(mclus_2,mclus_3,mclus_4)
    
    rm(mclus_2,mclus_3,mclus_4,mapper_matrix)
    gc()
    
    
  }
  list(mz=mz,best_clus_list = best_clus_list, m3c_list = m3c_list, plot_list = plot_list, img_lists = img_lists)
  
}

m3c_clustering_ita <- m3c_clustering(mz_input = 129.02)

m3c_clustering_taur <- m3c_clustering(mz_input = 124.01)


mz <- m3c_clustering_ita$mz

norm_peak_idx <- which.min(abs(normal_peak-mz))
norm_peak_idx_2 <- which.min(abs(normal_peak_2-mz))

idx_mz <- c(which.min(abs(mz[-c(norm_peak_idx,norm_peak_idx_2)]-129.0193)))
cell_lab_ita <- c()

repeated_prefixes <- c()

for(i in (indices)) {
  repeated_prefixes <- c(repeated_prefixes, rep(prefixes[i], times=sublist_lengths[i]))
}

combined_labels <- paste(repeated_prefixes, cell_lab, sep="_")

for(indo in indices){
  tmp_vec <- cell_info_list[[indo]][[1]][[o]][, c(idx_mz)]
  
  clus_num <- m3c_clustering_ita$best_clus_list[[indo]]
  
  if(clus_num == 1){
    lab_tmp <- rep("Ita-", length(tmp_vec))
    cell_lab_ita <- c(cell_lab_ita, lab_tmp)
    next
  }
  
  if(clus_num >= 5){
    lab_tmp <- rep("Ita-", length(tmp_vec))
    cell_lab_ita <- c(cell_lab_ita, lab_tmp)
    next
  }
  
  lab_tmp <- as.vector(m3c_clustering_ita$m3c_list[[indo]]$realdataresults[[clus_num]]$assignments)
  
  if(sum(is.na(lab_tmp)) == length(lab_tmp)){
    cell_lab_ita <- c(cell_lab_ita, lab_tmp)
  } else {
    mean_values <- sapply(1:clus_num, function(i) mean(tmp_vec[lab_tmp == i]))
    cluster_max <- which.max(mean_values)
    cluster_min <- which.min(mean_values)
    
    lab_tmp[lab_tmp != cluster_max & lab_tmp != cluster_min] <- "Ita-"
    lab_tmp[lab_tmp == cluster_min] <- "Ita-"
    lab_tmp[lab_tmp == cluster_max] <- "Ita+"
    
    cell_lab_ita <- c(cell_lab_ita, lab_tmp)
  }
}

idx_mz <- c(which.min(abs(mz[-c(norm_peak_idx,norm_peak_idx_2)]-124.007)))
cell_lab_taur <- c()

for(indo in indices){
  tmp_vec <- cell_info_list[[indo]][[1]][[o]][, c(idx_mz)]
  
  clus_num <- m3c_clustering_taur$best_clus_list[[indo]]
  
  if(clus_num == 1){
    lab_tmp <- rep("Taur-", length(tmp_vec))
    cell_lab_taur <- c(cell_lab_taur, lab_tmp)
    next
  }
  
  if(clus_num >= 5){
    lab_tmp <- rep("Taur-", length(tmp_vec))
    cell_lab_taur <- c(cell_lab_taur, lab_tmp)
    next
  }
  
  lab_tmp <- as.vector(m3c_clustering_taur$m3c_list[[indo]]$realdataresults[[clus_num]]$assignments)
  
  if(sum(is.na(lab_tmp)) == length(lab_tmp)){
    cell_lab_taur <- c(cell_lab_taur, lab_tmp)
  } else {
    mean_values <- sapply(1:clus_num, function(i) mean(tmp_vec[lab_tmp == i]))
    cluster_max <- which.max(mean_values)
    cluster_min <- which.min(mean_values)
    
    lab_tmp[lab_tmp != cluster_max & lab_tmp != cluster_min] <- "Taur-"
    lab_tmp[lab_tmp == cluster_min] <- "Taur-"
    lab_tmp[lab_tmp == cluster_max] <- "Taur+"
    
    cell_lab_taur <- c(cell_lab_taur, lab_tmp)
  }
}

cell_lab_comb <- paste(cell_lab_ita, cell_lab_taur,sep="")

#### tsne ####

idx_dataset <- c(1:nrow(table_raw))
idx_dataset <- which(grepl("HMDB", table_raw$moleculeIds) & !grepl("fticr", table_raw$datasetName))
mz_dataset <- sort(unique(round(table_raw[idx_dataset, ]$mz, 2)))
idx_mz_dataset_h <- which(round(mz, 2) %in% round(mz_dataset, 2))

idx_dataset <- c(1:nrow(table_raw))
idx_dataset <- which(grepl("C", table_raw$moleculeIds) & !grepl("fticr", table_raw$datasetName))
mz_dataset <- sort(unique(round(table_raw[idx_dataset, ]$mz, 2)))
idx_mz_dataset_k <- which(round(mz, 2) %in% round(mz_dataset, 2))

idx_mz_dataset <- intersect(idx_mz_dataset_k, idx_mz_dataset_h)

colz_clus <- colorRampPalette(c("gray","green","magenta","#804080"))(4)
colz <- colorRampPalette(c("darkgreen","darkblue","orange","brown"))(4)

idx_clus <- idx_mz_dataset

set.seed(42)
tsne_native <- M3C::tsne(t(cells_all_mz_mat)[idx_clus,],labels=as.factor(repeated_prefixes),seed=42,
                         colvec = colz, axistextsize = 18,
                         legendtextsize = 18, dotsize = 2, controlscale=T,scale=3)

tsne_clus <- M3C::tsne(t(cells_all_mz_mat)[idx_clus,],labels=as.factor(cell_lab_comb),seed=42,
                       colvec = colz_clus, axistextsize = 18,
                       legendtextsize = 18, dotsize = 2, controlscale=T,scale=3)

tsne_native
tsne_clus


#### volcano ####

volc_plot <- function(labs, str_a, str_b, ms_mat_cardinal, pad_thresh = 0.05, log2fc_thresh = 0.2, idx_mz_dataset, cohen_d_thresh = 0.2, color_neg = "green", color_pos = "magenta", color_neutral = "gray",samp_ratio=1) {
  idx_A <- which(grepl(str_a, labs))
  idx_B <- which(grepl(str_b, labs))
  
  ms_mat_cardinal[is.nan(ms_mat_cardinal)] <- 0
  idx_all <- c(idx_A, idx_B)
  
  log2fc_mat <- rep(0, ncol(ms_mat_cardinal))
  pval_mat <- rep(0, ncol(ms_mat_cardinal))
  padj_mat <- rep(0, ncol(ms_mat_cardinal))
  cohend_mat <- rep(0, ncol(ms_mat_cardinal))
  
  set.seed(42)
  
  for (g in 1:10) {
    
    idx_As <- sample(idx_A, size = length(idx_A) * samp_ratio)
    idx_Bs <- sample(idx_B, size = length(idx_B) * samp_ratio)
    
    idx_keep_volc <- c(idx_As, idx_Bs)
    
    pval <- c()
    padj <- c()
    cohen_d <- c()
    
    A_mat <- ms_mat_cardinal[idx_As, ]
    B_mat <- ms_mat_cardinal[idx_Bs, ]
    
    for (k in 1:length(mz)) {
      if (var(A_mat[, k]) == 0 || var(B_mat[, k]) == 0) {
        pval[k] <- 1
        padj[k] <- 1
        cohen_d[k] <- 0
        next
      }
      
      pval[k] <- t.test(x = A_mat[, k], y = B_mat[, k])$p.value
      padj[k] <- p.adjust(pval[k], method = "BH", n = length(mz))
      cohen_d[k] <- effsize::cohen.d(A_mat[, k], B_mat[, k], pooled = FALSE)$estimate
    }
    
    mean_A <- colMeans(A_mat)
    mean_B <- colMeans(B_mat)
    log2foldchange <- mean_B - mean_A
    
    log2foldchange[is.nan(log2foldchange)] <- 0
    padj[is.nan(padj)] <- 1
    padj[padj < 1e-100] <- 1e-100
    
    log2foldchange_tmp <- c()
    for (fc in log2foldchange) {
      if (!is.infinite(fc)) {
        log2foldchange_tmp <- c(log2foldchange_tmp, fc)
      } else {
        if (fc <= 0) {
          log2foldchange_tmp <- c(log2foldchange_tmp, -10)
        } else {
          log2foldchange_tmp <- c(log2foldchange_tmp, 10)
        }
      }
    }
    
    log2foldchange <- log2foldchange_tmp
    log2foldchange[log2foldchange > 5] <- 5
    log2foldchange[log2foldchange <= -5] <- -5
    
    log2fc_mat <- rbind(log2fc_mat, log2foldchange)
    pval_mat <- rbind(pval_mat, pval)
    padj_mat <- rbind(padj_mat, padj)
    cohend_mat <- rbind(cohend_mat, cohen_d)
  }
  
  log2fc_mat <- log2fc_mat[2:11, ]
  pval_mat <- pval_mat[2:11, ]
  cohend_mat <- cohend_mat[2:11, ]
  
  log2folgchange <- colMeans(log2fc_mat)
  padj <- colMeans(pval_mat)
  padj <- colMeans(padj_mat)
  cohen_d <- colMeans(cohend_mat) * -1
  
  s_volc <- as.data.frame(cbind(as.numeric(colnames(ms_mat_cardinal)), log2foldchange, pval, padj, cohen_d))
  s_volc$cohen_d[which(s_volc$cohen_d >= 1)] <- 1
  s_volc$cohen_d[which(s_volc$cohen_d <= -1)] <- -1
  s_volc <- s_volc[idx_mz_dataset, ]
  labelz <- round(mz, 2)[idx_mz_dataset]
  custcols <- rep(color_neutral, nrow(s_volc))
  
  custcols[which(s_volc$cohen_d <= -cohen_d_thresh & s_volc$padj <= pad_thresh)] <- color_neg
  custcols[which(s_volc$cohen_d >= cohen_d_thresh & s_volc$padj <= pad_thresh)] <- color_pos
  s_volc$custom_color <- custcols
  
  names(custcols) <- rownames(s_volc)
  
  volc_plot <- EnhancedVolcano::EnhancedVolcano(
    s_volc,
    lab = labelz,
    x = "cohen_d",
    y = "padj",
    xlab = "Cohen's D",
    pCutoff = (pad_thresh),
    FCcutoff = cohen_d_thresh,
    xlim = range(-1.2, 1.2),
    title = "Comparison Plot",
    subtitle = "",
    caption = paste0("total = ", nrow(s_volc), " m/z values"),
    pointSize = 4,
    legendPosition = "none",
    labSize = 8,
    colCustom = custcols
  )
  
  mz_a <- labelz[which(s_volc$cohen_d <= -cohen_d_thresh & s_volc$padj <= pad_thresh)]
  mz_b <- labelz[which(s_volc$cohen_d >= cohen_d_thresh & s_volc$padj <= pad_thresh)]
  
  list(volc_plot=volc_plot,s_volc = s_volc, mz_a = mz_a, mz_b = mz_b)
}

idx_s <- which(grepl(repeated_prefixes, pattern="SIMA9"))
idx_hipsc <- which(grepl(repeated_prefixes, pattern="hiPSC"))

volc_plot_itataur_s <- volc_plot(labs=cell_lab_comb[idx_s], "Ita-Taur\\+","Ita\\+Taur-", idx_mz_dataset = idx_mz_dataset, ms_mat_cardinal=cells_all_mz_mat[idx_s,],color_neg = "green", color_pos = "magenta", color_neutral = "gray")
volc_plot_native_s <- volc_plot(labs=cell_lab[idx_s], "VEH","LPS", idx_mz_dataset = idx_mz_dataset, ms_mat_cardinal=cells_all_mz_mat[idx_s,],color_neg = "darkblue", color_pos = "brown", color_neutral = "gray")


volc_plot_itataur_s$volc_plot
volc_plot_native_s$volc_plot


volc_plot_itataur_hipsc <- volc_plot(labs=cell_lab_comb[idx_hipsc], "Ita-Taur\\+","Ita\\+Taur-",idx_mz_dataset = idx_mz_dataset, ms_mat_cardinal=cells_all_mz_mat[idx_hipsc,],color_neg = "green", color_pos = "magenta", color_neutral = "gray")
volc_plot_native_hipsc <- volc_plot(labs=cell_lab[idx_hipsc], "VEH","LPS", idx_mz_dataset = idx_mz_dataset, ms_mat_cardinal=cells_all_mz_mat[idx_hipsc,],color_neg = "darkblue", color_pos = "brown", color_neutral = "gray")

#### Slice Culture Data ####

collect_SC_info <- function(path_SC,path_spotlist_SC,scaleing=TRUE){
  msData_SC <- Cardinal::readMSIData(path_SC, attach.only=TRUE, resolution = res_ppm, units = "ppm", mass.range = mass_range)
  
  msData_SC_fil <- msData_SC
  
  tmp_tab_SC <- read.table(path_spotlist_SC)
  simple_labels_spotlist_SC <- as.character(tmp_tab_SC$V4)
  
  simple_labels_SC <- as.factor(simple_labels_spotlist_SC)
  
  norm_peak_idx <- which.min(abs(normal_peak-msData_SC_fil@featureData@mz))#
  norm_peak_idx_2 <- which.min(abs(normal_peak_2-msData_SC_fil@featureData@mz))
  
  ####
  ms_mat_SC <- as.matrix(t(Cardinal::spectra(msData_SC_fil)))
  dim(ms_mat_SC)
  
  norm_peak_idx_SC <- which.min(abs(normal_peak-msData_SC_fil@featureData@mz))
  norm_peak_idx_2_SC <- which.min(abs(normal_peak_2-msData_SC_fil@featureData@mz))
  msData_SC_fil@featureData@mz[norm_peak_idx_SC]
  
  
  ms_mat_normed_SC <- ms_mat_SC/(ms_mat_SC[,norm_peak_idx_SC])
  
  
  ms_mat_normed_SC[is.nan(ms_mat_normed_SC)] <- 0
  ms_mat_normed_SC[is.infinite(ms_mat_normed_SC)] <- 0
  
  ms_mat_normed_SC <- ms_mat_normed_SC[,-c(norm_peak_idx_SC,norm_peak_idx_2_SC)]
  
  if(scaleing==TRUE){
    
    idx_scale <- which(!grepl("Brain", simple_labels_SC) & !grepl("CLO", simple_labels_SC))
    means <- colMeans(ms_mat_normed_SC[idx_scale,])
    sds <- apply(ms_mat_normed_SC[idx_scale,], 2, sd)
    
    print(length(idx_scale))
    print(dim(ms_mat_normed_SC))
    print(levels(simple_labels_SC))
    for(t in c(1:nrow(ms_mat_normed_SC))){
      ms_mat_normed_SC[t,] <- (ms_mat_normed_SC[t,]-means)/sds
    }
  }
  
  ms_mat_normed_SC[is.nan(ms_mat_normed_SC)] <- 0
  ms_mat_normed_SC[is.infinite(ms_mat_normed_SC)] <- 0
  
  return(list(ms_mat_normed_SC,simple_labels_SC,msData_SC_fil))
  
}

SC_info_1 <- collect_SC_info(path_SC = "U:/James/HMDBE_KEGG_10/230524_SC_FINAL.imzMl",path_spotlist_SC = "U:/James/spotlists/230524_SC_spotlist.txt",scaleing=T)
SC_info_2 <- collect_SC_info(path_SC = "U:/James/HMDBE_KEGG_10/230601_SC_FINAL.imzMl",path_spotlist_SC = "U:/James/spotlists/230601_SC_spotlist.txt",scaleing=T)
SC_info_3 <- collect_SC_info(path_SC = "U:/James/HMDBE_KEGG_10/230728_SC_FINAL.imzML",path_spotlist_SC = "U:/James/spotlists/230728_SC_spotlist.txt",scaleing=T)
SC_info_4 <- collect_SC_info(path_SC = "U:/James/HMDBE_KEGG_10/230802_SC_D_FINAL.imzMl",path_spotlist_SC = "U:/James/spotlists/230802_SC_D_spotlist.txt",scaleing=T)
SC_info_5 <- collect_SC_info(path_SC = "U:/James/HMDBE_KEGG_10/230802_SC_D_10_FINAL.imzMl",path_spotlist_SC = "U:/James/spotlists/230802_SC_D_10_spotlist.txt",scaleing=T)
SC_info_6 <- collect_SC_info(path_SC = "U:/James/HMDBE_KEGG_10/230908_Slide_F_FINAL.imzML",path_spotlist_SC = "U:/James/spotlists/230908_Slide_F_spotlist_2.txt",scaleing=T)

SC_info_list <- list(SC_info_1,SC_info_2,SC_info_3,SC_info_4,
                     SC_info_5,SC_info_6)

rm(SC_info_1,SC_info_2,SC_info_3,SC_info_4,SC_info_5,SC_info_6)

indices_SC <- which(sapply(SC_info_list, function(x) length(x) > 0))

indices_SC <- c(1,2,3,4,5,6)

simple_labels_SC_all <- as.factor(unlist(sapply(SC_info_list[indices_SC], function(x) as.character(x[[2]]))))
ms_mat_normed_SC_all <- do.call(rbind, lapply(SC_info_list[indices_SC], function(x) x[[1]]))

gc()

idx_of_i_clo_SC <- which(grepl("CLO", simple_labels_SC_all))
idx_of_i_veh_SC <- which(grepl("VEH", simple_labels_SC_all))
idx_of_i_lps_SC <-  which(grepl("LPS", simple_labels_SC_all))

set.seed(42)

idx_B <- sample(idx_of_i_lps_SC,size=10000)
idx_A <- c(sample(idx_of_i_veh_SC,size=10000))

volc_plot_SC <- volc_plot(labs=simple_labels_SC_all[c(idx_A,idx_B)], "VEH","LPS",idx_mz_dataset = idx_mz_dataset, ms_mat_cardinal=ms_mat_normed_SC_all[c(idx_A,idx_B),],color_neg = "#F8A1D1", color_pos = "brown", color_neutral = "gray",samp_ratio=0.1)

volc_plot_SC$volc_plot

mz_veh_SC <- volc_plot_SC$mz_a
mz_lps_SC <- volc_plot_SC$mz_b

set.seed(42)
idx_B_c <- sample(idx_of_i_clo_SC,size=10000) 
idx_A_c <- c(sample(idx_of_i_veh_SC,size=5000),sample(idx_of_i_lps_SC,size=5000))


simple_labels_MG <- as.vector(simple_labels_SC_all[c(idx_A_c,idx_B_c)])
simple_labels_MG[grepl(simple_labels_MG,pattern="VEH")] <- "MG"
simple_labels_MG[grepl(simple_labels_MG,pattern="LPS")] <- "MG"
simple_labels_MG <- as.factor(simple_labels_MG)

volc_plot_SC_clo <- volc_plot(labs=simple_labels_MG, "MG","CLO", idx_mz_dataset = idx_mz_dataset, ms_mat_cardinal=ms_mat_normed_SC_all[c(idx_A_c,idx_B_c),],color_neg = "darkblue", color_pos = "yellow", color_neutral = "gray",samp_ratio=0.1)

volc_plot_SC_clo$volc_plot

mz_mg_SC <- volc_plot_SC_clo$mz_a
mz_clo_SC <- volc_plot_SC_clo$mz_b


df_cohen <- as.data.frame(cbind("mz" = round(mz,2)[idx_mz_dataset], "Slice Culture" = volc_plot_SC$s_volc$cohen_d, "SIMA9" = volc_plot_itataur_s$s_volc$cohen_d, "hiPSC" = volc_plot_itataur_hipsc$s_volc$cohen_d))
df_cohen_melted <- reshape2::melt(df_cohen, id.vars = 'mz')

mz_oi <- sort(c(mz_veh_SC,mz_lps_SC))

df_cohen_filtered <- df_cohen_melted[df_cohen_melted$mz %in% mz_oi, ]

ggplot(df_cohen_filtered, aes(x = factor(mz, levels = (mz_oi)), y = variable)) + 
  geom_point(aes(fill = (value)), shape = 21, size = 20, stroke = 0.5, color = "grey") +
  scale_fill_gradientn(colors = c("green", "lightgreen", "white", "lightcoral", "red"), values = scales::rescale(c(-1, -0.5, 0, 0.5, 1)))+
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(angle = 0, hjust = 1),
    text = element_text(size = 25, face = "bold"),
    legend.position = "right"
  ) +
  labs(fill = "Higher in LPS")

data <- list("LPS & VEH" = mz_mg_SC,"CLO" = mz_clo_SC, "VEH" = mz_veh_SC, "LPS" = mz_lps_SC)
ggvenn::ggvenn(data, fill_color  = c("darkblue", "#DAA520", "#F8A1D1", "brown"), fill_alpha  = 0.3, label = "none", text_size = 15,show_percentage = F,set_name_color = "white",)

