
MiloGetIDs <- function(milo_object, milo_results, nhood_column = "Nhood", p_cutoff = 0.1, p_choice = "FDR"){

  milo_sig_results <- milo_results %>%
    dplyr::filter(!!sym(p_choice) < p_cutoff)
  if(nrow(milo_sig_results) < 1){
    stop("No statisticaly significant results in milo_results")
  }

  sig_nhdx_ix <- milo_sig_results[[nhood_column]]
  sig_nhds_id <- unlist(nhoodIndex(t_milo)[sig_nhdx_ix])

  milo_nhoods <- nhoods(milo_object)
  all_ids <- 1:nrow(milo_nhoods)
  nbhd_list <- lapply(seq_along(sig_nhds_id), function(i){
    temp_nbhd_id <- as.character(sig_nhds_id[i])
    chosen_indices <- all_ids[milo_nhoods[,temp_nbhd_id] == 1]
    return(chosen_indices)
  })
  names(nbhd_list) <- sig_nhds_id
  return(nbhd_list)

}


MiloWilcoxGenes <- function(seurat_object, milo_ids){
  milo_ids_names <- names(milo_ids)
  nbhd_genes <- lapply(seq_along(milo_ids), function(i){
    message("Finding marker genes for nhd ", i, " of ", length(milo_ids))
    temp_nbhd_ids <- milo_ids[[i]]
    seurat_object[["tempID"]] <- ifelse(colnames(seurat_object) %in% temp_nbhd_ids,
                                        yes = "in_cluster", no = "out_cluster")
    presto_results <- presto::wilcoxauc(seurat_object,
                                        group_by = "tempID",
                                        groups_use = c("in_cluster", "out_cluster")) %>%
      dplyr::filter(group == "in_cluster")
    presto_results$nhd_ix <- i
    presto_results$nhd_id <- milo_ids_names[i]
    return(presto_results)
  }) %>%
    do.call(rbind, .)
  return(nbhd_genes)
}


MiloHeatmap <- function(seurat_object, features, milo_ids, slot = "data", assay = "RNA", scale_rows = TRUE, max_zed = Inf, y_text_size = 10){
  features <- unique(features)
  non_milo_ids <- unique(unlist(milo_ids))
  milo_ids[["outgroup"]] <- non_milo_ids
  DefaultAssay(seurat_object) <- assay
  seurat_object <- seurat_object[features,]
  nbhd_expr <- lapply(seq_along(milo_ids), function(i){
    seurat_object[["tempID"]] <- ifelse(colnames(seurat_object) %in% milo_ids[[i]], "in_cluster", "out_cluster")
    my_data <- TrueAverageExpression(seurat_object, group.by = "tempID")
    my_data <- my_data[,"in_cluster", drop = FALSE]
    colnames(my_data) <- names(milo_ids)[[i]]
    return(my_data)
  }) %>% do.call(cbind, .)
  nbhd_expr <- as.data.frame(t(scale(t(nbhd_expr), center = TRUE, scale = TRUE)))
  nbhd_expr[nbhd_expr > max_zed] <- max_zed
  nbhd_expr[nbhd_expr < -max_zed] <- -max_zed
  hclust_results <- hclust(dist(nbhd_expr))
  feature_order <- hclust_results$labels[hclust_results$order]
  hclust_results_2 <- hclust(dist(t(nbhd_expr)))
  cluster_order <- hclust_results_2$labels[hclust_results_2$order]

  my_data_long <- nbhd_expr %>%
    as.data.frame %>%
    tibble::rownames_to_column(var = "feature") %>%
    pivot_longer(-feature, values_to = "Z_sc", names_to = "cluster") %>%
    dplyr::mutate(feature = factor(feature, levels = feature_order)) %>%
    dplyr::mutate(cluster = factor(cluster, levels = cluster_order))
  ggplot(my_data_long, aes(x = cluster, y = feature, fill = Z_sc))+
    geom_tile()+
    scale_fill_viridis_c()+
    ylab(NULL)+
    scale_y_discrete(position = "left", expand = c(0.01, 0.01))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = y_text_size),
          plot.margin = margin(l = 0.1),
          panel.background = element_rect(fill = "white"),
          axis.text.y = element_text(size = y_text_size),
          legend.text = element_text(size = y_text_size * 0.75),
          legend.title = element_text(size = y_text_size * 0.75),
          legend.key.width = unit(5, "pt"),
          legend.key.height = unit(10, "pt"))+
    ylab(NULL)+
    xlab(NULL)
}

SeuratAddNhds <- function(seurat_object, nhd_list, nhd_directions, milo_results_name = "milo_results"){
  nhd_df <- stack(nhd_list) %>%
    dplyr::mutate(nhd = as.numeric(as.character(ind)), cells = as.numeric(values)) %>%
    dplyr::select(nhd, cells)
  nhd_df$milo_results <- plyr::mapvalues(nhd_df$nhd, from = nhd_directions$nhd, to = nhd_directions$nhd_dir)
  nhd_dups <- nhd_df %>%
    dplyr::group_by(ind, directions) %>%
    dplyr::summarize()
  return(nhd_dups)
  #seurat_object[[milo_results_name]] <- ifelse(colnames(seurat_object) %in% nhd_df$)
}
