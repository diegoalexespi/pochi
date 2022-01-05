RunPAGA <- function(object,
                    fuzzy.graph,
                    distance.graph,
                    pythondir = "~/miniconda3/envs/scanpy-env/",
                    group.by = NULL){


  message("Loading reticulate")
  require(reticulate)
  reticulate::use_python(pythondir)
  message("Importing scipy")
  scipy <- reticulate::import("scipy")
  message("Importing scanpy")
  scanpy <- reticulate::import("scanpy")

  temporary.dir <- tempdir()
  message("Using temporary directory ", temporary.dir)
  fuzzy_dir <- paste0(temporary.dir, "/fuzzy.mtx")
  distance_dir <- paste0(temporary.dir, "/distance.mtx")

  Matrix::writeMM(obj = as(object[[fuzzy.graph]], "dgCMatrix"), file = fuzzy_dir)
  Matrix::writeMM(obj = as(object[[distance.graph]], "dgCMatrix"), file = distance_dir)

  message("D")
  connectivities = scipy$io$mmread(fuzzy_dir)
  connectivities = connectivities$to_csr()

  message("Ds")

  distances = scipy$io$mmread(distance_dir)
  distances = distances$to_csr()

  message("Dss")

  file.remove(fuzzy_dir)
  file.remove(distance_dir)

  return("AH")


#
#   graph_dir="/Volumes/SanDiskSSD/bar-or_lab/projects/tonsil_project/analysis/sct_seurat/graphs_gcpc/"
#   seurat_connectivities = mmread(graph_dir+"fuzzy.mtx").tocsr()
#   seurat_cosines = mmread(graph_dir+"cosine.mtx").tocsr()
#   Bint.uns['neighbors'] = {"params" : {"n_neighbors" : 30, "method" : "annoy", "metric" : "cosine", "use_rep" : "X_pca"},
#     "distances" : seurat_cosines, "connectivities" : seurat_connectivities}
#
#   # see PAGA notebook here
#   # can then continue with running UMAP again, this time initializing positions with PAGA
#   # you should be able to implement reticulate here...
#
#   paga_pos <- read.csv(
#     "/Volumes/SanDiskSSD/bar-or_lab/projects/tonsil_project/analysis/sct_harmony_nomtreg_310/merge_prior/scanpy/paga_0.8.txt")
#   paga_pos$seurat_cluster <- factor(paga_pos$seurat_cluster, unique(paga_pos$seurat_cluster))
#   paga_matrix <- data.frame("seurat_cluster" = B.combined$harmony_annoy_nn_snn_res.0.8) %>%
#     dplyr::left_join(y = paga_pos, by = "seurat_cluster") %>%
#     dplyr::select(X0, X1)
#
#   dims.to.use <- 1:35
#   B.combined <- BuildAnnoyUMAP(B.combined,
#                                reduction = "harmony",
#                                dims = dims.to.use,
#                                reduction.name = "harmony_UMAP",
#                                graph.key = "harmony_annoy",
#                                reduction.key = "harmonyumap_",
#                                init.pos = as.matrix(paga_matrix))

}

