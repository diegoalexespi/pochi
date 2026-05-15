#' Exclude cells over a certan (n) median absolute deviation
#'
#' @param object Seurat object
#' @param n.mad (n) median absolute deviation
#' @details Based on 
#' https://matthieuxmoreau.github.io/EarlyPallialNeurogenesis/html-Reports/Quality_Control.html 
#' code. Exclude cells over a certan (n) median absolute deviation
#'
#' @return Returns a Seurat object without outlined (identified by n*mad)
#'
#' @importFrom rlang %||%
#'
#' @export
QC.n.mad <- function(Seurat.object, n.mad=3) {
  Cell.QC.Stat <- Seurat.object@meta.data
  print(nrow(Cell.QC.Stat))
  # high and low median absolute deviation (mad) thresholds to exclude outlier cells
  max.mito.thr <- median(Cell.QC.Stat$percent.mt) + n.mad*mad(Cell.QC.Stat$percent.mt)
  min.mito.thr <- median(Cell.QC.Stat$percent.mt) - n.mad*mad(Cell.QC.Stat$percent.mt)
  
  #Plot
  p1 <- ggplot(Cell.QC.Stat, aes(x=nFeature_RNA, y=percent.mt)) +
    geom_point() +
    geom_hline(aes(yintercept = max.mito.thr), colour = "red", linetype = 2) +
    geom_hline(aes(yintercept = min.mito.thr), colour = "red", linetype = 2) +
    annotate(geom = "text", label = paste0(as.numeric(table(Cell.QC.Stat$percent.mt > max.mito.thr | Cell.QC.Stat$percent.mt < min.mito.thr)[2])," cells removed\n", as.numeric(table(Cell.QC.Stat$percent.mt > max.mito.thr | Cell.QC.Stat$percent.mt < min.mito.thr)[1])," cells remain"), x = 6000, y = 0.1)
  
  Cell.QC.Stat <- Cell.QC.Stat %>% filter(percent.mt < max.mito.thr) %>% filter(percent.mt > min.mito.thr)
  
  # Set low and hight thresholds on the number of detected genes
  min.Genes.thr <- median(log10(Cell.QC.Stat$nFeature_RNA)) - n.mad*mad(log10(Cell.QC.Stat$nFeature_RNA))
  max.Genes.thr <- median(log10(Cell.QC.Stat$nFeature_RNA)) + n.mad*mad(log10(Cell.QC.Stat$nFeature_RNA))
  # Set hight threshold on the number of transcripts
  max.nUMI.thr <- median(log10(Cell.QC.Stat$nCount_RNA)) + n.mad*mad(log10(Cell.QC.Stat$nCount_RNA))
  
  p2 <- ggplot(Cell.QC.Stat, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
    geom_point() +
    geom_smooth(method="lm") +
    geom_hline(aes(yintercept = min.Genes.thr), colour = "green", linetype = 2) +
    geom_hline(aes(yintercept = max.Genes.thr), colour = "green", linetype = 2) +
    geom_vline(aes(xintercept = max.nUMI.thr), colour = "red", linetype = 2)
  
  p3 <- p1 / p2
  #p1=ggExtra::ggMarginal(p1, type = "histogram", fill="lightgrey", bins=100)
  #p2=ggExtra::ggMarginal(p2, type = "histogram", fill="lightgrey")
  #print(typeof(p1))
  #p3 <- p1 / p2
  p3 <- p3 + plot_annotation(title = names(Seurat.object), theme = theme(plot.title = element_text(hjust = 0.5)))
  print(p3)
  
  # Filter cells based on these thresholds
  
  Cell.QC.Stat <- Cell.QC.Stat %>% filter(log10(nFeature_RNA) > min.Genes.thr) %>% filter(log10(nCount_RNA) < max.nUMI.thr)
  print(nrow(Cell.QC.Stat))
  print("########")
  Seurat.object <- subset(Seurat.object, cells = rownames(Cell.QC.Stat))
  return(Seurat.object)
}
