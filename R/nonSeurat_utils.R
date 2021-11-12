#' Determine overlap measure between a list of genes
#'
#' @param genelist Named list of gene vectors
#' @param distance.measure Which sort of distance measure to use (Jaccard by default)
#' @param as.distance Whether to transform similarities to distances
#' @param return.tidy Whether to return a tidy dataframe
#'
#' @return Returns a matrix of genelist distance measures
#'
#' @importFrom rlang %||%
#'
#' @export
FindJaccard <- function(genelist,
                        distance.measure = "jaccard",
                        as.distance = FALSE,
                        return.tidy = TRUE){
  genelist_data <- stack(genelist) %>%
    magrittr::set_colnames(c("gene", "geneset")) %>%
    dplyr::select(geneset, gene) %>%
    dplyr::mutate(occurrence = 1)

  genelist_wide <- genelist_data %>%
    tidyr::pivot_wider(id_cols = geneset, names_from = gene, values_from = occurrence, values_fill = 0) %>%
    as.data.frame() %>%
    tibble::column_to_rownames(var = "geneset")

  genelist_matrix <- proxy::dist(genelist_wide, method = distance.measure, diag = TRUE, convert_similarities = as.distance)

  if(return.tidy){
    genelist_matrix <- as.matrix(genelist_matrix)
    diag(genelist_matrix) <- 1
    genelist_matrix <- as.data.frame(genelist_matrix) %>% tibble::rownames_to_column(var = "genelist_A")
    genelist_matrix <- tidyr::pivot_longer(genelist_matrix, !genelist_A, names_to = "genelist_B")
    return(genelist_matrix)
  }

  return(genelist_matrix)
}




#' Read in a GMT file
#'
#' @param file_path path to .gmt file
#' @param file_sep which delimiter is used to separate entries in .gmt file
#'
#' @return Returns a data frame of the .gmt file terms and genes
#' @references code adapted from clusterProfiler::read.gmt
#'
#' @export
ReadGMT <- function(file_path = "", file_sep = ","){
  x <- readLines(file_path)
  res <- strsplit(x, file_sep)
  names(res) <- vapply(res, function(y) y[1], character(1))
  res <- lapply(res, "[", -c(1:2))
  ont2gene <- stack(res)
  ont2gene <- ont2gene[, c("ind", "values")]
  colnames(ont2gene) <- c("term", "gene")
  return(ont2gene)
}
