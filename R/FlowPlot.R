#' Plot a metadata selection in a biaxial FeatureScatter plot
#'
#' @param object Seurat object
#' @param feature1 The x-axis feature
#' @param feature2 The y-axis feature
#' @param assay Which assay
#' @param slot Which slot in the assay, typically data
#' @param pt.size Size of geom_point
#' @param noise.zero Whether to simulate values from 0 to a negative value,
#' done for visual purposes ONLY and for kernel density estimation when
#' overlaying a density visualization. Will simulate values using a normal
#' distribution centered at 0, take absolute values, and then scale from
#' noise.zero.bound to 0.
#' @param noise.zero.bound The negative value to which to simulate zeros to.
#' @param rasterize Whether to rasterize the resulting plot.
#' @details Scales the expression values for each feature in the specified assay
#' from 0 to 1, with an optional split.by argument
#'
#' @return Returns a biaxial plot utilizing ggpointdensity::geom_pointdensity
#'
#' @importFrom rlang %||%
#'
#' @export
FlowPlot <- function(object,
                     feature1,
                     feature2,
                     assay = "ADT",
                     slot = "data",
                     pt.size = 1,
                     noise.zero = FALSE,
                     noise.zero.bound = -0.01,
                     gg_bandwidth = 0.1,
                     rasterize = FALSE){
  DefaultAssay(object) <- assay
  plotting_data <- FetchData(object, vars = c(feature1, feature2), slot = slot)
  if(noise.zero){

    feature1_zeros <- sum(plotting_data[,feature1] == 0)
    feature1_sim_zeros <- -1 * abs(rnorm(feature1_zeros, 0, 1))
    feature1_sim_zeros <- scales::rescale(feature1_sim_zeros, to = c(noise.zero.bound, 0))
    plotting_data[,feature1][plotting_data[,feature1] == 0] <- feature1_sim_zeros

    feature2_zeros <- sum(plotting_data[,feature2] == 0)
    feature2_sim_zeros <- -1 * abs(rnorm(feature2_zeros, 0, 1))
    feature2_sim_zeros <- scales::rescale(feature2_sim_zeros, to = c(noise.zero.bound, 0))
    plotting_data[,feature2][plotting_data[,feature2] == 0] <- feature2_sim_zeros

  }

  avail_args <- names(as.list(args(ggpointdensity::geom_pointdensity)))

  if("method" %in% avail_args){
    g <- ggplot(plotting_data, aes(x = !!sym(feature1), y = !!sym(feature2)))+
      ggpointdensity::geom_pointdensity(size = pt.size,
                                        method = "kde2d",
                                        adjust = gg_bandwidth)+
      theme_classic()
  } else {
    message("geom_pointdensity devel version not installed. defaulting to geom_hex")
    g <- ggplot(plotting_data, aes(x = !!sym(feature1), y = !!sym(feature2)))+
      geom_hex()+
      theme_classic()
  }

  if(rasterize){
    ggrastr::rasterise(g)
  } else {
    g
  }


}
