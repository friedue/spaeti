#' Get the plot straight from an SCE object (wrapper fct)
#'
#' @description Wrapper around \code{\link{get_reducedDimPlot.sc}} and
#'  \code{\link{plt.DimRedPlot}}, which will enable direct plotting from an SCE
#'  object.
#'
#' @param object \code{SingleCellExperiment} object
#' @param which_reddim Deafult: "spots"
#' @param X directly supply the X coordinates, e.g. \code{pca_res$PC1}; if this
#' is used, \code{which_reddim} is ignored.
#' \code{X} should be of the same length (and order!) as \code{colData(object)}
#' (or \code{pData2.0(object)}).
#' @param Y directly supply the Y coordinates, e.g. \code{pca_res$PC2}.
#' Should be of the same length (and order!) as \code{colData(object)}.
#' \code{which_reddim} is ignored if \code{X} and  \code{Y} are specified.
#' @param image_data list containing the grobs of the PNGs and the max. height and
#' width of the images. It should also contain a column named "slide",
#' which should indicate a unique name per slide. This column will be used
#' to facet the final plot. Make sure that its values are actually part of
#' the sce.object's colData, too.
#' See \code{\link{reading_image_with_mods}} for details on how to generate an
#' appropriate list.
#' @param color_by entry of either \code{colData(object)}, \code{rowData(object)},
#' \code{rownames(object)}, which will be used to assign either discrete or
#' continous color schemes. Alternatively, you can supply a \code{list} with a
#' label stored in \code{list$title} and values stored in \code{list$result}
#' (\code{result} should be as long as \code{dim(object)[2]}. Default: \code{NULL}
#' @param shape_by the factor to assign the shape of the points to; max. 10 levels
#' are accepted for \code{shape_by}. If you want to specify the (same) shape for
#' \emph{all} points, use an integer (default shape is 19).
#' @param size_by the factor to assign the size of the points to;
#' default: \code{NULL}. If you want to specify the (same) size for \emph{all}
#' points, use one integer (the default size is 3).
#' @param label_by factor (part of colData(sce)) used to label groups of cells.
#' A typical example would be "cluster" or "condition". Default: NULL.
#' @param label_size Default: 4
#' @param circle_by specify the column of \code{drp_object$plot_data} whose
#' entires will be interpreted as factors and used to draw circles
#' around the points belonging to the same factor instance (e.g. "condition").
#' Alternatively, supply a list (!) where \code{list$cells} indicates the cell names
#' (= rownames of \code{drp_object$plot_data}) around which a circle (without fill)
#' will be drawn. Example: \code{list(cells = head(colnames(tiny_sce))}
#' @param circle_by Used for encircling points of interest; default: \code{NULL}.
#' You may want to define the quantiles of the points that will be encircled via
#' \code{which_quantile} (Default: 1). Alternatively, you can supply a list with
#' a single vector of cell names (\code{list$cells = head(colnames(sce))}), which
#' will circle the respective cells with an empty circle.
#' @param exprs_values types of expression values that will be used for plotting
#' if a gene name is being specified in either \code{color_by}, \code{shape_by},
#' \code{size_by}
#' @param add_cell_info if you want to add columns for certain entries of
#' \code{colData(object)}, specify their names here. Default: NULL
#' @param ... additional parameters for \code{\link{plt.SpotPlot}} such as
#' \code{color_by_factor}, \code{theme_size}, \code{legend}, \code{alpha},
#'  \code{remove_rug}, \code{set_colors}, \code{set_fill_colors}.
#'
#' @examples
#' data_dir <- "/Users/frd2007/Documents/Projects/2019-03_Sofia_SchietingerLab/Visium/data/"
#' samples <- c("pN3A","pN3B", "pN5A", "pN5B")
#'
#' imgs <- lapply(samples, function(x){
#'    reading_image_with_mods(
#'        paste0(data_dir, x, "/spatial/"),
#'        slide_name = x, testing_mods = FALSE)
#'    }) %>% do.call(rbind, .)
#'
#' # make the plot
#' plot_Spot_from_sce(sce.norm[, sce.norm$cell.labels == "cells" & sce.norm$N.id == i],
#'        color_by = "Ins1", exprs_values = "log1p_sctransform",
#'        plot_PNG = TRUE, image_data = imgs,
#'        size_by = 1,
#'        add_cell_info = "slide", shape_by = "cell.status"
#' @seealso \code{get_SpotPlotDRP} and \code{plt.SpotPlot}, for which this
#' function is a wrapper.
#'
#' @export
#'
 plot_Spot_from_sce <- function(sce.object, which_reddim = "spots",
  X = NULL, Y = NULL, image_data = NULL, plot_PNG = FALSE,
  color_by=NULL, shape_by=NULL, size_by=NULL, label_by = NULL, label_size = 4,
  circle_by = NULL, exprs_values=NULL,
  add_cell_info = NULL, ...){

  ## generate DRP object --------------------------------------
  if(is.numeric(size_by)){
    size_val1 <- NULL
  }else{
    size_val1 <- size_by
  }

  if(is.numeric(shape_by)){
    shape_val1 <- NULL
  }else{
    shape_val1 <- shape_by
  }

    ## determine whether the PNG can be added or not
   # if(is.null(image_data)){
  #      image_data <- sce.object@metadata$img_measures
  #      plot_PNG <- FALSE
  #  }
    if(!is.null(image_data)){
      if(!"grob" %in% names(image_data)){
        image_data <- NULL
        warning("No image data found.")
      }
    }

  drp <- get_SpotPlotDRP(sce.object,
      which_reddim = which_reddim,
      X = X, Y = Y,
      image_data = image_data,
      color_by = color_by,
      shape_by = shape_val1,
      size_by = size_val1,
      label_by = label_by,
      circle_by = circle_by,
      exprs_values = exprs_values,
      add_cell_info = add_cell_info)

   ## make the plot ----------------------------------------------
  P <- plt.SpotPlot(drp,
      color_by = drp$factors$color_by,
      shape_by = shape_by,
      size_by = size_by, label_by=label_by, label_size = label_size, circle_by = circle_by, ...)
  return(P)
}
