# Visualize tissue coordinates ===============================================
#' Visualizing tissue coordinates
#'
#' @description This will generate a ggplot starting from a \code{DimRedPlot} object.
#' See \code{\link{get_SpotPlot.sce}} to start from a \code{SingleCellExperiment}
#' object.
#'
#' @param drp_object a \code{\link{DimRedPlot}} object
#' @param X name of the drp_object column that contains the values for the x axis
#' @param Y name of the drp_object column that contains the values for the y axis
#' @param which_quantile default: 1
#' @param color_by specify the column of \code{drp_object$plot_data} whose
#' entries will be used to assign either discrete or continous color schemes.
#' @param color_by_factor if set to TRUE, the values supplied to \code{color_by}
#' will be interpreted as factors (as long as there are fewer than 14 individual
#' values). Default: FALSE
#' @param shape_by the factor to assign the shape of the points to; max. 10 levels
#' are accepted for \code{shape_by}. If you want to specify the (same) shape for
#' \emph{all} points, use an integer (default shape is 19).
#' @param size_by the factor to assign the size of the points to;
#' default: \code{NULL}. If you want to specify the (same) size for \emph{all}
#' points, use one integer (the default size is 3).
#' @param label_by factor (part of colData(sce)) used to label groups of cells.
#' A typical example would be "cluster" or "condition". Default: NULL.
#' @param label_size default: 4
#' @param circle_by specify the column of \code{drp_object$plot_data} whose
#' entires will be interpreted as factors and used to draw circles
#' around the points belonging to the same factor instance (e.g. "condition").
#' Alternatively, supply a list (!) where \code{list$cells} indicates the cell names
#' (= rownames of \code{drp_object$plot_data}) around which a circle (without fill)
#' will be drawn. Example: \code{list(cells = head(colnames(tiny_sce))}
#' @param ignore_drp_labels If set to TRUE, \emph{all} entries of
#' \code{drp_object$factors} will be ignored. If you want to ignore only selected
#' labels, indicate their name (e.g. c("color_by", "exprs_val_type")).
#' The default setting (\code{FALSE}) will make use of the information from
#' the drp_object \emph{unless} specified here via \code{color_by}, \code{shape_by},
#' \code{size_by}, \code{circle_by}.
#' @param theme_size font size for the plot
#' @param legend choose wether there should be no legend ("none"), a legend for
#' every factor or whether the legend will only show up for those factors that
#' have more than 1 level ("auto", default setting).
#' @param alpha define the opacity of the points, default: .65
#' @param remove_rug Define whether the bars at both axes that
#' represent the two 1d marginal distributions should be removed. Default: FALSE.
#' @param set_colors Set to \code{FALSE} if you want to add your own color scheme
#' to the resulting plot. The default behavior is to try to automatically assign
#' an optimized coloring scheme. Default: \code{TRUE}.
#' @param set_fill_colors Set to \code{FALSE} if you want to add your own color
#' scheme to the resulting plot \emph{for the background circles}. The default
#' behavior is to try to automatically assign an optimized coloring scheme.
#'
#' @seealso \code{\link{get_SpotPlot.sce}}, which generates the DimRedPlot
#' object that will be used here.
#'
#' @return \code{ggplot2} object (a plot)
#'
#' @export
#'
plt.SpotPlot <- function(drp_object, X="x_axs", Y="y_axs",
  color_by=NULL, color_by_factor = FALSE,
  shape_by=NULL, size_by=NULL, label_by = NULL,
  label_size = 4,
  circle_by = NULL, which_quantile = 1,
  ignore_drp_labels = FALSE,
  theme_size = 14, legend = "auto", alpha = .75,
  remove_rug = FALSE,
  set_colors = TRUE, set_fill_colors = TRUE){

  ##--set the defaults for the plot---------------------------------------------
  if(is.numeric(shape_by) & length(shape_by) == 1){
    shape_val <- shape_by
    shape_by <- NULL
    ignore_drp_labels <- c(ignore_drp_labels, "shape_by")
  }else{
    shape_val <- 19
  }

  if(is.numeric(size_by) & length(size_by) == 1 ){
    size_val <- size_by
    size_by <- NULL
    ignore_drp_labels <- c(ignore_drp_labels, "size_by")
  }else{
    size_val <- 3
  }

  ggplot2::update_geom_defaults("point", list(shape = shape_val, size = size_val))

  ##--extract basic plot information from DRP object----------------------------
  color_by <- scABC2:::fx.set_factors(drp_object, fctr = color_by,
    fctr_type = "color_by", ignore_drp_labels)
  shape_by <- scABC2:::fx.set_factors(drp_object, fctr = shape_by,
    fctr_type = "shape_by", ignore_drp_labels)
  size_by <- scABC2:::fx.set_factors(drp_object, fctr = size_by,
    fctr_type = "size_by", ignore_drp_labels)
  label_by <- scABC2:::fx.set_factors(drp_object, fctr = label_by,
    fctr_type = "label_by", ignore_drp_labels)

  if(!is.list(circle_by)){
    circle_by <- scABC2:::fx.set_factors(drp_object, fctr = circle_by,
      fctr_type = "circle_by", ignore_drp_labels)
  }

  ##--extract the data.frame needed for plotting--------------------------------
  df_to_plot <- drp_object$plot_data
  df_to_plot$key <- row.names(df_to_plot)

  ##--check legend argument ----------------------------------------------------
  legend <- match.arg(legend, c("auto", "none", "all"), several.ok = FALSE)

  ## if only one level for the variable, set legend to NULL
  if ( legend == "auto" ) {
    if (length(unique(df_to_plot$color_by)) == 1)
      color_by <- NULL
    if (length(unique(df_to_plot$shape_by)) == 1)
      shape_by  <- NULL
    if (length(unique(df_to_plot$size_by)) == 1)
      size_by <- NULL
    if (length(unique(df_to_plot$circle_by)) == 1)
      circle_by <- NULL
  }

  ##--generate base plot, apply colour_by, shape_by and size_by variables ------
  ## (NULL will simply be ignored by ggplot)
  if( !is.null(color_by) && !(color_by %in% names(df_to_plot)) ){color_by <- NULL}
  if( !is.null(shape_by) && !(shape_by %in% names(df_to_plot)) ){shape_by <- NULL}
  if( !is.null(size_by)  && !(size_by  %in% names(df_to_plot)) ){size_by <- NULL}

  ### if there's fewer than 10 individual numbers, use a discrete color scheme
  ### instead of a numeric one
  if(is.numeric(df_to_plot[,color_by]) & isTRUE(color_by_factor)){
    if( length(unique(df_to_plot[, color_by])) <= 14){
      df_to_plot[, color_by] <- as.character( df_to_plot[, color_by])
    }else{
      message("Although you indicated that you want to use the color_by option (", color_by, ") as factors, we've overruled that assessment because there are more than 14 individual values.")
    }
  }

  if(is.numeric(df_to_plot[,shape_by])){
    df_to_plot[, shape_by] <- factor( df_to_plot[, shape_by])
  }


  ##==BASE PLOT ================================================================
  plot_out <- ggplot(df_to_plot,
    aes_string(x = X,
      y = Y,
      key = "key",
      colour = ABCutilities::fx.parse_column_names(color_by),
      shape = ABCutilities::fx.parse_column_names(shape_by),
      size = ABCutilities::fx.parse_column_names(size_by)
    ))

  # return(plot_out)

  ##--rm x and y labels--------------------------------------------------------
   ld <- drp_object$label_data
  #
  # if( !is.null(ld$variance_pct) & !("variance_pct" %in% ignore_drp_labels) ){
  #   x_lab <- paste0(ld$x_lab, ": ", round(ld$variance_pct[1] * 100), "% variance")
  #   y_lab <- paste0(ld$y_lab, ": ", round(ld$variance_pct[2] * 100), "% variance")
  # }else{
  #   x_lab <- ld$x_lab
  #   y_lab <- ld$y_lab
  # }
  #
  plot_out <- plot_out + xlab("") +  ylab("")

  ##--add image----------------------------------------------------------------
  if(!is.null(drp_object$image_data)){
     img_data <- drp_object$image_data

      if(!"grob" %in% names(img_data)){
        warning("No grob in the image data object; will not print PNG.")
        }else{
          plot_out <- plot_out +
            facet_wrap(~slide) +
            geom_visium(data = img_data,
              aes(grob = grob), x = 0.5, y = 0.5) +
            coord_cartesian(expand=FALSE) + ## needed for alignment of the spots with the tissue
            xlim(0, max(img_data$width) )+
            ylim(max(img_data$height), 0)
        }
  }else{
       img_data <- NULL
       plot_out <- plot_out +
         xlim(min(df_to_plot[, X]), max(df_to_plot[, X])) +
         ylim(max(df_to_plot[, Y]), min(df_to_plot[, Y]))
  }

  ##--add background circles (do not move below geom_point!)--------------------
  if(!is.null(circle_by)){
    if(is.list(circle_by)){
      if(all(circle_by$cells %in% rownames(df_to_plot))){
        c.df <- df_to_plot[circle_by$cells,]
        plot_out <- plot_out + ggalt::geom_encircle(data = c.df, inherit.aes = FALSE,
          aes_string(x = X, y = Y))
      }else{
        warning("There will be no circles because not all of the cells supplied via circle_by$cells are part of the data.frame for plotting.")
      }
    }else{
      c.df <- scABC2:::fx.get_circle_df(df_to_plot,
        separate_by = circle_by,
        keep_quant= which_quantile)

      plot_out <- plot_out + ggalt::geom_encircle(data = c.df,
        aes_string(x = X,
          y = Y,
          fill = ABCutilities::fx.parse_column_names(circle_by)
        ),
        # color = NULL, shape = NULL, size = NULL)
        inherit.aes = FALSE,
        s_shape = 0.9, expand = 0.07, alpha = .2)
    }
  }

  ##--add points (leave this _after_ the background circles so that the points--
  ##--are plotted on top of the circles ----------------------------------------
  plot_out <- plot_out +  geom_point(alpha = alpha)

  ##--add labels ---------------------------------------------------------------
  if( !is.null(label_by)){
    plot_out <- plot_out +
      geom_text(aes_string(x = "x_center", y = "y_center", label = label_by),
        inherit.aes = FALSE, fontface = "bold", size = label_size)
  }

  ##--fix axes-------------------------------
  # if(!is.null(img_data)){
  #   xdims <- c(min(img_data$width), max(img_data$width))
  #  if(length(unique(xdims)) == 1){xdims[1] <- 0} ## assumes that the max of img_data$width is greater than 0
  #   ydims <- c(min(img_data$height), max(img_data$height))
  #   if(length(unique(ydims)) == 1){ydims[1] <- 0}
#
 #       plot_out <- plot_out + xlim( xdims[1], xdims[2] ) +
#              ylim(ydims[2] , ydims[1]) ## this flips the y axis to match the image
 # }

  ##--assign a sensible color scheme based on the features of color_by----------
  if(!is.null(color_by)){

    if(set_colors){
      if( !is.null(ld$exprs_val_type) & !("exprs_val_type" %in% ignore_drp_labels) ){
        color_label <- paste0(color_by, "\n(", ld$exprs_val_type, ")")
      }else{
        color_label <- color_by
      }

      plot_out <- ABCutilities::fx.resolve_plot_colors(plot_out, df_to_plot[,color_by],
        colour_by_name = color_label,
        fill = FALSE)
    }

    ## ignore alpha for the color legend if colors are not mapped to a continuous variable
    if( !is.numeric(df_to_plot[,color_by]) ){
      plot_out <- plot_out +  guides(color = guide_legend(override.aes = list(alpha = 1,
        size = 6,
        shape = 15),
        title = color_by
      ))
    }
  }

  ##--assign a sensible color scheme for the background circle------------------
  if( !is.null(circle_by) & !is.list(circle_by) & set_fill_colors){
    plot_out <- ABCutilities::fx.resolve_plot_colors(plot_out,
      df_to_plot[,circle_by],
      colour_by_name = circle_by,
      fill = TRUE)
    ## ignore alpha for the fill legend
    plot_out <- plot_out +
      guides(fill = guide_legend(override.aes = list(alpha = 1)))
  }

  ##--add sensible shapes-------------------------------------------------------
  if(!is.null(shape_by)){
    plot_out <- ABCutilities::fx.resolve_plot_shapes(plot_out,
      df_to_plot[,shape_by])
    ## ignore alpha for the scale legend
    plot_out <- plot_out +
      guides(shape = guide_legend(override.aes = list(alpha = 1, size = 6)))
  }

  ##--define plotting theme ----------------------------------------------------
  plot_out <- plot_out + theme_set(theme_bw(base_size = theme_size))

  plot_out <- plot_out +
          theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black"),
              #axis.text = element_blank(),
              #axis.ticks = element_blank()
          )

  ##--remove legend if so desired ----------------------------------------------
  if ( legend == "none" ){
    plot_out <- plot_out + theme(legend.position = "none")
  }

  ##--return plot---------------------------------------------------------------
  return(plot_out)
}

