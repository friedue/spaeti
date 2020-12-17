# Prepare DRP object for plotting tSNE or PCA ==================================
#' Prepare the "DimRedPlot" object for plotting spot information
#'
#' @description Wrapper function for extracting the X and Y coordinates for the
#' plot as well as the variables for coloring, shape assignment etc.
#' This function expects a SingleCellExperiment object.
#'
#' @details
#' \emph{which_reddim} The reducedDims slot of the SingleCellExperiment object
#' is a list (\code{\link{SimpleList}}. To know which coordinates to extract,
#' the accessor for the specific reducedDims list member must be given.
#'
#' @param object \code{SingleCellExperiment} object
#' @param which_reddim accessor to retrieve the list of choice for the reduced
#' dimension coordinates from \code{reducedDims(object)}. See \code{details}.
#' @param X directly supply the X coordinates, e.g. \code{pca_res$PC1}; if this
#' is used, \code{which_reddim} is ignored.
#' \code{X} should be of the same length (and order!) as \code{colData(object)}
#' (or \code{pData2.0(object)}).
#' @param Y directly supply the Y coordinates, e.g. \code{pca_res$PC2}.
#' Should be of the same length (and order!) as \code{colData(object)}.
#' \code{which_reddim} and \code{which_pcs} will be ignored if \code{X} and
#' \code{Y} are specified.
#' @param image_data list with the grob, height and width of the PNG images
#' of the corresponding tissue (usually the low resolution version). It should
#' also contain a column named "slide", which should indicate a unique name per
#' slide. This column will be used to facet the final plot. Make sure that its
#' values are actually part of the sce.object's colData, too.
#' See \code{\link{reading_image_with_mods}} for details on how to generate an
#' appropriate list.
#' @param color_by entry of either \code{colData(object)}, \code{rowData(object)},
#' \code{rownames(object)}, which will be used to assign either discrete or
#' continous color schemes. Typically, only one value is accepted here, except for
#' gene names, where, if multiple genes are indicated, the median expression across
#' all genes will be shown per cell. Alternatively, you can supply a \code{list}
#' with a label stored in \code{list$title} and values stored in \code{list$result}.
#' (\code{result} should be as long as \code{dim(object)[2]}. Default: \code{NULL}
#' @param shape_by default: \code{NULL}
#' @param size_by default: \code{NULL}
#' @param label_by default: \code{NULL}
#' @param circle_by Used for encircling points of interest; default: \code{NULL}.
#' If \code{color_by} is a list, this will be ignored here because the list will
#' be used to subset the final data.frame once the actual plot is being constructed.
#' @param exprs_values types of expression values that will be used for plotting if a gene name is being
#' specified in either \code{color_by}, \code{shape_by}, \code{size_by}
#' @param add_cell_info if you want to add columns for certain entries of
#' \code{colData(object)}, specify their names here. Default: NULL
#'
#' @return DimRedPlot object: a simple S3-based object, i.e. a list with:
#' \enumerate{
#' \item \code{plot_data}: data.frame fit for ggplot2-based plotting
#' \item \code{label_data}: additional stuff that may be used for labelling the
#' x and y axes of the final plot
#' \item \code{factors}: names of the plotting variables, e.g. which values
#' were used for defining the coloring etc.
#' }
#'
#' @seealso \code{\link{plt.SpotPlot}} for plotting the resulting DRP and
#' \code{\link{plot_Spot_from_sce}} for a wrapper to generate the DRP and the
#' plot in one go starting from an SCE and image_data
#'
#' @export
#'
#' @import scater
#'
get_SpotPlotDRP <- function(object, which_reddim = "spots",
    X = NULL, Y = NULL, image_data = NULL,
    color_by=NULL, shape_by=NULL, size_by=NULL,
    label_by=NULL, circle_by = NULL,
    exprs_values="counts",
    add_cell_info = NULL){

    if(is.null(which_reddim)  & ( is.null(X) | is.null(Y) ) ){
        stop("You must either provide the X and Y coordinates directly or specify the
         accessor for the reduced dimension coordinates via which_reddim.")
    }

    ## Define data.frame for plotting --------------------------------------------

    ## Either by using the cell names and the X and Y coordinates directly specified
    ## or by using the results stored within reducedDimension(object)
    if(!is.null(which_reddim)){
        use_reddim <- TRUE
    }

    if(!is.null(X) | !is.null(Y) ){
        if(any(is.null(X), is.null(Y))){
            # if only one of the axes is specified, we may resort to using which_reddim
            if(!is.null(which_reddim)){
                warning("You only supplied either X or Y. We'll be using which_reddim instead.")
            }else{
                stop("You must supply either which_reddim or both X and Y. Currently you're
             only providing either X or Y.")
            }
        }else{
            use_reddim <- FALSE
        }
    }

    ## let's define the data.frame using the specific X and Y coordinates---------
    if(!use_reddim){
        xlen <- scABC2:::fx.return_axis_length(X)
        ylen <- scABC2:::fx.return_axis_length(Y)
        if(xlen != length(colnames(object)) | ylen != length(colnames(object))){
            stop("X and Y must be of the same length (and order!) as colnames(object).")
        }else{

            ## make data.frame
            df_to_plot <- data.frame(x_axs = X, y_axs = Y, row.names = colnames(object))

            ## extract axes labels
            if( all(names(df_to_plot) == c("x_axs", "y_axs")) ){
                axs_labs <- NULL
            }else{
                axs_labs <- names(df_to_plot)
                names(df_to_plot) <- c("x_axs", "y_axs")
            }
            rd <- NULL
        }
    }else{ ## Extracting coordinates from the reducedDims() slot -----------------

        if(!any(class(object) %in% "SingleCellExperiment")){
            stop("The object is not of class SingleCellExperiment, therefore we cannot retrieve
           the reduced dimensions via reducedDims. Please provide X and Y directly.")
        }else{
            red_dim_all <- reducedDims(object)
        }

        ## check that there was actual data in that slot.................
        if(length(red_dim_all) == 0){
            stop("The object doesn't seem to have the reducedDim() slot filled.")
        }else{
            rd <- red_dim_all[[which_reddim]]
        }

        ## check the accessor....................
        if(is.null(rd)){
            stop("The accessor you provided for retrieving the coordinates from the
           reducedDims() list did not return any entries. Check that the name
           or index exists for reducedDims(object). Otherwise use X and Y to
           specify the vectors for the coordinates manually.")
        }

        ## finally, make the data.frame with coordinates for the components of interest...........
        ## subset the redDim data.frame
        which_pcs <- c(1,2)
        df_to_plot <- data.frame(rd[, which_pcs, drop=FALSE])

        ## retrieve the original column names
        axs_labs <- names(df_to_plot)

        ## assign the standardized column names
        names(df_to_plot) <- c("x_axs","y_axs")
    }

    ## flip the y axis to match the tissue picture ---------------------
    #if(is.null(image_data)){
    #    df_to_plot$y_axs <- max(df_to_plot$y_axs) - df_to_plot$y_axs + min(df_to_plot$y_axs)
    #} #### doesn't work yet #####

    ## add meta-data to the X and Y coordinates, which will be used to assign color, shape etc. ------------
    if( !is.null(color_by) | !is.null(shape_by)  | !is.null(size_by) | !is.null(circle_by) ){
        df_to_plot <- scABC2:::fx.add_factors_to_redDim_df(object, df_to_plot,
            color_by, shape_by,
            size_by, circle_by,
            exprs_values)
    }

    ## return the data.frame with more info than what would normally be added based on color_by etc.
    ## this is useful if the data.frame is to be used in various ways without having to re-generate
    ## it all the time (e.g., extracting the values for individual genes can be time-consuming)
    if( !is.null(add_cell_info) ){
        df_to_plot <- scABC2:::fx.add_cellInfo_to_redDim_df(object, df_to_plot, add_cell_info)
    }

    ## add centerX and centerY coordinates for labels of clusters
    if(!is.null(label_by)){

        if(!(label_by %in% names(df_to_plot))){
            df_to_plot <- scABC2:::fx.add_cellInfo_to_redDim_df(object, df_to_plot, label_by)
        }

        if(!(label_by %in% names(df_to_plot))){ warning("The label you defined is not part of the SCE object")
        }else{
            df_to_plot <- scABC2:::fx.calc_label_positions(df_to_plot, label_by)
        }
    }

    ##? Do I NEED TO CHECK IMAGE_DATA HERE?
    ## if PNG is to be added, we have to make sure that the sample info is part of the df because we're going to facet on it
    if(!is.null(image_data)){
        mm <- as.matrix(colData(object))
        which_col <- names(which(apply(mm, 2, function(x) all(x %in% image_data$slide))))
        if(length(which_col) == 0){warning(paste("Could not find the colData entry in the SCE object that contains the slide names:", unique(image_data$slide)) )
            }else{
            df_to_plot <- scABC2:::fx.add_cellInfo_to_redDim_df(object, df_to_plot, which_col)
            names(df_to_plot)[names(df) == which_col] <- "slide"
            }
    }

    ## make sure there are no duplicated entries, e.g. if color and size are both
    ## mapped to the same variable [ggplot2 cannot handle that and will complain]
    df_to_plot <- df_to_plot[,unique(names(df_to_plot))]

    ## gather data for labels in the plot etc.------------------------------------
    ld <- list()
    ld$x_lab <- axs_labs[1]
    ld$y_lab <- axs_labs[2]
    #ld$variance_pct <- fx.get_percent_var(rd,  pcs_select = which_pcs )
    ld$exprs_val_type <- exprs_values
    #ld$dim_red_type <- dim_red_type

    ## remember which factors were used for what type of plotting variable--------
    fc <- list()
    if(is.list(color_by)){
        fc$color_by <- color_by$title
    }else{
        fc$color_by <- paste(color_by, collapse = ",")
    }
    fc$size_by <- size_by
    fc$shape_by <- shape_by
    fc$label_by <- label_by
    if(!is.list(circle_by)){
        fc$circle_by <- circle_by
    }
    if(!is.null(image_data)){
        fc$facet_by <- "slide"
    }

    ## return an S3 object--------------------------------------------------------
    DRP <- list(plot_data = df_to_plot,
        image_data = image_data,
        label_data = ld,
        factors = fc)
    class(DRP) <- "DimRedPlot"

    return(DRP)
}
