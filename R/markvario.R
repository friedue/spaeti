### IDENTIFY GENES WHOSE EXPRESSION CHANGES WITH LOCATION ======================
#' Run the mark variogram computation
#'
#' @description Wraps the functionality of \code{\link{spatstat::markvario}}.
#' Inspired by [Trendsceek (Endsgard, Johnsson, Sandberg, 2018)](https://www.nature.com/articles/nmeth.4634),
#' which models spatial transcriptomics data as a mark point process and computes a ‘variogram’, which
#'  identifies genes whose expression level is dependent on their spatial location.
#'  The calculation of the r-metric was copied from Seurat v. 3.1 (\code{ComputeRMetric}).
#'
#' @details This process calculates gamma(r) values measuring the dependence between
#'  two spots that are a certain "r" distance apart.
#' Spatial data can be made up of point patterns (location) with attached values (e.g. expression).
#' The mark variogram indicates how similar two marks (here: genes) within a certain distance from each other are.
#'
#' @param sce.obj SCE object with XY coordinates in reducedDim slot and normalized expression values
#' that can be used for testing whether individual genes show expression patterns that correlate with
#' the spatial location of the respective cells
#' @param coords_accessor Indicate which reducedDim slot should be accessed to obtain a matrix of x and y coordinates
#' of the cells
#' @param exprs_values Which assay should be used from the SCE object. Default: "log1p_sctransform"
#' @param genes if you want to select a subset of genes for the analysis, specify their names here. Default: NULL will
#' use all genes in the SCE object.
#' @param r.metric r value at which to report the "trans" value of the mark variogram. Default: 5
#' @param verbose Default: TRUE
#' @param ... Additional Arguments passed to \code{\link{spatstat::markvario}}
#'
#' @references
#'
#' Daniel Edsgärd, Per Johnsson & Rickard Sandberg.
#' Identification of spatial expression trends in single-cell gene expression data
#' \emph{Nature Methods}, 15, 339–342 (2018).
#' doi: 10.1038/nmeth.4634
#'
#' Adrian Baddeley, Ege Rubak and Rolf Turner:
#' Spatial Point Patterns -- Methodology and Applications with R.
#' CRC Press, 2015. ISBN 978-1-4822-1020-0. 810 pp.
#' \url{http://www.spatstat.org/}
#'
#' @return The SCE object will be returned with additional entries in the rowData:
#'  "r.metric.5" and "spatially.variable_mv.rank"
#'
#' @importFrom spatstat markvario ppp
#' @import future
#'
#' @export
#'
run_markvario <- function(
  sce.obj,
  coords_accessor = "spots",
  exprs_values = "log1p_sctransform",
  genes = NULL,
  r.metric = 5,
  return_SCE = TRUE,
  verbose = TRUE,

    ...
) {

    ## get expression values ------------------------------
    if(is.null(genes)){
      genes <- rownames(sce.obj)
    }else{
      genes <- intersect(genes, rownames(sce.obj))
    }
    if(length(genes) == 0){stop("No genes were left after filtering for those defined via genes.")}

    if(!exprs_values %in% assayNames(sce.obj) ) {stop(paste0("Cannot find assay named", exprs_values))}
    exprs_data <- as.matrix(assay(sce.obj[genes,], exprs_values))

    ## generate mark object ------------------------------------
    pp <- make_pp(sce.obj, coords_accessor )

    if(verbose)print("For each gene, we're calculating the  mark variogram gamma(r), i.e. the dependence between the gene expression values of two cells a distance r apart")
    if (nbrOfWorkers() > 1) {
        chunks <- nbrOfWorkers()
        features <- rownames(x = exprs_data)
        features <- split(
            x = features,
            f = ceiling(x = seq_along(along.with = features) / (length(x = features) / chunks))
        )
        mv <- future_lapply(X = features, FUN = function(x) {
          # test only a couple of expression values per cell at a time
            pp[["marks"]] <- as.data.frame(x = t(x = exprs_data[x, ]))
            markvario(X = pp, normalise = TRUE, ...)
        })
        mv <- unlist(x = mv, recursive = FALSE)
        names(x = mv) <- rownames(x = exprs_data)
    } else {
        pp[["marks"]] <- as.data.frame(x = t(x = exprs_data))
        mv <- markvario(X = pp, normalise = TRUE, ...)
    }

    #return(mv)

    ## calculate the r metric at a given radius (r)
    if(verbose)print("Calculating r metric based on radius")
    mv.out <- ComputeRMetric(mv, r.metric) ## see below for details; taken from Seurat
    mv.out <- mv.out[order(mv.out[, 1]), , drop = FALSE] ## sort
  #  mv.out$spatially.variable_mv <- FALSE
  #  mv.out$spatially.variable_mv[1:(min(nrow(x = mv.out), nfeatures))] <- TRUE
    mv.out$spatially.variable_mv.rank <- 1:nrow(x = mv.out)

    if(return_SCE){
      ## add to SCE
      rd <- rowData(sce.obj)
      rowData(sce.obj) <- mv.out[rownames(rd),] %>% cbind(rd, .)
      return(sce.obj)
    }else{
      return(mv.out)
    }

}

#' Computes the metric at a given r (radius) value
#'
#' @details Taken as is from Seurat 3.1
#'
#' @param mv Results of running markvario
#' @param r.metric r value at which to report the "trans" value of the mark
#'  variogram
#'
#' @return Returns a data.frame with r.metric values
#'
ComputeRMetric <- function(mv, r.metric = 5) {
    r.metric.results <- unlist(x = lapply(
        X = mv,
        FUN = function(x) {
          # returns the edge-corrected difference of exprs values (trans)
          # for two points separated by the distance (r)  that is
          # closest to the user-defined r.metric
            x$trans[which.min(x = abs(x = x$r - r.metric))]
        }
    ))
    r.metric.results <- as.data.frame(x = r.metric.results)
    colnames(r.metric.results) <- paste0("r.metric.", r.metric)
    return(r.metric.results)
}


#' Generate a point pattern object
#'
#' @param sce.obj SCE object with XY coordinates in reducedDim slot and normalized expression values
#' that can be used for testing whether individual genes show expression patterns that correlate with
#' the spatial location of the respective cells
#' @param coords_accessor Indicate which reducedDim slot should be accessed to obtain a matrix of x
#' and y coordinates of the cells
#'
#' @importFrom spatstat ppp
make_pp <- function(sce.obj, coords_accessor ){
  if(!coords_accessor %in% reducedDimNames(sce.obj) ){stop("Cannot find the xy coordinates in the reducedDims slot.")}

  ## extract data
  spatial.coords <- reducedDim(sce.obj, coords_accessor)

  ## generate ppp object
  x.coord = spatial.coords[, 1]
  y.coord = spatial.coords[, 2]

  pp <- ppp(
    x = x.coord,
    y = y.coord,
    xrange = range(x.coord),
    yrange = range(y.coord)
  )

  return(pp)
}

#SpatiallyVariableFeatures.Assay <- function(object, selection.method = "markvariogram", decreasing = TRUE, #...) {
#    CheckDots(...)
#    vf <- SVFInfo(object = object, selection.method = selection.method, status = TRUE)
#    vf <- vf[rownames(x = vf)[which(x = vf[, "variable"][, 1])], ]
#    if (!is.null(x = decreasing)) {
#        vf <- vf[order(x = vf[, "rank"], decreasing = !decreasing), ]
#    }
#    return(rownames(x = vf)[which(x = vf[, "variable"][, 1])])
#}
