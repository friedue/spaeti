#' Assign a spot ID based on whether it overlaps with a selection of polygons
#'
#'
#' @param sce.object SingleCellExperiment Object
#' @param subset_by Default: "Sample"
#' @param coord_list Must contain one element per entry of "subset_by"
#' @param id_name Default: "sec"
#' @param plot_overlap Boolean. Default: FALSE. If set to TRUE, the function
#' will not change the sce.obj, but will simply plot the fit for visual inspection
#'
#' @return SCE object with additional column in the colData that indicates the
#' name of the polygon that a given spot overlapped with
#'
#' @seealso Wrapper around \code{\link{determine_overlap}},
#'  \code{\link{create_polygons_from_xy}}
#'
#' @examples \dontrun{
#' sec.coords <-list()
#' smp <- "pN3A"
#' sce.norm[, sce.norm$Sample == smp] %>% reducedDim(., "spots") %>% plot
#' sec.coords[[smp]]$secI <- locator(15)
#' sec.coords[[smp]]$secII <- locator(15)
#' }
#'
#' @export
#'
assign_ID_based_on_overlap <- function(sce.obj, subset_by = "Sample", coord_list, id_name = "sec", plot_overlap = FALSE){

    if(!subset_by %in% names(colData(sce.obj))){stop(paste(subset_by,"is not part of the SCE's colData"))}
    sets <- unique(colData(sce.obj)[, subset_by, drop = TRUE])
    if(plot_overlap){
        for(w in sets){
            cl <- coord_list[[w]]
            if(!is.null(cl)){
                determine_overlap(sce.obj[, colData(sce.obj)[, subset_by] == w],
                    cl, plot_overlap)
            }
        }
        return()
    }else{
        ids <- lapply(sets, function(w){
            cl <- coord_list[[w]]
            if(!is.null(cl)){
                determine_overlap(sce.obj[, colData(sce.obj)[, subset_by] == w],
                    cl, plot_overlap)
            }
        }) %>% unlist

        colData(sce.obj)[names(ids), id_name] <- ids
        return(sce.obj)
    }
}

#' Create polygon objects based on x and y coordinates
#'
#'
#' @param coord_list list of x and y coordinates per region of interest, where each
#' element of the list coorresponds to a region, containing x and y coordinates
#' as another list or data.frame. It should be a named list.
#' @return SpatialPolygons object
#'
create_polygons_from_xy <- function(coord_list){

    if(is.null(names(coord_list))){stop("The coord_list object must have names.")}

    ## create a list of polygons with the names from the coord_list
    ## as IDs
    poly.l <- lapply(seq_along(coord_list), function(i){
        id <- names(coord_list[i])
        X <- coord_list[[i]]
        if(!all(c("x","y") %in% names(X)) ){stop(paste("Element", id, "does not have x and y elements."))}
        pl.ind <- sp::Polygon(cbind(X$x, X$y))
        return(sp::Polygons(list(pl.ind), ID = id))
    }) %>% sp::SpatialPolygons()
}


#' Determine the overlap between spots from an SCE and user-defined polygons
#'
#' @param sce.obj SCE object containing the cells that should be assigned an
#' ID based on which polygon defined in coord_list they overlap with
#' @param coord_list list of x and y coordinates per region of interest, where each
#' element of the list corresponds to a region, containing x and y coordinates
#' as another list or data.frame. It should be a named list.
#' @param plot_overlap Boolean. Default: FALSE. If set to TRUE, the function
#' will not change the sce.obj, but will simply plot the fit for visual inspection
#'
#' @return if plot_overlap = FALSE, this will return an named vector of
#' points (cells) overlapping with the polygons
#'
#' @import sp
determine_overlap <- function(sce.obj, coord_list, plot_overlap = FALSE){

    polys <- create_polygons_from_xy(coord_list)
    pnts <- SpatialPoints(reducedDim(sce.obj, "spots"))
    ols <- over(pnts, polys)
    if(plot_overlap){
        plot(polys)
        plot(pnts, add = TRUE, pch = 16, col = ols)
        legend("topleft", names(polys), col = unique(ols), pch = 19)
        return()
    }else{
     ols.nm <- names(polys)[ols]
     names(ols.nm) <- names(ols)
     return(ols.nm)
    }
}


