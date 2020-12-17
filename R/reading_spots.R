#' Reading the coordinates of the spots/barcodes
#'
#' @description Needs the path to the "spatial/" directory that contains the files
#' "scalefactors_json.json" and "tissue_positions_list.csv".
#'
#' @param spatial_dir path to the SpaceRanger output directory that contains the
#' spatial information, e.g. "pN3A_IGO_06000_GL_1_trimmed/outs/spatial"
#' @param slide_name Default: "Sample.X"
#' @param rm_empty_spots Default: FALSE. If set to TRUE, the returned data.frame
#' will only contain those barcodes for which SpaceRanger determined overlap
#' with tissue.
#'
#' @return data.frame with the positions of each barcode
#'
#' @examples \dontrun{
#' bc_info <- read_spot_coordinates("Visium/data/pN3A/spatial/",
#'     slide_name = "pN3A", rm_empty_spots = FALSE)
#' }
#'
read_spot_coordinates <- function(spatial_dir, slide_name = "Sample.X", rm_empty_spots = FALSE){

    ## add '/' to the end of the path if it's missing
    spatial_dir <- fix_path(spatial_dir)

    sf_path <- paste0(spatial_dir, "scalefactors_json.json")
    pos_path <- paste0(spatial_dir, "tissue_positions_list.csv")

    if (!file.exists(sf_path)) stop(paste("Sorry,", sf_path, "does not exist"))
    if (!file.exists(pos_path)) stop(paste("Sorry,", pos_path, "does not exist"))

    ## scale factors
    scfc <- rjson::fromJSON(file = sf_path)

    ## read spot positions with cell barcodes
    bcs <- read.csv(pos_path,
        col.names=c("barcode","tissue","row","col","imagerow","imagecol"),
        header = FALSE)

    if(rm_empty_spots){
        ## remove all entries for spots that didn't overlap with tissue
        bcs <- bcs[bcs$tissue != 0, ]
    }

    ## adjust the values
    bcs$y.lowres <- bcs$imagerow * scfc$tissue_lowres_scalef
    bcs$x.lowres <- bcs$imagecol * scfc$tissue_lowres_scalef
    bcs$y.hires <- bcs$imagerow * scfc$tissue_hires_scalef
    bcs$x.hires <- bcs$imagecol * scfc$tissue_hires_scalef

    ## final formatting tweaks
    bcs$tissue <- factor(bcs$tissue)
    bcs$slide <- slide_name
    bcs$barcode <- as.character(bcs$barcode)

    ## make rownames in the same fashion as read10X.abc makes the cellnames
    cell.names <- gsub("-.*$", "", bcs$barcode)
    cell.names <- paste(slide_name, cell.names, sep = "_")
    rownames(bcs) <- cell.names

    return(bcs)
}


#' Loading PNG information for plotting with ggplot2
#'
#' @description This function loads the PNG image of the Visium slide and returns
#' a tibble that contains the corresponding information that can be used to plot
#' the image within the ggplot2 framework (see \code{\link{geom_visium}}).
#'
#' @param image_path path to the image to be loaded and modified
#' @param slide_name Default: "Sample.X"
#' @param which_res either one of "lowres" or "hires"
#' @seealso \code{\link{geom_visium}}
#'
#' @examples \dontrun{
#' img1 <- reading_image("data/pN3A/spatial/", slide_name = "pN3A",which_res = "lowres")
#' bc_info <- read_spot_coordinates("data/pN3A/spatial/", slide_name = "pN3A", rm_empty_spots = FALSE)
#'
#' }
#'
#' @return A tibble with the sample name and the PNG information (grob)
#' @import grid
#' @import tibble
reading_image <- function(image_path, which_res = c("lowres","hires"), slide_name = "Sample.X"){

    image_path <- fix_path(image_path)
    image_path <- paste0(image_path, "tissue_", which_res, "_image.png")
    if (!file.exists(image_path)) stop(paste("Sorry,", image_path, "does not exist."))

    ## load image
    img <- readbitmap::read.bitmap(image_path)
    ## determine height and width
    img_measures <- data.frame(
        height = nrow(img),
        width = ncol(img),
        slide = slide_name
        )

    ## convert to grob for compatibility with ggplot2
    #library(grid)
    img.grob <- rasterGrob(img, width = unit(1,"npc"), height = unit(1,"npc"))

    img.tib <- tibble(
        slide = factor(slide_name),
        grob = list(img.grob),
        height = img_measures$height,
        width = img_measures$width
        )

    return(img.tib)
}


#' Modify the image to improve plotting
#'
#' @description This function allows you to read in an image and change key parameters
#' such as brightness, saturation and hue via \code{\link{magick::image_modulate}}.
#'
#' @details Requires the magick R package. To determine the best settings, use
#' the function with testing_mods = TRUE, which will simply print the results of
#' the modification. Once you're satisfied with the outcome, set testing_mods =
#' FALSE, which will result in a tibble that's equivalent to the one generated
#' by \code{\link{reading_image}}.
#'
#' @param image_path path to the image to be loaded and modified
#' @param slide_name Default: "Sample.X"
#' @param which_res either one of "lowres" or "hires"
#' @param testing_mods Default: TRUE. To determine the best settings, use
#' the function with testing_mods = TRUE, which will simply print the results of
#' the modification. Once you're satisfied with the outcome, set testing_mods =
#' FALSE
#' @param brightness modulation of brightness as percentage of the current value
#'  (100 for no change)
#' @param saturation modulation of saturation as percentage of the current value
#'  (100 for no change)
#' @param hue modulation of hue is an absolute rotation of -180 degrees to +180
#' degrees from the current position corresponding to an argument range of 0 to
#' 200 (100 for no change)
#' @param ... additional parameters for \code{\link{magick::image_modulate}}
#'
#' @return If testing_mods is set to FALSE, then the function will return a
#' tibble that includes the modified glob of the original image.
#'
#' @seealso \code{\link{magick::image_modulate}}, \code{\link{reading_image}}
#'
#' @import grid
#' @import magick
#'
#' @export
reading_image_with_mods <- function(image_path, slide_name = "Sample.X",
    which_res = c("lowres","hires"), testing_mods = TRUE, brightness = 450,
    saturation = 120, hue = 90, ...){

    image_path <- fix_path(image_path)
    image_path <- paste0(image_path, "tissue_", which_res, "_image.png")
    if (!file.exists(image_path)){stop(paste("Sorry,", image_path, "does not exist."))}

    ori <- image_read(image_path)
    mod <- image_modulate(ori,
        brightness = brightness,
        saturation = saturation,
        hue = hue, ...)

    if(testing_mods == TRUE){
        print(mod)
        return()
    }else{
        ## extract grob for compatibility with ggplot2
        mod.grob <- rasterGrob(mod, width = unit(1,"npc"), height = unit(1,"npc"))

        img.tib <- tibble(
            slide = factor(slide_name),
            grob = list(mod.grob),
            height = dim(ori[[1]])[[3]],
            width = dim(ori[[1]])[[2]]
        )
        return(img.tib)
    }
}


#' Helper function for correcting file paths
#'
#' @description Ensures that every file path string has a "/" at the end.
#' @return Path with "/" at the end (string)
fix_path <- function(path){
    if (!grepl("\\/$", path)) {
        out <- paste(path, "/", sep = "")
    }else{
        out <- path
    }
    return(out)
}
