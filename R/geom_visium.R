#' Function for plotting tissue images
#'
#' @description Based on 10X Genomics code (https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/rkit)
#'
#' @examples \dontrun{
#' img1  <- reading_image(paste0(data_dir, sample, "/spatial/"),
#'     sample_name = sample, which_res = "lowres")
#' bc_info <- read_spot_coordinates(paste0(data_dir, sample, "/spatial/"),
#'     sample_name = sample, rm_empty_spots = FALSE)
#'
#' ## plot
#' ggplot(bc_info, aes(x = x.lowres, y = y.lowres)) +
#'     geom_visium(data = img1, aes(grob = grob), x = 0.5, y = 0.5) +
#'     geom_point(shape = 21, aes(colour = tissue), size = 1.75, stroke = 0.5)
#' }
#' @import ggplot2
geom_visium <-  function(mapping = NULL,
    data = NULL,
    stat = "identity",
    position = "identity",
    na.rm = FALSE,
    show.legend = NA,
    inherit.aes = FALSE,
    ...) {

    GeomCustom <- ggproto(
        "GeomCustom",
        Geom,
        setup_data = function(self, data, params) {
            data <- ggproto_parent(Geom, self)$setup_data(data, params)
            data
        },

        draw_group = function(data, panel_scales, coord) {
            vp <- grid::viewport(x=data$x, y=data$y)
            g <- grid::editGrob(data$grob[[1]], vp=vp)
            ggplot2:::ggname("geom_visium", g)
        },

        required_aes = c("grob","x","y")

    )

    layer(
        geom = GeomCustom,
        mapping = mapping,
        data = data,
        stat = stat,
        position = position,
        show.legend = show.legend,
        inherit.aes = inherit.aes,
        params = list(na.rm = na.rm, ...)
    )
}
