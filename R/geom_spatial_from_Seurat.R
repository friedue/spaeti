## For inspiration on how to implement cropping/zooming
#========================================================
# For plotting the tissue image
#' @importFrom ggplot2 ggproto Geom aes ggproto_parent alpha draw_key_point
#' @importFrom grid unit gpar editGrob pointsGrob viewport gTree addGrob grobName
#' @export
#'
# GeomSpatial <- ggproto(
#     "GeomSpatial",
#     Geom,
#     required_aes = c("x", "y"),
#     extra_params = c("na.rm", "image", "crop"),
#     default_aes = aes(
#         shape = 21,
#         colour = "black",
#         point.size.factor = 1.0,
#         fill = NA,
#         alpha = NA,
#         stroke = 0.25
#     ),
#     setup_data = function(self, data, params) {
#         data <- ggproto_parent(Geom, self)$setup_data(data, params)
#         # We need to flip the image as the Y coordinates are reversed
#         data$y = max(data$y) - data$y + min(data$y)
#         data
#     },
#     draw_key = draw_key_point,
#     draw_panel = function(data, panel_scales, coord, image, crop) {
#         # This should be in native units, where
#         # Locations and sizes are relative to the x- and yscales for the current viewport.
#         if (!crop) {
#             y.transform <- c(0, nrow(x = image)) - panel_scales$y.range
#             data$y <- data$y + sum(y.transform)
#             panel_scales$y.range <- c(0, nrow(x = image))
#             panel_scales$x.range <- c(0, ncol(x = image))
#         }
#
#         z = coord$transform(
#             data.frame(x = c(0, ncol(x = image)), y = c(0, nrow(x = image))),
#             panel_scales
#         )
#         # Flip Y axis for image
#         z$y = -rev(z$y) + 1
#         wdth = z$x[2] - z$x[1]
#         hgth = z$y[2] - z$y[1]
#         vp <- grid::viewport(
#             x = unit(x = z$x[1], units = "npc"),
#             y = unit(x = z$y[1], units = "npc"),
#             width = unit(x = wdth, units = "npc"),
#             height = unit(x = hgth, units = "npc"),
#             just = c("left", "bottom")
#         )
#         img.grob <- GetImage(object = image)
#         img <- editGrob(grob = img.grob, vp = vp)
#         spot.size <- slot(object = image, name = "spot.radius")
#         coords <- coord$transform(data, panel_scales)
#         pts <- pointsGrob(
#             x = coords$x,
#             y = coords$y,
#             pch = data$shape,
#             size = unit(spot.size, "npc") * data$point.size.factor,
#             gp = gpar(
#                 col = alpha(coords$colour, coords$alpha),
#                 fill = alpha(coords$fill, coords$alpha),
#                 lwd = coords$stroke)
#         )
#         vp <- viewport()
#         gt <- gTree(vp = vp)
#         gt <- addGrob(gTree = gt, child = img)
#         gt <- addGrob(gTree = gt, child = pts)
#         # Replacement for ggname
#         gt$name <- grobName(grob = gt, prefix = 'geom_spatial')
#         return(gt)
#         # ggplot2:::ggname("geom_spatial", gt)
#     }
# )
