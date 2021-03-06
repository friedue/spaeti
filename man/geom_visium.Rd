% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geom_visium.R
\name{geom_visium}
\alias{geom_visium}
\title{Function for plotting tissue images}
\usage{
geom_visium(
  mapping = NULL,
  data = NULL,
  stat = "identity",
  position = "identity",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = FALSE,
  ...
)
}
\description{
Based on 10X Genomics code (https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/rkit)
}
\examples{
\dontrun{
img1  <- reading_image(paste0(data_dir, sample, "/spatial/"),
    sample_name = sample, which_res = "lowres")
bc_info <- read_spot_coordinates(paste0(data_dir, sample, "/spatial/"),
    sample_name = sample, rm_empty_spots = FALSE)

## plot
ggplot(bc_info, aes(x = x.lowres, y = y.lowres)) +
    geom_visium(data = img1, aes(grob = grob), x = 0.5, y = 0.5) +
    geom_point(shape = 21, aes(colour = tissue), size = 1.75, stroke = 0.5)
}
}
