% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reading_spots.R
\name{read_spot_coordinates}
\alias{read_spot_coordinates}
\title{Reading the coordinates of the spots/barcodes}
\usage{
read_spot_coordinates(
  spatial_dir,
  slide_name = "Sample.X",
  rm_empty_spots = FALSE
)
}
\arguments{
\item{spatial_dir}{path to the SpaceRanger output directory that contains the
spatial information, e.g. "pN3A_IGO_06000_GL_1_trimmed/outs/spatial"}

\item{slide_name}{Default: "Sample.X"}

\item{rm_empty_spots}{Default: FALSE. If set to TRUE, the returned data.frame
will only contain those barcodes for which SpaceRanger determined overlap
with tissue.}
}
\value{
data.frame with the positions of each barcode
}
\description{
Needs the path to the "spatial/" directory that contains the files
"scalefactors_json.json" and "tissue_positions_list.csv".
}
\examples{
\dontrun{
bc_info <- read_spot_coordinates("Visium/data/pN3A/spatial/",
    slide_name = "pN3A", rm_empty_spots = FALSE)
}

}
