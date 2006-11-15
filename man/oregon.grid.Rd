\name{oregon.grid}
\alias{oregon.grid}
\title{Hexagon Grid Cell Polygons covering Oregon, USA}
\description{
 Polygon borders for 389 hexagonal grid cells covering Oregon, USA,
 in \code{\link{polygon}} format.
}
\usage{data (oregon.grid)}
\format{A data frame with 3112 rows and 2 columns (the components
  \code{"x"} and \code{"y"}).}
\details{
 The polygon format used for these grid cell boundaries is a slight
 variation from the standard R/S format.  Each cell polygon is
 described by seven coordinate pairs, the last repeating the first.
 Prior to the first coordinate pair of each cell is a row containing
 NA in the \code{"y"} column and, in the \code{"x"} column, an
 identifier for the cell.  The identifiers are the same as the
 row names in \code{\link{oregon.bird.dist}} and
 \code{\link{oregon.env.vars}}.  See \code{\link{map.groups}} for
 how the linkage is made in mapping.

 These grid cells are extracted from a larger set covering the
 conterminous United States and adjacent parts of Canada and
 Mexico, as described in White et al. (1992).  Only cells with
 at least 50 percent of their area contained within the state of
 Oregon are included.

 The map projection for the coordinates, as well as the point
 coordinates in \code{\link{oregon.env.vars}}, is the Lambert
 Conformal Conic with standard parallels at 33 and 45
 degrees North latitude, with the longitude of the central
 meridian at 120 degrees, 30 minutes West longitude,
 and with the projection origin latitude at 41 degrees,
 45 minutes North latitude.
}
\source{Denis White, \email{white.denis@epa.gov} }
\references{
 White, D., Kimerling, A.J., Overton, W.S.  (1992)  Cartographic
 and geometric components of a global sampling design for
 environmental monitoring, \emph{Cartography and Geographic
 Information Systems}, \bold{19}(1), 5-22.
}
\keyword{datasets}
