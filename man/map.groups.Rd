\name{map.groups}
\alias{map.groups}
\title{Map Groups of Observations}
\description{
  Draws maps of groups of observations created by clustering,
  classification or regression trees, or some other type of
  classification.
}
\usage{
  map.groups (pts, group, pch=19, size=2, col=NULL, border=NULL)
}
\arguments{
  \item{pts}{matrix or data frame with components \code{"x"},
             and \code{"y"} for each observation (see details).}
  \item{group}{vector of integer class numbers corresponding to
               \code{pts} (see details), and indexing colors
               in \code{col}.}
  \item{pch}{symbol number from \code{par("pch")} if < 100, otherwise
             parameter \code{n} for ngon.}
  \item{size}{size in cex units of shaded point symbol.}
  \item{col}{vector of fill colors from \code{\link{hsv}}, \code{\link{rgb}},
             etc, or if \code{NULL}, then use \code{\link{rainbow}}, or
             if \code{gray} or \code{grey}, then use \code{\link{gray}}.}
  \item{border}{vector of border colors from \code{\link{hsv}},
             \code{\link{rgb}},
             etc, or if \code{NULL}, then use \code{\link{rainbow}}, or
             if \code{gray} or \code{grey}, then use \code{\link{gray}}.}
}
\details{
  If the number of rows of \code{pts} is not equal to the length
  of \code{group}, then (1) \code{pts} are assumed to represent
  polygons and \code{\link{polygon}} is used, and (2) the fill and
  border colors indexed by \code{group} are matched to the polygons
  in \code{pts} through \code{names(group)} and \code{pts$x[is.na(pts$y)]}.
  Otherwise, \code{group} is assumed to parallel \code{pts}, and,
  if \code{pch < 100}, then \code{\link{points}} is used, otherwise
  \code{\link{ngon}}, to draw shaded polygon symbols for each
  observation in pts.

  This graphic function maintains the aspect ratio of the input coordinates
  in the output graphic by setting \code{par(pin=)}, thus the map does not
  scale dynamically in an on-screen window.
}
\value{The vector of fill colors supplied or generated.}
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link{ngon}}, 
  \code{\link{polygon}}, 
  \code{\link{group.clust}}, 
  \code{\link{group.tree}}, 
  \code{\link{map.key}}
}
\examples{
  data (oregon.bird.names, oregon.env.vars, oregon.bird.dist)
  data (oregon.border, oregon.grid)

  # range map for American Avocet
  spp <- match ("American avocet", oregon.bird.names[["common.name"]])
  group <- oregon.bird.dist[,spp] + 1
  names(group) <- row.names(oregon.bird.dist)
  map.groups (oregon.grid, group=group, col="gray")
  lines (oregon.border)

  # distribution of January temperatures
  cuts <- quantile (oregon.env.vars[["jan.temp"]], probs=seq(0,1,1/5))
  group <- cut (oregon.env.vars[["jan.temp"]], cuts, labels=FALSE, 
    include.lowest=TRUE)
  names(group) <- row.names(oregon.env.vars)
  map.groups (oregon.grid, group=group, col="gray")
  lines (oregon.border)

  # January temperatures using point symbols rather than polygons
  map.groups (oregon.env.vars, group, col="gray")
  lines (oregon.border)
}
\keyword{hplot}
\keyword{cluster}
\keyword{tree}
