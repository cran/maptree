\name{map.groups}
\alias{map.groups}
\title{Map Groups of Observations}
\description{
  Draws maps of groups of observations created by clustering,
  classification or regression trees, or some other type of
  classification.
}
\usage{
  map.groups (pts, group, pch=par("pch"), size=2, col=NULL, 
      border=NULL, new=TRUE)
}
\arguments{
  \item{pts}{matrix or data frame with components \code{"x"},
             and \code{"y"} for each observation (see details).}
  \item{group}{vector of integer class numbers corresponding to
               \code{pts} (see details), and indexing colors
               in \code{col}.}
  \item{pch}{symbol number from \code{par("pch")} if < 100, otherwise
             parameter \code{n} for ngon.}
  \item{size}{size in cex units of point symbol.}
  \item{col}{vector of fill colors from \code{\link{hsv}}, \code{\link{rgb}},
             etc, or if \code{NULL}, then use \code{\link{rainbow}}.}
  \item{border}{vector of border colors from \code{\link{hsv}},
             \code{\link{rgb}},
             etc, or if \code{NULL}, then use \code{\link{rainbow}}.}
  \item{new}{if \code{TRUE}, call \code{\link{plot.new}}.}
}
\details{
  If the number of rows of \code{pts} is not equal to the length
  of \code{group}, then (1) \code{pts} are assumed to represent
  polygons and \code{\link{polygon}} is used, (2) the identifiers in 
  \code{group} are matched to the polygons in \code{pts} through 
  \code{names(group)} and \code{pts$x[is.na(pts$y)]}, and (3) these
  identifiers are mapped to dense integers to reference colours.
  Otherwise, \code{group} is assumed to parallel \code{pts}, and,
  if \code{pch < 100}, then \code{\link{points}} is used, otherwise
  \code{\link{ngon}}, to draw shaded polygon symbols for each
  observation in pts.
}
\value{The vector of fill colors supplied or generated.}
\author{Denis White}
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
  kol <- gray (seq(0.8,0.2,length.out=length (table (group))))
  map.groups (oregon.grid, group=group, col=kol)
  lines (oregon.border)

  # distribution of January temperatures
  cuts <- quantile (oregon.env.vars[["jan.temp"]], probs=seq(0,1,1/5))
  group <- cut (oregon.env.vars[["jan.temp"]], cuts, labels=FALSE, 
    include.lowest=TRUE)
  names(group) <- row.names(oregon.env.vars)
  kol <- gray (seq(0.8,0.2,length.out=length (table (group))))
  map.groups (oregon.grid, group=group, col=kol)
  lines (oregon.border)

  # January temperatures using point symbols rather than polygons
  map.groups (oregon.env.vars, group, col=kol, pch=19)
  lines (oregon.border)
}
\keyword{hplot}
\keyword{cluster}
\keyword{tree}
