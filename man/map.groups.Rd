\name{map.groups}
\alias{map.groups}
\title{Map Groups of Observations}
\description{
  Draws maps of groups of observations created by clustering,
  classification or regression trees, or some other type of 
  classification.
}
\usage{
  map.groups (pts, size=5, col=NULL, gon=10, lns=NULL, 
      lcol=1, lwd=par("lwd"))
}
\arguments{
  \item{pts}{matrix or data frame with components \code{"x"}, 
             \code{"y"}, and \code{"group"}.}
  \item{size}{size in mm of shaded square symbol.}
  \item{col}{vector of colors from \code{\link{hsv}}, \code{\link{rgb}} 
             etc, or if NULL, then generated as equally spaced hues, or
             if "gray", then generated as gray tones.}
  \item{gon}{number of sides for polygon symbol (>8 => circle).}
  \item{lns}{matrix of S format \code{\link{lines}} coordinates, or 
             a list of several matrices.}
  \item{lcol}{index for line color, or a vector of several 
              indeces recycled if necessary.}
  \item{lwd}{line width, or a vector of several widths recycled
             if necessary.}
}
\details{
  Uses \code{\link{ngon}} to draw shaded polygon symbols for
  each observation in pts.
}
\value{The vector of colors supplied or generated.}
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link{ngon}}, 
  \code{\link{lines}}, 
  \code{\link{group.clust}}, 
  \code{\link{group.tree}}, 
  \code{\link{map.key}}
}
\examples{
  data (spp.names, hex.env, hex.spp, oregon)

  # range map for American Avocet
  spp <- match ("American avocet", spp.names[["common.name"]])
  map.groups (cbind (hex.env, group=hex.spp[,spp]+1), 
      col="gray", gon=6, lns=oregon)

  # distribution of January temperatures
  cuts <- quantile (hex.env[["jan.temp"]], probs=seq(0,1,1/5))
  group <- cut (hex.env[["jan.temp"]], cuts, labels=FALSE, 
      include.lowest=TRUE)
  map.groups (cbind (hex.env, group=group), 
      col="gray", gon=6, lns=oregon)

}
\keyword{hplot}
\keyword{cluster}
\keyword{tree}
