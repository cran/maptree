\name{draw.clust}
\alias{draw.clust}
\title{Graph a Hierarchical Cluster Tree}
\description{
  Graph a hierarchical cluster tree of class \code{"twins"} 
  or \code{"hclust"} using shaded square symbols at observations.
}
\usage{
  draw.clust (cluster, ps=par("ps"), size=10, col=NULL)
}
\arguments{
  \item{cluster}{object of class \code{"hclust"} or \code{"twins"}.}
  \item{ps}{pointsize of text, par parameter.}
  \item{size}{size in mm of shaded square symbol.}
  \item{col}{vector of colors from \code{\link{hsv}}, \code{\link{rgb}} 
             etc, or if NULL, then generated as equally spaced hues, or
             if "gray", then generated as gray tones.}
}
\details{
  Styled like plot.hclust.
}
\value{The vector of colors supplied or generated.}
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link{draw.tree}}, 
  \code{\link{map.groups}}
}
\examples{
  library (mva)
  data (hex.spp)

  draw.clust (prune.clust (hclust (dist (hex.spp)), k=6))
}
\keyword{hplot}
\keyword{cluster}
