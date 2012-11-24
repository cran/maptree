\name{draw.clust}
\alias{draw.clust}
\title{Graph a Hierarchical Cluster Tree}
\description{
  Graph a hierarchical cluster tree of class \code{twins}
  or \code{hclust} using colored symbols at observations.
}
\usage{
  draw.clust (cluster, data=NULL, cex=par("cex"), pch=par("pch"), size=2.5*cex, 
      col=NULL, nodeinfo=FALSE, cases="obs", new=TRUE)
}
\arguments{
  \item{cluster}{object of class \code{hclust} or \code{twins}.}
  \item{data}{clustered dataset for hclust application.}
  \item{cex}{size of text, par parameter.}
  \item{pch}{shape of symbol at leaves, par parameter.}
  \item{size}{size in cex units of symbol at leaves.}
  \item{col}{vector of colors from \code{\link{hsv}}, \code{\link{rgb}},
             etc, or if \code{NULL}, then use \code{\link{rainbow}}.}
  \item{nodeinfo}{if \code{TRUE}, add a line at each node with number
                  of observations included in each leaf.}
  \item{cases}{label for type of observations.}
  \item{new}{if \code{TRUE}, call \code{\link{plot.new}}.}
}
\details{
  An alternative to \code{\link{pltree}} and \code{\link{plot.hclust}}.
}
\value{The vector of colors supplied or generated.}
\author{Denis White}
\seealso{
  \code{\link[cluster]{agnes}}, 
  \code{\link[cluster]{diana}}, 
  \code{\link{hclust}}, 
  \code{\link{draw.tree}}, 
  \code{\link{map.groups}}
}
\examples{
  library (cluster)
  data (oregon.bird.dist)

  draw.clust (clip.clust (agnes (oregon.bird.dist), k=6))
}
\keyword{hplot}
\keyword{cluster}
