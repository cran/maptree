\name{clip.clust}
\alias{clip.clust}
\title{Prunes a Hierarchical Cluster Tree}
\description{
  Reduces a hierarchical cluster tree to a smaller tree either
  by pruning until a given number of observation groups remain,
  or by pruning tree splits below a given height.
}
\usage{
  clip.clust (cluster, data=NULL, k=NULL, h=NULL)
}
\arguments{
  \item{cluster}{object of class \code{hclust} or \code{twins}.}
  \item{data}{clustered dataset for hclust application.}
  \item{k}{desired number of groups.}
  \item{h}{height at which to prune for grouping.}

  At least one of \code{k} or \code{h} must be specified; \code{k} takes
  precedence if both are given.
}
\details{
  Used with \code{\link{draw.clust}}.  See example.
}
\value{Pruned cluster object of class \code{hclust}.}
\author{Denis White }
\seealso{
  \code{\link{hclust}}, 
  \code{\link[cluster]{twins.object}}, 
  \code{\link{cutree}}, 
  \code{\link{draw.clust}}
}
\examples{
  library (cluster)
  data (oregon.bird.dist)

  draw.clust (clip.clust (agnes (oregon.bird.dist), k=6))
}
\keyword{manip}
\keyword{cluster}
