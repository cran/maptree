\name{group.clust}
\alias{group.clust}
\title{Observation Groups for a Hierarchical Cluster Tree}
\description{
  Alternative to \code{\link{cutree}} that orders pruned groups from
  left to right in draw order.
}
\usage{
  group.clust (cluster, k=NULL, h=NULL)
}
\arguments{
  \item{cluster}{object of class \code{hclust} or \code{twins}.}
  \item{k}{desired number of groups.}
  \item{h}{height at which to prune for grouping.}

  At least one of \code{k} or \code{h} must be specified; \code{k} takes
  precedence if both are given.
}
\details{
  Normally used with \code{\link{map.groups}}.  See example.
}
\value{Vector of pruned cluster membership} 
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link{hclust}}, 
  \code{\link[cluster]{twins.object}}, 
  \code{\link{cutree}}, 
  \code{\link{map.groups}}
}
\examples{
  data (oregon.bird.dist, oregon.grid)

  group <- group.clust (hclust (dist (oregon.bird.dist)), k=6)
  names(group) <- row.names(oregon.bird.dist)
  map.groups (oregon.grid, group)
}
\keyword{manip}
\keyword{cluster}
