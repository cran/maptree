\name{group.clust}
\alias{group.clust}
\title{Observation Groups for a Hierarchical Cluster Tree}
\description{
  Alternative to \code{\link{cutree}} that orders pruned groups from left to
  right in draw order, and adds group membership of
  observations to pts matrix if provided.
}
\usage{
  group.clust (cluster, k=NULL, h=NULL, pts=NULL)
}
\arguments{
  \item{cluster}{object of class \code{"hclust"} or \code{"twins"}.}
  \item{k}{desired number of groups.}
  \item{h}{height at which to prune for grouping.}
  \item{pts}{vector, matrix, or data frame of length or nrow equal to
             \code{length(cluster[["order"]])} to which group vector is 
             joined with cbind.}

  At least one of k or h must be specified; k takes precedence
  if both are given.
}
\details{
  Normally used with \code{\link{map.groups}}.  See example.
}
\value{Vector of pruned cluster membership or the cbind of pts and
       this vector if pts != NULL.}
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link{hclust}}, 
  \code{\link{twins.object}}, 
  \code{\link{cutree}}, 
  \code{\link{map.groups}}
}
\examples{
  library (mva)
  data (hex.spp, hex.env, oregon)

  map.groups (group.clust (hclust (dist (hex.spp)), k=6,
      pts=hex.env), lns=oregon)
}
\keyword{manip}
\keyword{cluster}
