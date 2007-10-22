\name{kgs}
\alias{kgs}
\title{KGS Measure for Pruning Hierarchical Clusters}
\description{
  Computes the Kelley-Gardner-Sutcliffe penalty function for a 
  hierarchical cluster tree.
}
\usage{
  kgs (cluster, diss, alpha=1, maxclust=NULL)
}
\arguments{
  \item{cluster}{object of class \code{hclust} or \code{twins}.}
  \item{diss}{object of class \code{dissimilarity} or \code{dist}.}
  \item{alpha}{weight for number of clusters.}
  \item{maxclust}{maximum number of clusters for which to compute
                  measure.}
}
\details{
  Kelley et al. (see reference) proposed a method that can help decide 
  where to prune a hierarchical cluster tree.  At any level of the 
  tree the mean across all clusters of the mean within clusters of the 
  dissimilarity measure is calculated.  After normalizing, the number 
  of clusters times alpha is added.  The minimum of this function  
  corresponds to the suggested pruning size.

  The current implementation has complexity O(n*n*maxclust), thus 
  very slow with large n.  For improvements, at least it should only 
  calculate the spread for clusters that are split at each level, 
  rather than over again for all.

  Uses contributed function \code{\link[combinat]{combn}}.
}
\value{Vector of the penalty function for trees of size 2:maxclust. 
  The names of vector elements are the respective numbers of clusters.}
\references{
  Kelley, L.A., Gardner, S.P., Sutcliffe, M.J.  (1996)  An automated
  approach for clustering an ensemble of NMR-derived protein structures 
  into conformationally-related subfamilies, \emph{Protein Engineering}, 
  \bold{9}, 1063-1065.
}
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link[cluster]{twins.object}}, 
  \code{\link[cluster]{dissimilarity.object}}, 
  \code{\link{hclust}}, 
  \code{\link{dist}},
  \code{\link{clip.clust}},
  \code{\link[combinat]{combn}}
}
\examples{
  library (cluster)
  data (votes.repub)

  a <- agnes (votes.repub, method="ward")
  b <- kgs (a, a$diss, maxclust=20)
  plot (names (b), b, xlab="# clusters", ylab="penalty")
}
\keyword{manip}
\keyword{cluster}
