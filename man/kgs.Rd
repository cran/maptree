\name{kgs}
\alias{kgs}
\title{KGS Measure for Pruning Hierarchical Clusters}
\description{
  Computes the Kelley-Gardner-Sutcliffe penalty function for a 
  hierarchical cluster tree.
}
\usage{
  kgs (cluster, diss, maxclust=NULL)
}
\arguments{
  \item{cluster}{object of class \code{hclust} or \code{twins}.}
  \item{diss}{object of class \code{dissimilarity} or \code{dist}.}
  \item{maxclust}{maximum number of clusters for which to compute
                  measure.}
}
\details{
  Kelley et al. (see reference) proposed a method that can help decide 
  to where to prune a hierarchical cluster tree.  At any level of the 
  tree the mean across all clusters of the mean within clusters of the 
  dissimilarity measure is calculated.  After normalizing, the number 
  of clusters is added.  The minimum of this function is the suggested 
  pruning size.

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
  \code{\link[mva]{hclust}}, 
  \code{\link[mva]{dist}},
  \code{\link{prune.clust}},
  \code{\link[combinat]{combn}}
}
\examples{
  library (cluster)
  library (combinat)
  data (oregon.bird.dist)

  a <- agnes (oregon.bird.dist)
  plot (kgs (a, a$diss, maxclust=20))
}
\keyword{manip}
\keyword{cluster}
