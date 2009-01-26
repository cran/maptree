\name{twins.to.hclust}
\alias{twins.to.hclust}
\title{Converts agnes or diana object to hclust object}
\description{
  Alternative to \code{\link{as.hclust}} that retains cluster data.
}
\usage{
  twins.to.hclust (cluster)
}
\arguments{
  \item{cluster}{object of class \code{twins}.}

}
\details{
  Used internally in with \code{\link{clip.clust}} and 
  \code{\link{draw.clust}}.
}
\value{hclust object} 
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link{hclust}}, 
  \code{\link[cluster]{twins.object}}
}
\keyword{manip}
\keyword{cluster}
