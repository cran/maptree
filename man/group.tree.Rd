\name{group.tree}
\alias{group.tree}
\title{Observation Groups for Classification or Regression Tree}
\description{
  Alternative to \code{tree[["where"]]} that orders groups from left
  to right in draw order, and adds group membership of
  observations to pts matrix if provided.
}
\usage{
  group.tree (tree, pts=NULL)
}
\arguments{
  \item{tree}{object of class \code{"tree"}.}
  \item{pts}{vector, matrix, or data frame of length or nrow equal to 
             \code{length(tree[["where"]])} to which group vector 
             is joined with cbind.}
}
\details{
  Normally used with \code{\link{map.groups}}.  See example.
}
\value{Vector of rearranged \code{tree[["where"]]} or the cbind of pts and
       this vector if pts != NULL.}
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link{tree}}, 
  \code{\link{map.groups}}
}
\examples{
  library (tree)
  data (hex.env, oregon)

  map.groups (group.tree (prune.tree (tree (hex.env),
      best=6), pts=hex.env), lns=oregon)
}
\keyword{manip}
\keyword{tree}
