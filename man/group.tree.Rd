\name{group.tree}
\alias{group.tree}
\title{Observation Groups for Classification or Regression Tree}
\description{
  Alternative to \code{tree[["where"]]} that orders groups from left
  to right in draw order.
}
\usage{
  group.tree (tree)
}
\arguments{
  \item{tree}{object of class \code{tree} or \code{rpart}.}
}
\details{
  Normally used with \code{\link{map.groups}}.  See example.
}
\value{Vector of rearranged \code{tree[["where"]]}} 
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link[tree]{tree}}, 
  \code{\link[rpart]{rpart}}, 
  \code{\link{map.groups}}
}
\examples{
  library (tree)
  data (oregon.env.vars, oregon.grid)

  group <- group.tree (prune.tree (tree (oregon.env.vars), best=6))
  names(group) <- row.names(oregon.env.vars)
  map.groups (oregon.grid, group=group)
}
\keyword{manip}
\keyword{tree}
