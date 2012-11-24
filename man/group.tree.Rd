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
  \item{tree}{object of class \code{rpart} or \code{tree}.}
}
\details{
  Normally used with \code{\link{map.groups}}.  See example.
}
\value{Vector of rearranged \code{tree[["where"]]}} 
\author{Denis White}
\seealso{
  \code{\link[rpart]{rpart}}, 
  \code{\link{map.groups}}
}
\examples{
  library (rpart)
  data (oregon.env.vars, oregon.grid)

  group <- group.tree (clip.rpart (rpart (oregon.env.vars), best=7))
  names(group) <- row.names(oregon.env.vars)
  map.groups (oregon.grid, group=group)
}
\keyword{manip}
\keyword{tree}
