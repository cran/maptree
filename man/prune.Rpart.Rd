\name{prune.Rpart}
\alias{prune.Rpart}
\title{Prunes an Rpart Classification or Regression Tree}
\description{
  Reduces a prediction tree produced by \code{\link{rpart}} to a 
  smaller tree by specifying either a cost-complexity parameter,
  or a number of nodes to which to prune.
}
\usage{
  prune.Rpart (tree, cp=NULL, best=NULL)
}
\arguments{
  \item{tree}{object of class \code{rpart}.}
  \item{cp}{cost-complexity parameter.}
  \item{best}{number of nodes to which to prune.}

  If both \code{cp} and \code{best} are not \code{NULL}, then
  \code{cp} is used.
}
\details{
  A minor enhancement of the existing \code{\link{prune.rpart}} to
  incorporate the parameter \code{best} as it is used in  
  \code{\link{prune.tree}}.  See example.
}
\value{Pruned tree object of class \code{rpart}.}
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link{rpart}}, 
  \code{\link{prune.rpart}}
}
\examples{
  library (rpart)
  data (oregon.env.vars, oregon.border, oregon.grid)

  draw.tree (prune.Rpart (rpart (oregon.env.vars), best=6), 
    nodeinfo=TRUE, units="species", cases="cells")

  group <- group.tree (prune.Rpart (rpart (oregon.env.vars), best=6))
  names(group) <- row.names(oregon.env.vars)
  map.groups (oregon.grid, group)
  lines (oregon.border)
  map.key (0.05, 0.65, lables=as.character(seq(6)), 
    size=4, new=FALSE, sep=1, pch=19, head="node")
}
\keyword{manip}
\keyword{cluster}
