\name{clip.rpart}
\alias{clip.rpart}
\title{Prunes an Rpart Classification or Regression Tree}
\description{
  Reduces a prediction tree produced by \code{\link{rpart}} to a
  smaller tree by specifying either a cost-complexity parameter,
  or a number of nodes to which to prune.
}
\usage{
  clip.rpart (tree, cp=NULL, best=NULL)
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
  incorporate the parameter \code{best} as it is used in the 
  (now defunct)  \code{prune.tree} function in the old \pkg{tree}
  package.  See example.
}
\value{Pruned tree object of class \code{rpart}.}
\author{Denis White}
\seealso{
  \code{\link[rpart]{rpart}}, 
  \code{\link[rpart]{prune.rpart}}
}
\examples{
  library (rpart)
  data (oregon.env.vars, oregon.border, oregon.grid)

  draw.tree (clip.rpart (rpart (oregon.env.vars), best=7), 
    nodeinfo=TRUE, units="species", cases="cells", digits=0)

  group <- group.tree (clip.rpart (rpart (oregon.env.vars), best=7))
  names(group) <- row.names(oregon.env.vars)
  map.groups (oregon.grid, group)
  lines (oregon.border)
  map.key (0.05, 0.65, labels=as.character(seq(6)), 
    size=1, new=FALSE, sep=0.5, pch=19, head="node")
}
\keyword{manip}
\keyword{cluster}
