\name{draw.tree}
\alias{draw.tree}
\title{Graph a Classification or Regression Tree}
\description{
  Graph a classification or regression tree with a hierarchical
  tree diagram, optionally including shaded square symbols at
  leaves and additional info at intermediate nodes.
}
\usage{
  draw.tree (tree, cex=par("cex"), size=2, col=NULL, 
      nodeinfo=FALSE, units="units", cases="obs", digits=0)
}
\arguments{
  \item{tree}{object of class \code{tree}.}
  \item{cex}{size of text, par parameter.}
  \item{size}{if \code{size=0}, draw terminal symbol at leaves else a
              colored square with sides of length size in cex units.}
  \item{col}{vector of colors from \code{\link{hsv}}, \code{\link{rgb}},
             etc, or if \code{NULL}, then use \code{\link{rainbow}}, or
             if \code{gray} or \code{grey}, then use \code{\link{gray}}.}
  \item{nodeinfo}{if \code{TRUE}, add a line at each node with mean value
                  of response, number of observations, and percent
                  deviance explained (or classified correct).}
  \item{units}{label for units of mean value of response, if
               regression tree.}
  \item{cases}{label for type of observations.}
  \item{digits}{parameter for rounding mean value of response, if
                regression tree.}
}
\details{
  Styled like \code{post.tree}.  As in \code{plot.tree(,,type="u")}, each
  level has constant depth.  Specifying \code{nodeinfo=TRUE}, shows the
  deviance explained or the classification rate at each node.
}
\value{The vector of colors supplied or generated.}
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link[tree]{tree}}, 
  \code{\link[tree]{prune.tree}}, 
  \code{\link{draw.clust}}, 
  \code{\link{map.groups}}
}
\examples{
  library (tree)
  data (oregon.env.vars)

  draw.tree (prune.tree (tree (oregon.env.vars), best=6), 
      nodeinfo=TRUE, units="species", cases="cells")
}
\keyword{hplot}
\keyword{tree}
