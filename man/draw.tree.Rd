\name{draw.tree}
\alias{draw.tree}
\title{Graph a Classification or Regression Tree}
\description{
  Graph a classification or regression tree with a hierarchical
  tree diagram, optionally including shaded square symbols at   
  leaves and additional info at intermediate nodes.
}
\usage{
  draw.tree (tree, ps=par("ps"), size=10, col=NULL, 
      nodeinfo=FALSE, units="units", cases="obs", digits=0)
}
\arguments{
  \item{tree}{object of class \code{"tree"}.}
  \item{ps}{pointsize of text, par parameter.}
  \item{size}{if size=0, draw terminal symbol at leaves else a
              colored square with sides of length size in mm.}
  \item{col}{vector of colors from \code{\link{hsv}}, \code{\link{rgb}} 
             etc, or if NULL, then generated as equally spaced hues, or
             if "gray", then generated as gray tones.}
  \item{nodeinfo}{if TRUE, add a line at each node with mean value
                  of response, number of observations, and percent
                  deviance explained (or classified correct).}
  \item{units}{label for units of mean value of response, if
               regression tree.}
  \item{cases}{label for type of observations.}
  \item{digits}{parameter for rounding mean value of response, if
                regression tree.}
}
\details{
  Styled like S \code{post.tree}.  As in \code{plot.tree(,,type="u")}, each  
  level has constant depth.  With nodeinfo=TRUE, deviance explained
  at each node, or classification rate, is shown in text format.
}
\value{The vector of colors supplied or generated.}
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link{draw.clust}}, 
  \code{\link{map.groups}}
}
\examples{
  library (tree)
  data (hex.env)

  draw.tree (prune.tree (tree (hex.env), best=6), 
      nodeinfo=TRUE, units="species", cases="cells")
}
\keyword{hplot}
\keyword{tree}
