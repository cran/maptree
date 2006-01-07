\name{draw.tree}
\alias{draw.tree}
\title{Graph a Classification or Regression Tree}
\description{
  Graph a classification or regression tree with a hierarchical
  tree diagram, optionally including colored symbols at leaves 
  and additional info at intermediate nodes.
}
\usage{
  draw.tree (tree, cex=par("cex"), pch=par("pch"), size=2.5*cex, 
      col=NULL, nodeinfo=FALSE, units="", cases="obs", 
      digits=getOption("digits"), print.levels=TRUE, 
      new=TRUE)
}
\arguments{
  \item{tree}{object of class \code{rpart} or \code{tree}.}
  \item{cex}{size of text, par parameter.}
  \item{pch}{shape of symbol at leaves, par parameter.}
  \item{size}{if \code{size=0}, draw terminal symbol at leaves else a
              symbol of size in cex units.}
  \item{col}{vector of colors from \code{\link{hsv}}, \code{\link{rgb}},
             etc, or if \code{NULL}, then use \code{\link{rainbow}}.}
  \item{nodeinfo}{if \code{TRUE}, add a line at each node with mean value
                  of response, number of observations, and percent
                  deviance explained (or classified correct).}
  \item{units}{label for units of mean value of response, if
               regression tree.}
  \item{cases}{label for type of observations.}
  \item{digits}{number of digits to round mean value of response, if
                regression tree.}
  \item{print.levels}{if \code{TRUE}, print levels of factors at splits,
                      otherwise only the factor name.}
  \item{new}{if \code{TRUE}, call \code{\link{plot.new}}.}
}
\details{
  As in \code{plot.rpart(,uniform=TRUE)}, each level has constant depth.  
  Specifying \code{nodeinfo=TRUE}, shows the deviance explained or the 
  classification rate at each node.

  A split is shown, for numerical variables, as 
  \code{variable <> value} when the cases with lower values go left, or as 
  \code{variable >< value} when the cases with lower values go right.
  When the splitting variable is a factor, and print.levels=TRUE, 
  the split is shown as \code{levels = factor = levels} with the cases 
  on the left having factor levels equal to those on the left of the 
  factor name, and correspondingly for the right.
}
\value{The vector of colors supplied or generated.}
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link[rpart]{rpart}}, 
  \code{\link[tree]{tree}}, 
  \code{\link{draw.clust}}, 
  \code{\link{map.groups}}
}
\examples{
  library (rpart)
  data (oregon.env.vars)

  draw.tree (clip.rpart (rpart (oregon.env.vars), best=7), 
      nodeinfo=TRUE, units="species", cases="cells", digits=0)
}
\keyword{hplot}
\keyword{tree}
