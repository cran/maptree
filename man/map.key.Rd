\name{map.key}
\alias{map.key}
\title{Draw Key to accompany Map of Groups}
\description{
  Draws legends for maps of groups of observations.
}
\usage{
  map.key (x, y, lables=NULL, new=FALSE,  size=5, 
      ps=par("ps"), pch=19, head="", sep=1, col=NULL)
}
\arguments{
  \item{x, y}{coordinates of lower left position of key
              in proportional units (0-1) of plot.}
  \item{labels}{vector of labels for classes, or if \code{NULL},
                then integers \code{1:length(col)}, or \code{1}.}
  \item{new}{if \code{TRUE}, call \code{\link{plot}}.}
  \item{size}{size in mm of shaded square symbol.}
  \item{ps}{pointsize of text, \code{par} parameter.}
  \item{pch}{symbol number from \code{par("pch")} if < 100, otherwise
             parameter \code{n} for \code{\link{ngon}}.}
  \item{head}{text heading for key.}
  \item{sep}{separation in mm between adjacent symbols in key.
             If \code{sep=0}, assume a continuous scale, use square
             symbols, and put lables at breaks between squares.}
  \item{col}{vector of colors from \code{\link{hsv}}, \code{\link{rgb}},
             etc, or if \code{NULL}, then use \code{\link{rainbow}}, or
             if \code{gray} or \code{grey}, then use \code{\link{gray}}.}
}
\details{
  Uses \code{\link{points}} or \code{\link{ngon}}, depending on value of
  \code{pch}, to draw shaded polygon symbols for key.
}
\value{The vector of colors supplied or generated.}
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link{ngon}}, 
  \code{\link{map.groups}}
}
\examples{
  data (oregon.env.vars)

  # key for examples in help(map.groups)
  # range map for American Avocet
  map.key (0.2, 0.2, lables=c("absent","present"), 
    col="gray", size=4, new=TRUE, sep=1, pch=106, head="key")
  # distribution of January temperatures
  cuts <- quantile (oregon.env.vars[["jan.temp"]], probs=seq(0,1,1/5))
  map.key (0.2, 0.2, lables=as.character(round(cuts,0)), 
    col="gray", size=4, new=TRUE, sep=0, pch=106, head="key")

  # key for example in \code{help(group.tree)}
  map.key (0.2, 0.2, lables=as.character(seq(6)), 
    size=4, new=TRUE, sep=1, pch=15, head="node")
}
\keyword{hplot}
\keyword{aplot}
