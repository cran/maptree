\name{map.key}
\alias{map.key}
\title{Draw Key to accompany Map of Groups}
\description{
  Draws legends for maps of groups of observations.
}
\usage{
  map.key (x, y, lables=NULL, new=TRUE,  ps=par("ps"), 
      size=5, gon=10, head="", sep=1, col=NULL)
}
\arguments{
  \item{x, y}{coordinates of lower left position of key
              in proportional units (0-1) of plot.}
  \item{labels}{vector of labels for classes, or if NULL,
                then integers 1:length(col), or "1".}
  \item{new}{if TRUE, call \code{\link{plot}}.}
  \item{ps}{pointsize of text, par parameter.}
  \item{size}{size in mm of shaded square symbol.}
  \item{gon}{number of sides for polygon symbol (>8 => circle).}
  \item{head}{text heading for key.}
  \item{sep}{separation in mm between adjacent symbols in key.  
             If sep=0, assume a continuous scale, use gon=4,
             and put lables at breaks between squares.}
  \item{col}{vector of colors from \code{\link{hsv}}, \code{\link{rgb}} 
             etc, or if NULL, then generated as equally spaced hues, or
             if "gray", then generated as gray tones.}
}
\details{
  Uses \code{\link{ngon}} to draw shaded polygon symbols for
  key.
}
\value{The vector of colors supplied or generated.}
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link{ngon}}, 
  \code{\link{map.groups}}
}
\examples{
  data (hex.env)

  # key for examples in help(map.groups)
  # range map for American Avocet
  map.key (0.2, 0.2, lables=c("absent","present"), col="gray",
      size=3, gon=6, new=TRUE, sep=1, head="key")
  # distribution of January temperatures
  cuts <- quantile (hex.env[["jan.temp"]], probs=seq(0,1,1/5))
  map.key (0.2, 0.2, lables=as.character(round(cuts,0)), 
      col="gray", size=3, new=TRUE, sep=0, head="key")

  # key for example in help(group.tree)
  map.key (0.2, 0.2, lables=as.character(seq(6)), 
      size=3, gon=6, new=TRUE, sep=1, head="node")
}
\keyword{hplot}
\keyword{aplot}
