\name{map.key}
\alias{map.key}
\title{Draw Key to accompany Map of Groups}
\description{
  Draws legends for maps of groups of observations.
}
\usage{
  map.key (x, y, labels=NULL, cex=par("cex"), pch=par("pch"),
      size=2.5*cex, col=NULL, head="", sep=0.25*cex, new=FALSE)
}
\arguments{
  \item{x, y}{coordinates of lower left position of key
              in proportional units (0-1) of plot.}
  \item{labels}{vector of labels for classes, or if \code{NULL},
                then integers \code{1:length(col)}, or \code{1}.}
  \item{size}{size in cex units of shaded key symbol.}
  \item{pch}{symbol number from \code{par("pch")} if < 100, otherwise
             parameter \code{n} for \code{\link{ngon}}.}
  \item{cex}{pointsize of text, \code{par} parameter.}
  \item{head}{text heading for key.}
  \item{sep}{separation in cex units between adjacent symbols in key.
             If \code{sep=0}, assume a continuous scale, use square
             symbols, and put labels at breaks between squares.}
  \item{col}{vector of colors from \code{\link{hsv}}, \code{\link{rgb}},
             etc, or if \code{NULL}, then use \code{\link{rainbow}}.}
  \item{new}{if \code{TRUE}, call \code{\link{plot}}.}
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
  kol <- gray (seq(0.8,0.2,length.out=2))
  map.key (0.2, 0.2, labels=c("absent","present"), pch=106, 
    col=kol, head="key", new=TRUE)
  # distribution of January temperatures
  cuts <- quantile (oregon.env.vars[["jan.temp"]], probs=seq(0,1,1/5))
  kol <- gray (seq(0.8,0.2,length.out=5))
  map.key (0.2, 0.2, labels=as.character(round(cuts,0)), 
    col=kol, sep=0, head="key", new=TRUE)

  # key for example in \code{help(group.tree)}
  map.key (0.2, 0.2, labels=as.character(seq(6)), 
    pch=19, head="node", new=TRUE)
}
\keyword{hplot}
\keyword{aplot}
