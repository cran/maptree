\name{ngon}
\alias{ngon}
\title{Outline or Fill a Regular Polygon}
\description{
 Draws a regular polygon at specified coordinates as an
 outline or shaded.
}
\usage{
  ngon (xydc, n=4, angle=0, type=1)
}
\arguments{
  \item{xydc}{four element vector with \code{x} and \code{y}
              coordinates of center, \code{d} diameter in mm, and
              \code{c} color.}
  \item{n}{number of sides for polygon (>8 => circle).}
  \item{angle}{rotation angle of figure, in degrees.}
  \item{type}{\code{type=1} => interior filled, \code{type=2} => edge,
              \code{type=3} => both.}
}
\details{
 Uses \code{\link{polygon}} to draw shaded polygons and
 \code{\link{lines}} for outline.  If n is odd, there is
 a vertex at (0, d/2), otherwise the midpoint of a side is
 at (0, d/2).
}
\value{Invisible.}
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link{polygon}}, 
  \code{\link{lines}}, 
  \code{\link{map.key}}, 
  \code{\link{map.groups}}
}
\examples{
  plot (c(0,1), c(0,1), type="n")
  ngon (c(.5, .5, 10, "blue"), angle=30, n=3)
  apply (cbind (runif(8), runif(8), 6, 2), 1, ngon)
}
\keyword{aplot}
