\name{hex.env}
\alias{hex.env}
\title{Environmental Variables for Oregon, USA}
\description{
  Distributions of 10 environmental variables over 404 grid 
  cells in Oregon, USA.
}
\usage{data (hex.env)}
\format{A data frame with 404 rows and 10 columns.}
\details{
  Row labels are hexagon identifiers from White et al. (1992).
  Variables are 

  \tabular{ll}{
    bird.spp \tab number of native breeding bird species\cr
    x \tab x coordinate of center of grid cell\cr
    y \tab y coordinate of center of grid cell\cr
    jan.temp \tab mean minimum January temperature (C)\cr
    jul.temp \tab mean maximum July temperature (C)\cr
    rng.temp \tab mean difference between July and January temperatures (C)\cr
    ann.ppt \tab mean annual precipitation (mm)\cr
    min.elev \tab minimum elevation (m)\cr
    rng.elev \tab range of elevation (m)\cr
    max.slope \tab maximum slope (percent)
  }}
}
\source{Denis White, \email{white.denis@epa.gov} }
\references{
  White, D., Preston, E.M., Freemark, K.E., Kiester, A.R. (1999) 
  A hierarchical framework for conserving biodiversity, 
  \emph{Landscape ecological analysis: issues and applications}, 
  Klopatek, J.M., Gardner, R.H., editors, Springer-Verlag, pp. 
  127-153.

  White, D., Kimerling, A.J., Overton, W.S.  (1992)  Cartographic  
  and geometric components of a global sampling design for 
  environmental monitoring, \emph{Cartography and Geographic 
  Information Systems}, \bold{19}(1), 5-22.
}
\seealso{
  \code{\link{hex.spp}}, 
  \code{\link{oregon}}
}
\keyword{datasets}
