\name{oregon.bird.dist}
\alias{oregon.bird.dist}
\title{Presence/Absence of Bird Species in Oregon, USA}
\description{
  Binary matrix (1 = present) for distributions of 248 native
  breeding bird species for 389 grid cells in Oregon, USA.
}
\usage{data (oregon.bird.dist)}
\format{A data frame with 389 rows and 248 columns.}
\details{
  Row names are hexagon identifiers from White et al. (1992).
  Column names are species element codes developed by
  The Nature Conservancy (TNC),
  the Oregon Natural Heritage Program (ONHP), and
  the Association for Biodiversity Information (ABI).
}
\source{Denis White, \email{white.denis@epa.gov} }
\references{
  Master, L. (1996) Predicting distributions for vertebrate
  species: some observations, \emph{Gap Analysis: A Landscape
  Approach to Biodiversity Planning}, Scott, J.M., Tear, T.H.,
  and Davis, F.W., editors, American Society for Photogrammetry
  and Remote Sensing, Bethesda, MD, pp. 171-176.

  White, D., Preston, E.M., Freemark, K.E., Kiester, A.R. (1999)
  A hierarchical framework for conserving biodiversity,
  \emph{Landscape ecological analysis: issues and applications},
  Klopatek, J.M., Gardner, R.H., editors, Springer-Verlag, pp.
  127-153.

  White, D., Kimerling, A.J., Overton, W.S.  (1992)  Cartographic
  and geometric components of a global sampling design for
  environmental monitoring, \emph{Cartography and Geographic
  Information Systems}, \bold{19}(1), 5-22.

  TNC, \url{http://nature.org/}

  ONHP, \url{http://www.abi.org/nhp/us/or/}

  ABI, \url{http://www.abi.org/}

  ABI species data, \url{http://www.abi.org/datasets_zoo/overview.htm}
}
\seealso{
  \code{\link{oregon.env.vars}}, 
  \code{\link{oregon.bird.names}}, 
  \code{\link{oregon.grid}}, 
  \code{\link{oregon.border}}
}
\keyword{datasets}
