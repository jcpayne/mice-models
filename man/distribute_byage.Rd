% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MICE_initialization.R
\name{distribute_byage}
\alias{distribute_byage}
\title{Distribute initial recruitment by age.}
\usage{
distribute_byage(recruits, Mv, P)
}
\arguments{
\item{Mv}{Natural mortality (log scale) by age. Vector.}

\item{P}{Maximum age (age plus-group).  Scalar}

\item{recruits.}{Number of recruits in year 1 (scalar).}
}
\value{
Vector of number-at-age in timestep 1 (or year 1)
}
\description{
Distributes an initial recruitment of a species to all ages within timestep
1, by applying cumulative natural mortality.
}
\seealso{
Other initialization: \code{\link{calibrate}},
  \code{\link{distribute_bytimestep}},
  \code{\link{initialize_N}}, \code{\link{set_simParam}},
  \code{\link{solve_age0_survival}}
}
\keyword{internal}

