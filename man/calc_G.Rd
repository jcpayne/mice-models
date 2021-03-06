% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MICE_initialization.R
\name{calc_G}
\alias{calc_G}
\title{Create environmental driver, G}
\usage{
calc_G(Species, simParams)
}
\arguments{
\item{Species}{A Prey species object.}

\item{simParams}{A simulationParameters object}
}
\value{
A modified Prey object, with the G slot set.
}
\description{
Generates a unique time series for each prey species to simulate an
environment that switches between high (good reproduction) and low (poor
reproduction) states.  The function samples from a list of period lengths for
low and high periods that is read in from an internal file
(envir_periods.csv).  The period lengths and reproduction values for the
periods are balanced to match observed patterns.
}
\keyword{internal}

