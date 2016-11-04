### Create S4 Classes
#' An S4 class to hold global parameters and input data for simulations.
#'
#' Species-specific parameter values are held in two dataframes.  Global values
#' have their own slots.  This class is used to read input data from external
#' files, and is later used to initialize Species objects with species-specific
#' parameter values.
#' @slot name character. The object name (should match the name of the class).
#' @slot nYears numeric. Number of years for the simulation.
#' @slot nSteps numeric. Number of timesteps per year (may not be used by all
#'   species)
#' @slot scenario character. A name for the scenario, which must match a column
#'   name in input_params.csv.
#' @slot params A dataframe that holds parameter values for everything except
#'   diet preferences.
#' @slot diet_prefs A dataframe that holds diet preferences for predators.
#' @slot envir_periods A dataframe that holds a set of environmental driver
#'   period lengths for each prey species
#' @slot anchovy_residuals. A set of observed stock-recruit residuals for anchovies.
#'
#' @return A simulationParams object.
#' @export
#'
setClass("simulationParams",
  slots = list(
   name = "character",
   nYears = "numeric",
   nSteps = "numeric",
   scenario = "character",
   params = "data.frame",
   diet_prefs = "data.frame",
   envir_periods = "data.frame",
   anchovy_residuals = "data.frame"
  ),
  prototype = list(
   nYears = 100,
   nSteps = 4
  )
)

#' An S4 class to represent a generic species.
#'
#' It holds values that are used to track the population growth (including
#' reproduction, survival, and exploitation) of every species, and all other
#' species objects inherit from this class. The variable names that end in 'v'
#' are vectors; the names that end in 'mat' are matrices.  The rest are scalars,
#' except for a few that (as shown) are lists.  Lists are used to hold
#' parameters for a model.  Depending on the model, those parameters may be
#' scalars, vectors, or a combination.  Holding model parameters in lists makes
#' it easier to swap models later.
#'
#' @slot name character. Species name.  Should match the class name.
#' @slot uses_timestep logical. Is the species tracked by timestep within a year (or just by year)?
#' @slot P numeric. Age plus-group; a scalar.
#' @slot am numeric. Age at maturity; a scalar.
#' @slot maturev numeric. Proportion mature at age; a vector.
#' @slot Ninit numeric. Numbers at age in year 1 (or timestep 1 of year 1); a vector.
#' @slot Nmat matrix.  Population matrix, age x time (years, or timesteps x years).
#' @slot wv numeric. Weight at age; a vector.
#' @slot Mv numeric. Natural mortality rate (log scale) by age; a vector.
#' @slot fished logical.  Is the species fished?
#' @slot fishing_rule list. Parameters for the fishing rule model; a list.
#' @slot Sv numeric. Fishing selectivity at age; a vector.
#' @slot Fv numeric. Fishing mortality rate (log scale) by time; a vector.
#' @slot F_baseline numeric. Baseline fishing mortality rate; a scalar.
#' @slot gv numeric. Proportion of F that occurs per timestep; a vector.
#' @slot Cy numeric. Catch per year; a vector.
#' @slot Ct numeric. Catch per timestep; a vector.
#' @slot Zmat matrix. Total mortality rate (log scale) by time: a vector.
#' @slot R_baseyear numeric. Recruitment in a year 1, if supplied; a scalar.
#' @slot R numeric. Recruitment by year; a vector.
#' @slot lnSR_err numeric. Stock-recruitment errors (ln scale) by time; a vector.
#' @slot recruitment list. Parameters used in recruitment model (may vary by species).
#' @slot R_error list. Parameters that control stochastic recruitment error.
#'
#' @return Should not be instantiated; rather, instantiate the child classes.
#' @export
setClass("Species",
  slots = list(
    name="character",
    uses_timestep = "logical",
    P = "numeric",
    am = "numeric",
    maturev = "numeric",
    Ninit = "numeric",
    Nmat = "matrix",
    wv = "numeric",
    Mv = "numeric",
    fished = "logical",
    fishing_rule = "list",
    Sv = "numeric",
    Fv = "numeric",
    F_baseline = "numeric",
    gv = "numeric",
    Cy = "numeric",
    Ct = "numeric",
    Zmat = "matrix",
    R_baseyear = "numeric",
    R = "numeric",
    lnSR_err = "numeric",
    recruitment = "list",
    R_error = "list"
  ),
  prototype = list(
    fished = FALSE
  )
)

#' A S4 parent class for Prey.
#'
#' Inherits from the Species class, and has slots for values that are common to
#' prey, but not predators.
#' @slot B_sp numeric. Spawning biomass by time; a vector.
#' @slot B1P numeric. Age 1+ biomass by time; a vector.
#' @slot B1PperB0 numeric.  Age 1+ biomass relative to unfished age-1+ biomass; a vector (time).
#' @slot B0 numeric. Average age-1+ biomass in unfished condition.
#' @slot G numeric. Environmental driver for reproduction.
#'
#' @return Should not be instantiated; rather, instantiate the child classes.
#' @export
#' @seealso parent class \linkS4class{Species}
setClass("Prey",
  slots = list(
    B1P = "numeric",
    B1PperB0 = "numeric",
    B_sp = "numeric",
    B0 = "numeric",
    G = "numeric"
  ),
  prototype = list(
    uses_timestep = TRUE
  ),
  contains = "Species"
)

#' An S4 parent class for Predators.
#'
#' Inherits from Species and has slots for values that are common to predators but not prey.
#'
#' @slot sp_prefs list. Prey species preferences (one list element per prey species).
#' @slot age_prefs list. A vector of predator preferences for prey age classes; one list element for each prey species.
#' @slot omegas list. Parameters that scale available prey biomass per species so they sum (roughly) to 1
#' @slot prey_availability list. Parameters used in model of spatial availability
#'   (one set per prey species in list elements named for the prey species.  May also be used to store data for graphing availability.)
#' @slot prey_effect list. Parameters for the model that determines how much of an impact prey biomass has on predator recruitment
#' @slot D0 numeric. The available, scaled total prey biomass in unfished condition.  A scalar.
#' @slot Dy numeric. Available, scaled total prey biomass; a vector (time).
#' @slot Fmax numeric.  Maximum fecundity (population lambda, from data).
#' @slot bigPhi numeric.  Scaling parameter for reproduction.
#' @slot phi numeric.  Sensitivity of predator reproduction to prey biomass
#' @slot BHp numeric. Beverton-Holt p parameter, which is affected by prey; a vector (time).
#' @slot Mo numeric. Maximum number of offspring; a scalar.
#' @slot Nm numeric. Number of mature predators; a vector (time).
#' @slot N1P numeric. Number of age-1+ predators; a vector (time).
#' @slot K1P numeric.  Number of age-1+ predators when prey are unfished
#' @slot NperK numeric. Number of age-1+ predators per age-1+ predators when prey are unfished.
#' @slot RperS numeric. Recruits per spawner; a vector (time).
#' @slot RperS_Punt numeric. Recruits per spawner (Punt method); a vector (time).
#'
#' @return Should not be instantiated; rather, instantiate the child classes.
#' @export
#' @seealso parent class \linkS4class{Species}
setClass("Predator",
  slots = list(
    sp_prefs = "list",
    age_prefs = "list",
    omegas = "list",
    prey_availability = "list",
    prey_effect = "list",
    D0 = "numeric",
    Dy = "numeric",
    Fmax = "numeric",
    bigPhi = "numeric",
    phi = "numeric",
    R_Punt = "numeric",
    BHp = "numeric",
    Mo = "numeric",
    Nm = "numeric",
    N1P = "numeric",
    K1P = "numeric",
    NperK = "numeric",
    RperS = "numeric",
    RperS_Punt = "numeric"
  ),
  prototype = list(
    uses_timestep = FALSE
  ),
  contains = "Species"
)

#Define each of the species classes:
#' An S4 class for Anchovy.
#'
#' Inherits from Prey.
#'
#' @slot p_fail The probability of recruitment failure
#' @slot R_residuals A set of observed recruitment residuals.
#'
#' @return An object of the Anchovy class.
#' @export
#' @seealso parent classes \linkS4class{Prey} and \linkS4class{Species} (top-level).
#' @examples
#' anch<-new("Anchovy")
#' slotNames(anch) #see list of slots
#' anch@Zmat[,1:10] #show total mortality by age (Zmat) for timesteps 1 to 10
setClass("Anchovy",
  slots=list(
    p_fail="numeric",
    R_residuals = "list"
  ),
  prototype=list(
    name="Anchovy"),
  contains = "Prey"
)

#' An S4 class for Sardine
#'
#' Inherits from Prey class.  By default, this species is fished.
#' @return An object of Sardine class.
#' @export
#' @seealso parent classes \linkS4class{Prey} and \linkS4class{Species} (top-level).
#' @examples
#' sar<-new("Sardine")
#' slotNames(sar) #see list of slots
#' sar@G #Get values of environmental driver, G.

setClass("Sardine",
  prototype = list(
    name="Sardine",
    fished = T
  ),
  contains = "Prey"
)

#' An S4 class for Other_forage
#'
#' Inherits from Prey class.  By default, this species does not use timestep and
#' it has a fixed recruitment (R_base).
#' @return An object of Other_forage class.
#' @export
#' @family species_objects
#' @examples
#' other<-new("Other_forage")
#' #' oth<-new("Other_forage")
#' slotNames(oth) #see list of slots
#' oth@Zmat[,1:10] #show total mortality by age (Zmat) for timesteps 1 to 10
#' @seealso parent classes \linkS4class{Prey} and \linkS4class{Species} (top-level).
setClass("Other_forage",
  slots = list(
    R_base = "numeric"
  ),
  prototype = list(
    name = "Other_forage",
    uses_timestep = F
  ),
  contains = "Prey"
)

#' An S4 class for Sea lions
#'
#' Inherits from Predator class.
#' @return An object of Sealion class.
#' @export
#' @seealso parent classes \linkS4class{Predator} and \linkS4class{Species} (top-level).
#' @examples
#' lion<-new("Sealion")
#' slotNames(lion) #see list of slots
#' lion@Zmat[,1:10] #show total mortality by age (Zmat) for timesteps 1 to 10
setClass("Sealion",
  prototype = list(name = "Sealion"),
  contains = "Predator"
)

#' An S4 class for Pelicans.
#'
#' Inherits from the Predator class.
#' @return An object of the Pelican class.
#' @export
#' @examples
#' pel<-new("Pelican")
#' slotNames(pel) #see list of slots
#' pel@Zmat[,1:10] #show total mortality by age (Zmat) for timesteps 1 to 10
#' @seealso parent classes \linkS4class{Predator} and \linkS4class{Species} (top-level).
setClass("Pelican",
  prototype = list(name = "Pelican"),
  contains = "Predator"
)

