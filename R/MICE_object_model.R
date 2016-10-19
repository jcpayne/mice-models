### Create S4 Classes
#' An S4 class to hold global parameters and input data for simulations.
#'
#' Species-specific parameter values are held in two dataframes.  Global values
#' have their own slots.  This class is used to read input data from two external
#' files, and is later used to initialize Species objects with species-specific
#' parameter values.
#' @slot name character. The object name (should match the name of the class).
#' @slot nYears numeric. Number of years for the simulation.
#' @slot nSteps numeric. Number of timesteps per year (may not be used by all species)
#' @slot scenario character. A name for the scenario, which must match a column name in input_params.csv.
#' @slot params A dataframe that holds parameter values for everything except diet preferences.
#' @slot diet_prefs A dataframe that holds diet preferences for predators.
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
   diet_prefs = "data.frame"
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
#' @slot am numeric. Age at maturity; a vector.
#' @slot mv numeric. Proportion mature at age; a vector.
#' @slot Ninit numeric. Numbers at age in year 1 (or timestep 1 of year 1); a vector.
#' @slot Nmat matrix.  Population matrix, age x time (years, or timesteps x years).
#' @slot wv numeric. Weight at age; a vector.
#' @slot M numeric. Natural mortality rate (log scale); a scalar.
#' @slot fished logical.  Is the species fished?
#' @slot fishing_rule list. Parameters for the fishing rule model; a list.
#' @slot Sv numeric. Fishing selectivity at age; a vector.
#' @slot Fv numeric. Fishing mortality rate (log scale) by time; a vector.
#' @slot F_baseline numeric. Baseline fishing mortality rate; a scalar.
#' @slot gv numeric. Proportion of F that occurs per timestep; a vector.
#' @slot Cv numeric. Catch by time; a vector.
#' @slot Zv matrix. Total mortality rate (log scale) by time: a vector.
#' @slot R_baseyear numeric. Recruitment in a year 1, if supplied; a scalar.
#' @slot R numeric. Recruitment by year; a vector.
#' @slot lnSR_err numeric. Stock-recruitment errors (ln scale) by time; a vector.
#' @slot recruitment list. Parameters used in recruitment model (may vary by species)
#' @slot B1P numeric. Age 1+ biomass by time; a vector.
#'
#' @return Should not be instantiated; rather, instantiate the child classes.
#' @export
setClass("Species",
  slots = list(
    name="character",
    uses_timestep = "logical",
    P = "numeric",
    am = "numeric",
    mv = "numeric",
    Ninit = "numeric",
    Nmat = "matrix",
    wv = "numeric",
    M = "numeric",
    fished = "logical",
    fishing_rule = "list",
    Sv = "numeric",
    Fv = "numeric",
    F_baseline = "numeric",
    gv = "numeric",
    Cv = "numeric",
    Zv = "matrix",
    R_baseyear = "numeric",
    R = "numeric",
    lnSR_err = "numeric",
    recruitment = "list",
    B1P = "numeric"
  ),
  prototype = list(
    fished = FALSE
  )
)

#' A S4 parent class for Prey.
#'
#' Inherits from the Species class, and has slots for values that are common to
#' prey, but not predators.
#' @slot availability list. Parameters of the spatial availability function.
#' @slot observed_SR_err numeric. Observed stock-recruit residuals (ln scale), to resample from.
#'
#' @return Should not be instantiated; rather, instantiate the child classes.
#' @export
#' @seealso parent class \linkS4class{Species}
setClass("Prey",
  slots = list(
    B_sp = "numeric",
    availability = "list",
    observed_SR_err = "numeric"
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
#' @slot prey_effect list. Parameters for the model that determines how much of an impact prey biomass has on predator recruitment
#' @slot Dy numeric. Available, scaled total prey biomass; a vector (time).
#' @slot BHp numeric. Beverton-Holt p parameter, which is affected by prey; a vector (time).
#' @slot Mo numeric. Maximum number of offspring; a scalar.
#' @slot K numeric. Carrying capacity; a scalar.
#' @slot Nm numeric. Number of mature animals; a vector.
#'
#' @return Should not be instantiated; rather, instantiate the child classes.
#' @export
#' @seealso parent class \linkS4class{Species}
setClass("Predator",
  slots = list(
    sp_prefs = "list",
    age_prefs = "list",
    omegas = "list",
    prey_effect = "list",
    Dy = "numeric",
    BHp = "numeric",
    Mo = "numeric",
    K = "numeric",
    Nm = "numeric"
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
#' @return An object of the Anchovy class.
#' @export
#' @seealso parent classes \linkS4class{Prey} and \linkS4class{Species} (top-level).
#' @examples
#' anch<-new("Anchovy")
#' slotNames(anch) #see list of slots
#' anch@Zv[,1:10] #show total mortality by age (Zv) for timesteps 1 to 10
setClass("Anchovy",
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
#' sar@Zv[,1:10] #show total mortality by age (Zv) for timesteps 1 to 10

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
#' oth@Zv[,1:10] #show total mortality by age (Zv) for timesteps 1 to 10
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
#' lion@Zv[,1:10] #show total mortality by age (Zv) for timesteps 1 to 10
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
#' pel@Zv[,1:10] #show total mortality by age (Zv) for timesteps 1 to 10
#' @seealso parent classes \linkS4class{Predator} and \linkS4class{Species} (top-level).
setClass("Pelican",
  prototype = list(name = "Pelican"),
  contains = "Predator"
)

### Intialization functions
#These will load our objects with parameter values.

#First, we
#setOldClass("data.frame") #register the S3 "data.frame" class so it can be used in an S4 method signature.

#' Set parameter values in a simulationParams object.
#'
#' Reads two files (input_params.csv, which holds most of the input parameters,
#' and diet_prefs.csv, which holds predator by species dietary preferences) into
#' two data frames, and loads the data frames into an object of class simulationParams.
#' @param params  A simulationParams object.
#' @param scenario Character string; a scenario name (must match the name in input_params.csv)
#' @param nYears Integer.  Number of years to simulate.
#' @param nSteps Integer.  Number of timesteps per year.
#' @param ... Other arguments.
#' @return A modified simParams object, with values filled from the file.
#' @family parameter setting
#' @family initialization
#' @seealso \linkS4class{simulationParams}
#' @export
#' @keywords internal.
set_simParam<-function(params, inputdir,scenario, nYears, nSteps) {
  if(nSteps > 365) stop(paste("Are you sure that you want",nSteps,"timesteps within a year?"))

  #Read in and check the parameter and diet preferences files.
  sp<-params #the simulationParams object

  #Either get input filepath from function arguments, or look in package example data.
  if(missing(inputdir)){
    param_file<-system.file("extdata", "input_params.csv", package = "MICE")
    print("Species parameters are being loaded from the example input_params.csv file included with the MICE package")
    prefs_file<-system.file("extdata", "diet_prefs.csv", package = "MICE")
    print("inputdir missing: Diet preferences loaded from the example diet_prefs.csv file included with the MICE package")
  } else{
    param_file<-file.path(inputdir,"input_params.csv")
    prefs_file<-file.path(inputdir,"diet_prefs.csv")
  }
  paramdf<-read.csv(param_file,stringsAsFactors=FALSE,blank.lines.skip=T,header=T,sep=",")
  paramdf[paramdf==""]<-NA #Don't skip this step!
  #remove empty rows
  paramdf<-paramdf[rowSums(is.na(paramdf)) != ncol(paramdf),]
  prefs<-read.csv(prefs_file,stringsAsFactors=FALSE,header=T,sep=",")
  prefs[prefs==""]<-NA #Don't skip this step!
  if(missing(paramdf)||nrow(paramdf)==0) stop("Parameter file is missing or contains no data")
  if(missing(prefs)||nrow(prefs)==0) stop("Diet preference file is missing or contains no data")
  if(is.na(match(scenario,names(paramdf)))) stop("The scenario name does not match a column name in the input parameter file")

  #Since we passed all of the checks, assign the values to the slots of the object
  sp@nYears<-nYears
  sp@nSteps<-nSteps
  sp@scenario<-scenario
  sp@params<-paramdf
  sp@diet_prefs<-prefs
  sp
}

####Utility functions for initialization
#' Set diet preference slots in a Predator object.
#'
#' The function reads values from a simulationParams-class object and fills the
#' dietary preference slots in the Predator object.
#'
#' @param Predator A Predator object.
#' @param simParams A simulationParams object containing input data.
#'
#' @return A modified Predator object, with diet preference slots filled.
#' @export
#' @family parameter setting
#' @seealso \code{\link{get_params}} for reading in parameter values and
#'   \code{\link{set_species}} for calling both functions.
#' @keywords internal
get_prefs<-function(Predator,simParams){
  scenario<-simParams@scenario
  prefs<-simParams@diet_prefs
  predatorx<-Predator@name
  prey<-unique(prefs$prey[!is.null(prefs$scenario)]) #List the prey species in the file
  prey<-prey[order(prey)] #sort them alphabetically

  #Anticipating that the user may not want to change the diet preferences as often
  #as the other parameters, we revert to the "base" scenario if the scenario name
  #being run is not found, or all of the rows are blank under that scenario name.
  if(!(scenario %in% names(prefs)) || all(is.na(prefs[,scenario]))){
    if("base" %in% names(prefs)){
      print(paste("Diet preferences loaded from the base scenario because scenario",scenario,"is missing in diet_prefs."))
      scenario<-"base"
    } else{
      stop(paste("Neither scenario",scenario,"nor 'base' found in diet_prefs."))
    }
  }

  #Get the species preference for this predator, sorted in order of prey name (alphabetical). Returns a
  #named list of species preferences.
  for(i in 1:length(prey)){
    sp_prefx<-prefs[prefs$predator==predatorx & prefs$age == -1 & prefs$prey==prey[i],match(scenario,names(prefs))]
    if(i==1){
      sp_prefs<-list(sp_prefx)
    } else{
      sp_prefs<-c(sp_prefs,list(sp_prefx))
    }
  }
  names(sp_prefs)<-prey

  #Get the age preferences for each species in a list
  for(j in 1:length(prey)){
    temp<-prefs[prefs$predator==predatorx & prefs$age!= -1 & prefs$prey==prey[j],]
    temp<-temp[order(temp[,2]),] #sort on age, just in case
    age_pref_spx<-temp[,match(scenario,names(temp))]
    if(j==1){
      age_prefs<-list(age_pref_spx)
      names(age_prefs)[1]<-prey[j]
    } else{
      age_prefs<-c(age_prefs,list(age_pref_spx))
      names(age_prefs)[length(age_prefs)]<-prey[j]
    }
  }#for 1 to length(prey)

  #Multiply each age preference vector by the corresponding species preference
  #for (i in 1:length(prey)){
  #  age_prefs[[i]]<-age_prefs[[i]] * sp_pref[i]
  #}
  Predator@sp_prefs<-sp_prefs
  Predator@age_prefs<-age_prefs #fill the age_prefs slot in the Predator object with a matrix
  Predator
}#function get_prefs()

#' Set parameter values in a Species object
#'
#' Reads values from a simulationParams class object and sets them
#' (appropriately customized for each species) in a Species-class object.  The
#' parameters for certain models are encapsulated in lists, which makes it
#' easier to change them later.
#'
#' @param Species.  A Species object.
#' @param simParams A simulationParams object containing input data.
#'
#' @return A modified Species object, with parameters set from the input data.
#' @export
#'
#' @family parameter setting
#' @seealso \code{\link{get_prefs}} for reading in diet preferences,  and
#'   \code{\link{set_species}} for calling both functions.
#' @keywords internal
get_params<-function(Species,simParams){
  scenario<-simParams@scenario
  params<-simParams@params
  nSteps<-simParams@nSteps
  nYears<-simParams@nYears

  #First, set the scalar parameters that are in lists (param_list is not NA)
  temp<-params[params$species==Species@name & is.na(params$index) & !is.na(params$param_list),]
  paramlists<-unique(temp$param_list) #get the names of the param_lists
  if(length(paramlists)>0){
    #Step through the param_lists
    for(i in 1:length(paramlists)){

      #Put the elements of the param_list into a list and name them
      p.elements<-temp$parameter[temp$param_list==paramlists[i]]
      p.values<-temp[temp$param_list==paramlists[i],match(scenario,names(temp))]
      plist<-NA
      for(j in 1:length(p.elements)){
        if(j==1){
          plist<-list(p.values[j])
        } else{
          plist<-c(plist,list(p.values[j]))
        }
      }#for j
      names(plist)<-p.elements
      #Set the slot in the Species object
      slot(Species,as.character(paramlists[i]))<-plist
    }#for i in length(paramlists)
  }#paramlist length > 0


  #Next, set the scalar parameters that are not in lists
  #Get the params
  temp<-params[params$species==Species@name & is.na(params$index) & is.na(params$param_list),]
  if(nrow(temp) > 0){
    for(i in 1:nrow(temp)){
      val<-temp[i,scenario] #the parameter value
      slotname<-as.character(temp$parameter[i])
      slottype<-getSlots("Species")[match(slotname,names(getSlots("Species")))]
      if(slottype=="logical" && (val==1 || val==0)){
        val<-as.logical(val) #coerce 1/0 to logical if needed
      }
      slot(Species,slotname)<-val
    }
  }#nrow(temp) > 0

  #Save nAges for validation
  if(!is.na(match("P",slotNames(Species)))){
    nAges<-Species@P + 1
  } else{
    stop(paste("maxAge 'P' slot missing for",Species@name))
  }

  #Finally, deal with the by-age or by-timestep parameters, which have to be turned into vectors (possibly of uneven length)
  temp<-params[params$species==Species@name & !is.na(params$index),]
  if(nrow(temp) > 0){
    by_index_params<-unique(temp$parameter)
    for(i in 1:length(by_index_params)){
      #Create a three-column dataframe for index-parameter i (1=param_type, 2=index, 3=value for the given scenario)
      index_param<-temp[temp$parameter==by_index_params[i],c(match("param_type",names(temp)),match("index",names(temp)),match(scenario,names(temp)))]
      #Sort by index
      index_param<-index_param[order(index_param$index),]
      #Validation (check the length of by-age and by-timestep parameters)
      ptype<-unique(index_param$param_type) #age or timestep
      if(ptype=="age" && nrow(index_param) != nAges){
        stop(paste("Input parameter",by_index_params[i],"is",length(index_param),"long, but nAges is",nAges))
      }
      if(ptype=="timestep" && nrow(index_param) != nSteps){
        stop(paste("Input parameter",by_index_params[i],"is",length(index_param),"long, but nSteps is",nSteps))
      }
      #Copy it to the Species object
      slot(Species,as.character(by_index_params[i]))<-index_param[,3]
    }
  }#nrow(temp) > 0
  Species
}#function get_params()

#Two functions to spread an initial recruitment over ages and/or timesteps in year 1.

#' Distribute initial recruitment by age.
#'
#' Distributes an initial recruitment of a species to all ages within timestep
#' 1, by applying cumulative natural mortality.
#' @param recruits.  Number of recruits in year 1 (scalar).
#' @param M Natural mortality (log scale). Scalar.
#' @param P Maximum age (age plus-group).  Scalar
#'
#' @return Vector of number-at-age in timestep 1 (or year 1)
#' @export
#' @family initialization
#' @keywords internal
distribute_byage<- function(recruits,M,P) {
  Nv <- numeric(P+1) #Age classes go from 0 to P, but the age vector is 1-based, so goes from 1:(P+1)
  for (i in 1:(P+1)) {
    if(i == 1 ) {
      Nv[i] <- recruits
    } else {
      Nv[i] <- Nv[i-1] * exp(-M)
    }
  }
  Nv
}

#' Distribute numbers-at-age by timestep
#'
#'#Take an initial vector of population by age for a species and distribute it
#'across all timesteps within one year.  Apply cumulative natural mortality to each
#'step.

#' @param Nv Population size by age in timestep 1; a vector.
#' @param M Natural mortality rate (log scale); scalar.
#' @param L Number of timesteps per year; scalar.
#'
#' @return A small matrix of population size in year 1; dim=(nAges x L).
#' @export
#' @family initialization
#' @keywords internal
distribute_bytimestep<-function(Nv,M,L){
  Ny1mat<-matrix(data=NA,nrow=length(Nv),ncol=L)
  Ny1mat[,1]<-Nv
  for(i in 2:ncol(Ny1mat)){
    Ny1mat[,i]<-Ny1mat[,(i-1)] * exp(-M/L)
  }
  Ny1mat
}

#' Initialize a Species object.
#'
#' This function creates the population matrix, Nmat, for a Species, whose dimension is
#' either (nAges x nYears) if use_timestep is FALSE, or (nAges x (nYears *
#' nSteps)) if use_timestep is TRUE.  It fills Nmat[,year==1], based on values in the
#' simulationParams object.  It also dimensions empty vectors that will hold the
#' value of recruitment residuals, fishing mortality and catch (for fished
#' species), age 1+ biomass (for prey) and proportion mature-at-age, and calculates
#' Nm[1], and sets R[1] to zero.
#'
#' @param Species. A species object.
#' @param simulationParams.  A simulationParams object that holds global values and input data.
#'
#' @return A modified Species object, with arrays initialized.
#' @export
#' @keywords internal
#' @family parameter setting
#' @family initialization
initialize_N<-function(Species,simulationParams){
  L<-simulationParams@nSteps
  y<-simulationParams@nYears

  Ninit<-Species@Ninit #initial vector of numbers-at-age
  R_baseyear<-Species@R_baseyear #year-1 recruitment (an alternative to Ninit)
  M<-Species@M #natural mortality
  P<-Species@P #plus-group age
  am<-Species@am #age at maturity

  #Either read in the first year of numbers-at-age, or calculate it from an initial recruitment
  if(length(Ninit)==(P+1)){
    Nv<-Ninit
  } else{
    #Distribute initial recruits across age classes by applying cumulative natural mortality,
    #to create the by-age population vector for the first timestep
    Nv<-distribute_byage(R_baseyear,M,P)
  }

  #Now dimension the population matrix and distribute the first population vector across L timesteps if needed,
  #applying cumulative natural mortality
  if(Species@uses_timestep){
    #Dimension the Nmat matrix (P+1 age classes down, (y*L) timesteps across)
    partialmat<-matrix(data=NA,nrow=P+1,ncol=(L*(y-1))) #everything but year 1
    yr1mat<-distribute_bytimestep(Nv,M,L) #year 1
    Nmat<-cbind(yr1mat,partialmat)
    Species@Nmat<-Nmat
  }
  else {
    #Dimension Nmat as (P+1 age classes down, y years across), and fill the first column (for a species that doesn't use timestep)
    Species@Nmat<-matrix(data=NA,nrow=P+1,ncol=y)
    Species@Nmat[,1]<-Nv
  }

  nTimes<-dim(Species@Nmat)[2] #the time dimension; either years or years x timesteps

  #Dimension vectors to hold age 1+ biomass and spawning biomass for prey species
  if(is(Species,"Prey")){
    Species@B1P<-numeric(nTimes)
    Species@B_sp<-numeric(nTimes)
  }

  #Dimension an array to hold mortality rates by age (for debugging)
  #TODO: no reason to keep the whole by-age matrix; keep a total mortality vector only, once running correctly
  Species@Zv<-matrix(data=NA,nrow=P+1,ncol=nTimes)

  #Dimension the vector of recruitments (will be NA in timesteps where there is no recruitment)
  Species@R<-as.numeric(rep(NA,nTimes))

  #Dimension the vector of recruitment residuals
  Species@lnSR_err<-numeric(nTimes)

  #Dimension the vectors for prey biomass (Dy) and prey effect (BHp)
  if("Dy" %in% slotNames(Species)){
    Species@Dy<-as.numeric(rep(NA,nTimes))
  }
  if("BHp" %in% slotNames(Species)){
    Species@BHp<-as.numeric(rep(NA,nTimes))
  }

  #Generate the recruitment residuals vector (epsilon) for species that have stochastic recruitment error,
  #either by resampling or by random draws from a distribution.  epsilon is a vector nYears or (nSteps x nYears) long.
  if(Species@recruitment$R_stochastic){
    tmp<-Species@recruitment
    if(exists('resample_residuals',where=tmp) && Species@recruitment$resample_residuals){
      if(length(Species@observed_SR_err) > 0){
        Species<-resample_R_error(Species)
      } else{
        warning("Residuals will not be resampled because length(observed_SR_err) is zero.\n
                They will be drawn from a log-normal distribution instead.")
        Species<-calc_R_error(Species)
      }
    }else{
      Species<-calc_R_error(Species)
    }
    }#if(R_stochastic)

  #Set the proportion-mature-at-age vector, mv (for now, just a step function where mature=0 < am <= mature=1)
  mature<-numeric(P + 1)
  mature[(am+1):length(mature)]<-1
  Species@mv<-mature

  #Fill in the number mature, Nm
  if("Nm" %in% slotNames(Species)){
    Species@Nm[1]<-calc_Nm(Nv,Species@mv)
  }

  #Dimension the fishing mortality and catch vectors for species that are fished.
  if(Species@fished){
    if(length(Species@Sv)==0){
      stop(paste(Species@name,"is fished, but fishing selectivity Sv is missing from input_params.csv"))
    }
    Species@Fv<-numeric(nTimes) #Initialize as zeros (NAs would cause problems)
    ctch<-as.numeric(rep(NA,nTimes))
    ctch[1:L]<-0 #The first 1 year is not fished
    Species@Cv<-ctch
  }

  #Return the Species object
  Species
}

####Top-level functions to do initialization
#' Set initial parameter values in Species objects.
#'
#' @param species. A Species object.
#' @param params. A simParams object.
#' @param ... Other arguments (from child methods).
#'
#' @return A child method, depending on signature.
#' @export
#' @keywords internal
#' @family parameter setting
#' @examples speciesA<-set_species(speciesA)
#' @rdname set_species
setGeneric("set_species",
  def = function(species,params,...) standardGeneric("set_species")
)

#' The default method sets parameter values in a Species object.  It calls
#' functions to read in parameters from input files, and to initialize arrays.
#'
#' @describeIn set_species Default method.
setMethod("set_species",
  signature("Species","simulationParams"),
  function(species,params,...) {
    species<-get_params(species,params)
    species<-initialize_N(species,params) #Note: this must come after get_params
    species
})

#' The method for predators also reads in a predator's diet preferences (for
#' prey species and the age-classes of those prey).
#' @describeIn set_species Also reads in diet preferences.
setMethod("set_species",
  signature("Predator","simulationParams"),
  function(species,params,...) {
    species<-get_params(species,params)
    species<-get_prefs(species,params)
    species<-initialize_N(species,params) #Note: must come after get_params
    species
  })
