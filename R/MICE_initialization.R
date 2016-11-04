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

  #read in the environmental driver data from the internal file
  envt_file<-system.file("extdata","envir_driver.csv",package="MICE")
  envir_periods<-read.csv(envt_file,stringsAsFactors=FALSE,blank.lines.skip=T,header=T,sep=",")

  #Read in the anchovy residuals from the internal file
  anch_file<-system.file("extdata","anchovy_residuals.csv",package="MICE")
  anchovy_residuals<-read.csv(anch_file,stringsAsFactors=FALSE,blank.lines.skip=T,header=T,sep=",")

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
  sp@envir_periods<-envir_periods
  sp@anchovy_residuals<-anchovy_residuals
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
      print("Diet preferences loaded from the base scenario")
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
  #print("Diet preferences loaded")
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
      p.preycol<-temp$prey[temp$param_list==paramlists[i]]
      listdata<-data.frame(p.elements,p.values,p.preycol)
      names(listdata)<-(c("param","vals","preycol"))
      plist<-list()
      #Handle nested two-level list
      if(!all(is.na(p.preycol))){
        preyspp<-unique(p.preycol)

        for(j in 1:length(preyspp)){
          #subset the data for this prey species
          preyxdata<-listdata[listdata$preycol==preyspp[j],]
          #create a list for the prey species with named elements (e.g.: Anchovy$ac)
          spxlist<-list()
          for(k in 1:nrow(preyxdata)){
            spxlist<-c(spxlist,list(preyxdata$vals[k]))
            names(spxlist)[length(spxlist)]<-as.character(preyxdata$param[k])
          }
          #add the species list to the main list and name it
          plist<-c(plist,list(spxlist))
          names(plist)[length(plist)]<-as.character(preyspp[j])
        }#for each prey species
      }#if prey col not all missing

      else{
        for(j in 1:length(p.elements)){
          plist<-c(plist,list(p.values[j]))
          names(plist)[length(plist)]<-p.elements[j]
        }
      }#nothing in prey col (i.e., a 1-level list)

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
      slottype<-getSlots(class(Species))[match(slotname,names(getSlots(class(Species))))]
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
  print(paste(Species@name,"initialized."))
  Species
}#function get_params()

#' Load anchovy recruitment residuals
#'
#'Loads two sets of observed anchovy recruitment residuals into the anchovy
#'object (for spawning stock biomass below or above 500,000 tons).
#'
#' @param anchovy An Anchovy object.
#' @param simParams A simulationParameters object.
#'
#' @return A modified Anchovy object, with the recruitment residuals set in R_residuals.
#' @export
#'
#' @keywords internal
load_residuals<-function(anchovy,simParams){
  resids<-simParams@anchovy_residuals #get the data frame
  le500<-resids$residual[resids$residual_set=="SSB_LE_500K"]
  gt500<-resids$residual[resids$residual_set=="SSB_GT_500K"]
  R_residuals<-list(SSB_LE_500K=le500, SSB_GT_500K=gt500)
  anchovy@R_residuals<-R_residuals
  anchovy
}

#Two functions to spread an initial recruitment over ages and/or timesteps in year 1.

#' Distribute initial recruitment by age.
#'
#' Distributes an initial recruitment of a species to all ages within timestep
#' 1, by applying cumulative natural mortality.
#' @param recruits.  Number of recruits in year 1 (scalar).
#' @param Mv Natural mortality (log scale) by age. Vector.
#' @param P Maximum age (age plus-group).  Scalar
#'
#' @return Vector of number-at-age in timestep 1 (or year 1)
#' @export
#' @family initialization
#' @keywords internal
distribute_byage<- function(recruits,Mv,P) {
  Nv <- numeric(P+1) #Age classes go from 0 to P, but the age vector is 1-based, so goes from 1:(P+1)
  for (i in 1:(P+1)) {
    if(i == 1) {
      Nv[i] <- recruits
    }
    if(i > 1 && i <= P){
      Nv[i] <- Nv[i-1] * exp(-Mv[i-1])
    }
    #plus-group is calculated as N*survival/(1-survival), to generate a stable age structure
    if(i==(P+1)){
      Nv[i]<-(Nv[i-1] * exp(-Mv[i-1]))/(1-Mv[i-1])
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
#' @param Mv Natural mortality rate (log scale) by age; a vector.
#' @param L Number of timesteps per year; scalar.
#'
#' @return A small matrix of population size in year 1; dim=(nAges x L).
#' @export
#' @family initialization
#' @keywords internal
distribute_bytimestep<-function(Nv,Mv,L){
  Ny1mat<-matrix(data=NA,nrow=length(Nv),ncol=L)
  Ny1mat[,1]<-Nv
  for(i in 2:ncol(Ny1mat)){
    for(j in 1:nrow(Ny1mat)){
      Ny1mat[j,i]<-Ny1mat[j,(i-1)] * exp(-Mv[j]/L)
    }
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
initialize_N<-function(Species,simulationParams,preynames){
  L<-simulationParams@nSteps
  nYears<-simulationParams@nYears

  Ninit<-Species@Ninit #initial vector of numbers-at-age
  R_baseyear<-Species@R_baseyear #year-1 recruitment (an alternative to Ninit)
  M<-Species@Mv #natural mortality by age: we expect a scalar in the input parameters,
  # however, we re-dimension it as a by-age vector, below.
  P<-Species@P #plus-group age
  am<-Species@am #age at maturity
  wv<-Species@wv #weight at age

  #Re-dimension natural mortality as a by-age vector, if a scalar has been provided
  if(length(M)==1){
    Mv<-numeric(P+1)
    Mv[]<-M #set natural mortality rate to M in all age classes
    Species@Mv<-Mv
  }

  #Calculate age-0 mortality and solve for bigPhi, for Predators.
  #Must come before distributing numbers-at-age.
  if(is(Species,"Predator")){
    Species<-solve_age0_survival(Species) #Adjust age-0 survival for predators
    Mv<-Species@Mv #Copy the updated Mv for convenience
    Species<-solve_fecundity(Species)
  }

  #Either read in the first year of numbers-at-age, or calculate it from an initial recruitment
  if(length(Ninit)==(P+1)){
    Nv<-Ninit
  } else{
    #Distribute initial recruits across age classes by applying cumulative natural mortality,
    #to create the by-age population vector for the first timestep
    Nv<-distribute_byage(R_baseyear,Mv,P)
  }

  #Now dimension the population matrix and distribute the first population vector across L timesteps if needed,
  #applying cumulative natural mortality
  if(Species@uses_timestep){
    #Dimension the Nmat matrix (P+1 age classes down, (nYears*L) timesteps across)
    partialmat<-matrix(data=NA,nrow=P+1,ncol=(L*(nYears-1))) #everything but year 1
    yr1mat<-distribute_bytimestep(Nv,Mv,L) #year 1
    Nmat<-cbind(yr1mat,partialmat)
    Species@Nmat<-Nmat
  }
  else {
    #Dimension Nmat as (P+1 age classes down, nYears years across), and fill the first column (for a species that doesn't use timestep)
    Species@Nmat<-matrix(data=NA,nrow=P+1,ncol=nYears)
    Species@Nmat[,1]<-Nv
  }

  nTimes<-dim(Species@Nmat)[2] #the time dimension; either years or years x timesteps

  #Set the proportion-mature-at-age vector, maturev (for now, just a step function where mature=0 < am <= mature=1)
  maturev<-numeric(P + 1) #0 by default
  maturev[(am+1):length(maturev)]<-1
  Species@maturev<-maturev

  #Fill in the number mature, Nm
  if("Nm" %in% slotNames(Species)){
    Species@Nm[1]<-calc_Nm(Nv,maturev)
  }

  #Fill in K1P for predators: age-1+ carrying capacity.
  #Since Nmat only depends on natural mortality and R_baseyear, this gives
  #the predator population size at carrying capacity without any prey effect.
  if(is(Species,"Predator")){
    Species@K1P<-sum(Species@Nmat[2:(P+1),1])
  }

  #Dimension vectors to hold age 1+ biomass and spawning biomass for prey
  #species, and fill year 1 of both
  if("B1P" %in% slotNames(Species)){
    Species@B1P<-numeric(nTimes)
    Species@B_sp<-numeric(nTimes)
    Species@B1PperB0<-as.numeric(rep(NA,nTimes))
    if(!Species@uses_timestep){
      Species@B1P[1]<-calc_B1P(Nv,wv)
      Species@B_sp[1]<-calc_B_sp(wv, Nv,maturev)
    } else{
      for(i in 1:L){
        Species@B1P[i]<-calc_B1P(Nmat[,i],wv)
        Species@B_sp[i]<-calc_B_sp(wv,Nmat[,i],maturev)
      }
    }#uses timestep
  }#prey

  #Dimension an array to hold mortality rates by age (for debugging)
  #TODO: no reason to keep the whole by-age matrix; keep a total mortality vector only, once running correctly
  Species@Zmat<-matrix(data=NA,nrow=P+1,ncol=nTimes)

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

  #Dimension vectors for Punt prey effect, phi, and Punt recruits, R_Punt.
  if("phi" %in% slotNames(Species)){
    Species@phi<-as.numeric(rep(NA,nTimes))
  }
  if("R_Punt" %in% slotNames(Species)){
    Species@R_Punt<-as.numeric(rep(NA,nTimes))
  }
  if(is(Species,"Predator")){
    Species@RperS<-as.numeric(rep(NA,nTimes))
    Species@NperK<-as.numeric(rep(NA,nTimes))
  }

  #Create and dimension the vectors for saving prey availability data
  if(is(Species,"Predator")){
    for(i in 1:length(preynames)){
      preyxname<-preynames[[i]]
      #get existing prey_availability list
      splistx<-Species@prey_availability[[preyxname]]
      #Add two new numeric vectors to it, and name them B_avail and B_agepref
      splistx<-c(splistx,list(as.numeric(rep(NA,nTimes))),list(as.numeric(rep(NA,nTimes))))
      names(splistx)[length(splistx)]<-"B_avail"
      names(splistx)[length(splistx) - 1]<-"B_agepref"
      Species@prey_availability[[preyxname]]<-splistx
    }
  }

  #Generate the recruitment residuals vector (epsilon) for species that have stochastic recruitment error,
  #by random draws from a distribution.  epsilon is a vector nYears or (nSteps x nYears) long.
  #Note: we have to calculate anchovy recruitment error on the fly later, because it depends on R.
  if(Species@R_error$R_stochastic){
    tmp<-Species@R_error
    if((!exists('resample_residuals',where=tmp)) || (Species@R_error$resample_residuals==0)){
      Species<-calc_R_error(Species)
    }
  }#if(R_stochastic)

  #Dimension the fishing mortality and catch vectors for species that are fished.
  if(Species@fished){
    if(length(Species@Sv)==0){
      stop(paste(Species@name,"is fished, but fishing selectivity Sv is missing from input_params.csv"))
    }
    Species@Fv<-numeric(nTimes) #Initialize as zeros (NAs would cause problems)
    Species@Cy<-numeric(nYears)
    if(Species@uses_timestep){
      Species@Ct<-numeric(nTimes)
    }
  }

  #print("Matrices initialized.")
  #Return the Species object
  Species
}

#Solve for age-0 survival
#' Solve for age-0 survival
#'
#' Solves for age-0 survival of predators iteratively, by balancing reproduction
#' and mortality (ignores prey effect, so this is unfished prey condition).
#'
#' @param predator A Predator object.
#'
#' @return A predator with modified Mv slot (set age-0 survival, Mv[1]).
#' @export
#'
#' @keywords internal
#' @family initialization
solve_age0_survival<-function(predator){
  SJmin<-0
  SJmax<-10
  S1<-exp(-(predator@Mv[2])) #Get age-1 survival from mortality rate
  S0<-S1 #age-0 survival, initially set same as age-1
  am<-predator@am #age at maturity
  target<-1
    for(i in 1:50){
      mult<-(SJmin + SJmax)/2
      SJ<-S0 * mult
      predicted<-S1 + SJ * S1^(am-1)
      if(predicted < target){
        SJmin<-mult
      } else{
        SJmax<-mult
      }
    }
  S0<-S0 * mult
  #Set age-0 mortality rate (log scale) = -ln(S) in the predator object.
  predator@Mv[1]<--log(S0)
  predator
}

#' Solve for the scaling parameter for fecundity.
#'
#' Finds the value of a scaling parameter for predator reproduction (upper-case
#' Phi in Punt et al. 2016, equation 9 line 1).  Note: you must solve for age-0
#' survival before calling this function.
#'
#' @param predator A Predator object.
#'
#' @return A modified Predator object with bigPhi set in the prey_effect slot.
#' @export
#'
#' @keywords internal
solve_fecundity<-function(predator){
  SJmin<-0
  SJmax<-10
  S1<-exp(-(predator@Mv[2])) #age-1 survival
  S0<-exp(-(predator@Mv[1])) #age-0 survival
  am<-predator@am #age at maturity
  Fmax<-predator@Fmax
  target<-Fmax^am #(why? )
  for(i in 1:50){
    mult<-(SJmin + SJmax)/2
    SJ<-S0 * mult
    #I'm guessing this is something like how you calculate lambda for a Leslie matrix
    predicted<-Fmax^(am-1) * S1 + SJ * S1^(am-1)
    if(predicted < target){
      SJmin<-mult
    } else{
      SJmax<-mult
    }
  }
  bigPhi<-mult
  predator@prey_effect$bigPhi<-bigPhi
  predator
}

#' Create environmental driver, G
#'
#' Generates a unique time series for each prey species to simulate an
#' environment that switches between high (good reproduction) and low (poor
#' reproduction) states.  The function samples from a list of period lengths for
#' low and high periods that is read in from an internal file
#' (envir_periods.csv).  The period lengths and reproduction values for the
#' periods are balanced to match observed patterns.
#'
#' @param Species A Prey species object.
#' @param simParams A simulationParameters object
#'
#' @return A modified Prey object, with the G slot set.
#' @export
#'
#' @keywords internal
calc_G<-function(Species,simParams){
  spname<-Species@name
  nYears<-simParams@nYears
  periods<-simParams@envir_periods
  GG<-numeric(nYears)
  #Set G to all zeros if useG is false; otherwise proceed
  if(Species@recruitment$useG==0){
    Species@G<-numeric(nYears)
  } else{
    low_periods<-periods$period[periods$species==spname & periods$state=="low"]
    high_periods<-periods$period[periods$species==spname & periods$state=="high"]
    #Constants (multipliers for height (hm), and the random number (rm) that determines variance)
    #Note that high periods are less variable than low periods, and that Other and Anchovy are equal
    hm<-list(Sardine=list(low= -1.13,high = 2.31),Other_forage=list(low= -1.13,high = 2.31),Anchovy=list(low= -0.92,high = 1.47))
    rm<-list(Sardine=list(low= 0.932,high = 0.422),Other_forage=list(low= 0.932,high = 0.422),Anchovy=list(low= 0.35,high = 0.18))
    if(is(Species,"Sardine") || is(Species,"Other_forage")){
      state<-1
      year<-1
      while(year <= nYears){
        if(state==1){
          period<-sample(low_periods,1)
          rand<-rnorm(1,0,1) * rm[[spname]]$low
          height<- hm[[spname]]$low * (1 + rand)
          state<-2
        } else{
          period<-sample(high_periods,1)
          rand<-rnorm(1,0,1) * rm[[spname]]$high
          height<- hm[[spname]]$high * (1 + rand)
          state<-1
        }
        if(is(Species,"Sardine")){
          for(y in year:as.integer(min(nYears,year + period - 1))){
            GG[y]<-height
          }
        }
        if(is(Species,"Other_forage")){
          period<-100 #the sample function doesn't work as expected for a single number
          sigma<-Species@R_error$sigma #not certain this is the right parameter -- see CalcG in Punt code
          for(y in year:as.integer(min(nYears,year + period - 1))){
            GG[y]<-(height -.72)/1.72 * sigma - sigma^2/2
          }#for y in period
        }
        year<-year + period
      }#while year < nYears
    }#if species is Sardine or Other_forage
    if(is(Species,"Anchovy")){
      state<-1
      year<-1
      while(year <= nYears){
        if(state==1){
          period<-sample(low_periods,1)
          rand<-rnorm(1,0,1) * rm[[spname]]$low
          height<- hm[[spname]]$low * (1 - 1 + rand)
          state<-2
        } else{
          period<-sample(high_periods,1)
          rand<-rnorm(1,0,1) * rm[[spname]]$high
          height<- hm[[spname]]$high * (1 - 1 + rand)
          state<-1
        }
        for(y in year:as.integer(min(nYears,year + period - 1))){
          GG[y]<-height
        }
        year<-year + period
      }#while year < nYears
    }#if species is Anchovy

    #Set G in the Species object
    Species@G<-GG
  }#useG !=0

  Species
}#calc_G


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
            species<-calc_G(species,params)
            species<-initialize_N(species,params) #Note: this must come after get_params
            species
          })

#' The method for predators also reads in a predator's diet preferences (for
#' prey species and the age-classes of those prey).  Does not read in G.
#' @describeIn set_species Also reads in diet preferences.
setMethod("set_species",
          signature("Predator","simulationParams"),
          function(species,params,preynames) {
            species<-get_params(species,params)
            species<-get_prefs(species,params)
            species<-initialize_N(species,params,preynames) #Note: must come after get_params
            species
          })

#' The method for anchovy also reads in a set of recruitment residuals.
#' @describeIn set_species Also reads in recruitment residuals.
setMethod("set_species",
          signature("Anchovy","simulationParams"),
          function(species,params,preynames) {
            species<-get_params(species,params)
            species<-calc_G(species,params)
            species<-load_residuals(species,params)
            species<-initialize_N(species,params,preynames) #Note: must come after get_params
            species
          })

