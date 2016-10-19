#' Calculate total mortality rate
#'
#' A utility function used to calculate total mortality rate by age in species
#' that are tracked by timestep.
#'
#' @param M Natural mortality rate (scalar).
#' @param L Number of timesteps per year (scalar).
#' @param Ft Total fishing mortality rate (scalar).
#' @param Sv Fishing selectivity (vector).
#'
#' @return A vector of mortality rate (log scale) by age.
#' @export
#' @keywords internal
calc_mortality_t<-function(M,L,Ft,Sv){
  zv<-M/L + Ft * Sv
  zv
}

#
#' Calculate total mortality rate
#'
#' A utility function used to calculate total mortality rate by age in species
#' that are tracked by year.
#'
#' @param M Natural mortality rate (scalar).
#' @param Ft Total fishing mortality rate (scalar).
#' @param Sv Fishing selectivity (vector).
#'
#' @return A vector of mortality rate (log scale) by age.
#' @export
#' @keywords internal
calc_mortality_y<-function(M,Ft,Sv){
  zv<-M + Ft * Sv
  zv
}

#' Put survival rates into a Leslie matrix
#'
#' A utility function used in the calculation of numbers-at-age at time t.
#' Recruits are calculated separately, so fecundities are all set to 0 here.
#' Survival rates are in off-diagonal rows (survival from age 1 to age 2 is in
#' leslie[1,2], etc.). Note: age classes are 0:P, but vector and matrix row/col
#' indexes are 1-based, and go from 1:(P+1=nAges).
#' @param P A scalar: Maximum age (plus-group)
#' @param Zv A vector: Total mortality rate (log scale) at age.
#' @return A leslie matrix
#' @export
#' @keywords internal
set_lesliematrix<-function(P,Zv){
  nAges<-P+1
  leslie<-matrix(data=0,nrow=nAges,ncol=nAges)
  #Put survival rates in off-diagonal rows
  sv<-exp(-Zv)
  for (i in 1:(nAges-1)){
    leslie[i+1,i]<-sv[i]
  }
  leslie[nAges,nAges]<-sv[nAges] #Plus-group (age P) survival has this extra component
  leslie
}

#'Calculate Numbers at age for year y, for species tracked by year.
#'
#' Calculates numbers at age in year y, Nmat[,y], from Nmat[,y-1]. Sets the
#' following values in the Species object: Nmat[,y]; the previous year's total
#' mortality rate (Zv[y-1], log scale); the current year's fishing mortality
#' rate (Fv[y], log scale, fished species only) and the current year age 1+ biomass, B1P[y]
#' (prey species only).
#'
#' @param Species A Species object.
#' @param y The current year.
#'
#' @return A modified Species object; see Description.
#' @export
#' @keywords internal
#' @family Basic_dynamics
calc_Ny<-function(Species,y){
  fished<-Species@fished
  P<-Species@P
  M<-Species@M
  Nvprev<-Species@Nmat[,y-1] #Last year's N

  #Calculate total mortality rate for year y-1.
  if(fished){
    Sv<-Species@Sv
    Ftprev<-Species@Fv[y-1] #Last year's fishing mortality rate
    Zvprev<-calc_mortality_y(M,Ftprev,Sv) #last year's total mortality rate
  } else{
    #Set Z = natural mortality for all age classes
    Zvprev<-numeric(P+1) #a 1-d total mortality vector for year y-1
    Zvprev[]<-M
  }

  #For year 1, we assume only natural mortality (this function is first called in y=2, so the previous year's mortality is M.)
  if(y==2){
    Zvprev<-numeric(P+1)
    Zvprev[]<-M
  }
  Species@Zv[,y-1]<-Zvprev #TODO: change Species@Zv to a vector once running

  #Calculate N for this year
  leslie<-set_lesliematrix(P,Zvprev)
  Nv<-leslie %*% Nvprev #multiply numbers by survival for 0+ fish
  #Put Nv into the Species object
  Species@Nmat[,y]<-Nv #Numbers

  #If it is a prey species, calculate biomass of age-1+ animals and set that in the Species object
  if(is(Species,"Prey")){
    wv<-Species@wv
    Species@B1P[y]<-calc_B1P(Nv,wv)
  }

  #Calculate this year's total fishing mortality rate (a scalar) from the fishing_rule model.
  #It is done after setting Nv, because the fishing rule may depend on this year's N.
  #Set it in the Species object.
  if(fished){
    Ft<-calc_F_y(Species,y)
    Species@Fv[y]<-Ft #fishing mortality for year y
  }

  #Return the species object
  Species
}

###Utility functions for species with timesteps

#' Calculate Numbers at age at time t, for species tracked by timestep.
#'
#' Calculates numbers at age at time t, Nmat[,t], from Nmat[,t-1] and sets it in the
#' Species object.  Reproduction happens in timestep 1 each year, and during other
#' timesteps, only mortality is applied.  The function sets the following values
#' in the Species object: Nmat[,t]; the previous timestep's total mortality rate
#' (Zv[t-1], log scale); the current fishing mortality rate (Fv[t], log
#' scale, fished species only); and the current age 1+ biomass, B1P[t] (prey species
#' only).

#' @param Species A Species object.
#' @param ts Global timestep (ts = t * y)
#' @param t Within-year timestep
#' @param L Number of within-year timesteps per year.
#'
#' @return A modified Species object (see Description)
#' @export
#' @family Basic_dynamics
calc_Nt<-function(Species,ts,t,L) {
  #Get values we will need
  #browser()
  P<-Species@P
  nAges<-P+1
  Ntprev<-Species@Nmat[,ts-1] #Numbers-at-age from the previous timestep
  Nt<-as.numeric(rep(NA,length(Ntprev))) #Numbers at age from this timestep (NA, to start with)
  M<-Species@M
  fished<-Species@fished

  #Calculate mortality rate-at-age for the previous timestep
  if(fished){
    Sv<-Species@Sv #selectivity-at-age
    Ftprev<-Species@Fv[ts-1] #Retrieve last year's fishing mortality rate
    Zvprev<-calc_mortality_t(M,L,Ftprev,Sv)
  } else{
    #Set the mortality for each age to natural mortality divided by the number of steps per year.
    Zvprev<-numeric(nAges)
    Zvprev[]<-M/L #mortality/nSteps
  }
  Species@Zv[,ts-1]<-Zvprev #TODO: once running, just keep a vector of total mortality (not by-age)

  #Do the update (at t=1 we use the Leslie matrix; for later steps fish stay in the same age class)
  if(t==1){
    leslie<-set_lesliematrix(P,Zvprev) #Set up the Leslie matrix with the survivals
    Nt<-leslie %*% Ntprev
  }  else {
    Nt[1]<-Ntprev[1] #Mortality = 0 for recruits across all timesteps within a year
    Nt[2:nAges]<-Ntprev[2:nAges] * exp(-Zvprev[2:nAges]) #Fish stay in the same age classes but suffer mortality
  }

  #Set Nt in the Species object
  Species@Nmat[,ts]<-Nt

  #If it is a prey species, calculate biomass of age-1+ animals and set that in the Species object
  if(is(Species,"Prey")){
    wv<-Species@wv
    Species@B1P[ts]<-calc_B1P(Nt,wv)
  }

  #Calculate this year's total fishing mortality rate (a scalar) from the fishing_rule model.
  #It is done after setting Nt, because the fishing rule may depend on this year's N.
  #Set it in the Species object.
  if(fished){
    #Get the total fishing mortality rate (a scalar) for this timestep from the fishing_rule model.
    Ft<-calc_F_t(Species,ts,t)
    #Set Ft in the Species object
    Species@Fv[ts]<-Ft #fishing mortality for this timestep
  }

  #Return the Species object
  Species
}

#Spawner-recruit functions

###Set up generic and dispatch for alternative versions of calc_recruits.
#NOTE: if we start making Anchovy recruitment dependent on sardines, we'll need to add a new method for Anchovy
#Make calc_SR an S4 generic method

#' Calculate recruits.
#'
#' @param species. A Species object.
#' @param preylist. A list of Prey objects (only used for Predators).
#' @param ... Other arguments (from child methods).
#'
#' @return A child method, depending on signature.
#' @export
#' @keywords internal
#' @family recruitment
#' @rdname calc_recruits
setGeneric( "calc_recruits" , def=function(species,preylist,...) {
  standardGeneric( "calc_recruits" )
})

#' The default method calculates recruits in a Prey object.
#'
#' @describeIn calc_recruits Prey method.
setMethod("calc_recruits", signature("Prey"),function (species,preylist,ts,y,...) {
  calc_preyRecruits(species,preylist,ts,y,...)
})

#' The Other_forage method calculates recruits for Other_forage species.
#'
#' @describeIn calc_recruits Prey method.
setMethod("calc_recruits", signature("Other_forage"),function (species,y,...) {
  calc_other_R(species,y,...)
})

#' The Predator method calculates recruits for Predators.
#'
#' @describeIn calc_recruits Prey method.
setMethod("calc_recruits", signature("Predator","list"),function (species,preylist,y,ts,nSteps,...) {
  calc_predatorRecruits(species,preylist,y,ts,nSteps,...)
})

###Utility functions.

#' Calculate spawning biomass (tons)
#'
#' A utility function used in recruitment calculations, in some species.
#'
#' @param wv Weight at age (in kg)
#' @param Nv Numbers at age
#' @param mv maturity at age
#'
#' @return A scalar: biomass in tons
#' @export
#' @keywords internal
#' @family recruitment
calc_B_sp<-function(wv,Nv,mv){
  Bs<-sum(wv * Nv * mv)/1000
  Bs
}

#' Calculate age 1+ biomass
#'
#' Used for calculating prey biomass of a species in one timestep or year, where
#' prey size matters to the predator but prey maturity does not).
#'
#' @param Nv Numbers-at-age
#' @param wv Weight-at-age
#'
#' @return A scalar of total spawning biomass of fish aged 1+ (i.e., excludes age-0), in tons ('000 kg).
#' @export
#' @keywords internal
#' @family recruitment
calc_B1P<-function(Nv,wv){
  n_age<-length(Nv)
  B1P<-sum(wv[2:n_age] * Nv[2:n_age])/1000
  B1P
}


#' Calculate stochastic recruitment errors
#'
#' Calculate stochastic recruitment errors (either autocorrelated or not) for
#' all nYears for a Species.  Sigma and rho vary by species.  In this version,
#' the errors are drawn from a log-normal distribution
#' @param Species
#'
#' @return A modified Species object, with recruitment errors (1 x nTimes)
#'   vector stored in a slot called lnSR_err
#' @export
#' @family recruitment
#' @keywords internal
calc_R_error<-function(Species){
  sigma<-Species@recruitment$sigma
  rho<-Species@recruitment$rho
  epsilon<-Species@lnSR_err

  epsilon[1]<-Species@recruitment$epsilon_baseyear
  for(t in 2:length(epsilon)){
    epsR<-rnorm(1,mean=0,sd=1) * sigma
    epsilon[t]<-rho * epsilon[t-1] + sqrt((1-rho)^2) * epsR
  }

  #Set the value and return the object
  Species@lnSR_err<-epsilon
  Species
}

#Resample recruitment errors with replacement from values that were read in from the params.csv file. Put the resampled errors in the lnSR_err slot of the species object.
resample_R_error<-function(Species){
  epsilon<-numeric(dim(Species@Nmat)[2]) #a new vector to hold errors
  epsilon[1]<-Species@recruitment$epsilon_baseyear
  epsilon[2:length(epsilon)]<-sample(Species@observed_SR_err,length(epsilon)-1, replace=T)
  Species@lnSR_err<-epsilon
  Species
}

#' Calculate Prey recruits
#'
#' Calculate prey recruits at one timestep or year. As implemented, it
#' calculates the *number* of recruits from the *biomass* of spawners.  In
#' brief: 1) it calculates spawning biomass, 2) in the case of anchovy, it
#' calculates the biomass of sardines and adjusts the recruitment model
#' parameters in response, 3) it chooses a value of log recruitment error for
#' the year or timestep, and 4) calculates recruitment, using a Ricker model
#' with an environmental driver.
#'
#' @param Prey A prey object.
#' @param preylist A list of prey.  Used when the reproduction of one prey
#'   species depends on another (in this implementation, anchovy recruitment
#'   depends on the biomass of sardines, since sardines eat anchovies).
#' @param ts The current global timestep (ts = t + y * t; where t is the intra-year timestep)
#' @param y The current year (integer).
#' @return A modified Prey object, with recruits for the year set in the
#'   population matrix Nmat[1,t] (t may be a year or timestep).  Also sets R in
#'   R[t] and spawning biomass in B_sp[t], for convenience in graphing.
#' @export
#' @keywords internal
#' @family recruitment
calc_preyRecruits<-function(Prey,preylist,ts,y){
  if(Prey@uses_timestep){t<-ts} else {t<-y}
  alpha<-Prey@recruitment$alpha
  beta<-Prey@recruitment$beta
  #browser()
  gamma<-Prey@recruitment$gamma
  sigma2<-(Prey@recruitment$sigma)^2
  R_thresh<-Prey@recruitment$R_thresh
  G<-0 #Placeholder.  May want to generate environmental regimes with (-1 < G < +1).
  #Calculate spawning biomass
  Nv<-Prey@Nmat[,t]
  wv<-Prey@wv
  mv<-Prey@mv
  B_sp<-calc_B_sp(wv,Nv,mv)

  #Calculate recruits
  if(B_sp < R_thresh){
    R<-0
  } else{
    #If this is anchovy and there is a sardine effect on recruitment, get the corresponding alpha and beta
    if(is(Prey,"Anchovy") && Prey@recruitment$sardine_effect){
      #Get age-1+ biomass of sardine
      sardine_B1P<-preylist$sardine@B1P[t]
      #Set alpha and beta
      if(sardine_B1P < Prey@recruitment$sardine_thresh){
        alpha<-Prey@recruitment$alpha1
        beta<-Prey@recruitment$beta1
      } else{
        alpha<-Prey@recruitment$alpha2
        beta<-Prey@recruitment$beta2
      }
    }#if calculating sardine effect

    #Draw a spawner-recruit residual
    #epsilon will equal 0 if R_stochastic is false, and may be autocorrelated depending on the value of rho
    #(see calc_R_error, which is called at initialization).
    if(Prey@recruitment$R_stochastic){
      epsilon<-Prey@lnSR_err[t]
    } else{
      epsilon<-0
      sigma2<-0
    }
    #print(paste(t,Prey@name,alpha,beta))

    #Ricker with environmental driver G
    R<-calc_Ricker(alpha,beta,B_sp) #deterministic component
    R<-R * exp((gamma * G) + epsilon - sigma2/2) #variable component: envt plus error
    R
  }
  #Put the recruitment into the population matrix at the appropriate year or timestep;
  #also fill R and B_sp, and return the object
  Prey@B_sp[t]<-B_sp #sp
  Prey@Nmat[1,t]<-R
  Prey@R[t]<-R
  Prey
}

#' Calculate recruits per spawner (Ricker method)
#'
#' Note expansion by 1000 for this case.
#'
#' @param alpha Ricker alpha.
#' @param beta Ricker beta.
#' @param B_sp Spawning biomass (in tons)
#' @return Recruits (a scalar).
#' @export
#' @keywords internal
#' @family recruitment
calc_Ricker<-function(alpha,beta,B_sp){
  #R<-(1000 * alpha * B_sp * exp((-beta * B_sp) + (gamma * G) + epsilon - sigma2/2))
  R<-(1000 * alpha * B_sp * exp((-beta * B_sp)))
  R
}

#' Calculate recruits for Other_forage species
#'
#' As implemented, recruitment is constant at 10,000; however, log-scale
#' recruitment error is added as an autocorrelated random draw from a log-normal
#' distribution.
#'
#' @param Other An Other_forage species object.
#' @param y The current year (integer).
#' @return A modified Other_forage object, with recruitment set both in the population matrix (Nmat[1,y]) and in the R slot (R[y]).
#' @export
#' @keywords internal
#' @family recruitment
calc_other_R<-function(Other,y){
  epsilon<-Other@lnSR_err[y] #may be autocorrelated, depending on the value of rho
  sigma<-Other@recruitment$sigma
  #Get the fixed component of recruitment
  R_oth<-Other@recruitment$R_base
  #Add stochastic error if R_stochastic is true
  if(Other@recruitment$R_stochastic){
    R_oth<-R_oth * exp(epsilon -((sigma)^2)/2) #This year's recruitment
  }
  Other@Nmat[1,y]<-R_oth
  Other@R[y]<-R_oth
  #Set spawning biomass (only used for plotting)
  Other@B_sp[y]<-calc_B_sp(Other@wv,Other@Nmat[,y],Other@mv)
  Other
}

##Predators
###Utility functions
#These are for the various components of recruitment (prey biomass-at-age, predator preferences for species and age-classes, and spatial availability of the prey to the predator).
#Calculate the number of mature predators:
#Arguments: Current-year numbers at age (Nv) and proportion mature at age (mv) (both vectors).
#Returns a scalar (total number of mature predators).
calc_Nm<-function(Nv,mv){
  Nm<-Nv %*% mv
  Nm
}

#' Calculate prey Availability parameter
#'
#' The availability function is a surrogate for the spatial availability of prey
#' to predators (i.e., it answers the question "How much of the prey population
#' can the predators access?")
#' @param ap The p parameter of a Beverton-Holt function
#' @param ac The c parameter of a Beverton-Holt function
#' @param Bt The biomass of prey at time t
#' @return A value that represents some fraction of Bt.
#' @export
#' @keywords internal
#' @family recruitment
 calc_availability<-function(ap,ac,Bt){
   A<-Bt/((1/ap) + (1/ac) * Bt)
   A
 }

#' Calculate the age-preferred biomass of a prey species available to a predator
#'
#' A utility function that calculates the biomass (tons) of prey species X that
#' is available to predator species Y, as:  (biomass-at-age of prey) x (predator
#' preferences for prey age classes). This function returns an intermediate
#' product in the calculation of total available prey biomass. Note: this is
#' akin to Punt et al's use of biomass of age 1+ fish, but allows the predator
#' to eat some of the age-0 fish as well if their age-preferences indicate it.
#' @param predator A predator object.
#' @param prey A prey object.
#' @param t Time t.
#'
#' @return A biomass value.
#' @export
#' @keywords internal
#' @family recruitment
calc_Bprey_agepreferred<-function(predator,prey,t){
 age_prefs<-as.numeric(unlist(predator@age_prefs[prey@name])) #Predator preferences for this species' age classes: a vector of length nAges (prey ages)
 Nv<-prey@Nmat[,t] #prey numbers-at-age at this timestep
 wv<-prey@wv #prey weight-at-age
 B_preferred<-sum(Nv * wv * age_prefs)
 B_preferred<-B_preferred/1000 #convert to tons, to be consistent with B1P
 B_preferred
}

#' Calculate Dy.
#'
#' Calculates the scaled, total biomass of all prey, Dy, adjusted for species-
#' and age-class preferences and spatial availability, for a predator X.  Dy is
#' set so that Dy = 1 and mean diet proportions of prey species in the diet of
#' predator X = observed diet proportions when the prey are unfished.
#'
#' @param predator A Predator object.
#' @param Bprey A vector of available prey biomass for each species.
#' @param preylist A list of Prey objects.
#' @return A value for Dy (scalar).
#' @export
#' @keywords internal
#' @family recruitment
calc_Dy<-function(predator,Bprey, preylist){
 OP<-predator@sp_prefs$Other_prey #this predator's preference for Other_prey
 Dy<-OP
 for(i in 1: length(preylist)){
   preyname<-preylist[[i]]@name
   omega<-as.numeric(unlist(predator@omegas[preyname])) #get the omega
   Dy<-Dy + as.numeric(unlist(Bprey[preyname])) * omega
 }
 Dy
}

#' Calculate impact of prey biomass on predator reproduction
#'
#' Calculates the impact of available prey biomass on the maximum offspring per
#' predator, for predator species X.  The function has two parameters: K, which
#' controls how steep the logistic inflection is, and Dhat, which controls the
#' location of the inflection point on the X axis.  Those parameters and Mo, the
#' maximum offspring per spawner, are stored in the Predator object.
#'
#' @param predator A Predator object
#' @param Dy The total biomass of prey available to the predator.
#'
#' @return The Beverton-Holt p parameter value (a scalar), scaled to represent some
#'   fraction of the maximum offspring per predator.
#' @export
#' @keywords internal
#' @family recruitment
calc_p<-function(predator,Dy){
 k<-predator@prey_effect$k
 Dhat<-predator@prey_effect$Dhat
 Mo<-predator@Mo
 f<-1/(1 + exp(-k * (Dy-Dhat)))
 p<-Mo * f
 p
}

#Arguments: p, which defines max offspring per predator at low densities; Nm (Number of mature predators); and c (carrying capacity parameter).
#Returns a scalar (total recruits for predator species X)
#' Calculate recruits, using a Beverton-Holt model
#'
#' @param Nm Number of mature predators.
#' @param p A parameter that defines maximum offspring per predator at low densities.
#' @param c The parameter that defines carrying capacity.
#' @return Number of recruits (scalar).
#' @export
#' @keywords internal
#' @family recruitment
calc_Bev_Holt<-function(Nm,p,c){
 recruits<-Nm/((1/p) + (1/c) * Nm)
 recruits
}

#' Calculate Predator Recruits
#'
#' This function calls several utility functions to calculate predator recruits.
#' In brief: 1) it calculates the number of breeding predator adults; 2) it
#' calculates prey biomass; 3) it calculates the proportion of the prey biomass
#' that is available to the predator, given a) the predator's diet preferences
#' for prey species, b) the predator's preferences for certain age classes of
#' those species, and c) the predator's ability to access the prey (spatial
#' availability).
#'
#' @param Predator A Predator object.
#' @param preylist A list containing the Prey objects.
#' @param y The current year.
#' @param ts The current global timestep.
#' @param nSteps The number of timesteps per year.
#' @return A modified Predator object, with number of mature predators (Nm),
#'   predator recruits (Nmat[,t]) both set.  As a convience for later graphing,
#'   the number of recruits is also set in the predator's @R slot, and Dy (total,
#'   scaled prey biomass) and BHp (Beverton-Holt maximum growth parameter) are also set.
#' @export
#' @keywords internal
#' @family recruitment
calc_predatorRecruits<-function(Predator,preylist,y,ts,nSteps){
 #Set the time to year or timestep:
 #NOTE that each predator and/or prey species may be on a different time system.
 if(Predator@uses_timestep){t<-ts} else{t<-y}
 #Get variables
 nAges<-Predator@P + 1
 Nv<-Predator@Nmat[,t] #Numbers at age, time t (either year or timestep)
 mv<-Predator@mv #Proportion mature at age
 BHc<-Predator@recruitment$BHc #Beverton-Holt carrying capacity.  Steepness (BHp) is determined below by prey effect.
 R_thresh<-Predator@recruitment$R_thresh #recruitment threshold (minimum number of mature adults)

 #Calculate the number of mature predators
 Nm<-calc_Nm(Nv,mv)

 #Calculate prey biomass
 if(Nm < R_thresh){
   R<-0
 } else{
   Bprey<-list()
   #Calculate the available, age-preference-adjusted biomass of each prey species
   for(i in 1:length(preylist)){
     ap<-preylist[[i]]@availability$ap #the ap parameter for availabililty
     ac<-preylist[[i]]@availability$ac

     #Get prey biomass from the first timestep of year y (or just year y if not using timesteps)
     if(preylist[[i]]@uses_timestep){
       tt<-(((y-1) * nSteps) + 1)
     } else {
       tt<-y
     }
     #First get age-preferred biomass
     B_agepref<-calc_Bprey_agepreferred(Predator,preylist[[i]],tt)
     #Now adjust it for spatial availability
     B_avail<-calc_availability(ap,ac,B_agepref)
     #Add the biomass to a list and name it
     Bprey<-c(Bprey,list(B_avail))
     names(Bprey)[length(Bprey)]<-preylist[[i]]@name
   }

   #Calculate Dy, the total scaled available prey biomass
   Dy<-calc_Dy(Predator,Bprey,preylist)

   #Calculate the impact of the available biomass on the predator reproductive rate (affects the maximum
   #rate of increase, p, in the Beverton-Holt equation below)
   BHp<-calc_p(Predator,Dy)
   #browser()
   #Finally, calculate the recruits
   R<-calc_Bev_Holt(Nm,BHp,BHc)

   #Add recruitment error, if R_stochastic is true
   if(Predator@recruitment$R_stochastic){
     sigma2<-(Predator@recruitment$sigma)^2
     epsilon<-Predator@lnSR_err[t]
     R<-R * exp(epsilon - sigma2/2)
   }
 }#if Number of mature adults Nm is above recruitment threshold R_thresh

 #And put them into the predator population matrix
 Predator@Dy[t]<-Dy
 Predator@BHp[t]<-BHp
 Predator@Nm[t]<-Nm
 Predator@Nmat[1,t]<-R
 Predator@R[t]<-R #Set the R value (for convenience)
 Predator #return the Predator
}

#' Calculate fishing mortality according to a fishing rule
#'
#'There are two versions: one for species tracked by year, the other for species
#'tracked by timestep. Note: it would be easy to make this a generic and set up
#'alternative models for each species. This version is for species that are
#'tracked by timestep.  ts= t + (y * t).
#'
#' @param Species A Species object.
#' @param ts The global timestep ()
#' @param t The intra-year timestep.
#'
#' @return A scalar, Ft, the total fishing mortality rate (log scale) for timestep t.
#' @export
#' @family fishing
#' @keywords internal
calc_F_t<-function(Species,ts,t){
 Ft<-Species@F_baseline
 gv<-Species@gv #proportion fished per time period
 Bmin<-Species@fishing_rule$Bmin #Minimum biomass of age 1+ fish, below which fishing is stopped.

 #calculate the observed biomass of age 1+ fish
 Nv<-Species@Nmat[,ts]
 wv<-Species@wv
 B1P<-calc_B1P(Nv,wv)

 #The rule: if biomass < Bmin, set catch to zero
 if(B1P < Bmin){
   Ft<-0
 }
 Ft<-Ft * gv[t]
 Ft
}

#' Calculate fishing mortality according to a fishing rule
#'
#'This version is for species tracked by year.  Note: it would be easy to make
#'this a generic and set up alternative models for each species.
#' @param Species A Species object.
#' @param y The current year.
#' @param t The intra-year timestep.
#' @return A scalar, Fy: the total fishing mortality rate year (log scale) for year y.
#' @export
#' @family fishing
#' @keywords internal
calc_F_y<-function(Species,y){

 #Set fishing mortality at baseline (default)
 Fy<-Species@F_baseline
 Bmin<-Species@fishing_rule$Bmin #Minimum biomass of age 1+ fish, below which fishing is stopped.

 #calculate the observed biomass of age 1+ fish
 Nv<-Species@Nmat[,y]
 wv<-Species@wv
 B1P<-calc_B1P(Nv,wv)

 #The rule: if biomass < Bmin, set catch to zero
 if(B1P < Bmin){
   Fy<-0
 }
 Fy
}

#' Calculate catch (tons) for a timestep or year
#'
#' @param Species A Species object
#' @param t time period (may be year or global timestep)
#'
#' @return A Species object with catch set if fished=TRUE, or unset if fished=FALSE
#' @export
#' @family fishing
#' @keywords internal
calc_catch<-function(Species,t){
 if(Species@fished){
   Ft<-Species@Fv[t] #total fishing mortality Ft for time period t
   Nv<-Species@Nmat[,t] #numbers at age
   wv<-Species@wv #weight at age
   Sv<-Species@Sv #selectivity at age

   Fv<-Ft * Sv #Distribute fishing mortality rates by age (via selectivity)
   C<-sum(Nv * (1-exp(-Fv)) * wv)/1000
   #Put catch into the Species object and return it
   Species@Cv[t]<-C
 }
 Species
}

#' Do a calibration run.
#'
#' Do an initial run for prey species without fishing, to calculate parameters
#' in the Predator objects (omegas) that scale the unfished biomass of each prey species.
#'
#' @param predator A list of Predator objects
#' @param prey A list of Prey objects
#' @param nYears Number of years for the simulation
#' @param nSteps Number of timesteps per year.
#'
#' @return A modified list of predators in which the @omegas slot, which holds
#'   scaling parameters for the ratio of diet preference to unfished biomass for
#'   each prey species, has been set for each predator.  NOTE: the prey list is
#'   unchanged.
#' @export
#' @keywords internal
#' @family initialization
calc_omegas<-function(predator,prey,nYears,nSteps){
 #Run the main model without fishing or stochastic recruitment
 for(i in 1:length(prey)){
   prey[[i]]@fished<-FALSE
   prey[[i]]@recruitment$R_stochastic<-FALSE
 }
 #Simulate prey only
 for (y in 2:nYears){
   #Do prey calculations
   for(i in 1:length(prey)){
     if(!prey[[i]]@uses_timestep){
       prey[[i]]<-calc_Ny(prey[[i]],y)
       prey[[i]]<-calc_recruits(prey[[i]],prey,y)
       prey[[i]]<-calc_catch(prey[[i]],y)
     } else{
       for(t in 1:nSteps){
         ts<-t + ((y-1) * nSteps) #The step number, as a constantly increasing series.  Remember, year 1 has been initialized already.
         prey[[i]]<-calc_Nt(prey[[i]],ts,t,nSteps) #Calculate numbers
         if(t==1){
           prey[[i]]<-calc_recruits(species=prey[[i]],preylist=prey,ts=ts,y=y) #Calculate recruits
         }
         prey[[i]]<-calc_catch(prey[[i]],ts) #calculate catch
       }#for t 1:nSteps
     }#uses_timestep
   }#for preyx in prey
 }#for y in years

 #Set omega=(species preference/unfished age-preferred biomass) for each predator/prey combination.
  for(i in 1:length(predator)){
    omegas<-list()
    for(j in 1:length(prey)){
      preyname<-prey[[j]]@name
      spx_pref<-as.numeric(unlist(predator[[i]]@sp_prefs[preyname])) #the preference of predator i for prey species j
      #Calculate mean age-preferred biomass of prey j
      nmat<-prey[[j]]@Nmat #prey N matrix
      wv<-prey[[j]]@wv
      #get the predator's age prefs for this species
      ageprefsx<-as.numeric(unlist(predator[[i]]@age_prefs[preyname]))
      #Get total age-preferred biomass (a vector nTimes long)
      Btotals<-colSums(nmat * wv * ageprefsx)/1000 #in tons, to match B1P
      #Adjust it for availability
      ap<-prey[[j]]@availability$ap
      ac<-prey[[j]]@availability$ac
      Bavail<-sapply(Btotals,function(b) calc_availability(ap,ac,b))
      B_ap_unfished<-mean(Bavail)
      #Calculate omega = (pred preference for prey x) / (age-preferred biomass of prey x)
      omega<-spx_pref/B_ap_unfished
      #print(paste(predator[[i]]@name,prey[[j]]@name,B_ap_unfished,spx_pref,omega))
      #Add it to the list and name it
      omegas<-c(omegas,list(omega))
      names(omegas)[length(omegas)]<-preyname
    }#for prey
    predator[[i]]@omegas<-omegas
  }#for predators
  predator
}#calc_omegas()

#'Run simulations
#'
#' This function drives the model.  It creates and initializes objects by
#' reading two files (input_params.csv and diet_prefs.csv) into data frames,
#' replaces the "params" data frame if an alternative is passed in by the user,
#' and then steps through years and intra-year timesteps, calculating population
#' numbers-at-age, recruits, catch, etc. for each species. Species interactions
#' are possible because the calculations for each species are done at least once
#' per year.  It returns a list of data and plots for each species. Note: the order of
#' some operations matters. This function must initialize the simulationParams
#' object before Species objects; do prey calculations before predators;
#' calculate sardine numbers before anchovy (because anchovy reproduction
#' depends on sardine numbers); update population numbers before recruitment and
#' catch (because catch and recruitment may depend on population size).
#' @param params A data frame of parameters
#' @return A list with three components: $input, $data, and $plots.  The $input
#'   list contains the simulationParams object that was used to do the
#'   simulations.  The $data and $plots sub-lists both have named components for
#'   each species, so for example the plots for Anchovy can be accessed as
#'   <result>$plots$Anchovy, or the data frame with the numerical results for
#'   Anchovy can be accessed with <result>$data$Anchovy.
#' @export
#'
#' @examples
#' \dontrun{
#'  result<-simulate() #uses default input files
#'
#'  #Instead of using the default, modify some values in the parameter object
#'  my_params<-new("simulationParams") #create a new simulationParams object
#'  my_params<-set_simParam(my_params) #initialize the object with data from stored files
#'  p<-my_params@params #Get the input parameter data frame
#'  p[p$species=="Sardine",1:7] #Look at the parameters for one species.
#'  p$base[p$species=="Sardine" & p$parameter=="fished"]<-TRUE #Set fished to TRUE in the base scenario
#'  my_params@params<-p #put the changed value back into the simulationParams object.
#'
#'  #Re-run the model with the changed values (Alternatively, one could re-make the .csv input files).
#'  result<-runMICE(my_params) #uses a modified simulationParams object
#'
#'  #Look at the results
#'  result$plots$Pelican #plot pelican data (several plots are produced)
#'  anchovy.df<-result$data$Anchovy #save the Anchovy data to a data frame for further analysis
#' }
runMICE<-function(Params,inputdir,scenario="base",nYears=100,nSteps=4){

  plotsfile<-system.file("extdata", "plotlist.csv", package = "MICE")
  plotlist<-read.csv(plotsfile,stringsAsFactors=FALSE,header=T,sep=",")

 #Create the simulation object and set its values.  If the user has supplied a
 #new version of simParams, use that instead of the default
if(missing(Params)){
  simParams<-new("simulationParams")
  simParams<-set_simParam(simParams,inputdir,scenario, nYears, nSteps)
} else{
  simParams<-Params
}

 #Extract a few global parameters
 scenario<-simParams@scenario
 nYears<-simParams@nYears
 nSteps<-simParams@nSteps

 #Create and initialize the prey objects
 Anchovy<-new("Anchovy")
 Anchovy<-set_species(Anchovy,simParams)
 Sardine<-new("Sardine")
 Sardine<-set_species(Sardine,simParams)
 Other<-new("Other_forage")
 Other<-set_species(Other,simParams)

 #And the predators
 Sealion<-new("Sealion")
 Sealion<-set_species(Sealion,simParams)
 Pelican<-new("Pelican")
 Pelican<-set_species(Pelican,simParams)

 #Stick them in a list for convenience
 prey<-list(sardine=Sardine,anchovy=Anchovy,other=Other)
 predator<-list(sealion=Sealion,pelican=Pelican)

 #Run 1 simulation (nYears) with unfished prey to calculate mean biomasses.
 #Returns modified predator objects in the predator list with updated values for omegas,
 #which scale prey availability relative to unfished mean biomasses.  The prey objects are not affected.
 predator<-calc_omegas(predator,prey,nYears,nSteps)
 pb <- txtProgressBar(min = 0, max = nYears, style = 3)

 #Main loop
 for (y in 2:nYears){
   #Do prey calculations
   for(i in 1:length(prey)){
     if(!prey[[i]]@uses_timestep){
       prey[[i]]<-calc_Ny(prey[[i]],y)
       prey[[i]]<-calc_recruits(prey[[i]],prey,y)
       prey[[i]]<-calc_catch(prey[[i]],y)
     } else{
       for(t in 1:nSteps){
         ts<-t + ((y-1) * nSteps) #The step number, as a constantly increasing series.  Remember, year 1 has been initialized already.
         prey[[i]]<-calc_Nt(prey[[i]],ts,t,nSteps) #Calculate numbers
         if(t==1){
           prey[[i]]<-calc_recruits(species=prey[[i]],preylist=prey,ts=ts,y=y) #Calculate recruits
         }
         prey[[i]]<-calc_catch(prey[[i]],ts) #calculate catch
       }#for t 1:nSteps
     }#uses_timestep
   }#for preyx in prey
 } #TEMPORARY

 for (y in 2:nYears){
   #Do predator calculations
   for(i in 1:length(predator)){
     if(!predator[[i]]@uses_timestep){
       predator[[i]]<-calc_Ny(predator[[i]],y) #calculate Numbers at age
       predator[[i]]<-calc_recruits(predator[[i]],prey,y,ts=NA,nSteps) #pass in a list of prey objects because predator recruitment depends on prey numbers
       predator[[i]]<-calc_catch(predator[[i]],y) #calculate catch
     } else{
       for(t in 1:nSteps){
         ts<-t + ((y-1) * nSteps) #ts = the step number as a constantly increasing series.  Remember, year 1 has been initialized already.
         predator[[i]]<-calc_Nt(predator[[i]],ts,t,nSteps) #calculate numbers
         if(t==1){
           predator[[i]]<-calc_recruits(predator[[i]],prey,y,ts,nSteps)
         }
         predator[[i]]<-calc_catch(predator[[i]],ts)
       }#for t 1:nSteps
     }#uses_timestep
   }#for predator[[i]] in predators
   setTxtProgressBar(pb, y)
 }#for year

 #Prepare results
 res<-list(simParams)
 names(res)[1]<-"input_params"

 #Create a list of objects that will be added to the output
 objs<-list()

 #Calculate metrics and produce plots
 for(preyx in prey){
   #Add the prey object to the results
   objs<-c(objs,list(preyx))
   names(objs)[length(objs)]<-preyx@name
   #Compile results and plots
   preyresx<-compile_results(preyx,plotlist,nYears,nSteps)
   res<-c(res,list(preyresx))
   names(res)[length(res)]<-preyx@name
 }
 for(predx in predator){
   #Add the prey object to the results
   objs<-c(objs,list(predx))
   names(objs)[length(objs)]<-predx@name
   #Compile results and plots
   predresx<-compile_results(predx,plotlist,nYears,nSteps)
   res<-c(res,list(predresx))
   names(res)[length(res)]<-predx@name
 }

 res<-c(res,list(objs))
 names(res)[length(res)]<-"species"

 close(pb)
 system("say -v Vicki Finished!")
 res
}#end main()
