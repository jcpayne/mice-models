#' Calculate total mortality rate
#'
#' A utility function used to calculate total mortality rate by age in species
#' that are tracked by timestep.
#'
#' @param Mv Natural mortality rate (log scale) by age. Vector.
#' @param L Number of timesteps per year (scalar).
#' @param Ft Total fishing mortality rate (scalar).
#' @param Sv Fishing selectivity (vector).
#'
#' @return A vector of mortality rate (log scale) by age.
#' @export
#' @keywords internal
calc_mortality_t<-function(Mv,L,Ft,Sv){
  zv<-Mv/L + Ft * Sv
  zv
}

#
#' Calculate total mortality rate
#'
#' A utility function used to calculate total mortality rate by age in species
#' that are tracked by year.
#'
#' @param Mv Natural mortality rate (log scale) by age. Vector.
#' @param Ft Total fishing mortality rate (scalar).
#' @param Sv Fishing selectivity (vector).
#'
#' @return A vector of mortality rate (log scale) by age.
#' @export
#' @keywords internal
calc_mortality_y<-function(Mv,Ft,Sv){
  zv<-Mv + Ft * Sv
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
#' mortality rate (log scale) by age (Zmat[,y-1]); the current year's
#' fishing mortality rate (Fv[y], log scale, a scalar, fished species only) and
#' the current year age 1+ biomass, B1P[y]
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
  Mv<-Species@Mv #natural mortality at age
  Nvprev<-Species@Nmat[,y-1] #Last year's N
  wv<-Species@wv #weight at age

  #Calculate total mortality rate for year y-1.
  if(fished){
    Sv<-Species@Sv #fishing gear selectivity at age
    Ftprev<-Species@Fv[y-1] #Last year's fishing mortality rate
    Zvprev<-calc_mortality_y(Mv,Ftprev,Sv) #last year's total mortality rate
  } else{
    #Set Z = natural mortality for all age classes
    Zvprev<-numeric(P+1) #a 1-d total mortality vector for year y-1
    Zvprev<-Mv
  }

  #For year 1, we assume only natural mortality (this function is first called
  #in y=2, so the previous year's mortality is Mv.)
  if(y==2){
    Zvprev<-numeric(P+1)
    Zvprev<-Mv
  }
  Species@Zmat[,y-1]<-Zvprev

  #Calculate N for this year
  leslie<-set_lesliematrix(P,Zvprev)
  Nv<-leslie %*% Nvprev #multiply numbers by survival for 0+ fish
  #Put Nv into the Species object
  Species@Nmat[,y]<-Nv #Numbers

  #If it is a prey species, calculate biomass of age-1+ animals and set that in the Species object
  #Also set B1P/B0
  if(is(Species,"Prey")){
    Species@B1P[y]<-calc_B1P(Nv,wv)
  }

  #Calculate this year's total fishing mortality rate (a scalar) from the fishing_rule model.
  #It is done after setting Nv, because the fishing rule may depend on this year's N.
  #The fishing rule calculation is done in calc_Fy; calc_catch is just a utility function.
  if(fished){
    Fy<-calc_Fy(Species,y)
    Cy<-calc_catch(Nv,wv,Mv,Sv,Fy) #Catch

    #Set the fishing mortality rate and the catch in the Species object
    Species@Fv[y]<-Fy
    Species@Cy[y]<-Cy
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
#' (Zmat[,t-1], log scale); the current fishing mortality rate (Fv[t], log
#' scale) and catch; fished species only); and the current age 1+ biomass, B1P[t]
#' (prey species only).
#' @param Species A Species object.
#' @param y Simulation year
#' @param ts Global timestep (ts = t * y)
#' @param t Within-year timestep
#' @param L Number of within-year timesteps per year.
#'
#' @return A modified Species object (see Description)
#' @export
#' @family Basic_dynamics
#' @keywords internal
calc_Nt<-function(Species,y,ts,t,L) {
  #Get values we will need
  #browser()
  P<-Species@P
  nAges<-P+1
  Ntprev<-Species@Nmat[,ts-1] #Numbers-at-age from the previous timestep
  Nt<-as.numeric(rep(NA,length(Ntprev))) #Numbers at age from this timestep (NA, to start with)
  Mv<-Species@Mv
  wv<-Species@wv
  fished<-Species@fished

  #Calculate mortality rate-at-age for the previous timestep
  if(fished){
    Sv<-Species@Sv #selectivity-at-age
    Ftprev<-Species@Fv[ts-1] #Retrieve last year's fishing mortality rate
    Zvprev<-calc_mortality_t(Mv,L,Ftprev,Sv)
  } else{
    #Set the mortality for each age to natural mortality divided by the number of steps per year.
    Zvprev<-numeric(nAges)
    Zvprev<-Mv/L #mortality/nSteps
  }
  Species@Zmat[,ts-1]<-Zvprev

  #Do the update (at t=1 we use the Leslie matrix; for later steps fish stay in the same age class)
  if(t==1){
    leslie<-set_lesliematrix(P,Zvprev) #Set up the Leslie matrix with the survivals
    Nt<-leslie %*% Ntprev
  }  else {
    #Leaving this in for now, but NOTE: mortality between y-1 to y must be full Z rate (not Z/nSteps) if you do it this way
    #Nt[1]<-Ntprev[1] #Mortality = 0 for recruits across all timesteps within a year
    Nt[1:nAges]<-Ntprev[1:nAges] * exp(-Zvprev[1:nAges]) #Fish stay in the same age classes but suffer mortality
  }

  #Set Nt in the Species object
  Species@Nmat[,ts]<-Nt

  #If it is a prey species, calculate biomass of age-1+ animals and set that in the Species object
  if(is(Species,"Prey")){
    Species@B1P[ts]<-calc_B1P(Nt,wv)
  }

  #Calculate this year's total fishing mortality rate (a scalar) from the fishing_rule model.
  #It is done after setting Nt, because the fishing rule may depend on this year's N.
  #Set it in the Species object.
  if(fished){
    #Get the catch total for this year so far.
    Cy<-Species@Cy[y]

    #Calculate Ft, the total fishing mortality rate (a scalar) for this timestep from the fishing_rule model.
    Ft<-calc_Ft(Species,y,ts,t)

    #Calculate catch for the timestep and update the year's catch sum.
    Ct<-calc_catch(Nt,wv,Mv,Sv,Ft)
    Cy<-Cy + Ct

    #Set the fishing mortality rate and the catches in the Species object
    Species@Fv[ts]<-Ft
    Species@Ct[ts]<-Ct
    Species@Cy[y]<-Cy
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
#' @describeIn calc_recruits Other_forage method.
setMethod("calc_recruits", signature("Other_forage"),function (species,y,...) {
  calc_other_R(species,y,...)
})

#' The Predator method calculates recruits for Predators.
#'
#' @describeIn calc_recruits Predator method.
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
#' @param maturev maturity at age
#'
#' @return A scalar: biomass in tons
#' @export
#' @keywords internal
#' @family recruitment
calc_B_sp<-function(wv,Nv,maturev){
  Bs<-sum(wv * Nv * maturev)/1000
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
  useG<-Prey@recruitment$useG
  G<-Prey@G
  sigma2<-(Prey@recruitment$sigma)^2
  R_thresh<-Prey@recruitment$R_thresh

  #Calculate spawning biomass
  Nv<-Prey@Nmat[,t]
  wv<-Prey@wv
  maturev<-Prey@maturev
  B_sp<-calc_B_sp(wv,Nv,maturev)

  # if(is(Prey,"Anchovy")){
  #   browser()
  # }

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

    #Ricker with environmental driver G
    R<-calc_Ricker(alpha,beta,B_sp) #deterministic component

    #Multiply it by an environmental driver.  By default, anchovy does not use this (useG=0: see eqn 6.)
    R<-R * exp(useG * G[y]) #environmental component

    #Draw a spawner-recruit residual
    #epsilon will equal 0 if R_stochastic is false, and may be autocorrelated depending on the value of rho
    #(see calc_R_error, which is called at initialization).
    if(Prey@recruitment$R_stochastic){
      epsilon<-Prey@lnSR_err[t]

      #Overwrite it for anchovy if resample_residuals=T
      if(is(Prey,"Anchovy") && Prey@recruitment$resample_residuals){
        #Resample recruitment errors from values that were read in (2 series, depending on value of SSB)
        if(B_sp <= 500){
          epsilon<-sample(Prey@R_residuals$SSB_LE_500K,1, replace=T)
        } else{
          epsilon<-sample(Prey@R_residuals$SSB_LE_500K,1, replace=T)
        }
        #Put the resampled errors in the lnSR_err slot of the species object (since it can't be known earlier)
        Prey@lnSR_err[t]<-epsilon
      }#if resample residuals
    } else{
      epsilon<-0
      sigma2<-0
    }

    #Add the stochastic error to the recruitment
    R<- R * exp(epsilon - sigma2/2)

    #For anchovy (and possibly other prey), we force recruitment to fail a certain
    #percentage of the time regardless of all of the above calculations.
    if(exists('p_fail',where=Prey@recruitment)){
      p_fail<-Prey@recruitment$p_fail
      if(p_fail != 0){
          if(Prey@recruitment$R_stochastic){
            rnd<-runif(1) #draw a random uniform[0,1]
            if(rnd < p_fail){
              R<-0
            }
          } else {
            #The average value of the failure function is (1 - p_fail).  If it fails
            #25% of the time, then on average it succeeds 75% of the time.  We
            #multiply R by (1 - p_fail) so that whether or not it's switching on/off, we have
            #the same average value.
            R<-R * (1 - p_fail)
          }
        }#if p_fail != 0
      }#if p_fail exists
    }#B_sp is not below threshold

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
  R<-(1000 * alpha * B_sp * exp(-beta * B_sp))
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
  useG<-Other@recruitment$useG
  G<-Other@G

  #Get the fixed component of recruitment
  R_oth<-Other@R_baseyear

  #Add an environmental driver effect if useG=1 (kept separate from error so we
  #can plot just the environmental effect)
  if(Other@recruitment$useG==1){
    R_oth<-R_oth * exp(G[y])
  }
  #Add stochastic error if R_stochastic is true
  if(Other@recruitment$R_stochastic){
    R_oth<-R_oth * exp(epsilon -((sigma^2)/2)) #This year's recruitment
  }
  Other@Nmat[1,y]<-R_oth
  Other@R[y]<-R_oth
  #Set spawning biomass (only used for plotting)
  Other@B_sp[y]<-calc_B_sp(Other@wv,Other@Nmat[,y],Other@maturev)
  Other
}

##Predators
###Utility functions
#These are for the various components of recruitment (prey biomass-at-age, predator preferences for species and age-classes, and spatial availability of the prey to the predator).
#Calculate the number of mature predators:
#Arguments: Current-year numbers at age (Nv) and proportion mature at age (maturev) (both vectors).
#Returns a scalar (total number of mature predators).
calc_Nm<-function(Nv,maturev){
  Nm<-Nv %*% maturev
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
#' @param B0list A named list of available prey biomass for each species.
#' @param preylist A list of Prey objects.
#' @return A value for Dy (scalar).
#' @export
#' @keywords internal
#' @family recruitment
calc_Dy<-function(predator,B0list, preylist){
 OP<-predator@sp_prefs$Other_prey #this predator's preference for Other_prey
 Dy<-OP
 for(i in 1: length(preylist)){
   preyname<-preylist[[i]]@name
   omega<-as.numeric(unlist(predator@omegas[preyname])) #get the omega
   Dy<-Dy + as.numeric(unlist(B0list[preyname])) * omega
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
 D0<-predator@D0 #TODO: perhaps in the denominator, use Dy/D0 - Dhat instead.
 f<-1/(1 + exp(-k * (Dy-Dhat)))
 p<-Mo * f
 p
}

#' Calculate phi (impact of prey on predator reproduction)
#'
#' Calculates the impact of prey on predator reproduction; see equation 12 of
#' Punt et al. 2016
#'
#' @param predator A Predator object
#' @param Dy The relative (scaled) biomass of prey available to the predator.
#'   Note: we're using one set of thetas for all prey species, but we could
#'   potentially make them different for each prey (Punt et al. set them all equal
#'   in the base scenario.)
#'
#' @return The parameter phi
#' @export
#'
#' @keywords internal
#' @family recruitment
calc_phi<-function(predator,Dy){
  D0<-predator@D0 #the average biomass available in the absence of fishing
  t1<-predator@prey_effect$theta1
  t2<-predator@prey_effect$theta2
  t3<-predator@prey_effect$theta3
  f<-((1 - t1 - t2) * t3 * (Dy/D0 - t1))/(((1-t1) * t2 * (1 - t3)) + (t3 * (1 - t1) - t2) * (Dy/D0 - t1))
  phi<-max(0,f)
}

#' Calculate recruits using Punt method
#'
#' Calculates recruits using line 1 of equation 9 in Punt et al. 2016.
#'
#' @param predator A Predator object
#' @param Nm Number of mature predators (a scalar)
#' @param phi The prey impact parameter
#'
#' @return The number of recruits (a scalar)
#' @export
#'
#' @keywords internal
#' @family recruitment
calc_predRecruits_Punt<-function(predator,Nm,phi,N1P,K1P){
  bigPhi<-predator@prey_effect$bigPhi
  if(N1P/K1P <= 0){
    R<-(Nm * phi)
  }else{
    densDep<-max(0,(1 + (bigPhi - 1) * (1 - (phi * N1P/K1P)^2.39)))
    R<-(Nm * phi * densDep)
  }
  R
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
  maturev<-Predator@maturev #Proportion mature at age
  BHc<-Predator@recruitment$BHc #Beverton-Holt carrying capacity.  Steepness (BHp) is determined below by prey effect.
  R_thresh<-Predator@recruitment$R_thresh #recruitment threshold (minimum number of mature adults)

  #Calculate the number of mature predators and the number of age-1+ predators
  Nm<-calc_Nm(Nv,maturev)
  N1P<-sum(Nv[2:nAges])

  #Calculate prey biomass
  if(Nm < R_thresh){
    R<-0
  } else{
    Bprey<-list()
    #Calculate the available, age-preference-adjusted biomass of each prey species
    for(i in 1:length(preylist)){
      preyname<-preylist[[i]]@name
      ap<-Predator@prey_availability[[preyname]]$ap #the ap parameter for availabililty
      ac<-Predator@prey_availability[[preyname]]$ac

      #Get prey biomass from the first timestep of year y (or just year y if not using timesteps)
      #NOTE! The predator may not be on a timestep, in which case the passed pararmeter ts=NA.
      #So we re-calculate tt here because the prey *may* be on a timestep.
      if(preylist[[i]]@uses_timestep){
        tt<-(((y-1) * nSteps) + 1)
      } else {
        tt<-y
      }
      #First get age-preferred biomass
      B_agepref<-calc_Bprey_agepreferred(Predator,preylist[[i]],tt)
      #Now adjust it for spatial availability
      B_avail<-calc_availability(ap,ac,B_agepref)

      #Save the before-and-after-availability biomasses in the predator object for graphing
      Predator@prey_availability[[preyname]]$B_agepref[t]<-B_agepref
      Predator@prey_availability[[preyname]]$B_avail[t]<-B_avail

      #Add the biomass to a list and name it
      Bprey<-c(Bprey,list(B_avail))
      names(Bprey)[length(Bprey)]<-preylist[[i]]@name
    }#for i in length(preylist)

    #Calculate Dy, the total scaled available prey biomass
    Dy<-calc_Dy(Predator,Bprey,preylist)
    #debugging
    #print(paste("Bagepref:",B_agepref,"; B_avail",B_avail,"; Dy:",Dy))

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

    #Calculate phi, the Punt version of sensitivity of predator reproduction to prey
    phi<-calc_phi(Predator,Dy)

    #Calculate R_Punt, the Punt version of recruits
    K1P<-Predator@K1P #average number of age-1+ predators when prey are not fished
    R_Punt<-calc_predRecruits_Punt(Predator,Nm,phi,N1P,K1P)
    RperS_Punt<-R_Punt/Nm
  }#if Number of mature adults Nm is above recruitment threshold R_thresh

  #Set values in the Predator object
  Predator@Dy[t]<-Dy
  Predator@BHp[t]<-BHp
  Predator@Nm[t]<-Nm
  Predator@N1P[t]<-N1P
  Predator@NperK[t]<-N1P/Predator@K1P
  Predator@Nmat[1,t]<-R
  Predator@R[t]<-R #Set the R value (for convenience)
  Predator@RperS[t]<-R/Nm
  #Punt alternative method for reproduction
  Predator@phi[t]<-phi
  Predator@R_Punt[t]<-R_Punt
  Predator@RperS_Punt[t]<-RperS_Punt
  Predator #return the Predator
}#calcPredator_recruits()

#' Calculate fishing mortality according to a fishing rule
#'
#'There are two versions: one for species tracked by year, the other for species
#'tracked by timestep. Note: it would be easy to make this a generic and set up
#'alternative models for each species. This version is for species that are
#'tracked by timestep.  ts= t + (y * t).
#'
#' @param Species A Species object.
#' @param y The simulation year.
#' @param ts The global timestep ()
#' @param t The intra-year timestep.
#'
#' @return A scalar, Ft, the total fishing mortality rate (log scale) for timestep t.
#' @export
#' @family fishing
#' @keywords internal
calc_Ft<-function(Species,y,ts,t){
  Mv<-Species@Mv
  Ft<-Species@F_baseline
  gv<-Species@gv #proportion fished per time period
  Bmin<-Species@fishing_rule$Bmin #Minimum biomass of age 1+ fish, below which fishing is stopped.
  max_catch<-Species@fishing_rule$max_catch

  #calculate the observed biomass of age 1+ fish
  Nv<-Species@Nmat[,ts]
  wv<-Species@wv
  Sv<-Species@Sv
  B1P<-calc_B1P(Nv,wv)

  #The rule: if biomass < Bmin, set catch to zero
  if(B1P < Bmin){
   Ft<-0
  }
  Ft<-Ft * gv[t]

  #Calculate catch provisionally, and adjust F so that catch won't exceed max_catch.
  #Catch is not kept here.
  C_year<-Species@Cy[y]
  C_timestep<-calc_catch(Nv,wv,Mv,Sv,Ft)
  #check to see if the timestep catch plus the year catch is too large
  if((C_year + C_timestep) > max_catch){
    C_timestep<-(max_catch - C_year)
    Ft<-solve_for_Ft(C_timestep,Nv,wv,Mv,Sv,Ft)
  }
  Ft
}

#' Calculate fishing mortality according to a fishing rule
#'
#'This version is for species tracked by year.  Note: it would be easy to make
#'this a generic and set up alternative models for each species.
#' @param Species A Species object.
#' @param y The current year.
#' @param t The intra-year timestep.
#' @return A scalar, Fv: the total fishing mortality rate year (log scale) for year y.
#' @export
#' @family fishing
#' @keywords internal
calc_Fy<-function(Species,y){

  #Set fishing mortality at baseline (default)
  Fy<-Species@F_baseline
  Bmin<-Species@fishing_rule$Bmin_cutoff #Minimum biomass of age 1+ fish, below which fishing is stopped.
  max_catch<-Species@fishing_rule$max_catch
  Sv<-Species@Sv

  #calculate the observed biomass of age 1+ fish
  Nv<-Species@Nmat[,y]
  wv<-Species@wv
  B1P<-calc_B1P(Nv,wv)

  #The rule: if biomass < Bmin, set catch to zero
  if(B1P < Bmin){
    Fy<-0
  }

  #Calculate catch provisionally, and adjust F so that catch won't exceed max_catch.
  #Catch is not kept here.
  tempcatch<-calc_catch(Nv,wv,Mv,Sv,Fy)
  if(tempcatch > max_catch){
    tempcatch<-max_catch
    Fy<- solve_for_Ft(tempcatch,Nv,wv,Mv,Sv,Fy)
  }
  Fy
}

#' Calculate catch (biomass units)
#'
#' A utility function for calculating catch at time t.  Uses the Baranov equation,
#' with fishing selectivity by age.
#'
#' @param Species A Species object
#' @param t time period (year or timestep)
#' @param Ft Instantaneous fishing mortality rate (log scale); a scalar.
#' @param Nv Vector of numbers at age for a particular time t
#' @param wv Vector of weight at age
#' @param Sv Vector of fishing gear selectivity at age
#' @param Mv Vector of natural mortality rate at age
#'
#' @return Catch in biomass units.
#' @export
#' @family fishing
#' @keywords internal
calc_catch<-function(Nv,wv,Mv,Sv,Ft){
  Cv<-sapply(1:length(Nv),function(i) Nv[i]*(Sv[i] * Ft/Mv[i]) * (1-exp(-((Sv[i] * Ft) +  Mv[i]))))
  C<-sum(Cv)/1000
  C
}

#' Solve for F
#'
#'Iteratively solve for instantaneous fishing mortality, F, given some target
#'catch.  Assumes the target is between 0 and 2x the initial guess.
#'
#' @param C Target catch.
#' @param Nv Number at age vector.
#' @param wv Weight at age vector.
#' @param Mv Mortality at age vector.
#' @param Sv Selectivity at age vector.
#' @param Fstart Starting guess at F.
#'
#' @return Instantaneous fishing mortality rate, Ft (a scalar).
#' @export
#'
#' @keywords internal
solve_for_Ft<-function(C,Nv,wv,Mv,Sv,Fstart){
  Fmin<-0
  Fmax<-2
  target<-C
  for(i in 1:50){
    mult<-(Fmin + Fmax)/2
    Ftest<-Fstart * mult
    predicted<-calc_catch(Nv,wv,Mv,Sv,Ftest)
    if(predicted < target){
      Fmin<-mult
    } else{
      Fmax<-mult
    }
  }
  Ft<-Fstart * mult
  Ft
}

#' Do a calibration run.
#'
#' Do an initial run for prey species without fishing, to calculate parameters
#' in the Predator objects (omegas) that scale the unfished biomass of each prey species.
#'
#' @param predators A list of Predator objects
#' @param prey A list of Prey objects
#' @param nYears Number of years for the simulation
#' @param nSteps Number of timesteps per year.
#'
#' @return A modified list of predators and prey.  In the predators, the @omegas
#'   slot, which holds scaling parameters for the ratio of diet preference to
#'   unfished biomass for each prey species, has been set.  The function also
#'   sets D0, the scaled, available prey biomass at unfished levels.  In the
#'   prey objects, the age-1+ unfished biomass, B0, is set.
#' @export
#' @keywords internal
#' @family initialization
calibrate<-function(predators,prey,nYears,nSteps){
  print("Calibration run in progress...")
  starttime<-Sys.time()
  pb <- txtProgressBar(min = 0, max = nYears, style = 3)

  #Run the main model without fishing or stochastic recruitment
  for(i in 1:length(prey)){
    prey[[i]]@fished<-FALSE
  }
  #Simulate prey only
  for (y in 2:nYears){
    #Do prey calculations
    for(i in 1:length(prey)){
      if(!prey[[i]]@uses_timestep){
        prey[[i]]<-calc_Ny(prey[[i]],y)
        prey[[i]]<-calc_recruits(prey[[i]],prey,y)
        #prey[[i]]<-calc_catch(prey[[i]],y)
      } else{
        for(t in 1:nSteps){
          ts<-t + ((y-1) * nSteps) #The step number, as a constantly increasing series.  Remember, year 1 has been initialized already.
          prey[[i]]<-calc_Nt(prey[[i]],y,ts,t,nSteps) #Calculate numbers
          if(t==1){
            prey[[i]]<-calc_recruits(species=prey[[i]],preylist=prey,ts=ts,y=y) #Calculate recruits
          }
          #prey[[i]]<-calc_catch(prey[[i]],ts) #calculate catch
        }#for t 1:nSteps
      }#uses_timestep
    }#for preyx in prey
    setTxtProgressBar(pb, y)
  }#for y in years

  #Set the prey's mean age-1+ unfished biomass, B0.  This version is NOT
  #adjusted for age or species preferences; and it's only used for plotting
  #"Other_forage" B1+/B0.  Note that for species that use timesteps, we only
  #count biomass in step 1, because the biomass fluctuates through the year
  for(i in 1:length(prey)){
    nAges<-prey[[i]]@P + 1
    nmat<-prey[[i]]@Nmat
    wv<-prey[[i]]@wv
    B_coltotals<-colSums(nmat[2:nAges,] * wv[2:nAges])/1000 #in tons, to match B1P
    #Select only step 1, if the species uses timesteps
    if(prey[[i]]@uses_timestep){
      B_coltotals<-B_coltotals[seq(1,length(B_coltotals),nSteps)]
    }
    #Extract timestep 1 totals
    #Get the mean
    B0<-mean(B_coltotals,na.rm=T)

    #Set B0 value in the prey object
    prey[[i]]@B0<-B0
  }#for i in prey

 #Set omega=(species preference/unfished age-preferred biomass) for each predator/prey combination.
  for(i in 1:length(predators)){
    omegas<-list()
    B0list<-list() #Total prey biomass (adjusted for age and availability) for this predator (a named list)
    for(j in 1:length(prey)){
      preyname<-prey[[j]]@name

      spx_pref<-as.numeric(unlist(predators[[i]]@sp_prefs[preyname])) #the preference of predator i for prey species j
      #Calculate mean age-preferred biomass of prey j
      nmat<-prey[[j]]@Nmat #prey N matrix
      wv<-prey[[j]]@wv

      #get the predator's age prefs for this species
      ageprefsx<-as.numeric(unlist(predators[[i]]@age_prefs[preyname]))
      #Get total age-preferred biomass (a vector nTimes long)
      B_colsums<-colSums(nmat * wv * ageprefsx)/1000 #in tons, to match B1P
      #Take just the 1st timestep if it's a species tracked by timestep
      if(prey[[j]]@uses_timestep){
        B_agepreftotal<-B_colsums[seq(1,length(B_colsums),nSteps)]
      } else{
        B_agepreftotal<-B_colsums
      }
      #Adjust it for availability
      ap<-predators[[i]]@prey_availability[[preyname]]$ap
      ac<-predators[[i]]@prey_availability[[preyname]]$ac
      Bavail<-calc_availability(ap,ac,B_agepreftotal)
      B_ap_unfished<-mean(Bavail[50:length(Bavail)]) #discard the first 50 years to let age structure settle

      #Add the age- and availability-adjusted biomass of this prey species to the total
      B0list<-c(B0list,list(B_ap_unfished))
      names(B0list)[length(B0list)]<-prey[[j]]@name

      #Calculate omega = (pred preference for prey x) / (age-preferred biomass of prey x)
      omega<-spx_pref/B_ap_unfished
      #print(paste(predators[[i]]@name,prey[[j]]@name,B_ap_unfished,spx_pref,omega))
      #Add it to the list and name it
      omegas<-c(omegas,list(omega))
      names(omegas)[length(omegas)]<-preyname

    }#for prey

    #Set the omegas
    predators[[i]]@omegas<-omegas

    #Calc D0 (note: this uses the mean age-preferred avaiable biomass and omegas we just set)
    D0<-calc_Dy(predators[[i]],B0list,prey)
    predators[[i]]@D0<-D0
  }#for predators

  close(pb)

  print("Calibration finished.")

  #return predators and prey in a list
  specieslist<-list(predator=predators,prey=prey)
  specieslist
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
#' @param inputdir A directory where the input_params.csv and diet_prefs.csv files are.
#' @param scenario Character.  A scenario name (must match a column name in the input file.)
#' @param nYears Integer.  Number of years for the simulation.
#' @param nSteps Integer.  Number of timesteps within a year (ignored by some species).
#' @return A list with several components:
#'  $input_params: the simulationParams object that was used to make the run.
#'  $Anchovy: anchovy results.  This has two components:
#'     $data: a data frame of some of the most commonly-used results
#'     $plots: a set of plots, which may be accessed as a group or individually
#'  $Sardine (same as Anchovy), etc.
#'  $species.  The species objects that contain all values used.
#'     $Anchovy: the anchovy object
#'     $Sardine: the sardine object, etc.
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
#'  result$Pelican$plots #Show Pelican plots
#'  anchovy.df<-result$Anchovy$data #save the Anchovy data to a data frame for further analysis
#'
#'  #Look at a slot in one of the Species objects
#'  res$species$Anchovy@wv #weight at age
#' }
runMICE<-function(Params,inputdir,scenario="base",nYears=100,nSteps=4){

  starttime<-Sys.time()

  plotsfile<-system.file("extdata", "plotlist.csv", package = "MICE")
  plotlist<-read.csv(plotsfile,stringsAsFactors=FALSE,header=T,sep=",")

  print("Loading parameters and initializing objects...")
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
  preynames<-c("Anchovy","Sardine","Other_forage") #needed for various predator operations that iterate over prey

  #And the predators
  Sealion<-new("Sealion")
  Sealion<-set_species(Sealion,simParams,preynames)
  Pelican<-new("Pelican")
  Pelican<-set_species(Pelican,simParams,preynames)

  #Set environmental driver G equal for Anchovy and Other_forage if desired
  if(Other@recruitment$useAnchovyG){
    Other@G<-Anchovy@G
  }

  #Stick them in a list for convenience
  prey<-list(sardine=Sardine,anchovy=Anchovy,other=Other)
  predator<-list(sealion=Sealion,pelican=Pelican)

  #Run 1 simulation (nYears) with unfished prey to calculate mean biomasses.
  #Returns modified predator objects in the predator list with updated values for omegas,
  #which scale prey availability relative to unfished mean biomasses.  The prey objects are not affected.
  calibrated_list<-calibrate(predator,prey,nYears,nSteps)

  #Copy over the values we want (the intent is to NOT copy the rest of the
  #changes but if I switch to R6 classes the calibration will have to be run separately)
  calibrated_predlist<-calibrated_list$predator
  calibrated_preylist<-calibrated_list$prey
  for (i in 1:length(predator)){
    stopifnot(predator[[i]]@name == calibrated_predlist[[i]]@name)
    predator[[i]]@omegas<-calibrated_predlist[[i]]@omegas
    predator[[i]]@D0<-calibrated_predlist[[i]]@D0
  }
  for(i in 1:length(prey)){
    stopifnot(prey[[i]]@name == calibrated_preylist[[i]]@name)
    prey[[i]]@B0<-calibrated_preylist[[i]]@B0
  }

  #progress bar
  pb <- txtProgressBar(min = 0, max = (nYears * 2), style = 3)

  #Main loop
  for (y in 2:nYears){
   #Do prey calculations
   for(i in 1:length(prey)){
     if(!prey[[i]]@uses_timestep){
       prey[[i]]<-calc_Ny(prey[[i]],y)
       prey[[i]]<-calc_recruits(prey[[i]],prey,y)
       #prey[[i]]<-calc_catch(prey[[i]],y)
     } else{
       for(t in 1:nSteps){
         ts<-t + ((y-1) * nSteps) #The step number, as a constantly increasing series.  Remember, year 1 has been initialized already.
         prey[[i]]<-calc_Nt(prey[[i]],y,ts,t,nSteps) #Calculate numbers
         if(t==1){
           prey[[i]]<-calc_recruits(species=prey[[i]],preylist=prey,ts=ts,y=y) #Calculate recruits
         }
         #prey[[i]]<-calc_catch(prey[[i]],ts) #calculate catch
       }#for t 1:nSteps
     }#uses_timestep
   }#for preyx in prey
   setTxtProgressBar(pb, y)
  } #for y in nYears

  for (y in 2:nYears){
   #Do predator calculations
   for(i in 1:length(predator)){
     if(!predator[[i]]@uses_timestep){
       predator[[i]]<-calc_Ny(predator[[i]],y) #calculate Numbers at age
       predator[[i]]<-calc_recruits(predator[[i]],prey,y,ts=NA,nSteps) #pass in a list of prey objects because predator recruitment depends on prey numbers
       #predator[[i]]<-calc_catch(predator[[i]],y) #calculate catch
     } else{
       for(t in 1:nSteps){
         ts<-t + ((y-1) * nSteps) #ts = the step number as a constantly increasing series.  Remember, year 1 has been initialized already.
         predator[[i]]<-calc_Nt(predator[[i]],y,ts,t,nSteps) #calculate numbers
         if(t==1){
           predator[[i]]<-calc_recruits(predator[[i]],prey,y,ts,nSteps)
         }
         #predator[[i]]<-calc_catch(predator[[i]],ts)
       }#for t 1:nSteps
     }#uses_timestep
   }#for predator[[i]] in predators
    setTxtProgressBar(pb, (nYears + y))
  }#for year

  print("Preparing results...")
  #For graphing, we adjust B1P by B0 (can't be done earlier)
  for(i in 1:length(prey)){
    B1Px<-prey[[i]]@B1P
    #Only use biomass in step 1 of each year, because it fluctuates within the year
    if(prey[[i]]@uses_timestep){
      B1Px<-B1Px[seq(1,length(B1Px),nSteps)]
    }
    prey[[i]]@B1PperB0<-B1Px/prey[[i]]@B0
  }

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
    predresx<-compile_results(predx,plotlist,nYears,nSteps,preynames)
    res<-c(res,list(predresx))
    names(res)[length(res)]<-predx@name
  }

  #Add the species objects to the results list
  res<-c(res,list(objs))
  names(res)[length(res)]<-"objects"

  #Calculate stats and add to the results list
  stats<-getStats(res)
  res<-c(res,list(stats))
  names(res)[length(res)]<-"stats"

  #Create multi-species summary of biomass or numbers
  panelN<-make_panel_plot(res)
  res<-c(res,list(summaryplot2=panelN))
  #sard_effect_plot<-sardine_effect_plot(prey)

  close(pb)

  beep(10)
  #system("say -v Vicki Finished!")
  endtime<-Sys.time()
  print(paste("Elapsed time:",endtime - starttime,"seconds."))
  res
}#end main()
