#' Compile results and generate graphs for one species.
#'
#' Creates graphs by taking X and Y series from slots in a Species object, and
#' combines some of those data into a data frame that is returned with the
#' plots. Reads instructions from the plotlist.csv file to decide which plots
#' and series to produce (there are different sets of plots for prey and
#' predators).  Note that since a dataframe cannot be a ragged array, the number
#' of rows in the data frame that is returned by the function is either nYears or
#' nTimes, but certain items (catch, recruits) are *plotted* by year, even when
#' the species is tracked by timestep.
#'
#' @param species A Species object
#' @param slotlist A list of slots, from which graphs are to be made.
#' @param nYears Number of years in the simulation.
#' @param nSteps Number of timesteps per year (ignored for some species)
#'
#' @return A list with the values of certain vectors assembled into a data frame, and a list of plots.
#' @export
#' @family results
#' @keywords internal
compile_results<-function(species,plot_master,nYears,nSteps,preynames){
  if(!is(species,"Species")) stop("getresults() takes a Species object")

  pm<-plot_master
  pm[pm==""]<-NA #Don't skip this step!

  #Check to make sure we're not trying to plot a non-existent slot
  for (i in 1:nrow(pm)){
    if(is(species,pm$obj_type[i])){
      if(!is.na(pm$slot.x[i]) && !pm$slot.x[i] %in% slotNames(species)){
        stop(paste("Cannot plot",pm$slot.x[i],": slot not found in the",species@name,"object."))
      }
      if(!is.na(pm$slot.y[i]) && !pm$slot.y[i] %in% slotNames(species)){
        stop(paste("Cannot plot",pm$slot.y[i],": slot not found in the",species@name,"object."))
      }
    }
  }#for 1 to nrow(pm)

  speciesname<-species@name
  res<-list()
  plotlist<-list()

  #Set up two alternative X vectors (nYears, and nTimes) for the output dataframe
  #(NOTE: this x is only used for plotting if the xtype column of pm is "default".)
  x_ts<-1:(nSteps * nYears)
  x_yr<-1:nYears

  #Initialize the results dataframe (to which Y columns will be added).
  #For convenience, add a intra-year timestep to by-timestep dataframes
  #This dataframe is for output, not plotting.
  if(species@uses_timestep){
    data<-data.frame(x_ts,rep(1:nSteps,nYears))
    names(data)<-c("Timestep","Intr_year_step")
  } else{
    data<-data.frame(x_yr)
    names(data)<-"Year"
  }

  #Step through the plotmaster dataframe; make plots and add data to the result
  #dataframe as instructed.
  for(i in 1:nrow(pm)){
    yname<-NA
    pdata<-NA
    ptitle<-NA

    #Get x and y labels, and plot title if specified
    if(!is.na(pm$xlab[i])){
      xname<-pm$xlab[i]
    }
    if(!is.na(pm$title[i])){
      ptitle<-pm$title[i]
    }
    #Get the y range option
    y_range<-pm$y_range[i]

    #Get the predicted option
    predicted<-pm$predicted[i]

    series_type<-pm$series_type[i]

    #Make one set of plots for pm$obj_type=="Prey", a different set for Predators
    if(is(species,pm$obj_type[i])){
      if(!pm$plottype[i]=="availability"){

        #Get y axis label
        if(!is.na(pm$ylab[i])){
          yname<-pm$ylab[i]
        }else{
          yname<-pm$slot.y[i] #default; may be overridden
        }

        #Get y values
        vals<-slot(species,pm$slot.y[i])

        #Skip this row in pm if there are no data in the y slot.
        if(length(vals)==0){
          next
          #print(paste("There is no data for",yname,"in",speciesname)) #next
        }

        #Get x values if the x slot is specified, or create an x axis (time or years)
        if(!is.na(pm$slot.x[i])){
          px<-slot(species,pm$slot.x[i])
          xname<-pm$slot.x[i]
        } else{
          #Specify an x-axis for timesteps or years
          if(length(vals)==length(x_ts)){
            px<-x_ts
            xname<-"Timestep"
          } else {
            px<-x_yr
            xname<-"Year"
          }
        }

        #Make a small dataframe (pdata) for plotting
        if(length(vals) != length(px)){
          stop(paste("Cannot plot",yname,"; The length of",yname,"does not match the length of x (",length(px),")"))
        } else{
          pdata<-data.frame(px,vals)
        }

        #Subset it, in the case that the data are in timestep units but we want the plot in years,
        #and replace the x-axis with years.
        if(pm$xtype[i]=="step1only"){
          if(species@uses_timestep && nrow(pdata)==length(x_ts)){
            #First, subset the data
            pdata<-get_every_nth(nSteps,pdata,xname,yname) #subset the df
            #Don't change the x-axis or label unless the name is Timestep
            #because xy plots may have a different x-axis
            if(xname == "Timestep"){
              pdata<-data.frame(x_yr,pdata[,match(yname,names(pdata))])#replace the x column with year
              xname<-"Year"
            }
            names(pdata)<-c(xname,yname)
          }
        }

        #Add the (unshortened) Y timeseries to the results dataframe if in_dataframe=1 in plotmaster
        if(pm$in_dataframe[i]){
          data<-data.frame(data,vals)
          names(data)[ncol(data)]<-pm$slot.y[i] #name the column
        }
      }#plot type is not availability

      #Several plotting options:
      #Plot a time series
      if(pm$plottype[i]=="timeseries"){
        newplot<-makeplot(pdata,xname,yname,speciesname,ptitle,y_range,series_type)
      }

      #Plot a time series, summed by year
      if(pm$plottype[i]=="summed_by_year"){
        if(pm$slot.y[i]=="Ct" && species@uses_timestep && species@fished){
          xname<-"Year"
          newplot<-plot_summed_by_year(nYears,nSteps,pdata,xname,yname,speciesname,ptitle,y_range,series_type)
        }
      }

      #Plot two variables X and Y against each other
      if(pm$plottype[i]=="xy"){
        if(!is.na(predicted)){
          #get the predicted data
          fitted_data<-calc_predicted(species,predicted)
          names(fitted_data)<-names(pdata)
          #Send the observed and predicted data separately to ggplot
          newplot<-plot_fitted(fitted_data,pdata,xname,yname,speciesname,ptitle,y_range)
        } else{
          #pdata<-na.exclude(pdata) #remove rows with missing values
          pdata<-pdata[order(pdata[,1]),] #sort the data by x values.
          newplot<-makeplot(pdata,xname,yname,speciesname,ptitle,y_range,series_type)
        }
      }#plottype=xy

      #Plot just a function (no observed data)
      if(pm$plottype[i]=="function_only"){
        if(pm$predicted[i]=="availability"){
          fndata<-calc_predicted(species,predicted)
          newplot<-plot_function_only(fndata,xname,yname,speciesname,ptitle)
        }
      }

      #A special n-panel (1 per prey) plot of prey availability per predator
      if(pm$plottype[i]=="availability"){
        observed<-data.frame()
        predicted<-data.frame()
        for(k in 1:length(preynames)){
          #Put observed data in a dataframe
          x<-species@prey_availability[[preynames[[k]]]]$B_agepref
          y<-species@prey_availability[[preynames[[k]]]]$B_avail
          preysp<-rep(preynames[k],length(x))
          obs<-data.frame(preysp,x,y,rep("obs",length(x)))
          observed<-rbind(observed,obs)
        #Get predicted values
          pred<-calc_predicted(species,"availability",preynames[k])
          preysp<-rep(preynames[k],nrow(pred))
          pred<-data.frame(preysp,pred,rep("pred",nrow(pred))) #add the preyname column
          predicted<-rbind(predicted,pred)
        }#for prey species
        names(observed)<-c("species","x","y","dtype")
        names(predicted)<-c("species","x","y","dtype")
        pdata<-rbind(observed,predicted)
        yname<-pm$ylab[i]

        newplot<-plot_fitted_multipanel(pdata,xname,yname,speciesname,ptitle)
      }#plottype="availability

      if(pm$plottype[i]=="compare2xy" && pm$slot.y[i]=="R_Punt"){
        observed<-data.frame()
        predicted<-data.frame()
        #Put observed data in a dataframe
        xobs<-species@Nm
        yobs<-species@R
        xobs_Punt<-species@Nm
        yobs_Punt<-species@R_Punt
        method<-c(rep("Hilborn",length(x)),rep("Punt",length(x)))
        observed<-data.frame(method,c(xobs,xobs_Punt),c(yobs,yobs_Punt),rep("obs",length(x)))
        names(observed)<-c("method","x","y","dtype")

        #Get predicted values
        pred_R<-calc_predicted(species,"recruitment")
        pred_R<-data.frame(rep("Hilborn",nrow(pred_R)),pred_R,rep("pred",nrow(pred_R)))
        names(pred_R)<-c("method","x","y","dtype")
        pred_Punt<-calc_predicted(species,"R_Punt")
        pred_Punt<-data.frame(rep("Punt",nrow(pred_Punt)),pred_Punt,rep("pred",nrow(pred_Punt)))
        names(pred_Punt)<-c("method","x","y","dtype")
        predicted<-rbind(pred_R,pred_Punt)
        pdata<-rbind(observed,predicted)

        newplot<-plot_compare2fitted_multipanel(pdata,xname,yname,speciesname,ptitle)
      }#plottype=compare2xy

      #put the new plot in the growing plotlist
      plotlist<-c(plotlist,list(newplot))
      #name the list element
      names(plotlist)[length(plotlist)]<-pm$plotname[i]
    }#If the species is this obj_type
  }#for i in 1:nrow(pm)

  #Add the plot and dataframe to a list, and name them
  reslist<-list(data,plotlist)
  names(reslist)<-c("data","plots")
  reslist
} #compile_results()

#' Plot a time series (a utility function)
#'
#' @param data A data frame in which the X series to be plotted is in column 1.
#' @param xname A display name for the X-axis
#' @param yname A display name for the Y-axis
#' @param speciesname A display name for the species
#' @param ptitle A partial title
#' @param y_range Character; "include_0" to include zero, or "auto"
#'
#' @return A plot
#' @export
#' @keywords internal
#' @family results
makeplot<-function(pdata,xname,yname,speciesname,ptitle,y_range,series_type){
  if(y_range=="include_0"){
    if(series_type=="line"){
      tmp<-ggplot2::ggplot(pdata,aes(x=pdata[,1],y=pdata[,2])) +
        geom_line() + ylim(c(0,NA)) + xlab(xname) +
        ylab(yname) + ggtitle(paste(speciesname,ptitle))

    }
    if(series_type=="point"){
      tmp<-ggplot2::ggplot(pdata,aes(x=pdata[,1],y=pdata[,2])) +
        geom_point() + ylim(c(0,NA)) + xlab(xname) +
        ylab(yname) + ggtitle(paste(speciesname,ptitle))

    }
    if(series_type=="both" || series_type=="auto"){
      tmp<-ggplot2::ggplot(pdata,aes(x=pdata[,1],y=pdata[,2])) +
        geom_point() + geom_line() + ylim(c(0,NA)) + xlab(xname) +
        ylab(yname) + ggtitle(paste(speciesname,ptitle))
    }
  } else{
    tmp<-ggplot2::ggplot(pdata,aes(x=pdata[,1],y=pdata[,2])) +
      geom_point() + geom_line() + xlab(xname) + ylab(yname) +
      ggtitle(paste(speciesname,ptitle))
  }
  newname<-paste0("p.",yname)
  assign(newname,tmp)
  tmp
}

#' Sum a by-timestep series by years and plot it.
#'
#' A utility function for plotting.
#'
#' @param data A small data frame with the by-timestep series in it
#' @param item The name of the time series (matches a colname in data)
#' @param xname The name of the x series
#' @param speciesname The name of the species
#' @param nYears The number of years in the time series
#' @param nSteps The number of timesteps per year
#' @return A plot of values of a vector, summed by year.
#' @export
#' @keywords internal
#' @family results
plot_summed_by_year<-function(nYears,nSteps,pdata,xname,yname,...){
  #For species that use timesteps but where we want to plot values by years, make
  #a vector of the year number, e.g.,: 1 1 1 1 2 2 2 2 3 3 3 3 (here, nSteps=4).
  ts<-1:(nYears * nSteps)
  yrindex<-(ts + (nSteps-1)) %/% nSteps
  years<-1:nYears

  #Sum the y-values by year
  yrvals<-tapply(pdata[,2],yrindex,sum,na.rm=T)
  yrsumdata<-data.frame(years,yrvals)

  #Plot the new dataframe
  newplot<-makeplot(yrsumdata,xname,yname,...)
  newplot
}

#' Get every nth value of a time series
#'
#' A utility function that returns every nth intra-annual step from a time
#' series that includes intra-annual steps.
#'
#' @param nSteps The number of timesteps per year
#' @param pdata A small data frame with the by-timestep series in it
#' @param xname The name of the x series
#' @param yname The name of the y series
#' @return A data frame (x,y) of every nth value of a time series
#' @export
#' @keywords internal
#' @family results
get_every_nth<-function(nSteps,pdata,xname,yname){
  #Take every nth value of the item to plot (here, the 1st timestep of each year)
  x_nth<-pdata[seq(1,nrow(pdata),nSteps),1]
  y_nth<-pdata[seq(1,nrow(pdata),nSteps),2]

  #put it into a dataframe and return it
  nth_data<-data.frame(x_nth,y_nth)
  names(nth_data)<-c(xname,yname)
  nth_data
}

#' Plot a fitted function (observed vs. expected)
#'
#' @param fitted_data A dataframe of expected data
#' @param pdata A dataframe of observed data
#' @param xname X-axis label
#' @param yname Y-axis label
#' @param speciesname Species name (for title)
#' @param ptitle Partial title
#' @return A plot
#' @export
#' @keywords internal
#' @family results
plot_fitted<-function(fitted_data,pdata,xname,yname,speciesname,ptitle,...){
  tmp<-ggplot2::ggplot(pdata,aes(x=pdata[,1],y=pdata[,2])) +
    geom_point() + geom_line(data=fitted_data,aes(x=fitted_data[,1],y=fitted_data[,2])) +
    xlab(xname) + ylab(yname) + ggtitle(paste(speciesname,ptitle))
  tmp
}

#' Plot a multi-panel fitted function (observed vs. expected)
#'
#' @param observed.df A dataframe of observed data
#' @param predicted.df A dataframe of predicted data
#' @param xname X-axis label
#' @param yname Y-axis label
#' @param speciesname Species name (for title)
#' @param ptitle Partial title
#' @return A plot
#' @export
#' @keywords internal
#' @family results
plot_fitted_multipanel<-function(pdata,xname,yname,speciesname,ptitle,...){
  tmp<-ggplot2::ggplot(pdata,aes(x=x,y=y,group=dtype)) +
    geom_point(data=pdata[pdata$dtype=="obs",]) +
    geom_line(data=pdata[pdata$dtype=="pred",]) +
    facet_wrap(~species,scales="free") +
    xlab(xname) + ylab(yname) + ggtitle(paste(speciesname,ptitle))
  tmp
}

plot_compare2fitted_multipanel<-function(pdata,xname,yname,speciesname,ptitle,...){
  tmp<-ggplot2::ggplot(pdata,aes(x=x,y=y,colour=method)) +
    geom_point(data=pdata[pdata$dtype=="obs",]) +
    geom_line(data=pdata[pdata$dtype=="pred",]) +
    xlab(xname) + ylab(yname) + ggtitle(paste(speciesname,ptitle))
  tmp
}



#' Plot a function (expected values only)
#'
#' @param fndata A dataframe of expected data
#' @param xname X-axis label
#' @param yname Y-axis label
#' @param speciesname Species name (for title)
#' @param ptitle Partial title
#' @return A plot
#' @export
#' @keywords internal
#' @family results
plot_function_only<-function(fndata,xname,yname,speciesname,ptitle){
  tmp<-ggplot2::ggplot(fndata,aes(x=fndata[,1],y=fndata[,2])) +
    geom_line() + xlab(xname) + ylab(yname) + ggtitle(paste(speciesname,ptitle))
  tmp
}


#' Calculate predicted values
#'
#' A utility function for plotting, which calculates various kinds of predicted curves
#' for overlay in plots.
#' @param species A species object (whatever is being plotted)
#' @param predicted Character; the name of a predicted series
#' @param preyname Character; the name of a prey species (only used for availabililty)
#'
#' @return A data frame of predicted data (x and y; 100 rows)
#' @export
#' @keywords internal
#' @family results
calc_predicted<-function(species,predicted,preyname){
  if(predicted=="recruitment" && is(species,"Prey") && !is(species,"Other_forage")){
    alpha<-species@recruitment$alpha
    beta<-species@recruitment$beta
    bmax<-max(species@B_sp,na.rm=T)
    x<-seq(0,bmax,bmax/100)
    y<-sapply(x,function(b) calc_Ricker(alpha,beta,b))
  }
  if(predicted=="recruitment" && is(species,"Other_forage")){
    maxb<-max(species@B_sp,na.rm=T)
    x<-seq(0,maxb,maxb/100)
    y<-rep(species@R_baseyear,length(x))
  }
  if(predicted=="recruitment" && is(species,"Predator")){
    p<-species@recruitment$BHp
    c<-species@recruitment$BHc
    maxnm<-max(species@Nm,na.rm=T)
    x<-seq(0,maxnm,maxnm/100)
    y<-sapply(x,function(n) calc_Bev_Holt(n,p,c))
  }
  if(predicted=="prey_effect" && is(species,"Predator")){
    maxd<-max(species@Dy,na.rm=T)
    x<-seq(0,maxd,maxd/100)
    y<-sapply(x,function(d) calc_p(species,d))
  }
  #For plotting R_Punt, phi=1 (no effect of prey)
  if(predicted=="R_Punt" && is(species,"Predator")){
    maxd<-max(species@Nm,na.rm=T)
    x<-seq(0,maxd,maxd/100)
    y<-numeric(length(x))
    K1P<-species@K1P
    #Estimate N1P from Nm, via ratios of survival
    am<-species@am #age at maturity
    nAges<-species@P + 1
    S1<-exp(-(species@Mv[2])) #age-1 survival
    Sx<-sapply(1:(nAges-1),function(x) S1^(x-1))
    for(i in 1:length(x)){
      Nm<-x[i]
      Stot<-0
      Sam<-0
      for (j in 1:(nAges - 1)){
        Stot<-Stot + Sx[j]
        if(j >= am){
          Sam<-Sam + Sx[j]
        }
      }#for j
      N1P<-Nm * (Stot/Sam)
      y[i]<-calc_predRecruits_Punt(species,Nm,1,N1P,K1P)
    }
  }

  #For availability, create a data frame for the prey in question
  #(availability parameters vary by predator/prey combinations)
  if(predicted=="availability" && is(species,"Predator")){
    ap<-species@prey_availability[[preyname]]$ap
    ac<-species@prey_availability[[preyname]]$ac
    maxb<-max(species@prey_availability[[preyname]]$B_agepref,na.rm=T)
    x<-seq(0,maxb,maxb/100)
    y<-sapply(x,function(b) calc_availability(ap,ac,b))
  }

  #put the results in a data frame and add a column called "id" labeling the values as predicted
  df<-data.frame(x,y)
  df
}

#' Make a 4-panel plot
#'
#' Make a 4-panel plot that shows the abundance of the key species (spawning
#' biomass for prey; mature adults for predators) .
#'
#' @param prey A list containing the Prey objects
#' @param predator A list containing the Predator objects
#'
#' @return A single 4-panel plot showing the abundance of 4 key species
#' @export
#'
#' @keywords internal
make_4panel<-function(prey,predator,scenario,nSteps){
  pl_list<-c("Sardine","Anchovy","Other_forage","Sealion","Pelican")
  data4p<-data.frame()
  #Create a data frame of B_sp for prey, and Nm for predators
  for(preyx in prey){
    if(preyx@name %in% pl_list){
      names4p<-preyx@name
      #pick out every timestep 1 if necessary
      if(preyx@uses_timestep){
        xdata<-1:length(preyx@B_sp)
        ydata<-preyx@B_sp
        temp<-data.frame(xdata,ydata)
        temp<-get_every_nth(nSteps,temp,"Year","B_sp")
        newd<-data.frame(rep(preyx@name,nrow(temp)),temp)
      }else{
        newd<-data.frame(rep(preyx@name,length(preyx@B1PperB0)),1:length(preyx@B1PperB0),preyx@B1PperB0)
      }
      names(newd)<-c("species","year","spawners")

      if(nrow(data4p)==0){
        data4p<-newd
      } else{
        data4p<-rbind(data4p,newd)
      }
    }#preyname in pl_list
  }#preyx in prey

  for(predx in predator){
    if(predx@name %in% pl_list){
      names4p<-predx@name
      newd<-data.frame(rep(predx@name,length(predx@NperK)),1:length(predx@NperK),predx@NperK)
      names(newd)<-c("species","year","spawners")
      if(nrow(data4p)==0){
        data4p<-newd
      } else{
        data4p<-rbind(data4p,newd)
      }
    }#predname in pl_list
  }#predx in predator

  #Plot the data and return the plot
  p4<-ggplot(data4p,aes(x=year,y=spawners)) +
    geom_line() + facet_wrap(~species,scales="free",ncol=2) + ylim(c(0,NA)) +
    xlab("Year") + ylab("B1+, B1+/B0, N, or N/K") + ggtitle(paste("Scenario:",scenario))
  p4
}

#show all plots from a result object
#' Show all plots in a result object
#'
#' @param res A result object
#'
#' @return all the plots (in a list, but they display by default)
#' @export
#'
allplots<-function(res){
  species<-names(res$species)
  plist<-list()
  for(i in 1:length(species)){
    plist<-c(plist,list(res[[species[i]]]$plots))
  }
  plist
}

#' Show one type of plot for all species
#'
#' @param res Results object
#' @param plotname The name of a plot type (e.g., "B1P" or "R")
#'
#' @return Returns a list of plots (also displays them by default in RStudio)
#' @export
#'
#' @examples
#' \dontrun{
#'   justplot(res,"R") #plots the "R" (recruitment) plots for each species
#'   #To find plotnames, look at the plotnames for the species you're interested
#'   in, for example;
#'   names(res$Anchovy$plots)
#'   #or, to see what a particular plot looks like:
#'   res$Anchovy$plots$R
#' }
justplot<-function(res,plotname){
  species<-names(res$species)
  plist<-list()
  for(i in 1:length(species)){
    if(exists(plotname,where=res[[species[i]]]$plots)){
      plist<-c(plist,list(res[[species[i]]]$plots[[plotname]]))
    }
  }
  plist
}

make_panel_plot<-function(res){
  scenario<-res$input_params@scenario
  p1<-res$Anchovy$plots$B1P
  p1$labels$title<-"Anchovy"
  p2<-res$Sardine$plots$B1P
  p2$labels$title<-"Sardine"
  p3<-res$Other_forage$plots$B1PperB0
  p3$labels$title<-"Other_forage"
  p4<-res$Sealion$plots$NperK
  p4$labels$title<-"Sea Lion"
  p5<-res$Pelican$plots$NperK
  p5$labels$title<-"Pelican"

  grid.arrange(p1,p2,p3,p4,p5,ncol=2,top=paste("Scenario =",scenario))
}#splot

getStats<-function(res){
  scenario<-res$input_params@scenario
  nSteps<-res$input_params@nSteps
  C<-res$objects$Sardine@Cy
  #Get the B1P estimates for sardine and anchovy (they are by year)
  BsarperB0<-res$objects$Sardine@B1PperB0
  BanchperB0<-res$objects$Anchovy@B1PperB0
  #Get the actual sardine biomass, B1P, which is by step, so we select only step 1
  Bsar<-res$objects$Sardine@B1P
  Bsar<-Bsar[seq(1,length(Bsar),nSteps)]
  #Remove missing values (better to do this here, because you have to remove them from both numerator and denominator later)
  Bsar<-Bsar[!is.na(Bsar)]
  #And the data for predators
  NperK_pel<-res$objects$Pelican@NperK
  NperK_seal<-res$objects$Sealion@NperK
  #remove missing values
  NperK_pel<-NperK_pel[!is.na(NperK_pel)]
  NperK_seal<-NperK_seal[!is.na(NperK_seal)]


  #Assemble the pieces for the table
  meanC<-mean(C)
  pClt50<-length(C[C<50])/length(C)
  meanBsarperB0<-mean(BsarperB0,na.rm=T)
  meanBanchperB0<-mean(BanchperB0,na.rm=T)
  pBsargt400<-length(Bsar[Bsar > 400])/length(Bsar)
  pBsarlt150<-length(Bsar[Bsar < 150])/length(Bsar)
  meanPelNperK<-mean(NperK_pel)
  plt5_pel<-length(NperK_pel[NperK_pel < .5])/length(NperK_pel)
  plt1_pel<-length(NperK_pel[NperK_pel < .1])/length(NperK_pel)
  meanSealNperK<-mean(NperK_seal)
  plt5_seal<-length(NperK_seal[NperK_seal < .5])/length(NperK_seal)
  plt1_seal<-length(NperK_seal[NperK_seal < .1])/length(NperK_seal)

  #Stick them into a 1-row data frame
  nlabels<-c("C_mean","p(Catch < 50kt)","mean B/B0_sard","mean B/B0_anch","p(B_sard > 400kt)",
           "p(B_sard < 150kt)","mean N/K_pel","p(N_pel < .5K)","p(N_pel < .1K)",
           "mean N/K_seal","p(N_seal < .5K)","p(N_seal < .1K)")
  numbers<-data.frame(meanC,pClt50,meanBsarperB0,meanBanchperB0,pBsargt400,pBsarlt150,meanPelNperK,
                      plt5_pel,plt1_pel,meanSealNperK,plt5_seal,plt1_seal)

  stats<-data.frame(rbind(nlabels,numbers))
  names(stats)<-c("C","C50","B_sar","B_anch","Bsar400","Bsar150","Npel",
                  "Npel5","Npel1","Nseal","Nseal5","Nseal1")
  #return it
  stats
}

