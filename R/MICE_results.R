#' Compile results and generate graphs.
#'
#' Creates graphs of objects in Species slots, taking instructions from the
#' plotlist.csv file.  Note that since a dataframe cannot be a ragged array, the
#' number of rows in the data frame returned by the function is either nYears or
#' tTimes, but certain items (catch, recruits) are plotted by year, even for
#' species that are tracked by timestep.
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
compile_results<-function(species,plot_master,nYears,nSteps){
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

  #Set up two alternative X axes (nYears, and nTimes)
  if(species@uses_timestep){
    x<-1:(nSteps * nYears)
    xname<-"Timestep" #label for plot X-axis
  } else{
    x<-1:nYears
    xname<-"Year" #label for plot X-axis
  }

  #Initialize the dataframe (to which Y columns will be added)
  data<-data.frame(x)
  names(data)<-xname

  #Step through the plotmaster dataframe; make plots and add data to the result
  #dataframe as instructed.
  for(i in 1:nrow(pm)){
    yname<-NA
    pdata<-NA
    ptitle<-NA
    if(is(species,pm$obj_type[i])){
      #Get y values and name
      yname<-pm$slot.y[i] #default; may be overridden
      vals<-slot(species,pm$slot.y[i])
      #Get x values if specified (overwrite x)
      if(!is.na(pm$slot.x[i])){
        x<-slot(species,pm$slot.x[i])
        xname<-pm$slot.x[i]
      }
      if(length(vals)==0){
        next
        #print(paste("There is no data for",yname,"in",speciesname)) #next
      } else{
        if(length(vals) != length(x)){
          stop(paste("Cannot plot",yname,"; The length of",yname,"does not match the length of x (",length(x),")"))
        }

        #Make a small dataframe for plotting
        pdata<-data.frame(x,vals)

        #Add the y timeseries to the results dataframe if in_dataframe=1 in plotmaster
        if(pm$in_dataframe[i]){
          data<-data.frame(data,vals)
          names(data)[ncol(data)]<-pm$slot.y[i] #name the column
        }
        #Get x and y labels if specified
        if(!is.na(pm$xlab[i])){
          xname<-pm$xlab[i]
        }
        if(!is.na(pm$ylab[i])){
          yname<-pm$ylab[i]
        }
        if(!is.na(pm$title[i])){
          ptitle<-pm$title[i]
        }
        #Get the y range option
        y_range<-pm$y_range[i]

        #Get the predicted option
        predicted<-pm$predicted[i]

        #Several plotting options:
        #Plot a time series
        if(pm$plottype[i]=="timeseries"){
          newplot<-makeplot(pdata,xname,yname,speciesname,ptitle,y_range)
        }
        #Plot a time series, summed by year
        if(pm$plottype[i]=="summed_by_year"){
          if(pm$slot.y[i]=="Cv" && species@uses_timestep && species@fished){
            xname<-"Year"
            newplot<-plot_summed_by_year(nYears,nSteps,pdata,xname,yname,speciesname,ptitle,y_range)
          }
        }
        #Plot a time series, but only step 1 of each year.  This may be called in
        #some cases for species that don't use timesteps, in which case it defaults
        #to the plain timeseries plot.
        if(pm$plottype[i]=="step1only"){
          if(species@uses_timestep){
            xname<-"Timestep"
            newplot<-plot_every_nth(nSteps,pdata,xname,yname,speciesname,ptitle,y_range)
          } else{
            newplot<-makeplot(pdata,xname,yname,speciesname,ptitle,y_range)
          }
        }
        #Plot two variables X and Y against each other
        if(pm$plottype[i]=="xy"){
          if(!is.na(predicted)){
            #get the predicted data
            fitted_data<-calc_predicted(species,predicted)
            names(fitted_data)<-names(pdata)
            #bind the predicted data to the observed data (just the way ggplot likes it!)
            newplot<-plot_fitted(fitted_data,pdata,xname,yname,speciesname,ptitle,y_range)
          } else{
            pdata<-pdata[order(pdata[,1]),] #sort the data by col 1
            newplot<-makeplot(pdata,xname,yname,speciesname,ptitle,y_range)
          }
        }#plottype=xy
        #Plot just a function (no observed data)
        if(pm$plottype[i]=="function_only"){
          if(pm$predicted[i]=="availability"){
            fndata<-calc_predicted(species,predicted)
            newplot<-plot_function_only(fndata,xname,yname,speciesname,ptitle)
          }
        }

        #put the new plot in the growing plotlist
        plotlist<-c(plotlist,list(newplot))
        #name the list element
        names(plotlist)[length(plotlist)]<-pm$plotname[i]
      } #if(length(vals) > 0)
    }#If the species is this obj_type
  }#for i in 1:nrow(pm)

  #For convenience, add a intra-year timestep to by-timestep dataframes
  if(species@uses_timestep){
    step<-rep(1:nSteps,nYears)
    data<-data.frame(step,data)
  }
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
makeplot<-function(pdata,xname,yname,speciesname,ptitle,y_range){
  if(y_range=="include_0"){
    tmp<-ggplot2::ggplot(pdata,aes(x=pdata[,1],y=pdata[,2])) +
      geom_point() + geom_line() + ylim(c(0,NA)) + xlab(xname) +
      ylab(yname) + ggtitle(paste(speciesname,ptitle))
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

#' Plot every nth value.
#'
#' A utility function for plotting only the first timestep in each year.
#'
#' @param data A small data frame with the by-timestep series in it
#' @param item The name of the time series (matches a colname in data)
#' @param xname The name of the x series
#' @param speciesname The name of the species
#' @param nSteps The number of timesteps per year
#' @return A plot of every nth value of a vector
#' @export
#' @keywords internal
#' @family results
plot_every_nth<-function(nSteps,pdata,xname,yname,...){
  #Take every nth value of the item to plot (here, the 1st timestep of each year)
  x_nth<-pdata[seq(1,nrow(pdata),nSteps),1]
  y_nth<-pdata[seq(1,nrow(pdata),nSteps),2]

  #put it into a dataframe and call makeplot
  nth_data<-data.frame(x_nth,y_nth)
  names(nth_data)<-c("x",yname)

  #Plot the new dataframe
  newplot<-makeplot(nth_data,xname,yname,...)
  newplot
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
#'
#' @return A data frame of predicted data (x and y; 100 rows)
#' @export
#' @keywords internal
#' @family results
calc_predicted<-function(species,predicted){
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
    y<-rep(species@recruitment$R_base,length(x))
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
  #Here I am using B1P as a substitute for available prey biomass; it is slightly
  #different because a predator's age preferences may skew the biomass, but this
  #at least shows the availability function on the right scale.
  if(predicted=="availability" && is(species,"Prey")){
    ap<-species@availability$ap
    ac<-species@availability$ac
    maxb<-max(species@B1P,na.rm=T)
    x<-seq(0,maxb,maxb/100)
    y<-sapply(x,function(b) calc_availability(ap,ac,b))
  }

  #put the results in a data frame and add a column called "id" labeling the values as predicted
  df<-data.frame(x,y)
  df
}

