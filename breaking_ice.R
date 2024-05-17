library(IPDfromKM)
library(rms)
# using survminer https://rpkgs.datanovia.com/survminer/
library(survminer)
# library(plotly)


make_levels <- function(df){
  df$z <- factor(df$z)
  levels(df$z) <- c('control','experiment')
  df$z <- forcats::fct_relevel(df$z, 'experiment')
  return(df)
}

getOriginal <- function(){
  
  ## First step, go on digitilzese and for each arm (important, use 100 as the maximum in y axis), you export an csv files on your Desktop (EXP for the experimental arm, and CON for the control arm)
  
  ## first the experimental arm
  ##read the data
  E <- readxl::read_xlsx(path="./EXP.xlsx" ) #, sep = ",", header=T)
  # E <- EXP
  colnames(E) <-c("time", "survival probability")
  
  ## the same with the control arm (we called C)
  # C <- read.csv(file="/Users/kailash/Desktop/CON.csv", sep = ",", header=T)#sep = , I don't know why, but it's important
  C <- readxl::read_xlsx("./CON1.xlsx")
  # C <- CON
  colnames(C) <-c("time", "survival probability")
  
  ##time risk, separate by comma (no space allowed)
  trisk <- c(0,3,6,9,12,15,18,21,24,27,30)
  
  ## number at risk experimental arm
  nrisk.E <- c(200,135,89,45,17,8,4,2,0,0,0)
  
  ## number at risk control arm 
  nrisk.C <- c(200,98,57,35,16,8,6,5,3,1,0)
  
  return( processCurves(E, trisk, nrisk.E, C, nrisk.C) )
}

readCurvesFromExcel <- function(filename){
  sheets <- readxl::excel_sheets(filename)
  times <- readxl::read_xlsx(filename, sheet='Times')
  E <- readxl::read_xlsx(filename, sheet='EXP')
  C <- readxl::read_xlsx(filename, sheet='CON')
  
  return( processCurves(E, times$trisk, times$nrisk.E, C, times$nrisk.C) )
}

# st <- readExcel(filename='./Trial Data _ CONTACT-02-PFS.xlsx')

processCurves <- function(E, trisk, nrisk.E, C, nrisk.C) {
  ## preprocess the data
  pre_E <- preprocess(dat=E, trisk=trisk, nrisk=nrisk.E, maxy=100)
  pre_C <- preprocess(dat=C, trisk=trisk, nrisk=nrisk.C, maxy=100)
  
  ## then individual patient data = IPD, for experimental treat is 1, for control group treat is 0
  est_E_OS <- getIPD(prep=pre_E, armID=1, tot.events = NULL)
  est_C_OS <- getIPD(prep=pre_C, armID=0, tot.events = NULL)
  
  ### important, you can isolate the IPD part (time, status, treat) from est_E_OS and est_C_OS
  # est_E_OS$IPD
  # est_C_OS$IPD
  
  ### the summary function allow to have the data of events, censoring, etc. 
  # summary(est_E_OS)
  # summary (est_C_OS)
  
  ## here you can create a KM curve based on the novel IPD data generated based on the data that have been digitilzed
  ## you can also recapitulate the Cox analysis. 
  surv_data<- rbind(est_E_OS$IPD, est_C_OS$IPD)
  names(surv_data) <- c("time", "status", "arm")
  x <- surv_data$time
  y <- surv_data$status
  z <- surv_data$arm
  
  xyz <- data.frame(x=x,y=y,z=z)
  xyz <- make_levels(xyz)

  return(list(xyz=xyz, est_E_OS=est_E_OS, est_C_OS=est_C_OS, trisk=trisk))
}

original <- getOriginal()
# original <- readCurvesFromExcel(filename='./Trial Data _ CONTACT-02-PFS.xlsx')

originalFit <- survfit(Surv(x,y)~z, data = original$xyz)

simulateToxicity <- function(newdata, time_frame, perc) {
  # Parameters for the experimental arm (censored for additional toxicity, more likely to present the event)
  # here we will take the censored patients from 0 to 3 months
  # then randomly pick a percentage of those (15%)
  # and change their status between being censored to having the event
  # startTime <- 0 # start time
  # endTime <- 3 # end time
  # Eperc <- 15 # Percentage of "0" status to be changed to "1"
  
  # Subset rows with treat = 1 and time within [startTime, endTime]
  subset_treat_1 <- newdata[newdata$time >= time_frame[1] & newdata$time <= time_frame[2],]
  
  # Further subset those with status = 0
  subset_treat_1_status_0 <- subset_treat_1[subset_treat_1$status == 0,]
  
  # Calculate number of rows to change based on percentage
  e_change <- round(nrow(subset_treat_1_status_0) * (perc / 100.0))
  
  # Randomly select rows to change
  rows_to_change <- sample(nrow(subset_treat_1_status_0), e_change)
  
  # Change status to 1 for these rows
  subset_treat_1_status_0$status[rows_to_change] <- 1
  
  # Replace the original rows in newdata
  newdata[newdata$time >= time_frame[1] & newdata$time <= time_frame[2] & newdata$status == 0,] <- subset_treat_1_status_0
  
  return(newdata)
}

simulateDisappointment <- function(newdata, time_frame, perc, filter_treat=0) {
  # Parameters for the control arm (censored for patient disappointement, less likely to present the event)
  # here we will take the censored patients from 0 to 3 months
  # then randomly pick a percentage of those (15%)
  # and change the time of censoring to be the same as the last patient being censored in this arm
  # startTime <- 0 # start time
  # endTime <- 3 # end time
  # Cperc <- 15 # Percentage of times to be changed
  
  # Subset rows with treat = 0 and time within [startTime, endTime]
  subset_treat_0 <- newdata[newdata$time >= time_frame[1] & newdata$time <= time_frame[2],]
  
  # Further subset those with status = 0
  subset_treat_0_status_0 <- subset_treat_0[subset_treat_0$status == 0,]
  
  # Find the longest time among those with status = 0
  longest_time <- max(newdata$time[newdata$status == 0])

  # Calculate number of rows to change based on percentage
  c_change <- round(nrow(subset_treat_0_status_0) * (perc / 100.0))
  
  # Randomly select rows to change
  rows_to_change <- sample(nrow(subset_treat_0_status_0), c_change)
  
  # Change time for these rows
  subset_treat_0_status_0$time[rows_to_change] <- longest_time
  
  # Replace the original rows in newdata
  newdata[newdata$time >= time_frame[1] & newdata$time <= time_frame[2] & newdata$status == 0,] <- subset_treat_0_status_0
  
  return(newdata)
}

simulate <- function(df, modellingType, ...){
  if(modellingType == 'Toxicity'){
    return( simulateToxicity(df, ...) )
  }else if(modellingType == 'Disappointment'){
    return( simulateDisappointment(df, ...) )
  }else{
    error(paste0('modellingType not known: ',modellingType,' - has to be "Toxicity" or "Disappointment"'))
  }
}


newData <- function(time_frame_experimental=c(0,3),
                    time_frame_control=c(0,3),
                    Eperc=15,
                    Cperc=15,
                    modellingTypeExperimental='Toxicity',
                    modellingTypeControl='Disappointment',
                    seed=123,
                    original=getOriginal()
){

  # Set the seed for reproducibility
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  newdata <- rbind(
    simulate(original$est_E_OS$IPD, modellingTypeExperimental, time_frame_experimental, Eperc),
    simulate(original$est_C_OS$IPD, modellingTypeControl, time_frame_control, Cperc)
  )
  
  
  # ## the number of patients with data being change in each arm : 
  # c_change
  # e_change
  
  #### SURVIVAL ANALYIS
  ##preparation
  x <- newdata$time
  y <- newdata$status
  z <- newdata$treat

  xyz <- data.frame(x=x,y=y,z=z)
  xyz <- make_levels(xyz)
  
  return( list(xyz=xyz))
}

plotCurves <- function(original, originalFit, newXYZ, newFit, colours) {
  if(missing(colours)){
    colours <- c("#6699CC", "red", "green","orange")
  }
  df <- bind_rows(original=original$xyz, simulated=newXYZ, .id='group')
  print(levels(df$z))

  
  # g <- ggsurvplot(
  g <- ggsurvplot_combine(
    fit = list(original=originalFit, simulated=newFit),
    data = df,
    palette = colours,
      # c("#6699CC", "red", "green","orange"),    # custom color palettes
    legend.labs =
      c("Original - Experiment", "Original - Control", 'Sensitivity - Experiment', 'Sensitivity - Control'),
    # palette =
    #   c("red", "#6699CC","orange", "green"),    # custom color palettes
    # legend.labs =
    #   c("Original - Control", "Original - Experiment", 'Simulated - Control', 'Simulated - Experiment'),
      # c("Docetaxel - Control", "Sotorasib - Experiment", 'simulated - control', 'simulated - Experiment'),
    # originalFit,
    # data = original$xyz,
    # newFit(),
    # data = newXYZ(),
    size = 1,                 # change line size
    # palette =
    #   c("red", "#6699CC"),    # custom color palettes
    conf.int = FALSE,          # Add confidence interval
    pval = FALSE,              # Add p-value
    risk.table = TRUE,        # Add risk table
    # risk.table = "nrisk_cumcensor",        # Add risk table
    risk.table.col = "strata",  # Risk table color by groups
    # legend.labs =
    #   c("Docetaxel - Control", "Sotorasib - Experiment"),
    # Change legend labels
    risk.table.height = 0.25, # Useful to change when you have multiple groups
    xlab = "Time in months",   # customize X axis label.
    break.time.by = 3,     # break X axis in time intervals by 3.
    ggtheme = theme_minimal(),      # Change ggplot2 theme theme_light(), 
    ncensor.plot = FALSE,      # plot the number of censored subjects at time t
    ncensor.plot.height = 0.25
  )
}

calculateRiskMat <- function(xyz, trisk){
  
  xyz %>% mutate(
    t=cut(x, breaks=trisk, include.lowest=TRUE, labels=trisk[-length(trisk)])
    # nrisk=
    ) %>% 
    group_by(z, t) %>% summarise(censor.hat=sum(1-y), event.hat=sum(y)) %>% 
    arrange(t) %>% mutate(
      cumsum.censor=cumsum(coalesce(lag(censor.hat),0)),
      cumsum.event=cumsum(coalesce(lag(event.hat),0)),
      nrisk=sum(censor.hat)+sum(event.hat) - cumsum.censor - cumsum.event
    ) %>% ungroup() %>% arrange(z, t)
}
# calculateRiskMat(newXYZ$xyz, original$trisk)

calculateCensorPerc <- function(original, newData, roundTo=1){
  sim_riskmat <- calculateRiskMat(newXYZ$xyz, original$trisk)
  dataframes <- list(
    "Original - Experiment"=original$est_E_OS$riskmat,
    "Original - Control"=original$est_C_OS$riskmat,
    "Sensitivity - Experiment"=sim_riskmat %>% filter(z=='experiment'),
    "Sensitivity - Control"=sim_riskmat %>% filter(z=='control')
  )
  maxN <- max(sapply(dataframes, nrow))
  ret <- bind_cols(trisk=original$trisk[1:maxN], lapply(dataframes,
    function(df){
      res <- df %>% mutate(censor.perc=round(censor.hat / nrisk * 100,roundTo)) %>% pull(censor.perc)
      length(res) <- maxN
      # res[is.na(res)] <- max(res, na.rm=TRUE)
      res
    }
  ))
  return(ret)
}

# calculateCensorPerc(original, newData)

quickRun <- function(){
  print("Quick Run")
  newXYZ <- newData(
    modellingTypeExperimental = 'Disappointment',
    modellingTypeControl = 'Toxicity',
  )
  newFit <- survfit(Surv(x,y)~z, data=newXYZ$xyz)
  # tables <- ggsurvtable(newFit, data=newXYZ$xyz, bre)$risk.table$data
  g <- plotCurves(original, originalFit, newXYZ, newFit)
  g
}

## uncomment for quick development
print(quickRun())

# (simulateToxicity(original$est_C_OS$IPD, time_frame=c(0,3), perc=90 )) %>% head

