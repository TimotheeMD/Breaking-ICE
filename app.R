library(shiny)
library(gridlayout)
library(bslib)
library(IPDfromKM)
library(rms)
library(shinydashboard)
# using survminer https://rpkgs.datanovia.com/survminer/
library(survminer)
# library(plotly)

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
  
  return(list(xyz=list(x=x,y=y,z=z), est_E_OS=est_E_OS, est_C_OS=est_C_OS))

}

original <- getOriginal()
originalFit <- survfit(Surv(x,y)~z, data = original$xyz)

newData <- function(startTime=0,
                    endTime=3,
                    Eperc=15,
                    Cperc=15,
                    seed=123,
                    original=getOriginal()
                    ){
  newdata <- rbind(original$est_E_OS$IPD, original$est_C_OS$IPD)
  
  # Set the seed for reproducibility
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  # Parameters for the experimental arm (censored for additional toxicity, more likely to present the event)
  # here we will take the censored patients from 0 to 3 months
  # then randomly pick a percentage of those (15%)
  # and change their status between being censored to having the event
  # startTime <- 0 # start time
  # endTime <- 3 # end time
  # Eperc <- 15 # Percentage of "0" status to be changed to "1"
  
  # Parameters for the control arm (censored for patient disappointement, less likely to present the event)
  # here we will take the censored patients from 0 to 3 months
  # then randomly pick a percentage of those (15%)
  # and change the time of censoring to be the same as the last patient being censored in this arm
  # startTime <- 0 # start time
  # endTime <- 3 # end time
  # Cperc <- 15 # Percentage of times to be changed
  
  # Subset rows with treat = 1 and time within [startTime, endTime]
  subset_treat_1 <- newdata[newdata$treat == 1 & newdata$time >= startTime & newdata$time <= endTime,]
  
  # Further subset those with status = 0
  subset_treat_1_status_0 <- subset_treat_1[subset_treat_1$status == 0,]
  
  # Calculate number of rows to change based on percentage
  e_change <- round(nrow(subset_treat_1_status_0) * (Eperc / 100.0))
  
  # Randomly select rows to change
  rows_to_change <- sample(nrow(subset_treat_1_status_0), e_change)
  
  # Change status to 1 for these rows
  subset_treat_1_status_0$status[rows_to_change] <- 1
  
  # Replace the original rows in newdata
  newdata[newdata$treat == 1 & newdata$time >= startTime & newdata$time <= endTime & newdata$status == 0,] <- subset_treat_1_status_0
  
  # Subset rows with treat = 0 and time within [startTime, endTime]
  subset_treat_0 <- newdata[newdata$treat == 0 & newdata$time >= startTime & newdata$time <= endTime,]
  
  # Further subset those with status = 0
  subset_treat_0_status_0 <- subset_treat_0[subset_treat_0$status == 0,]
  
  # Find the longest time among those with status = 0 and treat = 0
  longest_time <- max(newdata$time[newdata$status == 0 & newdata$treat == 0])
  
  # Calculate number of rows to change based on percentage
  c_change <- round(nrow(subset_treat_0_status_0) * (Cperc / 100.0))
  
  # Randomly select rows to change
  rows_to_change <- sample(nrow(subset_treat_0_status_0), c_change)
  
  # Change time for these rows
  subset_treat_0_status_0$time[rows_to_change] <- longest_time
  
  # Replace the original rows in newdata
  newdata[newdata$treat == 0 & newdata$time >= startTime & newdata$time <= endTime & newdata$status == 0,] <- subset_treat_0_status_0
  
  ## the number of patients with data being change in each arm : 
  c_change
  e_change
  
  #### SURVIVAL ANALYIS
  ##preparation
  x <- newdata$time
  y <- newdata$status
  z <- newdata$treat
  
  return( list(xyz=list(x=x,y=y,z=z)))
}

# Define UI - version 0004

# Adding an image - code from Maria
# text_disclaimer_outdooractive = 
#   c("<a href='//www.outdooractive.com' title='This website uses technology and content from the Outdooractive Platform.'>
#     <img src='//res.oastatic.com/partner/outdooractive-black.png' alt='outdooractive' />
#     This website uses technology and content from the Outdooractive Platform.</a>",

ui <- fluidPage(

  titlePanel("Breaking ICE"),

      box(
        title="Trial Data",  br(),
        selectInput(
          inputId = "cut",
          label = "Select a trial",
          choices = list(
            "CONTACT-02-PFS" = "CONTACT-02-PFS",
            "CONTACT-02-OS" = "CONTACT-02-OS"
          ),
          selected = "CONTACT-02-PFS",
          width = "100%"
        ),
        textInput("text", "Or load your own trial data:")
      )
        ,

        sidebarLayout(

          position = c("right"),
          #fluid = TRUE,

        sidebarPanel(
          h3("Settings"),   "Change the parameters to see how it affects the survival analysis on the left",br(),
          sliderInput(
            inputId = "time_start",
            label = "Start Time",
            min = 0,
            max = ceiling(max(original$xyz$x)),
            value = 0,
            width = "100%"
          ),
          sliderInput(
            inputId = "time_end",
            label = "End Time",
            min = 0,
            max = ceiling(max(original$xyz$x)),
            value = 3,
            width = "100%"
          ),
          numericInput(
            inputId = "Eperc",
            label = "Eperc",
            value = 15,
            max = 100
          ),
          numericInput(
            inputId = "Cperc",
            label = "Cperc",
            value = 15,
            max = 100
          )
          #,
          #textInput("text", "Text input:")
        ),

        mainPanel(
          h3("KM curve"),
          plotOutput("plot"),
          # verbatimTextOutput("statistics")
          box(
            title="Results",  br(),
            "exp(z):", textOutput("statistics_z"),
            "p-level:", textOutput("statistics_plevel"),
            "interval:", textOutput("statistics_interval")
            # textInput("text", "Or load your own trial data:")
          )
        ),

      ),
)

server <- function(input, output) {
  
  newXYZ <- reactive({
    newData(input$time_start,input$time_end, input$Eperc, input$Cperc)$xyz
  })
  newFit <- reactive({
    survfit(Surv(x,y)~z, data=newXYZ())
  })
  newStatistics <- reactive({
    fit <- coxph(Surv(x,y)~z, data=newXYZ())
    s <- summary(fit)
    s
  })

  # output$plot <- renderPlot({
  #   ##then Kaplan-Meier Curve
  #   # par(mar = c(1, 1, 1, 1), xaxs = "i", yaxs = "i")
  #   plot(originalFit, xlim = c(0, 30), ylim = c(0, 1),
  #        col=c("#66CC66", "#60C060"), lwd=3, xaxt='n', bty = "l",
  #        las = 1, cex.axis = 1.5, tcl  = -0.5)
  #   axis(side=1, at=seq(0, 30, by=3), cex.axis = 1.5)
  # 
  #   ##then Kaplan-Meier Curve
  #   # par(mar = c(1, 1, 1, 1), xaxs = "i", yaxs = "i")
  #   lines(newFit(), xlim = c(0, 30), ylim = c(0, 1),
  #        col=c("red", "#6699CC"), lwd=3, xaxt='n', bty = "l",
  #        las = 1, cex.axis = 1.5, tcl  = -0.5)
  #   # axis(side=1, at=seq(0, 30, by=3), cex.axis = 1.5)
  #   abline(h = 0, v = 0)
  # 
  # 
  #   ## then Cox, HR and p-value
  #   # summary(coxph(Surv(x,y)~z))
  # 
  #   # legend
  #   legend("topright", inset=.02,
  #          legend=c("Original","Original", "Docetaxel", "Sotorasib"),
  #          col=c("#66CC66", "#60C060","red", "#6699CC"),
  #          title="Group",
  #          box.lty=0, lty=1:2, cex=0.8)
  # 
  # 
  # })


  newStatistics <- reactive({
    fit <- coxph(Surv(x,y)~z, data=newXYZ())
    s <- summary(fit)
    list(
      z=s$coefficients[2],
      plevel=s$coefficients[5],
      interval=s$conf.int[3:4]
    )
  })
  output$statistics_z <- renderText({newStatistics()$z})
  output$statistics_plevel <- renderText({newStatistics()$plevel})
  output$statistics_interval <- renderText({newStatistics()$interval})

  output$plot <- renderPlot({
    
    fit <- survfit(Surv(x,y)~z, data = original$xyz)

    # ggsurvplot(fit, data = original$xyz)


  g <- ggsurvplot(
    fit,
    data = original$xyz,
    size = 1,                 # change line size
    palette =
      c("red", "#6699CC"),    # custom color palettes
    conf.int = FALSE,          # Add confidence interval
    pval = TRUE,              # Add p-value
    risk.table = TRUE,        # Add risk table
    risk.table.col = "strata",  # Risk table color by groups
    legend.labs =
      c("Docetaxel - Control", "Sotorasib - Experiment"),
    # Change legend labels
    risk.table.height = 0.25, # Useful to change when you have multiple groups
    xlab = "Time in months",   # customize X axis label.
    break.time.by = 3,     # break X axis in time intervals by 3.
    ggtheme = theme_bw(),      # Change ggplot2 theme theme_light(), 
    ncensor.plot = FALSE,      # plot the number of censored subjects at time t
    ncensor.plot.height = 0.25
  )
  print(class(g))
  print(g)

  })

}

shinyApp(ui, server)
