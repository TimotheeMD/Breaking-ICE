library(shiny)
library(gridlayout)
library(bslib)
library(shinydashboard)

source("breaking_ice.R")

# Define UI - version 0004

# Adding an image - code from Maria
# text_disclaimer_outdooractive = 
#   c("<a href='//www.outdooractive.com' title='This website uses technology and content from the Outdooractive Platform.'>
#     <img src='//res.oastatic.com/partner/outdooractive-black.png' alt='outdooractive' />
#     This website uses technology and content from the Outdooractive Platform.</a>",

ui <- fluidPage(
  
  titlePanel("Breaking ICE"),

  fluidRow( 
    column (8, 
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
    fileInput("user_trial", "Or load your own trial data:", accept = c(".xlsx")),
  ))
  ,
     column (4, 
  box(
    title="Breaking-ICE",
    "Informative",br(),
    "Censoring",br(),
    "Exploration",br(),
  ))
  ),
  
  sidebarLayout(
    
    sidebarPanel(
      h3("Sensitivity analysis"),   
      "Change the parameters to see how it affects the survival analysis on the left",br(),
      # h4("Experimental arm"),
      h4(div("Experimental arm", style = "color: #6699CC")),
      sliderInput(
        inputId = "time_frame_experimental",
        label = "Time Frame",
        min = 0,
        max = ceiling(max(original$xyz$x)),
        value = c(0,3),
        width = "100%"
      ),
      numericInput(
        inputId = "Eperc",
        label = "Censoring (%)",
        value = 15,
        min = 0,
        max = 100
      ),
      radioButtons("modelling_type_experimental", label = ("Modelling type"), 
                  choices = list("Toxicity", "Disappointment"), 
                  selected = "Toxicity"),
      #h4("Control arm"),
      h4(div("Control arm", style = "color: red")),
      sliderInput(
        inputId = "time_frame_control",
        label = "Time Frame",
        min = 0,
        max = ceiling(max(original$xyz$x)),
        value = c(0,3),
        width = "100%"
      ),
      numericInput(
        inputId = "Cperc",
        label = "Censoring (%)",
        value = 15,
        min = 0,
        max = 100
      ),
      radioButtons("modelling_type_control", label = ("Modelling type"), 
                  choices = list("Toxicity", "Disappointment"), 
                  selected = "Disappointment")
    ),
    
    mainPanel(
      h3("Reconstructed Kaplan-Meier Curves"),
      plotOutput("plot"),
      # verbatimTextOutput("statistics"),
      h3("Statistical outputs (Cox)"),
      tags$style(HTML("
        #metrics {
          display: flex;
          justify-content: center;
        }
      ")),
      tableOutput("metrics")
      #box(
      #  title="Results",  br(),
      #  "exp(z):", textOutput("statistics_z"),
      #  "p-level:", textOutput("statistics_plevel"),
      #  "interval:", textOutput("statistics_interval")
      #)
    ),
    
  ),
)

simulateTest <- function(){
  newXYZ <<- newData(c(0,3),c(0,3), 50,50)$xyz
  newFit <<- survfit(Surv(x,y)~z, data=newXYZ)
}


server <- function(input, output) {
  
  userTrialData <- reactive({
    # do nothing if no file was uploaded
    if(! shiny::isTruthy(input$user_trial)){
      getOriginal()
    }else{
      # ext <- tools::file_ext(input$user_trial$name)
      # switch(ext,
      #        csv = vroom::vroom(input$user_trial$datapath, delim = ","),
      #        tsv = vroom::vroom(input$user_trial$datapath, delim = "\t"),
      #        validate("Invalid file; Please upload a .csv or .tsv file")
      # )
      readCurvesFromExcel(input$user_trial$datapath)
    }
  })
  
  userTrialFit <- reactive({
    survfit(Surv(x,y)~z, data=userTrialData()$xyz)
  })
  
  newXYZ <- reactive({
    newData(
      input$time_frame_experimental,input$time_frame_control,
      input$Eperc, input$Cperc,
      modellingTypeExperimental=input$modelling_type_experimental, modellingTypeControl=input$modelling_type_control,
      original=userTrialData()
    )$xyz
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
    bind_rows(lapply(
      list( Original=userTrialData()$xyz, Simulation=newXYZ()),
      function(xyz){
        fit <- coxph(Surv(x,y)~z, data=xyz)
        s <- summary(fit)
        data.frame(z = s$coefficients[2], "p-value" = s$coefficients[5],
                   "CI Low" = s$conf.int[3], "CI High" = s$conf.int[4]
        )
      }
    ), .id = "Model")
  })

  output$plot <- renderPlot({
    
    # fit <- survfit(Surv(x,y)~z, data = original$xyz)
    
    # fit <- survfit(Surv(x,y)~z+group, data = df)
    # fit <- survfit(Surv(x,y)~z, data = df)
    
    g <- plotCurves(userTrialData(), userTrialFit(), newXYZ(), newFit())
    
    print(g)
    
  })
  
  #output$statistics_z <- renderText({newStatistics()$z})
  #output$statistics_plevel <- renderText({newStatistics()$plevel})
  #output$statistics_interval <- renderText({newStatistics()$interval})
  
  output$metrics <- renderTable({
    newStatistics()
  })
  
}

shinyApp(ui, server)
