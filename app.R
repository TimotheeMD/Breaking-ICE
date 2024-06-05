library(shiny)
library(gridlayout)
library(bslib)
library(shinydashboard)
library(colourpicker)

source("breaking_ice.R")

TEMPLATE_DOWNLOAD_URL='https://raw.githubusercontent.com/TimotheeMD/SWOG1801_reanalysis/main/CON.xlsx'
TEMPLATE_DOWNLOAD_URL='template.xlsx'

# Define UI

tags$head(
  tags$style(HTML("
    body, html {
      height: 100%;
      margin: 0;
      padding: 0;
    }
    .container-fluid {
      height: 100%;
    }
    .row {
      height: 100%;
    }
  "))
)


ui <- fluidPage(
  fluidRow( 
    column (8,
  box(
    title="Trial Data",  br(),
    fluidRow(
      column(6,
        selectInput(
          inputId = "prepared_trial",
          label = "Select a trial",
          choices = names(allDatasets),
          selected = "CONTACT-02-PFS",
          width = "100%"
        ),
      ),
      column(6,
        fileInput("user_trial", tags$span("Or Load your own Trial Data", tags$a(href=TEMPLATE_DOWNLOAD_URL,"(Template)"), ":"), accept = c(".xlsx")),
      )
    ),
    width=12
  )),
  column (1, box()), # empty puffer
  column (3, 
  box(
    title="Breaking-ICE",br(),
    "Informative",br(),
    "Censoring",br(),
    "Exploration",br(),
    width=12
  ))
  ),
  
  sidebarLayout(
    
    sidebarPanel(
      h3("Sensitivity analysis"),   
      "Change the parameters to see how it affects the survival analysis",br(),
      "1- 'Toxicity' modelling modify censored patients into events",br(),
      "2- 'Disappointment' modelling modify the time of censoring into late censoring",br(),
      # h4("Experimental arm"),
      h4(div("Experimental arm", style = "color: #6699CC")),
      sliderInput(
        inputId = "time_frame_experimental",
        label = "Time Frame (when censoring is modified)",
        min = 0,
        max = ceiling(max(original$xyz$x)),
        value = c(0,3),
        width = "100%"
      ),
      numericInput(
        inputId = "Eperc",
        label = "Percentage of censored patients to be modified",
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
        label = "Time Frame (when censoring is modified)",
        min = 0,
        max = ceiling(max(original$xyz$x)),
        value = c(0,3),
        width = "100%"
      ),
      numericInput(
        inputId = "Cperc",
        label = "Percentage of censored patients to be modified",
        value = 15,
        min = 0,
        max = 100
      ),
      radioButtons("modelling_type_control", label = ("Modelling type"), 
                  choices = list("Toxicity", "Disappointment"), 
                  selected = "Disappointment"),
      actionButton("doReset", "Set All Parameters to Zero"),

      h3("Colours and Transparency"),
      # palette =
      #   c("#6699CC", "red", "green","orange"),    # custom color palettes
      # legend.labs =
      #   c("Original - Experiment", "Original - Control", 'Simulated - Experiment', 'Simulated - Control'),
      fluidRow(
        column(6,
          colourInput("col_original_exp", "Original - Experiment", "#6699CC", allowTransparent=TRUE, closeOnClick=TRUE),
          colourInput("col_simulated_exp", "Simulated - Experiment", "green", allowTransparent=TRUE, closeOnClick=TRUE),
        ),
        column(6,
          colourInput("col_original_con", "Original - Control", "red", allowTransparent=TRUE, closeOnClick=TRUE),
          colourInput("col_simulated_con", "Simulated - Control", "orange", allowTransparent=TRUE, closeOnClick=TRUE),
      ))
    ),
    
    mainPanel(
      h3("Reconstructed Kaplan-Meier Curves"),
      plotOutput("plot"),
      # verbatimTextOutput("statistics"),
      h3("Statistical outputs (Cox)"),
      tags$style(HTML("
        #metrics,#censor_perc_table {
          display: flex;
          justify-content: center;
        }
      ")),
      tableOutput("metrics"),
      h3("Censored patients between each time interval (%)"),
      tableOutput("censor_perc_table")
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


server <- function(input, output, session) {
  
  userTrialData <- reactive({
    # do nothing if no file was uploaded
    if(! shiny::isTruthy(input$user_trial)){
      trial <- input$prepared_trial
      if(trial %in% names(allDatasets)){
        allDatasets[[trial]] %>% glimpse
      }else{
        getOriginal()
      }
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
      list( "From Original"=userTrialData()$xyz, "Sensitivity Analysis"=newXYZ()),
      function(xyz){
        fit <- coxph(Surv(x,y)~z, data=xyz)
        s <- summary(fit)
        data.frame("HR" = 1/s$coefficients[2], "p-value" = s$coefficients[5],
                   "CI Low" = 1/s$conf.int[4], "CI High" = 1/s$conf.int[3]
        )
      }
    ), .id = "Model")
  })

  output$plot <- renderPlot({
    
    # fit <- survfit(Surv(x,y)~z, data = original$xyz)
    
    # fit <- survfit(Surv(x,y)~z+group, data = df)
    # fit <- survfit(Surv(x,y)~z, data = df)
    
    g <- plotCurves(userTrialData(), userTrialFit(), newXYZ(), newFit(),
            colours=c(input$col_original_exp, input$col_original_con, input$col_simulated_exp, input$col_simulated_con)
            )
    
    print(g)
    
  })
  
  #output$statistics_z <- renderText({newStatistics()$z})
  #output$statistics_plevel <- renderText({newStatistics()$plevel})
  #output$statistics_interval <- renderText({newStatistics()$interval})
  
  output$metrics <- renderTable({
    newStatistics()
  })
  
  output$censor_perc_table <- renderTable({
    calculateCensorPerc(userTrialData(), newXYZ())
  })
  
  observeEvent(input$doReset, {
    print("resetting")
    updateSliderInput(session, "time_frame_experimental", value=c(0,0))
    updateNumericInput(session, "Eperc", value=0)
    updateRadioButtons(session, "modelling_type_experimental", selected="Toxicity")
    updateSliderInput(session, "time_frame_control", value=c(0,0))
    updateNumericInput(session, "Cperc", value=0)
    updateRadioButtons(session, "modelling_type_control", selected="Disappointment")
  })
  
}

shinyApp(ui, server)
