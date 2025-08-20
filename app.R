library(shiny)
library(gridlayout)
library(bslib)
library(shinydashboard)
library(colourpicker)
library(DT)

source("breaking_ice.R")

TEMPLATE_DOWNLOAD_URL <- 'template.xlsx'
CUMULATIVE_BG_COLOR <- '#F0F0F0'

## Change here the defaults for the simulation parameters if they are not specified in the Excel 'Sensitivity' Sheet
DEFAULTS <- list(
  EXP=list(
    time_frame_max=0,
    modelling_type='Toxicity',
    percentage=0
  ),
  CON=list(
    time_frame_max=0,
    modelling_type='Disappointment',
    percentage=0
  )
)

updateSensitivityParameters <- function(session, trial_data, defaults = DEFAULTS) {
  # Get the maximum time from the trial data
  max_time <- max(trial_data$trisk, na.rm = TRUE)
  
  # Function to get parameter value, either from trial data or defaults
  getParameterValue <- function(param_name, arm) {
    if ('simulationDefaults' %in% names(trial_data)) {
      print('has simulationDefaults')
      str(trial_data$simulationDefaults)
      value <- (trial_data$simulationDefaults) %>% 
        filter(Parameter == param_name) %>% 
        pull(arm)
      if(length(value) > 0) return(value)
    }
    return(defaults[[arm]][[param_name]])
  }
  
  # Update all sensitivity parameters
  updateSliderInput(session, "time_frame_experimental", 
                    value = c(0, getParameterValue('time_frame_max', 'EXP')), 
                    max = max_time)
  updateNumericInput(session, "Eperc", 
                     value = getParameterValue('percentage', 'EXP'))
  updateRadioButtons(session, "modelling_type_experimental", 
                     selected = getParameterValue('modelling_type', 'EXP'))
  
  updateSliderInput(session, "time_frame_control", 
                    value = c(0, getParameterValue('time_frame_max', 'CON')), 
                    max = max_time)
  updateNumericInput(session, "Cperc", 
                     value = getParameterValue('percentage', 'CON'))
  updateRadioButtons(session, "modelling_type_control", 
                     selected = getParameterValue('modelling_type', 'CON'))
}

ui <- fluidPage(
  shinyjs::useShinyjs(),
  tags$head(
    tags$style(HTML(paste0("
    body, html {
      height: 100%;
      margin: 0;
      padding: 0;
    }
    .colourpicker-panel {
      bottom: 2.5em;
    }
    .td-cumulative {
      background: ",CUMULATIVE_BG_COLOR,";
    }
  ")))
  ),
  
  # Title section
  fluidRow(
    column(12, 
           tags$div(
             style = "text-align: center; margin-top: 40px; margin-bottom: 40px;",
             tags$h1(style = "font-weight: bold;", "BREAKING-ICE App©"),
             tags$h4(style = "font-style: italic;", "Informative Censoring Exploration")
           )
    )
  ),
  
  # Application content
  fluidRow( 
    column(10,
           box(
             fluidRow(
               column(5,
                      selectInput(
                        inputId = "prepared_trial",
                        label = tags$span("Select a Trial Data",
                                          ":"),
                        choices = names(allDatasets),
                        selected = "CONTACT-02-PFS",
                        width = "100%"
                      ),
                      uiOutput('trial_reference')
               ),
               column(4,
                      fileInput("user_trial", 
                                tags$span("Or Load your own Trial Data", 
                                          tags$a(href=TEMPLATE_DOWNLOAD_URL,"(Template)"), ":"),
                                accept = c(".xlsx")
                      )
               ),
               column(3,shinyjs::hidden(input_switch('do_user_file','Use uploaded file')))
             ),
             width = 12
           ))
  ),
  
  sidebarLayout(
    sidebarPanel(
      h3("Sensitivity analysis"),   
      "Change the parameters to see how it affects the survival analysis",br(),
      "1- 'Toxicity' modelling modify censored patients into events",br(),
      "2- 'Disappointment' modelling modify the time of censoring into late censoring",br(),
      h4(div("Experimental arm", style = "color: #6699CC")),
      sliderInput(
        inputId = "time_frame_experimental",
        label = "Time Frame (when censoring is modified)",
        min = 0,
        max = ceiling(max(original$xyz$x)),
        value = c(0, DEFAULTS$EXP$time_frame_max ),
        width = "100%"
      ),
      numericInput(
        inputId = "Eperc",
        label = "Percentage of censored patients to be modified",
        value = DEFAULTS$EXP$percentage,
        min = 0,
        max = 100
      ),
      radioButtons("modelling_type_experimental", label = ("Modelling type"), 
                   choices = list("Toxicity", "Disappointment"), 
                   selected = DEFAULTS$EXP$modelling_type),
      h4(div("Control arm", style = "color: red")),
      sliderInput(
        inputId = "time_frame_control",
        label = "Time Frame (when censoring is modified)",
        min = 0,
        max = ceiling(max(original$xyz$x)),
        value = c(0,DEFAULTS$CON$time_frame_max),
        width = "100%"
      ),
      numericInput(
        inputId = "Cperc",
        label = "Percentage of censored patients to be modified",
        value = DEFAULTS$CON$percentage,
        min = 0,
        max = 100
      ),
      radioButtons("modelling_type_control", label = ("Modelling type"), 
                   choices = list("Toxicity", "Disappointment"), 
                   selected = DEFAULTS$CON$modelling_type),
      actionButton("doReset", "Set All Parameters to Zero"),
      
      h3("Colours and Transparency"),
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
      
      h3("Statistical outputs (Cox)"),
      tags$style(HTML("
    #metrics, #censor_perc_table, #qual_table {
      display: flex;
      justify-content: flex-start;
    }
    .table-section {
      display: flex;
      justify-content: space-between;
      margin-bottom: 20px;
      align-items: center;
    }
    table:not(.dataTable) th, table:not(.dataTable) td {
      text-align: left !important;
      padding: 8px;
    }
  ")),
      
      tableOutput("metrics"),
      
      h3("Censored patients between each time interval (%)"),
      DTOutput('censor_perc_table_dt'),

      h3("Quality Control (percentage differences between original and digitilized curves)"),
      div(class = "table-section",
          div(style = "flex: 1;",
              tableOutput("qual_table")
          ),
          div(style = "margin-left: 20px;",
              tags$button(
                id = "switchToReverse", 
                onclick = "localStorage.setItem('selectedTrial', $('#prepared_trial').val()); window.open('https://www.timotheeolivier-research.com/reverse-km', '_blank');",
                "Go to Reverse Kaplan-Meier Analysis",
                class = "btn btn-primary",
                style = "font-size: 16px; padding: 10px 20px; background-color: #E9ECEF; color: #333; border: none; border-radius: 4px;"
              )
          )
      )
    )
  )
)

server <- function(input, output, session) {
  userTrialData <- reactive({
    trial_data <- if(!input$do_user_file || !shiny::isTruthy(input$user_trial)) {
      # Case 1: Using prepared trial data
      trial <- input$prepared_trial
      if(trial %in% names(allDatasets)) {
        print(paste("Loading dataset:", trial))
        ret <- getDataset(trial)
        ret
      } else {
        print("Loading original dataset")
        getOriginal()
      }
    } else {
      # Case 2: Using uploaded file
      print("Loading uploaded file")
      readCurvesFromExcel(input$user_trial$datapath)
    }
    
    # Update sensitivity parameters regardless of data source
    updateSensitivityParameters(session, trial_data)
    
    return(trial_data)
  })
  
  userTrialFit <- reactive({
    survfit(Surv(x,y)~z, data=userTrialData()$xyz)
  })
  
  newXYZ <- reactive({
    newData(
      input$time_frame_experimental,input$time_frame_control,
      input$Eperc, input$Cperc,
      modellingTypeExperimental=input$modelling_type_experimental, 
      modellingTypeControl=input$modelling_type_control,
      original=userTrialData()
    )$xyz
  })
  
  isAllParamsZero <- reactive({
    (input$Eperc == 0 || diff(input$time_frame_experimental)==0) &&
      (input$Cperc == 0 || diff(input$time_frame_control)==0)
  })
  
  newFit <- reactive({
    survfit(Surv(x,y)~z, data=newXYZ())
  })
  ##### change srpintf to 4f for 4 decimals, 5f for 5 decimals, etc.
  newStatistics <- reactive({
    res <- bind_rows(lapply(
      list("From Original"=userTrialData()$xyz, "Sensitivity Analysis"=newXYZ()),
      function(xyz){
        fit <- coxph(Surv(x,y)~z, data=xyz)
        s <- summary(fit)
        data.frame("HR" = 1/s$coefficients[2], "p-value" = sprintf("%.3f", s$coefficients[5]),
                   "CI Low" = 1/s$conf.int[4], "CI High" = 1/s$conf.int[3]
        )
      }
    ), .id = "Model")
    if(isAllParamsZero()){
      res[2,-1] <- NA
    }
    res
  })
  
  output$plot <- renderPlot({
    g <- if(isAllParamsZero()){
      plotCurvesWithoutSimulation(userTrialData(), userTrialFit(),
                                  colours=c(input$col_original_exp, input$col_original_con, 
                                            input$col_simulated_exp, input$col_simulated_con)
      )
    }else{
      plotCurves(userTrialData(), userTrialFit(), newXYZ(), newFit(),
                    colours=c(input$col_original_exp, input$col_original_con, 
                              input$col_simulated_exp, input$col_simulated_con)
      )
    }
    print(g)
  })
  
  output$number_at_risk_table <- renderTable({
    calculateNumberAtRisk(userTrialData(), newXYZ())
  }, rownames = TRUE, align='r', width='100%')
  
  output$metrics <- renderTable({
    newStatistics()
  })
  
  censor_perc_table <- reactive({
    data <- calculateCensorPerc(userTrialData(), newXYZ())
    data[] <- lapply(names(data), function(col_name) {
      if (col_name == "Time") {
        paste(data[[col_name]], ' - ', lead(data[[col_name]]))
      } else if (is.numeric(data[[col_name]])) {
        sprintf("%.1f", data[[col_name]])
      } else {
        data[[col_name]]
      }
    })
    if(isAllParamsZero()){
      data[,4:5] <- NA
    }
    return(data)
  })
  
  # output$censor_perc_table <- renderTable({
  #   censor_perc_table()
  # }, sanitize.text.function = function(x) x, align = "r")
  
  output$censor_perc_table_dt <- renderDT({
    df <- censor_perc_table()
    df <- bind_cols(df, calculateCensoredCumulative(userTrialData(), newXYZ())[,-1])
    if(isAllParamsZero()){
      df[,8:9] <- NA
    }
    
    sketch = htmltools::withTags(table(
      class = 'display',
      thead(
        tr(
          th(rowspan = 3, 'Time'),
          th(colspan = 4, 'Interval'),
          th(colspan = 4, 'Cumulative', class="td-cumulative")
        ),
        tr(
          th(colspan = 2, 'Original'),
          th(colspan = 2, 'Sensitivity'),
          th(colspan = 2, 'Original', class="td-cumulative"),
          th(colspan = 2, 'Sensitivity', class="td-cumulative")
        ),
        tr(
          lapply(rep(c('Experiment','Control'), 2), th),
          lapply(rep(c('Experiment','Control'), 2), th, class="td-cumulative")
        )
      )
    ))
    datatable(df, container = sketch, rownames = FALSE,
              # filter='none',selection='none',
              # pageLength=1000,
              # autoHideNavigation=TRUE,
              options=list(
                searching=FALSE, paging=FALSE, info=FALSE,ordering=FALSE,
                pageLength=nrow(df)
                ),
              ) %>%
      formatStyle(6:9, backgroundColor = CUMULATIVE_BG_COLOR) %>% 
      formatStyle(2:9, textAlign="right") %>% 
      formatRound(2:9, 1)
  })
  
  output$qual_table <- renderTable({
    trial_data <- userTrialData()
    print(str(trial_data))
    print(names(trial_data))
    
    if (!"Quality" %in% names(trial_data)) {
      print("Feuille 'Quality' absente")
      return(data.frame(
        "Metric" = c("HR", "Length of CI"),
        "Difference" = c("Not available", "Not available")
      ))
    }
    
    quality_data <- trial_data$Quality
    print(quality_data)
    
    HRPub <- quality_data[1, "HR"]
    CI_LowPub <- quality_data[1, "CI Low"]
    CI_HighPub <- quality_data[1, "CI High"]
    
    if (is.na(HRPub) | is.na(CI_LowPub) | is.na(CI_HighPub)) {
      print("Données manquantes dans Quality")
      return(data.frame(
        "Metric" = c("HR", "Length of CI"),
        "Difference" = c("Not available", "Not available")
      ))
    }
    
    CI_lengthPub <- as.numeric(CI_HighPub) - as.numeric(CI_LowPub)
    
    original_data <- trial_data$xyz
    fit <- coxph(Surv(x, y) ~ z, data = original_data)
    s <- summary(fit)
    
    HR <- 1 / s$coefficients[2]
    CI_Low <- 1 / s$conf.int[4]
    CI_High <- 1 / s$conf.int[3]
    CI_length <- CI_High - CI_Low
    
    HR_Qual <- abs((HR - as.numeric(HRPub)) / as.numeric(HRPub)) * 100
    CI_Qual <- abs((CI_length - CI_lengthPub) / CI_lengthPub) * 100
    
    data.frame(
      "Metric" = c("HR", "Length of CI"),
      "Difference" = c(sprintf("%.2f%%", HR_Qual), sprintf("%.2f%%", CI_Qual))
    )
  }, sanitize.text.function = function(x) x, align = "r")
  
  output$trial_reference <- renderUI({
    trial_data <- userTrialData()
    if("weblink" %in% names(trial_data)) {
      print(list(weblink=trial_data$weblink))
      if(length(trial_data$weblink) >= 2) {
        tags$div(
          "(",
          tags$a("Original Publication", href=trial_data$weblink[1], target="_blank"),
          "-",
          tags$a("Sensitivity Publication", href=trial_data$weblink[2], target="_blank"),
          ")"
        )
      } else {
        tags$a("(Original Publication)", href=trial_data$weblink[1], target="_blank")
      }
    } else {
      NULL
    }
  })
  
  
  observeEvent(input$doReset, {
    updateSliderInput(session, "time_frame_experimental", value=c(0,0))
    updateNumericInput(session, "Eperc", value=0)
    updateRadioButtons(session, "modelling_type_experimental", selected="Toxicity")
    updateSliderInput(session, "time_frame_control", value=c(0,0))
    updateNumericInput(session, "Cperc", value=0)
    updateRadioButtons(session, "modelling_type_control", selected="Disappointment")
  })
  
  observeEvent(input$user_trial,{
    shinyjs::show('do_user_file')
    update_switch('do_user_file', value=TRUE)
  })
  
  observeEvent(input$prepared_trial,{
    update_switch('do_user_file', value=FALSE)
  })
  
  output$trial_reference <- renderUI({
    trial_data <- userTrialData()
    if("weblink" %in% names(trial_data)) {
      if(length(trial_data$weblink) >= 2) {
        tags$div(
          "(",
          tags$a("Original Publication", href=trial_data$weblink[1], target="_blank"),
          "-",
          tags$a("Sensitivity Publication", href=trial_data$weblink[2], target="_blank"),
          ")"
        )
      } else {
        tags$a("(Original Publication)", href=trial_data$weblink[1], target="_blank")
      }
    } else {
      NULL
    }
  })
  
  observeEvent(input$switchToReverse, {
    current_trial <- input$prepared_trial
    session$sendCustomMessage(
      type = "storeTrialAndRedirect",
      message = list(
        trial = current_trial,
        url = "https://www.timotheeolivier-research.com/reverse-km"
      )
    )
  })
}

shinyApp(ui=ui, server=server)