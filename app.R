# Load packages ----
r <- getOption("repos") 
r <- append(BiocManager::repositories(), r)
options(repos = r)
library(shiny)
library(shinyjs)
library(rmarkdown)
library(HDMT)
# stored input data as an example ----
data("snp_input")

snp_input=snp_input[1:ceiling(nrow(snp_input)/10),]

# Source helper functions -----
source("helpers.R")

# User interface ----
ui <- fluidPage(
  titlePanel(h2("HDMT: a multiple-testing procedure for high-dimensional mediation analysis",align="center")),
  shinyjs::useShinyjs(),
  hr(),
  sidebarLayout(
    sidebarPanel(
      #tags$hr(),
      fileInput("file1", "Choose a CSV input File",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",".csv")),
      selectInput("var", 
                  label = "Or an example dataset for demo",
                  choices = c(" ", "meQTLdata"),
                  selected = " ",
                  multiple=FALSE),
      h6(helpText("An input dataset is a matrix containing two columns of p-values for candidate mediators.")),
      
      
      tags$hr(),
      radioButtons("radio1", "Choose an error measure",
                   choices = list("False discovery rate", "Family-wise error rate"),selected="False discovery rate"),
      
      #tags$hr(),
      sliderInput("range", 
                  label = "Choose a significance threshold",
                  min = 0, max = 0.2, step=0.05, value = c(0.05)),
      
      #tags$hr(),
      radioButtons("radio2", "Choose a correction method",
                   choices = list("Exact" = 1, "Approximation" = 0),selected = 1),
      br(),
      tags$hr(),
      br(),
      fluidRow(
        align = "center",
        actionButton("go", "Calculate"),
        shinyjs::hidden(p(id = "text1", "Processing..."))),
      br(),
      tags$head(tags$script(HTML('Shiny.addCustomMessageHandler("jsCode",function(message) {eval(message.value);});'))),
    ),
    
    mainPanel(
      
      # Output: Tabset w/ plot, summary, and table ----
      tabsetPanel(type = "tabs",
                  tabPanel("Results", 
                  
                  fluidRow(
                    align = "center",
                    plotOutput("plot")
                  ),
                  br(),
                  br(),
                  fluidRow(
                    align = "center",
                    h2(verbatimTextOutput("countsig"),align="center")
                  ),
                  br(),
                  #fluidRow(
                  #  align = "center",
                  #  tableOutput('table')
                  #),
                  br(),
                  #tags$hr(),
                  br(),
                  #uiOutput("downloadData1"),
                  fluidRow(
                    align = "center",
                    uiOutput("downloadData1"),
                  ),
                  
                  
                  ),
                  #rmarkdown::render("about.md") creat a html
                  tabPanel("More information",
                           fluidRow(includeHTML(rmarkdown::render("about.md")))
                           )
                  )
      
    )
  )
)

# Server logic ----
server <- function(input, output, session) {
  
  values <- reactiveValues()
  #fdr or fwer
  values$errmode <- NULL
  #flow control, refresh results if null 
  values$draw <- NULL
  #intermidiate results
  values$nullprop <- NULL
  values$result <- NULL
  #Are the results ready, indicator used to enable actionbutton
  values$progress <- NULL
  #input data
  values$data <- NULL
  shinyjs::hide("downloadData")
  # data1 <- reactive({req(input$var)
  #   "meQTLdata" = snp_input})
  data1 <- reactive({switch(input$var,
                            " " = NULL,
                            "meQTLdata" = snp_input)})
  
  data2 <- reactive({req(input$file1)
    df <- as.matrix(read.csv(input$file1$datapath,header = F)) })
  
  observeEvent(input$file1, {
    values$draw <- NULL
    values$data <- data2()
    shinyjs::enable("go")
    shinyjs::hide("text1")
  })
  
  observeEvent(input$var, {
    values$draw <- NULL
    values$data <- data1()
    shinyjs::enable("go")
    shinyjs::hide("text1")
  })
  
  observeEvent(input$range, {
    values$draw <- NULL
    #shinyjs::hide("downloadData")
    if (input$range == 0)
    {
      shinyjs::disable("go")
    }else
    {
      shinyjs::enable("go")
    }
  })
  
  observeEvent(input$radio1, {
    values$draw <- NULL
    values$errmode <- switch(input$radio1,
                              "False discovery rate" ="FDR",
                              "Family-wise error rate" = "FWER")
  })
  
  observeEvent(input$radio2, {
    values$draw <- NULL
  })
  
  observeEvent(input$go, {
    values$draw <- 'plot'
    #This is needed
    values$progress <- "Working"
    shinyjs::disable("go")
    shinyjs::show("text1")
    shinyjs::hide("downloadData")
    if (is.null(values$data))
    {
      js_string <- 'alert("Please provide an input dataset");'
      session$sendCustomMessage(type='jsCode', list(value = js_string))
      shinyjs::enable("go")
      shinyjs::hide("text1")
      values$draw <- NULL
      
    }
  })
  
  output$plot <- renderPlot({
    if (!is.null(values$draw))
    {
      
      nullprop <- get_nullprop(values$data)
      values$nullprop <- nullprop
      Corrected_qqplot(values$data,nullprop,input$radio2)

    }
  }, height = 400, width = 400)
  
  output$countsig <- renderText({
    if (!is.null(values$draw))
    {
      outtext <- NULL
      if (values$errmode =="FDR")
      {
        result <- get_fdr_est(values$data,values$nullprop,input$radio2,input$range[1])
        outtext <- paste0(result$nsig," mediation tests are signficant at FDR < ", input$range[1],"!")
      }

      if (values$errmode =="FWER")
      {
        result <- get_fwer_sig(values$data,values$nullprop,input$radio2,input$range[1])
        outtext <- paste0(result$nsig," mediation tests are signficant at FWER < ", input$range[1],"!")
      }
      
      values$result <- result
      values$progress <- "Ready"

      if (values$progress == "Ready")
      {
        shinyjs::enable("go")
        shinyjs::hide("text1")
        shinyjs::show("downloadData")     
      }
      
      outtext
    }
  })
  
  #output$table <- renderTable({
  #  if (!is.null(values$draw))
  #  {
  #    #withProgress(message = 'Making plot...', value = 1, {Sys.sleep(0.2)})
  #    result1 <- NULL
  #    if (!is.null(values$result))
  #    {
  #      if (values$result$nsig>0)
  #      {
  #        result1 <- values$result$output_p[1:min(values$result$nsig,10),]
  #      }
  #    }
  #    if (values$result$nsig==0)
  #    {
  #      shinyjs::disable("downloadData")
  #    }
  #    
  #    result1
  #  }
  #}, digits=-2, height = 200, width = 400)

  output$downloadData1 <- renderUI({
    
    req(values$draw)
    if (values$result$nsig>0) downloadButton('downloadData', label = 'Download significant tests') })
  
  output$downloadData <- downloadHandler(
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = "Significant_tests.txt",
    
    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(filename) {
      sep <-"\t"
      # Write to a file specified by the 'file' argument
      if (values$result$nsig>0)
      write.table(values$result$output_p, filename, sep = sep,
                  row.names = FALSE)

     },
    contentType = ".txt"
  )
  
}

# Run app ----
shinyApp(ui, server)