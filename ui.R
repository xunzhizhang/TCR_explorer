library(ggplot2)
library(Rcpp)
library(markdown)
library(Rtsne)
# library(cairo) for better fig in Linux
library(shiny)
options(stringsAsFactors = FALSE)
fluidPage(
  tags$h1(tags$strong("TCR Explorer"), align = "center",
          style = "font-family: Georgia; color: #000080;"),
  tags$h5(textOutput("time"), align = "right"),
  tags$hr(),
  sidebarLayout(
    sidebarPanel(
      selectInput("select_data", "Select demo data, user's data or raw TCR sequences",
                  c("Demo data" = "d",
                    "Visualize your data" = "u",
                    "Analyze your data" = "t")),
      tags$br(),
      selectInput("select_method","Select method to generate the coordinates",
                  c("PCA" = "p",
                    "t-SNE (PCA to obtain Top 50 elements)" = "s")),
      tags$br(),
      fileInput("tcr_file","Upload the .csv file (format as: tcr,x,y,color)",
                multiple=FALSE,
                accept=c(".csv")
      ),
      downloadButton('download_example', 'Example sequences'),
      downloadButton('download_result', 'Analyzed results'),
      tags$hr(),
      div(id = "tips", includeMarkdown("Tips.md")),
      width = 4
    ),
    
    mainPanel(
      # Interactive plot
      tags$h4("Brush on the left plot to zoom the target region"),
      fluidRow(
        column(
          width = 6,
          plotOutput(
            "plot_tcr",
            width = "100%", 
            click = "plot_click",
            brush = brushOpts(id = "plot_brush",
                              resetOnNew = TRUE)
          ),
          verbatimTextOutput("info")
        ),
        column(
          width = 6,
          plotOutput(
            "plot_tcr_z",
            width = "100%", 
            click = "plot_z_click"
          ),
          verbatimTextOutput("info_z")
        )
        
      )
    )
  )
  
)