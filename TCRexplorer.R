# This file is for test of Rcapp incooperation.
library(ggplot2)
library(Rcpp)
# library(cairo) for better fig in Linux
library(shiny)

ui <- fluidPage(
  tags$h1("TCR Sequences",
          style = "font-family: 'Imgact'; color: #00008B;"),
  tags$hr(),
  selectInput("select_data", "Select data demo, user's data or raw TCR sequences",
              c("Data demo" = "d",
                "User's data" = "u",
                "TCR sequences" = "t")),
  fileInput("tcr_file","Upload the .csv file (format as: tcr,x,y,color)",
            multiple=FALSE,
            accept=c(".csv")
  ),
  downloadButton('download_file', 'Download analyzed results'),
  
  # Interactive plot
  tags$h4("Brush to zoom the target region"),
  fluidRow(
    column(
      width = 6,
      plotOutput(
        "plot_tcr",
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
        click = "plot_z_click"
      ),
      verbatimTextOutput("info_z")
    )
    
  )
)

server <- function(input, output) {
  #Rcpp for test
  cppFunction(
    'CharacterVector colv(size_t len){
            CharacterVector color(len);
            for (size_t i=0; i < len; i++){
            color[i]="COLOR";
            }
            return color;
    }'
  )  
  
  sourceCpp(code = '
            #include <Rcpp.h>
            #include <iostream>
            #include <fstream>
            #include <sstream>
            #include <vector>
            #include <string>            
            using namespace Rcpp;

            int  BLOSUM62[128][128];     
            // [[Rcpp::export]]
            int loadBLOSUM62(){
                for (int i = 0; i < 128; ++i) {
                    for (int j = 0; j < 128; ++j) {
                       BLOSUM62[i][j] = 0;
                    }
                }
                std::ifstream file("blosum62.3col.txt");
                std::string   line;
            
                while(std::getline(file, line)){
                     std::stringstream   linestream(line);
                     std::string         a1;
                     std::string         a2;
                     int                 val;
                     std::getline(linestream, a1, \'\\t\');
                     std::getline(linestream, a2, \'\\t\');
                     linestream >> val;
                     BLOSUM62[(unsigned char) a1[0]][(unsigned char) a2[0]] = val;  
                }
              }
            
            // [[Rcpp::export]]
            int di(const std::string& v1, const std::string& v2){
              size_t len = v1.size();
              if (len > v2.size()) {
                  len = v2.size();
              }
              int d = 0;
              for (size_t i = 0 ; i != len; ++i) {
                  d += BLOSUM62[(unsigned char)v1[i]][(unsigned char)v2[i]];
              }
              return d;
            }
            // [[Rcpp::export]]
            NumericMatrix loc_m(CharacterVector s, size_t len){
                loadBLOSUM62();
                std::vector<std::string> seq = Rcpp::as <std::vector<std::string> > (s);
                NumericMatrix dmatrix(len, len);
                for (size_t i = 0; i < len; i++){
                    for (size_t j = 0; j <= i; j++){
                        dmatrix(i, j) = di(seq[i], seq[j]);
                    }
                }
                for (size_t i = 0; i < len; i++){
                    for (size_t j = 0; j < i; j++){
                        dmatrix(i, j) = 1.0 - 2.0 * di(seq[i], seq[j]) / (dmatrix(i,i) + dmatrix(j, j));
                        dmatrix(j, i) = dmatrix(i, j);
                    }
                }
                for (size_t i = 0; i < len; i++){
                    dmatrix(i, i) = 0.0;
                }
                
                return dmatrix;
            } 

  ')
  # Initilize ranges for axis
  ranges <- reactiveValues(rx = NULL, ry = NULL)
  ranges_z <- reactiveValues(rzx = NULL, rzy = NULL)
  # Read demo data/user data/TCR seq
  # Read the file as a data frame, 1st line as header, return empty data frame without input
  # Note: a .csv file with only seq column will be treated to generate entire file
  safe.read <- function(docu, slc) {
    if (slc == "d"){
      return(
        read.csv(
          "tcr.data.demo.csv",  # file under the same folder
          header = TRUE,
          sep = ",",
          stringsAsFactors = FALSE
        )
      )
    } else if (is.null(docu)) {
      empdata<-data.frame(tcr=c(NA),x=c(NA),y=c(NA),color=c(NA))
      return(empdata)
    } else if (slc == "u") {
      return(read.csv(
        docu$datapath,
        header = TRUE,
        sep = ",",
        stringsAsFactors = FALSE
      ))
    } else {
      tcr_seq <- read.csv(
        docu$datapath,
        header = TRUE,
        sep = ",",
        stringsAsFactors = FALSE
      )[, 1] # read a data frame with mere seq, get a char vector WITHOUT title
      # print(tcr_seq)
      # print(tcr_seq[1])
      # print(BLOSUM62)
      # analyze with C++, get a matrix and then transfer to dataframe
      slen <- length(tcr_seq)   # cal len of data in R
      # results <- r_result(tcr_seq, slen)
      # get dis matrix and then transfer to dataframe with PEA
      # sourceCpp("maintest.cpp")
      #dist_matrix <- as.matrix(dist(loc_m(tcr_seq, slen)))
      dist_matrix <- loc_m(tcr_seq, slen)
      save(dist_matrix, file = "saved.dist_matrix.Rdata")
      # print(dist_matrix[1:5, 1:5])
      # print(tcr_seq)
      # print(slen)
      # save.image(file = "saved.Rdata")
      h <- diag(slen) - matrix(1, slen, 1) %*% matrix(1, 1, slen) / slen
      ret <- eigen(h %*% (-0.5 * dist_matrix) %*% h, symmetric = TRUE)
      x <- ret$vectors [, 1:2] %*% diag(ret$values[1:2])[,1]
      y <- ret$vectors [, 1:2] %*% diag(ret$values[1:2])[,2]
      color <- colv(slen)
      results <- as.data.frame(cbind(tcr_seq,x,y,color))
      # print(results$x)
      # results$x <- as.numeric(as.character(results$x))
      # results$y <- as.numeric(as.character(results$y))
      # print(results$x)
      names(results) <- c("tcr","x","y","color")
      results$x <- as.numeric(as.character(results$x))
      results$y <- as.numeric(as.character(results$y))
      return(results)
    }
  }
  # set ranges for entire data
  # observeEvent(input$tcr_file,{
  #   point_result <- safe.read(input$tcr_file, input$select_data)
  #   if (!is.null(point_result$tcr)){
  #     ranges$rx <- c(min(point_result$x),max(point_result$x))
  #     ranges$ry <- c(min(point_result$y),max(point_result$y))
  #     print(ranges$rx)
  #   } else {
  #     ranges$rx <- NULL
  #     ranges$ry <- NULL
  # 
  #   }
  # })
  
  
  # Zoom after brush
  observeEvent(input$plot_brush,{
    brush <- input$plot_brush
    if (!is.null(brush)) {
      ranges_z$rzx <- c(brush$xmin, brush$xmax)
      ranges_z$rzy <- c(brush$ymin, brush$ymax)
      # print(ranges_z$rzx)
      # print(ranges_z$rzy)
    } else {
      ranges_z$rzx <- NULL
      ranges_z$rzy <- NULL
    }
  })
  
  # Construct plot with zooming function
  # Left plot with entire data
  output$plot_tcr <- renderPlot({
    ggplot(
      safe.read(input$tcr_file, input$select_data),
      aes(x = x, y = y)) + geom_point() + 
      coord_cartesian(xlim = ranges$rx, ylim = ranges$ry, expand = TRUE) +
      ggtitle("Entire data")
  })
  
  output$info <- renderPrint({
    nearPoints(safe.read(input$tcr_file, input$select_data),
               input$plot_click,
               xvar = "x",
               yvar = "y",
               threshold = 50,
               addDist = TRUE
    )
  })
  # Right plot for zooming
  output$plot_tcr_z <- renderPlot({
    ggplot(
      safe.read(input$tcr_file, input$select_data),
      aes(x = x, y = y)) + geom_point() + coord_cartesian(xlim = ranges_z$rzx, ylim = ranges_z$rzy, expand = TRUE) +
      ggtitle("Chosen region")
  })
  output$info_z <- renderPrint({
    nearPoints(safe.read(input$tcr_file, input$select_data),
               input$plot_z_click,
               xvar = "x",
               yvar = "y",
               threshold = 50,
               addDist = TRUE
    )
  })
  
  # Generate downloadable file with data frame generated by C++ or from input
  # Row names should be included in the data frame
  output$download_file <- downloadHandler(
    filename = function() {
      paste("TCR_analysis-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(safe.read(input$tcr_file, input$select_data), file, row.names=FALSE)
    }
  )
}


shinyApp(ui = ui, server = server)