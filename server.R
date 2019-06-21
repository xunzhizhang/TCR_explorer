library(ggplot2)
library(Rcpp)
library(markdown)
library(Rtsne)
# library(cairo) for better fig in Linux
library(shiny)
options(stringsAsFactors = FALSE)
dirtyFlag <- TRUE ## marked if the computation is needed or not
cacheResult <- list()

function(input, output) {
  sourceCpp("TCRfunctions.cpp")
  # Initilize ranges for axis
  ranges <- reactiveValues(rx = NULL, ry = NULL)
  ranges_z <- reactiveValues(rzx = NULL, rzy = NULL)
  # function to get first 500 rows of raw seq
  len.test <- function(raw_seq) {
    l <- length(raw_seq)
    if (l >500){
      cat("Warning: >500 rows detected! Only first 500 rows will be analyzed.\n")
      return(head(raw_seq, 500))
    } else {
      return(raw_seq)
    }
  }
  # Read demo data/user data/TCR seq
  # Read the file as a data frame, 1st line as header, return empty data frame without input
  # Note: a .csv file with only seq column will be treated to generate entire file
  safe.read <- function(docu, slc_data, slc_method) {
    if (slc_data == "d"){
      return(
        read.csv(
          "tcr.data.demo.csv",  # file under the same folder
          header = TRUE,
          sep = ",",
          stringsAsFactors = FALSE
        )
      )
    } else if (is.null(docu)) {
      empdata<-data.frame(tcr=c(NA),x=c(NA),y=c(NA),color=c(NA)) # color is required here
      return(empdata)
    } else if (slc_data == "u") {
      return(read.csv(
        docu$datapath,
        header = TRUE,
        sep = ",",
        stringsAsFactors = FALSE
      ))
    } else {
      results <- data.frame(tcr=c(NA),x=c(NA),y=c(NA),color=c(NA))
      tcr_seq <- read.csv(
        docu$datapath,
        header = TRUE,
        sep = ",",
        stringsAsFactors = FALSE
      )[, 1] # read a data frame with mere seq, get a char vector WITHOUT title
      # analyze with C++, get a matrix and then transfer to dataframe
      # Only first 500 rows are treated
      tcr_seq <- len.test(tcr_seq)
      slen <- length(tcr_seq) 
      dist_matrix <- loc_m(tcr_seq, slen)
      # calculate the coordinates
      if (slc_method == "p"){
        h <- diag(slen) - matrix(1, slen, 1) %*% matrix(1, 1, slen) / slen
        ret <- eigen(h %*% (-0.5 * dist_matrix) %*% h, symmetric = TRUE)
        x <- ret$vectors [, 1:2] %*% diag(ret$values[1:2])[,1]
        y <- ret$vectors [, 1:2] %*% diag(ret$values[1:2])[,2]
        cat("PCA\n")
      } else {
        # at least >3*perplexity+1 rows are required
        set.seed(1)
        tsne <- Rtsne(
          dist_matrix,
          dims = 2,
          perplexity = 30,
          max_iter = 1000,
          theta = 0.0
        )
        x <- tsne$Y[,1]
        y <- tsne$Y[,2]
        cat("t-SNE\n")
      }
      color <- colv(tcr_seq, slen)
      results <- as.data.frame(cbind(tcr_seq,x,y,color), stringsAsFactors = FALSE)
      names(results) <- c("tcr","x","y","color")
      results$tcr <- as.character(results$tcr)
      results$x <- as.numeric(as.character(results$x))
      results$y <- as.numeric(as.character(results$y))
      results$color <- as.character(results$color)
      return(results)
    }
  }
  
  # sort color with freq, and add freq to data
  color.dist <- function(frame) {
    sort_freq <- sort(table(frame$color), decreasing = T, na.last = NA)
    # avoid transfering to int with one cluster
    sort_freq <- as.table(sort_freq)
    # NA will be generated in table without color, which cannot be directly converted
    if (is.na(sort_freq[1])){
      sort_frame <- data.frame(color=c(NA),frequency=c(NA), plot_color=c("#ffffff"))
    } else {
      sort_frame <- as.data.frame(sort_freq, stringsAsFactors = FALSE)
      sort_frame$color <- rainbow(nrow(sort_frame))
    }
    names(sort_frame) <- c("color", "frequency", "plot_color")
    # match and add freq after the original data
    r <- data.frame(tcr=frame$tcr, 
                    x=frame$x, 
                    y=frame$y, 
                    color=frame$color,
                    frequency=sort_frame$frequency[match(frame$color, sort_frame$color)],
                    plot_color=sort_frame$plot_color[match(frame$color, sort_frame$color)]
    )
    names(r) <- c("tcr","x","y","color","frequency","plot_color")
    return(r)
  }
  
  # Zoom after brush
  observeEvent(input$plot_brush,{
    brush <- input$plot_brush
    if (!is.null(brush)) {
      ranges_z$rzx <- c(brush$xmin, brush$xmax)
      ranges_z$rzy <- c(brush$ymin, brush$ymax)
    } else {
      ranges_z$rzx <- NULL
      ranges_z$rzy <- NULL
    }
  })
  observeEvent(input$select_data, { dirtyFlag <<- TRUE })
  observeEvent(input$select_method, { dirtyFlag <<- TRUE })
  observeEvent(input$tcr_file, { dirtyFlag <<- TRUE })
  
  # Construct plot with zooming function
  docu <- function(){
    # cat("Dirty flag = ", dirtyFlag, "\n")
    if (dirtyFlag) {
      res <- safe.read(input$tcr_file, input$select_data, input$select_method)
      cacheResult <<- res
      dirtyFlag <<- FALSE
    } else {
      # do nothing
    }
    return(cacheResult)
  }
  data <- function(){
    color.dist(docu())
  }
  color_table <- function(){
    data()[, c("color", "frequency", "plot_color")]
  }
  color_table1 <- function(){
    color_table()[!duplicated(color_table()),]
  }                                            # extract unduplicated data
  color_table2 <- function(){
    color_table1()[order(-color_table1()$frequency),]
  } 
  # Left plot with entire data
  output$plot_tcr <- renderPlot({
    data2 <- data()
    data2$color <- factor(data2$color, levels = color_table2()$color)  ## control the order of color in the figure legend
    color_table_mapping <- color_table2()$plot_color
    names(color_table_mapping) <- color_table2()$color
    ggplot(
      data2,
      aes(x = x, y = y, color = color)) + 
      geom_point(aes(color=color)) +
      scale_color_manual(values = color_table_mapping) +
      coord_cartesian(xlim = ranges$rx, ylim = ranges$ry, expand = TRUE) +
      ggtitle("Entire data")
  })
  
  output$info <- renderPrint({
    nearPoints(docu(),
               input$plot_click,
               xvar = "x",
               yvar = "y",
               addDist = TRUE
    )
  })
  # Right plot for zooming
  output$plot_tcr_z <- renderPlot({
    data2 <- data()
    data2$color <- factor(data2$color, levels = color_table2()$color)  ## control the order of color in the figure legend
    color_table_mapping <- color_table2()$plot_color
    names(color_table_mapping) <- color_table2()$color
    ggplot(
      data2,
      aes(x = x, y = y)) + 
      geom_point(aes(col=color)) +
      scale_color_manual(values = color_table_mapping) +
      coord_cartesian(xlim = ranges_z$rzx, ylim = ranges_z$rzy, expand = TRUE) +
      ggtitle("Chosen region") 
  })
  output$info_z <- renderPrint({
    nearPoints(docu(),
               input$plot_z_click,
               xvar = "x",
               yvar = "y",
               addDist = TRUE
    )
  })
  # Downloadable file for raw TCR seq, as an example
  output$download_example <- downloadHandler(
    filename = function() {
      paste("TCR_seq_example", "csv", sep = ".")
    },
    content = function(file1) {
      file.copy("TCR_seq_example.csv", file1)
    },
    contentType = "text/csv"
  )
  # Generate downloadable file with data frame generated by C++ or from input
  # Row names should be included in the data frame
  output$download_result <- downloadHandler(
    filename = function() {
      paste("TCR_analysis-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file2) {
      write.csv(docu(), 
                file2, row.names=FALSE)
    },
    contentType = "text/csv"
  )
  # Current time
  output$time <- renderText({
    invalidateLater(1000)
    format(Sys.time())
  })
}