# SurfaceGenie_0.1/server.R
library(shiny)
library(plotly)
library(RColorBrewer) 
library(stringr)
#library(ggplot2)
library(svglite)
library(xlsx)
source("functions.R")
`%then%` <- shiny:::`%OR%`


function(input, output, session) {
  
  ##########          SurfaceGenie         ##########
  
  # Load and process data
  data_input <- reactive({
    withProgress(message = 'Reading Data', value = 0, {
      parts <- strsplit(input$file1$name, "\\.")[[1]]
      ext <- parts[length(parts)]
      if(ext == "csv"){
        df <- read.csv(input$file1$datapath, header=TRUE)
      } else if(ext=="txt" || ext=="tsv" || ext=="tab") {
        df <- read.delim(input$file1$datapath, header=TRUE)
      } else if(ext=="xlsx" || ext=="xls") {
        df <- read.xlsx(input$file1$datapath,  1)
      } else {
        validate(need(ext=="csv" || ext=="tsv" || ext=="xlsx" || ext=="xls" || ext=="tab" || ext=="txt", "You have the wrong file extension.  Please see the instructions for possible extensions and associated file types." ))
      }
      
      # validation
      validate(
        need(ncol(df)>1,"There is only one column of data.  This is probably caused by one of two things.  Either you used the wrong field delimiter - e.g. tabs instead of commas - or you need to enter one or more samples and their abundance data.")
      )

      # validation
      validate(
        need(ncol(df)>1,"There is only one column of data.  This is probably caused by one of two things.  Either you used the wrong field delimiter - e.g. tabs instead of commas - or you need to enter one or more samples and their abundance data."),
        need(length(which(as.vector(sapply(df[2:3],is.numeric))==FALSE))==0,"There is a non-numeric value in the abundance data columns.  The first column should be accessions and each column after that should be numbers representing the abundance of the protein for that sample.  Please refer to the instructions for more help.")
      )

      ## look for duplicates and for NA values.  These errors will not stop the program, but we will inform the user.
      # look for duplicates
      warningMsg <- ""
      e <- as.vector(df[duplicated(df),][,1])
      if( length(e)>0 ){
        warningMsg <- "You have duplicates in the list of accessions.  You may wish to address that."
      }
      # look for NA values.  This is not necessarily wrong, but want to make the user informed.
      if(length(which(is.na(df)))>0) {
        warningMsg <- c(warningMsg, "There are missing values in your data.  You may wish to address that.")
      }
      output$txtWarning<-renderText(warningMsg)

      input_size <- c(nrow(df), ncol(df))
      list(df, input_size)
    })
  })
  

  # Load annotation file
  annotation <- reactive({
#    print("Annotation")
    withProgress(message = 'Reading Annotations', value = 0, {
      if(input$species == "Human") {
        read.delim("ref/anno.tsv")
      } else if(input$species == "Rat") {
        read.delim("ref/annotation.rat.tsv")
      } else if(input$species == "Mouse") {
        read.delim("ref/annotation.mouse.tsv")
      } else if(input$species == "Other/Ignore") {
        read.delim("ref/annotation.none.tsv")
      }
    })
  })
    
  
  
  ##########  SurfaceGenie: Data Grouping  ##########
  
  data_output <- reactive({

    # setting up the progress meter
    progress<-shiny::Progress$new()
    progress$set(message = "Calculating Scores", value=0)
    on.exit(progress$close())
    
    # creating the progress function to pass to the Surface Genie function
    updateProgress <- function(value=NULL,detail=NULL) {
      if(is.null(value)){
        value <- progress$getValue()
        value <- value + (progress$getMax()-value)/10
      }
      progress$set(value=value, detail=detail)
    }

        
    if("grouping" %in% input$processing_opts){
      gtags <- c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5")
      groupcols <- list()
      if(input$numgroups >= 2){
        validate(
          need(input$group1, "Please indicate columns in Group 1") %then%
            need(laply(strsplit(input$group1, ","), as.integer), 
                 "Group 1 Error: Non-integer values have been entered.") %then%
            need(!(1 %in% strsplit(input$group1,",")[[1]]),
                 "Group 1 Error: Column 1 should contain accession numbers and can not be grouped" ),
          need(input$group2, "Please indicate columns in Group 2") %then%
            need(laply(strsplit(input$group2, ",")[[1]], as.integer), 
                 "Group 2 Error: Non-integer values have been entered.") %then%
            need(!(1 %in% strsplit(input$group2,",")[[1]]),
                 "Group 2 Error: Column 1 should contain accession numbers and can not be grouped" )
        )
        groupcols[1] <- input$group1
        groupcols[2] <- input$group2
      }
      if(input$numgroups >= 3){
        validate(
          need(input$group3, "Please indicate columns in Group 3") %then%
            need(laply(strsplit(input$group3, ",")[[1]], as.integer), 
                 "Group 3 Error: Non-integer values have been entered.") %then%
            need(!(1 %in% strsplit(input$group3,",")[[1]]),
                 "Group 3 Error: Column 1 should contain accession numbers and can not be grouped" )
        )
        groupcols[3] <- input$group3
      }
      if(input$numgroups >= 4){
        validate(
          need(input$group4, "Please indicate columns in Group 4") %then%
            need(laply(strsplit(input$group4, ",")[[1]], as.integer), 
                 "Group 4 Error: Non-integer values have been entered.") %then%
            need(!(1 %in% strsplit(input$group4,",")[[1]]),
                 "Group 4 Error: Column 1 should contain accession numbers and can not be grouped" )
        )
        groupcols[4] <- input$group4
      }
      if(input$numgroups >=5){
        validate(
          need(input$group5, "Please indicate columns in Group 5") %then%
            need(laply(strsplit(input$group5, ",")[[1]], as.integer), 
                 "Group 5 Error: Non-integer values have been entered.") %then%
            need(!(1 %in% strsplit(input$group5,",")[[1]]),
                 "Group 5 Error: Column 1 should contain accession numbers and can not be grouped" )
        )
        groupcols[5] <- input$group5
      }
      # call this if grouped
      sg<-SurfaceGenie(data_input()[[1]], input$processing_opts, 
                   input$groupmethod, input$numgroups, groupcols, annotation(), updateProgress)
      
    }
    else{
      # call this if not grouped
      sg<-SurfaceGenie(data_input()[[1]], input$processing_opts, 
                   groupmethod=NULL, numgroups=0, groupcols=NULL, annotation(), updateProgress)
    }
    if(length(which(as.vector(sapply(sg["SPC"],is.na))==TRUE))>0) {
      warningMsg <- "There are accessions that do not belong to the most recent swissprot version for the species you selected. They will be ignored for SurfaceGenie Score and IsoGenie Score. Download the csv file to find which ones."
      output$txtWarning<-renderText(warningMsg)
    }
    return(sg)
  })
  
  ##########  SurfaceGenie: Output Display  ##########
  
  # Apply export options
  data_export <- reactive({
    if(input$species == "Human") {
      df <- SG_export(data_output(), input$export_options1, input$export_options2h, input$scoring_opts)
    } else if(input$species == "Rat") {
      df <- SG_export(data_output(), input$export_options1, input$export_options2r, input$scoring_opts)
    } else if(input$species == "Mouse") {
      df <- SG_export(data_output(), input$export_options1, input$export_options2m, input$scoring_opts)
    } else if(input$species == "Other/Ignore") {
      df <- SG_export(data_output(), input$export_options1o, input$export_options2, input$scoring_opts_o)
    }
          
    output_size <- c(nrow(df), ncol(df))
    list(df, output_size)
  })
  
  # Example data to display
  example_data <- read.csv("ref/example_data.csv", header=TRUE)
  output$example_data <- renderTable(head(example_data, 5))
  
  # Display input data
  output$data_input <- renderTable({
    req(input$file1)
    head(data_input()[[1]], 10)
  })
  output$input_size <- renderText({
    req(input$file1)
    sprintf("[ %d rows x %d columns ]", data_input()[[2]][1], data_input()[[2]][2])
  })
  
  # Display output data
  output$data_output <- renderTable({
    req(input$file1)
    head(data_export()[[1]], 10)
  })
  output$output_size <- renderText({
    req(input$file1)
    sprintf("[ %d rows x %d columns ]", data_export()[[2]][1], data_export()[[2]][2])
  })

  # SG: SPC histogram
  output$SG_SPC_hist <- renderPlot({
    req(input$file1)
    SPC_hist(data_output())
  })
  
  output$SG_SPC_hist_PNGdl <- downloadHandler(
    filename = function() { 
      fname <- unlist(strsplit(as.character(input$file1), "[.]"))[1]
      paste(fname, '_SG_SPC_hist.png', sep='') 
    },
    content = function(file) {
      png(file,
          width=5*300,
          height=3*300,
          res=300)
      SPC_hist(data_output())
      dev.off()
    }
  )
  output$SG_SPC_hist_PNGdlbutton <- renderUI({
    req(input$file1)
    downloadButton("SG_SPC_hist_PNGdl", " .png", class="download_this")
  })
  output$SG_SPC_hist_SVGdl <- downloadHandler(
    filename = function() { 
      fname <- unlist(strsplit(as.character(input$file1), "[.]"))[1]
      paste(fname, '_SG_SPC_hist.svg', sep='') 
    },
    content = function(file) {
      svg(file,
          width=5,
          height=3,
          pointsize=10)
      SPC_hist(data_output())
      dev.off()
    }
  )
  output$SG_SPC_hist_SVGdlbutton <- renderUI({
    req(input$file1)
    downloadButton("SG_SPC_hist_SVGdl", " .svg", class="download_this")
  })

    # SG: Genie Score plot
  output$SG_dist <- renderPlotly({
    req(input$file1)
    SG_dist(data_output())
  })
  
  output$SG_dist_PNGdl <- downloadHandler(
    filename = function() { 
      fname <- unlist(strsplit(as.character(input$file1), "[.]"))[1]
      paste(fname, '_SG_dist.png', sep='') 
    },
    content = function(file) {
      png(file,
          width=5*300,
          height=4*300,
          res=300)
      #      SG_dist(data_output())
      SG_dist_export(data_output())
      dev.off()
    }
  )
  output$SG_dist_PNGdlbutton <- renderUI({
    req(input$file1)
    downloadButton("SG_dist_PNGdl", " .png", class="download_this")
  })
  output$SG_dist_SVGdl <- downloadHandler(
    filename = function() { 
      fname <- unlist(strsplit(as.character(input$file1), "[.]"))[1]
      paste(fname, '_SG_dist.svg', sep='') 
    },
    content = function(file) {
      svg(file,
          width=5,
          height=4,
          pointsize=10)
      SG_dist_export(data_output())
      dev.off()
    }
  )
  output$SG_dist_SVGdlbutton <- renderUI({
    req(input$file1)
    downloadButton("SG_dist_SVGdl", " .svg", class="download_this")
  })
  
  # IsoGenie: reverse Genie Score plot
  output$IsoGenie_dist <- renderPlotly({
    req(input$file1)
    IsoGenie_dist(data_output())
  })
  
  output$IsoGenie_dist_PNGdl <- downloadHandler(
    filename = function() { 
      fname <- unlist(strsplit(as.character(input$file1), "[.]"))[1]
      paste(fname, '_IsoGenie_dist.png', sep='') 
    },
    content = function(file) {
      png(file,
          width=5*300,
          height=4*300,
          res=300)
      #      IsoGenie_dist(data_output())
      IsoGenie_dist_export(data_output())
      dev.off()
    }
  )
  output$IsoGenie_dist_PNGdlbutton <- renderUI({
    req(input$file1)
    downloadButton("IsoGenie_dist_PNGdl", " .png", class="download_this")
  })
  output$IsoGenie_dist_SVGdl <- downloadHandler(
    filename = function() { 
      fname <- unlist(strsplit(as.character(input$file1), "[.]"))[1]
      paste(fname, '_IsoGenie_dist.svg', sep='') 
    },
    content = function(file) {
      svg(file,
          width=5,
          height=4,
          pointsize=10)
      IsoGenie_dist_export(data_output())
      dev.off()
    }
  )
  output$IsoGenie_dist_SVGdlbutton <- renderUI({
    req(input$file1)
    downloadButton("IsoGenie_dist_SVGdl", " .svg", class="download_this")
  })

  # OmniGenie: SurfaceGenie w/o SPC Score plot
  output$OmniGenie_dist <- renderPlotly({
    req(input$file1)
    OmniGenie_dist(data_output())
  })
  
  output$OmniGenie_dist_PNGdl <- downloadHandler(
    filename = function() { 
      fname <- unlist(strsplit(as.character(input$file1), "[.]"))[1]
      paste(fname, '_OmniGenie_dist.png', sep='') 
    },
    content = function(file) {
      png(file,
          width=5*300,
          height=4*300,
          res=300)
      #      OmniGenie_dist(data_output())
      OmniGenie_dist_export(data_output())
      dev.off()
    }
  )
  output$OmniGenie_dist_PNGdlbutton <- renderUI({
    req(input$file1)
    downloadButton("OmniGenie_dist_PNGdl", " .png", class="download_this")
  })
  output$OmniGenie_dist_SVGdl <- downloadHandler(
    filename = function() { 
      fname <- unlist(strsplit(as.character(input$file1), "[.]"))[1]
      paste(fname, '_OmniGenie_dist.svg', sep='') 
    },
    content = function(file) {
      svg(file,
          width=5,
          height=4,
          pointsize=10)
      OmniGenie_dist_export(data_output())
      dev.off()
    }
  )
  output$OmniGenie_dist_SVGdlbutton <- renderUI({
    req(input$file1)
    downloadButton("OmniGenie_dist_SVGdl", " .svg", class="download_this")
  })
  
  
  # IsoOmniGenie: Similarity w/o SPC score plot
  output$IsoOmniGenie_dist <- renderPlotly({
    req(input$file1)
    IsoOmniGenie_dist(data_output())
  })
  
  output$IsoOmniGenie_dist_PNGdl <- downloadHandler(
    filename = function() { 
      fname <- unlist(strsplit(as.character(input$file1), "[.]"))[1]
      paste(fname, '_IsoOmniGenie_dist.png', sep='') 
    },
    content = function(file) {
      png(file,
          width=5*300,
          height=4*300,
          res=300)
      #      IsoOmniGenie_dist(data_output())
      IsoOmniGenie_dist_export(data_output())
      dev.off()
    }
  )
  output$IsoOmniGenie_dist_PNGdlbutton <- renderUI({
    req(input$file1)
    downloadButton("IsoOmniGenie_dist_PNGdl", " .png", class="download_this")
  })
  output$IsoOmniGenie_dist_SVGdl <- downloadHandler(
    filename = function() { 
      fname <- unlist(strsplit(as.character(input$file1), "[.]"))[1]
      paste(fname, '_IsoOmniGenie_dist.svg', sep='') 
    },
    content = function(file) {
      svg(file,
          width=5,
          height=4,
          pointsize=10)
      IsoOmniGenie_dist_export(data_output())
      dev.off()
    }
  )
  output$IsoOmniGenie_dist_SVGdlbutton <- renderUI({
    req(input$file1)
    downloadButton("IsoOmniGenie_dist_SVGdl", " .svg", class="download_this")
  })
  

  # Downloadable csv of selected dataset
  output$csv_download <- downloadHandler(
    filename = function() {
      fname <- unlist(strsplit(as.character(input$file1), "[.]"))[1]
      paste(fname, "_SurfaceGenie.csv", sep = "")
    },
    content = function(filename) {
      write.csv(data_export()[[1]], filename, row.names = FALSE)
    }
  )
  output$csv_dlbutton <- renderUI({
    req(input$file1)
    downloadButton("csv_download", ".csv", class="download_this")
  })
  
  # Downloadable tsv of selected dataset
  output$tsv_download <- downloadHandler(
    filename = function() {
      fname <- unlist(strsplit(as.character(input$file1), "[.]"))[1]
      paste(fname, "_SurfaceGenie.tsv", sep = "")
    },
    content = function(filename) {
      write.table(data_export()[[1]], filename, sep="\t", row.names = FALSE)
    }
  )
  output$tsv_dlbutton <- renderUI({
    req(input$file1)
    downloadButton("tsv_download", ".tsv", class="download_this")
  })

  # Downloadable xlsx of selected dataset
  output$xlsx_download <- downloadHandler(
    filename = function() {
      fname <- unlist(strsplit(as.character(input$file1), "[.]"))[1]
      paste(fname, "_SurfaceGenie.xlsx", sep = "")
    },
      content = function(filename) {
      write.xlsx(data_export()[[1]], filename, row.names = FALSE)
    }
  )
  output$xlsx_dlbutton <- renderUI({
    req(input$file1)
    downloadButton("xlsx_download", ".xlsx", class="download_this")
  })
  
  ##########      SPC Quick Lookup      ##########
  
  # Load and process data
  SPC_quick_input <- reactive({
    numprots <- str_count(input$quicklookup, pattern="\\n") + 1
    validate(
      need(numprots <= 100, "Too many accession numbers! (Max: 100)")
    )
    proteins = data.frame(Accession=rep(0, numprots))
    for (p in 1:numprots)
    {
      proteins[p, "Accession"] <- strsplit(input$quicklookup, "\\n")[[1]][p]
    }
    return(proteins)
  })
  SPC_quick_output <- reactive({
    SPC_lookup(SPC_quick_input(), input$species2)
  })
  
  # Display data
  output$SPC_quick_output <- renderTable({
    req(input$quicklookup)
    SPC_quick_output()
  })
  output$SPC_quick_hist <- renderPlot({
    req(input$quicklookup)
    SPC_hist(SPC_quick_output())
  })
  
  
  ##########      SPC Bulk Lookup     ##########
  
  # Load and process data
  SPC_bulk_input <- reactive({
    withProgress(message = 'Reading Data', value = 0, {
      parts <- strsplit(input$file2$name, "\\.")[[1]]
      ext <- parts[length(parts)]
      if(ext == "csv"){
        df <- read.csv(input$file2$datapath, header=TRUE)
      } else if(ext=="txt" || ext=="tsv" || ext=="tab") {
        df <- read.delim(input$file2$datapath, header=TRUE)
      } else if(ext=="xlsx" || ext=="xls") {
        df <- read.xlsx(input$file2$datapath,  1)
      } else {
        validate(need(ext=="csv" || ext=="tsv" || ext=="xlsx" || ext=="xls" || ext=="tab" || ext=="txt", "You have the wrong file extension" ))
      }
      df
    })
  })
  
  SPC_bulk_output <- reactive({
    SPC_lookup(SPC_bulk_input(), input$species2)
  })
  SPC_bulk_output_for_export <- reactive({
    SPC_lookup_for_export(SPC_bulk_input(), input$species2)
  })
  
  # Display data
  output$SPC_bulk_output <- renderTable({
    req(input$file2)
    head(SPC_bulk_output(), 10)
  })
  output$SPC_bulk_hist <- renderPlot({
    req(input$file2)
    SPC_hist(SPC_bulk_output())
  })
  
  # Downloadable csv of selected dataset ----
  output$SPC_csv_download <- downloadHandler(
    filename = function() {
      fname <- unlist(strsplit(as.character(input$file2), "[.]"))[1]
      paste(fname, "_SPC.csv", sep = "")
    },
    content = function(filename) {
      write.csv(SPC_bulk_output_for_export(), filename, row.names = FALSE)
    }
  )
  output$SPC_csv_dlbutton <- renderUI({
    req(input$file2)
    downloadButton("SPC_csv_download", "Download .csv file output")
  })
  
}
