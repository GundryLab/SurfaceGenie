# SurfaceGenie_0.1/server.R
library(shiny)
library(plotly)
library(RColorBrewer) 
library(stringr)
#library(ggplot2)
library(svglite)
source("functions.R")
`%then%` <- shiny:::`%OR%`


function(input, output, session) {
  
  ##########          SurfaceGenie         ##########
  
  # Load and process data
  data_input <- reactive({
    df <- read.csv(input$file1$datapath, header=TRUE)
    input_size <- c(nrow(df), ncol(df))
    list(df, input_size)
  })
  
  ##########  SurfaceGenie: Data Grouping  ##########
  
  data_output <- reactive({
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
      SurfaceGenie(data_input()[[1]], input$processing_opts, 
                   input$groupmethod, input$numgroups, groupcols, input$species)
      
    }
    else{
      # call this if not grouped
      SurfaceGenie(data_input()[[1]], input$processing_opts, 
                   groupmethod=NULL, numgroups=0, groupcols=NULL, input$species)
    }
  })
  
  ##########  SurfaceGenie: Output Display  ##########
  
  myChoiceNames = list(
    "SPC score (SPC)",
    "Exclude HLA molecules",
    "CD molecules",
    "Number of CSPA experiments",
    "Gini coefficient (Gini)",
    "Signal strength (SS)",
#    "SurfaceGenie: Genie Score (GS)",
    "UniProt Linkout")
  myChoiceValues= list(
#    "SPC", "HLA", "CD", "CSPA #e", "Gini", "SS", "GS", "UniProt Linkout")
    "SPC", "HLA", "CD", "CSPA #e", "Gini", "SS", "UniProt Linkout")

  observe({
    updateCheckboxGroupInput(
      session, 'export_options', choiceNames=myChoiceNames, choiceValues = myChoiceValues, 
      selected = if (input$bar) myChoiceValues
    )
  })

  # Apply export options
  data_export <- reactive({
    df <- SG_export(data_output(), input$export_options, input$scoring_opts, input$species)
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
  
  # eineG: reverse Genie Score plot
  output$eineG_dist <- renderPlotly({
    req(input$file1)
    eineG_dist(data_output())
  })
  
  output$eineG_dist_PNGdl <- downloadHandler(
    filename = function() { 
      fname <- unlist(strsplit(as.character(input$file1), "[.]"))[1]
      paste(fname, '_eineG_dist.png', sep='') 
    },
    content = function(file) {
      png(file,
          width=5*300,
          height=4*300,
          res=300)
      #      eineG_dist(data_output())
      eineG_dist_export(data_output())
      dev.off()
    }
  )
  output$eineG_dist_PNGdlbutton <- renderUI({
    req(input$file1)
    downloadButton("eineG_dist_PNGdl", " .png", class="download_this")
  })
  output$eineG_dist_SVGdl <- downloadHandler(
    filename = function() { 
      fname <- unlist(strsplit(as.character(input$file1), "[.]"))[1]
      paste(fname, '_eineG_dist.svg', sep='') 
    },
    content = function(file) {
      svg(file,
          width=5,
          height=4,
          pointsize=10)
      eineG_dist_export(data_output())
      dev.off()
    }
  )
  output$eineG_dist_SVGdlbutton <- renderUI({
    req(input$file1)
    downloadButton("eineG_dist_SVGdl", " .svg", class="download_this")
  })

  # iGenie: SurfaceGenie w/o SPC Score plot
  output$iGenie_dist <- renderPlotly({
    req(input$file1)
    iGenie_dist(data_output())
  })
  
  output$iGenie_dist_PNGdl <- downloadHandler(
    filename = function() { 
      fname <- unlist(strsplit(as.character(input$file1), "[.]"))[1]
      paste(fname, '_iGenie_dist.png', sep='') 
    },
    content = function(file) {
      png(file,
          width=5*300,
          height=4*300,
          res=300)
      #      iGenie_dist(data_output())
      iGenie_dist_export(data_output())
      dev.off()
    }
  )
  output$iGenie_dist_PNGdlbutton <- renderUI({
    req(input$file1)
    downloadButton("iGenie_dist_PNGdl", " .png", class="download_this")
  })
  output$iGenie_dist_SVGdl <- downloadHandler(
    filename = function() { 
      fname <- unlist(strsplit(as.character(input$file1), "[.]"))[1]
      paste(fname, '_iGenie_dist.svg', sep='') 
    },
    content = function(file) {
      svg(file,
          width=5,
          height=4,
          pointsize=10)
      iGenie_dist_export(data_output())
      dev.off()
    }
  )
  output$iGenie_dist_SVGdlbutton <- renderUI({
    req(input$file1)
    downloadButton("iGenie_dist_SVGdl", " .svg", class="download_this")
  })
  
  
  # eineGi: Similarity w/o SPC score plot
  output$eineGi_dist <- renderPlotly({
    req(input$file1)
    eineGi_dist(data_output())
  })
  
  output$eineGi_dist_PNGdl <- downloadHandler(
    filename = function() { 
      fname <- unlist(strsplit(as.character(input$file1), "[.]"))[1]
      paste(fname, '_eineGi_dist.png', sep='') 
    },
    content = function(file) {
      png(file,
          width=5*300,
          height=4*300,
          res=300)
      #      eineGi_dist(data_output())
      eineGi_dist_export(data_output())
      dev.off()
    }
  )
  output$eineGi_dist_PNGdlbutton <- renderUI({
    req(input$file1)
    downloadButton("eineGi_dist_PNGdl", " .png", class="download_this")
  })
  output$eineGi_dist_SVGdl <- downloadHandler(
    filename = function() { 
      fname <- unlist(strsplit(as.character(input$file1), "[.]"))[1]
      paste(fname, '_eineGi_dist.svg', sep='') 
    },
    content = function(file) {
      svg(file,
          width=5,
          height=4,
          pointsize=10)
      eineGi_dist_export(data_output())
      dev.off()
    }
  )
  output$eineGi_dist_SVGdlbutton <- renderUI({
    req(input$file1)
    downloadButton("eineGi_dist_SVGdl", " .svg", class="download_this")
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
    #    downloadButton("csv_download", "Download .csv file output")
    downloadButton("csv_download", " .csv")
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
    SPC_lookup(SPC_quick_input())
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
    read.csv(input$file2$datapath, header=TRUE)
  })
  SPC_bulk_output <- reactive({
    SPC_lookup(SPC_bulk_input())
  })
  SPC_bulk_output_for_export <- reactive({
    SPC_lookup_for_export(SPC_bulk_input())
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
