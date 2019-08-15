# SurfaceGenie_0.1/functions.R
library(plyr)
library(stringr)
library(wordspace)
library(gplots)
library(RColorBrewer)
library(plotly)

##########  Genie Score Sub Functions  ##########

split_acc_iso <- function(protID) {
  return(unlist(strsplit(protID, "[-]"))[1])
}

get_SPC <- function(adata, Accession, species) {
  if(species=="human"){
    SPC_scores <- read.csv(file="ref/SPC.csv", header=TRUE)
  } else if(species=="rat"){
    SPC_scores <- read.csv(file="ref/Rat_SPC.csv", header=TRUE)
  } else if(species=="mouse") {
    SPC_scores <- read.csv(file="ref/Mouse_SPC.csv", header=TRUE)
  }
  noiso <- data.frame(Accession)
  noiso_SPC <- join(noiso, SPC_scores, by="Accession", type="left", match="first")
  noiso_SPC["SPC"][is.na(noiso_SPC["SPC"])]<-0
  adata["SPC"] <- noiso_SPC["SPC"]
  adata["noSPC"] <- matrix(rep(1, nrow(adata)))  
  return(adata)
}

get_Gini_coeff <- function(sdata, nsamps) {
  cmat <- matrix(rep(t(sdata), nsamps), ncol=nsamps, byrow=TRUE)
  sumdif <- sum(abs(cmat - t(cmat)))
  return(sumdif/(2*nsamps*sum(sdata)))
}

get_signal_strength <- function(sdata) {
  return(log10(max(sdata)+1))
}

get_dissimilar_score <- function(pdata, nsamps, spctype) {
  Gmax <- 1 - 1/nsamps
  if(spctype=="SPC"){
    return((pdata["Gini"]/Gmax)^2 * pdata["SPC"] * pdata["SS"])
  }else{
    return((pdata["Gini"]/Gmax)^2 * pdata["noSPC"] * pdata["SS"])
  }
}

get_similar_score <- function(pdata, nsamps, spctype) {
  Gmax <- 1 - 1/nsamps
  if(spctype=="SPC"){
    return(   ( ( 1-(pdata["Gini"]/Gmax)^2))   *   pdata["SPC"]   *   pdata["SS"]  )
  }else{
    return(   ( ( 1-(pdata["Gini"]/Gmax)^2))   *   pdata["noSPC"]   *   pdata["SS"]  )
  }
}

group_samples <- function(adata, groupmethod, groupcols){
  gtags <- c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5")
  numgroups <- length(groupcols)
  for(i in 1:numgroups){
    cols <- laply(strsplit(groupcols[[i]], ",")[[1]], as.integer)
    if(!(is.na(match("ave", groupmethod)))){
      adata[gtags[i]] = rowMeans(adata[cols])
    }
    if(!(is.na(match("med", groupmethod)))){
      for(j in 1:nrow(adata)) {
        v<-as.vector(t(adata[j,cols]))
        adata[gtags[i]] = median(v)
      }
    }
  }
  adata <- adata[c("Accession", gtags[1:numgroups])]
  return(adata)
}

append_UPL <- function(adata, Accession){
  baselink <- "https://www.uniprot.org/uniprot/"
  adata["UniProt Linkout"] <- laply(Accession, function(x) { paste(baselink, x, sep="") })
  return(adata)
}

get_CD <- function(adata, Accession, species) {
  if(species=="human") {
    CD <- read.csv("ref/Human_CD.csv", header=TRUE)
  } else if (species=="rat") {
    CD <- read.csv("ref/Rat_CD.csv", header=TRUE)
  } else if (species=="mouse") {
    CD <- read.csv("ref/Mouse_CD.csv", header=TRUE)
  }
  df <- data.frame(Accession)
  df <- join(df, CD, by="Accession", match="all")
  adata["CD"] <- df["CD"]
  return(adata)
}

get_numCSPA <- function(adata, Accession) {
  CSPA <- read.csv("ref/CSPA.csv", header=TRUE)
  df <- data.frame(Accession)
  df <- join(df, CSPA, by="Accession", match="all")
  adata["CSPA #e"] <- df["CSPA..e"]
  return(adata)
}

get_HLA <- function(adata, Accession, species) {
  if(species=="human") {
    HLA <- read.csv("ref/HLA.csv", header=TRUE)
  } else if (species=="rat") {
    HLA <- read.csv("ref/Rat_HLA.csv", header=TRUE)
  } else if (species=="mouse") {
    HLA <- read.csv("ref/Mouse_HLA.csv", header=TRUE)
  }
  df <- data.frame(Accession)
  df <- join(df, HLA, by="Accession", match="all")
  adata["HLA"] <- df["HLA"]
  return(adata)
}


get_geneName <- function(adata, Accession, species) {
  if(species=="human") {
    gn <- read.csv("ref/GeneName.csv", header=TRUE)
  } else if (species=="rat") {
    gn <- read.csv("ref/Rat_GeneName.csv", header=TRUE)
  } else if (species=="mouse") {
    gn <- read.csv("ref/Mouse_GeneName.csv", header=TRUE)
  }
  df <- data.frame(Accession)
  df <- join(df, gn, by="Accession", match="all")
  adata["geneName"] <- df["GeneName"]
  return(adata)
}

##########  Genie Score Main Function  ##########

SurfaceGenie <- function(adata, processing_opts, groupmethod, numgroups, groupcols, species, updateProgress) {
  accessions <- laply(laply(adata["Accession"], as.character), split_acc_iso)
  nsamps <- ncol(adata) - 1
  reqcols <- colnames(adata)
  # Sample grouping
  if("grouping" %in% processing_opts & numgroups > 1){
    adata <- group_samples(adata, groupmethod, groupcols)
    nsamps <- numgroups
    reqcols <- colnames(adata)
  }
  
  adata <- get_SPC(adata, accessions, species)
  adata <- get_CD(adata, accessions, species)
  adata <- get_HLA(adata, accessions, species)
  adata <- get_geneName(adata, accessions, species)

  # rather than pulling apart and inserting into a data frame and 
  # calling functions to calculate the GS, iGenie, etc scores, we 
  # will calculate them directly, drop them into a vector then
  # finally add the vectors to the data frame.  This has reduced 
  # the time for a 3700 line file from 14 secs to 3 (including the
  # annotation merge step)
  
  # set up vectors to eventually put in data frame
  Gmax <- 1 - 1/nsamps
  Ginis<-vector(mode="logical",length=0)
  SSs<-vector(mode="logical",length=0)
  GSs<-vector(mode="logical",length=0)
  eineGs<-vector(mode="logical",length=0)
  iGenies<-vector(mode="logical",length=0)
  eineGis<-vector(mode="logical",length=0)
  
  accessions <- laply(laply(adata["Accession"], as.character), split_acc_iso)  
  
  # loop through data set
  # update the progress meter
  rows = nrow(adata)
  for(irow in 1:rows){
    if(irow%%100==0){
      if (is.function(updateProgress)) {
        text <- paste0("row:", irow, "of ", rows)
        updateProgress(detail = text)
      }
    }
    
    #calculate the scores
    sdata <- adata[irow, 2:(nsamps + 1)]
    spc<-adata[irow,"SPC"]
    Gini <- get_Gini_coeff(sdata, nsamps)
    SS <- get_signal_strength(sdata)
    iGenie <-  (Gini/Gmax)^2 * SS
    GS <- iGenie * spc
    eineGi <- (1-(Gini/Gmax)^2) * SS
    eineG <-eineGi * spc
    
    # append the scores to each vector
    Ginis<-append(Ginis, Gini)
    SSs<-append(SSs, SS)
    iGenies<-append(iGenies, iGenie)
    GSs<-append(GSs, GS)
    eineGs<-append(eineGs, eineG)
    eineGis<-append(eineGis, eineGi)
  }
  
  # put the vectors into the data frame
  adata["Gini"]<-Ginis
  adata["SS"]<-SSs
  adata["GS"]<-GSs
  adata["eineG"]<-eineGs
  adata["iGenie"]<-iGenies
  adata["eineGi"]<-eineGis
  
  return(adata)
}

##########  SurfaceGenie Export  ##########

SG_export <- function(adata, exportvars1, exportvars2 , scoringvars, species) {
  
  accessions <- laply(laply(adata["Accession"], as.character), split_acc_iso)
  reqcols <- colnames(adata)[1:(ncol(adata)-11)]

  # Export option: append uniprot linkout column
  if("UniProt Linkout" %in% exportvars2){
    adata <- append_UPL(adata, accessions)
  }
  
  # Export option: append # CSPA experiments
  if("CSPA #e" %in% exportvars2){
    adata <- get_numCSPA(adata, accessions)
  }
  # Return data with export options as well as dataframe size
  return(adata[,c(reqcols, exportvars1, exportvars2, scoringvars)])
}

##########  SurfaceGenie Plots  ##########


SPC_hist <- function(adata) {
  scores <- adata[["SPC"]]
  counts <- count(scores)
  barplot(counts[["freq"]], names.arg=counts[["x"]], xlab="SPC Score", ylab="frequency", 
          main="SPC Score Histogram", col="#3498db", border="white")
}

SG_dist <- function(adata) {
  accessions <- laply(laply(adata["Accession"], as.character), split_acc_iso)
  adata <- adata[,c("Accession", "geneName", "GS", "CD")]
  adata <- adata[order(-adata$GS),]
  CD <- adata[,"CD"]
  df = data.frame(CD)
  df$CD<-as.character(df$CD)
  df$CD[!is.na(df$CD)]<-"CD"
  df$CD[is.na(df$CD)]<-"non-CD"
  adata["isCD"]<-df["CD"]
  fa<-list(family="Arial, sans-serif", size=12, color='black')
  ft<-list(family="Arial, sans-serif", size=14, color='black')
  plot_ly(data=adata,
          x=~1:nrow(adata),
          y=~GS, 
          type = 'scatter', 
          mode='markers', 
          hoverinfo = 'text', 
          hoverlabel = list(bgcolor='white'),
          text=paste("Gene Name: ", adata$geneName, "<br>Accession: ", adata$Accession, "<br>CD: ", adata$CD, "<br>GenieScore: ", adata$GS, "<br>Rank: ", 1:nrow(adata) ), 
          color=~isCD,
          colors=c("#3498db", "#c9c9d4") # blue, grey
  ) %>%
    layout(
      xaxis=list(title="rank", titlefont=fa, showgrid=FALSE),
      yaxis=list(title="Genie Score", titlefont=fa, showgrid=FALSE),
      title=list(text="<b>Genie Scores in Descending Order</b>", font=ft),
      legend=list(x=0.7,y=0.9) # controls the location on the plot of the legend
    )
}

SG_dist_export <- function(adata) {
  accessions <- laply(laply(adata["Accession"], as.character), split_acc_iso)
  adata <- adata[,c("GS", "CD")]
  adata <- adata[order(-adata$GS),]
  CD <- adata[,"GS"]
  CD[is.na(adata[,"CD"])] <- NA
  plot(1:nrow(adata), adata[,"GS"], xlab="rank", ylab="Genie Score", 
       main="Genie Scores in Descending Order",
       col="#C0C0C0")
  points(1:nrow(adata), CD, pch=16, col="#3498db")
  legend("topright", legend="CD molecules", col="#3498db", pch=16,
         inset=0.02, box.lwd=0)
}

eineG_dist <- function(adata) {
  accessions <- laply(laply(adata["Accession"], as.character), split_acc_iso)
  adata <- adata[,c("Accession", "geneName", "eineG", "CD")]
  adata <- adata[order(-adata$eineG),]
  CD <- adata[,"CD"]
  df = data.frame(CD)
  df$CD<-as.character(df$CD)
  df$CD[!is.na(df$CD)]<-"CD"
  df$CD[is.na(df$CD)]<-"non-CD"
  adata["isCD"]<-df["CD"]
  fa<-list(family="Arial, sans-serif", size=12, color='black')
  ft<-list(family="Arial, sans-serif", size=14, color='black')
  plot_ly(data=adata,
          x=~1:nrow(adata),
          y=~eineG, 
          type = 'scatter', 
          mode='markers', 
          hoverinfo = 'text', 
          hoverlabel = list(bgcolor='white'),
          text=paste("Gene Name: ", adata$geneName, "<br>Accession: ", adata$Accession, "<br>CD: ", adata$CD, "<br>IsoGenieScore: ", adata$eineG, "<br>Rank: ", 1:nrow(adata)), 
          color=~isCD,
          colors=c("#3498db", "#c9c9d4") # blue, grey
  ) %>%
    layout(
      title=list(text="<b>IsoGenie Scores in Descending Order</b>", font=ft),
      xaxis=list(title="rank", titlefont=fa, showgrid=FALSE),
      yaxis=list(title="IsoGenie Score", titlefont=fa, showgrid=FALSE),
      legend=list(x=0.7,y=0.9) # controls the location on the plot of the legend
    )
}

eineG_dist_export <- function(adata) {
  accessions <- laply(laply(adata["Accession"], as.character), split_acc_iso)
  adata <- adata[,c("eineG", "CD")]
  adata <- adata[order(-adata$eineG),]
  CD <- adata[,"eineG"]
  CD[is.na(adata[,"CD"])] <- NA
  plot(1:nrow(adata), adata[,"eineG"], xlab="rank", ylab="IsoGenie Score", 
       main="IsoGenie Scores in Descending Order",
       col="#C0C0C0")
  points(1:nrow(adata), CD, pch=16, col="#3498db")
  legend("topright", legend="CD molecules", col="#3498db", pch=16,
         inset=0.02, box.lwd=0)
}

iGenie_dist <- function(adata) {
  accessions <- laply(laply(adata["Accession"], as.character), split_acc_iso)
  adata <- adata[,c("Accession", "geneName", "iGenie", "CD")]
  adata <- adata[order(-adata$iGenie),]
  CD <- adata[,"CD"]
  df = data.frame(CD)
  df$CD<-as.character(df$CD)
  df$CD[!is.na(df$CD)]<-"CD"
  df$CD[is.na(df$CD)]<-"non-CD"
  adata["isCD"]<-df["CD"]
  fa<-list(family="Arial, sans-serif", size=12, color='black')
  ft<-list(family="Arial, sans-serif", size=14, color='black')
  plot_ly(data=adata,
          x=~1:nrow(adata),
          y=~iGenie, 
          type = 'scatter', 
          mode='markers', 
          hoverinfo = 'text', 
          hoverlabel = list(bgcolor='white'),
          text=paste("Gene Name: ", adata$geneName, "<br>Accession: ", adata$Accession, "<br>CD: ", adata$CD, "<br>OmniGenieScore: ", adata$iGenie, "<br>Rank: ", 1:nrow(adata)), 
          color=~isCD,
          colors=c("#3498db", "#c9c9d4") # blue, grey
  ) %>%
    layout(
      title=list(text="<b>OmniGenie Scores in Descending Order</b>", font=ft),
      xaxis=list(title="rank", titlefont=fa, showgrid=FALSE),
      yaxis=list(title="OmniGenie Score", titlefont=fa, showgrid=FALSE),
      legend=list(x=0.7,y=0.9) # controls the location on the plot of the legend
    )
}

iGenie_dist_export <- function(adata) {
  accessions <- laply(laply(adata["Accession"], as.character), split_acc_iso)
  adata <- adata[,c("iGenie", "CD")]
  adata <- adata[order(-adata$iGenie),]
  CD <- adata[,"iGenie"]
  CD[is.na(adata[,"CD"])] <- NA
  plot(1:nrow(adata), adata[,"iGenie"], xlab="rank", ylab="OmniGenie Score", 
       main="OmniGenie Scores in Descending Order",
       col="#C0C0C0")
  points(1:nrow(adata), CD, pch=16, col="#3498db")
  legend("topright", legend="CD molecules", col="#3498db", pch=16,
         inset=0.02, box.lwd=0)
}


eineGi_dist <- function(adata) {
  accessions <- laply(laply(adata["Accession"], as.character), split_acc_iso)
  adata <- adata[,c("Accession", "geneName", "eineGi", "CD")]
  adata <- adata[order(-adata$eineGi),]
  CD <- adata[,"CD"]
  df = data.frame(CD)
  df$CD<-as.character(df$CD)
  df$CD[!is.na(df$CD)]<-"CD"
  df$CD[is.na(df$CD)]<-"non-CD"
  adata["isCD"]<-df["CD"]
  fa<-list(family="Arial, sans-serif", size=12, color='black')
  ft<-list(family="Arial, sans-serif", size=14, color='black')
  plot_ly(data=adata,
          x=~1:nrow(adata),
          y=~eineGi, 
          type = 'scatter', 
          mode='markers', 
          hoverinfo = 'text', 
          hoverlabel = list(bgcolor='white'),
          text=paste("Gene Name: ", adata$geneName, "<br>Accession: ", adata$Accession, "<br>CD: ", adata$CD, "<br>IsoOmniGenieScore: ", adata$eineGi, "<br>Rank: ", 1:nrow(adata)), 
          color=~isCD,
          colors=c("#3498db", "#c9c9d4") # blue, grey
  ) %>%
    layout(
      title=list(text="<b>IsoOmniGenie Scores in Descending Order</b>", font=ft),
      xaxis=list(title="rank", titlefont=fa, showgrid=FALSE),
      yaxis=list(title="IsoOmniGenie Score", titlefont=fa, showgrid=FALSE),
      legend=list(x=0.7,y=0.9) # controls the location on the plot of the legend
    )
}

eineGi_dist_export <- function(adata) {
  accessions <- laply(laply(adata["Accession"], as.character), split_acc_iso)
  adata <- adata[,c("eineGi", "CD")]
  adata <- adata[order(-adata$eineGi),]
  CD <- adata[,"eineGi"]
  CD[is.na(adata[,"CD"])] <- NA
  plot(1:nrow(adata), adata[,"eineGi"], xlab="rank", ylab="IsoOmns Score", 
       main="IsoOmniGenie Scores in Descending Order",
       col="#C0C0C0")
  points(1:nrow(adata), CD, pch=16, col="#3498db")
  legend("topright", legend="CD molecules", col="#3498db", pch=16,
         inset=0.02, box.lwd=0)
}


##########  SPC Lookup  ##########

SPC_lookup <- function(sdata) {
  if("Accession" %in% colnames(sdata)){
    Accession <- laply(laply(sdata["Accession"], as.character), split_acc_iso)
  }
  else{
    Accession <- laply(laply(sdata[,1], as.character), split_acc_iso)
  }
  SPC_scores <- read.csv(file="ref/SPC_by_Source.csv", header=TRUE)
  noiso <- data.frame(Accession)
  noiso_SPC <- join(noiso, SPC_scores, by="Accession", match="first")
  sdata["SPC"] <- noiso_SPC["SPC"]
  sdata["SURFY"] <- noiso_SPC["SURFY"]
  sdata["Town"] <- noiso_SPC["Town"]
  sdata["Cunha"] <- noiso_SPC["Cunha"]
  sdata["Diaz-Ramos"] <- noiso_SPC["Diaz.Ramos"]
  sdata[,3:6][is.na(sdata[,3:6])] = " "
  sdata[,3:6][sdata[,3:6]>0] = "\u00A0\u00A0\u00A0\u00A0\u00A0\u2713"
  sdata[,3:6][sdata[,3:6]=="0"] = " "
  return(sdata)
}

SPC_lookup_for_export <- function(sdata) {
  if("Accession" %in% colnames(sdata)){
    Accession <- laply(laply(sdata["Accession"], as.character), split_acc_iso)
  }
  else{
    Accession <- laply(laply(sdata[,1], as.character), split_acc_iso)
  }
  SPC_scores <- read.csv(file="ref/SPC_by_Source.csv", header=TRUE)
  noiso <- data.frame(Accession)
  noiso_SPC <- join(noiso, SPC_scores, by="Accession", match="first")
  sdata["SPC"] <- noiso_SPC["SPC"]
  sdata["SURFY"] <- noiso_SPC["SURFY"]
  sdata["Town"] <- noiso_SPC["Town"]
  sdata["Cunha"] <- noiso_SPC["Cunha"]
  sdata["Diaz-Ramos"] <- noiso_SPC["Diaz.Ramos"]
  sdata[,3:6][is.na(sdata[,3:6])] = 0
  sdata[,3:6][sdata[,3:6]>0] = 1
  sdata[,3:6][sdata[,3:6]=="0"] = 0
  return(sdata)
}
