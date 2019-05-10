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

get_SPC <- function(adata, Accession) {
  SPC_scores <- read.csv(file="ref/SPC.csv", header=TRUE)
  noiso <- data.frame(Accession)
  noiso_SPC <- join(noiso, SPC_scores, by="Accession", type="left", match="first")
  noiso_SPC["SPC"][is.na(noiso_SPC["SPC"])]<-0
  adata["SPC"] <- noiso_SPC["SPC"]
  #  adata["SPCdisplay"] <- noiso_SPC["SPC"]
  adata["noSPC"] <- matrix(rep(1, nrow(adata)))  
  #  idx <- sapply(adata["SPC"] > 0, isTRUE)
  #  return(adata[idx,])
  return(adata)
}

filter_by_HLA <- function(adata, Accession) {
  HLA_molecs <- read.csv(file="ref/HLA.csv", header=TRUE)
  noiso <- data.frame(Accession)
  noiso_HLA <- join(noiso, HLA_molecs, by="Accession", match="all")
  idx <- sapply(is.na(noiso_HLA["HLA"]), isTRUE)
  return(adata[idx,])
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

get_CD <- function(adata, Accession) {
  CD <- read.csv("ref/CD.csv", header=TRUE)
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
get_geneName <- function(adata, Accession) {
  gn <- read.csv("ref/GeneName.csv", header=TRUE)
  df <- data.frame(Accession)
  df <- join(df, gn, by="Accession", match="all")
  adata["geneName"] <- df["GeneName"]
  return(adata)
}

##########  Genie Score Main Function  ##########

SurfaceGenie <- function(adata, processing_opts, groupmethod, numgroups, groupcols, markersample) {
  accessions <- laply(laply(adata["Accession"], as.character), split_acc_iso)
  nsamps <- ncol(adata) - 1
  reqcols <- colnames(adata)
  # Sample grouping
  if("grouping" %in% processing_opts & numgroups > 1){
    adata <- group_samples(adata, groupmethod, groupcols)
    nsamps <- numgroups
    reqcols <- colnames(adata)
  }
  
  # Get proteins where SPC score > 0
  
  adata <- get_SPC(adata, accessions)
  
  accessions <- laply(laply(adata["Accession"], as.character), split_acc_iso)  
  
  #  # probably move this into a loop for calc all four scores
  #  if(!("SPC" %in% processing_opts)){
  #    adata["SPC"] <- matrix(rep(1, nrow(adata)))
  #  }
  
  # Caluclate Gini coefficient and Signal Strength
  # Calculate all four scores
  for(irow in 1:nrow(adata)){
    sdata <- adata[irow, 2:(nsamps + 1)]
    adata[irow, "Gini"] <- get_Gini_coeff(sdata, nsamps)
    adata[irow, "SS"] <- get_signal_strength(sdata)
    adata[irow, "GS"] <- get_dissimilar_score(adata[irow,c("SPC","noSPC","Gini","SS")], nsamps, "SPC")
    adata[irow, "eineG"] <-get_similar_score(adata[irow,c("SPC","noSPC","Gini","SS")], nsamps, "SPC")
    adata[irow, "iGenie"] <-get_dissimilar_score(adata[irow,c("SPC","noSPC","Gini","SS")], nsamps,"noSPC")
    adata[irow, "eineGi"] <-get_similar_score(adata[irow,c("SPC","noSPC","Gini","SS")], nsamps,"noSPC")
  }
  # Return data with SPC score and GS  only
  return(adata)
  #return(adata[,c(reqcols, "SPC", "Gini", "SS", "GS")])
}

##########  SurfaceGenie Export  ##########

SG_export <- function(adata, exportvars, scoringvars) {
  print(colnames(adata))
  print(exportvars)
  print(scoringvars)
  accessions <- laply(laply(adata["Accession"], as.character), split_acc_iso)
  reqcols <- colnames(adata)[1:(ncol(adata)-8)]
  
  # Exclude HLA molecules
  if("HLA" %in% exportvars){
    adata <- filter_by_HLA(adata, accessions)
    accessions <- laply(laply(adata["Accession"], as.character), split_acc_iso)
    #need to remove HLA from exportvars because the other vars concern the number
    #of columns whereas HLA is the rows.  At the end of this function, the exportvars
    #gets passed and can't have any variables that concern rows.
    ev<-c("HLA")
    exportvars <- exportvars[!exportvars %in% ev]
  }
  
  # Export option: append uniprot linkout column
  if("UniProt Linkout" %in% exportvars){
    adata <- append_UPL(adata, accessions)
  }
  
  # Export option: append CD molecule info
  if("CD" %in% exportvars){
    adata <- get_CD(adata, accessions)
  }
  
  # Export option: append # CSPA experiments
  if("CSPA #e" %in% exportvars){
    adata <- get_numCSPA(adata, accessions)
  }
  
  # Return data with export options as well as dataframe size
  return(adata[,c(reqcols, exportvars, scoringvars)])
}

##########  SurfaceGenie Plots  ##########

SPC_hist <- function(adata) {
  scores <- adata[["SPC"]]
  bins <- seq(0, 4, length.out=5)
  hist(scores, breaks=bins, xlab="SPC Score", main="SPC Score Histogram",
       col="#3498db", border="white")
}

SG_dist <- function(adata) {
  accessions <- laply(laply(adata["Accession"], as.character), split_acc_iso)
  adata <- get_CD(adata, accessions)
  adata <- get_geneName(adata, accessions)
  adata <- adata[,c("Accession", "geneName", "GS", "CD")]
  adata <- adata[order(-adata$GS),]
  CD <- adata[,"CD"]
  df = data.frame(CD)
  df$CD<-as.character(df$CD)
  df$CD[!is.na(df$CD)]<-"CD"
  df$CD[is.na(df$CD)]<-"non-CD"
  adata["isCD"]<-df["CD"]
  fa<-list(family="Arial, sans-serif", size=12)
  ft<-list(family="Arial, sans-serif", size=14, color='black')
  plot_ly(data=adata,
          x=~1:nrow(adata),
          y=~GS, 
          type = 'scatter', 
          mode='markers', 
          hoverinfo = 'text', 
          hoverlabel = list(bgcolor='white'),
          text=paste("Gene Name: ", adata$geneName, "<br>Accession: ", adata$Accession, "<br>CD: ", adata$CD), 
          color=~isCD,
          colors=c("#3498db", "#c9c9d4") # blue, grey
  ) %>%
    layout(
      title="<b>Genie Scores in Descending Order</b>", titlefont=ft,
      xaxis=list(title="rank", titlefont=fa, showgrid=FALSE),
      yaxis=list(title="Genie Score", titlefont=fa, showgrid=FALSE),
      legend=list(x=0.7,y=0.9) # controls the location on the plot of the legend
    )
}

SG_dist_export <- function(adata) {
  accessions <- laply(laply(adata["Accession"], as.character), split_acc_iso)
  adata <- get_CD(adata, accessions)
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
  adata <- get_CD(adata, accessions)
  adata <- get_geneName(adata, accessions)
  adata <- adata[,c("Accession", "geneName", "eineG", "CD")]
  adata <- adata[order(-adata$eineG),]
  CD <- adata[,"CD"]
  df = data.frame(CD)
  df$CD<-as.character(df$CD)
  df$CD[!is.na(df$CD)]<-"CD"
  df$CD[is.na(df$CD)]<-"non-CD"
  adata["isCD"]<-df["CD"]
  fa<-list(family="Arial, sans-serif", size=12)
  ft<-list(family="Arial, sans-serif", size=14, color='black')
  plot_ly(data=adata,
          x=~1:nrow(adata),
          y=~eineG, 
          type = 'scatter', 
          mode='markers', 
          hoverinfo = 'text', 
          hoverlabel = list(bgcolor='white'),
          text=paste("Gene Name: ", adata$geneName, "<br>Accession: ", adata$Accession, "<br>CD: ", adata$CD), 
          color=~isCD,
          colors=c("#3498db", "#c9c9d4") # blue, grey
  ) %>%
    layout(
      title="<b>eineG Scores in Descending Order</b>", titlefont=ft,
      xaxis=list(title="rank", titlefont=fa, showgrid=FALSE),
      yaxis=list(title="eineG Score", titlefont=fa, showgrid=FALSE),
      legend=list(x=0.7,y=0.9) # controls the location on the plot of the legend
    )
}

eineG_dist_export <- function(adata) {
  accessions <- laply(laply(adata["Accession"], as.character), split_acc_iso)
  adata <- get_CD(adata, accessions)
  adata <- adata[,c("eineG", "CD")]
  adata <- adata[order(-adata$eineG),]
  CD <- adata[,"eineG"]
  CD[is.na(adata[,"CD"])] <- NA
  plot(1:nrow(adata), adata[,"eineG"], xlab="rank", ylab="eineG Score", 
       main="eineG Scores in Descending Order",
       col="#C0C0C0")
  points(1:nrow(adata), CD, pch=16, col="#3498db")
  legend("topright", legend="CD molecules", col="#3498db", pch=16,
         inset=0.02, box.lwd=0)
}

iGenie_dist <- function(adata) {
  accessions <- laply(laply(adata["Accession"], as.character), split_acc_iso)
  adata <- get_CD(adata, accessions)
  adata <- get_geneName(adata, accessions)
  adata <- adata[,c("Accession", "geneName", "iGenie", "CD")]
  adata <- adata[order(-adata$iGenie),]
  CD <- adata[,"CD"]
  df = data.frame(CD)
  df$CD<-as.character(df$CD)
  df$CD[!is.na(df$CD)]<-"CD"
  df$CD[is.na(df$CD)]<-"non-CD"
  adata["isCD"]<-df["CD"]
  fa<-list(family="Arial, sans-serif", size=12)
  ft<-list(family="Arial, sans-serif", size=14, color='black')
  plot_ly(data=adata,
          x=~1:nrow(adata),
          y=~iGenie, 
          type = 'scatter', 
          mode='markers', 
          hoverinfo = 'text', 
          hoverlabel = list(bgcolor='white'),
          text=paste("Gene Name: ", adata$geneName, "<br>Accession: ", adata$Accession, "<br>CD: ", adata$CD), 
          color=~isCD,
          colors=c("#3498db", "#c9c9d4") # blue, grey
  ) %>%
    layout(
      title="<b>iGenie Scores in Descending Order</b>", titlefont=ft,
      xaxis=list(title="rank", titlefont=fa, showgrid=FALSE),
      yaxis=list(title="iGenie Score", titlefont=fa, showgrid=FALSE),
      legend=list(x=0.7,y=0.9) # controls the location on the plot of the legend
    )
}

iGenie_dist_export <- function(adata) {
  accessions <- laply(laply(adata["Accession"], as.character), split_acc_iso)
  adata <- get_CD(adata, accessions)
  adata <- adata[,c("iGenie", "CD")]
  adata <- adata[order(-adata$iGenie),]
  CD <- adata[,"iGenie"]
  CD[is.na(adata[,"CD"])] <- NA
  plot(1:nrow(adata), adata[,"iGenie"], xlab="rank", ylab="iGenie Score", 
       main="iGenie Scores in Descending Order",
       col="#C0C0C0")
  points(1:nrow(adata), CD, pch=16, col="#3498db")
  legend("topright", legend="CD molecules", col="#3498db", pch=16,
         inset=0.02, box.lwd=0)
}


eineGi_dist <- function(adata) {
  accessions <- laply(laply(adata["Accession"], as.character), split_acc_iso)
  adata <- get_CD(adata, accessions)
  adata <- get_geneName(adata, accessions)
  adata <- adata[,c("Accession", "geneName", "eineGi", "CD")]
  adata <- adata[order(-adata$eineGi),]
  CD <- adata[,"CD"]
  df = data.frame(CD)
  df$CD<-as.character(df$CD)
  df$CD[!is.na(df$CD)]<-"CD"
  df$CD[is.na(df$CD)]<-"non-CD"
  adata["isCD"]<-df["CD"]
  fa<-list(family="Arial, sans-serif", size=12)
  ft<-list(family="Arial, sans-serif", size=14, color='black')
  plot_ly(data=adata,
          x=~1:nrow(adata),
          y=~eineGi, 
          type = 'scatter', 
          mode='markers', 
          hoverinfo = 'text', 
          hoverlabel = list(bgcolor='white'),
          text=paste("Gene Name: ", adata$geneName, "<br>Accession: ", adata$Accession, "<br>CD: ", adata$CD), 
          color=~isCD,
          colors=c("#3498db", "#c9c9d4") # blue, grey
  ) %>%
    layout(
      title="<b>eineGi Scores in Descending Order</b>", titlefont=ft,
      xaxis=list(title="rank", titlefont=fa, showgrid=FALSE),
      yaxis=list(title="eineGi Score", titlefont=fa, showgrid=FALSE),
      legend=list(x=0.7,y=0.9) # controls the location on the plot of the legend
    )
}

eineGi_dist_export <- function(adata) {
  accessions <- laply(laply(adata["Accession"], as.character), split_acc_iso)
  adata <- get_CD(adata, accessions)
  adata <- adata[,c("eineGi", "CD")]
  adata <- adata[order(-adata$eineGi),]
  CD <- adata[,"eineGi"]
  CD[is.na(adata[,"CD"])] <- NA
  plot(1:nrow(adata), adata[,"eineGi"], xlab="rank", ylab="eineGi Score", 
       main="eineGi Scores in Descending Order",
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
