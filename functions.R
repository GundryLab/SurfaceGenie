# SurfaceGenie_0.1/functions.R
library(plyr)
library(stringr)
library(wordspace)
library(gplots)
library(RColorBrewer) 

##########  Genie Score Sub Functions  ##########

split_acc_iso <- function(protID) {
  return(unlist(strsplit(protID, "[-]"))[1])
}

filter_by_SPC <- function(adata, Accession) {
  SPC_scores <- read.csv(file="ref/SPC.csv", header=TRUE)
  noiso <- data.frame(Accession)
  noiso_SPC <- join(noiso, SPC_scores, by="Accession", match="first")
  adata["SPC"] <- noiso_SPC["SPC"]
  idx <- sapply(adata["SPC"] > 0, isTRUE)
  return(adata[idx,])
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
  return(log10(max(sdata)))
}

get_SG_score <- function(pdata, nsamps) {
  Gmax <- 1 - 1/nsamps
  return((pdata["Gini"]/Gmax)^2 * pdata["SPC"] * pdata["SS"])
}

group_samples <- function(adata, groupmethod, groupcols){
  gtags <- c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5")
  numgroups <- length(groupcols)
  for(i in 1:numgroups){
    cols <- laply(strsplit(groupcols[[i]], ",")[[1]], as.integer)
    if(!(is.na(match("ave", groupmethod)))){
      adata[gtags[i]] = rowMeans(adata[cols])
    }
    if(!(is.na(match("sum", groupmethod)))){
      adata[gtags[i]] = rowSums(adata[cols])
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
  
  # Exclude HLA molecules
  if("HLA" %in% processing_opts){
    adata <- filter_by_HLA(adata, accessions)
    accessions <- laply(laply(adata["Accession"], as.character), split_acc_iso)
  }
  
  # Get proteins where SPC score > 0
  if("SPC" %in% processing_opts){
    adata <- filter_by_SPC(adata, accessions)
    accessions <- laply(laply(adata["Accession"], as.character), split_acc_iso)
  }
  else{
    adata["SPC"] <- matrix(rep(1, nrow(adata)))
  }
  
  # Exclude proteins where value = 0 in sample selected for markers
  if(!(is.null(markersample))){
    adata <- adata[(adata[,markersample] > 0),]
    accessions <- laply(laply(adata["Accession"], as.character), split_acc_iso)
  }

  # Caluclate Gini coefficient and Signal Strength
  for(irow in 1:nrow(adata)){
    sdata <- adata[irow, 2:(nsamps + 1)]
    adata[irow, "Gini"] <- get_Gini_coeff(sdata, nsamps)
    adata[irow, "SS"] <- get_signal_strength(sdata)
    adata[irow, "GS"] <- get_SG_score(adata[irow,c("SPC","Gini","SS")], nsamps)
  }
  
  # Return data with SPC score and GS  only
  return(adata[,c(reqcols, "SPC", "Gini", "SS", "GS")])
}

##########  SurfaceGenie Export  ##########

SG_export <- function(adata, exportvars) {
  accessions <- laply(laply(adata["Accession"], as.character), split_acc_iso)
  reqcols <- colnames(adata)[1:(ncol(adata)-4)]
  
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
  return(adata[,c(reqcols, exportvars)])
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
  adata <- adata[,c("GS", "CD")]
  adata <- adata[order(-adata$GS),]
  CD <- adata[,"GS"]
  CD[is.na(adata[,"CD"])] <- NA
  plot(1:nrow(adata), adata[,"GS"], xlab="Proteins", ylab="Genie Score", 
       main="Genie Scores in Descending Order",
       col="#C0C0C0")
  points(1:nrow(adata), CD, pch=16, col="#3498db")
  legend("topright", legend="CD molecules", col="#3498db", pch=16,
         inset=0.02, box.lwd=0)
}

SG_heatmap <- function(adata, hcopts) {
  sdata <- normalize.rows(data.matrix(adata[,2:ncol(adata)]))
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  rowhc <- "row" %in% hcopts
  colhc <- "column" %in% hcopts
  if("both" %in% hcopts){
    rowhc = TRUE
    colhc = TRUE
  }
  lhei <- c(1, 5)
  lwid <- c(1, 3)
  heatmap.2(sdata, Rowv=rowhc, Colv=colhc, dendrogram=hcopts, trace="none", col=hmcol,
            key.title=NA, key.xlab=NA, key.ylab=NA, lhei=lhei, lwid=lwid, labRow=FALSE,
            main="SurfaceGenie: Input Data Heatmap", ylab="Proteins")
}

##########  SPC Lookup  ##########

SPC_lookup <- function(sdata) {
  if("Accession" %in% colnames(sdata)){
    Accession <- laply(laply(sdata["Accession"], as.character), split_acc_iso)
  }
  else{
    Accession <- laply(laply(sdata[,1], as.character), split_acc_iso)
  }
  SPC_scores <- read.csv(file="ref/SPC.csv", header=TRUE)
  noiso <- data.frame(Accession)
  noiso_SPC <- join(noiso, SPC_scores, by="Accession", match="first")
  sdata["SPC"] <- noiso_SPC["SPC"]
  return(sdata)
}