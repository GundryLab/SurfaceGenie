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

# get_SPC <- function(adata, Accession, species) {
#   if(species=="human"){
#     SPC_scores <- read.csv(file="ref/SPC.csv", header=TRUE)
#   } else if(species=="rat"){
#     SPC_scores <- read.csv(file="ref/Rat_SPC.csv", header=TRUE)
#   } else if(species=="mouse") {
#     SPC_scores <- read.csv(file="ref/Mouse_SPC.csv", header=TRUE)
#   }
#   noiso <- data.frame(Accession)
#   noiso_SPC <- join(noiso, SPC_scores, by="Accession", type="left", match="first")
#   noiso_SPC["SPC"][is.na(noiso_SPC["SPC"])]<-0
#   adata["SPC"] <- noiso_SPC["SPC"]
#   adata["noSPC"] <- matrix(rep(1, nrow(adata)))  
#   return(adata)
# }

get_Gini_coeff <- function(sdata, nsamps) {
  cmat <- matrix(rep(t(sdata), nsamps), ncol=nsamps, byrow=TRUE)
  sumdif <- sum(abs(cmat - t(cmat)))
  return(sumdif/(2*nsamps*sum(sdata)))
}

get_signal_strength <- function(sdata) {
  return(log10(max(sdata)+1))
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
#  adata["UniProt Linkout"] <- laply(Accession, function(x) { paste(baselink, x, sep="") })
  adata["UniProt Linkout"] <- laply(adata["noiso"], function(x) { paste(baselink, x, sep="") })
  return(adata)
}



##########  Genie Score Main Function  ##########

SurfaceGenie <- function(adata, processing_opts, groupmethod, numgroups, groupcols, anno, updateProgress) {
  nsamps <- ncol(adata) - 1
  reqcols <- colnames(adata)
  adata["noiso"] <- laply(laply(adata["Accession"], as.character), split_acc_iso)
  final<-join(adata, anno, by="noiso", type="left", match="first")
  # Sample grouping
  if("grouping" %in% processing_opts & numgroups > 1){
    adata <- group_samples(adata, groupmethod, groupcols)
    nsamps <- numgroups
    reqcols <- colnames(adata)
  }

  # rather than pulling apart and inserting into a data frame and 
  # calling functions to calculate the GS, OmniGenie, etc scores, we 
  # will calculate them directly, drop them into a vector then
  # finally add the vectors to the data frame.  This has reduced 
  # the time for a 3700 line file from 14 secs to 3 (including the
  # annotation merge step)
  
  # set up vectors to eventually put in data frame
  Gmax <- 1 - 1/nsamps
  Ginis<-vector(mode="logical",length=0)
  SSs<-vector(mode="logical",length=0)
  GSs<-vector(mode="logical",length=0)
  IsoGenies<-vector(mode="logical",length=0)
  OmniGenies<-vector(mode="logical",length=0)
  IsoOmniGenies<-vector(mode="logical",length=0)
  
  # loop through data set
  # update the progress meter
  rows = nrow(final)
  for(irow in 1:rows){
    if(irow%%100==0){
      if (is.function(updateProgress)) {
        text <- paste0("row:", irow, "of ", rows)
        updateProgress(detail = text)
      }
    }
    
    #calculate the scores
    sdata <- final[irow, 3:(nsamps + 2)]
    spc <- final[irow,"SPC"]
    Gini <- get_Gini_coeff(sdata, nsamps)
    SS <- get_signal_strength(sdata)
    OmniGenie <-  (Gini/Gmax)^2 * SS
    IsoOmniGenie <- (1-(Gini/Gmax)^2) * SS
    GS <- OmniGenie * spc
    IsoGenie <-IsoOmniGenie * spc
    
    # append the scores to each vector
    Ginis<-append(Ginis, Gini)
    SSs<-append(SSs, SS)
    OmniGenies<-append(OmniGenies, OmniGenie)
    IsoOmniGenies<-append(IsoOmniGenies, IsoOmniGenie)
    GSs<-append(GSs, GS)
    IsoGenies<-append(IsoGenies, IsoGenie)
  }
  
  # put the vectors into the data frame
  final["Gini"]<-Ginis
  final["SS"]<-SSs
  final["OmniGenie"]<-OmniGenies
  final["IsoOmniGenie"]<-IsoOmniGenies
  final["GS"]<-GSs
  final["IsoGenie"]<-IsoGenies
  return(final)
}

##########  SurfaceGenie Export  ##########
# Most everything in here has been removed because rather than adding each 
# annotation as it was clicked, I just put it all in one annotation file and
# display it or not as it is clicked

SG_export <- function(adata, exportvars1, exportvars2 , scoringvars) {
  reqcols <- colnames(adata)[2:(ncol(adata)-13)]

  # Export option: append uniprot linkout column
  if("UniProt Linkout" %in% exportvars2){
    adata <- append_UPL(adata, accessions)
  }

  # Return data with export options as well as dataframe size
  return(adata[,c(reqcols, scoringvars, exportvars1, exportvars2)])
}

##########  SurfaceGenie Plots  ##########


SPC_hist <- function(adata) {
  scores <- adata[["SPC"]]
  scores[is.na(scores)]<-"NA"
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

IsoGenie_dist <- function(adata) {
  accessions <- laply(laply(adata["Accession"], as.character), split_acc_iso)
  adata <- adata[,c("Accession", "geneName", "IsoGenie", "CD")]
  adata <- adata[order(-adata$IsoGenie),]
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
          y=~IsoGenie, 
          type = 'scatter', 
          mode='markers', 
          hoverinfo = 'text', 
          hoverlabel = list(bgcolor='white'),
          text=paste("Gene Name: ", adata$geneName, "<br>Accession: ", adata$Accession, "<br>CD: ", adata$CD, "<br>IsoGenieScore: ", adata$IsoGenie, "<br>Rank: ", 1:nrow(adata)), 
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

IsoGenie_dist_export <- function(adata) {
  accessions <- laply(laply(adata["Accession"], as.character), split_acc_iso)
  adata <- adata[,c("IsoGenie", "CD")]
  adata <- adata[order(-adata$IsoGenie),]
  CD <- adata[,"IsoGenie"]
  CD[is.na(adata[,"CD"])] <- NA
  plot(1:nrow(adata), adata[,"IsoGenie"], xlab="rank", ylab="IsoGenie Score", 
       main="IsoGenie Scores in Descending Order",
       col="#C0C0C0")
  points(1:nrow(adata), CD, pch=16, col="#3498db")
  legend("topright", legend="CD molecules", col="#3498db", pch=16,
         inset=0.02, box.lwd=0)
}

OmniGenie_dist <- function(adata) {
  accessions <- laply(laply(adata["Accession"], as.character), split_acc_iso)
  adata <- adata[,c("Accession", "geneName", "OmniGenie", "CD")]
  adata <- adata[order(-adata$OmniGenie),]
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
          y=~OmniGenie, 
          type = 'scatter', 
          mode='markers', 
          hoverinfo = 'text', 
          hoverlabel = list(bgcolor='white'),
          text=paste("Gene Name: ", adata$geneName, "<br>Accession: ", adata$Accession, "<br>CD: ", adata$CD, "<br>OmniGenieScore: ", adata$OmniGenie, "<br>Rank: ", 1:nrow(adata)), 
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

OmniGenie_dist_export <- function(adata) {
  accessions <- laply(laply(adata["Accession"], as.character), split_acc_iso)
  adata <- adata[,c("OmniGenie", "CD")]
  adata <- adata[order(-adata$OmniGenie),]
  CD <- adata[,"OmniGenie"]
  CD[is.na(adata[,"CD"])] <- NA
  plot(1:nrow(adata), adata[,"OmniGenie"], xlab="rank", ylab="OmniGenie Score", 
       main="OmniGenie Scores in Descending Order",
       col="#C0C0C0")
  points(1:nrow(adata), CD, pch=16, col="#3498db")
  legend("topright", legend="CD molecules", col="#3498db", pch=16,
         inset=0.02, box.lwd=0)
}


IsoOmniGenie_dist <- function(adata) {
  accessions <- laply(laply(adata["Accession"], as.character), split_acc_iso)
  adata <- adata[,c("Accession", "geneName", "IsoOmniGenie", "CD")]
  adata <- adata[order(-adata$IsoOmniGenie),]
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
          y=~IsoOmniGenie, 
          type = 'scatter', 
          mode='markers', 
          hoverinfo = 'text', 
          hoverlabel = list(bgcolor='white'),
          text=paste("Gene Name: ", adata$geneName, "<br>Accession: ", adata$Accession, "<br>CD: ", adata$CD, "<br>IsoOmniGenieScore: ", adata$IsoOmniGenie, "<br>Rank: ", 1:nrow(adata)), 
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

IsoOmniGenie_dist_export <- function(adata) {
  accessions <- laply(laply(adata["Accession"], as.character), split_acc_iso)
  adata <- adata[,c("IsoOmniGenie", "CD")]
  adata <- adata[order(-adata$IsoOmniGenie),]
  CD <- adata[,"IsoOmniGenie"]
  CD[is.na(adata[,"CD"])] <- NA
  plot(1:nrow(adata), adata[,"IsoOmniGenie"], xlab="rank", ylab="IsoOmns Score", 
       main="IsoOmniGenie Scores in Descending Order",
       col="#C0C0C0")
  points(1:nrow(adata), CD, pch=16, col="#3498db")
  legend("topright", legend="CD molecules", col="#3498db", pch=16,
         inset=0.02, box.lwd=0)
}


##########  SPC Lookup  ##########

SPC_lookup <- function(sdata, species) {
  if( species == "Human") {
    if("Accession" %in% colnames(sdata)){
      Accession <- laply(laply(sdata["Accession"], as.character), split_acc_iso)
    }
    else{
      Accession <- laply(laply(sdata[,1], as.character), split_acc_iso)
    }
    SPC_scores <- read.csv(file="ref/SPC_by_Source_sprot.csv", header=TRUE)
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
  } else {
    if("Accession" %in% colnames(sdata)){
      Accession <- laply(laply(sdata["Accession"], as.character), split_acc_iso)
    }
    else{
      Accession <- laply(laply(sdata[,1], as.character), split_acc_iso)
    }
    if(species == "Rat"){
      SPC_scores <- read.csv(file="ref/Rat_SPC.csv", header=TRUE)
    } else if(species == "Mouse"){
      SPC_scores <- read.csv(file="ref/Mouse_SPC.csv", header=TRUE)
    }
    noiso <- data.frame(Accession)
    noiso_SPC <- join(noiso, SPC_scores, by="Accession", match="first")
    sdata["Human_Accession"] <- noiso_SPC["Human_Accession"]
    sdata["SPC"] <- noiso_SPC["SPC"]
#    print(sdata)
    return(sdata)
  }
}

SPC_lookup_for_export <- function(sdata, species) {
  if( species == "Human") {
    if("Accession" %in% colnames(sdata)){
      Accession <- laply(laply(sdata["Accession"], as.character), split_acc_iso)
    }
    else{
      Accession <- laply(laply(sdata[,1], as.character), split_acc_iso)
    }
    SPC_scores <- read.csv(file="ref/SPC_by_Source_sprot.csv", header=TRUE)
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
  } else {
    if("Accession" %in% colnames(sdata)){
      Accession <- laply(laply(sdata["Accession"], as.character), split_acc_iso)
    }
    else{
      Accession <- laply(laply(sdata[,1], as.character), split_acc_iso)
    }
    if(species == "Rat"){
      SPC_scores <- read.csv(file="ref/Rat_SPC.csv", header=TRUE)
    } else if(species == "Mouse"){
      SPC_scores <- read.csv(file="ref/Mouse_SPC.csv", header=TRUE)
    }
    noiso <- data.frame(Accession)
    noiso_SPC <- join(noiso, SPC_scores, by="Accession", match="first")
    sdata["Human_Accession"] <- noiso_SPC["Human_Accession"]
    sdata["SPC"] <- noiso_SPC["SPC"]
    #    print(sdata)
    return(sdata)
  }
}