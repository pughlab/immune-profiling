library(ggalluvial)
library(randomcoloR)


clontrack.fx <- function(datapath, plotpath, chain, patient_id, countfrac){
  
  if (!(countfrac %in% c("cloneFraction", "cloneCount"))) {
    stop("Error: unknown argument ", countfrac, ". Please provide either cloneFraction or cloneCount.")
    
  }  
  if (!(chain %in% c("TRA", "TRB", "TRD", "TRG"))) {
    stop("Error: unknown argument ", chain, ". Please provide one of the following: TRA, TRB, TRD, TRG.")
    
  }   
  
  flelst <- list.files(datapath, recursive = TRUE,
                       pattern = paste("CLONES", chain, sep = "_"))
  
  # subset to include only downsampled files
  ds_flelst <- flelst[grep("200000", flelst)]
  # subset to patient_id: in CHP_XXX format
  ds_flelst_pt <- ds_flelst[grepl(patient_id, ds_flelst)]
  
  message("list of available files for patient: ", patient_id)
  print(ds_flelst_pt)
  
  #Compile a big file with patient's mixcr files loaded in
  i <- 1
  for (f in ds_flelst_pt){
    mixcrfle <- read.table(paste(datapath, f, sep = ""), 
                           header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE,
                           na.strings = c("", "NA"))
    if(i == 1){
      compldfle <- mixcrfle[!duplicated(mixcrfle$aaSeqCDR3),]
      compldfle <- cbind(cloneno = row.names(compldfle), 
                         filename = f, 
                         compldfle)
      i <- i + 1   
    }
    else{
      compldfle1 <- mixcrfle[!duplicated(mixcrfle$aaSeqCDR3),]
      compldfle1 <- cbind(cloneno = row.names(compldfle1), 
                          filename = f, 
                          compldfle1)
      compldfle <- rbind(compldfle, compldfle1)
      rm(compldfle1)
    }
  }
  
  #If filenames require cleanup, uncomment the following and modify
  ## compldfle$samplename <- gsub(paste(".*",chain, sep = ""), "", compldfle$filename)
  ## compldfle$samplename <- gsub("-PBMC-DNA_2000000.txt", "", compldfle$samplename)  
  # Subset
  CDR3_fraction <- compldfle[, c("samplename","aaSeqCDR3","cloneFraction", "cloneCount")]
  # Number of samples
  mysamples <- unique(CDR3_fraction$samplename)
  #Assign colors to recurring clonotypes
  colorthem <- unique(CDR3_fraction$aaSeqCDR3[duplicated(CDR3_fraction$aaSeqCDR3)])
  norcolorthem <- CDR3_fraction$aaSeqCDR3[!CDR3_fraction$aaSeqCDR3 %in% colorthem]
  
  myColors <- distinctColorPalette(length(colorthem))
  myColors <- c(myColors, rep("white",length(norcolorthem)))
  names(myColors) <- c(colorthem, norcolorthem)
  
  message("the following are the clonotypes we color: ")  
  print(myColors[myColors != "white"])
  
  # Generate a row for each sample that doesnot have recurring clonotype
  ## This ensures colored alluvia are plotted
  
  for(c in colorthem){
    tmp <- CDR3_fraction[CDR3_fraction$aaSeqCDR3 == c,]
    nonexsiting <- mysamples[!mysamples %in% tmp$samplename]
    if(length(nonexsiting) > 0){
      newentries <- data.frame("samplename" = nonexsiting, "aaSeqCDR3" = c, 
                               "cloneFraction" = 0, "cloneCount" = 0)
      CDR3_fraction <- rbind(CDR3_fraction, newentries)
    }
  }
  
  p <-  ggplot(CDR3_fraction, aes(x = samplename, 
                                  y = eval(as.name(countfrac)),
                                  fill = aaSeqCDR3,
                                  stratum = aaSeqCDR3,
                                  alluvium = aaSeqCDR3,
                                  label = aaSeqCDR3))
  
  myp <- p + geom_alluvium(decreasing = FALSE) + 
    geom_stratum(decreasing = FALSE, stat = "alluvium") + 
    scale_fill_manual(values = myColors) +
    theme(axis.title.y = element_text(size = 50),
          axis.title.x = element_blank(),
          axis.line = element_line(color = "black"),
          axis.text = element_text(size = 50),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA),
          legend.key = element_rect(fill = "white", colour = "white"),
          legend.position = "none",
          plot.margin = unit(c(0.2,0,0,0),"cm")) + 
    labs(y = countfrac) 
  
  
  # return(myp)
  pdf(paste(plotpath, "clonetracking_", patient_id, 
            chain, countfrac, ".pdf", sep = ""),
      width = 15, 
      height = 20,
      useDingbats = FALSE,
      onefile = FALSE)       
  print(myp)  
  dev.off()       
  
}