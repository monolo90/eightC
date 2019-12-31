  prepareCisPairs_TF <- function(TF_name, TF_cell_type, TF_evidence="exp", maxdist=100){
    granges <- getTFpositions(TF_name, TF_cell_type, TF_evidence, maxdist)
    granges$score <- 1
    gc()
    return(prepareCisPairs(sort(granges), maxDist =10^6 ))
  }
  
  
  getTFpositions <- function(TF_name, TF_cell_type=NULL, TF_evidence, maxdist){
    TF_evidence <- tolower(TF_evidence)
    if(!TF_evidence %in% c("exp","motif","pred","custom")) stop("Argument 'TF_evidence' should be either 'exp','motif','pred' or 'custom'")
    if(TF_evidence == "exp")
    {granges <- getTFencode(TF_name,TF_cell_type, maxdist)}
    if(TF_evidence == "motif")
    {granges <- getTFmotif(TF_name)}# no cell type
    if(TF_evidence == "exp_motif")#HACER LAS DOS COMBINADAS
    {granges <- getTFpred(TF_name,TF_cell_type)}
    return(granges)
  }
  
  
  getTF.biosample_selection <- function(df){
  tabla <- as.data.frame(table(df$dataset_biosample_summary))
  if (length(tabla$Var1) > 1){
    
    elements <- as.vector((tabla[,1]))
    message(cat("WARNING HAVE THIS", length(tabla$Var1), "BIOSAMPLE DATASET: ","\n "),cat(paste(1:length(elements),elements, "\n"),fill = FALSE))
    my.variable1 <- as.numeric(readline(prompt="Select a number : "))
    df <- subset(df , dataset_biosample_summary ==tabla[my.variable1,1] )
  }
  return(df)
  }
  
  getTF.bed.selection.hg19.or.grh38 <- function(query_results){
    df.filtred <- subset(query_results ,grepl("[Hh][Gg]19",query_results$assembly))
  if(empty(df.filtred)){
    df.filtred <- subset(query_results ,grepl("[Gg][Rr][Cc][Hh]38",query_results$assembly))
      message("WARNING, use GRCH38")
    }
df.filtred2 <- subset(df.filtred , df.filtred$output_type == "optimal idr thresholded peaks")
  if(empty(df.filtred2)){
    df.filtred2 <-  subset(df.filtred , df.filtred$output_type == "conservative idr thresholded peaks")
    if(empty(df.filtred2)){
    df.filtred2 <- subset(df.filtred , df.filtred$output_type == "peaks and background as input for IDR")
      if(empty(df.filtred2)){
    df.filtred2 <- subset(df.filtred , df.filtred$output_type == "replicated peaks")
        if(empty(df.filtred2)){
    # df.filtred2 <- subset(df.filtred , df.filtred$output_type == "peaks") RESOLVER LO DE LAS REPLICAS.  
      }}}}
  if(unique(df.filtred2$file_type=="bed broadPeak")){
    message("WARNING, use bed broadPeak")
}
return(df.filtred2)
  }
  
  
  download.trycatch <- function(df.filtred) {
    out <- tryCatch(
      {
        message("Establish connection with ENCODE")
        downloadEncode(df.filtred)
      },
      
      error = function(e) {
        message(paste("No internet conection"))
        #message(cond)
        return(NULL)
      }
    )    
    return(out)
  }
  
  getTFencode <- function(TF_name, TF_cell_type, maxdist){
    if ("encode_df" %in% ls(envir = .GlobalEnv)) {
    } else {
      encode_df <- get_encode_df()
      encode_df <- subset(encode_df, encode_df$assembly=="GRCh38" | encode_df$assembly=="hg19")
      encode_df <<- subset(encode_df,  encode_df$file_format == "bed" & encode_df$assay=="ChIP-seq")#FILTRAR POR PICOS CORRECTOS Y VER QUE PONER EN LA SELECCION  AÑADIR CHIP-SEQ AQUI
    }
    query_results <- queryEncode(df=encode_df, biosample_name = TF_cell_type, organism = "Homo sapiens", target = TF_name, file_format = "bed", fixed = TRUE)
    df.filtred <- getTF.bed.selection.hg19.or.grh38(query_results)
    if(empty(df.filtred)){
      message(paste("NO exist",TF_name, "or", TF_cell_type))
    }
    query_results_biosample <- getTF.biosample_selection(df.filtred) 
    filtred.bed <- download.trycatch(query_results_biosample)
    if (is.null(filtred.bed)) {
      
      stop(return(message("NEED DOWNLOAD .BED FROM DATABASE")))
    }
    else {
      granges <- getTF.read.combined.resize(df.filtred, filtred.bed, maxdist)
      unlink(filtred.bed)
      return(granges)
    }
  }
  
  
  getTFmotif <- function(TF_name){
    Tf_jaspar <- paste(getwd(),"/TF_jaspar_corresp.csv",sep = "")
    table  <- read.csv(Tf_jaspar, header = TRUE, sep = ";")
    select_jaspar <- subset(table, table$DNA.binding.protein==TF_name)
    separate_jaspar <- as.vector(strsplit(as.character(select_jaspar$JASPAR.motifs.used),','))
    only.vector <- separate_jaspar[[1]]
    if (length(only.vector) > 1){
      message(cat("WARNING HAVE THIS", length(only.vector), "JASPAR motifs: ","\n "),cat(paste(1:length(only.vector),only.vector,"\n"),fill = FALSE))
      my.variable1 <- as.numeric(readline(prompt="Select a number : "))
      x <- only.vector[my.variable1]
    }
    if(length(only.vector)== 1){
      x <- only.vector  
      
    }
    extract <- str_extract(x, "[A-Z]+\\d+")
    fn <- paste(extract,"_top_30K_pVal_0.01.RData", sep = "")
    dir <- paste(getwd(),"/jaspar_RData/pVal_1e-02_top_30k/",fn,sep = "")  #find.file(fn, dir = dir, dirs = NULL)
    gr <- load(dir)
    return(get(gr))
  }  
  
  
  getTF.read.combined.resize <- function(df.filtred, filtred.bed, maxdist){
    granges.bed <- lapply(filtred.bed, readGeneric) #read one by one .bed genomicranges####
    granges.as.GRANGES <- GRangesList(granges.bed)
    granges.combined <- unlist(granges.as.GRANGES, recursive = TRUE, use.names = TRUE)
    chrom.names <<- paste("chr", c(1:22, "X", "Y", "M"), sep = "")
    granges.combined.filtred <- granges.combined[seqnames(granges.combined) %in% chrom.names]
    
    if (unique(df.filtred$assembly)%in% "hg19"){
     genome(granges.combined.filtred) <- "hg19"
     gr <- granges.combined.filtred
    }
    else {
      gr <- GRCH38.to.hg19(granges.combined.filtred)
    }
    granges.final <- resize(reduce(gr), fix = 'center', width = maxdist)
    return(granges.final)
  }

GRCH38.to.hg19 <- function(gr){
  path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
      ch = import.chain(path)
      seqlevelsStyle(gr) = "UCSC"  # necessary
      gr = liftOver(gr, ch)
      gr = unlist(gr)
      gr = new("GRanges", gr)
      genome(gr) = "hg19"
      return(gr)
}
  
######################################################################################
#############           NEGATIVE GENOMIC INTERACTIONS            ####################
#####################################################################################
  
gi_K562 <- prepareCisPairs_TF("YY1", "K562")
  gi_HCT116 <- prepareCisPairs_TF("YY1", "HCT116")
gi_JASPAR <- prepareCisPairs_TF("YY1", TF_evidence="motif")

gi_exp_K562 <-get(load("gi/gi_K562_YY1_hichip5kb_08.RData"))
gi_exp_HCT116 <-get(load("gi/gi_HCT116_YY1_hichip5kb_08.RData"))


No.Overlaping <- function(gi1, gi2){

ovAny = overlapsAny(gi1, gi2, maxgap = 1000)
no.overlaping <- gi1[!ovAny]
overlaping <- gi1[ovAny]

return(c(length(gi1),length(no.overlaping), length(overlaping)))
}

x = No.Overlaping(gi_HCT116, gi_exp_HCT116)
y = No.Overlaping(gi_K562, gi_exp_K562)
z = No.Overlaping(gi_JASPAR, gi_exp_HCT116)
u =No.Overlaping(gi_JASPAR, gi_exp_K562)


labels <-  c("No solapan","solapamiento")

x.piepercent<- round(c(x[2]/x[1]*100,x[3]/x[1]*100),2)
y.piepercent<- round(c(y[2]/y[1]*100,y[3]/y[1]*100),2)
z.piepercent<- round(c(z[2]/z[1]*100,z[3]/z[1]*100),2)
u.piepercent<- round(c(u[2]/u[1]*100,u[3]/u[1]*100),2)

# Give the chart file a name.
pdf(file = "overlaping.pdf")


# Plot the chart.
par(mfrow=c(2,2))
pie(c(x[2],x[3]), labels = paste(x.piepercent,"%",sep = " "), main = "Solapamiento en HCT116",col = rainbow(length(x.piepercent)))
legend("topright",c("No solapan","Solapamiento"), cex = 0.6,
       fill = rainbow(length(x.piepercent)))
pie(c(y[2],y[3]), labels = paste(y.piepercent,"%",sep = " "), main = "Solapamiento en K562",col = rainbow(length(y.piepercent)))
legend("topright",c("No solapan","Solapamiento"), cex = 0.6,
       fill = rainbow(length(y.piepercent)))
pie(c(z[2],z[3]), labels = paste(z.piepercent,"%",sep = " "), main = "Solapamiento en Jaspar/HCT116",col = rainbow(length(z.piepercent)))
legend("topright",c("No solapan","Solapamiento"), cex = 0.6,
       fill = rainbow(length(z.piepercent)))
pie(c(u[2],u[3]), labels = paste(u.piepercent,"%",sep = " "), main = "Solapamiento en Jaspar/K562",col = rainbow(length(u.piepercent)))
legend("topright",c("No solapan","Solapamiento"), cex = 0.6,
       fill = rainbow(length(u.piepercent)))
# Save the file.
dev.off()


######################################################################################
#############        NEGATIVE GENOMIC INTERACTIONS TO BED        ####################
#####################################################################################

suppressMessages(library(sevenC))
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicInteractions)
library(GenomicRanges)
library(InteractionSet)


No.Overlaping <- function(gi1, gi2){
  
  ovAny = overlapsAny(gi1, gi2, maxgap = 1000)
  no.overlaping <- gi1[!ovAny]
  
  return(no.overlaping)
}

gi_HCT116_Negative <- No.Overlaping(gi_HCT116, gi_exp_HCT116)
gi_K562_Negative <- No.Overlaping(gi_K562, gi_exp_K562)
gi_JASPAR_HCT116_Negative <- No.Overlaping(gi_JASPAR, gi_exp_HCT116)
gi_JASPAR_K562_Negative <- No.Overlaping(gi_JASPAR, gi_exp_K562)

###############################################################
######################## TO BED  ##############################


bs.genome <- BSgenome.Hsapiens.UCSC.hg19
library(gUtils)

  gi_to_bed_midpoint <- function(gi, name){
  # first anchor
  anchor1  <- anchors(gi, type="first")
  anchor2 <- anchors(gi, type="second")
  # add seqlengths, just in case
  seqlengths(anchor1) <- seqlengths(bs.genome)[names(seqlengths(anchor1))]
  seqlengths(anchor2) <- seqlengths(bs.genome)[names(seqlengths(anchor2))]
  
  # Delete chrY, bigwig problems
  anchor1.1  <- anchor1[!seqnames(anchor1) %in% "chrY"]
  anchor2.2 <- anchor2[!seqnames(anchor2) %in% "chrY"]
  # ORDER FOR JASPAR....
  
  anchor.start = anchor1.1[start(gr.mid(anchor1.1))  <= end(gr.mid(anchor2.2))]
  anchor.start2 = anchor2.2[start(gr.mid(anchor1.1)) > end(gr.mid(anchor2.2))]
  anchor.end = anchor2.2[start(gr.mid(anchor1.1))  <= end(gr.mid(anchor2.2))]
  anchor.end2 = anchor1.1[start(gr.mid(anchor1.1)) > end(gr.mid(anchor2.2))]
  
  start = append(anchor.start, anchor.start2)
  end =  append(anchor.end, anchor.end2)
  
  
  seq = as.character(seqnames(start))
  strand = as.character(c(rep("*", length(anchor1.1))))
  
  df <- data.frame(seqnames=seq,
                   starts= start(gr.mid(start)),
                   ends = end(gr.mid(end)) ,
               strands= strand)
              
  
  write.table(df, file = name, quote = F, row.names = F, col.names = F, sep = "\t")
  return(name)
  }

gi_to_bed_midpoint(gi_HCT116_Negative, "BED/gi_HCT116_Negative.bed")
gi_to_bed_midpoint(gi_K562_Negative, "BED/gi_K562_Negative.bed")
gi_to_bed_midpoint(gi_JASPAR_HCT116_Negative, "BED/gi_JASPAR_HCT116_Negative.bed")
gi_to_bed_midpoint(gi_JASPAR_K562_Negative, "BED/gi_JASPAR_K562_Negative.bed")


#sudo apt-get install libcurl4-openssl-dev libxml2-dev libssl-dev 
install.packages("RCurl")
install.packages("XML")
install.packages("curl")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("sevenC")

BiocManager::install(c("GenomicRanges","sevenC","GenomicInteractions","InteractionSet","ENCODExplorer","genomation","gwascat","BiocGenerics","rtracklayer","liftOver"))  
install.packages(c("DataCombine","knitr"))


library(sevenC)
library(GenomicRanges)
library(GenomicInteractions)
library(InteractionSet)
library(ENCODExplorer)  
library(genomation) 
library(DataCombine)
library(plyr)
library(gwascat)
library(BiocGenerics)
library(rtracklayer)
library(liftOver)
library(stringr)


