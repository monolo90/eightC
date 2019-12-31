  suppressMessages(library(sevenC))
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicInteractions)
library(GenomicRanges)
library(InteractionSet)
  
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")


No.Overlaping <- function(gi1, gi2){
  
  ovAny = overlapsAny(gi1, gi2, maxgap = 1000)
  no.overlaping <- gi1[!ovAny]
  
  return(no.overlaping)
}

    bs.genome <- BSgenome.Hsapiens.UCSC.hg19
    bigwig <- list.files(path="/home/manuel/Escritorio/TFM_FINAL/BIGWIG/") 
    gi <-get(load("gi/gi_K562_YY1_hichip5kb_08.RData"))
    gi2 <-get(load("gi/gi_HCT116_YY1_hichip5kb_08.RData"))
    gi_JASPAR_K562 <-  No.Overlaping(prepareCisPairs_TF("YY1", TF_evidence="motif"), gi)
    gi_JASPAR_HCT116 <-  No.Overlaping(prepareCisPairs_TF("YY1", TF_evidence="motif"), gi2)

########## GI SAMPLE 20.000  ##########
gi_K562 <-sample(gi[!seqnames(anchors(gi, type="first")) %in% "chrY" & !seqnames(anchors(gi, type="second")) %in% "chrY"] , 20000)
gi_HCT116 <-sample(gi2[!seqnames(anchors(gi2, type="first")) %in% "chrY" & !seqnames(anchors(gi2, type="second")) %in% "chrY"] , 20000)
gi_JASPAR_K562 <- sample(gi_JASPAR_K562[!seqnames(anchors(gi_JASPAR_K562, type="first")) %in% "chrY" & !seqnames(anchors(gi_JASPAR_K562, type="second")) %in% "chrY"] , 20000)
gi_JASPAR_HCT116 <- sample(gi_JASPAR_HCT116[!seqnames(anchors(gi_JASPAR_HCT116, type="first")) %in% "chrY" & !seqnames(anchors(gi_JASPAR_HCT116, type="second")) %in% "chrY"], 20000)

# gi_to_gr_midpoint <- function(gi){
#   # first anchor
#   anchor1  <- anchors(gi, type="first")
#   anchor2 <- anchors(gi, type="second")
#   # add seqlengths, just in case
#   seqlengths(anchor1) <- seqlengths(bs.genome)[names(seqlengths(anchor1))]
#   seqlengths(anchor2) <- seqlengths(bs.genome)[names(seqlengths(anchor2))]
#   
#   # Delete chrY, bigwig problems
#   anchor1.1  <- anchor1[!seqnames(anchor1) %in% "chrY"]
#   anchor2.2 <- anchor2[!seqnames(anchor2) %in% "chrY"]
#   # ORDER FOR JASPAR....
#   
#   anchor.start = anchor1.1[start(gr.mid(anchor1.1))  <= end(gr.mid(anchor2.2))]
#   anchor.start2 = anchor2.2[start(gr.mid(anchor1.1)) >= end(gr.mid(anchor2.2))]
#   anchor.end = anchor2.2[start(gr.mid(anchor1.1))  <= end(gr.mid(anchor2.2))]
#   anchor.end2 = anchor1.1[start(gr.mid(anchor1.1)) >= end(gr.mid(anchor2.2))]
#   
#   start = append(anchor.start, anchor.start2)
#   end =  append(anchor.end, anchor.end2)
#   
#   
#   seq = as.character(seqnames(start))
#   strand = as.character(c(rep("*", length(anchor1.1))))
#   
#   gr <- GRanges(seqnames=seq,
#                 IRanges(start = start(gr.mid(start)), 
#                         end=end(gr.mid(end))),
#                     strands= strand)
#   return(gr)
# }
# 
# library(gUtils)
# gi = gi_K562
# bigwig = bigwig.K562
# gi_K562<- gi_to_gr_midpoint(gi_K562)
# gi_JASPAR <-  gi_to_gr_midpoint(gi_JASPAR)

na.bigwig <- function(bigwig){
  for (i in 1:length(bigwig$target)) {
    if (is.na(bigwig$target[i])) {
      bigwig$target[i]=bigwig$assay[i]
    }
  }
  return(bigwig)
}


bigwig.HCT116 <- na.bigwig(subset(df_total, biosample_name %in% "HCT116"))
bigwig.K562 <- na.bigwig(subset(df_total, biosample_name %in% "K562"))

file.rename( from =  paste(bigwig.K562$file_accession,".bigWig",sep = ""), to = paste("K562/",na.rename(bigwig.K562),".bigWig",sep = ""))
file.rename( from =  paste(bigwig.HCT116$file_accession,".bigWig",sep = ""), to = paste("HCT116/",na.rename(bigwig.HCT116),".bigWig",sep = ""))

file.rename( from =   paste("K562/",na.rename(bigwig.K562),".bigWig",sep = ""), to =paste(bigwig.K562$file_accession,".bigWig",sep = ""))
file.rename( from =   paste("HCT116/",na.rename(bigwig.HCT116),".bigWig",sep = ""), to =paste(bigwig.HCT116$file_accession,".bigWig",sep = ""))



na.rename <- function(bigwig){
  for (i in 1:length(bigwig$target)) {
    if (is.na(bigwig$target[i])) {
      bigwig$target[i]=bigwig$assay[i]
    }
  }
  return(bigwig$target)
}


gi_bigwig.value <- function(gi, bigwig){
# first anchor
anchor1  <-  anchors(gi, type="first")
anchor2 <- anchors(gi, type="second")
# add seqlengths, just in case
seqlengths(anchor1) <- seqlengths(bs.genome)[names(seqlengths(anchor1))]
seqlengths(anchor2) <- seqlengths(bs.genome)[names(seqlengths(anchor2))]

# chr Y
anchor1.1  <- anchor1[!seqnames(anchor1) %in% "chrY"]
anchor2.2 <- anchor2[!seqnames(anchor2) %in% "chrY"]
#trycatch <- function(x){ addCovToGR (anchor1, x, window = 500, colname = "reads_500_size")}
#mclapply(bigwig,trycatch , mc.cores = 2)
my.dataframe = data.frame(anchor= "anchor_1",position=(1:length(anchor1.1)))
my.dataframe2 <- data.frame(anchor="anchor_2",position=(1:length(anchor2.2)))
  for(i in 1:nrow(bigwig.K562)){
  #setwd("~/Escritorio/TFM_FINAL/BIGWIG")
  GR_cov <- addCovToGR (anchor1.1, paste(bigwig[i,2],".bigWig",sep = ""), window = 5000, colname = "reads_500_size")
  GR_cov2 <- addCovToGR (anchor2.2, paste(bigwig[i,2],".bigWig",sep = ""), window = 5000, colname = "reads_500_size")
  GR_cov$ML_value <- mean(GR_cov$reads_500_size)
  GR_cov2$ML_value <- mean(GR_cov2$reads_500_size)
  VECTOR = unlist(GR_cov$ML_value, use.names=FALSE)
  VECTOR2 = unlist(GR_cov2$ML_value, use.names=FALSE)
  my.dataframe <- cbind(my.dataframe,  VECTOR)
  my.dataframe2 <- cbind(my.dataframe2,  VECTOR2)
    gc()
  } 
for (i in 1:length(bigwig$target)) {
    if (is.na(bigwig$target[i])) {
        bigwig$target[i]=bigwig$assay[i]
    }
  }
colnames(my.dataframe) <-c("ANCHOR", "POSITION", bigwig$target)
colnames(my.dataframe2) <-c("ANCHOR.2", "POSITION.2", paste(bigwig$target,".2",sep = ""))

complete.dataframe =cbind(my.dataframe,my.dataframe2)
return(complete.dataframe)
}

### K562 #####

YY1_LOOP_K562 = gi_bigwig.value(gi_K562, bigwig.K562)
JASPAR_NEGATIVE_K562 = gi_bigwig.value(gi_JASPAR_K562, bigwig.K562)


YY1_LOOP_K562$LOOP = factor("YES")
JASPAR_NEGATIVE_K562$LOOP = factor("NO")
YY1_K562_MATRIX = rbind(YY1_LOOP_K562, JASPAR_NEGATIVE_K562)
YY1_K562_MATRIX[ ,c("ANCHOR", "POSITION","ANCHOR.2", "POSITION.2")] <- list(NULL)
  
##### HCT116 #######

YY1_LOOP_HCT116 = gi_bigwig.value(gi_HCT116, bigwig.HCT116)
JASPAR_NEGATIVE_HCT116 = gi_bigwig.value(gi_JASPAR_HCT116, bigwig.HCT116)

YY1_LOOP_HCT116$LOOP = factor("YES")
JASPAR_NEGATIVE_HCT116$LOOP = factor("NO")
YY1_HTC116_MATRIX = rbind(YY1_LOOP_HCT116, JASPAR_NEGATIVE_HCT116)
YY1_HTC116_MATRIX[ ,c("ANCHOR", "POSITION","ANCHOR.2", "POSITION.2")] <- list(NULL)

write.csv(file="YY1_K562_MATRIX.csv", x=YY1_K562_MATRIX, row.names = FALSE)
write.csv(file="YY1_HTC116_MATRIX.csv", x=YY1_HTC116_MATRIX, row.names = FALSE)

# http://seqplots.ga/
#     
#   
#   
#   no HACER CLUSTER , TENEMOS QUE HACER ORDENADO
# 
# HACER UN .BED CON EL CENTRO DEL ANCLA 1 Y EL CENTRO DEL ANCLA 2/////////////// MIDDLE !!!! calculateDistances(gi, method="midpoint")  
# 
# mid(ANCLA1)---- mid(ANCLA2)
# bed con LAS DOS LINEAS CELULARES Y 46 - QUE HAY DOS Mask
# 
# SEQPLOT CARGAR LOS 2 BED Y 46 
#   
# 
# Batch operations
# 
# GENERAR UN CONJUNTO PARA CADA LINEA CELULAR UTILIZANDO BED Y JASPAR SOLO YY1, HTC116 Y K562....
# BUSCAR EL PORCENTAJE DE LOS QUE COINCIDEN 2 Y 1 DE GI InteractionSet(
    

# install.packages('devtools')
# install.packages('testthat')
# 

# install.packages("remotes")
# Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE)
#     
# remotes::install_github("mskilab/gUtils")


# 1. SEVENC PARA HACER YYY1 ENCODE Y JASPAR HTC116 Y K562
# OVERLAPS CON LOS LOOP QUE SON REALMENTE EXPERIMENTALES 

lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
library(gUtils)

seq = as.character(seqnames(anchor1.1))
strand = as.character(strand(anchor1.1))

gr <- GRanges(seqnames=seq,
      IRanges(start = start(gr.mid(anchor1)), 
      end=end(gr.mid(anchor2))))
df <- data.frame(seqnames=seq,
                 starts=start(gr.mid(anchor1.1)),
                 ends=end(gr.mid(anchor2.2)),
                 names=c(rep(".", length(anchor1.1))),
                 scores=c(rep(".", length(anchor1.1))),
                 strands=strand)

write.table(df, file = "HCT116_YY1_hichip5kb.bed", quote = F, row.names = F, col.names = F, sep = "\t")

seqlengths(gr) <- seqlengths(bs.genome)[names(seqlengths(gr))]


BiocManager::install("seqplots")





library(seqplots)


export.bed(gr, score = NULL, name="HCT116_YY1_hichip5kb", splitByChrom = FALSE)  

BiocManager::install("RIPSeeker")

library(RIPSeeker)
exportGRanges(gr, format="bed")
exportGRanges("example.txt", format="bed")

df <- data.frame(seqnames=seq,
                 starts=(if (start(gr.mid(anchor1.1)>end(gr.mid(anchor2.2)))){
                   start(gr.mid(anchor1.1))}
                   else {end(gr.mid(anchor2.2))}),
                 ends=(if (start(gr.mid(anchor1.1)>end(gr.mid(anchor2.2)))){
                   end(gr.mid(anchor2.2))}
                   else {start(gr.mid(anchor1.1))}),
                 strands=strand)

# 
# 

  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
    BiocManager::install("seqplots")
    
    bigwig.K562[,2]    
    
    #PRUEBA
    
GR_cov <- addCovToGR (anchor1.1[2], "ENCFF456JVK.bigWig", window = 500, colname = "reads_500_size")

mean(GR_cov$reads_500_size)

GR_cov <- addCovToGR (anchor1.1, paste(bigwig.K562[2,2],".bigWig",sep = ""), window = 500, colname = "reads_500_size")


anchor1  <-  anchors(gi_K562, type="first")
anchor2 <- anchors(gi_K562, type="second")