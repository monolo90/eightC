suppressMessages(library(InteractionSet))


###################################################################################
################################################################################
######################################################################
####################################################################
########################################################################
encode_df <- get_encode_df()
df.filtred <- subset(encode_df ,grepl("long range chromatin interactions",encode_df$output_type))
df.filtred <- subset(df.filtred, grepl("tsv",df.filtred$file_format))
df.filtred <- subset(df.filtred, grepl("ChIA-PET",df.filtred$assay))
df.filtred <- subset(df.filtred, grepl("2_1", df.filtred$technical_replicates))
accesion.download  <- df.filtred$file_accession
target <- df.filtred$target
cell.line <- df.filtred$biosample_name
LOOP.DATA<- data.frame((matrix(ncol = 6, nrow = 0)))
colnames(LOOP.DATA) <- c("Target","Cell_line","Prob_both","Prob_one","Prob_none","WIDTH") 
for(i in 1:length(accesion.download)) {
  tsv.path <-downloadEncode(accesion.download[i])
  # leemos la tabla y separamos cada rango en una tabla distinta
  tab <- read.csv(tsv.path, h = T, sep = "\t")
  t1 <- tab[,1:3]         ;           t2 <- tab[,4:6]
  # vamos a convertirlas a GRanges, y a GRanges no le gustan los nombres de las columnas que vienen dadas en el TSV, asi que cambiamos los nombres
  colnames(t1) <- colnames(t2) <- c("seqnames", "start", "end")
  # ahora sólo hay que pasar las tablas a GRanges y usar ambas como argumentos de la función GInteractions. Lo hacemos de una vez:
  anc1 <- GRanges(t1)
  anc2 <- GRanges(t2)
  gi<- GInteractions(anchor1 = GRanges(t1), anchor2 = GRanges(t2))
  gr <- prepareCisPairs_TF(target[i], cell.line[i])
  var.target <- target[i] ; var.cell.line <- cell.line[i]
WIDTH <- c(100,500,1000,2000) # ATNECION A LA RESOLUCION DE LOS LOOP
for(i in WIDTH) {
  anc1 <- resize(anc1, fix = "center", width = i)
  anc2 <- resize(anc2, fix = "center", width = i)
  #If overlap granges with first expermientas anchors or second anchors looping
  ov <- findOverlaps(gr, anc1, maxgap = 0)  
  ov2 <- findOverlaps(gr, anc2, maxgap = 0)
  #Named with number uniques anchors
  sbj <- unique(subjectHits(ov))
  sbj2 <- unique(subjectHits(ov2))
  #Ginteractions with both anchors, explain...
  which.both <- sbj[sbj %in% sbj2]
  one <- c(sbj[!sbj %in% which.both], sbj2[!sbj2 %in% which.both])
  both <- 100*length(which.both)/length(gi)
  one <-  100*length(only.one)/length(gi)
  none <-  100-(prob.one+prob.both)
  #par(mfrow=c(1,2))    # set the plotting area into a 1*2 array
  #pie(c(prob.one, prob.both, prob.none), c(prob.both, prob.one, prob.none),col = c("green", "yellow", "red"), xlab = c("WIDTH",i))
  #boxplot(list(prob.both, prob.one, prob.none), names = c("Both", "One", "None"), ylab = "Loop probability")
  LOOP.DATA[nrow(LOOP.DATA) + 1,] = c(var.target,var.cell.line,prob.both,prob.none,prob.one,i)
}
var.target <-rm()
var.cell.line <-rm()
}
saveRDS(LOOP.DATA, file="Loop_TF.rds")
#########################################
############################################################
##############################################
width(gi)