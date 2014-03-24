library("GenomicRanges")

###
# Filter Annotation File downloaded from UCSC Genome Browser Table Tools.
# Required BED Format: name chrom strand txStart txEnd name2 refLink.product
###
filterAnnotation <- function(annotationFile){
  filtered <- annotationFile 
  filtered <- filtered[grep("uncharacterized",filtered[,7],invert=TRUE),] # Remove 'Uncharacterized' Genes.
  filtered <- filtered[grep("hypothetical",filtered[,7],invert=TRUE),] # Remove 'Hypothetical' Genes.
  filtered <- filtered[grep("predicted",filtered[,7],invert=TRUE),] # Remove 'Predicted' Genes.
  
  #Filter out genes in chr_random.
  filtered <- filtered[grep("random",filtered[,2],invert=TRUE),]
  
  # Remove Gene with multiple TSS.
  pos <- filtered[filtered[,3]=="+",]
  neg <- filtered[filtered[,3]=="-",]
  #Correct TxStart and TxEnd for genes on "-" strand. 
  tempCol <- neg$mm9.refGene.txStart
  neg$mm9.refGene.txStart <- neg$mm9.refGene.txEnd
  neg$mm9.refGene.txEnd <- tempCol
  pos <- pos[!(duplicated(pos[c("mm9.refGene.txStart","mm9.refGene.name2")]) | duplicated(pos[c("mm9.refGene.txStart","mm9.refGene.name2")],fromLast=TRUE)),]
  neg <- neg[!(duplicated(neg[c("mm9.refGene.txStart","mm9.refGene.name2")]) | duplicated(neg[c("mm9.refGene.txStart","mm9.refGene.name2")],fromLast=TRUE)),]
  filtered <- rbind(pos,neg)
  
  return(filtered)
}

#Import mm9 RefSeq Annotation. Apply Filter
mm9.original <- read.csv("mm9.RefSeq.original",sep="\t")
mm9.Filtered <- filterAnnotation(mm9.original)

#Import B-cell & T-cell Rag1 binding genes lists sorted by descending RPKM
bcell.genes <- read.table("bcell.genes",header=F)
tcell.genes <- read.table("tcell.genes",header=F)
bcell.genes.filtered <- bcell.genes[which(bcell.genes$V1 %in% mm9.Filtered$mm9.refGene.name2),] #Remove genes not present in mm9.Filtered
tcell.genes.filtered <- tcell.genes[which(tcell.genes$V1 %in% mm9.Filtered$mm9.refGene.name2),] #Remove genes not present in mm9.Filtered

#Import Predicted RSS Data.
mm9.12rss <- read.table("mm9.12rss.filtered.sorted",header=F)
mm9.12rss <- cbind(mm9.12rss,rep("+",nrow(mm9.12rss)))
colnames(mm9.12rss) <-  c("chr","start","end","sequence","score","strand")
mm9.12rss <- with(mm9.12rss, GRanges(chr, IRanges(start, end),strand))

mm9.23rss <- read.table("mm9.23rss.filtered.sorted",header=F)
mm9.23rss <- cbind(mm9.23rss,rep("+",nrow(mm9.23rss)))
colnames(mm9.23rss) <-  c("chr","start","end","sequence","score","strand")
mm9.23rss <- with(mm9.23rss, GRanges(chr, IRanges(start, end),strand))


###
# Calculate RSS Concentration Tables.
###

addToPosTable <- function(group,table) {
  table <- rbind(table,c(i,length(group),sum(width(union(group,group))),length(union(group,group))
                         ,sum(width(union(group,group)))-(length(union(group,group))*27)
                         ,sum(countOverlaps(group,mm9.12rss,minoverlap=27L))
                         ,sum(countOverlaps(group,mm9.12rss,minoverlap=27L))/(sum(width(union(group,group)))-(length(union(group,group))*27))
                         ,sum(width(union(group,group)))-(length(union(group,group))*38)
                         ,sum(countOverlaps(group,mm9.23rss,minoverlap=38L))
                         ,sum(countOverlaps(group,mm9.23rss,minoverlap=38L))/(sum(width(union(group,group)))-(length(union(group,group))*38))))
  return(table)
}

addToNegTable <- function(groupPos,groupNeg,table) {
  overlap <- intersect(groupPos,groupNeg)
  overlapRSS12 <- sum(countOverlaps(overlap,mm9.12rss,minoverlap=27L))
  overlapRSS23 <- sum(countOverlaps(overlap,mm9.23rss,minoverlap=38L))
  table <- rbind(table,c(i,length(groupNeg),sum(width(union(groupNeg,groupNeg))),length(union(groupNeg,groupNeg))
                         ,sum(width(union(groupNeg,groupNeg)))-(length(union(groupNeg,groupNeg))*27)
                         ,sum(countOverlaps(groupNeg,mm9.12rss,minoverlap=27L))
                         ,(sum(countOverlaps(groupNeg,mm9.12rss,minoverlap=27L))-overlapRSS12)/(sum(width(union(groupNeg,groupNeg)))-sum(width(overlap))-(length(union(groupNeg,groupNeg))*27))
                         ,sum(width(union(groupNeg,groupNeg)))-(length(union(groupNeg,groupNeg))*38)
                         ,sum(countOverlaps(groupNeg,mm9.23rss,minoverlap=38L))
                         ,(sum(countOverlaps(groupNeg,mm9.23rss,minoverlap=38L))-overlapRSS23)/(sum(width(union(groupNeg,groupNeg)))-sum(width(overlap))-(length(union(groupNeg,groupNeg))*38))))
  return(table)
}


#Create Ranked BT Groups by combining top (1-200,201-400,401-600...3801-4000) from each b-cell & t-cell list. 
ranked <- seq(200,4000,200)

pos.1kb.Table <- NULL
neg.1kb.Table <- NULL
pos.2kb.Table <- NULL
neg.2kb.Table <- NULL
pos.asym.Table <- NULL
neg.asym.Table <- NULL

for (i in ranked){
  #Create BT Groups.
  b <- bcell.genes.filtered[(i-199):i]
  t <- tcell.genes.filtered[(i-199):i]
  bt <- c(as.matrix(b),as.matrix(t)) #Combine list
  bt <- bt[!duplicated(bt)] #Remove duplicate genes
  print(length(bt))
  
  #Retrieve gene information from mm9 Annotation. Genes that match bt list are classified into Rag1 Positive and the rest are grouped into Rag1 Negative.
  posGroup <- mm9.Filtered[which(mm9.Filtered$mm9.refGene.name2 %in% bt),] 
  posGroup <- posGroup[,c('mm9.refGene.chrom','mm9.refGene.txStart','mm9.refGene.strand','mm9.refGene.name2')]
  colnames(posGroup) <- c("chr","start","strand","id")
  posGroup <- with(posGroup, GRanges(chr, IRanges(start, start), strand, id=id))
  negGroup <- mm9.Filtered[which(!(mm9.Filtered$mm9.refGene.name2 %in% bt)),]
  negGroup <- negGroup[,c('mm9.refGene.chrom','mm9.refGene.txStart','mm9.refGene.strand','mm9.refGene.name2')]
  colnames(negGroup) <- c("chr","start","strand","id")
  negGroup <- with(negGroup, GRanges(chr, IRanges(start, start), strand, id=id))
  
  #1kb
  pos.1kb <- flank(posGroup,1000,both=TRUE)
  pos.1kb <- union(pos.1kb,pos.1kb)
  strand(pos.1kb) <- "+" #Set all strand to '+' for overlap with RSS.
  pos.1kb.Table <- addToPosTable(pos.1kb,pos.1kb.Table)
  
  neg.1kb <- flank(negGroup,1000,both=TRUE)
  neg.1kb <- union(neg.1kb,neg.1kb)
  strand(neg.1kb) <- "+" #Set all strand to '+' for overlap with RSS.
  neg.1kb.Table <- addToNegTable(pos.1kb,neg.1kb,neg.1kb.Table) 
  
  #2kb 
  pos.2kb <- flank(posGroup,2000,both=TRUE)
  pos.2kb <- union(pos.2kb,pos.2kb)
  strand(pos.2kb) <- "+" #Set all strand to '+' for overlap with RSS.
  pos.2kb.Table <- addToPosTable(pos.2kb,pos.2kb.Table)
  
  neg.2kb <- flank(negGroup,2000,both=TRUE)
  neg.2kb <- union(neg.2kb,neg.2kb)
  strand(neg.2kb) <- "+" #Set all strand to '+' for overlap with RSS.
  neg.2kb.Table <- addToNegTable(pos.2kb,neg.2kb,neg.2kb.Table)
  
  #Asymmetric
  pos.asym <- flank(posGroup[strand(posGroup)=="+",],1000,start=TRUE,both=FALSE)
  pos.asym <- flank(pos.asym[strand(pos.asym)=="+",],3000,start=FALSE,both=FALSE)
  pos.asym <- flank(pos.asym[strand(pos.asym)=="-",],3000,start=TRUE,both=FALSE)
  pos.asym <- flank(pos.asym[strand(pos.asym)=="-",],1000,start=FALSE,both=FALSE)
  pos.asym <- union(pos.asym,pos.asym)
  strand(pos.asym) <- "+" #Set all strand to '+' for overlap with RSS.
  pos.asym.Table <- addToPosTable(pos.asym,pos.asym.Table)
  
  neg.asym <- flank(negGroup[strand(negGroup)=="+",],1000,start=TRUE,both=FALSE)
  neg.asym <- flank(neg.asym[strand(neg.asym)=="+",],3000,start=FALSE,both=FALSE)
  neg.asym <- flank(neg.asym[strand(neg.asym)=="-",],3000,start=TRUE,both=FALSE)
  neg.asym <- flank(neg.asym[strand(neg.asym)=="-",],1000,start=FALSE,both=FALSE)
  neg.asym <- union(neg.asym,neg.asym)
  strand(neg.asym) <- "+" #Set all strand to '+' for overlap with RSS.
  neg.asym.Table <- addTonegTable(pos.asym,neg.asym,pos.asym.Table)
}



#Create Cumulative BT Groups by combining top (100,250,500,1000,2000,4000) from each b-cell & t-cell list. 
cumulative <- c(100,250,500,1000,2000,4000)