# Input and Output File Path
dataPath = "/home/mkim8/Data/RSS/R_Project/RSS/mm9"

####
## RIC Score Calculations
####
sequenceData = paste(dataPath,"/sample.fa",sep="")
source("RIC_Calculation.R")
runRIC(sequenceData)


#### 
## RSS Concentration Calculations.
####
mm9.original <- read.csv("mm9.RefSeq.original",sep="\t")
bcell.genes <- read.table("bcell.genes",header=F)
tcell.genes <- read.table("tcell.genes",header=F)
source("RSS.R")
calculate_RSS_Conc(mm9.orginal,bcell.genes,tcell.genes)


####
## Moving Window Calculations.
####

