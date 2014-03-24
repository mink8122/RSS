runRIC <- function(sequenceFile) {
  dir <- dirname(sequenceFile)
  system(paste(getwd(), "/bin/RIC.sh ", sequenceFile, sep="")) #Note: replace RIC.pl with C version.
  
  # Filter RSS.
  cutoff.RSS12 = -45
  cutoff.RSS23 = -55
  
#   system(paste("awk -F\"\t\" '{if($5>=-45) print $0}' ",dir,"/combined_Fasta.fa12.txt > ",dir,"/rss12.cutoff.txt",sep=""))
#   system(paste("awk -F\"\t\" '{if($5>=-65) print $0}' ",dir,"/combined_Fasta.fa23.txt > ",dir,"/rss23.cutoff.txt",sep=""))
  
  ##Correct for gene position.
  
}  