#Authors: Keith Dunaway and Charles Mordaunt
#setwd("~/Lab/Dup15/SingleCpG/")
setwd(dir = "/Users/kdunaway1/Lab/Dup15/SingleCpG/")
rm(list=ls()) # Clears memory
library(bsseq)
library(DSS)
library(permute)

#################################################
# Global Variables
#################################################
# Run perl Permeth_to_DSSformat.pl to make these files from PerMeth files
DSSprefix = c("DSSformat/JLKD001_",
              "DSSformat/JLKD002_",
              "DSSformat/JLKD003_",
              "DSSformat/JLKD004_",
              "DSSformat/JLKD005_",
              "DSSformat/JLKD014_",
              "DSSformat/JLKD040_",
              "DSSformat/JLKD041_",
              "DSSformat/JLKD042_",
              "DSSformat/JLKD026_",
              "DSSformat/JLKD028_",
              "DSSformat/JLDS017_",
              "DSSformat/JLDS018_",
              "DSSformat/JLDS019_",
              "DSSformat/JLKD006_",
              "DSSformat/JLKD007_",
              "DSSformat/JLKD008_",
              "DSSformat/JLKD021_",
              "DSSformat/JLKD022_",
              "DSSformat/JLKD023_")
numCtrl = 11
numExp = 9
genome = "hg38"
CGI_bedfile = "hg38_genome_CGI.bed"
GenesPlus_bedfile = "hg38_genes_plus.bed"
GenesMinus_bedfile = "hg38_genes_minus.bed"

nperm = 1000
mc.cores = 4
outprefix = "Dup15"
colorCtrl = "#999999"
colorExp = "#56B4E9"


#Pre-loop initialization
numMinCtrl = numCtrl - 1
numMinExp = numExp - 1
CTRLgroup = paste("C",1:numCtrl,sep="")
EXPgroup = paste("E",1:numExp,sep="")
dmrs_passed = NULL
dmrs_gold = NULL

# Create Matrices for Permutations
#idxMatrix <- permuteAll(nperm = nperm, design = numCtrl+numExp)
idxMatrix <- permuteAll(nperm = 10, design = numCtrl+numExp)
idxMatrixsub <- subsetByMatrix(c(CTRLgroup, EXPgroup), idxMatrix)
idxMatrix1 <- rbind(EXPgroup, idxMatrixsub[, 1:numExp])
idxMatrix2 <- rbind(CTRLgroup, idxMatrixsub[, (numExp+1):(numExp+numCtrl)])


#################################################
# Functions written by Kasper Daniel Hansen
# Downloaded from: https://github.com/kasperdanielhansen/bsseq/blob/master/R/permutations.R
# Paper: (Hansen, 2014)
#   "Large-scale hypomethylated blocks associated with Epstein-Barr virus–induced B-cell immortalization"
#################################################
permuteAll <- function(nperm, design) {
  # Creates an idxMatrix of permutations based on design (number of samples), and nperm (number of permutations)
  message(sprintf("[permuteAll] performing %d unrestricted permutations of the design matrix\n", nperm))
  CTRL <- how(nperm = nperm)
  idxMatrix <- shuffleSet(n = design, control = CTRL)
}
subsetByMatrix <- function(vec, mat) {
  # Assigns sample names to numbers in permutations
  apply(mat, 2, function(xx) vec[xx])
}
getNullDistribution_BSmooth.tstat <- function(BSseq, idxMatrix1, idxMatrix2,estimate.var = "same", local.correct = TRUE,cutoff, stat = "tstat.corrected", maxGap = 300, mc.cores = 1) {
  # For each permutation, creates tstat object and finds dmrs, returns a list
  stopifnot(nrow(idxMatrix1) == nrow(idxMatrix2))
  message(sprintf("[getNullDistribution_BSmooth.tstat] performing %d permutations\n", nrow(idxMatrix1)))
  nullDist <- lapply(1:nrow(idxMatrix1), function(ii) {
    ptime1 <- proc.time()
    BS.tstat <- BSmooth.tstat(BSseq, estimate.var = estimate.var,
                              group1 = idxMatrix1[ii,],
                              group2 = idxMatrix2[ii,],
                              local.correct = local.correct, maxGap = 10^8,
                              verbose = FALSE)
    dmrs0 <- dmrFinder(BS.tstat, stat = stat, cutoff = cutoff, maxGap = maxGap)
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    message(sprintf("[getNullDistribution_BSmooth.tstat] completing permutation %d in %.1f sec\n", ii, stime))
    dmrs0
  })
  nullDist
}
subsetDmrs <- function(xx) {
  if(is.null(xx) || is(xx, "try-error"))
    return(NULL)
  out <- xx[ xx[,"n"] >= 3 & abs(xx[, "meanDiff"]) > 0.1 & xx[, "invdensity"] <= 300, ]
  if(nrow(out) == 0)
    return(NULL)
  out
}
getFWER <- function(null) {
  # Calculates FWER from subsetted null DMR list
  reference <- null[[1]]
  null <- null[-1]
  null <- null[!sapply(null, is.null)]
  better <- sapply(1:nrow(reference), function(ii) {
    # meanDiff <- abs(reference$meanDiff[ii])
    areaStat <- abs(reference$areaStat[ii])
    n <- reference$n[ii]
    out <- sapply(null, function(nulldist) {
      # any(abs(nulldist$meanDiff) >= meanDiff &
      #     nulldist$n >= n)
      any(abs(nulldist$areaStat) >= areaStat &
            nulldist$n >= n)
    })
    sum(out)
  })
  better
}

#################################################
# Function written by Dave Tang
# Downloaded from: https://github.com/davetang/bedr/blob/master/R/bed_to_granges.R
#################################################
bed_to_granges <- function(file){
  # Converts bed file to GRanges object
  df <- read.table(file,
                   header=F,
                   stringsAsFactors=F)
  
  if(length(df) > 6){
    df <- df[,-c(7:length(df))]
  }
  
  if(length(df)<3){
    stop("File has less than 3 columns")
  }
  
  header <- c('chr','start','end','id','score','strand')
  names(df) <- header[1:length(names(df))]
  
  if('strand' %in% colnames(df)){
    df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
  }
  
  library("GenomicRanges")
  
  if(length(df)==3){
    gr <- with(df, GRanges(chr, IRanges(start, end)))
  } else if (length(df)==4){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id))
  } else if (length(df)==5){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score))
  } else if (length(df)==6){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
  }
  return(gr)
}

#################################################
# Functions 
#################################################
write.bed = function(df,file,header="",name="name"){
  #Generalized write bed function
  # Example Output:
  #  chr1    10469   10470   0.00-1  0       +       0       0       0,0,0
  #Area stat as name
  #color for direction (hyper = red, hypo = blue)
  colnum = 3
  if(name %in% colnames(df)) {namedat = df[,name];colnum = 4} else{namedat = rep(".",length(df[,1]))}
  if("score" %in% colnames(df)) {score = df$score;colnum = 5} else{score = rep("0",length(df[,1]))}
  if("strand" %in% colnames(df)) {strand = df$strand;colnum = 6} else{strand = rep("+",length(df[,1]))}
  if("thickStart" %in% colnames(df)) {thickStart = df$thickStart;colnum = 7} else{thickStart = rep("0",length(df[,1]))}
  if("thickEnd" %in% colnames(df)) {thickEnd = df$thickEnd;colnum = 8} else{thickEnd = rep("0",length(df[,1]))}
  if("itemRgb" %in% colnames(df)) {itemRgb = df$itemRgb;colnum = 9} else{itemRgb = rep("0,0,0",length(df[,1]))}
  if("blockCount" %in% colnames(df)) {blockCount = df$blockCount;colnum = 10} else{blockCount = rep("0",length(df[,1]))}
  if("blockSizes" %in% colnames(df)) {blockSizes = df$blockSizes;colnum = 11} else{blockSizes = rep("0",length(df[,1]))}
  if("blockStarts" %in% colnames(df)) {blockStarts = df$blockStarts;colnum = 12} else{blockStarts = rep("0",length(df[,1]))}
  
  subdf = cbind(df$chr,df$start,df$end,namedat,score,strand,thickStart,thickEnd,itemRgb,blockCount,blockSizes,blockStarts)
  
  if(header == ""){
    write.table(subdf[,1:colnum], file=file, quote=FALSE, sep='\t', col.names = FALSE, row.names=FALSE)
  }else{
    write(header,file=file)
    write.table(subdf[,1:colnum], file=file, quote=FALSE, sep='\t', col.names = FALSE, row.names=FALSE,append=TRUE)
  }
  cat(paste("Finished writing bedfile",file,"with",colnum,"columns\n"))
}
write.dmrs_bed = function(df,file,trackname,genome){
  #Area stat as name
  #color for direction (hyper = red, hypo = blue)
  df$itemRgb = ifelse(df$direction == "hypo","0,0,255","255,0,0")
  headerline = paste("track name=",trackname," description=",trackname," useScore=0 itemRgb=On genome=",genome,sep="")
  write.bed(df,file,header=headerline,name="areaStat")
}

#################################################
# DMR Finder Pipeline
#################################################

# Load Genes and CGI's
genome_CpG_islands <- bed_to_granges(CGI_bedfile)
genome_genes_plus <- bed_to_granges(GenesPlus_bedfile)
genome_genes_minus <- bed_to_granges(GenesMinus_bedfile)
genome_list <- list(genome_CpG_islands, genome_genes_plus, genome_genes_minus)
names(genome_list) <- c("CGI", "+", "-")

if(genome == "hg38"){chroms = paste("chr",1:22,sep=""); chroms = c(chroms,"chrX","chrY")
}else if(genome == "mm10"){chroms = paste("chr",1:19,sep="");chroms = c(chroms,"chrX","chrY")  
}else if(genome == "rn6"){chroms = paste("chr",1:20,sep="");chroms = c(chroms,"chrX","chrY")  
}else{cat(paste("Warning! Chromosome names unknown because genome is defined as ",genome,"\n"))}

# Loops through every chromosome and processes everything
for(n in 1:length(chroms)){
  chrom = chroms[n]
  cat("Loading",chrom,"...\n")
  DSSlist = list(read.table(paste(DSSprefix[1],chrom,".DSS.txt",sep=""), header=TRUE))
  for(i in 2:length(DSSprefix))
    DSSlist = c(DSSlist,list(read.table(paste(DSSprefix[i],chrom,".DSS.txt",sep=""), header=TRUE)))
  cat("Making BSobject\n")
  BSobj <- makeBSseqData(DSSlist,c(CTRLgroup,EXPgroup))
  BSobj_BSmooth <- BSmooth(BSobj, mc.cores = mc.cores)    
  cat("Completed smoothing\n")
  rm(DSSlist) #To free up memory
  pData <- pData(BSobj_BSmooth)
  pData$type <- as.character(c(rep(c("Ctrl"), numCtrl), rep(c("Exp") ,numExp)))
  pData$col <- c(rep(c(colorCtrl), numCtrl), rep(c(colorExp), numExp))
  pData(BSobj_BSmooth) <- pData
  BS.cov <- getCoverage(BSobj_BSmooth)
  keepLoci.ex <- which(rowSums(BS.cov[, BSobj_BSmooth$type == "Ctrl"] >= 1) >= numMinCtrl &rowSums(BS.cov[, BSobj_BSmooth$type == "Exp"] >= 1) >= numMinExp)
  BSobj_BSmooth_keep <- BSobj_BSmooth[keepLoci.ex,]    
  BSobj_BSmooth_tstat <- BSmooth.tstat(BSobj_BSmooth_keep,group1 = EXPgroup, group2 = CTRLgroup,estimate.var = "same", local.correct = TRUE, qSd = 0.75, k = 101)    
  cutoff <- qt(1-0.05/2, BSobj_BSmooth@colData@nrows-2)
  dmrs_all_BSobj_BSmooth_tstat <- dmrFinder(BSobj_BSmooth_tstat, cutoff = c(-cutoff, cutoff), maxGap = 300)    
  dmrs_select_BSobj_BSmooth_tstat <- subset(dmrs_all_BSobj_BSmooth_tstat, n >=3 & abs(meanDiff) > 0.1 & invdensity <= 300)
  
  null_dmrs <- getNullDistribution_BSmooth.tstat(BSobj_BSmooth_keep, idxMatrix1, idxMatrix2, cutoff = c(-cutoff, cutoff), mc.cores = mc.cores)

  null_dmrs_select <- lapply(null_dmrs, subsetDmrs)
  FWER <- getFWER(null_dmrs_select)
  dmrs_select_BSobj_BSmooth_tstat$FWER <- FWER
  dmrs_select_FWER <- subset(dmrs_select_BSobj_BSmooth_tstat, FWER < 0.05*nperm)
  
  dmrs_passed = rbind(dmrs_passed, dmrs_select_BSobj_BSmooth_tstat)
  dmrs_gold = rbind(dmrs_gold, dmrs_select_FWER)
  
  write.dmrs_bed(dmrs_select_BSobj_BSmooth_tstat,paste(outprefix,chrom,"DMRs.bed",sep="_"),paste(outprefix,chrom,"DMRs",sep="_"),genome)
  if(length(dmrs_select_FWER[,1]) > 0)
    write.dmrs_bed(dmrs_select_FWER,paste(outprefix,chrom,"DMRs_gold.bed",sep="_"),paste(outprefix,chrom,"DMRs_gold",sep="_"),genome)
  
  if(length(dmrs_select_BSobj_BSmooth_tstat[,1]) > 0){
    dmrfilename = paste(outprefix,chrom,"dmrs.pdf",sep="_")
    pdf(file = dmrfilename, width = 7.5, height = 3.67)  #Opens file for figures
    plotManyRegions(BSobj_BSmooth, dmrs_select_BSobj_BSmooth_tstat, extend = 5000, addRegions = dmrs_select_BSobj_BSmooth_tstat,lwd = rep(1.5, 22), BSseqTstat = BSobj_BSmooth_tstat, stat.ylim = c(-8, 8), stat.lwd = 1.5, annoTrack = genome_list)
    dev.off()
  }else{cat("Warning, no DMRs passed gold standard in" , chrom ,"\n")}
  rm(dmrs_all_BSobj_BSmooth_tstat,dmrs_select_BSobj_BSmooth_tstat,
     BSobj_BSmooth,pData,BS.cov,keepLoci.ex,BSobj_BSmooth_keep,
     BSobj_BSmooth_tstat,cutoff,dmrs_select_FWER) #cleanup memory
}

write.dmrs_bed(dmrs_passed,paste(outprefix,"_DMRs.bed",sep=""),paste(outprefix,"_DMRs",sep=""),genome)
write.dmrs_bed(dmrs_gold,paste(outprefix,"_DMRs_gold.bed",sep=""),paste(outprefix,"_DMRs_gold",sep=""),genome)

write.table(dmrs_passed, paste(outprefix,"_Select_DMRs.txt", sep = ""), sep = "\t")
write.table(dmrs_gold, paste(outprefix,"_Gold_DMRs.txt", sep = ""), sep = "\t")

