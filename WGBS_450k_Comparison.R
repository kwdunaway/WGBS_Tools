################################################################################
# Exphor: Roy Chu & Keith Dunaway
# Email: rgchu@ucdavis.edu & kwdunaway@ucdavis.edu
# Date: 4-11-2016
#
#
################################################################################


# This hg38 pipeline takes HM450 probe locations out of WGBS methylation data and
# analyzes these locations. There are some preprocessing steps that must be
# manually done.

# STEP 1
# Before this step, the WGBS data should have been converted to percentage
# methylation folders for each sample. For this step, run this script:
# AvgMeth.2col.pl [1] [2] [3] [4] [5] [6] [7]
# [1], [2], etc. represent inputs
# The inputs will be listed when you run AvgMeth.2col.pl with no inputs.
# [2] is the HM450_hg38.bed file.
# [3] is 3 (the fourth column in the HM450_hg38.bed file).
# [8/10/12/etc.] and [9/11/13/etc.] are additional folders/samples in the same
# format as [6] and [7]

# STEP 2
# Load the file onto your local machine and 
# *****************************************
# *****Fill in the following variables*****
# *****************************************

# Experimental sample names
# Example:
# Exp_samples <- c( "JLKD051" ,  "JLKD052" ,  "JLKD054"  ,  "JLKD055"  ,  "JLKD056"  ,  "JLKD057"  ,  "JLKD059"  ,  "JLKD061"  ,  "JLKD097"  ,  "JLKD009"  ,  "JLKD013" ,  "JLKD016"  ,  "JLKD018"  ,  "JLKD019"  ,  "JLKD020");
Exp_samples <- c();

# Control sample names
# Example:
# Ctl_samples <- c( "JLKD053"  ,  "JLKD062"  ,  "JLKD063"  ,  "JLKD064"  ,  "JLKD065"  ,  "JLKD066"  ,  "JLKD067"  ,  "JLKD068"  ,  "JLKD069"  ,  "JLKD070"  ,  "JLKD071"  ,  "JLKD080"  ,  "JLKD082"  ,  "JLKD083"  ,  "JLKD084"  ,  "JLKD098"  ,  "JLKD099"  ,  "JLKD002"  ,  "JLKD003"  ,  "JLKD005"  ,  "JLKD001"  ,  "JLKD014"  ,  "JLKD004"  ,  "JLKD040"  ,  "JLKD041"  ,  "JLKD042"  ,  "JLKD025"  ,  "JLKD027");
Ctl_samples <- c();

# Place the full file path
# Example: 
# WGBS_450k_filepath <- "C:/Users/Roy/Dropbox/R/WGBS_450k_Comparison_4_2016/URC_SingleCpG.txt.2col"
WGBS_450k_filepath <- ""



# STEP 3
# Run this R script

############################
# Call Necessary Libraries #
############################

# This statement will cycle through all necessary packages (in this case,
# just ggplot2) and install/call them. This is necessary to use them.

for (package in c('ggplot2', 'plyr', 'reshape2')) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}

####################
# Create Functions #
####################

# This is a function for calculating standard error
se <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))

# This is a function for detecting NaN (not the same as NA)
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))


####################
#   Read in file   #
####################

WGBS_450k_file <- read.delim(WGBS_450k_filepath)


######################################################
#   You can follow this process or create your own   #
######################################################

# _count is count of samples
# _meth is column for meth count data
# _total is column for CpG count data
# _PerMeth is percentage methylation

Exp_count <- 0
Exp_meth <- cbind(WGBS_450k_file[,1:4])
Exp_total <- cbind(WGBS_450k_file[,1:4])
Exp_PerMeth <- cbind(WGBS_450k_file[,1:4])
for(i in Exp_samples)  # for each row
{
  #cat(i)
  MethandTot <- grep(i,colnames(WGBS_450k_file))
  Exp_meth <- cbind(Exp_meth,WGBS_450k_file[MethandTot[1]])
  Exp_total <- cbind(Exp_total,WGBS_450k_file[MethandTot[2]])
  Exp_PerMeth <- cbind(Exp_PerMeth,WGBS_450k_file[MethandTot[1]]/WGBS_450k_file[MethandTot[2]])
  Exp_count <- Exp_count + 1
}

cat("Number of experimental samples")
Exp_count

Ctl_count <- 0
Ctl_meth <- cbind(WGBS_450k_file[,1:4])
Ctl_total <- cbind(WGBS_450k_file[,1:4])
Ctl_PerMeth <- cbind(WGBS_450k_file[,1:4])
for(i in Ctl_samples)  # for each row
{
  #cat(i)
  MethandTot <- grep(i,colnames(WGBS_450k_file))
  Ctl_meth <- cbind(Ctl_meth,WGBS_450k_file[MethandTot[1]])
  Ctl_total <- cbind(Ctl_total,WGBS_450k_file[MethandTot[2]])
  Ctl_PerMeth <- cbind(Ctl_PerMeth,WGBS_450k_file[MethandTot[1]]/WGBS_450k_file[MethandTot[2]])
  Ctl_count <- Ctl_count + 1
}

cat("Number of control samples")
Ctl_count

# Aggregate Samples and find percentage methylation

Exp_dim <- Exp_count+4
Ctl_dim <- Ctl_count+4

PerMeth_aggregated <- cbind(WGBS_450k_file[,1:4]) 

PerMeth_aggregated <- transform(PerMeth_aggregated, Exp_meth=rowSums(Exp_meth[,5:Exp_dim], na.rm=TRUE))
PerMeth_aggregated <- transform(PerMeth_aggregated, Exp_total=rowSums(Exp_total[,5:Exp_dim], na.rm=TRUE))
PerMeth_aggregated <- transform(PerMeth_aggregated, Exp_percent=Exp_meth/Exp_total)
PerMeth_aggregated[is.nan(PerMeth_aggregated)] <- NA

PerMeth_aggregated <- transform(PerMeth_aggregated, Ctl_meth=rowSums(Ctl_meth[,5:Ctl_dim], na.rm=TRUE))
PerMeth_aggregated <- transform(PerMeth_aggregated, Ctl_total=rowSums(Ctl_total[,5:Ctl_dim], na.rm=TRUE))
PerMeth_aggregated <- transform(PerMeth_aggregated, Ctl_percent=Ctl_meth/Ctl_total)
PerMeth_aggregated <- transform(PerMeth_aggregated, Avg_percent=(Exp_percent+Ctl_percent)/2)
PerMeth_aggregated <- transform(PerMeth_aggregated, Norm_percent = (Exp_percent - Ctl_percent))
PerMeth_aggregated[is.nan(PerMeth_aggregated)] <- NA



######################
# Hypothesis Testing #
######################

t.test(PerMeth_aggregated$Exp_percent, PerMeth_aggregated$Ctl_percent,paired=TRUE,alternative = "greater", na.rm=TRUE)
t.test(PerMeth_aggregated$Exp_percent, PerMeth_aggregated$Ctl_percent,paired=TRUE,alternative = "less", na.rm=TRUE)
t.test(PerMeth_aggregated$Exp_percent, PerMeth_aggregated$Ctl_percent,paired=TRUE,alternative = "two.sided", na.rm=TRUE)

############
# Graphing #
############

# Scatterplot

Exp_Ctl_scatter <- ggplot(data=PerMeth_aggregated) + labs(x="Control Percentage Methylation", y="Experimental Percentage Methylation")
Exp_Ctl_scatter <- Exp_Ctl_scatter + theme_bw() + theme(plot.title = element_text(size = rel(1.3), face = "bold", colour = "black", vjust = 2), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
Exp_Ctl_scatter + geom_point(data=PerMeth_aggregated, aes(y = PerMeth_aggregated$Exp_percent, x = PerMeth_aggregated$Ctl_percent, colour = "PerMeth")) + geom_abline() + scale_colour_manual(values=c("#56B4E9"))

# Bar Graph

PerMeth_bar_data <- data.frame(c("Exp","Ctl"))
colnames(PerMeth_bar_data) <- c("Variable")
PerMeth_bar_data$percent <- c(mean(PerMeth_aggregated$Exp_percent, na.rm=TRUE), mean(PerMeth_aggregated$Ctl_percent, na.rm=TRUE))
PerMeth_bar_data$se <- c(se(PerMeth_aggregated$Exp_percent), se(PerMeth_aggregated$Ctl_percent))

PerMeth_bar_data$Variable <- factor(PerMeth_bar_data$Variable,levels = c("Exp","Ctl"))
PerMeth_bar_data.m <- melt(PerMeth_bar_data, id.vars="Variable")

limits <- aes(x=Variable, ymax = percent + se, ymin=percent - se)
dodge <- position_dodge(width=0.9)

Exp_Ctl_bar <- ggplot(data=PerMeth_bar_data) + labs(x="Variable", y="Percentage Change", title = "Percentage Methylation")
Exp_Ctl_bar <- Exp_Ctl_bar + theme_bw() + theme(plot.title = element_text(size = rel(1.3), face = "bold", colour = "black", vjust = 2), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
# Adjust ylim to your liking and fill_manual to the colors you want
Exp_Ctl_bar + geom_bar(data=PerMeth_bar_data, aes(x=Variable, y=percent, fill = Variable), position = "dodge", stat = "identity") + scale_fill_manual(values=c("#56B4E9","#0072B2")) + coord_cartesian(ylim=c(0.4825, 0.4875)) + geom_hline(yintercept=0) + geom_errorbar(limits, position=dodge, width=0.25)

# Methylation Category

PerMeth_Low <- subset(PerMeth_aggregated, Avg_percent <= 0.10)
PerMeth_Low$methcat <- "low";
PerMeth_Partial <- subset(PerMeth_aggregated, Avg_percent > 0.10 & Avg_percent <= 0.75)
PerMeth_Partial$methcat <- "partial";
PerMeth_High <- subset(PerMeth_aggregated, Avg_percent > 0.75)
PerMeth_High$methcat <- "high";

t.test(PerMeth_Low$Exp_percent, PerMeth_Low$Ctl_percent,paired=TRUE,alternative = "greater")
t.test(PerMeth_Low$Exp_percent, PerMeth_Low$Ctl_percent,paired=TRUE,alternative = "less")
t.test(PerMeth_Low$Exp_percent, PerMeth_Low$Ctl_percent,paired=TRUE,alternative = "two.sided")

t.test(PerMeth_Partial$Exp_percent, PerMeth_Partial$Ctl_percent,paired=TRUE,alternative = "greater")
t.test(PerMeth_Partial$Exp_percent, PerMeth_Partial$Ctl_percent,paired=TRUE,alternative = "less")
t.test(PerMeth_Partial$Exp_percent, PerMeth_Partial$Ctl_percent,paired=TRUE,alternative = "two.sided")

t.test(PerMeth_High$Exp_percent, PerMeth_High$Ctl_percent,paired=TRUE,alternative = "greater")
t.test(PerMeth_High$Exp_percent, PerMeth_High$Ctl_percent,paired=TRUE,alternative = "less")
t.test(PerMeth_High$Exp_percent, PerMeth_High$Ctl_percent,paired=TRUE,alternative = "two.sided")


####
# Average Percentage Bar Graph
###
PerMeth_CatBar <- data.frame(c("low","partial","high"))
colnames(PerMeth_CatBar) <- c("methcat")
PerMeth_CatBar$Exp_percent <- c(mean(PerMeth_Low$Exp_percent), mean(PerMeth_Partial$Exp_percent), mean(PerMeth_High$Exp_percent))
PerMeth_CatBar$Ctl_percent <- c(mean(PerMeth_Low$Ctl_percent), mean(PerMeth_Partial$Ctl_percent), mean(PerMeth_High$Ctl_percent))

# Ordering
PerMeth_CatBar$methcat <- factor(PerMeth_CatBar$methcat,levels = c("low", "partial", "high"))
PerMeth_CatBar.m <- melt(PerMeth_CatBar, id.vars="methcat")
PerMeth_CatBar.m$se <- c(se(PerMeth_Low$Exp_percent), se(PerMeth_Partial$Exp_percent), se(PerMeth_High$Exp_percent), se(PerMeth_Low$Ctl_percent), se(PerMeth_Partial$Ctl_percent), se(PerMeth_High$Ctl_percent))

dodge <- position_dodge(width=0.9)

PerMeth_CatBar_graph <- ggplot(data=PerMeth_CatBar.m, aes(x=methcat, y=value, fill = variable)) + labs(x="Methylation Category", y="Average Percentage", title = "Average Percentage")
PerMeth_CatBar_graph <- PerMeth_CatBar_graph+ theme_bw() + theme(plot.title = element_text(size = rel(1.3), face = "bold", colour = "black", vjust = 2), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
# Adjust ylim to your liking and fill_manual to the colors you want
PerMeth_CatBar_graph + geom_bar(position = "dodge", stat = "identity") + scale_fill_manual(values=c("#56B4E9","#0072B2", "#9999CC")) + coord_cartesian(ylim=c(-1,1)) + geom_hline(yintercept=0) + geom_errorbar(aes(ymax = value + se, ymin = value - se), position=dodge, width=0.25)

####
# Normalized Bar Graph
###
PerMeth_CatBar <- data.frame(c("low","partial","high"))
colnames(PerMeth_CatBar) <- c("methcat")
PerMeth_CatBar$Norm_percent <- c(mean(PerMeth_Low$Exp_percent) - mean(PerMeth_Low$Ctl_percent), mean(PerMeth_Partial$Exp_percent) - mean(PerMeth_Partial$Ctl_percent), mean(PerMeth_High$Exp_percent - mean(PerMeth_High$Ctl_percent)))
PerMeth_CatBar$se <- c(se(PerMeth_Low$Norm_percent), se(PerMeth_Partial$Norm_percent), se(PerMeth_High$Norm_percent))

# Ordering
PerMeth_CatBar$methcat <- factor(PerMeth_CatBar$methcat,levels = c("low", "partial", "high"))
PerMeth_CatBar.m <- melt(PerMeth_CatBar, id.vars="methcat")

limits <- aes(x=methcat, ymax = Norm_percent + se, ymin=Norm_percent - se)
dodge <- position_dodge(width=0.9)

PerMeth_CatBar_graph <- ggplot(data=PerMeth_CatBar) + labs(x="Methylation Category", y="Normalized Percentage Change", title = "BA10 Hypermethylated")
PerMeth_CatBar_graph <- PerMeth_CatBar_graph+ theme_bw() + theme(plot.title = element_text(size = rel(1.3), face = "bold", colour = "black", vjust = 2), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
# Adjust ylim to your liking and fill_manual to the colors you want
PerMeth_CatBar_graph + geom_bar(data=PerMeth_CatBar, aes(x=methcat, y=Norm_percent, fill = methcat), position = "dodge", stat = "identity") + scale_fill_manual(values=c("#56B4E9","#0072B2", "#9999CC")) + coord_cartesian(ylim=c(-0.005,0.005)) + geom_hline(yintercept=0) + geom_errorbar(limits, position=dodge, width=0.25)


