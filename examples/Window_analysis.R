#Author: Keith Dunaway
library(ggplot2)

############################################
# Functions
############################################
windowprocess = function(df,countdf, exp_cols, ctrl_cols){
  # Function that processes 2 window % and count data.frames
  #   and makes them into a single data.frame with other info
  colnames(countdf) = paste("count", names(df), sep = "_")
  totalcols = c(exp_cols,ctrl_cols)
  df = cbind(df, countdf[,totalcols])
  
  df$methdif = rowMeans(df[,exp_cols]) - rowMeans(df[,ctrl_cols])
  df$countdif = rowMeans(countdf[,exp_cols]) - rowMeans(countdf[,ctrl_cols])
  direction_hyho = function(x){
    if(x > 0)
      return("hyper")
    else
      return("hypo")  
  }
  df$direction = sapply(1:length(df$methdif),function(x) direction_hyho(df$methdif[x]))  
  df$CpGs = round(rowMeans(countdf[,totalcols]))
  t.test.pval <- function(...) {
    obj<-try(t.test(...), silent=TRUE)
    if (is(obj, "try-error")) return(.999) else return(obj$p.value)
  }
  df$ttest = sapply(1:length(df[,1]), function(n) t.test.pval(x = df[n,exp_cols], y = df[n,ctrl_cols]))
  df$direction <- factor(df$direction, levels = c("hypo","hyper"))  
  df$countttest = sapply(1:length(df[,1]), function(n) t.test.pval(x = countdf[n,exp_cols], y = df[n,ctrl_cols]))  
  return (df)
}
WinPVal = function (sig_num, window_size, p_val = .05){
  # Function to determine the p value of finding multiple p values in a given window.
  q_val = 1-p_val
  if(sig_num < 1)
    return(1)
  total = q_val ^ window_size
  x=1
  while (x < sig_num) {
    result = choose(window_size,x) * q_val ^ (window_size-x) * p_val^x
    total = total + result
    x = x+1
  }
  return(1-total)
}
MultiHyp_WinPVal = function(sample_size,sig_num, window_size, p_val = .05){
  # Function to determine the p value of finding multiple p values in a given window given
  #  multiple sampling occuring.
  return(sample_size * WinPVal(sig_num, window_size, p_val))
}
GetRange_WinPVal = function(range, sample_size, pval = .05, hypval = .05){
  # Function that gets range of MultiHyp_WinPVal with multi-hyp values adjusted
  results = c()
  for(i in range){
    t = i - 1
    r = 0
    while(r < hypval){
      t=t+1
      r = MultiHyp_WinPVal(sample_size, i, t, pval)
    }
    results=c(results,i,t-1)
  }
  results = paste(results, collapse = ' ')
  return(results)
}
plotcor=function(x,y, rx=.8, ry=.15, pchval=16, cexval=.1, xlabel="x axis", ylabel="y axis", mainlabel = "Main"){
  cornums = paste("r=",round(cor(x,y), digits=2))
  plot(x,y,xlab=xlabel,ylab=ylabel,main=mainlabel,pch=pchval,cex=cexval) + text(rx,ry,cornums)
  abline(0,1)
}
makerandomwindows=function(windows=136484,samples=6){
  df = cbind(1:windows,1:windows,1:windows)
  dfnames=c("chr","start","stop")
  for (i in 1:samples ) {
    n = runif(windows, 0, 1)
    df = cbind(df,n)
    dfnames = c(dfnames,paste("samp_",i))
  }
  df = as.data.frame(df)
  names(df)<-dfnames
  return(df)
}

############################################
# Windowing and Clustering analyses
############################################
#Load and preprocess data
meth_20K = read.table("Window_20k20c1r12f.txt", header = TRUE, sep = "\t") 
count_20K = read.table("Window_20k20c1r12f.txt.count", header = TRUE, sep = "\t") 
Dup15_20K = windowprocess(meth_20K, count_20K, 4:9, 10:15)
rm(meth_20K,count_20K)

#P-value from ttest vs direction
Dup15_20K$direction <- factor(Dup15_20K$direction,levels=levels(Dup15_20K$direction)[order(levels(Dup15_20K$direction), decreasing = TRUE)])
GetRange_WinPVal(6:18,137363)
# Returns: 
# 6 7 7 11 8 16 9 21 10 26 11 33 12 39 13 47 14 54 15 62 16 71 17 79 18 88

#writing table to outfile
write.table(Dup15_20K, "Dup15_20K.txt", quote=F, sep = "\t", row.names = F)
# Run perl command window_cluster.pl using above parameters and file
# See examples/Window_Cluster.bash for example
hypowin_incluster =subset(Dup15_20K_sig, chr=="chr15" & start > 23020000 & end < 26060000 & direction == "hypo")
hyperwin_incluster =subset(Dup15_20K_sig, chr=="chr15" & start > 23020000 & end < 26060000 & direction == "hyper")


####################################################
# Graphing Windows over significant chr15 cluster
####################################################
Cluter_20Kwind <- subset(BrainNOCGI_20K, chr == "chr15" & start >= 23020000 & end <= 26040000)
Cluter_20Kwind$Ctrl = (Cluter_20Kwind$JLKD001 + Cluter_20Kwind$JLKD002 + Cluter_20Kwind$JLKD003 + Cluter_20Kwind$JLKD004 + Cluter_20Kwind$JLKD005 + Cluter_20Kwind$JLKD014 + Cluter_20Kwind$JLKD040 + Cluter_20Kwind$JLKD041 + Cluter_20Kwind$JLKD042 + Cluter_20Kwind$JLKD026 + Cluter_20Kwind$JLKD028)/11
Cluter_20Kwind$Dup15 = (Cluter_20Kwind$JLKD006 + Cluter_20Kwind$JLKD007 + Cluter_20Kwind$JLKD008 + Cluter_20Kwind$JLDS017 + Cluter_20Kwind$JLDS018 + Cluter_20Kwind$JLDS019 + Cluter_20Kwind$JLKD023 + Cluter_20Kwind$JLKD021+ Cluter_20Kwind$JLKD022)/9
Cluter_20Kwind$PWS = (Cluter_20Kwind$JLKD038 + Cluter_20Kwind$JLKD039 + Cluter_20Kwind$JLKD012 + Cluter_20Kwind$JLKD017 + Cluter_20Kwind$JLKD037) / 5
Cluter_20Kwind$AS = (Cluter_20Kwind$JLKD011 + Cluter_20Kwind$JLKD015 + Cluter_20Kwind$JLKD033)/3
Cluter_20Kwind$Downs = (Cluter_20Kwind$JLKD010 + Cluter_20Kwind$JLKD034 + Cluter_20Kwind$JLKD035 + Cluter_20Kwind$JLKD036)/4
Cluter_20Kwind$Idio = (Cluter_20Kwind$JLKD009 + Cluter_20Kwind$JLKD013 + Cluter_20Kwind$JLKD016 + Cluter_20Kwind$JLKD018 + Cluter_20Kwind$JLKD019 + Cluter_20Kwind$JLKD020)/6
                       
ggplot(Cluter_20Kwind) + 
  geom_line(aes(x=start, y=Ctrl/Ctrl), color = "black" , size = 1) + 
  geom_line(aes(x=start, y=Dup15/Ctrl), color = "blue" , size = 1) + 
  geom_line(aes(x=start, y=PWS/Ctrl), color = "red" , size = 1) + 
  geom_line(aes(x=start, y=AS/Ctrl), color = "green" , size = 1) + 
  geom_line(aes(x=start, y=Downs/Ctrl), color = "purple" , size = 1) + 
  geom_line(aes(x=start, y=Idio/Ctrl), color = "light blue" , size = 1) + 
  xlim(23020000, 26040000) +
  theme_bw()  
