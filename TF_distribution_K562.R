#   Title: Transcription Factors from K562 cell line
##  Description: This script makes a density graphic from number of transcription factors per promoter
######################################
#Set working Directory
setwd("")
######################################
#Read Files
eprom <- read.csv(file="k562_active_2.txt", header = T, sep = "\t")
control <-  read.csv(file="k562_control_2.txt", header = T, sep = "\t")
set_control <-  read.csv(file="k562_set_control_2.txt", header = T, sep = "\t")
both <-  read.csv(file="k562_active_both.txt", header = T, sep = "\t")
#####################################################################################
#ransform to marix data
mat_eprom <- (data.matrix(eprom[,2:ncol(eprom)]))
mat_control <- (data.matrix(control[,2:ncol(control)]))
mat_set_control <- (data.matrix(set_control[,2:ncol(set_control)]))
mat_both <- (data.matrix(both[,2:ncol(both)]))
#obtain sumatory row
sum_eprom <- rowSums(mat_eprom)
sum_control <- rowSums(mat_control)
sum_set_control <- rowSums(mat_set_control)
sum_both <- rowSums(mat_both)
#Make tables
breaks = seq(min(c(min(sum_eprom), min(sum_control), min(sum_set_control)))-0.1, max(c(max(sum_eprom), max(sum_control), max(sum_set_control)))+0.1, by=0.01) 
epromoter.cut = cut(sum_eprom, breaks, right=FALSE) 
control.cut = cut(sum_control, breaks, right=FALSE) 
set_control.cut = cut(sum_set_control, breaks, right=FALSE) 
both.cut = cut(sum_both, breaks, right=FALSE) 
#tables
epromoter.freq = as.data.frame(table(epromoter.cut))
control.freq = as.data.frame(table(control.cut))
set_control.freq = as.data.frame(table(set_control.cut))
both.freq = as.data.frame(table(both.cut))
#####################################################################################
breaks.names <- as.data.frame(breaks)
#epromoter
index1 <- as.data.frame(which(epromoter.freq$Freq != 0))
values1 <-  (as.data.frame( epromoter.freq$Freq[which(epromoter.freq$Freq != 0)]))/553
index1.match <- as.data.frame(breaks[index1$`which(epromoter.freq$Freq != 0)`])
colnames(index1.match) <- c("value")
#control
index2 <- as.data.frame(which(control.freq$Freq != 0))
values2 <-  (as.data.frame(control.freq$Freq[which(control.freq$Freq != 0)]))/12722
index2.match <- as.data.frame(breaks[index2$`which(control.freq$Freq != 0)`])
colnames(index2.match) <- c("value")
#set control
index3 <- as.data.frame(which(set_control.freq$Freq != 0))
values3 <-  (as.data.frame(set_control.freq$Freq[which(set_control.freq$Freq != 0)]))/649
index3.match <- as.data.frame(breaks[index3$`which(set_control.freq$Freq != 0)`])
colnames(index3.match) <- c("value")
#both
index4 <- as.data.frame(which(both.freq$Freq != 0))
values4 <-  (as.data.frame(both.freq$Freq[which(both.freq$Freq != 0)]))/145
index4.match <- as.data.frame(breaks[index4$`which(both.freq$Freq != 0)`])
colnames(index4.match) <- c("value")
####################################
#Add last row with cero
newarrow1 = c(76)
newarrow2 = c(0)
newarrow3 = c(0.254237288)

index1.match = rbind(index1.match, newarrow1)
values1 = rbind(values1, newarrow2)
#
index4.match = rbind(index4.match, newarrow1)
values4 = rbind(values4, newarrow2)
#
index3.match = rbind(index3.match, newarrow2)
values3 = rbind(values3, newarrow3)
#####################################################################################
#Smooth Data
smt <- 0.7  #smooth level  0 to 1
smooth.eprom = smooth.spline(index1.match$value, values1$`epromoter.freq$Freq[which(epromoter.freq$Freq != 0)]`, spar=smt)
#smooth.control = smooth.spline(index2.match$value, values2$`control.freq$Freq[which(control.freq$Freq != 0)]`, spar=0.6)
smooth.set_control = smooth.spline(index3.match$value, values3$`set_control.freq$Freq[which(set_control.freq$Freq != 0)]`, spar=smt)
smooth.both = smooth.spline(index4.match$value, values4$`both.freq$Freq[which(both.freq$Freq != 0)]`, spar=smt)
#########################################################
#Stats
a <-as.data.frame(sum_eprom)
b <-as.data.frame(sum_control)
c <-as.data.frame(sum_set_control)
d <-as.data.frame(sum_both) 
stat1 <- ks.test(a$sum_eprom, b$sum_control, alternative = "l")  
PVALUE1 = as.numeric(stat1$p.value)
PVALUE1 <- signif(PVALUE1, digits=4)
stat2 <- ks.test(a$sum_eprom, c$sum_set_control, alternative = "l")  
PVALUE2 = as.numeric(stat2$p.value)
PVALUE2 <- signif(PVALUE2, digits=4)
#
stat3 <- ks.test(d$sum_both, c$sum_set_control, alternative = "l")  
PVALUE3 = as.numeric(stat3$p.value)
PVALUE3 <- signif(PVALUE3, digits=4)
stat4 <- ks.test(d$sum_both, b$sum_control, alternative = "l")  
PVALUE4 = as.numeric(stat4$p.value)
PVALUE4 <- signif(PVALUE4, digits=4)
####################################################################################
#Plot Epromoters, Epromoters (hela and k562), set control
cairo_pdf("TF_dist_K562.pdf", 8, 6, bg="transparent")
plot(smooth.set_control, col = "darkblue", type = "l", lwd=9, ylim=c(0.0055, 0.20) , xlim=c(2.5, 73.8), xlab = "Number of transcription factors per promoter", ylab = "Density", main = "K562 transcription factors distribution")
lines(smooth.both, col = "chocolate", lwd=9)
lines(smooth.eprom, col = "firebrick", lwd=9)
legend(x=35, y=0.18, bty = "n", y.intersp= 0.6, cex = 1, legend=c("EPromoters", "EPromoters(K562 & HeLa)", "Set Control"), fill=c("firebrick3", "chocolate", "darkblue"))
legend(25,0.11,legend=paste("p-val (Kolmogorov, EProm > Set Control):", PVALUE2, sep = ""), cex = 0.8, bty = "n")
legend(25,0.10,legend=paste("p-val (Kolmogorov, EProm (common)> Set Control):", PVALUE3, sep = ""), cex = 0.8, bty = "n")
dev.off()