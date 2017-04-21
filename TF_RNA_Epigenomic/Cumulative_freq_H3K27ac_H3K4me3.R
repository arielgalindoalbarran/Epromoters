#   Title: Cumulative frequency H3K27ac/H3K4me3
##  Description: This script makes a graphic from H3K27ac and H3K4me3 fold cahnge ammounts in Epromoters and control (same expression)
### in K562 and HeLa cell lines.

######################################################################################
#Set working directory
setwd(" ")
######################################################################################
#Read coverage area from Ip files:
##K562 cell line
k562_H3k4me3 <- read.csv(file="ENCFF001FWX.tab", header = FALSE, sep = "\t") #H3k4me3
k562_H3k27ac <- read.csv(file="ENCFF000BWY.tab", header = FALSE, sep = "\t") #H3k27ac
k562_H3k4me1 <- read.csv(file="ENCFF000BXQ.tab", header = FALSE, sep = "\t") #H3k4me1
##HeLa cell line
hela_H3k4me3 <- read.csv(file="ENCFF000BCP.tab", header = FALSE, sep = "\t") #H3k4me3
hela_H3k27ac <- read.csv(file="ENCFF000BBR.tab", header = FALSE, sep = "\t") #H3k27ac
hela_H3k4me1 <- read.csv(file="ENCFF000BBF.tab", header = FALSE, sep = "\t") #H3k4me1

#Read Epromoters and control gene list: 
##K562  cell line
control_k562 <- read.csv(file="k562_control_list.txt", header = FALSE, sep = "\t")     #Epromoter list
reference_k562 <- read.csv(file="k562_reference_list.txt", header = FALSE, sep = "\t") #control (same exp) list
##heLa cell line
control_hela <- read.csv(file="hela_control_list.txt", header = FALSE, sep = "\t")     #Epromoter list
reference_hela <- read.csv(file="hela_reference_list.txt", header = FALSE, sep = "\t") #control (same exp) list

######################################################################################
#Block 1: Process data and export tables
######################################################################################
#Order
k562_H3k4me3 <- k562_H3k4me3[order(k562_H3k4me3$V1),] 
k562_H3k27ac <- k562_H3k27ac[order(k562_H3k27ac$V1),] 
k562_H3k4me1 <- k562_H3k4me1[order(k562_H3k4me1$V1),] 

hela_H3k4me3 <- hela_H3k4me3[order(hela_H3k4me3$V1),] 
hela_H3k27ac <- hela_H3k27ac[order(hela_H3k27ac$V1),] 
hela_H3k4me1 <- hela_H3k4me1[order(hela_H3k4me1$V1),] 

#Join Files in one
k562_histonemarks_all <- cbind(k562_H3k4me3, k562_H3k27ac$V2, k562_H3k4me1$V2)
hela_histonemarks_all <- cbind(hela_H3k4me3, hela_H3k27ac$V2, hela_H3k4me1$V2)

#delate rownames
rownames(k562_histonemarks_all) <- NULL
rownames(hela_histonemarks_all) <- NULL

#Change column names
colnames(k562_histonemarks_all) <- c("gene", "H3K4me3", "H3K27ac", "H3k4me1")  
colnames(hela_histonemarks_all) <- c("gene", "H3K4me3", "H3K27ac", "H3k4me1")

######################################################################################
#Export like tables
##k562 cell line
write.table(k562_histonemarks_all, file = "k562_histonemarks_all.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
##heLa cell line
write.table(hela_histonemarks_all, file = "hela_histonemarks_all.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

#Important note: Some gene names was restored to original name because excell autoformat (sep15 for example). 
######################################################################################
#Read last files:
k562_histonemarks_all <- read.csv(file="k562_histonemarks_all.txt", header = TRUE, sep = "\t")
hela_histonemarks_all <- read.csv(file="hela_histonemarks_all.txt", header = TRUE, sep = "\t")

#Make an consolidate by meanvof repeated data.  
##K562 cell line
consol.k562_histonemarks_all <- as.data.frame(t(sapply(split(k562_histonemarks_all,k562_histonemarks_all$gene), function(y) sapply(y,mean))))
consol.k562_histonemarks_all$gene <- NULL
gene <- rownames(consol.k562_histonemarks_all)
rownames(consol.k562_histonemarks_all) <- NULL
consol.k562_histonemarks_all <- cbind(gene,consol.k562_histonemarks_all)
remove (gene)
##heLa cell line
consol.hela_histonemarks_all <- as.data.frame(t(sapply(split(hela_histonemarks_all,hela_histonemarks_all$gene), function(y) sapply(y,mean))))
consol.hela_histonemarks_all$gene <- NULL
gene <- rownames(consol.hela_histonemarks_all)
rownames(consol.hela_histonemarks_all) <- NULL
consol.hela_histonemarks_all <- cbind(gene,consol.hela_histonemarks_all)
remove (gene)

#Get gene names from each list
k562_referencenames <- reference_k562$V1
k562_controlnames <- control_k562$V1
hela_referencenames <- reference_hela$V1
hela_controlnames <- control_hela$V1
#Note, if genes are the same, use alldatasetnames, if not:
alldatasetnamesk562 <- consol.k562_histonemarks_all$gene
alldatasetnameshela <- consol.hela_histonemarks_all$gene

#Compare and get positions from data file above:
##k562 cell line
k562_control_values <- intersect(alldatasetnamesk562, k562_controlnames)
k562_reference_values <- intersect(alldatasetnamesk562, k562_referencenames) 
##heLa cell line
hela_control_values <- intersect(alldatasetnameshela,hela_controlnames)
hela_reference_values <- intersect(alldatasetnameshela,hela_referencenames)

#Get positions to one variable:
##k562 cell line
selected_k562_control <- c(k562_control_values)
selected_k562_reference <- c(k562_reference_values) 
##heLa cell line
selected_hela_control <- c(hela_control_values)
selected_hela_reference <- c(hela_reference_values)

#generate index
##k562 cell line
index_k562_control <- unlist(sapply (selected_k562_control, function(x){ which(consol.k562_histonemarks_all$gene == x)} ))
index_k562_control.names <- names(index_k562_control)  
index_k562_reference <- unlist(sapply (selected_k562_reference, function(x){ which(consol.k562_histonemarks_all$gene == x)} )) #con problema
index_k562_reference.names <- names(index_k562_reference)                 
##heLa cell line
index_hela_control <- unlist(sapply (selected_hela_control, function(x){ which(consol.hela_histonemarks_all$gene == x)} ))
index_hela_control.names <- names(index_hela_control)
index_hela_reference <- unlist(sapply (selected_hela_reference, function(x){ which(consol.hela_histonemarks_all$gene == x)} ))
index_hela_reference.names <- names(index_hela_reference)


#Generate export files:
#k562 cell line
selected.genes.table.k562.control <- consol.k562_histonemarks_all[index_k562_control, ]
rownames(selected.genes.table.k562.control) <- index_k562_control.names
selected.genes.table.k562.reference <- consol.k562_histonemarks_all[index_k562_reference, ]
rownames(selected.genes.table.k562.reference) <- index_k562_reference.names
#heLa cell line
selected.genes.table.hela.control <- consol.hela_histonemarks_all[index_hela_control, ]
rownames(selected.genes.table.hela.control) <- index_hela_control.names
selected.genes.table.hela.reference <- consol.hela_histonemarks_all[index_hela_reference, ]
rownames(selected.genes.table.hela.reference) <- index_hela_reference.names

#Generate tables
##k562
selected.genes.table.k562.control.file <- "k562_control_TF_list.tab"
selected.genes.table.k562.reference.file <- "k562_control_TF_list.tab"
##hela
selected.genes.table.hela.control.file <- "k562_control_TF_list.tab"
selected.genes.table.hela.reference.file <- "k562_control_TF_list.tab"


#Export like tables:
##k562 cell line
write.table(selected.genes.table.k562.control, file = "Histone_marks_k562_control.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(selected.genes.table.k562.reference, file = "Histone_marks_k562_reference.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
##heLa cell line
write.table(selected.genes.table.hela.control, file = "Histone_marks_hela_control.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(selected.genes.table.hela.reference, file = "Histone_marks_hela_reference.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


######################################################################################
#Block 2: Plot data
######################################################################################

#Read files
#K562 cell line
k562_control <- read.csv(file="Histone_marks_k562_control.txt", header = T, sep = "\t")
k562_reference <- read.csv(file="Histone_marks_k562_reference.txt", header = T, sep = "\t")
#HeLa cell line
hela_control <- read.csv(file="Histone_marks_hela_control.txt", header = T, sep = "\t")
hela_reference <- read.csv(file="Histone_marks_hela_reference.txt", header = T, sep = "\t")
#All
k562_all <- read.csv(file="k562_histonemarks_inactives.txt", header = T, sep = "\t")
hela_all <- read.csv(file="hela_histonemarks_inactives.txt", header = T, sep = "\t")
######################################################################################
#Make the cumulative sum:

##K562 cell line
x <- (k562_control$H3K27ac+1)/(k562_control$H3K4me3+1)
y <- (k562_reference$H3K27ac+1)/(k562_reference$H3K4me3+1)
###ranks
breaks_x_y = seq(min(c(min(x), min(y)))-0.01, max(c(max(x), max(y)))+0.01, by=0.1) 
x.cut = cut(x, breaks_x_y, right=FALSE) 
y.cut = cut(y, breaks_x_y, right=FALSE) 
###tables
x.freq = table(x.cut)
y.freq = table(y.cut)
###cumulatves frequences
cum_x <- cumsum(x.freq)
cum_y <- cumsum(y.freq)
###Control random:
k562_control_all <- (k562_all$H3K27ac+1)/(k562_all$H3K4me3+1)
k562_random <- sample(k562_control_all, 528)
k562_random.cut = cut(k562_random, breaks_x_y, right=FALSE) 
k562_random.freq = table(k562_random.cut)
cum_k562_random <- cumsum(k562_random.freq)


##HeLa cell line
k <- (hela_control$H3K27ac+1)/(hela_control$H3K4me3+1)
l <- (hela_reference$H3K27ac+1)/(hela_reference$H3K4me3+1)
#ranks
breaks_k_l = seq(min(c(min(k), min(l)))-0.01, max(c(max(k), max(l)))+0.01, by=0.1) 
k.cut = cut(k, breaks_k_l, right=FALSE) 
l.cut = cut(l, breaks_k_l, right=FALSE) 
#tables
k.freq = table(k.cut)
l.freq = table(l.cut)
#cumulatves frequences
cum_k <- cumsum(k.freq)
cum_l <- cumsum(l.freq)
#Control random:
hela_control_all <- (hela_all$H3K27ac+1)/(hela_all$H3K4me3+1)
hela_random <- sample(hela_control_all, 462)
hela_random.cut = cut(hela_random, breaks_k_l, right=FALSE) 
hela_random.freq = table(hela_random.cut)
cum_hela_random <- cumsum(hela_random.freq)  
######################################################################################
#Statistics
##K562 cell line
stat.k562 <- wilcox.test(x,y, paired = TRUE) 
PVALUE.k562 = as.numeric(stat.k562$p.value)
PVALUE.k562 <- signif(PVALUE.k562, digits=4)
###Random
stat.k562.random <- wilcox.test(k562_random, y, paired = TRUE) 
PVALUE.k562.random = as.numeric(stat.k562.random$p.value)
PVALUE.k562.random <- signif(PVALUE.k562.random, digits=4)
##HeLa cell line
stat.hela <- wilcox.test(k,l, paired = TRUE) 
PVALUE.hela = as.numeric(stat.hela$p.value)
PVALUE.hela <- signif(PVALUE.hela, digits=4)
###Random
stat.hela.random <- wilcox.test(hela_random, l, paired = TRUE) 
PVALUE.hela.random = as.numeric(stat.hela.random$p.value)
PVALUE.hela.random <- signif(PVALUE.hela.random, digits=4)
######################################################################################
#Plot data

##K562 cell line
###Add "cero" to the vector
cum_x0 = c(0, cum_x)
cum_y0 = c(0, cum_y)
cum_k562_random0 = c(0, cum_k562_random)
breaks_x_y_log=log2(breaks_x_y)
###Make soft lines
k562_smooth.control = smooth.spline(breaks_x_y_log, cum_x0, spar=1)
k562_smooth.reference = smooth.spline(breaks_x_y_log, cum_y0, spar=1)
k562_smooth.random = smooth.spline(breaks_x_y_log, cum_k562_random0, spar=1)
###graph
plot(breaks_x_y_log, cum_x0, main="K562", 
     xlab="H3K27ac / H3K4me3", ylab="Cumulative frequency", col=c("white"), xlim = c(-3.5,3.5))    
lines(k562_smooth.random, col=("gray47"), lwd=6)
lines(k562_smooth.control, col=("blue"), lwd=6)
lines(k562_smooth.reference, col=("firebrick"), lwd=6)
legend(x=0.5, y=400, bty = "n", y.intersp= 0.7, cex = 0.8, legend=c("EPromoters", "Control promoters (same expression)", "Random promoters (same number)"), fill=c("firebrick3", "blue", "black"))
legend(1.5,200,legend=paste(PVALUE.k562, sep = ""), cex = 0.8, bty = "n")
legend(-1,200,"Wilcox paired test Epromoter-Control promoter:", cex = 0.8, bty = "n")
legend(1.5,150,legend=paste(PVALUE.k562.random, sep = ""), cex = 0.8, bty = "n")
legend(-1,150,"Wilcox paired test Epromoter-Random promoter:", cex = 0.8, bty = "n")

##heLa cell line
###Add "cero" to the vector
cum_k0 = c(0, cum_k)
cum_l0 = c(0, cum_l)
cum_hela_random0 = c(0, cum_hela_random)
breaks_k_l_log=log2(breaks_k_l)
###Make soft lines
hela_smooth.control = smooth.spline(breaks_k_l_log, cum_k0, spar=1)
hela_smooth.reference = smooth.spline(breaks_k_l_log, cum_l0, spar=1)
hela_smooth.random = smooth.spline(breaks_k_l_log, cum_hela_random0, spar=1)
##graph
plot(breaks_k_l_log, cum_k0, main="HELA", 
     xlab="H3K27ac / H3K4me3", ylab="Cumulative frequency", col=c("white"), xlim = c(-3,3))    
lines(hela_smooth.random, col=("gray47"), lwd=6)
lines(hela_smooth.control, col=("blue"), lwd=6)
lines(hela_smooth.reference, col=("firebrick"), lwd=6)
legend(x=0.5, y=350, bty = "n", y.intersp= 0.7, cex = 0.8, legend=c("EPromoters", "Control promoters (same expression)", "Random promoters (same number)"), fill=c("firebrick3", "blue", "black"))
legend(1.5,150,legend=paste(PVALUE.hela, sep = ""), cex = 0.8, bty = "n")
legend(-1,150,"Wilcox paired test Epromoter-Control promoter:", cex = 0.8, bty = "n")
legend(1.5,100,legend=paste(PVALUE.hela.random, sep = ""), cex = 0.8, bty = "n")
legend(-1,100,"Wilcox paired test Epromoter-Random promoter:", cex = 0.8, bty = "n")
