######################################
#Title: Comparation k562 and hela RNAseq and EPromoters
#Description: Analyzing FPKM from RNAseq in different promoters. 
######################################
#Working Directory
setwd("")
##################################################################################
#Read Files
k562_only <- read.csv(file="k562_epromoters_RPKM.txt", header = TRUE, sep = "\t")
hela_only <- read.csv(file="hela_epromoters_RPKM.txt", header = TRUE, sep = "\t")
both <- read.csv(file="both_eprom_RPKM.txt", header = TRUE, sep = "\t")
inactive <-  read.csv(file="inactives_promoters_RPKM.txt", header = TRUE, sep = "\t")
#Variables
k562_only_k562_log2 <- na.omit(log2(k562_only$RPKM_k562_mean))
k562_only_hela_log2 <- na.omit(log2(k562_only$RPKM_hela_.mean))
hela_only_k562_log2 <- na.omit(log2(hela_only$RPKM_k562_mean))
hela_only_hela_log2 <- na.omit(log2(hela_only$RPKM_hela_.mean))
both_only_k562_log2 <- na.omit(log2(both$RPKM_k562_mean))
both_only_hela_log2 <- na.omit(log2(both$RPKM_hela_.mean))
inactive_both_k562_log2 <- na.omit(log2(inactive$RPKM_k562_mean))
inactive_both_hela_log2 <- na.omit(log2(inactive$RPKM_hela_.mean))
#################################################################################
##Statistics
#Compared to Control:
STAT1 <- wilcox.test(k562_only$RPKM_k562_mean, inactive$RPKM_k562_mean, alternative = "g", paired = FALSE)
STAT2 <- wilcox.test(k562_only$RPKM_hela_.mean, inactive$RPKM_hela_.mean, alternative = "g", paired = FALSE)
STAT3 <- wilcox.test(hela_only$RPKM_k562_mean, inactive$RPKM_k562_mean, alternative = "g", paired = FALSE)
STAT4 <- wilcox.test(hela_only$RPKM_hela_.mean, inactive$RPKM_hela_.mean, alternative = "g", paired = FALSE)
STAT5 <- wilcox.test(both$RPKM_k562_mean, inactive$RPKM_k562_mean, alternative = "g", paired = FALSE)
STAT6 <- wilcox.test(both$RPKM_hela_.mean, inactive$RPKM_hela_.mean, alternative = "g", paired = FALSE)
#Compared between cell lines exclusive Epromoters:
STAT7 <- wilcox.test(k562_only$RPKM_k562_mean, k562_only$RPKM_hela_.mean, paired = TRUE)
STAT8 <- wilcox.test(hela_only$RPKM_k562_mean, hela_only$RPKM_hela_.mean, paired = TRUE)
STAT9 <- wilcox.test(both$RPKM_k562_mean, both$RPKM_hela_.mean, paired = TRUE)
#Obtain only the pvalue and the significatives:
PVALUE1 = signif(STAT1$p.value,digits=4)
PVALUE2 = signif(STAT2$p.value,digits=4)
PVALUE3 = signif(STAT3$p.value,digits=4)
PVALUE4 = signif(STAT4$p.value,digits=4)
PVALUE5 = signif(STAT5$p.value,digits=4)
PVALUE6 = signif(STAT6$p.value,digits=4)
PVALUE7 = signif(STAT7$p.value,digits=4)
PVALUE8 = signif(STAT8$p.value,digits=4)
PVALUE9 = signif(STAT9$p.value,digits=4)
################################################################################
#Plot
cairo_pdf("Tranacription_and_enhancer_act_K562.pdf", 8, 6, bg="transparent")
boxplot(k562_only_k562_log2, k562_only_hela_log2, hela_only_k562_log2, hela_only_hela_log2, both_only_k562_log2, both_only_hela_log2, las = 0,   
        at =c(1,2,4,5,7,8), 
        col =c("firebrick3", "blue", "firebrick3", "blue", "firebrick3", "blue"),
        main ="RNAseq K562 and Hela EPromoters",
        ylab ="RPKM (Log2)", 
        xlab = "EPromoters only in K562         EPromoters Only in Hela        EPromoters in common",
        ylim=c(-15,10))
legend(x=6, y=-13, bty = "n", y.intersp= 0.8, legend=c("K562", "Hela"), fill=c("firebrick3", "blue"))
legend(0,11.6,legend=paste("p-val:  ", PVALUE7, sep = ""), cex = 1, bty = "n")
legend(0,10.6,legend=paste("K562 vs Hela"), cex = 1, bty = "n")
legend(4,10.6,legend=paste(PVALUE8, sep = ""), cex = 1, bty = "n")
legend(7,10.6,legend=paste(PVALUE9, sep = ""), cex = 1, bty = "n")
dev.off()