#!/usr/bin/Rscript 

##------------------------------------------------------------------------------
## Guillaume Charbonnier
## Created: 2016-06-02 17h35
##------------------------------------------------------------------------------

"A script to produce boxplot for Salva's article.

Usage: 
	boxplot_ratio.R --tab_k562_k562 <input.tab_k562_k562> --tab_k562_hela <input.tab_k562_hela> --tab_k562_common <input.tab_k562_common> --tab_hela_k562 <input.tab_hela_k562> --tab_hela_hela <input.tab_hela_hela> --tab_hela_common <input.tab_hela_common> --tab_k562_all <input.tab_k562_all> --tab_hela_all <input.tab_hela_all> --output <output.pdf>
                
Options:
	--tab_k562_k562 <tab_k562_k562> 	{input.tab_k562_k562}
	--tab_k562_hela <tab_k562_hela> 	{input.tab_k562_hela}
	--tab_k562_common <tab_k562_common> 	{input.tab_k562_common}
	--tab_hela_k562 <tab_hela_k562> 	{input.tab_hela_k562}
	--tab_hela_hela <tab_hela_hela> 	{input.tab_hela_hela}
	--tab_hela_common <tab_hela_common> 	{input.tab_hela_common}
	--tab_k562_all <input.tab_k562_all> 	{input.tab_k562_all}
	--tab_hela_all <input.tab_hela_all> 	{input.tab_hela_all}
	--output <output> 	{output.pdf}
" -> doc

## -----------------------------------------------------------------------------
## Libraries
## -----------------------------------------------------------------------------
# http://r.789695.n4.nabble.com/Install-package-automatically-if-not-there-td2267532.html

load.fun <- function(x) {
        x <- as.character(substitute(x))
        if(isTRUE(x %in% .packages(all.available=TRUE))) {
                suppressWarnings(eval(parse(text=paste("require(", x, ")", sep=""))))
        } else {
                ## Modified part and not really tested...
                # Should add here the part for cran packages    
                ## maybe split function in two part, one with cran install method and the other with bioclite method.
                ###Bioconductor part:
                ## try http:// if https:// URLs are not supported
                #source("http://bioconductor.org/biocLite.R")
                #eval(parse(text=paste("biocLite('", x, "')", sep=""))) 
        }
}

## -----------------------------------------------------------------------------
## Arg parsing
## -----------------------------------------------------------------------------
suppressMessages(load.fun("docopt"))
my_opts <- docopt(doc)

# A control inspired by Denis code but here we have 5 arguments and length=15 so this is not the correct way to do it.
#if(length(commandArgs()) != 5 ){
#       print(commandArgs())
#       cat("Use -h for more informations\n")
#       q(status=1)
#}


## -----------------------------------------------------------------------------
## Parsing args
## -----------------------------------------------------------------------------

tab_k562_k562 <- my_opts$"--tab_k562_k562"
tab_k562_hela <- my_opts$"--tab_k562_hela"
tab_k562_common <- my_opts$"--tab_k562_common"
tab_hela_k562 <- my_opts$"--tab_hela_k562"
tab_hela_hela <- my_opts$"--tab_hela_hela"
tab_hela_common <- my_opts$"--tab_hela_common"
tab_k562_all <- my_opts$"--tab_k562_all"
tab_hela_all <- my_opts$"--tab_hela_all"
output <- my_opts$"--output"

## -----------------------------------------------------------------------------
## Loading libs
## -----------------------------------------------------------------------------
#suppressMessages(load.fun(biomaRt))

## -----------------------------------------------------------------------------
## Defining function
## -----------------------------------------------------------------------------

print("Loading data")
d_tab_k562_k562 <- read.table(file=tab_k562_k562, header=TRUE)
d_tab_k562_hela <- read.table(file=tab_k562_hela, header=TRUE)
d_tab_k562_common <- read.table(file=tab_k562_common, header=TRUE)
d_tab_hela_k562 <- read.table(file=tab_hela_k562, header=TRUE)
d_tab_hela_hela <- read.table(file=tab_hela_hela, header=TRUE)
d_tab_hela_common <- read.table(file=tab_hela_common, header=TRUE)
d_tab_k562_all <- read.table(file=tab_k562_all, header=TRUE)
d_tab_hela_all <- read.table(file=tab_hela_all, header=TRUE)

med_k562 <- median(d_tab_k562_all$transcript_cov)
med_hela <- median(d_tab_hela_all$transcript_cov)

list_norm_by_med <- list(
	k562_k562=d_tab_k562_k562$transcript_cov - med_k562,
	k562_common=d_tab_k562_common$transcript_cov - med_k562,
	k562_hela=d_tab_k562_hela$transcript_cov - med_k562,
	hela_hela=d_tab_hela_hela$transcript_cov - med_hela,
	hela_common=d_tab_hela_common$transcript_cov - med_hela,
	hela_k562=d_tab_hela_k562$transcript_cov - med_hela)

list_norm_by_med
print("first word after dollar is the signal strain, second word is class promoter strain")
wilcox.test(list_norm_by_med$k562_k562, list_norm_by_med$k562_hela, alternative="greater")
wilcox.test(list_norm_by_med$k562_common, list_norm_by_med$k562_hela, alternative="greater")
wilcox.test(list_norm_by_med$hela_hela, list_norm_by_med$hela_k562, alternative="greater")
wilcox.test(list_norm_by_med$hela_common, list_norm_by_med$hela_k562, alternative="greater")

print("Do not forget we did median-normalization and that greatly affect conclusion we can have with the following comparisons...")
wilcox.test(list_norm_by_med$k562_k562, list_norm_by_med$hela_k562, alternative="greater")
wilcox.test(list_norm_by_med$hela_hela, list_norm_by_med$k562_hela, alternative="greater")

print("Paired-test")
wilcox.test(list_norm_by_med$k562_k562, list_norm_by_med$hela_k562, alternative="greater", paired=TRUE)
wilcox.test(list_norm_by_med$hela_hela, list_norm_by_med$k562_hela, alternative="greater", paired=TRUE)



#pdf(output)
## First test:
#boxplot(
#	d_tab_k562_k562$transcript_cov,
#	d_tab_hela_k562$transcript_cov,
#	d_tab_k562_hela$transcript_cov,
#	d_tab_hela_hela$transcript_cov,
#	d_tab_k562_common$transcript_cov,
#	d_tab_hela_common$transcript_cov,
#	notch=TRUE,ylab="log2(H3K27ac/H3K4me3)",
#	cex.axis=0.8,
#	names=c(
#		"K562 signal\nin K562\nspecific\nEProm",
#		"Hela signal\nin K562\nspecific\nEProm",
#		"K562 signal\nin Hela\nspecific\nEProm",
#		"Hela signal\nin Hela\nspecific\nEProm",
#		"K562 signal\nin common\nEprom",
#		"Hela signal\nin common\nEProm"),
#	las=2
#	)
#dev.off()
#
#
## Second test
#pdf(paste0(output,"2.pdf"))
#
#par(mfrow=c(1,2))
#
#boxplot(
#	d_tab_k562_k562$transcript_cov,
#	d_tab_k562_hela$transcript_cov,
#	d_tab_k562_common$transcript_cov,
#	notch=TRUE,
#	#ylab="log2(H3K27ac/H3K4me3)",
#	main="K562 signal",
#	cex.axis=0.8,
#	names=c(
#		"K562\nspecific\nEProm",
#		"Hela\nspecific\nEProm",
#		"common\nEprom"),
#	las=2
#	)
#
#boxplot(
#	d_tab_hela_k562$transcript_cov,
#	d_tab_hela_hela$transcript_cov,
#	d_tab_hela_common$transcript_cov,
#	notch=TRUE,ylab="log2(H3K27ac/H3K4me3)",
#	main="Hela signal",
#	cex.axis=0.8,
#	names=c(
#		"K562\nspecific\nEProm",
#		"Hela\nspecific\nEProm",
#		"common\nEProm"),
#	las=2
#	)
#
#dev.off()
#
#
#
#
## Third test: reordering the classes to keep the same reading logic between K5562 and Hela.
#
#pdf(paste0(output,"3.pdf"))
#
#par(mfrow=c(1,2))
#
#boxplot(
#	d_tab_k562_k562$transcript_cov,
#	d_tab_k562_common$transcript_cov,
#	d_tab_k562_hela$transcript_cov,
#	notch=TRUE,
#	#ylab="log2(H3K27ac/H3K4me3)",
#	main="K562 signal",
#	cex.axis=0.8,
#	names=c(
#		"K562\nspecific\nEProm",
#		"common\nEprom",
#		"Hela\nspecific\nEProm"),
#	las=2
#	)
#
#boxplot(
#	d_tab_hela_hela$transcript_cov,
#	d_tab_hela_common$transcript_cov,
#	d_tab_hela_k562$transcript_cov,
#	notch=TRUE,ylab="log2(H3K27ac/H3K4me3)",
#	main="Hela signal",
#	cex.axis=0.8,
#	names=c(
#		"Hela\nspecific\nEProm",
#		"common\nEProm",
#		"K562\nspecific\nEProm"),
#	las=2
#	)
#
#dev.off()
#
#
#
#
## Fourth test: same yscale
#pdf(paste0(output,"4.pdf"))
#
#par(mfrow=c(1,2))
#
#boxplot(
#	d_tab_k562_k562$transcript_cov,
#	d_tab_k562_common$transcript_cov,
#	d_tab_k562_hela$transcript_cov,
#	notch=TRUE,
#	#ylab="log2(H3K27ac/H3K4me3)",
#	main="K562 signal",
#	cex.axis=0.8,
#	ylim=c(-4,4),
#	names=c(
#		"K562\nspecific\nEProm",
#		"common\nEprom",
#		"Hela\nspecific\nEProm"),
#	las=2
#	)
#
#boxplot(
#	d_tab_hela_hela$transcript_cov,
#	d_tab_hela_common$transcript_cov,
#	d_tab_hela_k562$transcript_cov,
#	notch=TRUE,ylab="log2(H3K27ac/H3K4me3)",
#	main="Hela signal",
#	cex.axis=0.8,
#	ylim=c(-4,4),
#	names=c(
#		"Hela\nspecific\nEProm",
#		"common\nEProm",
#		"K562\nspecific\nEProm"),
#	las=2
#	)
#
#dev.off()
#
#
## Fifth test: same yscale, normalized by median
#pdf(paste0(output,"5.pdf"))
#
#par(mfrow=c(1,2))
#
#boxplot(d_norm_by_med[,c(1:3)],
#	notch=TRUE,
#	#ylab="log2(H3K27ac/H3K4me3)",
#	main="K562 signal",
#	cex.axis=0.8,
#	ylim=c(-4,6),
#	names=c(
#		"K562\nspecific\nEProm",
#		"common\nEprom",
#		"Hela\nspecific\nEProm"),
#	las=2
#	)
#
#boxplot(d_norm_by_med[,c(4:6)],
#	notch=TRUE,ylab="Median-normalized log2(H3K27ac/H3K4me3)",
#	main="Hela signal",
#	cex.axis=0.8,
#	ylim=c(-4,6),
#	names=c(
#		"Hela\nspecific\nEProm",
#		"common\nEProm",
#		"K562\nspecific\nEProm"),
#	las=2
#	)
#dev.off()
#

# Sixth test: same yscale, normalized by median, improved aspect
pdf(output)

op <- par(mfrow=c(1,2),
	oma = c(4.5,4,0,0) + 0.1,
	mar = c(0,0,1.5,1) + 0.1)

#par(op)

#boxplot(d_norm_by_med[,c(1:3)],
boxplot(list_norm_by_med$k562_k562,list_norm_by_med$k562_common,list_norm_by_med$k562_hela,
	notch=TRUE,
	#ylab="log2(H3K27ac/H3K4me3)",
	main="K562 signal",
	#cex.axis=0.8,
	ylim=c(-4,6),
	names=c(
		"K562\nspecific\nEProm",
		"common\nEprom",
		"Hela\nspecific\nEProm"),
	las=2
	)

#boxplot(d_norm_by_med[,c(4:6)],
boxplot(list_norm_by_med$hela_hela,list_norm_by_med$hela_common,list_norm_by_med$hela_k562,
	notch=TRUE,ylab="Median-normalized log2(H3K27ac/H3K4me3)",
	main="Hela signal",
	#cex.axis=0.8,
	axes=FALSE,
	ylim=c(-4,6),
	names=c(
		"Hela\nspecific\nEProm",
		"common\nEProm",
		"K562\nspecific\nEProm"),
	las=2
	)
axis(side=1, at=c(1,2,3), las=2, labels=c(
	"Hela\nspecific\nEProm",
	"common\nEProm",
	"K562\nspecific\nEProm"))

axis(side=2,labels=FALSE)
box()

title(ylab = "Median-normalized log2(H3K27ac/H3K4me3)",
      outer = TRUE, line = 3)

dev.off()


