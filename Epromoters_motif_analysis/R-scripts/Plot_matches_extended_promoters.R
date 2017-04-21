#############################
## How to run this script:
## # cat Plot_matches_extended_promoters.R | /usr/bin/R --slave --no-save --no-restore --no-environ --args " prefix = 'testing'; parsed.inactive.file = 'CapStarrseq_Active_Prom_HELA_merge_IP_extended_matrix_scan_pval_1e-3_Human_bg_mkv_2_PARSED.tab'; parsed.active.file = 'CapStarrseq_Active_Prom_K562_merge_IP_extended_matrix_scan_pval_1e-3_Human_bg_mkv_2_PARSED.tab';"


## Required libraries
suppressPackageStartupMessages(library("IRanges", warn.conflicts=FALSE, character.only = TRUE))
suppressPackageStartupMessages(library("RColorBrewer", warn.conflicts=FALSE, character.only = TRUE))
suppressPackageStartupMessages(library("gplots", warn.conflicts=FALSE, character.only = TRUE))
suppressPackageStartupMessages(library("jpeg", warn.conflicts=FALSE, character.only = TRUE))


###########################################
## Read arguments from the command line.
##
## Arguments passed on the command line
## will over-write the default arguments
## specified above.
message("Reading arguments from command-line")
args <- commandArgs(trailingOnly=TRUE);
if (length(args >= 1)) {
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

message("Checking mandatory arguments")
if (!exists("parsed.active.file")) {
  stop("Missing mandatory argument (Active promoters set): parsed.active.file ")
} else if (!exists("parsed.inactive.file")) {
  stop("Missing mandatory argument (Inactive promoters set): parsed.inactive.file ")
} else if (!exists("prefix")) {
  stop("Missing mandatory argument (prefix): prefix ")
} 

if (!exists("p.val")) {
  p.val <- 1e-3
} else if (!exists("bin")) {
  bin <- 50
}


# cat Plot_matches_extended_promoters.R | /usr/bin/R --slave --no-save --no-restore --no-environ --args " prefix = 'testing'; parsed.inactive.file = 'CapStarrseq_Active_Prom_HELA_merge_IP_extended_matrix_scan_pval_1e-3_Human_bg_mkv_2_PARSED.tab'; parsed.active.file = 'CapStarrseq_Active_Prom_K562_merge_IP_extended_matrix_scan_pval_1e-3_Human_bg_mkv_2_PARSED.tab';"

# prefix <- "/home/jaimicore/Documents/PhD/Human_promoters_project/bin/enrichment_by_scan/Extended_promoters_analysis"
# parsed.inactive.file <- "/home/jaimicore/Documents/PhD/Human_promoters_project/bin/enrichment_by_scan/parsed_files_bk/CapStarrseq_InactiveProm_FDR95_All_samples_extended_matrix_scan_pval_1e-3_Human_bg_mkv_2_PARSED.tab"
# parsed.active.file <- "/home/jaimicore/Documents/PhD/Human_promoters_project/bin/enrichment_by_scan/parsed_files_bk/CapStarrseq_Active_Prom_K562_merge_IP_extended_matrix_scan_pval_1e-3_Human_bg_mkv_2_PARSED.tab"

# parsed.active.file <- "/home/jaimicore/Documents/PhD/Human_promoters_project/bin/enrichment_by_scan/parsed_files/CapStarrseq_Active_Prom_K562_merge_IP_extended_matrix_scan_pval_1e-3_Human_bg_mkv_2_PARSED.tab"
# grep -v '^;' CapStarrseq_InactiveProm_FDR95_All_samples_extended_matrix_scan_pval_1e-3_Human_bg_mkv_2.tab | awk -F '\t'  ' $2!="limit" && ($11 >= 4) {print $1"\t"$3"\t"($6+$5)/2"\t"$9} ' | more

##################################################################
## Parse the matrix-scan input
## Create a temporary file with those columns that will be used 
## to plot the TFBSs.
## This is done because reduce the time to read the files.
# message("Parsing matrix-scan tables")
# 
# parsed.active.file <- paste(prefix, "_active_set_parsed.tab", sep = "")
# parsed.inactive.file <- paste(prefix, "_inactive_set_parsed.tab", sep = "")
#   
# cat.file <- paste(prefix, "_temporal_cat_commands.txt", sep = "")
# cat(paste("grep -v '^;' ", matrix.scan.active, " | awk -F '\\t' '$2!=\"limit\" && ($11 >= 4) {print $1\"\\t\"$3\"\\t\"($6+$5)/2\"\\t\"$9} ' > ", parsed.active.file, " ;\n", sep = ""), file = cat.file, append = TRUE)
# cat(paste("grep -v '^;' ", matrix.scan.inactive, " | awk -F '\\t' '$2!=\"limit\" && ($11 >= 4) {print $1\"\\t\"$3\"\\t\"($6+$5)/2\"\\t\"$9} ' > ", parsed.inactive.file, " ;\n", sep = ""), file = cat.file, append = TRUE)

# system(paste("bash", cat.file))

# parsed.active.file <- paste(prefix, "_active_set_parsed.tab", sep = "")
# parsed.inactive.file <- paste(prefix, "_inactive_set_parsed.tab", sep = "")

# parsed.inactive.file = 'parsed_files/CapStarrseq_Active_Prom_HELA_merge_IP_extended_matrix_scan_pval_1e-3_Human_bg_mkv_2_PARSED.tab'
# parsed.active.file = 'parsed_files/CapStarrseq_Active_Prom_K562_merge_IP_extended_matrix_scan_pval_1e-3_Human_bg_mkv_2_PARSED.tab'

#############################################
## Read matrix-scan table Active Promoters
message("Reading active promoters data")
scan.results.active <- read.csv(file = parsed.active.file, sep = "\t", header = TRUE, comment.char = ";")
colnames(scan.results.active) <- c("seq_id", "ft_name", "bspos", "Pval")


#############################################
## Read matrix-scan table Inactive Promoters
message("Reading inactive promoters data")
scan.results.inactive <- read.csv(file = parsed.inactive.file, sep = "\t", header = TRUE, comment.char = ";")
colnames(scan.results.inactive) <- c("seq_id", "ft_name", "bspos", "Pval")

# save(file = "Actives.RData", scan.results.active)
# save(scan.results.inactive, file = "Inactives.RData")

# copy.active <- scan.results.active
# copy.inactive <- scan.results.inactive 

#################
## Set p-value
# p.val <- 1e-4
p.val <- as.numeric(p.val)

#############################
## Get the matrices name's
matrix.names.active <- unique(as.vector(scan.results.active$ft_name))
matrix.names.inactive <- unique(as.vector(scan.results.inactive$ft_name))

#############################
## Get the sequences name's
seq.id.active <- unique(as.vector(scan.results.active$seq_id))
seq.id.inactive <- unique(as.vector(scan.results.inactive$seq_id))

bin <- 50
message(paste("Setting bins of size", bin))
bin <- as.numeric(bin)

windows <- IRanges(start = seq(from = -500, to = 300 - bin + 1, by = bin), width = bin)

matrix.names.active <- unique(as.vector(scan.results.active$ft_name))

profiles <- NULL
pvalues <- NULL
windows.labels <- NULL
mean.fold.change <- NULL
median.fold.change <- NULL

message("Print the profiles in a PDF file")
pdf.file.name <- paste(prefix, "_TFBSs_distribution_in_extended_promoters_scaled_TEST_BIASED_NORMALIZATION.pdf", sep = "")
pdf(pdf.file.name)

#################################################
## Normalize the nmber of hits before plotting ##
#################################################
## c(96, 122) # Interesting cases
sapply(1:length(matrix.names.active), function(m){
# sapply(c(122,165,171,199), function(m){


  print(m)
  
  matrix.query <- matrix.names.active[m]
  matrix.query.selection.active <- scan.results.active[scan.results.active$ft_name == matrix.query,]
  matrix.query.selection.inactive <- scan.results.inactive[scan.results.inactive$ft_name == matrix.query,]
  
  matrix.query.selection.active$bspos <- matrix.query.selection.active$bspos + 1000
  matrix.query.selection.inactive$bspos <- matrix.query.selection.inactive$bspos + 1000
  
  
  active.IR <- IRanges(start = matrix.query.selection.active$bspos, end = matrix.query.selection.active$bspos)
  inactive.IR <- IRanges(start = matrix.query.selection.inactive$bspos, end = matrix.query.selection.inactive$bspos)
  
  counts.per.range.active <- countOverlaps(windows, active.IR)
  counts.per.range.inactive <- countOverlaps(windows, inactive.IR)
  
  
  ## Chi-square test to calculate a p-value
  ## H0 = the number of hits (TFBSs) follow the same  distribution in the
  ## Active and Inactive promoters
  ## The test is restricted to the core promoter region (-200 , +50)
  counts.per.range.matrix <- matrix(c(counts.per.range.active[6:11], counts.per.range.inactive[6:11]), ncol = 2)

  chi <- chisq.test(counts.per.range.active[6:11], p = counts.per.range.inactive[6:11]/sum(counts.per.range.inactive[6:11]))

  TFBSs.enrichment.pval <- round(as.numeric(chi[[3]]), digits = 20)
#   TFBSs.enrichment.eval <- TFBSs.enrichment.pval * length(matrix.names.active)
  TFBSs.enrichment.eval <- TFBSs.enrichment.pval
  TFBSs.enrichment.eval <- prettyNum(TFBSs.enrichment.eval, scientific=TRUE, digits = 3)
  
  ## Expected counts
#   counts.per.range.inactive <- round(counts.per.range.inactive * (sum(counts.per.range.active)/sum(counts.per.range.inactive)))
#   profiles <<- cbind(profiles,(counts.per.range.active - counts.per.range.inactive))
  profiles <<- cbind(profiles,(counts.per.range.active - round(counts.per.range.inactive * (sum(counts.per.range.active)/sum(counts.per.range.inactive)))))

  counts.per.range.active <- counts.per.range.active/sum(counts.per.range.active)
  counts.per.range.inactive <- counts.per.range.inactive * (sum(counts.per.range.active)/sum(counts.per.range.inactive))
#   counts.per.range.inactive <- counts.per.range.inactive/sum(counts.per.range.inactive[16:21])

  fold.change <- log2(counts.per.range.active/counts.per.range.inactive)
  mean.fold.change <<- cbind(mean.fold.change, mean(fold.change[6:11]))
  median.fold.change <<- cbind(median.fold.change, median(fold.change[6:11]))

  total.hits.active <- sum(counts.per.range.active)
  total.hits.inactive <- sum(counts.per.range.inactive)
  
  windows.df <- data.frame(windows)

  names(counts.per.range.active) <- windows.df$start
  names(counts.per.range.active)[names(counts.per.range.active) == "0"] <- "1"

  names(counts.per.range.inactive) <- windows.df$start
  names(counts.per.range.inactive)[names(counts.per.range.inactive) == "0"] <- "1"

  y.val.active <- c(0, counts.per.range.active, 0)
  x.val.active <- as.numeric(c(-500, names(counts.per.range.active), 300))

  y.val.inactive <- c(0, counts.per.range.inactive, 0)
  x.val.inactive <- as.numeric(c(-500, names(counts.per.range.inactive), 300))

  pvalues <<- append(pvalues, TFBSs.enrichment.eval)
  windows.labels <<- names(counts.per.range.active)
  
  plot(x = c(-250,-250,50,50),
       y = c(-0.1,-0.4,-0.4,-0.1),
       type = "l",
#        ylim = c(0, max(c(counts.per.range.inactive, counts.per.range.inactive))+50),
        ylim = c(0, 0.25),
       xlim = c(-500,300),
       col = "#ffeda0",
       lwd = 1,
       ## Labels
       main = paste("Motif:", matrix.query),
       xlab = "Distance to TSS (bp)",
       ylab = "Normalized Nb Hits",
       ## Hide x-axis
       xaxt='n', 
  )  
#   polygon(x = c(-250,-250, 50, 50), y = c(0, 1000, 1000, 0), col="#ffeda0", border = NA, lty = 0, )
polygon(x = c(-250,-250, 50, 50), y = c(0, 1, 1, 0), col="#ffeda0", border = NA, lty = 0, )
    
  ## Draw the grid
  abline(v=(x.val.active), col="lightgray", lty="dotted")
  abline(h=(seq(from = 0, to = 100, by = 5)), col="lightgray", lty="dotted")

  ## Draw the TSS (+1) position
  abline(v = 1, col="#045a8d", lwd = 2, lty = 2)
  
  ## Set x-axis values 
  axis(side = c(1,2,3,4), at = as.character(x.val.active), labels = as.character(x.val.active))

  ## Draw the lines for the active promoters 
  lines(x = x.val.active, y = y.val.active, type = "l", col = "#00BFC4", lty = 1, lwd = 3)
  ## "#006837"

  ## Draw the lines for the inactive promoters 
  lines(x = x.val.inactive, y = y.val.inactive, type = "l", col = "#F8766D", lty = 1, lwd = 3)
  ## "#bd0026"

  ## Draw the legend
  legend("topleft", legend = c("Epromoters", "Non-Epromoters", "Core Promoter", "TSS"), fill = c("#00BFC4", "#F8766D", "#ffeda0", "#045a8d"), bty="o", bg="white")
  legend("topright", legend = paste("E-value:", TFBSs.enrichment.eval), bty="o", bg="white")

  # logo.file <- paste("logos/", matrix.query, "_logo.jpeg", sep = "")
  # logo <- readJPEG(logo.file)
  # rasterImage(logo, 
  #           xleft = 60,
  #           xright = 325, 
  #           ybottom = 0.22,
  #           ytop = 0.235)
})
dev.off()

# save(mean.fold.change, file="mean_FC.RData")
# save(median.fold.change, file="median_FC.RData")


## To export the table only
merged.TFs.names <- read.csv("Root_motifs_names.tab", sep = "\t", header = FALSE)
vector <- sapply(matrix.names.active, function(n){ which(n == as.vector(merged.TFs.names$V1)) })
TF.names <- as.vector(merged.TFs.names[as.vector(vector),2])

sorted.table <- data.frame(as.vector(mean.fold.change), as.numeric(pvalues), matrix.names.active, TF.names) 
sorted.table <- sorted.table[rev(order(sorted.table[,1])),]
colnames(sorted.table) <- c("Mean_Fold_Change", "E-value", "Merge", "TFs")
write.table(sorted.table, file = paste(prefix, "_enriched_TFs.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)

ordered.merge <- as.vector(sorted.table$Merge)

##########################
## Draw Profiles heatmap

# a <- profiles

# profiles <- data.frame(t(profiles))
# rownames(profiles) <- matrix.names.active
# colnames(profiles) <- as.character(data.frame(windows)$start)
# profiles <- as.matrix(profiles)


# save(profiles, paste(prefix, "_profiles.RData", sep = ""))

# ## Color palette
# rgb.palette <- colorRampPalette(brewer.pal(7, "RdBu"), space="Lab")
# 
# ## Color code to represent the core promoter 
# ## in the heatmap
# windows.bin <- as.numeric(data.frame(windows)$start)
# windows.bin.colors <- sapply(bin, function(b){
#   if(bin < -250){
#     return("#fdd49e")
#   } else if (bin >= 100){
#     return("#fdd49e")
#   } else {
#     return("#b30000")
#   }
# })
# 
# # profiles.k562 <- profiles
# 
# 
# 
# ## Heatmap
# pdf("heatmap_profiles.pdf")
# heatmap.2(
#   
#   profiles[ordered.merge,],
# #   profiles[,order(profiles[18,], profiles[19,])][,11:31],
#   
#           ## Dendrogram control
#           dendrogram = c("none"),
#           Rowv = FALSE,
#           Colv = FALSE,
#           
# #            hclustfun = function(d){hclust(d, method="ward")},
#           
#           ## Color
#           col = rgb.palette,
#           
#           ## Trace
#           trace = "none",
#           
#           ## Key control
#           key = TRUE,
#           keysize = 2,
#           density.info = "none",
#           key.xlab = "",
#           key.ylab = "",
#           key.title = "",
#           
#           offsetCol = 0.25,
#           cexRow = 0.25,
#           
# #           ColSideColors = windows.bin.colors,
#         )
# dev.off()

#cat /home/jcastro/Documents/JaimeCastro/PhD/Human_promoters_project/bin/enrichment_by_scan/Plot_matches_extended_promoters.R | /usr/bin/R --slave --no-save --no-restore --no-environ --args " matrix.scan.active = '/home/jcastro/Documents/JaimeCastro/PhD/Human_promoters_project/test_metrics_with_yeast_data/CapStarrseq_Active_Prom_K562_merge_IP_extended_matrix_scan_pval_1e-3_HOCOMOCO_bg_mkv_2.tab'; matrix.scan.inactive = '/home/jcastro/Documents/JaimeCastro/PhD/Human_promoters_project/test_metrics_with_yeast_data/CapStarrseq_Active_Prom_K562_merge_IP_extended_matrix_scan_pval_1e-3_HOCOMOCO_bg_mkv_2.tab'; p.val = '1e-4'; bin = '50'; pdf.file = './test.pdf'"



# a <- 1:20
# b <<- NULL
# for (i in 1:100){
#   b <<- append(b, print(sample(a)[1]))
# }
# c <- heatmap(matrix(b, ncol = 10))
# 
# d <- dist(matrix(b, ncol = 10))
# e <- hclust(d)
# plot(e)
