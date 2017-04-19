# Formatting Chiapet tsv files to bed12 in order to see interactions
# See bed12 reference: https://genome.ucsc.edu/FAQ/FAQformat.html#format1

# Loading tsv chiapet
input <- read.delim("DATA/ENCFF002ENS_B2_T1_SNYD_RAD21.tsv")

# Select lines with same chr in column 1 and 4 (chrom1 and chrom2)
data <- input[input$chrom1 == input$chrom2,]
# nFiltered <- dim(input)[1] - dim(data)[1]
# print(
#   paste(
#     nFiltered,'chiapet inter-chromosome interactions removed.', dim(data)[1], 'intra-chromosome interactions kept.'
#     )
# )
# [1] "22 chiapet inter-chromosome interactions removed. 11111 intra-chromosome interactions kept."

# Define the region inluding the 2 blocks
data[,'startRegion'] <- apply(X = data[,c('start1','start2')], MARGIN = 1, FUN = min)
data[,'stopRegion'] <- apply(X = data[,c('stop1','stop2')], MARGIN = 1, FUN = max)

# Compute block1 and block2 size
data$block1Size <- data$stop1 - data$start1
data$block2Size <- data$stop2 - data$start2

# Get Blocksizes
## Need to check which block is the first
data[data$startRegion == data$start1,'blockSizes'] <- paste(block1Size, block2Size, sep=",")[data$startRegion == data$start1]
data[data$startRegion == data$start2,'blockSizes'] <- paste(block2Size, block1Size, sep=",")[data$startRegion == data$start2]

# Get Blockstart
data[data$startRegion == data$start1,'blockStarts'] <- paste(0, data$start2 - data$start1, sep=",")[data$startRegion == data$start1]
data[data$startRegion == data$start2,'blockStarts'] <- paste(0, data$start1 - data$start2, sep=",")[data$startRegion == data$start2]

# Item RGB
data[,'itemRgb'] <- '255,0,0'
data[,'blockCount'] <- 2

# Creating bed12 output
output <- data[,c('chrom1','startRegion', 'stopRegion', 'pair_name', 'IAB', 'strand1', 'startRegion', 'stopRegion', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts')]

# We write the same output to two different bed files because IGV does not display interactions in the same we if there is (or is not) a "junctions" before ".bed"
write.table(x = output, file="testBed12.bed",sep="\t", col.names = F, row.names = F, quote = F)
write.table(x = output, file="testBed12_junctions.bed",sep="\t", col.names = F, row.names = F, quote = F)