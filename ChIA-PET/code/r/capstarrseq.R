# est ce que les enhProm interagisse plus souvent avec des promoteurs de gènes.
# est ce que les enhProm interagisse ont plus d'interaction que les autres promoteurs.
# A quelle distance sont localisée ces interactions.

# Load parameters

nb_sd <-  2
############################################
# Load data about promoter enhancer activity
############################################
setwd("PROGS/starrseq_r_code/")
  all <- read.table("../../DATA/allResults_FC_unix.txt", head=T)
  
  # Select pics located in promoter regions:
  all.prox <- subset(all, location =="P")
  
  # Select pics called with inflexion point
  
  all.prox.sign <- subset(all.prox, K562_X_replicates_group_inflexionPoint == "Active")
  
  # load control regions
  
  ctrl <- read.table("../../DATA/ControlRegions_FC0.8-0.9_All_samples.txt", head=F)
  
  
  # Select the same number of ctrl and sign.
  
  nb.sign <- nrow(all.prox.sign)
  ctrl.selected.pos <- sample(1:nrow(ctrl), size = nb.sign)
  ctrl.select <- ctrl[ctrl.selected.pos,]
  ctrl.select.cordinate <- paste(paste(ctrl.select[,1], ctrl.select[,2], sep=":"), ctrl.select[,3], sep="-")
  pos <- all$cordinate %in% ctrl.select.cordinate
  
  ctrl.select.full <- all[pos,]
  
  write.table(ctrl.select.full, "../../RESULTS/CONTROL/control.tsv", col.names=NA, sep="\t", quote=F)
  write.table(ctrl.select.full, "../../RESULTS/ENH_PROM/enh_prom_k562.tsv", col.names=NA, sep="\t", quote=F)
  
  # now we should work with all.prox.sign (Promoter with enhancer activity) and ctrl.select.full (Promoter without enhancer activity).

############################################
# Load data about long distance 
# interaction of promoters
############################################
# first replicate
# ppi for promoter-promoter interactions

ppi <- read.table("../../DATA/ENCFF002ENM_B2_T2_SNYD_POLR2A.tsv", header=TRUE)
br <- seq (10,20, len=100)
hist(log2(ppi$stop1-ppi$start1), main="toot")
summary(ppi$stop1-ppi$start1)
ppi$  
# Compute a threshold based on left and right anchor length.

tresh <- mad(log2(c(ppi.1$left_anchor_len, ppi.1$right_anchor_len))) * nb_sd
med <- median(log2(c(ppi.1$left_anchor_len, ppi.1$right_anchor_len)))
up.bound <- med+tresh
low.bound <- med -tresh

anchor.len.sort  <- sort(c(ppi.1$left_anchor_len, ppi.1$right_anchor_len))
indexes <- 1:(nrow(ppi.1)*2)


plot(indexes, log2(anchor.len.sort) , pch=".", ylim=c(8,20))
abline(h=med, col="green")
abline(h=up.bound, col="blue")
abline(h=low.bound)

pos.deleted <- log2(anchor.len.sort) > (med + tresh) | log2(anchor.len.sort) < (med - tresh)
points(indexes[pos.deleted],  log2(anchor.len.sort)[pos.deleted] , col="red", pch=".")

print(paste("Selected range for pic length : [", round(2^(low.bound)), ",",  round(2^(up.bound)), "]", sep=""))

# Now apply this threshold and select only fragment with left and right pair

left.is.ok <- log2(ppi.1$left_anchor_len) < (med + tresh) & log2(ppi.1$left_anchor_len) > (med - tresh)
right.is.ok <- log2(ppi.1$right_anchor_len) < (med + tresh) & log2(ppi.1$right_anchor_len) > (med - tresh)
ppi.1.sel <- ppi.1[left.is.ok & left.is.ok, ]
dim(ppi.1.sel)

dim(ppi.1)

file.show("../../DATA/table_s3g_K562_saturated_IXN_1_sel.txt")

print(paste("% of feature kept in IXN_1 (both anchor): ", nrow(ppi.1.sel)/nrow(ppi.1)*100 ))


################################################
# Select by nb_count (anchor)
################################################
hist(log2(ppi.1.sel$nb_count+1), br=100)
summary(ppi.1.sel$nb_count)
ppi <- ppi[ppi.1.sel$nb_count > ]
