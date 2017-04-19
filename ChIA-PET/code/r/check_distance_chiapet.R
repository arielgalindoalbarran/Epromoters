
# Draft code run in DATA directory

file_list <- list.files(pattern="ENC.*tsv")

for (file_name in file_list){
	d <- read.table(file=file_name, header=TRUE)
	d$midpoint_distance <- abs((d$start1+d$stop1)/2-(d$start2+d$stop2)/2)
	pdf(paste0("hist/",file_name,".pdf"))
	hist(d$midpoint_distance, labels=TRUE, breaks=seq(0,225545000,by=5000), xlim=c(0,50000))
	dev.off()
	}

d <- read.table("Hela_TableS2_Kuznetsova_GB2016_POLII.csv", header=TRUE)
d$midpoint_distance <- abs((d$Head_start+d$Head_end)/2-(d$Tail_start+d$Tail_end)/2)
pdf(paste0("hist/Hela_TableS2_Kuznetsova_GB2016_POLII.csv.pdf"))
hist(d$midpoint_distance, labels=TRUE, breaks=seq(0,225545000,by=5000), xlim=c(0,50000))
dev.off()


d <- read.table("Hela_TANG2015_GSM1872889_PET_clusters_RNAPII.txt")
d$midpoint_distance <- abs((d$V2+d$V3)/2-(d$V5+d$V6)/2)
pdf(paste0("hist/Hela_TANG2015_GSM1872889_PET_clusters_RNAPII.txt.pdf"))
hist(d$midpoint_distance, labels=TRUE, breaks=seq(0,247910000,by=5000), xlim=c(0,50000))
dev.off()



# Draft to run in RESULTS/FILTER_ON_CHIAPET/ directory

file_list <- list.files(pattern="fltr.txt$")



#

d <- read.table("K562_MERGE_BALL_TALL_SNYD_CHIA_POLR2A_fltr_pp_chiapet_overlap_prom.txt") 
dim(d[d$chrom1!=d$chrom2,]))




