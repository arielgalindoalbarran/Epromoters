################
## Recover all genes in CaptStarrSeq relation

## 1) Command line retrieve of genes in CapStarrSeq that came from salvatores annotation
##

###
#library
library(GenomicRanges)
library(biomaRt)
library(org.Hs.eg.db)


# use sql to get alias table and gene_info table (contains the symbols)
# first open the database connection
dbCon <- org.Hs.eg_dbconn()
# write your SQL query
sqlQuery <- 'SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'
# execute the query on the database
aliasSymbol <- dbGetQuery(dbCon, sqlQuery)

version <- "20160422_extended_regiones" 

main.dir <- "./"

gtex.version <- "v6"
eqtl.data <- paste(data.dir,paste("eQTL/gtex",gtex.version,sep="_"), sep="/")

r.scripts <- paste(main.dir, "R-scripts", sep="/")


cap.starr.seq.folder <- paste(main.dir,"CapStarrSeq", sep="/")

## All genes associated by Salvatore column 8
## beware of SEPT genes that are edited as sept- moth date mars-o to MARS2
all.gene.names.file <- paste(cap.starr.seq.folder,"results/capStarrSeq/eQTL_compare/20160422_extended_regiones/All_genes.txt", sep="/" )
print (all.gene.names.file)
all.gene.names <- read.table(file=all.gene.names.file, sep="\t", stringsAsFactors=FALSE,header=FALSE)[,1]
## Remove gene names with conflicts 
all.gene.names<-all.gene.names[-which(all.gene.names %in% c("NAA38", "MYDGF", "GCOM1"))]

## start mart
archive.ensembl <- "feb2014.archive.ensembl.org"
##
##archive.ensembl <- "dec2013.archive.ensembl.org"
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host=archive.ensembl )
#listDatasets(mart)
ensembl = useDataset("hsapiens_gene_ensembl",mart=mart)

## Retrieve as much as possible gene coordinates using ID as HGNC symbol
all.gene.coord <- getBM(attributes=c("chromosome_name", "start_position","end_position","hgnc_symbol"), filters = "hgnc_symbol", values = unique(c(all.gene.names)) , mart=ensembl )

all.gene.coord <- all.gene.coord[!grepl("HSCHR",all.gene.coord$chromosome_name),]
all.gene.coord <- all.gene.coord[!grepl("HG", all.gene.coord$chromosome_name),]
all.gene.coord <- all.gene.coord[!grepl("LRG", all.gene.coord$chromosome_name),]

print ("Gene liste cleaned")
################
## Recover coordinates for missing genes

##load missing genes fileinformation for genes with annotation incosnsitencies 
cap.starr.missing.genesfile <- paste(cap.starr.seq.folder, "data/genomes/missing_gene_ids.txt", sep="/")
cap.starr.missing.genes.table <- read.table(file=cap.starr.missing.genesfile, sep="\t", header=TRUE, stringsAsFactor=FALSE)
print (cap.starr.missing.genes.table$gene_id_cap)
## Select genes with missing data from the full list
missing.genes <- all.gene.names[which(!(all.gene.names%in%unique(all.gene.coord$hgnc_symbol)))]

## Keep only genes that do not have pre-curated information
print (grep("NAA38", missing.genes))
genes.with.curated.data <- unique(unlist(lapply(cap.starr.missing.genes.table$gene_id_cap , grep, missing.genes)))
missing.genes <- missing.genes[-genes.with.curated.data]
print (grep("NAA38", missing.genes))
print ("Missing genes selected")

one.gene.distances <- function(x){
    gene.name.aux <- x
    #gene.name.aux <-  "ACKR1"
    print(x)
    gene.name <- unique(unlist(strsplit( gene.name.aux, split=" ")))
    gene.hgnc <- aliasSymbol[which(aliasSymbol[,2] %in% gene.name ),5]
    gene.syn.hgnc <- aliasSymbol[which(aliasSymbol[,5] %in% gene.name ),2]
    all.names <- unique(c(gene.hgnc,gene.syn.hgnc,gene.name))
    all.names <-sub ("\'","",all.names ) ## remove weird characters

    gene.coord <- getBM(attributes=c("chromosome_name", "start_position","end_position","hgnc_symbol"), filters = c("hgnc_symbol"), values = all.names  , mart=ensembl )

    print (all.names)
    print (gene.coord)

    
    if (all.names%in%c("NAA38", "MYDGF") || gene.coord$hgnc_symbol%in%c("NAA38","MYDGF")   ){
        print (all.names)
        print (gene.coord)
    }
    ## correct annotation modifications that creates problems
    if (gene.name.aux =="LSM8"){
        gene.coord$hgnc_symbol<-"LSM8"		  
    }
    
    ## print ( gene.name.aux )
    
    ## remove none canonical chromosomes
    gene.coord <- gene.coord[!grepl("HSCHR",gene.coord$chromosome_name),]
    gene.coord <- gene.coord[!grepl("HG",gene.coord$chromosome_name),]
    gene.coord <- gene.coord[!grepl("LRG",gene.coord$chromosome_name),]
    
    if (length(unique(gene.coord$chromosome_name))>1){
        print (gene.name)
        print (gene.coord)
    }
        
    if (dim(gene.coord)[1]==0){
        return(data.frame(chromosome_name="<NA>",start_position="<NA>",end_position="<NA>",hgnc_symbol="<NA>", gene_name=gene.name.aux  ))
    }else {
        gene.coord$gene_name <-   gene.name.aux 
        return(as.data.frame(gene.coord))
    }
}

missing.genes.information <- lapply(missing.genes,one.gene.distances )
missing.genes.information.df <- do.call("rbind", missing.genes.information)

## Line to check if I lost genes.

if ( length( missing.genes.information.df[missing.genes.information.df$chromosome_name=="<NA>","gene_name"])>=1){
    print ("Some genes did not have any information")
}


## Now put all genes together in one dataframe to be used in an other program for quick access to genes.

## Get gene information for the ones that required manual annotation
cap.starr.missing.genes.t.p <- data.frame(chromosome_name=cap.starr.missing.genes.table$chromosome_name,start_position=cap.starr.missing.genes.table$start_position, end_position=cap.starr.missing.genes.table$end_position, hgnc_symbol=cap.starr.missing.genes.table$gene_id_cap, gene_name=cap.starr.missing.genes.table$gene_id_cap)


## Prepare files to be rbound 
all.gene.coord.comp <-  all.gene.coord
all.gene.coord.comp$gene_name <- all.gene.coord.comp$hgnc_symbol

all.gene.coord.comp <- rbind(all.gene.coord.comp,cap.starr.missing.genes.t.p, missing.genes.information.df )

## save data in text file
all.gene.coord.file <- paste(cap.starr.seq.folder,"results",version,"AllGenes_transcript_coordinates.txt", sep="/" )

print (all.gene.coord.file)
write.table(all.gene.coord.comp, file=all.gene.coord.file, quote=FALSE, row.names=FALSE, sep="\t")

