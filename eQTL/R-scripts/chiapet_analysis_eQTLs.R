################
## Recover gene list and reannotate it to correct errors
options(stringsAsFactors = FALSE)
library(ggplot2)
library(plyr)
library("extrafont")
## Read in chiapet table data


main.dir <- "./"

data.dir <- file.path(main.dir,"data")
starr.dir <- file.path(data.dir,"starr_seq_promoters")
results.dir <- file.path(main.dir,"results","chia_pet20160601")

results.file <- file.path(results.dir,"chia_pet_interactions_eQTL_support.txt")

chiapet.file <- file.path(starr.dir, "ChIA-PET_All_20160420.txt")

chiapet.table <- read.table(file=chiapet.file, header=TRUE, stringsAsFactors=FALSE, sep="\t")

################
## Read in table to correct gene annotation coordiante inconsistencies 

correct.annotations.file <- paste(data.dir,"genomes","ChiaPet_annotation_errors.txt" ,sep="/")
correct.annotation <- read.table(file=correct.annotations.file , header=TRUE, stringsAsFactors=FALSE, sep="\t")


################
## Recover gene names that were modified by excel
for (error in correct.annotation$error){
    print (error)
    gene <- correct.annotation[correct.annotation$error==error, "gene_name"]
    print (gene)
    chiapet.table[chiapet.table$Promoter_1_gene_id %in% error, "Promoter_1_gene_id"] <- gene
    chiapet.table[chiapet.table$Promoter_2_gene_id %in% error, "Promoter_2_gene_id"] <- gene   
}

################
## Anotate unique pairs

## Vectors to store information abou pairs

## Unique pairs within the list
unique.pairs <- c()

## The unique pairID for each pair in the array to re-assign it to the chiapet table to reduce duplicated lines
all.pairs <- c()

## For each pair of chia-pet interactions
## 1) Report unique pair if it has not been reproted in a previous pass of the loop
## 2) If it is a new pair report it in unique pairs and in the all.pairs vectors
## 3) If it was reported before, save the previously reported id for the pair in the all.pairs vector

for (i in 1:dim(chiapet.table)[1] ){

    ## Split the chiapet pair
    gene1 <- chiapet.table[i,c("Promoter_1_gene_id")]
    gene2 <- chiapet.table[i,c("Promoter_2_gene_id")]

    ## For the first pair
    ## Initialize the vectors
    if(length(unique.pairs)<1){
        pair <-  paste( chiapet.table[i,c("Promoter_1_gene_id","Promoter_2_gene_id") ], collapse="_")
        unique.pairs <- c(unique.pairs, pair)
        all.pairs <-c(all.pairs , pair)
        next
    }

    ## Create the two possible pairs with both gene IDs
    pair1 <- paste( chiapet.table[i,c("Promoter_1_gene_id","Promoter_2_gene_id") ], collapse="_")
    pair2 <-  paste( chiapet.table[i,c("Promoter_2_gene_id","Promoter_1_gene_id") ], collapse="_")

    ## If the pair was reported before
    if ( (pair1 %in% unique.pairs) | (pair2 %in% unique.pairs)) {

        ## Add the pair that was reported to the all.pairs vector
        if (pair1 %in% unique.pairs){
            
            all.pairs <-c(all.pairs , pair1)
            
        } else if (pair2 %in% unique.pairs) {
            
            all.pairs <-c(all.pairs , pair2)

        }
        next

        ## If the pair was never reported before
        ## Report the pair in the unique pairs vector
        ## and the all.pairs vector
    } else{
       pair <-  paste( chiapet.table[i,c("Promoter_1_gene_id","Promoter_2_gene_id") ], collapse="_")
       unique.pairs <- c(unique.pairs, pair)
       all.pairs <-c(all.pairs , pair)


   }
}

## Report the unique pairs ID into the main table
## by adding as a column the all.pairs vector
chiapet.table$uniquePair <- all.pairs


################
## Clean repeted pairs from chiapet interaction table
## Remove same gene interactions

## Data.frame to store the cleaned chiapet data
clean.chiapet.table <- data.frame( Promoter_1_gene_id = character(),
                                  Promoter_1_Type = character() ,
                                  ## Cell_line_interaction = character (),
                                  Promoter_2_gene_id = character () ,
                                  Promoter_2_Type   = character (),
                                  ## direction = character () ,
                                  uniquePair = character()  )

## Column names  for data.frame to store the cleaned chiapet data
cpnames <- c( "Promoter_1_gene_id", "Promoter_1_Type",
             ##"Cell_line_interaction",
             "Promoter_2_gene_id" , "Promoter_2_Type",
             ## "Direction",
             "uniquePair")

## For each of the unique pairs
for (pair in unique.pairs){
    print (pair)
    ## Separate the genes in the pair
    genes <- unlist(strsplit(pair,"_"))

    ## Remove same gene interaction
    if (genes[1]==genes[2]){
        next
    }

    ## Select from the original table the pairs with the same uniquepair ID
    chia.pairs <- chiapet.table[chiapet.table$uniquePair %in% pair,]

    ## Empty data frame to store sorted pairs
    ## The genes will be sorted in the same order as the pair ID
    sorted.chia.pairs <- data.frame( Promoter_1_gene_id = character(),
                                    Promoter_1_Type = character() ,
                                    #Cell_line_interaction = character (),
                                    Promoter_2_gene_id = character () ,
                                    Promoter_2_Type   = character (),
                                    #direction = character () ,
                                    uniquePair = character()  )

    ## For each pair
    for (p in 1:dim(chia.pairs)[1]){
        ## print(p)

        ## pair
        cp <- chia.pairs[p,]

        ## Get the gene pair
        ## Get the order in the uniquePair ID
        gene1 <- genes[1]
        numgen1 <- which (cp %in% gene1)
        promoter1type <- cp [,numgen1+1]
        gene2 <- genes[2]
        numgen2 <- which (cp %in% gene2)
        promoter2type <- cp [,numgen2+1]

        ## Cell interaction and direction will be ignored
        ##cell.interaction <- cp[,3]
        ##direction <- paste (cp[,1] ,cp[,4] ,sep="_")

        ## Add the sorted chiapet information into the sorted table
        sorted.chia.pairs <- rbind( sorted.chia.pairs, c(gene1,promoter1type,
                                                         ## cell.interaction ,
                                                         gene2,promoter2type ,
                                                         ## direction,
                                                         cp$uniquePair ))
        colnames(sorted.chia.pairs ) <- cpnames
    }

    ## Remove duplicates from the sorted.chia.pairs
    unique.info <- unique(sorted.chia.pairs)

    colnames(unique.info) <-cpnames

    ## Add the cleaned pair set into the clean chiapet.table
    clean.chiapet.table <- rbind (clean.chiapet.table , unique.info)
}

print ("Chiapet interactions cleaned")


################################################################
## For each ChIA-PET pair recover the eQTLs that fall within the
## promoters of the pair, and that have a regulatory effect on any
## of the genes in the pair

## eQTLs from gtex were previously overlaped with the capstarr-seq promoters that weres
## extenden 1.5Kb to each side
eQTL.match.version <- "20160422_extended_regiones"

eQTL.match.dir <- paste(main.dir, "results",eQTL.match.version,sep="/")
eQTL.match <- paste(eQTL.match.dir,"/*","_eQTL_","*", sep="")

## header for eQTL overlaps table 
header <- c("capsstarr_cell_type","chr","start","end","width","strand","coordinate","category","gene_id_cap","transcript_id_cap","chr_eqtl","start_eqtl","end_eqtl","width_eqtl","strand_eqtl","snp_eqtl","gene_eqtl","beta_eqtl","t_stat_eqtl","se_eqtl","p_value_eqtl","nom_thresh_eqtl","min_p_eqtl","gene_emp_p_eqtl","k_eqtl","n_eqtl","gene_q_value_eqtl","beta_noNorm_eqtl","snp_chrom_eqtl","snp_pos_eqtl","minor_allele_samples_eqtl","minor_allele_count_eqtl","maf_eqtl","ref_factor_eqtl","ref_eqtl","alt_eqtl","snp_id_1kg_project_phaseI_v3_eqtl","rs_id_dbSNP142_GRCh37p13_eqtl","num_alt_per_site_eqtl","has_best_p_eqtl","is_chosen_snp_eqtl","gene_name_eqtl","gene_source_eqtl","gene_chr_eqtl","gene_start_eqtl","gene_stop_eqtl","orientation_eqtl","tss_position_eqtl","gene_type_eqtl","gencode_attributes_eqtl","tss_distance_eqtl","same_gene","gtex_tissue")

## Columns to be reported in the final table
selected.columns <- c("gene1","gene2",
                      ##"chiapet_cell_type",
                      "promType1", "promType2","uniquePair",
                      "capsstarr_cell_type","coordinate","gene_id_cap",
                      "chr_eqtl","start_eqtl","end_eqtl",
                      "snp_eqtl","ref_eqtl", "alt_eqtl","snp_id_1kg_project_phaseI_v3_eqtl",
                      "rs_id_dbSNP142_GRCh37p13_eqtl","tss_distance_eqtl","gene_name_eqtl",
                      "beta_eqtl","p_value_eqtl",
                      "gtex_tissue" )

## List to store pairs with their affecting eQTLs
results.list <-list() 

## For each chiapet interaction 
for (i in c(1:dim( clean.chiapet.table )[1])){
    ## i<-8
    one.chiapet <- clean.chiapet.table[i,]
    
    ## get genes on chia-pet interaction
    gene1 <- one.chiapet$Promoter_1_gene_id
    gene2 <- one.chiapet$Promoter_2_gene_id
    
    ## get their promoter types
    promType1 <- one.chiapet$Promoter_1_Type
    promType2 <- one.chiapet$Promoter_2_Type
    
    ## get the uniquePair ID
    uniquePair <- one.chiapet$uniquePair
    
    ## From all files containing eQTL and capstarr seq exptende promoters overlaps
    cell.files <- paste( eQTL.match.dir,"/*overlap*.txt",sep="" )   
    
     
    ## Grep command to get for both genes in a pair all eQTLs that could be related to them
    ## either falling whithin the promoter of the gene or having the gene as a target.
    grep.genes.command.cell <- paste("grep -E ", paste("\'",gene1,"|",gene2 ,"\'",sep=""), cell.files," | perl -pe \' s/.+CapStarrSeq_// ; s/:chr/\tchr/  \' " )
    
    gene.promoter.eQTL.cell.v <- lapply(system(  grep.genes.command.cell , intern = TRUE) , strsplit, split="\t")
    
    gene.promoter.eQTL <- do.call("rbind.data.frame",lapply(gene.promoter.eQTL.cell.v ,unlist) )
    
    ## If the pair of chiapet connected genes are reported in the capstarr overlaped by eQTL data sets report them
    all.gene.prom <- dim(gene.promoter.eQTL)[1]
    
    if( all.gene.prom >0 ){
        ## Add header and remove undesire things left by the grep command
        colnames(gene.promoter.eQTL ) <- header
        gene.promoter.eQTL$cell.type <- gsub ("_.+", "", gene.promoter.eQTL$capsstarr_cell_type)
        gene.promoter.eQTL$capsstarr_cell_type <- sub("_overlap_eQTL.+","" , gene.promoter.eQTL$capsstarr_cell_type , perl=TRUE)
        ##results$capsstarr_cell_type <- sub("_overlap_eQTL.+","" , results$capsstarr_cell_type , perl=TRUE)
        
        
        ## Select only the eQTL interactions for the genes in the pair
        ## Possible interactions, all other interactions will be discarded
        ## 1)gene1 to gene1
        ## 2)gene1 to gene2
        ## 3)gene2 to gene1
        ## 4)gene2 to gene2
        
        gene.promoter.eQTL <- gene.promoter.eQTL[ ( ( gene.promoter.eQTL$gene_id_cap == gene1 & gene.promoter.eQTL$gene_name_eqtl ==gene1 ) |
                                                       (    gene.promoter.eQTL$gene_id_cap == gene1 & gene.promoter.eQTL$gene_name_eqtl == gene2 ) |
                                                           (  gene.promoter.eQTL$gene_id_cap == gene2 & gene.promoter.eQTL$gene_name_eqtl ==gene1 ) |
                                                               (    gene.promoter.eQTL$gene_id_cap == gene2 & gene.promoter.eQTL$gene_name_eqtl == gene2 ) ), ]
        
        
        ## Test ot check marhc7 baz2b pair
        ##t <- "2_160472950_C_T_b37"
        ## gene.promoter.eQTL[gene.promoter.eQTL$snp_eqtl==t,]
        ## gene.promoter.eQTL2[gene.promoter.eQTL2$snp_eqtl==t,]
        
        ## If there were none eQTL interactiosn that affected the genes in the pair
        ## report empty values
        if (!(dim(gene.promoter.eQTL)[1]>0)){
            print (i)
            no.recovery <-  t(as.data.frame( c(gene1, gene2,
                                                promType1,promType2,uniquePair,
                                           rep("NA", length (selected.columns )-4 ) )) )
            colnames(no.recovery ) <- c(selected.columns,"ePromOvlp.eQTL")
            results.list[[paste("r",i,sep="" )]] <- no.recovery
                     
            next
        }

        ## If here were eQTLs interacting with the genes in the pair
        ## annoate them based on the type of promoter they overlaped
        ## hpromoter or epromoter
        
        print ("Found")
        print (i)

        
        gene.promoter.eQTL$gene1 <- gene1
        gene.promoter.eQTL$gene2 <- gene2

        gene.promoter.eQTL$promType1 <- promType1
        gene.promoter.eQTL$promType2 <- promType2
        
        ## gene.promoter.eQTL$chiapet_cell_type <- cell.lines

        gene.promoter.eQTL$uniquePair <- uniquePair 

        gene.promoter.eQTL <-   gene.promoter.eQTL[, selected.columns]

        ## Vector to store the types of promoters
        ePromOvlp.eQTL <-c()

        ## For each eQTL overlaping promoters 
        for (eq in 1:dim(gene.promoter.eQTL )[1]){

            one.eqtl <- gene.promoter.eQTL[eq,]

            ## Is any of the promoters an ePromoter
            eproms <- which (grepl ("EPromoter" , one.eqtl) %in% TRUE)

            ## None of the promoters where epromoter
            ## but there is an eQTL
            if (length(eproms)==0){ 
                
                ePromOvlp.eQTL <- c(ePromOvlp.eQTL, -1)                    

            ## Both promoters are gene promoters    
            }else if(length(eproms)==2){
                
                ePromOvlp.eQTL <- c(ePromOvlp.eQTL, 2) ## both genes have an Epromoter
                next

                ## The pair is assymetrical
                ## One is hpromoter and the other one is epromoter
            }else if (length(eproms)==1){ ## One of the genes was epromoter

                ## Check if the eQTL fell in the one that is an ePromoter
                ## Get the gene name that had the epromoter and check if in this line it corresponds to the promoter where the eQTL fell
                if (eproms==3){
                    geneE <- one.eqtl$gene1
                    
                }else if (eproms == 4) {
                    
                    geneE <- one.eqtl$gene2

                    ## All options are contempleted report error if the value is an other one
                }else {
                    print(eq)
                    print ("Error" )
                    stop()
                    
                }
                ## There is one ePromoter and the eQTL falls within the extention of it
                if (geneE == one.eqtl$gene_id_cap ){
                    
                    ePromOvlp.eQTL <- c(ePromOvlp.eQTL, 1)
                    
                }else {
                    ## There was en ePromoter but the eQTL fell within the hPromoter
                    ePromOvlp.eQTL <- c(ePromOvlp.eQTL, 0)
                    
                }
                
                
            }else {
                print ("Error")
            }
            
            
            
        }
        gene.promoter.eQTL$ePromOvlp.eQTL <- ePromOvlp.eQTL
        
        results.list[[paste("r",i,sep="" )]] <- gene.promoter.eQTL

        

       
        
        ## Else if there was no existance of a capstarr overlaping eQTL data report it as NA
    } else {
        print (i)
        no.recovery <-  t(as.data.frame( c(gene1, gene2,
                                           promType1,promType2,uniquePair,
                                           rep("NA", length (selected.columns )-4 ) )) )
        colnames(no.recovery ) <- c(selected.columns,"ePromOvlp.eQTL")
        results.list[[paste("r",i,sep="" )]] <- no.recovery
        
    
    }
    
    
}

## Bind results stored in the list to on data.frame
results <- do.call("rbind.data.frame", results.list )

## Print data.frame into table
write.table(results, file=results.file, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")



################
## Analyise ChiaPet pairs interactions and their promoters

## results <- read.table(file=results.file, sep="\t", stringsAsFactors=FALSE,header=TRUE)

## total counts of epromoters and promoter
all.total.number.h.promoter <- 20719 ## number stracted from the file data/genomes/All_hpromoterRegions_salvatore20160419.txt
total.capstarr.control <-12742
total.epromoter <- 981
repeated.epromoter <- 145 ## epromoters in both cell lines k562 and hela

## total number of reported epromoter and the number minus the ones that are in both cell lines
total.capstarr.k562 <- 632
norep.capstarr.k562 <- total.capstarr.k562 - repeated.epromoter
    
total.capstarr.hela <- 494
norep.capstarr.hela <- total.capstarr.hela - repeated.epromoter 

## total number of promoters that went in to the capstarr seq assay
total.num.capstarr.promoter <- total.capstarr.control + norep.capstarr.hela + norep.capstarr.k562 + repeated.epromoter 

## frequence of promoters in the epromoters set
freq.not.capstarr <-  (all.total.number.h.promoter -  total.num.capstarr.promoter)/ all.total.number.h.promoter
freq.capstarr.control <- total.capstarr.control /  all.total.number.h.promoter 
freq.capstarr.epromoter <- (norep.capstarr.k562 + norep.capstarr.hela + repeated.epromoter ) / all.total.number.h.promoter 


## Get only the pairs that have eQTL information
chia.eqtl.support <- na.omit(results)

## Cathegorize the eQTL promoter interaction into hProm and eProm
chia.eqtl.support$ePromOvlp.eQTL2 <- NA

chia.eqtl.support[chia.eqtl.support$ePromOvlp.eQTL==-1,  "ePromOvlp.eQTL2"] <- "hProm" ## None of the genes is an ePromoter
chia.eqtl.support[chia.eqtl.support$ePromOvlp.eQTL==0,  "ePromOvlp.eQTL2"] <- "hProm" ## One of the genes is an ePromoter but the eQTL is in hPromoterc
chia.eqtl.support[chia.eqtl.support$ePromOvlp.eQTL==1,  "ePromOvlp.eQTL2"] <- "eProm" ## One promoter is a ePromoter and it contains the eQTL
chia.eqtl.support[chia.eqtl.support$ePromOvlp.eQTL==2,  "ePromOvlp.eQTL2"] <- "eProm" ## Both genes in pair are ePromoters

################
## Classify interactions if the eQTL falls within the capstarrSeq

chia.eqtl.support$ePromOvlp.eQTLgene <- NA

chia.eqtl.support[chia.eqtl.support$gene_id_cap==chia.eqtl.support$gene_name_eqtl,"ePromOvlp.eQTLgene"] <- "CloseGene"  ## The eQTL falls in the CapStarrGene Close Gene
chia.eqtl.support[!(chia.eqtl.support$gene_id_cap==chia.eqtl.support$gene_name_eqtl),"ePromOvlp.eQTLgene"] <- "DistalGene"  ## The eQTL falls in the CapStarrGene Close Gene
 

##########
## Select reciprocal eQTLs
## eQTLs that affect both genes in a pair

## Emtpy data.frame to repor the eQTLs that affect both genes of a pair 
eQTLs.affect.both.gene <- data.frame("eQTL_id"=character(),
                                     "uniquePair"=character()
                                     )

## For each of the unique reported pairs
for (pair in unique.pairs ){

    ## Retrieve all promoter eQTL interactions for the pair
    chia.eqtl.support.pair <- chia.eqtl.support[chia.eqtl.support$uniquePair==pair,]

    ## Separete the gens in the pair
    genes <- unlist(strsplit(pair, "_"))

    ## Get the eQTLs affecting each of the genes in the pair
    eQTL.gene1 <- chia.eqtl.support.pair[chia.eqtl.support.pair$gene_name_eqtl==genes[1],"snp_eqtl"]
    eQTL.gene2 <- chia.eqtl.support.pair[chia.eqtl.support.pair$gene_name_eqtl==genes[2],"snp_eqtl"]

    ## Add flag, the number should always add up !!!!!

    
    ## Select all eQTLS affecting gene1 and gene2 , plus the ones affecting gene2 and gene1
    eQTL_id <- unique(c(eQTL.gene1[eQTL.gene1%in% eQTL.gene2 ], eQTL.gene2[eQTL.gene2 %in% eQTL.gene1]))

    ## if there are zero eQTLs affecting both genes the eQTL_id objet will have length zero
    ## and then the uniquePair one will have zero as well
    ## so no line will be added to the data.frame
    uniquePair <- rep(pair, length(eQTL_id))
    eQTLs.affect.both.gene <- rbind(eQTLs.affect.both.gene, cbind( eQTL_id , uniquePair ))
    
}

## retrieve all the eQTL information for the pairs that were affected by them
## reciprocally

## Save the index of the pair_eQTL information line that correspond to a
## reciprocal eQTL interaction
bieffect.index <- c()

## Check the summary of the effect for both genes
eQTLs.affect.both.genes.eqtl.summary <- data.frame(
    "num.effect_gene1"=numeric(),
    "num.effect_gene2"=numeric(),
    "shared.tissues"=numeric(),
    "total.tissues"=numeric(),
    "meanBeta_gene1"=numeric(),
    "meanBeta_gene2"=numeric(),
    "num.effect_distal"=numeric(),
    "num.effect_close"=numeric(),
    "shared.tissues.distance"=numeric(),
    "total.tissues.distance"=numeric(),
    "meanBeta_distal"=numeric(),
    "meanBeta_close"=numeric(),
    "SignDiff" = logical (),
    "PromType"=character()
    )


## for each eQTL affecting both genes of a pair
for (i in c(1:dim(eQTLs.affect.both.gene)[1])){

    ## Get the index of the  pair_eQTL information line that corresponds to the reciprocal
    ## set
    index <- which(chia.eqtl.support$snp_eqtl==eQTLs.affect.both.gene[i,"eQTL_id"] & chia.eqtl.support$uniquePair==eQTLs.affect.both.gene[i,"uniquePair"]   )

    ## Get the index
    bieffect.index <- c(bieffect.index , index)

    ## retrieve data and get summary for eqtl effects
    local.effects2 <- chia.eqtl.support[index , ]

    
    ## Type
    types <- unique(local.effects2$ePromOvlp.eQTL2)

    for (type in types){
        local.effects <- local.effects2[local.effects2$ePromOvlp.eQTL2==type,]
        ## Count for genes 1 and 2
        l.gene1 <-  unique(local.effects$gene1)
        l.gene2 <- unique(local.effects$gene2)
        
        beta.g1 <- local.effects[local.effects$gene_name_eqtl==l.gene1,"beta_eqtl"]
        beta.g2 <- local.effects[local.effects$gene_name_eqtl==l.gene2,"beta_eqtl"]
        
        num.e.g1 <- length(beta.g1)
        num.e.g2 <- length(beta.g2)
        
        t.g1 <- unique(local.effects[local.effects$gene_name_eqtl==l.gene1,"gtex_tissue"])
        t.g2 <- unique(local.effects[local.effects$gene_name_eqtl==l.gene2,"gtex_tissue"])
        
        s.t <- length(unique(t.g1[t.g1 %in% t.g2 ], t.g2[t.g2 %in% t.g1 ])) ## shared tissues
        
        t.s <- length(unique(t.g1, t.g2)) ## total tissues
        
        mean.e.g1 <- mean(beta.g1)
        mean.e.g2 <- mean(beta.g2)
        
        ## counts for distal and close genes
        beta.d <- local.effects[local.effects$ePromOvlp.eQTLgene=="DistalGene","beta_eqtl"]
        beta.c <- local.effects[local.effects$ePromOvlp.eQTLgene=="CloseGene","beta_eqtl"]
        
        num.e.d <- length(beta.d)
        num.e.c <- length(beta.c)
        
        t.d <- unique(local.effects[local.effects$ePromOvlp.eQTLgene=="DistalGene","gtex_tissue"])
        t.c <- unique(local.effects[local.effects$ePromOvlp.eQTLgene=="CloseGene","gtex_tissue"])
        
        s.t.dist <- length(unique(t.d[t.d %in% t.c ], t.c[t.d %in% t.d ])) ## shared tissues
        
        t.s.dist <- length(unique(t.d, t.c)) ## total tissues
        
        mean.e.d <- mean(beta.d)
        mean.e.c <- mean(beta.c)
        
        signdif.dist <- (mean.e.d<0 & mean.e.c>0 ) | (mean.e.d>0 & mean.e.c<0 )
        
        
        eQTLs.affect.both.genes.eqtl.summary <- rbind(eQTLs.affect.both.genes.eqtl.summary, cbind(eQTLs.affect.both.gene[i,], num.e.g1,  num.e.g2,  s.t, t.s, mean.e.g1,  mean.e.g2,   num.e.d,  num.e.c , s.t.dist,  t.s.dist,mean.e.d, mean.e.c ,signdif.dist, type ) )
    }
}

## Select the eQTLs that show the reciprocal effect in the
## promoters were they act on both genes of a pair
chia.eqtl.support.eqtl.bimodal <- chia.eqtl.support[bieffect.index , ]
require(reshape2)
## Save in a table
table.bieffect.results.file <- file.path(results.dir,"chia_pet_interactions_eQTL_support_reciprocal_effect_in_pairs_strict.txt")
write.table(chia.eqtl.support.eqtl.bimodal,file=table.bieffect.results.file , quote=FALSE, sep="\t", row.names=FALSE)

## Report analysis of bimodal eQTLS

counts.eqtl.names <- c("eQTL_id", "uniquePair", "num.effect_gene1","num.effect_gene2","shared.tissues","total.tissues","meanBeta_gene1","meanBeta_gene2","num.effect_distal", "num.effect_close",
                       "shared.tissues.distance", "total.tissues.distance",  "meanBeta_distal","meanBeta_close","SignDiff","PromType" )
eQTLs.affect.both.genes.counts <-  eQTLs.affect.both.genes.eqtl.summary

colnames(eQTLs.affect.both.genes.counts) <- counts.eqtl.names

table.bieffect.eqtl.results.file <- file.path(results.dir,"bimodal_eQTLs_effect_in_pairs_strict_meanbetas.txt")
write.table(eQTLs.affect.both.genes.counts,file=table.bieffect.eqtl.results.file  , quote=FALSE, sep="\t", row.names=FALSE)
graph.eqtls <- melt(eQTLs.affect.both.genes.counts[,c("eQTL_id","uniquePair", "meanBeta_distal", "meanBeta_close", "PromType","SignDiff")], id=c("eQTL_id","uniquePair","PromType","SignDiff"), variable.name="Distance",  value.name="mean_beta_eqtl")


################
## Calculate the  mean difference of beta effect between samples
## Use none parametric test since the distributions are bimodal 

## For all eQTLS
mean.diff.w.test <- wilcox.test(chia.eqtl.support[chia.eqtl.support$ePromOvlp.eQTL2=="hProm","beta_eqtl" ] , chia.eqtl.support[chia.eqtl.support$ePromOvlp.eQTL2=="eProm", "beta_eqtl"], paired=FALSE, conf.int=TRUE)

## All distal eQTLs
chia.eqtl.support.distal <- chia.eqtl.support[chia.eqtl.support$ePromOvlp.eQTLgene=="DistalGene",]

table.distal.results.file <- file.path(results.dir,"distal_effect_eQTLs.txt")
write.table( chia.eqtl.support.distal ,file=table.distal.results.file  , quote=FALSE, sep="\t", row.names=FALSE)

distal.beta.hprom <- chia.eqtl.support.distal[chia.eqtl.support.distal$ePromOvlp.eQTL2=="hProm","beta_eqtl" ]
distal.beta.eprom <-  chia.eqtl.support.distal[chia.eqtl.support.distal$ePromOvlp.eQTL2=="eProm","beta_eqtl" ]

mean.diff.distal.w.test <- wilcox.test( distal.beta.hprom , distal.beta.eprom , conf.int=TRUE)

### Distal by effect
## split distribution into two
## Split done using the mixtools library.

library(mixtools)

## Code from  http://stats.stackexchange.com/questions/57993/how-to-explain-how-i-divided-a-bimodal-distribution-based-on-kernel-density-esti/78397#78397cite

## Function to find the cutoff based on the adjusted model
find.cutoff <- function(proba=0.5, i=index.lower) {
    ## Cutoff such that Pr[drawn from bad component] == proba
    f <- function(x) {
        proba - (model$lambda[i]*dnorm(x, model$mu[i], model$sigma[i]) /
                     (model$lambda[1]*dnorm(x, model$mu[1], model$sigma[1]) + model$lambda[2]*dnorm(x, model$mu[2], model$sigma[2])))
        }
        return(uniroot(f=f, lower=-2, upper=2, extendInt="yes")$root)  # Careful with division by zero if changing lower and upper
}

## Fit a two-component mixture model
x <- chia.eqtl.support.distal$beta_eqtl
model <- normalmixEM(x=x, k=2)
index.lower <- which.min(model$mu)  # Index of component with lower mean



cutoffs <- c(find.cutoff(proba=0.5), find.cutoff(proba=0.75))  # Around c(1.8, 1.5)

hist(x)
abline(v=cutoffs, col=c("red", "blue"), lty=2)

cutoff0.5 <- cutoffs[1]

chia.eqtl.support.distal.pos <- chia.eqtl.support.distal[chia.eqtl.support.distal$beta_eqtl>cutoff0.5 ,]
chia.eqtl.support.distal.neg <- chia.eqtl.support.distal[chia.eqtl.support.distal$beta_eqtl<=cutoff0.5 ,]


mean.diff.distal.pos.eqtl <- wilcox.test (chia.eqtl.support.distal.pos[chia.eqtl.support.distal.pos$ePromOvlp.eQTL2=="hProm","beta_eqtl"], chia.eqtl.support.distal.pos[chia.eqtl.support.distal.pos$ePromOvlp.eQTL2=="eProm","beta_eqtl"], alternative="less" ) ## beta in hProm is less than in eProm
mean.diff.distal.neg.eqtl <- wilcox.test (chia.eqtl.support.distal.neg[chia.eqtl.support.distal.neg$ePromOvlp.eQTL2=="hProm","beta_eqtl"], chia.eqtl.support.distal.neg[chia.eqtl.support.distal.neg$ePromOvlp.eQTL2=="eProm","beta_eqtl"] , alternative="greater") ## beta in hProm is higher than in eProm, negative values. 

## Comparison for asymmetrical pairs
chia.eqtl.support.asymmetrical <- chia.eqtl.support[chia.eqtl.support$ePromOvlp.eQTL==0 | chia.eqtl.support$ePromOvlp.eQTL==1, ]


mean.diff.w.test.assym <- wilcox.test( as.vector(chia.eqtl.support.asymmetrical[chia.eqtl.support.asymmetrical$ePromOvlp.eQTL2=="hProm","beta_eqtl" ] ), as.vector( chia.eqtl.support.asymmetrical[chia.eqtl.support.asymmetrical$ePromOvlp.eQTL2=="eProm", "beta_eqtl"] ), paired=FALSE, conf.int=TRUE)

## Only eQTLs acting in both genes of a pair
mean.diff.w.test.bimodal <- wilcox.test(chia.eqtl.support.eqtl.bimodal[chia.eqtl.support.eqtl.bimodal$ePromOvlp.eQTL2=="hProm","beta_eqtl" ] , chia.eqtl.support.eqtl.bimodal[chia.eqtl.support.eqtl.bimodal$ePromOvlp.eQTL2=="eProm", "beta_eqtl"], paired=FALSE, conf.int=TRUE )

################
## Test diference between bimodals for distal and close genes
#test <- chia.eqtl.support.eqtl.bimodal[  chia.eqtl.support.eqtl.bimodal$ePromOvlp.eQTL2=="eProm" ,]

dist.eprom <- chia.eqtl.support.eqtl.bimodal[ (chia.eqtl.support.eqtl.bimodal$ePromOvlp.eQTLgene=="DistalGene" & chia.eqtl.support.eqtl.bimodal$ePromOvlp.eQTL2=="eProm" ),]

dist.hprom <- chia.eqtl.support.eqtl.bimodal[ (chia.eqtl.support.eqtl.bimodal$ePromOvlp.eQTLgene=="DistalGene" & chia.eqtl.support.eqtl.bimodal$ePromOvlp.eQTL2=="hProm" ),]

close.eprom <- chia.eqtl.support.eqtl.bimodal[ (chia.eqtl.support.eqtl.bimodal$ePromOvlp.eQTLgene=="CloseGene" & chia.eqtl.support.eqtl.bimodal$ePromOvlp.eQTL2=="eProm" ),]

close.hprom <- chia.eqtl.support.eqtl.bimodal[ (chia.eqtl.support.eqtl.bimodal$ePromOvlp.eQTLgene=="CloseGene" & chia.eqtl.support.eqtl.bimodal$ePromOvlp.eQTL2=="hProm" ),]


## Compare distance and promoter betavalues

dist.epromVSdist.hprom <- wilcox.test(dist.eprom[,"beta_eqtl" ] , dist.hprom[, "beta_eqtl"], paired=FALSE )

dist.epromVSclose.eprom <- wilcox.test(dist.eprom[,"beta_eqtl" ] , close.eprom[, "beta_eqtl"], paired=FALSE )

dist.epromVSclose.hprom <- wilcox.test(dist.eprom[,"beta_eqtl" ] , close.hprom[, "beta_eqtl"], paired=FALSE )

dist.hpromVSclose.eprom <- wilcox.test(dist.hprom[,"beta_eqtl" ] , close.eprom[, "beta_eqtl"], paired=FALSE )

dist.hpromVSclose.hprom <- wilcox.test(dist.hprom[,"beta_eqtl" ] , close.hprom[, "beta_eqtl"], paired=FALSE )

close.epromVSclose.hprom <- wilcox.test(close.eprom[,"beta_eqtl" ] , close.hprom[, "beta_eqtl"], paired=FALSE )
    


## adjust pvalue

adjusted.pvalues <- p.adjust(c (mean.diff.w.test$p.value ,
                                mean.diff.w.test.assym$p.value,
                                mean.diff.w.test.bimodal$p.value ,
                                dist.epromVSdist.hprom$p.value ,
                                dist.epromVSclose.eprom$p.value ,
                                dist.epromVSclose.hprom$p.value ,
                                dist.hpromVSclose.eprom$p.value ,
                                dist.hpromVSclose.hprom$p.value ,
                                close.epromVSclose.hprom$p.value,
                                mean.diff.distal.w.test$p.value,
                                mean.diff.distal.pos.eqtl$p.value,
                                mean.diff.distal.neg.eqtl$p.value
                                ) , method="BH")

mean.diff.test.pval <-  adjusted.pvalues [1]
mean.diff.test.pval.assym <- adjusted.pvalues[2]
mean.diff.test.pval.bimodal<- adjusted.pvalues[3]

dist.epromVSdist.hprom.pval <- adjusted.pvalues[4]
dist.epromVSclose.eprom.pval <- adjusted.pvalues[5]
dist.epromVSclose.hprom.pval <- adjusted.pvalues[6]
dist.hpromVSclose.eprom.pval <- adjusted.pvalues[7]
dist.hpromVSclose.hprom.pval <- adjusted.pvalues[8]
close.epromVSclose.hprom.pval <- adjusted.pvalues[9]
mean.diff.distal.w.test.pval <- adjusted.pvalues[10]
mean.diff.distal.pos.eqtl.p.val <- adjusted.pvalues[11]
mean.diff.distal.neg.eqtl.p.val <- adjusted.pvalues[12]


rowNames <- c(
    "All_ePromVShProm", 
    "Assymetrical_ePromVShProm", 
    "Bimodal_ePromVShProm",
    "Distal_eProm_VS_Distal_hProm" ,
    "Distal_eProm_VS_Close_eProm" ,
    "Distal_eProm_VS_Close_hProm" ,
    "Distal_hProm_VS_Close_eProm" ,
    "Distal_hProm_VS_Close_hProm" ,
    "Close_eProm_VS_Close_hProm" ,
    "AllDistal_ePromVShPRom",
    "AllDistalPos_ePromVShPRom",
    "AllDistalNeg_ePromVShPRom"
   
    )

pvals <- c(
    mean.diff.w.test$p.value ,
    mean.diff.w.test.assym$p.value,
    mean.diff.w.test.bimodal$p.value ,
    dist.epromVSdist.hprom$p.value ,
    dist.epromVSclose.eprom$p.value ,
    dist.epromVSclose.hprom$p.value ,
    dist.hpromVSclose.eprom$p.value ,
    dist.hpromVSclose.hprom$p.value ,
    close.epromVSclose.hprom$p.value,
    mean.diff.distal.w.test$p.value,
    mean.diff.distal.pos.eqtl$p.value,
    mean.diff.distal.neg.eqtl$p.value
    )

ad.pvals <- c(
    mean.diff.test.pval ,
    mean.diff.test.pval.assym ,
    mean.diff.test.pval.bimodal ,
    
    dist.epromVSdist.hprom.pval ,
    dist.epromVSclose.eprom.pval ,
    dist.epromVSclose.hprom.pval ,
    dist.hpromVSclose.eprom.pval ,
    dist.hpromVSclose.hprom.pval ,
    close.epromVSclose.hprom.pval ,
    mean.diff.distal.w.test.pval ,
    mean.diff.distal.pos.eqtl.p.val ,
    mean.diff.distal.neg.eqtl.p.val 
    
    )


comp.dist.pvalues <- data.frame(rowNames, pvals,ad.pvals     )


file.pvals <- file.path(results.dir,"Pvalues_beta_comparisons.txt")
write.table( comp.dist.pvalues , file=file.pvals , quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

################
## Report summary statistics

eQTLs.beta.effect.summary.n <- data.frame("eQTL_selection"=character(),
                                          "Min." = numeric(),
                                          "1st Qu." = numeric(),
                                          "Median" = numeric(),
                                          "Mean" = numeric(),
                                       "3rd Qu." = numeric(),
                                          "Max." = numeric()
                                          )

eQTLs.beta.effect.summary <- rbind (
    c("All_eqtl",summary(chia.eqtl.support$beta_eqtl) ) ,
    c("Bimodal_eqtl", summary(chia.eqtl.support.eqtl.bimodal$beta_eqtl) ) ,
    
    c("hProm_eqtl",summary( chia.eqtl.support[chia.eqtl.support$ePromOvlp.eQTL2=="hProm","beta_eqtl" ] ) ) ,
    c("eProm_eqtl",summary( chia.eqtl.support[chia.eqtl.support$ePromOvlp.eQTL2=="eProm","beta_eqtl" ] ) ) ,

    c("Distal_hProm_eqtl",summary( chia.eqtl.support.distal[chia.eqtl.support.distal$ePromOvlp.eQTL2=="hProm","beta_eqtl" ] ) ) ,
    c("Distal_eProm_eqtl",summary( chia.eqtl.support.distal[chia.eqtl.support.distal$ePromOvlp.eQTL2=="eProm","beta_eqtl" ] ) ) ,

    c("Distal_hProm_Pos_eqtl",summary( chia.eqtl.support.distal.pos[chia.eqtl.support.distal.pos$ePromOvlp.eQTL2=="hProm","beta_eqtl"] ) ) ,
    c("Distal_eProm_Pos_eqtl",summary( chia.eqtl.support.distal.pos[chia.eqtl.support.distal.pos$ePromOvlp.eQTL2=="eProm","beta_eqtl"] ) ) ,
    
    c("Distal_hProm_Neg_eqtl",summary( chia.eqtl.support.distal.neg[chia.eqtl.support.distal.neg$ePromOvlp.eQTL2=="hProm","beta_eqtl"] ) ) ,
    c("Distal_eProm_Neg_eqtl",summary( chia.eqtl.support.distal.neg[chia.eqtl.support.distal.neg$ePromOvlp.eQTL2=="eProm","beta_eqtl"] ) ) ,
    
    
    c("Asymmetrical_eqtl",summary( chia.eqtl.support.asymmetrical$beta_eqtl ) ) ,
    
    c("Asymmetrical_hProm_eqtl",summary( chia.eqtl.support.asymmetrical[chia.eqtl.support.asymmetrical$ePromOvlp.eQTL2=="hProm","beta_eqtl" ] )) ,
    c("Asymmetrical_eProm_eqtl",summary( chia.eqtl.support.asymmetrical[chia.eqtl.support.asymmetrical$ePromOvlp.eQTL2=="eProm", "beta_eqtl"] )) ,
    
    c("Bimodal_hProm_eqtl",summary( chia.eqtl.support.eqtl.bimodal[chia.eqtl.support.eqtl.bimodal$ePromOvlp.eQTL2=="hProm","beta_eqtl" ] )) ,
    c("Bimodal_eProm_eqtl",summary( chia.eqtl.support.eqtl.bimodal[chia.eqtl.support.eqtl.bimodal$ePromOvlp.eQTL2=="eProm", "beta_eqtl"] )) ,
    
    c("Bimodal_distal_eProm_eqtl", summary(dist.eprom$beta_eqtl)) ,
    c("Bimodal_distal_hProm_eqtl", summary(dist.hprom$beta_eqtl)),
    
    c("Bimodal_close_hProm_eqtl", summary(close.hprom$beta_eqtl)) ,
    c("Bimodal_close_eProm_eqtl", summary(close.eprom$beta_eqtl))
    )

eQTLs.beta.effect.summary <- rbind(eQTLs.beta.effect.summary.n, eQTLs.beta.effect.summary)
colnames(eQTLs.beta.effect.summary) <-  c("eQTLs.beta.effect.summary" , "Min." , "1st Qu." ,  "Median" , "Mean" , "3rd Qu." , "Max.")

file.summary <- file.path(results.dir,"chia_pet_eQTLs_summary_statistics2.txt")
write.table(eQTLs.beta.effect.summary , file=file.summary, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

################
## Draw boxplots

font <- "Arial Narrow"
## draw significance line
df1 <- data.frame(a = c(1, 1,2,2), b = c(2, 2.2, 2.2, 2))


## All eQTLS comparisons
bp.all <- ggplot (chia.eqtl.support, aes(factor( ePromOvlp.eQTL2), beta_eqtl )) + theme(text=element_text(family=font, size=12)) +
    geom_violin(aes(color=factor( ePromOvlp.eQTL2 )) )  +  coord_flip() +
        ggtitle ("Beta effect comparison, all eQTLS") +
            xlab("eQTL possition") + ylab ("Beta effect") +
                geom_boxplot(width=0.1) +
                    scale_colour_discrete( name="Promoter Type", labels=c("ePromoter", "hPromoter"))  +
                    geom_line(data = df1, aes(x = a, y = b)) + annotate("text", x = 1.5, y = 1.5, label=paste("AdjustedPval=",format(mean.diff.test.pval, scientific=TRUE, digits=3)) , size = 4, family=font)

file.bp.all <- file.path(results.dir,"chia_pet_interactions_Olvp_eQTL_betaEffect_epromVShprom_V2.pdf")
ggsave(file.bp.all, bp.all) 

## Plot only distal

bp.distal <- ggplot (chia.eqtl.support.distal, aes(factor( ePromOvlp.eQTL2), beta_eqtl )) + theme(text=element_text(family=font, size=12)) +
    geom_violin(aes(color=factor( ePromOvlp.eQTL2 )) )  +  coord_flip() +
        ggtitle ("Beta effect comparison, distal eQTLS") +
            xlab("eQTL possition") + ylab ("Beta effect") +
                geom_boxplot(width=0.1) +
                    scale_colour_discrete( name="Promoter Type", labels=c("ePromoter", "hPromoter"))  +
                    geom_line(data = df1, aes(x = a, y = b)) + annotate("text", x = 1.5, y = 1.5, label=paste("AdjustedPval= NS",format(mean.diff.distal.w.test.pval, scientific=TRUE, digits=3)) , size = 4, family=font)

file.bp.distal <- file.path(results.dir,"chia_pet_interactions_Olvp_Distal_eQTL_betaEffect_epromVShprom_V2.pdf")
ggsave(file.bp.distal, bp.distal) 

## Plot Distal separated by positive and negative effect

## Positive effect
bp.distal.pos <- ggplot (chia.eqtl.support.distal.pos, aes(factor( ePromOvlp.eQTL2), beta_eqtl )) + theme(text=element_text(family=font, size=12)) +
    geom_violin(aes(color=factor( ePromOvlp.eQTL2 )) )  +  coord_flip() +
        ggtitle ("Beta effect comparison, distal eQTLS") +
            xlab("eQTL possition") + ylab ("Beta effect") +
                geom_boxplot(width=0.1) +
                    scale_colour_discrete( name="Promoter Type", labels=c("ePromoter", "hPromoter"))  +
                    geom_line(data = df1, aes(x = a, y = b)) + annotate("text", x = 1.5, y = 1.5, label=paste("AdjustedPval= ",format( mean.diff.distal.pos.eqtl.p.val , scientific=TRUE, digits=3)) , size = 4, family=font)

file.bp.distal.pos <- file.path(results.dir,"chia_pet_interactions_Olvp_Distal_eQTL_betaEffect_epromVShprom_OnlyPositiveEffect_V2.pdf")
ggsave(file.bp.distal.pos, bp.distal.pos) 


## Negative effect
bp.distal.neg <- ggplot (chia.eqtl.support.distal.neg, aes(factor( ePromOvlp.eQTL2), beta_eqtl )) + theme(text=element_text(family=font, size=12)) +
    geom_violin(aes(color=factor( ePromOvlp.eQTL2 )) )  +  coord_flip() +
        ggtitle ("Beta effect comparison, distal eQTLS") +
            xlab("eQTL possition") + ylab ("Beta effect") +
                geom_boxplot(width=0.1) +
                    scale_colour_discrete( name="Promoter Type", labels=c("ePromoter", "hPromoter"))  +
                    geom_line(data = df1, aes(x = a, y = b)) + annotate("text", x = 1.5, y = 1.5, label=paste("AdjustedPval= ",format( mean.diff.distal.neg.eqtl.p.val , scientific=TRUE, digits=3)) , size = 4, family=font)

file.bp.distal.neg <- file.path(results.dir,"chia_pet_interactions_Olvp_Distal_eQTL_betaEffect_epromVShprom_OnlyNegativeEffect_V2.pdf")
ggsave(file.bp.distal.neg, bp.distal.neg) 



## Only assimetrical pairs hPromoter pair with ePromoter

bp.assym <- ggplot (chia.eqtl.support.asymmetrical, aes(factor( ePromOvlp.eQTL2), beta_eqtl )) + theme(text=element_text(family=font , size=12)) +
    geom_violin(aes(color=factor( ePromOvlp.eQTL2 )) )  + coord_flip() +
        ggtitle ("Beta effect comparison, eQTLS in asymmetrical promoters") +
            xlab("eQTL possition") + ylab ("Beta effect") +
                scale_colour_discrete( name="Promoter Type", labels=c("ePromoter", "hPromoter"))  +
                    geom_boxplot(width=0.1) +
                    geom_line(data = df1, aes(x = a, y = b)) + annotate("text", x = 1.5, y = 1.5, label=paste("AdjustedPval=",format(mean.diff.test.pval.assym,  scientific=TRUE, digits=3  )) , size = 4, family=font)

file.bp.assym  <- file.path(results.dir,"chia_pet_interactions_Olvp_eQTL_betaEffect_epromVShprom_asymmetrical_pairs_V2.pdf")
ggsave(file.bp.assym , bp.assym ) 

## Strict reciprocal eQTL affecting both genes in a pair

bp.bimodal <- ggplot ( chia.eqtl.support.eqtl.bimodal , aes(factor( ePromOvlp.eQTL2), beta_eqtl )) + theme(text=element_text(family=font, size=12)) +
    geom_violin(aes(color=factor( ePromOvlp.eQTL2 )) )  + coord_flip() +
        ggtitle ("Beta effect comparison between eQTLs affecting both genes") +
            xlab("eQTL possition") + ylab ("Beta effect") +
                scale_colour_discrete( name="Promoter type", labels=c("ePromoter", "hPromoter")) +
                    geom_boxplot(width=0.1) +
                        geom_line(data = df1, aes(x = a, y = b)) + annotate("text", x = 1.5, y = 1.5, label=paste("AdjustedPval=",format(mean.diff.test.pval.bimodal,scientific=TRUE, digits=3)) , size = 4, family=font)

file.bp.bimodal <- file.path(results.dir,"chia_pet_interactions_Olvp_eQTL_betaEffect_epromVShprom_eTQL_effect2Bothgenes_strict_V2.pdf")
ggsave(file.bp.bimodal, bp.bimodal)

## Strict by close or distal interaction


## draw significance line
dfa <- data.frame(a = c(1, 1,2,2), b = c(2, 2.2, 2.2, 2))


bp.bimodal.dist <-ggplot ( chia.eqtl.support.eqtl.bimodal , aes(factor( ePromOvlp.eQTLgene), beta_eqtl )) + coord_flip() +
    geom_violin(aes(color=factor( ePromOvlp.eQTL2 )) ) + facet_grid(.~ ePromOvlp.eQTL2)+ theme(text=element_text(family=font, size=12)) +
        ggtitle ("Beta effect comparison between eQTLs affecting both genes") + 
            xlab("eQTL possition") + ylab ("Beta effect") +
                scale_colour_discrete( name="Promoter type", labels=c("ePromoter", "hPromoter")) +
                    geom_boxplot(width=0.1) +
                        geom_line(data = dfa, aes(x = a, y = b)) + annotate("text", x = 1.5, y = 1, label=paste("AdjustedPval=",c( format (dist.epromVSdist.hprom.pval ,  scientific=TRUE, digits=3), format(close.epromVSclose.hprom.pval,  scientific=TRUE, digits=3))),  , size = 2.6, family=font)

file.bp.bimodal.dist <- file.path(results.dir,"chia_pet_interactions_Olvp_eQTL_betaEffect_epromVShprom_eTQL_effect2Bothgenes_strict_byDistanceV2.pdf")
ggsave(file.bp.bimodal.dist, bp.bimodal.dist)


## Draw eQTL sing change box plot
## mean beta values
     
bp.bimodal.sign.change.dist <- ggplot ( graph.eqtls , aes(factor(Distance), mean_beta_eqtl )) + coord_flip() +
    geom_violin(aes(color=factor( SignDiff )) ) + facet_grid(.~ SignDiff  )+ theme(text=element_text(family=font, size=12))+
        ggtitle ("Beta effect comparison between eQTLs affecting both genes, change in effect") + 
            xlab("eQTL effect change") + ylab ("Beta effect") +
                scale_colour_discrete( name="Change in effect", labels=c("NoChange", "Change")) +
                    geom_boxplot(width=0.1)

#+
 #                       geom_line(data = dfa, aes(x = a, y = b))

#+ annotate("text", x = 1.5, y = 1, label=paste("AdjustedPval=",c( format (dist.epromVSdist.hprom.pval ,  scientific=TRUE, digits=3), format(close.epromVSclose.hprom.pval,  scientific=TRUE, digits=3))),  , size = 2.6, family=font)

file.bp.bimodal.sign.change.dist <- file.path(results.dir,"eQTL_meanbetaEffect_signChange_byDistanceV2.pdf")
ggsave(file.bp.bimodal.sign.change.dist , bp.bimodal.sign.change.dist )

## Mean beta values between long and distal effect

bp.bimodal.meanbeta.dist <- ggplot ( graph.eqtls , aes(factor(Distance), mean_beta_eqtl )) + coord_flip() +
    geom_violin(aes(color=factor( PromType )) ) + facet_grid(.~ PromType  )+ theme(text=element_text(family=font, size=12))+
        ggtitle ("Beta effect comparison between eQTLs affecting both genes, by Promoter Type") + 
            xlab("eQTL effect change") + ylab ("Beta effect") +
                scale_colour_discrete( name="Change in effect", labels=c("ePromoter", "hPromoter")) +
                    geom_boxplot(width=0.1)
file.bp.bimodal.meanbeta.dist<- file.path(results.dir,"eQTL_meanbetaEffect_byProm_byDistanceV2.pdf")
ggsave(file.bp.bimodal.meanbeta.dist , bp.bimodal.meanbeta.dist )



## Mean beta values between long and distal effect
graph.eqtls.singChange <- graph.eqtls[graph.eqtls$SignDiff==TRUE,]

bp.bimodal.meanbeta.dist.signChange <- ggplot ( graph.eqtls.singChange  , aes(factor(Distance), mean_beta_eqtl )) + coord_flip() +
    geom_violin(aes(color=factor( PromType )) ) + facet_grid(.~ PromType  )+ theme(text=element_text(family=font, size=12))+
        ggtitle ("Beta effect comparison between eQTLs affecting both genes, by Promoter Type, EffectSignChange") + 
            xlab("eQTL effect change") + ylab ("Beta effect") +
                scale_colour_discrete( name="Change in effect", labels=c("ePromoter", "hPromoter")) +
                    geom_boxplot(width=0.1)
file.bp.bimodal.meanbeta.dist.signChange<- file.path(results.dir,"eQTL_meanbetaEffect_byProm_byDistanceV2_effectSignChange.pdf")
ggsave(file.bp.bimodal.meanbeta.dist.signChange , bp.bimodal.meanbeta.dist.signChange )

################
## correlation analysis of eQTLs affecting gene pairs
plotCorrelation <- function (df,fileName, mainTitle){
    correlationB <- cor.test(df$mean.e.d, df$mean.e.c)
    lmvalue <- lm (mean.e.d~mean.e.c,df )
    pdf(fileName)
    plot(df$mean.e.d, df$mean.e.c, main=mainTitle, xlab="Mean Beta Effect Close Gene", ylab="Mean Beta Effect Distal Gene")
    abline(lmvalue, col="red")
    legend(x="topright",paste("Pearson's Correlation",correlationB$estimate ))
    dev.off()
}

file.corrAll<- file.path(results.dir,"eQTL_meanbetaEffect_CorrelationDistance.pdf")
plotCorrelation(df=eQTLs.affect.both.genes.eqtl.summary , fileName=file.corrAll, mainTitle="Reciprocal eQTLs \n Mean Beta Effect, All eQTLS")

file.correprom<- file.path(results.dir,"eQTL_meanbetaEffect_CorrelationDistance_hProm.pdf")
eQTLs.affect.both.genes.eqtl.summary.hprom <- eQTLs.affect.both.genes.eqtl.summary[eQTLs.affect.both.genes.eqtl.summary$type=="hProm",]
plotCorrelation(df=eQTLs.affect.both.genes.eqtl.summary.hprom , fileName=file.correprom, mainTitle="Reciprocal eQTLs \n Mean Beta Effect, hProm")

file.corrhprom<- file.path(results.dir,"eQTL_meanbetaEffect_CorrelationDistance_eProm.pdf")
eQTLs.affect.both.genes.eqtl.summary.eprom <- eQTLs.affect.both.genes.eqtl.summary[eQTLs.affect.both.genes.eqtl.summary$type=="eProm",]
plotCorrelation(df=eQTLs.affect.both.genes.eqtl.summary.eprom , fileName=file.corrhprom, mainTitle="Reciprocal eQTLs \n Mean Beta Effect,eProm")





################
## Report table

## Total number of promoters with ChIAPET interaction
genes.inChiaPet <- length(unique(c(results$gene1, results$gene2)))

## Total number of Epromoters with ChIAPET interaction

genes1.eProm <- results[ which (grepl ("EPromoter" ,results$promType1 ) %in% TRUE) , c("gene1")] 
genes2.eProm <- results[ which (grepl ("EPromoter" ,results$promType2 ) %in% TRUE) , c("gene2")] 

genes.eProm <- unique(c(genes1.eProm , genes2.eProm ))

genes.in.eProm.chiaPet <- length (genes.eProm)

## Total number of Hpromoters with ChIAPET interaction

genes1.hprom <- results[ which (grepl ("hPromoter" ,results$promType1 ) %in% TRUE) , c("gene1")] 
genes2.hprom <- results[ which (grepl ("hPromoter" ,results$promType2 ) %in% TRUE) , c("gene2")] 

genes.hprom <- unique(c(genes1.hprom , genes2.hprom ))

genes.in.hProm.chiaPet <- length (genes.hprom)


## Total number of promoter with ChiaPet interaction that contain an eQTL regulating distal genes
chia.eqtl.support.distal <- chia.eqtl.support[chia.eqtl.support$ePromOvlp.eQTLgene=="DistalGene",]
genes.chiaPet.eQTL.distal <- length(unique( chia.eqtl.support.distal$gene_id_cap))

## Total number of ePromoter with ChiaPet that contain an eQTL regulating distal genes

chia.eqtl.support.distal.eprom <- chia.eqtl.support[chia.eqtl.support$ePromOvlp.eQTLgene=="DistalGene" & chia.eqtl.support$ePromOvlp.eQTL2=="eProm" ,]
genes.chiaPet.eQTL.eProm.distal <- length(unique (chia.eqtl.support.distal.eprom$gene_id_cap ))

## Total number of hpromoter with ChiaPet that contain an eQTL regulating distal genes

chia.eqtl.support.distal.hprom <- chia.eqtl.support[chia.eqtl.support$ePromOvlp.eQTLgene=="DistalGene" & chia.eqtl.support$ePromOvlp.eQTL2=="hProm" ,]
genes.chiaPet.eQTL.hProm.distal <- length(unique (chia.eqtl.support.distal.hprom$gene_id_cap ))


## Total number of promoter with ChiaPet interaction that contain an eQTL regulating close genes
chia.eqtl.support.close <- chia.eqtl.support[chia.eqtl.support$ePromOvlp.eQTLgene=="CloseGene",]
genes.chiaPet.eQTL.close <- length(unique( chia.eqtl.support.close$gene_id_cap))

## Total number of ePromoter with ChiaPet that contain an eQTL regulating close genes

chia.eqtl.support.close.eprom <- chia.eqtl.support[chia.eqtl.support$ePromOvlp.eQTLgene=="CloseGene" & chia.eqtl.support$ePromOvlp.eQTL2=="eProm" ,]
genes.chiaPet.eQTL.eProm.close <- length(unique (chia.eqtl.support.close.eprom$gene_id_cap ))

## Total number of hpromoter with ChiaPet that contain an eQTL regulating close genes

chia.eqtl.support.close.hprom <- chia.eqtl.support[chia.eqtl.support$ePromOvlp.eQTLgene=="CloseGene" & chia.eqtl.support$ePromOvlp.eQTL2=="hProm" ,]
genes.chiaPet.eQTL.hProm.close <- length(unique (chia.eqtl.support.close.hprom$gene_id_cap ))


## Report
promoter.count.eqtl.d <- data.frame(
    PromoterType=character(),
    Number=numeric()
    )


promoter.count.eqtl.names <- c("Promoter Type", "Number of Promoters in type")

promoter.count.eqtl <- rbind (
    c("Num. Proms in ChiaPet Interactions" , genes.inChiaPet ),
    c("Num. eProms in ChiaPet Interactions" , genes.in.eProm.chiaPet ),
    c("Num. hProms in ChiaPet Interactions" , genes.in.hProm.chiaPet ),    
    c("Num. Proms in ChiaPet with eQTL Distal Target" ,genes.chiaPet.eQTL.distal ),
    c("Num. eProms in ChiaPet with eQTL Distal Target", genes.chiaPet.eQTL.eProm.distal ),
    c("Num. hProms in ChiaPet with eQTL Distal Target", genes.chiaPet.eQTL.hProm.distal ),
    c("Num. Proms in ChiaPet with eQTL Close Target" ,genes.chiaPet.eQTL.close ),
    c("Num. eProms in ChiaPet with eQTL Close Target", genes.chiaPet.eQTL.eProm.close ),
    c("Num. hProms in ChiaPet with eQTL Close Target", genes.chiaPet.eQTL.hProm.close )
    )

promoter.count.eqtl <- rbind (promoter.count.eqtl.d  ,promoter.count.eqtl )
colnames(promoter.count.eqtl) <- promoter.count.eqtl.names

file.prom.counts <- file.path(results.dir,"promoters_in_chia_pet_w_eQTLs_counts2.txt")
write.table( promoter.count.eqtl , file=file.prom.counts , quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)


comp.all.distal.tab <- matrix( c(genes.in.eProm.chiaPet, genes.in.hProm.chiaPet,  genes.chiaPet.eQTL.eProm.distal , genes.chiaPet.eQTL.hProm.distal), byrow=TRUE, ncol=2 )
colnames(comp.all.distal.tab) <- c("eProm","hProm")
rownames(comp.all.distal.tab) <- c("All Prom","Prom w eQTL Distal")
prop.test (comp.all.distal.tab)



################
## From ChIaPet pairs perspetive
## count:
## The number of eQTLs falling in ePromoter
## The number of eQTLs falling in hPromoter
## The number of eQTLs that are reciprocal and fall in ePromoter
## The number of eQTLs that are reciprocal and fall in hPromoter

chiaPet.pair.report <- data.frame(uniquePair = character(),
                                  Promoter_1_Type = character() ,
                                  Promoter_2_Type   = character (),
                                  Num.eQTLs = numeric(),
                                  Num.eQTLs.eProm = numeric(),
                                  Num.eQTLs.hProm = numeric(),
                                  Num.eQTLs.Re.eProm = numeric(),
                                  Num.eQTLs.Re.hProm = numeric(),
                                  Frac.eQTLs.in.Re.eProm = numeric(),
                                  Frac.eQTLs.in.Re.hProm = numeric()
                                  )

chiaPet.pair.report.names <- c( "uniquePair", "Promoter_1_Type" , "Promoter_2_Type",
                               "Num.eQTLs","Num.eQTLs.eProm","Num.eQTLs.hProm",
                               "Num.eQTLs.Re.eProm", "Num.eQTLs.Re.hProm",
                               "Frac.eQTLs.in.Re.eProm", "Frac.eQTLs.in.Re.hProm" )

for (pair in unique.pairs ){
    pair.info <- clean.chiapet.table [clean.chiapet.table$uniquePair==pair, ]

    if (dim(pair.info)[1] >1){
        print ("Error, several pair report")
        print (pair)
        stop()
    }
    prom1 <- pair.info$Promoter_1_Type
    prom2 <- pair.info$Promoter_2_Type

    num.eqtl <- length(unique(chia.eqtl.support[chia.eqtl.support$uniquePair==pair,c("snp_eqtl")]))
    
    num.eqtl.eprom <- length(unique(chia.eqtl.support[  ( chia.eqtl.support$uniquePair==pair & chia.eqtl.support$ePromOvlp.eQTL2=="eProm" )  ,c("snp_eqtl")]))
    num.eqtl.hprom <- length(unique(chia.eqtl.support[  ( chia.eqtl.support$uniquePair==pair & chia.eqtl.support$ePromOvlp.eQTL2=="hProm" )  ,c("snp_eqtl")]))

    num.eqtl.re.eprom <- length(unique(chia.eqtl.support.eqtl.bimodal[  ( chia.eqtl.support.eqtl.bimodal$uniquePair==pair & chia.eqtl.support.eqtl.bimodal$ePromOvlp.eQTL2=="eProm" )  ,c("snp_eqtl")]))
    num.eqtl.re.hprom <- length(unique(chia.eqtl.support.eqtl.bimodal[  ( chia.eqtl.support.eqtl.bimodal$uniquePair==pair & chia.eqtl.support.eqtl.bimodal$ePromOvlp.eQTL2=="hProm" )  ,c("snp_eqtl")]))

    fraction.eqtls.recip.eprom <- num.eqtl.re.eprom / num.eqtl.eprom 
    fraction.eqtls.recip.hprom <- num.eqtl.re.hprom / num.eqtl.hprom 
    
}







