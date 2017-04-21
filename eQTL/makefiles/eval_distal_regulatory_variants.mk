################################################################
## Get eQTL data regulatory variants 

## The file used as input in this analisis resulted from the R-script
## R-scripts/chiapet_analysis_eQTLs.R 
## Process with R to mark the ones that kill or create a site.
## Requires RSAT to run (http://rsat.eu/)

################
## directories
include ${RSAT}/makefiles/util.mk
MAIN_DIR=./
MAKEFILE=${MAIN_DIR}/makefiles/eval_distal_regulatory_variants.mk


RESULTS_DIR=${MAIN_DIR}/results
ANYLISIS=chia_pet
VERSION=20160601

SOURCE_DIR=${RESULTS_DIR}/${ANYLISIS}${VERSION}
SOURCE_FILE_NAME=distal_effect_eQTLs
SOURCE_FILE=${SOURCE_DIR}/${SOURCE_FILE_NAME}.txt
VAR_SCAN_DIR=${SOURCE_DIR}/variation-scan
VARBED_FILE=${VAR_SCAN_DIR}/${SOURCE_FILE_NAME}.varBed
VARSEQ_FILE=${VAR_SCAN_DIR}/${SOURCE_FILE_NAME}.varSeq

VARBED_HEADER="\#chr	start	end	strand	id	ref	alt	so_term	validate	minor_allele_freq	is_supvar	in_supvar"


RETRIEVE_VAR_CMD=retrieve-variation-seq -i ${VARBED_FILE}  -format varBed -species Homo_sapiens -release 79 -assembly GRCh37 -mml 30  -o ${VARSEQ_FILE}
VAR_SCAN_CMD=variation-scan -i ${ONE_TISSUE_varSeq} -m ${MOTIF_COLLECTION} -uth pval ${PVAL} -lth pval_ratio ${PVAL_RATIO} \
	-bg ${BG_MODEL} -o ${ONE_TISSUE_varresult}

create_variation_bed_file:
	@echo ${SOURCE_FILE}
	@mkdir -p ${VAR_SCAN_DIR}
	@echo ${VARBED_HEADER} > ${VARBED_FILE}
	@grep -v "^gene1" ${SOURCE_FILE}| perl -pe 's/ +//g'| awk '{print $$9"\t"$$10"\t"$$11"\t+\t"$$12"\t"$$13"\t"$$14"\tSNV\t0\tNA\t0\t0"} ' | perl -pe 's/chr//g' |sort -u >> ${VARBED_FILE}
	@echo ${VARBED_FILE}


retrieve_variation_seq:
	${RETRIEVE_VAR_CMD}
	@echo ${VARSEQ_FILE}



DATA_FOLDER=${MAIN_DIR}/data
BG_MODEL=${DATA_FOLDER}/bg_model/All_promoterRegions_BG_model_mkv_2.oligos
MOTIF_FOLDER=${DATA_FOLDER}/motifs
MOTIF_COLLECTION=${MOTIF_FOLDER}/Hocomoco_Jaspar_nonredundant_motifs_named_clusters_20151217.tf
MOTIF_suffix=allClusteredMotifs


PVAL=1e-3
PVAL_RATIO=10
VARSEQ_RESULTS_FILE=${VAR_SCAN_DIR}/${SOURCE_FILE_NAME}_varscan_pval${PVAL}_pratio${PVAL_RATIO}_${MOTIF_suffix}.txt


VAR_SCAN_CMD=variation-scan -i ${VARSEQ_FILE} -m ${MOTIF_COLLECTION} -uth pval ${PVAL} -lth pval_ratio ${PVAL_RATIO} \
	-bg ${BG_MODEL} -o ${VARSEQ_RESULTS_FILE}

variation_scan:
	${VAR_SCAN_CMD}


## Get only variants affecting specific pairs.
PAIRS="MARCH7_BAZ2B|CSDE1_SIKE1|CCNYL1_METTL21A"

SOURCE_FILE_SELECTED_NAME=${SOURCE_FILE_NAME}_${MOTIF_suffix}_pairs_crispr
SOURCE_FILE_SELECTED=${SOURCE_DIR}/${SOURCE_FILE_SELECTED_NAME}.txt

select_pairs:
	grep -E ${PAIRS} ${SOURCE_FILE} | grep eProm > ${SOURCE_FILE_SELECTED}

run_retrieve_var_selected_pairs:
	${MAKE}  create_variation_bed_file retrieve_variation_seq SOURCE_FILE_NAME=${SOURCE_FILE_SELECTED_NAME}


MOTIF_COLLECTION_SELECTED=${MOTIF_FOLDER}/Hocomoco_Jaspar_nonredundant_motifs_named_clusters_20151217.tf
MOTIF_suffix_SELECTED=allClusteredMotifs
run_var_scan_pairs:
	${MAKE}  variation_scan SOURCE_FILE_NAME=${SOURCE_FILE_SELECTED_NAME}

