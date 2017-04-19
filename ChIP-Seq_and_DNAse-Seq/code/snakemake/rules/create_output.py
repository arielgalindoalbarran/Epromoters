rule create_output_table: 
    input: expand("RESULTS/PROM_PROM/{cond}-{exp}_fltr_pp_chiapet_overlap_prom.stat", exp=EXP, cond=COND)
    output: "RESULTS/OUTPUT/table_summary_stats.txt"
    params: ppn="nodes=1:ppn=1"
    shell: """
    for i in `ls RESULTS/PROM_PROM/*log RESULTS/PROM_PROM/*stat` ; do awk 'BEGIN{{FS=OFS="\t"}}{{print $1,$2,FILENAME}}' $i; done >  {output}
    """
    