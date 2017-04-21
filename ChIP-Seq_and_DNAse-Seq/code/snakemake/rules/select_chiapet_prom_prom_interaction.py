rule select_chiapet_prom_prom_interaction:
    input: txt="RESULTS/FILTER_ON_CHIAPET/{exp}_fltr.txt", \
        all_prom="DATA/allResults_FC_unix_1k_1k.bed",\
        enhProm="RESULTS/ENH_PROM_AND_CTRL/{cond}.bed"
    output: txt="RESULTS/PROM_PROM/{cond}-{exp}_fltr_pp_chiapet_overlap_prom.txt", \
            bed="RESULTS/PROM_PROM/{cond}-{exp}_fltr_pp_chiapet_overlap_prom.bed", \
            gtf="RESULTS/PROM_PROM/{cond}-{exp}_fltr_pp_chiapet_overlap_prom.gtf", \
            log="RESULTS/PROM_PROM/{cond}-{exp}_fltr_pp_chiapet_overlap_prom.log", \
            ip="RESULTS/PROM_PROM/{cond}-{exp}_interacting_promoters.bed",\
            eip="RESULTS/PROM_PROM/{cond}-{exp}_enhProm_in_interacting_promoters.bed"
    threads:
        1
    shell: """
    # Select R1 fragment (line num as bed name)
    awk 'BEGIN{{FS=OFS="\\t"}}{{print $2,$3,$4,$8}}' {input.txt} | perl -ne 'print if($. > 1)' > {input.txt}_R1_tmp.bed
    touch {input.txt}_R1_tmp.bed # if not created
    # Select R2 fragment (line num as bed name)
    awk 'BEGIN{{FS=OFS="\\t"}}{{print $5,$6,$7,$8}}' {input.txt} | perl -ne 'print if($. > 1)' > {input.txt}_R2_tmp.bed
    touch {input.txt}_R1_tmp.bed # if not created
    
    # Check intersection of R1 and R2 ends (use uniq counting, -u).
    bedtools intersect -u -a {input.txt}_R1_tmp.bed -b {input.all_prom} > {output.txt}_chiapet_overlap_prom_R1.bed
    bedtools intersect -u -a {input.txt}_R2_tmp.bed -b {input.all_prom} > {output.txt}_chiapet_overlap_prom_R2.bed
    
    # Check that a particular interaction id was found in both bedtools results
    (cut -f4 {output.txt}_chiapet_overlap_prom_R2.bed ; cut -f4 {output.txt}_chiapet_overlap_prom_R1.bed)  \
        | sort -g  | uniq -c  | grep -w 2 | sed 's/ *2 //'> {output.txt}.id_list 

    #retrieve the list of chiapet fragments
    awk -F '\\t' 'NR==FNR {{id[$1]; next}} $8 in id || FNR == 1' {output.txt}.id_list {input.txt} > {output.txt}

    # Convert the list of chiapet fragments into GTF format
    awk 'BEGIN{{FS=OFS="\t"}}{{if(NR > 2){{print $2,"SOURCE", "exon", $3,$4,"0","+",0,"transcript_id \\"FRAG_"FNR"\\";\\n"$5,"SOURCE", "exon",$6,$7,"0","+",0,"transcript_id \\"FRAG_"FNR"\\";"}}}}' {output.txt} > {output.gtf}
    
    
    # Get the list of interacting promoters
    cut -f2,3,4 {output.txt} | perl -ne 'print if($. > 1)' > {output.bed}
    cut -f5,6,7 {output.txt} | perl -ne 'print if($. > 1)' >> {output.bed}
    bedtools intersect -u -a {input.all_prom} -b {output.bed} > {output.ip}
    
    
    # basic stats
    echo -ne "Number_of_filtered_interactions\\t" >> {output.log}
    wc -l  {input.txt}  |perl -npe 's/ +/\t/' >> {output.log}
    echo -ne "Number_of_promoters\\t" >> {output.log}
    wc -l  {input.all_prom}  |perl -npe 's/ +/\t/' >> {output.log}
    echo -ne "Number_of_interacting_promoters\\t" >> {output.log} 
    wc -l  {output.txt}  |perl -npe 's/ +/\t/' >> {output.log}
    echo -ne "Number_of_enhProm_promoters\\t" >> {output.log} 
    wc -l  {input.enhProm}  |perl -npe 's/ +/\t/' >> {output.log}

    # Compute the number of enhProm falling into interacting promoters
    bedtools intersect -u -a {input.enhProm} -b {output.ip} > {output.eip}
    
    echo -ne "Number_of_enhProm_falling_in_interacting_promoters\\t" >> {output.log}
    wc -l  {output.eip}  |perl -npe 's/ +/\t/'     >> {output.log}    
    
    """
