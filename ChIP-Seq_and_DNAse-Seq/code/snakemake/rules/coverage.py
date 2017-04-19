"""
Rules to compute coverage around TSS with different settings for FT and Broad domains.
"""

rule coverage_tf:
    """
    Get coverage in promoter region for FT.
    
    Usage: 
    expand('RESULTS/coverage_ft/{mark}_{strain}.tab', \
            mark=['POLR2A', 'RAD21'], \
            strain=['k562', "hela"])
    """
    input:
        gtf='data/All_hpromoterRegions_m200p50_geneID.gtf',
        bw='data/Encode/{strain}/{mark}/{accession_id}.bigWig',
        chrominfo='annotation/ChromInfo.txt',
        gtftoolkit_coverage='code/python/gtftoolkit/bw_coverage.py',
        gtftoolkit_convert='code/python/gtftoolkit/convert.py',
        gtftoolkit_col_from_tab='code/python/gtftoolkit/col_from_tab.py'
    output:
        tab='results/coverage_tf/{strain}/{mark}/{accession_id}.tab'
    params:
        ppn="nodes=1:ppn=1"
    shell:
        """
        {input.gtftoolkit} --verbose \
            --infile {input.gtf} \
            --type transcript \
            --coveragefileList {input.bw} \
            --chromInfo {input.chrominfo} \
            --score b0 --nbThreads 1 | \
            {input.gtftoolkit_convert} \
                --format tx_tab --verbose | \
                {input.gtftoolkit_col_from_tab} \
                    --columns transcript_id,transcript_cov > \
                    {output.tab}
        """


rule coverage_broad:
    """
    Get coverage in promoter region for broad domain.
    The difference is that we use only top25% on the region -1000+1000 to define the coverage.
    """
    input:
        gtf='data/All_hpromoterRegions_m1000p1000_geneID.gtf',
        bw='data/Encode/{strain}/{mark}/{accession_id}.bigWig',
        chrominfo='annotation/ChromInfo.txt',
        gtftoolkit_coverage='code/python/gtftoolkit/bw_coverage.py',
        gtftoolkit_convert='code/python/gtftoolkit/convert.py',
        gtftoolkit_col_from_tab='code/python/gtftoolkit/col_from_tab.py'
    output:
        tab='results/coverage_broad/{strain}/{mark}/{accession_id}.tab'
    params:
        ppn="nodes=1:ppn=1"
    shell:
        """
        {input.gtftoolkit_coverage} --verbose \
            --infile {input.gtf} \
            --type transcript \
            --coveragefileList {input.bw} \
            --chromInfo {input.chrominfo} \
            --binNb 20 --nHighest 5 \
            --score b0 --nbThreads 1 | \
            {input.gtftoolkit_convert} \
                --format tx_tab --verbose | \
                {input.gtftoolkit_col_from_tab} \
                    --columns transcript_id,transcript_cov > \
                    {output.tab}
        """

rule coverage_summary:
    """
    Create the summary table requested by Salva
    """
    input:
        expand('RESULTS/coverage_{ft_or_broad}/{mark}_{strain}.tab',\
            zip, ft_or_broad=["ft", "ft", "broad", "broad"], \
            mark=["RAD21", "POLR2A", "H3K27ac", "H3K4me3"], \
            strain=["K562"]*4)
    output:
        'RESULTS/coverage_summary/coverage_summary.tsv'
    params:
        ppn="nodes=1:ppn=1"
    shell:"""
    paste \
        RESULTS/coverage_ft/RAD21_K562.tab \
        RESULTS/coverage_ft/POLR2A_K562.tab \
        RESULTS/coverage_broad/H3K4me3_K562.tab \
        RESULTS/coverage_broad/H3K27ac_K562.tab > {output}
    """


