"""
Rules to get control list from coverage files and salva's list of inactive promoters
"""
#
#rule get_inactive_prom:
#    """
#    A rule to convert table containing inactive promoters into list to use with gtftoolkit control_list.
#
#    Usage: Used by rule gtftoolkit_control_list
#    """
#    input: "DATA/SalvaWork/CapStarrseq_InactiveProm_FDR95_All_samples.txt"
#    output: "DATA/SalvaWork/CapStarrseq_InactiveProm_FDR95_All_samples_onlyFirstTranscriptId.txt"
#    params: ppn="nodes=1:ppn=1"
#    shell:"""
#    cut --fields 8 {input} | cut --fields 1 --delimiter " " | \
#            sed 's/\r//g' > {output}
#    """
#
#rule get_active_prom:
#    """
#    A rule to convert table containing active promoters into list to use with gtftoolkit control_list.
#
#    Usage: Used by rule gtftoolkit_control_list
#    """
#    input: "DATA/SalvaWork/CapStarrseq_Active_Prom_K562_merge_IP.txt"
#    output: "DATA/SalvaWork/CapStarrseq_Active_Prom_K562_merge_IP_onlyFirstTranscriptId.txt"
#    params: ppn="nodes=1:ppn=1"
#    shell:"""
#    cut --fields 8 {input} | cut --fields 1 --delimiter " " | \
#            sed 's/\r//g' > {output}
#    """
#

rule gtftoolkit_control_list:
    """
    Get control listcoverage in promoter region.
    
    Usage:
    expand('RESULTS/control_list_{ft_or_broad}/{mark}_{strain}/done',
            zip, ft_or_broad=["ft", "ft", "broad", "broad"],
            mark=["RAD21", "POLR2A", "H3K27ac", "H3K4me3"],
            strain=["K562"]*4),
    """
    input:
        infile='results/coverage_{tf_or_broad}/{strain}/{mark}/{accession_id}.tab',
        ref='data/reference_list_{strain}.txt',
        gtftoolkit='code/python/gtftoolkit/control_list.py'
    output:
        control_list='results/control_list_{tf_or_broad}/{strain}/{mark}/{accession_id}/control_list.txt',
        path="results/control_list_{tf_or_broad}/{strain}/{mark}/{accession_id}"
    wildcard_constraints:
        tf_or_broad='tf|broad',
        strain='hela|K562',
        mark="RAD21|POLR2A|H3K27ac|H3K4me3|tfactors|dnase|faire|epimarks"
    conda:
        "../envs/gtftoolkit.yaml"
    threads:
        1
    shell:"""
    # We remove the content in the output directory first because gtftoolkit control_list do not want to overwrite files.
    rm -rf {output.path}

    python {input.gtftoolkit} \
            --in-file {input.infile} \
            --referenceGeneFile {input.ref} \
            --out-dir {output.path}

    # gtftoolkit produces output files with random ID so we copy the file of interest to a name compliant with snakemake requirement.
    cp {output.path}/control_list_*.txt {output.control_list}
    """
