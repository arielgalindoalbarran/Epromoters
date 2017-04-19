rule get_prom:
    input: bed="DATA/allResults_FC_unix.bed", chrinfo="DATA/ChromInfo.txt"
    output: "DATA/allResults_FC_unix_1k_1k.bed"
    params: ppn="nodes=1:ppn=1"
    shell: """
    sleep 3
    bedtools slop -s -l 750 -r 950 -i {input.bed} -g {input.chrinfo} > {output}
    """