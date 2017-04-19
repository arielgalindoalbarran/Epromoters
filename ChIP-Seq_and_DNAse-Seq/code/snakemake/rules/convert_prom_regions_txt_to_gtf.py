rule convert_prom_regions_txt_to_gtf:
    """
    Aim:
        Convert promoter regions from Salva tab file into gtf for use with gtftoolkit coverage.
    Notes:
        Quantify marks with the list of all hPro.
        PolII et Rad21: -200 + 50
        Ne pas utiliser le gtf mais all_hpromoterRegions.
        We select only the first transcript as ID:
        sed is used because return characters are conserved in All_hpromoterRegions_bed4.txt
    """
    input:
        txt='data/All_hpromoterRegions.txt',
        chrominfo='annotation/ChromInfo.txt',
        gtftoolkit='code/python/gtftoolkit/bed_to_gtf.py'
    output:
        gtfFT='data/All_hpromoterRegions_m200p50_geneID.gtf',
        gtfBroad='data/All_hpromoterRegions_m1000p1000_geneID.gtf'
    params:
        ppn="nodes=1:ppn=1"
    shell:
        """
        cut --fields 7 {input.txt} > {input.txt}_bed4.txt

        cut --fields 1-3 {input.txt} > {input.txt}_bed1-3.txt

        cut --fields 5-6 {input.txt} > {input.txt}_bed5-6.txt

        # We add the first transcript ID as bed names.
        # We create Bed from all_hpromoterRegions.
        paste --delimiters='\t' \
                {input.txt}_bed1-3.txt \
                {input.txt}_bed4.txt \
                {input.txt}_bed5-6.txt | \
                tail -n +2 > \
                {input.txt}.bed
                
        {input.gtftoolkit} \
                --infile {input.txt}.bed \
                --outfile {output.gtfFT}

        bedtools slop \
            -i {input.txt}.bed \
            -g {input.chrominfo} \
            -l 800 -r 950 > \
            {input.txt}_m1000p1000.bed

        {input.gtftoolkit} \
            --infile {input.txt}_m1000p1000.bed \
            --outfile {output.gtfBroad}
        """
