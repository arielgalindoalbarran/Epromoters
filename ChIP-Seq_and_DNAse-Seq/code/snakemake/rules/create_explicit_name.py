"""
"""

rule create_explicit_name_part1_r:
    """
    """
    input:  tsv='data/Encode/{strain}/epimarks/metadata.tsv',\
            path='data/Encode/{strain}/epimarks/'
    output: ln_command='data/Encode/{strain}/epimarks/ln_command_to_eval.txt'
    params: ppn="nodes=1:ppn=1"
    run:
        R("""
        metadata <- read.delim("{input.tsv}")
        
        encode_name <- paste(
            metadata$File.accession,
            metadata$File.format,
            sep=".")

        path_encode_name <- paste(
            "{input.path}",
            encode_name,
            sep="")

        explicit_name <- paste(
            metadata$Assay,
            metadata$Biosample.term.name,
            metadata$Experiment.target,
            metadata$File.accession,
            sep="_")

        explicit_name <- paste(
            explicit_name,
            metadata$File.format,
            sep=".")

        path_explicit_name <- paste(
            "{input.path}",
            explicit_name,
            sep="")

        ln_command <- 'ln -s'
        merge_names <- cbind(ln_command,encode_name,path_explicit_name)

        write.table(x = merge_names, file = "{output.ln_command}", col.names = F, row.names = F, quote = F, sep = " ")
    """)

rule create_explicit_name_part2_soft_link:
    """
    """
    input:  tsv='data/Encode/{strain}/epimarks/metadata_names.tsv'
    output: done='data/Encode/{strain}/epimarks/soft_link.done'
    params: ppn="nodes=1:ppn=1"
    shell:"""
    touch {output.done} 
    """


