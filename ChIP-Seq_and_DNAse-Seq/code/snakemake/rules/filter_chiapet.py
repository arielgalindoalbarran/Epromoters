rule filter_chiapet:
    input: "DATA/{exp}.tsv"
    output: "RESULTS/FILTER_ON_CHIAPET/{exp}_fltr.txt"
    threads:
        1
    run: R("""
    d <- read.table("{input}", header=TRUE)
    write.table(d, "{output}", sep="\t", quote=F, col.names=NA)
    """)
