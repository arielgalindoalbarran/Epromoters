rule compute_proba_prom_prom:
    input: "RESULTS/PROM_PROM/{cond}-{exp}_fltr_pp_chiapet_overlap_prom.log"
    output: "RESULTS/PROM_PROM/{cond}-{exp}_fltr_pp_chiapet_overlap_prom.stat"
    params: ppn="nodes=1:ppn=1"
    run: R("""
    Sys.sleep(6)
    data <- read.table('{input}', head=F, row=1)
    nbp <- data["Number_of_promoters",1]
    nbpi <- data["Number_of_interacting_promoters",1]
    p <- nbpi / nbp
    nep <- data["Number_of_enhProm_promoters",1]
    nepi <- data["Number_of_enhProm_falling_in_interacting_promoters",1]
    res <- pbinom(lower.tail=FALSE, prob= p, size=nep, q=nepi-1)
    cat('{output}\n', file='{output}')
    cat(paste("Binomial_test\\t", res, "\n", sep=""), file = "{output}", append=T)
    """)