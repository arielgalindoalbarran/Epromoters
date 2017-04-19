rule stats_on_chiapet:
    input: "DATA/{exp}.tsv"
    output: "RESULTS/STAT_ON_CHIAPET/{exp}.pdf"
    params: ppn="nodes=1:ppn=1"
    run: R("""
    d <- read.table("{input}", header=TRUE)
    len.1 <- log2(d$stop1-d$start1)
    len.2 <- log2(d$stop2-d$start2)
    
    pdf("{output}")
    
    par(mfrow=c(2,3))
    hist(log2(len.1), main="Distribution of log2(fragment len) (1)")
    hist(log2(len.2), main="Distribution of log2(fragment len) (2)")
    
    hist(d$Z, main="Distribution of log2(Z) values")

    hist(d$IAB, main="Distribution of IAB values")
    hist(log2(d$IAB), main="Distribution of log2(IAB) values")
     
    plot(0)
    
    hist(d$CA, main="Distribution of CA values")
    hist(log2(d$CA), main="Distribution of log2(CA) values")
    
    hist(log2(d$CB), main="Distribution of log2(CB) values")
    hist(d$CB, main="Distribution of log2(CB) values")
    dev.null <- dev.off()
    """)