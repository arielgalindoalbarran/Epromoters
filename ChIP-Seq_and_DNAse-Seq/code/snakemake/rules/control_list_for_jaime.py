rule control_list_for_jaime:
    """
    Rules to get control list for e-promoters active both in Hela and K562, as requested by Jaime. 
    2015-12-09 15:45

    cut -f1 data/reference_and_control/hela_control_list.txt data/reference_and_control/K562_control_list.txt
    """
    input:  control_list_hela="data/reference_and_control/hela_control_list.txt",\
            control_list_k562="data/reference_and_control/K562_control_list.txt",\
            reference_list_hela="data/reference_and_control/hela_reference_list.txt",\
            reference_list_k562="data/reference_and_control/K562_reference_list.txt"
    output: reference_common_epromoters='results/control_list_for_jaime/reference_common_epromoters.txt', \
            control_based_on_k562='results/control_list_for_jaime/control_based_on_k562.txt',\
            control_based_on_hela='results/control_list_for_jaime/control_based_on_hela.txt',\
            reference_based_on_k562='results/control_list_for_jaime/reference_based_on_k562.txt',\
            reference_based_on_hela='results/control_list_for_jaime/reference_based_on_hela.txt'
    params: ppn="nodes=1:ppn=1"
    shell:"""
    cut -f1 {input.reference_list_hela} {input.reference_list_k562} | sort | uniq -d > {output.reference_common_epromoters}
    
    rm -f linesToSelectInRefListHela linesToSelectInRefListK562 \
            {output.control_based_on_k562} {output.control_based_on_hela} \
            {output.reference_based_on_k562} {output.reference_based_on_hela}

    for i in `cat {output.reference_common_epromoters}`
        do
            # -n argument add the line number of match. It is separated from the content of the line with ":".
            grep -n -E "^${{i}}\t" {input.reference_list_hela} | cut -f1 --delimiter ":" >> linesToSelectInRefListHela
            grep -n -E "^${{i}}\t" {input.reference_list_k562} | cut -f1 --delimiter ":" >> linesToSelectInRefListK562
        done

    for i in `cat linesToSelectInRefListHela`
        do
            head -$i {input.reference_list_hela} | tail -1 >> {output.reference_based_on_hela}
            head -$i {input.control_list_hela} | tail -1 >> {output.control_based_on_hela}
        done

    for i in `cat linesToSelectInRefListK562`
        do
            head -$i {input.reference_list_k562} | tail -1 >> {output.reference_based_on_k562}
            head -$i {input.control_list_k562} | tail -1 >> {output.control_based_on_k562}
        done

    rm -f linesToSelectInRefListHela linesToSelectInRefListK562
    """

rule control_list_for_jaime_plots:
    """
    Some plots to check what is produced by "control_list_for_jaime" rule.
    """
    input:  control_based_on_k562='results/control_list_for_jaime/control_based_on_k562.txt',\
            control_based_on_hela='results/control_list_for_jaime/control_based_on_hela.txt',\
            reference_based_on_k562='results/control_list_for_jaime/reference_based_on_k562.txt',\
            reference_based_on_hela='results/control_list_for_jaime/reference_based_on_hela.txt'
    output: superposition_points_ctr_and_ref='results/control_list_for_jaime/superposition_points_ctr_and_ref.pdf',\
            superposition_lines_ctr_and_ref='results/control_list_for_jaime/superposition_lines_ctr_and_ref.pdf',\
            correlation_k562_and_hela='results/control_list_for_jaime/correlation_k562_and_hela.pdf'
    params: ppn="nodes=1:ppn=1"
    run:R("""
    control_based_on_k562 <- read.delim("{input.control_based_on_k562}", header=FALSE)

    reference_based_on_k562 <- read.delim("{input.reference_based_on_k562}", header=FALSE)

    control_based_on_hela <- read.delim("{input.control_based_on_hela}", header=FALSE)

    reference_based_on_hela <- read.delim("{input.reference_based_on_hela}", header=FALSE)

    pdf("{output.superposition_points_ctr_and_ref}")
    plot(control_based_on_k562$V2, pch=1, col=1,xlab="Index (e-promoter or common to Hela and K562 or control)",ylab="Expression (ask Ariel for details)")
    points(reference_based_on_k562$V2,pch=2, col=2)
    points(control_based_on_hela$V2,pch=3, col=3)
    points(reference_based_on_hela$V2,pch=4, col=4)
    legend("bottomright",legend=c("ctr K562","ref K562","ctr Hela","ref Hela"), pch=seq(1,4), col=seq(1,4))
    out <- dev.off()

    pdf("{output.superposition_lines_ctr_and_ref}")
    plot(control_based_on_k562$V2, pch=1, col=1, type="l", xlab="Index (e-promoter common to Hela and K562 or control)",ylab="Expression (ask Ariel for details)")
    lines(reference_based_on_k562$V2,pch=2, col=2)
    lines(control_based_on_hela$V2,pch=3, col=3)
    lines(reference_based_on_hela$V2,pch=4, col=4)
    legend("bottomright",legend=c("ctr K562","ref K562","ctr Hela","ref Hela"), pch=seq(1,4), col=seq(1,4))
    out <- dev.off()

    pdf("{output.correlation_k562_and_hela}")
    plot(x=control_based_on_k562$V2, y=control_based_on_hela$V2, pch=1, col=1)
    out <- dev.off()
    """)

