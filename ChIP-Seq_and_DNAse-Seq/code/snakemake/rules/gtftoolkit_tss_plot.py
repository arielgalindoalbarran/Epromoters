rule control_list_to_gtf:
    """
    This rule creates gtf file for control and reference files produced by 'gtftoolkit control_list' in order to use them for 'gtftoolkit tss_plot'.

    Note that we use the gtf made for broad domain to get coordinates because they are already centered (-1kb+1kb) around TSS.
    """
    input:
        control_list="data/reference_and_control/{strain}_{control_or_reference}_list.txt",
        gtf="annotation/RefSeqGenes_hg19_mainChr.gtf"
    output:
        gtf="results/control_list_to_gtf/{strain}_{control_or_reference}_list.gtf"
    params:
        ppn="nodes=1:ppn=1"
    shell:
        """
        rm -f {output.gtf}
        for i in `cut -f1 {input.control_list}`
        do
            grep $i {input.gtf} >> {output.gtf}
        done
        """

rule gtftoolkit_tss_plot_v3:
    """
    Rule to produce TSS plots for Salva and Ariel project.
    Now you can specify upstream, downstream et window size to make -1kb+1kb plots as requested by Salva.

    Usage:
    expand("result/tss_plot/{broad_or_tf}/{strain}/{mark}/{accession_id}/{control_or_reference}/transform_{transform}", broad_or_tf ....)
    u5000d5000w200
    u1000d1000w40
    """
    input:
        gtf="results/control_list_to_gtf/{strain}_{control_or_reference}_list.gtf",
        chrominfo="annotation/ChromInfo.txt",
        bw="data/Encode/{strain}/{mark}/{accession_id}.bigWig"
    output:
        output_dir="results/tss_plot_v3/{strain}/{mark}/{accession_id}/{control_or_reference}/transform_{transform}/u{upstream}d{downstream}w{window}"
    params:
        ppn="nodes=1:ppn=15"
    shell:
        """
        gtftoolkit tss_plot \
            -k 15 -i {input.gtf} -c {input.chrominfo} -v -n 4 -T {wildcards.transform} \
            -l {input.bw} \
            --labels {wildcards.strain}_{wildcards.mark}_{wildcards.accession_id}_{wildcards.control_or_reference} \
            -o {output.output_dir} \
            -u {wildcards.upstream} -d {wildcards.downstream} -w {wildcards.window}
        """

rule gtftoolkit_tss_plot:
    """
    Rule to produce TSS plots for Salva and Ariel project.

    Usage:
    expand("result/tss_plot/{broad_or_tf}/{strain}/{mark}/{accession_id}/{control_or_reference}/transform_{transform}", broad_or_tf ....)
    """
    input:
        gtf="results/control_list_to_gtf/{broad_or_tf}/{strain}/{mark}/{accession_id}/{control_or_reference}.gtf",
        chrominfo="annotation/ChromInfo.txt",
        bw="data/Encode/{strain}/{mark}/{accession_id}.bigWig" 
    output:
        output_dir="results/tss_plot/{broad_or_tf}/{strain}/{mark}/{accession_id}/{control_or_reference}/transform_{transform}"
    params:
        ppn="nodes=1:ppn=15"
    shell:
        """
        gtftoolkit tss_plot \
            -k 15 -i {input.gtf} -c {input.chrominfo} -v -n 4 -T {wildcards.transform} \
            -l {input.bw} \
            --labels {wildcards.strain}_{wildcards.mark}_{wildcards.accession_id}_{wildcards.control_or_reference} \
            -o {output.output_dir} \
            -u 5000 -d 5000 -w 200
        """

rule merge_control_and_reference_tss_plot:
    input:
        input_dir_ctr="results/tss_plot_v3/{strain}/{mark}/{accession_id}/control/transform_{transform}/u{upstream}d{downstream}w{window}",
        input_dir_ref="results/tss_plot_v3/{strain}/{mark}/{accession_id}/reference/transform_{transform}/u{upstream}d{downstream}w{window}"
    output:
        merged_tss_plot="results/merge_control_and_reference_tss_plot/{strain}_{mark}_{accession_id}_merged_tss_plot_transform_{transform}_u{upstream}d{downstream}w{window}.pdf"
    params: ppn="nodes=1:ppn=1"
    run: 
        R("""
        load.fun <- function(x) {{
            x <- as.character(substitute(x)) 
            if(isTRUE(x %in% .packages(all.available=TRUE))) {{ 
                eval(parse(text=paste("require(", x, ")", sep=""))) 
            }} else {{ 
                eval(parse(text=paste("install.packages('", x, "', 
                repos = 'http://cran.us.r-project.org')", sep=""))) 
                eval(parse(text=paste("require(", x, ")", sep=""))) 
            }} 
        }}
        
        load.bioc <- function(x) {{ 
            x <- as.character(substitute(x)) 
            if(isTRUE(x %in% .packages(all.available=TRUE))) {{ 
                eval(parse(text=paste("require(", x, ")", sep=""))) 
            }} else {{ 
                eval(parse(text="source('http://bioconductor.org/biocLite.R'"))
                eval(parse(text=paste("biocLite('", x, "')", sep="")))
            }} 
        }} 
        
        
        suppressWarnings(suppressMessages(load.fun("reshape2")))
        suppressWarnings(suppressMessages(load.fun("ggplot2")))
        suppressWarnings(suppressMessages(load.fun("colorRamps")))
        suppressWarnings(suppressMessages(load.fun("data.table")))
    
    
        ########################################
        # Function declaration
        ########################################
        
        transformData <- function(expr, mark, add.min=F, to.log=True){{
            mark <- as.character(mark)
            for(k in unique(mark)){{
                if(to.log){{
                    tmp <- log2(expr[mark==k])
                }}else{{
                    tmp <- expr[mark==k]
                }}
                med <- median(tmp)
                tmp.cent <- tmp-med
                mad.v <- mad(tmp)
                tmp.cent.mad <- tmp.cent/mad.v
                if(add.min){{
                    min.v <- min(tmp.cent.mad)
                    expr[mark==k] <- tmp.cent.mad - min.v
                }}else{{
                    expr[mark==k] <- tmp.cent.mad
                }}
                
            }}
            return(expr)
        }}
        
        message <- function(msg){{
            cat(paste("    |--- ", msg, "\n", sep=""))
        }}
        
        ########################################
        # Load dataset
        ########################################
        
        
        ### Reference
        path_ref <- "{input.input_dir_ref}"
        file_ref <- list.files(path=path_ref, pattern="coverage_matrix*")
        d_ref <- as.data.frame(
            fread(
                paste(
                    path_ref,
                    file_ref,
                    sep="/"),
                sep='\t', 
                header=F)
            )
        
        
        ### Control
        path_ctr <- "{input.input_dir_ctr}"
        file_ctr <- list.files(path=path_ctr, pattern="coverage_matrix*")
        d_ctr <- as.data.frame(
            fread(
                paste(
                    path_ctr,
                    file_ctr,
                    sep="/"),
                sep='\t', 
                header=F)
            )
        
        ########################################
        # Variables
        ########################################
        
        nb.samples <- length(table(d_ref[,1]))
        nb.regions <- table(d_ref[,1])[1]
        col <-c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf')[1:nb.samples]
        pdf_file <- "{output.merged_tss_plot}"
        from = -5000
        to   <- 5000
        bin_nb <- 200
        bin_nb_dws <-  0
        bin_nb_ups <-  0
        bin_nb_total <- bin_nb_ups + bin_nb + bin_nb_dws
        nb_class <- 4
        ft_type <- 'promoter'
        transform <- 'none'
        split_mode <- 'exprs'
        distance <- 'pearson'
        q.show <- 'False'
        bw_to_split <- c('H4K5ac')
        if(bw_to_split == ''){{
                bw_to_split <- c()
        }}
        
           
        color.ramp <- colorRampPalette(c('#4575b4','#74add1','#abd9e9','#e0f3f8','#fee090','#fdae61','#f46d43','#d73027'))(10)
        
        ########################################
        # Overlaid profiles
        ########################################
        
        # Rename the columns
        colnames(d_ref)[1:6] <- c("file", "chrom", "start", "end", "gene", "strand")
        colnames(d_ctr)[1:6] <- c("file", "chrom", "start", "end", "gene", "strand")
        
        
        # Reshape
        suppressMessages(dm_ref <- melt(d_ref, id=colnames(d_ref)[1:6]))
        colnames(dm_ref) <- c("mark", "chrom", "start", "end", "gene", "strand", "pos", "exprs")
        dm_ref$group <- 'Active'
        
        suppressMessages(dm_ctr <- melt(d_ctr, id=colnames(d_ctr)[1:6]))
        colnames(dm_ctr) <- c("mark", "chrom", "start", "end", "gene", "strand", "pos", "exprs")
        dm_ctr$group <- 'Control'
        
        dm_merge <- rbind(dm_ref,dm_ctr)
    
    
        from <- -5000
        to   <- 5000
        bin_nb <- 200
    
        levels(dm_merge$pos) <- seq(from=from,
            to=to,
            length.out=bin_nb)
            dm_merge$pos <- as.numeric(as.character(dm_merge$pos))
    
    
        write.table(dm_merge, "testTable.txt") 
        pdf("justAtest.pdf")
        plot(x=c(1,2,3))
        dev.off()
    
    
        df <- data.frame(gp = factor(rep(letters[1:3], each = 10)),
                     y = rnorm(30))
                     p <- ggplot(df, aes(x = gp, y = y)) +
                        geom_point()
        pdf("anotherTest.pdf")
        print(p)
        dev.off()
    
        p <- ggplot(dm_merge, aes(x=pos,y=exprs, group=group, colour = group))
        p <- p + stat_summary(fun.y = mean, geom="line", size=0.8)
        #p <- p + scale_color_manual(values=levels(as.factor(col)))
        #p <- p + theme(legend.title=element_blank())
        p <- p + theme_bw()
        # p <- p + theme(legend.text=element_text(size=6))
        #p <- p + theme(legend.title=element_blank())
        #p <- p + theme(legend.position = "bottom",
        #    legend.key = element_rect(colour = "white"))
        #p <- p + theme(axis.text.x = element_text(colour="grey20",
        #    size=5,angle=90,
        #    face="plain"))
        #p <- p + theme(axis.title.x = element_text(colour="grey20",
        #    size=12,angle=0,
        #    face="plain"))
        #p <- p + theme(axis.title.y = element_text(colour="grey20",
        #    size=10,angle=90,
        #    hjust=.5,vjust=1,
        #    face="plain"))
        ## p <- p + theme(axis.title.x=element_blank())
        ## p <- p + theme(axis.text.x = element_blank())
        #p <- p + ylab(ylab)
        #p <- p + xlab(xlab)
        #p <- p + ggtitle("Classes of gene expression defined by quartiles")
        #p <- p + theme(axis.text.x = element_text(angle = 35, size=8))
        
    
        pdf(pdf_file)
        print(p)
        dev.off()
        
        """)

