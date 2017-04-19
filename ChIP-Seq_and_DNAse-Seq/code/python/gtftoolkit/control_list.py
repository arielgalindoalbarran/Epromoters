from __future__ import print_function

"""
    The ``control_list`` script
    ======================

    Based on a reference gene list, returns a list of genes matched for
    epigenetic signal in promoter, transcript body or tss. Note that
    it also works for any kind of signal (e.g RNA-Seq). In this case,
    the expression/signal values should be provided for all genes through
    the signalFile argument.
"""

import sys
import os
import time
import logging
import argparse

parser = argparse.ArgumentParser(add_help=True)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('more optional arguments')
        

def make_parser():
    """The program parser."""

    required.add_argument(
        '--in-file', '-i',
        metavar='File',
        help='A two columns tab-file with a non redondant list of ids '
             '(col 1, e.g transcript_id) and expression/signal values in column 2. No header.'
             '. This file should include the expression/signal for reference IDs '
             '(see --referenceGeneFile).'
             ' Any unwanted ID (e.g transcript sharing gene_id with '
             'the reference) should be discarded.',
        default=None,
        type=argparse.FileType('r'),
        required=True)
    
    required.add_argument(
        '--referenceGeneFile', '-r',
        metavar='File',
        help='The file containing the reference gene list (1 column,'
             ' transcript ids).'
             ' No header.',
        default=None,
        type=argparse.FileType('r'),
        required=True)
    
    optional.add_argument(
        '--out-dir', '-o',
        help='Name of the output directory.',
        type=str,
        default="control_list",
        required=False)
    
    optional.add_argument(
        '--log2', '-l',
        help='If selected, data will be log transformed.',
        action="store_true",
        required=False)

    optional.add_argument(
        '--pseudo-count', '-p',
        help='The value for a pseudo-count to be added.',
        type=float,
        default=1,
        required=False)

    optional.add_argument('-v', '--verbose',
                            help="Verbose if true.",
                            action='store_true')

    return parser

def command_log(log_file):

    logger = logging.getLogger(__name__)

    logging.basicConfig(filename=log_file,
                        level=logging.INFO,
                        format='%(asctime)s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    logger.info(" ".join(sys.argv))



def make_outdir_and_file(out_dir=None, alist=None, add_date=True):
    """Create output directory and a set of associated files (alist). Return a list of file handler (write mode) for the requested files."""

    if out_dir is not None:
        # Working directory
        if os.path.exists(out_dir):
            sys.stderr.write("Aborting. Directory already exist.")
            sys.exit()
        else:
            os.makedirs(out_dir)

    out_fh_list = list()

    timestr = time.strftime("%Y%m%d-%H%M%S")

    for fn in alist:

        fn_split = os.path.splitext(fn)
        if add_date:
            fn = fn_split[0] + "_" + timestr + fn_split[1]
        else:
            fn = fn_split[0] + fn_split[1]
        if out_dir is not None:
            out_fh_list.append(open(os.path.join(out_dir, fn), "w"))
        else:
            out_fh_list.append(open(fn, "w"))

    return out_fh_list


def control_list(in_file=None,
                 out_dir=None,
                 referenceGeneFile=None,
                 log2=False,
                 pseudo_count=1,
                 verbose=None):

    for p, line in enumerate(in_file):

        line = line.split("\t")

        try:
            float(line[1])
        except:
            print("It seems that column 2 of input file contains non numeric values.")
            print("Check that no header is present and that "
                  "columns are ordered properly.")
            sys.exit()

    if verbose:
        sys.stderr.write("Selecting genes (R call).")


    # Preparing pdf file name
    file_out_list = make_outdir_and_file(out_dir, ["control_list.txt",
                                                   "reference_list.txt",
                                                   "diagnostic_diagrams.pdf",
                                                   "R_code_control_list.R",
                                                   "command.log"])

    control_file, reference_file_out, pdf_file, r_code_file, log_file = file_out_list

    # Store the user command
    command_log(log_file.name)

    code_body = """

            ########################################
            # Load and/or install packages
            ########################################

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

        

            suppressWarnings(suppressMessages(load.fun("beanplot")))
            suppressWarnings(suppressMessages(load.fun("data.table")))
            suppressWarnings(suppressMessages(load.fun("ggplot2")))
            suppressWarnings(suppressMessages(load.fun("reshape2")))

            ########################################
            # Function declaration
            ########################################

            message <- function(msg){{
                cat(paste("    |--- ", msg, "\\n", sep=""), file=stderr())
            }}

            ########################################
            ## Get the list of reference genes
            ########################################
            
            reference_genes <- as.character(read.table('{reference_file}',
                                            header=F)[,1])

            ########################################
            ## Get expression data
            ########################################
            
            exp_data <- read.table('{signal_file}', sep="\\t", head =F)
            
            exp_data_vec <- exp_data[,2] + {pseudo_count}
            
            ########################################
            ## Log transformation
            ########################################            
            
            to.log <- '{log2}'
            
            if(to.log == 'True'){{
                if(length(exp_data_vec[exp_data_vec == 0])){{
                    message("Can't use log transformation on zero or negative values. Use -p. Exiting.")
                    q("no", 1, FALSE)
                }}
                exp_data_vec <- log2(exp_data_vec)
            }}
            
            names(exp_data_vec) <- exp_data[,1]

            # Now we have sorted expression data with gene/tx names.
            exp_data_vec <- sort(exp_data_vec)

            # T/F vector indicating which in the 
            # expression data list are found in reference_gene
            which_reference_genes <- names(exp_data_vec) %in% reference_genes            

            # convert the T/F vector to positions indicating wich position in
            # expression data is a reference gene/tx

            which_reference_genes <- which(which_reference_genes)

            message(paste("Found ", length(which_reference_genes),
                    " genes of the reference in the provided signal file", sep=""))
            
            not_found <- !(reference_genes %in% names(exp_data_vec))
            if(length(reference_genes[not_found]) > 0){{
                message(paste("List of reference genes not found :", reference_genes[not_found]))
            }}else{{
                message("All reference genes were found.")
            }}

            ########################################
            ## Search for gene with matched signal
            ########################################
            
            control_list<- c()

            nb.candidate.left <- length(exp_data_vec) -  length(which_reference_genes)
            
            if(nb.candidate.left < length(which_reference_genes) ){{
                message("Not enough element to perform selection. Exiting")
                q("no", 1, FALSE)
            }}

            cpt <- 1
            
            candidates <- exp_data_vec

            for(i in which_reference_genes){{
                p <- i
                not_candidate_pos <- unique(c(which_reference_genes, control_list))
                candidates[not_candidate_pos] <- NA
                diff <- abs(exp_data_vec[p] - candidates)
                control_list[cpt] <- which.min(diff)
                cpt <- cpt + 1
                
            }}


            write.table(cbind(exp_data_vec[which_reference_genes]),
                    "{reference_file_out}",
                    sep="\t", 
                    quote=F, 
                    col.names=NA)

            write.table(cbind(exp_data_vec[control_list]),
                    "{control_file}",
                    sep="\t", 
                    quote=F, 
                    col.names=NA)
                    

            message("Preparing diagnostic plots.")
            m <- as.data.frame(cbind(
                                exp_data_vec[control_list], 
                                exp_data_vec[which_reference_genes]
                )
            )
                            
            colnames(m) <- c("Control", "Reference")
            
            
            pdf("{pdf_file}")
            
            # Preparing beanplot (side=both)
            beanplot(m, 
                    col = list("blue", "darkgrey"), 
                    border="white", 
                    log="",
                    ll=0.13,
                    what=c(0,1,0,1),
                    side="both"
            )
            grid(col="grey", nx=0, ny=NULL)
            beanplot(m, 
                    col = list("blue", "darkgrey"), 
                    border="white", 
                    log="", main="Beanplots of Control and Reference values",
                    ll=0.13,
                    what=c(0,1,0,1),
                    add=T,
                    side="both"
            )
            
            # Preparing beanplot (side=default)
            beanplot(m, 
                    col = list("blue", "darkgrey"), 
                    border="white", 
                    log="",
                    ll=0.13,
                    what=c(0,1,0,1)
            )
            grid(col="grey", nx=0, ny=NULL)
            beanplot(    m, 
                    col = list("blue", "darkgrey"), 
                    border="white", 
                    log="", main="Beanplots of Control and Reference values",
                    ll=0.13,
                    what=c(0,1,0,1),
                    add=T
            )
            
            # Preparing boxplot
            boxplot(m, col=c("blue","darkgrey"))
            grid(col="grey", nx=0, ny=NULL)

            # Preparing qqplot
            plot(    sort(m$Control),
                    sort(m$Reference), 
                    pch=16, xlab="Control", 
                    ylab="Reference", 
                    main="QQplot of control and reference values",  
                    panel.first=grid()
            )
            
            # Preparing histograms
            
            m.m <- melt(m, id.vars = c(NULL))
            col <- m.m$variable
            levels(col) <- c("blue","darkgrey")
            p <- ggplot(m.m, 
                        aes(    x=value, 
                                color=variable, 
                                fill=variable))
            p <- p + geom_bar(position="dodge")
            p <- p + scale_fill_manual(values=levels(col))
            p <- p + scale_color_manual(values=levels(col))
            p <- p + theme(legend.title=element_blank())
            
            suppressMessages(print(p))
            
            out <- dev.off()
    """.format(signal_file=in_file.name,
               reference_file=referenceGeneFile.name,
               log2=log2,
               pseudo_count=str(pseudo_count),
               pdf_file=pdf_file.name,
               control_file=control_file.name,
               reference_file_out=reference_file_out.name)

    if verbose:
        
        sys.stderr.write("Printing R code to: " + r_code_file.name)

    r_code_file.write(code_body)
    r_code_file.close()

    if verbose:
        sys.stderr.write("Executing R code.")

    # Execute R code.
    os.system("cat " + r_code_file.name + "| R --slave")


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    control_list(**args)

if __name__ == '__main__':
    main()