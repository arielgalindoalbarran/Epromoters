from __future__ import print_function
from __future__ import division
import os
import re
import sys
import time
import argparse
import pandas as pd
from tempfile import NamedTemporaryFile
import multiprocessing
import shutil
from collections import defaultdict
from bx.bbi.bigwig_file import BigWigFile
from pybedtools import BedTool
from itertools import repeat
import numpy as np
import logging
import shutil


def silentremove(filename):
    import os
    import errno
    try:
        os.remove(filename)
    except OSError as e:
        pass


def make_tmp_file(prefix='', store=True, suffix=''):
    """ Create a temporary file. This file will be deleted after a call to the gtftoolkit command.
    """
    tmp_file = NamedTemporaryFile(delete=False, prefix=prefix, suffix=suffix)
    if store:
        TMP_FILE_LIST.append(tmp_file.name)
    return tmp_file

TMP_FILE_LIST = []

"""
Make a coverage plot (e.g around TSS) with one or several bigWig.
-------------------------------------------------------

Use bx-python as interface to the kent utilities.

"""


def intervals(l, n):
    """Returns a list of tuple that correspond to n ~ equally spaced intervals
    between 0 and l.
    >>> l = range(10)
    >>> intervals(l, 3)
    >>> [(0, 3), (3, 6), (6, 10)]
"""
    result = list()

    def chunks(l, n):
        """ Yield n successive chunks from l.
        """
        newn = int(len(l) / n)
        for i in xrange(0, n - 1):
            yield l[i * newn:i * newn + newn + 1]
        yield l[n * newn - newn:]

    for i in chunks(l, n):
        result.append((i[0], i[len(i) - 1]))

    result[len(result) - 1] = (result[len(result) - 1][0],
                               result[len(result) - 1][1] + 1)
    return result


def message(msg, nl=True):
    msg = str(msg)
    """Send a message to sdterr."""
    if nl:
        sys.stderr.write("    |--- " + msg + "\n")
    else:
        sys.stderr.write("    |--- " + msg)


def big_wig_summary_worker(xxx_todo_changeme):
    """
    This function compute bigwig coverage. 'span' is a tuple that correspond to a
    fraction (from, to) of the bedfile to be processed. The bedfile is produced through
    the bw_coverage function (promoters, tx,...). Each worker will process all bigwig files
    but it will only process a fraction (span) of the bed file regions

    * span: the fraction (lines) of the bed file [from, to] to be processed.
    * bw_list: the list of bigWig files to be processed.
    * region_bed_file_name: the bed file containing the region for which coverage is to be computed.
    * bin_nb: The number of bin into which the region should be splitted.
    * If the number of nucleotides is < nbBin. Estimate the values with the mean.
    * Add a pseudo-count.
    * n_highest: compute the score based on the n highest values in the bins.
    * profile: don't return a single value but a vector (profile) for each regions.
    * verbose: run in verbose mode.
    """
    (span,
     bw_list,
     region_bed_file_name,
     bin_nb,
     pseudo_count,
     n_highest,
     profile,
     verbose) = xxx_todo_changeme
    if n_highest is None:
        n_highest = bin_nb

    results = list()

    # The result will be a list of file if profile
    # is True
    if profile:
        matrix_file_list = list()

    for big_wig, cpt in zip(bw_list, range(len(bw_list))):

        bigwig = BigWigFile(open(big_wig, "r"))

        if verbose:

            message(
                "Computing coverage for file: " +
                big_wig +
                " [" +
                str(multiprocessing.current_process()) +
                "], " +
                str(span[1] - span[0]) + " chunks to process.")

        # Load the regions for which the coverage is to be processed.

        tx_bed = BedTool(region_bed_file_name)

        # There is no  current chrom at the moment
        chr_cur = None

        # The fraction of bed file
        from_here = span[0]
        to_here = span[1]

        if profile:
            matrix_file = make_tmp_file(prefix=os.path.basename(big_wig),
                                        store=False)
            if verbose:
                message("Writing to: " + matrix_file.name)

        for i in tx_bed[slice(from_here, to_here)]:

            if chr_cur is None:
                chr_cur = i.chrom
            else:
                if i.chrom != chr_cur:
                    chr_cur = i.chrom

            # Note: bigWig is zero-based/half open as bed.

            # If end - start is greater than bin_nb
            # bigwig.summarize will return a set of Zero.
            # No way to do it.
            if (i.end - i.start) < bin_nb:
                message(i.name +
                        " : the number of nucleotides in feature (" +
                        str(i.end - i.start) +
                        " is below the window length (" + str(bin_nb) +
                        "). " +
                        "Skipping")
                continue

            bw_summary = bigwig.summarize(i.chrom,
                                          i.start,
                                          i.end,
                                          bin_nb)

            if bw_summary is None:

                if verbose:
                    message(i.name + ": no coverage info found (check if chrom, " +
                            i.chrom + " , is " +
                            "declared in bigWig file.")

                # Should be replaced with NA
                bw_summary = np.repeat(pseudo_count, bin_nb)

            else:
                bw_summary = bw_summary.sum_data / ((i.end - i.start) / bin_nb)
                bw_summary = bw_summary + pseudo_count

            # Replace nan by pseudo_count
            # Should be replaced with NA
            nan_pos = np.isnan(bw_summary)

            bw_summary[nan_pos] = pseudo_count

            # Data should be oriented in 5' -> 3'
            if i.strand == '-':
                bw_summary = bw_summary[::-1]

            if profile:

                bw_name_short = os.path.splitext(os.path.basename(big_wig))[0]
                out_text = [bw_name_short,
                            i.chrom,
                            str(i.start),
                            str(i.end),
                            str(i.name),
                            i.strand]
                out_text = out_text + [str(x) for x in bw_summary]
                out_text = "\t".join(out_text)
                matrix_file.write(out_text + "\n")

            else:

                if bw_summary is None:

                    out = 0

                else:

                    out = sorted(bw_summary, reverse=True)
                    out = out[0:n_highest]
                    out = np.mean(out)

                results.append((i.name, "b", str(cpt), float(out)))

        # Add the matrix to the list (only if profile is true)
        if profile:
            matrix_file_list.append(matrix_file.name)
            matrix_file.close()

    if profile:

        return matrix_file_list

    else:
        return results


def command_log(log_file):

    logger = logging.getLogger(__name__)

    logging.basicConfig(filename=log_file,
                        level=logging.INFO,
                        format='%(asctime)s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    logger.info(" ".join(sys.argv))


def flatten_list(alist):
    """flatten a list of lists (with string...)."""
    import itertools
    return list(itertools.chain.from_iterable(
        itertools.repeat(x, 1) if isinstance(x, str) else x for x in alist))


def make_outdir_and_file(out_dir=None, alist=None, add_date=True):
    """Create output directory and a set of associated files (alist). Return a list of file handler (write mode) for the requested files."""

    if out_dir is not None:
        # Working directory
        if os.path.exists(out_dir):
            message("Aborting. Directory already exist.")
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

    return(out_fh_list)


def make_parser():
    """The program parser."""

    parser = argparse.ArgumentParser(add_help=True)

    parser.add_argument(
        '--infile', '-i',
        metavar='File',
        help='The input gtf file (default to STDIN)',
        default=sys.stdin,
        dest='in_file',
        nargs='?',
        type=argparse.FileType('r'),
        required=False)

    parser.add_argument(
        '--chromInfo', '-c',
        help='Chromosome information. A tabulated two-columns file with '
        'chromosomes as column 1 and sizes as column 2.',
        type=argparse.FileType('r'),
        dest='chrom_file',
        required=True)

    parser.add_argument(
        '--coveragefileList', '-l',
        type=str,
        default=None,
        help='Coverage file (bigWig) or a comma separated list of coverage'
        ' files.',
        dest='bw_list',
        required=True)

    parser.add_argument(
        '--outdir',
        '-o',
        default='tss_cov_results',
        dest='out_dir',
        help="The output directory where results will be stored.",
        required=False)

    parser.add_argument(
        '--upstream', '-u',
        type=int,
        default=1000,
        help="Extend the region in 5' by a given value (int). Used to define the region around the TSS searching for promoters or regions around TTS. Default 0.")

    parser.add_argument(
        '--downstream', '-d',
        type=int,
        default=1000,
        help="Extend the region in 3' by a given value (int). Used to define the region around the TSS searching for promoters or regions around TTS. Default 0.")

    parser.add_argument(
        '--binNb', '-w',
        type=int,
        dest='bin_nb',
        default=100,
        help='Split the region into w bins (see -n).',
        required=False)

    parser.add_argument(
        '--nbThreads', '-k',
        type=int,
        default=1,
        dest='nb_proc',
        help='Use this many threads to compute coverage.',
        required=False)

    parser.add_argument(
        '--labels',
        '-L',
        type=str,
        default=None,
        help="A comma separated list of label"
        " for each bigwigs (e.g short version for plotting).",
        required=False)

    parser.add_argument("--verbose",
                        "-v",
                        action='store_true',
                        help="Display processing and warning messages.",
                        required=False)

    parser.add_argument("--tmp-dir",
                        "-K",
                        type=str,
                        default=None,
                        help="Keep all temporary files in this directory.",
                        required=False)

    return parser


def cov_around_tss(
        in_file=None,
        out_dir=None,
        bw_list=None,
        upstream=1000,
        downstream=1000,
        chrom_file=None,
        tmp_file=None,
        bin_nb=100,
        nb_proc=None,
        labels=None,
        color_ramp=None,
        profile_colors=None,
        nb_palette_col=10,
        lower_limit=0,
        upper_limit=1,
        tx_order=None,
        heatmap_width=None,
        heatmap_height=None,
        show_rownames=False,
        color_overlay=None,
        tmp_dir=None,
        verbose=False):
    """
    Create a coverage matrix around TSS.
    """

    pseudo_count = 1

    # The input bigWigs
    bw_list = bw_list.split(",")
    bw_tmp_list = list()

    # Compute the list of bigWigs (including those infered from wildcards).
    for bw_file in bw_list:
        import glob
        if re.search("\*", bw_file):
            bw_tmp_list += sorted(glob.glob(bw_file))
        else:
            bw_tmp_list += [bw_file]

    bw_list = bw_tmp_list

    # Create a list of labels for the diagrams.
    # Take user input in account
    add_label_to_file_name = False
    if labels is not None:
        add_label_to_file_name = True
        if re.search(",", labels):
            labels = labels.split(",")
        else:
            labels = [labels]
        user_entered_labels = True
    else:
        labels = [os.path.splitext(os.path.basename(x))[0] for x in bw_list]
        user_entered_labels = False

    if len(labels) != len(bw_list):
        message("ERROR: The number of labels should be the same"
                " as the number of bigWig.")
        message("Labels: " + " ".join(labels))
        message("bigWig: " + " ".join(bw_list))
        sys.exit()

    bw_name_short = [os.path.splitext(os.path.basename(x))[0] for x in bw_list]

    pos = list()

    """
    # Get labels of reference bigWig ('to split')
    for bwts in bw_to_split:
        pos.append(bw_name_short.index(bwts))

    bw_to_split = [labels[p] for p in pos]
    """

    # Print file names and associated labels
    # and create a dict (bigWig to label)
    bw_to_labels = dict()
    for bw, lab, short_name in zip(bw_list, labels, bw_name_short):
        message(bw + " has label : " + lab)
        bw_to_labels[short_name] = lab

    # Prepare output  files

    if add_label_to_file_name:

        file_label = "_" + "-".join(labels)
    else:
        file_label = ""

    file_out_list = make_outdir_and_file(out_dir,
                                         ["command.log",
                                          "coverage_matrix.txt"
                                          ],
                                         add_date=False)

    log_file, coverage_matrix_file = file_out_list

    # Store the user command
    command_log(log_file.name)

    if chrom_file is None:
        message("Please provide a file with chromosome length.")
        exit(0)

    if verbose:
        message("Getting promoter regions [-" +
                str(upstream) +
                ",+" +
                str(downstream) +
                "].")

    # Where to store the region coord
    region_bed_file = make_tmp_file(prefix="tss_", suffix=".bed")
    region_bed_file_slop = make_tmp_file(prefix="tss_slop_", suffix=".bed")

    # Search tss coordinates

    region_coord = BedTool(in_file)
    tx_tss = defaultdict(lambda: 0)
    tx_chrom = defaultdict(lambda: 0)
    tx_strand = defaultdict(lambda: 0)

    for i in BedTool(in_file):

        if i[2] == 'exon':
            tx_id = re.search('transcript_id "(.*?)"', i[8]).group(1)
            if tx_id in tx_tss:
                if i.strand == '-':
                    if i.end > tx_tss[tx_id]:
                        tx_tss[tx_id] = i.end
                elif i.strand == '+':
                    if i.start < tx_tss[tx_id]:
                        tx_tss[tx_id] = i.start
            else:
                tx_strand[tx_id] = i.strand
                tx_chrom[tx_id] = i.chrom
                if i.strand == '-':
                    tx_tss[tx_id] = i.end
                elif i.strand == '+':
                    tx_tss[tx_id] = i.start

    for tx_id, tss in tx_tss.iteritems():
        if tx_strand[tx_id] == '-':
            out_list = [tx_chrom[tx_id],
                        str(int(tx_tss[tx_id]) - 1),
                        str(int(tx_tss[tx_id])),
                        tx_id,
                        "0",
                        tx_strand[tx_id]]
        elif tx_strand[tx_id] == '+':
            out_list = [tx_chrom[tx_id],
                        str(tx_tss[tx_id]),
                        str(int(tx_tss[tx_id]) + 1),
                        tx_id,
                        "0",
                        tx_strand[tx_id]]
        region_bed_file.write("\t".join(out_list) + "\n")
    region_bed_file.close()

    bed_obj = BedTool(region_bed_file.name)
    bed_obj = bed_obj.slop(s=True,
                           l=upstream,
                           r=downstream,
                           g=chrom_file.name)

    tokens = intervals(range(len(bed_obj)), nb_proc)
    bed_obj.saveas(region_bed_file_slop.name)

    # Computing coverage of features.
    # Each worker will send a file
    pool = multiprocessing.Pool(nb_proc)
    matrix_files = pool.map_async(big_wig_summary_worker,
                                  zip(tokens,
                                      repeat(bw_list),
                                      repeat(region_bed_file_slop.name),
                                      repeat(bin_nb),
                                      repeat(pseudo_count),
                                      repeat(None),
                                      repeat(True),
                                      repeat(verbose))).get(9999999)

    # If their was a problem processing the regions
    # the workers return False
    if False in matrix_files:
        message("Aborting...")
        sys.exit()

    matrix_files = flatten_list(matrix_files)

    # merging worker's files.
    #matrix_merge_file = make_tmp_file(prefix='matrix_merge_file', suffix='.txt')

    with coverage_matrix_file as outfile:
        for fname in matrix_files:
            with open(fname) as infile:
                outfile.write(infile.read())
    if verbose:
        message("Merged file (" +
                'promoter' + ") " +
                coverage_matrix_file.name +
                " ready.")

    for i in flatten_list(TMP_FILE_LIST):
        if tmp_dir is not None:
            if not os.path.exists(tmp_dir):
                msg = "Creating directory {d}."
                message(msg.format(d=tmp_dir))
                os.mkdir(tmp_dir)
            if not os.path.isdir(tmp_dir):
                msg = "{d} is not a directory."
                message(msg.format(d=tmp_dir))
                exit(0)
            if verbose:
                message("Keeping temporary file :" + i)
                base_name_i = os.path.basename(i)
                shutil.move(i, os.path.join(tmp_dir, base_name_i))
        else:
            if verbose:
                message("Deleting temporary file :" + i)
            silentremove(i)


def main():
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    cov_around_tss(**args)

if __name__ == '__main__':
    main()
