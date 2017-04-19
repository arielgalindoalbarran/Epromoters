"""
Compute transcript coverage over bed feature with one or several bigWig.
------------------------------------------------------------------------

Use bx-python as interface to the kent utilities.

Author: D. Puthier
"""
from __future__ import print_function
from __future__ import division
import re
import sys
import argparse
import multiprocessing
from tempfile import NamedTemporaryFile
from itertools import repeat
from pybedtools import BedTool
import pandas as pd
from bx.bbi.bigwig_file import BigWigFile
import numpy as np

#-------------------------------------------------------------------------------
# A set of basic functions
#-------------------------------------------------------------------------------

def close_properly(*args):
    """Close a set of file if they are not None."""
    for afile in args:
        if afile is not None:
            if afile != sys.stdout:
                afile.close()

def write_properly(string, afile):
    """Write a string to a file. If file is None, write string to stdout."""
    if afile is not None:
        afile.write(string + "\n")
    else:
        sys.stdout.write(string + "\n")

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

#-------------------------------------------------------------------------------
# The main parser
#-------------------------------------------------------------------------------

def make_parser():

    parser = argparse.ArgumentParser(add_help=True)

    input_gtf = parser.add_argument_group('Input gtf file')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('more optional arguments')


    # required but not required...
    input_gtf.add_argument(
        '--inputfile', '-i',
        metavar='File',
        help='The input bed file (default to STDIN)',
        default=sys.stdin,
        nargs='?',
        type=argparse.FileType('r'),
        required=False)

    required.add_argument(
        '--bw-list', '-l',
        type=str,
        default=None,
        help='Coverage file (bigWig) or a comma separated list of coverage'
        ' files.',
        required=True)

    required.add_argument(
        '--score',
        '-s',
        type=str,
        help="""The score to be computed. Each bigwig file is identified
        by b0, b1, b2,... 
        The following operations are accepted: +, -, *, /, **.
        Each element in the formula corresponds to the average signal value in a given region.
        E.g : (b0+b1)/b2*0.5.""",
        required=True)

    optional.add_argument(
        '--bin-nb', '-w',
        type=int,
        default=1,
        help='Split the region into w bins (see -n).',
        required=False)

    optional.add_argument(
        '--nb-proc', '-k',
        type=int,
        default=1,
        help='Use this many threads to compute coverage.',
        required=False)

    optional.add_argument(
        '--pseudo-count',
        '-p',
        type=int,
        default=1,
        help='A pseudo-count to add in case count is equal to 0.')


    optional.add_argument(
        '--n-highest', '-n',
        type=int,
        default=None,
        help='If the region is splitted into several bins (see -w) '
        'compute the region score for a given bigWig (b0, b1,...)'
        ' based on the average of the n bin with highest coverage values.    ',
        required=False)

    optional.add_argument(
        '--out-file', '-o',
        help='The output file name.',
        type=argparse.FileType('w'),
        default=sys.stdout,
        nargs='?',
        required=False)

    optional.add_argument(
        '--verbose', '-v',
        help='Ask for verbosity.',
        action="store_true")

    return parser

#-------------------------------------------------------------------------------
# The function  performing coverage computation on bwig
#-------------------------------------------------------------------------------

# Avoid some warnings that come from bx python
# Although there are no problem at all with the regions
# for which coverage is to be computed.
# Seems erratic as calling twice time with the same
# regions give or not a runtime warning.
np.seterr(divide='ignore', invalid='ignore')


def big_wig_summary_worker((span,
                            bw_list,
                            region_bed_file_name,
                            bin_nb,
                            pseudo_count,
                            n_highest,
                            verbose)):
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
    * verbose: run in verbose mode.
    """

    if n_highest is None:
        n_highest = bin_nb

    results = list()

    for big_wig, cpt in zip(bw_list, range(len(bw_list))):

        bigwig = BigWigFile(open(big_wig, "r"))

        if verbose:
            mesg = "Computing coverage for file: {bw} [{proc}], {chunk} chunks to process.\n"
            sys.stderr.write(mesg.format(bw=big_wig,
                                         proc=str(multiprocessing.current_process()),
                                         chunk=str(span[1] - span[0])))

        # Load the regions for which the coverage is to be processed.

        tx_bed = BedTool(region_bed_file_name)

        # There is no  current chrom at the moment
        chr_cur = None

        # The fraction of bed file
        from_here = span[0]
        to_here = span[1]

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
                mesg = "{name} : the number of nucleotides in feature ({nb})" + \
                       "is below the window length ({wl}). Skipping."

                sys.stderr.write(mesg.format(name=i.name,
                                             nb=str(i.end - i.start),
                                             wl=str(bin_nb)))
                continue

            bw_summary = bigwig.summarize(i.chrom,
                                          i.start,
                                          i.end,
                                          bin_nb)
            if bw_summary is None:

                if verbose:
                    mesg = "{name} : no coverage info found (check if chrom," +\
                           "{chr} is declared in bigWig file"


                # Should be replaced with NA
                bw_summary = np.repeat(pseudo_count, bin_nb)

            else:
                bw_summary = bw_summary.sum_data / ((i.end - i.start) / bin_nb)
                bw_summary = bw_summary + pseudo_count

            # Replace nan by pseudo_count
            # Should be replaced with NA
            nan_pos = np.isnan(bw_summary)

            bw_summary[nan_pos] = pseudo_count


            if bw_summary is None:

                out = 0

            else:

                out = sorted(bw_summary, reverse=True)
                out = out[0:n_highest]
                out = np.mean(out)

            results.append((i.name, "b", str(cpt), float(out)))

    return results



#-------------------------------------------------------------------------------
# The main function
#-------------------------------------------------------------------------------

def bw_coverage(
        inputfile=None,
        out_file=None,
        bw_list=None,
        pseudo_count=1,
        score=None,
        bin_nb=1,
        n_highest=None,
        nb_proc=1,
        verbose=True):
    """
    Compute transcript coverage with one or several bigWig.
    -------------------------------------------------------
    Uses bx-python as interface to kent utilities.
    """

    # Check if the score is well written

    if not re.search(r"^[b\d\/\*\+\-\(\)\.]+$", score):
        sys.stderr.write(
            "Score should contain the following characters: "
            "b0, b1 (...) and operators +, ., -, *, /, **, (, ).")
        sys.exit(0)

    # Check if the score to compute fits with
    # The number of input bigWigs
    bw_list = bw_list.split(",")
    bwig_in_score = re.finditer(r"b\d+", score)
    bwig_expected_in_score = ["b" + str(x) for x in range(len(bw_list))]

    for i in bwig_in_score:

        if i.group(0) not in bwig_expected_in_score:
            sys.stderr.write("The indicated column (" +
                             i.group(0) +
                             ") was not found.")
            sys.exit(0)

    # Check the number of windows
    if n_highest is None:
        n_highest = bin_nb

    if verbose:
        sys.stderr.write("Number of bins: " + str(bin_nb) + "\n")
        sys.stderr.write("N highest values: " + str(n_highest) + "\n")

    if n_highest > bin_nb:
        sys.stderr.write("The number of window used for computing the score"
                         " (-n) can not be greater than the number of"
                         " windows (-w)")
        sys.exit()


    # Check input file is in bed6 format

    region_bo = BedTool(inputfile.name)


    if region_bo.field_count() != 6:
        sys.stderr.write("Bed file should should be in Bed6 format. Use '.' if strand is undefined.\n")
        sys.exit()

    tokens = intervals(range(len(BedTool(inputfile.name))), nb_proc)

    pool = multiprocessing.Pool(nb_proc)
    coverage_list = pool.map_async(big_wig_summary_worker,
                                   zip(tokens,
                                       repeat(bw_list),
                                       repeat(inputfile.name),
                                       repeat(bin_nb),
                                       repeat(pseudo_count),
                                       repeat(n_highest),
                                       repeat(verbose))).get(9999999)


    if False in coverage_list:
        sys.stderr.write("Aborting...")
        sys.exit()

    # Unlist the list of list
    coverage_list = [item for sublist in coverage_list for item in sublist]


    # Prepare a data.frame to collect the results
    dataframe = pd.DataFrame(columns=None)

    if verbose:
        sys.stderr.write("Retrieving results.\n")

    for i in coverage_list:
        dataframe.ix[i[0], i[1] + str(i[2])] = float(i[3])

    if verbose:
        sys.stderr.write("Computing score.\n")


    dataframe = dataframe.eval(score)

    dataframe.to_csv(out_file,
                     sep="\t",
                     header=False)

    close_properly(inputfile, out_file)


#-------------------------------------------------------------------------------
# Call  to main
#-------------------------------------------------------------------------------


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    bw_coverage(**args)

if __name__ == '__main__':
    main()
