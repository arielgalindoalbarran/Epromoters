#!/usr/bin/env python
"""
 Convert a bed file to a gtf (but with lots of empty fields...). May be helpful sometimes...
"""

from __future__ import print_function
import sys
import argparse
from pybedtools import BedTool
from tempfile import NamedTemporaryFile


def close_properly(*args):
    """Close a set of file if they are not None."""
    for afile in args:
        if afile is not None:
            if afile != sys.stdout:
                afile.close()

def chomp(string):
    """Delete \r\n from end of string."""
    string = string.rstrip('\r\n')
    return string

def write_properly(string, afile):
    """Write a string to a file. If file is None, write string to stdout."""
    if afile is not None:
        afile.write(string + "\n")
    else:
        sys.stdout.write(string + "\n")


def make_parser():
    """The program parser."""
    parser = argparse.ArgumentParser(add_help=True)

    parser_grp = parser.add_argument_group('Arguments')

    parser_grp.add_argument('-i', '--inputfile',
                            help="Path to BED file you would like to "
                            "behave as if it was a GTF.",
                            default=sys.stdin,
                            metavar="BED",
                            type=argparse.FileType("r"))

    parser_grp.add_argument('-o', '--outputfile',
                            help="output file",
                            default=sys.stdout,
                            metavar="GTF",
                            type=argparse.FileType('w'))

    parser_grp.add_argument('-t', '--ft-type',
                            help="The type of features you are trying to "
                            "immitate...",
                            default='transcript',
                            type=str,
                            required=False)

    parser_grp.add_argument('-s', '--source',
                            help="The source of annotation.",
                            default='Unknown',
                            type=str,
                            required=False)

    parser_grp.add_argument('-v', '--verbose',
                            help="Verbose if true.",
                            action='store_true')
    
    return parser


def bed_to_gtf(
        inputfile=None,
        outputfile=None,
        ft_type="transcript",
        source="Unknown",
        verbose=False):
    """
 Convert a bed file to a gtf. 
    """

    if verbose:
        sys.stderr.write("Converting the bed file into GTF file\n.")

    if inputfile.name == '<stdin>':
        tmp_file = NamedTemporaryFile()
        for i in inputfile:
            write_properly(chomp(str(i)), tmp_file)

        tmp_file.close()
        inputfile.close()

        bed_obj = BedTool(tmp_file.name)
    else:
        bed_obj = BedTool(inputfile.name)

    n = 1
    for i in bed_obj:

        if i.strand == "":
            i.strand = "."
        if i.name == "":
            i.name = str("feature_" + str(n))
        if i.score == "":
            i.score = "0"

        if ft_type == "exon":
            key_value = "gene_id \"" + i.name + "\"; " + \
                        "transcript_id \"" + i.name + "\"; " + \
                        "exon_id \"" + i.name + "\";"
        elif ft_type == "gene":
            key_value = "gene_id \"" + i.name + "\";"
        else:
            key_value = "gene_id \"" + i.name + "\"; " + \
                        "transcript_id \"" + i.name + "\";"

        list_out = [i.chrom,
                    source,
                    ft_type,
                    str(i.start + 1),
                    str(i.end),
                    str(i.score),
                    i.strand,
                    ".",
                    key_value]

        write_properly("\t".join(list_out), outputfile)

        n += 1

    close_properly(outputfile)


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    bed_to_gtf(**args)

if __name__ == '__main__':
    main()
