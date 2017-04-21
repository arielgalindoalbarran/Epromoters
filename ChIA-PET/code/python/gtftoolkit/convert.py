#!/usr/bin/env python
"""
 Description: Convert a GTF to various format.
"""
from __future__ import print_function
import sys
import argparse
import re


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

def is_comment(string):
    """Check wether a string is a comment (starts with #...).

    :param string: The string to be tested.
    :Example:

    >>> from gtftk.utils import is_comment
    >>> assert is_comment("#bla")

    """

    return string.startswith("#")

def make_parser():
    """The program parser."""
    parser = argparse.ArgumentParser(add_help=True, 
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser_grp = parser.add_argument_group('Arguments')

    parser_grp.add_argument('-i', '--inputfile',
                            help="Path to the GTF file. Default to STDIN",
                            default=sys.stdin,
                            metavar="GTF",
                            type=argparse.FileType("r"))

    parser_grp.add_argument('-o', '--outputfile',
                            help="output file",
                            default=sys.stdout,
                            metavar="GTF",
                            type=argparse.FileType('w'))

    parser_grp.add_argument('-n', '--names',
                            help="The key(s) that should be used as name.",
                            default="gene_id,transcript_id",
                            metavar="NAME",
                            type=str)

    parser_grp.add_argument('-s', '--separator',
                            help="The separator to be used for separating name elements (see -n).",
                            default="|",
                            metavar="SEP",
                            type=str)

    parser_grp.add_argument('-f', '--format',
                            help='Currently one of bed3, bed/bed6',
                            type=str,
                            choices=('bed', 'bed3', 'bed6'),
                            default='bed6',
                            metavar="FORMAT",
                            required=False)

    parser_grp.add_argument('-t', '--feature-type',
                            help='Which feature to get: transcript, gene or exon.',
                            type=str,
                            choices=('transcript', 'gene', 'exon'),
                            default='transcript',
                            metavar="FEATURE_TYPE",
                            required=False)   

    parser_grp.add_argument('-c', '--concat',
                            help='Controls whether the NAME column should'
                            ' contain: transcript ID (t), gene ID (g) or both'
                            '(b) separated.',
                            type=str,
                            choices=('t', 'g', 'b'),
                            default='t',
                            metavar="CONCAT",
                            required=False) 

    return parser


def convert(inputfile=None,
            outputfile=None,
            format="bed",
            names="gene_id,transcript_id",
            separator="|",
            feature_type="transcript",
            concat="t",
            verbosity=0):
    """Convert a GTF to various format."""


    if feature_type in ['transcript', 'gene']:

        feature_starts = dict()
        feature_ends = dict()
        feature_chrom = dict()
        feature_strand = dict()

        for line in inputfile:
        
            columns = line.split("\t")
            
            if columns[2] == "exon":

                start_current = int(columns[3])
                end_current = int(columns[4])
            
                if feature_type == "transcript":

                    tx_id = re.search('.*transcript_id "(.*?)"', 
                                            columns[8]).group(1)

                    gene_id = re.search('.*gene_id "(.*?)"', 
                                             columns[8]).group(1)

                    feature_name =   tx_id + separator + gene_id             
                                            
                else:
                    feature_name = re.search('.*gene_id "(.*?)"', 
                                             columns[8]).group(1)

                if feature_name not in feature_starts:
                    feature_strand[feature_name] = columns[6]
                    feature_chrom[feature_name] = columns[0]
                    feature_starts[feature_name] = start_current
                    feature_ends[feature_name] = end_current
            
                else:

                    if start_current < feature_starts[feature_name]:
            
                        feature_starts[feature_name] = start_current
            
                    if end_current > feature_ends[feature_name]:
            
                        feature_ends[feature_name] = end_current
            

        for i in feature_starts:

            if format == "bed3":
                out = [feature_chrom[i],
                       str(feature_starts[i] -1),
                       str(feature_ends[i])]

            elif format in ["bed","bed6"]:
                name = i
                if feature_type == "transcript":
                    if concat == "t":
                        name = i.split(separator)[0]
                    elif concat == "b":
                        name = i
                    elif concat == "g":
                        name = i.split(separator)[1]
                    
                out = [feature_chrom[i],
                       str(feature_starts[i] -1),
                       str(feature_ends[i]),
                       name,
                       "0",
                       feature_strand[i]]
            outputfile.write("\t".join(out) + "\n")
               
    else:
        for line in inputfile:
        
            columns = line.split("\t")
            
            if columns[2] == "exon":

                start_current = int(columns[3])
                end_current = int(columns[4])

                tx_id = re.search('.*transcript_id "(.*?)"', 
                                  columns[8]).group(1)
                gene_id = re.search('.*gene_id "(.*?)"', 
                                    columns[8]).group(1)
                
                if concat == "t":
                    name = tx_id
                elif concat == "b":
                    name = tx_id + separator + gene_id
                elif concat == "g":
                    name = gene_id
                    
                out = [columns[0],
                       str(start_current - 1),
                       str(end_current),
                       name,
                       "0",
                       columns[6]]
                
                outputfile.write("\t".join(out) + "\n")
                               
    close_properly(outputfile, inputfile)
    

def main():
    """The main program."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    convert(**args)

if __name__ == '__main__':
    main()

