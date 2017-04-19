"""
    The ``col_from_tab`` module
    ======================

    Select columns from a tabulated file based on their names.

"""

import re
import sys
import argparse

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
    required = parser.add_argument_group('required arguments')
    input_gtf = parser.add_argument_group('Input gtf file')
    optional = parser.add_argument_group('more optional arguments')


    # required but not required...
    required.add_argument(
        '--infile', '-i',
        metavar='File',
        help='The input tabulated file (default to STDIN)',
        default=sys.stdin,
        type=argparse.FileType('r'),
        required=False)

    required.add_argument(
        '--columns', '-c',
        type=str,
        default=None,
        help='Comma separated list of column names.',
        required=True)

    optional.add_argument(
        '--outfile', '-o',
        help='The output file name.',
        type=argparse.FileType('w'),
        default=sys.stdout,
        nargs='?',
        required=False)

    optional.add_argument(
        '--sep', '-s',
        default="\t",
        help="The input and output columns separator.",
        required=False)

    optional.add_argument(
        '--invert', '-V',
        action="store_true",
        help="Invert match. All columns except those provided.",
        required=False)

    optional.add_argument(
        '--keep', '-k',
        action="store_true",
        help="Print the header.",
        required=False)
     
    optional.add_argument('-v', '--verbose',
                           help="Verbose if true.",
                           action='store_true')
       

    return parser

def col_from_tab(infile=None,
                 outfile=None,
                 columns=None,
                 invert=False,
                 sep=None,
                 keep=False,
                 verbose=False):
    
    if re.search(",", columns):
        columns = columns.split(",")
    else:
        columns = [columns]



    for p,line in enumerate(infile):
    
        line = chomp(line)
        line = line.split(sep)

        if p == 0:
            
            if not invert:

                pos_list = list()
                
                for i in range(len(columns)):
                    
                    pos = line.index(columns[i]) if columns[i]  in line else -1
    
                    if pos > -1:
                        pos_list.append(pos)
                    else:
                        sys.stderr.write("Column " + columns[i] + " not found")
                        sys.exit()
            else:

                pos_list = range(len(line))

                for i in range(len(columns)):

                    pos = line.index(columns[i]) if columns[i]  in line else -1

                    if pos > -1:
                        pos_list.remove(pos)
                    else:
                        sys.stderr.write("Column " + columns[i] + " not found")
                        sys.exit()

            if keep:
                out = sep.join([line[k] for k in pos_list])
                write_properly(out, outfile)
        else:
            out = sep.join([line[k] for k in pos_list])
            write_properly(out, outfile)
            
    close_properly(infile, outfile)
            
def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    col_from_tab(**args)

if __name__ == '__main__':
    main()            
                    
                
                
                
    
    