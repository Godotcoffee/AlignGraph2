import argparse
import sys

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='AlignGraph2', formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--version', action='version', version='%(prog)s 1.0beta')
    parser.add_argument('-r', '--read', required=True, type=str, default=argparse.SUPPRESS,
                        help='read path')
    parser.add_argument('-c', '--contig', required=True, type=str, default=argparse.SUPPRESS,
                        help='contig path')
    parser.add_argument('-R', '--ref', required=True, type=str, default=argparse.SUPPRESS,
                        help='reference path')
    parser.add_argument('-o', '--output', required=True, type=str, default=argparse.SUPPRESS,
                        help='output directory')
    parser.add_argument('-k', required=False, type=float, default=14,
                        help='size of k-mer')
    parser.add_argument('--ratio', required=False, type=float, default=0.2,
                        help='threshold of solid k-mer set')
    parser.add_argument('-t', '--thread', required=False, type=int, default=16,
                        help='thread number')
    parser.add_argument('--clean', required=False, action='store_true',
                        help='clean file after running')

    if len(sys.argv) <= 1:
        parser.print_help(file=sys.stderr)
        parser.exit()
