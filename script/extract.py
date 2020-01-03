import os
import sys
import argparse

from Bio import SeqIO


def action(seq_path, assembly_dir, config_path, out_dir):
    exclude_set = set()

    with open(config_path) as f:
        for item in f.read().splitlines():
            exclude_set.add(item)

    include_rec = list()
    exclude_rec = list()

    for seq_record in SeqIO.parse(seq_path, 'fasta'):
        if seq_record.id in exclude_set:
            exclude_rec.append(seq_record)
        else:
            include_rec.append(seq_record)

    add_rec = list()

    for fasta_path in (p for p in os.listdir(assembly_dir)
                       if os.path.isfile(os.path.join(assembly_dir, p)) and os.path.splitext(p)[1] == '.fasta'):
        for seq_record in SeqIO.parse(os.path.join(assembly_dir, fasta_path), 'fasta'):
            add_rec.append(seq_record)

    SeqIO.write(include_rec + add_rec, os.path.join(out_dir, 'all.fasta'), 'fasta')
    SeqIO.write(exclude_rec, os.path.join(out_dir, 'exclude.fasta'), 'fasta')
    SeqIO.write(add_rec, os.path.join(out_dir, 'add.fasta'), 'fasta')
    SeqIO.write(include_rec, os.path.join(out_dir, 'include.fasta'), 'fasta')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Extract', formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--version', action='version', version='%(prog)s 1.10beta')
    parser.add_argument('-s', '--sequence', required=True, type=str, default=argparse.SUPPRESS,
                        help='Sequence path')
    parser.add_argument('-c', '--config', required=True, type=str, default=argparse.SUPPRESS,
                        help='Configuration path')
    parser.add_argument('-a', '--assembly', required=True, type=str, default=argparse.SUPPRESS,
                        help='Assembly directory')
    parser.add_argument('-o', '--output', required=True, type=str, default=argparse.SUPPRESS,
                        help='Output directory')

    if len(sys.argv) <= 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()

    seq_path = os.path.abspath(args.sequence)
    assembly_dir = os.path.abspath(args.assembly)
    config_path = os.path.abspath(args.config)
    out_dir = os.path.abspath(args.output)

    action(seq_path, assembly_dir, config_path, out_dir)
