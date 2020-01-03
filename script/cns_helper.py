import os
from Bio import SeqIO


def merge_files(output_path, *input_paths):
    buf_size = 65536

    with open(output_path, 'wb') as out_file:
        for input_path in input_paths:
            with open(input_path, 'rb') as in_file:
                while True:
                    data = in_file.read(buf_size)
                    if data:
                        out_file.write(data)
                    else:
                        break


def merge_fasta(output_fasta, *input_fasta):
    merge_files(output_fasta, *input_fasta)


def split_fasta(input_fasta, output_dir):
    idx = 0
    mapping = dict()

    for seq_record in SeqIO.parse(input_fasta, 'fasta'):
        SeqIO.write([seq_record], os.path.join(output_dir, str(idx)), 'fasta')
        mapping[seq_record.id] = str(idx)
        idx += 1

    return mapping


def split_ref(ref_path, output_dir, mapping):

    ctg_to_file = dict()

    try:
        for ctg, idx in mapping.items():
            ctg_to_file[ctg] = open(os.path.join(output_dir, '{}.ref'.format(idx)), 'w')

        with open(ref_path) as f:
            while True:
                line_one = f.readline()
                if not line_one:
                    break

                line_two = f.readline()
                line_three = f.readline()

                fasta_name = line_one.split()[1]

                out_f = ctg_to_file.get(fasta_name)
                if out_f:
                    out_f.write(line_one)
                    out_f.write(line_two)
                    out_f.write(line_three)

    finally:
        for f in ctg_to_file.values():
            f.close()

