from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO
import sys
import os
import subprocess

def long2ref(mecat_cmd, ctg_path, ref_path, wrk_dir, thread_num, out_path):
    split_len = 10000
    seq_cnt = 0
    
    id_len = dict()
    idx_id = list()

    new_seq_list = list()

    for seq in SeqIO.parse(ctg_path, 'fasta'):
        seq_str = str(seq.seq)
        num_part = (len(seq) + split_len - 1) //split_len
        seq_size = len(seq.seq)

        for i in range(0, num_part):
            pseq = SeqRecord(Seq(
                seq_str[i * split_len:min(seq_size, (i + 1) * split_len)],
                IUPAC.protein),
                id = '{}_{}'.format(seq_cnt, i), description='')
            pseq.letter_annotations["phred_quality"] = [0] * (min(seq_size, (i + 1) * split_len) - i * split_len)
            new_seq_list.append(pseq)

        idx_id.append(seq.id)
        id_len[seq.id] = seq_size

        seq_cnt +=1
    
    tmp_read_path = os.path.join(wrk_dir, 'read.fq')

    SeqIO.write(new_seq_list, tmp_read_path, 'fastq')

    tmp_ref_path = os.path.join(wrk_dir, 'mecat.ref')
    tmp_ref2_path = os.path.join(wrk_dir, 'mecat.x.ref')

    ret = subprocess.run([mecat_cmd,
                    '-t', str(thread_num),
                    '-d', tmp_read_path,
                    '-r', ref_path,
                    '-b', str(1),
                    '-w', './wrk',
                    '-o', 'dummy',
                    '-p', tmp_ref_path],
                   cwd=wrk_dir)

    if ret.returncode == 0:
            import filter
            filter.filter_error(tmp_ref_path, tmp_ref2_path)
    
    with open(tmp_ref2_path) as in_ref, open(out_path, 'w') as out_ref:
        for line in in_ref:
            if line and line[0].isdigit():
                sp = line.split()
                o_id = new_seq_list[int(sp[0])].id
                sp2 = o_id.split('_')
                oseq_id = idx_id[int(sp2[0])]
                oseq_len = id_len[oseq_id]
                offset = int(sp2[1]) * split_len
            
                out_ref.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    oseq_id, sp[1], sp[2], sp[3], int(sp[4]) + offset, int(sp[5]) + offset, oseq_len, sp[7], sp[8], sp[9]
                ))
            else:
                out_ref.write(line)


if __name__ == "__main__":
    long2ref('/rhome/huangs/Github/AlignGraph2/mecat_plus/MECAT-master_1/Linux-amd64/bin/mecat2ref',
         sys.argv[1], sys.argv[2], sys.argv[3], 16, sys.argv[4])
