from Bio import SeqIO
import sys
import os

def split_pre_process(ctg_path, ref_path, aln_path, dir_path, out_path):
    cfg_path = os.path.join(dir_path, 'config.txt')

    ctg = SeqIO.to_dict(SeqIO.parse(ctg_path, 'fasta'))
    ref = SeqIO.to_dict(SeqIO.parse(ref_path, 'fasta'))

    cnt = 0

    line_cnt = 0
    ref_id = ''
    file_list = list()
    ctg_list = list()

    aln_map = dict()
    out_flist = list()

    with open(cfg_path, 'r') as cf:
        for line in cf:
            if line_cnt == 0:
                ref_id = line[:-1]
                line_cnt += 1
            elif line_cnt < 4:
                file_list.append(line[:-1])
                line_cnt += 1
            elif line != '\n':
                ctg_list.append(line[:-1])
                line_cnt += 1
            else:
                to_dir = os.path.join(out_path, str(cnt))
                os.makedirs(to_dir, exist_ok=True)

                aln_map[ref_id] = cnt

                SeqIO.write(ref[ref_id], os.path.join(to_dir, 'ref.fasta'), 'fasta')
                ctg_seq = list()
                for i in range(0, len(ctg_list), 2):
                    ctg_seq.append(ctg[ctg_list[i]])
                SeqIO.write(ctg_seq, os.path.join(to_dir, 'ctg.fasta'), 'fasta')
                for f in file_list:
                    if os.path.exists(os.path.join(to_dir, f)):
                        os.remove(os.path.join(to_dir, f))
                    os.symlink(os.path.abspath(os.path.join(dir_path, f)), os.path.join(to_dir, f))
                with open(os.path.join(to_dir, 'config.txt'), 'w') as co:
                    co.write('{}\n'.format(ref_id))
                    co.write('{}\n'.format(file_list[0]))
                    co.write('{}\n'.format(file_list[1]))
                    co.write('{}\n'.format(file_list[2]))
                    for c in ctg_list:
                        co.write('{}\n'.format(c))
                out_flist.append(open(os.path.join(to_dir, 'aln'), 'w'))
                cnt += 1
                line_cnt = 0
                file_list.clear()
                ctg_list.clear()

    with open(aln_path) as ai:
        line_cnt = 0
        line1 = ''
        line2 = ''
        line3 = ''
        for line in ai:
            if line_cnt % 3 == 0:
                line1 = line
            elif line_cnt % 3 == 1:
                line2 = line
            else:
                line3 = line
            
            if line_cnt % 3 == 2:
                sp = line1.split('\t')
                if sp[1] in aln_map:
                    out_flist[aln_map[sp[1]]].write('{}{}{}'.format(line1, line2, line3))

            line_cnt += 1
    
    for f in out_flist:
        f.close()


def merge_out(sp_dir, mg_dir):
    with open(os.path.join(mg_dir, 'contig.txt'), 'w') as co, open(os.path.join(mg_dir, 'coninfo'), 'w') as ci:
        for d in os.listdir(sp_dir):
            if os.path.isfile(os.path.join(sp_dir, d, 'DONE')):
                with open(os.path.join(sp_dir, d, 'contig.txt')) as one_c:
                    co.writelines(one_c.readlines())
                for ctg in (f for f in os.listdir(os.path.join(sp_dir, d)) if os.path.splitext(f)[1] == '.fasta'):
                    seq = next(SeqIO.parse(os.path.join(sp_dir, d, ctg), 'fasta'))
                    seq.id = d + '_' +  seq.id
                    SeqIO.write(seq, os.path.join(mg_dir, '{}_{}'.format(d, ctg)), 'fasta')
                for con in (f for f in os.listdir(os.path.join(sp_dir, d)) if os.path.splitext(f)[1] == '.con'):
                    with open(os.path.join(sp_dir, d, con)) as cii:
                        ci.writelines(cii.readlines())
                        ci.write('\n')
