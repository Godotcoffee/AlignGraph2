from Bio import SeqIO

import sys
import threading

def rv(ch):
    c = ch.upper()
    if c == 'A':
        return 'T'
    elif c == 'C':
        return 'G'
    elif c == 'G':
        return 'T'
    elif c == 'T':
        return 'A'
    else:
        return 'N'


def paf2aln(ctg_path, ref_path, paf_path, out_path, thread_num=16):
    ctg_dict = SeqIO.to_dict(SeqIO.parse(ctg_path, "fasta"))
    ref_dict = SeqIO.to_dict(SeqIO.parse(ref_path, "fasta"))

    test_cnt = 0

    all_lines = open(paf_path).readlines()

    print_lock = threading.Lock()

    def thread_fun(start):
        global test_cnt
        with open(out_path + '_{}'.format(start), 'w') as aln:
            for i in range(start, len(all_lines), thread_num):
                line = all_lines[i]
                sp = line.split('\t')
                ctg_name = sp[0]
                ref_name = sp[5]
                forward = 'F' if sp[4] == '+' else 'R'
                ctg_len = sp[1]
                ctg_start = sp[2]
                ctg_end = sp[3]
                ref_len = sp[6]
                ref_start = sp[7]
                ref_end = sp[8]
                cigar_str = sp[13]

                aln.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(ctg_name, ref_name, forward, 'NULL', ctg_start, ctg_end, ctg_len, ref_start, ref_end, ref_len))

                cb = int(ctg_start) if forward == 'F' else int(ctg_end) - 1
                cs = 1 if forward == 'F' else -1

                rb = int(ref_start)

                stb1, stb2 = '', ''
                nb = ''

                query = ctg_dict[ctg_name]
                ref = ref_dict[ref_name]
                #query = 'A' * 100000

                for ch in cigar_str[5:-1]:
                    if ch.isdigit():
                        nb += ch
                    else:
                        gap = int(nb)
                        nb = ''
                        for j in range(gap):
                            if ch == 'M':
                                stb1 += query[cb] if forward == 'F' else rv(query[cb])
                                cb += cs
                                stb2 += ref[rb]
                                rb += 1
                            elif ch == 'D':
                                stb1 += '-'
                                stb2 += ref[rb]
                                rb += 1
                            elif ch == 'I':
                                stb1 += query[cb] if forward == 'F' else rv(query[cb])
                                cb += cs
                                stb2 += '-'
                
                aln.write('{}\n{}\n'.format(stb1, stb2))

                with print_lock:
                    test_cnt += 1
                    if test_cnt % 100 == 0:
                        print(test_cnt)

    thread_list = list()

    for i in range(0, thread_num):
        thread_list.append(threading.Thread(target=thread_fun, args=(i, )))
    
    for thread in thread_list:
        thread.start()

    for thread in thread_list:
        thread.join()
