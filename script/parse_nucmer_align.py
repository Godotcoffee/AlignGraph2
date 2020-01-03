#!/bin/env python
import os
import subprocess
from subprocess import PIPE
import sys
import threading


def get_fasta_inf(fasta_path):
    name_list = list()
    length_list = list()

    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>') or line.startswith(';'):
                name_list.append(line[1:].split()[0])
                length_list.append(0)
            else:
                length_list[-1] += len(line.strip())

    return list(zip(name_list, length_list))


def parse_mummer_delta(delta_path):
    align_set = set()

    with open(delta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                sp = line[1:].split()
                align_set.add(((sp[1], int(sp[3])), (sp[0], int(sp[2]))))

    return align_set


def parse_align(line_iter):
    begin_tag = '-- BEGIN alignment '
    end_tag = '--   END alignment '

    is_begin, in_count, ref_align, query_align = False, 0, list(), list()
    ref_start, ref_end, query_start, query_end = 0, 0, 0, 0

    align = list()

    for line in line_iter:
        if not line or line[0].isspace():
            continue
        if line.startswith(begin_tag):
            is_begin, in_count, ref_align, query_align = True, 0, list(), list()

            sp = line[len(begin_tag) + 1:-2].split()
            ref_start, ref_end = int(sp[1]), int(sp[3])
            query_start, query_end = int(sp[6]), int(sp[8])
        elif line.startswith(end_tag):
            is_begin = False

            align.append((query_start, query_end, ref_start, ref_end,
                          ''.join(query_align).replace('.', '-'),
                          ''.join(ref_align).replace('.', '-')))
        elif is_begin:
            if in_count == 0:
                ref_align.append(line.split()[1])
                in_count = 1
            elif in_count == 1:
                query_align.append(line.split()[1])
                in_count = 0

    return align


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


def main(argv):
    show_align_exec = argv[0]
    delta_path = argv[1]
    out_path = argv[2]
    wrk_dir = argv[3]
    thread_num = int(argv[4])

    print_lock = threading.Lock()

    align_pair = list(parse_mummer_delta(delta_path))

    main.print_cnt = 0

    def thread_fun(start):
        tmp_path = os.path.join(wrk_dir, str(start))

        with open(tmp_path, 'w') as f:
            for i in range(start, len(align_pair), thread_num):
                (query_name, query_len), (ref_name, ref_len) = align_pair[i]
                cp = subprocess.run([show_align_exec, delta_path, ref_name, query_name], stdout=PIPE, stderr=PIPE)
                aligns = parse_align(cp.stdout.decode('ascii').split('\n'))
                with print_lock:
                    main.print_cnt += 1
                    print('\'{}\' ref to \'{}\' total={} ({}/{})'.format(query_name, ref_name, len(aligns), main.print_cnt, len(align_pair)))
                for align in aligns:
                    is_forward = align[0] <= align[1]
                    f.write('{}\t{}\t'.format(query_name, ref_name))
                    f.write('{}\t'.format('F' if is_forward else 'R'))
                    f.write('{}\t'.format('NULL'))
                    f.write('{}\t'.format(align[0] - 1 if is_forward else align[1] - 1))
                    f.write('{}\t'.format(align[1] if is_forward else align[0]))
                    f.write('{}\t'.format(query_len))
                    f.write('{}\t{}\t{}\n'.format(align[2] - 1, align[3], ref_len))
                    f.write('{}\n{}\n'.format(align[4], align[5]))

    thread_list = list()
    merge_list = list()

    for i in range(0, thread_num):
        thread_list.append(threading.Thread(target=thread_fun, args=(i, )))
        merge_list.append(os.path.join(wrk_dir, str(i)))

    for thread in thread_list:
        thread.start()

    for thread in thread_list:
        thread.join()

    merge_files(out_path, *merge_list)


if __name__ == '__main__':
    if len(sys.argv) < 6:
        print('Usage: parse_nucmer_align.py show_align_exec delat_file output_file wrk_dir thread')
    exit()
    main(sys.argv[1:])
