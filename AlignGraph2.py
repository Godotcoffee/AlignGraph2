import argparse
import shutil
import sys
import os
import subprocess
import time

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
    parser.add_argument('-k', required=False, type=int, default=14,
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

    args = parser.parse_args()

    start_time = time.time()

    # Input
    read_path = os.path.realpath(args.read)
    ctg_path = os.path.realpath(args.contig)
    ref_path = os.path.realpath(args.ref)
    out_dir = os.path.realpath(args.output)
    k = args.k
    ratio = args.ratio
    thread_num = args.thread

    # Get root
    root_dir = os.path.dirname(os.path.realpath(__file__))
    # PAGraph module
    kmer_count_cmd = os.path.join(root_dir, 'PAGraph', 'build', 'kmer_counter')
    pre_process_cmd = os.path.join(root_dir, 'PAGraph', 'build', 'pre_process')
    pagraph_cmd = os.path.join(root_dir, 'PAGraph', 'build', 'pagraph')
    pa_cns_cmd = os.path.join(root_dir, 'PAGraph', 'build', 'pa_cns')
    # MECAT
    uname = os.uname()
    mecat_ref_cmd = os.path.join(root_dir, 'thirdparty', 'mecat',
                                 '{}-{}'.format(uname.sysname, 'amd64' if uname.machine == 'x86_64' else uname.machine),
                                 'bin', 'mecat2ref')
    # MUMmer
    nucmer_cmd = os.path.join(root_dir, 'thirdparty', 'mummer', 'nucmer')
    delta_fileter_cmd = os.path.join(root_dir, 'thirdparty', 'mummer', 'delta-filter')

    # Test cmd
    subprocess.run([kmer_count_cmd], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    subprocess.run([pre_process_cmd], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    subprocess.run([pagraph_cmd], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    subprocess.run([pa_cns_cmd], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    subprocess.run([mecat_ref_cmd], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    subprocess.run([nucmer_cmd], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    wrk_dir = os.path.join(out_dir, 'working_dir')
    mecat_dir = os.path.join(wrk_dir, 'mecat')
    mecat_ctg_dir = os.path.join(mecat_dir, 'ctg')
    mecat_ref_dir = os.path.join(mecat_dir, 'ref')
    mummer_dir = os.path.join(wrk_dir, 'mummer')
    mummer_tmp_dir = os.path.join(mummer_dir, 'tmp')
    input_dir = os.path.join(wrk_dir, 'input')
    pagraph_dir = os.path.join(wrk_dir, 'pagraph')
    cns_dir = os.path.join(wrk_dir, 'cns')
    cns_in_dir = os.path.join(cns_dir, 'input')
    cns_out_dir = os.path.join(cns_dir, 'output')
    cns_wrk_dir = os.path.join(cns_dir, 'wrk')
    cns_wrk_mecat_dir = os.path.join(cns_wrk_dir, 'mecat')
    cns_wrk_split_dir = os.path.join(cns_wrk_dir, 'split')
    cns_wrk_ref_dir = os.path.join(cns_wrk_dir, 'ref')

    os.makedirs(wrk_dir, exist_ok=True)
    os.makedirs(mecat_ctg_dir, exist_ok=True)
    os.makedirs(mecat_ref_dir, exist_ok=True)
    os.makedirs(mummer_dir, exist_ok=True)
    os.makedirs(mummer_tmp_dir, exist_ok=True)
    os.makedirs(input_dir, exist_ok=True)
    os.makedirs(pagraph_dir, exist_ok=True)
    os.makedirs(cns_dir, exist_ok=True)
    os.makedirs(cns_in_dir, exist_ok=True)
    os.makedirs(cns_out_dir, exist_ok=True)
    os.makedirs(cns_wrk_dir, exist_ok=True)
    os.makedirs(cns_wrk_mecat_dir, exist_ok=True)
    os.makedirs(cns_wrk_split_dir, exist_ok=True)
    os.makedirs(cns_wrk_ref_dir, exist_ok=True)

    solid_kmer_path = os.path.join(wrk_dir, 'solid_kmer_set.bin')
    read_to_ctg_path = os.path.join(mecat_ctg_dir, 'read_to_contig.txt')
    read_to_ref_path = os.path.join(mecat_ref_dir, 'read_to_ref.txt')
    ctg_to_ref_path = os.path.join(mummer_dir, 'ctg_to_ref.txt')

    # K-mer count
    print('K-Mer counting...')
    subprocess.run([kmer_count_cmd,
                    '-t', str(thread_num),
                    '-i', read_path,
                    '-o', solid_kmer_path,
                    '-k', str(k),
                    '-m', str(ratio)])
    print('Done')

    # Read to Contig
    print('Read to Contig...')
    subprocess.run([mecat_ref_cmd,
                    '-t', str(thread_num),
                    '-d', read_path,
                    '-r', ctg_path,
                    '-b', str(1),
                    '-w', './wrk',
                    '-o', read_to_ctg_path],
                   cwd=mecat_ctg_dir)
    print('Done')

    # Read to Ref
    print('Read to Ref...')
    subprocess.run([mecat_ref_cmd,
                    '-t', str(thread_num),
                    '-d', read_path,
                    '-r', ref_path,
                    '-b', str(1),
                    '-w', './wrk',
                    '-o', read_to_ref_path],
                   cwd=mecat_ref_dir)
    print('Done')
    
    # Contig to Ref
    print('Contig to Ref...')

    # nucmer
    subprocess.run([nucmer_cmd,
                    '-t', str(thread_num),
                    ref_path,
                    ctg_path],
                   cwd=mummer_dir)

    # delta-filter
    out_filter_path = os.path.join(mummer_dir, 'out.filter.delta')
    with open(out_filter_path, 'w') as filter_f:
        subprocess.run([delta_fileter_cmd,
                        '-q',
                        '-r',
                        os.path.join(mummer_dir, 'out.delta')],
                       cwd=mummer_dir,
                       stdout=filter_f)
    
    # convert
    show_align_cmd = os.path.join(root_dir, 'thirdparty', 'mummer', 'show-aligns')
    import script.parse_nucmer_align
    script.parse_nucmer_align.main([show_align_cmd, out_filter_path, ctg_to_ref_path, mummer_tmp_dir, str(thread_num)])

    print('Done')

    # pre_process
    print('Pre process...')

    subprocess.run([pre_process_cmd, 
                    '-r', read_path,
                    '-c', ctg_path,
                    '-x', read_to_ctg_path,
                    '-y', read_to_ref_path,
                    '-z', ctg_to_ref_path,
                    '-o', input_dir])

    print('Done')
    

    # pagraph
    print('PAGraph...')

    subprocess.run([pagraph_cmd,
                    '-t', str(thread_num),
                    '-k', solid_kmer_path,
                    '-r', read_path,
                    '-c', ctg_path,
                    '-R', ref_path,
                    '-p', input_dir,
                    '-a', ctg_to_ref_path,
                    '-o', pagraph_dir])

    print('Done')

    import script.cns_helper as cns_helper

    # extract
    print('Extract and split...')

    import script.extract
    script.extract.action(ctg_path, pagraph_dir, os.path.join(pagraph_dir, 'contig.txt'), cns_in_dir)
    mapping = cns_helper.split_fasta(os.path.join(cns_in_dir, 'add.fasta'), cns_wrk_split_dir)

    print('Done')

    # align
    print('Align and split...')

    subprocess.run([mecat_ref_cmd,
                    '-t', str(thread_num),
                    '-d', read_path,
                    '-r', os.path.join(cns_in_dir, 'all.fasta'),
                    '-b', str(1),
                    '-w', './wrk',
                    '-o', os.path.join(cns_wrk_mecat_dir, 'merge.ref')],
                   cwd=cns_wrk_mecat_dir)

    cns_helper.split_ref(os.path.join(cns_wrk_mecat_dir, 'merge.ref'), cns_wrk_ref_dir, mapping)

    print('Done')

    # Correct
    print('Correct...')
    part_len = 5000
    top_k = 3000

    fasta_cor_path_list = list()

    for ctg, idx in mapping.items():
        print('\tcorrecting {}'.format(ctg))

        subprocess.run([pa_cns_cmd, '-i', os.path.join(cns_wrk_split_dir, idx),
                        '-a', os.path.join(cns_wrk_ref_dir, '{}.ref'.format(idx)),
                        '-o', os.path.join(cns_out_dir, '{}.new'.format(idx)),
                        '-l', str(part_len),
                        '-t', str(thread_num),
                        '-k', str(top_k)])

        fasta_cor_path_list.append(os.path.join(cns_out_dir, '{}.new'.format(idx)))

    f_out = os.path.join(cns_out_dir, 'cor.fasta')
    cns_helper.merge_fasta(f_out, *fasta_cor_path_list)
    print('Out file: {}'.format(f_out))

    print('Done')

    final_out = os.path.join(out_dir, 'final.fasta')
    cns_helper.merge_fasta(final_out, os.path.join(cns_in_dir, 'include.fasta'), os.path.join(cns_out_dir, 'cor.fasta'))
    shutil.copyfile(os.path.join(cns_in_dir, 'include.fasta'), os.path.join(out_dir, 'remainder.fasta'))
    shutil.copyfile(os.path.join(cns_in_dir, 'exclude.fasta'), os.path.join(out_dir, 'exclude.fasta'))
    shutil.copyfile(os.path.join(cns_out_dir, 'cor.fasta'), os.path.join(out_dir, 'add.fasta'))

    print('Final output: {}'.format(final_out))

    print('Time used: {:.3f} seconds'.format(time.time() - start_time))
