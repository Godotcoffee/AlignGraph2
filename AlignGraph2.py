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
    parser.add_argument('-r', '--read', metavar='[fastq]', required=True, type=str, default=argparse.SUPPRESS,
                        help='read path')
    parser.add_argument('-c', '--contig', metavar='[fasta]', required=True, type=str, default=argparse.SUPPRESS,
                        help='contig path')
    parser.add_argument('-g', '--genome', metavar='[fasta]', required=True, type=str, default=argparse.SUPPRESS,
                        help='reference path')
    parser.add_argument('-o', '--output', metavar='[dir]', required=True, type=str, default=argparse.SUPPRESS,
                        help='output directory')
    parser.add_argument('-m', required=False, action='store_true', default=False,
                        help='customized alignment algorithm mecat2ref+')
    parser.add_argument('-b', metavar='[int]', required=False, type=int, default=200,
                        help='size of similar genome blocks for mecat2ref+ [50-1000]')
    parser.add_argument('--alpha', metavar='[real]', required=False, type=float, default=0.5,
                        help='lower bound of k-mer scoring function for mecat2ref+ [0-1]')
    parser.add_argument('--beta', metavar='[real]', required=False, type=float, default=2.0,
                        help='upper bound of k-mer scoring function for mecat2ref+ [1-infinity]')
    parser.add_argument('--delta', metavar='[real]', required=False, type=float, default=0.9,
                        help='threshold for alignment scoring [0-1]')
    parser.add_argument('-k', metavar='[int]', required=False, type=int, default=14,
                        help='size of k-mer [4-15]')
    #parser.add_argument('--ratio', required=False, type=float, default=0.2,
    #                    help='threshold of solid k-mer set')
    parser.add_argument('--epsilon', metavar='[int]', required=False, type=int, default=10,
                        help='distance to join two vertices in A-Bruijn graph [5-100]')
    parser.add_argument('-l', metavar='[int]', required=False, type=int, default=50,
                        help='minimum path length for graph traversal [0-infinity]')
    parser.add_argument('-a', metavar='[int]', required=False, type=int, default=10000,
                        help='size of long read blocks for consensus [100-100000]')
    parser.add_argument('-t', metavar='[int]', required=False, type=int, default=16,
                        help='thread number')
    #parser.add_argument('--clean', required=False, action='store_true',
    #                    help='clean file after running')

    if len(sys.argv) <= 1:
        parser.print_help(file=sys.stderr)
        parser.exit()

    args = parser.parse_args()

    start_time = time.time()

    aligned_mode = 'MECAT' if args.m else 'NUCMER'
    #print(aligned_mode)

    # Input
    read_path = os.path.realpath(args.read)
    ctg_path = os.path.realpath(args.contig)
    ref_path = os.path.realpath(args.genome)
    out_dir = os.path.realpath(args.output)
    k = args.k
    block1 = args.b
    lower_score = args.alpha
    uppser_score = args.beta
    filter_score = args.delta
    err_dist = args.epsilon
    min_len = args.l
    block2 = args.a
    #ratio = args.ratio
    thread_num = args.t

    # Check ranges of arguments
    if not 4 <= k <= 15:
        print('Size of k-mer must be [4-15]', file=sys.stderr)
        exit(1)
    if not 50 <= block1 <= 1000:
        print('Size of similar genome blocks must be [50-1000]', file=sys.stderr)
        exit(1)
    if not 0.0 <= lower_score <= 1.0:
        print('Lower bound of k-mer scoring must be [0-1]', file=sys.stderr)
        exit(1)
    if not 1.0 <= uppser_score:
        print('Upper bound of k-mer scoring must be >= 1', file=sys.stderr)
        exit(1)
    if not 0.0 <= filter_score <= 1.0:
        print('threshold for alignment scoring must be [0-1]', file=sys.stderr)
        exit(1)
    if not 5 <= err_dist <= 100:
        print('Distance to join two vertices must be [5-100]', file=sys.stderr)
        exit(1)
    if not 0 <= min_len:
        print('Minimum path length must not be negative', file=sys.stderr)
        exit(1)
    if not 100 <= block2 <= 100000:
        print('Size of long read blocks must be [100-100000]', file=sys.stderr)
        exit(1)
    if not 0 <= thread_num:
        print('Thread number must not be negative', file=sys.stderr)
        exit(1)

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
    mecat_ref2_cmd = os.path.join(root_dir, 'mecat_plus', 'MECAT-master_1',
                                 '{}-{}'.format(uname.sysname, 'amd64' if uname.machine == 'x86_64' else uname.machine),
                                 'bin', 'mecat2ref')
    # MUMmer
    nucmer_cmd = os.path.join(root_dir, 'thirdparty', 'mummer', 'nucmer')
    delta_fileter_cmd = os.path.join(root_dir, 'thirdparty', 'mummer', 'delta-filter')

    # k8-linux
    k8_cmd = os.path.join(root_dir, 'thirdparty', 'k8-linux')

    # paftools
    paf_tools = os.path.join(root_dir, 'thirdparty', 'minimap2', 'paftools.js')

    # Test cmd
    subprocess.run([kmer_count_cmd], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    subprocess.run([pre_process_cmd], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    subprocess.run([pagraph_cmd], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    subprocess.run([pa_cns_cmd], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    subprocess.run([mecat_ref_cmd], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    subprocess.run([nucmer_cmd], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    subprocess.run([k8_cmd], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    wrk_dir = os.path.join(out_dir, 'working_dir')
    mecat_dir = os.path.join(wrk_dir, 'mecat')
    mecat_ctg_dir = os.path.join(mecat_dir, 'ctg')
    mecat_ref_dir = os.path.join(mecat_dir, 'ref')
    mummer_dir = os.path.join(wrk_dir, 'mummer')
    mummer_tmp_dir = os.path.join(mummer_dir, 'tmp')
    input_dir = os.path.join(wrk_dir, 'input', 'p')
    sp_input_dir = os.path.join(wrk_dir, 'input', 's')
    pagraph_dir = os.path.join(wrk_dir, 'pagraph')
    pagraph_m_dir = os.path.join(wrk_dir, 'pagraph2')
    cns_dir = os.path.join(wrk_dir, 'cns')
    cns_in_dir = os.path.join(cns_dir, 'input')
    cns_out_dir = os.path.join(cns_dir, 'output')
    cns_wrk_dir = os.path.join(cns_dir, 'wrk')
    cns_wrk_mecat_dir = os.path.join(cns_wrk_dir, 'mecat')
    cns_wrk_split_dir = os.path.join(cns_wrk_dir, 'split')
    cns_wrk_ref_dir = os.path.join(cns_wrk_dir, 'ref')
    out_paf_path = os.path.join(mummer_dir, 'out.paf')

    os.makedirs(wrk_dir, exist_ok=True)
    os.makedirs(mecat_ctg_dir, exist_ok=True)
    os.makedirs(mecat_ref_dir, exist_ok=True)
    os.makedirs(mummer_dir, exist_ok=True)
    os.makedirs(mummer_tmp_dir, exist_ok=True)
    try:
        if os.path.exists(sp_input_dir):
            shutil.rmtree(sp_input_dir)
        if os.path.exists(pagraph_dir):
            shutil.rmtree(pagraph_dir)
        if os.path.exists(pagraph_m_dir):
            shutil.rmtree(pagraph_m_dir)
        if os.path.exists(cns_dir):
            shutil.rmtree(cns_dir)
    except:
        pass
    os.makedirs(input_dir, exist_ok=True)
    os.makedirs(sp_input_dir, exist_ok=True)
    os.makedirs(pagraph_dir, exist_ok=True)
    os.makedirs(pagraph_m_dir, exist_ok=True)
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
    read_to_ref2_path = os.path.join(mecat_ref_dir, 'read_to_ref2.txt')
    ctg_to_ref_path = os.path.join(mummer_dir, 'ctg_to_ref.txt')
    
    # K-mer count
    print('K-Mer counting...')
    subprocess.run([kmer_count_cmd,
                    '-t', str(thread_num),
                    '-i', read_path,
                    '-o', solid_kmer_path,
                    '-k', str(k)
                    #'-m', str(ratio)
                    ])
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
    okok = False
    
    if aligned_mode == 'MECAT':
        try:
            #if os.path.exists(read_to_ref_path):
            #    os.remove(read_to_ref_path)
            #if os.path.exists(read_to_ref2_path):
            #    os.remove(read_to_ref2_path)
            #print(mecat_ref2_cmd)
            ret = subprocess.run([mecat_ref2_cmd,
                        '-t', str(thread_num),
                        '-d', read_path,
                        '-r', ref_path,
                        '-b', str(1),
                        '-w', './wrk',
                        '-o', 'dummy',
                        '-p', read_to_ref2_path,
                        '-l', str(lower_score),
                        '-u', str(uppser_score),
                        '-z', str(block2),
                        '-y', str(filter_score)],
                    cwd=mecat_ref_dir)
            #print(ret.returncode)
                    
            if ret.returncode == 0:
                import script.filter
                script.filter.filter_error(read_to_ref2_path, read_to_ref_path)
                okok = True
        except:
            pass

    if not okok:
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
    
    #if aligned_mode != 'MECAT' or not os.path.exists(mecat_ref2_cmd):
    okok2 = False

   
    try:
        import script.long2ref
        script.long2ref.long2ref(mecat_ref2_cmd, mecat_ref_cmd, ctg_path, ref_path, mummer_dir, thread_num, ctg_to_ref_path)
        okok2 = True
    except:
        pass

    if okok2 == False:
        # nucmer
        subprocess.run([nucmer_cmd,
                        '-t', str(thread_num),
                        #'-l', '100',
                        #'-c', '1000',
                        '-g', '100',
                        '--banded',
                        ref_path,
                        ctg_path],
                    cwd=mummer_dir)

        # delta-filter
        #out_filter_path = os.path.join(mummer_dir, 'out.filter.delta')
        #with open(out_filter_path, 'w') as filter_f:
        #    subprocess.run([delta_fileter_cmd,
        #                    '-q',
        #                    '-r',
        #                    os.path.join(mummer_dir, 'out.delta')],
        #                   cwd=mummer_dir,
        #                   stdout=filter_f)
        
        # convert
        #show_align_cmd = os.path.join(root_dir, 'thirdparty', 'mummer', 'show-aligns')
        #import script.parse_nucmer_align
        #script.parse_nucmer_align.main([show_align_cmd, out_filter_path, ctg_to_ref_path, mummer_tmp_dir, str(thread_num)])
        
        # convert
        with open(out_paf_path, 'w') as paf_f:
            subprocess.run([k8_cmd,
                            paf_tools,
                            'delta2paf',
                            os.path.join(mummer_dir, 'out.delta')],
                        cwd=mummer_dir,
                        stdout=paf_f)
        
        import script.paf2aln

        script.paf2aln.paf2aln(ctg_path, ref_path, out_paf_path, ctg_to_ref_path, thread_num)
    #else:
    #    import script.long2ref
    #    script.long2ref.long2ref(mecat_ref2_cmd, ctg_path, ref_path, mummer_dir, thread_num, ctg_to_ref_path)

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
    
    # split
    import script.split_helper
    script.split_helper.split_pre_process(ctg_path, ref_path, ctg_to_ref_path, input_dir, sp_input_dir)

    # pagraph
    print('PAGraph...')

    for d in os.listdir(sp_input_dir):
        in_dir = os.path.join(sp_input_dir, d)
        tmp_out_dir = os.path.join(pagraph_dir, d)
        os.makedirs(tmp_out_dir, exist_ok=True)

        #if os.path.isfile(os.path.join(out_dir, 'DONE')):
        #    print('Ignore {}'.format(d))
        #    continue

        p_ctg_path = os.path.join(in_dir, 'ctg.fasta')
        p_ref_path = os.path.join(in_dir, 'ref.fasta')
        p_aln_path = os.path.join(in_dir, 'aln')
        
        subprocess.run([
            pagraph_cmd,
            '-t', str(thread_num),
            '-r', 'dummy',
            '-k', solid_kmer_path,
            '-c', p_ctg_path,
            '-R', p_ref_path,
            '-p', in_dir,
            '-a', p_aln_path,
            '-o', tmp_out_dir,
            '-r', str(min_len),
            '--epsilon', str(err_dist)
        ])

        with open(os.path.join(tmp_out_dir, 'DONE'), 'w'):
            pass

    #subprocess.run([pagraph_cmd,
    #                '-t', str(thread_num),
    #                '-k', solid_kmer_path,
    #                '-r', read_path,
    #                '-c', ctg_path,
    #                '-R', ref_path,
    #                '-p', input_dir,
    #                '-a', ctg_to_ref_path,
    #                '-o', pagraph_dir])
    
    # merge
    script.split_helper.merge_out(pagraph_dir, pagraph_m_dir)

    print('Done')

    import script.cns_helper as cns_helper
    
    # extract
    print('Extract and split...')

    import script.extract
    script.extract.action(ctg_path, pagraph_m_dir, os.path.join(pagraph_m_dir, 'contig.txt'), cns_in_dir)
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
    #part_len = 10000
    part_len = block2
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


