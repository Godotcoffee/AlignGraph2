
import os
import shutil
import subprocess
import time
import stat
import tarfile
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Installation of AlignGraph2', formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--version', action='version', version='%(prog)s 1.0beta')
    parser.add_argument('-t', '--thread', metavar='[int]', required=False, type=int, default=4,
                        help='thread number')
    #parser.add_argument('--clean', required=False, action='store_true',
    #                    help='clean file after running')

    args = parser.parse_args()

    total_start_time = time.time()
    root_dir = os.path.dirname(os.path.realpath(__file__))
    log_path = os.path.join(root_dir, 'build.log')
    build_thread = args.thread

    uname = os.uname()
    mecat_ref_cmd = os.path.join(root_dir, 'thirdparty', 'mecat',
                                '{}-{}'.format(uname.sysname, 'amd64' if uname.machine == 'x86_64' else uname.machine),
                                'bin', 'mecat2ref')
    mecat_ref2_cmd = os.path.join(root_dir, 'mecat_plus', 'MECAT-master_1',
                                 '{}-{}'.format(uname.sysname, 'amd64' if uname.machine == 'x86_64' else uname.machine),
                                 'bin', 'mecat2ref')

    try:
        subprocess.run(['cmake', '--version'],
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.STDOUT)
    except FileNotFoundError:
        print('Error: Couldn\'t find CMake')
        exit(1)

    print('Build with {} threads'.format(build_thread))
    print(r'Number of threads can be modified by adding \'-t {thread_num}\'')

    with open(log_path, 'w') as log_f:
        pass    # Clear file
    print('Log file is located in {}'.format(log_path))

    with open(log_path, 'a') as log_f:
        # Build PAGraph
        print('Building PAGraph...')
        start_time = time.time()

        pagraph_dir = os.path.join(root_dir, 'PAGraph')
        pagraph_build = os.path.join(pagraph_dir, 'build')

        os.makedirs(pagraph_build, exist_ok=True)

        print('#', file=log_f, flush=True)
        print('# Begin of Building PAGraph', file=log_f, flush=True)
        print('#', file=log_f, flush=True)

        print('#', file=log_f, flush=True)
        print('# Being of CMake', file=log_f, flush=True)
        print('#', file=log_f, flush=True)

        ret = subprocess.run(['cmake', '..'],
                             stdout=log_f,
                             stderr=subprocess.STDOUT,
                             cwd=pagraph_build)

        print('#', file=log_f, flush=True)
        print('# End of CMake', file=log_f, flush=True)

        if ret.returncode != 0:
            print('Failed to build PAGraph, please check {}'.format(log_f.name))
            exit(1)

        print('#', file=log_f, flush=True)
        print('# Begin of make', file=log_f, flush=True)
        print('#', file=log_f, flush=True)

        ret = subprocess.run(['make', '-j{}'.format(build_thread)],
                             stdout=log_f,
                             stderr=subprocess.STDOUT,
                             cwd=os.path.join(pagraph_dir, 'build'))

        print('# End of make', file=log_f, flush=True)
        print('#', file=log_f, flush=True)

        if ret.returncode != 0:
            print('Failed to build pagraph, please check {}'.format(log_f.name))
            exit(1)

        print('#', file=log_f, flush=True)
        print('# End of Building PAGraph', file=log_f, flush=True)
        print('#', file=log_f, flush=True)

        print('Done. Took {:.3f} seconds'.format(time.time() - start_time))

        # Build MECAT
        print('Building MECAT...')
        start_time = time.time()

        mecat_dir = os.path.join(root_dir, 'thirdparty', 'mecat')
        mecat_p_dir = os.path.join(root_dir, 'mecat_plus', 'MECAT-master_1')

        print('#', file=log_f, flush=True)
        print('# Begin of Building MECAT+', file=log_f, flush=True)
        print('#', file=log_f, flush=True)

        print('#', file=log_f, flush=True)
        print('# Begin of make', file=log_f, flush=True)
        print('#', file=log_f, flush=True)

        ret = subprocess.run(['make', '-j{}'.format(build_thread)],
                             stdout=log_f,
                             stderr=subprocess.STDOUT,
                             cwd=mecat_dir)

        if not os.path.exists(mecat_ref_cmd):
            print('Failed to build mecat2ref, please check {}'.format(log_f.name))
            exit(1)
        
        try:
            ret = subprocess.run(['make', '-j{}'.format(build_thread)],
                                stdout=log_f,
                                stderr=subprocess.STDOUT,
                                cwd=mecat_p_dir)
        except:
            pass

        if not os.path.exists(mecat_ref2_cmd):
            print('Warning: failed to build mecat+')

        print('# End of make', file=log_f, flush=True)
        print('#', file=log_f, flush=True)

        print('# End of Building MECAT', file=log_f, flush=True)
        print('#', file=log_f, flush=True)

        print('Done. Took {:.3f} seconds'.format(time.time() - start_time))

        # Build MUMmer
        print('Building MUMmer...')
        start_time = time.time()

        mummer_dir = os.path.join(root_dir, 'thirdparty', 'mummer')
        mummer_tar = os.path.join(root_dir, 'thirdparty', 'mummer-4.0.0beta2.tar.gz')

        if os.path.exists(mummer_dir) == False:
            t = tarfile.open(mummer_tar)
            t.extractall(path=os.path.join(root_dir, 'thirdparty'))
            os.rename(os.path.join(root_dir, 'thirdparty', 'mummer-4.0.0beta2'), mummer_dir)

        print('#', file=log_f, flush=True)
        print('# Begin of Building MUMmer', file=log_f, flush=True)
        print('#', file=log_f, flush=True)

        print('#', file=log_f, flush=True)
        print('# Begin of configure', file=log_f, flush=True)
        print('#', file=log_f, flush=True)

        ret = subprocess.run(['./configure'],
                             stdout=log_f,
                             stderr=subprocess.STDOUT,
                             cwd=mummer_dir)

        if ret.returncode != 0:
            print('Failed to build mummer, please check {}'.format(log_f.name))
            exit(1)

        print('# End of configure', file=log_f, flush=True)
        print('#', file=log_f, flush=True)

        print('#', file=log_f, flush=True)
        print('# Begin of make', file=log_f, flush=True)
        print('#', file=log_f, flush=True)

        ret = subprocess.run(['make', '-j{}'.format(build_thread)],
                             stdout=log_f,
                             stderr=subprocess.STDOUT,
                             cwd=mummer_dir)

        if ret.returncode != 0:
            print('Failed to build mummer, please check {}'.format(log_f.name))
            exit(1)

        print('# End of make', file=log_f, flush=True)
        print('#', file=log_f, flush=True)

        print('# End of Building MUMmer', file=log_f, flush=True)
        print('#', file=log_f, flush=True)

        print('Done. Took {:.3f} seconds'.format(time.time() - start_time))

    print('All Done. Total Took {:.3f} seconds'.format(time.time() - total_start_time))

