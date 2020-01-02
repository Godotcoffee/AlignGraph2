
import os
import shutil
import subprocess
import time

if __name__ == '__main__':
    total_start_time = time.time()
    root_dir = os.path.dirname(os.path.realpath(__file__))
    log_path = os.path.join(root_dir, 'build.log')
    with open(log_path, 'w') as log_f:
        pass    # Clear file
    print('Log file is located in {}'.format(log_path))

    with open(log_path, 'a') as log_f:
        # Build PAGraph
        print('Building PAGraph...')
        start_time = time.time()

        pagraph_dir = os.path.join(root_dir, 'PAGraph')
        pagraph_build = os.path.join(pagraph_dir, 'build')
        if os.path.isdir(pagraph_build):
            shutil.rmtree(pagraph_build)
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
                             cwd=os.path.join(pagraph_dir, 'build'))
        print('#', file=log_f, flush=True)
        print('# End of CMake', file=log_f, flush=True)

        if ret.returncode != 0:
            print('Failed to build PAGraph, please check {}'.format(log_f.name))
            exit(1)

        print('#', file=log_f, flush=True)
        print('# Begin of make', file=log_f, flush=True)
        print('#', file=log_f, flush=True)
        ret = subprocess.run(['make', '-j4'],
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

        # Build MECAT+
        print('Building MECAT+...')
        start_time = time.time()

        print('#', file=log_f, flush=True)
        print('# Begin of Building MECAT+', file=log_f, flush=True)
        print('#', file=log_f, flush=True)

        # TODO

        print('#', file=log_f, flush=True)
        print('# End of Building MECAT+', file=log_f, flush=True)
        print('#', file=log_f, flush=True)

        print('Done. Took {:.3f} seconds'.format(time.time() - start_time))

    # Build MECAT+

