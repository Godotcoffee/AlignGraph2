import sys
import os

def filter_error(ref_path, out_path):
    with open(ref_path, 'r', errors='ignore') as in_file, open(out_path, 'w') as o_file:
        line_cnt = 0

        line1, line2, line3 = '', '', ''

        for line in in_file:
            if line and line[0].isdigit():
                line_cnt = 0
                line1 = line
                line_cnt += 1
            elif line_cnt == 1:
                line2 = line
                line_cnt += 1
            elif line_cnt == 2:
                line3 = line
                line_cnt += 1

            if line_cnt == 3 and len(line2) == len(line3):
                o_file.write('{}{}{}'.format(line1, line2, line3))
                line_cnt = 0


