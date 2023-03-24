import os
import re
import argparse
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser('Extract the junction coordinates of gapped-reads with 4 gaps for each segment in each cell/sample from sam file. Note: the given boundaries are corresponding to the last base at the 3prime end of the first proportion of the gapped reads and the 5prime end of the last proportion of the gapped reads.')

parser.add_argument(
    'files', nargs='+', help='input sam file for each cell')
parser.add_argument(
    '-r', '--ref', required=True, help='a list of ref sequence id text file (e.g. PR8_ref_seq_id.txt)')
parser.add_argument(
    '-m', '--min_length', required=True, type=int, help='the minimum length of parts of reads mapped to 5 and 3 end separately')
parser.add_argument(
    '-sl', '--skip_length', required=True, type=int, help='the minimum length of the skip region in cigar')
parser.add_argument(
    '-o', '--output', required=True, help='the output file name')

args = parser.parse_args()

uppercase = re.compile(r'([A-Z])')


def is_header(line):
    return line.startswith('@')


def in_range(value, range_5, range_3):
    return range_5 <= value <= range_3


def extract_align_info(split_line):
    sposition1 = int(split_line[3])
    cigar = {'cha': [], 'len': []}
    last = None
    for i, s in enumerate(re.split(uppercase, split_line[5])):
        if s:
            if i % 2 == 0:
                last = int(s)
            else:
                cigar['cha'].append(s)
                cigar['len'].append(last)
    return sposition1, cigar


def dip_filter_record(split_line, sposition1, cigar, min_length, skip_length):
    pattern = ''.join(cigar['cha'])
    if pattern == 'MNMNMNMNM':
        # filter reads based on the length of parts of mapped reads
        if cigar['len'][0] >= min_length and cigar['len'][2] >= min_length and cigar['len'][4] >= min_length and cigar['len'][6] >= min_length and cigar['len'][8] >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][1] >= skip_length and cigar['len'][3] >= skip_length and cigar['len'][5] >= skip_length and cigar['len'][7] >= skip_length:
                readid = split_line[0]
                segment = split_line[2]
                boundary5 = sposition1 + cigar['len'][0] - 1
                secondm_5start = sposition1 + cigar['len'][0] + cigar['len'][1]
                secondm_3stop = sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] - 1
                thirdm_5start = sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3]
                thirdm_3stop = sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] - 1
                fourthm_5start = sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5]
                fourthm_3stop = sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] - 1
                boundary3 = sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7]
                dip_filtered_boundary = (segment, readid, boundary5, secondm_5start, secondm_3stop, thirdm_5start, thirdm_3stop, fourthm_5start, fourthm_3stop, boundary3, np.nan, np.nan, np.nan, np.nan)
                return dip_filtered_boundary
    if pattern == 'SMNMNMNMNM':
        # filter reads based on the length of parts of mapped reads
        if cigar['len'][1] >= min_length and cigar['len'][3] >= min_length and cigar['len'][5] >= min_length and cigar['len'][7] >= min_length and cigar['len'][9] >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][2] >= skip_length and cigar['len'][4] >= skip_length and cigar['len'][6] >= skip_length and cigar['len'][8] >= skip_length:
                readid = split_line[0]
                segment = split_line[2]
                boundary5 = sposition1 + cigar['len'][1] - 1
                secondm_5start = sposition1 + cigar['len'][1] + cigar['len'][2]
                secondm_3stop = sposition1 + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] - 1
                thirdm_5start = sposition1 + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4]
                thirdm_3stop = sposition1 + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] - 1
                fourthm_5start = sposition1 + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6]
                fourthm_3stop = sposition1 + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] - 1
                boundary3 = sposition1 + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] + cigar['len'][8]
                dip_filtered_boundary = (segment, readid, boundary5, secondm_5start, secondm_3stop, thirdm_5start, thirdm_3stop, fourthm_5start, fourthm_3stop, boundary3, np.nan, np.nan, np.nan, np.nan)
                return dip_filtered_boundary


def main():
    segments = set()
    with open(args.ref) as f:
        for line in f:
            split_line = line.split('\r')
            segments.add(split_line[0])
    # segment_boundary = {'CY098877.1', 'CY098876.1', 'CY098875.1', 'CY098870.1', 'CY098873.1', 'CY098872.1', 'CY098871.1', 'CY098874.1'}

    results = []
    for file_ in args.files:
        print 'Processing file: ' + file_
        sid = os.path.basename(file_).split('.')[0]
        header_only_file = True
        with open(file_) as f:
            for line in f:
                if not is_header(line):
                    header_only_file = False
                    split_line = line.split('\t')
                    segment = split_line[2]
                    if segment in segments:
                        sposition1, cigar = extract_align_info(split_line)
                        dip_filtered_boundary = dip_filter_record(split_line, sposition1, cigar, args.min_length, args.skip_length)
                        if dip_filtered_boundary is not None:
                            results.append((sid,) + dip_filtered_boundary)
        if header_only_file:
            for segment in segments:
                results.append((sid,) + (segment, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan))

    results = pd.DataFrame(results)
    results.to_csv(args.output)

if __name__ == '__main__':
    main()
