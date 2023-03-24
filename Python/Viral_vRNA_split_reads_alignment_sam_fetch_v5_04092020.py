import re
import argparse
import collections


parser = argparse.ArgumentParser('Extract all the alignment fallen in the coding seqeunce regions for each viral segment. Note any alignment with any number of Ns will be considered')

parser.add_argument(
    '--ref_CDS_position', required=True, help='reference sequence coding region positions (e.g., a tab-delimited file containing information such as AF389115.1 28 2307')

parser.add_argument(
    '--input_sam_file', required=True, help='input sam file')

parser.add_argument(
    '--output_sam_file', required=True, help='output sam file')

args = parser.parse_args()

UPPERCASE = re.compile(r'([A-Z])')

PATTERNS = {
    'MNM',
    'SMNM',
    'MNMS',
    'MNMDM',
    'SMNMS',
    'MDMNM',
    'MNMIM',
    'MIMNM',
    'SMDMNM',
    'MNMDMS',
    'SMNMDM',
    'SMNMIM',
    'SMIMNM',
    'MDMNMS',
    'MIMNMS',
    'MNMIMS',
    'MDMIMNM',
    'SMDMNMS',
    'SMNMIMS',
    'SMIMNMS',
    'MDMDMNM',
    'MNMIMIM',
    'SMNMDMS',
    'MNMDMDM',
    'MNMDMIM',
    'MIMDMNM',
    'MIMIMNM',
    'MNMIMDM',
    'MDMNMIMS',
    'SMIMNMDM',
    'MIMNMIM',
    'MDMNMIM',
    'MIMNMDM',
    'MDMNMDM',
    'MDMDMNMS',
    'MNMNM',
    'MNMNMS',
    'SMNMNM',
    'MDMNMNM',
    'MNMNMNM',
    'SMNMNMS',
    'MNMDMNM',
    'MNMDMNMS',
    'MNMNMIM',
    'MNMNMDMS',
    'MNMNMDM',
    'MNMNMIMS',
    'SMNMNMDM',
    'MIMNMNMS',
    'SMNMDMNM',
    'SMDMNMNM',
    'MNMIMNM',
    'SMNMDMDM',
    'SMDMDMNM',
    'MIMNMNM',
    'MDMNMDMS',
    'MDMNMDMNM',
    'SMDMNMDM',
    'MDMNMNMS',
    'SMNMIMNM',
    'SMIMNMNM',
    'MDMNMDMNMS',
    'MNMDMDMS',
    'SMNMNMNM',
    'MNMNMNMS',
    'MNMNMNMDM',
    'MDMNMNMNM',
    'MNMNMDMNM',
    'MNMDMNMNM',
    'MNMNMNMIM',
    'MNMNMNMNM',
    'SMNMNMNMS',
    'MNMIMDMS',
    'SMNMIMDM',
    'SMDMNMNMS',
    'MDMDMNMNM',
    'MNMNMDMDM',
    'MNMDMNMDM',
    'SMNMNMIM',
    'SMNMNMDMS',
    'SMNMDMNMS',
    'MDMDMDMNM',
    'SMIMDMNM',
    'MDMDMNMDM',
    'MNMDMNMNMS',
    'SMDMNMIM',
    'SMDMDMDMNM',
    'MNMDMDMDM',
    'SMDMIMNM',
    'SMNMNMNMNM',
    'MDMNMDMDM',
    'MDMNMNMDM',
    'SMIMNMNMS',
    'SMNMDMDMS',
    'SMDMDMNMNM',
    'SMDMNMDMS',
    'SMDMDMNMS',
    'MNMIMNMS',
    'MDMDMIMNM',
    'SMNMNMNMDM',
    'MNMNMIMNM',
    'MNMDMDMNM',
    'MNMNMDMDMS',
    'MNMDMDMNMS',
    'MNMDMIMNM',
    'MDMIMNMS',
    'MNMNMDMDMNM',
    'MNMIMNMNM',
    'MNMIMNMDM',
    'MNMIMDMNM',
    'SMNMDMDMNM',
    'MNMDMNMDM',
    'MIMNMDMDM'
}


def is_header(line):
    return line.startswith('@')


def is_in_range(pattern, cigar, split_line, cds):
    segment = split_line[2]
    info = cds[segment]
    firstbase = int(split_line[3])

    if pattern == 'MNM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] - 1

    if pattern == 'SMNM':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] - 1

    if pattern == 'MNMS':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] - 1

    if pattern == 'MNMDM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] - 1

    if pattern == 'SMNMS':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] - 1

    if pattern == 'MDMNM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] - 1

    if pattern == 'MNMIM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][4] - 1

    if pattern == 'MIMNM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] - 1

    if pattern == 'SMDMNM':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] - 1

    if pattern == 'MNMDMS':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] - 1

    if pattern == 'SMNMDM':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] - 1

    if pattern == 'SMNMIM':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][5] - 1

    if pattern == 'SMIMNM':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] - 1

    if pattern == 'MDMNMS':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] - 1

    if pattern == 'MIMNMS':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] - 1

    if pattern == 'MNMIMS':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][4] - 1

    if pattern == 'MDMIMNM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] - 1

    if pattern == 'SMDMNMS':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] - 1

    if pattern == 'SMNMIMS':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][5] - 1

    if pattern == 'SMIMNMS':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] - 1

    if pattern == 'MDMDMNM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] - 1

    if pattern == 'MNMIMIM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][4] + cigar['len'][6] - 1

    if pattern == 'SMNMDMS':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] - 1

    if pattern == 'MNMDMDM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] - 1

    if pattern == 'MNMDMIM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][6] - 1

    if pattern == 'MIMDMNM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] - 1

    if pattern == 'MIMIMNM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][2] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] - 1

    if pattern == 'MNMIMDM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] - 1

    if pattern == 'MDMNMIMS':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][6] - 1

    if pattern == 'SMIMNMDM':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] - 1

    if pattern == 'MIMNMIM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][6] - 1

    if pattern == 'MDMNMIM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][6] - 1

    if pattern == 'MIMNMDM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] - 1

    if pattern == 'MDMNMDM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] - 1

    if pattern == 'MDMDMNMS':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] - 1

    if pattern == 'MNMNM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] - 1

    if pattern == 'MNMNMS':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] - 1

    if pattern == 'SMNMNM':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] - 1

    if pattern == 'MDMNMNM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] - 1

    if pattern == 'MNMNMNM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] - 1

    if pattern == 'SMNMNMS':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] - 1

    if pattern == 'MNMDMNM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] - 1

    if pattern == 'MNMDMNMS':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] - 1

    if pattern == 'MNMNMIM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][6] - 1

    if pattern == 'MNMNMDMS':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] - 1

    if pattern == 'MNMNMDM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] - 1

    if pattern == 'MNMNMIMS':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][6] - 1

    if pattern == 'SMNMNMDM':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] - 1

    if pattern == 'MIMNMNMS':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] - 1

    if pattern == 'SMNMDMNM':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] - 1

    if pattern == 'SMDMNMNM':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] - 1

    if pattern == 'MNMIMNM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] - 1

    if pattern == 'SMNMDMDM':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] - 1

    if pattern == 'SMDMDMNM':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] - 1

    if pattern == 'MIMNMNM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] - 1

    if pattern == 'MDMNMDMS':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] - 1

    if pattern == 'MDMNMDMNM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] + cigar['len'][8] - 1

    if pattern == 'SMDMNMDM':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] - 1

    if pattern == 'MDMNMNMS':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] - 1

    if pattern == 'SMNMIMNM':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] - 1

    if pattern == 'SMIMNMNM':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] - 1

    if pattern == 'MDMNMDMNMS':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] + cigar['len'][8] - 1

    if pattern == 'MNMDMDMS':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] - 1

    if pattern == 'SMNMNMNM':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] - 1

    if pattern == 'MNMNMNMS':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] - 1

    if pattern == 'MNMNMNMDM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] + cigar['len'][8] - 1

    if pattern == 'MDMNMNMNM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] + cigar['len'][8] - 1

    if pattern == 'MNMNMDMNM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] + cigar['len'][8] - 1

    if pattern == 'MNMDMNMNM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] + cigar['len'][8] - 1

    if pattern == 'MNMNMNMIM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][8] - 1

    if pattern == 'MNMNMNMNM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] + cigar['len'][8] - 1

    if pattern == 'SMNMNMNMS':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] - 1

    if pattern == 'MNMIMDMS':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] - 1

    if pattern == 'SMNMIMDM':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] - 1

    if pattern == 'SMDMNMNMS':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] - 1

    if pattern == 'MDMDMNMNM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] + cigar['len'][8] - 1

    if pattern == 'MNMNMDMDM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] + cigar['len'][8] - 1

    if pattern == 'MNMDMNMDM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] + cigar['len'][8] - 1

    if pattern == 'SMNMNMIM':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][7] - 1

    if pattern == 'SMNMNMDMS':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] - 1

    if pattern == 'SMNMDMNMS':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] - 1

    if pattern == 'MDMDMDMNM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] + cigar['len'][8] - 1

    if pattern == 'SMIMDMNM':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] - 1

    if pattern == 'MDMDMNMDM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] + cigar['len'][8] - 1

    if pattern == 'MNMDMNMNMS':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] + cigar['len'][8] - 1

    if pattern == 'SMDMNMIM':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][7] - 1

    if pattern == 'SMDMDMDMNM':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] + cigar['len'][8] + cigar['len'][9] - 1

    if pattern == 'MNMDMDMDM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] + cigar['len'][8] - 1

    if pattern == 'SMDMIMNM':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] - 1

    if pattern == 'SMNMNMNMNM':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] + cigar['len'][8] + cigar['len'][9] - 1

    if pattern == 'MDMNMDMDM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] + cigar['len'][8] - 1

    if pattern == 'MDMNMNMDM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] + cigar['len'][8] - 1

    if pattern == 'SMIMNMNMS':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] - 1

    if pattern == 'SMNMDMDMS':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] - 1

    if pattern == 'SMDMDMNMNM':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] + cigar['len'][8] + cigar['len'][9] - 1

    if pattern == 'SMDMNMDMS':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] - 1

    if pattern == 'SMDMDMNMS':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] - 1

    if pattern == 'MNMIMNMS':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] - 1

    if pattern == 'MDMDMIMNM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][6] + cigar['len'][7] + cigar['len'][8] - 1

    if pattern == 'SMNMNMNMDM':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] + cigar['len'][8] + cigar['len'][9] - 1

    if pattern == 'MNMNMIMNM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][6] + cigar['len'][7] + cigar['len'][8] - 1

    if pattern == 'MNMDMDMNM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] + cigar['len'][8] - 1

    if pattern == 'MNMNMDMDMS':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] + cigar['len'][8] - 1

    if pattern == 'MNMDMDMNMS':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] + cigar['len'][8] - 1

    if pattern == 'MNMDMIMNM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][6] + cigar['len'][7] + cigar['len'][8] - 1

    if pattern == 'MDMIMNMS':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] - 1

    if pattern == 'MNMNMDMDMNM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] + cigar['len'][8] + cigar['len'][9] + cigar['len'][10] - 1

    if pattern == 'MNMIMNMNM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] + cigar['len'][8] - 1

    if pattern == 'MNMIMNMDM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] + cigar['len'][8] - 1

    if pattern == 'MNMIMDMNM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] + cigar['len'][8] - 1

    if pattern == 'SMNMDMDMNM':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] + cigar['len'][8] + cigar['len'][9] - 1

    if pattern == 'MNMDMNMDM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] + cigar['len'][8] - 1

    if pattern == 'MIMNMDMDM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] + cigar['len'][8] - 1

    return (info['min'] <= firstbase) and (lastbase <= info['max'])


def extract_align_info(split_line):
    sposition1 = int(split_line[3])
    cigar = {'cha': [], 'len': []}
    last = None
    for i, s in enumerate(re.split(UPPERCASE, split_line[5])):
        if s:
            if i % 2 == 0:
                last = int(s)
            else:
                cigar['cha'].append(s)
                cigar['len'].append(last)
    return sposition1, cigar


def main():
    cds = {}
    with open(args.ref_CDS_position) as f:
        for line in f:
            split_line = line.split('\t')
            cds[split_line[0]] = {
                'min': int(split_line[1]),
                'max': int(split_line[2])
            }

    headers = []
    with open(args.input_sam_file) as f:
        segs = collections.defaultdict(list)
        for line in f:
            if not is_header(line):
                split_line = line.split('\t')
                segs[split_line[2]].append(split_line)
            else:
                headers.append(line)

    segs = collections.OrderedDict(sorted(segs.items(), key=lambda t: t[0]))
    filtered_segs = collections.OrderedDict()
    for seg, split_lines in segs.iteritems():
        filtered_split_lines = []
        for split_line in split_lines:
            sposition1, cigar = extract_align_info(split_line)
            pattern = ''.join(cigar['cha'])
            if pattern in PATTERNS:
                if is_in_range(pattern, cigar, split_line, cds):
                    filtered_split_lines.append(split_line)
        filtered_segs[seg] = filtered_split_lines

    with open(args.output_sam_file, 'w') as f:
        f.writelines(headers)
        for seg, filtered_split_lines in filtered_segs.iteritems():
            f.writelines(map(lambda x: '\t'.join(x), filtered_split_lines))


if __name__ == '__main__':
    main()
