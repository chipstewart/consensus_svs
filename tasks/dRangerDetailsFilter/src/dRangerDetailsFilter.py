__author__ = 'stewart'

import numpy as np
import pandas as pd
import sys,subprocess,os,argparse
if not (sys.version_info[0] == 2  and sys.version_info[1] in [7]):
    raise "Must use Python 2.7.x"

def parseOptions():
    description = '''
    Filter dRanger somatic.details.txt file to output file
    '''

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', '--pair_id', metavar='pair_id', type=str, help='pair_id.',default='sample')
    parser.add_argument('-d', '--dRanger_input_file', metavar='dRanger_input_file', type=str, help='details.txt.')
    parser.add_argument('-t', '--min_tumor_count_threshold', metavar='min_tumor_count_threshold', type=str, help='min_tumor_count_threshold.',default='4')
    parser.add_argument('-n', '--max_normal_count_threshold', metavar='max_normal_count_threshold', type=str, help='max_normal_count_threshold.',default='1')
    parser.add_argument('-o', '--output', metavar='output', type=str, help='output area', default='./')

    args = parser.parse_args()

    return args

if __name__ == "__main__":

    args = parseOptions()
    pair_id = args.pair_id
    dRanger_input_file = args.dRanger_input_file
    min_tumor_count_threshold = int(args.min_tumor_count_threshold)
    max_normal_count_threshold = int(args.max_normal_count_threshold)
    output = args.output
    if not os.path.exists(dRanger_input_file):
        exit(1)

    df = pd.read_csv(dRanger_input_file,sep='\t')
    df['individual'] = pair_id
    df1=df.loc[df['VCF_NALT']<=max_normal_count_threshold]
    df1=df1.loc[df1['VCF_TALT']>=min_tumor_count_threshold]

    outfile=output + '/' + pair_id + '.somatic.sv.detail.pass.txt'
    df1.to_csv(path_or_buf=outfile, sep='\t', index=False)

    outfile = output + '/' + pair_id + '.somatic.sv.detail.all.txt'
    df.to_csv(path_or_buf=outfile, sep='\t', index=False)
