#!/usr/bin/env python3

import tompltools
import tompytools
from Bio import SearchIO


# print a line of tab-separated values from the hsp
def write_tsv_lines(hsp, read_length):
    return('%s\t%s\t%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%f\t%f' %
          (hsp.query_id,         # 1
           hsp.hit_id,           # 2
           hsp.hit_description,  # 3
           read_length,          # 4
           hsp.query_span,       # 5
           hsp.query_start,      # 6
           hsp.query_end,        # 7
           hsp.hit_span,         # 8
           hsp.hit_start,        # 9
           hsp.hit_end,          # 10
           hsp.evalue,           # 11
           hsp.bitscore))         # 12


def main():
    # parse io files
    parsed_arguments = tompltools.parse_cli_arguments()
    other_input = parsed_arguments.other_input[0]
    output_table = parsed_arguments.output_table[0]
    # for testing
    # other_input='test/2125-01-11-1_R2.xml'
    # output_table='test/2125-01-11-1_R2.table'

    # header for tsv output
    tsv_header='\t'.join(
        ['query_id',
         'hit_id',
         'hit_description',
         'read_length',
         'query_span',
         'query_start',
         'query_end',
         'hit_span',
         'hit_start',
         'hit_end',
         'evalue',
         'bitscore'])
    tsv_lines = [tsv_header]

    # get a searchIO generator
    tompytools.generate_message('Reading BLAST results from %s' % other_input)
    qresults = SearchIO.parse(other_input, 'blast-xml')

    # loop over the qresults
    tompytools.generate_message('Generating tab-separated output')
    for qresult in qresults:
        read_length = qresult.seq_len
        for hit in qresult:
            for hsp in hit:
                tsv_lines.append(write_tsv_lines(hsp, read_length))

    tompytools.generate_message('Writing tsv to %s' % output_table)
    with open(output_table, 'w') as f:
        f.write('\n'.join(tsv_lines))
        f.write('\n')

if __name__ == '__main__':
    main()
