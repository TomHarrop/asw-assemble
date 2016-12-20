#!/usr/bin/env python3

import tompltools
import tompytools
import gzip
from Bio import SeqIO
from Bio.Blast import NCBIWWW
import tempfile

def main():
    
    # parse the io files
    parsed_arguments = tompltools.parse_cli_arguments()
    input_fq = parsed_arguments.input_fq[0]
    other_output = parsed_arguments.other_output[0]

    # unzip the file transparently
    tompytools.generate_message('Reading fastq input from %s' % input_fq)
    with gzip.open(input_fq, 'rt') as fastq_file:
        # read the fastq
        input_fq_iterator = SeqIO.parse(fastq_file, "fastq")
        fasta_tempfile = tempfile.mkstemp()[1]

        tompytools.generate_message('Writing fasta to temporary file %s'
                                    % fasta_tempfile)
        with open(fasta_tempfile, 'w') as tmp_fasta:
            count = SeqIO.write(input_fq_iterator, tmp_fasta, "fasta")

    # read in as fasta
    fasta_string = open(fasta_tempfile).read()

    # run BLAST query
    tompytools.generate_message('Running BLAST for %s fasta sequences' % count)
    result_handle = NCBIWWW.qblast(
        'blastn',
        'nt',
        fasta_string,
        entrez_query='(all[filter] NOT predicted[title]) AND txid50557[ORGN]',
        #entrez_query='txid50557[ORGN]',
        expect=10, hitlist_size=10)

    # write output
    tompytools.generate_message('Writing XML to %s' % other_output)
    with open(other_output, "w") as f:
        f.write(result_handle.read())

if __name__ == '__main__':
    main()
