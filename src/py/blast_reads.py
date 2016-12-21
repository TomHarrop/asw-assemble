#!/usr/bin/env python3

import tompltools
import tompytools
import gzip
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
import tempfile
import os


def main():

    # make sure blastdb location is set
    if os.getenv('BLASTDB'):
        print('Using BLAST database folder: %s' % os.getenv('BLASTDB'))
    else:
        raise EnvironmentError("BLASTDB environment variable not set")

    # how many cpus are available?
    if os.getenv('SLURM_JOB_CPUS_PER_NODE'):
        max_cpus = int(os.getenv('SLURM_JOB_CPUS_PER_NODE'))
    else:
        max_cpus = 1
    tompytools.generate_message('Using %s CPU(s)' % max_cpus)

    # parse the io files
    parsed_arguments = tompltools.parse_cli_arguments()
    input_fq = parsed_arguments.input_fq[0]
    other_output = parsed_arguments.other_output[0]

    # unzip the file transparently
    tompytools.generate_message('Reading fastq input from %s' % input_fq)
    with gzip.open(input_fq, 'rt') as fastq_file:
        # read the fastq
        input_fq_iterator = SeqIO.parse(fastq_file, "fastq")
        fasta_tempfile = tempfile.mkstemp(suffix=".fasta")[1]

        tompytools.generate_message('Writing fasta to temporary file %s'
                                    % fasta_tempfile)
        with open(fasta_tempfile, 'w') as tmp_fasta:
            count = SeqIO.write(input_fq_iterator, tmp_fasta, "fasta")

    # read in as fasta
    # fasta_string = open(fasta_tempfile).read()

    # run BLAST query
    tompytools.generate_message('Running BLAST for %s fasta sequences' % count)
    blastn_cline = NcbiblastnCommandline(
        query=fasta_tempfile,
        out=other_output,
        outfmt=5,
        db='nt',
        evalue=1,
        max_target_seqs=10,
        num_threads=max_cpus,
        task='blastn')
    print(blastn_cline)
    blastn_cline()

    # write output
    # tompytools.generate_message('Writing XML to %s' % other_output)
    # with open(other_output, "w") as f:
    #     f.write(result_handle.read())

if __name__ == '__main__':
    main()
