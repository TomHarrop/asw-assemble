#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

################################
# ASW genome assembly pipeline #
################################

################
# Requirements #
################

import tompltools
import tompytools
import ruffus
import os

############
# Pipeline #
############

def main():

    #########
    # SETUP #
    #########

    # test function for checking input/output passed to job_script and parsing
    # by src/sh/io_parser
    test_job_function = tompltools.generate_job_function(
        job_script='src/sh/io_parser',
        job_name='test',
        verbose=True)

    # parse email etc. here?
    parser = ruffus.cmdline.get_argparse(
        description='ASW genome assembly pipeline.')
    # parser.add_argument('--email', '-e',
    #                     help='Logon email address for JGI',
    #                     type=str,
    #                     dest='jgi_logon')
    # parser.add_argument('--password', '-p',
    #                     help='JGI password',
    #                     type=str,
    #                     dest='jgi_password')
    options = parser.parse_args()
    # jgi_logon = options.jgi_logon
    # jgi_password = options.jgi_password

    # initialise pipeline
    main_pipeline = ruffus.Pipeline.pipelines['main']

    # find fastq.gz files
    dir_listing = [x[0] for x in os.walk(top='data', followlinks=True)]
    fastq_file_list = []
    for directory in dir_listing:
        file_list = os.scandir(directory)
        fastq_file_list.append([x.path for x in file_list
                               if (x.name.endswith('fastq.gz')
                                   or x.name.endswith('.fastq'))
                               and x.is_file()])

    fastq_files = list(tompytools.flatten_list(fastq_file_list))

    # extract only ASW gDNA fastq data, i.e.
    # 2125-01-11-1 = ASW PE
    # 2125-01-06-1 = ASW MP
    active_fq_files = [x for x in fastq_files
                       if ('2125-01-11-1' in x
                           or '2125-01-06-1' in x)]

    # load files into ruffus 
    raw_fq_files = main_pipeline.originate(
        name='raw_fq_files',
        task_func=os.path.isfile,
        output=active_fq_files)

    # merge libraries
    merged_fq_files = main_pipeline.transform(
        name='merged_fq_files',
        task_func=tompltools.generate_job_function(
            job_script='src/sh/merge_fq',
            job_name='merge_fq'),
        input=raw_fq_files,
        filter=ruffus.formatter(
            r'data/NZGL02125/.*/'
            '[^-]+-(?P<LIB>[^_]+).+_R(?P<RN>\d)_.*.fastq.gz'),
        output=[r'output/fq_merged/{LIB[0]}_R{RN[0]}_merged.fastq.gz'])

    # make pairs and send to cutadapt for trimming external adaptors
    trim_cutadapt = main_pipeline.collate(
        name='trim_cutadapt',
        task_func=tompltools.generate_job_function(
            job_script='src/sh/cutadapt_pe',
            job_name='cutadapt'),
        input=merged_fq_files,
        filter=ruffus.formatter(
            r'.+/(?P<LIB>[^_]+)_R(?P<RN>\d)_merged.fastq.gz'),
        output=[['output/cutadapt/{LIB[0]}_R1_trimmed.fastq.gz',
                'output/cutadapt/{LIB[0]}_R2_trimmed.fastq.gz']])

    # send trimmed reads to splitnextera
    mp_splitnextera = main_pipeline.subdivide(
        name='splitnextera',
        task_func=tompltools.generate_job_function(
            job_script='src/sh/splitnextera',
            job_name='splitnextera'),
        input=trim_cutadapt,
        filter=ruffus.regex(
            r'.+?/2125-01-06-1_R(?P<RN>\d)_trimmed.fastq.gz'),
        output=['output/splitnextera/2125-01-06-1.pe.fastq.gz',
                'output/splitnextera/2125-01-06-1.se.fastq.gz',
                'output/splitnextera/2125-01-06-1.mp.fastq.gz',
                'output/splitnextera/2125-01-06-1.unknown.fastq.gz'])

    # decontaminate PhiX (other?) sequences
    decon_mp = main_pipeline.transform(
        name='decon_mp',
        task_func=tompltools.generate_job_function(
            job_script='src/sh/decon',
            job_name='decon_mp'),
        input=mp_splitnextera,
        filter=ruffus.formatter(
            r'.+/2125-01-06-1\.(?P<VL>[^.]+)\.fastq.gz'),
        output=['output/decon/2125-01-06-1_{VL[0]}.fastq.gz'])

    decon_pe = main_pipeline.collate(
        name='decon_pe',
        task_func=tompltools.generate_job_function(
            job_script='src/sh/decon',
            job_name='decon_pe'),
        input=trim_cutadapt,
        filter=ruffus.regex(
            r'.+?/2125-01-11-1_R(?P<RN>\d)_trimmed.fastq.gz'),
        output=[r'output/decon/2125-01-11-1.fastq.gz'])

    decon = [decon_mp, decon_pe]

    # digital normalisation and error correction w/ bbnorm
    bbnorm = main_pipeline.subdivide(
        name='bbnorm',
        task_func=tompltools.generate_job_function(
            job_script='src/sh/bbnorm',
            job_name='bbnorm',
            mem_per_cpu=7000,
            cpus_per_task=8),
        input=decon,
        filter=ruffus.formatter(r'.+/(?P<LN>[^(_|.)]+)(?P<VL>_?\w*).fastq.gz'),
        output=[r'output/bbnorm/{LN[0]}{VL[0]}.fastq.gz'])

    # print raw and normalised kmer distribution plots
    main_pipeline.merge(
        name='kmer_distribution_plots',
        task_func=tompltools.generate_job_function(
            job_script='src/r/kmer_distribution_plots.R',
            job_name='kmer_distribution_plots'),
        input=bbnorm,
        output=['output/bbnorm/plots.pdf', 'output/bbnorm/plot_data.Rds'])

    # run fastqc on decontaminated and normalised libraries
    main_pipeline.transform(
        name='fastqc',
        task_func=tompltools.generate_job_function(
            job_script='src/sh/fastqc',
            job_name='fastqc',
            cpus_per_task=1),
        input=bbnorm,
        filter=ruffus.formatter(r'.+/(?P<LN>[^(_|.)]+)(?P<VL>_?\w*).fastq.gz'),
        output=[r'output/fastqc/{LN[0]}{VL[0]}_fastqc.html'])

    # trim reads to 100 bp for edena?
    clip_to_100b = main_pipeline.subdivide(
        name='clip_to_100b',
        task_func=tompltools.generate_job_function(
            job_script='src/sh/clip_to_100b',
            job_name='clip_to_100b'),
        input=bbnorm,
        # next line is for ALL libraries
        # filter=ruffus.formatter(r'.+/(?P<LN>[^(_|.)]+)(?P<VL>_?\w*).fastq.gz'),
        # next line is for ONLY PE ThruPLEX library
        filter=ruffus.regex(r'.+/2125-01-11-1.fastq.gz'),
        output=[r'output/trunc_100/2125-01-11-1_R1.fastq.gz',
                r'output/trunc_100/2125-01-11-1_R2.fastq.gz'])

    # subsample reads, blast with biopython and parse results
    fq_subsample = main_pipeline.subdivide(
        name='fq_subsample',
        task_func=tompltools.generate_job_function(
            job_script='src/sh/fq_subsample',
            job_name='fq_subsample'),
        input=bbnorm,
        filter=ruffus.formatter(r'.+/(?P<LN>[^(_|.)]+)(?P<VL>_?\w*).fastq.gz'),
        output=[r'output/subsample/{LN[0]}{VL[0]}_R1.fastq.gz',
                r'output/subsample/{LN[0]}{VL[0]}_R2.fastq.gz'])
    blast_reads = main_pipeline.transform(
        name='blast_reads',
        task_func=test_job_function,
        input=fq_subsample,
        filter=ruffus.suffix('.fastq.gz'),
        output=['.xml'])
    main_pipeline.transform(
        name='parse_blast_results',
        task_func=test_job_function,
        input=blast_reads,
        filter=ruffus.suffix('.xml'),
        output=['.table'])

    # overlap step with edena
    edena_overlaps = main_pipeline.collate(
        name='edena_overlaps',
        task_func=tompltools.generate_job_function(
            job_script='src/sh/edena_overlaps',
            job_name='edena_overlaps'),
        input=clip_to_100b,
        filter=ruffus.formatter(r'.+/(?P<LN>[^_]+)_R\d.fastq.gz'),
        output=[r'output/edena/{LN[0]}.ovc'])

    # prepare files with velveth
    # set threads for velvet to 1 !!!
    min_kmer = 71
    max_kmer = 87
    step = 8
    kmer_lengths = [x for x in range(min_kmer, max_kmer + 1, step)]
    velveth_output = list(
        tompytools.flatten_list(
            [('output/velveth_' + str(x) + '/Sequences')
             for x in kmer_lengths]))
    velveth = main_pipeline.merge(
        name='hash_files',
        task_func=test_job_function,
        # task_func=tompltools.generate_job_function(
        #     job_script='src/sh/hash_files',
        #     job_name='hash_files'),
        input=bbnorm,
        output=velveth_output)

    # # run velvetg on hash files
    # velvet_assembly = main_pipeline.subdivide(
    #     name='velvet_assembly',
    #     task_func=test_job_function,
    #     input=velveth_output,
    #     filter=ruffus.suffix("Sequences"),
    #     output=["Log"])\
    #     .follows(hash_files)

    ###################
    # RUFFUS COMMANDS #
    ###################

    # print the flowchart
    ruffus.pipeline_printout_graph(
        "ruffus/flowchart.pdf", "pdf",
        pipeline_name="ASW genome assembly pipeline")

    # run the pipeline
    ruffus.cmdline.run(options, multithread=8)

if __name__ == "__main__":
    main()
