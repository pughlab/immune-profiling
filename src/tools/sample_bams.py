from __future__ import division
from multiprocessing import Process
import subprocess
from multiprocessing import Pool
import time
import os
import sys
import argparse
import glob

def get_count(bam, max_workers):
    """
    Count total number of paired reads for a given bam file
    :param bam: input bam files
    :param max_workers: number of threads
    :return: total paired reads from input bam file
    """
    print ("Count total number of paired reads in %s ..."%bam)
    cmd = ['samtools','view','-c','-f', '3','-@',str(max_workers),bam]
    out, err = subprocess.Popen(cmd, stdin = subprocess.PIPE, stdout=subprocess.PIPE).communicate()
    return int(out.split()[0])

def subset(seed, bam, read, output, count, max_workers):
    """
    Down sampling bam file for a given rate
    :param seed: random sampling rate
    :param bam: input bam file
    :param read: dilution read
    :param output: output path
    :param count: total number of paired reads of the bam file
    :param max_workers: number of threads
    :return: None
    """
    bam_name = os.path.basename(bam).split(".")[0]
    sample_rate = round((seed + read/count), 8)
    print (sample_rate)
    sorted_bam = os.path.join(output,"%s_%s_%s_%s.bam"%(bam_name,seed, str(sample_rate).split(".")[1],read))
    cmd = 'samtools view -s %s -f 3 -@ %s -b %s | samtools sort -n > %s'%(sample_rate, max_workers, bam, sorted_bam)

    print (cmd)
    if os.system(cmd) == 0:
        bam2fq_cmd = 'bedtools bamtofastq -i %s -fq %s -fq2 %s'%(sorted_bam,
                                                                 sorted_bam.replace(".bam","-1.fastq"),
                                                                 sorted_bam.replace(".bam", "-2.fastq"))
        print (bam2fq_cmd)
        if os.system(bam2fq_cmd) == 0:
            return sorted_bam
        else:
            print ("Fail to convert sorted bam %s to fastq" %sorted_bam)
            sys.exit(1)
    else:
        print("Fail to downsample bam %s" %bam)
        sys.exit(1)

def execute(seed, bams, reads, output, max_workers, counts):
    """
    Down sampling bam files for given rates and then convert sorted bam files to paired fastq files
    :param seed: random sampling seed
    :param bams: input bam files
    :param reads: comma separated dilution reads
    :param output: directory to store fastq files
    :param max_workers: number of concurrent processes
    :param counts: total number of paired reads in input bam files
    :return: None
    """
    print ("Convert sorted bam to fastq...")
    pool = Pool(processes=max_workers)
    results = [ pool.apply_async(subset, [int(seed), str(bam), int(read), str(output), count, max_workers])
                for bam, read, count in zip(bams,reads, counts) ]
    sub_bams = []
    for async_result in results:
        try:
            bam = async_result.get()
            while bam is None:
                time.sleep(1)
            sub_bams.append(bam)
        except ValueError as e:
           print(e)

    sorted_bams = glob.glob(os.path.join(output, "%s*.bam"%bams[0].split(".")[0]))
    sorted_bams_1 = glob.glob(os.path.join(output, "%s*.bam" % bams[1].split(".")[0]))

    for bam1, bam2 in zip(sorted_bams, sorted_bams_1):
        fastq_name = "{0}_{1}".format(os.path.basename(bam1).split(".")[0],os.path.basename(bam2).split(".")[1])
        cat_cmd = "cat {0} {1} > {2}".format(bam1.replace(".bam","-%s.fastq"%seed),
                                 bam2.replace(".bam","-%s.fastq"%seed),
                                 os.path.join(output, "%s-%s.fastq" %(fastq_name,seed)))
        print (cat_cmd)
        if os.sytem(cat_cmd) != 0:
            print ("Failed to create fasstq {0}.fastq".format(fastq_name))
            sys.exit(1)
    print ("Done!")

def get_options():
    """
    get input options from command line
    :return: command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bams", type=str, required=True,
                        help="Comma separated bam files including path")
    parser.add_argument("-d", "--dilutions", type=str, required=True,
                        help="Comma separated dilutions used for random sampling, e.g 100000000,99000000,98000000,97000000")
    parser.add_argument("-j", "--jobs", type=int, required=False, default=2,
                        help="Max Workers for threading")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Directory for storing output fastq files")
    parser.add_argument('-r', '--reads', type=int, required=False, default=100000000,
                        help="Total reads in final merged bam")
    parser.add_argument("-s", "--seeds", type=str, required=False, default="1,101",
                        help="Range for random sampling, eg. 1,101 is 1-100")
    parser.add_argument("-del", "--delete", action='store_true',
                        help="If set to true, all intermediate sorted Bam files will be removed after fastq files generated")
    return parser.parse_args()

if __name__ == '__main__':
    options = get_options()
    print (options)
    dilutions = [ int(read) for read in options.dilutions.split(",") ]
    start,end = [ int(x) for x in options.seeds.split(",") ]
    bams = [ str(bam) for bam in options.bams.split(",") ]
    counts = [ get_count(bam, options.jobs) for bam in bams ]
    for read in dilutions:
        reads = [read, options.reads - read]
        print ("Reads need to be extracted: %s" %reads)
        for seed in range(start, end):
           print ("Seed : %s" %seed)
           p = Process(target=execute, args=(seed, bams, reads, options.output, options.jobs, counts))
           p.start()


