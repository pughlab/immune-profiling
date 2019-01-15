from __future__ import division
from multiprocessing import Process
import subprocess
from multiprocessing import Pool, TimeoutError
import time
import os
import argparse

def info(title):
    print (title)
    print ('module name:', __name__)
    if hasattr(os, 'getppid'):  # only available on Unix
        print ('parent process:', os.getppid())
    print ('process id:', os.getpid())

def get_count(bam, max_workers):
    print ("count total number of reads in %s ..."%bam)
    cmd = ['sambamba','view','-c','-t',str(max_workers),bam]
    out, err = subprocess.Popen(cmd, stdin = subprocess.PIPE, stdout=subprocess.PIPE).communicate()
    return int(out.split()[0])

def subset(seed, bam, read, output, count, max_workers):
    bam_name = os.path.basename(bam).split(".")[0]
    sample_rate = round((seed + read/count), 8)
    print (sample_rate)
    output = os.path.join(output,"%s_%s_%s_%s.bam"%(bam_name,seed, str(sample_rate).split(".")[1],read))
    cmd = 'samtools view -s %s -b %s | samtools sort > %s'%(sample_rate, bam, output) 

    print (cmd)
    if os.system(cmd) == 0:
        return output
    else:
        return None

def execute(seed, bams, reads, output, max_workers, counts):
    print ("total number of reads in the bams: %s" %counts)
    pool = Pool(processes=max_workers)
    results = [ pool.apply_async(subset, [int(seed), str(bam), int(read), str(output), count, max_workers]) for bam, read, count in zip(bams,reads, counts) ]
    sub_bams = []
    for async_result in results:
        try:
            bam = async_result.get()
            while bam is None:
                time.sleep(1)
            sub_bams.append(bam)
        except ValueError as e:
           print(e)

    bam_names = [os.path.basename(one_bam).split(".")[0] for one_bam in sub_bams] 
    merged_bam = os.path.join(os.path.dirname(sub_bams[0]),"%s.bam"%"_".join(bam_names))
    merge_cmd = 'java -jar $picard_dir/picard.jar MergeSamFiles I=%s I=%s O=%s'%(sub_bams[0], sub_bams[1], merged_bam)
    if os.system(merge_cmd) == 0:
        #delete subset bams
        for one_bam in sub_bams:
            os.system("rm -f %s"%one_bam)

        bam2fq_cmd = 'bedtools bamtofastq -i %s -fq %s -fq2 %s'%(merged_bam,
                                                                 merged_bam.replace(".bam","-1.fastq"),
                                                                 merged_bam.replace(".bam", "-2.fastq"))
        os.system(bam2fq_cmd)

def get_options():

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
    return parser.parse_args()

if __name__ == '__main__':
    options = get_options()
    print (options)
    dilutions = [ int(read) for read in options.dilutions.split(",") ]
    start,end = [ int(x) for x in options.seeds.split(",") ]
    bams = [ str(bam) for bam in options.bams.split(",") ]
    counts = [ get_count(bam, options.jobs) for bam in bams ]
    # dilutions = [[100000000, 0], [99000000, 1000000],[98000000, 2000000], [97000000, 3000000],
    #                 [96000000, 4000000], [95000000, 5000000], [94000000, 6000000],[93000000, 7000000],
    #                 [92000000, 8000000], [91000000, 9000000], [90000000, 10000000], [8000000, 12000000],
    #                 [6000000, 14000000], [4000000, 16000000],[2000000, 18000000],[0, 100000000]]
    for read in dilutions:
        reads = [read, options.reads - read]
        print ("reads to be extracted: %s" %reads)
        for seed in range(start, end):
           print ("seed : %s" %seed)
           p = Process(target=execute, args=(seed, bams, reads, options.output, options.jobs, counts))
           p.start()


