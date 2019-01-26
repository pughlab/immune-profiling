from __future__ import division
from multiprocessing.pool import ThreadPool
import subprocess
import os
import argparse
import re
import time

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

def create_script(sh_file, cmds, max_workers, num_nodes=1):
    """
    generate a cluster shell script for downsampling bam file
    :param sh_file: qsub script name
    :param cmds: commands to be executed on cluster
    :param max_workers: maximum number of processes for a job
    :param num_nodes: number of cluster nodes required
    :return: None
    """
    output = os.path.dirname(sh_file)
    job_name = os.path.splitext(os.path.basename(sh_file))[0]
    err_file = os.path.join(output,"{0}.error".format(job_name))
    complete_file = os.path.join(output, "{0}.complete".format(job_name))
    with open(sh_file, 'w') as of:
        of.write("#!/bin/bash\n")
        of.write("#PBS -N {0}\n".format(job_name))
        of.write("#PBS -l nodes={0}:ppn={1}\n".format(num_nodes,max_workers))
        of.write("#PBS -l walltime=2:30:00\n")
        of.write("#PBS -l vmem=8g\n")
        of.write("#PBS -j eo\n")
        of.write("#PBS Join_Path={0}\n".format(os.path.join(output,"%s.err"%job_name)))
        of.write("module load samtools/1.9\n")
        of.write("module load bedtools/2.27.1\n")
        of.write("{0}\n".format(cmds[0]))
        of.write("if [ $? -ne 0 ]; then \n\ttouch {0};exit 1 \nfi\n".format(err_file))
        of.write("{0}\n".format(cmds[1]))
        of.write("if [ $? -ne 0 ]; then \n\ttouch {0}\nelse\n\ttouch {1} \nfi\n".format(err_file, complete_file))
    os.system("chmod 755 %s" % sh_file)

def create_master_script(master_file, cmds, num_files, max_workers, num_nodes, remove_files=False):
    """

    :param master_file: master script to monitor all downsampling jobs
    :param cmds: commands to be executed on cluster
    :param num_files: number of expected fastq files
    :param max_workers: number of processes for a job
    :param num_nodes: number of cluster nodes
    :param remove_files: if set, all intermediate bam files will be removed after fastq files generated
    :return: None
    """
    output = os.path.dirname(master_file)
    job_name = os.path.splitext(os.path.basename(master_file))[0]
    with open(master_file, 'w') as of:
        of.write("#!/bin/bash\n")
        of.write("#PBS -N {0}\n".format(job_name))
        of.write("#PBS -l nodes={0}:ppn={1}\n".format(num_nodes,max_workers))
        of.write("#PBS -l walltime=1:30:00\n")
        of.write("#PBS -l vmem=8g\n")
        of.write("#PBS -o {0}\n".format(os.path.join(output,"%s.out"%job_name)))
        of.write("#PBS -e {0}\n".format(os.path.join(output, "%s.err" %job_name)))
        of.write('find %s -name "*.sh" -exec qsub {} \;\n'%(output))
        of.write("fail=$(ls -l %s/*.error |wc -l) || fail=0\n"%output)
        of.write("complete=$(ls -l {0}/*.complete |wc -l) || complte=0\n".format(output))
        of.write("if [ $fail -ne 0 ]; then\n exit 1 \nfi\n")
        of.write("while [ $complete -lt %s ];do\n sleep 5\ndone\n"%num_files)
        for i in range(0, len(cmds)):
            of.write("{0}\n".format(cmds[i]))
        if remove_files:
            of.write("rm {0}".format(os.path.join(output,"*.bam")))
    os.system("chmod 755 %s" % master_file)

def subset(seed, bam, read, output, count, max_workers, num_nodes, qsub_dir):
    """
    Down sampling bam file for a given sampling rate
    :param seed: seed for calculate random sampling rate
    :param bam: input bam file
    :param read: dilution read
    :param output: output path
    :param count: total number of paired reads of the bam file
    :param max_workers: number of threads for subset a bam
    :param num_nodes: number of cluster nodes
    :return: name sorted bam name
    """
    bam_name = os.path.basename(bam).split(".")[0]
    sample_rate = round((seed + read/count), 8)
    sorted_bam = os.path.join(output,"%s_%s_%s_%s.bam"%(bam_name,seed, str(sample_rate).split(".")[1],read))
    cmds = list()
    cmds.append('samtools view -s {0} -f 3 -@ {1} -b {2} | samtools sort -n -T {3} > {4}'.format(sample_rate,
                                                                                max_workers, bam, output, sorted_bam))
    cmds.append('bedtools bamtofastq -i {0} -fq {1} -fq2 {2}'.format(sorted_bam,
                                                                 sorted_bam.replace(".bam","-1.fastq"),
                                                                 sorted_bam.replace(".bam", "-2.fastq")))
    # creating qsub script
    create_script(os.path.join(qsub_dir,os.path.basename(sorted_bam).replace("bam","sh")),
                  cmds, max_workers, num_nodes)

    return os.path.basename(sorted_bam)

def check_qsub_path(script_dir):
    """
    create a new directory if not already exists
    :param script_dir: directory for storing shell scripts
    :return: None
    """
    if not os.path.exists(script_dir):
        try:
            os.mkdir(script_dir, 0o755)
        except OSError as e:
            print("Error:Script directory exists")

def get_paired_bam(sorted_bams, one_bam, paired_bam_name, read):
    """
    find paired sorted bam file
    :param sorted_bams: name sorted bam files
    :param one_bam: bam file
    :param paired_bam_name: paired bam file name
    :param read: dilution reads
    :return: paired bam name
    """
    try:
        print (sorted_bams)
        seed, rate, dilution = (os.path.splitext(one_bam))[0].split("_")[-3:]
        pattern = "{0}_{1}.*{2}".format(paired_bam_name,seed, str(read-int(dilution)))
        print ("search pattern: {0}".format(pattern))
        res = re.compile(pattern)
        paired_bams = list(filter(res.match, sorted_bams))
        return paired_bams[0]
    except:
        return None

def get_cat_cmd(bam, paired_bam, output, fq_name, remove_files=False):
    cat_cmds = []
    for i in range(1, 3):
        cat_cmds.append("cat {0} {1} > {2}".format(os.path.join(output,bam.replace(".bam","-%s.fastq"%i)),
                             os.path.join(output,paired_bam.replace(".bam","-%s.fastq"%i)),
                             os.path.join(output, "%s-%s.fastq" %(fq_name,i))))

        if remove_files:
            cat_cmds.append("if [ $? -eq 0 ]; then\n\trm {0} {1}\nfi\n".format(bam.replace(".bam","-%s.fastq"%i),
                                                paired_bam.replace(".bam","-%s.fastq"%i)))
    return cat_cmds

def get_bam_name_parts(one_bam):
    one_bam_parts = os.path.splitext(one_bam)[0]
    one_bam_parts = one_bam_parts.split("_")
    seed, rate, dilution = one_bam_parts[-3:]
    del one_bam_parts[-3:]
    if len(one_bam_parts) > 1:
        one_bam_name = "_".join(one_bam_parts)
    else:
        one_bam_parts = one_bam_parts[0]

    return one_bam_name, seed, rate, dilution

def get_fastq_name(one_bam, paired_bam, total_reads):
    one_bam_name, seed, rate1, dilution1 = get_bam_name_parts(one_bam)
    paired_bam_name, seed, rate2, dilution2 = get_bam_name_parts(paired_bam)
    return  ("{0}_{1}_{2}_{3}_{4}".format(one_bam_name, paired_bam_name, seed, dilution1, total_reads),
            seed, "%s,%s"%(dilution1,dilution2), "%s,%s"%(rate1,rate2))

def get_options():
    """
    get input options from command line
    :return: commandline arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bams", type=str, required=True,
                        help="Comma separated bam files including path")
    parser.add_argument("-d", "--dilutions", type=str, required=True,
                        help="Comma separated dilutions used for random sampling, e.g 100000000,99000000,98000000,97000000")
    parser.add_argument("-j", "--jobs", type=int, required=False, default=8,
                        help="Max Workers for threading")
    parser.add_argument("-n", "--nodes", type=int, required=False, default=1,
                        help="Number of cluster nodes for a job")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Directory for storing output fastq files")
    parser.add_argument("-q", "--qsub", type=str, required=False, default="qsub_script",
                        help="Directory for storing qsub scripts")
    parser.add_argument('-r', '--reads', type=int, required=False, default=100000000,
                        help="Total reads in final merged bam")
    parser.add_argument("-s", "--seeds", type=str, required=False, default="1,101",
                        help="Range for random sampling, eg. 1,101 is 1-100")
    parser.add_argument("-rf", "--remove_files", action='store_true',
                        help="If set, all intermediate sorted Bam files will be removed after fastq files generated")
    return parser.parse_args()


if __name__ == '__main__':
    #collect arguments
    options = get_options()
    print ("user options:{0}".format(options))
    dilutions = [ int(read) for read in options.dilutions.split(",") ]
    start,end = [ int(x) for x in options.seeds.split(",") ]
    bams = [ str(bam) for bam in options.bams.split(",") ]
    counts = [ get_count(bam, options.jobs) for bam in bams ]

    #create a new directory for storing qsub scripts if not already exists
    script_dir = os.path.join(options.output, options.qsub)
    check_qsub_path(script_dir)

    #generate qsub scripts...
    pool = ThreadPool(processes=options.jobs)
    workers = list()
    for read in dilutions:
        reads = [read, options.reads - read]
        print ("Reads need to be extracted: %s" %reads)
        for seed in range(start, end):
            print ("Seed : %s" %seed)
            for one_bam, read, count in zip(bams, reads, counts):
                workers.append(pool.apply_async(subset, (seed,one_bam,read,options.output, count,
                                                         options.jobs,options.nodes,script_dir)))

    sorted_bams = []
    bam_names = [ list(os.path.splitext(os.path.basename(b)))[0] for b in bams ]

    for t in workers:
        while (t.get() == None):
            time.sleep(5)
        sorted_bams.append(t.get())
    cmds = []

    # create log file
    log_file = os.path.join(options.output,"log.txt")
    with open(log_file,'w') as of:
        of.write("merge_fq\tbams\tseed\tdilutions\ttotal_reads\tsample_rates\n")
        for one_bam in sorted_bams:
            if bam_names[0] not in one_bam: continue
            paired_bam = get_paired_bam(sorted_bams, one_bam, bam_names[1], options.reads)
            print (one_bam, paired_bam)
            if paired_bam:
                fq_name, seed_str, dil_str, rate_str = get_fastq_name(one_bam, paired_bam, options.reads)
                print ("concat fastq name: {}".format(fq_name))
                cmds += get_cat_cmd(one_bam, paired_bam, options.output, fq_name, options.remove_files)

                of.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(
                                    fq_name,
                                    '%s,%s'%(os.path.basename(one_bam),os.path.basename(paired_bam)),
                                    seed_str,
                                    dil_str,
                                    options.reads,
                                    rate_str))

    #create master qsub script
    master_file = os.path.join(options.output, "{0}-{1}.sh".format(bam_names[0],bam_names[1]))
    print (master_file)
    create_master_script(master_file, cmds, len(sorted_bams), options.jobs, options.nodes, options.remove_files)
