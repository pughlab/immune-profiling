import argparse
import glob
import os
from threading import Thread
from multiprocessing import Queue

def get_options():
    """
        collect options from the user
        :return: list of command line arguments
        :rtype: list
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="path to normalized TPM files")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="output file for expression matrix")
    return parser.parse_args()

def process_expression_file(file):
    """
        save normalized expression value to dictionary
        :param file_name: output file name
        :type file_name: str
    """
    file_name = (os.path.basename(file)).split(".")[0]
    exp_dict = {}
    with open(file) as f:
        file_name = os.path.basename(file).split(".")[0]
        #print("file name: %s"%file_name)
        exp_dict[file_name] = {}
        for line in f.readlines()[1:]:
            items = line.strip().split("\t")
            try:
                exp_dict[file_name][items[0].split(".")[0]] = items[5]
            except:
                continue
        return exp_dict

def write_to_CSV(file_name, sample_ids, tpm_matrix):
    """
        Write normalized expression value to a CSV file
        :param file_name: output file name
        :param sample_ids: sample IDs
        :param tpm_matrix: TPM matrix
        :type file_name: str
        :type sample_ids: list
        :type tpm_matrix: dictionary
    """
    # Write to CSV file
    with open(file_name, 'w') as of:
        of.write("%s\n"%",".join(sample_ids))
        for k, v in tpm_matrix.items():
            of.write("%s,%s\n"%(k, ",".join(v)))

def main():
    args = get_options()
    print (args)

    # Create Thread to parse each TPM file
    que = Queue()
    threads_list = list()

    tpm_files = glob.glob(os.path.join(args.input, "*.txt"))
    for file in tpm_files:
        print ('process %s...' %file)
        t = Thread(target=lambda q, arg1: q.put(process_expression_file(file)), args=(que, file))
        t.start()
        threads_list.append(t)

    # Join all the threads
    for t in threads_list:
        t.join()

    # Process results
    sample_ids = ['esemblID']
    tpm_matrix = {}
    while not que.empty():
        result = que.get()
        for sample_id, tpm_dict in result.items():
            sample_ids.append(sample_id)
            for k, v in tpm_dict.items():
                if k not in tpm_matrix.keys():
                    tpm_matrix[k] = [v]
                else:
                    tpm_matrix[k].append(v)

    # Write to a CSV file
    write_to_CSV(args.output, sample_ids, tpm_matrix)

if __name__ == "__main__":
    main()