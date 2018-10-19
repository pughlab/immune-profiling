import fileinput
from optparse import OptionParser
import glob
import os
from threading import Thread

def get_options():
    """
        collect options from the user
        :return: list of command line arguments
        :rtype: list
    """
    parser = OptionParser()
    parser.add_option('-i', '--input',
                      dest="input",
                      type="string",
                      help="mixcr log file path"
                      )
    parser.add_option('-o', '--output',
                      dest="output",
                      type="string",
                      help="directory to save alignment or assembly stats from log"
                      )
    parser.add_option('-s', '--alignments',
                      dest="alignment",
                      default="alginment.txt",
                      type="string",
                      help="aligment stats file name without path"
                      )
    parser.add_option('-a', '--assembly',
                      dest="assembly",
                      default="assembly.txt",
                      type="string",
                      help="assembly stats file name without path"
                      )
    (options, args) = parser.parse_args()
    return options

def process_alignment_stats(log_path):
    """
        parse alignment log file to collect alignment stats
        :param log_path: log file directory
        :type log_path: str
    """
    align_files = glob.glob(os.path.join(log_path, "LOG_ALIGN*.txt"))
    stats_dict = {}
    with fileinput.input(files=align_files) as f:
        for line in f:
            if line.startswith("Analysis Date") or \
                line.startswith("Overlapped") or \
                line.startswith("Analysis time") or \
                line.startswith("Command line arguments") or \
                line.startswith("=="): continue
            parts = line.strip().split(":")
            key, val = parts[0], parts[1]
            if key == 'Input file(s)':
                files = [os.path.basename(f) for f in val.split(",")]
                val = ",".join(files)
            elif key == 'Output file':
                val = os.path.basename(val)
            elif key == 'Version':
                val = val.split(";")[0]
            elif "%)" in val:
                val = val.split()[0]
            print ("%s=>%s" %(key, val))

            if key not in stats_dict:
                stats_dict[key] = [val]
            else:
                stats_dict[key].append(val)
    print (stats_dict)

def process_assembly_stats(files):
    """
        parse assembly log file to collect assembly stats
        :param log_path: log file directory
        :type log_path: str
    """
    align_files = glob.glob(os.path.join(log_path, "LOG_ASSEMBLE*.txt"))
    stats_dict = {}
    with fileinput.input(files=align_files) as f:
        for line in f:
            if line.startswith("Analysis Date") or \
                line.startswith("Analysis time") or \
                line.startswith("Command line arguments") or \
                line.startswith("=="): continue
            parts = line.strip().split(":")
            key, val = parts[0], parts[1]
            if key == 'Input file(s)':
                files = [os.path.basename(f) for f in val.split(",")]
                val = ",".join(files)
            elif key == 'Output file':
                val = os.path.basename(val)
            elif key == 'Version':
                val = val.split(";")[0]
            elif "%)" in val:
                val = val.split()[0]
            print ("%s=>%s" %(key, val))

            if key not in stats_dict:
                stats_dict[key] = [val]
            else:
                stats_dict[key].append(val)
    print (stats_dict)


def main():
    args = get_options()
    print (args)
    t = Thread(target=process_alignment_stats, args=[args.input])
    t.start()
    t = Thread(target=process_assembly_stats, args=[args.input])
    t.start()

if __name__ == "__main__":
    main()