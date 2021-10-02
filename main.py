import argparse, logging, os
from subprocess import run
from pathlib import Path
from datetime import datetime
from contextlib import contextmanager

CURR_TIME = datetime.now().strftime("%d-%m-%Y_%H-%M-%S")
ENDING = ['.fq.gz', '.fastq.gz', '.fastq', '.fq', '.fasta']

@contextmanager
def working_directory(directory):
    owd = os.getcwd()
    try:
        os.chdir(directory)
        yield directory
    finally:
        os.chdir(owd)

def check_index(path_index: str):
    if Path(path_index).is_file():
        pass
    else:
        exit(f"No index file found on {path_index}")

    return 

def add_fwd_slash(args: argparse.Namespace):
    if any(end in args.index[0] for end in ENDING):
        pass
    elif not any(end in args.index[0] for end in ENDING) \
        and args.index[0][-1] != "/":
        args.index[0] += "/"
    
    if args.dir[0][-1] != "/":
        args.dir[0] += "/"

    return args

def build_directory(args: argparse.Namespace):
    run(["mkdir", "-p", "results_" + CURR_TIME + "/{1_quality_control,2_trimmed_output,3_kallisto_results"])
    d = run(["pwd"], capture_output=True, text=True)[0] + "/results_" + CURR_TIME + "/"
    args.output[0] = d

    return args

def create_index(args: argparse.Namespace):
    if "/" in args.transcript[0]:
        idx_name = args.transcript[0].split("/")[-1].split(".")[0]
    else:
        idx_name = args.transcript[0].split(".")[0]
        
    idx = run(["kallisto", "index", args.index + idx_name + ".idx", \
        args.transcript], capture_output=True, text=True)

    if idx.stdout:
        print(idx.stdout)
        logging.info(idx.stdout)
    if idx.stderr:
        print(idx.stderr)
        logging.info(idx.stderr)
    
    return args.index + idx_name + ".idx"

def run_qctk(args: argparse.Namespace):
    if not args.single:
        for sample in args.samples:
            qc = run(["fastqc", "-o", args.output[0] + "1_quality_control", "--no-extract",
                args.dir[0] + sample + args.complement[0] + ".fq.gz",
                args.dir[0] + sample + args.complement[1] + ".fq.gz"],
                    capture_output=True, text=True)
            if qc.stdout:
                print(qc.stdout)
                logging.info(qc.stdout)
            if qc.stderr:
                print(qc.stderr)
                logging.info(qc.stderr)

            trim = run(["trim_galore", "--quality", "20", "--fastqc", "--length", "25", "--paired",
                "--output-dir", args.output[0] + "2_trimmed_output",
                    args.dir[0] + sample + args.complement[0] + ".fq.gz",
                    args.dir[0] + sample + args.complement[1] + ".fq.gz"],
                        capture_output=True, text=True)
            if trim.stdout:
                print(trim.stdout)
                logging.info(trim.stdout)
            if trim.stderr:
                print(trim.stderr)
                logging.info(trim.stderr)

            kall = run(["kallisto", "quant", "-t", args.threads[0], "-b", args.bootstrap[0],
                "-i", args.index, "-o", args.output[0] + "3_kallisto_output/" + sample,
                args.output[0] + "2_trimmed_output/" + sample + args.complement[0] + "_val_1.fq.gz",
                args.output[0] + "2_trimmed_output/" + sample + args.complement[1] + "_val_2.fq.gz"],
                capture_output=True, text=True)
            if kall.stdout: 
                print(kall.stdout)
                logging.info(kall.stdout)
            if kall.stderr: 
                print(kall.stderr)
                logging.info(kall.stderr)

    elif args.single:
        for sample in args.samples:
            qc = run(["fastqc", "-o", args.output[0] + "1_quality_control", "--no-extract",
                args.dir[0] + sample + ".fq.gz"],
                    capture_output=True, text=True)
            if qc.stdout:
                print(qc.stdout)
                logging.info(qc.stdout)
            if qc.stderr:
                print(qc.stderr)
                logging.info(qc.stderr)

            trim = run(["trim_galore", "--quality", "20", "--fastqc", "--length", "25",
                "--output-dir", args.output[0] + "2_trimmed_output",
                    args.dir[0] + sample + ".fq.gz"],
                        capture_output=True, text=True)
            if trim.stdout:
                print(trim.stdout)
                logging.info(trim.stdout)
            if trim.stderr:
                print(trim.stderr)
                logging.info(trim.stderr)

            kall = run(["kallisto", "quant", "-t", args.threads[0], "-b", args.bootstrap[0],
                "--single", "-i", args.index, "-o", args.output[0] + "3_kallisto_output/" + sample,
                args.output[0] + "2_trimmed_output/" + sample + "_val_1.fq.gz"],
                capture_output=True, text=True)
            if kall.stdout: 
                print(kall.stdout)
                logging.info(kall.stdout)
            if kall.stderr: 
                print(kall.stderr)
                logging.info(kall.stderr)

    return print("Finished pseudoalignment!")

def read_samples(args: argparse.Namespace) -> dict:
    if args.simple & args.complement is not None:
        logging.info("Single-ended analysis does not contain complements. \
            Complements are for paired-ended (e.g. sample1_R1.fq.gz and sample1_R2.fq.gz)")
        print("Single-ended analysis does not contain complements. \
            Complements are for paired-ended (e.g. sample1_R1.fq.gz and sample1_R2.fq.gz)")
        exit()
    elif not args.simple & args.complement is not None:
        if len(args.complement) > 2:
            logging.info("Complement (-c or --complement) argument has to be maximum \
                of 2 for paired-ended reads.")
            print("Complement (-c or --complement) argument has to be maximum \
                of 2 for paired-ended reads.")
            exit()

    if not Path(args.output[0]).is_dir():
        logging.info("Directory {0} does not exist.".format(args.output[0]))
        print("Directory {0} does not exist.".format(args.output[0]))
        exit()
    
    if not Path(args.dir[0]).is_dir():
        logging.info("Directory {0} does not exist.".format(args.dir[0]))
        print("Directory {0} does not exist.".format(args.dir[0]))
        exit()

    if not args.simple:
        for file in args.samples:
            if Path(args.dir[0] + file + args.complement[0] + ".fq.gz").is_file() \
                & Path(args.dir[0] + file + args.complement[1] + ".fq.gz").is_file():
                    pass
            else:
                logging.info("File {0} does not exist in {1}.".format(file, str(args.dir[0])))
                print("File {0} does not exist in {1}.".format(file, str(args.dir[0])))
                exit()
    elif args.simple:
        for file in args.samples:
            if Path(args.dir[0] + file + ".fq.gz").is_file():
                pass
            else:
                logging.info("File {0} does not exist in {1}.".format(file, str(args.dir[0])))
                print("File {0} does not exist in {1}.".format(file, str(args.dir[0])))
                exit()
    
    print("All files exists. Continuing the analysis.")

    return args

def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Kallisto for every samples given.")
    parser.add_argument("-s", "--samples", nargs="+", required=True, 
        help="<Required> List of samples to iterate over", type=str)
    parser.add_argument("-c", "--complement", nargs="+", required=False, 
        help="<Optional> Complementary for paired-ended", type=str)
    parser.add_argument("-d", "--dir", nargs=1, required=True, 
        help="<Required> Directory to files", type=str, metavar="dir/to/files")
    parser.add_argument("-i", "--index", nargs=1, required=True, 
        help="<Required> Directory to index already built or path to build the index file", type=str)
    parser.add_argument("-t", "--transcript", nargs=1, required=False, 
        help="<Optional> Path to transcript file, has to be passed along with -i/--index")
    parser.add_argument("--threads", nargs=1, required=False,
        help="<Optional> Number of threads to be used in quantification for Kallisto. Default: 1.")
    parser.add_argument("-b", "--bootstrap", nargs=1, required=False,
        help="<Optional> Number of bootstrap samples. Default: 0.")
    parser.add_argument("--single", action='store_true', required=False,
        help="<Optional> Flag to indicate single-ended quantification without complements.")

    args = parser.parse_args()

    # Creating and logging info for the current run
    logging.basicConfig(filename= CURR_TIME + ".log", filemode="a", \
        format='%(asctime)s %(message)s', datefmt='%d/%m/%Y %H:%M:%S', \
        level=logging.DEBUG)
    logging.info("Samples used: {0}".format(args.samples))
    logging.info("Complements: {0}".format(args.complement))
    logging.info("Working directory for samples: {0}".format(args.dir))

    print("Samples used: {0}".format(args.samples))
    print("Complements: {0}".format(args.complement))
    print("Working directory for samples: {0}".format(args.dir))
    
    r = False
    while r not in ['y', 'n']:
        r = str(input("\nIs this correct? [y/n]\n")).lower()
        if r in ['y', 'n']:
            if r == 'y':
                pass
            elif r == 'n':
                print("Arguments not right. Exiting.")
                exit()
        else:
            print("\nType `y` for Yes or `n` for No")
    
    print("Every argument has been checked. Continuing the analysis.\n")
    
    return args

if __name__ == '__main__':
    arguments = arguments()
    arguments = add_fwd_slash(arguments)
    arguments = build_directory(arguments)
    fargs = read_samples(args=arguments)
    if fargs.transcript and fargs.index:
        idx = create_index(fargs)
        check_index(idx)
        fargs.index = idx
    elif fargs.index and fargs.transcript is None:
        check_index(fargs.index)
    else:
        exit("Neither index or transcript file has been passed")

    run_qctk(fargs)