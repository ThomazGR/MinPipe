import argparse
from subprocess import run
from pathlib import Path
from datetime import datetime
import logging

def run_kallisto(args: argparse.Namespace):
    for sample in args.samples:
        o = run(["kallisto", "quant", "-t", "4", "-b", "100","-i", "MMUM37.idx", "-o", args.output[0] + "/" + sample, \
            args.dir[0] + "/" + sample + args.complement[0] + ".fq.gz", \
            args.dir[0] + "/" + sample + args.complement[1] + ".fq.gz"] \
            , capture_output=True, text=True)
        if o.stdout: 
            print(o.stdout)
            logging.info(o.stdout)
        if o.stderr: 
            print(o.stderr)
            logging.info(o.stderr)

    return print("Finished!")

def read_samples(args: argparse.Namespace) -> dict:
    if len(args.complement) > 2:
        logging.info("Complement (-c or --complement) argument has to be maximum of 2 for paired-ended reads.")
        print("Complement (-c or --complement) argument has to be maximum of 2 for paired-ended reads.")
        exit()

    if not Path(args.output[0]).is_dir():
        logging.info("Directory {0} does not exist.".format(args.output[0]))
        print("Directory {0} does not exist.".format(args.output[0]))
        exit()
    
    if not Path(args.dir[0]).is_dir():
        logging.info("Directory {0} does not exist.".format(args.dir[0]))
        print("Directory {0} does not exist.".format(args.dir[0]))
        exit()

    for file in args.samples:
        if Path(args.dir[0] + "/" + file + args.complement[0] + ".fq.gz").is_file() \
            & Path(args.dir[0] + "/" + file + args.complement[1] + ".fq.gz").is_file():
                pass
        else:
            logging.info("File {0} does not exist in {1}.".format(file, str(args.dir[0])))
            print("File {0} does not exist in {1}.".format(file, str(args.dir[0])))
            exit()
    
    print("All files exists. Continuing the analysis.")

    return args

def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Kallisto for every samples given.")
    parser.add_argument("-s", "--samples", nargs="+", required=True, help="<Required> List of samples to iterate over", type=str)
    parser.add_argument("-c", "--complement", nargs="+", required=False, help="<Optional> Complementary for paired-ended", type=str)
    parser.add_argument("-o", "--output", nargs=1, required=True,help="<Required> Directory to output", type=str, metavar="dir/to/folder/output")
    parser.add_argument("-d", "--dir", nargs=1, required=True, help="<Required> Directory to files", type=str, metavar="dir/to/files")
    args = parser.parse_args()

    # Creating and logging info for the current run
    logging.basicConfig(filename=datetime.now().strftime("%d-%m-%Y_%H-%M-%S") + ".log", filemode="a", \
        format='%(asctime)s %(message)s', datefmt='%d/%m/%Y %H:%M:%S', \
        level=logging.DEBUG)
    logging.info("Samples used: {0}".format(args.samples))
    logging.info("Complements: {0}".format(args.complement))
    logging.info("Working directory for samples: {0}".format(args.dir))
    logging.info("Output directory for results: {0}".format(args.output))

    print("Samples used: {0}".format(args.samples))
    print("Complements: {0}".format(args.complement))
    print("Working directory for samples: {0}".format(args.dir))
    print("Output directory for results: {0}".format(args.output))
    
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
    fargs = read_samples(args=arguments)
    run_kallisto(fargs)