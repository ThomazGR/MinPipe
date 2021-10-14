import argparse, logging
from os import cpu_count
from subprocess import run
from pathlib import Path
from datetime import datetime
from src.utility import working_directory as wd, disk_usage as du, decide_format as decide

CURR_TIME = datetime.now().strftime("%d-%m-%Y_%H-%M-%S")

def picard_qc(args: argparse.Namespace):
    logger = logging.getLogger("main.logger")

    for sample in args.samples:
        qsd = run(["picard", "QualityScoreDistribution", 
        "-I", args.output + "3_kallisto_results/" + sample + "/pseudoalignments.bam",
        "-O", args.output + "4_picard_qc/" + sample + ".txt",
        "-CHART", args.output + "4_picard_qc/"  + sample + ".pdf"],
        capture_output=True, text=True)
        if qsd.stdout:
            logger.info(qsd.stdout)
        if qsd.stderr:
            logger.info(qsd.stderr)

        
    return

def check_index(path_index: str):
    if Path(path_index).is_file():
        pass
    else:
        exit(f"No index file found on {path_index}")

    return 

def build_directory(args: argparse.Namespace) -> argparse.Namespace:
    run(["mkdir", "-p", "results_" + CURR_TIME + "/1_quality_control"])
    run(["mkdir", "-p", "results_" + CURR_TIME + "/2_trimmed_output"])
    run(["mkdir", "-p", "results_" + CURR_TIME + "/3_kallisto_results"])
    run(["mkdir", "-p", "results_" + CURR_TIME + "/4_picard_qc"])

    d = "/results_" + CURR_TIME + "/"

    args.output = d
    return args

def create_index(args: argparse.Namespace) -> argparse.Namespace:
    logger = logging.getLogger("main.logger")
    if "/" in args.transcript[0]:
        idx_name = args.transcript[0].split("/")[-1].split(".")[0]
    else:
        idx_name = args.transcript[0].split(".")[0]
        
    idx = run(["kallisto", "index", "index/" + idx_name + ".idx",
        "index/" + args.transcript[0]], capture_output=True, text=True)

    if idx.stdout:
        logger.info(idx.stdout)
    if idx.stderr:
        logger.info(idx.stderr)
    
    args.index = "index/" + idx_name + ".idx"

    return args

def run_qctk(args: argparse.Namespace):
    logger = logging.getLogger("main.logger")
    if not args.single:
        for sample in args.samples:
            qc = run(["fastqc", "-o", args.output + "1_quality_control", "--no-extract",
                "input/" + sample + args.complement[0] + args.format,
                "input/" + sample + args.complement[1] + args.format],
                    capture_output=True, text=True)
            if qc.stdout:
                logger.info(qc.stdout)
            if qc.stderr:
                logger.info(qc.stderr)

            trim = run(["trim_galore", "--quality", "20", "--fastqc", "--length", "25", "--paired",
                "-o", args.output + "2_trimmed_output",
                    "input/" + sample + args.complement[0] + args.format,
                    "input/" + sample + args.complement[1] + args.format],
                        capture_output=True, text=True)
            if trim.stdout:
                logger.info(trim.stdout)
            if trim.stderr:
                logger.info(trim.stderr)

            kall = run(["kallisto", "quant", "-t", args.threads[0], "-b", args.bootstrap[0],
                "-i", args.index, "-o", args.output + "3_kallisto_output/" + sample,
                "--pseudobam",
                args.output + "2_trimmed_output/" + sample + args.complement[0] + "_val_1" + args.format,
                args.output + "2_trimmed_output/" + sample + args.complement[1] + "_val_2" + args.format],
                capture_output=True, text=True)
            if kall.stdout: 
                logger.info(kall.stdout)
            if kall.stderr: 
                logger.info(kall.stderr)

    elif args.single:
        for sample in args.samples:
            qc = run(["fastqc", "-o", args.output + "1_quality_control", "--no-extract",
                "input/" + sample + args.format],
                    capture_output=True, text=True)
            if qc.stdout:
                logger.info(qc.stdout)
            if qc.stderr:
                logger.info(qc.stderr)

            trim = run(["trim_galore", "--quality", "20", "--fastqc", "--length", "25",
                "-o", args.output + "2_trimmed_output",
                    "input/" + sample + args.format],
                        capture_output=True, text=True)
            if trim.stdout:
                logger.info(trim.stdout)
            if trim.stderr:
                logger.info(trim.stderr)

            kall = run(["kallisto", "quant", "-t", args.threads[0], "-b", args.bootstrap[0],
                "--pseudobam",
                "--single", "-i", args.index, "-o", args.output + "3_kallisto_output/" + sample,
                args.output + "2_trimmed_output/" + sample + "_trimmed" + args.format],
                capture_output=True, text=True)
            if kall.stdout: 
                logger.info(kall.stdout)
            if kall.stderr: 
                logger.info(kall.stderr)

    return print("Finished pseudoalignment!")

def read_samples(args: argparse.Namespace) -> dict:
    logger = logging.getLogger("main.logger")
    if args.single and args.complement is not None:
        logger.info("Single-ended analysis does not contain complements. \
            Complements are for paired-ended (e.g. sample1_R1.fastq.gz and sample1_R2.fastq.gz)")
        exit()
    elif not args.single and args.complement is not None:
        if len(args.complement) > 2:
            logger.info("Complement (-c or --complement) argument has to be maximum \
                of 2 for paired-ended reads.")
            exit()

    if not args.single:
        for file in args.samples:
            if Path("input/" + file + args.complement[0] + args.format).is_file() \
                and Path("input/" + file + args.complement[1] + args.format).is_file():
                    pass
            else:
                logger.info(f"File {file} does not exist in input/ folder.")
                exit()
    elif args.single:
        for file in args.samples:
            if Path("input/" + file + args.format).is_file():
                pass
            else:
                logger.info(f"File {file} does not exist in input/ folder.")
                exit()
    
    print("All files exists. Continuing the analysis.")

    return args

def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Kallisto for every samples given.")
    parser.add_argument("-s", "--samples", nargs="+", required=True, 
        help="<Required> List of samples to iterate over", type=str)
    parser.add_argument("-c", "--complement", nargs="+", required=False, 
        help="<Optional> Complementary for paired-ended", type=str)
    parser.add_argument("-i", "--index", nargs=1, required=False,
        help="<Optional> Name of the index file to be used in pseudoalignment. If no file has been passed \
            it has to have a UNIQUE .idx file in `index` folder or it will raise an error.")
    parser.add_argument("-t", "--transcript", nargs=1, required=False, 
        help="<Optional> Name of the transcript file to be indexed.")
    parser.add_argument("--threads", nargs=1, required=False, default=["1"],
        help="<Optional> Number of threads to be used in quantification for Kallisto. Default: 1.")
    parser.add_argument("-b", "--bootstrap", nargs=1, required=False, default=["100"],
        help="<Optional> Number of bootstrap samples. Default: 100.")
    parser.add_argument("--single", action='store_true', required=False,
        help="<Optional> Flag to indicate single-ended quantification without complements.")

    args = parser.parse_args()

    if args.index is None and args.transcript is None:
        lst = run(["ls", "index/"], capture_output=True, text=True).stdout.split("\n")
        n = [k for k in lst if ".idx" in k]
        if n == 0:
            print("A index file `.idx` needs to be inserted in the `index` folder or a transcript file name \
                needs to be passed to be consumed from the index folder.")
            exit()
        elif n > 1:
            print("No specific `.idx` file has been passed in the -i/--index parameters but many `.idx` \
                files has been found on the `index` folder. Please specify the file.")
            exit()
        elif n == 1:
            args.index = "index/" + n[0]
    elif args.index is None and args.transcript:
        try:
            args = create_index(args.transcript)
            print("Index created")
            check_index(args.index)
        except Exception as ex:
            print(ex)
            exit()

    # Creating and logging info for the current run
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%d/%m/%Y %H:%M:%S")
    logger = logging.getLogger("main.logger")
    logger.addHandler(logging.FileHandler(CURR_TIME + ".log", "a"))
    logger.info("Samples used: {0}".format(args.samples))
    logger.info("Complements: {0}".format(args.complement))
    logger.info("Working directory for samples: {0}".format(args.dir))
    
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
    arguments = build_directory(arguments)
    arguments = decide(arguments)
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