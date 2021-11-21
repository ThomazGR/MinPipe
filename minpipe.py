import argparse
import logging
from subprocess import run
from pathlib import Path
from datetime import datetime
from src.json.support import json_parse
from src.yaml.support import yaml_parse
from src.utility import (
    decide_format as decide,
    build_directory,
    check_index
)

CURR_TIME = datetime.now().strftime("%d-%m-%Y_%H-%M-%S")


def download_hsa_transcript(args: argparse.Namespace) -> argparse.Namespace:
    logger = logging.getLogger("main.logger")

    try:
        run([
            "wget", "-P", "index/", "-c",
            "http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz",
            "-O", "index/homo_sapiens_GRCh38_cdna.fa.gz"
        ])
    except Exception as exc:
        logger.info(exc)
        quit()

    logger.info("Homo sapiens transcript has been downloaded!")
    args.transcript = ["homo_sapiens_GRCh38_cdna.fa.gz"]
    return args


def download_mmu_trscript(args: argparse.Namespace) -> argparse.Namespace:
    logger = logging.getLogger("main.logger")

    try:
        run([
            "wget", "-P", "index/", "-c",
            "http://ftp.ensembl.org/pub/release-104/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz",
            "-O", "index/mus_musculus_GRCm39_cdna.fa.gz"
        ])
    except Exception as exc:
        logger.info(exc)
        quit()

    logger.info("Mus musculus transcript has been downloaded!")
    args.transcript = ["mus_musculus_GRCm39_cdna.fa.gz"]
    return args


def picard_qc(args: argparse.Namespace):
    logger = logging.getLogger("main.logger")

    for sample in args.samples:
        qsd = run(["picard", "QualityScoreDistribution",
                   "-I", args.output + "3_kallisto_results/" + sample + "/pseudoalignments.bam",
                   "-O", args.output + "4_picard_qc/" + sample + ".txt",
                   "-CHART", args.output + "4_picard_qc/" + sample + ".pdf"],
                  capture_output=True, text=True)
        if qsd.stdout:
            logger.info(qsd.stdout)
        if qsd.stderr:
            logger.info(qsd.stderr)

    logger.info("Finished Quality Score Distribution for all files!")

    return True


def create_index(args: argparse.Namespace) -> argparse.Namespace:
    logger = logging.getLogger("main.logger")
    if "/" in args.transcript[0]:
        idx_name = args.transcript[0].split("/")[-1].split(".")[0]
    else:
        idx_name = args.transcript[0].split(".")[0]

    idx = run(["kallisto", "index", "-i", "index/" + idx_name + ".idx",
               "index/" + args.transcript[0]], capture_output=True, text=True)

    if idx.stdout:
        logger.info(idx.stdout)
    if idx.stderr:
        logger.info(idx.stderr)

    args.index = ["index/" + idx_name + ".idx"]

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

            run(["mkdir", "-p", f"results_{CURR_TIME}/3_kallisto_results/{sample}"])

            kall = run(["kallisto", "quant", "-t", args.threads[0], "-b", args.bootstrap[0],
                        "-i", "index/" + args.index[0], "-o", args.output + "3_kallisto_results/" + sample,
                        "--pseudobam",
                        args.output + "2_trimmed_output/" + sample + args.complement[0] + "_val_1.fq.gz",
                        args.output + "2_trimmed_output/" + sample + args.complement[1] + "_val_2.fq.gz"],
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

            run(["mkdir", "-p", f"results_{CURR_TIME}/3_kallisto_results/{sample}"])

            kall = run(["kallisto", "quant", "-t", args.threads[0], "-b", args.bootstrap[0],
                        "--pseudobam",
                        "--single", "-i", "index/" + args.index[0], "-o", args.output + "3_kallisto_results/" + sample,
                        args.output + "2_trimmed_output/" + sample + "_trimmed.fq.gz"],
                       capture_output=True, text=True)
            if kall.stdout:
                logger.info(kall.stdout)
            if kall.stderr:
                logger.info(kall.stderr)

    logger.info("Finished pseudoalignment!")

    return True


def read_samples(args: argparse.Namespace) -> argparse.Namespace:
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

    logger.info("All files exists. Continuing the analysis.")

    return args


def check_idx_trans(args: argparse.Namespace) -> argparse.Namespace:
    logger = logging.getLogger("main.logger")
    if args.index is None and args.transcript is None:
        logger.info("No index or transcript has been passed")
        exit()
    elif args.index is None and args.transcript:
        if any(fmt in args.transcript[0] for fmt in ['.fa', '.fa.gz',
                                                     '.fastq', '.fastq.gz',
                                                     '.fq', '.fq.gz']):
            try:
                args = create_index(args)
                logger.info("Index created")
                check_index(args.index)
            except Exception as ex:
                logger.info(ex)
                exit()
        else:
            if args.transcript[0].lower() == "mmu":
                args = download_mmu_trscript(args)
                args = create_index(args)
                logger.info("Mmu transcript downloaded and index created.")
                check_index(args.index)
            elif args.transcript[0].lower() == "hsa":
                args = download_hsa_transcript(args)
                args = create_index(args)
                logger.info("Hsa transcript downloaded and index created.")
                check_index(args.index)
            else:
                logger.info("Species or format not supported. \
                    Select hsa or mmu, or download your own transcript.")
                exit()
    elif args.index and args.transcript is None:
        try:
            check_index("index/" + args.index[0])
        except Exception as ex:
            logger.info(ex)
            exit()
    else:
        logger.info("You can only pass `--index` or `--transcript` argument. \
            If both are passed the pipeline don't know if needs to build another index with \
                the transcript passed or use the index without building a new one.")
        exit()

    return args


def arguments() -> argparse.Namespace:
    # # # # # # # # # # # # # # # # # #
    # Building argparse
    # # # # # # # # # # # # # # # # # #
    parser = argparse.ArgumentParser(description="Run Kallisto for every samples given.")
    parser.add_argument("-s", "--samples", nargs="+", required=False,
                        help="<Optional> List of samples to iterate over. All arguments can be passed as Json.",
                        type=str)
    parser.add_argument("-c", "--complement", nargs="+", required=False,
                        help="<Optional> Complementary for paired-ended", type=str)
    parser.add_argument("-i", "--index", nargs=1, required=False,
                        help="<Optional> Name of the index file to be used in pseudoalignment. Either `index` or \
            `transcript`has to be passed.")
    parser.add_argument("-t", "--transcript", nargs=1, required=False,
                        help="<Optional> Name of the transcript file to be indexed. `mmu` or `hsa` can be passed so \
            the transcript will be downloaded automatically and index will be built.")
    parser.add_argument("--threads", nargs=1, required=False, default=["1"],
                        help="<Optional> Number of threads to be used in quantification for Kallisto. Default: 1.")
    parser.add_argument("-b", "--bootstrap", nargs=1, required=False, default=["100"],
                        help="<Optional> Number of bootstrap samples. Default: 100.")
    parser.add_argument("--single", action='store_true', required=False,
                        help="<Optional> Flag to indicate single-ended quantification without complements.")
    parser.add_argument("--ext-qc", action="store_true", required=False,
                        help="<Optional> Flag to indicate that will have extensive QC. **MAY NEED MORE FILES**")
    parser.add_argument("--json", nargs=1, required=False,
                        help="<Optional> Use a Json file to pass all the arguments instead of command line interface. \
            See the parameters.json inside the examples folder to understand how to pass args.")
    parser.add_argument("--yaml", nargs=1, required=False,
                        help="<Optional> Use a YAML/YML file to pass all the arguments instead of \
                        command line interface. \
                        See the parameters.yml inside the examples folder to understand how to pass args.")

    args = parser.parse_args()

    # Creating and logging info for the current run
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%d/%m/%Y %H:%M:%S")
    logger = logging.getLogger("main.logger")
    logger.addHandler(logging.FileHandler(CURR_TIME + ".log", "a"))

    # # # # # # # # # # # # # # # # # #
    # Handling YAML/JSON/CLI parsing
    # # # # # # # # # # # # # # # # # #
    if args.json is not None and args.yaml is not None:
        logger.info("Both YAML and Json detected. Please pass only 1 file through argument handling.")
        exit()
    elif args.json is not None:
        try:
            Path(f"input/{args.json[0]}").is_file()
        except FileNotFoundError as fnfe:
            logger.info(fnfe)
            quit()
        logger.info("Json file found. Proceeding with Json parser for arguments.")
        lst_args = json_parse(args.json[0])
        args = parser.parse_args(lst_args)
    elif args.yaml is not None:
        try:
            Path(f"input/{args.yaml[0]}").is_file()
        except FileNotFoundError as fnfe:
            logger.info(fnfe)
            quit()
        logger.info("YAML file found. Proceeding with YAML parser for arguments.")
        lst_args = yaml_parse(args.yaml[0])
        args = parser.parse_args(lst_args)
    else:
        pass

    logger.info(f"Samples used: {args.samples}")
    logger.info(f"Complements: {args.complement}")
    logger.info(f"Index: {args.index}")
    logger.info(f"Transcript: {args.transcript}")
    logger.info(f"Threads number: {args.threads[0]}")
    logger.info(f"Bootstrap number: {args.bootstrap[0]}")
    logger.info(f"Single ended: {args.single}")
    logger.info(f"Extensive Quality Control: {args.ext_qc}")

    # # # # # # # # # # # # # # # # # #
    # Argument checking
    # # # # # # # # # # # # # # # # # #
    r = False
    while r not in ['y', 'n']:
        r = str(input("\nIs this correct? [y/n]\n")).lower()
        if r in ['y', 'n']:
            if r == 'y':
                pass
            elif r == 'n':
                logger.info("Arguments not right. Exiting.")
                exit()
        else:
            print("\nType `y` for Yes or `n` for No")

    logger.info("Every argument has been checked. Continuing the analysis.\n")

    return args


if __name__ == '__main__':
    arguments = arguments()
    arguments = check_idx_trans(arguments)
    arguments = build_directory(arguments, CURR_TIME)
    arguments = decide(arguments)
    fargs = read_samples(args=arguments)
    run_qctk(fargs)

    if fargs.ext_qc:
        picard_qc(fargs)
    else:
        pass
