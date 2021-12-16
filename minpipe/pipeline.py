from subprocess import run
from datetime import datetime
import logging
import re

class PipelineCreator():
    def __init__(self, single: bool, complement: list, samples: list, format: str, output: str, index: str, input: str,
                logger = None, threads: str = "4", bootstrap: str = "100", min_len: str = "25", quality: str = "20") -> None:
        self.single = single
        self.complement = complement
        self.samples = samples
        self.format = format
        self.output = output
        self.input = input
        self.index = index
        self.threads = threads
        self.bootstrap = bootstrap
        self.min_len = min_len
        self.quality = quality
        self.curr_time = datetime.now().strftime("%d-%m-%Y_%H-%M-%S")
        #self.format, if format is passed then no decide_format needed
        #self.input where is passed input path to sample files
        #if format in sample names no need to decide_format and get format from sample name
        #if path/in/sample/names no need to get input path
        #if both format and path in sample names then get both from samples and insert to self.input and self.format
        
        find_format_path = re.compile("[a-zA-Z0-9-_]*/[a-zA-Z0-9-]*(\.fa\.gz|\.fq\.gz|\.fastq\.gz|\.fasta\.gz)")

        if all(
            bool(re.search(find_format_path, sample))
            for sample in self.samples
            ):
            self.format = format
            self.input = input



        self.logger = logger
        if logger is None:
            logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%d/%m/%Y %H:%M:%S")
            self.logger = logging.getLogger("main.logger")
            self.logger.addHandler(logging.FileHandler(f"{self.curr_time}.log", "a"))

        pass

    def __run_paired(self) -> None:
        for sample in self.samples:
            qc = run(["fastqc", "-o", f"{self.output}1_quality_control", "--no-extract",
                    f"input/{sample}{self.complement[0]}{self.format}",
                    f"input/{sample}{self.complement[1]}{self.format}"],
                    capture_output=True, text=True)
            self.logger.info(qc.stdout)
            self.logger.info(qc.stderr)

            trim = run(["trim_galore", "--quality", self.quality, "--fastqc", "--length", self.min_len, "--paired",
                        "-o", f"{self.output}2_trimmed_output",
                        f"input/{sample}{self.complement[0]}{self.format}",
                        f"input/{sample}{self.complement[1]}{self.format}"],
                    capture_output=True, text=True)
            self.logger.info(trim.stdout)
            self.logger.info(trim.stderr)

            run(["mkdir", "-p", f"results_{self.curr_time}/3_kallisto_results/{sample}"])

            kall = run(["kallisto", "quant", "-t", self.threads[0], "-b", self.bootstrap[0],
                        "-i", f"index/{self.index[0]}", "-o", f"{self.output}3_kallisto_results/{sample}",
                        "--pseudobam",
                        f"{self.output}2_trimmed_output/{sample}{self.complement[0]}_val_1.fq.gz",
                        f"{self.output}2_trimmed_output/{sample}{self.complement[1]}_val_2.fq.gz"],
                    capture_output=True, text=True)
            self.logger.info(kall.stdout)
            self.logger.info(kall.stderr)
        
        pass

    def __run_single(self) -> None:
        for sample in self.samples:
            qc = run(["fastqc", "-o", f"{self.output}1_quality_control", "--no-extract",
                    f"input/{sample}{self.format}"],
                    capture_output=True, text=True)
            self.logger.info(qc.stdout)
            self.logger.info(qc.stderr)

            trim = run(["trim_galore", "--quality", self.quality, "--fastqc", "--length", self.min_len,
                        "-o", f"{self.output}2_trimmed_output",
                        f"input/{sample}{self.format}"],
                    capture_output=True, text=True)
            self.logger.info(trim.stdout)
            self.logger.info(trim.stderr)

            run(["mkdir", "-p", f"results_{self.curr_time}/3_kallisto_results/{sample}"])

            kall = run(["kallisto", "quant", "-t", self.threads[0], "-b", self.bootstrap[0],
                        "--pseudobam",
                        "--single", "-i", f"index/{self.index[0]}", "-o", f"{self.output}3_kallisto_results/{sample}",
                        f"{self.output}2_trimmed_output/{sample}_trimmed.fq.gz"],
                    capture_output=True, text=True)
            self.logger.info(kall.stdout)
            self.logger.info(kall.stderr)
        
        pass

    def run_full(self) -> None:
        if not self.single:
            self.__run_paired()
        elif self.single:
            self.__run_single()
            
        self.logger.info("Finished pseudoalignment!")

        pass