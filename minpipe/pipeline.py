from subprocess import run
from datetime import datetime
from pathlib import Path
from os import makedirs
import logging
import warnings

from check import TestIndexTranscript, TestSamples
from libinst import CheckLibs
from quality import ExtensiveQC
# import re


class PipelineCreator:
    def __init__(self, samples: list, single: bool = False, complement: list = None, file_format: str = None,
                 output_path: str = None, index: str = None, input_path: str = "input/", logger: logging.Logger = None,
                 transcript: str = None, threads: int = 4, bootstrap: int = 100, min_len: int = 25, quality: int = 20,
                 ext_qc: bool = False) -> None:
        """
        Construct the PipelineCreator object to run full pipeline writing results to parameter/default folder.

        :type single: bool
        :type complement: list
        :type samples: list
        :type file_format: str
        :type index: str
        :type transcript: str
        :type logger: logging.Logger
        :type threads: int
        :type bootstrap: int
        :type min_len: int
        :type quality: int
        :type ext_qc: bool
        :type output_path: str
        :type input_path: str
        """
        self.single = single
        self.complement = complement
        self.samples = samples
        self.format = file_format
        self.output = output_path
        self.input = input_path
        self.index = index
        self.transcript = transcript
        self.threads = str(threads)
        self.bootstrap = str(bootstrap)
        self.min_len = str(min_len)
        self.quality = str(quality)
        self.ext_qc = ext_qc
        self.logger = logger
        self.curr_time = str(datetime.now().strftime("%d-%m-%Y_%H-%M-%S"))
        # self.format, if format is passed then no decide_format needed
        # self.input where is passed input path to sample files
        # if format in sample names no need to decide_format and get format from sample name
        # if path/in/sample/names no need to get input path
        # if both format and path in sample names then get both from samples and insert to self.input and self.format

        pass

    def __enter__(self):
        # find_format_path = re.compile(r"[a-zA-Z0-9-_]*/[a-zA-Z0-9-]*(\.fa\.gz|\.fq\.gz|\.fastq\.gz|\.fasta\.gz)")
        if self.format is None:
            self.__decide_format()

        if self.input[-1] != "/":
            self.input += "/"
            assert Path(self.input).is_dir() is True, f"Input path should be a valid path. Passed `{self.input}`"
        else:
            assert Path(self.input).is_dir() is True, f"Input path should be a valid path. Passed `{self.input}`"

        if self.output is None:
            self.output = f"results_{self.curr_time}/"
        else:
            if self.output[-1] == "/":
                assert Path(self.output).is_dir() is True, f"Output path should be a valid path. Passed `{self.output}`"
            else:
                self.output += "/"
                assert Path(self.output).is_dir() is True, f"Output path should be a valid path. Passed `{self.output}`"

        # # # # # # # # # # # # # # # # # #
        # Argument checking
        # # # # # # # # # # # # # # # # # #
        print("# # # # # # # # # # # # # # #")
        print("# Current parameters")
        print("# # # # # # # # # # # # # # #")
        print(f"Single-ended analysis? {self.single}")
        print(f"Complements: {self.complement}")
        print(f"Samples: {self.samples}")
        print(f"File format: {self.format}")
        print(f"Output path: {self.output}")
        print(f"Input path: {self.input}")
        print(f"Index used: {self.index}")
        print(f"Threads used: {self.threads}")
        print(f"Quantification bootstrap: {self.bootstrap}")
        print(f"Minimum length of trimmage: {self.min_len}")
        print(f"Minimum quality of trimmage: {self.quality}")
        print(f"Logging object: {bool(self.logger)}")
        print(f"Time of start: {self.curr_time}")

        response = False
        while response not in ['y', 'n']:
            response = str(input("\nIs this correct? [y/n]\n")).lower()
            if response not in ['y', 'n']:
                print("\nType `y` for Yes or `n` for No")
            else:
                if response == 'y':
                    pass
                elif response == 'n':
                    print("Arguments not right. Exiting.")
                    exit()

        self.__build_directory()

        if self.logger is None:
            logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%d/%m/%Y %H:%M:%S")
            self.logger = logging.getLogger("main.logger")
            self.logger.addHandler(logging.FileHandler(f"{self.curr_time}.log", "a"))

        self.logger.info(f"Samples used: {self.samples}")
        self.logger.info(f"Complements: {self.complement}")
        self.logger.info(f"Index: {self.index}")
        self.logger.info(f"Transcript: {self.transcript}")
        self.logger.info(f"Threads number: {self.threads}")
        self.logger.info(f"Bootstrap number: {self.bootstrap}")
        self.logger.info(f"Single ended: {self.single}")
        self.logger.info(f"Extensive Quality Control: {self.ext_qc}")
        self.logger.info(f"Minimum quality for trimmage: {self.quality}")
        self.logger.info(f"Minimum length for trimmage: {self.min_len}")
        self.logger.info(f"Input path: {self.input}")
        self.logger.info(f"Output path: {self.output}")

        lib_is_installed = CheckLibs(self.logger)
        lib_is_installed.check_all()
        
        test_index_transc = TestIndexTranscript(self.logger, transcript=self.transcript, index=self.index)
        index = test_index_transc.check_idx_trans()
        if ".idx" in index: self.index = index

        test_samples = TestSamples(self.logger, single=self.single, complement=self.complement, 
                                    samples=self.samples, file_format=self.format)
        test_samples.read_samples()

        return self

    def __start_log(self) -> None:
        if self.logger is None:
            logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%d/%m/%Y %H:%M:%S")
            self.logger = logging.getLogger("main.logger")
            self.logger.addHandler(logging.FileHandler(f"{self.curr_time}.log", "a"))
            
            self.logger.info(f"Samples used: {self.samples}")
            self.logger.info(f"Complements: {self.complement}")
            self.logger.info(f"Index: {self.index}")
            self.logger.info(f"Transcript: {self.transcript}")
            self.logger.info(f"Threads number: {self.threads}")
            self.logger.info(f"Bootstrap number: {self.bootstrap}")
            self.logger.info(f"Single ended: {self.single}")
            self.logger.info(f"Extensive Quality Control: {self.ext_qc}")
            self.logger.info(f"Minimum quality for trimmage: {self.quality}")
            self.logger.info(f"Minimum length for trimmage: {self.min_len}")
            self.logger.info(f"Input path: {self.input}")
            self.logger.info(f"Output path: {self.output}")
        
        pass

    def __run_paired(self) -> None:
        """
        Run paired-ended analysis using current tools
        :return: Writes quality control, trimmed plus quality control and kallisto abundance/BAM results
        """
        for sample in self.samples:
            qc = run(["fastqc", "-o", f"{self.output}1_quality_control", "--no-extract",
                      f"{self.input}{sample}{self.complement[0]}{self.format}",
                      f"{self.input}{sample}{self.complement[1]}{self.format}"],
                     capture_output=True, text=True)
            self.logger.info(qc.stdout)
            self.logger.info(qc.stderr)

            trim = run(["trim_galore", "--quality", self.quality, "--fastqc", "--length", self.min_len, "--paired",
                        "-o", f"{self.output}2_trimmed_output",
                        f"{self.input}{sample}{self.complement[0]}{self.format}",
                        f"{self.input}{sample}{self.complement[1]}{self.format}"],
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
        """
        Run single-ended analysis using current tools
        :return: Writes quality control, trimmed plus quality control and kallisto abundance/BAM results
        """
        for sample in self.samples:
            qc = run(["fastqc", "-o", f"{self.output}1_quality_control", "--no-extract",
                      f"{self.input}{sample}{self.format}"],
                     capture_output=True, text=True)
            self.logger.info(qc.stdout)
            self.logger.info(qc.stderr)

            trim = run(["trim_galore", "--quality", self.quality, "--fastqc", "--length", self.min_len,
                        "-o", f"{self.output}2_trimmed_output",
                        f"{self.input}{sample}{self.format}"],
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

    def __decide_format(self) -> None:
        results = {}
        for file_format in ['.fq.gz', '.fastq.gz', '.fastq', '.fq']:
            value = sum(file_format in file for file in run(["ls", f"{self.input}/"],
                                                            capture_output=True, text=True).stdout.split("\n"))
            results.update({format: value})

        self.format = str(max(results, key=results.get))

        if max(results.values()) <= len(self.samples):
            warnings.warn("Samples may not have the same file extension or it has not been detected.",
                          category=UserWarning)

        pass

    def __build_directory(self) -> None:
        makedirs(f"{self.output}1_quality_control")
        makedirs(f"{self.output}2_trimmed_output")
        makedirs(f"{self.output}3_kallisto_results")
        makedirs(f"{self.output}4_picard_qc")

        pass

    def run_full(self) -> None:
        """
        Run full pipeline choosing between single or paired-ended type
        :return: None
        """
        self.__start_log()

        if not self.single:
            self.__run_paired()
        elif self.single:
            self.__run_single()

        if self.ext_qc:
            qc = ExtensiveQC(samples=self.samples, output=self.output, logger=self.logger)
            qc.run_all()

        self.logger.info("Finished pseudoalignment!")

        pass
