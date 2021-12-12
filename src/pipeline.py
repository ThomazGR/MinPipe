from subprocess import run
from pathlib import Path
import logging
from utility import decide_format as decide, check_index, build_directory
from datetime import datetime

class TestSamples():
    def __init__(self, logger, single, complement, samples, format) -> None:
        self.samples = samples
        self.complement = complement
        self.single = single
        self.logger = logger
        pass

    def read_samples(self) -> None:
        if self.single and self.complement is not None:
            self.logger.info("Single-ended analysis does not contain complements. \
                Complements are for paired-ended (e.g. sample1_R1.fastq.gz and sample1_R2.fastq.gz)")
            exit()
        elif not self.single and self.complement is not None:
            if len(self.complement) > 2:
                self.logger.info("Complement (-c or --complement) argument has to be maximum \
                    of 2 for paired-ended reads.")
                exit()

        if not self.single:
            for file in self.samples:
                if Path("input/" + file + self.complement[0] + self.format).is_file() \
                        and Path("input/" + file + self.complement[1] + self.format).is_file():
                    pass
                else:
                    self.logger.info(f"File {file} does not exist in input/ folder.")
                    exit()
        elif self.single:
            for file in self.samples:
                if Path("input/" + file + self.format).is_file():
                    pass
                else:
                    self.logger.info(f"File {file} does not exist in input/ folder.")
                    exit()

        self.logger.info("All files exists. Continuing the analysis.")

        pass

    def __exit__(self, type, value, traceback):
        try:
            run(args=[], text=True, capture_output=True)
        except Exception as exc:
            print(exc)

class TestIndexTranscript():
    def __init__(self, logger, transcript, index) -> None:
        self.logger = logger
        self.transcript = transcript
        self.index = index
        pass

    def __download_hsa_transcript(self) -> None:
        try:
            run([
                "wget", "-P", "index/", "-c",
                "http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz",
                "-O", "index/homo_sapiens_GRCh38_cdna.fa.gz"
            ])
        except Exception as exc:
            self.logger.info(exc)
            quit()

        self.logger.info("Homo sapiens transcript has been downloaded!")
        self.transcript = ["homo_sapiens_GRCh38_cdna.fa.gz"]
        
        pass


    def __download_mmu_trscript(self) -> None:
        try:
            run([
                "wget", "-P", "index/", "-c",
                "http://ftp.ensembl.org/pub/release-104/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz",
                "-O", "index/mus_musculus_GRCm39_cdna.fa.gz"
            ])
        except Exception as exc:
            self.logger.info(exc)
            quit()

        self.logger.info("Mus musculus transcript has been downloaded!")
        self.transcript = ["mus_musculus_GRCm39_cdna.fa.gz"]
        
        pass

    def create_index(self) -> None:
        if "/" in self.transcript[0]:
            idx_name = self.transcript[0].split("/")[-1].split(".")[0]
        else:
            idx_name = self.transcript[0].split(".")[0]

        idx = run(["kallisto", "index", "-i", "index/" + idx_name + ".idx",
                "index/" + self.transcript[0]], capture_output=True, text=True)
        self.logger.info(idx.stdout)
        self.logger.info(idx.stderr)

        self.index = ["index/" + idx_name + ".idx"]

        pass

    def check_idx_trans(self) -> None:
        if self.index is None and self.transcript is None:
            self.logger.info("No index or transcript has been passed")
            exit()
        elif self.index is None and self.transcript:
            if any(fmt in self.transcript[0] for fmt in ['.fa', '.fa.gz',
                                                        '.fastq', '.fastq.gz',
                                                        '.fq', '.fq.gz']):
                try:
                    self.create_index(self)
                    self.logger.info("Index created")
                    check_index(self.index)
                except Exception as ex:
                    self.logger.info(ex)
                    exit()
            else:
                if self.transcript[0].lower() == "mmu":
                    self.__download_mmu_trscript(self)
                    self.create_index(self)
                    self.logger.info("Mmu transcript downloaded and index created.")
                    check_index(self.index)
                elif self.transcript[0].lower() == "hsa":
                    self.__download_hsa_transcript(self)
                    self.create_index(self)
                    self.logger.info("Hsa transcript downloaded and index created.")
                    check_index(self.index)
                else:
                    self.logger.info("Species or format not supported. \
                        Select hsa or mmu, or download your own transcript.")
                    exit()
        elif self.index and self.transcript is None:
            try:
                check_index("index/" + self.index[0])
            except Exception as ex:
                self.logger.info(ex)
                exit()
        else:
            self.logger.info("You can only pass `--index` or `--transcript` argument. \
                If both are passed the pipeline don't know if needs to build another index with \
                    the transcript passed or use the index without building a new one.")
            exit()

        pass

class PipelineCreator():
    def __init__(self, single: bool, complement: list, samples: list, format: str, 
                output: str, index: str, threads: str, bootstrap: str) -> None:
        self.single = single
        self.complement = complement
        self.samples = samples
        self.format = format
        self.output = output
        self.index = index
        self.threads = threads
        self.bootstrap = bootstrap
        self.curr_time = datetime.now().strftime("%d-%m-%Y_%H-%M-%S")

        pass

    def __run_paired(self) -> None:
        for sample in self.samples:
            qc = run(["fastqc", "-o", self.output + "1_quality_control", "--no-extract",
                    "input/" + sample + self.complement[1] + self.format],
                    "input/" + sample + self.complement[0] + self.format,
                    capture_output=True, text=True)
            self.logger.info(qc.stdout)
            self.logger.info(qc.stderr)

            trim = run(["trim_galore", "--quality", "20", "--fastqc", "--length", "25", "--paired",
                        "-o", self.output + "2_trimmed_output",
                        "input/" + sample + self.complement[0] + self.format,
                        "input/" + sample + self.complement[1] + self.format],
                    capture_output=True, text=True)
            self.logger.info(trim.stdout)
            self.logger.info(trim.stderr)

            run(["mkdir", "-p", f"results_{self.curr_time}/3_kallisto_results/{sample}"])

            kall = run(["kallisto", "quant", "-t", self.threads[0], "-b", self.bootstrap[0],
                        "-i", "index/" + self.index[0], "-o", self.output + "3_kallisto_results/" + sample,
                        "--pseudobam",
                        self.output + "2_trimmed_output/" + sample + self.complement[0] + "_val_1.fq.gz",
                        self.output + "2_trimmed_output/" + sample + self.complement[1] + "_val_2.fq.gz"],
                    capture_output=True, text=True)
            self.logger.info(kall.stdout)
            self.logger.info(kall.stderr)
        
        pass

    def __run_single(self) -> None:
        for sample in self.samples:
            qc = run(["fastqc", "-o", self.output + "1_quality_control", "--no-extract",
                    "input/" + sample + self.format],
                    capture_output=True, text=True)
            self.logger.info(qc.stdout)
            self.logger.info(qc.stderr)

            trim = run(["trim_galore", "--quality", "20", "--fastqc", "--length", "25",
                        "-o", self.output + "2_trimmed_output",
                        "input/" + sample + self.format],
                    capture_output=True, text=True)
            self.logger.info(trim.stdout)
            self.logger.info(trim.stderr)

            run(["mkdir", "-p", f"results_{self.curr_time}/3_kallisto_results/{sample}"])

            kall = run(["kallisto", "quant", "-t", self.threads[0], "-b", self.bootstrap[0],
                        "--pseudobam",
                        "--single", "-i", "index/" + self.index[0], "-o", self.output + "3_kallisto_results/" + sample,
                        self.output + "2_trimmed_output/" + sample + "_trimmed.fq.gz"],
                    capture_output=True, text=True)
            self.logger.info(kall.stdout)
            self.logger.info(kall.stderr)
        
        pass

    def run_full(self) -> None:
        if not self.single:
            self.__run_paired(self)
        elif self.single:
            self.__run_single(self)
            
        self.logger.info("Finished pseudoalignment!")

        pass