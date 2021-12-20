from subprocess import run
import logging

#TODO: variant calling for each sample
class ExtensiveQC:
    def __init__(self, samples: list, output: str, logger: logging.Logger = None) -> None:
        self.samples = samples
        self.output = output
        self.logger = logger
        pass

    def QualityScoreDist(self):
        for sample in self.samples:
            qsd = run(["picard", "QualityScoreDistribution",
                   "-I", f"{self.output}3_kallisto_results/{sample}/pseudoalignments.bam",
                   "-O", f"{self.output}4_picard_qc/{sample}.txt",
                   "-CHART", f"{self.output}4_picard_qc/{sample}.pdf"],
                  capture_output=True, text=True)
            self.logger.info(qsd.stderr)
            self.logger.info(qsd.stdout)

    def run_all(self):
        self.QualityScoreDist()