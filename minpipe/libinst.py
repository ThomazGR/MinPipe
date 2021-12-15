from subprocess import run
from pathlib import Path

class CheckLibs():
    def __init__(self, logger) -> None:
        self.logger = logger
        pass

    def __check_kallisto(self) -> None:
        kall = run(["which", "kallisto"])
        if all(contains in kall.stdout for contains in ["kallisto", "/"]) and Path(kall.stdout).is_dir():
            self.logger.info("Contains Kallisto? Yes")
        else:
            self.logger.info("Contains Kallisto? No")
            self.logger.info("EXITING")
            exit()
    
    def __check_trim_galore(self) -> None:
        tg = run(["which", "trim_galore"])
        if all(contains in tg.stdout for contains in ["trim_galore", "/"]) and Path(tg.stdout).is_dir():
            self.logger.info("Contains Trim Galore? Yes")
        else:
            self.logger.info("Contains Trim Galore? No")
            self.logger.info("EXITING")
            exit()
    
    def __check_fastqc(self) -> None:
        fqc = run(["which", "fastqc"])
        if all(contains in fqc.stdout for contains in ["fastqc", "/"]) and Path(fqc.stdout).is_dir():
            self.logger.info("Contains FastQC? Yes")
        else:
            self.logger.info("Contains FastQC? No")
            self.logger.info("EXITING")
            exit()

    def check_all(self) -> None:
        self.__check_fastqc(self)
        self.__check_trim_galore(self)
        self.__check_kallisto(self)