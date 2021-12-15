from contextlib import contextmanager
import argparse
from subprocess import run
from pathlib import Path
from os import makedirs, getcwd, chdir


def disk_usage(path: str = "/") -> dict:
    import shutil
    total, used, free = [i / 1073741824 for i in shutil.disk_usage(path)]
    return {"total": total, "used": used, "free": free,
            "used_p": round((used / total) * 100, 2), "free_p": round((free / total) * 100, 2)}


@contextmanager
def working_directory(directory):
    owd = getcwd()
    chdir(directory)
        


def decide_format(input_folder) -> str:
    results = {}
    for format in ['.fq.gz', '.fastq.gz', '.fastq', '.fq']:
        value = sum(format in s for s in run(["ls", f"{input_folder}/"], capture_output=True, text=True).stdout.split("\n"))
        results.update({format: value})

    return str(max(results, key=results.get))


def build_directory(curr_time: str) -> None:
    makedirs(f"results_{curr_time}/1_quality_control")
    makedirs(f"results_{curr_time}/2_trimmed_output")
    makedirs(f"results_{curr_time}/3_kallisto_results")
    makedirs(f"results_{curr_time}/4_picard_qc")

    pass