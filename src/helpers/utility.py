from contextlib import contextmanager
import argparse
from subprocess import run
from pathlib import Path


def disk_usage(path: str = "/") -> dict:
    import shutil
    total, used, free = [i / 1073741824 for i in shutil.disk_usage(path)]
    return {"total": total, "used": used, "free": free,
            "used_p": round((used / total) * 100, 2), "free_p": round((free / total) * 100, 2)}


@contextmanager
def working_directory(directory):
    import os
    owd = os.getcwd()
    try:
        os.chdir(directory)
        yield directory
    finally:
        os.chdir(owd)


def decide_format(args: argparse.Namespace) -> argparse.Namespace:
    from subprocess import run
    results = {}
    for format in ['.fq.gz', '.fastq.gz', '.fastq', '.fq']:
        value = sum(format in s for s in run(["ls", "input/"], capture_output=True, text=True).stdout.split("\n"))
        results.update({format: value})

    args.format = max(results, key=results.get)

    return args


def build_directory(args: argparse.Namespace, curr_time: str) -> argparse.Namespace:
    run(["mkdir", "-p", "results_" + curr_time + "/1_quality_control"])
    run(["mkdir", "-p", "results_" + curr_time + "/2_trimmed_output"])
    run(["mkdir", "-p", "results_" + curr_time + "/3_kallisto_results"])
    run(["mkdir", "-p", "results_" + curr_time + "/4_picard_qc"])

    d = "results_" + curr_time + "/"

    args.output = d
    return args
