from contextlib import contextmanager
import argparse

def disk_usage(path: str = "/") -> dict:
    import shutil
    total, used, free = [i/1073741824 for i in shutil.disk_usage(path)]
    return {"total": total, "used":used, "free":free, 
    "used_p":round((used / total) * 100, 2), "free_p":round((free / total) * 100, 2)}

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
        results.update({format:value})

    args.format = max(results, key=results.get)

    return args