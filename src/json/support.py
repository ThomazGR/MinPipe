import json
import argparse
from itertools import chain

def json_parse(file: str = None) -> list:
    try:
        with open(f"examples/{file}.json") as fd:
            data = json.load(fd)
    except FileNotFoundError as fnfe:
        raise Exception(fnfe)

    vals = []
    for k in data.keys():
        vals.append(f"--{str(k)}")
        if type(data[k]) == list:
            for i in data[k]:
                vals.append(str(i))
        elif type(data[k]) == bool:
            pass
        else:
            vals.append(str(data[k]))

    return vals