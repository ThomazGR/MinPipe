import json

def json_parse(file: str = None) -> list:
    try:
        with open(f"input/{file}") as fd:
            data = json.load(fd)
    except FileNotFoundError as fnfe:
        raise Exception(fnfe)

    try:
        vals = []
        for k in data.keys():
            if type(data[k]) == bool and data[k] is True:
                vals.append(f"--{str(k)}")
            elif type(data[k]) == bool and data[k] is False:
                pass
            elif type(data[k]) == list:
                vals.append(f"--{str(k)}")
                for i in data[k]:
                    vals.append(str(i))
            else:
                vals.append(f"--{str(k)}")
                vals.append(str(data[k]))
    except Exception as exc:
        raise Exception(exc)

    return vals