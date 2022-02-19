from itertools import chain


def yaml_parse(yml_file: str = None):
    with open(f"input/{yml_file}") as yml:
        comp = list(chain.from_iterable([i.split("\n") for i in yml.readlines()]))

    comp = list(chain.from_iterable([i.split(":") for i in comp]))
    med = [st.strip() if st != "- params" else '' for st in comp]
    fltrd = list(filter(None, med))
    fnl = [st.replace("- ", "") if "- " in st else st for st in fltrd]

    tbp = []
    for index, value in enumerate(fnl):
        if value in ["samples", "complement", "index", "transcript", \
                     "threads", "bootstrap", "single", "ext-qc"]:
            fnl[index] = f"--{value}"

        if value == "true":
            tbp.append(index)
        elif value == "false":
            tbp.append(index - 1)
            tbp.append(index)

    for i in sorted(tbp, reverse=True): fnl.pop(i)

    return fnl
