from itertools import chain

with open("examples/parameters.yml") as yml:
    comp = list(chain.from_iterable([i.split("\n") for i in yml.readlines()]))
    comp = list(chain.from_iterable([i.split(":") for i in comp]))
    med = [st.strip() if st != "- params" else '' for st in comp]
    fltrd = list(filter(None, med))
    print([st.replace("- ", "") if "- " in st else st for st in fltrd])