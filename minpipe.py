import argparse
from pathlib import Path

from minpipe.parser import json, yaml
from minpipe.pipeline import PipelineCreator

def main():
    # # # # # # # # # # # # # # # # # #
    # Building argparse
    # # # # # # # # # # # # # # # # # #
    parser = argparse.ArgumentParser(description="Run full pipeline for every samples given.")
    parser.add_argument("-s", "--samples", nargs="+", required=False,
                        help="<Optional> List of samples to iterate over. All arguments can be passed as Json.",
                        type=str)
    parser.add_argument("-c", "--complement", nargs="+", required=False,
                        help="<Optional> Complementary for paired-ended", type=str)
    parser.add_argument("-i", "--index", nargs="?", required=False,
                        help="<Optional> Name of the index file to be used in pseudoalignment. Either `index` or \
            `transcript`has to be passed.")
    parser.add_argument("-t", "--transcript", nargs="?", required=False,
                        help="<Optional> Name of the transcript file to be indexed. `mmu` or `hsa` can be passed so \
            the transcript will be downloaded automatically and index will be built.")
    parser.add_argument("-f", "--format", nargs="?", required=False, help="<Optional> File format to be used in the \
            end of file names.")
    parser.add_argument("--threads", nargs="?", required=False, default="1",
                        help="<Optional> Number of threads to be used in quantification for Kallisto. Default: 1.")
    parser.add_argument("-b", "--bootstrap", nargs="?", required=False, default="100",
                        help="<Optional> Number of bootstrap samples. Default: 100.")
    parser.add_argument("--single", action='store_true', required=False,
                        help="<Optional> Flag to indicate single-ended quantification without complements.")
    parser.add_argument("--ext-qc", action="store_true", required=False,
                        help="<Optional> Flag to indicate that will have extensive QC. **MAY NEED MORE FILES**")
    parser.add_argument("--json", nargs=1, required=False,
                        help="<Optional> Use a Json file to pass all the arguments instead of command line interface. \
            See the parameters.json inside the examples folder to understand how to pass args.")
    parser.add_argument("--yaml", nargs=1, required=False,
                        help="<Optional> Use a YAML/YML file to pass all the arguments instead of \
                        command line interface. \
                        See the parameters.yml inside the examples folder to understand how to pass args.")

    args = parser.parse_args()

    # # # # # # # # # # # # # # # # # #
    # Handling YAML/JSON/CLI parsing
    # # # # # # # # # # # # # # # # # #
    if args.json is not None and args.yaml is not None:
        print("Both YAML and Json detected. Please pass only 1 file through argument handling.")
        exit()
    elif args.json is not None:
        try:
            Path(f"input/{args.json[0]}").is_file()
        except FileNotFoundError as fnfe:
            print(fnfe)
            quit()
        print("Json file found. Proceeding with Json parser for arguments.")
        lst_args = json.json_parse(args.json[0])
        args = parser.parse_args(lst_args)
    elif args.yaml is not None:
        try:
            Path(f"input/{args.yaml[0]}").is_file()
        except FileNotFoundError as fnfe:
            print(fnfe)
            quit()
        print("YAML file found. Proceeding with YAML parser for arguments.")
        lst_args = yaml.yaml_parse(args.yaml[0])
        args = parser.parse_args(lst_args)
    else:
        pass

    pipe = PipelineCreator(samples=args.samples, single=args.single, complement=args.complement,
        file_format=args.file_format, output_path=args.output_path, input_path=args.input_path,
        index=args.index, transcript=args.transcript, threads=args.threads, bootstrap=args.bootstrap,
        min_len=args.min_len, quality=args.quality, ext_qc=args.ext_qc)
    pipe.run_full() #min_len, quality, ext_qc, bootstrap, threads

    return


if __name__ == '__main__':
    main()
