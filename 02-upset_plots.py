import glob
import os
import subprocess
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import upsetplot

ROOT_DIR = "results/4_5_hours_regular/03_peak_annotations/subsetted_peaks/Promoter"
OUTPUT_DIR = "UpSet_Plots/promoters"
PREFIX_ORDER = ["DMSO", "CCT", "VEH", "INH"]
BLACKLIST_FACTORS = []
SHELL_CMD = """
bedtools multiinter -i {files} -header -names {names} > {output_path}
"""
COLUMNS_TO_IGNORE = ["chrom", "start", "end", "num", "list"]


def parse_name(filename):
    """
    Parse the file name to extract relevant information.

    Parameters
    ----------
    filename : str
        Name of the file.

    Returns
    -------
    str
        Parsed name in the format Prefix-Name.
    """
    stem = Path(filename).stem
    return f"{(part:=stem.split('_'))[0]}_{part[1]}"


def multi_intersect(files, names, output_path="./tmp.bed"):
    """
    Perform multi-intersection of BED files and save the output.

    Parameters
    ----------
    files : list
        List of file paths to be processed.
    names : list
        List of names corresponding to files.
    output_path : str, optional
        Output path for the multi-intersection file. Default is "./tmp.bed".

    Returns
    -------
    None
    """
    formatted_cmd = SHELL_CMD.format(
        files=" ".join(files), names=" ".join(names), output_path=output_path
    )
    result = subprocess.run(formatted_cmd, shell=True, capture_output=True)
    if result.returncode != 0:
        print(f"Error executing command for directory {dir}:")
        print(result.stderr.decode())
    else:
        print(f"Output saved for {dir} in {output_path}")


def visualize_upset(multi_intersect_output="./tmp.bed", output=None):
    """
    Visualize the multi-intersection output using an upset plot.

    Parameters
    ----------
    multi_intersect_output : str, optional
        Path to the multi-intersection output file. Default is "./tmp.bed".

    Returns
    -------
    None
    """
    df = pd.read_csv(multi_intersect_output, sep="\t")
    matrix = df.drop(columns=COLUMNS_TO_IGNORE)
    data = matrix.groupby(by=matrix.columns.to_list()).size()
    upsetplot.plot(data, sort_by="cardinality")

    if output is None:
        plt.show()
    else:
        plt.savefig(output, bbox_inches="tight", transparent=True)


def fetch_bed_files(root_dir="."):
    """
    Fetch all BED files from a directory.

    Parameters
    ----------
    root_dir : str, optional
        Root directory from which to start the search. Default is ".".

    Returns
    -------
    list
        List of paths to BED files.
    """
    return glob.glob(os.path.join(root_dir, "**/*.bed"), recursive=True)


def group_files(files):
    """
    Group files by prefixes defined in PREFIX_ORDER.

    Parameters
    ----------
    files : list
        List of file paths.

    Returns
    -------
    dict
        Dictionary with keys as prefixes and values as lists of file paths.
    """
    grouped_files = dict()
    for prefix in PREFIX_ORDER:
        grouped_files[prefix] = [
            file
            for file in files
            if Path(file).stem.startswith(prefix)
            and not any(factor in file for factor in BLACKLIST_FACTORS)
        ]
    return grouped_files


if __name__ == "__main__":
    bed_files = fetch_bed_files(ROOT_DIR)
    grouped_files = group_files(bed_files)
    for group, files in grouped_files.items():
        if len(files) == 0:
            continue
        names = [parse_name(filename) for filename in files]
        multi_intersect(files, names, output_path="./tmp.bed")

        output_filename = os.path.join(OUTPUT_DIR, group + "_upset.svg")
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        visualize_upset(multi_intersect_output="./tmp.bed", output=output_filename)
