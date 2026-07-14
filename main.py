import subprocess
from pathlib import Path
from Bio import AlignIO
from score import score_all_alignments
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

REFERENCE_DIR = Path("References")

# CONVERTING MSF REFERENCES TO FASTA
def convert_reference(msf_file: Path, overwrite: bool = False) -> Path:
    """
    Converts a single MSF alignment to FASTA.

    Parameters
    ----------
    msf_file : Path
        Path to the .msf file.

    overwrite : bool
        If True, overwrite an existing FASTA file.

    Returns
    -------
    Path
        Path to the generated FASTA file.
    """

    fasta_file = msf_file.with_suffix(".fasta")

    if fasta_file.exists() and not overwrite:
        print(f"Skipping {msf_file.name} (already converted)")
        return fasta_file

    alignment = AlignIO.read(msf_file, "msf")
    AlignIO.write(alignment, fasta_file, "fasta")

    print(f"Converted {msf_file.name} -> {fasta_file.name}")

    return fasta_file


def convert_all_references(overwrite: bool = False):
    """
    Converts every .msf file in the References folder to FASTA.
    """

    if not REFERENCE_DIR.exists():
        raise FileNotFoundError(
            f"Reference directory '{REFERENCE_DIR}' does not exist."
        )

    msf_files = sorted(REFERENCE_DIR.glob("*.msf"))

    if not msf_files:
        print("No MSF files found.")
        return

    for msf_file in msf_files:
        convert_reference(msf_file, overwrite)

    print(f"\nConverted {len(msf_files)} reference alignments.")

# CONVERTING CLUSTAL ALIGNMENTS TO FASTA
def convert_clustal_to_fasta(clustal_file):
    alignment = AlignIO.read(clustal_file, "clustal")
    fasta_file = clustal_file.replace(".aln-clustal", "_clustal.fasta")
    fasta_file = fasta_file.replace("Clustal_raw", "Outputs")
    AlignIO.write(alignment, fasta_file, "fasta")

# MAKES A LIST OF THE TFA FILES IN THE DATA FOLDER
def get_tfa_files():
    """
    Returns a list of all .tfa files in the Data folder
    located in the current working directory.
    """
    data_dir = Path.cwd() / "Data"
    return list(data_dir.glob("*.tfa"))

def make_tfa_list(folder_name="Data"):
    lst = []
    for input_file in get_tfa_files():
        lst.append(folder_name + "/" + input_file.name)
    return lst

# MAFFT
def run_mafft(input_file):
    output_file = input_file.replace(".tfa", "_mafft.fasta")
    output_file = output_file.replace("Data", "Outputs")
    subprocess.run(
        ["mafft.bat", input_file],
        stdout=open(output_file, "w"),
        check=True
    )
    return output_file

# MUSCLE
def run_muscle(input_file):
    output_file = input_file.replace(".tfa", "_muscle.fasta")
    output_file = output_file.replace("Data", "Outputs")
    subprocess.run(
        ["muscle", "-align", input_file, "-output", output_file],
        check=True
    )
    return output_file

def main():

    # Load the results
    df = pd.read_csv("Results/alignment_scores.csv")

    sns.boxplot(data=df, x="Aligner", y="SP")

    plt.title("SP Scores by Aligner")
    plt.xlabel("Aligner")
    plt.ylabel("SP Score")

    plt.show()

    # summary = (
    #     df.groupby("Aligner")
    #     .agg(
    #         Mean_SP=("SP", "mean"),
    #         Median_SP=("SP", "median"),
    #         StdDev_SP=("SP", "std"),
    #         Mean_TC=("TC", "mean"),
    #         Median_TC=("TC", "median"),
    #         StdDev_TC=("TC", "std")
    #     )
    #     .round(4)
    # )

    # summary.to_csv("Results/summary_statistics.csv")

    """One time runs

    tfa_list = make_tfa_list()
    score_all_alignments()

    Convert all MSF reference files to FASTA
    convert_all_references()
    
    Convert all Clustal alignments to FASTA
    for tfa_file in tfa_list:
        edit_file = tfa_file.replace("Data", "Clustal_raw").replace(".tfa", ".aln-clustal")
        convert_clustal_to_fasta(edit_file)    
    
    Run MAFFT and MUSCLE on all TFA files
    for tfa_file in tfa_list:
        print(f"Processing file: {tfa_file}")
        mafft_output = run_mafft(tfa_file)
        print(f"MAFFT output saved to: {mafft_output}")

        muscle_output = run_muscle(tfa_file)
        print(f"MUSCLE output saved to: {muscle_output}")"""

if __name__ == "__main__":
    main()