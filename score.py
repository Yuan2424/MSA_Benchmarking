from pathlib import Path
import subprocess
import csv
import re
import math
import shutil

# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------

REFERENCE_DIR = Path("References")
OUTPUT_DIR = Path("Outputs")
RESULT_DIR = Path("Results")

FASTSP_JAR = Path("FastSP.jar")

RESULT_DIR.mkdir(exist_ok=True)


# ---------------------------------------------------------------------
# FastSP scoring
# ---------------------------------------------------------------------

def score_alignment(reference: Path, prediction: Path):
    """
    Scores a predicted alignment against a reference alignment
    using FastSP.

    Parameters
    ----------
    reference : Path
        Reference FASTA alignment.

    prediction : Path
        Predicted FASTA alignment.

    Returns
    -------
    dict
        {
            "SP": float,
            "TC": float
        }
    """

    print("Java:", shutil.which("java"))
    result = subprocess.run(
        [
            "java",
            "-jar",
            str(FASTSP_JAR),
            "-r",
            str(reference),
            "-e",
            str(prediction)
        ],
        capture_output=True,
        text=True,
        check=True
    )

    output = result.stdout

    sp_match = re.search(r"SP-Score\s+(\S+)", output)
    tc_match = re.search(r"TC\s+(\S+)", output)

    if sp_match is None or tc_match is None:
        raise RuntimeError(
            "Could not parse FastSP output:\n\n" + output
        )

    sp = float(sp_match.group(1))
    tc = float(tc_match.group(1))

    return {
        "SP": sp,
        "TC": tc
    }


# ---------------------------------------------------------------------
# Score every alignment
# ---------------------------------------------------------------------

def score_all_alignments():
    """
    Scores every alignment in the Outputs folder.

    Results are written to:
        Results/alignment_scores.csv
    """

    results = []

    predictions = sorted(OUTPUT_DIR.glob("*.fasta"))

    if not predictions:
        print("No prediction alignments found.")
        return

    for prediction in predictions:

        # Example:
        # BB11035_mafft.fasta
        #
        # dataset = BB11035
        # aligner = mafft

        dataset, aligner = prediction.stem.rsplit("_", 1)

        reference = REFERENCE_DIR / f"{dataset}.fasta"

        if not reference.exists():
            print(f"Reference not found for {prediction.name}")
            continue

        scores = score_alignment(reference, prediction)

        results.append({
            "Dataset": dataset,
            "Aligner": aligner,
            "SP": scores["SP"],
            "TC": scores["TC"]
        })

        sp_text = (
            "NaN"
            if math.isnan(scores["SP"])
            else f"{scores['SP']:.4f}"
        )

        tc_text = (
            "NaN"
            if math.isnan(scores["TC"])
            else f"{scores['TC']:.4f}"
        )

        print(
            f"{dataset:<10}"
            f"{aligner:<10}"
            f"SP={sp_text:<8}"
            f"TC={tc_text}"
        )


    csv_path = RESULT_DIR / "alignment_scores.csv"

    with open(csv_path, "w", newline="") as f:

        writer = csv.DictWriter(
            f,
            fieldnames=[
                "Dataset",
                "Aligner",
                "SP",
                "TC"
            ]
        )

        writer.writeheader()
        writer.writerows(results)

    print(f"\nResults written to {csv_path}")