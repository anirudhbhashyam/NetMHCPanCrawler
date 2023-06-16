import argparse

import asyncio

from datetime import datetime

import multiprocessing

from netmhcpan import utils
from netmhcpan import crawler
from netmhcpan import peptide

import pandas as pd

from pathlib import Path

import subprocess

import typing


CPD = Path(__file__).parent
PEPTIDES_DATA_PATH = CPD / "peptide_spaces"
ALLELES_DATA_PATH = CPD / "allele_sets"
RESULTS_PATH = CPD / "results"

def process_args() -> argparse.Namespace:
    processor = argparse.ArgumentParser()

    processor.add_argument(
        "--allele_data_filepaths",
        type = str,
        nargs = "+",
        help = "The filepath to the alleles.",
        required = True,
    )

    processor.add_argument(
        "--mhc_class",
        "-mhc",
        type = str,
        choices = ["I", "II"],
        default = "I",
        help = "The MHC class to use.",
    )

    return processor.parse_args()


def run_netmhc(peptides_world_filepath: Path, allele_data_filepath: Path, i: int, mhc_class: str) -> None:
    subprocess.check_call(
        [
            "netmhc",
            "--peptides_filepath",
            peptides_world_filepath,
            "--alleles_filepath",
            allele_data_filepath,
            "--save_filepath",
            RESULTS_PATH / f"{peptides_world_filepath.stem}_predictions_{i}.csv",
            "--mhc_class",
            mhc_class,
        ]
    )


def main(args: argparse.Namespace) -> int:
    seqs = [
        "MWLSYFVASFRLFARTRSMWSFNPETNILLNVP",
        "MWLSYFIASFRLFARTRSMWSFNPETNILLNVP",
        "MWISYFVQSIRLFMRTGSWWSFNPETNCLLNVP",
        "IWILYFVNSIRLFIRTGSWWSFNPETNNLMCID",
        "MWIVYFVNSIRLFIRTGSWWSFNPETNNLMCID",
        "MWIVYFVNSIRLFIRTGSFWSFNPETNNLMCID",
    ]
    peptides = [peptide.Peptide(peptide.ConsensusSequence(seq)) for seq in seqs]

    peptides_world_filepath = PEPTIDES_DATA_PATH / f"peptides_{datetime.now().strftime('%Y%m%d%H')}_class_{args.mhc_class}.txt" 
    world_gen = peptide.create_peptides_world(peptides, 9) if args.mhc_class == "II" else peptide.create_peptides_world(peptides, 8, 10)
    with open(peptides_world_filepath, "w") as f:
        for pep in world_gen:
            f.write(pep + "\n")

    # with multiprocessing.Pool(4) as pool:
    #     pool.starmap(
    #         run_netmhc, 
    #         [
    #             (peptides_world_filepath, allele_data_filepath)
    #             for allele_data_filepath in args.allele_data_filepaths
    #         ]
    #     )

    for i, allele_data_filepath in enumerate(args.allele_data_filepaths):
        run_netmhc(peptides_world_filepath, allele_data_filepath, i, args.mhc_class)

    # Collate the prediction data.
    df = pd.concat(
        [
            pd.read_csv(filepath, index_col = 0)
            for filepath in RESULTS_PATH.glob(f"{peptides_world_filepath.stem}_*.csv")
        ],
        axis = 0,
    )
    df.to_csv(RESULTS_PATH / f"{peptides_world_filepath.stem}_predictions.csv")
        
    return 0


if __name__ == "__main__":
    args = process_args()
    raise SystemExit(main(args))