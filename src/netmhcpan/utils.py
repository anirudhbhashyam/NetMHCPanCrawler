from dataclasses import dataclass

import enum

import itertools

import numpy as np

import os

import pandas as pd

from pathlib import Path

from selenium import webdriver

import typing


PathType: typing.TypeAlias = str | os.PathLike
WebDriverType: typing.TypeAlias = webdriver.Chrome


HEADER_SCHEMA_CLASS_I = {
    "Pos": str,
    "MHC": str,
    "Peptide": np.int32,
    "Core": str,
    "Of": np.int32,
    "Gp": np.int32,
    "Gl": np.int32,
    "Ip": np.int32,
    "Il": np.int32,
    "Icore": str,
    "Identity": np.float32,
    "Score_EL": np.float32,
    "%Rank_EL": np.float32,
    "Score_BA": np.float32,
    "%Rank_BA": np.float32,
    "Aff(nM)": np.float32,
    "BindLevel": str,
}
HEADER_SCHEMA_CLASS_II = {
    "Pos": str,
    "MHC": str,
    "Peptide": str,
    "Of": np.int32,
    "Core": str,
    "Core_Rel": str,
    "Identity": np.float32,
    "Score_EL": np.float32,
    "%Rank_EL": np.float32,
    "Exp_Bind": np.float32,
    "Score_BA": np.float32,
    "Affinity(nM)": np.float32,
    "%Rank_BA": np.float32,
    "BindLevel": str,
}
MHC_CLASS_DATA = {
    "I": {
        "header_schema": HEADER_SCHEMA_CLASS_I,
        "job_url": "https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/",
        "results_url": "https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi",
    },
    "II": {
        "header_schema": HEADER_SCHEMA_CLASS_II,
        "job_url": "https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.0/",
        "results_url": "https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi",
    },
}


class MHCClassData(typing.NamedTuple):
    class_: str
    header_schema: dict[str, typing.Any]
    job_url: str
    results_url: str
    ba_checkbox_name: str


class NetMHCPanCrawlerArgs(typing.NamedTuple):
    peptides_filepath: PathType
    alleles_filepath: PathType
    browser_binary_filepath: PathType = "/Applications/Brave Browser.app/Contents/MacOS/Brave Browser"
    driver_filepath: PathType = "./chromedriver"
    mhc_class: str = "I"


# class Browser(enum.Enum):
#     CHROME = "chrome"
#     FIREFOX = "firefox"
#     EDGE = "edge"
#     SAFARI = "safari"


# @dataclass
# class SeleniumBrowser:
#     browser: Browser
#     binary_filepath: PathType
#     driver_filepath: PathType

#     def get_options() -> Options:
#         from selenium.webdriver.chrome.options import Options as ChromeOptions
#         from selenium.webdriver.firefox.options import Options as FirefoxOptions
#         from selenium.webdriver.edge.options import Options as EdgeOptions
#         from selenium.webdriver.safari.options import Options as SafariOptions

#         if self.browser == Browser.CHROME:
#             options = ChromeOptions()
#         elif self.browser == Browser.FIREFOX:
#             options = FirefoxOptions()
#         elif self.browser == Browser.EDGE:
#             options = EdgeOptions()
#         elif self.browser == Browser.SAFARI:
#             options = SafariOptions()
#         else:
#             raise ValueError("Browser must be one of chrome, firefox, edge, or safari.")


def read_txt(filepath: PathType) -> list[str] | None:
    try:
        with open(filepath, "r") as f:
            return f.read().splitlines()
    except FileNotFoundError:
        return


def save_data(df: pd.DataFrame, filepath: PathType) -> None:
    if not isinstance(filepath, Path):
        filepath = Path(filepath)

    if not filepath.exists() and not df.empty:
        df.to_csv(filepath)


def process_hla_prediction_data(
    peptide_hla_prediction_results: pd.DataFrame, 
    threshold: float, 
    full_peptide_seqs: typing.Iterable[str],
    mhc_class: str = "I"
) -> pd.DataFrame:
    if mhc_class not in ("I", "II"):
        raise ValueError("MHC class must be either I or II")
        
    if mhc_class == "I":
        aff_col = "Aff(nM)"
    else:
        aff_col = "Affinity(nM)"

    peptide_hla_prediction_hits = peptide_hla_prediction_results[
        peptide_hla_prediction_results[aff_col].astype(float) <= threshold
    ].copy().reset_index(drop = True)

    peptide_hla_prediction_hits["full_sequence"] = peptide_hla_prediction_hits["Peptide"].apply(
        lambda seq: next((x for x in full_peptide_seqs if seq in x), None)
    )

    new_index = pd.MultiIndex.from_frame(
        peptide_hla_prediction_hits[["full_sequence", "Peptide"]], 
        names = ("full_sequence", "Peptide"),
    )

    peptide_hla_prediction_hits = peptide_hla_prediction_hits\
        .set_index(new_index)\
        .drop(columns = ["full_sequence", "Peptide"])

    return peptide_hla_prediction_hits[["MHC", "Score_BA", aff_col]]


def n_gram_split(seq: str, n: int = 8) -> typing.Iterator[str]:
    for i in range(len(seq) - n + 1):
        yield seq[i : i + n]


def parse_pre_text(pre_text: str, row_length: int) -> typing.Iterator[list[str]]:
    for line in pre_text.split("\n"):
        if line.startswith(" " * 3):
            row = [x for x in line.split() if x != " " and (x.isalpha() or re.match(DECIMAL_RE, x))]
            if len(row) < row_length:
                row += ["None"] * (row_length - len(row))
            elif len(row) > row_length:
                row = row[: row_length]
            yield row


def split_sequence(optional_seq: str, option_delim: str = "/") -> typing.Iterator[str]:
    s = optional_seq.split()
    mask_s = {i: m.replace(option_delim, "") for i, m in enumerate(s) if option_delim in m}
    for option_seq in itertools.product(*mask_s.values()):
        for i, pos in enumerate(mask_s):
            s[pos] = option_seq[i]
        yield "".join(s)
