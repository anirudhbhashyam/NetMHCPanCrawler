from dataclasses import dataclass

import enum

import os

import numpy as np

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


async def _save_dtu_mhc_data(df: pd.DataFrame, filepath: PathType, overwrite: bool = False) -> None:
    if not isinstance(filepath, Path):
        filepath = Path(filepath)

    if (not filepath.exists() and not df.empty) or overwrite:
        df.to_csv(filepath)

    await asyncio.sleep(0.1)



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

