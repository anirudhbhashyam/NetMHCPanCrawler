import argparse

import asyncio

from bs4 import BeautifulSoup

from dataclasses import dataclass

import httpx

import numpy as np

import os

import pandas as pd

from pathlib import Path

import re

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.common import exceptions
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import Select

import time

import typing

PathType = typing.TypeVar("PathType", str, os.PathLike)


MAIN_DATA_PATH = Path("../data").resolve()
DATA_PATH = Path("./dtu_mhc_data").resolve()
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


@dataclass
class NetMHCPanCrawler:
    driver: webdriver.Chrome
    mhc_data: MHCClassData
    peptides: str = None
    alleles: str = None
    
    _jobid_regex: typing.ClassVar[str] = re.compile(r"(?<=jobid=)([A-Z0-9]+?&)")

    def set_peptides(self, peptides: typing.Iterable[str]) -> None:
        if len(peptides) > 5000:
            raise ValueError("The number of peptides must be less than 5000.")
        self.peptides = "\n".join(peptides)
    
    def set_alleles(self, alleles: typing.Iterable[str]) -> None:
        if len(alleles) > 20:
            raise ValueError("The number of alleles must be less than 50.")
        self.alleles = ",".join(alleles)

    def _connect(self, url: str) -> None:
        # set a timeout for the driver.
        self.driver.set_page_load_timeout(200)
        try:
            self.driver.get(url)
        except exceptions.TimeoutException:
            self.driver.execute_script("window.stop();")
            print("Timeout for page load. Try again later.") 
            quit()
    
    def submit_job(self) -> str:
        self._connect(self.mhc_data.job_url)

        select_tag = self.driver.find_element(By.NAME, "inp")
        select_input_format = Select(select_tag)
        select_input_format.select_by_value("1")

        peptide_input = self.driver.find_element(By.NAME, "PEPPASTE")
        peptide_input.send_keys(self.peptides)

        allele_input = self.driver.find_element(By.NAME, "allele")
        allele_input.send_keys(self.alleles)

        ba_checkbox = self.driver.find_element(By.NAME, self.mhc_data.ba_checkbox_name)

        self.driver.execute_script("arguments[0].scrollIntoView();", ba_checkbox)
        self.driver.execute_script("arguments[0].click();", ba_checkbox)
        self.driver.execute_script("document.querySelector('input[type=\"submit\"]').click();")
        return re.search(self._jobid_regex, self.driver.current_url).group(1).replace("&", "")
    
    async def get_data(self, job_id: str, header_schema: typing.Iterable[str], data_path: str = DATA_PATH, overwrite: bool = False) -> pd.DataFrame | None:
        url_to_query = f"{self.mhc_data.results_url}?jobid={job_id}&wait=20"
        filepath = data_path / f"{job_id}.csv"
        
        if filepath.exists() and not overwrite:
            return pd.read_csv(filepath)
    
        async with httpx.AsyncClient(timeout = 20) as c:
            request_data = await c.get(url_to_query)

        soup = BeautifulSoup(request_data.text, "html.parser")

        pre_html_content = soup.find("pre")

        if soup is None:
            return 

        if pre_html_content is None:
            return 
        
        pre_text = pre_html_content.text
        data_rows = list(self._parse_pre_text(pre_text, row_length = len(header_schema)))
        df = pd.DataFrame(data_rows, columns = header_schema)
        await self._save_dtu_mhc_data(df, filepath, overwrite = overwrite)
        return df
    
    @staticmethod
    def _parse_pre_text(pre_text: str, row_length: int) -> typing.Iterator[list[str]]:
        decimal_re = re.compile(r"[0-9]*\.?[0-9]*")
        for line in pre_text.split("\n"):
            if line.startswith(" " * 3):
                row = [x for x in line.split() if x != " " and (x.isalpha() or re.match(decimal_re, x))]
                if len(row) < row_length:
                    row += ["None"] * (row_length - len(row))
                elif len(row) > row_length:
                    row = row[: row_length]
                yield row

    @staticmethod
    async def _save_dtu_mhc_data(df: pd.DataFrame, filepath: Path, overwrite: bool = False) -> None:
        if (not filepath.exists() and not df.empty) or overwrite:
            df.to_csv(filepath)
        await asyncio.sleep(0.1)


def read_txt(filepath: Path) -> list[str] | None:
    try:
        with open(filepath, "r") as f:
            return f.read().splitlines()
    except FileNotFoundError:
        return 
    
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


def process_args() -> argparse.Namespace:
    processor = argparse.ArgumentParser()

    processor.add_argument(
        "-p",
        "--peptides_filepath",
        required = True,
        help = "The peptides to predict.",
    )

    processor.add_argument(
        "-a",
        "--alleles_filepath",
        required = True,
        help = "The alleles to predict.",
    )

    processor.add_argument(
        "--browser_binary_filepath",
        "-bbf",
        type = str,
        default = "/Applications/Brave Browser.app/Contents/MacOS/Brave Browser",
        help = "The path to the browser binary.",
    )

    processor.add_argument(
        "--driver_filepath",
        "-df",
        type = str,
        default = "./chromedriver",
        help = "The path to the chromedriver.",
    )

    processor.add_argument(
        "--mhc_class",
        "-mhc",
        type = str,
        choices = ["I", "II"],
        default = "I",
        help = "The MHC class for which prediction is done.",
    )

    return processor.parse_args()


def init_selenium(browser_binary_filepath: Path, driver_filepath: Path) -> webdriver.Chrome:
    options = Options()
    # options.add_argument("--headless")
    options.page_load_strategy = "eager"
    options.binary_location = str(browser_binary_filepath)
    return webdriver.Chrome(str(driver_filepath), options = options)


async def run(args: NetMHCPanCrawlerArgs) -> pd.DataFrame:
    driver = init_selenium(
        Path(args.browser_binary_filepath),
        Path(args.driver_filepath)
    )

    mhc_data = MHCClassData(
        class_ = args.mhc_class,
        job_url = MHC_CLASS_DATA[args.mhc_class]["job_url"],
        results_url = MHC_CLASS_DATA[args.mhc_class]["results_url"],
        header_schema = MHC_CLASS_DATA[args.mhc_class]["header_schema"],
        ba_checkbox_name = "BApred" if args.mhc_class == "I" else "BA",
    )

    crawler = NetMHCPanCrawler(
        driver,
        mhc_data,
        MHC_CLASS_DATA[args.mhc_class]["job_url"],
        MHC_CLASS_DATA[args.mhc_class]["results_url"],
    )

    peptides = read_txt(Path(args.peptides_filepath))
    alleles = read_txt(Path(args.alleles_filepath))

    crawler.set_peptides(peptides)
    crawler.set_alleles(alleles)

    job_id = crawler.submit_job()

    print(f"Submitted job: {job_id}. Waiting for results...")

    while True:
        data = await crawler.get_data(job_id, mhc_data.header_schema)
        if data is not None:
            break
        await asyncio.sleep(5)

    return data


async def main(args: argparse.Namespace) -> int:
    data = await run(args)
    print(data)
    return 0


if __name__ == "__main__":
    args = process_args()
    raise SystemExit(asyncio.run(main(args)))