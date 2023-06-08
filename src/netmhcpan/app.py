import argparse

import asyncio

from netmhcpan import crawler
from netmhcpan import utils

import pandas as pd

from pathlib import Path

from selenium import webdriver
from selenium.webdriver.chrome.options import Options as ChromeOptions

import sys

import time


def process_args() -> argparse.Namespace:
    processor = argparse.ArgumentParser()

    processor.add_argument(
        "-p",
        "--peptides_filepath",
        required = True,
        help = "The peptides to predict. Must be a txt file with new line separated peptides.",
    )

    processor.add_argument(
        "-a",
        "--alleles_filepath",
        required = True,
        help = "The alleles to predict. Must be a txt file with new line separated alleles. "
        "Additionally, the alleles must conform to the format and naming scheme adopted by NetMHCpan.",
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
        help = "The path to the browser driver.",
    )

    processor.add_argument(
        "--mhc_class",
        "-mhc",
        type = str,
        choices = ["I", "II"],
        default = "I",
        help = "The MHC class for which prediction should be performed.",
    )

    processor.add_argument(
        "--save_filepath",
        "-sf",
        type = str,
        default = None,
        help = "The filepath to save the results to. If not provided, the results will not be saved and will be written to stdout.",
    )

    return processor.parse_args()


def init_selenium(browser_binary_filepath: utils.PathType, driver_filepath: utils.PathType) -> utils.WebDriverType:
    options = ChromeOptions()
    options.add_argument("--headless")
    options.page_load_strategy = "eager"
    options.binary_location = str(browser_binary_filepath)
    return webdriver.Chrome(str(driver_filepath), options = options)


async def run(args: utils.NetMHCPanCrawlerArgs | argparse.Namespace) -> pd.DataFrame | None:
    driver = init_selenium(
        Path(args.browser_binary_filepath),
        Path(args.driver_filepath)
    )

    mhc_data = utils.MHCClassData(
        class_ = args.mhc_class,
        job_url = utils.MHC_CLASS_DATA[args.mhc_class]["job_url"],
        results_url = utils.MHC_CLASS_DATA[args.mhc_class]["results_url"],
        header_schema = utils.MHC_CLASS_DATA[args.mhc_class]["header_schema"],
        ba_checkbox_name = "BApred" if args.mhc_class == "I" else "BA",
    )

    craw = crawler.NetMHCPanCrawler(
        driver,
        mhc_data,
    )

    peptides = utils.read_txt(Path(args.peptides_filepath))
    alleles = utils.read_txt(Path(args.alleles_filepath))

    craw.set_peptides(peptides)
    craw.set_alleles(alleles)

    job_id = craw.submit_job()

    print(f"Submitted job: {job_id}. Waiting for results...")

    data = await craw.query_job(job_id)

    return data


async def amain(args: argparse.Namespace) -> pd.DataFrame:
    data = await run(args)
    return data


def main() -> int:
    args = process_args()
    data = asyncio.run(amain(args))
    if not args.save_filepath:
        print(data)
    else:
        utils.save_data(data, args.save_filepath)
    return 0