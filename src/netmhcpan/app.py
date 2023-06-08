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
        help = "The path to the browser driver.",
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


def init_selenium(browser_binary_filepath: utils.PathType, driver_filepath: utils.PathType) -> utils.WebDriverType:
    options = ChromeOptions()
    options.add_argument("--headless")
    options.page_load_strategy = "eager"
    options.binary_location = str(browser_binary_filepath)
    return webdriver.Chrome(str(driver_filepath), options = options)


async def run(args: utils.NetMHCPanCrawlerArgs | argparse.Namespace) -> pd.DataFrame:
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

    wait_start_time = time.time()
    elapsed_seconds = 0.0
    max_time = 60 * 4
    while elapsed_seconds < max_time:
        data = await craw.get_data(job_id)
        if data is not None:
            break
        await asyncio.sleep(5)
        elapsed_seconds += time.time() - wait_start_time

    return data


async def amain(args: argparse.Namespace) -> pd.DataFrame:
    data = await run(args)
    return data


def main() -> int:
    args = process_args()
    data = asyncio.run(amain(args))
    print(data)
    return 0