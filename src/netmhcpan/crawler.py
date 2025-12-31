import asyncio

from bs4 import BeautifulSoup

from dataclasses import dataclass

import httpx

from netmhcpan import utils

import pandas as pd

import re

from selenium import webdriver
from selenium.webdriver import Chrome
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.common import exceptions
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import Select

import time

import typing


@dataclass
class NetMHCPanCrawler:
    mhc_data: utils.MHCClassData
    driver: utils.WebDriverType = None
    peptides: str = None
    alleles: str = None
    
    _jobid_regex: typing.ClassVar[str] = re.compile(r"(?<=jobid=)([A-Z0-9]+?)&")
    _decimal_re: typing.ClassVar[str] = re.compile(r"[0-9]*\.?[0-9]*")

    def set_peptides(self, peptides: typing.Iterable[str]) -> None:
        if len(peptides) > 5000:
            raise ValueError("The number of peptides must be less than 5000.")
        self.peptides = "\n".join(peptides)
    
    def set_alleles(self, alleles: typing.Iterable[str]) -> None:
        if len(alleles) > 20:
            raise ValueError("The number of alleles must be less than 20.")
        self.alleles = ",".join(alleles)

    def _connect(self, url: str) -> None:
        if self.driver is None:
            raise ValueError("Driver is not set.")
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
        time.sleep(5)
        return re.search(self._jobid_regex, self.driver.current_url).group(1)

    async def query_job(self, job_id: str) -> pd.DataFrame | None:
        wait_start_time = time.monotonic()
        elapsed_seconds = 0.0
        max_time = 60 * 500
        while elapsed_seconds < max_time:
            data = await self.get_data(job_id)
            if data is not None:
                break
            await asyncio.sleep(5)
            elapsed_seconds += time.monotonic() - wait_start_time
            print(f"Waited: {elapsed_seconds}s.")
        return data
    
    async def get_data(self, job_id: str, data_path: str = None) -> pd.DataFrame | None:
        url_to_query = f"{self.mhc_data.results_url}?jobid={job_id}&wait=20"

        if data_path is not None:
            filepath = data_path / f"{job_id}.csv"
            if filepath.exists():
                return pd.read_csv(filepath)
        
        async with httpx.AsyncClient(timeout = 20) as c:
            request_data = await c.get(url_to_query)

        soup = BeautifulSoup(request_data.text, "html.parser")

        if soup is None:
            return 

        pre_html_content = soup.find("pre")

        if pre_html_content is None:
            return 
        
        pre_text = pre_html_content.text
        data_rows = list(self._parse_pre_text(pre_text, row_length = len(self.mhc_data.header_schema)))
        df = pd.DataFrame(data_rows, columns = self.mhc_data.header_schema)
        return pd.DataFrame(data_rows, columns = self.mhc_data.header_schema)
    
    def _parse_pre_text(self, pre_text: str, row_length: int) -> typing.Iterator[list[str]]:
        for line in pre_text.split("\n"):
            if line.startswith("-") or line.startswith("#") or line == "":
                continue
            if not line.strip()[0].isdigit():
                continue

            row = [x for x in line.split() if x != " " and (x.isalpha() or re.match(self._decimal_re, x))]
            if len(row) < row_length:
                row += ["None"] * (row_length - len(row))
            elif len(row) > row_length:
                row = row[: row_length]
            yield row

            