# NetMHCpan Crawler

A tool to run the NetMHCpan server(s) locally using selenium.


# Usage

```sh
git clone https://github.com/anirudhbhashyam/NetMHCPanCrawler
pip install .
netmhc --peptides_filepath <filepath> --alleles_filepath <filepath> --mhc_class <class>
```
The `mhc_class` can either be `I` or `II`. It is by default `I`. Try `netmhc --help` for more information.
