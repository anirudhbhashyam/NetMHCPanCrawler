# NetMHCpan Crawler

A tool to run the NetMHCpan server(s) locally using selenium.
Currently the *Chrome* driver is supported. So binaries for browsers using *Chromium* maybe used. 


# Usage

```sh
git clone https://github.com/anirudhbhashyam/NetMHCPanCrawler
pip install .
netmhc --peptides_filepath <filepath> --alleles_filepath <filepath> --mhc_class [-mhc] <class>
```
The `mhc_class` can either be `I` or `II`. It is by default `I`. Try `netmhc --help` for more information. The command above will display the results of the prediction to stdout.
To save the data add an additional option.

```sh
netmhc --peptides_filepath <filepath> --alleles_filepath <filepath> --mhc_class [-mhc] <class> --save_filepath [-sf] <filepath>
```
The current behaviour supports saving to CSV.
