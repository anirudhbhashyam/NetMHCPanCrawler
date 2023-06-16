from dataclasses import dataclass

import itertools

from netmhcpan import utils

import typing


@dataclass
class ConsensusSequence:
    seq: str
    delim: str = ""

    def __post_init__(self) -> None:
        if self.delim not in self.seq:
            raise ValueError("Delimiter not found in sequence.")
    
    def get_possibilities(self) -> typing.Iterator[str]:
        if self.delim == "":
            yield self.seq
            return
        s = self.seq.split()
        mask_s = {i: m.replace(self.delim, "") for i, m in enumerate(s) if self.delim in m}
        for option_seq in itertools.product(*mask_s.values()):
            for i, pos in enumerate(mask_s):
                s[pos] = option_seq[i]
            yield "".join(s)


class Peptide:
    def __init__(self, c_seq: ConsensusSequence) -> None:
        self.c_seq = c_seq

    def get_peptide_space(self, n_low: int, n_high: int = None) -> typing.Iterator[str]:
        for peptide in self.c_seq.get_possibilities():
            if n_high is None:
                n_high = len(peptide)
            for n in range(n_low, n_high + 1):
                yield from utils.n_gram_split(peptide, n)

    def __str__(self) -> str:
        return self.c_seq.seq

    def __repr__(self) -> str:
        return f"Peptide({self.c_seq.seq!r})"


def create_peptides_world(peptides: typing.Iterable[Peptide], n_low: int, n_high: int = None) -> typing.Iterator[str]:
    for peptide in peptides:
        yield from peptide.get_peptide_space(n_low, n_high)


def main() -> int:
    c_seq = ConsensusSequence("RLF I/M/C RTGSWWSFNPET N/C/M L", "/")
    # print(list(c_seq.get_possibilities()))
    peptide = Peptide(c_seq)
    list(peptide.get_peptide_space(8, 10))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())