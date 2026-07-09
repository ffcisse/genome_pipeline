"""Shared FASTA parsing and sequence-cleaning helpers.

Ported from the exploratory notebooks (01_proteome_overview_FINAL.ipynb for
protein sequences, dna_seq_analysis.ipynb for CDS sequences) so every pipeline
stage reads FASTA the same validated way.
"""

import gzip
import io
import tarfile
from pathlib import Path

from Bio import SeqIO

VALID_AAS = set("ACDEFGHIKLMNPQRSTVWY")


def clean_sequence(seq: str) -> str:
    """Remove stop codons and any non-standard residues (e.g. X, B, Z, U)."""
    seq = seq.replace("*", "")
    return "".join(aa for aa in seq if aa in VALID_AAS)


def iter_fasta_records(path):
    """Yield Bio.SeqRecord objects from a FASTA file that may be plain,
    gzipped, or a gzipped tar containing exactly one FASTA member -- the
    three formats found across the real Mycocosm downloads."""
    path = Path(path)
    name = path.name
    if name.endswith(".tar.gz") or name.endswith(".tgz"):
        with tarfile.open(path, "r:gz") as tar:
            member = next(
                m for m in tar.getmembers() if ".fasta" in m.name or m.name.endswith(".fa")
            )
            handle = io.TextIOWrapper(tar.extractfile(member), encoding="utf-8")
            yield from SeqIO.parse(handle, "fasta")
    elif name.endswith(".gz"):
        with gzip.open(path, "rt") as handle:
            yield from SeqIO.parse(handle, "fasta")
    else:
        with open(path) as handle:
            yield from SeqIO.parse(handle, "fasta")
