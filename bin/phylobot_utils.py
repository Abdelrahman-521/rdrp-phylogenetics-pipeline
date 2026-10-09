"""Shared helpers for the PhyloBot pipeline scripts.

Kept dependency-free (standard library only) so every step can import it.
"""

from __future__ import annotations

import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

EUTILS_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
SEQUENCE_CHARS = set("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz*-")


def read_list(path: str | Path) -> list[str]:
    """Read one item per line, skipping blank lines and # comments."""
    items = []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            items.append(line)
    return items


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    """Read a tab-separated file with a header row; # lines are comments."""
    rows: list[dict[str, str]] = []
    header: list[str] | None = None
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if header is None:
            header = fields
            continue
        rows.append(dict(zip(header, fields)))
    return rows


def parse_fasta_text(text: str) -> tuple[list[tuple[str, str]], list[str]]:
    """Split raw text into well-formed FASTA records and junk lines.

    NCBI sometimes mixes error messages into a download. Anything that is not a
    header or a sequence line belonging to a header is returned as junk instead
    of silently ending up in the alignment.
    """
    records: list[tuple[str, str]] = []
    junk: list[str] = []
    header: str | None = None
    seq: list[str] = []

    def flush() -> None:
        if header is not None:
            if seq:
                records.append((header, "".join(seq)))
            else:
                junk.append(">" + header)

    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            flush()
            header, seq = line[1:].strip(), []
        elif header is not None and set(line) <= SEQUENCE_CHARS:
            seq.append(line)
        else:
            junk.append(line)
    flush()
    return records, junk


def record_id(header: str) -> str:
    """First word of a FASTA header, as Biopython would report it."""
    return header.split()[0] if header.split() else header


def source_accession(rec_id: str) -> str:
    """Accession a record came from, without its version.

    'CAA25029.1' -> 'CAA25029'
    'lcl|MG833031.1_prot_AWC26954.1_1' -> 'MG833031'  (translated CDS from nuccore)
    """
    rid = rec_id.split("|", 1)[1] if "|" in rec_id else rec_id
    rid = rid.split("_prot_", 1)[0]
    return rid.split(".", 1)[0]


def accession_base(accession: str) -> str:
    return accession.split(".", 1)[0]


def missing_accessions(records: list[tuple[str, str]], accessions: list[str]) -> list[str]:
    present = {source_accession(record_id(h)) for h, _ in records}
    return [acc for acc in accessions if accession_base(acc) not in present]


def write_fasta(records: list[tuple[str, str]], path: str | Path, width: int = 70) -> None:
    with open(path, "w") as out:
        for header, seq in records:
            out.write(f">{header}\n")
            for i in range(0, len(seq), width):
                out.write(seq[i : i + width] + "\n")


def efetch(db: str, accession: str, rettype: str, *, email: str | None = None,
           api_key: str | None = None, retries: int = 4, pause: float = 0.4,
           timeout: float = 30.0) -> str:
    """Fetch one record from NCBI E-utilities with exponential backoff.

    Returns the response text ('' if NCBI answered but had nothing). Raises the
    last error if every attempt fails at the network level.
    """
    params = {"db": db, "id": accession, "rettype": rettype, "retmode": "text"}
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    url = EUTILS_EFETCH + "?" + urllib.parse.urlencode(params)
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as response:
                text = response.read().decode("utf-8", errors="replace")
            time.sleep(pause)  # stay under NCBI's 3 requests/second limit
            return text
        except urllib.error.HTTPError as err:
            # 400 means "this id is not in this database": not worth retrying.
            if err.code == 400:
                time.sleep(pause)
                return ""
            last_error = err
        except (urllib.error.URLError, TimeoutError, ConnectionError) as err:
            last_error = err
        time.sleep(pause * (2 ** (attempt + 1)))
    assert last_error is not None
    raise last_error
