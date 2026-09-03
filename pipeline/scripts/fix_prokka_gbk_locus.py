#!/usr/bin/env python3
"""Rewrite Prokka GenBank LOCUS lines so Biopython (and PhiSpy 5) can parse them.

Prokka often emits a two-line LOCUS header (name+bp on line 1, date on line 2)
with a long contig id; Biopython's parser rejects that layout. Every record in a
multi-contig .gbk can start this way; this script fixes each LOCUS block in order.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

DATE_ONLY = re.compile(r"^\d{2}-[A-Z]{3}-\d{4}\s*$")


def rewrite_locus_line(l0: str) -> str:
    name_m = re.match(r"LOCUS\s+(\S+)", l0)
    raw = name_m.group(1) if name_m else "contig"

    # Prokka contig ids often look like NODE_1_length_419635_cov_12.34 — do not
    # truncate the whole token to 16 chars (that yields NODE_1_length_41).
    ml = re.match(r"(.+)_length_(\d+)_cov", raw)
    if ml:
        base = ml.group(1)
        length = ml.group(2)
    else:
        m_short = re.match(r"(.+)_length_(\d+)$", raw)
        if m_short:
            base = m_short.group(1)
            length = m_short.group(2)
        else:
            mlen = re.search(r"length_(\d+)_cov", raw)
            length = mlen.group(1) if mlen else ""
            if not length:
                mbp = re.search(r"(\d+)\s+bp", l0)
                length = mbp.group(1) if mbp else "0"
            base = raw.split("_length_")[0] if "_length_" in raw else raw

    short = re.sub(r"[^a-zA-Z0-9_]", "_", base)[:16] or "contig"
    return f"LOCUS       {short}              {length} bp"


def main() -> None:
    path = Path(sys.argv[1])
    lines = path.read_text().splitlines()
    out: list[str] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("LOCUS"):
            new_locus = rewrite_locus_line(line)
            out.append(new_locus)
            i += 1
            if i < len(lines) and DATE_ONLY.match(lines[i]):
                i += 1
            continue
        out.append(line)
        i += 1
    text = "\n".join(out)
    if not text.endswith("\n"):
        text += "\n"
    sys.stdout.write(text)


if __name__ == "__main__":
    main()
