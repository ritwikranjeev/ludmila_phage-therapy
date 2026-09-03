#!/usr/bin/env python3
"""Extract FASTA records whose description looks like a hypothetical protein."""

import argparse
import sys
from pathlib import Path
from typing import Iterator, Tuple


def is_hypothetical_description(desc: str) -> bool:
    low = desc.lower()
    if "hypothetical" in low and "protein" in low:
        return True
    if "uncharacterized" in low and "protein" in low:
        return True
    return False


def iter_faa_records(path: Path) -> Iterator[Tuple[str, str]]:
    with path.open() as handle:
        header = None
        seq_chunks = []
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line.rstrip("\n")
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            yield header, "".join(seq_chunks)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input_faa", type=Path, help="Prokka .faa")
    ap.add_argument("output_faa", type=Path, help="Filtered .faa")
    args = ap.parse_args()

    n_in = n_out = 0
    args.output_faa.parent.mkdir(parents=True, exist_ok=True)
    with args.output_faa.open("w") as fout:
        for header, seq in iter_faa_records(args.input_faa):
            n_in += 1
            parts = header[1:].split(None, 1)
            rest = parts[1] if len(parts) > 1 else ""
            if is_hypothetical_description(rest):
                fout.write(header + "\n" + seq)
                n_out += 1

    print("Records scanned: {}, hypothetical written: {}".format(n_in, n_out), file=sys.stderr)


if __name__ == "__main__":
    main()
