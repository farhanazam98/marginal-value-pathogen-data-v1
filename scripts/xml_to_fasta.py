#!/usr/bin/env python3
"""Stream-convert a UniRef100 XML archive (.xml.gz) to FASTA.

Every historical UniProt release ships UniRef100 as XML only (the standalone
uniref100.fasta.gz exists solely for whichever release is currently live), so
every year in the temporal study has to be built from XML. This script is the
calibration target for that conversion cost.

Design constraints, in order of importance:

1. Streaming. Decompress, parse, and write in one pass. Decompressed XML never
   lands on disk, so peak disk is compressed input plus FASTA output, nothing
   more. The 2010 file is 2.7 GB compressed and roughly 20 GB decompressed, so
   this is not a micro-optimization.
2. No DOM parse. UniRef100 XML is flat and predictable (verified empirically:
   one <entry>, one <representativeMember>, one <sequence> per cluster), so a
   line-oriented scan is both faster and lower-memory than lxml.iterparse with
   element clearing, and needs no clearing discipline to stay flat in RAM.
3. One representative sequence per cluster. Verified against the 2010 file that
   <member> blocks carry no <sequence> element (941,729 entries and 941,729
   <representativeMember> against 390,145 <member> blocks yielded exactly
   941,729 <sequence> tags), so taking every <sequence> is correct and does not
   double-count multi-member clusters.

Memory stays flat at roughly the size of the single largest sequence, because
the only accumulator is the current sequence's line buffer.

Two <sequence> layouts occur: 2010-2018 wrap residues over several lines with
</sequence> on its own line; 2020+ put the whole sequence on one line. Both
are handled. A sequence body that hits another tag or EOF before closing now
raises, rather than latching `in_seq` on and buffering the rest of the file --
the unbounded-RAM failure of the multi-line-only version on 2020+ input.
"""

import argparse
import gzip
import json
import re
import sys
import time

# The XML declares ISO-8859-1, not UTF-8. Organism names carry accented
# characters, so decoding as UTF-8 would raise or silently mangle bytes.
XML_ENCODING = "ISO-8859-1"

LEN_ATTR = re.compile(r'length="(\d+)"')


def convert(xml_gz_path, fasta_path, progress_every=25_000_000):
    """Convert one .xml.gz to FASTA. Returns a stats dict."""
    n_entries = 0          # <entry> tags seen
    n_seqs = 0             # sequences written
    n_residues = 0         # total residues written
    n_len_mismatch = 0     # declared length != observed length
    max_seq_len = 0
    fasta_bytes = 0

    acc = None             # current cluster id, e.g. UniRef100_Q197F8
    desc = None            # current cluster name, e.g. "Cluster: ..."
    declared_len = None
    in_seq = False
    buf = []

    t0 = time.time()
    last_report = 0

    with gzip.open(xml_gz_path, "rt", encoding=XML_ENCODING) as fin, \
         open(fasta_path, "w", encoding="ascii", newline="\n") as fout:

        for line in fin:
            t = line.strip()
            flush = False

            if in_seq:
                # Multi-line layout: residue lines (letters only) until a bare
                # </sequence>. Any other tag here means it never closed -- raise
                # rather than latch on and buffer the rest of the file.
                if t == "</sequence>":
                    in_seq = False
                    flush = True
                elif t.startswith("<"):
                    raise RuntimeError(
                        f"sequence for {acc!r} not closed before {t[:40]!r}"
                    )
                else:
                    buf.append(t)
            # Ordered by descending frequency in the file to keep the branch
            # chain short on the common path.
            elif t.startswith("<sequence "):
                m = LEN_ATTR.search(t)
                declared_len = int(m.group(1)) if m else None
                gt = t.find(">")
                end = t.find("</sequence>", gt + 1)
                if end == -1:
                    in_seq = True  # residues continue on following lines
                else:
                    # Single-line layout (2020+): emit residues as one FASTA
                    # line, no rewrap.
                    buf.append(t[gt + 1:end])
                    flush = True
            elif t.startswith("<entry "):
                n_entries += 1
                # id="UniRef100_Q197F8"
                i = t.find('id="')
                if i != -1:
                    j = t.find('"', i + 4)
                    acc = t[i + 4:j]
            elif t.startswith("<name>"):
                # Only the entry-level <name> appears at this nesting; the
                # representative's protein name is a <property>, not a <name>.
                if desc is None:
                    desc = t[6:-7] if t.endswith("</name>") else None

            if flush:
                seq_len = sum(len(chunk) for chunk in buf)
                header = f">{acc} {desc}\n" if desc else f">{acc}\n"
                body = "\n".join(buf) + "\n"
                fout.write(header)
                fout.write(body)

                fasta_bytes += len(header) + len(body)
                n_seqs += 1
                n_residues += seq_len
                if seq_len > max_seq_len:
                    max_seq_len = seq_len
                if declared_len is not None and declared_len != seq_len:
                    n_len_mismatch += 1

                buf = []
                acc = None
                desc = None
                declared_len = None

            if progress_every and n_seqs and n_seqs - last_report >= progress_every:
                last_report = n_seqs
                print(f"  ... {n_seqs:,} sequences, {n_residues:,} residues",
                      file=sys.stderr, flush=True)

    if in_seq:
        raise RuntimeError(
            f"reached EOF inside an unterminated sequence for {acc!r} -- "
            "truncated input or XML line-layout assumption violated"
        )

    elapsed = time.time() - t0

    return {
        "xml_gz_path": xml_gz_path,
        "fasta_path": fasta_path,
        "n_entries": n_entries,
        "n_seqs": n_seqs,
        "n_residues": n_residues,
        "n_len_mismatch": n_len_mismatch,
        "max_seq_len": max_seq_len,
        "fasta_bytes": fasta_bytes,
        "wall_seconds": round(elapsed, 2),
        "mean_seq_len": round(n_residues / n_seqs, 2) if n_seqs else None,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("xml_gz")
    ap.add_argument("fasta_out")
    ap.add_argument("--stats-json", default=None,
                    help="write the stats dict here as JSON")
    args = ap.parse_args()

    stats = convert(args.xml_gz, args.fasta_out)

    for k, v in stats.items():
        print(f"{k}: {v}")

    if args.stats_json:
        with open(args.stats_json, "w") as f:
            json.dump(stats, f, indent=2)
        print(f"wrote {args.stats_json}")


if __name__ == "__main__":
    main()
