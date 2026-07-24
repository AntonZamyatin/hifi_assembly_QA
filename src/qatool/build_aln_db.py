"""
Project: hifi assembly QA
Module: bulid_aln_db.py
Description: 

Author: Anton Zamyatin
Affiliation: Pennsylvania State University
Contact: zamyatin.anton@gmail.com

Created: 2/20/26

License: MIT License
"""
import argparse
import pysam
import sqlite3
from qatool.aln_utils import cigar_stats

parser = argparse.ArgumentParser()
parser.add_argument("bam", help="Path to BAM file")
parser.add_argument("db", help="Path to output SQLite DB")
parser.add_argument("asm", help="Assembly name")

SCHEMA_SQL = """
PRAGMA journal_mode = WAL;
PRAGMA synchronous = NORMAL;
PRAGMA temp_store = MEMORY;

CREATE TABLE IF NOT EXISTS aln (
    assembly    TEXT NOT NULL,
    read_name   TEXT NOT NULL,
    contig      TEXT,
    start       INTEGER,
    end         INTEGER,
    strand      INTEGER,     -- 1 forward, 0 reverse, NULL unmapped
    cigar       TEXT,
    pid         REAL,
    soft_left   INTEGER NOT NULL,
    soft_right  INTEGER NOT NULL,
    hard_left   INTEGER NOT NULL,
    hard_right  INTEGER NOT NULL,
    q_aln_len   INTEGER NOT NULL,
    q_len       INTEGER NOT NULL,
    ed          INTEGER,
    plen        REAL,
    aln_type    TEXT NOT NULL
);

-- Fast lookup by (assembly, read_name)
CREATE INDEX IF NOT EXISTS aln_by_read ON aln(assembly, read_name);

-- Optional: helpful for some downstream queries
CREATE INDEX IF NOT EXISTS aln_by_locus ON aln(assembly, contig, start, end);
"""

def get_aln_type(aln: pysam.AlignedSegment) -> str:
    if aln.is_supplementary:
        return "supplementary"
    if aln.is_secondary:
        return "secondary"
    return "primary_wsup" if aln.has_tag("SA") else "primary"

def build_alignment_db_from_bam(
    bam_path: str,
    db_path: str,
    assembly: str,
    *,
    batch_size: int = 50_000,
) -> None:
    """
    Build a persistent SQLite DB for random access by read_name.
    """
    con = sqlite3.connect(db_path)
    try:
        con.executescript(SCHEMA_SQL)
        cur = con.cursor()

        insert_sql = """
        INSERT INTO aln (
            assembly, read_name, contig, start, end, strand, cigar, pid,
            soft_left, soft_right, hard_left, hard_right,
            q_aln_len, q_len, ed, plen, aln_type
        ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        """

        rows = []
        n_written = 0

        with pysam.AlignmentFile(bam_path, "rb") as bam:
            # until_eof=True ensures we iterate the whole file regardless of coordinate fetch
            for aln in bam.fetch(until_eof=True):
                if aln.is_unmapped:
                    continue
                atype = get_aln_type(aln)
                contig = bam.get_reference_name(aln.reference_id)
                start = int(aln.reference_start)
                end = int(aln.reference_end)
                strand = 1 if aln.is_forward else 0
                cigar = aln.cigarstring if aln.cigarstring is not None else None
                soft_left, soft_right, hard_left, hard_right, q_aln_len, q_len, plen, pid, ed = cigar_stats(aln)

                rows.append((
                    assembly, aln.query_name,
                    contig, start, end, strand, cigar, pid,
                    int(soft_left), int(soft_right), int(hard_left), int(hard_right),
                    int(q_aln_len), int(q_len), ed,
                    (None if plen is None else float(plen)),
                    atype
                ))

                if len(rows) >= batch_size:
                    cur.executemany(insert_sql, rows)
                    con.commit()
                    n_written += len(rows)
                    rows.clear()

            if rows:
                cur.executemany(insert_sql, rows)
                con.commit()
                n_written += len(rows)
                rows.clear()

        # Analyze helps SQLite planner
        con.execute("ANALYZE;")
        con.commit()
        print(f"[build_alignment_db_from_bam] wrote {n_written:,} rows for assembly={assembly} into {db_path}")
    finally:
        con.close()

if __name__ == "__main__":
    args = parser.parse_args()
    build_alignment_db_from_bam(args.bam, args.db, args.asm)
