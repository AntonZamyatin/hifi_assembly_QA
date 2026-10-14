"""qatool.build_aln_db: an assembly loads whole or not at all, so "it has rows"
means "it is loaded" (prepare_data.sh skips an assembly that has rows).

Run: venv/bin/python tests/test_build_aln_db.py
"""

import os
import sqlite3
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fixtures import write_bam

import qatool.build_aln_db as build

_failures = []


def check(cond, msg):
    if not cond:
        _failures.append(msg)
        print(f"  FAIL: {msg}")
    return cond


def _rows(db, asm):
    con = sqlite3.connect(db)
    try:
        return con.execute("SELECT COUNT(*) FROM aln WHERE assembly=?", (asm,)).fetchone()[0]
    finally:
        con.close()


def test_load_is_all_or_nothing(tmp):
    bam = write_bam(os.path.join(tmp, "a.bam"),
                    [(f"r{i}", "ctgA", 100 * i, 50, 0, 0) for i in range(10)])
    db = os.path.join(tmp, "aln.sqlite")
    build.build_alignment_db_from_bam(bam, db, "A", batch_size=3)
    check(_rows(db, "A") == 10, f"a whole load: {_rows(db, 'A')}")

    real, n = build.cigar_stats, [0]

    def dies(aln):          # interrupted after two batches
        n[0] += 1
        if n[0] > 7:
            raise KeyboardInterrupt
        return real(aln)
    build.cigar_stats = dies
    try:
        build.build_alignment_db_from_bam(bam, db, "B", batch_size=3)
    except KeyboardInterrupt:
        pass
    finally:
        build.cigar_stats = real
    check(_rows(db, "B") == 0, f"an interrupted load leaves no rows: {_rows(db, 'B')}")
    check(_rows(db, "A") == 10, "the other assembly is untouched")
    con = sqlite3.connect(db)
    idx = {r[0] for r in con.execute("SELECT name FROM sqlite_master WHERE type='index'")}
    mode = con.execute("PRAGMA journal_mode").fetchone()[0]
    con.close()
    check({"aln_by_read", "aln_by_locus"} <= idx, f"indexes made: {idx}")
    check(mode == "wal", f"left in WAL mode: {mode}")


def main():
    tmp = tempfile.mkdtemp(prefix="qaview_alndb_tests_")
    tests = [test_load_is_all_or_nothing]
    for fn in tests:
        print(f"{fn.__name__} ...")
        fn(tmp)
    print()
    if _failures:
        print(f"{len(_failures)} FAILURE(S)")
        return 1
    print(f"all {len(tests)} tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
