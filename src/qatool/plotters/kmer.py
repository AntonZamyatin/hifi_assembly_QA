"""K-mer dot plots: a single-pane version and a per-segment-pair mosaic version."""

import math
import time

import numpy as np
import matplotlib.pyplot as plt


def dotplot_kmer(seq_a: str, seq_b: str, k: int = 15, binsize: int = 1000,
                 ax=None, max_cells: int = 10_000_000,
                 mask_ambiguous: bool = True, log_scale: bool = False,
                 max_kmer_freq: int = 100, verbose: bool = False):
    """
    Generate a Gepard-like k-mer dot plot using binned k-mer matches.

    Parameters:
    - seq_a: First sequence string.
    - seq_b: Second sequence string.
    - k: k-mer size.
    - binsize: Window size for binning positions.
    - max_cells: Max allowed dense matrix cells (bins_a * bins_b).
    - max_kmer_freq: k-mers appearing more than this many times in seq_a are
        excluded from the index (repeat / low-complexity filter). 0 = no filter.
        Default 100 is appropriate for genome-scale sequences; increase if the
        plot is missing signal in repetitive regions.
    - verbose: print k-mer counts and per-step timing.
    """
    if k <= 0:
        raise ValueError("k must be positive")
    if binsize <= 0:
        raise ValueError("binsize must be positive")

    seq_a = seq_a.upper()
    seq_b = seq_b.upper()
    len_a = len(seq_a)
    len_b = len(seq_b)
    if len_a < k or len_b < k:
        raise ValueError("Sequences must be at least k long")

    bins_a = math.ceil(len_a / binsize)
    bins_b = math.ceil(len_b / binsize)
    if bins_a * bins_b > max_cells:
        raise ValueError(
            f"dotplot matrix would be {bins_a}x{bins_b} cells; "
            "increase binsize or max_cells"
        )

    # Byte-level lookup: ASCII -> 2-bit base (A=0,C=1,G=2,T=3), ambiguous -> 255.
    enc_table = np.full(256, 255, dtype=np.uint8)
    for _b, _v in zip(b'ACGT', [0, 1, 2, 3]):
        enc_table[_b] = _v

    def _kmer_hashes(enc: np.ndarray, do_rc: bool) -> tuple:
        """
        Compute k-mer hashes for all valid positions using a k-iteration vectorised loop.
        Each iteration adds one shifted base, consuming only O(n) memory (no (n,k) array).
        Returns (positions, fwd_hashes, rc_hashes|None).
        """
        n = len(enc)
        n_pos = n - k + 1
        if n_pos <= 0:
            empty = np.array([], dtype=np.int64)
            return np.array([], dtype=np.intp), empty, empty if do_rc else None

        ph = 4 ** np.arange(k - 1, -1, -1, dtype=np.int64)  # fwd: 4^(k-1)..4^0
        pl = 4 ** np.arange(k,           dtype=np.int64)     # rc:  4^0..4^(k-1)
        ei = enc.astype(np.int64)

        fwd = np.zeros(n_pos, dtype=np.int64)
        rci = np.zeros(n_pos, dtype=np.int64) if do_rc else None
        for j in range(k):
            sl = ei[j: j + n_pos]
            fwd += sl * ph[j]
            if do_rc:
                rci += sl * pl[j]

        if mask_ambiguous:
            # Mark windows that contain any ambiguous base (255 in encoded form).
            inv_pos = np.where(enc == 255)[0]
            valid = np.ones(n_pos, dtype=bool)
            for p in inv_pos:
                valid[max(0, p - k + 1): min(n_pos, p + 1)] = False
            pos = np.where(valid)[0]
            fwd = fwd[valid]
            rc  = (np.int64(4 ** k - 1) - rci[valid]) if do_rc else None
        else:
            pos = np.arange(n_pos, dtype=np.intp)
            rc  = (np.int64(4 ** k - 1) - rci) if do_rc else None

        return pos, fwd, rc

    enc_a = enc_table[np.frombuffer(seq_a.encode('ascii'), dtype=np.uint8)]
    enc_b = enc_table[np.frombuffer(seq_b.encode('ascii'), dtype=np.uint8)]
    pos_a, fwd_a, rc_a = _kmer_hashes(enc_a, do_rc=True)
    pos_b, fwd_b, _    = _kmer_hashes(enc_b, do_rc=False)  # only fwd for seq_b

    bin_a = (pos_a // binsize).astype(np.int64)
    bin_b = (pos_b // binsize).astype(np.int64)

    # Pack (hash << SHIFT | bin) into a single int64 so np.sort (no index array)
    # replaces argsort + fancy-index.  Sorting the b-array also makes searchsorted
    # queries sequential → cache-friendly (avoids random jumps through a_packed).
    # SHIFT bits hold the bin index; hash occupies the upper bits.
    # Constraint: 2*k + SHIFT < 63 (int64 must not overflow).
    # For k<=20 SHIFT=21 supports up to 2M bins (~2Gbp at binsize=1000).
    SHIFT = max(int(max(bins_a, bins_b) - 1).bit_length(), 1)
    if 2 * k + SHIFT >= 63:
        # Fall back to argsort for unusually large k or tiny binsize.
        so_a    = np.argsort(np.concatenate([fwd_a, rc_a]), kind='stable')
        a_h_all = np.concatenate([fwd_a, rc_a])[so_a]
        a_b_all = np.concatenate([bin_a, bin_a])[so_a]
        so_b    = np.argsort(fwd_b, kind='stable')
        fwd_b_s = fwd_b[so_b]
        bin_b_s = bin_b[so_b]
        lo = np.searchsorted(a_h_all, fwd_b_s, side='left')
        hi = np.searchsorted(a_h_all, fwd_b_s, side='right')
        MASK = None
        def _get_ba_rep(lo_m, cnt_m):
            grp = np.empty(len(lo_m), dtype=np.int64); grp[0] = 0
            if len(lo_m) > 1: grp[1:] = np.cumsum(cnt_m[:-1])
            diffs = np.ones(int(cnt_m.sum()), dtype=np.int64); diffs[0] = lo_m[0]
            if len(lo_m) > 1: diffs[grp[1:]] = lo_m[1:] - lo_m[:-1] - cnt_m[:-1] + 1
            return a_b_all[np.cumsum(diffs)]
    else:
        MASK  = (1 << SHIFT) - 1
        a_h_c = np.concatenate([fwd_a, rc_a])
        a_b_c = np.concatenate([bin_a, bin_a])
        a_pk  = (a_h_c << SHIFT) | a_b_c;  a_pk.sort()

        # Repeat-filter: drop k-mers whose hash appears > max_kmer_freq times in
        # seq_a (satellite DNA, transposons, …).  These would explode total matches
        # and produce only uninformative noise in the plot.
        if max_kmer_freq > 0 and len(a_pk) > 0:
            a_hash_part = a_pk >> SHIFT
            run_ends    = np.concatenate([np.where(np.diff(a_hash_part))[0] + 1, [len(a_pk)]])
            run_starts  = np.concatenate([[0], run_ends[:-1]])
            run_lens    = run_ends - run_starts
            valid_runs  = run_lens <= max_kmer_freq
            keep        = np.zeros(len(a_pk) + 1, dtype=np.int8)
            keep[run_starts[valid_runs]] += 1
            keep[run_ends[valid_runs]]   -= 1
            a_pk = a_pk[np.cumsum(keep)[:-1].astype(bool)]
            n_filtered = valid_runs.size - valid_runs.sum()
            if verbose and n_filtered:
                print(f"  [dotplot_kmer] filtered {n_filtered} high-freq k-mer hashes "
                      f"(>{max_kmer_freq} occurrences)")

        # Collapse duplicate (hash, bin_a) pairs: multiple copies of the same k-mer
        # in the same bin (tandem repeats, satellite DNA) reduce to one entry.
        # a_pk is already sorted so np.unique is a O(n) linear scan.
        # Without this, a k-mer appearing N times in bin_a and M times in bin_b
        # would contribute N*M identical flat indices — pure wasted work.
        a_pk = np.unique(a_pk)

        b_pk  = (fwd_b << SHIFT) | bin_b;   b_pk.sort()
        b_pk  = np.unique(b_pk)
        fwd_b_s = b_pk >> SHIFT
        bin_b_s = b_pk & MASK
        lo = np.searchsorted(a_pk, fwd_b_s << SHIFT,       side='left')
        hi = np.searchsorted(a_pk, (fwd_b_s + 1) << SHIFT, side='left')
        def _get_ba_rep(lo_m, cnt_m):
            grp = np.empty(len(lo_m), dtype=np.int64); grp[0] = 0
            if len(lo_m) > 1: grp[1:] = np.cumsum(cnt_m[:-1])
            total = int(cnt_m.sum())
            diffs = np.ones(total, dtype=np.int64); diffs[0] = lo_m[0]
            if len(lo_m) > 1: diffs[grp[1:]] = lo_m[1:] - lo_m[:-1] - cnt_m[:-1] + 1
            return a_pk[np.cumsum(diffs)] & MASK

    cnt = hi - lo
    has_match = cnt > 0
    mat = np.zeros((bins_a, bins_b), dtype=np.uint32)
    _tt = time.perf_counter()
    if has_match.any():
        lo_m  = lo[has_match]
        cnt_m = cnt[has_match]
        bb_m  = bin_b_s[has_match]
        if verbose:
            print(f"  [dotplot_kmer] total k-mer matches to accumulate: {int(cnt_m.sum()):,}")
        ba_rep = _get_ba_rep(lo_m, cnt_m)
        bb_rep = np.repeat(bb_m, cnt_m)
        flat   = ba_rep * bins_b + bb_rep
        mat    = np.bincount(flat, minlength=bins_a * bins_b).reshape(bins_a, bins_b).astype(np.uint32)
        if verbose:
            print(f"  [dotplot_kmer]   accumulate: {time.perf_counter()-_tt:.2f}s")

    if ax is None:
        _, ax = plt.subplots(figsize=(16, 5))

    plot_mat = np.log1p(mat) if log_scale else mat
    nz = plot_mat[plot_mat > 0]
    vmin = 0.0
    vmax = float(np.percentile(nz, 99.5)) if nz.size > 0 else 1.0
    ax.imshow(plot_mat, origin='lower', aspect='equal', cmap='Greys', vmin=vmin, vmax=vmax)
    ax.set_xlabel("Seq B bins")
    ax.set_ylabel("Seq A bins")
    ax.set_title(f"K-mer dotplot (k={k}, binsize={binsize})")
    return ax


def dotplot_kmer_mosaic(seqs_a: list, seqs_b: list,
                        k: int = 15, binsize: int = 1000, ax=None,
                        max_cells: int = 10_000_000,
                        mask_ambiguous: bool = True, log_scale: bool = False,
                        max_kmer_freq: int = 100, verbose: bool = False):
    """
    K-mer dotplot rendered as a mosaic of independent per-segment-pair sub-dotplots.

    For each pair (seqs_a[i], seqs_b[j]) a separate k-mer match matrix is computed
    and placed in the corresponding cell of the combined grid.  Because every segment
    starts on a bin boundary (offset = cumulative ceil(seg_len/binsize)*binsize), no
    bin ever straddles two segments — eliminating the "mixed bin" artefact.

    Coordinate space: x ∈ [0, Σ_j ceil(len_j/binsize)*binsize],
                      y ∈ [0, Σ_i ceil(len_i/binsize)*binsize].
    Each pixel is exactly binsize × binsize bp.
    """
    if k <= 0 or binsize <= 0:
        raise ValueError("k and binsize must be positive")
    if not seqs_a or not seqs_b:
        raise ValueError("seqs_a and seqs_b must be non-empty")

    enc_table = np.full(256, 255, dtype=np.uint8)
    for _b, _v in zip(b'ACGT', [0, 1, 2, 3]):
        enc_table[_b] = _v

    bins_a = [max(math.ceil(len(s) / binsize), 1) for s in seqs_a]
    bins_b = [max(math.ceil(len(s) / binsize), 1) for s in seqs_b]
    total_bins_a = sum(bins_a)
    total_bins_b = sum(bins_b)

    if total_bins_a * total_bins_b > max_cells:
        raise ValueError(
            f"Mosaic matrix {total_bins_a}×{total_bins_b} = {total_bins_a*total_bins_b} cells; "
            "increase binsize or max_cells"
        )

    max_bins_any = max(max(bins_a), max(bins_b))
    SHIFT = max(int(max_bins_any - 1).bit_length(), 1)
    if 2 * k + SHIFT >= 63:
        raise ValueError(
            f"k={k} too large for binsize={binsize} (SHIFT={SHIFT}); reduce k or increase binsize"
        )
    MASK = (1 << SHIFT) - 1

    def _seg_pk(seq: str, do_rc: bool) -> np.ndarray:
        """Packed sorted (hash << SHIFT | bin) array for one segment sequence."""
        seq = seq.upper()
        if len(seq) < k:
            return np.array([], dtype=np.int64)
        enc = enc_table[np.frombuffer(seq.encode('ascii'), dtype=np.uint8)]
        n_pos = len(enc) - k + 1
        ph = 4 ** np.arange(k - 1, -1, -1, dtype=np.int64)
        pl = 4 ** np.arange(k,           dtype=np.int64) if do_rc else None
        ei = enc.astype(np.int64)
        fwd = np.zeros(n_pos, dtype=np.int64)
        rci = np.zeros(n_pos, dtype=np.int64) if do_rc else None
        for j in range(k):
            sl = ei[j: j + n_pos]
            fwd += sl * ph[j]
            if do_rc:
                rci += sl * pl[j]
        if mask_ambiguous:
            inv_pos = np.where(enc == 255)[0]
            valid = np.ones(n_pos, dtype=bool)
            for p in inv_pos:
                valid[max(0, p - k + 1): min(n_pos, p + 1)] = False
            pos = np.where(valid)[0]
            fwd = fwd[valid]
            rc  = (np.int64(4 ** k - 1) - rci[valid]) if do_rc else None
        else:
            pos = np.arange(n_pos, dtype=np.intp)
            rc  = (np.int64(4 ** k - 1) - rci) if do_rc else None
        seg_bins = (pos // binsize).astype(np.int64)
        if do_rc:
            h_all = np.concatenate([fwd, rc])
            b_all = np.concatenate([seg_bins, seg_bins])
        else:
            h_all = fwd
            b_all = seg_bins
        pk = (h_all << SHIFT) | b_all
        pk.sort()
        # Repeat filter on a-side
        if do_rc and max_kmer_freq > 0 and len(pk) > 0:
            hp = pk >> SHIFT
            re = np.concatenate([np.where(np.diff(hp))[0] + 1, [len(pk)]])
            rs = np.concatenate([[0], re[:-1]])
            rl = re - rs
            vr = rl <= max_kmer_freq
            keep = np.zeros(len(pk) + 1, dtype=np.int8)
            keep[rs[vr]] += 1
            keep[re[vr]] -= 1
            pk = pk[np.cumsum(keep)[:-1].astype(bool)]
        return np.unique(pk)   # collapse duplicate (hash, bin) pairs

    # Hash each segment once
    pks_a = [_seg_pk(s, do_rc=True)  for s in seqs_a]   # fwd + RC
    pks_b = [_seg_pk(s, do_rc=False) for s in seqs_b]   # fwd only

    cum_a = [0] + list(np.cumsum(bins_a))
    cum_b = [0] + list(np.cumsum(bins_b))

    mat = np.zeros((total_bins_a, total_bins_b), dtype=np.uint32)

    _tt = time.perf_counter()
    for i, a_pk in enumerate(pks_a):
        if len(a_pk) == 0:
            continue
        for j, b_pk in enumerate(pks_b):
            if len(b_pk) == 0:
                continue
            b_fwd_s = b_pk >> SHIFT
            bin_b_s = b_pk & MASK
            lo = np.searchsorted(a_pk, b_fwd_s << SHIFT,       side='left')
            hi = np.searchsorted(a_pk, (b_fwd_s + 1) << SHIFT, side='left')
            cnt = hi - lo
            has_match = cnt > 0
            if not has_match.any():
                continue
            lo_m  = lo[has_match]
            cnt_m = cnt[has_match]
            bb_m  = bin_b_s[has_match]
            total = int(cnt_m.sum())
            # Expand a_pk bin indices (ragged-range cumsum)
            grp = np.empty(len(lo_m), dtype=np.int64); grp[0] = 0
            if len(lo_m) > 1: grp[1:] = np.cumsum(cnt_m[:-1])
            diffs = np.ones(total, dtype=np.int64); diffs[0] = lo_m[0]
            if len(lo_m) > 1:
                diffs[grp[1:]] = lo_m[1:] - lo_m[:-1] - cnt_m[:-1] + 1
            ba_rep = a_pk[np.cumsum(diffs)] & MASK
            bb_rep = np.repeat(bb_m, cnt_m)
            flat   = ba_rep * bins_b[j] + bb_rep
            sub    = np.bincount(flat, minlength=bins_a[i] * bins_b[j])
            mat[cum_a[i]:cum_a[i+1], cum_b[j]:cum_b[j+1]] = \
                sub.reshape(bins_a[i], bins_b[j]).astype(np.uint32)
    if verbose:
        print(f"  [dotplot_kmer_mosaic] accumulate {len(pks_a)}x{len(pks_b)} cells: "
              f"{time.perf_counter()-_tt:.2f}s")

    if ax is None:
        _, ax = plt.subplots(figsize=(16, 5))

    plot_mat = np.log1p(mat) if log_scale else mat.astype(float)
    nz = plot_mat[plot_mat > 0]
    vmin = 0.0
    vmax = float(np.percentile(nz, 99.5)) if nz.size > 0 else 1.0
    ax.imshow(plot_mat, origin='lower', aspect='equal', cmap='Greys',
              vmin=vmin, vmax=vmax,
              extent=(0, total_bins_b * binsize, 0, total_bins_a * binsize))
    ax.set_xlabel("Seq B bins")
    ax.set_ylabel("Seq A bins")
    ax.set_title(f"K-mer dotplot mosaic (k={k}, binsize={binsize})")
    return ax
