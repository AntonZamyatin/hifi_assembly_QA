"""Run the pair HMM forwards: draw reads from a template (stage D of docs/phmm.md).

The model is generative, so it can be walked in either direction. `dp.py` asks
"how probable is this read"; this module asks "what read would this process
produce". Having both is what makes the consistency test possible: sample many
reads and the observed frequency of each must equal what the scorer assigns it.
Nothing else in the project can check the scorer against something that is not
itself a rearrangement of the same recursion.

Two ways to stop, and the choice decides what is being sampled.

`tau > 0` -- the chain terminates on its own, from M or I. This is the joint
model over (read, length), and the frequency of a read is exactly
`tau * exp(forward(...).logp)`.

`max_len` with `tau = 0` -- the chain never terminates, so the walk stops after
that many read bases. This is the length-conditioned model of §5.6, the one the
project actually uses, and the frequency of a read is *proportional* to
`exp(forward(...).logp)` with `1/Z(l)` as the constant. That proportionality is
testable without knowing Z, and inverting it estimates Z empirically.

Qualities are conditioned on, never generated: the model has nothing to say about
where a Phred score comes from, so one is supplied and used for both the error
rate and the returned array.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from qamodel.phmm.params import PhmmParams, eps_from_phred

_BASES = np.frombuffer(b"ACGT", dtype=np.uint8)
_OTHER = {b: np.array([x for x in _BASES if x != b], dtype=np.uint8)
          for b in _BASES}


@dataclass(frozen=True)
class Sample:
    """One read, with the truth that produced it.

    `path` is the state string, so `path.count('M')` matched columns and so on.
    `end` is the last template index consumed, `origin - 1` if the walk consumed
    none. `truncated` means the walk wanted template the window did not have --
    the sample is then not a draw from the model and must be discarded, not used.
    """

    read: str
    qual: np.ndarray
    path: str
    origin: int
    end: int
    truncated: bool


def sample_read(template, params: PhmmParams, rng: np.random.Generator, *,
                origin: int = 0, qual: int = 40,
                max_len: int | None = None) -> Sample:
    """Draw one read from `template` starting at `origin`.

    Stops when the chain enters the end state (needs `tau > 0`) or when
    `max_len` read bases have been emitted, whichever happens first. With
    `tau = 0` and no `max_len` the walk cannot terminate, so that is rejected.
    """
    if params.tau <= 0.0 and max_len is None:
        raise ValueError("tau = 0 never terminates; pass max_len to stop the walk")

    t = np.frombuffer(template.upper().encode("ascii"), dtype=np.uint8) \
        if isinstance(template, str) else np.asarray(template, dtype=np.uint8)
    T = t.shape[0]
    if not 0 <= origin < T:
        raise ValueError(f"origin {origin} outside template of length {T}")

    eps = float(eps_from_phred(qual))
    d_i, d_d = params.delta_i, params.delta_d
    g_i, g_d, tau = params.gamma_i, params.gamma_d, params.tau

    read, path = [], []
    v = origin
    truncated = False

    state = "I" if rng.random() < d_i else "M"          # B -> I or B -> M

    while True:
        if state == "M":
            if v >= T:
                truncated = True
                break
            base = t[v]
            if rng.random() < eps:
                base = _OTHER[base][rng.integers(3)]     # a specific wrong base
            read.append(base)
            path.append("M")
            v += 1
        elif state == "I":
            read.append(_BASES[rng.integers(4)])
            path.append("I")
        else:                                            # D emits nothing
            if v >= T:
                truncated = True
                break
            path.append("D")
            v += 1

        # The stop test sits here on purpose: it can only fire just after an
        # emission, so a walk never ends in D -- which is the D -> E exclusion
        # showing up in the sampler exactly as it does in the recursion.
        if max_len is not None and len(read) >= max_len:
            break

        r = rng.random()
        if state == "M":
            if r < tau:
                break
            elif r < tau + d_i:
                state = "I"
            elif r < tau + d_i + d_d:
                state = "D"
            else:
                state = "M"
        elif state == "I":
            if r < tau:
                break
            elif r < tau + g_i:
                state = "I"
            else:
                state = "M"
        else:
            state = "D" if r < g_d else "M"              # D cannot terminate

    return Sample(
        read=np.array(read, dtype=np.uint8).tobytes().decode("ascii") if read else "",
        qual=np.full(len(read), qual, dtype=np.int64),
        path="".join(path),
        origin=origin,
        end=v - 1,
        truncated=truncated,
    )


def sample_reads(template, params: PhmmParams, rng: np.random.Generator,
                 n: int, **kw) -> list[Sample]:
    """`n` independent draws. Truncated ones are returned, not filtered."""
    return [sample_read(template, params, rng, **kw) for _ in range(n)]


def run_lengths(path: str, state: str) -> list[int]:
    """Lengths of the maximal runs of `state` in a path, for checking gap models."""
    out, run = [], 0
    for c in path:
        if c == state:
            run += 1
        elif run:
            out.append(run)
            run = 0
    if run:
        out.append(run)
    return out
