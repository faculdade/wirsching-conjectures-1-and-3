# Reproducibility material: Wirsching's Conjectures 1 and 3

Code for R. A. Tavares, *Wirsching's Positive-Predecessor-Density
Program: Proofs of Conjectures 1 and 3*.

Archived on Zenodo. **Cite the version DOI**, not the concept one:

| | DOI | resolves to |
|---|---|---|
| `v1.0.2` | minted at release, recorded in the paper | this exact snapshot |
| `v1.0.1` | [10.5281/zenodo.22116346](https://doi.org/10.5281/zenodo.22116346) | commit `c0339d1` |
| `v1.0.0` | [10.5281/zenodo.22116226](https://doi.org/10.5281/zenodo.22116226) | commit `40e378f` |
| all versions | [10.5281/zenodo.22116225](https://doi.org/10.5281/zenodo.22116225) | whatever the latest release is |

The paper cites the version DOI, because a reader reproducing a printed
number needs the code that produced it, not the newest code. Each new
release mints a new version DOI and leaves the older archives untouched;
the concept DOI never changes.

`paper/` carries the paper as it stood at the release that archived it.
It can lag the paper of record by one revision, and does whenever the
paper is corrected after a release: adding this DOI to the paper was
itself such a correction. The paper of record is the preprint.

The paper proves two of the three conjectures G. J. Wirsching left open
in *On the problem of positive predecessor density in 3n+1 dynamics*,
Discrete Contin. Dyn. Syst. 9(3) (2003), 771-787,
[doi:10.3934/dcds.2003.9.771](https://doi.org/10.3934/dcds.2003.9.771).
Conjecture 1 and Conjecture 3 are proved, condition `(*4)` follows at
every window radius with `mu = 1/3`, and Wirsching's chain reduces to
the single condition `(*3)`. Conjecture 2 is his route to `(*3)` and
remains open.

## Layout

One folder per section of the paper, numbered as the paper numbers them.
Sections 1, 2, 3 and 5 carry no computation.

| Folder | Paper section | Backs |
|---|---|---|
| `section-04-conjecture-1/` | 4, Wirsching's Conjecture 1 | the generating-function cancellation, the partition bound (5), and the explicit constants of Lemmas 4.1 and 4.2 |
| `section-06-periodic-correction/` | 6, the periodic correction | the certified enclosure of `osc(H)` and the derivative bounds (10) |
| `section-07-saddlepoint/` | 7, the uniform saddlepoint approximation | the constant chain behind Theorem 7.5 and the threshold `N >= 19` |
| `section-08-envelope-lemma/` | 8, the envelope lemma | Lemma 8.1 and Proposition 8.2 |
| `section-09-berg-kruppel-identity/` | 9, the Berg-Kruppel identity | Proposition 9.1 and the closed-form constant `C_P` |
| `section-10-conjecture-3/` | 10, Conjecture 3 | numerical evidence for Theorem 10.1 along Wirsching's own sequence |
| `section-11-condition-star4/` | 11, the chain closed to a single condition | the two steps the source states rather than carries out, `mu = 1/3`, and the residual in Remark 11.3 |
| `section-12-numerics/` | 12, Discussion | the derivation of `e^{H(0)}` and its cross-check against section 6's certified enclosure; also the depth-500 evaluation of `phi`, which preceded the theorem and which the current paper no longer quotes (see that folder's README) |
| `paper/` | | the paper itself, source and PDF |

Each folder has its own README saying what its scripts verify, how to
run them, and what to expect.

## What is certified and what is not

Only two kinds of claim in the paper are computer-assisted, and both
live in `section-06-periodic-correction/`: the enclosure of `osc(H)` and
the bounds on `H'` and `H''`. Those use Arb ball arithmetic through
python-flint, at a working precision stated at the top of the script,
and the script prints its input balls and final enclosure, so the
printed run output is the certificate.

Non-constancy of `H` is **not** among them. It is proved analytically,
from the classical zero-free theorem for zeta on `Re s = 1`, with no
computation.

Everything else here is high-precision floating arithmetic through
mpmath, or exact integer and rational arithmetic. Where a script claims
a rigorous inequality it asserts it and exits non-zero on failure.

## Running

```
python3 -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
```

Then, in any folder, `python3 <script>.py`. Or all of them:

```
set -e
for d in section-*/; do
  echo "=== $d ==="
  ( cd "$d" && for f in *.py; do echo "--- $f"; python3 "$f" || exit 1; done ) || exit 1
done
```

Every script is self-contained: no imports from elsewhere in this
repository, no hardcoded paths, no network. `section-12-numerics/experiment_conjecture3.py`
is the long one; the rest are minutes or less.

## Relation to the two earlier repositories

Two repositories preceded this one and remain where they are, with their
DOIs valid:

- `faculdade/wirsching-conjecture3-proof`, archived at
  [doi:10.5281/zenodo.21854549](https://doi.org/10.5281/zenodo.21854549),
  commit `f8248c3`, for the standalone Conjecture 3 preprint. Its
  scripts are here under sections 6 to 10, renumbered to this paper's
  sections. In that archive the results are numbered as they were first
  published: Theorem 10.1 is Theorem 1, Theorem 10.2 is Theorem 2,
  Corollary 10.3 is Corollary 3, Theorem 6.1 is Theorem 4,
  Proposition 6.3 is Proposition 6, Proposition 6.5 is Proposition 8,
  Theorem 7.5 is Theorem 13, Lemma 8.1 is Lemma 16, Proposition 8.2 is
  Proposition 17, Proposition 9.1 is Proposition 18, and
  Proposition 10.6 is Proposition 20.
- `faculdade/collatz-wirsching-2003`, for the unsubmitted manuscript
  that carried Conjecture 1 together with material on Conjecture 2. Its
  Conjecture 1 scripts are here under section 4 and its evaluation of
  `phi` under section 12. Its microcanonical experiments are not here:
  that material is not in this paper.

This repository is the one the paper cites.
