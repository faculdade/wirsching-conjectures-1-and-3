---
built: 2026-08-27
source-published: Wirsching 2003, Discrete Contin. Dyn. Syst. 9(3), 771-787
source-preliminary: the Eichstaett preprint of the same paper
---

# Source concordance: preliminary versus published Wirsching 2003

The paper's footnote in its section on Wirsching's chain says that every
section, equation and result reference to the 2003 paper is to the
published article and was checked against it. This file is that check.

Both versions are now held. The published article is *Discrete and
Continuous Dynamical Systems* **9** (2003), no. 3, 771-787,
doi:10.3934/dcds.2003.9.771, 17 pages. The preliminary version is the
one circulated from Katholische Universitaet Eichstaett, 19 pages,
self-paginated.

Every claim the paper draws from the source is below, with the text as
it appears in each version. Text is quoted from `pdftotext -layout`
output, so ligatures and some subscripts are lost; differences of that
kind are marked and are extraction artifacts, not textual differences.

**Result: identical throughout. No numbering shift, and no difference of
content in any item the paper uses.** The published version is the one
the bibliography now cites, and the page on which each section and named
result begins is printed there.

## Sections

| Section | Title | Published p. | Preliminary |
|---|---|---|---|
| 1 | A first condition for positive density | 773 | present, same title |
| 2 | Generators for Elka functions | 776 | present, same title |
| 3 | Normalization and digital topology | 778 | present, same title |
| 4 | Transition operators | 780 | present, same title |
| 5 | The limiting transition operator S_infty | 780 | present, same title |
| 6 | A strongly stable Markov chain | 783 | present, same title |
| 7 | The asymptotics of a quotient | 784 | present, same title |

## Named results

| Item | Published p. | Identical? | Note |
|---|---|---|---|
| Definition 1 | 772 | yes | "if there is a real constant c > 0 such that liminf ... >= c/a for each a in A". This is the wording the paper's quantifier reading rests on, and it is verbatim the same in both |
| Theorem 1 | 775 | yes | (*1) => uniform positive predecessor density on the non-cyclic a not= 0 mod 3 |
| Conjecture 1 | 778 | yes | (*2) => (*1) |
| Theorem 2 | 780 | yes | (*3) => (*2). Stated at the END of section 3, which is where the paper cites it. Both versions carry the same typo, "and in index l_0" for "an index" |
| Theorem 3 | 781 | yes | strong operator convergence of the transition operators, uniform on bounded equi-continuous families. Cited in this paper's section on the invariant density |
| Theorem 6 | 783 | yes | properties (a), (b), (c) of W_3 |
| Corollary 7 | 783 | yes | existence and uniqueness of phi |
| Corollary 8 | 783 | yes | geometric convergence of the W_3 iterates, same constant 2^{-n+1} |
| Conjecture 2 | 784 | yes | (*4) => (*3). Stated at the END of section 6 |
| Conjecture 3 | 787 | yes | (*5), with the parenthesis "(Note that the results given in section 9 of [1] suggest that conjecture 3 is true.)" present in both |

## Numbered equations the paper cites

| Equation | Published p. | Identical? | What it is |
|---|---|---|---|
| (1.4) | 774 | yes | the 3-adic average of the Elka functions. Extraction artifact only: `pdftotext` drops the slash in "a not= 0 mod 3" in the published file |
| (1.5) | 775 | yes | the class A_delta |
| (2.1) | 777 | yes | the recursion for the generators g_l |
| (2.2) | 777 | yes | the generating function for p_l(m) |
| (2.3) | 777 | yes | e_l = p_l * g_l |
| (2.4) | 778 | yes | the averaged form. Both versions say "each urn U_j has capacity c_j", which is one more than the maximum occupancy of the urn source; the paper already flags this |
| (3.2) | 779 | yes | the class A-tilde_delta |
| (7.3) | 785 | yes | K_l(x,0) = K_l(x + 3^-(l+1), 1) |
| (7.4) | 785 | yes | K_l(x,0) > phi(x) > K_l(x,1) for 0 < x < 1/3 |
| (7.5) | 786 | yes | the sufficient estimate for (*4) |
| (7.7) | 786 | yes | phi'(x) = (9/2) phi(3x) on [0, 2/3] |
| (7.9) | 786 | yes | the untruncated equation |
| (7.10) | 786 | yes | the truncated equation |
| (7.11) | 786 | yes | the closed form of phi_0 |
| (7.12) | 786 | yes | the asymptotic-solution property |
| (7.13) | 786 | yes | the limit 2/3, uniformly on the class |
| (7.14) | 787 | yes | the class lifting |

## How this was checked

`pdftotext -layout` on both files, page by page, then a scripted
comparison that printed, for each of the thirty-four locators above, the
surrounding text from each version side by side. That printed output was
read, and the quotations in the tables above are taken from it. The
extraction and the comparison are reproducible from the commands
recorded in this repository's history. What was not done: no
character-level diff of the two full texts was run, so this certifies
the thirty-four locators the paper actually uses (seven sections,
ten named results, seventeen numbered equations) and nothing beyond
them.

## What changed in the paper as a result

Nothing mathematical. The bibliography's caveat is gone and the entry
cites the published article; a footnote in the paper's section on
Wirsching's chain states that every locator is to that article and
prints the page on which each section and named result begins.
