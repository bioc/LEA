# Frozen .lfmm reader compatibility contract

Baseline: published e9b4998186416322bdb729dfd2b02e4cc074029a, R 4.6.1.
Fixtures and expected.rds were captured before changing R/readLfmm.R. Never
regenerate expected outcomes from the experimental reader to make tests pass.
The 35 fixtures are stored as named raw vectors in fixtures.rds so Git cannot
normalize intentional CRLF/control bytes. Their original per-file SHA-256
entries remain in baseline-manifest.json; serialization changed no fixture
bytes. Tests write raw bytes to temporary .lfmm paths on every platform.

**A — documented numeric records.** One individual per nonblank record and one
locus per column, with uniform widths. Preserve real dosage values as well as
integer genotypes; spaces/tabs, LF/CRLF, blank lines, and absent final newline.
Preserve 9 and -9 as numbers and NA as R missingness. Preserve dimensions,
orientation, V1…Vp column labels, automatic NULL row names, and legacy storage
type: integer when every nonmissing column is inferred integer, double when any
column requires double, logical for wholly missing matrices. Decimal/exponent
spellings such as 1.0/1e0 require double even if numerically integral. Single-row
and single-column matrices retain two dimensions. Keep existing warnings.

**B — historical/nonstandard or ambiguous input.** Preserve read.table behavior
via fallback where necessary: comments, quotes, header/row-name inference,
logical/complex/character tokens, BOM, nonfinite numeric forms and unusual
encoding. Numeric extremes follow the original conversion/rounding semantics.
Bad text tokens historically can produce a character matrix; this audit does
not silently redefine that behavior as numeric input or change downstream
statistical validation. Uncertainty means legacy parsing, not speculative
numeric interpretation. Small files may also use legacy parsing.

**C — malformed structure.** Empty and comment-only inputs and ragged rows that
the published parser rejects must remain rejected. In particular balanced
short/long records must never be flattened/repartitioned. Cases triggering
historical header heuristics stay in B rather than gaining new interpretations.
Extra fields must not be truncated. A fast path validates complete record widths
before scanning and verifies all consumed dimensions; fallback delegates the
same original parser and cannot invent a new accepted numeric shape.

No change to public arguments, extension rules, missing-code semantics, LFMM2
calculations, formula preprocessing or dependencies. Files must remain stable
while being read; changed size/timestamps trigger legacy fallback. No reader can
guarantee an atomic snapshot against undetectable concurrent in-place edits.
