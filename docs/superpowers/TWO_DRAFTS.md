# Two V3 drafts, for comparison

Written overnight 2026-08-25 into 08-26. Nothing is merged, pushed or deployed. Neither branch
exists on the remote; `git ls-remote --heads origin` lists five branches and neither of these
is among them.

---

## Read this first

**Do not run `python3 GPCompaReports_v2/generate_site.py` without an `--output` argument.**

It writes into `GPCompaReports_v2/output/`, which holds the pre-branch build that the freeze
test compares against. Overwriting it makes the test compare the branch to itself, so it passes
forever and stops protecting the PI-approved report pages. `bash scripts/build_demo.sh` is safe:
it writes only to a draft directory and reads `output/data/` without writing there.

---

## What to open

Both builds are already on disk. You do not need to build anything to look.

| | primary draft | alternative draft |
|---|---|---|
| branch | `feature/v3-demo` (53 commits) | `feature/v3-demo-alt` (55 commits) |
| landing page | `GPCompaReports_v2/output_v3_demo/index.html` | `GPCompaReports_v2/output_v3_alt/index.html` |
| screenshots | `docs/superpowers/screenshots/` | `docs/superpowers/screenshots-alt/` |
| full handoff | `docs/superpowers/V3_DRAFT_HANDOFF.md` | same, plus `.superpowers/sdd/alt-draft-report.md` |

The alternative branches **from** the primary, so it carries every one of the primary's 53
commits. **The only difference between the two drafts is the landing page.** Report pages, the
Contact page, the Downloads page, the statistics page, the tests and all other code are identical.

---

## How they differ

### Primary: linear document, eight sections

Hero with a factual paragraph, then a dataset table, then the material in teaching order.

1. Hero: name, what the database is, two calls to action, search
2. Dataset and analysis: the frozen counts, then the limitation statement
3. How to read a GPCompaRe contact: the real PAR2 worked example
4. Core Functional Residues across Class A GPCRs: the manuscript's topology figure
5. What each receptor report contains
6. Analysing your own structures
7. Methods and data sources
8. Version, downloads and citation

Reads top to bottom. A reader who scrolls once has met the concept, the caveats and the sources.

### Alternative: receptor-first console, six sections

The search field is the primary object. Beside it sits one complete real receptor record,
rendered at first paint, so a visitor with no receptor in mind learns the unit of content by
reading one instead of reading a description of one. A ligand-class facet row gives a second
way in.

1. Lookup console: search, ten ligand-class chips with live counts, and the PAR2 record card
2. Dataset and analysis
3. How to read one contact
4. Inside a receptor report
5. Recurrent positions across the family
6. Reference: methods, sources, software, release and citation

Optimised for the visitor who arrived wanting one receptor. The explanatory material sits
around the search rather than ahead of it.

### The judged comparison

Four directions were proposed and scored by three lenses: a sceptical journal reviewer, a bench
scientist looking up one receptor, and the lab member who must maintain this for five years
without a front-end toolchain.

| direction | mean score |
|---|---|
| receptor-first console (built as the alternative) | 7.0 |
| data-repository index | 6.7 |
| methods document with a table of contents | 6.3 |
| lead with the scientific finding | 5.7 |

The two rejected directions are described in `.superpowers/sdd/` if you want to see them.

---

## What is true of both

Every one of these was verified against the data, not asserted.

- **The freeze holds.** Approved report sections are byte-identical to the pre-branch build
  across all five demo receptors, and the injected section nav and toolbar were confirmed at
  runtime to land outside every approved section.
- **Every frozen count reproduces** from the batch CSVs: 283 receptors, 60 families, 566 models,
  213,456 contact-pair records, 23,025 threshold-passing changes, 4 excluded, 279 with variant data.
- **The PAR2 worked example is exact**: M159 `3.36x36` against F300 `6.48x48`, active 0.67,
  inactive 5.66, ΔRRCS −4.99, threshold 2.60, 27th of 776. Matches Figure 4d of the manuscript.
- **The analysis layer from the first round is entirely gone**, not hidden.
- 20 tests pass on both branches.

---

## Claims removed during review, so you know what was nearly shipped

Each of these rendered to a reader before being caught:

- **"M159 is the methionine of the xMY motif."** Your own fact-check rejects that label: the
  adjacent Tyr3.37 fails magnitude and presence under two aggregations, and the supported label
  is M3.36, itself a flagged reconstruction. The copy now credits only the canonical aromatic
  6.48 toggle-switch position.
- **"CFR top 50"** badges on residues ranked as low as 356, because the rank map was built from
  the whole table rather than the published subset. On PAR2 that was 37 of 100 rows.
- **A hero claiming every report links gnomAD variation.** Four receptors have none. The same
  claim had already been scoped correctly further down the same page.
- **"Receptor-level Core Functional Residues"** listed as report content, when what a report
  shows is the cross-family set mapped onto that receptor.
- **A statistics tooltip** stating a CFR score formula the code does not compute.
- Six further unsupported claims in the analysis layer that was subsequently deleted outright.

---

## Decisions waiting for you

In rough priority order.

1. **Which landing direction**, or a hybrid. They are independent branches, so keeping both a
   while costs nothing.
2. **The CFR cutoff inconsistency, which is a scientific call and was deliberately not made.**
   The statistics page's chart and table use the published top 30, but the CFR contact-pair
   network and the variant-enrichment analysis still use the top 50. So the table headed
   "Top CFR Contact Pairs (both residues are CFRs)" includes residues outside the top 30.
   Aligning them on 30 would change the published enrichment percentages and the chi-squared
   result, which may be quoted in the manuscript. That is why it was left alone.
3. **Pin the freeze baseline** before this stops being a draft. Right now it is honest only
   because the file on disk happens to predate the branch.
4. **The "Complete RRCS results" heading**, capped at 1,000 rows and therefore incomplete for
   7 of 283 receptors. Proposed correction, pending your PI's approval and deliberately not
   applied: "RRCS results, up to 1,000 contact pairs".
5. **Four receptors absent from the variant dataset. Traced 2026-08-26, and it is two
   different causes, not one.** `npy6r_human` and `taar3_human` resolve to annotated pseudogene
   loci (`NPY6RP`, `TAAR3P`): gnomAD carries variants there but **zero missense**, so nothing can
   be placed on a residue and zero is correct. `gp182_human` and `p2ry8_human` are a
   data-collection limitation: gnomAD does hold missense variants for both (953 and 526), and
   they are absent because their identifiers were not resolved when the dataset was collected.
   P2RY8 is protein coding and sits in the pseudoautosomal region, where UniProt lists a second
   Ensembl ID that gnomAD does not recognise. **Known release limitation, not an open website
   bug**: the site keeps the submitted 279-receptor dataset and states both causes separately on
   the statistics page. Repairing the two would change frozen report pages, the 279 count and
   Figure 3's denominators, so it is recorded for the later manuscript revision.
6. **The database download must become a real archive before publication.** The software
   placeholder may persist.
7. **External setup nobody can do but you**: Formspree account routed to `tggrlab@gmail.com`,
   reCAPTCHA keys, domain restriction, named contacts, licences, Zenodo deposit. Listed in
   `docs/superpowers/EXTERNAL_SETUP.md`.

---

## Two findings that reach beyond the website

**Generic-number coverage is much lower than the figure in circulation.** Two different
coverage measures had been conflated. Verified across all 283 receptors:

| coverage of above-threshold contact positions | min | median | max |
|---|---|---|---|
| assigned to a named segment | 45.8% | 83.9% | 100% (20 receptors) |
| carrying a GPCRdb generic number | 8.2% | 63.8% | 87.7% (none reach 100%) |

The CFR analysis depends on the second row, because mapping a receptor-level CFR to a generic
position requires `display_number`. So CFR recurrence counts, the Top 30 table, the statistics
page's frequency column and Figure 2a's shading are all bounded by roughly two thirds coverage,
not five sixths. Both landing pages now state that recurrence counts are lower bounds. **Whether
the manuscript needs the same framing is your call. Nothing in the manuscript was touched.**
Note also that a figure quoted as 45.9% is 45.8% on recomputation.

**Y160 is an annotation-coverage gap, not a numbering conflict. Closed 2026-08-26.** Y160 is
`3.37x37` in the manuscript, in the full GPCRdb residue map and in the report snake plot. All
three agree. The original sparse per-report annotation CSV has no row for residue 160, so
report-level contact data may not display that generic number. This is the same annotation-coverage
limitation that makes the site's own recurrence count 356 against the submitted 368, not a
disagreement about the position. It does not affect the selected M159-F300 explainer. No report or
annotation change was made, because the report content is frozen.

---

## One process gotcha worth knowing

Draft build directories are gitignored, so they persist across branch checkouts. While the
alternative was being built, `output_v3_demo/` from the primary was still on disk, and tests
that read it ran against the *other* branch's output. Nothing was harmed, and the alternative's
own build was verified separately by hand, but a test that reads a gitignored directory can
silently validate the wrong branch.
