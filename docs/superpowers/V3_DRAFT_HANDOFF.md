# GPCompaRe V3 review draft: handoff

## STOP: a full site build destroys this branch's safeguard

The freeze test, `GPCompaReports_v2/tests/test_freeze.py`, is the safeguard this
whole draft exists to protect: it proves the approved report pages have not
changed. It works by comparing a fresh build against
`GPCompaReports_v2/output/reports/par2_human.html`, the last full build the PI
approved. That directory is gitignored (`.gitignore:17`), so it is not version
controlled, and its contents currently predate this branch.

**Running `python3 GPCompaReports_v2/generate_site.py` with no `--output`
argument writes into `GPCompaReports_v2/output/` and overwrites that
baseline.** `--output` defaults to exactly that directory
(`generate_site.py:25-26`), the same one the freeze test reads from. Once it is
overwritten, the freeze test is comparing this branch's build against itself.
It will pass forever after that, and the property this entire draft protects
stops being checked, silently, with no error anywhere.

**Do this instead.** `bash scripts/build_demo.sh` is safe: it writes only to
`GPCompaReports_v2/output_v3_demo/`, and the only thing it does with
`GPCompaReports_v2/output/` is read from `output/data/`, to seed the
conservation and AlphaMissense caches. It never writes into `output/`.

The final whole-branch review judged this effectively blocking, because
rebuilding the full site is a plausible first action for anyone opening this
draft.

---

Branch `feature/v3-demo`. This document is the record for the project owner
of what this draft contains, what was verified and how, and what is not
resolved yet. It is a draft review artifact, not a release.

Plan: `docs/superpowers/plans/2026-08-25-v3-review-draft.md`.
Running ledger kept during the work: `.superpowers/sdd/progress.md`.

---

## 0. Read this first: this is a five-receptor demo, not the full site

**Scope.** `scripts/build_demo.sh` builds only five receptors
(`par2_human`, `par1_human`, `adrb2_human`, `5ht2a_human`, `cxcr4_human`),
so `GPCompaReports_v2/output_v3_demo/reports/` contains only those five
report pages. The landing hero, the browse page and the search index are
unchanged production content: they still list and index all 283
receptors. **278 of the 283 receptor links in this draft do not
resolve.** `scripts/build_demo.sh` already prints this warning at the end
of its own output; this handoff did not previously say so, and it is the
most likely way a reader misjudges the scope of what they are looking at.

**How to open it.** `output_v3_demo/` is gitignored, so a fresh clone has
nothing to open until it is built:

```bash
bash scripts/build_demo.sh
```

Then open directly as local files, no server needed:

- `GPCompaReports_v2/output_v3_demo/index.html` (landing page)
- `GPCompaReports_v2/output_v3_demo/reports/par2_human.html`
- `GPCompaReports_v2/output_v3_demo/reports/par1_human.html`
- `GPCompaReports_v2/output_v3_demo/reports/adrb2_human.html`
- `GPCompaReports_v2/output_v3_demo/reports/5ht2a_human.html`
- `GPCompaReports_v2/output_v3_demo/reports/cxcr4_human.html`
- `GPCompaReports_v2/output_v3_demo/contact.html`
- `GPCompaReports_v2/output_v3_demo/downloads.html`

## 1. Nothing merged, deployed or pushed

Verified directly before writing this document, not assumed:

- `git status --short` on `feature/v3-demo` shows no uncommitted changes
  under `GPCompaReports_v2/`, `docs/`, or any tracked path touched by this
  work. (A handful of pre-existing untracked files at the repo root, e.g.
  `PROJECT_SUMMARY.md`, `all_receptors.txt`, are unrelated leftovers from
  other sessions and were not touched.)
- `git log origin/master..feature/v3-demo` and `git ls-remote --heads
  origin` (network reachable, confirmed) show `feature/v3-demo` is not
  among the remote branches. It has never been pushed.
- `feature/v3-demo` has not been merged into `master`,
  `feature/site-updates-2026-04-30` (the live site source per the
  project's `CLAUDE.md`), or any other branch.
- The `gh-pages` branch, both locally and at `origin/gh-pages`, is still at
  its 2026-05-07 deploy (`site-v2-released-47-g3a10a53`), which is the exact
  commit this branch's own base (`3a10a53`) matches. Nothing from this draft
  has been deployed.

## 2. The selected PAR2 contact, with its full verification table

The worked example on the landing page is the PAR2 (F2RL1) M159 to F300
contact. All figures below were re-derived independently for this handoff
from the underlying batch CSV,
`The_batch_RRCS_analyzer/batch_analysis_full/batch_analysis_20260202_151051/csv_data/par2_human_rrcs_delta.csv`
(776 rows for PAR2), and cross-checked against the built
`GPCompaReports_v2/output_v3_demo/reports/par2_human.html`.

| Field | Value | How it was checked |
|---|---|---|
| Residue 1 | M159, generic number `3.36x36`, TM3 | matches `par2_human.html`'s own "Largest changes" table row |
| Residue 2 | F300, generic number `6.48x48`, TM6 | same row |
| Active RRCS | 0.6651 (shown rounded as 0.67) | CSV `active_rrcs` for res1=159, res2=300; report shows 0.665 |
| Inactive RRCS | 5.6580 (shown rounded as 5.66) | CSV `inactive_rrcs`; report shows 5.658 |
| ΔRRCS | −4.9928 (shown rounded as −4.99) | CSV `delta_rrcs` = active − inactive; report shows −4.993 |
| PAR2 threshold | 2.6016 (shown rounded as 2.60) | recomputed as `max(mean(abs_delta) + std(abs_delta), 0.2)` over all 776 PAR2 rows |
| Contact passes threshold | Yes (4.9928 > 2.6016) | direct comparison |
| Rank in PAR2 | 27 of 776, by \|ΔRRCS\| | the 776-row CSV is already sorted descending by `abs_delta`; the M159/F300 row is at 0-based index 26 |

A second row from the same CSV matters for item 10 below: M159 to Y160
(res1=159, res2=160) has ΔRRCS = −1.7826 (\|ΔRRCS\| = 1.78), which is below
the 2.6016 threshold and does not pass. That is the numeric basis for the
manuscript's Figure 4e and for the Y160 concern recorded in item 10.

## 3. Every scientific count on the landing page, and its exact source

All counts were confirmed present, verbatim, in the built
`GPCompaReports_v2/output_v3_demo/index.html` and are pinned by
`GPCompaReports_v2/tests/test_landing_counts.py`.

| Count as shown | Where it appears | Source |
|---|---|---|
| 283 | Hero tagline, "Browse 283 GPCRs", search placeholder, dataset table "Receptors" row | `landing_page.py`: `total_gpcrs = len(store.gpcr_ids)` |
| 60 | Dataset table, "Receptor families" | `landing_page.py`: `info_df['receptor_family'].nunique()` |
| 566 | Dataset table, "Structural models" | `landing_page.py`: `n_models = 2 * len(store.gpcr_ids)` |
| 213,456 | Dataset table, "Contact-pair records" | `landing_page.py`: sum of `len(delta_df)` across all 283 receptors' delta tables |
| 23,025 | Dataset table, "Threshold-passing changes" | `landing_page.py`: sum, per receptor, of rows where `abs_delta >=` that receptor's own `_calc_significance_threshold` (`gpcr_report_helpers.py`) |
| 4 | Dataset table, "Excluded"; Methods section names them (GPR26, LGR5, OPSX, RXFP1) | `landing_page.py`: `len(store.metadata) - len(store.gpcr_ids)` |
| 0.67 / 5.66 / −4.99 / 2.60 / 27 of 776 | PAR2 worked-example table | hand-authored prose in `landing.html`, sourced from `par2_human.html`'s own rendered table; independently re-derived for item 2 above |
| 30 | "See the ranked Top 30 recurrent positions" and the CFR topology figure caption | `cfr_analysis.py`, dot-plot now uses `.head(30)` (fixed in commit `52afa57` to match the pre-existing "Top 30" heading on the statistics page) |
| 58% / 9% to 88% | "generic numbering reaches a median of 58% of the distinct above-threshold contact positions in a receptor, ranging from 9% to 88%" | hand-authored prose in `landing.html`; added in commit `6b6fcaa` as an occurrence-weighted 64% (range 8% to 88%), then corrected in commit `608362c` to the more conservative position-based 58% (range 9% to 88%); the underlying number is still not independently re-derived by any code in this branch, see item 12 |
| 1,000 | "up to 1,000 rows" in the report-contents list | hardcoded `n=1000` default in `_get_complete_rrcs()`, `gpcr_report_helpers.py` |
| 279 | "for the 279 receptors that have variant data" | derived count: 283 `*_rrcs_delta.csv` files minus 4 receptors with no matching `*_variants.csv` (see item 10); hardcoded prose in `landing.html` |
| 2 | "two AlphaFold-Multistate models per receptor" | restates `n_models = 2 * total_gpcrs`; constant in the prose |

## 4. Complete user-facing copy added or changed

Everything below is new copy introduced on this branch. No copy inside an
existing approved report section was changed (see item 8).

### Landing page (`GPCompaReports_v2/templates/landing.html`)

**Hero**
> Paired active and inactive AlphaFold-Multistate models for 283 human
> Class A GPCRs, compared using residue-residue contact scores (RRCS). Each
> report links the structural result to GPCRdb generic numbering,
> conservation, AlphaMissense predictions, cross-family Core Functional
> Residue positions, and gnomAD variation for the 279 receptors that have
> it.

The hero previously claimed every report links gnomAD variation. Four
receptors (`gp182_human`, `npy6r_human`, `p2ry8_human`, `taar3_human`) have
no variants section at all (item 10). The correct scoping to 279 receptors
already existed further down the same page; the hero was never carried up
to match it. Fixed in commit `608362c`, which also credits cross-family CFR
positions here instead of the unqualified "Core Functional Residues" the
hero previously implied every report computed on its own.

Buttons: "Browse 283 GPCRs", "Database statistics". Search placeholder:
"Search 283 GPCRs. Try ADRB2, 5HT2A, CXCR4..."

**Dataset and analysis** table rows: Receptors, Receptor families,
Structural models, Contact-pair records, Threshold-passing changes,
Excluded, with the notes shown in item 3's table above, plus the limitation
statement:
> Positions and contact pairs highlighted by GPCompaRe should be
> interpreted as model-derived structural candidates, not direct
> experimental evidence of function.

**How to read a GPCompaRe contact** (the PAR2 worked example): full prose
explaining RRCS as a contact score computed in both models, ΔRRCS defined
as active minus inactive, the verification table from item 2, a statement
that this contact is stronger in the inactive model (blue) versus active
(red), a note that F300 occupies the aromatic 6.48 toggle-switch position
(one of the canonical activation landmarks) and that both residues sit in
transmembrane helices, and a link to the PAR2 snake plot. This replaces an
earlier sentence describing M159 as "the methionine of the xMY motif,"
removed in commit `2445658`; see item 12 for why that label did not hold up.

**Core Functional Residues across Class A GPCRs**: defines a receptor-level
CFR, explains recurrence across at least three receptors, states that
recurrence counts are lower bounds because generic numbering reaches only
a median 58% (range 9% to 88%) of a receptor's distinct above-threshold
contact positions, states that no motif was picked in advance and the
analysis recovers DRY/NPxxY/CWxP machinery on its own, and links to the
ranked Top 30 table. The lower-bounds caveat was added in commit `6b6fcaa`
(originally an occurrence-weighted 64%, range 8% to 88%) and corrected in
commit `608362c` to the position-based 58% (range 9% to 88%); see item 12.

**What each receptor report contains**: a nine-item list (interactive snake
plot; ranked contact pairs up to 1,000 rows; receptor threshold;
cross-family CFR positions mapped onto the receptor; generic numbering;
conservation; gnomAD variants for the 279 receptors that have them;
AlphaMissense; CSV export), plus example-report links to ADRB2 and PAR2.
The CFR bullet previously read "Receptor-level Core Functional Residues,"
which named the wrong set: a report shows the family-scope CFR list mapped
onto that one receptor, not something computed at the receptor level.
Fixed in commit `608362c` to "Cross-family Core Functional Residue
positions mapped onto this receptor."

**Analysing your own structures**: states that GPCompaRe covers analysis
of user-supplied paired active and inactive structures using the same
RRCS comparison, that a releasable version of that program is not yet
available, and that local analyses stay local and are not added to the
public database, with a link to Downloads. This replaces an earlier
sentence saying GPCompaRe "also includes" the program, which contradicted
the Downloads page's own statement that no releasable version exists;
resolved in commit `38c5d1b`, see item 12.

**Methods and data sources**: attributes structures to AlphaFold-Multistate,
scoring to RRCS, generic numbering to GPCRdb, variation to gnomAD v4,
pathogenicity to AlphaMissense, conservation to ProtVar and UniProt; states
the per-receptor threshold formula; states the 566-model, 4-excluded count
with the four receptor names.

**Version, downloads and citation**: "GPCompaRe database release 2026.08",
"DOI pending release", links to Downloads and Contact, and a statement that
corrections and substantial updates will be documented as new releases with
no fixed schedule.

### Contact page (`GPCompaReports_v2/templates/contact.html`, new page)

Heading "Contact", subtitle inviting data/annotation issues, website
problems, receptor questions, or custom-analysis enquiries. Form fields:
Name, Email, Institution, Type of enquiry (six options), Receptor or gene,
Relevant page URL, Message. Two standing notes: "GPCompaRe is an academic
resource maintained on a best-effort basis." and "Do not submit
patient-identifiable or other sensitive information through this form."
While the Formspree endpoint is a placeholder, a visible warning reads
"This form is not connected yet. See the setup notes before publishing this
page." and the Send button is disabled. Lab panel: TGGR Laboratory /
Translational Genetics and Genomics Research Lab / Genetics Institute, Tel
Aviv Sourasky Medical Center, `tggrlab@gmail.com`, and a placeholder
sentence: "Named scientific and technical contacts will be listed here once
assigned."


> **Superseded 2026-08-26.** The form is now live: it posts as plain HTML to `https://formspree.io/f/mljrqazl`, the submit button and all visible fields ship enabled, and the "not connected" notice no longer renders. What remains is dashboard-only and listed under "Contact form: manual Formspree dashboard checks" in `RELEASE_CHECKLIST.md`.

### Downloads page (`GPCompaReports_v2/templates/downloads.html`, new page)

Heading "Downloads", subtitle stating neither archive is packaged yet.
Release metadata: "GPCompaRe database release 2026.08", "DOI pending
release". Two panels, each listing planned archive contents:

- **Database release** (status "Not yet available"): lists the canonical
  receptor-level files, receptor summary, cross-receptor results, README,
  data dictionary, version/date, changelog, provenance, licenses, file
  manifest, checksums. Note: "Packaging is pending a licensing review of
  the third-party annotations." Callout: "This placeholder is provisional.
  The public release of this database will include a complete data
  archive, not this placeholder."
- **Analysis software** ("GPCompaRe software v1.0.0", status "Not yet
  available"): lists the program, install instructions, supported Python
  version, dependencies, one redistributable example, an example command,
  expected outputs, known limitations, software license, validation
  command. Note: "This placeholder may remain in place after publication.
  No fixed schedule is set for releasing a user-facing analysis program."

Footer note: "The website generator is available in the repository for
reproducibility. It is not the analysis program."

## 5. Files changed

Derived directly from `git diff --stat 3a10a53..HEAD`, where `3a10a53` is
this branch's actual base (the merge-base with
`feature/site-updates-2026-04-30`, the live site branch, and the same
commit the brief's own verification command uses). Regenerated for this
fix wave: the table below no longer lists
`GPCompaReports_v2/scripts/build_par_dossier_csvs.py`, `poster_figures.py`
and `preview_par2.py` (665 lines combined). Commit `03c94a3` untracked all
three after an earlier over-broad `git add -A GPCompaReports_v2/` had
swept them into history by mistake. They remain on disk unchanged; only
their git tracking changed. Also regenerated to add
`GPCompaReports_v2/templates/statistics.html`, touched for the first time
by commits `608362c` and `cc8ae9b` (the CFR tooltip fix and the lower-bounds
callout, item 12). File count moved from 39 to 40, insertions from 4,678 to
4,697, deletions from 431 to 436. Verified that the branch now tracks
nothing beyond its own work:

```
$ comm -13 <(git ls-tree -r --name-only 3a10a53 | sort) <(git ls-files | sort)
```

returns exactly this branch's own added and changed files, none of the
three scripts among them.

```
$ git diff --stat 3a10a53..HEAD -- . ':!docs/superpowers/V3_DRAFT_HANDOFF.md' ':!docs/superpowers/screenshots'
 .gitignore                                         |    3 +
 GPCompaReports_v2/analysis/cfr_analysis.py         |   14 +-
 GPCompaReports_v2/generate_site.py                 |    3 +
 GPCompaReports_v2/static/css/landing.css           |  281 +---
 GPCompaReports_v2/static/css/site.css              |  139 +-
 GPCompaReports_v2/static/css/v3.css                |   55 +
 GPCompaReports_v2/static/img/fig-cfr-topology.png  |  Bin 0 -> 406319 bytes
 .../static/img/par2-explainer-snake.png            |  Bin 0 -> 257868 bytes
 GPCompaReports_v2/static/js/landing.js             |   38 +-
 GPCompaReports_v2/static/js/v3-deeplink.js         |  100 ++
 GPCompaReports_v2/static/js/v3-nav.js              |  117 ++
 GPCompaReports_v2/templates/_partials/footer.html  |    2 +-
 GPCompaReports_v2/templates/_partials/topbar.html  |   12 +-
 GPCompaReports_v2/templates/base.html              |    7 +-
 GPCompaReports_v2/templates/contact.html           |  102 ++
 GPCompaReports_v2/templates/downloads.html         |   90 ++
 GPCompaReports_v2/templates/gpcr_report.html       |   56 +-
 GPCompaReports_v2/templates/landing.html           |  304 ++--
 GPCompaReports_v2/templates/statistics.html        |   12 +-
 GPCompaReports_v2/tests/__init__.py                |    0
 GPCompaReports_v2/tests/test_brand.py              |   23 +
 GPCompaReports_v2/tests/test_freeze.py             |   61 +
 GPCompaReports_v2/tests/test_landing_counts.py     |   56 +
 GPCompaReports_v2/tests/test_new_pages.py          |   57 +
 GPCompaReports_v2/tests/test_selection.py          |   27 +
 GPCompaReports_v2/website/brand.py                 |    8 +
 .../website/page_generators/contact_page.py        |   28 +
 .../website/page_generators/downloads_page.py      |   36 +
 .../website/page_generators/gpcr_index.py          |    4 +-
 .../website/page_generators/gpcr_report_page.py    |    8 +-
 .../website/page_generators/landing_page.py        |   55 +-
 .../website/page_generators/statistics_page.py     |    4 +-
 GPCompaReports_v2/website/site_generator.py        |   38 +-
 docs/superpowers/EXTERNAL_SETUP.md                 |   76 +
 .../plans/2026-08-24-gpcompare-v3-demo.md          | 1566 ++++++++++++++++++++
 .../plans/2026-08-25-v3-review-draft.md            | 1079 ++++++++++++++
 .../specs/2026-08-24-gpcompare-v3-design.md        |  162 ++
 .../specs/2026-08-24-v3-report-layouts.html        |  308 ++++
 .../specs/2026-08-25-v3-open-questions-brief.md    |  136 ++
 scripts/build_demo.sh                              |   66 +
 40 files changed, 4697 insertions(+), 436 deletions(-)
```

This excludes `docs/superpowers/V3_DRAFT_HANDOFF.md` itself and the ten
files under `docs/superpowers/screenshots/`. Every revision of this document,
including this one, is staged and committed as a docs-only change on top of
the diff above, the same convention this item used when the document was
first written.

Deleted in this range and not shown above because they no longer exist to
diff: `GPCompaReports_v2/analysis/receptor_profile.py`,
`GPCompaReports_v2/static/js/v3-analysis.js`,
`GPCompaReports_v2/tests/test_receptor_profile.py` (the removed analysis
layer, see item 7).

## 6. Test results

> **Release gate (added after this section was written).** The counts below
> are a snapshot from an earlier commit and are no longer current: the suite
> is now 23 Python tests plus a 15-test browser suite for `v3-nav.js` and
> `v3-deeplink.js`. More importantly, the concern this section raises about
> skipped tests is now written down and enforced. See
> `docs/superpowers/RELEASE_CHECKLIST.md`. `bash scripts/run_browser_tests.sh`
> is a required pre-release and pre-deployment check, it runs automatically
> inside `scripts/deploy_pages.sh` against the freshly built site before
> anything is pushed, and a skipped run is treated as a failure rather than
> a pass.


```
$ python3 -m pytest GPCompaReports_v2/tests -v
...
GPCompaReports_v2/tests/test_brand.py::test_brand_constants PASSED
GPCompaReports_v2/tests/test_brand.py::test_no_old_brand_text_left_in_templates PASSED
GPCompaReports_v2/tests/test_brand.py::test_urls_and_hosts_are_untouched PASSED
GPCompaReports_v2/tests/test_freeze.py::test_report_sections_are_unchanged_with_features_off PASSED
GPCompaReports_v2/tests/test_freeze.py::test_no_v3_markup_inside_approved_sections PASSED
GPCompaReports_v2/tests/test_landing_counts.py::test_frozen_counts_appear_verbatim PASSED
GPCompaReports_v2/tests/test_landing_counts.py::test_the_promotional_stat_band_is_gone PASSED
GPCompaReports_v2/tests/test_landing_counts.py::test_the_limitation_statement_is_present PASSED
GPCompaReports_v2/tests/test_landing_counts.py::test_no_count_up_animation_remains PASSED
GPCompaReports_v2/tests/test_new_pages.py::test_contact_form_has_every_required_field PASSED
GPCompaReports_v2/tests/test_new_pages.py::test_contact_form_accepts_no_file_uploads PASSED
GPCompaReports_v2/tests/test_new_pages.py::test_contact_page_invents_no_credentials PASSED
GPCompaReports_v2/tests/test_new_pages.py::test_downloads_page_offers_no_dead_links PASSED
GPCompaReports_v2/tests/test_new_pages.py::test_downloads_page_invents_no_doi PASSED
GPCompaReports_v2/tests/test_new_pages.py::test_both_pages_are_reachable_from_the_nav PASSED
GPCompaReports_v2/tests/test_selection.py::test_no_filters_returns_everything PASSED
GPCompaReports_v2/tests/test_selection.py::test_limit_takes_a_prefix PASSED
GPCompaReports_v2/tests/test_selection.py::test_only_selects_named_receptors_in_given_order PASSED
GPCompaReports_v2/tests/test_selection.py::test_only_beats_limit PASSED
GPCompaReports_v2/tests/test_selection.py::test_unknown_id_raises_with_the_bad_name_in_the_message PASSED

============================== 20 passed in 0.76s ==============================
```

Re-run for this fix wave: same 20 tests, same names, all still pass. Test count
is unchanged; `test_frozen_counts_appear_verbatim` (`test_landing_counts.py`)
was strengthened in commit `608362c` rather than replaced, since the family
count of 60 could no longer be pinned with a bare substring check (item 12).

Caveat: most of these tests are skipped, not run, unless
`GPCompaReports_v2/output_v3_demo/` exists on disk (built by
`scripts/build_demo.sh`). That build was present when this run was taken.
Verified directly what happens without one: with `output_v3_demo/` moved
aside, 9 of these 20 tests still run and pass. `test_brand.py` (3 tests)
and `test_selection.py` (5 tests) carry no skip guard at all, and
`test_landing_counts.py::test_no_count_up_animation_remains` reads
`static/js/landing.js` directly rather than the build, so it runs too.
The other 11, all of `test_freeze.py` and `test_new_pages.py` plus the
three remaining `test_landing_counts.py` tests, skip. There is no CI for
this project, so those 11 only protect against regressions when someone
remembers to build the demo first. See item 12.

## 7. The analysis layer is completely removed

"The analysis layer" refers to the opt-in per-receptor segment fingerprint
and CFR badge feature (`GPCompaReports_v2/analysis/receptor_profile.py`,
`GPCompaReports_v2/static/js/v3-analysis.js`,
`GPCompaReports_v2/tests/test_receptor_profile.py`), deleted in commit
`96accd9` because, per that commit's message, "it needed three rounds of
fixes, its comparison had to be renormalised onto a different base, and it
never produced a defensible result."

```
$ ls GPCompaReports_v2/analysis/receptor_profile.py \
     GPCompaReports_v2/static/js/v3-analysis.js \
     GPCompaReports_v2/tests/test_receptor_profile.py
ls: cannot access '...receptor_profile.py': No such file or directory
ls: cannot access '...v3-analysis.js': No such file or directory
ls: cannot access '...test_receptor_profile.py': No such file or directory

$ git grep -niE "receptor_profile|v3-analysis|buildFingerprint|fingerprint|cfr.?badge" \
    -- GPCompaReports_v2
(no output, exit code 1: zero matches in any tracked file)

$ grep -rliE "receptor_profile|v3-analysis|fingerprint" GPCompaReports_v2/output_v3_demo/
(no output, exit code 1: zero matches in the built demo)
```

Sticky section nav, Compact view, and the snake-plot deep link (`sec=`,
`view=`, etc.) are unrelated features and were kept.

## 8. Approved report content is unchanged by default

`GPCompaReports_v2/tests/test_freeze.py` does not build anything itself:
it reads whatever PAR2 report happens to be sitting at
`GPCompaReports_v2/output_v3_demo/reports/par2_human.html` (built with
every V3 toggle off by `scripts/build_demo.sh`) and compares it, section
by section, to `GPCompaReports_v2/output/reports/par2_human.html`, the
last full build the PI approved. Two assertions:

- `test_report_sections_are_unchanged_with_features_off`: strips tags from
  every `<section class="report-section">` block and requires the text to
  match the baseline exactly, section for section.
- `test_no_v3_markup_inside_approved_sections`: requires no `v3-` prefixed
  class or id token inside any approved section, catching markup that
  carries no visible text (an empty mount point, a hidden span).

Both pass (item 6). Independent confirmation from the diff itself: `git
diff 3a10a53..HEAD -- GPCompaReports_v2/templates/gpcr_report.html` shows
every change inside the template is structural, not textual: an `id=` and
`data-section-title=` attribute added to each `<section>` tag (for the new
sticky nav and deep links), a `v3.css` stylesheet link, and a
`window.SnakeAPI` object plus two new `<script>` tags appended after the
existing inline script. No heading, label, or paragraph text inside any
`report-section` was touched.

**Checked beyond the suite, by the final whole-branch review.** The
automated test only reads PAR2. The freeze property was also checked by
hand across all five demo receptors (`par2_human`, `par1_human`,
`adrb2_human`, `5ht2a_human`, `cxcr4_human`), using the same
section-by-section text comparison the test itself runs, and every one
came back byte-identical to its pre-branch baseline. Re-confirmed directly
while writing this revision of the document: comparing
`GPCompaReports_v2/output/reports/*.html` against
`GPCompaReports_v2/output_v3_demo/reports/*.html` for all five gives zero
section mismatches. Separately, the injected sticky nav and toolbar
(`v3-nav.js`) were confirmed at runtime, in headless Chromium, to render
outside every `<section class="report-section">` on all five report
pages, closing part of the runtime gap the next paragraph describes. That
check was manual; it added no automated coverage, and `test_freeze.py`
itself was not changed to run a browser.

**What this does not prove.** Both assertions are static-markup checks
against the built HTML; neither one loads a page or runs a browser, so
the freeze test itself is silent on runtime JavaScript behaviour (the
manual check above covers only where the injected markup renders, not
what it does). This branch adds two script files that run on every
report page, `v3-nav.js` (builds the sticky section nav and the
Compact-view toggle from `[data-section-title]`) and `v3-deeplink.js`
(reads and writes snake-plot view state in the URL hash). Neither is
covered by this test or any other (see item 12). Separately, "unchanged
by default" is only true for a first-time visitor: `v3-nav.js` persists
its Compact-view toggle to `localStorage` under the key
`gpcompare-compact`, so a reader who once switched Compact view on will
see it on again on a later visit to any report page, with no toggle
interaction in between.

See item 10 for a caveat on how this baseline file itself is sourced.

## 9. Screenshots

Ten PNGs under `docs/superpowers/screenshots/`, light theme forced, one
pair per target. Desktop is 1440px wide and mobile 390px wide throughout;
the six landing/PAR2/CFR captures use a fixed 900px/844px-tall viewport,
and the four Contact/Downloads captures (re-taken in this fix wave) use a
full-page capture trimmed to actual content height, since both pages are
short enough to show whole rather than cropped.

| File | Captures |
|---|---|
| `landing-hero-desktop.png` / `-mobile.png` | Landing page top: hero, nav, dataset table start |
| `par2-worked-example-desktop.png` / `-mobile.png` | "How to read a GPCompaRe contact" section, full |
| `cfr-section-desktop.png` / `-mobile.png` | "Core Functional Residues across Class A GPCRs" section, full |
| `contact-page-desktop.png` / `-mobile.png` | Contact page, full: heading, all seven form fields, both standing notes, the "not connected yet" warning, the disabled Send button, and the Lab panel |
| `downloads-page-desktop.png` / `-mobile.png` | Downloads page, full: release metadata, both archive panels including the complete provisional callout, and the footer note |

**Four of the ten were clipped in an earlier pass and have been
re-captured.** The Contact and Downloads pairs were first captured at a
fixed viewport height (900px desktop, 844px mobile) with no scrolling, so
each image stopped at the fold: `contact-page-desktop.png` cut off the
disabled Send button and the "not connected" warning, the two most
decision-relevant things on that page; `contact-page-mobile.png` stopped
mid-form with no Message field, standing notes, warning, Send button, or
Lab panel; `downloads-page-desktop.png` cut mid-sentence through the
provisional callout; `downloads-page-mobile.png` showed only the top
third of the Database release panel, with the Analysis software panel
entirely off-frame. All four were re-captured as full-page images (a tall
off-screen viewport, auto-trimmed to the actual content height) and
re-verified by eye; the row descriptions above reflect what they now
show, not what an earlier draft of this table claimed.

The re-capture also needed one flag not in the original recipe:
`--allow-file-access-from-files`. Without it, headless Chromium loads the
page's own stylesheets as if the `file://` origin were unset (this
appears to be specific to this machine's snap-confined Chromium build)
and the screenshot comes back as unstyled default-browser HTML, no error
printed. This was caught by eye, not by tooling; every re-capture in this
pass was taken with the flag and confirmed styled correctly.

All ten screenshots were opened and visually reviewed as PNGs before
being accepted for this fix wave.

**A capture bug was found and worked around, not just used as given.** The
proven recipe from `.superpowers/sdd/progress.md` (forcing light theme with
`--blink-settings=preferredColorScheme=1`) works correctly for a page URL
with no fragment. It does not work for a `file://` URL with a `#fragment`:
Chromium headless's `--screenshot` renders those as a blank page. This was
first caught because the initial worked-example and CFR captures
(`index.html#worked-heading`, `index.html#cfr-heading`) came back
byte-identical to each other and an order of magnitude smaller than the
hero capture. It was confirmed with a minimal two-line reproduction (a
tall test page with a coloured anchor target): the DOM loads correctly
(`--dump-dom` shows the full page) but the screenshot itself paints solid
background colour, at any window size, with or without
`--virtual-time-budget` or `--run-all-compositing-stages-before-draw`.

The fix used: render the whole landing page as one tall screenshot (no
fragment), locate the two target sections by scanning the image for the
page's own 1px section-divider rule (a distinct border colour, found at
consistent y-offsets), and crop to each section with Pillow. Both sections'
crops were re-verified by eye afterward. This is a capture-tooling issue
only; it does not reflect anything about the site's own anchor or deep-link
behaviour, which is untouched by this task.

**Dark theme, closed.** All ten committed screenshots were forced to
light theme (`--blink-settings=preferredColorScheme=1`); none of them show
the site's dark palette (`primitives-dark.css`), and that has not changed
in this fix wave. The open question this originally raised, whether dark
theme even renders correctly, is now closed: the final whole-branch review
captured and checked dark theme separately, on four page types, and it
renders correctly with no defect found. That check was a throwaway
capture, not saved to `docs/superpowers/screenshots/`, so there is still
no persisted dark-theme image in this repository, only the confirmation
that it works. See item 12; the "still open" list there no longer carries
this item.

## 10. Concerns affecting scientific accuracy, reproducibility, or the manuscript-to-site match

**Y160: annotation coverage, not a numbering conflict. CLOSED 2026-08-26.**
Y160 is `3.37x37` in the manuscript, in the full GPCRdb residue map, and in
the report snake plot (verified in the SVG: a residue circle titled
`Y160 3.37x37`, and `3.37x37` is present in the built report). All three
sources agree on the position. What the original sparse per-report annotation
CSV lacks is a row for residue 160 at all: it carries 159 (M, `3.36x36`) and
then skips to the next annotated residue. So some report-level contact data
may not display that generic number.

This is an annotation-coverage limitation, the same one that makes a site
recomputation give 356 recurrent positions where the submitted full-GPCRdb-map
analysis gives 368. It is not a disagreement about which number Y160 carries,
and there is no authoritative-source decision to make. Two earlier framings
here were wrong and are superseded: first that there was "no entry for
position 160 anywhere", then that this was "two numbering sources
disagreeing".

It does not affect the selected M159-F300 explainer, which uses `3.36x36` and
`6.48x48`, both present. Consistently, the M159-Y160 contact falls below
PAR2's own threshold (\|ΔRRCS\| 1.78 against 2.60, item 2), matching Figure 4e.
No report or annotation change was made, because the report content is
frozen.

**Four receptors have contact data but are absent from the variant dataset,
for two different reasons.** The batch run has 283 `*_rrcs_delta.csv` files
(in `The_batch_RRCS_analyzer/batch_analysis_full/batch_analysis_20260202_151051/csv_data/`)
but only 279 `*_variants.csv` files (in that same run's `variants/`
subdirectory, not `csv_data/`). `gp182_human`, `npy6r_human`, `p2ry8_human`
and `taar3_human` have no variants file, so those four report pages render
no variants section at all (confirmed during Task 6 of this effort: zero
hits for "gnomAD missense variants located" on any of the four built pages).

An earlier version of this paragraph was headed "no gnomAD variant data".
That is wrong for two of the four, and the trace of 2026-08-26 replaced it:

- `npy6r_human` and `taar3_human` resolve to annotated pseudogene loci
  (`NPY6RP` `ENSG00000226306`, `TAAR3P` `ENSG00000179073`). gnomAD carries
  1,567 and 568 variants there and **zero missense** in each, so no variant
  can be placed on a residue. Their absence is expected and correct.
- `gp182_human` and `p2ry8_human` **do** have gnomAD missense variants:
  `ENSG00000166856` (indexed by gnomAD as ACKR5) has 953, and
  `ENSG00000182162` has 526. They are absent because their identifiers were
  not resolved during data collection. The pipeline queried GP182 as `ACKR5`,
  a symbol HGNC withdrew, and P2RY8 is protein coding in the pseudoautosomal
  region where UniProt lists a second Ensembl ID, `ENSG00000292333`, that
  gnomAD does not recognise. The pipeline's log cannot distinguish these from
  a genuine zero: an Ensembl symbol-lookup failure and an empty gnomAD result
  both return `[]` and log "No gnomAD variants found".

The statistics page states both causes separately, and a test rejects any
wording that collapses the four into one "no data available" line. The site
keeps the submitted 279-receptor dataset unchanged: repairing the two would
add variants sections to frozen report pages, move the 279 count away from
submitted Supplementary Table S1, which carries the identical gap, and change
Figure 3's denominators. **Known release limitation carried to the later
manuscript revision, not an open website bug.**

**"Complete RRCS results" is capped at 1,000 rows.** The report section
heading calls this "Complete," but `_get_complete_rrcs()` truncates at
`n=1000`. Seven of 283 receptors have more contact pairs than that in their
delta table (verified directly against the batch CSVs: `lgr4_human`
2,490 rows, `lgr6_human` 2,382, `tshr_human` 1,864, `rxfp2_human` 1,767,
`lshr_human` 1,704, `fshr_human` 1,688, `gp149_human` 1,139), so the
heading is inaccurate for those seven. That said, the section's own
subtitle already reads "Top 1000 contact pairs ranked by |ΔRRCS|"
(`gpcr_report.html:329`), which softens the case: a reader who continues
past the heading to the subtitle is told the true scope, so no reader is
left with a wrong number, only a heading that overstates it in isolation.
The proposed correction, "RRCS results, up to 1,000 contact pairs," is
pending PI approval and is deliberately not applied in this draft.

**The freeze test's baseline is unpinned.** (Restated as the blocking
warning at the very top of this document; full detail here.) `test_freeze.py` compares the
new build against `GPCompaReports_v2/output/reports/par2_human.html`,
whichever build happens to be sitting in that directory. `output/` is
gitignored (`.gitignore:17`), so that baseline is not itself version
controlled; it happens to predate this branch right now, but the first
full 283-receptor build run from `feature/v3-demo` would overwrite it. If
that happens, the freeze test would start comparing the branch to itself
and pass forever, silently losing the one property this whole draft is
built to preserve. This needs pinning, for example a checked-in fixture or
a build tagged and preserved outside `output/`, before this stops being a
draft.

**"CFR" denotes two different sets on the statistics page.** The ranked
CFR table and its dot plot both take the top 30 of the full ranked table:
`cfr_analysis.make_cfr_dotplot` (`cfr_analysis.py:112`) and
`statistics_page.generate_statistics_page` (`statistics_page.py:43`) both
call `.head(30)`, matching the "Top 30 Core Functional Residue Positions"
heading and the landing page's "ranked Top 30 recurrent positions" link.
But `cfr_analysis.build_cfr_network` (`cfr_analysis.py:138`) and
`variant_correlation._get_cfr_position_map` (`variant_correlation.py:41`)
still call `.head(50)`. Two consequences follow directly: the table headed
"Top CFR Contact Pairs (both residues are CFRs)" (`statistics.html:129`)
can list, and does list, pairs where a residue ranks between 31 and 50, so
some rows are not actually pairs of Top 30 CFRs despite the page's own
heading; and the AlphaMissense pathogenicity enrichment percentages and
the chi-squared statistic shown in the pathogenicity chart on the same
page are computed against the top 50 CFR positions, not the top 30.

This split was deliberately not changed. Aligning both functions on 30
would change the published enrichment percentages and the chi-squared
p-value, and those numbers may already be quoted in the manuscript;
changing the underlying population is a scientific call, not a copy fix,
and is not something a documentation or bug-fix pass should decide on its
own. Two options for the owner:

- **Align on 30.** Change `build_cfr_network` and `_get_cfr_position_map`
  to `.head(30)`, so every use of "CFR" on the statistics page means the
  same 30 positions the landing page links to. This changes the
  contact-pair table's contents, the pathogenicity enrichment percentages,
  and the chi-squared p-value; any manuscript citation of the current
  numbers would need to be re-verified against the new ones.
- **Leave the populations different, but say so.** Keep `.head(50)` where
  it is and add a line to the statistics page stating that the contact-pair
  table and the pathogenicity statistics are computed over the top 50 CFR
  positions, a broader set than the top 30 shown in the ranked table above.
  This preserves the numbers as currently quoted, at the cost of "CFR"
  meaning two different things on one page unless a reader notices the new
  disclaimer.

No default is proposed here; this is a decision being surfaced, not one
resolved in this draft. See item 13.

**The database download must be a real archive before public release.**
`downloads.html` currently marks both the database archive and the
software archive "Not yet available," and that is a safe *draft* state
because it invents no file, link, or DOI (verified by
`test_downloads_page_invents_no_doi` and
`test_downloads_page_offers_no_dead_links`). It is not a safe *publication*
state for the database half.
Only the software placeholder may persist past publication, since there is
no releasable user-facing analysis program yet. The database placeholder
must be replaced with a real, licensed archive (see item 13, decision 4)
before the site is published for the paper. `downloads_page.py`'s own
module docstring already states this ordering; this item exists so it is
not missed at review time.

## 11. External setup still required

Nothing below is an oversight; each is a gap the implementation left
because the account, credential, or license decision did not exist to
invent. The full checklist, with exact steps and file locations, is
`docs/superpowers/EXTERNAL_SETUP.md`. Summary:

- **Contact page**: a Formspree account and endpoint (recipient must be
  `tggrlab@gmail.com` or forward there), a reCAPTCHA v3 key pair, a domain
  restriction once the site has a production domain, and named scientific
  and technical contacts to replace the current placeholder paragraph. The
  Send button re-enables itself automatically once the endpoint is real; no
  further template change is needed for that part.
- **Downloads and data release**: a data license decision, a software
  license decision, a Zenodo deposit and DOI, and (per item 10 above)
  replacing the database "Not yet available" block with a real archive.
  The software archive may remain a placeholder past publication.
- **What is already safe to ship without these**: both the Contact and
  Downloads pages can be built and deployed as they stand. The Contact
  form is disabled with a visible notice rather than silently dropping
  submissions, and the Downloads page invents no file, link, or DOI.
  Publishing either page with real credentials, named contacts, or a real
  data archive is a separate, later release step.


> **Superseded 2026-08-26.** The form is now live: it posts as plain HTML to `https://formspree.io/f/mljrqazl`, the submit button and all visible fields ship enabled, and the "not connected" notice no longer renders. What remains is dashboard-only and listed under "Contact form: manual Formspree dashboard checks" in `RELEASE_CHECKLIST.md`.

## 12. Other findings from the working ledger

Carried forward from `.superpowers/sdd/progress.md` for completeness. None
of these blocks this draft; they are recorded so they are not lost.

**Fixed during this branch, kept for the record:**

- A controller error swept two untracked preview build directories
  (`output_par2_preview/`, `output_tm_preview/`, 23 and 14 files) into git
  history via an over-broad `git add -A GPCompaReports_v2/` for what was
  meant to be a one-line fix. Both directories were untracked again and
  added to `.gitignore`. The plan was updated to forbid broad `git add`,
  and this task followed that rule (explicit paths only, see item 5's
  closing note).
- A second, later instance of the same mistake: an over-broad
  `git add -A GPCompaReports_v2/` swept three local tooling scripts
  (`scripts/build_par_dossier_csvs.py`, `poster_figures.py`,
  `preview_par2.py`, 665 lines) into history. Untracked again in commit
  `03c94a3`; they remain on disk, unaffected. See item 5 for the
  verification command and the regenerated files-changed table.
- The landing page's worked-example prose called M159 "the methionine of
  the xMY motif." The project's own fact-check rejects that label: the
  adjacent Tyr3.37 fails both magnitude and presence under two separate
  aggregations, so it cannot anchor an "xMY" motif read off this data, and
  the label this dataset does support, M3.36, is itself flagged elsewhere
  as a reconstruction, not a directly observed motif. Stating "xMY" on the
  public landing page asserted a specific residue identity beyond what the
  underlying analysis defends. Fixed in commit `2445658`: the sentence now
  describes F300 as occupying the 6.48 toggle-switch position instead,
  which is uncontested, and states only that both residues sit in
  transmembrane helices. See item 4 for the copy as it now reads.
- The landing page's "Analysing your own structures" section said
  GPCompaRe "also includes" a program for analysing user-supplied
  structures, while the Downloads page says no releasable version exists.
  A visitor reading only the landing page would expect a working tool the
  Downloads page has never offered. Fixed in commit `38c5d1b`: the landing
  page now says a releasable version is not yet available, matching
  Downloads. See item 4.
- The Downloads page's `.downloads-block-head` row (heading, status badge,
  and on the software panel also the version string) did not wrap. At a
  390px viewport its combined minimum width exceeded the available space,
  and because both panels share one CSS grid track, the overflow applied
  to both, not only the panel with the extra child. Found while
  re-capturing `downloads-page-mobile.png` for this fix wave (the earlier
  clipped capture had cut off the affected area, hiding it). Fixed by
  adding `flex-wrap: wrap` to `.downloads-block-head` in `site.css`;
  confirmed by eye in the re-captured screenshot. This did not touch
  `.badge`, the shared badge class that also appears in `gpcr_report.html`;
  only the Downloads-specific head row changed, so the approved report
  pages are unaffected.
- The statistics page carried both a "Top 50" dot-plot title and a "Top 30"
  table heading for the same Core Functional Residue ranking. The dot plot
  was changed to `.head(30)` and retitled to match, since the new landing
  CFR section links to it saying "Top 30." (`cfr_analysis.py:138`, a
  different function, `build_cfr_network`, still uses `.head(50)`, but it
  renders no "Top 50" claim on any built page, so it is not a live
  inconsistency today.)
- A stale docstring on `identify_cfrs` claimed it truncated to the top 50
  CFR positions; it never truncated at all. Fixed.
- A fifth landing-page change arrived after this document's item 4 was
  first drafted and independently of the fix wave that first recorded it
  here: commit `6b6fcaa` added a paragraph to the CFR section stating that
  recurrence counts are lower bounds, because GPCRdb generic numbering
  reaches only a median 64% (range 8% to 88%) of a receptor's
  above-threshold contact positions. That figure counted contact-position
  *occurrences*, weighting a position once per receptor where it recurs;
  commit `608362c` corrected it to a median 58% (range 9% to 88%) of
  *distinct* above-threshold contact positions, the more conservative,
  position-based count. Commit `cc8ae9b` then added the same lower-bounds
  caveat, at 58%, to the statistics page's own CFR callout
  (`statistics.html:66-73`), which previously claimed every residue is
  mapped "across all 283 GPCRs," contradicting the landing page. Item 4's
  CFR summary and item 3's count table now carry the corrected 58%/9%/88%
  figure. That figure is still hand-authored prose, not something any code
  in this branch independently re-derives.
- The Contact page's placeholder text originally pointed readers at
  `docs/superpowers/EXTERNAL_SETUP.md`, an internal repository path, and
  would have leaked the internal working-directory name onto a page meant
  for public visitors. Replaced with plain prose. Verified: no built page
  contains the strings "superpowers" or "EXTERNAL_SETUP".
- `downloads_page.py` passed `software_name` to the template but the
  template never rendered it. Fixed. Two tests were weaker than their
  names suggested (`test_both_pages_are_reachable_from_the_nav` only
  checked that a string appeared in markup, true even before either page
  existed; the archive-extension check missed `.tgz`/`.7z` and case
  variants). Both strengthened, and a new
  `test_downloads_page_invents_no_doi` was added. Suite grew from 19 to 20
  tests.
- One editorial aside ("A paper companion that offers no data is not a
  companion") was trimmed from the Downloads placeholder note as
  out-of-register phrasing that had leaked from a work-in-progress draft
  into public copy; the substantive commitment sentence was kept.
- The final whole-branch review forced a further wave of fixes, commit
  `608362c`, that this document did not yet reflect until this revision:
  - The hero claimed every report links gnomAD variation. Four receptors
    have none (item 10). The correct scoping already existed further down
    the same page; the hero was never brought in line with it. Item 4 has
    the corrected hero text.
  - "Receptor-level Core Functional Residues," in the report-contents
    list, named the wrong set: a report shows the family-scope CFR list
    mapped onto that one receptor, not something computed at the receptor
    level. Reworded to "Cross-family Core Functional Residue positions
    mapped onto this receptor." Item 4 has the corrected bullet.
  - The statistics page's CFR score tooltip said "frequency × mean
    absolute delta." The code (`cfr_analysis.py:91-96`) computes the mean
    of normalised frequency and normalised mean delta, not a product.
    Tooltip corrected to describe the actual formula
    (`statistics.html:105`).
  - Compact view folded a report section's body but left the snake plot
    rendered while hiding its legend and controls, so an approved figure
    appeared without its colour key. Cause: `#snake-plot-container`
    carries `display: flex` from an id selector in `site.css`, which
    outranks the `.report-section.v3-folded > *` attribute rule used to
    hide folded content. Fixed with an id-matching rule in `v3.css` that
    hides the plot and its legend specifically when folded.
  - The dataset table's "Contact-pair records" note said "residue pairs
    scored in both states." Only 60.1% of the 213,456 rows are actually
    non-zero in both states; the rest are zero in one. Reworded to
    "residue pairs detected in either state," matching the wording the
    approved browse page already uses for the same population.
  - `test_landing_counts.py::test_frozen_counts_appear_verbatim` pinned
    the receptor-family count with a bare `'60' in html` check. The
    string "60" occurs many times on the landing page, so a rebuild that
    changed the family count to 61 would still have passed. Anchored to
    the count's own table row instead, and proven able to fail (item 6).
- Commit `cc8ae9b` fixed one further false claim the same review caught:
  the statistics page's CFR callout said every residue is mapped "across
  all 283 GPCRs." Generic numbering does not cover every residue; the
  callout now states the same 58% median lower-bound caveat the landing
  page carries (`statistics.html:66-73`, item 3, item 4).

**Verified independently by the final whole-branch review, recorded here
for the first time:**

- The freeze property holds under a stronger check than `test_freeze.py`
  itself runs: byte-identical approved sections across all five demo
  receptors, not only PAR2. Re-confirmed directly while writing this
  revision; see item 8.
- The injected sticky nav and toolbar were confirmed at runtime, in
  headless Chromium, to render outside every approved report section on
  all five report pages. A manual check, not automated coverage; see
  item 8.
- Every frozen figure and the entire PAR2 worked example reproduce
  exactly from the underlying batch CSVs, corroborating item 2's
  independently re-derived verification table.
- Dark theme was captured and checked on four page types and renders
  correctly with no defect found. See item 9, which previously listed
  this as unreviewed and now records it as closed.
- The DRY, NPxxY and CWxP recovery claim on the landing page holds under
  direct inspection of the built Top 30 table
  (`GPCompaReports_v2/output_v3_demo/statistics.html`): `3.50x50` (the DRY
  arginine) ranks 1, `7.53x53` (the NPxxY tyrosine) ranks 2, and `6.48x48`
  (the CWxP tryptophan) ranks 3. Re-confirmed directly for this revision.
  The hedge word "parts" in that landing-page sentence is accurate:
  `N7.49` (NPxxY) and both `C6.47` and `P6.50` (CWxP) do not appear
  anywhere in the Top 30, so no one motif is recovered in full, only
  parts of each.

**A pattern worth naming before the next review pass.** Three of the
false claims the final review caught, the hero's unscoped gnomAD claim,
the "Receptor-level Core Functional Residues" bullet, and the
occurrence-weighted 64% coverage figure, were sentences lifted verbatim
from the implementation plan and never checked against the pipeline that
actually produces the numbers. Each earlier review in this branch's
history fixed only the single instance directly in front of it, not the
pattern behind all three. Recommendation for whoever runs the next pass:
before reviewing prose, walk the plan's copy blocks against the built
site as a checklist, one sentence-to-code check per claim, rather than
reacting to whichever sentence a reviewer happens to notice first.

**Still open, not blocking, worth carrying forward:**

- `_prepare_output` (`GPCompaReports_v2/website/site_generator.py:107-116`)
  never cleans a pre-existing output directory before copying `static/`
  into it, so a deleted static asset persists in that directory forever.
  `scripts/build_demo.sh` works around this by clearing `$OUT/static`
  first; a durable fix belongs in the generator itself, as separate work.
- Three of the four new landing-count tests are skipped rather than run
  when `GPCompaReports_v2/output_v3_demo/` does not exist on disk, and so
  is every test in `test_freeze.py` and `test_new_pages.py`, 11 tests in
  total. `test_brand.py` (3 tests), `test_selection.py` (5 tests), and
  `test_landing_counts.py::test_no_count_up_animation_remains` carry no
  such guard and run regardless, so 9 of 20 tests run and pass with no
  demo build present (verified directly, see item 6). This project has no
  CI, so in practice a regression among the skipped 11 is only caught when
  someone has run `scripts/build_demo.sh` first. A CI run without a
  preceding build reports 9 passed and 11 skipped, not the whole suite
  failed, and not the whole suite skipped either; a reviewer who does not
  check which 11 skipped could still wave that through as clean.
- Neither `v3-nav.js` nor `v3-deeplink.js`, the two script files this
  branch adds that run on every report page (item 8), has any test
  coverage, and this project has no JavaScript test harness at all, only
  the Python/pytest suite under `GPCompaReports_v2/tests/`. Adding one is
  out of scope for this draft but is the natural next step before relying
  on either script's behaviour for a release.
- Round 2 of this effort reused round 1's task-brief and task-report
  filenames under `.superpowers/sdd/`, so round 1's working files were
  overwritten as round 2 proceeded. Round 1 is complete and was committed
  before that happened; the durable record of it is git history plus
  `.superpowers/sdd/progress.md`, not those overwritten files.
- The general `.btn:disabled` rule added to `site.css` for the Contact
  page's disabled Send button applies site-wide, not just on that page. It
  was confirmed not to affect anything else (it does not match the
  pagination buttons, which are not `.btn`, and never matches an anchor),
  but a one-line comment noting the scope would help a future author.
- `.callout-head` is used inline as a run-in phrase in the Downloads
  database panel, where every other use in the codebase is block-level.
  The two Downloads panels also render at slightly different heights
  because only the database panel carries the extra callout box.
- For anyone touching the Plotly figures later: Plotly 6 encodes numeric
  trace arrays as base64 `bdata` dictionaries, so `len(trace['x'])`
  reports the dictionary's key count, not the number of data points. Count
  points via the `text` array instead.

Full blow-by-blow detail, including exact image-crop parameters used for
the CFR topology and PAR2 snake-plot figures, is in
`.superpowers/sdd/progress.md`.

## 13. Decisions needed from the owner

Everything above is a finding. This is what those findings are asking the
owner to decide, in priority order, with what each decision blocks.

1. **Which numbering source is authoritative for PAR2 position 160.** The
   manuscript's Figure 4 and the report's snake plot both use `3.37x37`
   for Y160; the annotation CSV that feeds the report's tables has no
   entry for that position (item 10). Blocks correcting the annotation
   file and closing out the manuscript-to-site match for this receptor.
2. **Whether "Complete RRCS results" should be relabelled** to "RRCS
   results, up to 1,000 contact pairs" (item 10). A drafted correction
   exists and is deliberately not applied. Blocks a one-line template
   change; nothing else depends on it.
3. ~~**Whether the four receptors missing gnomAD variant data should have
   it**~~ **RESOLVED for this release, 2026-08-26.** Traced to two distinct
   causes. `npy6r_human` and `taar3_human` are annotated pseudogene loci with
   zero missense variants in gnomAD, so their absence is expected and correct.
   `gp182_human` and `p2ry8_human` do have gnomAD missense variants (953 and
   526) and are absent because their identifiers were not resolved during data
   collection. The website records both causes separately on the statistics
   page and keeps the submitted 279-receptor dataset. The two identifier
   failures are a **known release limitation carried to the later manuscript
   revision**, not an open website bug: fixing them would add variants sections
   to two frozen report pages, move the 279 count and change Figure 3's
   denominators. No website work is blocked.
4. **Data licensing, software licensing, and a Zenodo deposit and DOI**
   (item 11). Blocks replacing the Downloads page's "Not yet available"
   database placeholder with a real archive, which item 10 identifies as
   required before public release. The software placeholder may remain
   past release regardless of this decision.
5. **A Formspree account, reCAPTCHA keys, a domain restriction, and named
   scientific/technical contacts** for the Contact page (item 11). Blocks
   enabling the Send button and replacing the current placeholder contact
   paragraph. The Contact page is otherwise safe to ship as-is in the
   meantime (item 11).
6. **Whether to pin `test_freeze.py`'s baseline** rather than leave it as
   whatever build happens to sit in the gitignored `output/` directory
   (item 10). Blocks keeping the one property this whole draft protects,
   the frozen approved report content, reliable across future builds
   without anyone noticing if it silently stopped being checked. The final
   whole-branch review judged this effectively blocking already: see the
   warning at the very top of this document for what happens if a full
   build runs before this is pinned.
7. **Whether to scale this draft to the full 283-receptor build** for the
   next review pass, rather than the five receptors it currently covers
   (item 0). Blocks giving the owner a look at the real site instead of a
   five-receptor sample; everything else in this document was verified
   against that five-receptor sample and would need re-verifying at full
   scale.
8. **Whether to align the statistics page's two CFR populations, top 30
   versus top 50, or document the split instead** (item 10). The ranked
   table and dot plot use the top 30; the "Top CFR Contact Pairs" table
   and the pathogenicity enrichment/chi-squared statistics still use the
   top 50. Aligning on 30 would change the published enrichment
   percentages and the chi-squared result, which may already be quoted in
   the manuscript, so this was deliberately left as a decision rather than
   a fix applied in this draft. Blocks nothing else in this document; it
   is scoped to the statistics page alone.
