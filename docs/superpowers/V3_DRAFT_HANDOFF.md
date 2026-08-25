# GPCompaRe V3 review draft: handoff

Branch `feature/v3-demo`. This document is the record for the project owner
of what this draft contains, what was verified and how, and what is not
resolved yet. It is a draft review artifact, not a release.

Plan: `docs/superpowers/plans/2026-08-25-v3-review-draft.md`.
Running ledger kept during the work: `.superpowers/sdd/progress.md`.

---

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
> report links the structural result to GPCRdb generic numbering, Core
> Functional Residues, conservation, gnomAD variation and AlphaMissense
> predictions.

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
(red), a note that M159 is the methionine of the xMY motif and F300 is the
6.48 toggle-switch aromatic, and a link to the PAR2 snake plot.

**Core Functional Residues across Class A GPCRs**: defines a receptor-level
CFR, explains recurrence across at least three receptors, states that no
motif was picked in advance and the analysis recovers DRY/NPxxY/CWxP
machinery on its own, and links to the ranked Top 30 table.

**What each receptor report contains**: a nine-item list (interactive snake
plot; ranked contact pairs up to 1,000 rows; receptor threshold; CFRs;
generic numbering; conservation; gnomAD variants for the 279 receptors that
have them; AlphaMissense; CSV export), plus example-report links to ADRB2
and PAR2.

**Analysing your own structures**: states GPCompaRe includes a program for
analysing user-supplied paired structures, that local analyses stay local
and are not added to the public database, with a link to Downloads.

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
commit the brief's own verification command uses).

```
$ git diff --stat 3a10a53..HEAD
 .gitignore                                                    |    3 +
 GPCompaReports_v2/analysis/cfr_analysis.py                    |   14 +-
 GPCompaReports_v2/generate_site.py                             |    3 +
 GPCompaReports_v2/scripts/build_par_dossier_csvs.py            |  257 ++++
 GPCompaReports_v2/scripts/poster_figures.py                    |  319 ++++
 GPCompaReports_v2/scripts/preview_par2.py                      |   89 ++
 GPCompaReports_v2/static/css/landing.css                       |  281 +---
 GPCompaReports_v2/static/css/site.css                          |  138 +-
 GPCompaReports_v2/static/css/v3.css                            |   49 +
 GPCompaReports_v2/static/img/fig-cfr-topology.png              |  Bin 0 -> 406319 bytes
 GPCompaReports_v2/static/img/par2-explainer-snake.png          |  Bin 0 -> 257868 bytes
 GPCompaReports_v2/static/js/landing.js                         |   38 +-
 GPCompaReports_v2/static/js/v3-deeplink.js                     |  100 ++
 GPCompaReports_v2/static/js/v3-nav.js                          |  117 ++
 GPCompaReports_v2/templates/_partials/footer.html              |    2 +-
 GPCompaReports_v2/templates/_partials/topbar.html               |   12 +-
 GPCompaReports_v2/templates/base.html                          |    7 +-
 GPCompaReports_v2/templates/contact.html                       |  102 ++
 GPCompaReports_v2/templates/downloads.html                     |   90 ++
 GPCompaReports_v2/templates/gpcr_report.html                   |   56 +-
 GPCompaReports_v2/templates/landing.html                       |  297 ++--
 GPCompaReports_v2/tests/__init__.py                            |    0
 GPCompaReports_v2/tests/test_brand.py                          |   23 +
 GPCompaReports_v2/tests/test_freeze.py                         |   61 +
 GPCompaReports_v2/tests/test_landing_counts.py                 |   51 +
 GPCompaReports_v2/tests/test_new_pages.py                      |   57 +
 GPCompaReports_v2/tests/test_selection.py                      |   27 +
 GPCompaReports_v2/website/brand.py                              |    8 +
 GPCompaReports_v2/website/page_generators/contact_page.py       |   28 +
 GPCompaReports_v2/website/page_generators/downloads_page.py     |   36 +
 GPCompaReports_v2/website/page_generators/gpcr_index.py         |    4 +-
 GPCompaReports_v2/website/page_generators/gpcr_report_page.py   |    8 +-
 GPCompaReports_v2/website/page_generators/landing_page.py       |   55 +-
 GPCompaReports_v2/website/page_generators/statistics_page.py    |    4 +-
 GPCompaReports_v2/website/site_generator.py                     |   38 +-
 docs/superpowers/EXTERNAL_SETUP.md                              |   76 +
 docs/superpowers/plans/2026-08-24-gpcompare-v3-demo.md          | 1566 ++++++++
 docs/superpowers/plans/2026-08-25-v3-review-draft.md            | 1079 +++++
 docs/superpowers/specs/2026-08-24-gpcompare-v3-design.md        |  162 ++
 docs/superpowers/specs/2026-08-24-v3-report-layouts.html        |  308 ++
 docs/superpowers/specs/2026-08-25-v3-open-questions-brief.md    |  136 ++
 scripts/build_demo.sh                                           |   66 +
 42 files changed, 5336 insertions(+), 431 deletions(-)
```

This does not yet include the two files created by this task,
`docs/superpowers/V3_DRAFT_HANDOFF.md` and the ten files under
`docs/superpowers/screenshots/`, which are staged and committed as a
separate change on top of the above.

Deleted in this range and not shown above because they no longer exist to
diff: `GPCompaReports_v2/analysis/receptor_profile.py`,
`GPCompaReports_v2/static/js/v3-analysis.js`,
`GPCompaReports_v2/tests/test_receptor_profile.py` (the removed analysis
layer, see item 7).

## 6. Test results

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

Caveat: most of these tests are skipped, not run, unless
`GPCompaReports_v2/output_v3_demo/` exists on disk (built by
`scripts/build_demo.sh`). That build was present when this run was taken.
There is no CI for this project, so this suite only protects against
regressions when someone remembers to build the demo first. See item 12.

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

`GPCompaReports_v2/tests/test_freeze.py` builds the PAR2 report with every
V3 toggle off and compares it, section by section, to
`GPCompaReports_v2/output/reports/par2_human.html`, the last full build the
PI approved. Two assertions:

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

See item 10 for a caveat on how this baseline file itself is sourced.

## 9. Screenshots

Ten PNGs under `docs/superpowers/screenshots/`, light theme forced, desktop
1440x900 and mobile 390x844, one pair per target:

| File | Captures |
|---|---|
| `landing-hero-desktop.png` / `-mobile.png` | Landing page top: hero, nav, dataset table start |
| `par2-worked-example-desktop.png` / `-mobile.png` | "How to read a GPCompaRe contact" section, full |
| `cfr-section-desktop.png` / `-mobile.png` | "Core Functional Residues across Class A GPCRs" section, full |
| `contact-page-desktop.png` / `-mobile.png` | Contact page, form and lab panel |
| `downloads-page-desktop.png` / `-mobile.png` | Downloads page, both archive panels |

All ten were opened and visually reviewed as PNGs before being accepted;
none is blank, clipped, or showing the wrong content.

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

## 10. Concerns affecting scientific accuracy, reproducibility, or the manuscript-to-site match

**The Y160 numbering mismatch, corrected form.** The manuscript's Figure 4
gives Y160 the generic number `3.37x37`. The snake plot on the PAR2 report
*does* carry that number (verified directly in the SVG: a residue circle
titled `Y160 3.37x37`), so the site can render it. What lacks an entry for
position 160 is the annotation CSV that feeds the report's tables, which is
a different source from the one the snake plot draws on. This is two
numbering sources disagreeing, not one missing number. (An earlier version
of this finding, recorded during this branch's own prep work, stated it
more broadly as "no entry for position 160 anywhere"; that was wrong and
has been corrected here.) Separately and consistently: the M159-Y160
contact itself falls below PAR2's own threshold (\|ΔRRCS\| 1.78 against
2.60, item 2), which matches Figure 4e, so the two sources agree on the
science even though they disagree on the label. Flag for investigation
before release. No change was made to the approved report, the manuscript,
or the annotation file to resolve this. It needs a human decision about
which numbering source is authoritative.

**Four receptors have contact data but no gnomAD variant data.** The batch
run has 283 `*_rrcs_delta.csv` files but only 279 `*_variants.csv` files.
`gp182_human`, `npy6r_human`, `p2ry8_human` and `taar3_human` have no
variants file, so those four report pages render no variants section at
all (confirmed during Task 6 of this effort: zero hits for "gnomAD missense
variants located" on any of the four built pages). The landing page's copy
already reflects this ("for the 279 receptors that have variant data").
Whether those four receptors should have variant data, and if so why the
pipeline did not produce it, is a pre-release question for the lab, not
something resolved by this draft.

**"Complete RRCS results" is capped at 1,000 rows.** The report section
heading calls this "Complete," but `_get_complete_rrcs()` truncates at
`n=1000`. Seven of 283 receptors have more contact pairs than that in their
delta table (verified directly against the batch CSVs: `lgr4_human`
2,490 rows, `lgr6_human` 2,382, `tshr_human` 1,864, `rxfp2_human` 1,767,
`lshr_human` 1,704, `fshr_human` 1,688, `gp149_human` 1,139), so the
heading is inaccurate for those seven. The proposed correction, "RRCS
results, up to 1,000 contact pairs," is pending PI approval and is
deliberately not applied in this draft.

**The freeze test's baseline is unpinned.** `test_freeze.py` compares the
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

**The database download must be a real archive before public release.**
`downloads.html` currently marks both the database archive and the
software archive "Not yet available," and that is a safe *draft* state
because it invents no file, link, or DOI (verified by
`test_downloads_page_invents_no_doi` and
`test_downloads_page_offers_no_dead_links`). It is not a safe *publication*
state for the database half.
Only the software placeholder may persist past publication, since there is
no releasable user-facing analysis program yet. The database placeholder
must be replaced with a real, licensed archive (see item 11, items 6 and
8) before the site is published for the paper. `downloads_page.py`'s own
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
- The statistics page carried both a "Top 50" dot-plot title and a "Top 30"
  table heading for the same Core Functional Residue ranking. The dot plot
  was changed to `.head(30)` and retitled to match, since the new landing
  CFR section links to it saying "Top 30." (`cfr_analysis.py:132`, a
  different function, `build_cfr_network`, still uses `.head(50)`, but it
  renders no "Top 50" claim on any built page, so it is not a live
  inconsistency today.)
- A stale docstring on `identify_cfrs` claimed it truncated to the top 50
  CFR positions; it never truncated at all. Fixed.
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

**Still open, not blocking, worth carrying forward:**

- `_prepare_output` (`GPCompaReports_v2/website/site_generator.py:99-108`)
  never cleans a pre-existing output directory before copying `static/`
  into it, so a deleted static asset persists in that directory forever.
  `scripts/build_demo.sh` works around this by clearing `$OUT/static`
  first; a durable fix belongs in the generator itself, as separate work.
- Three of the four new landing-count tests, and in fact every test module
  in this suite, are skipped rather than run when
  `GPCompaReports_v2/output_v3_demo/` does not exist on disk. This project
  has no CI, so in practice a regression these tests would catch is only
  caught when someone has run `scripts/build_demo.sh` first. A CI run
  without a preceding build reports the whole suite as skipped, not
  failed, which could look deceptively clean.
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
