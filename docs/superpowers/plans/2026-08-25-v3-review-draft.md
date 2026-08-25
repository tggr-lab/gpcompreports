# GPCompaRe V3 Review Draft Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Turn the V3 demo branch into a review draft: remove the analysis layer, rebuild the landing page as an academic resource page grounded in real data, and add Contact and Downloads pages.

**Architecture:** Same Python + Jinja2 generator, same hand-written CSS and vanilla JS, no build step. The landing page becomes eleven ordered sections rendered from real values. Two new standalone page generators join the existing four. Nothing is invented: every number comes from the build, every visual from the paper or the live data.

**Tech Stack:** Python 3.13, Jinja2 3.1.6, pandas 3.0.0, Plotly 6.5.2, pytest 9.0.2, vanilla ES5-style JS, CSS custom properties, poppler-utils for one figure extraction.

## Global Constraints

Every task's requirements implicitly include this section.

- **Draft only. Do not merge, do not deploy, do not push.** Work on `feature/v3-demo`.
- **Visible product name is `GPCompaRe database`**, short form `GPCompaRe`. Every lowercase `gpcompreports` URL and host string stays exactly as it is: the GitHub remote, `https://tggr-lab.github.io/gpcompreports/`, `https://gpcompreports.goatcounter.com/count`.
- **No npm, no bundler, no framework, no server, no new build system.**
- **PI-approved report content and its default presentation must not change.** `GPCompaReports_v2/tests/test_freeze.py` enforces this and must stay green.
- **No generated prose on individual report pages.**
- **No unsupported biological claims.** No claim of experimental function or causality.
- **No em dashes in user-facing copy.** Use a colon, a comma, or a full stop.
- **Use real artifacts, values and visual language. Never invent an explanatory number or diagram.**
- **Terminology is exact:** CFR always means Core Functional Residues. Receptor-level CFRs and recurrent CFR positions are different scopes of the same concept. A recurrent CFR position is one identified as a CFR in at least three receptors. Ranked references use the published **Top 30**. The string "Top 50" must not appear anywhere in user-facing copy.
- **GPCRdb numbering is read from `display_number`.** `generic_number` is empty in all 283 annotation files.
- **Colour convention:** red (`--semantic-active`) is stronger in the active state, blue (`--semantic-inactive`) is stronger in the inactive state.
- **The TGGR logo asset must not be altered or removed.** It may be rendered smaller.
- No marketing register. Banned: "unlock", "powerful platform", "comprehensive insights", "cutting-edge", and anything of that family.
- Do not promise response times, update schedules, feature development, an API, user submissions, or clinical interpretation.
- **Do not modify the system Python installation.** Do not run `pip3 install --break-system-packages`. Use the existing environment, or a project-local virtual environment if something is genuinely missing. Nothing in this plan requires a new dependency: jinja2, pandas, plotly, scipy, pytest and poppler-utils are all already present.
- **The frozen dataset counts, verified against the build. Use these exact figures and no others:** 283 receptors, 60 receptor families, 213,456 contact-pair records, 23,025 threshold-passing changes, 566 models (two per receptor), 4 receptors excluded for lacking a complete active and inactive pair (`gpr26`, `lgr5`, `opsx`, `rxfp1`, from 287 in `class_A_all.csv`).
- **This limitation statement must appear somewhere prominent on the landing page, verbatim:** "Positions and contact pairs highlighted by GPCompaRe should be interpreted as model-derived structural candidates, not direct experimental evidence of function."
- All paths below are relative to `The Ultimate RRCS database project/alphafold_multistate_gpcr/`.

---

## File Structure

**Deleted:**
- `GPCompaReports_v2/static/js/v3-analysis.js`
- `GPCompaReports_v2/analysis/receptor_profile.py`
- `GPCompaReports_v2/tests/test_receptor_profile.py`
- The `.v3-layer`, `.v3-strip`, `.v3-tile*`, `.v3-fp*`, `.v3-badge*`, `.v3-seg-lab`, `.v3-gn` blocks in `GPCompaReports_v2/static/css/v3.css`

**Created:**
- `GPCompaReports_v2/templates/contact.html`
- `GPCompaReports_v2/templates/downloads.html`
- `GPCompaReports_v2/website/page_generators/contact_page.py`
- `GPCompaReports_v2/website/page_generators/downloads_page.py`
- `GPCompaReports_v2/static/img/fig-cfr-topology.png` (extracted from the paper)
- `GPCompaReports_v2/static/img/par2-explainer-snake.png` (cropped from the built PAR2 snake plot)
- `GPCompaReports_v2/tests/test_landing_counts.py`
- `GPCompaReports_v2/tests/test_new_pages.py`
- `docs/superpowers/EXTERNAL_SETUP.md`

**Modified:**
- `GPCompaReports_v2/templates/landing.html` (rebuilt body)
- `GPCompaReports_v2/templates/_partials/topbar.html` and `_partials/footer.html` (navigation)
- `GPCompaReports_v2/static/css/landing.css` (new section styles)
- `GPCompaReports_v2/static/js/landing.js` (remove the count-up animation)
- `GPCompaReports_v2/website/site_generator.py` (drop the median loop, add two generators)
- `GPCompaReports_v2/website/page_generators/gpcr_report_page.py` (drop the analysis payload)
- `GPCompaReports_v2/website/page_generators/landing_page.py` (new context values)
- `GPCompaReports_v2/templates/gpcr_report.html` (drop two analysis tags)

---

### Task 1: Remove the analysis layer

It did not produce a defensible result. Remove it, do not hide it. Everything else V3 added to report pages stays.

**Files:**
- Delete: `GPCompaReports_v2/static/js/v3-analysis.js`, `GPCompaReports_v2/analysis/receptor_profile.py`, `GPCompaReports_v2/tests/test_receptor_profile.py`
- Modify: `GPCompaReports_v2/static/css/v3.css`, `GPCompaReports_v2/website/site_generator.py`, `GPCompaReports_v2/website/page_generators/gpcr_report_page.py`, `GPCompaReports_v2/templates/gpcr_report.html`

**Interfaces:**
- Produces: report pages with sticky nav, Compact view and deep links intact, and no analysis-layer toggle, payload, or markup anywhere.
- The `V3Nav` toolbar still exists because Compact view mounts into it.

- [ ] **Step 1: Confirm what depends on the module before deleting**

Run: `grep -rn "receptor_profile\|v3-analysis\|v3_analysis_json\|median_profile\|segment_profile\|segment_coverage\|cfr_ranks\|low_confidence_segments" GPCompaReports_v2/ --include=*.py --include=*.html --include=*.js --include=*.css`

Expected: hits only in the files this task deletes or modifies. If anything else references them, stop and report it: something outside the layer grew a dependency on it.

- [ ] **Step 2: Delete the three files**

```bash
git rm GPCompaReports_v2/static/js/v3-analysis.js \
       GPCompaReports_v2/analysis/receptor_profile.py \
       GPCompaReports_v2/tests/test_receptor_profile.py
```

- [ ] **Step 3: Remove the payload from the report generator**

In `GPCompaReports_v2/website/page_generators/gpcr_report_page.py`, delete the `cfr_ranks` derivation, the `v3_payload` dict, and the `v3_analysis_json=` keyword in the `template.render(...)` call. Delete the `from ...analysis import receptor_profile as rprofile` import.

- [ ] **Step 4: Remove the median computation from the orchestrator**

In `GPCompaReports_v2/website/site_generator.py`, delete the block that prints "Computing database median segment profile...", its loop, and both `self.analysis_results['median_profile']` and `self.analysis_results['median_coverage']` assignments. Delete the now-unused `receptor_profile` and `gpcr_report_helpers` imports if nothing else in that file uses them. Check before deleting.

- [ ] **Step 5: Remove the two tags from the report template**

In `GPCompaReports_v2/templates/gpcr_report.html`, in `{% block scripts %}`, delete the `<script type="application/json" id="v3-analysis-data">` tag and the `<script src="...v3-analysis.js" defer></script>` tag. Leave the `v3-nav.js` and `v3-deeplink.js` tags alone.

- [ ] **Step 6: Remove the layer's styles**

In `GPCompaReports_v2/static/css/v3.css`, delete every rule whose selector begins `.v3-layer`, `.v3-strip`, `.v3-tile`, `.v3-fp`, `.v3-legend`, `.v3-seg-lab`, `.v3-badge`, `.v3-badges`, `.v3-gn`, and the `.v3-analysis-on` variants. Keep `.v3-secnav*`, `.v3-toolbar`, `.v3-toggle`, `.v3-folded`, the `.report-section { scroll-margin-top }` rule, and the `@media print` block.

- [ ] **Step 7: Verify nothing references the removed names**

Run: `grep -rn "v3-analysis\|v3-badge\|v3-strip\|v3-fp\|v3-layer\|receptor_profile\|gpcompare-analysis" GPCompaReports_v2/ --include=*.py --include=*.html --include=*.js --include=*.css`
Expected: no output.

- [ ] **Step 8: Rebuild and confirm the report still works**

Run: `python3 GPCompaReports_v2/generate_site.py --output /tmp/v3rm --only par2_human`
Then: `grep -c "v3-secnav\|v3-nav.js\|v3-deeplink.js" /tmp/v3rm/reports/par2_human.html`
Expected: the nav and deep-link assets are still referenced. And `grep -c "v3-analysis" /tmp/v3rm/reports/par2_human.html` returns 0.

- [ ] **Step 9: Run the suite**

Run: `python3 -m pytest GPCompaReports_v2/tests -v`
Expected: green. The total drops by however many tests lived in `test_receptor_profile.py`; record the before and after counts in your report rather than assuming a number. `test_freeze.py` must still pass, which is the one that matters here: rebuild the demo first with `bash scripts/build_demo.sh` so it compares against a current build.

- [ ] **Step 10: Commit**

```bash
git add -A GPCompaReports_v2/
git commit -m "report: remove the opt-in analysis layer"
```

---

### Task 2: Landing hero, name, logo, and a dataset table replacing the stat band

**Files:**
- Modify: `GPCompaReports_v2/templates/landing.html` (hero and stat band)
- Modify: `GPCompaReports_v2/static/css/landing.css`
- Modify: `GPCompaReports_v2/static/js/landing.js`
- Modify: `GPCompaReports_v2/website/page_generators/landing_page.py`
- Create: `GPCompaReports_v2/tests/test_landing_counts.py`

**Interfaces:**
- Produces: a hero carrying the name, the factual description and the receptor scope; the search box; and a compact unanimated "Dataset and analysis" table in place of the stat band.
- Produces landing context values `n_contact_records=213456` and `n_threshold_changes=23025`, computed at build time, not hardcoded in the template.

**The stat band goes away entirely.** Do not replace it with another four-tile display. The four animated tiles were visual promotion; what replaces them is a documentation table whose purpose is cohort transparency. Ligand types leave the headline counts altogether.

**Background on the counts, so you fix the right thing.** The live page appears to "flash" wrong values: 51 GPCRs, 11 families, 2 ligand types, 38K contacts, before settling on 283, 60, 10, 213K. That is not stale data and there is no second source. `GPCompaReports_v2/static/js/landing.js` animates each stat from zero to its `data-target` over 1200ms with an easeOutQuart curve, and those four numbers are all about 18% of their targets, which is one frame of that animation about 58ms in. The values were never wrong.

It is still worth removing. Animating a scientific count means the page deliberately displays incorrect numbers for a second, and a screenshot taken for a talk can catch one. The code already honours `prefers-reduced-motion`, so some readers never see the effect at all, which makes it decorative rather than functional.

- [ ] **Step 1: Write the failing test**

Create `GPCompaReports_v2/tests/test_landing_counts.py`:

```python
"""The landing page's scientific counts must be the frozen dataset's, rendered
server side, never animated from a placeholder.

Every intermediate frame of a count-up animation is a wrong number on a page whose
whole purpose is to be quoted, so the animation is gone and this pins it.
"""
from pathlib import Path

import pytest

BASE = Path(__file__).resolve().parent.parent
INDEX = BASE / 'output_v3_demo' / 'index.html'
LANDING_JS = BASE / 'static' / 'js' / 'landing.js'

# Verified against the build on 2026-08-25. If a rebuild changes any of these,
# the dataset changed and that is a release event, not a test to update.
FROZEN = {
    '283': 'receptors',
    '60': 'receptor families',
    '213,456': 'contact-pair records',
    '23,025': 'threshold-passing changes',
    '566': 'models',
}


@pytest.mark.skipif(not INDEX.exists(), reason='needs a demo build on disk')
def test_frozen_counts_appear_verbatim():
    html = INDEX.read_text(encoding='utf-8')
    for value, what in FROZEN.items():
        assert value in html, 'missing %s count: %s' % (what, value)


@pytest.mark.skipif(not INDEX.exists(), reason='needs a demo build on disk')
def test_the_promotional_stat_band_is_gone():
    html = INDEX.read_text(encoding='utf-8')
    assert 'stat-band' not in html, 'the stat band survived'
    assert 'stat-num' not in html, 'stat tiles survived'
    assert 'data-target' not in html, 'an animated counter survived'


@pytest.mark.skipif(not INDEX.exists(), reason='needs a demo build on disk')
def test_the_limitation_statement_is_present():
    html = INDEX.read_text(encoding='utf-8')
    assert 'model-derived structural candidates' in html
    assert 'not direct experimental evidence of function' in html


def test_no_count_up_animation_remains():
    js = LANDING_JS.read_text(encoding='utf-8')
    assert 'data-target' not in js, 'the count-up animation is still wired up'
    assert 'animateCount' not in js, 'the count-up animation function survives'
```

- [ ] **Step 2: Run it and watch it fail**

Run: `python3 -m pytest GPCompaReports_v2/tests/test_landing_counts.py -v`
Expected: `test_no_count_up_animation_remains` FAILS on `data-target`.

- [ ] **Step 3: Compute the two new counts at build time**

In `GPCompaReports_v2/website/page_generators/landing_page.py`, add to `generate_landing_page`, before the render call:

```python
    from ..page_generators import gpcr_report_helpers as helpers

    n_contact_records = 0
    n_threshold_changes = 0
    for gid in store.gpcr_ids:
        delta_df = store.delta_data.get(gid)
        if delta_df is None or delta_df.empty:
            continue
        n_contact_records += len(delta_df)
        thr = helpers._calc_significance_threshold(delta_df)
        n_threshold_changes += int((delta_df['abs_delta'] >= thr).sum())
```

Pass `n_contact_records=n_contact_records`, `n_threshold_changes=n_threshold_changes`,
`n_models=2 * len(store.gpcr_ids)`, and `n_excluded=len(store.metadata) - len(store.gpcr_ids)`
into `template.render(...)`. Drop `total_contacts_k` and `n_ligand_types` if nothing else uses them.

These must be computed, not hardcoded. If a future rebuild changes them, the page should change
with the data rather than silently disagree with it.

- [ ] **Step 4: Delete the stat band and add the dataset table**

In `GPCompaReports_v2/templates/landing.html`, delete the entire `<section class="stat-band" ...>` block. Replace it with:

```html
<section class="dataset" aria-labelledby="dataset-heading">
  <div class="section-inner">
    <header class="section-head">
      <h2 id="dataset-heading">Dataset and analysis</h2>
    </header>
    <table class="dataset-table">
      <tbody>
        <tr><th scope="row">Receptors</th><td>{{ total_gpcrs }}</td>
            <td class="dataset-note">human Class A GPCRs with a complete active and inactive pair</td></tr>
        <tr><th scope="row">Receptor families</th><td>{{ n_families }}</td>
            <td class="dataset-note"></td></tr>
        <tr><th scope="row">Structural models</th><td>{{ n_models }}</td>
            <td class="dataset-note">two AlphaFold-Multistate models per receptor</td></tr>
        <tr><th scope="row">Contact-pair records</th><td>{{ "{:,}".format(n_contact_records) }}</td>
            <td class="dataset-note">residue pairs scored in both states</td></tr>
        <tr><th scope="row">Threshold-passing changes</th><td>{{ "{:,}".format(n_threshold_changes) }}</td>
            <td class="dataset-note">contacts exceeding their own receptor's threshold</td></tr>
        <tr><th scope="row">Excluded</th><td>{{ n_excluded }}</td>
            <td class="dataset-note">receptors without a complete active and inactive pair</td></tr>
      </tbody>
    </table>
    <p class="dataset-limitation">
      Positions and contact pairs highlighted by GPCompaRe should be interpreted as
      model-derived structural candidates, not direct experimental evidence of function.
    </p>
  </div>
</section>
```

The limitation statement sits here, directly under the counts, which is the most prominent
place a reader will meet it.

- [ ] **Step 5: Delete the animation**

In `GPCompaReports_v2/static/js/landing.js`, delete the entire "Stat count-up on scroll-into-view" block: `statNums`, `easeOutQuart`, `animateCount`, and the `IntersectionObserver` that drives them, including the non-observer fallback. Leave the search teaser code untouched. If `reduceMotion` is now unused, delete it too; check first.

- [ ] **Step 6: Style the dataset table**

Append to `GPCompaReports_v2/static/css/landing.css`:

```css
/* ----------- Dataset and analysis ----------- */
.dataset { padding: var(--s-6) 0; border-bottom: 1px solid var(--borderColor-muted); }
.dataset-table { border-collapse: collapse; font-size: 14.5px; margin: 0 0 var(--s-4); }
.dataset-table th { text-align: left; font-weight: 600; padding: 6px var(--s-5) 6px 0; white-space: nowrap; }
.dataset-table td { padding: 6px var(--s-5) 6px 0; font-variant-numeric: tabular-nums; }
.dataset-table .dataset-note { color: var(--fgColor-muted); font-size: 13px; font-variant-numeric: normal; }
.dataset-limitation { max-width: 70ch; color: var(--fgColor-muted); font-size: 14px; margin: 0; }
@media (max-width: 700px) { .dataset-table .dataset-note { display: none; } }
```

- [ ] **Step 7: Rewrite the hero**

In `GPCompaReports_v2/templates/landing.html`, replace the `.hero-copy` block's title and tagline with:

```html
      <h1 class="hero-title"><span class="hero-title-serif">GPCompaRe database</span></h1>
      <p class="hero-tagline">
        Paired active and inactive AlphaFold-Multistate models for {{ total_gpcrs }} human
        Class&nbsp;A GPCRs, compared with residue-residue contact scores (RRCS). Each report
        links the structural result to GPCRdb generic numbering, Core Functional Residues,
        conservation, gnomAD variation and AlphaMissense predictions.
      </p>
```

The hero carries the receptor scope only. Every other count now lives in the dataset table.
Keep the two call-to-action buttons and the search block exactly as they are. The search must stay prominent.

- [ ] **Step 8: Make the logo smaller**

In `GPCompaReports_v2/static/css/landing.css`, find the `.hero-emblem` rule and reduce the rendered size so the title, tagline and search carry more weight. Do not touch `GPCompaReports_v2/static/img/lablogo.png` itself. Add a comment saying the asset is unchanged and only its rendered size differs.

- [ ] **Step 9: Rebuild and check every number against the frozen dataset**

Run: `bash scripts/build_demo.sh`, then:

```bash
grep -o '283\|60\|566\|213,456\|23,025\|>4<' GPCompaReports_v2/output_v3_demo/index.html | sort | uniq -c
grep -c 'stat-band\|stat-num\|data-target' GPCompaReports_v2/output_v3_demo/index.html
```

Expected: all six frozen figures present, and 0 hits for the stat band and any animated counter.
If a computed count disagrees with 213,456 or 23,025, stop and report it: the dataset changed and
that is a release event, not a number to quietly update.

- [ ] **Step 10: Run the tests**

Run: `python3 -m pytest GPCompaReports_v2/tests -v`
Expected: green, including all five new tests.

- [ ] **Step 11: Commit**

```bash
git add -A GPCompaReports_v2/
git commit -m "landing: factual hero, dataset table replacing the stat band"
```

---

### Task 3: Remove the ranking, the schematic, and the old About and Data Sources sections

**Files:**
- Modify: `GPCompaReports_v2/templates/landing.html`
- Modify: `GPCompaReports_v2/static/css/landing.css`
- Modify: `GPCompaReports_v2/website/page_generators/landing_page.py`

Do not replace the ranking with a normalised one. A new ranking is a new analysis and would need its own validation.

- [ ] **Step 1: Delete the featured section**

In `GPCompaReports_v2/templates/landing.html`, delete the entire `<section class="featured" ...>` block, headed "Top Conformational Changes", including the card loop and the sparkline SVG.

- [ ] **Step 2: Delete the explainer schematic**

Delete the `<section class="explainer" ...>` block added earlier, including its two-circle SVG with the invented `RRCS 0.4`, `RRCS 2.9` and `ΔRRCS +2.5` values. Task 4 replaces it with a real worked example.

- [ ] **Step 3: Delete the old About and Data Sources sections**

`GPCompaReports_v2/templates/landing.html` still carries a `<section class="about" id="about" ...>` block headed "About GPCompaRe" and a `<section class="sources" ...>` block headed "Data Sources" with a grid of source cards. Task 6 adds a Methods and data sources section that covers the same ground properly.

Delete both old sections now. The new sections must **replace** them, not sit alongside them: two accounts of the same thing on one page is how a page starts contradicting itself.

Then check nothing links to the removed anchor:

```bash
grep -rn 'href="#about"\|#about' GPCompaReports_v2/templates/
```

If the topbar or footer links to `#about`, remove that link too.

- [ ] **Step 4: Drop the now-unused context**

In `GPCompaReports_v2/website/page_generators/landing_page.py`, remove `top5`, `sparkline_bins`, `delta_clip`, the `bin_signed_delta` call, and the `SPARKLINE_BINS` and `DELTA_CLIP` constants if nothing else uses them. Keep `bin_signed_delta` itself only if another module imports it; check with grep before deleting.

- [ ] **Step 5: Remove the dead styles**

In `GPCompaReports_v2/static/css/landing.css`, delete the `.featured*`, `.feature-card*`, `.feature-gene`, `.feature-meta`, `.feature-delta*`, `.sparkline`, `.feature-footer`, `.feature-stat`, `.feature-open` rules, the `.explainer*` and `.ex-*` rules, and the `.about*` and `.sources*` / `.source-card*` / `.source-icon` rules.

Keep `.ex-swatch`, `.ex-swatch.is-active` and `.ex-swatch.is-inactive` if you have already added Task 4's styles, since the worked example reuses those class names for its colour key. If Task 4 has not run yet, delete them and let Task 4 reintroduce them.

- [ ] **Step 6: Rebuild and confirm**

Run: `bash scripts/build_demo.sh`
Then: `grep -c 'class="featured"\|sparkline\|class="about"\|class="sources"' output_v3_demo/index.html`
Expected: 0.

- [ ] **Step 7: Run the suite and commit**

```bash
python3 -m pytest GPCompaReports_v2/tests -v
git add -A GPCompaReports_v2/
git commit -m "landing: drop the ranking, the schematic, and the old About and Sources sections"
```

---

### Task 4: Real PAR2 contact explainer

**Files:**
- Create: `GPCompaReports_v2/static/img/par2-explainer-snake.png`
- Modify: `GPCompaReports_v2/templates/landing.html`, `GPCompaReports_v2/static/css/landing.css`

**The contact is already selected and verified.** Do not choose a different one and do not recompute the values from scratch. These figures were read out of `par2_human_rrcs_delta.csv` and match the paper's Figure 4d exactly:

| field | value |
|---|---|
| Residue 1 | M159, `3.36x36`, TM3 |
| Residue 2 | F300, `6.48x48`, TM6 |
| Active RRCS | 0.6651 |
| Inactive RRCS | 5.6580 |
| ΔRRCS | −4.9928 |
| PAR2 threshold | 2.602, so this contact passes |
| Rank by \|Δ\| | 27 of 776, so it is not the largest change |

M159 is the methionine of the xMY motif. F300 occupies the corresponding aromatic 6.48 toggle-switch position. Both sit in transmembrane helices. Display values rounded to two decimals: active 0.67, inactive 5.66, ΔRRCS −4.99.

**Wording that must not drift.** Do not write that F300 "is part of CWxP" or "is the CWxP toggle". State only that it occupies the 6.48 toggle-switch position. Naming the motif asserts a sequence identity this page has not established for PAR2.

**The deep link cannot select this contact, and the copy must not imply it does.** The hash grammar supports `view`, `min`, `max`, `dir` and `sec` only. Per-residue selection exists inside the snake plot but was never exposed to the URL, and a magnitude filter of `min=4.9` matches 29 PAR2 contacts rather than one. So the link opens the snake plot at the ΔRRCS view and is labelled exactly that: "Open the PAR2 snake plot". Do not label it "open this contact" or "see this pair".

- [ ] **Step 1: Verify the numbers yourself before rendering them**

Run:

```bash
python3 - <<'PY'
import pandas as pd
b='The_batch_RRCS_analyzer/batch_analysis_full/batch_analysis_20260202_151051/csv_data'
d=pd.read_csv(f'{b}/par2_human_rrcs_delta.csv')
r=d[((d.res1==159)&(d.res2==300))|((d.res1==300)&(d.res2==159))]
print(r[['res1','res2','active_rrcs','inactive_rrcs','delta_rrcs','abs_delta']].to_string(index=False))
PY
```

Expected: active 0.6651, inactive 5.6580, delta −4.9928. If any figure differs, stop and report it rather than publishing a number that does not match the database.

- [ ] **Step 2: Crop the snake plot fragment**

The built PAR2 report contains the snake plot as inline SVG. Render the region around residues 159, 160 and 300 to a PNG for the landing page:

```bash
python3 GPCompaReports_v2/generate_site.py --output /tmp/par2crop --only par2_human
```

Then extract the `<svg>` from `/tmp/par2crop/reports/par2_human.html` (it lives inside `#snake-plot-container`), write it to a standalone `.svg` file, and rasterise a crop containing TM3 and TM6 with `chromium --headless --screenshot`. Chromium here is a snap package and cannot read `/tmp`, so stage the file inside the project directory, for example `GPCompaReports_v2/static/img/_work/`, and delete the staging directory afterwards.

Save the result as `GPCompaReports_v2/static/img/par2-explainer-snake.png`. The crop must show M159 and F300 and be legible at roughly 480px wide.

If the crop cannot be made legible, stop and report it. Do not draw a substitute.

- [ ] **Step 3: Add the section**

In `GPCompaReports_v2/templates/landing.html`, immediately after the hero and stat band, insert:

```html
<section class="worked" aria-labelledby="worked-heading">
  <div class="section-inner">
    <header class="section-head">
      <h2 id="worked-heading">How to read a GPCompaRe contact</h2>
      <p class="section-sub">One real contact from the PAR2 report, from the topology diagram to the number.</p>
    </header>

    <div class="worked-grid">
      <figure class="worked-fig">
        <img src="{{ static_prefix }}static/img/par2-explainer-snake.png"
             alt="Fragment of the PAR2 snake plot showing residues M159 in TM3 and F300 in TM6">
        <figcaption>PAR2 (F2RL1), transmembrane helices 3 and 6.</figcaption>
      </figure>

      <div class="worked-copy">
        <p>
          RRCS assigns a contact score to a residue pair from the interatomic distances within a
          structure. The same pair is scored in the active model and in the inactive model, and
          the difference is ΔRRCS: <strong>active minus inactive</strong>.
        </p>

        <table class="worked-table">
          <caption>M159 and F300 as they appear in the PAR2 report</caption>
          <tbody>
            <tr><th scope="row">Residue 1</th><td>M159 <code>3.36x36</code>, TM3</td></tr>
            <tr><th scope="row">Residue 2</th><td>F300 <code>6.48x48</code>, TM6</td></tr>
            <tr><th scope="row">Active RRCS</th><td>0.67</td></tr>
            <tr><th scope="row">Inactive RRCS</th><td>5.66</td></tr>
            <tr><th scope="row">ΔRRCS</th><td class="is-inactive">−4.99</td></tr>
            <tr><th scope="row">PAR2 threshold</th><td>2.60, so this contact passes</td></tr>
          </tbody>
        </table>

        <p class="worked-key">
          The contact is stronger in the inactive model, so it is
          <span class="ex-swatch is-inactive"></span> <strong>blue</strong>. Contacts stronger in
          the active model are <span class="ex-swatch is-active"></span> <strong>red</strong>.
        </p>

        <p>
          M159 is the methionine of the xMY motif. F300 occupies the corresponding aromatic
          6.48 toggle-switch position. Both sit in transmembrane helices.
        </p>

        <p>
          <a href="{{ static_prefix }}reports/par2_human.html#view=delta_rrcs&amp;sec=sec-snake">
            Open the PAR2 snake plot</a>
        </p>
      </div>
    </div>
  </div>
</section>
```

- [ ] **Step 4: Style it**

Append to `GPCompaReports_v2/static/css/landing.css`:

```css
/* ----------- Worked example ----------- */
.worked { padding: var(--s-7) 0 var(--s-6); border-bottom: 1px solid var(--borderColor-muted); }
.worked-grid { display: grid; grid-template-columns: minmax(0, 480px) 1fr; gap: var(--s-6); align-items: start; }
.worked-fig { margin: 0; }
.worked-fig img { width: 100%; height: auto; border: 1px solid var(--borderColor-muted); border-radius: var(--r-md); }
.worked-fig figcaption { font-size: 12.5px; color: var(--fgColor-muted); margin-top: var(--s-2); }
.worked-copy p { margin: 0 0 var(--s-3); max-width: 60ch; }
.worked-table { border-collapse: collapse; margin: 0 0 var(--s-4); font-size: 14px; }
.worked-table caption { text-align: left; font-size: 12.5px; color: var(--fgColor-muted); padding-bottom: var(--s-2); }
.worked-table th { text-align: left; font-weight: 600; padding: 5px 14px 5px 0; color: var(--fgColor-muted); white-space: nowrap; }
.worked-table td { padding: 5px 0; font-variant-numeric: tabular-nums; }
.worked-table code { font-family: var(--font-mono); font-size: 12.5px; background: var(--bgColor-neutral-muted); border-radius: 4px; padding: 1px 5px; }
.worked-table .is-inactive { color: var(--semantic-inactive); font-weight: 700; }
.worked-key { display: flex; align-items: center; flex-wrap: wrap; gap: 6px; }
.ex-swatch { display: inline-block; width: 12px; height: 12px; border-radius: 3px; }
.ex-swatch.is-active { background: var(--semantic-active); }
.ex-swatch.is-inactive { background: var(--semantic-inactive); }
@media (max-width: 820px) { .worked-grid { grid-template-columns: 1fr; } }
```

- [ ] **Step 5: Verify the deep link resolves**

Run: `bash scripts/build_demo.sh`, then confirm `output_v3_demo/reports/par2_human.html` exists and that `sec-snake` is a real id in it: `grep -c 'id="sec-snake"' output_v3_demo/reports/par2_human.html` returns 1.

- [ ] **Step 6: Run the suite and commit**

```bash
python3 -m pytest GPCompaReports_v2/tests -v
git add -A GPCompaReports_v2/
git commit -m "landing: real PAR2 worked example replacing the schematic"
```

---

### Task 5: Core Functional Residues section

**Files:**
- Create: `GPCompaReports_v2/static/img/fig-cfr-topology.png`
- Modify: `GPCompaReports_v2/templates/landing.html`, `GPCompaReports_v2/static/css/landing.css`

The visual is panel **a** of Figure 2 from the manuscript's own figure set, at `/home/yamir/.claude/uploads/7fdf1d09-11cb-4b9c-957a-047980866073/c637be9e-Draft17_Main_Figures.pdf`, page 2. It is a consensus Class A topology coloured by how many of the 283 receptors show a threshold-passing change at each generic position, with DRY, NPxxY, CWxP, PIF and D2.50x50 circled as landmarks. It is the lab's own figure, so there is no licensing question, and using it means the site and the paper show the reader the same picture.

- [ ] **Step 1: Extract the figure**

```bash
pdftoppm -f 2 -l 2 -r 200 -png \
  "/home/yamir/.claude/uploads/7fdf1d09-11cb-4b9c-957a-047980866073/c637be9e-Draft17_Main_Figures.pdf" \
  GPCompaReports_v2/static/img/_fig2
```

That writes `_fig2-2.png`. Crop panel **a** (the topology diagram and its colour legend, excluding panels b and c) and save as `GPCompaReports_v2/static/img/fig-cfr-topology.png`. Delete the intermediate file. Keep the legend: without it the colours mean nothing.

- [ ] **Step 2: Add the section**

Insert after the worked-example section in `GPCompaReports_v2/templates/landing.html`:

```html
<section class="cfr-intro" aria-labelledby="cfr-heading">
  <div class="section-inner">
    <header class="section-head">
      <h2 id="cfr-heading">Core Functional Residues across Class A GPCRs</h2>
    </header>

    <div class="cfr-grid">
      <figure class="cfr-fig">
        <img src="{{ static_prefix }}static/img/fig-cfr-topology.png"
             alt="Consensus Class A topology with each generic position shaded by the number of receptors showing a threshold-passing contact change, and the DRY, NPxxY, CWxP, PIF and D2.50x50 landmarks marked">
        <figcaption>
          Each position is shaded by how many of the {{ total_gpcrs }} receptors show a
          threshold-passing change there. Lettered and outlined positions mark canonical
          activation landmarks.
        </figcaption>
      </figure>

      <div class="cfr-copy">
        <p>
          A receptor-level Core Functional Residue (CFR) is a residue that participates in at
          least one contact whose absolute ΔRRCS passes that receptor's own threshold.
        </p>
        <p>
          Mapping those receptor-level CFRs onto GPCRdb generic positions gives recurrent CFR
          positions across the family. A position is counted as recurrent when it is identified as
          a CFR in at least three receptors.
        </p>
        <p>
          No motif was selected in advance. The cross-receptor analysis recovers established parts
          of the DRY, NPxxY and CWxP activation machinery on its own.
        </p>
        <p>
          <a href="{{ static_prefix }}statistics.html">See the ranked Top 30 recurrent positions</a>
        </p>
      </div>
    </div>
  </div>
</section>
```

- [ ] **Step 3: Style it**

Append to `GPCompaReports_v2/static/css/landing.css`:

```css
/* ----------- CFR intro ----------- */
.cfr-intro { padding: var(--s-7) 0 var(--s-6); border-bottom: 1px solid var(--borderColor-muted); }
.cfr-grid { display: grid; grid-template-columns: minmax(0, 1fr) minmax(0, 380px); gap: var(--s-6); align-items: start; }
.cfr-fig { margin: 0; }
.cfr-fig img { width: 100%; height: auto; }
.cfr-fig figcaption { font-size: 12.5px; color: var(--fgColor-muted); margin-top: var(--s-2); max-width: 70ch; }
.cfr-copy p { margin: 0 0 var(--s-3); max-width: 60ch; }
@media (max-width: 820px) { .cfr-grid { grid-template-columns: 1fr; } }
```

- [ ] **Step 4: Check the terminology**

Run: `grep -rn "Top 50" GPCompaReports_v2/templates GPCompaReports_v2/static`
Expected: no output. The string must not appear in user-facing copy.

- [ ] **Step 5: Rebuild, run the suite, commit**

```bash
bash scripts/build_demo.sh
python3 -m pytest GPCompaReports_v2/tests -v
git add -A GPCompaReports_v2/
git commit -m "landing: Core Functional Residues section using the paper's topology figure"
```

---

### Task 6: Reports, custom analysis, methods and citation sections

**Files:**
- Modify: `GPCompaReports_v2/templates/landing.html`, `GPCompaReports_v2/static/css/landing.css`

- [ ] **Step 1: Add the three remaining sections**

Insert after the CFR section, in this order:

```html
<section class="contains" aria-labelledby="contains-heading">
  <div class="section-inner">
    <header class="section-head">
      <h2 id="contains-heading">What each receptor report contains</h2>
    </header>
    <ul class="contains-list">
      <li>An interactive snake plot with selectable colour views and contact arcs</li>
      <li>Ranked contact-pair results with active RRCS, inactive RRCS and ΔRRCS, up to 1,000 ranked contact pairs</li>
      <li>The receptor's own threshold, and which contacts pass it</li>
      <li>Receptor-level Core Functional Residues</li>
      <li>GPCRdb generic numbering</li>
      <li>Sequence conservation</li>
      <li>gnomAD missense variation at contact-forming residues</li>
      <li>AlphaMissense pathogenicity predictions</li>
      <li>Ranked tables with CSV export</li>
    </ul>
    <p class="contains-note">
      Use the search above to open a receptor. Example reports:
      <a href="{{ static_prefix }}reports/adrb2_human.html">ADRB2</a>,
      <a href="{{ static_prefix }}reports/par2_human.html">PAR2</a>.
    </p>
  </div>
</section>

<section class="custom" aria-labelledby="custom-heading">
  <div class="section-inner">
    <header class="section-head">
      <h2 id="custom-heading">Analysing your own structures</h2>
    </header>
    <p>
      GPCompaRe also includes a program for analysing your own paired active and inactive
      structures. It applies the same RRCS comparison used to build this database and writes
      standalone HTML and CSV output.
    </p>
    <p>
      Analyses you run locally stay local. They are not added to the public database.
    </p>
    <p><a href="{{ static_prefix }}downloads.html">Downloads</a></p>
  </div>
</section>

<section class="methods-sources" aria-labelledby="methods-heading">
  <div class="section-inner">
    <header class="section-head">
      <h2 id="methods-heading">Methods and data sources</h2>
    </header>
    <p>
      Structures are paired active-state and inactive-state models from
      <a href="https://github.com/huhlim/alphafold-multistate" target="_blank" rel="noopener">AlphaFold-Multistate</a>.
      Contacts are scored with
      <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC6954041/" target="_blank" rel="noopener">RRCS</a>.
      Generic numbering and segment boundaries come from
      <a href="https://gpcrdb.org/" target="_blank" rel="noopener">GPCRdb</a>.
      Population variation comes from
      <a href="https://gnomad.broadinstitute.org/" target="_blank" rel="noopener">gnomAD v4</a>,
      pathogenicity predictions from
      <a href="https://alphamissense.hegelab.org/" target="_blank" rel="noopener">AlphaMissense</a>,
      and conservation from
      <a href="https://www.ebi.ac.uk/ProtVar/" target="_blank" rel="noopener">ProtVar</a>.
    </p>
    <p>
      Each receptor has its own ΔRRCS threshold, computed as max(mean(|Δ|) + σ, 0.2).
    </p>
    <p>
      {{ n_models }} models were analysed, two per receptor. Four further receptors listed in the
      Class&nbsp;A reference set are excluded because they lack a complete active and inactive
      pair: GPR26, LGR5, OPSX and RXFP1.
    </p>
    <p>
      Positions and contact pairs highlighted by GPCompaRe should be interpreted as model-derived
      structural candidates, not direct experimental evidence of function.
    </p>
  </div>
</section>

<section class="release" aria-labelledby="release-heading">
  <div class="section-inner">
    <header class="section-head">
      <h2 id="release-heading">Version, downloads and citation</h2>
    </header>
    <dl class="release-list">
      <dt>Release</dt><dd>GPCompaRe database release 2026.08</dd>
      <dt>DOI</dt><dd>DOI pending release</dd>
      <dt>Downloads</dt><dd><a href="{{ static_prefix }}downloads.html">Database release and analysis software</a></dd>
      <dt>Contact</dt><dd><a href="{{ static_prefix }}contact.html">Contact the lab</a></dd>
    </dl>
    <p class="release-scope">
      GPCompaRe was developed as the companion resource for the associated study. Corrections and
      substantial updates will be documented as new releases. No fixed update schedule is
      currently planned.
    </p>
  </div>
</section>
```

- [ ] **Step 2: Style them**

Append to `GPCompaReports_v2/static/css/landing.css`:

```css
/* ----------- Contains, custom, methods, release ----------- */
.contains, .custom, .methods-sources, .release { padding: var(--s-7) 0 var(--s-6); border-bottom: 1px solid var(--borderColor-muted); }
.release { border-bottom: 0; }
.contains-list { columns: 2; column-gap: var(--s-6); margin: 0 0 var(--s-4); padding-left: 1.1rem; }
.contains-list li { margin-bottom: var(--s-2); break-inside: avoid; }
.contains-note, .custom p, .methods-sources p { max-width: 70ch; margin: 0 0 var(--s-3); }
.release-list { display: grid; grid-template-columns: max-content 1fr; gap: var(--s-2) var(--s-5); margin: 0 0 var(--s-4); }
.release-list dt { font-weight: 600; color: var(--fgColor-muted); }
.release-list dd { margin: 0; }
.release-scope { max-width: 70ch; color: var(--fgColor-muted); font-size: 14px; }
@media (max-width: 820px) { .contains-list { columns: 1; } }
```

- [ ] **Step 3: Confirm the section order matches the spec**

Run:

```bash
grep -o 'class="\(hero\|dataset\|worked\|cfr-intro\|contains\|custom\|methods-sources\|release\|featured\|explainer\|about\|sources\|stat-band\)"' \
  GPCompaReports_v2/templates/landing.html
```

Expected, in this order and nothing else: `hero`, `dataset`, `worked`, `cfr-intro`, `contains`, `custom`, `methods-sources`, `release`.

If `featured`, `explainer`, `about`, `sources` or `stat-band` still appears, Task 2 or Task 3 left something behind. Fix that before continuing rather than adding sections on top of it.

- [ ] **Step 4: Rebuild, run the suite, commit**

```bash
bash scripts/build_demo.sh
python3 -m pytest GPCompaReports_v2/tests -v
git add -A GPCompaReports_v2/
git commit -m "landing: reports, custom analysis, methods and release sections"
```

---

### Task 7: Contact page

**Files:**
- Create: `GPCompaReports_v2/templates/contact.html`, `GPCompaReports_v2/website/page_generators/contact_page.py`
- Create: `docs/superpowers/EXTERNAL_SETUP.md`
- Modify: `GPCompaReports_v2/website/site_generator.py`, `GPCompaReports_v2/templates/_partials/topbar.html`

**Interfaces:**
- Produces: `generate_contact_page(env, store, output_dir)` writing `output/contact.html`.

**Do not invent** a Formspree endpoint, account id, reCAPTCHA site key, secret, or any personal name, role or email beyond `tggrlab@gmail.com`. Use a clearly marked placeholder and document the setup.

- [ ] **Step 1: Write the generator**

Create `GPCompaReports_v2/website/page_generators/contact_page.py`:

```python
"""Generate the Contact page.

The form posts to Formspree. The endpoint is a placeholder until the lab creates
the account, see docs/superpowers/EXTERNAL_SETUP.md. Nothing here invents a
credential.
"""

FORMSPREE_ENDPOINT_PLACEHOLDER = 'FORMSPREE_ENDPOINT_NOT_YET_CONFIGURED'
CONTACT_EMAIL = 'tggrlab@gmail.com'


def generate_contact_page(env, store, output_dir):
    template = env.get_template('contact.html')
    html = template.render(
        static_prefix='',
        active_page='contact',
        nav_home_url='index.html',
        nav_browse_url='browse/index.html',
        nav_stats_url='statistics.html',
        nav_downloads_url='downloads.html',
        nav_contact_url='contact.html',
        page_title='Contact · GPCompaRe',
        total_gpcrs=len(store.gpcr_ids),
        formspree_endpoint=FORMSPREE_ENDPOINT_PLACEHOLDER,
        contact_email=CONTACT_EMAIL,
    )
    (output_dir / 'contact.html').write_text(html, encoding='utf-8')
    print("  Generated: contact.html")
```

- [ ] **Step 2: Write the template**

Create `GPCompaReports_v2/templates/contact.html` extending `base.html`, with a `content` block containing a form whose `action` is `{{ formspree_endpoint }}` and `method="POST"`. Required fields: Name (`name`), Email (`email`, `type="email"`), Institution (`institution`), Type of enquiry (`enquiry_type`, a `<select>`), Message (`message`, a `<textarea>`). Optional: Receptor or gene (`receptor`), Relevant page URL (`page_url`, `type="url"`). Every field needs a real `<label for>`. Required fields carry `required`. No `<input type="file">` anywhere.

The enquiry options, exactly: Data or annotation issue; Website problem; Custom analysis program; Scientific question; Data, download or citation question; Other.

Include, as visible copy:

```html
      <p class="contact-note">
        GPCompaRe is an academic resource maintained on a best-effort basis.
      </p>
      <p class="contact-note">
        Do not submit patient-identifiable or other sensitive information through this form.
      </p>
```

Include a lab block naming the TGGR Laboratory and `{{ contact_email }}`, with a clearly marked gap for the scientific and technical contact names. Do not fill those in.

When the endpoint is unconfigured, show a notice **and disable the submit button**, so no placeholder submission is possible:

```html
      {% set unconfigured = formspree_endpoint.startswith('FORMSPREE_') %}
      {% if unconfigured %}
      <p class="contact-warning" role="status">
        This form is not connected yet. See the setup notes before publishing this page.
      </p>
      {% endif %}
      <button type="submit" class="btn btn-primary"{% if unconfigured %} disabled aria-disabled="true"{% endif %}>
        Send
      </button>
```

A disabled button is not a security control, it is honesty: an enabled button on an unwired form invites someone to type a real question and lose it.

- [ ] **Step 3: Wire it into the pipeline and the nav**

In `GPCompaReports_v2/website/site_generator.py`, import `generate_contact_page` and call it after the statistics page. In `GPCompaReports_v2/templates/_partials/topbar.html`, add Downloads and Contact links to both the desktop nav and the mobile drawer, using `nav_downloads_url` and `nav_contact_url`. Every existing page generator must now pass those two variables; update `landing_page.py`, `gpcr_index.py`, `statistics_page.py` and `gpcr_report_page.py` accordingly, remembering that report pages need the `../` prefix.

**On the shared topbar and the freeze rule, already checked.** The topbar is shared with report pages, so this adds two links to a PI-approved page. It does not break the freeze test: that test inspects only the content inside `<section class="report-section">` elements, and the topbar sits well before the first of them, verified in the built PAR2 page. Top-level navigation is wayfinding, the same category as the sticky section nav this branch already added to report pages.

Confirm it rather than trusting this note. After Step 5, run `python3 -m pytest GPCompaReports_v2/tests/test_freeze.py -v` and check it is green. **If the freeze test fails, do not relax it.** Stop, revert the topbar change, and report it: the fallback is to leave the report topbar as it was and surface Downloads and Contact only on non-report pages and in the footer.

- [ ] **Step 4: Document the external setup**

Create `docs/superpowers/EXTERNAL_SETUP.md` listing what a human must do before this page works: create the Formspree account and paste the endpoint into `contact_page.py`; obtain reCAPTCHA v3 site and secret keys; restrict the form to the GPCompaRe domain; supply the scientific and technical contact names and roles; confirm the data and software licenses; create the Zenodo deposit and record the DOI. State plainly that none of these are configured.

- [ ] **Step 5: Verify**

Run: `bash scripts/build_demo.sh`, then confirm `output_v3_demo/contact.html` exists, contains all five required fields and both optional ones, contains no `<input type="file">`, and shows the not-connected notice.

- [ ] **Step 6: Run the suite and commit**

```bash
python3 -m pytest GPCompaReports_v2/tests -v
git add -A GPCompaReports_v2/ docs/
git commit -m "site: Contact page with an unconfigured Formspree form"
```

---

### Task 8: Downloads page

**Files:**
- Create: `GPCompaReports_v2/templates/downloads.html`, `GPCompaReports_v2/website/page_generators/downloads_page.py`
- Create: `GPCompaReports_v2/tests/test_new_pages.py`
- Modify: `GPCompaReports_v2/website/site_generator.py`

**Do not create fake download links.** Neither archive is ready. The database half is not packaged and its third-party annotation licensing is unresolved; the software half has no releasable user-facing program yet, confirmed by the project owner. Both are clearly labelled placeholders in this draft.

**The two placeholders are not equivalent, and the page must say so.** The software placeholder may stand until a user-facing program is genuinely releasable, which may be well after publication. The database placeholder may not: **the final public release must contain a real database archive**, because a paper companion that offers no data is not a companion. Record that distinction on the page and again in the handoff, so the database archive cannot quietly ship as a placeholder.

- [ ] **Step 1: Write the generator**

Create `GPCompaReports_v2/website/page_generators/downloads_page.py` with the signature
`generate_downloads_page(env, store, output_dir)`, matching `generate_contact_page` exactly so
`site_generator.py` calls both the same way. It renders `downloads.html` with the same nav
context variables as the contact page (`static_prefix`, `active_page='downloads'`,
`nav_home_url`, `nav_browse_url`, `nav_stats_url`, `nav_downloads_url`, `nav_contact_url`,
`total_gpcrs`), `page_title='Downloads · GPCompaRe'`, plus
`release_name='GPCompaRe database release 2026.08'`,
`software_name='GPCompaRe software v1.0.0'`, and `doi_status='DOI pending release'`.
It writes `output_dir / 'downloads.html'` and prints `"  Generated: downloads.html"`.

- [ ] **Step 2: Write the template**

Two clearly separated blocks.

**A, the database release.** State that it will contain: the canonical receptor-level files used to build the website; a receptor summary; cross-receptor results; a README; a data dictionary; the database version and release date; a changelog; source and provenance information; licenses and third-party notices; a file manifest; and checksums. Mark it "Not yet available" with no link. Add one line: packaging is pending a licensing review of the third-party annotations.

**B, the analysis software.** State that it will contain: the user-facing analysis program; installation instructions; the supported Python version; dependencies; one redistributable example; an example command; expected HTML and CSV outputs; known limitations; a software license; and a validation command. Mark it "Not yet available" with no link.

Add: "The website generator is available in the repository for reproducibility. It is not the analysis program."

Include the release name and `DOI pending release`.

- [ ] **Step 3: Write the test**

Create `GPCompaReports_v2/tests/test_new_pages.py`:

```python
"""The Contact and Downloads pages must not ship invented credentials or dead links."""
import re
from pathlib import Path

import pytest

BASE = Path(__file__).resolve().parent.parent
OUT = BASE / 'output_v3_demo'
CONTACT = OUT / 'contact.html'
DOWNLOADS = OUT / 'downloads.html'

pytestmark = pytest.mark.skipif(not OUT.exists(), reason='needs a demo build on disk')


def test_contact_form_has_every_required_field():
    html = CONTACT.read_text(encoding='utf-8')
    for field in ('name', 'email', 'institution', 'enquiry_type', 'message'):
        assert 'name="%s"' % field in html, 'missing field: %s' % field


def test_contact_form_accepts_no_file_uploads():
    html = CONTACT.read_text(encoding='utf-8')
    assert 'type="file"' not in html


def test_contact_page_invents_no_credentials():
    html = CONTACT.read_text(encoding='utf-8')
    assert 'formspree.io/f/' not in html, 'a Formspree endpoint was invented'
    assert not re.search(r'6L[\w-]{20,}', html), 'a reCAPTCHA site key was invented'


def test_downloads_page_offers_no_dead_links():
    html = DOWNLOADS.read_text(encoding='utf-8')
    assert '.zip' not in html and '.tar.gz' not in html, 'a download link was invented'
    assert 'DOI pending release' in html
    assert 'zenodo.org/record' not in html, 'a Zenodo DOI was invented'


def test_both_pages_are_reachable_from_the_nav():
    index = (OUT / 'index.html').read_text(encoding='utf-8')
    assert 'contact.html' in index
    assert 'downloads.html' in index
```

- [ ] **Step 4: Run it and watch it fail, then wire the generator in**

Run: `python3 -m pytest GPCompaReports_v2/tests/test_new_pages.py -v`
Expected: failures on the missing `downloads.html`. Then add `generate_downloads_page` to `site_generator.py` next to the contact page.

- [ ] **Step 5: Rebuild and confirm green**

```bash
bash scripts/build_demo.sh
python3 -m pytest GPCompaReports_v2/tests -v
```

- [ ] **Step 6: Commit**

```bash
git add -A GPCompaReports_v2/
git commit -m "site: Downloads page with labelled placeholders for both archives"
```

---

### Task 9: Screenshots and handoff

**Files:**
- Create: `docs/superpowers/V3_DRAFT_HANDOFF.md`
- Create: screenshots under `docs/superpowers/screenshots/`

Chromium and Firefox are installed, both as snap packages. Snap confinement cannot read `/tmp`, so render from inside the project directory.

- [ ] **Step 1: Capture the screenshots**

For each of the landing hero, the PAR2 worked example, the CFR section, the Contact page and the Downloads page, capture desktop (1440x900) and mobile (390x844):

```bash
chromium --headless --disable-gpu --hide-scrollbars \
  --window-size=1440,900 \
  --screenshot="docs/superpowers/screenshots/landing-desktop.png" \
  "file://$(pwd)/GPCompaReports_v2/output_v3_demo/index.html"
```

Repeat with `--window-size=390,844` for mobile and the appropriate `#` anchors for the sections. If chromium fails, try firefox with `--screenshot`. If neither works, report that rather than shipping the draft without visual evidence.

- [ ] **Step 2: Write the handoff**

Create `docs/superpowers/V3_DRAFT_HANDOFF.md` containing, in order:

1. Confirmation that nothing was merged, deployed or pushed
2. The selected PAR2 contact with its full verification table
3. A table of every scientific count on the landing page and the exact source of each
4. The complete user-facing copy added or changed
5. The list of files changed
6. Test results
7. Confirmation that the analysis layer is completely removed, with the grep that proves it
8. Confirmation that approved report content is unchanged by default, citing the freeze test
9. External setup still required, referencing `EXTERNAL_SETUP.md`
10. Any concern affecting scientific accuracy, reproducibility, or the manuscript-to-site match

Under item 10, record this as a **pre-release investigation item**: the manuscript's Figure 4 gives Y160 the generic number `3.37x37`, but PAR2's annotation file has no entry for position 160 at all, so the site cannot render that numbering. The M159-Y160 contact also falls below PAR2's threshold at 1.78 against 2.602, which matches Figure 4e.

Flag it for investigation before release. **Do not change the approved report or the manuscript to resolve it**, and do not patch the annotation file. It needs a human decision about which source is right.

Also under item 10: state that the database download must be a real archive before public release, and that only the software download may remain a placeholder past that point.

Also under item 10: the report section headed "Complete RRCS results" is capped at 1000 rows and is therefore incomplete for 7 of 283 receptors. The proposed correction, "RRCS results, up to 1,000 contact pairs", is **pending PI approval and deliberately not applied in this draft**.

- [ ] **Step 3: Final verification**

```bash
python3 -m pytest GPCompaReports_v2/tests -v
git status --short
git log --oneline 3a10a53..HEAD | head -20
```

Confirm: suite green, no uncommitted tracked changes, and no push has occurred.

- [ ] **Step 4: Commit**

```bash
git add -A docs/
git commit -m "docs: V3 draft handoff and screenshots"
```

---

## Out of scope for this draft

- Merging, deploying or pushing anything
- Changing the "Complete RRCS results" heading, which needs PI approval first
- Creating a Zenodo record or any DOI
- Any new normalised receptor ranking to replace the removed one
- Any change to report page content or its default presentation
- Re-adding the analysis layer in any form
