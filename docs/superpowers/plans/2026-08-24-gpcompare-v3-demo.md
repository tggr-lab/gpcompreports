# GPCompaRe V3 Demo Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a clickable V3 demo of the GPCompaRe site on five receptors, implementing the real V3 code rather than a throwaway mockup.

**Architecture:** Keep the existing Python + Jinja2 generator. Add one Python analysis module for the opt-in analysis layer, three vanilla-JS modules under `static/js/`, and a receptor selector so the generator can build a named subset. Every new reader-facing feature defaults to off and is remembered in `localStorage`, so a report with all toggles off renders the content the PI already approved.

**Tech Stack:** Python 3.13, Jinja2 3.1.6, pandas 3.0.0, Plotly 6.5.2, pytest 9.0.2, vanilla ES5-style JS (matching `site.js`), CSS custom properties.

## Global Constraints

Every task's requirements implicitly include this section.

- **Brand display text is `GPCompaRe database`** (short form `GPCompaRe`). Applies to templates, page titles, footer, and docstrings.
- **Never change any lowercase `gpcompreports` URL or host string.** Intentional survivors: the GitHub remote, `https://tggr-lab.github.io/gpcompreports/`, `https://gpcompreports.goatcounter.com/count` in `base.html`, and the `deploy_pages.sh` echo. The Python package directory stays `GPCompaReports_v2/`.
- **No em dashes in user-facing copy.** Use a colon, a comma, or a full stop.
- **Freeze rule:** with every V3 toggle off, the text content inside `<section class="report-section">` elements must be unchanged from the current build. Task 9 tests this.
- **No npm, no bundler, no framework.** New JS goes in `GPCompaReports_v2/static/js/` and is loaded with `defer` from a template block.
- **GPCRdb numbering lives in `display_number`.** `generic_number` and `gpcrdb_number` are empty in all 283 annotation files. Never read them without the `display_number` fallback.
- **Colour convention:** red = active-favoring (`--semantic-active`), blue = inactive-favoring (`--semantic-inactive`).
- **Install with** `pip3 install --break-system-packages` (no venv on this machine).
- **Branch:** cut `feature/v3-demo` from `feature/site-updates-2026-04-30` at `3a10a53`, not from `master` (master is 30 commits behind the deployed branch).
- **Demo receptors:** `par2_human`, `par1_human`, `adrb2_human`, `5ht2a_human`, `cxcr4_human`.
- All paths below are relative to `The Ultimate RRCS database project/alphafold_multistate_gpcr/`.

---

## File Structure

**New files:**
- `GPCompaReports_v2/website/brand.py` — brand constants, single source of the display name.
- `GPCompaReports_v2/analysis/receptor_profile.py` — segment profile, database median profile, key numbers, badge flags. Pure functions over DataFrames, no Jinja2, no I/O.
- `GPCompaReports_v2/static/js/v3-nav.js` — sticky section nav and compact view.
- `GPCompaReports_v2/static/js/v3-deeplink.js` — URL hash state for snake plot view and filters.
- `GPCompaReports_v2/static/js/v3-analysis.js` — analysis layer toggle, fingerprint SVG, table badges.
- `GPCompaReports_v2/static/css/v3.css` — styles for all of the above.
- `GPCompaReports_v2/tests/` — pytest suite.
- `scripts/build_demo.sh` — one command to build the demo.

**Modified files:**
- `GPCompaReports_v2/generate_site.py` — add `--only`.
- `GPCompaReports_v2/website/site_generator.py` — receptor selection, brand, median profile.
- `GPCompaReports_v2/website/page_generators/gpcr_report_page.py` — pass analysis-layer data.
- `GPCompaReports_v2/templates/` — `base.html`, `_partials/topbar.html`, `_partials/footer.html`, `landing.html`, `gpcr_report.html`.

---

### Task 1: Test harness and receptor selector

The demo needs to build five named receptors. `--limit N` takes the first N alphabetically, which does not include PAR2.

**Files:**
- Create: `GPCompaReports_v2/tests/__init__.py` (empty)
- Create: `GPCompaReports_v2/tests/test_selection.py`
- Modify: `GPCompaReports_v2/website/site_generator.py`
- Modify: `GPCompaReports_v2/generate_site.py:23-33`

**Interfaces:**
- Produces: `select_gpcr_ids(all_ids: list[str], limit: int|None = None, only: list[str]|None = None) -> list[str]` in `website/site_generator.py`. `only` wins over `limit`. Unknown ids raise `ValueError`.
- Produces: `SiteGenerator(..., only=None)` keyword argument.

- [ ] **Step 1: Write the failing test**

Create `GPCompaReports_v2/tests/test_selection.py`:

```python
import pytest

from GPCompaReports_v2.website.site_generator import select_gpcr_ids

ALL = ['5ht2a_human', 'adrb2_human', 'cxcr4_human', 'par1_human', 'par2_human']


def test_no_filters_returns_everything():
    assert select_gpcr_ids(ALL) == ALL


def test_limit_takes_a_prefix():
    assert select_gpcr_ids(ALL, limit=2) == ['5ht2a_human', 'adrb2_human']


def test_only_selects_named_receptors_in_given_order():
    assert select_gpcr_ids(ALL, only=['par2_human', 'adrb2_human']) == [
        'par2_human', 'adrb2_human']


def test_only_beats_limit():
    assert select_gpcr_ids(ALL, limit=1, only=['par2_human']) == ['par2_human']


def test_unknown_id_raises_with_the_bad_name_in_the_message():
    with pytest.raises(ValueError, match='nosuch_human'):
        select_gpcr_ids(ALL, only=['nosuch_human'])
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd "The Ultimate RRCS database project/alphafold_multistate_gpcr" && python3 -m pytest GPCompaReports_v2/tests/test_selection.py -v`
Expected: FAIL, `ImportError: cannot import name 'select_gpcr_ids'`

- [ ] **Step 3: Write minimal implementation**

Add to `GPCompaReports_v2/website/site_generator.py`, above the `SiteGenerator` class:

```python
def select_gpcr_ids(all_ids, limit=None, only=None):
    """Pick which receptors get report pages.

    `only` is an explicit ordered list of filesystem ids and wins over `limit`.
    `limit` takes the first N in the order the store discovered them.
    """
    if only:
        known = set(all_ids)
        missing = [g for g in only if g not in known]
        if missing:
            raise ValueError(
                "unknown gpcr id(s): %s" % ", ".join(sorted(missing)))
        return list(only)
    if limit:
        return list(all_ids[:limit])
    return list(all_ids)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python3 -m pytest GPCompaReports_v2/tests/test_selection.py -v`
Expected: 5 passed

- [ ] **Step 5: Wire it into the generator**

In `GPCompaReports_v2/website/site_generator.py`, change `__init__` to accept `only=None` and store `self.only = only`. In `run()`, replace the report-generation block (currently `if self.limit: ... generate_all_reports(...)`) with:

```python
        selected = select_gpcr_ids(self.store.gpcr_ids, self.limit, self.only)
        print(f"  Individual reports ({len(selected)} pages)...")
        generate_all_reports(env, self.store, self.output_dir,
                             analysis_results=self.analysis_results,
                             only=selected)
```

In `GPCompaReports_v2/website/page_generators/gpcr_report_page.py`, change the signature of `generate_all_reports` from `(env, store, output_dir, analysis_results=None, limit=None)` to `(env, store, output_dir, analysis_results=None, only=None)` and replace the loop's id source with `gpcr_ids = only if only is not None else store.gpcr_ids`. Delete the `limit` slicing.

In `GPCompaReports_v2/generate_site.py`, add after the `--limit` argument:

```python
    parser.add_argument('--only', type=str, default=None,
                        help='Comma-separated gpcr ids to build, e.g. par2_human,adrb2_human')
```

and pass it through:

```python
    generator = SiteGenerator(
        batch_dir=args.batch_dir,
        metadata_csv=args.metadata,
        output_dir=args.output,
        limit=args.limit,
        only=[s.strip() for s in args.only.split(',')] if args.only else None,
    )
```

- [ ] **Step 6: Verify end to end**

Run: `python3 GPCompaReports_v2/generate_site.py --output /tmp/v3sel --only par2_human,adrb2_human`
Expected: build completes, and `ls /tmp/v3sel/reports/` prints exactly `adrb2_human.html` and `par2_human.html`.

- [ ] **Step 7: Commit**

```bash
git add GPCompaReports_v2/tests GPCompaReports_v2/website/site_generator.py GPCompaReports_v2/generate_site.py GPCompaReports_v2/website/page_generators/gpcr_report_page.py
git commit -m "build: --only selects named receptors"
```

---

### Task 2: Brand sweep to GPCompaRe database

**Files:**
- Create: `GPCompaReports_v2/website/brand.py`
- Create: `GPCompaReports_v2/tests/test_brand.py`
- Modify: `GPCompaReports_v2/templates/base.html:28`, `_partials/topbar.html:12`, `_partials/footer.html:13`, `landing.html:8`
- Modify: the four page generators' `page_title` strings

**Interfaces:**
- Produces: `BRAND_NAME = 'GPCompaRe database'` and `BRAND_SHORT = 'GPCompaRe'` in `website/brand.py`.

- [ ] **Step 1: Write the failing test**

Create `GPCompaReports_v2/tests/test_brand.py`:

```python
from pathlib import Path

from GPCompaReports_v2.website.brand import BRAND_NAME, BRAND_SHORT

TEMPLATES = Path(__file__).resolve().parent.parent / 'templates'
BASE = Path(__file__).resolve().parent.parent.parent


def test_brand_constants():
    assert BRAND_NAME == 'GPCompaRe database'
    assert BRAND_SHORT == 'GPCompaRe'


def test_no_old_brand_text_left_in_templates():
    for path in TEMPLATES.rglob('*.html'):
        assert 'GPCompaReports' not in path.read_text(), path


def test_urls_and_hosts_are_untouched():
    base_html = (TEMPLATES / 'base.html').read_text()
    assert 'gpcompreports.goatcounter.com' in base_html
    deploy = (BASE / 'scripts' / 'deploy_pages.sh').read_text()
    assert 'tggr-lab.github.io/gpcompreports/' in deploy
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python3 -m pytest GPCompaReports_v2/tests/test_brand.py -v`
Expected: FAIL, `ModuleNotFoundError: No module named 'GPCompaReports_v2.website.brand'`

- [ ] **Step 3: Write minimal implementation**

Create `GPCompaReports_v2/website/brand.py`:

```python
"""Single source of the site's display name.

Display text only. Every URL and host string stays lowercase `gpcompreports`,
see docs/superpowers/specs/2026-08-24-gpcompare-v3-design.md section 5.
"""

BRAND_NAME = 'GPCompaRe database'
BRAND_SHORT = 'GPCompaRe'
```

- [ ] **Step 4: Replace the display strings**

In `templates/_partials/topbar.html:12`, change the wordmark span to `<span class="wordmark-text">GPCompaRe</span>` and the `aria-label` to `"GPCompaRe database home"`.

In `templates/_partials/footer.html:13`, change `<strong>GPCompaReports</strong>` to `<strong>GPCompaRe database</strong>`.

In `templates/base.html:28`, change the title fallback to:
`{{ page_title or 'GPCompaRe database: active-inactive contact changes across human Class A GPCRs' }}`

In `templates/landing.html:8`, change the hero title span to `<span class="hero-title-serif">GPCompaRe</span>` and in the About heading (line 99) change `About GPCompaReports` to `About GPCompaRe`.

In the four page generators, replace every `· GPCompaReports` suffix in `page_title` with `· GPCompaRe`, and in `landing_page.py` change the `page_title` to `'GPCompaRe database: active-inactive contact changes across human Class A GPCRs'`.

- [ ] **Step 5: Run test to verify it passes**

Run: `python3 -m pytest GPCompaReports_v2/tests/test_brand.py -v`
Expected: 3 passed

- [ ] **Step 6: Confirm no URL churn**

Run: `grep -rn "gpcompreports" GPCompaReports_v2/templates scripts/deploy_pages.sh | wc -l`
Expected: at least 2 (the GoatCounter host and the deploy echo), unchanged from before the task.

- [ ] **Step 7: Commit**

```bash
git add GPCompaReports_v2/website/brand.py GPCompaReports_v2/tests/test_brand.py GPCompaReports_v2/templates GPCompaReports_v2/website/page_generators
git commit -m "brand: GPCompaReports -> GPCompaRe database (display text only)"
```

---

### Task 3: Section anchors and sticky section nav

**Files:**
- Create: `GPCompaReports_v2/static/js/v3-nav.js`
- Create: `GPCompaReports_v2/static/css/v3.css`
- Modify: `GPCompaReports_v2/templates/gpcr_report.html` (add `id` and `data-section-title` to the nine sections, load the new assets)

**Interfaces:**
- Produces: each `<section class="report-section">` carries `id="sec-<slug>"` and `data-section-title="<short label>"`.
- Produces: `window.V3Nav = { build: function() {...} }`.

- [ ] **Step 1: Add anchors to the nine sections**

In `templates/gpcr_report.html`, add an id and a short title to each `<section class="report-section">`, in document order:

| line (approx) | id | data-section-title |
|------|----|--------------------|
| 64 | `sec-snake` | Snake plot |
| 122 | `sec-changes` | Largest changes |
| 167 | `sec-tm` | TM analysis |
| 231 | `sec-charts` | Visualizations |
| 274 | `sec-variants` | Variants |
| 325 | `sec-rrcs` | Complete RRCS |
| 372 | `sec-distribution` | Distribution |
| 390 | `sec-methods` | Methods |
| 423 | `sec-understanding` | Understanding |

Example for the first one:

```html
    <section class="report-section" id="sec-snake" data-section-title="Snake plot">
```

- [ ] **Step 2: Write the nav module**

Create `GPCompaReports_v2/static/js/v3-nav.js`:

```javascript
/* v3-nav.js — sticky section nav for report pages.
   Builds itself from [data-section-title] so adding a section needs no JS change.
   Navigation is wayfinding, not content: this is always on. */
(function () {
  'use strict';
  var V3Nav = (window.V3Nav = window.V3Nav || {});

  function build() {
    var sections = document.querySelectorAll('.report-section[data-section-title]');
    if (sections.length < 2) return;

    var nav = document.createElement('nav');
    nav.className = 'v3-secnav';
    nav.setAttribute('aria-label', 'Report sections');
    var list = document.createElement('div');
    list.className = 'v3-secnav-list';
    nav.appendChild(list);

    var links = [];
    Array.prototype.forEach.call(sections, function (sec) {
      var a = document.createElement('a');
      a.href = '#' + sec.id;
      a.textContent = sec.dataset.sectionTitle;
      a.dataset.target = sec.id;
      list.appendChild(a);
      links.push(a);
    });

    var host = document.querySelector('.report-page .section-inner');
    var anchor = document.querySelector('.report-stats');
    if (!host || !anchor) return;
    anchor.parentNode.insertBefore(nav, anchor.nextSibling);

    if (!('IntersectionObserver' in window)) return;
    var seen = {};
    var obs = new IntersectionObserver(function (entries) {
      entries.forEach(function (e) { seen[e.target.id] = e.isIntersecting; });
      var activeId = null;
      Array.prototype.forEach.call(sections, function (sec) {
        if (!activeId && seen[sec.id]) activeId = sec.id;
      });
      links.forEach(function (a) {
        a.classList.toggle('on', a.dataset.target === activeId);
      });
    }, { rootMargin: '-80px 0px -70% 0px' });
    Array.prototype.forEach.call(sections, function (s) { obs.observe(s); });
  }

  V3Nav.build = build;
  document.addEventListener('DOMContentLoaded', build);
})();
```

- [ ] **Step 3: Write the styles**

Create `GPCompaReports_v2/static/css/v3.css`:

```css
/* v3.css — V3 additions. Loaded only on pages that opt in via extra_css. */

.v3-secnav {
  position: sticky;
  top: var(--topbar-height, 68px);
  z-index: 20;
  background: var(--bgColor-default);
  border-bottom: 1px solid var(--borderColor-muted);
  margin: 0 0 var(--s-5);
}
.v3-secnav-list {
  display: flex;
  gap: var(--s-1);
  overflow-x: auto;
  padding: var(--s-2) 0;
  scrollbar-width: thin;
}
.v3-secnav a {
  flex: 0 0 auto;
  font-size: 13px;
  color: var(--fgColor-muted);
  text-decoration: none;
  padding: 4px 10px;
  border-radius: var(--r-pill);
  white-space: nowrap;
}
.v3-secnav a:hover { background: var(--bgColor-neutral-muted); color: var(--fgColor-default); }
.v3-secnav a.on { background: var(--bgColor-accent-muted); color: var(--fgColor-accent); font-weight: 600; }

@media print {
  .v3-secnav, .v3-toolbar, .topbar, .mobile-drawer { display: none !important; }
}
```

- [ ] **Step 4: Load the assets on report pages**

In `templates/gpcr_report.html`, inside `{% block head %}` after the Plotly script, add:

```html
<link rel="stylesheet" href="{{ static_prefix }}static/css/v3.css">
```

At the end of `{% block scripts %}`, before `{% endblock %}`, add:

```html
<script src="{{ static_prefix }}static/js/v3-nav.js" defer></script>
```

- [ ] **Step 5: Verify in a browser**

Run: `python3 GPCompaReports_v2/generate_site.py --output /tmp/v3demo --only par2_human`
Then open `/tmp/v3demo/reports/par2_human.html`.
Expected: a pill nav sits under the four stat tiles, sticks below the topbar on scroll, and the active pill tracks the section in view.

- [ ] **Step 6: Commit**

```bash
git add GPCompaReports_v2/static/js/v3-nav.js GPCompaReports_v2/static/css/v3.css GPCompaReports_v2/templates/gpcr_report.html
git commit -m "report: sticky section nav"
```

---

### Task 4: Compact view

Off by default. Folds each section to its heading. Remembered per reader.

**Files:**
- Modify: `GPCompaReports_v2/static/js/v3-nav.js` (append)
- Modify: `GPCompaReports_v2/static/css/v3.css` (append)
- Modify: `GPCompaReports_v2/templates/gpcr_report.html` (toolbar mount point)

**Interfaces:**
- Consumes: `.report-section[data-section-title]` from Task 3.
- Produces: `localStorage['gpcompare-compact']`, values `'1'` or absent.
- Produces: a `<div class="v3-toolbar">` inserted before the nav, which Task 7 also mounts into.

- [ ] **Step 1: Add the toolbar mount and compact logic**

Append to `GPCompaReports_v2/static/js/v3-nav.js`, inside the IIFE before `V3Nav.build = build;`:

```javascript
  var COMPACT_KEY = 'gpcompare-compact';

  function readFlag(key) {
    try { return localStorage.getItem(key) === '1'; } catch (e) { return false; }
  }
  function writeFlag(key, on) {
    try { on ? localStorage.setItem(key, '1') : localStorage.removeItem(key); } catch (e) {}
  }
  V3Nav.readFlag = readFlag;
  V3Nav.writeFlag = writeFlag;

  function ensureToolbar() {
    var existing = document.querySelector('.v3-toolbar');
    if (existing) return existing;
    var anchor = document.querySelector('.report-stats');
    if (!anchor) return null;
    var bar = document.createElement('div');
    bar.className = 'v3-toolbar';
    anchor.parentNode.insertBefore(bar, anchor.nextSibling);
    return bar;
  }
  V3Nav.ensureToolbar = ensureToolbar;

  function applyCompact(on) {
    document.querySelectorAll('.report-section[data-section-title]').forEach(function (sec) {
      sec.classList.toggle('v3-folded', on);
    });
  }

  function buildCompactToggle() {
    var bar = ensureToolbar();
    if (!bar) return;
    var btn = document.createElement('button');
    btn.type = 'button';
    btn.className = 'btn btn-ghost btn-sm v3-toggle';
    btn.id = 'v3-compact-toggle';
    var on = readFlag(COMPACT_KEY);
    function paint() {
      btn.setAttribute('aria-pressed', on ? 'true' : 'false');
      btn.classList.toggle('on', on);
      btn.textContent = on ? 'Compact view: on' : 'Compact view';
    }
    btn.addEventListener('click', function () {
      on = !on;
      writeFlag(COMPACT_KEY, on);
      applyCompact(on);
      paint();
    });
    bar.appendChild(btn);
    paint();
    applyCompact(on);
  }
```

Then inside `build()`, as its first statement, call `buildCompactToggle();`.

Note the insertion order: `ensureToolbar` inserts the toolbar directly after `.report-stats`, and `build()` inserts the nav after the same anchor afterwards, so the nav ends up above the toolbar. That is intentional: nav first, controls under it.

- [ ] **Step 2: Add the styles**

Append to `GPCompaReports_v2/static/css/v3.css`:

```css
.v3-toolbar { display: flex; gap: var(--s-2); flex-wrap: wrap; margin: 0 0 var(--s-5); }
.v3-toggle.on { background: var(--bgColor-accent-muted); color: var(--fgColor-accent); border-color: var(--fgColor-accent); }

.report-section.v3-folded > *:not(.report-section-head) { display: none; }
.report-section.v3-folded { padding-bottom: var(--s-2); }
.report-section.v3-folded .report-section-head { margin-bottom: 0; cursor: pointer; }
.report-section.v3-folded .report-section-head::after {
  content: '▸'; margin-left: var(--s-2); color: var(--fgColor-muted); font-size: 12px;
}
```

- [ ] **Step 3: Let a folded heading unfold itself**

Append inside the IIFE in `v3-nav.js`:

```javascript
  document.addEventListener('click', function (ev) {
    var head = ev.target.closest ? ev.target.closest('.report-section.v3-folded .report-section-head') : null;
    if (!head) return;
    head.parentNode.classList.remove('v3-folded');
  });
```

- [ ] **Step 4: Verify**

Rebuild with the Task 3 command and open the page.
Expected, in order: with a clean profile the page renders fully expanded and the button reads "Compact view". Click it, every section folds to its heading and the label reads "Compact view: on". Reload, it is still folded. Click a folded heading, that one section opens. Click the button again, everything expands and reloading keeps it expanded.

- [ ] **Step 5: Commit**

```bash
git add GPCompaReports_v2/static/js/v3-nav.js GPCompaReports_v2/static/css/v3.css
git commit -m "report: compact view toggle, off by default"
```

---

### Task 5: Deep links for snake plot state

**Files:**
- Create: `GPCompaReports_v2/static/js/v3-deeplink.js`
- Modify: `GPCompaReports_v2/templates/gpcr_report.html` (export the snake internals, load the module)

**Interfaces:**
- Produces: `window.SnakeAPI = { setView, getView, setMagRange, setDirection, getState }`, exported from the existing snake IIFE.
- Produces: hash grammar `#view=<name>&min=<float>&max=<float>&dir=<both|active|inactive>&sec=<section-id>`. Every key is optional.

- [ ] **Step 1: Export the snake internals**

In `templates/gpcr_report.html`, replace the three existing export lines near the end of the snake IIFE (`window.switchSnakeView = switchSnakeView;` through `window.toggleLoop = toggleLoop;`) with:

```javascript
  window.switchSnakeView = switchSnakeView;
  window.toggleSnakeLinks = toggleSnakeLinks;
  window.toggleLoop = toggleLoop;
  window.SnakeAPI = {
    setView: switchSnakeView,
    getView: function () {
      var b = document.querySelector('.snake-btn.active');
      return b ? b.dataset.view : 'delta_rrcs';
    },
    setMagRange: function (lo, hi) {
      var minSlider = document.getElementById('snake-mag-min');
      var maxSlider = document.getElementById('snake-mag-max');
      if (typeof lo === 'number' && !isNaN(lo)) {
        filterState.magMin = Math.max(minAbsDelta, Math.min(lo, maxAbsDelta));
        if (minSlider) minSlider.value = String(filterState.magMin);
      }
      if (typeof hi === 'number' && !isNaN(hi)) {
        filterState.magMax = Math.max(minAbsDelta, Math.min(hi, maxAbsDelta));
        if (maxSlider) maxSlider.value = String(filterState.magMax);
      }
      if (filterState.magMin > filterState.magMax) filterState.magMin = filterState.magMax;
      updateMagDisplays();
      redrawSnakeLinks();
    },
    setDirection: setDirection,
    setLinksVisible: toggleSnakeLinks,
    getState: function () {
      return {
        view: window.SnakeAPI.getView(),
        min: filterState.magMin,
        max: filterState.magMax,
        dir: filterState.direction,
        links: filterState.visible
      };
    },
    bounds: function () { return { lo: minAbsDelta, hi: maxAbsDelta }; }
  };
  document.dispatchEvent(new CustomEvent('snake-ready'));
```

- [ ] **Step 2: Write the deep-link module**

Create `GPCompaReports_v2/static/js/v3-deeplink.js`:

```javascript
/* v3-deeplink.js — put snake-plot state in the URL hash so a link can point at a view.
   Grammar: #view=cfr&min=0.5&max=3&dir=active&sec=sec-changes
   Reading the hash changes what is displayed but adds nothing to the page. */
(function () {
  'use strict';

  function parseHash() {
    var raw = (window.location.hash || '').replace(/^#/, '');
    var out = {};
    if (!raw) return out;
    raw.split('&').forEach(function (pair) {
      var kv = pair.split('=');
      if (kv.length === 2 && kv[0]) out[decodeURIComponent(kv[0])] = decodeURIComponent(kv[1]);
    });
    return out;
  }

  function writeHash(state) {
    var parts = [];
    if (state.view && state.view !== 'delta_rrcs') parts.push('view=' + state.view);
    var b = window.SnakeAPI.bounds();
    if (typeof state.min === 'number' && state.min > b.lo) parts.push('min=' + state.min.toFixed(2));
    if (typeof state.max === 'number' && state.max < b.hi) parts.push('max=' + state.max.toFixed(2));
    if (state.dir && state.dir !== 'both') parts.push('dir=' + state.dir);
    var hash = parts.length ? '#' + parts.join('&') : '';
    if (hash !== window.location.hash) {
      history.replaceState(null, '', window.location.pathname + window.location.search + hash);
    }
  }

  function applyHash() {
    var h = parseHash();
    if (!window.SnakeAPI) return;
    if (h.view) window.SnakeAPI.setView(h.view);
    if (h.min || h.max) {
      window.SnakeAPI.setLinksVisible(true);
      var cb = document.getElementById('snake-links-cb');
      if (cb) cb.checked = true;
      window.SnakeAPI.setMagRange(
        h.min ? parseFloat(h.min) : undefined,
        h.max ? parseFloat(h.max) : undefined);
    }
    if (h.dir) window.SnakeAPI.setDirection(h.dir);
    if (h.sec) {
      var target = document.getElementById(h.sec);
      if (target) target.scrollIntoView();
    }
  }

  function watch() {
    var sync = function () { writeHash(window.SnakeAPI.getState()); };
    document.querySelectorAll('.snake-btn').forEach(function (b) {
      b.addEventListener('click', function () { setTimeout(sync, 0); });
    });
    document.querySelectorAll('.snake-dir-btn').forEach(function (b) {
      b.addEventListener('click', function () { setTimeout(sync, 0); });
    });
    ['snake-mag-min', 'snake-mag-max'].forEach(function (id) {
      var el = document.getElementById(id);
      if (el) el.addEventListener('change', sync);
    });
  }

  function start() {
    applyHash();
    watch();
    window.addEventListener('hashchange', applyHash);
  }

  // Ordering: the snake plot's code is an INLINE script, so it runs at parse
  // time. This module is deferred, so it runs after parsing, by which point
  // SnakeAPI already exists and `snake-ready` has already fired. Check for the
  // object first and only fall back to the event, otherwise a listener attached
  // here would wait for an event that already happened and never apply a hash.
  if (window.SnakeAPI) {
    start();
  } else {
    document.addEventListener('snake-ready', start);
  }
})();
```

- [ ] **Step 3: Load it**

In `templates/gpcr_report.html`, in `{% block scripts %}` next to the nav script, add:

```html
<script src="{{ static_prefix }}static/js/v3-deeplink.js" defer></script>
```

The snake IIFE is inline and therefore runs at parse time, before any deferred script. So by the time `v3-deeplink.js` executes, `SnakeAPI` already exists and `snake-ready` has already fired. The module checks for the object first and treats the event only as a fallback. Listening for the event alone would wait forever.

Verify this specifically: a built page loaded at `#view=cfr` must come up on the Core Functional view. If it comes up on the default view, the ordering assumption is wrong and the whole feature is inert.

- [ ] **Step 4: Verify**

Rebuild, then open `reports/par2_human.html#view=cfr&dir=active`.
Expected: the page loads with the Core Functional view selected and the Active-favoring direction button active. Then click through a few snake views and watch the address bar update. Copy the URL, open it in a new tab, and confirm the same state comes back.

- [ ] **Step 5: Commit**

```bash
git add GPCompaReports_v2/static/js/v3-deeplink.js GPCompaReports_v2/templates/gpcr_report.html
git commit -m "report: deep-linkable snake plot state via URL hash"
```

---

### Task 6: Analysis layer data

Pure Python. No template changes. This task produces the numbers Task 7 renders.

**Files:**
- Create: `GPCompaReports_v2/analysis/receptor_profile.py`
- Create: `GPCompaReports_v2/tests/test_receptor_profile.py`

**Interfaces:**
- Produces: `SEGMENTS: list[str]` — the 16 segment names in N-to-C order.
- Produces: `segment_profile(delta_df, annot_map, threshold) -> dict[str, float]` — percent of above-threshold contact endpoints per segment, summing to 100 across the segments present.
- Produces: `median_profile(profiles: list[dict]) -> dict[str, float]`.
- Produces: `is_low_confidence(segment: str) -> bool` — True for terminus and loop segments.
- Produces: `key_numbers(delta_df, annot_map, threshold, cfr_ranks) -> dict` with keys `total_contacts`, `above_threshold`, `above_threshold_pct`, `threshold`, `top_segments` (list of `(segment, pct)`, length 3), `largest_structured` (dict or None), `cfr_top_movers` (int), `top_mover_count` (int).

- [ ] **Step 1: Write the failing test**

Create `GPCompaReports_v2/tests/test_receptor_profile.py`:

```python
import pandas as pd

from GPCompaReports_v2.analysis import receptor_profile as rp


def _delta(rows):
    df = pd.DataFrame(rows, columns=['res1', 'res2', 'delta_rrcs'])
    df['abs_delta'] = df['delta_rrcs'].abs()
    return df


ANNOT = {
    1: {'position': 1, 'amino_acid': 'A', 'protein_segment': 'TM6', 'display_number': '6.48x48'},
    2: {'position': 2, 'amino_acid': 'C', 'protein_segment': 'TM3', 'display_number': '3.50x50'},
    3: {'position': 3, 'amino_acid': 'D', 'protein_segment': 'N-term', 'display_number': ''},
    4: {'position': 4, 'amino_acid': 'E', 'protein_segment': 'ECL3', 'display_number': ''},
}


def test_segment_profile_counts_both_endpoints_of_each_contact():
    df = _delta([(1, 2, 5.0), (3, 4, 9.0)])
    prof = rp.segment_profile(df, ANNOT, threshold=1.0)
    assert prof['TM6'] == 25.0
    assert prof['TM3'] == 25.0
    assert prof['N-term'] == 25.0
    assert prof['ECL3'] == 25.0


def test_segment_profile_ignores_below_threshold_contacts():
    df = _delta([(1, 2, 5.0), (3, 4, 0.5)])
    prof = rp.segment_profile(df, ANNOT, threshold=1.0)
    assert prof['TM6'] == 50.0
    assert prof['N-term'] == 0.0


def test_segment_profile_of_empty_frame_is_all_zero():
    prof = rp.segment_profile(_delta([]), ANNOT, threshold=1.0)
    assert set(prof) == set(rp.SEGMENTS)
    assert sum(prof.values()) == 0.0


def test_median_profile_is_per_segment():
    a = {s: 0.0 for s in rp.SEGMENTS}
    b = {s: 0.0 for s in rp.SEGMENTS}
    c = {s: 0.0 for s in rp.SEGMENTS}
    a['TM6'], b['TM6'], c['TM6'] = 10.0, 20.0, 30.0
    assert rp.median_profile([a, b, c])['TM6'] == 20.0


def test_low_confidence_covers_termini_and_loops_but_not_helices():
    assert rp.is_low_confidence('N-term')
    assert rp.is_low_confidence('C-term')
    assert rp.is_low_confidence('ECL3')
    assert rp.is_low_confidence('ICL1')
    assert not rp.is_low_confidence('TM6')
    assert not rp.is_low_confidence('H8')


def test_largest_structured_skips_the_terminus_artifact():
    # The 9.0 contact is N-term to ECL3 and must not win.
    df = _delta([(1, 2, 5.0), (3, 4, 9.0)])
    kn = rp.key_numbers(df, ANNOT, threshold=1.0, cfr_ranks={})
    assert kn['largest_structured']['abs_delta'] == 5.0
    assert kn['largest_structured']['seg1'] == 'TM6'


def test_largest_structured_is_none_when_every_contact_is_low_confidence():
    df = _delta([(3, 4, 9.0)])
    kn = rp.key_numbers(df, ANNOT, threshold=1.0, cfr_ranks={})
    assert kn['largest_structured'] is None


def test_key_numbers_counts_cfr_top_movers_by_display_number():
    df = _delta([(1, 2, 5.0)])
    kn = rp.key_numbers(df, ANNOT, threshold=1.0, cfr_ranks={'3.50x50': 4})
    assert kn['cfr_top_movers'] == 1
    assert kn['above_threshold'] == 1
    assert kn['total_contacts'] == 1
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python3 -m pytest GPCompaReports_v2/tests/test_receptor_profile.py -v`
Expected: FAIL, `ModuleNotFoundError: No module named 'GPCompaReports_v2.analysis.receptor_profile'`

- [ ] **Step 3: Write the implementation**

Create `GPCompaReports_v2/analysis/receptor_profile.py`:

```python
"""Per-receptor derived numbers for the opt-in V3 analysis layer.

Pure functions over the frames the store already holds. Nothing here is
rendered unless the reader turns the analysis layer on, see
docs/superpowers/specs/2026-08-24-gpcompare-v3-design.md section 5.

GPCRdb numbering is read from `display_number`. The `generic_number` column
exists in every annotation CSV and is empty in all 283 of them.
"""

import statistics

SEGMENTS = [
    'N-term', 'TM1', 'ICL1', 'TM2', 'ECL1', 'TM3', 'ICL2', 'TM4',
    'ECL2', 'TM5', 'ICL3', 'TM6', 'ECL3', 'TM7', 'H8', 'C-term',
]

# AlphaFold-multistate is least reliable in termini and loops, and the single
# largest delta in a receptor usually lands there. Rows in these segments are
# marked low confidence and excluded from "largest structured change".
STRUCTURED = {'TM1', 'TM2', 'TM3', 'TM4', 'TM5', 'TM6', 'TM7', 'H8'}


def is_low_confidence(segment):
    """True for termini and loops, False for the seven helices and H8."""
    return segment not in STRUCTURED


def _seg(annot_map, position):
    entry = annot_map.get(int(position))
    if not entry:
        return 'unassigned'
    return entry.get('protein_segment') or 'unassigned'


def _num(annot_map, position):
    entry = annot_map.get(int(position))
    if not entry:
        return ''
    return entry.get('display_number') or entry.get('generic_number') or ''


def _aa(annot_map, position):
    entry = annot_map.get(int(position))
    if not entry:
        return ''
    return entry.get('amino_acid') or ''


def segment_profile(delta_df, annot_map, threshold):
    """Percent of above-threshold contact endpoints falling in each segment.

    Both ends of every above-threshold contact are counted, so a TM6-TM3
    contact contributes one to each. Returns a value for all 16 segments.
    """
    profile = {s: 0.0 for s in SEGMENTS}
    if delta_df is None or delta_df.empty:
        return profile
    sig = delta_df[delta_df['abs_delta'] >= threshold]
    if sig.empty:
        return profile
    counts = {}
    total = 0
    for _, row in sig.iterrows():
        for col in ('res1', 'res2'):
            seg = _seg(annot_map, row[col])
            counts[seg] = counts.get(seg, 0) + 1
            total += 1
    if not total:
        return profile
    for seg in SEGMENTS:
        profile[seg] = round(100.0 * counts.get(seg, 0) / total, 1)
    return profile


def median_profile(profiles):
    """Per-segment median across many receptor profiles."""
    out = {}
    for seg in SEGMENTS:
        values = [p.get(seg, 0.0) for p in profiles]
        out[seg] = round(statistics.median(values), 1) if values else 0.0
    return out


def key_numbers(delta_df, annot_map, threshold, cfr_ranks, top_n=6):
    """Numbers for the analysis-layer strip.

    `cfr_ranks` maps a display_number to its cross-receptor CFR rank.
    `largest_structured` is the biggest |delta| whose two endpoints are both
    in a helix or H8, or None if the receptor has no such contact.
    """
    total = 0 if delta_df is None else len(delta_df)
    if delta_df is None or delta_df.empty:
        return {
            'total_contacts': 0, 'above_threshold': 0, 'above_threshold_pct': 0.0,
            'threshold': round(float(threshold), 2), 'top_segments': [],
            'largest_structured': None, 'cfr_top_movers': 0, 'top_mover_count': 0,
        }

    sig = delta_df[delta_df['abs_delta'] >= threshold]
    profile = segment_profile(delta_df, annot_map, threshold)
    ranked = sorted(profile.items(), key=lambda kv: -kv[1])
    top_segments = [(s, p) for s, p in ranked if p > 0][:3]

    largest = None
    for _, row in sig.sort_values('abs_delta', ascending=False).iterrows():
        s1, s2 = _seg(annot_map, row['res1']), _seg(annot_map, row['res2'])
        if s1 in STRUCTURED and s2 in STRUCTURED:
            largest = {
                'abs_delta': round(float(row['abs_delta']), 2),
                'delta': round(float(row['delta_rrcs']), 2),
                'seg1': s1, 'seg2': s2,
                'num1': _num(annot_map, row['res1']),
                'num2': _num(annot_map, row['res2']),
                'label1': '%s%d' % (_aa(annot_map, row['res1']), int(row['res1'])),
                'label2': '%s%d' % (_aa(annot_map, row['res2']), int(row['res2'])),
            }
            break

    movers = sig.sort_values('abs_delta', ascending=False).head(top_n)
    mover_numbers = set()
    for _, row in movers.iterrows():
        for col in ('res1', 'res2'):
            num = _num(annot_map, row[col])
            if num:
                mover_numbers.add(num)
    cfr_hits = len([n for n in mover_numbers if n in cfr_ranks])

    return {
        'total_contacts': total,
        'above_threshold': int(len(sig)),
        'above_threshold_pct': round(100.0 * len(sig) / total, 1) if total else 0.0,
        'threshold': round(float(threshold), 2),
        'top_segments': top_segments,
        'largest_structured': largest,
        'cfr_top_movers': cfr_hits,
        'top_mover_count': int(len(movers)),
    }
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python3 -m pytest GPCompaReports_v2/tests/test_receptor_profile.py -v`
Expected: 8 passed

- [ ] **Step 5: Commit**

```bash
git add GPCompaReports_v2/analysis/receptor_profile.py GPCompaReports_v2/tests/test_receptor_profile.py
git commit -m "analysis: per-receptor segment profile and key numbers"
```

---

### Task 7: Analysis layer UI

**Files:**
- Modify: `GPCompaReports_v2/website/site_generator.py` (compute the median profile once)
- Modify: `GPCompaReports_v2/website/page_generators/gpcr_report_page.py` (per-receptor payload)
- Modify: `GPCompaReports_v2/templates/gpcr_report.html` (payload script tag)
- Create: `GPCompaReports_v2/static/js/v3-analysis.js`
- Modify: `GPCompaReports_v2/static/css/v3.css` (append)

**Interfaces:**
- Consumes: `receptor_profile.key_numbers`, `segment_profile`, `median_profile`, `is_low_confidence` from Task 6; `V3Nav.ensureToolbar`, `V3Nav.readFlag`, `V3Nav.writeFlag` from Task 4.
- Produces: `<script type="application/json" id="v3-analysis-data">` containing `{key_numbers, profile, median_profile, segments, cfr_ranks, low_confidence_segments}`.
- Produces: `localStorage['gpcompare-analysis']`.

- [ ] **Step 1: Compute the median profile once per build**

In `GPCompaReports_v2/website/site_generator.py`, import at the top:

```python
from ..analysis import receptor_profile as rprofile
from .page_generators import gpcr_report_helpers as helpers
```

`site_generator.py` already lives in `website/`, so the helpers import is a single-dot
relative import. `..analysis` reaches up to `GPCompaReports_v2.analysis`.

Then in `run()`, immediately after the CFR and variant analysis lines, add:

```python
        print("  Computing database median segment profile...")
        profiles = []
        for gid in self.store.gpcr_ids:
            delta_df = self.store.delta_data.get(gid)
            if delta_df is None or delta_df.empty:
                continue
            annot_map = self.store.get_annotation_map(gid)
            thr = helpers._calc_significance_threshold(delta_df)
            profiles.append(rprofile.segment_profile(delta_df, annot_map, thr))
        self.analysis_results['median_profile'] = rprofile.median_profile(profiles)
```

This loop runs over all 283 receptors even when only five get pages, because the median is the comparison that makes one receptor readable.

Use `store.get_annotation_map(gid)` rather than building the map from the annotation DataFrame directly. It routes every field through `_safe_str`, so a missing segment arrives as `''`. Building it with `row.to_dict()` instead leaves `NaN` floats, and `NaN or 'unassigned'` evaluates to `NaN` because NaN is truthy, which would poison the profile keys.

- [ ] **Step 2: Build the per-receptor payload**

In `GPCompaReports_v2/website/page_generators/gpcr_report_page.py`, add near the other imports:

```python
from ...analysis import receptor_profile as rprofile
```

Inside `generate_all_reports`, before the `template.render(` call, add:

```python
        cfr_ranks = {}
        cfr_table = (analysis_results or {}).get('cfr', {}).get('cfr_table')
        if cfr_table is not None and not cfr_table.empty:
            cfr_ranks = dict(zip(cfr_table['generic_number'], cfr_table['rank']))

        v3_payload = {
            'key_numbers': rprofile.key_numbers(delta_df, annot_map, sig_threshold, cfr_ranks),
            'profile': rprofile.segment_profile(delta_df, annot_map, sig_threshold),
            'median_profile': (analysis_results or {}).get('median_profile', {}),
            'segments': rprofile.SEGMENTS,
            'cfr_ranks': {k: int(v) for k, v in cfr_ranks.items()},
            'low_confidence_segments': [s for s in rprofile.SEGMENTS
                                        if rprofile.is_low_confidence(s)],
        }
```

and add `v3_analysis_json=json.dumps(v3_payload, separators=(',', ':')),` to the `template.render(...)` kwargs.

- [ ] **Step 3: Emit the payload and load the module**

In `templates/gpcr_report.html`, in `{% block scripts %}`, next to the other payload script tags add:

```html
<script type="application/json" id="v3-analysis-data">{{ v3_analysis_json | safe }}</script>
<script src="{{ static_prefix }}static/js/v3-analysis.js" defer></script>
```

- [ ] **Step 4: Write the analysis layer module**

Create `GPCompaReports_v2/static/js/v3-analysis.js`:

```javascript
/* v3-analysis.js — opt-in analysis layer.
   Off by default. When on it ADDS a key-number strip, a segment fingerprint,
   and row badges. It never removes or rewords anything already on the page. */
(function () {
  'use strict';
  var KEY = 'gpcompare-analysis';
  var data = null;

  function load() {
    var el = document.getElementById('v3-analysis-data');
    if (!el) return null;
    try { return JSON.parse(el.textContent); } catch (e) { return null; }
  }

  function tile(label, value, sub) {
    var d = document.createElement('div');
    d.className = 'v3-tile';
    var l = document.createElement('div'); l.className = 'v3-tile-lab'; l.textContent = label;
    var v = document.createElement('div'); v.className = 'v3-tile-val'; v.innerHTML = value;
    d.appendChild(l); d.appendChild(v);
    if (sub) { var s = document.createElement('div'); s.className = 'v3-tile-sub'; s.textContent = sub; d.appendChild(s); }
    return d;
  }

  function buildStrip(k) {
    var wrap = document.createElement('div');
    wrap.className = 'v3-strip';
    wrap.id = 'v3-strip';

    var segs = (k.top_segments || []).map(function (p) {
      return p[0] + ' ' + p[1].toFixed(1) + '%';
    }).join(' · ') || 'no above-threshold contacts';
    wrap.appendChild(tile('Change concentrates in', segs));

    if (k.largest_structured) {
      var L = k.largest_structured;
      var val = '<span class="v3-gn">' + (L.num1 || L.label1) + '</span> ↔ ' +
                '<span class="v3-gn">' + (L.num2 || L.label2) + '</span>';
      wrap.appendChild(tile('Largest change, structured region', val,
        L.label1 + ' (' + L.seg1 + ') / ' + L.label2 + ' (' + L.seg2 + ')  |Δ| ' + L.abs_delta));
    } else {
      wrap.appendChild(tile('Largest change, structured region', 'none',
        'every above-threshold contact touches a terminus or loop'));
    }

    wrap.appendChild(tile('Top movers that are CFRs',
      k.cfr_top_movers + ' of ' + k.top_mover_count,
      'cross-receptor Core Functional Residues'));
    return wrap;
  }

  function buildFingerprint(profile, median, segments) {
    var W = 640, H = 116, padB = 20, plot = H - padB;
    var max = 0;
    segments.forEach(function (s) {
      max = Math.max(max, profile[s] || 0, median[s] || 0);
    });
    max = Math.max(max, 1) * 1.15;
    var bw = W / segments.length;
    var NS = 'http://www.w3.org/2000/svg';
    var svg = document.createElementNS(NS, 'svg');
    svg.setAttribute('viewBox', '0 0 ' + W + ' ' + H);
    svg.setAttribute('role', 'img');
    svg.setAttribute('aria-label', 'Above-threshold contact distribution by segment, this receptor against the median of the database');

    segments.forEach(function (s, i) {
      var x = i * bw, v = profile[s] || 0, m = median[s] || 0;
      var h = v / max * plot, mh = m / max * plot;
      var structured = /^TM[1-7]$/.test(s) || s === 'H8';
      var bar = document.createElementNS(NS, 'rect');
      bar.setAttribute('x', x + bw * 0.18); bar.setAttribute('y', plot - h);
      bar.setAttribute('width', bw * 0.5); bar.setAttribute('height', Math.max(h, 0.8));
      bar.setAttribute('rx', '1.5');
      bar.setAttribute('fill', structured ? 'var(--brand-primary)' : 'var(--fgColor-muted)');
      bar.setAttribute('opacity', structured ? '1' : '0.45');
      var t = document.createElementNS(NS, 'title');
      t.textContent = s + ': ' + v.toFixed(1) + '% here, ' + m.toFixed(1) + '% median';
      bar.appendChild(t);
      svg.appendChild(bar);

      var line = document.createElementNS(NS, 'line');
      line.setAttribute('x1', x + bw * 0.10); line.setAttribute('x2', x + bw * 0.78);
      line.setAttribute('y1', plot - mh); line.setAttribute('y2', plot - mh);
      line.setAttribute('stroke', 'var(--fgColor-muted)');
      line.setAttribute('stroke-width', '1.4');
      line.setAttribute('stroke-dasharray', '3 2');
      svg.appendChild(line);

      var lab = document.createElementNS(NS, 'text');
      lab.setAttribute('x', x + bw * 0.43); lab.setAttribute('y', H - 6);
      lab.setAttribute('text-anchor', 'middle'); lab.setAttribute('class', 'v3-seg-lab');
      lab.textContent = s;
      svg.appendChild(lab);
    });

    var box = document.createElement('div');
    box.className = 'v3-fp';
    box.id = 'v3-fp';
    var head = document.createElement('div');
    head.className = 'v3-fp-head';
    head.innerHTML = '<span>Where above-threshold contacts sit, against the median of the database</span>' +
      '<span class="v3-legend"><i class="bar"></i>this receptor <i class="med"></i>median</span>';
    box.appendChild(head);
    box.appendChild(svg);
    return box;
  }

  function badge(text, cls) {
    var s = document.createElement('span');
    s.className = 'v3-badge ' + cls;
    s.textContent = text;
    return s;
  }

  function decorateChangesTable() {
    var table = document.getElementById('top-changes-table');
    if (!table) return;
    var lowSegs = data.low_confidence_segments || [];
    Array.prototype.forEach.call(table.tBodies[0].rows, function (tr) {
      if (tr.querySelector('.v3-badge')) return;
      var cells = tr.cells;
      if (cells.length < 3) return;
      var text = (cells[1].textContent + ' ' + cells[2].textContent);
      var holder = document.createElement('span');
      holder.className = 'v3-badges';

      var nums = text.match(/\d+\.\d+x\d+/g) || [];
      var best = null;
      nums.forEach(function (n) {
        var r = data.cfr_ranks[n];
        if (r && (best === null || r < best)) best = r;
      });
      if (best !== null) {
        holder.appendChild(badge(best <= 30 ? 'CFR #' + best : 'CFR top 50', 'v3-b-cfr'));
      }
      // A residue with no GPCRdb generic number is outside the helices, which
      // means a terminus or a loop: exactly where the model is least reliable.
      // lowSegs is carried in the payload for the tooltip text below.
      if (!nums.length) {
        var b = badge('low confidence', 'v3-b-low');
        b.title = 'Outside the helices (' + lowSegs.join(', ') + '). ' +
                  'Predicted conformations are least reliable there.';
        holder.appendChild(b);
      }
      if (holder.children.length) cells[cells.length - 1].appendChild(holder);
    });
  }

  function apply(on) {
    document.documentElement.classList.toggle('v3-analysis-on', on);
    if (on) decorateChangesTable();
  }

  function build() {
    data = load();
    if (!data || !window.V3Nav) return;
    var bar = window.V3Nav.ensureToolbar();
    if (!bar) return;

    var strip = buildStrip(data.key_numbers);
    var fp = buildFingerprint(data.profile, data.median_profile, data.segments);
    var host = document.createElement('div');
    host.className = 'v3-layer';
    host.appendChild(strip);
    host.appendChild(fp);
    // Mount after the toolbar, NOT after .report-stats. Tasks 3 and 4 both
    // insert directly after .report-stats, so anything inserted there later
    // jumps ahead of them. Final order must be: stats, nav, toolbar, layer.
    bar.parentNode.insertBefore(host, bar.nextSibling);

    var btn = document.createElement('button');
    btn.type = 'button';
    btn.className = 'btn btn-ghost btn-sm v3-toggle';
    btn.id = 'v3-analysis-toggle';
    var on = window.V3Nav.readFlag(KEY);
    function paint() {
      btn.setAttribute('aria-pressed', on ? 'true' : 'false');
      btn.classList.toggle('on', on);
      btn.textContent = on ? 'Analysis layer: on' : 'Analysis layer';
    }
    btn.addEventListener('click', function () {
      on = !on;
      window.V3Nav.writeFlag(KEY, on);
      apply(on);
      paint();
    });
    bar.appendChild(btn);
    paint();
    apply(on);
  }

  document.addEventListener('DOMContentLoaded', function () { setTimeout(build, 0); });
})();
```

- [ ] **Step 5: Add the styles**

Append to `GPCompaReports_v2/static/css/v3.css`:

```css
.v3-layer { display: none; margin: 0 0 var(--s-5); }
.v3-analysis-on .v3-layer { display: block; }
.v3-badges { display: none; }
.v3-analysis-on .v3-badges { display: inline-flex; gap: 4px; margin-left: 6px; }

.v3-strip { display: grid; grid-template-columns: repeat(3, 1fr); gap: var(--s-2); margin-bottom: var(--s-3); }
.v3-tile { border: 1px solid var(--borderColor-muted); border-radius: var(--r-md); padding: 10px 12px; }
.v3-tile-lab { font-size: 10.5px; text-transform: uppercase; letter-spacing: .06em;
               color: var(--fgColor-muted); font-weight: 600; margin-bottom: 4px; }
.v3-tile-val { font-size: 15px; font-weight: 700; }
.v3-tile-sub { font-size: 12px; color: var(--fgColor-muted); margin-top: 3px; }
.v3-gn { font-family: var(--font-mono); font-size: 12.5px;
         background: var(--bgColor-neutral-muted); border-radius: 4px; padding: 1px 5px; }

.v3-fp { border: 1px solid var(--borderColor-muted); border-radius: var(--r-md); padding: 12px 12px 4px; }
.v3-fp-head { display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 8px;
              font-size: 10.5px; text-transform: uppercase; letter-spacing: .06em;
              color: var(--fgColor-muted); font-weight: 600; margin-bottom: 8px; }
.v3-legend { text-transform: none; letter-spacing: 0; font-weight: 400; }
.v3-legend i { display: inline-block; width: 9px; height: 9px; border-radius: 2px; margin: 0 4px 0 10px; vertical-align: -1px; }
.v3-legend i.bar { background: var(--brand-primary); }
.v3-legend i.med { background: var(--fgColor-muted); }
.v3-seg-lab { font-size: 8.5px; fill: var(--fgColor-muted); }

.v3-badge { font-size: 10px; font-weight: 700; border-radius: var(--r-pill); padding: 1px 7px; white-space: nowrap; }
.v3-b-cfr { background: var(--bgColor-accent-muted); color: var(--fgColor-accent); }
.v3-b-low { background: var(--bgColor-neutral-muted); color: var(--fgColor-muted); }

@media (max-width: 820px) { .v3-strip { grid-template-columns: 1fr; } }
```

- [ ] **Step 6: Verify**

Run: `python3 GPCompaReports_v2/generate_site.py --output /tmp/v3demo --only par2_human,adrb2_human`
Open `/tmp/v3demo/reports/par2_human.html`.
Expected with a clean profile: no strip, no fingerprint, no badges, and a button reading "Analysis layer". Click it: three tiles appear, the fingerprint draws with PAR2's TM6 bar near 17.5% against a dashed median near 11.7%, and CFR badges appear in the largest-changes table. The R36/Q319 row gets a "low confidence" badge and is not what the "largest structured" tile names. Reload and the layer is still on.

Open `adrb2_human.html` with the layer on: it must render without error despite 35% of its endpoints being `unassigned`.

- [ ] **Step 7: Commit**

```bash
git add GPCompaReports_v2/static/js/v3-analysis.js GPCompaReports_v2/static/css/v3.css GPCompaReports_v2/website GPCompaReports_v2/templates/gpcr_report.html
git commit -m "report: opt-in analysis layer with segment fingerprint and CFR badges"
```

---

### Task 8: Front door

**Files:**
- Modify: `GPCompaReports_v2/templates/landing.html` (new explainer section above the stat band)
- Modify: `GPCompaReports_v2/static/css/landing.css` (append)

**Interfaces:**
- Consumes: `total_gpcrs` from `landing_page.py`, already passed.

- [ ] **Step 1: Add the explainer**

In `templates/landing.html`, insert between the `</section>` that closes `.hero` and the `<section class="stat-band">`:

```html
<section class="explainer" aria-labelledby="explainer-heading">
  <div class="section-inner">
    <header class="section-head">
      <h2 id="explainer-heading">What this database measures</h2>
      <p class="section-sub">Every receptor here is compared in two predicted conformations. The number we report is how much each pair of residues touches in one state versus the other.</p>
    </header>

    <div class="explainer-grid">
      <figure class="explainer-fig">
        <svg viewBox="0 0 320 150" role="img" aria-labelledby="explainer-svg-title">
          <title id="explainer-svg-title">Two residues in contact, drawn in the inactive and the active conformation</title>
          <g>
            <text x="78" y="18" text-anchor="middle" class="ex-cap">Inactive</text>
            <circle cx="52" cy="66" r="16" class="ex-res"/>
            <circle cx="104" cy="66" r="16" class="ex-res"/>
            <line x1="68" y1="66" x2="88" y2="66" class="ex-bond ex-bond-weak"/>
            <text x="52" y="71" text-anchor="middle" class="ex-lab">A</text>
            <text x="104" y="71" text-anchor="middle" class="ex-lab">B</text>
            <text x="78" y="108" text-anchor="middle" class="ex-num">RRCS 0.4</text>
          </g>
          <g>
            <text x="242" y="18" text-anchor="middle" class="ex-cap">Active</text>
            <circle cx="222" cy="66" r="16" class="ex-res"/>
            <circle cx="262" cy="66" r="16" class="ex-res"/>
            <line x1="238" y1="66" x2="246" y2="66" class="ex-bond ex-bond-strong"/>
            <text x="222" y="71" text-anchor="middle" class="ex-lab">A</text>
            <text x="262" y="71" text-anchor="middle" class="ex-lab">B</text>
            <text x="242" y="108" text-anchor="middle" class="ex-num">RRCS 2.9</text>
          </g>
          <text x="160" y="138" text-anchor="middle" class="ex-delta">ΔRRCS +2.5, this contact is stronger in the active state</text>
        </svg>
      </figure>

      <div class="explainer-copy">
        <p><strong>RRCS</strong> scores how strongly two residues are in contact in one structure. Comparing the same pair across an active and an inactive model gives <strong>ΔRRCS</strong>, the change in that contact between states.</p>
        <p class="explainer-key">
          <span class="ex-swatch is-active"></span> <strong>Red</strong> marks contacts that are stronger in the active state.
          <span class="ex-swatch is-inactive"></span> <strong>Blue</strong> marks contacts that are stronger in the inactive state.
        </p>
        <p>Each receptor page reports every contact pair we detected, which of them change more than that receptor's own threshold, and where those changes sit in the structure.</p>
        <p><a href="{{ static_prefix }}statistics.html">See the patterns across all {{ total_gpcrs }} receptors</a></p>
      </div>
    </div>
  </div>
</section>
```

- [ ] **Step 2: Style it**

Append to `GPCompaReports_v2/static/css/landing.css`:

```css
/* ----------- Explainer ----------- */
.explainer { padding: var(--s-7) 0 var(--s-6); border-bottom: 1px solid var(--borderColor-muted); }
.explainer-grid { display: grid; grid-template-columns: 340px 1fr; gap: var(--s-6); align-items: center; }
.explainer-fig { margin: 0; }
.explainer-fig svg { width: 100%; height: auto; }
.ex-res { fill: var(--bgColor-neutral-muted); stroke: var(--borderColor-default); stroke-width: 1.5; }
.ex-lab { font-size: 13px; font-weight: 700; fill: var(--fgColor-default); font-family: var(--font-ui); }
.ex-cap { font-size: 11px; font-weight: 600; fill: var(--fgColor-muted);
          text-transform: uppercase; letter-spacing: .07em; font-family: var(--font-ui); }
.ex-num { font-size: 12px; fill: var(--fgColor-muted); font-family: var(--font-mono); }
.ex-delta { font-size: 11.5px; fill: var(--fgColor-default); font-family: var(--font-ui); }
.ex-bond { stroke-linecap: round; }
.ex-bond-weak { stroke: var(--fgColor-muted); stroke-width: 2; stroke-dasharray: 3 3; }
.ex-bond-strong { stroke: var(--semantic-active); stroke-width: 6; }
.explainer-copy p { margin: 0 0 var(--s-3); max-width: 60ch; }
.explainer-key { display: flex; align-items: center; flex-wrap: wrap; gap: 6px; }
.ex-swatch { display: inline-block; width: 12px; height: 12px; border-radius: 3px; margin-left: 10px; }
.ex-swatch.is-active { background: var(--semantic-active); margin-left: 0; }
.ex-swatch.is-inactive { background: var(--semantic-inactive); }
@media (max-width: 820px) { .explainer-grid { grid-template-columns: 1fr; } }
```

- [ ] **Step 3: Verify**

Rebuild and open `/tmp/v3demo/index.html`.
Expected: below the hero, a two-column block with the contact diagram on the left and the explanation on the right. The strong bond is red. In dark mode the diagram stays legible because every colour is a token.

- [ ] **Step 4: Commit**

```bash
git add GPCompaReports_v2/templates/landing.html GPCompaReports_v2/static/css/landing.css
git commit -m "landing: teach RRCS before counting records"
```

---

### Task 9: Demo build and freeze verification

The freeze rule is the one thing in this plan that a reviewer cannot check by eye. This task proves it.

**Files:**
- Create: `scripts/build_demo.sh`
- Create: `GPCompaReports_v2/tests/test_freeze.py`

**Interfaces:**
- Consumes: everything above.
- Produces: `GPCompaReports_v2/output_v3_demo/` containing the five-receptor demo.

- [ ] **Step 1: Write the demo build script**

Create `scripts/build_demo.sh`:

```bash
#!/usr/bin/env bash
# Build the five-receptor V3 demo.
#
# Conservation and AlphaMissense caches live under the OUTPUT dir, so a fresh
# output directory silently loses the full snake-plot colour views. Seed them
# from the main build before generating.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$ROOT/GPCompaReports_v2/output_v3_demo"
SRC_DATA="$ROOT/GPCompaReports_v2/output/data"
RECEPTORS="par2_human,par1_human,adrb2_human,5ht2a_human,cxcr4_human"

mkdir -p "$OUT/data"
if [ -d "$SRC_DATA" ]; then
  echo "==> Seeding conservation / AlphaMissense caches..."
  for gid in ${RECEPTORS//,/ }; do
    cp -f "$SRC_DATA/conservation_${gid}.json"  "$OUT/data/" 2>/dev/null || true
    cp -f "$SRC_DATA/alphamissense_${gid}.json" "$OUT/data/" 2>/dev/null || true
  done
fi

echo "==> Building demo..."
python3 "$ROOT/GPCompaReports_v2/generate_site.py" --output "$OUT" --only "$RECEPTORS"

echo ""
echo "==> Done. Open:"
echo "    $OUT/index.html"
```

Then: `chmod +x scripts/build_demo.sh`

- [ ] **Step 2: Write the freeze test**

Create `GPCompaReports_v2/tests/test_freeze.py`:

```python
"""The PI approved the v2 report pages. With every V3 toggle off, the text
inside the report sections must be what it was. This test compares a freshly
built demo report against the last full build in GPCompaReports_v2/output/.
"""
import re
from pathlib import Path

import pytest

BASE = Path(__file__).resolve().parent.parent
OLD = BASE / 'output' / 'reports' / 'par2_human.html'
NEW = BASE / 'output_v3_demo' / 'reports' / 'par2_human.html'

SECTION_RE = re.compile(
    r'<section class="report-section"[^>]*>(.*?)</section>', re.S)
TAG_RE = re.compile(r'<[^>]+>')
WS_RE = re.compile(r'\s+')


def _section_text(html):
    out = []
    for body in SECTION_RE.findall(html):
        text = WS_RE.sub(' ', TAG_RE.sub(' ', body)).strip()
        out.append(text)
    return out


@pytest.mark.skipif(not OLD.exists() or not NEW.exists(),
                    reason='needs both a full build and a demo build on disk')
def test_report_sections_are_unchanged_with_features_off():
    old = _section_text(OLD.read_text(encoding='utf-8'))
    new = _section_text(NEW.read_text(encoding='utf-8'))
    assert old, 'baseline produced no sections, the regex or the fixture is wrong'
    assert len(new) == len(old), 'section count changed'
    for i, (a, b) in enumerate(zip(old, new)):
        assert a == b, 'section %d text changed' % i


V3_TOKEN_RE = re.compile(r'v3-[a-z-]+')


@pytest.mark.skipif(not NEW.exists(), reason='needs a demo build on disk')
def test_no_v3_markup_inside_approved_sections():
    """Structural half of the freeze rule.

    The text comparison above catches added words. This catches added markup
    that carries no text yet, which is how a V3 feature would most likely
    leak: an empty mount point or a hidden span inside approved content.
    """
    html = NEW.read_text(encoding='utf-8')
    for i, body in enumerate(SECTION_RE.findall(html)):
        found = V3_TOKEN_RE.findall(body)
        assert not found, 'section %d contains V3 markup: %s' % (i, sorted(set(found)))
```

- [ ] **Step 3: Build the demo**

Run: `bash scripts/build_demo.sh`
Expected: the build reports 5 report pages and prints the path to `output_v3_demo/index.html`.

- [ ] **Step 4: Run the freeze test**

Run: `python3 -m pytest GPCompaReports_v2/tests/test_freeze.py -v`
Expected: PASS.

If it fails, the diff it reports is a section whose text changed. Fix by moving whatever was added out of the `<section class="report-section">` elements, not by relaxing the test. The analysis layer mounts before `.report-stats`, and its badges live inside `.v3-badges` spans, so a failure here means something leaked into approved content.

- [ ] **Step 5: Run the whole suite**

Run: `python3 -m pytest GPCompaReports_v2/tests -v`
Expected: all tests pass.

- [ ] **Step 6: Commit**

```bash
git add scripts/build_demo.sh GPCompaReports_v2/tests/test_freeze.py
git commit -m "demo: five-receptor V3 build plus freeze verification"
```

---

## Out of scope for the demo

Deliberately excluded, from spec sections 5 and 8. Do not build these:

- Provenance pages (Methods, Downloads, Cite this, changelog, version stamp).
- Machine readability (Open Graph, `schema.org` Dataset JSON-LD, sitemap).
- Generator cleanup: dead v1 `generate_all_reports`, duplicated TM segment lists, the wrong CFR-score tooltip in `statistics.html:104`, the reversed colour docstring in `landing_page.py:19`.
- The URL move. The site stays at `tggr-lab.github.io/gpcompreports/`.
- Surfacing the snake-plot filter panel. Task 5 opens it when a link carries a magnitude
  range, but the panel stays behind its checkbox otherwise.
- Lazy chart rendering, the full accessibility pass, and a dedicated mobile pass. The demo
  ships the responsive breakpoints written into Tasks 7 and 8 and the print rules in Task 3,
  nothing more.
- Deploying anything. The demo is local, and deploys are always manual.
