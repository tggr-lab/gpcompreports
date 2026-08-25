# GPCompaRe V3 design

Date: 2026-08-24
Status: design agreed, not started
Site: `GPCompaReports_v2/` on branch `feature/site-updates-2026-04-30` (head `3a10a53`)

## 1. What V3 is for

The site serves one reader: someone who arrives from the GPCompaRe manuscript. That reader
needs to grasp what RRCS is, trust the numbers, and cite the resource. V3 targets those three
jobs and nothing else.

V3 does not reproduce the paper's figures or tables. The site stays a browsable database of
283 human Class A GPCRs, and the paper cites it.

## 2. The governing constraint

The PI approved the v2 report pages. Yamir does not want what a report shows or says to
change. Every V3 decision follows from one rule:

> The default view of an approved page is exactly what the PI signed off.
> Anything new is opt-in, off by default, and remembered per reader.

Two consequences worth stating, because both killed earlier proposals:

- **No generated prose on reports.** An auto-written per-receptor summary would put 283
  sentences in front of a referee, each one a place to object.
- **Nothing collapses by default.** Hiding approved content is a subtraction, so the Compact
  view described in section 5 has to start off.

## 3. Weak spots V3 fixes

Yamir picked three out of four candidates. Page weight was explicitly not one of them.

1. **No orientation.** The landing page sells a database (283 receptors, N families, K contacts)
   without teaching what RRCS is, what a positive delta means, or why red and blue.
2. **The report page is a wall.** Nine stacked sections with no way to move between them.
   Note the tables are already paginated at 25 rows, so row count is not the problem.
3. **No provenance or citation.** Nothing states the data version, the build date, the pipeline,
   or how to cite. A reviewer cannot pin the numbers to anything permanent.

## 4. Engineering envelope

New vanilla-JS modules alongside `static/js/site.js`. No bundler, no npm, no framework. One
Python command still builds the site. A generator cleanup pass rides along, invisible to readers.

## 5. Scope

### V3.0 Foundation

- **Brand sweep.** Display text becomes `GPCompaRe database` in templates, page titles, footer,
  docstrings. **URLs and hosts do not move.** Yamir chose on 2026-08-24 to keep
  `tggr-lab.github.io/gpcompreports/` and may revisit later, so the brand-versus-URL split set in
  May extends to the new name. Leave the GitHub remote, the GoatCounter host, and every lowercase
  `gpcompreports` string alone.
- **Generator cleanup.** Delete the dead v1 copy of `generate_all_reports`
  (`gpcr_report_helpers.py:134-263`). Share the TM segment list instead of defining it twice
  (`gpcr_report_helpers.py:45-49` and `tm_domain_analysis.py:10-11`). Fix the CFR-score tooltip in
  `statistics.html:104`, which describes `frequency x mean |delta|` while `cfr_analysis.py:94`
  computes `(norm_frequency + norm_delta) / 2`. Fix the reversed colour docstring at
  `landing_page.py:19`. Document or retire the always-empty `generic_number` column. Add tests
  covering the threshold formula, `_compute_residue_metrics`, and the `display_number` fallback,
  since those three feed every residue-level view on the site.
- **Provenance.** A data version stamp on every page. A Methods page covering the pipeline from
  AlphaFold-multistate models through RRCS to the per-receptor threshold. A Downloads page for the
  dataset. A "Cite this" block. A changelog.
- **Machine readability.** Open Graph tags, `schema.org` Dataset JSON-LD, `sitemap.xml`,
  canonical links.

### V3.1 Front door

- Rebuild the landing page so it teaches before it counts. A figure-first explainer showing one
  receptor in two states, a handful of contacts, and the red/blue convention rendered rather than
  described. Inline SVG, in the same visual language as the snake plots and sparklines. Counters
  move below the explainer.
- Mobile work lands here: the snake plot and the wide tables are the two things that break on a
  phone, and both need a horizontal-scroll container rather than a shrunk layout.

Sequenced ahead of navigation at Yamir's request: first impression outranks depth for a reader
arriving from a paper. It depends on the Methods page from V3.0 so the landing page can link out
instead of explaining everything itself.

### V3.2 Navigation and deep links

- **Sticky section nav** on report pages, always on. Wayfinding is not content, so this needs no
  approval.
- **Compact view**, off by default, remembered in `localStorage`. Folds the nine sections to
  headers for a returning reader. The PI always opens the approved page.
- **Deep links.** Snake plot view, magnitude range, direction filter, and section anchor encoded in
  the URL hash, so `reports/par2_human.html#view=cfr&min=0.5` restores that state. Nothing new is
  displayed. This is the feature that most earns "paper companion": the manuscript can point at a
  view without the site hosting a figure.
- **Surface the filter panel.** The magnitude slider and direction toggle already exist behind a
  control most readers never press.
- Print stylesheet and accessibility pass ride along.

### V3.3 Analysis layer

One toggle, off by default, remembered per reader. When on it reveals:

- A key-number strip: where change concentrates by segment, the largest change within a structured
  region, how many top movers are CFRs, how many contact-site variants AlphaMissense calls pathogenic.
- A segment fingerprint: this receptor's above-threshold contact distribution against the median of
  283. The comparison carries the information. A bare per-receptor chart shows TM6 highest almost
  everywhere and says nothing.
- Row badges on existing tables: CFR membership, low confidence on terminus and loop rows,
  pathogenic variant present.

Nothing is removed when the layer is on. Nothing is added when it is off.

Last in the sequence because it is the only piece that puts new material on an approved page, so it
is the only piece needing the PI. Everything else ships while that conversation happens, and a no
costs nothing else in V3.

Mockup with real PAR2 numbers, three header layouts and the badge treatment:
`2026-08-24-v3-report-layouts.html`, beside this file. Open it in a browser. The three layouts there
were drawn before the freeze was known, so treat them as the content of the layer rather than as
proposals for the default header.

## 6. Data notes that constrain the build

Verified on disk on 2026-08-24. Anything in V3.3 depends on all three.

- **GPCRdb numbering lives in `display_number`.** The `generic_number` and `gpcrdb_number` columns
  are present and empty in all 283 annotation files. Every existing module falls back with
  `display_number or generic_number`. New code must read `display_number`.
- **The largest change in a receptor is usually an artifact.** PAR2's biggest is R36 (N-term)
  against Q319 (ECL3) at 12.88. ADRB2's is H390 against Q391, both C-terminal. These sit in
  disordered regions where AlphaFold-multistate is least reliable. The key-number strip must rank
  within structured segments, and the low-confidence badge exists to mark these rows.
- **Segment coverage is uneven.** ADRB2 puts 35% of its above-threshold contact endpoints in a
  segment labelled `unassigned`. Every derived view degrades gracefully or it breaks on receptors
  like that one.

Worked example, PAR2 against the database median, showing why the comparison matters:

| segment | PAR2 | median of 283 |
|---------|------|---------------|
| TM2     | 14.9% | 5.4% |
| TM6     | 17.5% | 11.7% |
| TM7     | 14.9% | 10.2% |
| TM5     | 4.5%  | 8.5% |
| ECL2    | 2.6%  | 6.0% |

## 7. Decisions taken, with reasons

| Decision | Reason |
|----------|--------|
| Paper companion, not a community database | Chosen audience. Optimise first-visit clarity and citability over depth features. |
| Harden the front door rather than mirror the paper | The site stays a database. The paper cites it. |
| No prose on reports | 283 generated sentences are 283 chances for a referee to object. |
| Rank computed but not printed | ADRB2 ranks 252 of 283. A reviewer who sees their receptor called 252nd attacks the metric instead of reading on. |
| No family-relative claims | The PAR family has n=4. "Above the family median" means little at that size. |
| Keep the legacy URL | Yamir's call, 2026-08-24. Revisit before the manuscript commits an address to print. |
| No new toolchain | The site is mostly finished and needs better explanation, not better ergonomics. A second build system would cost every future editor an `npm install`. |

## 8. Open questions

1. Does the opt-in analysis layer need the PI's approval before shipping, given it is off by default?
2. If the URL is revisited, the move has to happen before the manuscript prints an address.
   After that the old path is permanent.
3. What goes on the Downloads page: the batch CSVs as-is, or a curated release archive?
