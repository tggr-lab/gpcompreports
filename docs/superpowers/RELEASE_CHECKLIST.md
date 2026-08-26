# Release checklist

Required checks before publishing a build, whether to gh-pages or to any
other public location. This list is short on purpose. Everything on it is
enforced by `scripts/deploy_pages.sh` except where noted.

---

## The rule about skipped tests

**A skipped test is not a passing test.**

Most of this project's suite is guarded: tests skip when the build they read
is not on disk, and the browser tests skip when Playwright is not installed
in the interpreter running them. That is correct for day-to-day work, where a
developer without a demo build should not be blocked by 11 red tests.

It is not correct for a release. At release time a skip means the property was
never checked, which is indistinguishable from the property being broken. So
for the checks below, a skip is a failure.

Concretely: `python3 -m pytest GPCompaReports_v2/tests` reporting
`23 passed, 2 skipped` is **not** sufficient to release. One of those skips is
the browser suite, and it has to actually run.

---

## 1. Full suite, with a build present

```bash
bash scripts/build_demo.sh          # or a full build
python3 -m pytest GPCompaReports_v2/tests -rs
```

Read the `-rs` skip report. Every skip must be one you can explain. The two
expected at demo scale are the browser suite (covered by check 2, which runs
it properly) and `test_full_build_covers_every_manifest_receptor`, which
cannot assert full coverage against a five-receptor demo.

Before a real release, run the freeze test against the full build so that
second skip disappears:

```bash
FREEZE_BUILD_DIR=output python3 -m pytest GPCompaReports_v2/tests/test_freeze.py -rs
```

## 2. Browser smoke tests (required, must not be skipped)

```bash
bash scripts/run_browser_tests.sh
```

This covers `v3-nav.js` and `v3-deeplink.js`, the only JavaScript the V3
branch adds to a report page. Neither can be checked by the Python suite: the
Compact-view defect that reached a built page was a CSS specificity problem,
invisible to any test that does not render the page.

`run_browser_tests.sh` **cannot** silently pass. It exits non-zero when the
venv is absent, when Playwright is not installed in it, when there is no build
to test, and when pytest collects no tests. Treat any non-zero exit as a
release blocker, not as "not applicable".

To test exactly what is about to be published rather than the demo:

```bash
SMOKE_BUILD_DIR=output bash scripts/run_browser_tests.sh
```

`scripts/deploy_pages.sh` runs this automatically against the freshly built
site, after the build and before anything is pushed.

## 3. Freeze baseline intact

The PI approved the v2 report content. `tests/fixtures/freeze_manifest.json`
pins it. The freeze test must pass against the build being released, and the
manifest must not have been regenerated as part of the release. Regenerating
it is a separate, deliberate act with manual review, never a release step.

## 4. Nothing unresolved from the handoff

`docs/superpowers/V3_DRAFT_HANDOFF.md` section 13 lists decisions the owner
must make, and `docs/superpowers/EXTERNAL_SETUP.md` lists setup only a human
can do (Formspree, reCAPTCHA, licensing, Zenodo). The Downloads page must not
ship a dead link or an invented DOI.

---

## First-time setup for the browser tests

Needed once per machine. The system Python is externally managed (PEP 668),
so this goes in a project venv and does not touch the system install.

```bash
python3 -m venv .venv-test
.venv-test/bin/python -m pip install -r GPCompaReports_v2/requirements-test.txt
.venv-test/bin/python -m playwright install chromium
```

The pinned version is `playwright==1.62.0`, which is what these tests were
written and verified against. `playwright install chromium` is a separate
step and is not optional: the pip package alone has no browser to drive. It
downloads a private Chromium (about 400 MB) into `~/.cache/ms-playwright`.

That browser is deliberately not the snap Chromium. Snap's `home` interface
silently denies access to hidden directories, which has cost this project
time more than once.

If you upgrade Playwright, re-run `playwright install chromium` for the new
version and update the pin in `GPCompaReports_v2/requirements-test.txt`.
