#!/usr/bin/env python3
"""Fetch the per-residue AlphaMissense LANDSCAPE from ProtVar for each GPCR.

Why this exists
---------------
The batch analysis pipeline (single_gpcr_processor.py) calls ProtVar's
`/score/{acc}/{pos}` endpoint, which returns AM scores for ALL 19 possible
substitutions at the queried position. But the batch keeps only the AM
score for the gnomAD-observed substitution and discards the rest. As a
result, the snake plot's AlphaMissense view colors a position by the
*observed*-variant landscape, not the *theoretical* landscape.

For positions where gnomAD has no pathogenic-scored variant but a
non-observed substitution would be pathogenic, the current view says
"ambiguous" or "benign" when the truer answer is "could be pathogenic".

This script fetches the full per-residue AM landscape and caches it next
to the conservation cache, so the snake plot can color by max AM score
across all 19 substitutions.

Behavior
--------
- Idempotent: skips a receptor if its cache already exists, unless --force.
- Resumable per-position: existing scores in a partial cache are preserved.
- Polite: 6 workers, exponential backoff on 429/5xx.

Cache format (mirrors fetch_conservation.py):
    {
      "gpcr_id":          "par2_human",
      "uniprot_id":       "P55085",
      "sequence_length":  397,
      "max_am_score":     {"159": 0.99, ...},    # max across all 19 mt
      "mean_am_score":    {"159": 0.43, ...},    # mean across all 19 mt
      "n_substitutions":  {"159": 19, ...},      # count of mt with a score
      "am_class":         {"159": "pathogenic"|"ambiguous"|"benign"|null}
                                                 # bucket of mean_am_score
    }
The mean is the default the snake plot uses: max-across-19 saturates
(every TM position has at least one pathogenic substitution), so it's not
informative. Mean spreads the signal so positions where MOST substitutions
break things stand out from positions where only 1-2 would.

Usage:
  python3 scripts/fetch_alphamissense.py par2_human
  python3 scripts/fetch_alphamissense.py                # all 283 GPCRs
  python3 scripts/fetch_alphamissense.py --limit 3
"""

import argparse
import concurrent.futures as cf
import json
import random
import sys
import time
from pathlib import Path

import requests


REPO_ROOT = Path(__file__).resolve().parent.parent
GPCOMP_V2_ROOT = REPO_ROOT / 'GPCompReports_v2'
BATCH_ROOT = REPO_ROOT / 'The_batch_RRCS_analyzer'

PROTVAR_BASE = 'https://www.ebi.ac.uk/ProtVar/api'
UNIPROT_BASE = 'https://rest.uniprot.org/uniprotkb'

DEFAULT_OUT_DIR = GPCOMP_V2_ROOT / 'output' / 'data'
DEFAULT_WORKERS = 6
REQUEST_TIMEOUT = 15
MAX_RETRIES = 5

# AlphaMissense classification cutoffs (the values AM ships with)
AM_PATHOGENIC_CUTOFF = 0.564
AM_BENIGN_CUTOFF = 0.34


def _bucket_class(score):
    if score is None:
        return None
    if score >= AM_PATHOGENIC_CUTOFF:
        return 'pathogenic'
    if score >= AM_BENIGN_CUTOFF:
        return 'ambiguous'
    return 'benign'


# ---------------------------------------------------------------------------
# GPCR discovery (mirrors fetch_conservation.py)
# ---------------------------------------------------------------------------

def _discover_gpcr_ids():
    batch_dir = next(BATCH_ROOT.glob('batch_analysis_full/batch_analysis_*'), None)
    if batch_dir is None:
        return []
    csv_dir = batch_dir / 'csv_data'
    return sorted(p.stem.replace('_rrcs_delta', '') for p in csv_dir.glob('*_rrcs_delta.csv'))


def _load_uniprot_map():
    batch_dir = next(BATCH_ROOT.glob('batch_analysis_full/batch_analysis_*'), None)
    if batch_dir is None:
        return {}
    summary = batch_dir / 'summary' / 'processing_results.csv'
    if not summary.exists():
        return {}
    import ast, csv
    result = {}
    with summary.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            gid = row.get('gpcr', '').strip()
            try:
                parsed = ast.literal_eval(row.get('summary', '{}'))
            except (ValueError, SyntaxError):
                parsed = {}
            accession = parsed.get('uniprot_id') or ''
            if gid and accession:
                result[gid] = accession.strip().upper()
    return result


def fetch_sequence_length(session, accession):
    url = f'{UNIPROT_BASE}/{accession}.json'
    r = session.get(url, timeout=REQUEST_TIMEOUT)
    r.raise_for_status()
    return int(r.json()['sequence']['length'])


# ---------------------------------------------------------------------------
# ProtVar AM extraction
# ---------------------------------------------------------------------------

def fetch_am_aggregates_for_position(session, accession, position):
    """Return (max_score, mean_score, n_substitutions) for one position.

    A ProtVar response payload is a list of `{name, mt, ...}` entries; the
    AM rows have name == 'AM' and the substitution residue is in `mt`. We
    aggregate amPathogenicity across all such rows.

    Returns (None, None, 0) on miss/error.
    """
    url = f'{PROTVAR_BASE}/score/{accession}/{position}'
    delay = 0.4
    for attempt in range(MAX_RETRIES):
        try:
            r = session.get(url, timeout=REQUEST_TIMEOUT)
        except (requests.ConnectionError, requests.Timeout):
            if attempt >= MAX_RETRIES - 1:
                return None, None, 0
            time.sleep(delay + random.random() * 0.1)
            delay *= 2
            continue
        if r.status_code == 200:
            try:
                payload = r.json()
            except ValueError:
                return None, None, 0
            if not isinstance(payload, list):
                return None, None, 0
            best = None
            total = 0.0
            n = 0
            for item in payload:
                if item.get('name') != 'AM':
                    continue
                v = item.get('amPathogenicity')
                if v is None:
                    continue
                try:
                    fv = float(v)
                except (TypeError, ValueError):
                    continue
                n += 1
                total += fv
                if best is None or fv > best:
                    best = fv
            mean = (total / n) if n else None
            return best, mean, n
        if r.status_code == 404:
            return None, None, 0
        if r.status_code in (429, 500, 502, 503, 504):
            time.sleep(delay + random.random() * 0.2)
            delay = min(delay * 2, 8.0)
            continue
        return None, None, 0
    return None, None, 0


def fetch_gpcr_alphamissense(gpcr_id, accession, out_dir, force=False, workers=DEFAULT_WORKERS):
    out_path = out_dir / f'alphamissense_{gpcr_id}.json'
    existing_max, existing_mean, existing_n, existing_class = {}, {}, {}, {}
    existing_length = None
    if out_path.exists() and not force:
        try:
            obj = json.loads(out_path.read_text())
            existing_max = obj.get('max_am_score', {}) or {}
            existing_mean = obj.get('mean_am_score', {}) or {}
            existing_n = obj.get('n_substitutions', {}) or {}
            existing_class = obj.get('am_class', {}) or {}
            existing_length = obj.get('sequence_length')
        except (ValueError, OSError):
            existing_max, existing_mean, existing_n, existing_class = {}, {}, {}, {}

    session = requests.Session()
    session.headers.update({'Accept': 'application/json'})

    try:
        length = existing_length or fetch_sequence_length(session, accession)
    except Exception as e:
        print(f'  [{gpcr_id}] failed to fetch sequence length for {accession}: {e}', flush=True)
        return None

    # If max is cached but mean isn't (older cache format), refetch everything
    # so we can populate the mean column. Otherwise resume per-position.
    if existing_max and not existing_mean:
        print(f'  [{gpcr_id}] cache lacks mean_am_score column; refetching all {length} positions', flush=True)
        existing_max, existing_mean, existing_n, existing_class = {}, {}, {}, {}

    missing = [p for p in range(1, length + 1) if str(p) not in existing_max]
    if not missing:
        print(f'  [{gpcr_id}] cached ({length} positions), skipping', flush=True)
        return out_path

    print(f'  [{gpcr_id}] {accession} length={length}, fetching {len(missing)} positions…', flush=True)
    new_max = dict(existing_max)
    new_mean = dict(existing_mean)
    new_n = dict(existing_n)
    new_class = dict(existing_class)

    def _one(pos):
        return (pos, *fetch_am_aggregates_for_position(session, accession, pos))

    t0 = time.time()
    done = 0
    checkpoint_every = 50
    with cf.ThreadPoolExecutor(max_workers=workers) as ex:
        for pos, mx, mn, n_subs in ex.map(_one, missing):
            key = str(pos)
            new_max[key] = mx
            new_mean[key] = mn
            new_n[key] = n_subs
            # Class is bucketed from the MEAN — that's what the snake plot uses
            new_class[key] = _bucket_class(mn)
            done += 1
            if done % checkpoint_every == 0:
                _write_cache(out_path, gpcr_id, accession, length,
                             new_max, new_mean, new_n, new_class)
    _write_cache(out_path, gpcr_id, accession, length,
                 new_max, new_mean, new_n, new_class)

    covered = sum(1 for v in new_mean.values() if v is not None)
    n_path = sum(1 for v in new_class.values() if v == 'pathogenic')
    n_amb = sum(1 for v in new_class.values() if v == 'ambiguous')
    n_ben = sum(1 for v in new_class.values() if v == 'benign')
    elapsed = time.time() - t0
    print(
        f'  [{gpcr_id}] saved {out_path.name}: '
        f'{covered}/{length} pos with AM (mean-based path/amb/ben = '
        f'{n_path}/{n_amb}/{n_ben}) ({elapsed:.1f}s)', flush=True
    )
    return out_path


def _write_cache(out_path, gpcr_id, accession, length, max_scores, mean_scores,
                 n_subs, classes):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps({
        'gpcr_id': gpcr_id,
        'uniprot_id': accession,
        'sequence_length': length,
        'max_am_score': max_scores,
        'mean_am_score': mean_scores,
        'n_substitutions': n_subs,
        'am_class': classes,
    }, separators=(',', ':')))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('gpcrs', nargs='*', help='GPCR ids to fetch (default: all)')
    ap.add_argument('--limit', type=int, default=None, help='Limit to first N from default list')
    ap.add_argument('--out', type=Path, default=DEFAULT_OUT_DIR, help='Output dir for alphamissense_*.json cache files')
    ap.add_argument('--workers', type=int, default=DEFAULT_WORKERS, help='Concurrent requests per GPCR')
    ap.add_argument('--force', action='store_true', help='Re-fetch even if cached')
    args = ap.parse_args()

    uniprot_map = _load_uniprot_map()
    if not uniprot_map:
        print('ERROR: could not load UniProt id map from batch summary', file=sys.stderr)
        return 2

    targets = list(args.gpcrs) if args.gpcrs else _discover_gpcr_ids()
    if args.limit:
        targets = targets[:args.limit]

    missing_acc = [g for g in targets if g not in uniprot_map]
    if missing_acc:
        print(f'WARNING: no UniProt id for: {missing_acc[:5]}{"..." if len(missing_acc) > 5 else ""}', file=sys.stderr)
        targets = [g for g in targets if g in uniprot_map]

    args.out.mkdir(parents=True, exist_ok=True)
    print(f'Cache dir: {args.out}')
    print(f'Targets: {len(targets)} GPCRs')

    t0 = time.time()
    for gid in targets:
        try:
            fetch_gpcr_alphamissense(gid, uniprot_map[gid], args.out, force=args.force, workers=args.workers)
        except KeyboardInterrupt:
            print('\nInterrupted. Partial cache retained; rerun to resume.', file=sys.stderr)
            return 130
        except Exception as e:
            print(f'  [{gid}] error: {e}', file=sys.stderr)
    print(f'Done in {time.time() - t0:.1f}s')
    return 0


if __name__ == '__main__':
    sys.exit(main())
