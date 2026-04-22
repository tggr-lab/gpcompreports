#!/usr/bin/env python3
"""Fetch per-residue conservation scores from ProtVar for each GPCR and cache
them to disk.

The batch analysis pipeline only fetches ProtVar scores for variant or
significant-RRCS positions, so most residues have no conservation data in
`*_variants.csv`. The v2 snake plot wants a conservation grade for every
residue in the receptor.

This script fills that gap:
  1. For each requested GPCR, look up its UniProt accession from
     `class_A_all.csv` (or the batch processing summary).
  2. Query UniProt for the canonical sequence length.
  3. Call `GET https://www.ebi.ac.uk/ProtVar/api/score/{acc}/{pos}` for every
     residue 1..length, extract the CONSERV entry's score.
  4. Write `GPCompReports/output/data/conservation_{gpcr_id}.json`:
       {"gpcr_id": "5ht1a_human", "uniprot_id": "P08908",
        "sequence_length": 422, "scores": {"1": 0.42, "2": 0.88, ...}}

Behavior:
- Idempotent: if the cache file already exists, skipped unless --force.
- Resumable at per-position level: if a partial cache is present, we keep
  its existing scores and fill the missing positions only.
- Polite: thread pool of 6 workers, exponential backoff on 429/5xx, short
  timeouts, tolerates missing CONSERV entries (stored as null).

Usage:
  python3 scripts/fetch_conservation.py                      # all 283 GPCRs
  python3 scripts/fetch_conservation.py 5ht1a_human 5ht1b_human
  python3 scripts/fetch_conservation.py --limit 3
  python3 scripts/fetch_conservation.py --force 5ht1a_human  # re-fetch
"""

import argparse
import concurrent.futures as cf
import json
import random
import re
import sys
import time
from pathlib import Path

import requests


REPO_ROOT = Path(__file__).resolve().parent.parent
GPCOMP_V2_ROOT = REPO_ROOT / 'GPCompReports_v2'
BATCH_ROOT = REPO_ROOT / 'The_batch_RRCS_analyzer'
METADATA_CSV = REPO_ROOT / 'class_A_all.csv'

PROTVAR_BASE = 'https://www.ebi.ac.uk/ProtVar/api'
UNIPROT_BASE = 'https://rest.uniprot.org/uniprotkb'

# The v2 report generator reads this cache directory.
DEFAULT_OUT_DIR = GPCOMP_V2_ROOT / 'output' / 'data'
DEFAULT_WORKERS = 6
REQUEST_TIMEOUT = 15
MAX_RETRIES = 5


# ---------------------------------------------------------------------------
# GPCR discovery
# ---------------------------------------------------------------------------

def _discover_gpcr_ids():
    """All gpcr_ids the batch pipeline produced (one per delta CSV)."""
    batch_dir = next(BATCH_ROOT.glob('batch_analysis_full/batch_analysis_*'), None)
    if batch_dir is None:
        return []
    csv_dir = batch_dir / 'csv_data'
    return sorted(p.stem.replace('_rrcs_delta', '') for p in csv_dir.glob('*_rrcs_delta.csv'))


def _load_uniprot_map():
    """Map filesystem gpcr_id -> UniProt accession (e.g. '5ht1a_human' -> 'P08908').

    Read from the processing summary emitted by the batch pipeline.
    """
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
            parsed = {}
            try:
                parsed = ast.literal_eval(row.get('summary', '{}'))
            except (ValueError, SyntaxError):
                pass
            accession = parsed.get('uniprot_id') or ''
            if gid and accession:
                result[gid] = accession.strip().upper()
    return result


# ---------------------------------------------------------------------------
# UniProt sequence length
# ---------------------------------------------------------------------------

def fetch_sequence_length(session, accession):
    url = f'{UNIPROT_BASE}/{accession}.json'
    r = session.get(url, timeout=REQUEST_TIMEOUT)
    r.raise_for_status()
    return int(r.json()['sequence']['length'])


# ---------------------------------------------------------------------------
# ProtVar per-residue CONSERV
# ---------------------------------------------------------------------------

def fetch_conservation_for_position(session, accession, position):
    """Return the CONSERV score (float) for one residue, or None on miss.

    Retries on transient failures (429, 5xx) with exponential backoff.
    """
    url = f'{PROTVAR_BASE}/score/{accession}/{position}'
    delay = 0.4
    for attempt in range(MAX_RETRIES):
        try:
            r = session.get(url, timeout=REQUEST_TIMEOUT)
        except (requests.ConnectionError, requests.Timeout):
            if attempt >= MAX_RETRIES - 1:
                return None
            time.sleep(delay + random.random() * 0.1)
            delay *= 2
            continue
        if r.status_code == 200:
            try:
                payload = r.json()
            except ValueError:
                return None
            if not isinstance(payload, list):
                return None
            for item in payload:
                if item.get('name') == 'CONSERV':
                    score = item.get('score')
                    if score is None:
                        return None
                    try:
                        return float(score)
                    except (TypeError, ValueError):
                        return None
            return None  # no CONSERV entry for this position
        if r.status_code == 404:
            return None
        if r.status_code in (429, 500, 502, 503, 504):
            time.sleep(delay + random.random() * 0.2)
            delay = min(delay * 2, 8.0)
            continue
        # Other 4xx: skip silently
        return None
    return None


def fetch_gpcr_conservation(gpcr_id, accession, out_dir, force=False, workers=DEFAULT_WORKERS):
    """Fetch (or resume fetching) conservation scores for one GPCR."""
    out_path = out_dir / f'conservation_{gpcr_id}.json'
    existing = {}
    existing_length = None
    if out_path.exists() and not force:
        try:
            obj = json.loads(out_path.read_text())
            existing = obj.get('scores', {}) or {}
            existing_length = obj.get('sequence_length')
        except (ValueError, OSError):
            existing = {}

    session = requests.Session()
    session.headers.update({'Accept': 'application/json'})

    try:
        length = existing_length or fetch_sequence_length(session, accession)
    except Exception as e:
        print(f'  [{gpcr_id}] failed to fetch sequence length for {accession}: {e}', flush=True)
        return None

    # Only fetch positions we don't already have
    missing = [p for p in range(1, length + 1) if str(p) not in existing]
    if not missing:
        print(f'  [{gpcr_id}] cached ({length} positions), skipping', flush=True)
        return out_path

    print(f'  [{gpcr_id}] {accession} length={length}, fetching {len(missing)} positions…', flush=True)
    new_scores = dict(existing)

    def _one(pos):
        return pos, fetch_conservation_for_position(session, accession, pos)

    t0 = time.time()
    done = 0
    checkpoint_every = 50
    with cf.ThreadPoolExecutor(max_workers=workers) as ex:
        for pos, score in ex.map(_one, missing):
            new_scores[str(pos)] = score
            done += 1
            if done % checkpoint_every == 0:
                _write_cache(out_path, gpcr_id, accession, length, new_scores)
    _write_cache(out_path, gpcr_id, accession, length, new_scores)

    covered = sum(1 for v in new_scores.values() if v is not None)
    elapsed = time.time() - t0
    print(
        f'  [{gpcr_id}] saved {out_path.name}: '
        f'{covered}/{length} positions with scores ({elapsed:.1f}s)', flush=True
    )
    return out_path


def _write_cache(out_path, gpcr_id, accession, length, scores):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps({
        'gpcr_id': gpcr_id,
        'uniprot_id': accession,
        'sequence_length': length,
        'scores': scores,
    }, separators=(',', ':')))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('gpcrs', nargs='*', help='GPCR ids to fetch (default: all)')
    ap.add_argument('--limit', type=int, default=None, help='Limit to first N from default list')
    ap.add_argument('--out', type=Path, default=DEFAULT_OUT_DIR, help='Output dir for conservation_*.json cache files')
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
            fetch_gpcr_conservation(gid, uniprot_map[gid], args.out, force=args.force, workers=args.workers)
        except KeyboardInterrupt:
            print('\nInterrupted. Partial cache retained; rerun to resume.', file=sys.stderr)
            return 130
        except Exception as e:
            print(f'  [{gid}] error: {e}', file=sys.stderr)
    print(f'Done in {time.time() - t0:.1f}s')
    return 0


if __name__ == '__main__':
    sys.exit(main())
