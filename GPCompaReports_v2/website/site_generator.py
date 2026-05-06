"""Site generator for GPCompReports v2 — Primer + teal/orange design system.

Self-contained v2 pipeline. The v1 site lives in ../GPCompReports/ and is
unchanged by this build; the two projects share only the upstream RRCS batch
data (../The_batch_RRCS_analyzer/batch_analysis_full/).
"""

import shutil
from pathlib import Path

from jinja2 import Environment, FileSystemLoader

from ..analysis.data_loader import GPCRDataStore
from ..analysis.cross_gpcr_analysis import run_cross_gpcr_analysis
from ..analysis.tm_domain_analysis import run_tm_domain_analysis
from ..analysis.cfr_analysis import run_cfr_analysis
from ..analysis.variant_correlation import run_variant_analysis
from .page_generators.landing_page import generate_landing_page
from .page_generators.gpcr_index import generate_gpcr_index
from .page_generators.gpcr_report_page import generate_all_reports
from .page_generators.statistics_page import generate_statistics_page


class SiteGenerator:
    """Orchestrates the v2 pipeline: load -> analyze -> generate pages."""

    def __init__(self, batch_dir=None, metadata_csv=None, output_dir=None, limit=None):
        self.base_dir = Path(__file__).resolve().parent.parent
        self.output_dir = Path(output_dir) if output_dir else self.base_dir / 'output'
        self.template_dir = self.base_dir / 'templates'
        self.static_dir = self.base_dir / 'static'
        self.limit = limit

        self.store = GPCRDataStore(batch_dir=batch_dir, metadata_csv=metadata_csv)
        self.analysis_results = {}

    def run(self):
        print("=" * 60)
        print("GPCompReports v2 — Site Generator")
        print("=" * 60)

        print("\n[1/5] Loading data...")
        self.store.load_all()

        print("\n[2/5] Running cross-GPCR analysis...")
        self.analysis_results['cross_gpcr'] = run_cross_gpcr_analysis(self.store)

        print("\n[3/5] Running TM domain, CFR, and variant analysis...")
        self.analysis_results['tm_domain'] = run_tm_domain_analysis(self.store)
        self.analysis_results['cfr'] = run_cfr_analysis(self.store)
        cfr_table = self.analysis_results['cfr'].get('cfr_table')
        self.analysis_results['variant'] = run_variant_analysis(self.store, cfr_table)

        print("\n[4/5] Preparing output directory...")
        self._prepare_output()

        print("\n[5/5] Generating website pages...")
        env = Environment(loader=FileSystemLoader(str(self.template_dir)))

        print("  Landing page...")
        generate_landing_page(env, self.store, self.output_dir)

        print("  GPCR browser...")
        generate_gpcr_index(env, self.store, self.output_dir)

        print("  Statistics page...")
        generate_statistics_page(env, self.store, self.analysis_results, self.output_dir)

        if self.limit:
            print(f"  Individual reports (limit: {self.limit} pages)...")
        else:
            print(f"  Individual reports ({len(self.store.gpcr_ids)} pages)...")
        generate_all_reports(env, self.store, self.output_dir,
                             analysis_results=self.analysis_results, limit=self.limit)

        print("\n" + "=" * 60)
        print("Site generation complete!")
        print(f"Open: {self.output_dir / 'index.html'}")
        print("=" * 60)

    def _prepare_output(self):
        for subdir in ['reports', 'browse', 'data', 'static/css', 'static/js', 'static/img']:
            (self.output_dir / subdir).mkdir(parents=True, exist_ok=True)
        for src in self.static_dir.rglob('*'):
            if src.is_file():
                rel = src.relative_to(self.static_dir)
                dst = self.output_dir / 'static' / rel
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, dst)
        print(f"  Output directory: {self.output_dir}")
