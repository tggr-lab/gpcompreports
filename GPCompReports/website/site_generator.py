"""Site generator orchestrator: data -> analysis -> pages."""

import shutil
from pathlib import Path

from jinja2 import Environment, FileSystemLoader

from ..analysis.data_loader import GPCRDataStore
from ..analysis.cross_gpcr_analysis import run_cross_gpcr_analysis
from ..analysis.tm_domain_analysis import run_tm_domain_analysis
from ..analysis.cfr_analysis import run_cfr_analysis
from ..analysis.variant_correlation import run_variant_analysis
from .page_generators.landing_page import generate_landing_page
from .page_generators.landing_page_v2 import generate_landing_page_v2
from .page_generators.gpcr_index import generate_gpcr_index
from .page_generators.gpcr_report_page import generate_all_reports
from .page_generators.statistics_page import generate_statistics_page


class SiteGenerator:
    """Orchestrates the complete site generation pipeline."""

    def __init__(self, batch_dir=None, metadata_csv=None, output_dir=None, limit=None):
        self.base_dir = Path(__file__).resolve().parent.parent
        self.output_dir = Path(output_dir) if output_dir else self.base_dir / 'output'
        self.template_dir = self.base_dir / 'templates'
        self.static_dir = self.base_dir / 'static'
        self.limit = limit

        self.store = GPCRDataStore(batch_dir=batch_dir, metadata_csv=metadata_csv)
        self.analysis_results = {}

    def run(self):
        """Execute the full pipeline: load -> analyze -> generate."""
        print("=" * 60)
        print("GPCompReports Site Generator")
        print("=" * 60)

        # Step 1: Load data
        print("\n[1/5] Loading data...")
        self.store.load_all()

        # Step 2: Run analysis
        print("\n[2/5] Running cross-GPCR analysis...")
        self.analysis_results['cross_gpcr'] = run_cross_gpcr_analysis(self.store)

        print("\n[3/5] Running TM domain, CFR, and variant analysis...")
        self.analysis_results['tm_domain'] = run_tm_domain_analysis(self.store)
        self.analysis_results['cfr'] = run_cfr_analysis(self.store)

        cfr_table = self.analysis_results['cfr'].get('cfr_table')
        self.analysis_results['variant'] = run_variant_analysis(self.store, cfr_table)

        # Step 3: Prepare output directory
        print("\n[4/5] Preparing output directory...")
        self._prepare_output()

        # Step 4: Generate pages
        print("\n[5/5] Generating website pages...")
        env = Environment(loader=FileSystemLoader(str(self.template_dir)))

        print("  Landing page...")
        generate_landing_page(env, self.store, self.output_dir)
        generate_landing_page_v2(env, self.store, self.output_dir)

        print("  GPCR browser...")
        generate_gpcr_index(env, self.store, self.output_dir)

        print("  Statistics page...")
        generate_statistics_page(env, self.store, self.analysis_results, self.output_dir)

        if self.limit:
            print(f"  Individual reports (limit: {self.limit} pages)...")
        else:
            print("  Individual reports (283 pages)...")
        generate_all_reports(env, self.store, self.output_dir, limit=self.limit)

        print("\n" + "=" * 60)
        print("Site generation complete!")
        print(f"Open: {self.output_dir / 'index.html'}")
        print("=" * 60)

    def _prepare_output(self):
        """Create output directory structure and copy static assets."""
        # Create directories
        for subdir in ['reports', 'browse', 'data', 'static/css', 'static/js', 'static/img']:
            (self.output_dir / subdir).mkdir(parents=True, exist_ok=True)

        # Copy static assets
        for src in self.static_dir.rglob('*'):
            if src.is_file():
                rel = src.relative_to(self.static_dir)
                dst = self.output_dir / 'static' / rel
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, dst)

        print(f"  Output directory: {self.output_dir}")
