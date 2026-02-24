"""Generate the landing page."""

from jinja2 import Environment


def generate_landing_page(env: Environment, store, output_dir):
    """Generate index.html landing page."""
    info_df = store.get_all_info_df()
    top5 = info_df.nlargest(5, 'sum_abs_delta').to_dict('records')

    total_contacts = int(info_df['total_contacts'].sum())
    n_families = info_df['receptor_family'].nunique()
    n_ligand_types = info_df['ligand_type'].nunique()

    template = env.get_template('landing.html')
    html = template.render(
        static_prefix='',
        active_page='home',
        total_gpcrs=len(store.gpcr_ids),
        n_families=n_families,
        n_ligand_types=n_ligand_types,
        total_contacts_k=f"{total_contacts // 1000}",
        top5=top5,
    )

    out_path = output_dir / 'index.html'
    out_path.write_text(html, encoding='utf-8')
    print(f"  Generated: {out_path.name}")
