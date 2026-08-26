# External setup required before release

None of the items below are configured. They are gaps left deliberately by the
implementation, not oversights, because the account, the credential, or the
name did not exist to invent. This file is the checklist for a human on the
TGGR team to close them before the Contact and Downloads pages go live.

## Contact page

1. **Formspree account and endpoint.** Create a Formspree account for the lab
   The account's notification recipient MUST be `tggrlab@gmail.com`, or forward to it.
   That is the only address the Contact page shows visitors, so an account created under
   any other address would silently send submissions somewhere nobody is watching.

   (or the corresponding author) and add a new form for the GPCompaRe contact
   page. Replace `FORMSPREE_ENDPOINT_PLACEHOLDER` in
   `GPCompaReports_v2/website/page_generators/contact_page.py` with the real
   form endpoint URL (`https://formspree.io/f/<form-id>`).
2. **reCAPTCHA v3.** Obtain a site key and secret key from Google reCAPTCHA
   v3, configure them in the Formspree form settings (or add the widget to
   `contact.html` directly), and wire the secret into Formspree's spam
   filtering. No key of any kind is present anywhere in this codebase today.
3. **Domain restriction.** Once the site has a production domain, restrict
   the Formspree form (and the reCAPTCHA key) to that domain so the endpoint
   cannot be used from an unrelated page.
4. **Scientific and technical contacts.** `contact.html` currently states
   that named contacts are not yet assigned. Supply the corresponding
   author's name and role (scientific contact) and, if different, a
   technical contact for data or website issues, then fill in the marked
   gap in `GPCompaReports_v2/templates/contact.html` (the `.contact-gap`
   paragraph). Do not guess these from the manuscript author list without
   asking the people involved whether they want to be listed publicly with
   an email address.
5. **Re-enable the submit button.** `contact.html` disables the Send button
   whenever `formspree_endpoint` starts with `FORMSPREE_` (see
   `contact_page.py`). Once step 1 is done and the real endpoint is in
   place, the button re-enables itself automatically; no template change is
   needed for that part.

## Downloads page and data release

6. **Data license.** Confirm what license the RRCS dataset (CSV/matrix
   outputs) is released under before it is offered for download.
7. **Software license.** Confirm what license the analysis pipeline and
   site generator are released under, if the code itself is to be
   downloadable.
8. **Zenodo deposit and DOI.** Create a Zenodo deposit for the dataset
   release (and optionally the software), obtain a DOI, and record it on
   the landing page (`Version, downloads and citation` section currently
   says "DOI pending release") and on the Downloads page once it exists.
9. **The two Downloads placeholders are not equivalent.** `downloads.html`
   (`GPCompaReports_v2/website/page_generators/downloads_page.py`) currently
   marks both the database archive and the software archive "Not yet
   available", with no `.zip`/`.tar.gz` link and no Zenodo record. Only the
   software half may ship that way at publication: there is no releasable
   user-facing analysis program yet, and that placeholder may stand until
   there is, possibly well after publication. The database half must not:
   a paper companion that offers no data is not a companion, so before the
   final public release, item 6 (data license) and item 8 (Zenodo DOI) above
   must be closed and the "Not yet available" database block replaced with a
   real archive link. Do not let the database placeholder ship as final
   because the software placeholder was allowed to.

## What is safe to ship without these

The Contact page can be built and deployed with the Formspree endpoint left
unconfigured: the submit button is disabled and a visible notice tells
visitors the form is not connected, so no message can be silently lost.


> **Superseded 2026-08-26.** The form is now live: it posts as plain HTML to `https://formspree.io/f/mljrqazl`, the submit button and all visible fields ship enabled, and the "not connected" notice no longer renders. What remains is dashboard-only and listed under "Contact form: manual Formspree dashboard checks" in `RELEASE_CHECKLIST.md`.
Shipping the page in this state is a documentation gap, not a functional
bug. Publishing it with real Formspree/reCAPTCHA credentials and named
contacts is a separate release step gated on the items above.

The Downloads page can likewise be built and deployed with both archives
marked "Not yet available": it invents no file, link, or DOI, so nothing is
promised that does not exist. That is a safe *draft* state, not a safe
*publication* state for the database half; see item 9.
