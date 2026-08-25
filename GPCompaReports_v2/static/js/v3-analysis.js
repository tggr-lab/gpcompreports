/* v3-analysis.js — opt-in analysis layer.
   Off by default. When on it ADDS a key-number strip, a segment fingerprint,
   and row badges. It never removes or rewords anything already on the page. */
(function () {
  'use strict';
  var KEY = 'gpcompare-analysis';
  var data = null;
  // 30, matching the site's own published Top 30 (Statistics page and the
  // report's Core Functional snake-plot view both use 30); a badge citing a
  // rank past this would point a reader to a residue neither view shows.
  var CFR_RANK_CUTOFF = 30;

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
    }).join(' · ') || 'no above-threshold contacts with a resolved segment';
    // "resolved to a named segment" not "annotated": an endpoint can be
    // annotated with a segment outside the 16 the profile knows, and those
    // are excluded from the base. Keep this wording in step with the
    // fingerprint header below.
    wrap.appendChild(tile('Largest segment shares', segs,
      (k.top_segments || []).length ?
        'shares of endpoints that resolved to a named segment, not of all above-threshold contacts' :
        ''));

    if (k.largest_structured) {
      var L = k.largest_structured;
      var val = '<span class="v3-gn">' + (L.num1 || L.label1) + '</span> ↔ ' +
                '<span class="v3-gn">' + (L.num2 || L.label2) + '</span>';
      wrap.appendChild(tile('Largest change, structured region', val,
        L.label1 + ' (' + L.seg1 + ') / ' + L.label2 + ' (' + L.seg2 + ')  |Δ| ' + L.abs_delta));
    } else {
      wrap.appendChild(tile('Largest change, structured region', 'none',
        'every above-threshold contact has an endpoint in a terminus, a loop, or without a resolved segment'));
    }

    wrap.appendChild(tile('Top movers that are CFRs',
      k.cfr_top_movers + ' of ' + k.top_mover_count,
      'cross-receptor Core Functional Residues'));
    return wrap;
  }

  function buildFingerprint(profile, median, segments, coverage, medianCoverage) {
    var W = 640, H = 116, padB = 20, plot = H - padB;
    // A median profile that is absent or empty (e.g. a preview build whose
    // analysis dict never computed one) is not the same thing as a database
    // median measured at zero in every segment. Treat it as absent and omit
    // the comparator entirely, rather than drawing it flat along the axis.
    var hasMedian = !!median && Object.keys(median).length > 0;
    var max = 0;
    segments.forEach(function (s) {
      max = Math.max(max, profile[s] || 0, hasMedian ? (median[s] || 0) : 0);
    });
    max = Math.max(max, 1) * 1.15;
    var bw = W / segments.length;
    var NS = 'http://www.w3.org/2000/svg';
    var svg = document.createElementNS(NS, 'svg');
    svg.setAttribute('viewBox', '0 0 ' + W + ' ' + H);
    svg.setAttribute('role', 'img');
    svg.setAttribute('aria-label', hasMedian ?
      'Above-threshold contact distribution by segment, this receptor against the per-segment median across receptors' :
      'Above-threshold contact distribution by segment, this receptor');

    segments.forEach(function (s, i) {
      var x = i * bw, v = profile[s] || 0, m = hasMedian ? (median[s] || 0) : 0;
      var h = v / max * plot, mh = m / max * plot;
      var structured = /^TM[1-7]$/.test(s) || s === 'H8';
      var bar = document.createElementNS(NS, 'rect');
      bar.setAttribute('x', x + bw * 0.18); bar.setAttribute('y', plot - h);
      bar.setAttribute('width', bw * 0.5); bar.setAttribute('height', Math.max(h, 0.8));
      bar.setAttribute('rx', '1.5');
      bar.setAttribute('fill', structured ? 'var(--brand-primary)' : 'var(--fgColor-muted)');
      bar.setAttribute('opacity', structured ? '1' : '0.45');
      var t = document.createElementNS(NS, 'title');
      t.textContent = hasMedian ?
        (s + ': ' + v.toFixed(1) + '% here, ' + m.toFixed(1) + '% median') :
        (s + ': ' + v.toFixed(1) + '% here');
      bar.appendChild(t);
      svg.appendChild(bar);

      if (hasMedian) {
        var line = document.createElementNS(NS, 'line');
        line.setAttribute('x1', x + bw * 0.10); line.setAttribute('x2', x + bw * 0.78);
        line.setAttribute('y1', plot - mh); line.setAttribute('y2', plot - mh);
        line.setAttribute('stroke', 'var(--fgColor-muted)');
        line.setAttribute('stroke-width', '1.4');
        line.setAttribute('stroke-dasharray', '3 2');
        svg.appendChild(line);
      }

      var lab = document.createElementNS(NS, 'text');
      lab.setAttribute('x', x + bw * 0.43); lab.setAttribute('y', H - 6);
      lab.setAttribute('text-anchor', 'middle'); lab.setAttribute('class', 'v3-seg-lab');
      lab.textContent = s;
      svg.appendChild(lab);
    });

    // `profile` and `median` are both shares of endpoints that resolved to a
    // named segment, so the two are on a common, comparable base (see
    // receptor_profile.segment_profile). `coverage` and `medianCoverage`,
    // computed server side, say how much of each side's above-threshold
    // contacts that base actually covers. Disclose both, factually, only
    // when either falls short of full coverage. The comparator is a
    // per-segment median across receptors, not a receptor itself, so it has
    // no coverage of its own: medianCoverage is worded as a median across
    // receptors, never as a property of the dashed line.
    coverage = coverage || 0;
    medianCoverage = medianCoverage || 0;

    var box = document.createElement('div');
    box.className = 'v3-fp';
    box.id = 'v3-fp';
    var head = document.createElement('div');
    head.className = 'v3-fp-head';
    var headText = hasMedian ?
      'Where above-threshold contacts sit, against the per-segment median across receptors (not itself a profile)' :
      'Where above-threshold contacts sit';
    var partial = hasMedian ? (coverage < 100 || medianCoverage < 100) : coverage < 100;
    if (partial) {
      headText += hasMedian ?
        (' (this receptor: ' + coverage.toFixed(1) + '% of endpoints resolved to a named segment, ' +
          'median across receptors: ' + medianCoverage.toFixed(1) + '%)') :
        (' (this receptor: ' + coverage.toFixed(1) + '% of endpoints resolved to a named segment)');
    }
    var headLabel = document.createElement('span');
    headLabel.textContent = headText;
    if (partial) {
      headLabel.title = hasMedian ?
        ('Bars are shares of above-threshold contact endpoints that resolved to a named segment, not of all ' +
          'endpoints. Coverage is the percentage that resolved: ' + coverage.toFixed(1) + '% for this receptor. ' +
          'Median across receptors: ' + medianCoverage.toFixed(1) + '%.') :
        ('Bars are shares of above-threshold contact endpoints that resolved to a named segment, not of all ' +
          'endpoints. Coverage is the percentage that resolved: ' + coverage.toFixed(1) + '% for this receptor.');
    }
    var legend = document.createElement('span');
    legend.className = 'v3-legend';
    legend.innerHTML = hasMedian ?
      '<i class="bar"></i>this receptor <i class="med"></i>median' :
      '<i class="bar"></i>this receptor';
    head.appendChild(headLabel);
    head.appendChild(legend);
    box.appendChild(head);
    box.appendChild(svg);
    return box;
  }

  function badge(text, cls) {
    var s = document.createElement('span');
    s.className = 'v3-badge ' + cls;
    // The visible label lives in CSS (`::after { content: attr(data-label) }`
    // in v3.css), never in textContent: site.js reads textContent for CSV
    // export, search, and column sort, and this label must not contaminate
    // any of those. aria-label keeps the label available to assistive tech.
    s.setAttribute('data-label', text);
    s.setAttribute('aria-label', text);
    // aria-label is ignored on an element whose implicit role is `generic`
    // (a bare <span>); role="img" gives it an accessible-name-bearing role
    // so screen readers announce the label.
    s.setAttribute('role', 'img');
    return s;
  }

  function decorateChangesTable() {
    var table = document.getElementById('top-changes-table');
    if (!table) return;
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
      if (best !== null && best <= CFR_RANK_CUTOFF) {
        holder.appendChild(badge('CFR #' + best, 'v3-b-cfr'));
      }
      // Neither endpoint has a GPCRdb generic number in this receptor's
      // annotation. That alone does not tell us where the residue sits: the
      // annotation file simply may not cover it. State only what is known.
      if (!nums.length) {
        var b = badge('no GPCRdb number', 'v3-b-low');
        b.title = 'Neither endpoint has a GPCRdb generic number, so its ' +
                  'position cannot be placed in the reference numbering.';
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
    var fp = buildFingerprint(data.profile, data.median_profile, data.segments,
      data.coverage, data.median_coverage);
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
