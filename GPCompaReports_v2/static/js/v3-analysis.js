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

    // The 16 named segments do not necessarily account for every above-
    // threshold contact endpoint: an endpoint whose residue has no resolved
    // segment (no GPCRdb generic number) counts toward the receptor's total
    // but draws no bar. That share can be large (ADRB2 is about 35%), and
    // without a note the bars read as a complete picture when they are not.
    // Surface it as a small factual addition to the header, not an
    // interpretation, only when it is more than rounding noise.
    var accounted = 0;
    segments.forEach(function (s) { accounted += profile[s] || 0; });
    var unaccounted = Math.round(Math.max(0, 100 - accounted) * 10) / 10;

    var box = document.createElement('div');
    box.className = 'v3-fp';
    box.id = 'v3-fp';
    var head = document.createElement('div');
    head.className = 'v3-fp-head';
    var headText = 'Where above-threshold contacts sit, against the median of the database';
    if (unaccounted >= 1) {
      headText += ' (' + unaccounted.toFixed(1) + '% of endpoints are unassigned, not shown)';
    }
    var headLabel = document.createElement('span');
    headLabel.textContent = headText;
    if (unaccounted >= 1) {
      headLabel.title = unaccounted.toFixed(1) + '% of this receptor\'s above-threshold ' +
        'contact endpoints have no resolved segment (no GPCRdb generic number), so they are ' +
        'counted in the totals but do not appear as a bar.';
    }
    var legend = document.createElement('span');
    legend.className = 'v3-legend';
    legend.innerHTML = '<i class="bar"></i>this receptor <i class="med"></i>median';
    head.appendChild(headLabel);
    head.appendChild(legend);
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
