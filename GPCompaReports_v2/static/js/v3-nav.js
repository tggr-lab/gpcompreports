/* v3-nav.js — sticky section nav for report pages.
   Builds itself from [data-section-title] so adding a section needs no JS change.
   Navigation is wayfinding, not content: this is always on. */
(function () {
  'use strict';
  var V3Nav = (window.V3Nav = window.V3Nav || {});

  function build() {
    buildCompactToggle();
    var sections = document.querySelectorAll('.report-section[data-section-title]');
    if (sections.length < 2) return;

    var nav = document.createElement('nav');
    nav.className = 'v3-secnav';
    nav.setAttribute('aria-label', 'Report sections');
    var list = document.createElement('div');
    list.className = 'v3-secnav-list';
    nav.appendChild(list);

    var links = [];
    Array.prototype.forEach.call(sections, function (sec) {
      var a = document.createElement('a');
      a.href = '#' + sec.id;
      a.textContent = sec.dataset.sectionTitle;
      a.dataset.target = sec.id;
      list.appendChild(a);
      links.push(a);
    });

    var host = document.querySelector('.report-page .section-inner');
    var anchor = document.querySelector('.report-stats');
    if (!host || !anchor) return;
    anchor.parentNode.insertBefore(nav, anchor.nextSibling);

    if (!('IntersectionObserver' in window)) return;
    var seen = {};
    var obs = new IntersectionObserver(function (entries) {
      entries.forEach(function (e) { seen[e.target.id] = e.isIntersecting; });
      var activeId = null;
      Array.prototype.forEach.call(sections, function (sec) {
        if (!activeId && seen[sec.id]) activeId = sec.id;
      });
      links.forEach(function (a) {
        a.classList.toggle('on', a.dataset.target === activeId);
      });
    }, { rootMargin: '-80px 0px -70% 0px' });
    Array.prototype.forEach.call(sections, function (s) { obs.observe(s); });
  }

  var COMPACT_KEY = 'gpcompare-compact';

  function readFlag(key) {
    try { return localStorage.getItem(key) === '1'; } catch (e) { return false; }
  }
  function writeFlag(key, on) {
    try { on ? localStorage.setItem(key, '1') : localStorage.removeItem(key); } catch (e) {}
  }
  V3Nav.readFlag = readFlag;
  V3Nav.writeFlag = writeFlag;

  function ensureToolbar() {
    var existing = document.querySelector('.v3-toolbar');
    if (existing) return existing;
    var anchor = document.querySelector('.report-stats');
    if (!anchor) return null;
    var bar = document.createElement('div');
    bar.className = 'v3-toolbar';
    anchor.parentNode.insertBefore(bar, anchor.nextSibling);
    return bar;
  }
  V3Nav.ensureToolbar = ensureToolbar;

  function applyCompact(on) {
    document.querySelectorAll('.report-section[data-section-title]').forEach(function (sec) {
      sec.classList.toggle('v3-folded', on);
    });
  }

  function buildCompactToggle() {
    var bar = ensureToolbar();
    if (!bar) return;
    var btn = document.createElement('button');
    btn.type = 'button';
    btn.className = 'btn btn-ghost btn-sm v3-toggle';
    btn.id = 'v3-compact-toggle';
    var on = readFlag(COMPACT_KEY);
    function paint() {
      btn.setAttribute('aria-pressed', on ? 'true' : 'false');
      btn.classList.toggle('on', on);
      btn.textContent = on ? 'Compact view: on' : 'Compact view';
    }
    btn.addEventListener('click', function () {
      on = !on;
      writeFlag(COMPACT_KEY, on);
      applyCompact(on);
      paint();
    });
    bar.appendChild(btn);
    paint();
    applyCompact(on);
  }

  document.addEventListener('click', function (ev) {
    var head = ev.target.closest ? ev.target.closest('.report-section.v3-folded .report-section-head') : null;
    if (!head) return;
    head.parentNode.classList.remove('v3-folded');
  });

  V3Nav.build = build;
  document.addEventListener('DOMContentLoaded', build);
})();
